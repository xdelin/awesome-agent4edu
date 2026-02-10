"""XMILE XML generation and parsing for Stella .stmx files."""

import math
import re
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from typing import Optional
import uuid
from html import escape


# XML namespaces
XMILE_NS = "http://docs.oasis-open.org/xmile/ns/XMILE/v1.0"
ISEE_NS = "http://iseesystems.com/XMILE"

# Stella/XMILE built-in functions (not variable names)
STELLA_FUNCTIONS = {
    'IF', 'THEN', 'ELSE', 'AND', 'OR', 'NOT',
    'MIN', 'MAX', 'ABS', 'SIN', 'COS', 'TAN',
    'EXP', 'LN', 'LOG', 'LOG10', 'SQRT', 'INT',
    'ROUND', 'MOD', 'TIME', 'DT', 'STARTTIME', 'STOPTIME',
    'DELAY', 'DELAY1', 'DELAY3', 'DELAYN',
    'SMOOTH', 'SMOOTH3', 'SMOOTHN', 'SMTH1', 'SMTH3', 'SMTHN',
    'TREND', 'FORCST', 'PULSE', 'STEP', 'RAMP',
    'RANDOM', 'NORMAL', 'POISSON', 'EXPRND',
    'PREVIOUS', 'INIT', 'SELF', 'SUM', 'MEAN',
    'GRAPH', 'LOOKUP', 'INTERPOLATE', 'HISTORY',
    'SAFEDIV', 'NPV', 'IRR', 'COUNTER',
    'TRUE', 'FALSE', 'PI', 'E', 'INF', 'NAN',
}

# Layout constants
AUX_RADIUS = 18  # Default aux circle radius in pixels


# =============================================================================
# Geometry Primitives for Layout Collision/Crossing Detection
# =============================================================================

@dataclass
class BoundingBox:
    """Axis-aligned bounding box for collision detection."""
    x: float  # center x
    y: float  # center y
    width: float
    height: float

    def intersects(self, other: "BoundingBox", margin: float = 0) -> bool:
        """Check if two boxes overlap (with optional margin)."""
        return (abs(self.x - other.x) < (self.width + other.width) / 2 + margin and
                abs(self.y - other.y) < (self.height + other.height) / 2 + margin)


def _ccw(a: tuple[float, float], b: tuple[float, float], c: tuple[float, float]) -> bool:
    """Check if three points are in counter-clockwise order."""
    return (c[1] - a[1]) * (b[0] - a[0]) > (b[1] - a[1]) * (c[0] - a[0])


def segments_intersect(
    p1: tuple[float, float], p2: tuple[float, float],
    p3: tuple[float, float], p4: tuple[float, float]
) -> bool:
    """Check if line segment p1-p2 intersects segment p3-p4.

    Uses the counter-clockwise orientation test for robust intersection detection.
    """
    return (_ccw(p1, p3, p4) != _ccw(p2, p3, p4) and
            _ccw(p1, p2, p3) != _ccw(p1, p2, p4))


def segment_intersects_box(
    p1: tuple[float, float], p2: tuple[float, float],
    box: BoundingBox
) -> bool:
    """Check if line segment p1-p2 intersects a bounding box.

    Tests intersection with all four edges of the box.
    """
    # Box corners
    left = box.x - box.width / 2
    right = box.x + box.width / 2
    top = box.y - box.height / 2
    bottom = box.y + box.height / 2

    # Check if segment is entirely inside box
    if (left <= p1[0] <= right and top <= p1[1] <= bottom and
        left <= p2[0] <= right and top <= p2[1] <= bottom):
        return True

    # Check intersection with each edge
    edges = [
        ((left, top), (right, top)),      # top edge
        ((right, top), (right, bottom)),  # right edge
        ((right, bottom), (left, bottom)),  # bottom edge
        ((left, bottom), (left, top)),    # left edge
    ]

    for edge_p1, edge_p2 in edges:
        if segments_intersect(p1, p2, edge_p1, edge_p2):
            return True

    return False


@dataclass
class Stock:
    """Represents a stock (reservoir) in the model."""
    name: str
    initial_value: str
    units: str = ""
    inflows: list[str] = field(default_factory=list)
    outflows: list[str] = field(default_factory=list)
    non_negative: bool = True
    x: Optional[float] = None  # None means auto-position
    y: Optional[float] = None  # None means auto-position
    width: int = 45  # Default stock width
    height: int = 35  # Default stock height


@dataclass
class Flow:
    """Represents a flow between stocks."""
    name: str
    equation: str
    units: str = ""
    from_stock: Optional[str] = None  # None means external source
    to_stock: Optional[str] = None    # None means external sink
    non_negative: bool = True
    x: Optional[float] = None  # None means auto-position
    y: Optional[float] = None  # None means auto-position
    points: list[tuple[float, float]] = field(default_factory=list)
    graphical_function: Optional["GraphicalFunction"] = None


@dataclass
class Aux:
    """Represents an auxiliary variable."""
    name: str
    equation: str
    units: str = ""
    x: Optional[float] = None  # None means auto-position
    y: Optional[float] = None  # None means auto-position
    graphical_function: Optional["GraphicalFunction"] = None


@dataclass
class GraphicalFunction:
    """Represents a graphical function (lookup table) definition."""
    ypts: list[float]
    xscale: Optional[tuple[float, float]] = None
    xpts: Optional[list[float]] = None
    yscale: Optional[tuple[float, float]] = None
    gf_type: Optional[str] = None


@dataclass
class Connector:
    """Represents a dependency connector between variables."""
    uid: int
    from_var: str
    to_var: str
    angle: float = 0


@dataclass
class SimSpecs:
    """Simulation specifications."""
    start: float = 0
    stop: float = 100
    dt: float = 0.25
    method: str = "Euler"
    time_units: str = "Years"


class StellaModel:
    """Represents a complete Stella system dynamics model."""

    def __init__(self, name: str = "Untitled"):
        self.name = name
        self.uuid = str(uuid.uuid4())
        self.sim_specs = SimSpecs()
        self.stocks: dict[str, Stock] = {}
        self.flows: dict[str, Flow] = {}
        self.auxs: dict[str, Aux] = {}
        self.connectors: list[Connector] = []
        self._connector_uid = 0

    def _next_connector_uid(self) -> int:
        """Get the next unique connector ID."""
        self._connector_uid += 1
        return self._connector_uid

    def _normalize_name(self, name: str) -> str:
        """Convert display name to internal name (spaces to underscores)."""
        return name.replace(" ", "_")

    def _display_name(self, name: str) -> str:
        """Convert internal name to display name (underscores to spaces)."""
        return name.replace("_", " ")

    def _extract_variable_refs(self, equation: str) -> set[str]:
        """Extract variable names referenced in an equation.

        Returns normalized variable names (spaces converted to underscores).
        Filters out Stella built-in functions and keywords.
        """
        if not equation:
            return set()

        # Extract potential variable names (alphanumeric with underscores only)
        # Note: Stella allows spaces in variable names, but in equations they should
        # be written with underscores or quoted. We extract standard identifiers.
        tokens = re.findall(r'\b([A-Za-z_][A-Za-z0-9_]*)\b', equation)

        refs = set()
        for token in tokens:
            # Check if it's a function or keyword (case-insensitive)
            if token.upper() not in STELLA_FUNCTIONS:
                # Try to filter out pure numbers
                try:
                    float(token)
                except ValueError:
                    # Normalize (in case there are any spaces, though regex won't match them)
                    refs.add(self._normalize_name(token))

        return refs

    def _build_dependency_graph(self) -> tuple[dict[str, set[str]], dict[str, set[str]]]:
        """Build bidirectional adjacency lists from connectors and flow-stock relationships.

        Returns:
            (outgoing, incoming) where:
            - outgoing[node] = set of nodes this node connects TO
            - incoming[node] = set of nodes that connect TO this node
        """
        from collections import defaultdict

        outgoing: dict[str, set[str]] = defaultdict(set)
        incoming: dict[str, set[str]] = defaultdict(set)

        # Initialize all elements as nodes
        for name in self.stocks:
            outgoing.setdefault(name, set())
            incoming.setdefault(name, set())
        for name in self.flows:
            outgoing.setdefault(name, set())
            incoming.setdefault(name, set())
        for name in self.auxs:
            outgoing.setdefault(name, set())
            incoming.setdefault(name, set())

        # Add connector edges
        for conn in self.connectors:
            from_var = conn.from_var
            to_var = conn.to_var
            if from_var in outgoing and to_var in incoming:
                outgoing[from_var].add(to_var)
                incoming[to_var].add(from_var)

        # Add implicit flow-stock edges
        for name, flow in self.flows.items():
            if flow.from_stock and flow.from_stock in self.stocks:
                outgoing[flow.from_stock].add(name)
                incoming[name].add(flow.from_stock)
            if flow.to_stock and flow.to_stock in self.stocks:
                outgoing[name].add(flow.to_stock)
                incoming[flow.to_stock].add(name)

        return dict(outgoing), dict(incoming)

    def _find_subsystems(self, outgoing: dict[str, set[str]], incoming: dict[str, set[str]]) -> list[set[str]]:
        """Find connected components (subsystems) in the graph.

        Returns list of node sets, sorted by size (largest first).
        """
        all_nodes = set(self.stocks) | set(self.flows) | set(self.auxs)
        visited: set[str] = set()
        subsystems: list[set[str]] = []

        # Build undirected graph for component detection
        undirected: dict[str, set[str]] = {node: set() for node in all_nodes}
        for node, neighbors in outgoing.items():
            for neighbor in neighbors:
                if neighbor in undirected:
                    undirected[node].add(neighbor)
                    undirected[neighbor].add(node)
        for node, neighbors in incoming.items():
            for neighbor in neighbors:
                if neighbor in undirected:
                    undirected[node].add(neighbor)
                    undirected[neighbor].add(node)

        def dfs(node: str, component: set[str]):
            if node in visited:
                return
            visited.add(node)
            component.add(node)
            for neighbor in undirected.get(node, set()):
                dfs(neighbor, component)

        for node in sorted(all_nodes):  # Sorted for determinism
            if node not in visited:
                component: set[str] = set()
                dfs(node, component)
                if component:
                    subsystems.append(component)

        return sorted(subsystems, key=len, reverse=True)

    def _find_stock_chains(self, subsystem: set[str]) -> list[list[str]]:
        """Find linear chains of stocks connected by flows within a subsystem.

        Returns list of chains, where each chain is ordered by flow direction.
        Stocks are ordered by incoming flow count (true sources first).
        """
        stocks_in_subsystem = [s for s in subsystem if s in self.stocks]
        if not stocks_in_subsystem:
            return []

        # Build directed stock adjacency from flows
        stock_adj: dict[str, list[str]] = {s: [] for s in stocks_in_subsystem}
        incoming_count: dict[str, int] = {s: 0 for s in stocks_in_subsystem}

        for flow in self.flows.values():
            if (flow.from_stock in stocks_in_subsystem
                    and flow.to_stock in stocks_in_subsystem):
                if flow.to_stock not in stock_adj[flow.from_stock]:
                    stock_adj[flow.from_stock].append(flow.to_stock)
                    incoming_count[flow.to_stock] += 1

        # Find true sources: sort by incoming_count ASC, then alphabetically for determinism
        # This ensures stocks with 0 incoming flows come first (true sources)
        start_candidates = sorted(
            stocks_in_subsystem,
            key=lambda s: (incoming_count[s], s)
        )

        # Trace chains starting from true sources
        chains: list[list[str]] = []
        visited: set[str] = set()

        for start in start_candidates:
            if start in visited:
                continue

            chain = []
            current: Optional[str] = start

            while current and current not in visited:
                visited.add(current)
                chain.append(current)
                # Follow flow direction to next stock (sorted for determinism)
                next_stocks = sorted(stock_adj.get(current, []))
                unvisited_next = [s for s in next_stocks if s not in visited]
                current = unvisited_next[0] if unvisited_next else None

            if chain:
                chains.append(chain)

        return chains

    def _find_aux_target_position(
        self,
        aux_name: str,
        outgoing: dict[str, set[str]]
    ) -> Optional[tuple[float, float]]:
        """Find the position where an aux should be placed (near its targets).

        Returns average position of all targets, or None if no targets have positions.
        """
        target_positions: list[tuple[float, float]] = []

        # Find what this aux connects to via connectors
        for conn in self.connectors:
            if conn.from_var == aux_name:
                target = conn.to_var
                if target in self.stocks:
                    stock = self.stocks[target]
                    if stock.x is not None and stock.y is not None:
                        target_positions.append((stock.x, stock.y))
                elif target in self.flows:
                    flow = self.flows[target]
                    if flow.x is not None and flow.y is not None:
                        target_positions.append((flow.x, flow.y))
                elif target in self.auxs:
                    aux = self.auxs[target]
                    if aux.x is not None and aux.y is not None:
                        target_positions.append((aux.x, aux.y))

        if target_positions:
            avg_x = sum(p[0] for p in target_positions) / len(target_positions)
            avg_y = sum(p[1] for p in target_positions) / len(target_positions)
            return (avg_x, avg_y)
        return None

    def _position_subsystem(
        self,
        subsystem: set[str],
        outgoing: dict[str, set[str]],
        incoming: dict[str, set[str]],
        start_x: float,
        stock_y: float,
        aux_y: float,
        stock_spacing: float,
        aux_spacing: float
    ) -> tuple[float, float, float, float]:
        """Position all elements in a subsystem. Returns bounding box (min_x, min_y, max_x, max_y)."""

        # Find stock chains
        chains = self._find_stock_chains(subsystem)

        # Position stocks in chains horizontally
        x = start_x
        for chain in chains:
            for stock_name in chain:
                stock = self.stocks[stock_name]
                if stock.x is None or stock.y is None:
                    stock.x = x
                    stock.y = stock_y
                x += stock_spacing

        # Position flows between their stocks
        flows_in_subsystem = [f for f in self.flows if f in subsystem]
        for flow_name in flows_in_subsystem:
            flow = self.flows[flow_name]
            if flow.x is not None and flow.y is not None:
                continue

            from_stock = self.stocks.get(flow.from_stock) if flow.from_stock else None
            to_stock = self.stocks.get(flow.to_stock) if flow.to_stock else None

            if from_stock and to_stock:
                from_x = from_stock.x if from_stock.x is not None else start_x
                to_x = to_stock.x if to_stock.x is not None else start_x
                from_y = from_stock.y if from_stock.y is not None else stock_y
                to_y = to_stock.y if to_stock.y is not None else stock_y
                flow.x = (from_x + to_x) / 2
                flow.y = (from_y + to_y) / 2
            elif from_stock:
                flow.x = (from_stock.x or start_x) + 90
                flow.y = from_stock.y or stock_y
            elif to_stock:
                flow.x = (to_stock.x or start_x) - 90
                flow.y = to_stock.y or stock_y

        # Position auxs near their targets with horizontal spread for same-target overlap
        auxs_in_subsystem = [a for a in sorted(subsystem) if a in self.auxs]
        positioned_auxs: set[str] = set()
        AUX_SPREAD_SPACING = 70  # pixels between auxs targeting the same element

        # First pass: group auxs by their primary connector target
        target_to_auxs: dict[str, list[str]] = {}
        for aux_name in auxs_in_subsystem:
            aux = self.auxs[aux_name]
            if aux.x is not None and aux.y is not None:
                positioned_auxs.add(aux_name)
                continue

            # Find targets via connectors
            targets = [c.to_var for c in self.connectors if c.from_var == aux_name]
            if targets:
                # Use first target (sorted for determinism)
                primary_target = sorted(targets)[0]
                if primary_target not in target_to_auxs:
                    target_to_auxs[primary_target] = []
                target_to_auxs[primary_target].append(aux_name)

        # Position each group with horizontal spread centered on target
        for target, aux_names in target_to_auxs.items():
            # Get target position
            target_x: Optional[float] = None
            target_y: Optional[float] = None
            if target in self.stocks:
                target_x = self.stocks[target].x
                target_y = self.stocks[target].y
            elif target in self.flows:
                target_x = self.flows[target].x
                target_y = self.flows[target].y
            elif target in self.auxs:
                target_x = self.auxs[target].x
                target_y = self.auxs[target].y

            if target_x is None:
                continue

            # Sort aux names for determinism
            aux_names = sorted(aux_names)
            n = len(aux_names)

            # Center the group horizontally above target
            # Offsets: for n=3 -> [-70, 0, +70], for n=2 -> [-35, +35]
            start_offset = -((n - 1) * AUX_SPREAD_SPACING) / 2

            for i, aux_name in enumerate(aux_names):
                aux = self.auxs[aux_name]
                aux.x = target_x + start_offset + (i * AUX_SPREAD_SPACING)
                # Position above flows (80px), at aux_y for stocks
                if target in self.flows and target_y is not None:
                    aux.y = target_y - 80
                else:
                    aux.y = aux_y
                positioned_auxs.add(aux_name)

        # Second pass: position remaining auxs (those without connector targets)
        # These may be referenced in flow equations
        aux_x_offset = 0
        for aux_name in auxs_in_subsystem:
            if aux_name in positioned_auxs:
                continue

            aux = self.auxs[aux_name]
            if aux.x is not None and aux.y is not None:
                continue

            # Check if this aux is referenced by any flow (it's a parameter for that flow)
            for flow_name, flow in self.flows.items():
                if flow_name not in subsystem:
                    continue
                flow_refs = self._extract_variable_refs(flow.equation)
                if aux_name in flow_refs:
                    # Position near the flow
                    if flow.x is not None:
                        aux.x = flow.x + aux_x_offset
                        aux.y = aux_y
                        aux_x_offset += aux_spacing
                        positioned_auxs.add(aux_name)
                        break

        # Third pass: any remaining auxs go in a row at start_x
        remaining_x = start_x
        for aux_name in auxs_in_subsystem:
            if aux_name in positioned_auxs:
                continue

            aux = self.auxs[aux_name]
            if aux.x is None or aux.y is None:
                aux.x = remaining_x
                aux.y = aux_y - 60  # Put orphan auxs above the main aux row
                remaining_x += aux_spacing

        # Calculate bounding box
        all_x: list[float] = []
        all_y: list[float] = []

        for name in subsystem:
            if name in self.stocks and self.stocks[name].x is not None:
                all_x.append(self.stocks[name].x)  # type: ignore
                all_y.append(self.stocks[name].y)  # type: ignore
            if name in self.flows and self.flows[name].x is not None:
                all_x.append(self.flows[name].x)  # type: ignore
                all_y.append(self.flows[name].y)  # type: ignore
            if name in self.auxs and self.auxs[name].x is not None:
                all_x.append(self.auxs[name].x)  # type: ignore
                all_y.append(self.auxs[name].y)  # type: ignore

        if all_x and all_y:
            return (min(all_x), min(all_y), max(all_x), max(all_y))
        return (start_x, aux_y, start_x + stock_spacing, stock_y)

    def _arrange_subsystems(
        self,
        subsystems: list[set[str]],
        bounds: list[tuple[float, float, float, float]],
        gap: float
    ):
        """Arrange subsystems: largest stays in place, smaller ones offset to the right."""
        if len(subsystems) <= 1:
            return

        # First subsystem (largest) stays in place
        # Offset subsequent subsystems to the right
        current_x = bounds[0][2] + gap  # max_x of first + gap

        for i, subsystem in enumerate(subsystems[1:], start=1):
            min_x = bounds[i][0]
            max_x = bounds[i][2]
            offset_x = current_x - min_x

            # Shift all elements in this subsystem
            for name in subsystem:
                if name in self.stocks and self.stocks[name].x is not None:
                    self.stocks[name].x += offset_x
                if name in self.flows and self.flows[name].x is not None:
                    self.flows[name].x += offset_x
                if name in self.auxs and self.auxs[name].x is not None:
                    self.auxs[name].x += offset_x

            current_x = current_x + (max_x - min_x) + gap

    def add_stock(
        self,
        name: str,
        initial_value: str,
        units: str = "",
        inflows: Optional[list[str]] = None,
        outflows: Optional[list[str]] = None,
        non_negative: bool = True,
        x: Optional[float] = None,
        y: Optional[float] = None
    ) -> Stock:
        """Add a stock to the model."""
        stock = Stock(
            name=name,
            initial_value=initial_value,
            units=units,
            inflows=[self._normalize_name(f) for f in (inflows or [])],
            outflows=[self._normalize_name(f) for f in (outflows or [])],
            non_negative=non_negative,
            x=x,
            y=y
        )
        self.stocks[self._normalize_name(name)] = stock
        return stock

    def add_flow(
        self,
        name: str,
        equation: str,
        units: str = "",
        from_stock: Optional[str] = None,
        to_stock: Optional[str] = None,
        non_negative: bool = True,
        x: Optional[float] = None,
        y: Optional[float] = None,
        graphical_function: Optional[GraphicalFunction] = None
    ) -> Flow:
        """Add a flow to the model."""
        flow = Flow(
            name=name,
            equation=equation,
            units=units,
            from_stock=self._normalize_name(from_stock) if from_stock else None,
            to_stock=self._normalize_name(to_stock) if to_stock else None,
            non_negative=non_negative,
            x=x,
            y=y,
            graphical_function=graphical_function,
        )
        self.flows[self._normalize_name(name)] = flow

        # Update stock inflows/outflows
        if from_stock:
            from_key = self._normalize_name(from_stock)
            if from_key in self.stocks:
                flow_key = self._normalize_name(name)
                if flow_key not in self.stocks[from_key].outflows:
                    self.stocks[from_key].outflows.append(flow_key)

        if to_stock:
            to_key = self._normalize_name(to_stock)
            if to_key in self.stocks:
                flow_key = self._normalize_name(name)
                if flow_key not in self.stocks[to_key].inflows:
                    self.stocks[to_key].inflows.append(flow_key)

        return flow

    def add_aux(
        self,
        name: str,
        equation: str,
        units: str = "",
        x: Optional[float] = None,
        y: Optional[float] = None,
        graphical_function: Optional[GraphicalFunction] = None
    ) -> Aux:
        """Add an auxiliary variable to the model."""
        aux = Aux(
            name=name,
            equation=equation,
            units=units,
            x=x,
            y=y,
            graphical_function=graphical_function,
        )
        self.auxs[self._normalize_name(name)] = aux
        return aux

    def add_connector(self, from_var: str, to_var: str) -> Connector:
        """Add a connector (dependency) between variables."""
        connector = Connector(
            uid=self._next_connector_uid(),
            from_var=self._normalize_name(from_var),
            to_var=self._normalize_name(to_var)
        )
        self.connectors.append(connector)
        return connector

    def _calculate_stock_sizes(self):
        """Calculate appropriate width/height for each stock based on connectivity.

        Stocks with more flows get larger to allow visual separation of flow attachments.
        Maintains a pleasing aspect ratio (roughly 1.3:1 width:height).
        """
        MIN_WIDTH = 45
        MAX_WIDTH = 120
        MIN_HEIGHT = 35
        MAX_HEIGHT = 90
        FLOW_WIDTH_CONTRIBUTION = 15  # Extra width per flow beyond 2
        ASPECT_RATIO = 1.3  # width:height ratio

        for stock in self.stocks.values():
            num_flows = len(stock.inflows) + len(stock.outflows)

            # Start at minimum, add width for extra flows
            width = MIN_WIDTH + max(0, num_flows - 2) * FLOW_WIDTH_CONTRIBUTION
            width = min(width, MAX_WIDTH)

            # Scale height to maintain aspect ratio
            height = int(width / ASPECT_RATIO)
            height = max(MIN_HEIGHT, min(height, MAX_HEIGHT))

            stock.width = width
            stock.height = height

    def _calculate_stock_spacing(self) -> int:
        """Calculate appropriate spacing between stocks based on their sizes.

        Returns the minimum spacing needed to prevent overlap and allow
        room for flows, valves, and labels.
        """
        if not self.stocks:
            return 200  # Default

        # Find the maximum stock width
        max_width = max(s.width for s in self.stocks.values())

        # Spacing = max stock width + gap for flows/valves/labels
        # Gap should be at least 100px for flow valve + labels
        MIN_GAP = 100

        return max_width + MIN_GAP

    def _auto_layout(self):
        """Auto-arrange visual positions using graph-based hierarchical layout.

        Uses connector relationships to position elements:
        1. Calculates stock sizes based on connectivity
        2. Calculates dynamic spacing based on stock sizes
        3. Builds dependency graph from connectors
        4. Detects subsystems (connected components)
        5. Positions stocks in chains following flow direction
        6. Positions auxs near their connector targets
        7. Separates independent subsystems visually

        Falls back to simple row layout if no connectors exist.
        Always recalculates flow.points to ensure flows connect to stocks correctly.
        """
        # Calculate stock sizes first (affects spacing and flow attachment)
        self._calculate_stock_sizes()

        # Calculate dynamic spacing based on stock sizes
        STOCK_SPACING = self._calculate_stock_spacing()

        # Layout constants
        AUX_SPACING = 80
        STOCK_Y = 300
        AUX_Y = 150
        START_X = 200
        SUBSYSTEM_GAP = 250

        # Build dependency graph from connectors
        outgoing, incoming = self._build_dependency_graph()

        # Find subsystems (connected components)
        subsystems = self._find_subsystems(outgoing, incoming)

        # Position each subsystem
        subsystem_bounds: list[tuple[float, float, float, float]] = []

        for subsystem in subsystems:
            bounds = self._position_subsystem(
                subsystem, outgoing, incoming,
                START_X, STOCK_Y, AUX_Y, STOCK_SPACING, AUX_SPACING
            )
            subsystem_bounds.append(bounds)

        # Arrange subsystems relative to each other (largest centered, others offset)
        if len(subsystems) > 1 and len(subsystem_bounds) > 1:
            self._arrange_subsystems(subsystems, subsystem_bounds, SUBSYSTEM_GAP)

        # Always recalculate flow points to connect stocks at their actual positions
        self._recalculate_flow_points(START_X, STOCK_Y)

        # Calculate connector angles based on final positions
        self._calculate_connector_angles()

        # Resolve any layout violations (collisions, crossings)
        self._resolve_layout_violations()

    def _calculate_flow_offset(self, index: int, total: int) -> float:
        """Calculate vertical offset for flow attachment point.

        When multiple flows share a stock endpoint, offset them vertically
        to prevent overlap. Centers the group around the stock center.

        Args:
            index: This flow's index in the group (0-based)
            total: Total number of flows in the group

        Returns:
            Vertical offset in pixels (positive = down, negative = up)
        """
        if total <= 1:
            return 0.0
        # Center the group: e.g., 3 flows -> offsets of -20, 0, +20
        return (index - (total - 1) / 2) * 20.0

    def _recalculate_flow_points(self, start_x: float, stock_y: float):
        """Recalculate flow.points to connect stocks at their actual positions.

        Uses orthogonal (L-shaped/U-shaped) routing for multiple flows from the
        same stock to prevent visual overlap. Routes flows in the direction of
        their actual destination (up for destinations above, down for below).
        """
        ROUTE_OFFSET = 40  # Vertical offset between routing paths (pixels)

        # Group flows by their source and destination stocks
        outflows_by_stock: dict[str, list[str]] = {}
        inflows_by_stock: dict[str, list[str]] = {}

        for name, flow in self.flows.items():
            if flow.from_stock:
                if flow.from_stock not in outflows_by_stock:
                    outflows_by_stock[flow.from_stock] = []
                outflows_by_stock[flow.from_stock].append(name)
            if flow.to_stock:
                if flow.to_stock not in inflows_by_stock:
                    inflows_by_stock[flow.to_stock] = []
                inflows_by_stock[flow.to_stock].append(name)

        # Sort flow lists for determinism
        for stock_name in outflows_by_stock:
            outflows_by_stock[stock_name].sort()
        for stock_name in inflows_by_stock:
            inflows_by_stock[stock_name].sort()

        # Calculate flow points with orthogonal routing
        for name, flow in self.flows.items():
            from_stock = self.stocks.get(flow.from_stock) if flow.from_stock else None
            to_stock = self.stocks.get(flow.to_stock) if flow.to_stock else None

            if from_stock and to_stock:
                from_x = from_stock.x if from_stock.x is not None else start_x
                from_y = from_stock.y if from_stock.y is not None else stock_y
                to_x = to_stock.x if to_stock.x is not None else start_x
                to_y = to_stock.y if to_stock.y is not None else stock_y

                # Use actual stock dimensions for attachment points
                from_half_width = from_stock.width / 2
                from_half_height = from_stock.height / 2
                to_half_width = to_stock.width / 2
                to_half_height = to_stock.height / 2

                # Get flow's position in outflow group
                outflows = outflows_by_stock.get(flow.from_stock, [name])
                total = len(outflows)

                # Calculate entry offset for inflows
                to_offset = 0.0
                if flow.to_stock and flow.to_stock in inflows_by_stock:
                    inflows = inflows_by_stock[flow.to_stock]
                    if len(inflows) > 1:
                        inflow_index = inflows.index(name)
                        to_offset = self._calculate_flow_offset(inflow_index, len(inflows))

                if total == 1:
                    # Single flow: straight 2-point path
                    flow.points = [
                        (from_x + from_half_width, from_y),
                        (to_x - to_half_width, to_y + to_offset)
                    ]
                else:
                    # Multiple flows: route based on destination position
                    # Group flows by direction (up/same/down) relative to source
                    flows_up = []
                    flows_same = []
                    flows_down = []

                    for flow_name in outflows:
                        f = self.flows[flow_name]
                        dest = self.stocks.get(f.to_stock)
                        if dest and dest.y is not None:
                            dest_y = dest.y
                            if dest_y < from_y - 20:
                                flows_up.append(flow_name)
                            elif dest_y > from_y + 20:
                                flows_down.append(flow_name)
                            else:
                                flows_same.append(flow_name)
                        else:
                            flows_same.append(flow_name)

                    # Determine this flow's direction and index within that group
                    if name in flows_up:
                        go_up = True
                        group_index = flows_up.index(name)
                        group_size = len(flows_up)
                    elif name in flows_down:
                        go_up = False
                        group_index = flows_down.index(name)
                        group_size = len(flows_down)
                    else:
                        # Same level: alternate up/down
                        same_index = flows_same.index(name)
                        if same_index == 0:
                            # First same-level flow: straight
                            flow.points = [
                                (from_x + from_half_width, from_y),
                                (to_x - to_half_width, to_y + to_offset)
                            ]
                            continue
                        go_up = (same_index % 2 == 1)
                        group_index = (same_index - 1) // 2
                        group_size = (len(flows_same) - 1) // 2 + 1

                    # Calculate routing offset based on position in group
                    offset = (group_index + 1) * ROUTE_OFFSET
                    route_y = from_y - from_half_height - offset if go_up else from_y + from_half_height + offset

                    # Offset entry X based on position in direction group to prevent overlap
                    ENTRY_X_OFFSET = 30  # Horizontal spacing between entry points
                    exit_x = from_x + from_half_width
                    exit_y = from_y - from_half_height if go_up else from_y + from_half_height
                    entry_x = to_x - to_half_width - (group_index * ENTRY_X_OFFSET)
                    entry_y = to_y - to_half_height if go_up else to_y + to_half_height

                    flow.points = [
                        (exit_x, exit_y),       # Exit from stock edge
                        (exit_x, route_y),      # Go up/down
                        (entry_x, route_y),     # Go across
                        (entry_x, entry_y)      # Go to destination edge
                    ]

            elif from_stock:
                # Source-only flow (external sink)
                from_x = from_stock.x if from_stock.x is not None else start_x
                from_y = from_stock.y if from_stock.y is not None else stock_y
                from_half_width = from_stock.width / 2

                # Calculate offset for this flow among other outflows
                from_offset = 0.0
                if flow.from_stock and flow.from_stock in outflows_by_stock:
                    outflows = outflows_by_stock[flow.from_stock]
                    if len(outflows) > 1:
                        index = outflows.index(name)
                        from_offset = self._calculate_flow_offset(index, len(outflows))

                flow.points = [
                    (from_x + from_half_width, from_y + from_offset),
                    (from_x + 160, from_y + from_offset)
                ]

            elif to_stock:
                # Sink-only flow (external source)
                to_x = to_stock.x if to_stock.x is not None else start_x
                to_y = to_stock.y if to_stock.y is not None else stock_y
                to_half_width = to_stock.width / 2

                # Calculate offset for this flow among other inflows
                to_offset = 0.0
                if flow.to_stock and flow.to_stock in inflows_by_stock:
                    inflows = inflows_by_stock[flow.to_stock]
                    if len(inflows) > 1:
                        index = inflows.index(name)
                        to_offset = self._calculate_flow_offset(index, len(inflows))

                flow.points = [
                    (to_x - 160, to_y + to_offset),
                    (to_x - to_half_width, to_y + to_offset)
                ]

    def _calculate_connector_angles(self):
        """Calculate connector angles based on source and target positions.

        Uses atan2 to compute the angle from source to target.
        Convention: degrees, 0 = right, counter-clockwise positive.
        Note: -dy because screen y-coordinates increase downward.
        """
        # Build position lookup for all elements
        positions: dict[str, tuple[float, float]] = {}
        for name, stock in self.stocks.items():
            if stock.x is not None and stock.y is not None:
                positions[name] = (stock.x, stock.y)
        for name, flow in self.flows.items():
            if flow.x is not None and flow.y is not None:
                positions[name] = (flow.x, flow.y)
        for name, aux in self.auxs.items():
            if aux.x is not None and aux.y is not None:
                positions[name] = (aux.x, aux.y)

        for conn in self.connectors:
            from_pos = positions.get(conn.from_var)
            to_pos = positions.get(conn.to_var)

            if from_pos and to_pos:
                dx = to_pos[0] - from_pos[0]
                dy = to_pos[1] - from_pos[1]

                # Handle zero distance (same position)
                if abs(dx) < 0.001 and abs(dy) < 0.001:
                    conn.angle = 0
                else:
                    # -dy because y increases downward in screen coordinates
                    conn.angle = math.degrees(math.atan2(-dy, dx))

    # =========================================================================
    # Layout Collision/Crossing Detection and Resolution
    # =========================================================================

    def _get_element_box(self, name: str) -> Optional[BoundingBox]:
        """Get bounding box for any model element."""
        if name in self.stocks:
            stock = self.stocks[name]
            if stock.x is not None and stock.y is not None:
                return BoundingBox(stock.x, stock.y, stock.width, stock.height)
        elif name in self.flows:
            flow = self.flows[name]
            if flow.x is not None and flow.y is not None:
                # Flow valve is roughly 20x20
                return BoundingBox(flow.x, flow.y, 20, 20)
        elif name in self.auxs:
            aux = self.auxs[name]
            if aux.x is not None and aux.y is not None:
                return BoundingBox(aux.x, aux.y, AUX_RADIUS * 2, AUX_RADIUS * 2)
        return None

    def _get_all_bounding_boxes(self) -> dict[str, BoundingBox]:
        """Get bounding boxes for all positioned elements."""
        boxes: dict[str, BoundingBox] = {}
        for name in self.stocks:
            box = self._get_element_box(name)
            if box:
                boxes[name] = box
        for name in self.auxs:
            box = self._get_element_box(name)
            if box:
                boxes[name] = box
        # Note: flows are not included as their position is the valve,
        # and flow lines are handled separately
        return boxes

    def _get_connector_segments(self) -> dict[int, tuple[tuple[float, float], tuple[float, float]]]:
        """Get line segments for all connectors (from source to target position)."""
        segments: dict[int, tuple[tuple[float, float], tuple[float, float]]] = {}

        # Build position lookup
        positions: dict[str, tuple[float, float]] = {}
        for name, stock in self.stocks.items():
            if stock.x is not None and stock.y is not None:
                positions[name] = (stock.x, stock.y)
        for name, flow in self.flows.items():
            if flow.x is not None and flow.y is not None:
                positions[name] = (flow.x, flow.y)
        for name, aux in self.auxs.items():
            if aux.x is not None and aux.y is not None:
                positions[name] = (aux.x, aux.y)

        for conn in self.connectors:
            from_pos = positions.get(conn.from_var)
            to_pos = positions.get(conn.to_var)
            if from_pos and to_pos:
                segments[conn.uid] = (from_pos, to_pos)

        return segments

    def _get_flow_segments(self) -> dict[str, list[tuple[tuple[float, float], tuple[float, float]]]]:
        """Get line segments for all flow paths."""
        segments: dict[str, list[tuple[tuple[float, float], tuple[float, float]]]] = {}

        for name, flow in self.flows.items():
            if flow.points and len(flow.points) >= 2:
                flow_segs: list[tuple[tuple[float, float], tuple[float, float]]] = []
                for i in range(len(flow.points) - 1):
                    flow_segs.append((flow.points[i], flow.points[i + 1]))
                segments[name] = flow_segs

        return segments

    def _detect_aux_collisions(self) -> list[tuple[str, str]]:
        """Detect pairs of auxs that overlap."""
        collisions: list[tuple[str, str]] = []
        aux_names = list(self.auxs.keys())

        for i, name1 in enumerate(aux_names):
            box1 = self._get_element_box(name1)
            if not box1:
                continue
            for name2 in aux_names[i + 1:]:
                box2 = self._get_element_box(name2)
                if box2 and box1.intersects(box2, margin=5):
                    collisions.append((name1, name2))

        return collisions

    def _detect_connector_flow_crossings(self) -> list[tuple[int, str]]:
        """Detect connectors that cross flow lines. Returns (connector_uid, flow_name) pairs.

        Note: A connector is expected to touch its target flow, so we skip checking
        if a connector crosses the flow it's connected TO.
        """
        crossings: list[tuple[int, str]] = []

        connector_segs = self._get_connector_segments()
        flow_segs = self._get_flow_segments()

        # Build map of connector uid -> target name
        conn_targets: dict[int, str] = {}
        for conn in self.connectors:
            conn_targets[conn.uid] = conn.to_var

        for conn_uid, (cp1, cp2) in connector_segs.items():
            target = conn_targets.get(conn_uid)
            for flow_name, segments in flow_segs.items():
                # Skip if this is the connector's target flow
                if flow_name == target:
                    continue
                for fp1, fp2 in segments:
                    if segments_intersect(cp1, cp2, fp1, fp2):
                        crossings.append((conn_uid, flow_name))
                        break  # One crossing per connector-flow pair is enough

        return crossings

    def _detect_flow_stock_crossings(self) -> list[tuple[str, str]]:
        """Detect flows that pass through stocks (not their source/dest). Returns (flow_name, stock_name) pairs."""
        crossings: list[tuple[str, str]] = []

        flow_segs = self._get_flow_segments()

        for flow_name, segments in flow_segs.items():
            flow = self.flows[flow_name]
            for stock_name, stock in self.stocks.items():
                # Skip source and destination stocks
                if stock_name in (flow.from_stock, flow.to_stock):
                    continue

                box = self._get_element_box(stock_name)
                if not box:
                    continue

                for p1, p2 in segments:
                    if segment_intersects_box(p1, p2, box):
                        crossings.append((flow_name, stock_name))
                        break

        return crossings

    def _detect_connector_stock_crossings(self) -> list[tuple[int, str]]:
        """Detect connectors that pass through stocks. Returns (connector_uid, stock_name) pairs.

        Skips stocks that are the source or target of the connector.
        """
        crossings: list[tuple[int, str]] = []

        connector_segs = self._get_connector_segments()

        # Build map of connector uid -> (from_var, to_var)
        conn_endpoints: dict[int, tuple[str, str]] = {}
        for conn in self.connectors:
            conn_endpoints[conn.uid] = (conn.from_var, conn.to_var)

        for conn_uid, (cp1, cp2) in connector_segs.items():
            from_var, to_var = conn_endpoints.get(conn_uid, ("", ""))
            for stock_name, stock in self.stocks.items():
                # Skip if this stock is the source or target of the connector
                if stock_name in (from_var, to_var):
                    continue

                box = self._get_element_box(stock_name)
                if not box:
                    continue

                if segment_intersects_box(cp1, cp2, box):
                    crossings.append((conn_uid, stock_name))

        return crossings

    def _separate_auxs(self, name1: str, name2: str):
        """Push two overlapping auxs apart."""
        aux1 = self.auxs.get(name1)
        aux2 = self.auxs.get(name2)

        if not aux1 or not aux2:
            return
        if aux1.x is None or aux1.y is None or aux2.x is None or aux2.y is None:
            return

        # Direction from aux1 to aux2
        dx = aux2.x - aux1.x
        dy = aux2.y - aux1.y
        dist = math.sqrt(dx * dx + dy * dy)

        if dist < 0.001:
            # Same position - push horizontally
            dx, dy = 1.0, 0.0
            dist = 1.0

        # Minimum distance is 2 * radius + margin
        min_dist = AUX_RADIUS * 2 + 10
        if dist >= min_dist:
            return  # Already separated

        # Push each aux half the needed distance
        push = (min_dist - dist) / 2 + 2
        aux1.x -= push * dx / dist
        aux1.y -= push * dy / dist
        aux2.x += push * dx / dist
        aux2.y += push * dy / dist

    def _reposition_aux_to_avoid_crossing(self, conn_uid: int, obstacle_name: str, obstacle_type: str = "flow"):
        """Move aux so its connector doesn't cross the specified obstacle.

        Args:
            conn_uid: The connector's unique ID
            obstacle_name: Name of the flow or stock to avoid
            obstacle_type: Either "flow" or "stock"
        """
        # Find the connector and its source aux
        conn = None
        for c in self.connectors:
            if c.uid == conn_uid:
                conn = c
                break

        if not conn:
            return

        aux = self.auxs.get(conn.from_var)
        if not aux or aux.x is None or aux.y is None:
            return

        # Get target position and size for proportional offsets
        target_pos: Optional[tuple[float, float]] = None
        target_size = 45  # Default size for offset calculation
        if conn.to_var in self.stocks:
            stock = self.stocks[conn.to_var]
            if stock.x is not None and stock.y is not None:
                target_pos = (stock.x, stock.y)
                target_size = max(stock.width, stock.height)
        elif conn.to_var in self.flows:
            flow = self.flows[conn.to_var]
            if flow.x is not None and flow.y is not None:
                target_pos = (flow.x, flow.y)
                target_size = 20  # Flow valve size
        elif conn.to_var in self.auxs:
            other_aux = self.auxs[conn.to_var]
            if other_aux.x is not None and other_aux.y is not None:
                target_pos = (other_aux.x, other_aux.y)
                target_size = AUX_RADIUS * 2

        if not target_pos:
            return

        # Build obstacle check function based on type
        if obstacle_type == "flow":
            flow_segs = self._get_flow_segments().get(obstacle_name, [])
            if not flow_segs:
                return

            def crosses_obstacle(candidate: tuple[float, float]) -> bool:
                for fp1, fp2 in flow_segs:
                    if segments_intersect(candidate, target_pos, fp1, fp2):
                        return True
                return False
        else:  # stock
            stock_box = self._get_element_box(obstacle_name)
            if not stock_box:
                return

            def crosses_obstacle(candidate: tuple[float, float]) -> bool:
                return segment_intersects_box(candidate, target_pos, stock_box)

        # Calculate proportional offsets based on target element size
        base_offset = target_size + AUX_RADIUS + 20
        diag_offset = int(base_offset * 0.75)
        far_offset = int(base_offset * 1.5)

        offsets = [
            (0, -base_offset), (0, base_offset),  # above, below
            (-base_offset, 0), (base_offset, 0),  # left, right
            (-diag_offset, -diag_offset), (diag_offset, -diag_offset),  # diagonal up
            (-diag_offset, diag_offset), (diag_offset, diag_offset),  # diagonal down
            (0, -far_offset), (0, far_offset),  # further above/below
            (-far_offset, 0), (far_offset, 0),  # further left/right
        ]

        for dx, dy in offsets:
            candidate = (target_pos[0] + dx, target_pos[1] + dy)

            if not crosses_obstacle(candidate):
                aux.x, aux.y = candidate
                return

        # Fallback: keep current position (crossing unavoidable with simple repositioning)

    def _reroute_flow_around_stock(self, flow_name: str, stock_name: str):
        """Add waypoints to route flow around a stock it currently crosses.

        This is a best-effort fix - if the flow has already been modified multiple
        times (>8 points), we skip to avoid infinite loops.
        """
        flow = self.flows.get(flow_name)
        stock = self.stocks.get(stock_name)

        if not flow or not stock or not flow.points or len(flow.points) < 2:
            return
        if stock.x is None or stock.y is None:
            return

        # Guard against infinite rerouting - if flow already has many points, skip
        if len(flow.points) > 8:
            return

        box = BoundingBox(stock.x, stock.y, stock.width, stock.height)
        # Clearance proportional to stock size (minimum 20px, plus half the larger dimension)
        clearance = 20 + max(stock.width, stock.height) / 2

        # Find which segment intersects and modify
        for i in range(len(flow.points) - 1):
            p1, p2 = flow.points[i], flow.points[i + 1]

            if not segment_intersects_box(p1, p2, box):
                continue

            # Determine if this is a horizontal or vertical segment
            is_horizontal = abs(p1[1] - p2[1]) < abs(p1[0] - p2[0])

            if is_horizontal:
                # Route above or below the stock - pick the side further from current Y
                dist_above = abs(p1[1] - (stock.y - stock.height / 2))
                dist_below = abs(p1[1] - (stock.y + stock.height / 2))

                if dist_above > dist_below:
                    route_y = stock.y - stock.height / 2 - clearance
                else:
                    route_y = stock.y + stock.height / 2 + clearance

                # Insert waypoints to go around
                new_points = list(flow.points[:i + 1])
                new_points.append((p1[0], route_y))
                new_points.append((p2[0], route_y))
                new_points.extend(flow.points[i + 1:])
                flow.points = [(float(x), float(y)) for x, y in new_points]
            else:
                # Vertical segment - route left or right
                dist_left = abs(p1[0] - (stock.x - stock.width / 2))
                dist_right = abs(p1[0] - (stock.x + stock.width / 2))

                if dist_left > dist_right:
                    route_x = stock.x - stock.width / 2 - clearance
                else:
                    route_x = stock.x + stock.width / 2 + clearance

                new_points = list(flow.points[:i + 1])
                new_points.append((route_x, p1[1]))
                new_points.append((route_x, p2[1]))
                new_points.extend(flow.points[i + 1:])
                flow.points = [(float(x), float(y)) for x, y in new_points]

            return  # Only fix one intersection per call

    def _resolve_layout_violations(self, max_iterations: int = 10):
        """Iteratively resolve collisions and crossings in the layout."""
        # Track what we've already tried to fix to avoid infinite loops
        processed_flow_stock: set[tuple[str, str]] = set()
        processed_connector_flow: set[tuple[int, str]] = set()
        processed_connector_stock: set[tuple[int, str]] = set()

        for iteration in range(max_iterations):
            # Detect all violations
            aux_collisions = self._detect_aux_collisions()
            connector_flow_crossings = self._detect_connector_flow_crossings()
            connector_stock_crossings = self._detect_connector_stock_crossings()
            flow_stock_crossings = self._detect_flow_stock_crossings()

            # Filter out already-processed items
            new_connector_flow = [c for c in connector_flow_crossings if c not in processed_connector_flow]
            new_connector_stock = [c for c in connector_stock_crossings if c not in processed_connector_stock]
            new_flow_stock = [c for c in flow_stock_crossings if c not in processed_flow_stock]

            if not aux_collisions and not new_connector_flow and not new_connector_stock and not new_flow_stock:
                return  # Layout is valid (or as good as we can make it)

            # Resolve aux collisions first (simplest)
            for name1, name2 in aux_collisions:
                self._separate_auxs(name1, name2)

            # Resolve connector-flow crossings by repositioning auxs
            for conn_uid, flow_name in new_connector_flow:
                self._reposition_aux_to_avoid_crossing(conn_uid, flow_name, "flow")
                processed_connector_flow.add((conn_uid, flow_name))

            # Resolve connector-stock crossings by repositioning auxs
            for conn_uid, stock_name in new_connector_stock:
                self._reposition_aux_to_avoid_crossing(conn_uid, stock_name, "stock")
                processed_connector_stock.add((conn_uid, stock_name))

            # Resolve flow-stock crossings by rerouting flows
            for flow_name, stock_name in new_flow_stock:
                self._reroute_flow_around_stock(flow_name, stock_name)
                processed_flow_stock.add((flow_name, stock_name))

            # Recalculate connector angles after moving auxs
            self._calculate_connector_angles()

        # If we hit max iterations, layout is best-effort (some violations may remain)

    def to_xml(self) -> str:
        """Generate XMILE XML string for the model."""
        self._auto_layout()

        lines = []
        lines.append('<?xml version="1.0" encoding="utf-8"?>')
        lines.append(f'<xmile version="1.0" xmlns="{XMILE_NS}" xmlns:isee="{ISEE_NS}">')

        # Header
        lines.append('\t<header>')
        lines.append('\t\t<smile version="1.0" namespace="std, isee"/>')
        lines.append(f'\t\t<name>{escape(self.name)}</name>')
        lines.append(f'\t\t<uuid>{self.uuid}</uuid>')
        lines.append('\t\t<vendor>isee systems, inc.</vendor>')
        lines.append('\t\t<product version="1.9.3" isee:build_number="1954" isee:saved_by_v1="true" lang="en">Stella Professional</product>')
        lines.append('\t</header>')

        # Sim specs
        if self.sim_specs.dt < 1:
            dt_str = f'<dt reciprocal="true">{int(1/self.sim_specs.dt)}</dt>'
        else:
            dt_str = f'<dt>{self.sim_specs.dt}</dt>'
        lines.append(f'\t<sim_specs isee:sim_duration="1.5" isee:simulation_delay="0.0015" isee:restore_on_start="false" method="{self.sim_specs.method}" time_units="{self.sim_specs.time_units}" isee:instantaneous_flows="false">')
        lines.append(f'\t\t<start>{self.sim_specs.start}</start>')
        lines.append(f'\t\t<stop>{self.sim_specs.stop}</stop>')
        lines.append(f'\t\t{dt_str}')
        lines.append('\t</sim_specs>')

        # Preferences
        lines.append('\t<isee:prefs show_module_prefix="true" live_update_on_drag="true" show_restore_buttons="false" layer="model" interface_scale_ui="true" interface_max_page_width="10000" interface_max_page_height="10000" interface_min_page_width="0" interface_min_page_height="0" saved_runs="5" keep="false" rifp="true"/>')

        # Model
        lines.append('\t<model>')
        lines.append('\t\t<variables>')

        # Stocks
        for name, stock in self.stocks.items():
            display = escape(self._display_name(stock.name))
            lines.append(f'\t\t\t<stock name="{display}">')
            lines.append(f'\t\t\t\t<eqn>{escape(stock.initial_value)}</eqn>')
            for inflow in stock.inflows:
                lines.append(f'\t\t\t\t<inflow>{inflow}</inflow>')
            for outflow in stock.outflows:
                lines.append(f'\t\t\t\t<outflow>{outflow}</outflow>')
            if stock.non_negative:
                lines.append('\t\t\t\t<non_negative/>')
            if stock.units:
                lines.append(f'\t\t\t\t<units>{escape(stock.units)}</units>')
            lines.append('\t\t\t</stock>')

        # Flows
        for name, flow in self.flows.items():
            display = escape(self._display_name(flow.name))
            lines.append(f'\t\t\t<flow name="{display}">')
            lines.append(f'\t\t\t\t<eqn>{escape(flow.equation)}</eqn>')
            if flow.graphical_function is not None:
                self._add_graphical_function_str(lines, flow.graphical_function)
            if flow.non_negative:
                lines.append('\t\t\t\t<non_negative/>')
            if flow.units:
                lines.append(f'\t\t\t\t<units>{escape(flow.units)}</units>')
            lines.append('\t\t\t</flow>')

        # Auxiliaries
        for name, aux in self.auxs.items():
            display = escape(self._display_name(aux.name))
            lines.append(f'\t\t\t<aux name="{display}">')
            lines.append(f'\t\t\t\t<eqn>{escape(aux.equation)}</eqn>')
            if aux.graphical_function is not None:
                self._add_graphical_function_str(lines, aux.graphical_function)
            if aux.units:
                lines.append(f'\t\t\t\t<units>{escape(aux.units)}</units>')
            lines.append('\t\t\t</aux>')

        lines.append('\t\t</variables>')

        # Views
        lines.append('\t\t<views>')
        self._add_view_styles_str(lines)

        # Main view
        lines.append('\t\t\t<view isee:show_pages="false" background="white" page_width="768" page_height="596" isee:page_cols="2" isee:page_rows="2" isee:popup_graphs_are_comparative="true" type="stock_flow">')
        self._add_inner_view_styles_str(lines)

        # Stock visuals (positions guaranteed by _auto_layout)
        for name, stock in self.stocks.items():
            display = escape(self._display_name(stock.name))
            sx = int(stock.x) if stock.x is not None else 0
            sy = int(stock.y) if stock.y is not None else 0
            lines.append(f'\t\t\t\t<stock x="{sx}" y="{sy}" width="{stock.width}" height="{stock.height}" name="{display}"/>')

        # Flow visuals (positions guaranteed by _auto_layout)
        for name, flow in self.flows.items():
            display = escape(self._display_name(flow.name))
            fx = flow.x if flow.x is not None else 0
            fy = int(flow.y) if flow.y is not None else 0
            if flow.points:
                lines.append(f'\t\t\t\t<flow x="{fx}" y="{fy}" name="{display}">')
                lines.append('\t\t\t\t\t<pts>')
                for px, py in flow.points:
                    lines.append(f'\t\t\t\t\t\t<pt x="{px}" y="{py}"/>')
                lines.append('\t\t\t\t\t</pts>')
                lines.append('\t\t\t\t</flow>')
            else:
                lines.append(f'\t\t\t\t<flow x="{fx}" y="{fy}" name="{display}"/>')

        # Aux visuals (positions guaranteed by _auto_layout)
        for name, aux in self.auxs.items():
            display = escape(self._display_name(aux.name))
            ax = int(aux.x) if aux.x is not None else 0
            ay = int(aux.y) if aux.y is not None else 0
            lines.append(f'\t\t\t\t<aux x="{ax}" y="{ay}" name="{display}"/>')

        # Connector visuals
        for conn in self.connectors:
            lines.append(f'\t\t\t\t<connector uid="{conn.uid}" angle="{conn.angle}">')
            lines.append(f'\t\t\t\t\t<from>{conn.from_var}</from>')
            lines.append(f'\t\t\t\t\t<to>{conn.to_var}</to>')
            lines.append('\t\t\t\t</connector>')

        lines.append('\t\t\t</view>')
        lines.append('\t\t</views>')
        lines.append('\t</model>')
        lines.append('</xmile>')

        return '\n'.join(lines)

    def _add_view_styles_str(self, lines: list[str]):
        """Add the default view styles as strings."""
        lines.append('\t\t\t<style color="black" background="white" font_style="normal" font_weight="normal" text_decoration="none" text_align="center" vertical_text_align="center" font_color="black" font_family="Arial" font_size="10pt" padding="2" border_color="black" border_width="thin" border_style="none">')
        lines.append('\t\t\t\t<text_box color="black" background="white" text_align="left" vertical_text_align="top" font_size="12pt"/>')
        lines.append('\t\t\t</style>')

    def _add_inner_view_styles_str(self, lines: list[str]):
        """Add the inner view styles as strings."""
        lines.append('\t\t\t\t<style color="black" background="white" font_style="normal" font_weight="normal" text_decoration="none" text_align="center" vertical_text_align="center" font_color="black" font_family="Arial" font_size="10pt" padding="2" border_color="black" border_width="thin" border_style="none">')
        lines.append('\t\t\t\t\t<stock color="blue" background="white" font_color="blue" font_size="9pt" label_side="top">')
        lines.append('\t\t\t\t\t\t<shape type="rectangle" width="45" height="35"/>')
        lines.append('\t\t\t\t\t</stock>')
        lines.append('\t\t\t\t\t<flow color="blue" background="white" font_color="blue" font_size="9pt" label_side="bottom"/>')
        lines.append('\t\t\t\t\t<aux color="blue" background="white" font_color="blue" font_size="9pt" label_side="bottom">')
        lines.append('\t\t\t\t\t\t<shape type="circle" radius="18"/>')
        lines.append('\t\t\t\t\t</aux>')
        lines.append('\t\t\t\t\t<connector color="#FF007F" background="white" font_color="#FF007F" font_size="9pt" isee:thickness="1"/>')
        lines.append('\t\t\t\t</style>')

    def _format_point_list(self, points: list[float]) -> str:
        return " ".join(f"{p:g}" for p in points)

    def _add_graphical_function_str(self, lines: list[str], gf: GraphicalFunction):
        attrs = f' type="{escape(gf.gf_type)}"' if gf.gf_type else ""
        lines.append(f'\t\t\t\t<gf{attrs}>')
        if gf.xpts is not None:
            lines.append(f'\t\t\t\t\t<xpts>{self._format_point_list(gf.xpts)}</xpts>')
        elif gf.xscale is not None:
            lines.append(f'\t\t\t\t\t<xscale min="{gf.xscale[0]:g}" max="{gf.xscale[1]:g}"/>')
        if gf.yscale is not None:
            lines.append(f'\t\t\t\t\t<yscale min="{gf.yscale[0]:g}" max="{gf.yscale[1]:g}"/>')
        lines.append(f'\t\t\t\t\t<ypts>{self._format_point_list(gf.ypts)}</ypts>')
        lines.append('\t\t\t\t</gf>')


def parse_stmx(filepath: str) -> StellaModel:
    """Parse an existing .stmx file and return a StellaModel."""
    tree = ET.parse(filepath)
    root = tree.getroot()

    # Handle namespaces with full Clark notation
    xmile = f"{{{XMILE_NS}}}"
    isee = f"{{{ISEE_NS}}}"

    def find_elem(parent, *tags):
        """Find element trying both namespaced and non-namespaced tags."""
        for tag in tags:
            # Try with XMILE namespace
            elem = parent.find(f".//{xmile}{tag}")
            if elem is not None:
                return elem
            # Try without namespace
            elem = parent.find(f".//{tag}")
            if elem is not None:
                return elem
        return None

    def find_child(parent, tag):
        """Find direct child element."""
        elem = parent.find(f"{xmile}{tag}")
        if elem is None:
            elem = parent.find(tag)
        return elem

    def findall_children(parent, tag):
        """Find all direct children with given tag."""
        elems = parent.findall(f"{xmile}{tag}")
        if not elems:
            elems = parent.findall(tag)
        return elems

    def parse_point_list(text: Optional[str]) -> list[float]:
        if not text:
            return []
        return [float(val) for val in text.split()]

    def parse_gf(elem: Optional[ET.Element]) -> Optional[GraphicalFunction]:
        if elem is None:
            return None
        gf_type = elem.get("type")
        xpts_elem = find_child(elem, "xpts")
        xscale_elem = find_child(elem, "xscale")
        yscale_elem = find_child(elem, "yscale")
        ypts_elem = find_child(elem, "ypts")

        xpts = parse_point_list(xpts_elem.text) if xpts_elem is not None else None
        xscale = None
        if xscale_elem is not None:
            xscale = (
                float(xscale_elem.get("min", "0")),
                float(xscale_elem.get("max", "0")),
            )
        yscale = None
        if yscale_elem is not None:
            yscale = (
                float(yscale_elem.get("min", "0")),
                float(yscale_elem.get("max", "0")),
            )
        ypts = parse_point_list(ypts_elem.text if ypts_elem is not None else None)
        if not ypts:
            return None
        return GraphicalFunction(
            ypts=ypts,
            xscale=xscale,
            xpts=xpts,
            yscale=yscale,
            gf_type=gf_type if gf_type else None,
        )

    # Get model name
    header = find_child(root, "header")
    name_elem = find_child(header, "name") if header is not None else None
    model_name = name_elem.text if name_elem is not None else "Untitled"
    model = StellaModel(name=model_name)

    # Parse sim_specs
    sim_specs = find_child(root, "sim_specs")
    if sim_specs is not None:
        start = find_child(sim_specs, "start")
        if start is not None and start.text:
            model.sim_specs.start = float(start.text)

        stop = find_child(sim_specs, "stop")
        if stop is not None and stop.text:
            model.sim_specs.stop = float(stop.text)

        dt = find_child(sim_specs, "dt")
        if dt is not None and dt.text:
            if dt.get("reciprocal") == "true":
                model.sim_specs.dt = 1.0 / float(dt.text)
            else:
                model.sim_specs.dt = float(dt.text)

        method = sim_specs.get("method")
        if method:
            model.sim_specs.method = method

        time_units = sim_specs.get("time_units")
        if time_units:
            model.sim_specs.time_units = time_units

    # Find variables section
    model_elem = find_child(root, "model")
    variables = find_child(model_elem, "variables") if model_elem is not None else None

    if variables is not None:
        # Parse stocks
        for stock_elem in findall_children(variables, "stock"):
            name = stock_elem.get("name")
            eqn = find_child(stock_elem, "eqn")
            initial_value = eqn.text if eqn is not None else "0"

            units_elem = find_child(stock_elem, "units")
            units = units_elem.text if units_elem is not None else ""

            inflows = [inf.text for inf in findall_children(stock_elem, "inflow") if inf.text]
            outflows = [outf.text for outf in findall_children(stock_elem, "outflow") if outf.text]

            non_negative = find_child(stock_elem, "non_negative") is not None

            stock = Stock(
                name=name,
                initial_value=initial_value,
                units=units,
                inflows=inflows,
                outflows=outflows,
                non_negative=non_negative
            )
            model.stocks[model._normalize_name(name)] = stock

        # Parse flows
        for flow_elem in findall_children(variables, "flow"):
            name = flow_elem.get("name")
            eqn = find_child(flow_elem, "eqn")
            equation = eqn.text if eqn is not None else "0"
            gf = parse_gf(find_child(flow_elem, "gf"))

            units_elem = find_child(flow_elem, "units")
            units = units_elem.text if units_elem is not None else ""

            non_negative = find_child(flow_elem, "non_negative") is not None

            flow = Flow(
                name=name,
                equation=equation,
                units=units,
                non_negative=non_negative,
                graphical_function=gf,
            )
            model.flows[model._normalize_name(name)] = flow

        # Parse auxiliaries
        for aux_elem in findall_children(variables, "aux"):
            name = aux_elem.get("name")
            eqn = find_child(aux_elem, "eqn")
            equation = eqn.text if eqn is not None else "0"
            gf = parse_gf(find_child(aux_elem, "gf"))

            units_elem = find_child(aux_elem, "units")
            units = units_elem.text if units_elem is not None else ""

            aux = Aux(
                name=name,
                equation=equation,
                units=units,
                graphical_function=gf,
            )
            model.auxs[model._normalize_name(name)] = aux

    # Determine flow from/to stocks based on stock inflows/outflows
    for stock_name, stock in model.stocks.items():
        for inflow in stock.inflows:
            norm_inflow = model._normalize_name(inflow)
            if norm_inflow in model.flows:
                model.flows[norm_inflow].to_stock = stock_name
        for outflow in stock.outflows:
            norm_outflow = model._normalize_name(outflow)
            if norm_outflow in model.flows:
                model.flows[norm_outflow].from_stock = stock_name

    # Parse visual positions and connectors from views
    views = find_child(model_elem, "views") if model_elem is not None else None
    view = find_child(views, "view") if views is not None else None

    if view is not None:
        # Extract stock positions from view
        for stock_elem in findall_children(view, "stock"):
            name = stock_elem.get("name")
            x_attr = stock_elem.get("x")
            y_attr = stock_elem.get("y")
            if name:
                norm_name = model._normalize_name(name)
                if norm_name in model.stocks:
                    if x_attr is not None:
                        model.stocks[norm_name].x = float(x_attr)
                    if y_attr is not None:
                        model.stocks[norm_name].y = float(y_attr)

        # Extract flow positions from view
        for flow_elem in findall_children(view, "flow"):
            name = flow_elem.get("name")
            x_attr = flow_elem.get("x")
            y_attr = flow_elem.get("y")
            if name:
                norm_name = model._normalize_name(name)
                if norm_name in model.flows:
                    if x_attr is not None:
                        model.flows[norm_name].x = float(x_attr)
                    if y_attr is not None:
                        model.flows[norm_name].y = float(y_attr)

        # Extract aux positions from view
        for aux_elem in findall_children(view, "aux"):
            name = aux_elem.get("name")
            x_attr = aux_elem.get("x")
            y_attr = aux_elem.get("y")
            if name:
                norm_name = model._normalize_name(name)
                if norm_name in model.auxs:
                    if x_attr is not None:
                        model.auxs[norm_name].x = float(x_attr)
                    if y_attr is not None:
                        model.auxs[norm_name].y = float(y_attr)

        # Extract connectors
        for conn_elem in findall_children(view, "connector"):
            uid = int(conn_elem.get("uid", 0))
            angle = float(conn_elem.get("angle", 0))

            from_elem = find_child(conn_elem, "from")
            to_elem = find_child(conn_elem, "to")

            if from_elem is not None and to_elem is not None and from_elem.text and to_elem.text:
                connector = Connector(
                    uid=uid,
                    from_var=from_elem.text,
                    to_var=to_elem.text,
                    angle=angle
                )
                model.connectors.append(connector)
                model._connector_uid = max(model._connector_uid, uid)

    return model
