"""Geometry primitives and force-directed layout for Stella model diagrams.

Phase 2 TODO (visual tuning after Stella Professional testing):
- Re-enable _resolve_layout_violations checks if connector/flow lines cross unrelated stocks
- Add L-R post-processing if linear chains need left-to-right ordering
- Further tune edge weights (FLOW_WEIGHT, CONNECTOR_WEIGHT) based on real models
- Profile to_xml() for larger (30+ element) models
"""

import math
from collections.abc import Mapping
from dataclasses import dataclass


# =============================================================================
# Geometry Primitives
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
    """Check if line segment p1-p2 intersects segment p3-p4."""
    return (_ccw(p1, p3, p4) != _ccw(p2, p3, p4) and
            _ccw(p1, p2, p3) != _ccw(p1, p2, p4))


def segment_intersects_box(
    p1: tuple[float, float], p2: tuple[float, float],
    box: BoundingBox
) -> bool:
    """Check if line segment p1-p2 intersects a bounding box."""
    left = box.x - box.width / 2
    right = box.x + box.width / 2
    top = box.y - box.height / 2
    bottom = box.y + box.height / 2

    # Segment entirely inside box
    if (left <= p1[0] <= right and top <= p1[1] <= bottom and
        left <= p2[0] <= right and top <= p2[1] <= bottom):
        return True

    edges = [
        ((left, top), (right, top)),
        ((right, top), (right, bottom)),
        ((right, bottom), (left, bottom)),
        ((left, bottom), (left, top)),
    ]

    for edge_p1, edge_p2 in edges:
        if segments_intersect(p1, p2, edge_p1, edge_p2):
            return True

    return False


# =============================================================================
# Fruchterman-Reingold Force-Directed Layout
# =============================================================================

def force_directed_layout(
    nodes: list[str],
    edges: list[tuple[str, str, float]],
    fixed_positions: Mapping[str, tuple[float, float]],
    canvas_width: float = 1600.0,
    canvas_height: float = 1000.0,
    iterations: int = 300,
    min_separation: float = 50.0,
) -> dict[str, tuple[float, float]]:
    """Compute node positions using Fruchterman-Reingold force-directed layout.

    Pure function. Deterministic (sorted iteration + circular initial placement).

    Args:
        nodes: List of node names.
        edges: List of (source, target, weight) tuples.
        fixed_positions: Nodes with user-specified positions (pinned, exert forces but don't move).
        canvas_width: Canvas width in pixels.
        canvas_height: Canvas height in pixels.
        iterations: Number of simulation steps.
        min_separation: Minimum distance between node centers after layout.

    Returns:
        Dict mapping node name to (x, y) position.
    """
    if not nodes:
        return {}

    n = len(nodes)

    if n == 1:
        name = nodes[0]
        if name in fixed_positions:
            return {name: fixed_positions[name]}
        return {name: (canvas_width / 2, canvas_height / 2)}

    # Ideal edge length
    area = canvas_width * canvas_height
    k = math.sqrt(area / n)

    # Initial placement: fixed nodes at their positions, free nodes on a circle
    pos: dict[str, tuple[float, float]] = {}
    cx, cy = canvas_width / 2, canvas_height / 2
    radius = min(canvas_width, canvas_height) * 0.35

    free_nodes = sorted([nd for nd in nodes if nd not in fixed_positions])
    for name in fixed_positions:
        if name in nodes:
            pos[name] = fixed_positions[name]

    for i, name in enumerate(free_nodes):
        angle = 2 * math.pi * i / max(len(free_nodes), 1)
        pos[name] = (cx + radius * math.cos(angle), cy + radius * math.sin(angle))

    # Build adjacency for attractive forces
    adj: dict[str, list[tuple[str, float]]] = {nd: [] for nd in nodes}
    for src, tgt, w in edges:
        if src in adj and tgt in adj:
            adj[src].append((tgt, w))
            adj[tgt].append((src, w))

    # Simulation with linear cooling
    t = max(canvas_width, canvas_height) * 0.1  # initial temperature
    dt = t / (iterations + 1)

    for _ in range(iterations):
        disp: dict[str, tuple[float, float]] = {nd: (0.0, 0.0) for nd in nodes}

        # Repulsive forces between all pairs
        sorted_nodes = sorted(nodes)
        for i, u in enumerate(sorted_nodes):
            ux, uy = pos[u]
            dx_acc, dy_acc = 0.0, 0.0
            for v in sorted_nodes[i + 1:]:
                vx, vy = pos[v]
                dx = ux - vx
                dy = uy - vy
                dist = math.sqrt(dx * dx + dy * dy)
                if dist < 0.01:
                    dist = 0.01
                force = (k * k) / dist
                fx = (dx / dist) * force
                fy = (dy / dist) * force
                dx_acc += fx
                dy_acc += fy
                disp[v] = (disp[v][0] - fx, disp[v][1] - fy)
            disp[u] = (disp[u][0] + dx_acc, disp[u][1] + dy_acc)

        # Attractive forces along edges
        for src, tgt, w in edges:
            if src not in pos or tgt not in pos:
                continue
            sx, sy = pos[src]
            tx, ty = pos[tgt]
            dx = sx - tx
            dy = sy - ty
            dist = math.sqrt(dx * dx + dy * dy)
            if dist < 0.01:
                continue
            force = w * (dist * dist) / k
            fx = (dx / dist) * force
            fy = (dy / dist) * force
            disp[src] = (disp[src][0] - fx, disp[src][1] - fy)
            disp[tgt] = (disp[tgt][0] + fx, disp[tgt][1] + fy)

        # Apply displacements (skip fixed nodes), clamp by temperature
        for name in sorted(nodes):
            if name in fixed_positions:
                continue
            dx, dy = disp[name]
            dist = math.sqrt(dx * dx + dy * dy)
            if dist < 0.01:
                continue
            scale = min(dist, t) / dist
            nx = pos[name][0] + dx * scale
            ny = pos[name][1] + dy * scale
            pos[name] = (nx, ny)

        t -= dt

    # Post-processing: enforce minimum separation
    sorted_nodes = sorted(nodes)
    for _ in range(50):
        moved = False
        for i, u in enumerate(sorted_nodes):
            if u in fixed_positions:
                continue
            ux, uy = pos[u]
            for v in sorted_nodes[i + 1:]:
                vx, vy = pos[v]
                dx = ux - vx
                dy = uy - vy
                dist = math.sqrt(dx * dx + dy * dy)
                if dist < min_separation:
                    if dist < 0.01:
                        dx, dy, dist = 1.0, 0.0, 1.0
                    push = (min_separation - dist) / 2 + 1
                    nx = (dx / dist) * push
                    ny = (dy / dist) * push
                    if u not in fixed_positions:
                        pos[u] = (ux + nx, uy + ny)
                    if v not in fixed_positions:
                        pos[v] = (vx - nx, vy - ny)
                    moved = True
        if not moved:
            break

    # Rescale free nodes to fit within canvas, expanding only if needed
    padding = 50.0
    free = [nd for nd in nodes if nd not in fixed_positions]
    if not free:
        return pos

    free_x = [pos[nd][0] for nd in free]
    free_y = [pos[nd][1] for nd in free]
    min_x, max_x = min(free_x), max(free_x)
    min_y, max_y = min(free_y), max(free_y)
    span_x = max_x - min_x
    span_y = max_y - min_y

    # Available space (use canvas, expand if layout truly needs more)
    avail_w = canvas_width - 2 * padding
    avail_h = canvas_height - 2 * padding

    # Scale down if layout exceeds canvas, but don't upscale small layouts
    if span_x > avail_w:
        scale_x = avail_w / span_x
    else:
        scale_x = 1.0
    if span_y > avail_h:
        scale_y = avail_h / span_y
    else:
        scale_y = 1.0
    scale = min(scale_x, scale_y)

    # Center of current free-node positions
    center_x = (min_x + max_x) / 2
    center_y = (min_y + max_y) / 2

    # Target center
    target_cx = canvas_width / 2
    target_cy = canvas_height / 2

    for name in free:
        x, y = pos[name]
        x = (x - center_x) * scale + target_cx
        y = (y - center_y) * scale + target_cy
        # Hard clamp to padding
        x = max(padding, x)
        y = max(padding, y)
        pos[name] = (x, y)

    return pos
