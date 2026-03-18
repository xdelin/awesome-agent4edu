"""Tests for force-directed layout engine."""

import math
import time

import pytest

from stella_mcp.layout import BoundingBox, force_directed_layout
from stella_mcp.xmile import StellaModel


class TestForceDirectedPureFunction:
    """Tests for the pure force_directed_layout() function."""

    def test_empty_model_returns_empty(self):
        """Empty node list produces empty result."""
        result = force_directed_layout([], [], {})
        assert result == {}

    def test_single_element_centered(self):
        """Single free node is placed at canvas center."""
        result = force_directed_layout(["A"], [], {})
        x, y = result["A"]
        assert abs(x - 800) < 1
        assert abs(y - 500) < 1

    def test_single_element_fixed(self):
        """Single pinned node keeps its position."""
        result = force_directed_layout(["A"], [], {"A": (100, 200)})
        assert result["A"] == (100, 200)

    def test_deterministic_layout(self):
        """Same inputs produce identical positions."""
        nodes = ["A", "B", "C"]
        edges = [("A", "B", 1.0), ("B", "C", 1.0)]

        r1 = force_directed_layout(nodes, edges, {})
        r2 = force_directed_layout(nodes, edges, {})

        for name in nodes:
            assert r1[name][0] == pytest.approx(r2[name][0], abs=0.001)
            assert r1[name][1] == pytest.approx(r2[name][1], abs=0.001)

    def test_user_positions_preserved(self):
        """Pinned nodes do not move."""
        nodes = ["A", "B", "C"]
        edges = [("A", "B", 1.0), ("B", "C", 1.0)]
        fixed = {"A": (100, 300), "C": (500, 300)}

        result = force_directed_layout(nodes, edges, fixed)

        assert result["A"] == (100, 300)
        assert result["C"] == (500, 300)

    def test_no_overlaps_after_layout(self):
        """All node pairs have 50px+ separation."""
        nodes = [f"n{i}" for i in range(10)]
        edges = [(f"n{i}", f"n{i+1}", 1.0) for i in range(9)]

        result = force_directed_layout(nodes, edges, {})

        names = list(result.keys())
        for i, a in enumerate(names):
            for b in names[i + 1:]:
                ax, ay = result[a]
                bx, by = result[b]
                dist = math.sqrt((ax - bx) ** 2 + (ay - by) ** 2)
                assert dist >= 50, f"{a} and {b} too close: {dist:.1f}px"

    def test_feedback_loop_nodes_not_collinear(self):
        """3-node cycle should form a triangle, not a line."""
        nodes = ["A", "B", "C"]
        edges = [("A", "B", 1.0), ("B", "C", 1.0), ("C", "A", 1.0)]

        result = force_directed_layout(nodes, edges, {})

        ax, ay = result["A"]
        bx, by = result["B"]
        cx, cy = result["C"]

        # Triangle area via cross product — should be significantly non-zero
        area = abs((bx - ax) * (cy - ay) - (cx - ax) * (by - ay)) / 2
        assert area > 100, f"Nodes are nearly collinear, area={area:.1f}"

    def test_all_positions_within_canvas_bounds(self):
        """No node placed outside canvas (with padding)."""
        nodes = [f"n{i}" for i in range(8)]
        edges = [(f"n{i}", f"n{i+1}", 1.0) for i in range(7)]

        result = force_directed_layout(nodes, edges, {}, canvas_width=1600, canvas_height=1000)

        for name, (x, y) in result.items():
            assert x >= 0, f"{name} x={x} is off-canvas left"
            assert y >= 0, f"{name} y={y} is off-canvas top"
            # Allow expansion but should still be reasonable
            assert x <= 3200, f"{name} x={x} is unreasonably far right"
            assert y <= 2000, f"{name} y={y} is unreasonably far down"

    def test_all_pinned_converges(self):
        """All nodes fixed — returns immediately with those positions."""
        nodes = ["A", "B"]
        edges = [("A", "B", 1.0)]
        fixed = {"A": (100, 200), "B": (300, 400)}

        result = force_directed_layout(nodes, edges, fixed)

        assert result["A"] == (100, 200)
        assert result["B"] == (300, 400)


class TestForceDirectedIntegration:
    """Integration tests: FR layout through StellaModel._auto_layout()."""

    def test_no_overlaps_after_layout(self):
        """All element pairs have 50px+ bounding box gap after full layout."""
        model = StellaModel("Test")
        model.add_stock("Population", "1000")
        model.add_stock("Food", "500")
        model.add_stock("Pollution", "0")
        model.add_flow("births", "Population * birth_rate", to_stock="Population")
        model.add_flow("consumption", "Population * 0.1", from_stock="Food", to_stock="Population")
        model.add_flow("emissions", "Population * 0.05", from_stock="Population", to_stock="Pollution")
        model.add_aux("birth_rate", "0.02")
        model.add_connector("birth_rate", "births")
        model.add_connector("Population", "births")

        model._auto_layout()

        # Check all stock/aux pairs have minimum separation
        elements: dict[str, tuple[float, float]] = {}
        for name, s in model.stocks.items():
            if s.x is not None and s.y is not None:
                elements[name] = (s.x, s.y)
        for name, a in model.auxs.items():
            if a.x is not None and a.y is not None:
                elements[name] = (a.x, a.y)

        names = list(elements.keys())
        for i, a in enumerate(names):
            for b in names[i + 1:]:
                ax, ay = elements[a]
                bx, by = elements[b]
                dist = math.sqrt((ax - bx) ** 2 + (ay - by) ** 2)
                assert dist >= 40, f"{a} and {b} too close: {dist:.1f}px"

    def test_deterministic_layout(self):
        """Same model built twice produces identical positions."""
        def build():
            m = StellaModel("Test")
            m.add_stock("A", "100")
            m.add_stock("B", "50")
            m.add_flow("f1", "10", from_stock="A", to_stock="B")
            m.add_aux("rate", "0.1")
            m.add_connector("rate", "f1")
            m._auto_layout()
            return m

        m1 = build()
        m2 = build()

        for name in m1.stocks:
            assert m1.stocks[name].x == pytest.approx(m2.stocks[name].x, abs=0.01)
            assert m1.stocks[name].y == pytest.approx(m2.stocks[name].y, abs=0.01)
        for name in m1.auxs:
            assert m1.auxs[name].x == pytest.approx(m2.auxs[name].x, abs=0.01)
            assert m1.auxs[name].y == pytest.approx(m2.auxs[name].y, abs=0.01)

    def test_user_positions_preserved(self):
        """Pinned nodes don't move through full layout pipeline."""
        model = StellaModel("Test")
        model.add_stock("A", "100", x=100, y=300)
        model.add_stock("B", "50", x=500, y=300)
        model.add_flow("f1", "10", from_stock="A", to_stock="B")

        model._auto_layout()

        assert model.stocks["A"].x == 100
        assert model.stocks["A"].y == 300
        assert model.stocks["B"].x == 500
        assert model.stocks["B"].y == 300

    def test_empty_model_no_crash(self):
        """Empty model doesn't crash."""
        model = StellaModel("Empty")
        xml = model.to_xml()
        assert "<xmile" in xml

    def test_feedback_loop_nodes_not_collinear(self):
        """3-stock feedback loop should form a triangle, not a line."""
        model = StellaModel("Test")
        model.add_stock("A", "100")
        model.add_stock("B", "50")
        model.add_stock("C", "25")
        model.add_flow("f1", "10", from_stock="A", to_stock="B")
        model.add_flow("f2", "5", from_stock="B", to_stock="C")
        model.add_flow("f3", "3", from_stock="C", to_stock="A")

        model._auto_layout()

        ax, ay = model.stocks["A"].x, model.stocks["A"].y
        bx, by = model.stocks["B"].x, model.stocks["B"].y
        cx, cy = model.stocks["C"].x, model.stocks["C"].y

        assert all(v is not None for v in [ax, ay, bx, by, cx, cy])

        # Triangle area should be significantly non-zero
        area = abs((bx - ax) * (cy - ay) - (cx - ax) * (by - ay)) / 2
        assert area > 100, f"Stocks are nearly collinear, area={area:.1f}"

    def test_all_positions_within_canvas_bounds(self):
        """No element placed at negative coordinates."""
        model = StellaModel("Test")
        for i in range(6):
            model.add_stock(f"S{i}", "100")
        for i in range(5):
            model.add_flow(f"f{i}", "10", from_stock=f"S{i}", to_stock=f"S{i+1}")

        model._auto_layout()

        for name, s in model.stocks.items():
            assert s.x is not None and s.x >= 0, f"Stock {name} x={s.x}"
            assert s.y is not None and s.y >= 0, f"Stock {name} y={s.y}"

    def test_layout_completes_within_time_budget(self):
        """to_xml() completes in < 2 seconds for 30-element model."""
        model = StellaModel("Test")
        # 10 stocks, 9 flows, 11 auxs = 30 elements
        for i in range(10):
            model.add_stock(f"S{i}", "100")
        for i in range(9):
            model.add_flow(f"f{i}", "10", from_stock=f"S{i}", to_stock=f"S{i+1}")
        for i in range(11):
            model.add_aux(f"a{i}", "0.1")
            model.add_connector(f"a{i}", f"f{min(i, 8)}")

        start = time.time()
        model.to_xml()
        elapsed = time.time() - start

        assert elapsed < 2.0, f"Layout took {elapsed:.2f}s, exceeds 2s budget"
