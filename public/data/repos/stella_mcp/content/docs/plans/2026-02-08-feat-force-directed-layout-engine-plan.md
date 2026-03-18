---
title: "feat: Force-directed layout engine"
type: feat
date: 2026-02-08
---

# Force-Directed Layout Engine

## Overview

Replace the rule-based layout in `xmile.py` with Fruchterman-Reingold force-directed layout. The current system places elements in fixed bands (stocks at Y=300, auxs at Y=150) and nudges collisions apart. This fails for feedback loops, circular topologies, and dense models. FR uses physics simulation — repulsive forces between all nodes, attractive forces along edges — to produce naturally-spaced, readable diagrams. Pure Python, zero new dependencies, ~80 lines for the core algorithm.

## Key Decisions

| Decision | Choice | Rationale |
|----------|--------|-----------|
| Algorithm | Fruchterman-Reingold | Best for small graphs with feedback loops |
| Dependencies | None (pure Python) | Only `mcp` is currently required |
| Canvas | Fixed 1600x1000 | Close to XMILE page geometry (1536x1192); adjust if needed |
| Flow positioning | Derived from stocks (midpoint) | Flows sit on the pipe, not in simulation |
| User positions | Pinned as fixed nodes | Exert forces but don't move |
| Determinism | Sorted iteration + circular initial placement | Same model = same layout every time |
| L-R ordering | Not guaranteed | FR does not enforce left-to-right for linear chains; evaluate in Phase 2 |

## Architecture

Minimal extraction — only the pure algorithm and geometry move to a new file.

**New file `stella_mcp/layout.py` contains:**
- Geometry primitives (`BoundingBox`, `_ccw`, `segments_intersect`, `segment_intersects_box`)
- `force_directed_layout()` — pure function, no knowledge of `StellaModel`

**Stays on `StellaModel` in `xmile.py`:**
- `_auto_layout()` — existing entry point, calls pipeline methods as before
- `_position_subsystem()` — modified to call `layout.force_directed_layout()` internally
- `_build_dependency_graph()`, `_find_subsystems()`, `_arrange_subsystems()` — unchanged
- `_recalculate_flow_points()`, `_calculate_connector_angles()` — unchanged initially; flow routing may need direction-aware fixes if FR positions break assumptions (see Known Risks)
- Data model classes, XML serialization/parsing

**Delete (dead code after FR):**
- `_find_stock_chains` — chain-based seeding biases against cycles
- `_find_aux_target_position` — aux positioning handled by FR simulation
- `_calculate_stock_spacing` — FR uses forces, not fixed spacing

**Defer (do not delete yet):** `_resolve_layout_violations` and its helpers. FR prevents node-node overlaps but does not prevent connector/flow lines crossing through unrelated stocks. We accept crossing regressions in Phase 1. In Phase 2, if Stella testing reveals crossings, re-enable specific checks (e.g., `_reroute_flow_around_stock`). If FR handles everything, delete in Phase 2.

## Core Algorithm

```python
def force_directed_layout(
    nodes: list[str],
    edges: list[tuple[str, str, float]],  # (source, target, weight)
    fixed_positions: dict[str, tuple[float, float]],
    canvas_width: float = 1600.0,
    canvas_height: float = 1000.0,
) -> dict[str, tuple[float, float]]:
```

Pure function. Deterministic (sorted iteration). ~80 lines.

**Force equation:** Attractive force along edges = `weight * (d² / k)` where `k = sqrt(area / n)`. Repulsive force between all pairs = `k² / d`. Standard FR with linear cooling.

**Edge weights (initial values, tune in Phase 2):**
- Flow-stock edges: weight 2.0 (strong — stocks sharing flows should be close)
- Connector edges: weight 0.5 (weaker — auxs nearby but not tightly bound)

**Initial placement:** Deterministic circle around canvas center, sorted by name. Fixed nodes use their specified positions.

**Edge cases:** Empty model → empty dict. Single element → canvas center. All pinned → converges immediately.

## Known Risks

1. **Flow routing:** `_recalculate_flow_points()` assumes roughly horizontal stock-to-stock layout (flow exit points use `from_half_width`, implying left-to-right). FR places stocks at arbitrary angles. Attachment points may need to become direction-aware — exit from the edge closest to the destination. Will surface during Phase 1 testing.

2. **Crossing regressions:** Without `_resolve_layout_violations`, connector and flow lines may cross through unrelated stocks. Accepted for Phase 1; fix targeted checks in Phase 2 if needed.

3. **L-R ordering:** Linear chains (A→B→C) may not render left-to-right. Evaluate in Phase 2 whether a post-processing step is needed.

## Implementation

### Phase 1: Core Algorithm + Tests

**New file: `stella_mcp/layout.py`**
- Geometry primitives extracted from `xmile.py`
- `force_directed_layout()` pure function

**Modified: `stella_mcp/xmile.py`**
- `_position_subsystem()` calls `layout.force_directed_layout()` instead of chain-based placement
- Delete `_find_stock_chains`, `_find_aux_target_position`, `_calculate_stock_spacing`
- Skip `_resolve_layout_violations` call in `_auto_layout()` (keep the code for now)
- Update geometry imports to use `layout` module

**New tests (write first):**
- `test_no_overlaps_after_layout` — all element pairs have 50px+ bounding box gap
- `test_deterministic_layout` — same model twice → identical positions
- `test_user_positions_preserved` — pinned nodes don't move
- `test_empty_model_no_crash`
- `test_single_element_centered`
- `test_feedback_loop_nodes_not_collinear` — 3-node cycle forms a triangle, not a line
- `test_all_positions_within_canvas_bounds` — no node placed outside canvas
- `test_layout_completes_within_time_budget` — `to_xml()` < 2s for 30-element model

**Existing tests to update** (absolute → range/distance assertions):
- `test_stock_chain_horizontal` (line 243): stocks at similar Y within 100px; remove L-R assertions
- `test_stock_chain_follows_flow_direction` (line 310): same treatment — relax L-R assertions
- `test_aux_near_flow_via_connector` (line 171): aux within 200px of connector target
- `test_auxs_dont_overlap_after_layout` (line 655): should still pass with FR alone
- `test_dense_aux_network_no_overlaps` (line 737): should still pass with FR alone
- `test_feedback_loop_doesnt_cross_itself` (line 717): may need relaxed crossing assertion without `_resolve_layout_violations`
- Update geometry imports (`xmile` → `layout`)

### Phase 2: Tune and Validate

Open generated `.stmx` files in Stella Professional. If flow routing breaks, make attachment points direction-aware. If crossings are unacceptable, re-enable targeted checks from `_resolve_layout_violations`. If L-R ordering matters, add a post-processing step. Tune edge weights. Profile `to_xml()` for 5/15/30 element models.

## Acceptance Criteria

- [x] Models with 5-15 elements open in Stella with no overlaps
- [x] Feedback loops render as circles or grids, not tangles
- [x] 50px+ breathing room between all element bounding boxes
- [x] User-specified positions are always preserved
- [x] Deterministic: same model → same positions
- [x] All positions within canvas bounds
- [x] `to_xml()` < 2 seconds for 30-element models
- [x] All existing tests pass (with updated assertions)
- [x] Zero new dependencies

## References

- [Fruchterman-Reingold Paper](https://doi.org/10.1002/spe.4380211102)
- [NetworkX spring_layout](https://networkx.org/documentation/stable/reference/generated/networkx.drawing.layout.spring_layout.html)
