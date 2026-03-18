# Example Project: The Memory Forge

This folder contains a **complete example** of a novel developed using the Crucible Writing System. It demonstrates the full workflow from initial premise through planning, outlining, drafting, and editing.

## About This Example

**Title:** The Memory Forge
**Genre:** Epic Fantasy
**Target Length:** 100,000 words
**Current Progress:** 41,950 words (10 chapters, 12 scenes)

**Premise:** A protagonist with memory-forging powers must protect the last memory-keepers while resisting the transformation her gift inflicts, ultimately discovering that identity is choice, not origin.

---

## Folder Structure

```
Example Project/
|
|-- CLAUDE.md              # Root memory document (project state + commands)
|-- style-profile.md       # Author voice and style guidelines
|-- readme.md              # This file
|
|-- planning/              # Phase 1: Foundation documents
|   |-- CLAUDE.md
|   |-- crucible-thesis.md         # Core story argument and structure
|   |-- world-forge.md             # Worldbuilding reference
|   |-- constellation-bible.md     # Character profiles and arcs
|   |-- dark-mirror-profile.md     # Antagonist development
|   |-- mercy-ledger.md            # Cost/consequence tracking
|   |
|   |-- strand-maps/               # The three narrative strands
|   |   |-- quest-strand-map.md    # External plot (Burden)
|   |   |-- fire-strand-map.md     # Internal transformation
|   |   |-- constellation-strand-map.md  # Relationship dynamics
|   |
|   |-- forge-points/              # Major turning points
|       |-- fp0-ignition.md        # Story catalyst
|       |-- fp1-first-forge-point.md
|       |-- fp2-second-forge-point.md
|       |-- fp3-third-forge-point.md
|       |-- apex.md                # Climax structure
|
|-- outline/               # Phase 2: Scene-by-scene breakdown
|   |-- master-outline.md          # Full chapter/scene overview
|   |-- foreshadowing-tracker.md   # Setup/payoff management
|   |
|   |-- chapters/                  # Individual chapter outlines
|       |-- chapter-01.md through chapter-10.md
|
|-- draft/                 # Phase 3: Prose output
    |-- edit-report-ch1-3.md       # Editing feedback (chapters 1-3)
    |-- edit-report-complete.md    # Full manuscript review
    |
    |-- chapters/                  # Written chapters
        |-- chapter-01.md through chapter-10.md
```

---

## How to Use This Example

### As a Learning Reference

1. **Start with `planning/crucible-thesis.md`** - See how a premise becomes a structured story argument
2. **Review the strand maps** - Understand how three narrative threads interweave
3. **Compare outlines to drafts** - See how scene beats translate to prose
4. **Check edit reports** - Learn what Crucible's review process catches

### As a Template

Copy this folder structure to start your own project:

```bash
cp -r "Example Project" "My Novel"
```

Then clear the content files and run `/crucible-plan` with your premise.

### Key Files to Study

| File | What It Demonstrates |
|------|---------------------|
| `CLAUDE.md` (root) | How to track project state across sessions |
| `style-profile.md` | Defining author voice for consistent prose |
| `crucible-thesis.md` | The 36-beat structure in action |
| `forge-points/` | Major turning points with emotional beats |
| `edit-report-complete.md` | The review/revision workflow |

---

## The Story at a Glance

**Protagonist:** Kira, a memory-forger whose gift threatens to consume her humanity

**Core Bond:** Tam, a childhood friend who anchors her to her former self

**Antagonist:** Lord Sethren, a cult leader who shares her power but chose domination over protection. He views memories as chains; she views them as sacred.

**Central Conflict:** Each use of memory-forging makes Kira more powerful but less human. She must decide what she's willing to sacrifice to protect the last memory-keepers.

**Resolution:** Kira surrenders her gift entirely to defeat the cult, becoming a guardian rather than a wielder. Identity is defined by choice, not ability.

---

## The Three Strands

The Crucible Structure weaves three narrative threads:

1. **Quest Strand (Burden)** - Protecting the living memory-keepers; the external plot weight Kira carries

2. **Fire Strand (Transformation)** - Each use of power risks her humanity; the internal cost of her journey

3. **Constellation Strand (Bonds)** - Tam and other relationships that anchor her; what she's fighting to preserve

These strands converge at four **Forge Points** where the story's direction fundamentally shifts.

---

## Workflow Phases Demonstrated

| Phase | Commands Used | Output |
|-------|--------------|--------|
| Planning | `/crucible-plan` | 7 foundation documents |
| Outlining | `/crucible-outline` | Chapter breakdowns + trackers |
| Writing | `/crucible-write` | Scene-by-scene prose |
| Editing | `/crucible-edit` | Review reports + revisions |

---

## Notes

- The `CLAUDE.md` files use a **lazy-loading architecture** - subdirectory files only load when you work in that phase, preserving context window
- The `style-profile.md` was derived from sample prose to ensure consistent voice
- Edit reports show the bi-chapter review process that catches continuity and craft issues

For full documentation on the Crucible Writing System, see the [main README](../README.md).
