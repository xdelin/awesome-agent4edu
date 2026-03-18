---
allowed-tools: Read, Write, Edit, Glob, Grep, Bash
argument-hint: [book number] | "continue"
description: Create chapter-by-chapter outlines from Crucible planning documents. Requires completed planning phase.
---

# /crucible-outline

Generate detailed chapter outlines from your Crucible planning documents.

## Execution Instructions

**IMPORTANT:** When this command is invoked, you MUST:

1. **Invoke the `crucible-suite:crucible-outliner` skill** using the Skill tool
2. The skill will guide you through the complete outlining workflow
3. Verify planning documents are complete before outlining
4. Follow the skill's beat-to-chapter mapping process

## Usage

- `/crucible-suite:crucible-outline 1` - Outline Book 1
- `/crucible-suite:crucible-outline 2` - Outline Book 2 (for series)
- `/crucible-suite:crucible-outline continue` - Resume outlining session

## What This Does

1. Activates the crucible-outliner skill (you must follow its instructions)
2. Loads your Crucible planning documents
3. Maps the 36 beats to chapters
4. Guides you through scene-by-scene outlining
5. Tracks foreshadowing plants and payoffs
6. Generates complete outline package:
   - Master outline (full book structure)
   - Chapter summaries (quick reference)
   - Scene breakdown (detailed scene list)
   - Foreshadowing tracker (plant/payoff tracking)
   - Individual chapter outlines

## Prerequisites

Requires completed planning documents from `/crucible-suite:crucible-plan`:
- Crucible Thesis
- Strand Maps (Quest, Fire, Constellation)
- Forge Point Blueprints
- Dark Mirror Profile
- Constellation Bible
- Mercy Ledger
- World Forge

## The 36-Beat Structure

Your outline will follow the Crucible's five movements:

| Movement | Beats | Approx. Chapters |
|----------|-------|------------------|
| I. Ignition | 1-6 | 1-5 |
| II. First Tempering | 7-14 | 6-10 |
| III. Scattering | 15-22 | 11-16 |
| IV. Brightest Burning | 23-28 | 17-21 |
| V. Final Forging | 29-34 | 22-24 |
| Coda | 35-36 | 25 |

## What Happens Next

After outlining completes, you can:
- `/crucible-suite:crucible-write` - Begin drafting prose from your outline
- Adjust specific chapter outlines
- Add more scene detail
- Outline additional books in a series
