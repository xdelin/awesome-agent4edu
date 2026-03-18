---
allowed-tools: Read, Write, Edit, Glob, Grep, Bash
argument-hint: [chapter number] | "all"
description: Revision and editing workflows for Crucible manuscripts. Handles developmental, line, copy editing, and polish.
---

# /crucible-edit

Revise and edit your Crucible manuscript at multiple levels.

## Execution Instructions

**IMPORTANT:** When this command is invoked, you MUST:

1. **Invoke the `crucible-suite:crucible-editor` skill** using the Skill tool
2. The skill will guide you through the complete editing workflow
3. Follow the skill's assessment and multi-level editing protocol
4. If bi-chapter review reports exist, prioritize the issues they flagged

## Usage

- `/crucible-suite:crucible-edit 1` - Edit chapter 1
- `/crucible-suite:crucible-edit 5-10` - Edit chapters 5 through 10
- `/crucible-suite:crucible-edit all` - Full manuscript editing pass

## What This Does

1. Activates the crucible-editor skill (you must follow its instructions)
2. Assesses the draft to determine editing needs
3. Works through four editing levels (as needed):
   - **Developmental** - Structure, pacing, character arcs
   - **Line editing** - Prose quality, voice consistency
   - **Copy editing** - Grammar, consistency, style
   - **Polish** - Final refinement, word choice, rhythm
4. Generates detailed change reports
5. Preserves original versions for comparison

## Editing Levels

### Developmental Editing
- Structural issues (scene order, chapter breaks)
- Pacing problems (sections that drag or rush)
- Character arc verification
- Plot coherence checks

### Line Editing
- Show vs tell improvements
- Filter word removal
- Weak verb strengthening
- Dialogue refinement
- Sensory detail balance

### Copy Editing
- Grammar and punctuation
- Tense and POV consistency
- Name/term spelling consistency
- Timeline reference accuracy

### Polish
- Opening/closing sentence impact
- Sentence rhythm refinement
- Word choice precision
- Final read-aloud quality

## Prerequisites

Requires draft chapters from `/crucible-suite:crucible-write` or imported manuscript.

Optional but recommended:
- Original chapter outlines (for developmental editing)
- Style profile (for voice consistency)
- Story bible (for continuity checking)
- Bi-chapter review reports (prioritizes flagged issues)

## What Happens Next

After editing completes:
- Review the diff report showing all changes
- Compile the final edited manuscript
- Return to specific chapters for further refinement
