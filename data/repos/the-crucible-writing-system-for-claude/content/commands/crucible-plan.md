---
allowed-tools: Read, Write, Edit, Glob, Grep, Bash
argument-hint: [premise text] | "continue"
description: Start planning an epic fantasy novel using the Crucible Structure. Provide a premise/synopsis to begin, or use "continue" to resume an in-progress planning session.
---

# /crucible-plan

Start or continue planning an epic fantasy novel using the Crucible Structure.

## Execution Instructions

**IMPORTANT:** When this command is invoked, you MUST:

1. **Invoke the `crucible-suite:crucible-planner` skill** using the Skill tool
2. The skill will guide you through the complete planning workflow
3. Follow the skill's questioning protocol with multi-choice questions

## Usage

- `/crucible-suite:crucible-plan [your premise]` - Start new planning with a premise
- `/crucible-suite:crucible-plan continue` - Resume existing planning session

## What This Does

1. Activates the crucible-planner skill (you must follow its instructions)
2. Guides you through interactive questions for all three strands (Quest, Fire, Constellation)
3. Generates 7 comprehensive planning documents:
   - Crucible Thesis (philosophical core)
   - Strand Maps (3 separate arc maps)
   - Forge Point Blueprints (5 convergence crises)
   - Dark Mirror Profile (antagonist design)
   - Constellation Bible (character relationships)
   - Mercy Ledger (mercy/payoff tracking)
   - World Forge (world-building)
4. Saves state for resumable sessions

## Prerequisites

None - this is the first step in the Crucible workflow.

## Example

```
/crucible-suite:crucible-plan A young blacksmith discovers she can forge weapons that steal the memories of those they cut. When her village is destroyed by a memory-hunting cult, she must master her forbidden gift to save the last people who remember the old ways.
```

## What Happens Next

After planning completes, you can:
- `/crucible-suite:crucible-outline` - Create chapter outlines from your planning documents
- Review and adjust any planning document
- Add more detail to specific elements
