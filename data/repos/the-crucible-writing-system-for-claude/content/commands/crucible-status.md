---
allowed-tools: Read, Glob, Grep, Bash
argument-hint:
description: Show current project progress across all Crucible phases (planning, outlining, writing, editing).
model: claude-haiku-4-5-20251001
---

# /crucible-status

Display comprehensive status of your Crucible writing project.

## Execution Instructions

**IMPORTANT:** When this command is invoked, you MUST:

1. **Run the status reporter script** to gather project data:
```bash
python "${CLAUDE_PLUGIN_ROOT}/scripts/status_reporter.py" "." "text"
```

2. If the script is not available or fails, manually gather status by reading:
   - `.crucible/state/planning-state.json`
   - `.crucible/state/outline-state.json`
   - `.crucible/state/draft-state.json`
   - `.crucible/state/edit-state.json`

3. **Present the formatted output** to the user

## Usage

- `/crucible-suite:crucible-status` - Show full project status

## What This Shows

### Project Overview
- Project title and target word count
- Current phase (planning/outlining/writing/editing)
- Overall completion percentage

### Planning Status
- Which planning documents are complete
- Missing or incomplete sections
- Last modified dates

### Outline Status
- Chapters outlined vs. total
- Beat coverage verification
- Foreshadowing thread tracking

### Writing Status
- Chapters written vs. outlined
- Current word count vs. target
- Scenes completed in current chapter
- Next bi-chapter review trigger

### Editing Status
- Chapters through each editing level
- Issues found and resolved
- Word count changes

### Recent Activity
- Last 5 actions taken
- Files recently modified
- Backup status

## Example Output

```
CRUCIBLE PROJECT STATUS
═══════════════════════════════════════

📚 The Memory Forge
   Phase: WRITING
   Progress: 42% complete

PLANNING    ████████████████████ 100%
├─ Thesis       ✓ complete
├─ Strand Maps  ✓ complete (3/3)
├─ Forge Points ✓ complete (5/5)
├─ Dark Mirror  ✓ complete
├─ Constellation ✓ complete
├─ Mercy Ledger ✓ complete
└─ World Forge  ✓ complete

OUTLINING   ████████████████████ 100%
├─ Chapters outlined: 25/25
├─ Foreshadowing threads: 12 tracked
└─ All beats mapped: ✓

WRITING     ████████░░░░░░░░░░░░  42%
├─ Chapters written: 10/25
├─ Current chapter: 11 (scene 2/4)
├─ Word count: 63,450 / 150,000
└─ Next review: After chapter 12

EDITING     ░░░░░░░░░░░░░░░░░░░░   0%
└─ Not yet started

Last backup: 2 hours ago
```

## When to Use

- At the start of a session to see where you left off
- Before `/crucible-suite:crucible-continue` to understand project state
- When you're unsure what to do next
- To verify all planning is complete before outlining
