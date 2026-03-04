---
allowed-tools: Read, Write, Glob, Bash
argument-hint: [timestamp] | "latest"
description: Restore your project from a backup. Lists available backups and allows selective restoration.
model: claude-haiku-4-5-20251001
---

# /crucible-restore

Restore your Crucible project from a previous backup.

## Execution Instructions

**IMPORTANT:** When presenting backup options to the user, ALWAYS use the AskUserQuestion tool (NOT plain text A/B/C options).

### Step 1: List Available Backups

Run the restore script to list backups:

```bash
bash "${CLAUDE_PLUGIN_ROOT}/scripts/run_python.sh" "${CLAUDE_PLUGIN_ROOT}/scripts/restore_project.py" --list
```

This returns JSON with all available backups including timestamps, sizes, and file counts.

### Step 2: Present Options to User

Use AskUserQuestion to let user select which backup to restore and what scope.

### Step 3: Perform Restoration

For full restoration:
```bash
bash "${CLAUDE_PLUGIN_ROOT}/scripts/run_python.sh" "${CLAUDE_PLUGIN_ROOT}/scripts/restore_project.py" --restore "TIMESTAMP" --scope full
```

For selective restoration:
```bash
bash "${CLAUDE_PLUGIN_ROOT}/scripts/run_python.sh" "${CLAUDE_PLUGIN_ROOT}/scripts/restore_project.py" --restore "TIMESTAMP" --scope chapters
bash "${CLAUDE_PLUGIN_ROOT}/scripts/run_python.sh" "${CLAUDE_PLUGIN_ROOT}/scripts/restore_project.py" --restore "TIMESTAMP" --scope planning
bash "${CLAUDE_PLUGIN_ROOT}/scripts/run_python.sh" "${CLAUDE_PLUGIN_ROOT}/scripts/restore_project.py" --restore "TIMESTAMP" --scope state
bash "${CLAUDE_PLUGIN_ROOT}/scripts/run_python.sh" "${CLAUDE_PLUGIN_ROOT}/scripts/restore_project.py" --restore "TIMESTAMP" --scope story-bible
```

For incremental backup restoration:
```bash
bash "${CLAUDE_PLUGIN_ROOT}/scripts/run_python.sh" "${CLAUDE_PLUGIN_ROOT}/scripts/restore_project.py" --incremental --timestamp "YYYYMMDD"
```

Add `--dry-run` to preview what would be restored without making changes.

## Usage

- `/crucible-suite:crucible-restore` - List available backups
- `/crucible-suite:crucible-restore latest` - Restore from most recent backup
- `/crucible-suite:crucible-restore 2024-01-15-1432` - Restore from specific timestamp

## What This Does

1. Lists all available backups with timestamps
2. Shows what changed since each backup
3. Allows selective or full restoration
4. Creates a pre-restore backup (safety net)

## Backup Types

### Automatic Backups
Created automatically when you:
- Complete a chapter
- Finish a planning document
- Complete a bi-chapter review
- Use Write or Edit tools on project files

### Manual Backups
Can be triggered anytime with the backup script.

## Restoration Options

When restoring, you can choose:

### Full Restore
Restores entire project to backup state:
- All chapters
- Planning documents
- Story bible
- State files

### Selective Restore
Restore specific components:
- Single chapter
- Planning documents only
- Story bible only
- State files only

## Example

```
/crucible-suite:crucible-restore

Available Backups for "The Memory Forge"
═══════════════════════════════════════

1. 2024-01-15-1432 (2 hours ago)
   └─ Chapter 11 partial, 63,450 words

2. 2024-01-15-1030 (6 hours ago)
   └─ Chapter 10 complete, 60,200 words

3. 2024-01-14-2145 (yesterday)
   └─ Chapter 10 partial, 58,100 words

4. 2024-01-14-1500 (yesterday)
   └─ Bi-chapter review complete (Ch 8-9)

```

Then use AskUserQuestion to let user select:
```json
{
  "questions": [
    {
      "header": "Restore",
      "question": "Which backup would you like to restore?",
      "options": [
        {"label": "Latest (2 hrs ago)", "description": "Chapter 11 partial, 63,450 words"},
        {"label": "6 hours ago", "description": "Chapter 10 complete, 60,200 words"},
        {"label": "Yesterday evening", "description": "Chapter 10 partial, 58,100 words"},
        {"label": "Enter timestamp", "description": "Specify a different backup timestamp"}
      ],
      "multiSelect": false
    }
  ]
}
```

## Safety Features

- **Pre-restore backup**: Current state is always saved before restoration
- **Confirmation required**: Shows exactly what will change
- **Selective restore**: Can restore just what you need
- **Undo restore**: Can restore the pre-restore backup if needed

## When to Use

- After a session went wrong
- When you prefer an earlier version
- After accidental deletions
- To recover from corruption
- To compare with previous versions
