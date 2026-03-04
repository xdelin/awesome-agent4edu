#!/usr/bin/env python3
"""
Script: backup_on_change.py
Purpose: Incremental backup triggered by PostToolUse hook on Write|Edit
         Also checks for chapter completion and reminds about bi-chapter reviews
Crucible Suite Plugin
"""

import sys
import json
import shutil
import re
from pathlib import Path
from datetime import datetime

# Ensure Python 3.8+
if sys.version_info < (3, 8):
    print("Error: Python 3.8+ required", file=sys.stderr)
    sys.exit(1)

# Import shared utilities from cross_platform
from cross_platform import (
    find_crucible_project_with_type,
    get_backup_directory,
    encode_path_b64
)


def get_timestamp():
    return datetime.now().strftime("%Y%m%d-%H%M%S")


def should_backup_file(file_path: str) -> bool:
    """Determine if a file change should trigger a backup."""
    if not file_path:
        return False

    path = Path(file_path)
    path_str = str(path).lower()

    # Never backup backup files
    if "backup" in path_str:
        return False

    # Backup if it's in .crucible directory
    if ".crucible" in path.parts:
        return True

    # Backup CLAUDE.md
    if path.name == "CLAUDE.md":
        return True

    # Backup markdown files in common Crucible directories
    crucible_dirs = ["draft", "chapters", "planning", "manuscript"]
    if path.suffix == ".md":
        for dir_name in crucible_dirs:
            if dir_name in path.parts:
                return True

    # Backup key JSON files at project root
    key_json_files = ["story-bible.json", "style-profile.json"]
    if path.name in key_json_files:
        return True

    # Backup JSON state files
    if path.suffix == ".json" and "state" in path.parts:
        return True

    # Backup any file in a crucible-project directory structure
    # Check if any parent directory contains .crucible
    try:
        for parent in path.parents:
            if (parent / ".crucible").exists():
                # This file is in a Crucible project
                # Backup markdown and JSON files
                if path.suffix in [".md", ".json"]:
                    return True
                break
    except (OSError, PermissionError):
        pass

    return False


def is_chapter_file(file_path: str) -> tuple:
    """
    Check if a file is a chapter file and extract chapter number.

    Returns:
        tuple: (is_chapter: bool, chapter_number: int or None)
    """
    if not file_path:
        return False, None

    path = Path(file_path)
    name = path.name.lower()

    # Common chapter file patterns
    patterns = [
        r'chapter[_-]?(\d+)',
        r'ch[_-]?(\d+)',
        r'(\d+)[_-]?chapter',
    ]

    for pattern in patterns:
        match = re.search(pattern, name, re.IGNORECASE)
        if match:
            return True, int(match.group(1))

    # Check if in a chapters/draft directory
    if any(part in ['chapters', 'draft', 'manuscript'] for part in path.parts):
        # Try to find a number in the filename
        numbers = re.findall(r'\d+', name)
        if numbers:
            return True, int(numbers[0])

    return False, None


def check_review_status(project_root: Path) -> dict:
    """
    Check if a bi-chapter review is due.

    NOTE: This logic must stay in sync with sync_draft_state() in
    skills/crucible-writer/scripts/update_story_bible.py, which is the
    CANONICAL SOURCE for bi-chapter review logic. That function is used
    when chapters are completed via save_draft.py. This function is used
    by the backup hook to remind users during active writing.

    The review trigger rule: A review is needed when 2 or more chapters
    have been completed since the last review.

    Returns:
        dict with review status information
    """
    if project_root is None:
        return {"review_due": False}

    state_dir = project_root / ".crucible" / "state"
    draft_state_file = state_dir / "draft-state.json"

    if not draft_state_file.exists():
        return {"review_due": False}

    try:
        with open(draft_state_file, "r", encoding="utf-8") as f:
            state = json.load(f)
    except (json.JSONDecodeError, OSError):
        return {"review_due": False}

    chapters_complete = state.get("chapters_complete", 0)
    last_review_chapter = state.get("last_review_at_chapter", 0)

    # Check if review is due: trigger when 2+ chapters since last review
    # This matches the logic in update_story_bible.py:sync_draft_state()
    chapters_since_review = chapters_complete - last_review_chapter

    if chapters_since_review >= 2:
        return {
            "review_due": True,
            "chapters_complete": chapters_complete,
            "review_start": last_review_chapter + 1,
            "review_end": chapters_complete
        }

    return {
        "review_due": False,
        "chapters_complete": chapters_complete
    }


def incremental_backup(file_path: str, project_root: Path = None, structure_type: str = None) -> dict:
    """Create an incremental backup of a changed file."""

    if not should_backup_file(file_path):
        return {
            "success": True,
            "action": "skipped",
            "reason": "File not in backup scope"
        }

    source = Path(file_path)
    if not source.exists():
        return {
            "success": True,
            "action": "skipped",
            "reason": "File does not exist (may be new)"
        }

    # Find project root and structure type using shared detection logic
    if project_root is None or structure_type is None:
        detected_root, detected_type = find_crucible_project_with_type(source.parent)
        if project_root is None:
            project_root = detected_root
        if structure_type is None:
            structure_type = detected_type

    if project_root is None:
        return {
            "success": False,
            "error": "Could not find Crucible project root"
        }

    # Default to dotcrucible if type detection failed
    if structure_type is None:
        structure_type = "dotcrucible"

    # Set up incremental backup directory based on structure type
    backup_base = get_backup_directory(project_root, structure_type)
    backup_dir = backup_base / "incremental"
    backup_dir.mkdir(parents=True, exist_ok=True)

    # Create backup with timestamp
    timestamp = get_timestamp()
    # Use try/except for Python 3.8 compatibility (is_relative_to added in 3.9)
    try:
        relative_path = source.relative_to(project_root)
    except ValueError:
        relative_path = Path(source.name)
    # Use base64url encoding for reversible, unambiguous path encoding
    # This replaces the old underscore-based flattening which was ambiguous
    safe_name = encode_path_b64(str(relative_path))
    backup_name = f"{timestamp}-{safe_name}"
    backup_path = backup_dir / backup_name

    try:
        shutil.copy2(source, backup_path)

        # Clean up old incremental backups (keep last 100)
        # Wrap in try/except to handle permission errors or file-in-use on Windows
        try:
            all_backups = sorted(backup_dir.glob("*"), key=lambda x: x.stat().st_mtime)
            if len(all_backups) > 100:
                for old_backup in all_backups[:-100]:
                    try:
                        old_backup.unlink()
                    except OSError:
                        # Skip files that can't be deleted (in use, permission denied)
                        pass
        except OSError:
            # If we can't list/stat backups, skip cleanup but don't fail the backup
            pass

        return {
            "success": True,
            "action": "backed_up",
            "source": str(source),
            "backup": str(backup_path),
            "timestamp": timestamp
        }

    except Exception as e:
        return {
            "success": False,
            "error": f"Backup failed: {str(e)}"
        }


def main():
    """
    Main entry point for PostToolUse hook.

    Output format follows official Claude Code hooks specification (hooks.md):
    - Exit 0 with hookSpecificOutput.additionalContext when context needed
    - Exit 0 with empty JSON {} when no context needed
    - Exit 2 with stderr message to block (not used here - backups are non-blocking)
    """
    # Read hook input from stdin
    try:
        input_data = json.load(sys.stdin)
    except (json.JSONDecodeError, OSError):
        # Handle both JSON parse errors and stdin read errors
        input_data = {}

    # Extract file path from hook data
    # PostToolUse hook provides tool_input which may contain file_path
    file_path = None

    tool_input = input_data.get("tool_input", {})
    if isinstance(tool_input, dict):
        file_path = tool_input.get("file_path")

    # Also check direct input
    if not file_path:
        file_path = input_data.get("file_path")

    if not file_path:
        # No file path - nothing to do, exit silently with empty output
        print("{}")
        sys.exit(0)

    # Perform backup (errors are logged internally, non-blocking)
    result = incremental_backup(file_path)

    # Check if this was a chapter file and if review is now due
    try:
        is_chapter, chapter_num = is_chapter_file(file_path)
        file_parent = Path(file_path).parent if file_path else None
        project_root, _ = find_crucible_project_with_type(file_parent)
        review_status = check_review_status(project_root)
    except OSError:
        # Handle filesystem errors gracefully
        is_chapter = False
        chapter_num = None
        review_status = {"review_due": False}

    # Build spec-compliant output: ONLY hookSpecificOutput when context is needed
    # Per hooks.md, PostToolUse additionalContext adds info for Claude to consider
    output = {}

    # If a chapter file was modified and review is due, add context for Claude
    if is_chapter and review_status.get("review_due"):
        review_start = review_status["review_start"]
        review_end = review_status["review_end"]

        output["hookSpecificOutput"] = {
            "hookEventName": "PostToolUse",
            "additionalContext": (
                f"[CRUCIBLE REMINDER] Chapter {chapter_num} was just modified. "
                f"You have now completed {review_status['chapters_complete']} chapters. "
                f"A bi-chapter review is due for chapters {review_start}-{review_end}. "
                f"Please run /crucible-suite:crucible-review {review_start}-{review_end} "
                f"before continuing to the next chapter."
            )
        }

    print(json.dumps(output, indent=2) if output else "{}")
    sys.exit(0)


if __name__ == "__main__":
    main()
