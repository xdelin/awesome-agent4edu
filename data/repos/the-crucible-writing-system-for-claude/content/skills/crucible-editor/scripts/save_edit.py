#!/usr/bin/env python3
"""
Script: save_edit.py
Purpose: Save editing progress and track changes
Crucible Suite Plugin
"""

import sys
import json
import os
from pathlib import Path
from datetime import datetime
import shutil

# Ensure Python 3.8+
if sys.version_info < (3, 8):
    print("Error: Python 3.8+ required", file=sys.stderr)
    sys.exit(1)


def find_state_file(project_dir: Path) -> Path:
    """Find the edit state file, checking new location first then legacy."""
    # New location: .crucible/state/edit-state.json
    new_path = project_dir / ".crucible" / "state" / "edit-state.json"
    if new_path.exists():
        return new_path

    # Legacy location: edit-state.json at project root
    legacy_path = project_dir / "edit-state.json"
    if legacy_path.exists():
        return legacy_path

    return None


def save_editing_progress(project_path: str, chapter_num: int = None,
                         edit_type: str = None, changes: list = None):
    """Save editing progress and update tracking."""
    project_dir = Path(project_path)

    # Load current state
    state_file = find_state_file(project_dir)
    if not state_file:
        return {"success": False, "error": "Project not initialized"}

    with open(state_file, "r", encoding="utf-8") as f:
        state = json.load(f)

    # Update state
    state["last_modified"] = datetime.now().isoformat()

    if chapter_num is not None:
        chapter_key = f"chapter_{chapter_num}"
        if chapter_key not in state["chapters"]:
            state["chapters"][chapter_key] = {
                "developmental": False,
                "line_edit": False,
                "copy_edit": False,
                "polish": False,
                "changes": []
            }

        if edit_type:
            state["chapters"][chapter_key][edit_type] = True

        if changes:
            state["chapters"][chapter_key]["changes"].extend(changes)

    # Record in history
    state["edit_history"].append({
        "timestamp": datetime.now().isoformat(),
        "chapter": chapter_num,
        "edit_type": edit_type,
        "change_count": len(changes) if changes else 0
    })

    # Save state
    with open(state_file, "w", encoding="utf-8") as f:
        json.dump(state, f, indent=2)

    # Update tracking
    tracking_file = project_dir / "edit-tracking.json"
    if tracking_file.exists():
        with open(tracking_file, "r", encoding="utf-8") as f:
            tracking = json.load(f)

        if chapter_num is not None:
            if chapter_num not in tracking["chapters_edited"]:
                tracking["chapters_edited"].append(chapter_num)

        with open(tracking_file, "w", encoding="utf-8") as f:
            json.dump(tracking, f, indent=2)

    return {
        "success": True,
        "message": f"Progress saved for chapter {chapter_num}" if chapter_num else "Progress saved",
        "timestamp": datetime.now().isoformat()
    }


def main():
    # Read input
    if len(sys.argv) >= 2:
        project_path = sys.argv[1]
        chapter_num = int(sys.argv[2]) if len(sys.argv) > 2 else None
        edit_type = sys.argv[3] if len(sys.argv) > 3 else None
    else:
        try:
            input_data = json.load(sys.stdin)
            project_path = input_data.get("project_path", "./edit-project")
            chapter_num = input_data.get("chapter_num")
            edit_type = input_data.get("edit_type")
        except json.JSONDecodeError:
            print(json.dumps({"success": False, "error": "Invalid input"}))
            sys.exit(1)

    result = save_editing_progress(project_path, chapter_num, edit_type)
    print(json.dumps(result, indent=2))
    sys.exit(0 if result["success"] else 1)


if __name__ == "__main__":
    main()
