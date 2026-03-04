#!/usr/bin/env python3
"""
Script: init_edit.py
Purpose: Initialize an editing project for Crucible manuscripts
Crucible Suite Plugin
"""

import sys
import json
import os
from pathlib import Path
from datetime import datetime

# Ensure Python 3.8+
if sys.version_info < (3, 8):
    print("Error: Python 3.8+ required", file=sys.stderr)
    sys.exit(1)


def init_editing_project(project_path: str, title: str, source_draft: str = None):
    """Initialize a new editing project structure."""
    project_dir = Path(project_path)

    # Create directory structure
    directories = [
        project_dir / "original",      # Original drafts (backup)
        project_dir / "working",       # Current editing versions
        project_dir / "edited",        # Completed edited versions
        project_dir / "reports",       # Edit reports and change logs
        project_dir / ".crucible" / "state",  # State files
    ]

    for directory in directories:
        directory.mkdir(parents=True, exist_ok=True)

    # Initialize project state
    state = {
        "title": title,
        "created": datetime.now().isoformat(),
        "last_modified": datetime.now().isoformat(),
        "source_draft": source_draft,
        "current_phase": "assessment",
        "chapters": {},
        "global_issues": [],
        "edit_history": []
    }

    # Save state to .crucible/state/edit-state.json
    state_file = project_dir / ".crucible" / "state" / "edit-state.json"
    with open(state_file, "w", encoding="utf-8") as f:
        json.dump(state, f, indent=2)

    # Create edit tracking file
    tracking = {
        "developmental_complete": False,
        "line_edit_complete": False,
        "copy_edit_complete": False,
        "polish_complete": False,
        "chapters_edited": [],
        "word_count_original": 0,
        "word_count_edited": 0
    }

    tracking_file = project_dir / "edit-tracking.json"
    with open(tracking_file, "w", encoding="utf-8") as f:
        json.dump(tracking, f, indent=2)

    return {
        "success": True,
        "project_path": str(project_dir),
        "state_file": str(state_file),
        "message": f"Editing project initialized: {title}"
    }


def main():
    # Read input from stdin or command line
    if len(sys.argv) >= 3:
        project_path = sys.argv[1]
        title = sys.argv[2]
        source_draft = sys.argv[3] if len(sys.argv) > 3 else None
    else:
        try:
            input_data = json.load(sys.stdin)
            project_path = input_data.get("project_path", "./edit-project")
            title = input_data.get("title", "Untitled")
            source_draft = input_data.get("source_draft")
        except json.JSONDecodeError:
            project_path = "./edit-project"
            title = "Untitled"
            source_draft = None

    result = init_editing_project(project_path, title, source_draft)
    print(json.dumps(result, indent=2))
    sys.exit(0 if result["success"] else 1)


if __name__ == "__main__":
    main()
