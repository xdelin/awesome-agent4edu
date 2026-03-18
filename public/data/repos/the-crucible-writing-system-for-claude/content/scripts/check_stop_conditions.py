#!/usr/bin/env python3
"""
Script: check_stop_conditions.py
Purpose: Stop hook to enforce bi-chapter reviews before session end
Crucible Suite Plugin

This script runs when Claude tries to stop and blocks if a bi-chapter
review is due but hasn't been performed.

Supports all project structure types:
- dotcrucible: .crucible/state/draft-state.json (standard structure)
- rootlevel: project-state.json at root (detected by state.json presence)
- legacy: project-state.json at root (detected by planning/ dir or story-bible.json)
"""

import sys
import json
from pathlib import Path
from typing import Optional

# Ensure Python 3.8+
if sys.version_info < (3, 8):
    print("Error: Python 3.8+ required", file=sys.stderr)
    sys.exit(1)

# Import shared project detection from cross_platform
from cross_platform import find_crucible_project_with_type


def get_draft_state_path(project_root: Path, structure_type: str) -> Optional[Path]:
    """
    Get the draft state file path based on project structure type.

    Args:
        project_root: Path to the project root
        structure_type: One of "dotcrucible", "rootlevel", or "legacy"

    Returns:
        Path to the draft state file if found, None otherwise
    """
    # Paths to check in order of preference
    paths_to_check = []

    if structure_type == "dotcrucible":
        # Standard structure: .crucible/state/draft-state.json
        paths_to_check = [
            project_root / ".crucible" / "state" / "draft-state.json",
        ]
    elif structure_type == "rootlevel":
        # Planner-created structure: check both locations
        paths_to_check = [
            project_root / ".crucible" / "state" / "draft-state.json",
            project_root / "project-state.json",
        ]
    else:
        # Legacy structure: check multiple possible locations
        paths_to_check = [
            project_root / "project-state.json",
            project_root / ".crucible" / "state" / "draft-state.json",
        ]

    # Return first existing path
    for path in paths_to_check:
        if path.exists():
            return path

    return None


def check_review_needed(project_root: Path, structure_type: str = None) -> dict:
    """
    Check if a bi-chapter review is needed.

    Supports all project structure types for finding draft state.

    Args:
        project_root: Path to the project root
        structure_type: Project structure type (dotcrucible, rootlevel, legacy)

    Returns:
        dict with keys:
        - review_needed: bool
        - chapters_to_review: tuple (start, end) or None
        - reason: str
    """
    if project_root is None:
        return {
            "review_needed": False,
            "chapters_to_review": None,
            "reason": "No Crucible project found"
        }

    # Default to dotcrucible if not specified
    if structure_type is None:
        structure_type = "dotcrucible"

    # Get draft state path based on structure type
    draft_state_file = get_draft_state_path(project_root, structure_type)

    if draft_state_file is None:
        return {
            "review_needed": False,
            "chapters_to_review": None,
            "reason": "No draft state found (not in writing phase)"
        }

    try:
        with open(draft_state_file, "r", encoding="utf-8") as f:
            state = json.load(f)
    except (json.JSONDecodeError, OSError) as e:
        return {
            "review_needed": False,
            "chapters_to_review": None,
            "reason": f"Could not read draft state: {e}"
        }

    chapters_complete = state.get("chapters_complete", 0)
    last_review_at_chapter = state.get("last_review_at_chapter", 0)

    # Bi-chapter review is due when 2+ chapters completed since last review.
    # REVIEW LOGIC: Must match sync_draft_state() in update_story_bible.py
    # Also mirrored in: backup_on_change.py, load_project_context.py, update_draft_state.py
    chapters_since_review = chapters_complete - last_review_at_chapter

    if chapters_since_review >= 2:
        review_start = last_review_at_chapter + 1
        review_end = chapters_complete
        return {
            "review_needed": True,
            "chapters_to_review": (review_start, review_end),
            "chapters_complete": chapters_complete,
            "last_review_at_chapter": last_review_at_chapter,
            "reason": f"Bi-chapter review required for chapters {review_start}-{review_end}"
        }

    return {
        "review_needed": False,
        "chapters_to_review": None,
        "chapters_complete": chapters_complete,
        "last_review_at_chapter": last_review_at_chapter,
        "reason": "No review currently due"
    }


def main():
    # Read hook input from stdin
    try:
        input_data = json.load(sys.stdin)
    except json.JSONDecodeError:
        input_data = {}

    # Check if stop hook is already active (prevent infinite loops)
    stop_hook_active = input_data.get("stop_hook_active", False)
    if stop_hook_active:
        # Don't block if we're already in a stop hook continuation
        # This prevents infinite loops
        sys.exit(0)

    # Find project and check conditions using structure-aware detection
    project_root, structure_type = find_crucible_project_with_type()
    result = check_review_needed(project_root, structure_type)

    if result["review_needed"]:
        start, end = result["chapters_to_review"]
        # Exit code 2 blocks the stop and sends stderr to Claude
        error_message = (
            f"STOP BLOCKED: Bi-chapter review required.\n\n"
            f"You have completed {result['chapters_complete']} chapters but the last review "
            f"only covered up to chapter {result['last_review_at_chapter']}.\n\n"
            f"ACTION REQUIRED: Run the bi-chapter review for chapters {start}-{end} using:\n"
            f"  /crucible-suite:crucible-review {start}-{end}\n\n"
            f"After completing the review, you may stop the session."
        )
        print(error_message, file=sys.stderr)
        sys.exit(2)
    else:
        # Exit code 0 allows the stop
        sys.exit(0)


if __name__ == "__main__":
    main()
