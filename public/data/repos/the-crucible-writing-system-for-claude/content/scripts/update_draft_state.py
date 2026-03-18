#!/usr/bin/env python3
"""
Script: update_draft_state.py
Purpose: Update draft state for chapter tracking and bi-chapter review triggers
Crucible Suite Plugin
"""

import sys
import json
import os
from pathlib import Path
from datetime import datetime
import argparse

# Ensure Python 3.8+
if sys.version_info < (3, 8):
    print("Error: Python 3.8+ required", file=sys.stderr)
    sys.exit(1)

# Import shared project detection from cross_platform
from cross_platform import find_crucible_project


def _load_story_bible(project_root: Path) -> dict:
    """Load story bible if it exists, return None otherwise."""
    if project_root is None:
        return None
    bible_path = project_root / "story-bible.json"
    if not bible_path.exists():
        return None
    try:
        with open(bible_path, "r", encoding="utf-8") as f:
            return json.load(f)
    except (json.JSONDecodeError, OSError):
        return None


def _sync_with_story_bible(project_root: Path, state: dict) -> None:
    """Sync draft state with story bible to prevent desynchronization."""
    bible = _load_story_bible(project_root)
    if bible is None:
        return

    if "progress" not in bible:
        return

    # Sync progress fields
    bible["progress"]["current_chapter"] = state.get("current_chapter", 1)
    bible["progress"]["current_scene"] = state.get("current_scene", 1)
    bible["progress"]["chapters_complete"] = state.get("chapters_complete", 0)

    if "meta" in bible:
        bible["meta"]["updated"] = datetime.now().isoformat()

    # Save updated bible
    bible_path = project_root / "story-bible.json"
    try:
        with open(bible_path, "w", encoding="utf-8") as f:
            json.dump(bible, f, indent=2)
    except OSError:
        pass  # Best effort sync


def load_draft_state(project_root: Path) -> dict:
    """Load the current draft state or create default."""
    state_file = project_root / ".crucible" / "state" / "draft-state.json"

    default_state = {
        "chapters_complete": 0,
        "last_review_at_chapter": 0,
        "review_pending": False,
        "current_chapter": 1,
        "current_scene": 1,
        "total_chapters": 25,
        "target_words": 150000,
        "last_updated": None
    }

    if state_file.exists():
        try:
            with open(state_file, "r", encoding="utf-8") as f:
                state = json.load(f)
                # Merge with defaults to ensure all fields exist
                for key, value in default_state.items():
                    if key not in state:
                        state[key] = value
                return state
        except (json.JSONDecodeError, OSError):
            return default_state

    return default_state


def save_draft_state(project_root: Path, state: dict) -> bool:
    """Save the draft state to disk."""
    state_dir = project_root / ".crucible" / "state"
    state_dir.mkdir(parents=True, exist_ok=True)

    state_file = state_dir / "draft-state.json"
    state["last_updated"] = datetime.now().isoformat()

    try:
        with open(state_file, "w", encoding="utf-8") as f:
            json.dump(state, f, indent=2)
        return True
    except OSError as e:
        print(f"Error saving state: {e}", file=sys.stderr)
        return False


def update_chapter_complete(project_root: Path, chapter_num: int) -> dict:
    """Mark a chapter as complete and check for review trigger."""
    state = load_draft_state(project_root)

    # Update chapter count
    state["chapters_complete"] = max(state["chapters_complete"], chapter_num)
    state["current_chapter"] = chapter_num + 1
    state["current_scene"] = 1

    # Check if bi-chapter review is needed
    # REVIEW LOGIC: Trigger when 2+ chapters completed since last review.
    # Canonical source: sync_draft_state() in update_story_bible.py
    # Also mirrored in: backup_on_change.py, load_project_context.py
    chapters_since_review = state["chapters_complete"] - state["last_review_at_chapter"]

    if chapters_since_review >= 2:
        state["review_pending"] = True
        review_needed = True
        review_chapters = (state["last_review_at_chapter"] + 1, state["chapters_complete"])
    else:
        review_needed = False
        review_chapters = None

    save_draft_state(project_root, state)

    # Sync with story bible to prevent desynchronization
    _sync_with_story_bible(project_root, state)

    return {
        "success": True,
        "chapters_complete": state["chapters_complete"],
        "review_pending": state["review_pending"],
        "review_needed": review_needed,
        "review_chapters": review_chapters,
        "next_chapter": state["current_chapter"]
    }


def update_review_complete(project_root: Path, reviewed_up_to: int = None) -> dict:
    """Mark a bi-chapter review as complete."""
    state = load_draft_state(project_root)

    if reviewed_up_to is None:
        reviewed_up_to = state["chapters_complete"]

    state["last_review_at_chapter"] = reviewed_up_to
    state["review_pending"] = False
    state["last_review_timestamp"] = datetime.now().isoformat()

    save_draft_state(project_root, state)

    return {
        "success": True,
        "last_review_at_chapter": state["last_review_at_chapter"],
        "review_pending": False
    }


def update_current_position(project_root: Path, chapter: int, scene: int) -> dict:
    """Update the current writing position."""
    state = load_draft_state(project_root)

    state["current_chapter"] = chapter
    state["current_scene"] = scene

    save_draft_state(project_root, state)

    return {
        "success": True,
        "current_chapter": chapter,
        "current_scene": scene
    }


def get_status(project_root: Path) -> dict:
    """Get the current draft state status."""
    state = load_draft_state(project_root)

    chapters_since_review = state["chapters_complete"] - state["last_review_at_chapter"]
    chapters_until_review = max(0, 2 - chapters_since_review)

    return {
        "success": True,
        "chapters_complete": state["chapters_complete"],
        "current_chapter": state["current_chapter"],
        "current_scene": state["current_scene"],
        "total_chapters": state["total_chapters"],
        "last_review_at_chapter": state["last_review_at_chapter"],
        "review_pending": state["review_pending"],
        "chapters_until_review": chapters_until_review,
        "last_updated": state["last_updated"]
    }


def main():
    parser = argparse.ArgumentParser(description="Update Crucible draft state")
    parser.add_argument("project_path", nargs="?", default=".",
                        help="Path to Crucible project")
    parser.add_argument("--chapter-complete", type=int,
                        help="Mark chapter N as complete")
    parser.add_argument("--review-complete", type=int, nargs="?", const=0,
                        help="Mark review as complete (optionally specify chapter, 0=use current)")
    parser.add_argument("--set-position", nargs=2, type=int, metavar=("CHAPTER", "SCENE"),
                        help="Set current writing position")
    parser.add_argument("--status", action="store_true",
                        help="Get current draft status")
    parser.add_argument("--json", action="store_true",
                        help="Output as JSON")

    args = parser.parse_args()

    # Find project
    project_root = find_crucible_project(Path(args.project_path))

    if project_root is None:
        result = {"success": False, "error": "No Crucible project found"}
        if args.json:
            print(json.dumps(result, indent=2))
        else:
            print(f"Error: {result['error']}", file=sys.stderr)
        sys.exit(1)

    # Execute requested action
    result = None

    if args.chapter_complete is not None:
        result = update_chapter_complete(project_root, args.chapter_complete)
        if not args.json and result["review_needed"]:
            print(f"Chapter {args.chapter_complete} marked complete.")
            print(f"[WARN] BI-CHAPTER REVIEW NEEDED for chapters {result['review_chapters'][0]}-{result['review_chapters'][1]}")

    elif args.review_complete is not None:
        # 0 means use current chapter count, positive means specific chapter
        chapter = args.review_complete if args.review_complete > 0 else None
        result = update_review_complete(project_root, chapter)
        if not args.json:
            print(f"Review marked complete through chapter {result['last_review_at_chapter']}.")

    elif args.set_position is not None:
        chapter, scene = args.set_position
        result = update_current_position(project_root, chapter, scene)
        if not args.json:
            print(f"Position set to Chapter {chapter}, Scene {scene}.")

    elif args.status:
        result = get_status(project_root)
        if not args.json:
            print(f"Chapters complete: {result['chapters_complete']}/{result['total_chapters']}")
            print(f"Current position: Chapter {result['current_chapter']}, Scene {result['current_scene']}")
            print(f"Last review at chapter: {result['last_review_at_chapter']}")
            if result['review_pending']:
                print("[WARN] REVIEW PENDING")
            else:
                print(f"Chapters until next review: {result['chapters_until_review']}")

    else:
        # Default: show status
        result = get_status(project_root)
        if not args.json:
            print(f"Chapters complete: {result['chapters_complete']}/{result['total_chapters']}")
            print(f"Review pending: {result['review_pending']}")

    if args.json and result:
        print(json.dumps(result, indent=2))

    sys.exit(0 if result and result.get("success") else 1)


if __name__ == "__main__":
    main()
