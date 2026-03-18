#!/usr/bin/env python3
"""
Script: migrate_project.py
Purpose: Auto-migrate Crucible projects between phases and legacy formats
Crucible Suite Plugin

This script handles:
1. Migrating legacy state.json -> .crucible/state/planning-state.json
2. Migrating legacy outline state.json -> .crucible/state/outline-state.json (when outline/ exists)
3. Migrating legacy project-state.json -> .crucible/state/draft-state.json
4. Ensuring .crucible/state/ directory exists
5. Creating story-bible.json when transitioning to writing phase
"""

import json
import os
import sys
import shutil
from pathlib import Path
from datetime import datetime


def ensure_crucible_structure(project_path: Path) -> None:
    """Ensure .crucible/state/ directory exists."""
    state_dir = project_path / ".crucible" / "state"
    state_dir.mkdir(parents=True, exist_ok=True)

    backups_dir = project_path / ".crucible" / "backups"
    backups_dir.mkdir(parents=True, exist_ok=True)


def migrate_planning_state(project_path: Path) -> bool:
    """
    Migrate legacy state.json to .crucible/state/planning-state.json.
    Returns True if migration occurred.
    """
    legacy_path = project_path / "state.json"
    new_path = project_path / ".crucible" / "state" / "planning-state.json"

    if legacy_path.exists() and not new_path.exists():
        ensure_crucible_structure(project_path)

        # Read legacy state
        with open(legacy_path, "r", encoding="utf-8") as f:
            state = json.load(f)

        # Write to new location
        with open(new_path, "w", encoding="utf-8") as f:
            json.dump(state, f, indent=2)

        # Rename legacy file to indicate migration
        backup_path = project_path / "state.json.migrated"
        shutil.move(str(legacy_path), str(backup_path))

        print(f"[OK] Migrated planning state: state.json -> .crucible/state/planning-state.json")
        return True

    return False


def migrate_draft_state(project_path: Path) -> bool:
    """
    Migrate legacy project-state.json to .crucible/state/draft-state.json.
    Returns True if migration occurred.
    """
    legacy_path = project_path / "project-state.json"
    new_path = project_path / ".crucible" / "state" / "draft-state.json"

    if legacy_path.exists() and not new_path.exists():
        ensure_crucible_structure(project_path)

        # Read legacy state
        with open(legacy_path, "r", encoding="utf-8") as f:
            state = json.load(f)

        # Write to new location
        with open(new_path, "w", encoding="utf-8") as f:
            json.dump(state, f, indent=2)

        # Rename legacy file to indicate migration
        backup_path = project_path / "project-state.json.migrated"
        shutil.move(str(legacy_path), str(backup_path))

        print(f"[OK] Migrated draft state: project-state.json -> .crucible/state/draft-state.json")
        return True

    return False


def migrate_outline_state(project_path: Path) -> bool:
    """
    Migrate legacy outline state.json to .crucible/state/outline-state.json.
    Note: Outline projects use state.json at project root (same as planning).
    This only migrates if there's an outline/ directory present (to distinguish from planning).
    Returns True if migration occurred.
    """
    legacy_path = project_path / "state.json"
    new_path = project_path / ".crucible" / "state" / "outline-state.json"
    outline_dir = project_path / "outline"

    # Only migrate if outline directory exists (indicates this is an outline project)
    if legacy_path.exists() and outline_dir.exists() and not new_path.exists():
        ensure_crucible_structure(project_path)

        # Read legacy state
        with open(legacy_path, "r", encoding="utf-8") as f:
            state = json.load(f)

        # Check if this looks like an outline state (has chapters array, not answers dict)
        if "chapters" in state and isinstance(state.get("chapters"), list):
            # Write to new location
            with open(new_path, "w", encoding="utf-8") as f:
                json.dump(state, f, indent=2)

            # Rename legacy file to indicate migration
            backup_path = project_path / "state.json.migrated"
            shutil.move(str(legacy_path), str(backup_path))

            print(f"[OK] Migrated outline state: state.json -> .crucible/state/outline-state.json")
            return True

    return False


def create_story_bible_from_planning(project_path: Path, title: str = None, target_words: int = 100000, target_chapters: int = 25) -> bool:
    """
    Create story-bible.json when transitioning from planning to writing.
    Returns True if created.
    """
    bible_path = project_path / "story-bible.json"

    if bible_path.exists():
        print(f"[INFO] story-bible.json already exists")
        return False

    # Try to get title from planning state
    if title is None:
        planning_state = project_path / ".crucible" / "state" / "planning-state.json"
        legacy_state = project_path / "state.json"

        state_to_check = planning_state if planning_state.exists() else (legacy_state if legacy_state.exists() else None)

        if state_to_check:
            try:
                with open(state_to_check, "r", encoding="utf-8") as f:
                    state = json.load(f)
                    title = state.get("project", {}).get("title", "Untitled")
            except (json.JSONDecodeError, OSError):
                title = "Untitled"
        else:
            title = "Untitled"

    # Calculate words per chapter
    words_per_chapter = target_words // target_chapters

    # Create story bible structure
    story_bible = {
        "meta": {
            "title": title,
            "series": None,
            "book_number": 1,
            "target_words": target_words,
            "target_chapters": target_chapters,
            "words_per_chapter": words_per_chapter,
            "created": datetime.now().isoformat(),
            "updated": datetime.now().isoformat()
        },
        "progress": {
            "current_chapter": 1,
            "current_scene": 1,
            "total_words": 0,
            "chapters_complete": 0,
            "status": "initialized"
        },
        "chapters": {},
        "character_states": {},
        "established_facts": [],
        "foreshadowing": {
            "planted": [],
            "paid_off": []
        },
        "timeline": {},
        "invented_details": [],
        "locations": {},
        "relationships": {},
        "mercy_engine": {"mercy_acts": [], "mercy_refused": [], "mercy_balance": 0}
    }

    with open(bible_path, "w", encoding="utf-8") as f:
        json.dump(story_bible, f, indent=2)

    print(f"[OK] Created story-bible.json for: {title}")
    return True


def create_style_profile(project_path: Path) -> bool:
    """
    Create style-profile.json if it doesn't exist.
    Returns True if created.
    """
    style_path = project_path / "style-profile.json"

    if style_path.exists():
        print(f"[INFO] style-profile.json already exists")
        return False

    style_profile = {
        "meta": {
            "captured": False,
            "sample_source": None,
            "updated": datetime.now().isoformat()
        },
        "sentences": {
            "average_length": None,
            "variation": None,
            "structure_mix": None,
            "rhythm": None
        },
        "vocabulary": {
            "register": None,
            "density": None,
            "signature_words": [],
            "avoid_words": []
        },
        "dialogue": {
            "attribution_style": None,
            "character_distinction": None,
            "subtext_level": None
        },
        "description": {
            "sensory_focus": None,
            "metaphor_density": None,
            "setting_integration": None
        },
        "interiority": {
            "pov_depth": None,
            "emotion_approach": None,
            "reflection_frequency": None
        },
        "pacing": {
            "transition_style": None,
            "white_space": None,
            "tension_technique": None
        },
        "signature_elements": [],
        "genre_conventions": {
            "genre": None,
            "subgenre": None,
            "specific_conventions": []
        }
    }

    with open(style_path, "w", encoding="utf-8") as f:
        json.dump(style_profile, f, indent=2)

    print(f"[OK] Created style-profile.json")
    return True


def transition_to_writing(project_path: Path, title: str = None, target_words: int = 100000, target_chapters: int = 25) -> dict:
    """
    Transition a project from planning to writing phase.
    Creates all necessary files and directories.
    """
    results = {
        "structure_created": False,
        "planning_migrated": False,
        "story_bible_created": False,
        "style_profile_created": False,
        "draft_directories_created": False
    }

    # Ensure .crucible structure
    ensure_crucible_structure(project_path)
    results["structure_created"] = True

    # Migrate planning state if needed
    results["planning_migrated"] = migrate_planning_state(project_path)

    # Create story bible
    results["story_bible_created"] = create_story_bible_from_planning(
        project_path, title, target_words, target_chapters
    )

    # Create style profile
    results["style_profile_created"] = create_style_profile(project_path)

    # Create draft directories
    draft_dir = project_path / "draft" / "chapters"
    manuscript_dir = project_path / "manuscript"

    if not draft_dir.exists():
        draft_dir.mkdir(parents=True, exist_ok=True)
        results["draft_directories_created"] = True
        print(f"[OK] Created draft/chapters/ directory")

    if not manuscript_dir.exists():
        manuscript_dir.mkdir(parents=True, exist_ok=True)
        print(f"[OK] Created manuscript/ directory")

    return results


def migrate_story_bible(project_path: Path) -> bool:
    """
    Migrate legacy story-bible.json to include all current schema fields.
    Returns True if migration occurred, False otherwise.
    """
    bible_path = project_path / "story-bible.json"

    if not bible_path.exists():
        return False

    try:
        with open(bible_path, "r", encoding="utf-8") as f:
            bible = json.load(f)
    except (json.JSONDecodeError, OSError) as e:
        print(f"[ERROR] Failed to read story-bible.json: {e}")
        return False

    # Inline schema defaults (no import dependency)
    schema_defaults = {
        "meta": {
            "title": "Untitled",
            "series": None,
            "book_number": 1,
            "target_words": 100000,
            "target_chapters": 25,
            "words_per_chapter": 4000,
            "created": datetime.now().isoformat(),
            "updated": datetime.now().isoformat()
        },
        "progress": {
            "current_chapter": 1,
            "current_scene": 1,
            "total_words": 0,
            "chapters_complete": 0,
            "status": "initialized"
        },
        "chapters": {},
        "character_states": {},
        "established_facts": [],
        "foreshadowing": {
            "planted": [],
            "paid_off": []
        },
        "timeline": {},
        "invented_details": [],
        "locations": {},
        "relationships": {},
        "mercy_engine": {
            "mercy_acts": [],
            "mercy_refused": [],
            "mercy_balance": 0
        }
    }

    migrated = False

    # Check for missing top-level keys
    for key, default_value in schema_defaults.items():
        if key not in bible:
            bible[key] = default_value
            migrated = True

    # Check for missing nested keys in structural sections
    nested_sections = ["foreshadowing", "mercy_engine", "meta", "progress"]
    for section in nested_sections:
        if section in bible and isinstance(bible[section], dict):
            for nested_key, nested_default in schema_defaults[section].items():
                if nested_key not in bible[section]:
                    bible[section][nested_key] = nested_default
                    migrated = True

    if migrated:
        try:
            # Update the updated timestamp
            if "meta" in bible:
                bible["meta"]["updated"] = datetime.now().isoformat()

            with open(bible_path, "w", encoding="utf-8") as f:
                json.dump(bible, f, indent=2)
            print(f"[OK] Migrated story-bible.json: added missing schema fields")
            return True
        except OSError as e:
            print(f"[ERROR] Failed to save migrated story-bible.json: {e}")
            return False
    else:
        print(f"[INFO] story-bible.json already up to date")
        return False


def auto_migrate(project_path: Path) -> dict:
    """
    Automatically detect and perform any needed migrations.
    """
    results = {
        "planning_migrated": False,
        "outline_migrated": False,
        "draft_migrated": False,
        "story_bible_migrated": False,
        "structure_ensured": False
    }

    # Ensure structure exists
    ensure_crucible_structure(project_path)
    results["structure_ensured"] = True

    # Check and migrate planning state
    results["planning_migrated"] = migrate_planning_state(project_path)

    # Check and migrate outline state
    results["outline_migrated"] = migrate_outline_state(project_path)

    # Check and migrate draft state
    results["draft_migrated"] = migrate_draft_state(project_path)

    # Check and migrate story bible
    results["story_bible_migrated"] = migrate_story_bible(project_path)

    if not any([results["planning_migrated"], results["outline_migrated"], results["draft_migrated"], results["story_bible_migrated"]]):
        print(f"[INFO] No migrations needed")
    else:
        print(f"[OK] Migration complete")

    return results


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Migrate Crucible projects between phases and formats"
    )
    parser.add_argument(
        "project_path",
        nargs="?",
        default=".",
        help="Path to the project directory"
    )
    parser.add_argument(
        "--to-writing",
        action="store_true",
        help="Transition project to writing phase (creates story-bible, etc.)"
    )
    parser.add_argument(
        "--title",
        help="Book title (for --to-writing)"
    )
    parser.add_argument(
        "--target-words",
        type=int,
        default=100000,
        help="Target word count (default: 100000)"
    )
    parser.add_argument(
        "--target-chapters",
        type=int,
        default=25,
        help="Target chapter count (default: 25)"
    )
    parser.add_argument(
        "--auto",
        action="store_true",
        help="Auto-detect and migrate legacy formats"
    )

    args = parser.parse_args()
    project_path = Path(args.project_path).resolve()

    if not project_path.exists():
        print(f"[ERROR] Project path does not exist: {project_path}")
        sys.exit(1)

    if args.to_writing:
        results = transition_to_writing(
            project_path,
            title=args.title,
            target_words=args.target_words,
            target_chapters=args.target_chapters
        )
        print(f"\nTransition results: {results}")
    elif args.auto:
        results = auto_migrate(project_path)
        print(f"\nMigration results: {results}")
    else:
        # Default: auto-migrate
        results = auto_migrate(project_path)
        print(f"\nMigration results: {results}")


if __name__ == "__main__":
    main()
