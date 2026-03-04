#!/usr/bin/env python3
"""
Script: detect_project.py
Purpose: Detect Crucible project state for /crucible-suite:crucible-continue
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

# Import shared project detection from cross_platform
from cross_platform import find_crucible_project_with_type


def count_words_in_files(directory: Path, pattern: str = "*.md") -> int:
    """Count total words in files matching pattern."""
    total = 0
    if directory.exists():
        for file_path in directory.rglob(pattern):
            try:
                with open(file_path, "r", encoding="utf-8") as f:
                    total += len(f.read().split())
            except (OSError, UnicodeDecodeError):
                pass
    return total




def planning_files_exist(planning_dir: Path) -> bool:
    """
    Check if planning documents have actually been compiled.

    The Q&A phase marks documents_complete in state, but files aren't created
    until compile_documents.py runs in Phase 3. This function verifies the
    actual files exist.
    """
    return (
        planning_dir.exists() and
        (planning_dir / "crucible-thesis.md").exists()
    )

def detect_project_state(project_root: Path, structure_type: str = None) -> dict:
    """Detect the current state of a Crucible project."""

    if project_root is None:
        return {
            "found": False,
            "error": "No Crucible project found"
        }

    # Determine paths based on structure type
    # New unified structure: all content at project root, state in .crucible/state/
    crucible_dir = project_root / ".crucible"
    planning_dir = project_root / "planning"
    outline_dir = project_root / "outline"
    draft_dir = project_root / "draft"
    story_bible_file = project_root / "story-bible.json"

    if structure_type == "dotcrucible":
        state_dir = crucible_dir / "state"
        root_state_file = None
        # Also check for legacy state.json that might have been migrated
        legacy_state_file = crucible_dir / "state" / "planning-state.json"
    else:
        # rootlevel structure (planner-created, legacy)
        state_dir = project_root
        root_state_file = project_root / "state.json"
        legacy_state_file = None

    result = {
        "found": True,
        "project_root": str(project_root),
        "structure_type": structure_type,
        "title": "Untitled",
        "phase": "not_started",
        "progress": {}
    }

    # Load state.json for rootlevel structure
    project_state = {}
    if root_state_file and root_state_file.exists():
        try:
            with open(root_state_file, "r", encoding="utf-8") as f:
                project_state = json.load(f)
                # Get title from state
                if "project" in project_state and "title" in project_state["project"]:
                    result["title"] = project_state["project"]["title"]
        except (json.JSONDecodeError, OSError):
            pass

    # Try to get title from CLAUDE.md (fallback or override)
    claude_md = project_root / "CLAUDE.md"
    if claude_md.exists():
        try:
            with open(claude_md, "r", encoding="utf-8") as f:
                content = f.read()
                for line in content.split("\n"):
                    if "Book Title:" in line:
                        result["title"] = line.split(":", 1)[1].strip()
                        break
        except (OSError, UnicodeDecodeError):
            pass

    # Check planning state
    if structure_type in ("rootlevel", "legacy") and project_state:
        # Use state.json progress for rootlevel projects
        progress = project_state.get("progress", {})
        documents_complete = progress.get("documents_complete", [])
        current_doc = progress.get("current_document", 1)

        # Planning is only complete if BOTH:
        # 1. All 9 document Q&A sessions are done (state tracking)
        # 2. The actual files have been generated (compile_documents.py ran)
        qa_complete = len(documents_complete) >= 9
        files_exist = planning_files_exist(planning_dir)
        actually_complete = qa_complete and files_exist

        result["progress"]["planning"] = {
            "status": "complete" if actually_complete else "in_progress",
            "current_document": current_doc,
            "documents_complete": len(documents_complete),
            "documents_total": 9,
            "qa_complete": qa_complete,
            "files_generated": files_exist
        }
        if not actually_complete:
            result["phase"] = "planning"
    elif structure_type == "dotcrucible":
        planning_state_file = state_dir / "planning-state.json"
        if planning_state_file.exists():
            try:
                with open(planning_state_file, "r", encoding="utf-8") as f:
                    planning_state = json.load(f)
                    # Check both state AND actual files
                    progress_data = planning_state.get("progress", {})
                    documents_complete = progress_data.get("documents_complete", [])
                    qa_complete = len(documents_complete) >= 9 if isinstance(documents_complete, list) else False
                    files_exist = planning_files_exist(planning_dir)
                    actually_complete = qa_complete and files_exist

                    result["progress"]["planning"] = {
                        "status": "complete" if actually_complete else "in_progress",
                        "current_document": progress_data.get("current_document"),
                        "documents_complete": len(documents_complete) if isinstance(documents_complete, list) else 0,
                        "documents_total": 9,
                        "qa_complete": qa_complete,
                        "files_generated": files_exist
                    }
                    if not actually_complete:
                        result["phase"] = "planning"
            except (json.JSONDecodeError, OSError):
                pass
        elif planning_dir.exists():
            # Check for planning documents (file-based fallback)
            docs = list(planning_dir.glob("*.md")) + list(planning_dir.glob("**/*.md"))
            if docs:
                result["progress"]["planning"] = {
                    "status": "complete",
                    "documents_complete": len(docs),
                    "qa_complete": True,  # Assumed if files exist
                    "files_generated": True
                }
    else:
        # Generic file-based check
        if planning_dir.exists():
            docs = list(planning_dir.glob("*.md")) + list(planning_dir.glob("**/*.md"))
            if docs:
                result["progress"]["planning"] = {
                    "status": "complete",
                    "documents_complete": len(docs),
                    "qa_complete": True,  # Assumed if files exist
                    "files_generated": True
                }

    # Check outline state
    if structure_type == "dotcrucible":
        outline_state_file = state_dir / "outline-state.json"
        if outline_state_file.exists():
            try:
                with open(outline_state_file, "r", encoding="utf-8") as f:
                    outline_state = json.load(f)
                    result["progress"]["outline"] = {
                        "status": outline_state.get("status", "in_progress"),
                        "current_chapter": outline_state.get("current_chapter"),
                        "chapters_complete": outline_state.get("chapters_complete", 0),
                        "total_chapters": outline_state.get("total_chapters")
                    }
                    if outline_state.get("status") != "complete":
                        result["phase"] = "outlining"
            except (json.JSONDecodeError, OSError):
                pass

    # File-based outline detection - check for outline/ directory at project root
    if "outline" not in result["progress"] and outline_dir.exists():
        master_outline = outline_dir / "master-outline.md"
        result["progress"]["outline"] = {
            "status": "complete" if master_outline.exists() else "in_progress",
            "outline_dir": str(outline_dir)
        }

    # Check draft state
    if structure_type == "dotcrucible":
        draft_state_file = state_dir / "draft-state.json"
    else:
        draft_state_file = None

    if draft_state_file and draft_state_file.exists():
        try:
            with open(draft_state_file, "r", encoding="utf-8") as f:
                draft_state = json.load(f)
                result["progress"]["writing"] = {
                    "status": draft_state.get("status", "in_progress"),
                    "current_chapter": draft_state.get("current_chapter"),
                    "current_scene": draft_state.get("current_scene"),
                    "chapters_complete": draft_state.get("chapters_complete", 0),
                    "total_chapters": draft_state.get("total_chapters"),
                    "word_count": draft_state.get("word_count", 0),
                    "target_words": draft_state.get("target_words")
                }
                if draft_state.get("status") != "complete":
                    result["phase"] = "writing"
        except (json.JSONDecodeError, OSError):
            pass
    # File-based draft detection for both structures
    chapters_dir = draft_dir / "chapters"
    if "writing" not in result["progress"] and (draft_dir.exists() or chapters_dir.exists()):
        # Check for chapter files in draft/chapters/ (ch01.md, ch02.md, etc.)
        search_dir = chapters_dir if chapters_dir.exists() else draft_dir
        chapters = list(search_dir.glob("ch*.md")) + list(search_dir.glob("chapter*.md"))
        word_count = count_words_in_files(search_dir, "*.md")
        if chapters or word_count > 0:
            result["progress"]["writing"] = {
                "status": "in_progress" if chapters else "not_started",
                "chapters_complete": len(chapters),
                "word_count": word_count
            }
            if chapters:
                result["phase"] = "writing"

    # Check edit state (only for dotcrucible structure)
    if structure_type == "dotcrucible":
        edit_state_file = state_dir / "edit-state.json"
        if edit_state_file.exists():
            try:
                with open(edit_state_file, "r", encoding="utf-8") as f:
                    edit_state = json.load(f)
                    result["progress"]["editing"] = {
                        "status": edit_state.get("current_phase", "in_progress"),
                        "chapters_edited": len(edit_state.get("chapters", {}))
                    }
                    result["phase"] = "editing"
            except (json.JSONDecodeError, OSError):
                pass

    # Determine resume point
    result["resume"] = get_resume_point(result)

    return result


def get_resume_point(state: dict) -> dict:
    """Determine the best point to resume from."""
    phase = state.get("phase", "not_started")
    progress = state.get("progress", {})

    if phase == "not_started":
        return {
            "action": "start_planning",
            "message": "No progress found. Start with /crucible-suite:crucible-plan"
        }

    if phase == "planning":
        planning = progress.get("planning", {})
        # Check if Q&A is done but files aren't generated
        if planning.get("qa_complete") and not planning.get("files_generated"):
            return {
                "action": "compile_planning",
                "current_document": planning.get("current_document"),
                "message": "Q&A complete (9/9 docs). Run compile_documents.py to generate planning files."
            }
        return {
            "action": "continue_planning",
            "current_document": planning.get("current_document"),
            "message": f"Resume planning: {planning.get('documents_complete', 0)}/9 documents complete"
        }

    if phase == "outlining":
        outline = progress.get("outline", {})
        return {
            "action": "continue_outlining",
            "current_chapter": outline.get("current_chapter"),
            "message": f"Resume outlining: Chapter {outline.get('current_chapter', '?')}"
        }

    if phase == "writing":
        writing = progress.get("writing", {})
        return {
            "action": "continue_writing",
            "current_chapter": writing.get("current_chapter"),
            "current_scene": writing.get("current_scene"),
            "word_count": writing.get("word_count", 0),
            "message": f"Resume writing: Chapter {writing.get('current_chapter', '?')}, Scene {writing.get('current_scene', '?')}"
        }

    if phase == "editing":
        editing = progress.get("editing", {})
        return {
            "action": "continue_editing",
            "message": f"Resume editing: {editing.get('chapters_edited', 0)} chapters edited"
        }

    return {
        "action": "unknown",
        "message": "Unable to determine resume point"
    }


def main():
    # Read input
    start_path = None

    if len(sys.argv) >= 2:
        start_path = Path(sys.argv[1])
    else:
        try:
            input_data = json.load(sys.stdin)
            if "path" in input_data:
                start_path = Path(input_data["path"])
        except json.JSONDecodeError:
            pass

    project_root, structure_type = find_crucible_project_with_type(start_path)
    result = detect_project_state(project_root, structure_type)

    print(json.dumps(result, indent=2))
    sys.exit(0 if result.get("found") else 1)


if __name__ == "__main__":
    main()
