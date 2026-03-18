#!/usr/bin/env python3
"""
Script: load_project_context.py
Purpose: Load project context on session start (SessionStart hook)
Crucible Suite Plugin

This script runs at the start of a Claude Code session and provides
context about the current Crucible project state.

Note: The hook's matcher in hooks.json is set to 'startup|resume|clear|compact',
which includes 'compact' events. This is intentional - during context window
compaction, earlier context may be summarized or truncated, so re-injecting
project context ensures Crucible project awareness is maintained.
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


def get_project_context(project_root: Path, structure_type: str = None) -> str:
    """Generate context string for the session."""

    if project_root is None:
        return ""  # No project found, no context to inject

    crucible_dir = project_root / ".crucible"
    state_dir = crucible_dir / "state"

    context_parts = [
        "=== CRUCIBLE PROJECT DETECTED ===",
        ""
    ]

    # Get title from CLAUDE.md
    title = "Untitled Project"
    claude_md = project_root / "CLAUDE.md"
    if claude_md.exists():
        try:
            with open(claude_md, "r", encoding="utf-8") as f:
                content = f.read()
                for line in content.split("\n"):
                    if "Book Title:" in line:
                        title = line.split(":", 1)[1].strip()
                        break
        except (OSError, UnicodeDecodeError):
            pass

    context_parts.append(f"Project: {title}")
    context_parts.append(f"Location: {project_root}")
    context_parts.append("")

    # Determine phase and provide relevant context
    phase = "not_started"
    phase_context = []

    # Check for edit state (highest priority)
    edit_state_file = state_dir / "edit-state.json"
    if edit_state_file.exists():
        phase = "editing"
        try:
            with open(edit_state_file, "r", encoding="utf-8") as f:
                state = json.load(f)
                phase_context.append(f"Currently in: EDITING phase")
                phase_context.append(f"Edit level: {state.get('current_phase', 'unknown')}")
                chapters_edited = len(state.get("chapters", {}))
                phase_context.append(f"Chapters edited: {chapters_edited}")
        except (json.JSONDecodeError, OSError):
            pass

    # Check for draft state (new location first, then legacy)
    draft_state_file = state_dir / "draft-state.json"
    legacy_draft_state = project_root / "project-state.json"

    draft_state_to_use = None
    if draft_state_file.exists():
        draft_state_to_use = draft_state_file
    elif legacy_draft_state.exists():
        draft_state_to_use = legacy_draft_state

    if not phase_context and draft_state_to_use:
        phase = "writing"
        try:
            with open(draft_state_to_use, "r", encoding="utf-8") as f:
                state = json.load(f)
                phase_context.append(f"Currently in: WRITING phase")
                phase_context.append(f"Current chapter: {state.get('current_chapter', '?')}")
                phase_context.append(f"Current scene: {state.get('current_scene', '?')}")
                phase_context.append(f"Word count: {state.get('word_count', 0):,}")
                phase_context.append(f"Target: {state.get('target_words', 150000):,}")

                # Check if bi-chapter review is due
                # Logic must match sync_draft_state() in update_story_bible.py
                # Review is needed when 2+ chapters completed since last review
                chapters_complete = state.get('chapters_complete', 0)
                last_review = state.get('last_review_at_chapter', 0)
                chapters_since_review = chapters_complete - last_review
                if chapters_since_review >= 2:
                    review_start = last_review + 1
                    review_end = chapters_complete
                    phase_context.append("")
                    phase_context.append("[WARN] BI-CHAPTER REVIEW DUE")
                    phase_context.append(f"Review chapters {review_start}-{review_end} before continuing")
        except (json.JSONDecodeError, OSError):
            pass

    # Check for outline state
    outline_state_file = state_dir / "outline-state.json"
    if not phase_context and outline_state_file.exists():
        phase = "outlining"
        try:
            with open(outline_state_file, "r", encoding="utf-8") as f:
                state = json.load(f)
                phase_context.append(f"Currently in: OUTLINING phase")
                phase_context.append(f"Current chapter: {state.get('current_chapter', '?')}")
                phase_context.append(f"Chapters complete: {state.get('chapters_complete', 0)}/{state.get('total_chapters', '?')}")
        except (json.JSONDecodeError, OSError):
            pass

    # Check for planning state (new location first, then legacy)
    planning_state_file = state_dir / "planning-state.json"
    legacy_planning_state = project_root / "state.json"

    if not phase_context:
        state_file_to_use = None
        if planning_state_file.exists():
            state_file_to_use = planning_state_file
        elif legacy_planning_state.exists():
            state_file_to_use = legacy_planning_state

        if state_file_to_use:
            phase = "planning"
            try:
                with open(state_file_to_use, "r", encoding="utf-8") as f:
                    state = json.load(f)
                    # Check if this looks like a planning state file
                    if "progress" in state and "current_document" in state.get("progress", {}):
                        phase_context.append(f"Currently in: PLANNING phase")
                        progress = state.get("progress", {})
                        phase_context.append(f"Current document: {progress.get('current_document', '?')}")
                        docs_complete = progress.get('documents_complete', [])
                        if isinstance(docs_complete, list):
                            phase_context.append(f"Documents complete: {len(docs_complete)}/9")
                        else:
                            phase_context.append(f"Documents complete: {docs_complete}/9")
            except (json.JSONDecodeError, OSError):
                pass

    # Check if we have planning docs but no state file (at project root)
    planning_dir = project_root / "planning"
    if not phase_context and planning_dir.exists():
        docs = list(planning_dir.glob("*.md")) + list(planning_dir.glob("**/*.md"))
        if docs:
            phase_context.append("Planning documents found")
            phase_context.append(f"Planning status: Complete ({len(docs)} documents)")

    if phase_context:
        context_parts.extend(phase_context)
    else:
        context_parts.append("No active session state found.")
        context_parts.append("Use /crucible-suite:crucible-plan to start a new project")
        context_parts.append("or /crucible-suite:crucible-continue to resume an existing one.")

    context_parts.append("")
    context_parts.append("Available commands:")
    context_parts.append("  /crucible-suite:crucible-status   - Show full project status")
    context_parts.append("  /crucible-suite:crucible-continue - Resume from current point")
    context_parts.append("")
    context_parts.append("=== END CRUCIBLE CONTEXT ===")

    return "\n".join(context_parts)


def main():
    """Main entry point with top-level error handling for hook stability."""
    try:
        # Read hook input from stdin
        try:
            input_data = json.load(sys.stdin)
        except json.JSONDecodeError:
            input_data = {}

        # Find project using shared detection logic
        project_root, structure_type = find_crucible_project_with_type()

        # Generate context
        context = get_project_context(project_root, structure_type)

        # Output for SessionStart hook
        # This format allows the context to be injected into the session
        output = {
            "hookSpecificOutput": {
                "hookEventName": "SessionStart",
                "additionalContext": context
            }
        }

        print(json.dumps(output, indent=2))
        sys.exit(0)

    except Exception as e:
        # Graceful degradation: return valid hook output even on error
        # Log the error to stderr for debugging
        print(f"Error in load_project_context: {e}", file=sys.stderr)

        # Return valid but empty hook output
        output = {
            "hookSpecificOutput": {
                "hookEventName": "SessionStart",
                "additionalContext": ""
            }
        }
        print(json.dumps(output, indent=2))
        sys.exit(0)  # Exit 0 so hook doesn't block session start


if __name__ == "__main__":
    main()
