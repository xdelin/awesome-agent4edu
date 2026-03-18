#!/usr/bin/env python3
"""
Save state for a Crucible Planner project.

Usage:
    # Save a single answer with full context
    python save_state.py <project_path> --answer <document> <key> <question> <answer> <description>

    # Update progress position
    python save_state.py <project_path> --progress <doc_num> <question_num>

    # Mark document complete
    python save_state.py <project_path> --complete <doc_num>

    # Set project scope
    python save_state.py <project_path> --scope <target_length> <complexity>

    # Set project metadata (title, premise)
    python save_state.py <project_path> --set <field> <value>

Can also be imported and used programmatically:
    from save_state import save_state, update_answer, mark_document_complete

Answer storage format:
    Each answer is stored as:
    {
        "question": "The full question text",
        "answer": "The selected option label",
        "description": "The option description"
    }
"""

import argparse
import json
import os
import sys
from datetime import datetime


def find_state_file(project_path: str) -> str:
    """Find the state file, checking new location first, then legacy."""
    # New location: .crucible/state/planning-state.json
    new_path = os.path.join(project_path, ".crucible", "state", "planning-state.json")
    if os.path.exists(new_path):
        return new_path

    # Legacy location: state.json at root
    legacy_path = os.path.join(project_path, "state.json")
    if os.path.exists(legacy_path):
        return legacy_path

    return None


def get_state_path(project_path: str) -> str:
    """Get the path where state should be saved (prefers new location)."""
    # If new location exists, use it
    new_path = os.path.join(project_path, ".crucible", "state", "planning-state.json")
    if os.path.exists(new_path):
        return new_path

    # If legacy exists but new doesn't, use legacy (don't migrate automatically)
    legacy_path = os.path.join(project_path, "state.json")
    if os.path.exists(legacy_path):
        return legacy_path

    # Neither exists - use new location and create directory
    state_dir = os.path.join(project_path, ".crucible", "state")
    os.makedirs(state_dir, exist_ok=True)
    return new_path


def load_state(project_path: str) -> dict:
    """Load state from project directory."""
    state_path = find_state_file(project_path)
    if state_path is None:
        raise FileNotFoundError(
            f"No state file found. Checked:\n"
            f"  - {os.path.join(project_path, '.crucible', 'state', 'planning-state.json')}\n"
            f"  - {os.path.join(project_path, 'state.json')}"
        )

    with open(state_path, 'r') as f:
        return json.load(f)


def save_state(project_path: str, state: dict) -> None:
    """Save state to project directory."""
    state["project"]["last_updated"] = datetime.now().isoformat()

    state_path = get_state_path(project_path)
    with open(state_path, 'w') as f:
        json.dump(state, f, indent=2)

    print(f"[OK] State saved: {state_path}")


def update_answer(project_path: str, document: str, key: str, question_text: str, answer: str, description: str) -> dict:
    """
    Update a single answer in the state with full context.

    Args:
        project_path: Path to the project directory
        document: Document key (e.g., "doc1_crucible_thesis" or "doc5_forge_points.fp0_ignition")
        key: The state key (e.g., "burden_type")
        question_text: The full question that was asked
        answer: The selected option label
        description: The description of the selected option

    Stores as:
        {
            "question": question_text,
            "answer": answer,
            "description": description
        }
    """
    state = load_state(project_path)

    # Create the answer object with full context
    answer_obj = {
        "question": question_text,
        "answer": answer,
        "description": description
    }

    # Navigate to the right place in state
    if document == "scope":
        # Handle scope as a special case - store at top level
        # Map common keys to scope fields
        scope_key_map = {
            "target_length": "target_length",
            "novel_length": "target_length",
            "length": "target_length",
            "complexity": "complexity",
            "narrative_complexity": "complexity",
            "pov": "complexity"
        }
        actual_key = scope_key_map.get(key, key)

        # For scope, we store just the answer value, not the full object
        state["scope"][actual_key] = answer

        # If target_length is set, calculate chapters
        if actual_key == "target_length":
            chapters_map = {
                "Standard": 22,
                "standard": 22,
                "Epic": 30,
                "epic": 30,
                "Extended/Series": 40,
                "Extended": 40,
                "extended": 40
            }
            state["scope"]["chapters"] = chapters_map.get(answer, 25)

        # Don't increment total_questions_answered for scope
        save_state(project_path, state)
        return state

    elif document.startswith("doc5_forge_points"):
        # Handle nested forge points
        parts = document.split(".")
        if len(parts) == 2:
            state["answers"]["doc5_forge_points"][parts[1]][key] = answer_obj
        else:
            state["answers"]["doc5_forge_points"][key] = answer_obj
    elif document.startswith("doc7_constellation_bible"):
        parts = document.split(".")
        if len(parts) == 2:
            if parts[1] == "protagonist":
                state["answers"]["doc7_constellation_bible"]["protagonist"][key] = answer_obj
            else:
                # Handle character additions
                state["answers"]["doc7_constellation_bible"]["characters"].append({
                    "name": parts[1],
                    key: answer_obj
                })
        else:
            state["answers"]["doc7_constellation_bible"][key] = answer_obj
    elif document.startswith("doc8_mercy_ledger"):
        parts = document.split(".")
        if len(parts) == 2:
            state["answers"]["doc8_mercy_ledger"][parts[1]][key] = answer_obj
        else:
            state["answers"]["doc8_mercy_ledger"][key] = answer_obj
    else:
        state["answers"][document][key] = answer_obj

    # Update progress
    state["progress"]["total_questions_answered"] += 1

    save_state(project_path, state)
    return state


def get_answer_value(answer_data) -> str:
    """
    Extract the answer value from answer data.
    Handles both old format (string) and new format (dict with 'answer' key).

    Args:
        answer_data: Either a string (old format) or dict (new format)

    Returns:
        The answer string value
    """
    if isinstance(answer_data, dict):
        return answer_data.get("answer", "")
    return answer_data if answer_data else ""


def update_progress(project_path: str, document: int, question: int) -> dict:
    """Update current progress position."""
    state = load_state(project_path)
    state["progress"]["current_document"] = document
    state["progress"]["current_question"] = question
    save_state(project_path, state)
    return state


def mark_document_complete(project_path: str, document_num: int) -> dict:
    """Mark a document as complete."""
    state = load_state(project_path)
    if document_num not in state["progress"]["documents_complete"]:
        state["progress"]["documents_complete"].append(document_num)
    state["progress"]["current_document"] = document_num + 1
    state["progress"]["current_question"] = 1
    save_state(project_path, state)
    print(f"[OK] Document {document_num} marked complete")
    return state


def set_scope(project_path: str, target_length: str, complexity: str) -> dict:
    """Set project scope."""
    state = load_state(project_path)
    state["scope"]["target_length"] = target_length
    state["scope"]["complexity"] = complexity

    # Calculate approximate chapters
    chapters_map = {
        "standard": 22,
        "epic": 30,
        "extended": 40
    }
    state["scope"]["chapters"] = chapters_map.get(target_length, 25)

    save_state(project_path, state)
    return state


def set_field(project_path: str, field: str, value: str) -> dict:
    """
    Set a project field (title, premise, or other metadata).

    Args:
        project_path: Path to the project directory
        field: Field name (title, premise, etc.)
        value: Value to set

    Returns:
        Updated state dict
    """
    state = load_state(project_path)

    # Map common field names to state locations
    if field in ("title", "premise"):
        state["project"][field] = value
    elif field in ("target_length", "complexity", "chapters"):
        state["scope"][field] = value
        # Recalculate chapters if length changed
        if field == "target_length":
            chapters_map = {
                "standard": 22, "Standard": 22,
                "epic": 30, "Epic": 30,
                "extended": 40, "Extended": 40
            }
            state["scope"]["chapters"] = chapters_map.get(value, 25)
    else:
        # Store in project section by default
        state["project"][field] = value

    save_state(project_path, state)
    return state


def main():
    parser = argparse.ArgumentParser(
        description='Save state for a Crucible Planner project.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Save a single answer with full context
    python save_state.py ./my-project --answer doc1_crucible_thesis burden_type "What form does the external burden take?" "Physical object" "An artifact that must be destroyed or protected"

    # Save nested answer (forge points, mercies)
    python save_state.py ./my-project --answer doc5_forge_points.fp0_ignition quest_crisis "What Quest crisis demands attention?" "Artifact stolen" "The burden is taken by enemy forces"

    # Update progress position
    python save_state.py ./my-project --progress 1 3

    # Mark document complete
    python save_state.py ./my-project --complete 1

    # Set project scope
    python save_state.py ./my-project --scope epic dual

    # Set project metadata (title, premise)
    python save_state.py ./my-project --set title "The Memory Forge"
    python save_state.py ./my-project --set premise "A blacksmith who forges weapons that steal memories"

    # Set multiple fields at once
    python save_state.py ./my-project --set title "My Novel" --set premise "The premise"
        """
    )

    parser.add_argument('project_path', help='Path to the Crucible project directory')

    # Mutually exclusive operation modes
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        '--answer',
        nargs=5,
        metavar=('DOCUMENT', 'KEY', 'QUESTION', 'ANSWER', 'DESCRIPTION'),
        help='Save answer with full context: --answer <document> <key> <question_text> <answer> <description>'
    )
    group.add_argument(
        '--progress',
        nargs=2,
        metavar=('DOC_NUM', 'QUESTION_NUM'),
        type=int,
        help='Update progress position: --progress <doc_num> <question_num>'
    )
    group.add_argument(
        '--complete',
        type=int,
        metavar='DOC_NUM',
        help='Mark a document as complete: --complete <doc_num>'
    )
    group.add_argument(
        '--scope',
        nargs=2,
        metavar=('LENGTH', 'COMPLEXITY'),
        help='Set project scope: --scope <standard|epic|extended> <single|dual|ensemble>'
    )
    group.add_argument(
        '--set',
        nargs=2,
        metavar=('FIELD', 'VALUE'),
        action='append',
        help='Set project field: --set <field> <value> (can be used multiple times)'
    )

    args = parser.parse_args()

    try:
        if args.answer:
            document, key, question_text, answer, description = args.answer
            update_answer(args.project_path, document, key, question_text, answer, description)
            print(f"[OK] Saved: {document}.{key}")
            print(f"   Q: {question_text[:60]}{'...' if len(question_text) > 60 else ''}")
            print(f"   A: {answer} - {description[:40]}{'...' if len(description) > 40 else ''}")

        elif args.progress:
            doc_num, question_num = args.progress
            update_progress(args.project_path, doc_num, question_num)
            print(f"[OK] Progress updated: Document {doc_num}, Question {question_num}")

        elif args.complete:
            mark_document_complete(args.project_path, args.complete)

        elif args.scope:
            target_length, complexity = args.scope
            set_scope(args.project_path, target_length, complexity)
            print(f"[OK] Scope set: {target_length}, {complexity}")

        elif args.set:
            # Handle multiple --set arguments
            for field, value in args.set:
                set_field(args.project_path, field, value)
                print(f"[OK] Set {field}: {value[:60]}{'...' if len(value) > 60 else ''}")

        else:
            # Default: just touch the timestamp
            state = load_state(args.project_path)
            save_state(args.project_path, state)

    except FileNotFoundError as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        sys.exit(1)
    except KeyError as e:
        print(f"[ERROR] Invalid document or key: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
