#!/usr/bin/env python3
"""
Load and display state for a Crucible Planner project.

Usage:
    python load_state.py <project_path>

Displays current progress and can be used to resume interrupted sessions.
"""

import json
import os
import sys
from datetime import datetime


DOCUMENT_NAMES = {
    1: "Crucible Thesis",
    2: "Quest Strand Map",
    3: "Fire Strand Map",
    4: "Constellation Strand Map",
    5: "Forge Point Blueprints",
    6: "Dark Mirror Profile",
    7: "Constellation Bible",
    8: "Mercy Ledger",
    9: "World Forge"
}

DOCUMENT_QUESTIONS = {
    1: 10,
    2: 7,
    3: 7,
    4: 7,
    5: 20,  # 4 questions × 5 forge points
    6: 9,
    7: 12,
    8: 16,  # 4 questions × 4 mercies
    9: 9
}


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


def display_state(state: dict) -> None:
    """Display state in a readable format."""
    project = state["project"]
    progress = state["progress"]
    scope = state["scope"]
    
    print("=" * 60)
    print(f"CRUCIBLE PROJECT: {project['title']}")
    print("=" * 60)

    print(f"\nPremise:")
    print(f"   {project['premise']}")

    print(f"\nCreated: {project['created'][:10]}")
    print(f"   Updated: {project['last_updated'][:10]}")

    if scope["target_length"]:
        print(f"\nScope:")
        print(f"   Length: {scope['target_length'].title()}")
        print(f"   Complexity: {scope['complexity'].title() if scope['complexity'] else 'Not set'}")
        print(f"   Target Chapters: {scope['chapters']}")

    print(f"\nProgress:")
    
    # Show document completion
    total_docs = 9
    complete = len(progress["documents_complete"])
    print(f"   Documents Complete: {complete}/{total_docs}")
    
    # Progress bar
    bar_length = 30
    filled = int(bar_length * complete / total_docs)
    bar = "#" * filled + "-" * (bar_length - filled)
    print(f"   [{bar}] {complete * 100 // total_docs}%")
    
    # Current position
    current_doc = progress["current_document"]
    current_q = progress["current_question"]
    
    if current_doc <= 9:
        doc_name = DOCUMENT_NAMES.get(current_doc, f"Document {current_doc}")
        total_q = DOCUMENT_QUESTIONS.get(current_doc, "?")
        print(f"\nCurrent Position:")
        print(f"   Document {current_doc}: {doc_name}")
        print(f"   Question {current_q} of {total_q}")
    else:
        print(f"\n[OK] All documents complete!")

    # Document status
    print(f"\nDocument Status:")
    for doc_num in range(1, 10):
        doc_name = DOCUMENT_NAMES.get(doc_num, f"Document {doc_num}")
        if doc_num in progress["documents_complete"]:
            status = "[x]"
        elif doc_num == current_doc:
            status = "[>]"
        else:
            status = "[ ]"
        print(f"   {status} {doc_num}. {doc_name}")

    print(f"\nTotal Questions Answered: {progress['total_questions_answered']}")
    print("=" * 60)


def get_resume_info(state: dict) -> dict:
    """Get information needed to resume the session."""
    progress = state["progress"]
    
    return {
        "document": progress["current_document"],
        "question": progress["current_question"],
        "document_name": DOCUMENT_NAMES.get(progress["current_document"], "Unknown"),
        "documents_complete": progress["documents_complete"],
        "answers": state["answers"]
    }


def main():
    if len(sys.argv) < 2:
        print("Usage: python load_state.py <project_path>")
        sys.exit(1)
    
    project_path = sys.argv[1]
    
    try:
        state = load_state(project_path)
        display_state(state)
        
        # Output resume info
        resume = get_resume_info(state)
        print(f"\nTo Resume:")
        print(f"   Continue from Document {resume['document']}: {resume['document_name']}")
        print(f"   Starting at Question {resume['question']}")

    except FileNotFoundError as e:
        print(f"[ERROR] {e}")
        print("   Have you initialized this project?")
        print("   Run: python init_project.py <path> <title> <premise>")
        sys.exit(1)


if __name__ == "__main__":
    main()
