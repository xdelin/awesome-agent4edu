#!/usr/bin/env python3
"""
Initialize a new Crucible Planner project.

Usage:
    python init_project.py <project_path> <title> <premise>

Example:
    python init_project.py "./my-novel" "The Ashen Crown" "A healer discovers her power comes from stealing life"
"""

import json
import os
import sys
from datetime import datetime


def get_template_dir() -> str:
    """Get the path to the plugin's templates directory."""
    # Script is at: skills/crucible-planner/scripts/init_project.py
    # Templates are at: templates/
    script_dir = os.path.dirname(os.path.abspath(__file__))
    plugin_root = os.path.dirname(os.path.dirname(os.path.dirname(script_dir)))
    return os.path.join(plugin_root, "templates")


def copy_template(template_path: str, dest_path: str, substitutions: dict = None) -> bool:
    """
    Copy a template file to destination with optional placeholder substitution.

    Args:
        template_path: Path to the template file
        dest_path: Destination path
        substitutions: Dict of {placeholder: value} for substitution

    Returns:
        True if copied successfully, False otherwise
    """
    if not os.path.exists(template_path):
        return False

    try:
        with open(template_path, "r", encoding="utf-8") as f:
            content = f.read()

        if substitutions:
            for placeholder, value in substitutions.items():
                content = content.replace(placeholder, value)

        with open(dest_path, "w", encoding="utf-8") as f:
            f.write(content)
        return True
    except (OSError, UnicodeDecodeError):
        return False


def init_project(project_path: str, title: str, premise: str) -> dict:
    """Initialize a new Crucible planning project."""

    # Create directory structure
    directories = [
        project_path,
        os.path.join(project_path, ".crucible"),
        os.path.join(project_path, ".crucible", "state"),
        os.path.join(project_path, "planning"),
        os.path.join(project_path, "planning", "strand-maps"),
        os.path.join(project_path, "planning", "forge-points"),
    ]

    for directory in directories:
        os.makedirs(directory, exist_ok=True)

    # Initialize state
    state = {
        "project": {
            "title": title,
            "premise": premise,
            "created": datetime.now().isoformat(),
            "last_updated": datetime.now().isoformat()
        },
        "scope": {
            "target_length": None,  # "standard", "epic", "extended"
            "complexity": None,     # "single", "dual", "ensemble"
            "chapters": None
        },
        "progress": {
            "current_document": 1,
            "current_question": 1,
            "documents_complete": [],
            "total_questions_answered": 0
        },
        "answers": {
            "doc1_crucible_thesis": {},
            "doc2_quest_strand": {},
            "doc3_fire_strand": {},
            "doc4_constellation_strand": {},
            "doc5_forge_points": {
                "fp0_ignition": {},
                "fp1_first": {},
                "fp2_second": {},
                "fp3_third": {},
                "apex": {}
            },
            "doc6_dark_mirror": {},
            "doc7_constellation_bible": {
                "protagonist": {},
                "characters": []
            },
            "doc8_mercy_ledger": {
                "mercy_1": {},
                "mercy_2": {},
                "mercy_3": {},
                "mercy_4": {}
            },
            "doc9_world_forge": {}
        }
    }

    # Save state to .crucible/state/planning-state.json
    state_dir = os.path.join(project_path, ".crucible", "state")
    state_path = os.path.join(state_dir, "planning-state.json")
    with open(state_path, 'w') as f:
        json.dump(state, f, indent=2)

    # Copy CLAUDE.md templates for context management
    template_dir = get_template_dir()
    templates_copied = []

    # Root CLAUDE.md with placeholder substitution
    root_template = os.path.join(template_dir, "CLAUDE.md")
    root_dest = os.path.join(project_path, "CLAUDE.md")
    if copy_template(root_template, root_dest, {
        "[TITLE]": title,
        "[SERIES NAME] Book [#]": "Standalone",
        "[TARGET]": "100,000",
        "[CURRENT]": "0",
        "[planning|outlining|writing|editing]": "planning",
        "[#]": "1",
        "[DATE]": datetime.now().strftime("%Y-%m-%d"),
        "[Decision 1]": premise[:100] + "..." if len(premise) > 100 else premise,
        "[Decision 2]": "",
        "[Question needing author input]": "Review premise and begin planning",
    }):
        templates_copied.append("CLAUDE.md")

    # Planning directory CLAUDE.md
    planning_template = os.path.join(template_dir, "project-structure", "planning", "CLAUDE.md")
    planning_dest = os.path.join(project_path, "planning", "CLAUDE.md")
    if copy_template(planning_template, planning_dest):
        templates_copied.append("planning/CLAUDE.md")

    print(f"[OK] Initialized Crucible project: {title}")
    print(f"     Location: {os.path.abspath(project_path)}")
    print(f"     State file: {state_path}")
    if templates_copied:
        print(f"     Context templates: {', '.join(templates_copied)}")
    print(f"\nDirectory structure created:")
    print(f"     {project_path}/")
    print(f"     +-- CLAUDE.md              (project memory)")
    print(f"     +-- .crucible/")
    print(f"     |   +-- state/")
    print(f"     |       +-- planning-state.json")
    print(f"     +-- planning/")
    print(f"         +-- CLAUDE.md          (planning context)")
    print(f"         +-- strand-maps/")
    print(f"         +-- forge-points/")

    return state


def main():
    if len(sys.argv) < 4:
        print("Usage: python init_project.py <project_path> <title> <premise>")
        print("Example: python init_project.py ./my-novel \"The Ashen Crown\" \"A healer discovers...\"")
        sys.exit(1)
    
    project_path = sys.argv[1]
    title = sys.argv[2]
    premise = " ".join(sys.argv[3:])
    
    init_project(project_path, title, premise)


if __name__ == "__main__":
    main()
