#!/usr/bin/env python3
"""
Script: compile_edited.py
Purpose: Compile edited chapters into final manuscript
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


def compile_edited_manuscript(project_path: str, output_format: str = "md"):
    """Compile all edited chapters into a single manuscript."""
    project_dir = Path(project_path)
    edited_dir = project_dir / "edited"

    if not edited_dir.exists():
        return {"success": False, "error": "No edited directory found"}

    # Load state for title
    state_file = find_state_file(project_dir)
    title = "Untitled"
    if state_file and state_file.exists():
        with open(state_file, "r", encoding="utf-8") as f:
            state = json.load(f)
            title = state.get("title", "Untitled")

    # Find all chapter files (support both ch*.md and chapter-*.md patterns)
    chapter_files = sorted(edited_dir.glob("ch*.md"),
                          key=lambda x: int(x.stem.replace("ch", "")))
    if not chapter_files:
        # Fallback to legacy chapter-*.md pattern
        chapter_files = sorted(edited_dir.glob("chapter-*.md"),
                              key=lambda x: int(x.stem.split("-")[1]))

    if not chapter_files:
        return {"success": False, "error": "No edited chapters found"}

    # Compile manuscript
    manuscript_content = f"# {title}\n\n"
    manuscript_content += f"*Edited: {datetime.now().strftime('%Y-%m-%d')}*\n\n"
    manuscript_content += "---\n\n"

    total_words = 0
    chapter_count = 0

    for chapter_file in chapter_files:
        with open(chapter_file, "r", encoding="utf-8") as f:
            content = f.read()

        manuscript_content += content
        manuscript_content += "\n\n---\n\n"

        word_count = len(content.split())
        total_words += word_count
        chapter_count += 1

    # Write compiled manuscript
    output_file = project_dir / f"edited-manuscript.{output_format}"
    with open(output_file, "w", encoding="utf-8") as f:
        f.write(manuscript_content)

    # Generate summary
    summary = {
        "success": True,
        "title": title,
        "chapters_compiled": chapter_count,
        "total_words": total_words,
        "output_file": str(output_file),
        "compiled_at": datetime.now().isoformat()
    }

    # Write summary file
    summary_file = project_dir / "reports" / "compilation-summary.json"
    summary_file.parent.mkdir(parents=True, exist_ok=True)
    with open(summary_file, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    return summary


def main():
    # Read input
    if len(sys.argv) >= 2:
        project_path = sys.argv[1]
        output_format = sys.argv[2] if len(sys.argv) > 2 else "md"
    else:
        try:
            input_data = json.load(sys.stdin)
            project_path = input_data.get("project_path", "./edit-project")
            output_format = input_data.get("output_format", "md")
        except json.JSONDecodeError:
            project_path = "./edit-project"
            output_format = "md"

    result = compile_edited_manuscript(project_path, output_format)
    print(json.dumps(result, indent=2))
    sys.exit(0 if result["success"] else 1)


if __name__ == "__main__":
    main()
