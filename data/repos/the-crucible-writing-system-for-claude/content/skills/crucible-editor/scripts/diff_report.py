#!/usr/bin/env python3
"""
Script: diff_report.py
Purpose: Generate change report comparing original and edited versions
Crucible Suite Plugin
"""

import sys
import json
from pathlib import Path
from datetime import datetime
import difflib

# Ensure Python 3.8+
if sys.version_info < (3, 8):
    print("Error: Python 3.8+ required", file=sys.stderr)
    sys.exit(1)


def generate_diff_report(project_path: str, chapter_num: int = None):
    """Generate a diff report between original and edited versions."""
    project_dir = Path(project_path)
    original_dir = project_dir / "original"
    edited_dir = project_dir / "edited"
    reports_dir = project_dir / "reports"

    reports_dir.mkdir(parents=True, exist_ok=True)

    if chapter_num:
        # Single chapter comparison (try ch*.md first, then chapter-*.md)
        original_file = original_dir / f"ch{chapter_num:02d}.md"
        edited_file = edited_dir / f"ch{chapter_num:02d}.md"

        if not original_file.exists() or not edited_file.exists():
            # Fallback to legacy pattern
            original_file = original_dir / f"chapter-{chapter_num:02d}.md"
            edited_file = edited_dir / f"chapter-{chapter_num:02d}.md"

        if not original_file.exists() or not edited_file.exists():
            return {"success": False, "error": f"Chapter {chapter_num} files not found"}

        chapters = [(chapter_num, original_file, edited_file)]
    else:
        # All chapters (try ch*.md first, then chapter-*.md)
        chapters = []
        edited_files = sorted(edited_dir.glob("ch*.md"))
        if not edited_files:
            edited_files = sorted(edited_dir.glob("chapter-*.md"))

        for edited_file in edited_files:
            if edited_file.stem.startswith("ch") and not edited_file.stem.startswith("chapter"):
                num = int(edited_file.stem.replace("ch", ""))
            else:
                num = int(edited_file.stem.split("-")[1])
            original_file = original_dir / edited_file.name
            if original_file.exists():
                chapters.append((num, original_file, edited_file))

    report_content = "# Edit Diff Report\n\n"
    report_content += f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}\n\n"

    total_additions = 0
    total_deletions = 0
    total_changes = 0

    for num, original_file, edited_file in chapters:
        try:
            with open(original_file, "r", encoding="utf-8") as f:
                original_lines = f.readlines()
        except (FileNotFoundError, UnicodeDecodeError) as e:
            print(f"Warning: Could not read {original_file}: {e}", file=sys.stderr)
            continue

        try:
            with open(edited_file, "r", encoding="utf-8") as f:
                edited_lines = f.readlines()
        except (FileNotFoundError, UnicodeDecodeError) as e:
            print(f"Warning: Could not read {edited_file}: {e}", file=sys.stderr)
            continue

        # Generate diff (use actual file names)
        diff = list(difflib.unified_diff(
            original_lines, edited_lines,
            fromfile=f"original/{original_file.name}",
            tofile=f"edited/{edited_file.name}",
            lineterm=""
        ))

        # Count changes
        additions = sum(1 for line in diff if line.startswith("+") and not line.startswith("+++"))
        deletions = sum(1 for line in diff if line.startswith("-") and not line.startswith("---"))

        total_additions += additions
        total_deletions += deletions
        total_changes += len(diff)

        # Word count comparison
        original_words = sum(len(line.split()) for line in original_lines)
        edited_words = sum(len(line.split()) for line in edited_lines)
        word_diff = edited_words - original_words

        report_content += f"## Chapter {num}\n\n"
        report_content += f"- **Lines added:** {additions}\n"
        report_content += f"- **Lines removed:** {deletions}\n"
        report_content += f"- **Word count:** {original_words} -> {edited_words} ({word_diff:+d})\n\n"

        if diff:
            report_content += "### Changes\n\n```diff\n"
            report_content += "\n".join(diff[:50])  # Limit to first 50 lines
            if len(diff) > 50:
                report_content += f"\n... and {len(diff) - 50} more lines\n"
            report_content += "\n```\n\n"

    # Summary
    report_content = report_content.replace(
        "# Edit Diff Report\n\n",
        f"# Edit Diff Report\n\n**Summary:** {total_additions} additions, {total_deletions} deletions across {len(chapters)} chapters\n\n"
    )

    # Write report
    report_file = reports_dir / "diff-report.md"
    with open(report_file, "w", encoding="utf-8") as f:
        f.write(report_content)

    return {
        "success": True,
        "report_file": str(report_file),
        "chapters_compared": len(chapters),
        "total_additions": total_additions,
        "total_deletions": total_deletions
    }


def main():
    # Read input
    if len(sys.argv) >= 2:
        project_path = sys.argv[1]
        chapter_num = int(sys.argv[2]) if len(sys.argv) > 2 else None
    else:
        try:
            input_data = json.load(sys.stdin)
            project_path = input_data.get("project_path", "./edit-project")
            chapter_num = input_data.get("chapter_num")
        except json.JSONDecodeError:
            project_path = "./edit-project"
            chapter_num = None

    result = generate_diff_report(project_path, chapter_num)
    print(json.dumps(result, indent=2))
    sys.exit(0 if result["success"] else 1)


if __name__ == "__main__":
    main()
