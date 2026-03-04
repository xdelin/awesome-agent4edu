#!/usr/bin/env python3
"""
Script: status_reporter.py
Purpose: Generate comprehensive project status report for /crucible-suite:crucible-status
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


def count_words_in_file(file_path: Path) -> int:
    """Count words in a single file."""
    try:
        with open(file_path, "r", encoding="utf-8") as f:
            return len(f.read().split())
    except (OSError, UnicodeDecodeError):
        return 0


def get_last_modified(directory: Path) -> str:
    """Get the most recent modification time from files in a directory."""
    latest = None
    if directory.exists():
        for file_path in directory.rglob("*"):
            if file_path.is_file():
                mtime = file_path.stat().st_mtime
                if latest is None or mtime > latest:
                    latest = mtime

    if latest:
        return datetime.fromtimestamp(latest).strftime("%Y-%m-%d %H:%M")
    return "Unknown"


def generate_status_report(project_root: Path, structure_type: str = None) -> dict:
    """Generate a comprehensive status report."""

    if project_root is None:
        return {
            "success": False,
            "error": "No Crucible project found"
        }

    # Determine paths based on structure type
    # New unified structure: all content at project root, state in .crucible/state/
    crucible_dir = project_root / ".crucible"
    planning_dir = project_root / "planning"
    outline_dir = project_root / "outline"
    draft_dir = project_root / "draft"
    story_bible_file = project_root / "story-bible.json"
    backup_dir = crucible_dir / "backups"

    if structure_type == "dotcrucible":
        state_dir = crucible_dir / "state"
        state_file = None  # Uses individual state files in state_dir
    else:
        # rootlevel structure (planner-created, legacy)
        state_dir = project_root  # state.json at root
        state_file = project_root / "state.json"

    report = {
        "success": True,
        "project_root": str(project_root),
        "structure_type": structure_type,
        "title": "Untitled",
        "phase": "not_started",
        "overall_progress": 0,
        "planning": {},
        "outline": {},
        "writing": {},
        "editing": {},
        "backups": {},
        "recent_files": []
    }

    # Load state.json for rootlevel structure
    project_state = {}
    if state_file and state_file.exists():
        try:
            with open(state_file, "r", encoding="utf-8") as f:
                project_state = json.load(f)
                # Get title from state
                if "project" in project_state and "title" in project_state["project"]:
                    report["title"] = project_state["project"]["title"]
        except (json.JSONDecodeError, OSError):
            pass

    # Get title from CLAUDE.md (fallback or override)
    claude_md = project_root / "CLAUDE.md"
    if claude_md.exists():
        try:
            with open(claude_md, "r", encoding="utf-8") as f:
                content = f.read()
                for line in content.split("\n"):
                    if "Book Title:" in line:
                        report["title"] = line.split(":", 1)[1].strip()
                        break
        except (OSError, UnicodeDecodeError):
            pass

    # Planning status
    # Map document numbers to display names (9 documents total)
    planning_doc_map = {
        1: "Thesis",
        2: "Quest Strand",
        3: "Fire Strand",
        4: "Constellation Strand",
        5: "Forge Points",
        6: "Dark Mirror",
        7: "Constellation Bible",
        8: "Mercy Ledger",
        9: "World Forge"
    }

    # File-based detection (for compiled documents)
    file_docs = {
        "crucible-thesis.md": "Thesis",
        "strand-maps": "Strand Maps",
        "forge-points": "Forge Points",
        "dark-mirror-profile.md": "Dark Mirror",
        "constellation-bible.md": "Constellation Bible",
        "mercy-ledger.md": "Mercy Ledger",
        "world-forge.md": "World Forge"
    }

    planning_complete = 0
    planning_total = 9  # Total planning documents
    planning_status = {}

    # Check for rootlevel state.json progress
    if structure_type in ("rootlevel", "legacy") and project_state:
        progress = project_state.get("progress", {})
        documents_complete = progress.get("documents_complete", [])
        current_doc = progress.get("current_document", 1)

        # Map document numbers to state.json document keys
        doc_num_to_key = {
            1: "doc1_crucible_thesis",
            2: "doc2_quest_strand",
            3: "doc3_fire_strand",
            4: "doc4_constellation_strand",
            5: "doc5_forge_points",
            6: "doc6_dark_mirror",
            7: "doc7_constellation_bible",
            8: "doc8_mercy_ledger",
            9: "doc9_world_forge"
        }

        for doc_num, name in planning_doc_map.items():
            doc_key = doc_num_to_key.get(doc_num)
            if doc_key and doc_key in documents_complete:
                planning_status[name] = "complete"
                planning_complete += 1
            elif doc_num == current_doc:
                planning_status[name] = "in_progress"
            else:
                planning_status[name] = "pending"
    else:
        # File-based detection
        for doc, name in file_docs.items():
            doc_path = planning_dir / doc
            if doc_path.exists():
                planning_status[name] = "complete"
                planning_complete += 1
            else:
                planning_status[name] = "pending"
        planning_total = len(file_docs)

    report["planning"] = {
        "complete": planning_complete,
        "total": planning_total,
        "percentage": int((planning_complete / planning_total) * 100) if planning_total > 0 else 0,
        "documents": planning_status,
        "last_modified": get_last_modified(planning_dir),
        "current_document": project_state.get("progress", {}).get("current_document") if structure_type in ("rootlevel", "legacy") else None
    }

    # Outline status - check for outline/ directory at project root
    master_outline = outline_dir / "master-outline.md"
    outline_exists = outline_dir.exists() and master_outline.exists()

    # Try to get total from state or default
    total_chapters = 25  # Default
    if structure_type in ("rootlevel", "legacy"):
        # Get from project_state scope
        total_chapters = project_state.get("scope", {}).get("chapters", 25) or 25
    else:
        outline_state_file = state_dir / "outline-state.json"
        if outline_state_file.exists():
            try:
                with open(outline_state_file, "r", encoding="utf-8") as f:
                    state = json.load(f)
                    total_chapters = state.get("total_chapters", 25)
            except (json.JSONDecodeError, OSError):
                pass

    report["outline"] = {
        "complete": total_chapters if outline_exists else 0,
        "total": total_chapters,
        "percentage": 100 if outline_exists else 0,
        "last_modified": get_last_modified(outline_dir) if outline_exists else "N/A"
    }

    # Writing status
    chapters_dir = draft_dir / "chapters"

    draft_chapters = []
    if chapters_dir.exists():
        # Look for ch01.md, ch02.md, etc. pattern
        draft_chapters = sorted(chapters_dir.glob("ch*.md"))
        if not draft_chapters:
            draft_chapters = sorted(chapters_dir.glob("chapter*.md"))
    elif draft_dir.exists():
        draft_chapters = sorted(draft_dir.glob("ch*.md"))
        if not draft_chapters:
            draft_chapters = sorted(draft_dir.glob("chapter*.md"))

    # Prefer story-bible.json for word count and chapter count (authoritative source)
    total_words = 0
    story_bible_word_count = 0
    chapters_complete_from_bible = 0
    if story_bible_file.exists():
        try:
            with open(story_bible_file, "r", encoding="utf-8") as f:
                story_bible = json.load(f)
                story_bible_word_count = story_bible.get("progress", {}).get("total_words", 0)
                # Count completed chapters from story bible
                chapters = story_bible.get("chapters", {})
                chapters_complete_from_bible = sum(
                    1 for k, ch in chapters.items()
                    if not k.startswith("_") and ch.get("completed", False)
                )
        except (json.JSONDecodeError, OSError):
            pass

    if story_bible_word_count > 0:
        total_words = story_bible_word_count
    else:
        for chapter in draft_chapters:
            total_words += count_words_in_file(chapter)

    # Use story bible chapter count if available, otherwise file count
    chapters_complete = chapters_complete_from_bible if chapters_complete_from_bible > 0 else len(draft_chapters)

    # Get target words from state
    target_words = 150000  # Default
    current_chapter = None
    current_scene = None

    if structure_type == "dotcrucible":
        draft_state_file = state_dir / "draft-state.json"
    else:
        draft_state_file = None  # rootlevel doesn't use separate draft state file yet

    if draft_state_file and draft_state_file.exists():
        try:
            with open(draft_state_file, "r", encoding="utf-8") as f:
                state = json.load(f)
                target_words = state.get("target_words", 150000)
                current_chapter = state.get("current_chapter")
                current_scene = state.get("current_scene")
        except (json.JSONDecodeError, OSError):
            pass

    report["writing"] = {
        "chapters_complete": chapters_complete,
        "total_chapters": total_chapters,
        "current_chapter": current_chapter,
        "current_scene": current_scene,
        "word_count": total_words,
        "target_words": target_words,
        "percentage": int((total_words / target_words) * 100) if target_words > 0 else 0,
        "last_modified": get_last_modified(draft_dir)
    }

    # Editing status
    editing_info = {
        "started": False,
        "chapters_edited": 0,
        "current_phase": "not_started"
    }

    if structure_type == "dotcrucible":
        edit_state_file = state_dir / "edit-state.json"
        if edit_state_file.exists():
            try:
                with open(edit_state_file, "r", encoding="utf-8") as f:
                    state = json.load(f)
                    editing_info["started"] = True
                    editing_info["chapters_edited"] = len(state.get("chapters", {}))
                    editing_info["current_phase"] = state.get("current_phase", "assessment")
            except (json.JSONDecodeError, OSError):
                pass

    report["editing"] = editing_info

    # Backup status
    manifest_file = backup_dir / "backup-manifest.json"

    backup_info = {
        "count": 0,
        "latest": None
    }

    if manifest_file.exists():
        try:
            with open(manifest_file, "r", encoding="utf-8") as f:
                manifest = json.load(f)
                backup_info["count"] = len(manifest)
                if manifest:
                    backup_info["latest"] = manifest[-1].get("timestamp")
        except (json.JSONDecodeError, OSError):
            pass

    report["backups"] = backup_info

    # Determine current phase
    if report["editing"]["started"]:
        report["phase"] = "editing"
    elif report["writing"]["word_count"] > 0:
        report["phase"] = "writing"
    elif report["outline"]["complete"] > 0:
        report["phase"] = "outlining"
    elif report["planning"]["complete"] > 0 or report["planning"].get("current_document"):
        report["phase"] = "planning"
    else:
        report["phase"] = "not_started"

    # Calculate overall progress (weighted)
    # Planning: 10%, Outline: 15%, Writing: 60%, Editing: 15%
    overall = (
        (report["planning"]["percentage"] * 0.10) +
        (report["outline"]["percentage"] * 0.15) +
        (report["writing"]["percentage"] * 0.60) +
        (report["editing"].get("percentage", 0) * 0.15)
    )
    report["overall_progress"] = int(overall)

    return report


def format_report_text(report: dict) -> str:
    """Format the report as beautiful text output."""
    if not report.get("success"):
        return f"Error: {report.get('error', 'Unknown error')}"

    # ASCII-safe box drawing characters
    H, V = "-", "|"
    CORNER = "+"

    width = 58
    inner = width - 4

    def box_line(char="-"):
        return f"+{char * (width - 2)}+"

    def text_line(text, align="left"):
        if align == "center":
            padded = text.center(inner)
        else:
            padded = text.ljust(inner)
        return f"| {padded} |"

    def progress_bar(percent, bar_width=24):
        filled = int(bar_width * percent / 100)
        empty = bar_width - filled
        return f"[{'#' * filled}{'-' * empty}]"

    def section_header(title):
        pad_total = width - len(title) - 4
        pad_left = pad_total // 2
        pad_right = pad_total - pad_left
        return f"+{'-' * pad_left} {title} {'-' * pad_right}+"

    # Build the report
    lines = []

    # Title box
    lines.append(box_line("="))
    lines.append(text_line(""))
    lines.append(text_line("~ THE CRUCIBLE ~", "center"))
    lines.append(text_line(report['title'], "center"))
    lines.append(text_line(""))
    lines.append(box_line("="))

    # Phase and progress
    phase_display = report['phase'].replace('_', ' ').upper()
    lines.append(text_line(""))
    lines.append(text_line(f"  Phase:    {phase_display}"))
    lines.append(text_line(f"  Progress: {progress_bar(report['overall_progress'])} {report['overall_progress']:3d}%"))
    lines.append(text_line(""))

    # Planning section
    lines.append(section_header("PLANNING"))
    lines.append(text_line(""))
    plan_pct = report['planning']['percentage']
    lines.append(text_line(f"  {progress_bar(plan_pct)} {plan_pct:3d}%  ({report['planning']['complete']}/{report['planning']['total']})"))
    lines.append(text_line(""))

    for doc, status in report["planning"]["documents"].items():
        if status == "complete":
            icon = "[x]"
        elif status == "in_progress":
            icon = "[>]"
        else:
            icon = "[ ]"
        lines.append(text_line(f"    {icon} {doc}"))

    lines.append(text_line(""))

    # Outline section
    lines.append(section_header("OUTLINE"))
    lines.append(text_line(""))
    out_pct = report['outline']['percentage']
    lines.append(text_line(f"  {progress_bar(out_pct)} {out_pct:3d}%  ({report['outline']['complete']}/{report['outline']['total']} ch)"))
    lines.append(text_line(""))

    # Writing section
    lines.append(section_header("WRITING"))
    lines.append(text_line(""))
    write_pct = report['writing']['percentage']
    lines.append(text_line(f"  {progress_bar(write_pct)} {write_pct:3d}%"))
    lines.append(text_line(""))
    lines.append(text_line(f"    Chapters:  {report['writing']['chapters_complete']}/{report['writing']['total_chapters']}"))
    lines.append(text_line(f"    Words:     {report['writing']['word_count']:,} / {report['writing']['target_words']:,}"))

    if report["writing"]["current_chapter"]:
        lines.append(text_line(f"    Current:   Ch {report['writing']['current_chapter']}, Scene {report['writing'].get('current_scene', '?')}"))
    lines.append(text_line(""))

    # Editing section
    lines.append(section_header("EDITING"))
    lines.append(text_line(""))
    if report['editing']['started']:
        lines.append(text_line(f"    Status:    {report['editing']['current_phase'].replace('_', ' ').title()}"))
        lines.append(text_line(f"    Chapters:  {report['editing']['chapters_edited']} edited"))
    else:
        lines.append(text_line("    Status:    Not started"))
    lines.append(text_line(""))

    # Footer
    lines.append(section_header("BACKUP"))
    lines.append(text_line(""))
    backup_text = report['backups']['latest'] or 'Never'
    lines.append(text_line(f"    Last:      {backup_text}"))
    lines.append(text_line(""))
    lines.append(box_line("="))

    return "\n".join(lines)


def main():
    # Read input
    start_path = None
    output_format = "json"

    if len(sys.argv) >= 2:
        start_path = Path(sys.argv[1])
        if len(sys.argv) >= 3:
            output_format = sys.argv[2]
    else:
        try:
            input_data = json.load(sys.stdin)
            if "path" in input_data:
                start_path = Path(input_data["path"])
            output_format = input_data.get("format", "json")
        except json.JSONDecodeError:
            pass

    project_root, structure_type = find_crucible_project_with_type(start_path)
    report = generate_status_report(project_root, structure_type)

    if output_format == "text":
        print(format_report_text(report))
    else:
        print(json.dumps(report, indent=2))

    sys.exit(0 if report.get("success") else 1)


if __name__ == "__main__":
    main()
