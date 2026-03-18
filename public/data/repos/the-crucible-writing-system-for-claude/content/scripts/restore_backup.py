#!/usr/bin/env python3
"""
Script: restore_backup.py
Purpose: Restore Crucible project from backups
Crucible Suite Plugin

Supports:
- Listing available backups with timestamps and sizes
- Restoring from full zip backups
- Restoring individual files from incremental backups
- Creating safety backup before any restore
"""

import sys
import json
import os
import shutil
import zipfile
from pathlib import Path
from datetime import datetime
from typing import Optional, List, Dict, Any

# Ensure Python 3.8+
if sys.version_info < (3, 8):
    print("Error: Python 3.8+ required", file=sys.stderr)
    sys.exit(1)

# Import shared utilities
try:
    from cross_platform import (
        find_crucible_project_with_type, ensure_directory, get_timestamp,
        safe_read_json, format_output, get_backup_directory,
        get_backup_paths_for_structure
    )
except ImportError:
    # Fallback if not in path
    def find_crucible_project_with_type(start_path=None):
        current = Path.cwd() if start_path is None else Path(start_path)
        for directory in [current] + list(current.parents):
            if (directory / ".crucible").exists():
                return directory, "dotcrucible"
            if (directory / "state.json").exists():
                return directory, "rootlevel"
            if (directory / "story-bible.json").exists():
                return directory, "legacy"
            if (directory / "planning").exists():
                return directory, "legacy"
        return None, None

    def ensure_directory(path):
        path.mkdir(parents=True, exist_ok=True)
        return True

    def get_timestamp():
        return datetime.now().strftime("%Y%m%d-%H%M%S")

    def safe_read_json(path):
        try:
            with open(path, "r", encoding="utf-8") as f:
                return json.load(f)
        except (json.JSONDecodeError, FileNotFoundError, OSError):
            return None

    def format_output(data, for_hook=False):
        return json.dumps(data, indent=2)

    def get_backup_directory(project_root, structure_type):
        if structure_type == "dotcrucible":
            backup_dir = project_root / ".crucible" / "backups"
        else:
            backup_dir = project_root / ".crucible-backups"
        ensure_directory(backup_dir)
        return backup_dir

    def get_backup_paths_for_structure(project_root, structure_type):
        dirs = []
        files = []
        for dirname in ["planning", "outline", "draft", "manuscript", "chapters"]:
            dir_path = project_root / dirname
            if dir_path.exists():
                dirs.append(dir_path)
        for filename in ["CLAUDE.md", "story-bible.json", "style-profile.json"]:
            file_path = project_root / filename
            if file_path.exists():
                files.append(file_path)
        if structure_type == "dotcrucible":
            crucible_dir = project_root / ".crucible"
            if crucible_dir.exists():
                state_dir = crucible_dir / "state"
                if state_dir.exists():
                    dirs.append(state_dir)
        elif structure_type == "rootlevel":
            state_file = project_root / "state.json"
            if state_file.exists():
                files.append(state_file)
        return {"dirs": dirs, "files": files}


def list_backups(project_root: Path, structure_type: str = None) -> List[Dict[str, Any]]:
    """
    List all available backups for a project.

    Args:
        project_root: Path to the Crucible project root
        structure_type: Project structure type (dotcrucible, rootlevel, legacy)

    Returns:
        List of backup info dictionaries with timestamp, path, size, type
    """
    backups = []

    # Auto-detect structure type if not provided
    if structure_type is None:
        _, structure_type = find_crucible_project_with_type(project_root)
        if structure_type is None:
            structure_type = "dotcrucible"

    backup_dir = get_backup_directory(project_root, structure_type)
    incremental_dir = backup_dir / "incremental"

    # List full backups (zip files)
    if backup_dir.exists():
        for zip_file in sorted(backup_dir.glob("crucible-backup-*.zip"), reverse=True):
            try:
                stat = zip_file.stat()
                # Extract timestamp from filename: crucible-backup-YYYYMMDD-HHMMSS.zip
                name_parts = zip_file.stem.split("-")
                if len(name_parts) >= 4:
                    timestamp_str = f"{name_parts[2]}-{name_parts[3]}"
                else:
                    timestamp_str = "unknown"

                backups.append({
                    "type": "full",
                    "path": str(zip_file),
                    "filename": zip_file.name,
                    "timestamp": timestamp_str,
                    "size_bytes": stat.st_size,
                    "size_human": _format_size(stat.st_size),
                    "modified": datetime.fromtimestamp(stat.st_mtime).isoformat()
                })
            except OSError:
                continue

    # List incremental backups
    if incremental_dir.exists():
        for backup_file in sorted(incremental_dir.glob("*"), reverse=True):
            if backup_file.is_file():
                try:
                    stat = backup_file.stat()
                    backups.append({
                        "type": "incremental",
                        "path": str(backup_file),
                        "filename": backup_file.name,
                        "timestamp": datetime.fromtimestamp(stat.st_mtime).strftime("%Y%m%d-%H%M%S"),
                        "size_bytes": stat.st_size,
                        "size_human": _format_size(stat.st_size),
                        "modified": datetime.fromtimestamp(stat.st_mtime).isoformat()
                    })
                except OSError:
                    continue

    return backups


def _format_size(size_bytes: int) -> str:
    """Format byte size as human-readable string."""
    for unit in ["B", "KB", "MB", "GB"]:
        if size_bytes < 1024:
            return f"{size_bytes:.1f} {unit}"
        size_bytes /= 1024
    return f"{size_bytes:.1f} TB"


def create_pre_restore_backup(project_root: Path, structure_type: str = None) -> Dict[str, Any]:
    """
    Create a safety backup before performing a restore.

    Args:
        project_root: Path to the Crucible project root
        structure_type: Project structure type (dotcrucible, rootlevel, legacy)

    Returns:
        Result dictionary with success status and backup path
    """
    # Auto-detect structure type if not provided
    if structure_type is None:
        _, structure_type = find_crucible_project_with_type(project_root)
        if structure_type is None:
            structure_type = "dotcrucible"

    base_backup_dir = get_backup_directory(project_root, structure_type)
    backup_dir = base_backup_dir / "pre-restore"
    ensure_directory(backup_dir)

    timestamp = get_timestamp()
    backup_name = f"pre-restore-{timestamp}.zip"
    backup_path = backup_dir / backup_name

    # Get directories and files to backup based on structure type
    backup_paths = get_backup_paths_for_structure(project_root, structure_type)
    backup_targets = backup_paths["dirs"]
    backup_files = backup_paths["files"]

    files_backed_up = 0
    total_size = 0

    try:
        with zipfile.ZipFile(backup_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
            for target in backup_targets:
                if target.exists():
                    for file_path in target.rglob("*"):
                        if file_path.is_file():
                            arcname = file_path.relative_to(project_root)
                            zipf.write(file_path, arcname)
                            files_backed_up += 1
                            total_size += file_path.stat().st_size

            for file_path in backup_files:
                if file_path.exists():
                    arcname = file_path.relative_to(project_root)
                    zipf.write(file_path, arcname)
                    files_backed_up += 1
                    total_size += file_path.stat().st_size

        return {
            "success": True,
            "backup_path": str(backup_path),
            "files_backed_up": files_backed_up,
            "size_bytes": total_size,
            "message": f"Pre-restore backup created: {backup_name}"
        }

    except Exception as e:
        return {
            "success": False,
            "error": f"Failed to create pre-restore backup: {str(e)}"
        }


def restore_backup(backup_path: Path, project_root: Path, dry_run: bool = False) -> Dict[str, Any]:
    """
    Restore a project from a full zip backup.

    Args:
        backup_path: Path to the backup zip file
        project_root: Path to the Crucible project root
        dry_run: If True, only show what would be restored

    Returns:
        Result dictionary with success status and details
    """
    backup_path = Path(backup_path)
    project_root = Path(project_root)

    if not backup_path.exists():
        return {
            "success": False,
            "error": f"Backup file not found: {backup_path}"
        }

    if not zipfile.is_zipfile(backup_path):
        return {
            "success": False,
            "error": f"Not a valid zip file: {backup_path}"
        }

    files_to_restore = []
    try:
        with zipfile.ZipFile(backup_path, 'r') as zipf:
            files_to_restore = zipf.namelist()

            if dry_run:
                return {
                    "success": True,
                    "dry_run": True,
                    "files_count": len(files_to_restore),
                    "files": files_to_restore[:50],  # Limit preview to 50 files
                    "truncated": len(files_to_restore) > 50,
                    "message": f"Would restore {len(files_to_restore)} files from {backup_path.name}"
                }

            # Create pre-restore backup for safety
            pre_backup = create_pre_restore_backup(project_root)
            if not pre_backup["success"]:
                return {
                    "success": False,
                    "error": f"Failed to create safety backup: {pre_backup.get('error', 'Unknown error')}"
                }

            # Extract all files
            zipf.extractall(project_root)

        return {
            "success": True,
            "restored_from": str(backup_path),
            "files_restored": len(files_to_restore),
            "pre_restore_backup": pre_backup["backup_path"],
            "message": f"Restored {len(files_to_restore)} files from {backup_path.name}"
        }

    except zipfile.BadZipFile as e:
        return {
            "success": False,
            "error": f"Corrupted backup file: {str(e)}"
        }
    except Exception as e:
        return {
            "success": False,
            "error": f"Restore failed: {str(e)}"
        }


def restore_incremental(
    file_pattern: str,
    project_root: Path,
    structure_type: str = None,
    dry_run: bool = False
) -> Dict[str, Any]:
    """
    Restore specific files from incremental backups.

    Args:
        file_pattern: Glob pattern to match backup files (e.g., "chapter-01*")
        project_root: Path to the Crucible project root
        structure_type: Project structure type (dotcrucible, rootlevel, legacy)
        dry_run: If True, only show what would be restored

    Returns:
        Result dictionary with success status and details
    """
    project_root = Path(project_root)

    # Auto-detect structure type if not provided
    if structure_type is None:
        _, structure_type = find_crucible_project_with_type(project_root)
        if structure_type is None:
            structure_type = "dotcrucible"

    backup_base = get_backup_directory(project_root, structure_type)
    incremental_dir = backup_base / "incremental"

    if not incremental_dir.exists():
        return {
            "success": False,
            "error": "No incremental backups found"
        }

    # Find matching backup files
    matching_files = list(incremental_dir.glob(file_pattern))
    if not matching_files:
        return {
            "success": False,
            "error": f"No backups matching pattern: {file_pattern}"
        }

    # Sort by modification time (newest first)
    matching_files.sort(key=lambda p: p.stat().st_mtime, reverse=True)

    if dry_run:
        return {
            "success": True,
            "dry_run": True,
            "matching_files": [
                {
                    "path": str(f),
                    "filename": f.name,
                    "size": _format_size(f.stat().st_size),
                    "modified": datetime.fromtimestamp(f.stat().st_mtime).isoformat()
                }
                for f in matching_files[:20]
            ],
            "count": len(matching_files),
            "message": f"Found {len(matching_files)} matching backups"
        }

    # Get the most recent matching file
    most_recent = matching_files[0]

    # Determine restore destination based on filename
    # Incremental backups are named like: chapter-01-scene-02-20241218-120000.md
    # We need to figure out where to restore it
    filename = most_recent.name
    # Remove timestamp suffix (last two dash-separated segments before extension)
    parts = filename.rsplit("-", 2)
    if len(parts) >= 3:
        original_name = parts[0] + most_recent.suffix
    else:
        original_name = filename

    # Common restore locations
    possible_destinations = [
        project_root / "draft" / original_name,
        project_root / "outline" / original_name,
        project_root / "manuscript" / original_name,
    ]

    # Create pre-restore backup
    pre_backup = create_pre_restore_backup(project_root)
    if not pre_backup["success"]:
        return {
            "success": False,
            "error": f"Failed to create safety backup: {pre_backup.get('error', 'Unknown error')}"
        }

    # Try to find existing file to replace, or use draft as default
    restore_dest = project_root / "draft" / original_name
    for dest in possible_destinations:
        if dest.exists():
            restore_dest = dest
            break

    try:
        ensure_directory(restore_dest.parent)
        shutil.copy2(most_recent, restore_dest)

        return {
            "success": True,
            "restored_from": str(most_recent),
            "restored_to": str(restore_dest),
            "pre_restore_backup": pre_backup["backup_path"],
            "alternatives_available": len(matching_files) - 1,
            "message": f"Restored {original_name} from incremental backup"
        }

    except Exception as e:
        return {
            "success": False,
            "error": f"Failed to restore file: {str(e)}"
        }


def selective_restore(
    backup_path: Path,
    project_root: Path,
    scope: str,
    structure_type: str = None,
    dry_run: bool = False
) -> Dict[str, Any]:
    """
    Selectively restore specific types of files from a backup.

    Args:
        backup_path: Path to the backup zip file
        project_root: Path to the Crucible project root
        scope: What to restore - "chapter", "planning", "story-bible", "state", or "all"
        structure_type: Project structure type (dotcrucible, rootlevel, legacy)
        dry_run: If True, only show what would be restored

    Returns:
        Result dictionary with success status and details
    """
    backup_path = Path(backup_path)
    project_root = Path(project_root)

    if not backup_path.exists():
        return {
            "success": False,
            "error": f"Backup file not found: {backup_path}"
        }

    if not zipfile.is_zipfile(backup_path):
        return {
            "success": False,
            "error": f"Not a valid zip file: {backup_path}"
        }

    # Auto-detect structure type if not provided
    if structure_type is None:
        _, structure_type = find_crucible_project_with_type(project_root)
        if structure_type is None:
            structure_type = "dotcrucible"

    # Define scope filters
    scope_filters = {
        "chapter": ["draft/", "chapters/", "manuscript/"],
        "planning": ["planning/"],
        "story-bible": ["story-bible.json", ".crucible/story-bible/", "story-bible/"],
        "state": [".crucible/state/", "state.json", "state/"],
        "outline": ["outline/"],
        "all": None  # No filter, restore everything
    }

    if scope not in scope_filters:
        return {
            "success": False,
            "error": f"Invalid scope: {scope}. Valid options: {', '.join(scope_filters.keys())}"
        }

    filter_patterns = scope_filters[scope]

    try:
        with zipfile.ZipFile(backup_path, 'r') as zipf:
            all_files = zipf.namelist()

            # Filter files based on scope
            if filter_patterns is None:
                files_to_restore = all_files
            else:
                files_to_restore = []
                for f in all_files:
                    for pattern in filter_patterns:
                        if f.startswith(pattern) or f == pattern.rstrip('/'):
                            files_to_restore.append(f)
                            break

            if not files_to_restore:
                return {
                    "success": False,
                    "error": f"No files matching scope '{scope}' found in backup"
                }

            if dry_run:
                return {
                    "success": True,
                    "dry_run": True,
                    "scope": scope,
                    "files_count": len(files_to_restore),
                    "files": files_to_restore[:50],
                    "truncated": len(files_to_restore) > 50,
                    "message": f"Would restore {len(files_to_restore)} {scope} files from {backup_path.name}"
                }

            # Create pre-restore backup for safety
            pre_backup = create_pre_restore_backup(project_root, structure_type)
            if not pre_backup["success"]:
                return {
                    "success": False,
                    "error": f"Failed to create safety backup: {pre_backup.get('error', 'Unknown error')}"
                }

            # Extract only the selected files
            for file_name in files_to_restore:
                zipf.extract(file_name, project_root)

        return {
            "success": True,
            "scope": scope,
            "restored_from": str(backup_path),
            "files_restored": len(files_to_restore),
            "pre_restore_backup": pre_backup["backup_path"],
            "message": f"Restored {len(files_to_restore)} {scope} files from {backup_path.name}"
        }

    except zipfile.BadZipFile as e:
        return {
            "success": False,
            "error": f"Corrupted backup file: {str(e)}"
        }
    except Exception as e:
        return {
            "success": False,
            "error": f"Selective restore failed: {str(e)}"
        }


def get_latest_backup(project_root: Path, structure_type: str = None) -> Optional[Path]:
    """
    Get the path to the most recent full backup.

    Args:
        project_root: Path to the Crucible project root
        structure_type: Project structure type

    Returns:
        Path to the latest backup or None if no backups exist
    """
    backups = list_backups(project_root, structure_type)
    full_backups = [b for b in backups if b["type"] == "full"]

    if not full_backups:
        return None

    # Backups are already sorted newest first
    return Path(full_backups[0]["path"])


def main():
    """Main entry point with CLI argument parsing."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Restore Crucible project from backups",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  List all backups:
    python restore_backup.py --list

  Preview a restore (dry run):
    python restore_backup.py --restore backup.zip --dry-run

  Restore from the most recent backup:
    python restore_backup.py --restore latest

  Restore from a full backup:
    python restore_backup.py --restore crucible-backup-20241218-120000.zip

  Restore specific file from incremental:
    python restore_backup.py --incremental "chapter-01*"

  Selectively restore only chapters:
    python restore_backup.py --restore latest --scope chapter

  Selectively restore only state files:
    python restore_backup.py --restore backup.zip --scope state
        """
    )

    parser.add_argument("--list", action="store_true",
                        help="List available backups")
    parser.add_argument("--restore", type=str, metavar="BACKUP",
                        help="Restore from a full backup zip file (use 'latest' for most recent)")
    parser.add_argument("--incremental", type=str, metavar="PATTERN",
                        help="Restore from incremental backups matching pattern")
    parser.add_argument("--scope", type=str, metavar="SCOPE",
                        choices=["chapter", "planning", "story-bible", "state", "outline", "all"],
                        help="Selectively restore only specific types of files")
    parser.add_argument("--dry-run", action="store_true",
                        help="Preview changes without executing")
    parser.add_argument("--project", type=str, metavar="PATH",
                        help="Project root path (default: auto-detect)")
    parser.add_argument("--json", action="store_true",
                        help="Output results as JSON")

    args = parser.parse_args()

    # Determine project root and structure type
    if args.project:
        project_root = Path(args.project)
        _, structure_type = find_crucible_project_with_type(project_root)
    else:
        project_root, structure_type = find_crucible_project_with_type()
        if project_root is None:
            result = {
                "success": False,
                "error": "No Crucible project found. Run from project directory or specify --project"
            }
            print(format_output(result) if args.json else result["error"], file=sys.stderr)
            sys.exit(1)

    project_root = Path(project_root)
    if structure_type is None:
        structure_type = "dotcrucible"

    # Execute requested action
    if args.list:
        backups = list_backups(project_root, structure_type)
        if args.json:
            print(format_output({"success": True, "backups": backups, "structure_type": structure_type}))
        else:
            if not backups:
                print("No backups found.")
            else:
                print(f"Found {len(backups)} backup(s):\n")
                for b in backups:
                    print(f"  [{b['type'].upper():11}] {b['filename']}")
                    print(f"                Size: {b['size_human']}, Modified: {b['modified']}")
                    print()
        sys.exit(0)

    elif args.restore:
        # Handle "latest" keyword
        if args.restore.lower() == "latest":
            backup_path = get_latest_backup(project_root, structure_type)
            if backup_path is None:
                result = {
                    "success": False,
                    "error": "No backups found to restore from"
                }
                print(format_output(result) if args.json else result["error"], file=sys.stderr)
                sys.exit(1)
        else:
            backup_path = Path(args.restore)
            # If not absolute, check in backups directory
            if not backup_path.is_absolute():
                backup_dir = get_backup_directory(project_root, structure_type)
                if (backup_dir / backup_path).exists():
                    backup_path = backup_dir / backup_path

        # Use selective restore if scope is specified
        if args.scope:
            result = selective_restore(
                backup_path, project_root, args.scope,
                structure_type=structure_type, dry_run=args.dry_run
            )
        else:
            result = restore_backup(backup_path, project_root, dry_run=args.dry_run)

        print(format_output(result) if args.json else result.get("message", result.get("error", "Unknown result")))
        sys.exit(0 if result["success"] else 1)

    elif args.incremental:
        result = restore_incremental(args.incremental, project_root, structure_type=structure_type, dry_run=args.dry_run)
        print(format_output(result) if args.json else result.get("message", result.get("error", "Unknown result")))
        sys.exit(0 if result["success"] else 1)

    else:
        # No action specified - default to list
        parser.print_help()
        sys.exit(0)


if __name__ == "__main__":
    main()
