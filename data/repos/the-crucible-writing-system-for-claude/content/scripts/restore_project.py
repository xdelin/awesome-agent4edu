#!/usr/bin/env python3
"""
Script: restore_project.py
Purpose: Restore a Crucible project from backup
Crucible Suite Plugin

Supports:
- Listing available backups with metadata
- Creating pre-restore safety backup
- Full restoration from ZIP backup
- Selective restoration (specific chapters, planning, etc.)
- All project structure types (dotcrucible, rootlevel, legacy)
"""

import sys
import json
import zipfile
import shutil
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
        find_crucible_project_with_type,
        ensure_directory,
        get_backup_directory,
        get_timestamp,
        format_output,
        safe_read_json,
        decode_path_b64,
        is_base64_encoded_path
    )
except ImportError:
    # Fallback for decode_path_b64 and is_base64_encoded_path
    import base64
    import re

    def decode_path_b64(encoded: str):
        try:
            padding = 4 - (len(encoded) % 4)
            if padding != 4:
                encoded = encoded + ("=" * padding)
            return base64.urlsafe_b64decode(encoded.encode("ascii")).decode("utf-8")
        except (ValueError, UnicodeDecodeError):
            return None

    def is_base64_encoded_path(s: str) -> bool:
        if re.search(r'\.(md|json|txt|py)$', s, re.IGNORECASE):
            return False
        return bool(re.match(r'^[A-Za-z0-9_-]+$', s))
    # Fallback implementations
    def find_crucible_project_with_type(start_path=None):
        current = Path.cwd() if start_path is None else Path(start_path)
        for directory in [current] + list(current.parents):
            if (directory / ".crucible").exists():
                return directory, "dotcrucible"
            if (directory / "state.json").exists():
                return directory, "rootlevel"
            if (directory / "story-bible.json").exists():
                return directory, "legacy"
        return None, None

    def ensure_directory(path):
        path.mkdir(parents=True, exist_ok=True)
        return True

    def get_backup_directory(project_root, structure_type):
        if structure_type == "dotcrucible":
            return project_root / ".crucible" / "backups"
        return project_root / ".crucible-backups"

    def get_timestamp():
        return datetime.now().strftime("%Y%m%d-%H%M%S")

    def format_output(data, for_hook=False):
        return json.dumps(data, indent=2)

    def safe_read_json(path):
        try:
            with open(path, "r", encoding="utf-8") as f:
                return json.load(f)
        except (json.JSONDecodeError, FileNotFoundError, OSError):
            return None


def list_backups(project_root: Path = None, structure_type: str = None) -> Dict[str, Any]:
    """
    List all available backups for the project.

    Returns:
        Dict with:
        - success: bool
        - backups: List of backup info dicts
        - backup_dir: Path to backup directory
    """
    if project_root is None:
        project_root, structure_type = find_crucible_project_with_type()

    if project_root is None:
        return {
            "success": False,
            "error": "No Crucible project found"
        }

    if structure_type is None:
        structure_type = "dotcrucible"

    backup_dir = get_backup_directory(project_root, structure_type)

    if not backup_dir.exists():
        return {
            "success": True,
            "backups": [],
            "backup_dir": str(backup_dir),
            "message": "No backups found"
        }

    backups = []

    # List full backups (ZIP files)
    for zip_file in sorted(backup_dir.glob("crucible-backup-*.zip"), reverse=True):
        try:
            stat = zip_file.stat()
            # Extract timestamp from filename
            timestamp = zip_file.stem.replace("crucible-backup-", "")

            backup_info = {
                "type": "full",
                "file": zip_file.name,
                "path": str(zip_file),
                "timestamp": timestamp,
                "size_bytes": stat.st_size,
                "size_human": format_size(stat.st_size),
                "modified": datetime.fromtimestamp(stat.st_mtime).isoformat()
            }

            # Try to get file count from ZIP
            try:
                with zipfile.ZipFile(zip_file, 'r') as zf:
                    backup_info["files_count"] = len(zf.namelist())
            except zipfile.BadZipFile:
                backup_info["files_count"] = "unknown"

            backups.append(backup_info)
        except OSError:
            continue

    # Check for incremental backups
    incremental_dir = backup_dir / "incremental"
    if incremental_dir.exists():
        incremental_files = list(incremental_dir.glob("*"))
        if incremental_files:
            # Group by date
            dates = set()
            for f in incremental_files:
                # Extract date from filename (YYYYMMDD-HHMMSS-filename)
                parts = f.name.split("-")
                if len(parts) >= 2 and parts[0].isdigit():
                    dates.add(parts[0])

            backups.append({
                "type": "incremental",
                "path": str(incremental_dir),
                "count": len(incremental_files),
                "dates": sorted(dates, reverse=True)[:5]  # Show last 5 dates
            })

    # Load manifest for additional metadata
    manifest_path = backup_dir / "backup-manifest.json"
    manifest = safe_read_json(manifest_path)

    return {
        "success": True,
        "backups": backups,
        "backup_dir": str(backup_dir),
        "manifest": manifest,
        "project_root": str(project_root),
        "structure_type": structure_type
    }


def format_size(size_bytes: int) -> str:
    """Format bytes to human-readable size."""
    for unit in ['B', 'KB', 'MB', 'GB']:
        if size_bytes < 1024:
            return f"{size_bytes:.1f} {unit}"
        size_bytes /= 1024
    return f"{size_bytes:.1f} TB"


def create_pre_restore_backup(project_root: Path, structure_type: str) -> Dict[str, Any]:
    """
    Create a safety backup before restoration.

    This allows undoing the restore if needed.
    """
    backup_dir = get_backup_directory(project_root, structure_type)
    ensure_directory(backup_dir)

    timestamp = get_timestamp()
    backup_name = f"pre-restore-{timestamp}.zip"
    backup_path = backup_dir / backup_name

    files_backed_up = 0
    total_size = 0

    try:
        with zipfile.ZipFile(backup_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
            # Backup key directories
            dirs_to_backup = ["planning", "outline", "draft", "manuscript", "chapters"]
            for dirname in dirs_to_backup:
                dir_path = project_root / dirname
                if dir_path.exists():
                    for file_path in dir_path.rglob("*"):
                        if file_path.is_file():
                            arcname = file_path.relative_to(project_root)
                            zipf.write(file_path, arcname)
                            files_backed_up += 1
                            total_size += file_path.stat().st_size

            # Backup key files
            files_to_backup = ["CLAUDE.md", "story-bible.json", "style-profile.json"]
            for filename in files_to_backup:
                file_path = project_root / filename
                if file_path.exists():
                    zipf.write(file_path, filename)
                    files_backed_up += 1
                    total_size += file_path.stat().st_size

            # Backup state directory for dotcrucible
            if structure_type == "dotcrucible":
                state_dir = project_root / ".crucible" / "state"
                if state_dir.exists():
                    for file_path in state_dir.rglob("*"):
                        if file_path.is_file():
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
            "error": f"Failed to create pre-restore backup: {e}"
        }


def restore_from_backup(
    backup_path: Path,
    project_root: Path,
    structure_type: str,
    scope: str = "full",
    dry_run: bool = False
) -> Dict[str, Any]:
    """
    Restore project from a backup ZIP file.

    Args:
        backup_path: Path to the backup ZIP file
        project_root: Project root directory
        structure_type: Project structure type
        scope: What to restore - "full", "chapters", "planning", "state", "story-bible"
        dry_run: If True, only report what would be restored

    Returns:
        Dict with restoration results
    """
    if not backup_path.exists():
        return {
            "success": False,
            "error": f"Backup file not found: {backup_path}"
        }

    try:
        with zipfile.ZipFile(backup_path, 'r') as zf:
            all_files = zf.namelist()

            # Filter files based on scope
            files_to_restore = []

            for filename in all_files:
                should_restore = False

                if scope == "full":
                    should_restore = True
                elif scope == "chapters":
                    should_restore = any(d in filename for d in ["draft/", "chapters/", "manuscript/"])
                elif scope == "planning":
                    should_restore = "planning/" in filename
                elif scope == "state":
                    should_restore = "state/" in filename or filename in ["state.json", "project-state.json"]
                elif scope == "story-bible":
                    should_restore = "story-bible" in filename.lower()

                if should_restore:
                    files_to_restore.append(filename)

            if dry_run:
                return {
                    "success": True,
                    "dry_run": True,
                    "files_to_restore": files_to_restore,
                    "count": len(files_to_restore),
                    "scope": scope
                }

            # Create pre-restore backup first
            pre_backup = create_pre_restore_backup(project_root, structure_type)
            if not pre_backup["success"]:
                return {
                    "success": False,
                    "error": f"Failed to create safety backup: {pre_backup.get('error')}"
                }

            # Perform restoration
            restored_count = 0
            errors = []

            for filename in files_to_restore:
                try:
                    target_path = project_root / filename
                    ensure_directory(target_path.parent)

                    # Extract file
                    with zf.open(filename) as src:
                        with open(target_path, 'wb') as dst:
                            dst.write(src.read())

                    restored_count += 1
                except Exception as e:
                    errors.append(f"{filename}: {e}")

            return {
                "success": True,
                "restored_count": restored_count,
                "scope": scope,
                "backup_used": str(backup_path),
                "pre_restore_backup": pre_backup["backup_path"],
                "errors": errors if errors else None,
                "message": f"Restored {restored_count} files from backup"
            }

    except zipfile.BadZipFile:
        return {
            "success": False,
            "error": f"Invalid or corrupted backup file: {backup_path}"
        }
    except Exception as e:
        return {
            "success": False,
            "error": f"Restoration failed: {e}"
        }


def find_original_path(safe_name: str, project_root: Path) -> Optional[Path]:
    """
    Reconstruct the original file path from a backup safe_name.

    Supports two encoding formats:
    1. Base64url encoding (new, v1.0.5+): Unambiguous, reversible encoding
    2. Underscore encoding (legacy): Ambiguous but handled with heuristics

    Args:
        safe_name: The encoded filename (e.g., "ZHJhZnQvY2hhcHRlci0xLm1k" or "draft_chapter_1.md")
        project_root: Project root to search for matching files

    Returns:
        Path to the original file, or None if not found
    """
    # First, try base64url decoding (new format from v1.0.5+)
    # This is unambiguous and preferred
    if is_base64_encoded_path(safe_name):
        decoded = decode_path_b64(safe_name)
        if decoded:
            # Successfully decoded - construct path
            candidate = project_root / decoded
            # Return the path even if file doesn't exist (for restore purposes)
            return candidate

    # Fall back to legacy underscore-based heuristics for backward compatibility
    # Common directory prefixes used in Crucible projects
    common_prefixes = [
        ".crucible_state_",
        ".crucible_story-bible_",
        ".crucible_",
        "draft_",
        "chapters_",
        "manuscript_",
        "planning_",
        "outline_",
    ]

    # Try each prefix
    for prefix in common_prefixes:
        if safe_name.startswith(prefix):
            # Reconstruct the directory path
            dir_part = prefix.rstrip("_").replace("_", "/")
            file_part = safe_name[len(prefix):]
            candidate = project_root / dir_part / file_part

            # Check if parent directory exists (file might not exist yet)
            if candidate.parent.exists():
                return candidate

    # Try to find the file by basename in common directories
    basename = safe_name.split("_")[-1] if "_" in safe_name else safe_name
    search_dirs = ["draft", "chapters", "manuscript", "planning", "outline", ".crucible/state"]

    for search_dir in search_dirs:
        dir_path = project_root / search_dir
        if dir_path.exists():
            # Look for files with matching basename
            for file_path in dir_path.rglob(f"*{basename}"):
                if file_path.is_file():
                    return file_path

    # Check if it's a root-level file
    root_files = ["CLAUDE.md", "story-bible.json", "style-profile.json", "state.json", "project-state.json"]
    for root_file in root_files:
        if safe_name == root_file or safe_name.endswith(f"_{root_file}"):
            return project_root / root_file

    # Last resort: try direct underscore-to-slash conversion
    # This will work for simple paths without underscores in filenames
    naive_path = project_root / safe_name.replace("_", "/")
    if naive_path.parent.exists():
        return naive_path

    return None


def restore_incremental(
    project_root: Path,
    structure_type: str,
    target_file: str = None,
    timestamp: str = None
) -> Dict[str, Any]:
    """
    Restore from incremental backups.

    Args:
        project_root: Project root directory
        structure_type: Project structure type
        target_file: Specific file to restore (optional)
        timestamp: Restore to specific timestamp (optional)

    Returns:
        Dict with restoration results
    """
    backup_dir = get_backup_directory(project_root, structure_type)
    incremental_dir = backup_dir / "incremental"

    if not incremental_dir.exists():
        return {
            "success": False,
            "error": "No incremental backups found"
        }

    # Get all incremental backups
    backups = sorted(incremental_dir.glob("*"), key=lambda x: x.name, reverse=True)

    if not backups:
        return {
            "success": False,
            "error": "No incremental backup files found"
        }

    # Filter by timestamp if specified
    if timestamp:
        backups = [b for b in backups if b.name.startswith(timestamp)]
        if not backups:
            return {
                "success": False,
                "error": f"No backups found for timestamp: {timestamp}"
            }

    # Filter by target file if specified
    if target_file:
        target_normalized = target_file.replace("/", "_").replace("\\", "_")
        backups = [b for b in backups if target_normalized in b.name]
        if not backups:
            return {
                "success": False,
                "error": f"No backups found for file: {target_file}"
            }

    # Get the most recent backup for each original file
    # Filename format: YYYYMMDD-HHMMSS-safe_name
    latest_backups = {}
    for backup in backups:
        # Parse filename: YYYYMMDD-HHMMSS-safe_name
        parts = backup.name.split("-", 2)
        if len(parts) >= 3:
            safe_name = parts[2]
            if safe_name not in latest_backups:
                latest_backups[safe_name] = backup

    if not latest_backups:
        return {
            "success": False,
            "error": "Could not parse incremental backup filenames"
        }

    # Create pre-restore backup
    pre_backup = create_pre_restore_backup(project_root, structure_type)
    if not pre_backup["success"]:
        return {
            "success": False,
            "error": f"Failed to create safety backup: {pre_backup.get('error')}"
        }

    # Restore files
    restored = []
    errors = []
    unresolved = []

    for safe_name, backup_file in latest_backups.items():
        # Find the original path using smart reconstruction
        target = find_original_path(safe_name, project_root)

        if target is None:
            unresolved.append({
                "safe_name": safe_name,
                "backup_file": str(backup_file),
                "suggestion": f"Could not determine original path for: {safe_name}"
            })
            continue

        try:
            ensure_directory(target.parent)
            shutil.copy2(backup_file, target)
            restored.append(str(target))
        except Exception as e:
            errors.append(f"{safe_name}: {e}")

    result = {
        "success": True,
        "restored": restored,
        "restored_count": len(restored),
        "pre_restore_backup": pre_backup["backup_path"],
        "message": f"Restored {len(restored)} files from incremental backups"
    }

    if errors:
        result["errors"] = errors
    if unresolved:
        result["unresolved"] = unresolved
        result["message"] += f" ({len(unresolved)} files could not be mapped to original paths)"

    return result


def main():
    """Main entry point for CLI usage."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Restore Crucible project from backup"
    )
    parser.add_argument(
        "--list", "-l",
        action="store_true",
        help="List available backups"
    )
    parser.add_argument(
        "--restore", "-r",
        metavar="BACKUP",
        help="Restore from specified backup file or timestamp"
    )
    parser.add_argument(
        "--incremental", "-i",
        action="store_true",
        help="Restore from incremental backups"
    )
    parser.add_argument(
        "--scope", "-s",
        choices=["full", "chapters", "planning", "state", "story-bible"],
        default="full",
        help="What to restore (default: full)"
    )
    parser.add_argument(
        "--file", "-f",
        metavar="FILE",
        help="Specific file to restore (for incremental)"
    )
    parser.add_argument(
        "--timestamp", "-t",
        metavar="TIMESTAMP",
        help="Target timestamp (YYYYMMDD or YYYYMMDD-HHMMSS)"
    )
    parser.add_argument(
        "--dry-run", "-n",
        action="store_true",
        help="Show what would be restored without doing it"
    )
    parser.add_argument(
        "--project", "-p",
        metavar="PATH",
        help="Project root path (default: auto-detect)"
    )

    args = parser.parse_args()

    # Find project
    if args.project:
        project_root = Path(args.project)
        _, structure_type = find_crucible_project_with_type(project_root)
    else:
        project_root, structure_type = find_crucible_project_with_type()

    if project_root is None:
        print(format_output({
            "success": False,
            "error": "No Crucible project found"
        }))
        sys.exit(1)

    # Handle actions
    if args.list:
        result = list_backups(project_root, structure_type)

    elif args.restore:
        backup_dir = get_backup_directory(project_root, structure_type)

        # Find the backup file
        if args.restore.endswith(".zip"):
            backup_path = Path(args.restore)
            if not backup_path.is_absolute():
                backup_path = backup_dir / args.restore
        else:
            # Search by timestamp
            matches = list(backup_dir.glob(f"*{args.restore}*.zip"))
            if not matches:
                print(format_output({
                    "success": False,
                    "error": f"No backup found matching: {args.restore}"
                }))
                sys.exit(1)
            backup_path = matches[0]

        result = restore_from_backup(
            backup_path,
            project_root,
            structure_type,
            scope=args.scope,
            dry_run=args.dry_run
        )

    elif args.incremental:
        result = restore_incremental(
            project_root,
            structure_type,
            target_file=args.file,
            timestamp=args.timestamp
        )

    else:
        # Default: list backups
        result = list_backups(project_root, structure_type)

    print(format_output(result))
    sys.exit(0 if result.get("success") else 1)


if __name__ == "__main__":
    main()
