#!/usr/bin/env python3
"""
Script: backup_project.py
Purpose: Create a full backup of a Crucible project
Crucible Suite Plugin
"""

import sys
import json
import os
import shutil
import zipfile
from pathlib import Path
from datetime import datetime

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
                incremental_dir = crucible_dir / "backups" / "incremental"
                if incremental_dir.exists():
                    dirs.append(incremental_dir)
        elif structure_type == "rootlevel":
            state_file = project_root / "state.json"
            if state_file.exists():
                files.append(state_file)
        return {"dirs": dirs, "files": files}


def cleanup_old_backups(backup_dir: Path, keep_count: int = 10) -> dict:
    """
    Delete old full backups, keeping only the most recent ones.

    Args:
        backup_dir: Path to backup directory
        keep_count: Number of backups to keep (default: 10)

    Returns:
        Dict with deleted count and freed bytes
    """
    if not backup_dir.exists():
        return {"deleted": 0, "freed_bytes": 0}

    # Find all full backup zip files
    backup_files = sorted(
        backup_dir.glob("crucible-backup-*.zip"),
        key=lambda p: p.stat().st_mtime,
        reverse=True  # Newest first
    )

    if len(backup_files) <= keep_count:
        return {"deleted": 0, "freed_bytes": 0}

    # Delete oldest backups beyond keep_count
    to_delete = backup_files[keep_count:]
    deleted = 0
    freed_bytes = 0

    for old_backup in to_delete:
        try:
            size = old_backup.stat().st_size
            old_backup.unlink()
            deleted += 1
            freed_bytes += size
        except OSError:
            continue

    return {"deleted": deleted, "freed_bytes": freed_bytes}


def create_full_backup(project_root: Path = None, backup_dir: Path = None) -> dict:
    """
    Create a complete backup of the Crucible project.

    Supports all project structure types:
    - dotcrucible: Standard .crucible/ directory structure
    - rootlevel: state.json at project root (planner-created)
    - legacy: planning/ directory or story-bible.json at root
    """
    # Find project root and structure type
    if project_root is None:
        project_root, structure_type = find_crucible_project_with_type()
        if project_root is None:
            return {
                "success": False,
                "error": "No Crucible project found. Run from project directory or specify path."
            }
    else:
        project_root = Path(project_root)
        # Detect structure type for specified project
        _, structure_type = find_crucible_project_with_type(project_root)
        if structure_type is None:
            structure_type = "legacy"  # Default fallback

    # Get backup directory based on structure
    if backup_dir is None:
        backup_dir = get_backup_directory(project_root, structure_type)
    else:
        backup_dir = Path(backup_dir)
        ensure_directory(backup_dir)

    # Create timestamped backup
    timestamp = get_timestamp()
    backup_name = f"crucible-backup-{timestamp}"
    backup_path = backup_dir / f"{backup_name}.zip"

    # Get paths to backup based on structure type
    backup_paths = get_backup_paths_for_structure(project_root, structure_type)
    backup_dirs = backup_paths["dirs"]
    backup_files = backup_paths["files"]

    files_backed_up = 0
    total_size = 0

    try:
        with zipfile.ZipFile(backup_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
            # Backup directories
            for target in backup_dirs:
                if target.exists():
                    for file_path in target.rglob("*"):
                        if file_path.is_file():
                            arcname = file_path.relative_to(project_root)
                            zipf.write(file_path, arcname)
                            files_backed_up += 1
                            total_size += file_path.stat().st_size

            # Backup individual files
            for file_path in backup_files:
                if file_path.exists():
                    arcname = file_path.relative_to(project_root)
                    zipf.write(file_path, arcname)
                    files_backed_up += 1
                    total_size += file_path.stat().st_size

        # Record backup in manifest
        manifest_path = backup_dir / "backup-manifest.json"
        manifest = []
        if manifest_path.exists():
            try:
                with open(manifest_path, "r", encoding="utf-8") as f:
                    manifest = json.load(f)
            except json.JSONDecodeError:
                manifest = []

        manifest.append({
            "timestamp": timestamp,
            "file": backup_path.name,
            "files_count": files_backed_up,
            "size_bytes": total_size,
            "structure_type": structure_type,
            "created": datetime.now().isoformat()
        })

        # Keep only last 10 backups in manifest (aligned with disk retention)
        manifest = manifest[-10:]

        with open(manifest_path, "w", encoding="utf-8") as f:
            json.dump(manifest, f, indent=2)

        # Clean up old backups (keep last 10 full backups)
        cleanup_result = cleanup_old_backups(backup_dir, keep_count=10)

        result = {
            "success": True,
            "backup_path": str(backup_path),
            "files_backed_up": files_backed_up,
            "size_bytes": total_size,
            "timestamp": timestamp,
            "structure_type": structure_type,
            "message": f"Backup created: {backup_path.name}"
        }

        if cleanup_result["deleted"] > 0:
            result["cleanup"] = {
                "old_backups_deleted": cleanup_result["deleted"],
                "space_freed_bytes": cleanup_result["freed_bytes"]
            }

        return result

    except Exception as e:
        return {
            "success": False,
            "error": f"Backup failed: {str(e)}"
        }


def main():
    # Read input
    project_root = None
    backup_dir = None

    if len(sys.argv) >= 2:
        project_root = Path(sys.argv[1])
        if len(sys.argv) >= 3:
            backup_dir = Path(sys.argv[2])
    else:
        try:
            input_data = json.load(sys.stdin)
            if "project_root" in input_data:
                project_root = Path(input_data["project_root"])
            if "backup_dir" in input_data:
                backup_dir = Path(input_data["backup_dir"])
        except json.JSONDecodeError:
            pass

    result = create_full_backup(project_root, backup_dir)
    print(format_output(result))
    sys.exit(0 if result["success"] else 1)


if __name__ == "__main__":
    main()
