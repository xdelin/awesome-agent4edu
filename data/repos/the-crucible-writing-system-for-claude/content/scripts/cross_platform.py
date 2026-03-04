#!/usr/bin/env python3
"""
Script: cross_platform.py
Purpose: Shared utilities for cross-platform compatibility
Crucible Suite Plugin
"""

import sys
import os
import json
import shutil
import base64
from pathlib import Path
from datetime import datetime
from typing import Optional, Dict, Any, List, Tuple

# Ensure Python 3.8+
if sys.version_info < (3, 8):
    print("Error: Python 3.8+ required", file=sys.stderr)
    sys.exit(1)


def _check_directory_for_markers(directory: Path) -> Tuple[Optional[Path], Optional[str]]:
    """
    Check a single directory for Crucible project markers.

    This helper function checks for project markers in priority order
    and returns the directory and structure type if found.

    Priority order:
    1. .crucible/ directory -> "dotcrucible" structure
    2. state.json file -> "rootlevel" structure
    3. story-bible.json file -> "legacy" structure
    4. planning/ directory -> "legacy" structure

    Args:
        directory: The directory to check for markers

    Returns:
        Tuple of (directory, structure_type) if markers found,
        (None, None) if no markers found
    """
    # Priority 1: .crucible directory (standard structure)
    crucible_dir = directory / ".crucible"
    if crucible_dir.exists() and crucible_dir.is_dir():
        return directory, "dotcrucible"

    # Priority 2: state.json at root (planner-created structure)
    state_file = directory / "state.json"
    if state_file.exists():
        return directory, "rootlevel"

    # Priority 3: story-bible.json (legacy/transition structure)
    story_bible = directory / "story-bible.json"
    if story_bible.exists():
        return directory, "legacy"

    # Priority 4: planning/ directory (legacy structure)
    planning_dir = directory / "planning"
    if planning_dir.exists() and planning_dir.is_dir():
        return directory, "legacy"

    return None, None


def find_crucible_project_with_type(start_path: Optional[Path] = None) -> Tuple[Optional[Path], Optional[str]]:
    """
    Find a Crucible project by looking for project markers.

    This is the CANONICAL project detection function with structure type info.
    All Crucible scripts should use this function for consistent project detection.

    Search order:
    1. Current directory and parent directories (upward search)
    2. Immediate subdirectories of start_path (one level deep, fallback)

    Detection priority (first match wins within each directory):
    1. .crucible/ directory -> "dotcrucible" structure (standard)
    2. state.json at root -> "rootlevel" structure (planner-created)
    3. story-bible.json at root -> "legacy" structure
    4. planning/ directory -> "legacy" structure

    Args:
        start_path: Directory to start searching from (default: cwd)

    Returns:
        Tuple of (project_root, structure_type) where structure_type is:
        - "dotcrucible": .crucible/ directory structure (standard)
        - "rootlevel": state.json at project root (planner-created projects)
        - "legacy": story-bible.json or planning/ directory (older projects)
        - None if no project found (both elements will be None)
    """
    if start_path is None:
        env_project_dir = os.environ.get("CLAUDE_PROJECT_DIR")
        current = Path(env_project_dir) if env_project_dir else Path.cwd()
    else:
        current = Path(start_path)

    # Check current directory and parents (upward search)
    for directory in [current] + list(current.parents):
        result = _check_directory_for_markers(directory)
        if result[0] is not None:
            return result

    # Fallback: Search immediate subdirectories (one level deep)
    try:
        for subdir in current.iterdir():
            if not subdir.is_dir() or subdir.name.startswith('.'):
                continue
            result = _check_directory_for_markers(subdir)
            if result[0] is not None:
                return result
    except (PermissionError, OSError):
        # Handle case where directory listing fails
        pass

    return None, None


def find_crucible_project(start_path: Optional[Path] = None) -> Optional[Path]:
    """
    Find a Crucible project by looking for project markers.

    This is a convenience wrapper that returns only the path.
    For structure type information, use find_crucible_project_with_type().

    Args:
        start_path: Directory to start searching from (default: cwd)

    Returns:
        Path to project root or None if not found
    """
    project_root, _ = find_crucible_project_with_type(start_path)
    return project_root


# Alias for backward compatibility
def get_crucible_root() -> Optional[Path]:
    """Find the Crucible project root. Alias for find_crucible_project()."""
    return find_crucible_project()


def get_plugin_root() -> Optional[Path]:
    """Get the plugin root directory from environment or detection."""
    # Check environment variable first
    env_root = os.environ.get("CLAUDE_PLUGIN_ROOT")
    if env_root:
        return Path(env_root)

    # Try to detect from script location
    script_path = Path(__file__).resolve()
    # scripts/ is inside plugin root
    if script_path.parent.name == "scripts":
        return script_path.parent.parent

    return None


def ensure_directory(path: Path) -> bool:
    """Ensure a directory exists, creating it if necessary."""
    try:
        path.mkdir(parents=True, exist_ok=True)
        return True
    except OSError as e:
        print(f"Error creating directory {path}: {e}", file=sys.stderr)
        return False


def safe_read_json(path: Path) -> Optional[Dict[str, Any]]:
    """Safely read a JSON file, returning None on error."""
    try:
        with open(path, "r", encoding="utf-8") as f:
            return json.load(f)
    except (json.JSONDecodeError, FileNotFoundError, OSError) as e:
        print(f"Error reading {path}: {e}", file=sys.stderr)
        return None


def safe_write_json(path: Path, data: Dict[str, Any]) -> bool:
    """
    Safely write a JSON file with proper error handling.

    Uses atomic write pattern: write to temp file, then rename.
    Handles Windows-specific issues where rename fails if target exists.
    """
    try:
        # Ensure parent directory exists
        path.parent.mkdir(parents=True, exist_ok=True)

        # Write to temp file first, then rename (atomic on most systems)
        temp_path = path.with_suffix(".tmp")
        with open(temp_path, "w", encoding="utf-8") as f:
            json.dump(data, f, indent=2, ensure_ascii=False)

        # Rename temp to final - handle Windows atomicity issues
        try:
            temp_path.replace(path)
        except OSError:
            # On Windows, replace() fails if target exists and is locked
            # Fall back to delete-then-move pattern
            if sys.platform == "win32":
                try:
                    path.unlink(missing_ok=True)
                except OSError:
                    pass
                shutil.move(str(temp_path), str(path))
            else:
                raise
        return True
    except OSError as e:
        print(f"Error writing {path}: {e}", file=sys.stderr)
        # Clean up temp file if it exists
        try:
            temp_path.unlink(missing_ok=True)
        except OSError:
            pass
        return False


def get_timestamp() -> str:
    """Get a filesystem-safe timestamp string."""
    return datetime.now().strftime("%Y%m%d-%H%M%S")


def get_iso_timestamp() -> str:
    """Get an ISO format timestamp."""
    return datetime.now().isoformat()


def count_words(text: str) -> int:
    """Count words in a text string."""
    return len(text.split())


def extract_title_from_claude_md(project_root: Path) -> Optional[str]:
    """
    Extract the book title from CLAUDE.md file.

    Args:
        project_root: Path to project root containing CLAUDE.md

    Returns:
        Title string or None if not found
    """
    claude_md = project_root / "CLAUDE.md"
    if not claude_md.exists():
        return None

    try:
        with open(claude_md, "r", encoding="utf-8") as f:
            for line in f:
                if "Book Title:" in line:
                    return line.split(":", 1)[1].strip()
    except (OSError, UnicodeDecodeError):
        pass

    return None


def count_words_in_file(path: Path) -> int:
    """Count words in a file."""
    try:
        with open(path, "r", encoding="utf-8") as f:
            return count_words(f.read())
    except (FileNotFoundError, OSError):
        return 0


def find_files(directory: Path, pattern: str) -> List[Path]:
    """Find files matching a glob pattern."""
    if not directory.exists():
        return []
    return sorted(directory.glob(pattern))


def backup_file(source: Path, backup_dir: Path) -> Optional[Path]:
    """Create a timestamped backup of a file."""
    if not source.exists():
        return None

    ensure_directory(backup_dir)

    timestamp = get_timestamp()
    backup_name = f"{source.stem}-{timestamp}{source.suffix}"
    backup_path = backup_dir / backup_name

    try:
        shutil.copy2(source, backup_path)
        return backup_path
    except OSError as e:
        print(f"Error backing up {source}: {e}", file=sys.stderr)
        return None


def encode_path_b64(path: str) -> str:
    """
    Encode a relative path using URL-safe base64.

    This creates a reversible, filesystem-safe encoding for backup filenames.
    Unlike underscore-based flattening, this encoding is unambiguous and
    can always be decoded back to the original path.

    Args:
        path: A relative path string (e.g., "draft/chapter-1.md")

    Returns:
        URL-safe base64 encoded string without padding (e.g., "ZHJhZnQvY2hhcHRlci0xLm1k")

    Example:
        >>> encode_path_b64("draft/chapter-1.md")
        'ZHJhZnQvY2hhcHRlci0xLm1k'
    """
    return base64.urlsafe_b64encode(path.encode("utf-8")).decode("ascii").rstrip("=")


def decode_path_b64(encoded: str) -> Optional[str]:
    """
    Decode a URL-safe base64 encoded path back to the original.

    Args:
        encoded: URL-safe base64 encoded string (with or without padding)

    Returns:
        Original path string, or None if decoding fails

    Example:
        >>> decode_path_b64("ZHJhZnQvY2hhcHRlci0xLm1k")
        'draft/chapter-1.md'
    """
    try:
        # Add back padding if needed (base64 requires length to be multiple of 4)
        padding = 4 - (len(encoded) % 4)
        if padding != 4:
            encoded = encoded + ("=" * padding)
        return base64.urlsafe_b64decode(encoded.encode("ascii")).decode("utf-8")
    except (ValueError, UnicodeDecodeError):
        return None


def is_base64_encoded_path(s: str) -> bool:
    """
    Check if a string looks like a base64-encoded path.

    Used to distinguish between legacy underscore-encoded backup filenames
    and new base64-encoded filenames for backward compatibility.

    Args:
        s: String to check (typically a backup filename without timestamp prefix)

    Returns:
        True if the string matches base64url pattern (alphanumeric + - + _ only,
        no path separators or dots except possibly at the very end)
    """
    # Base64url uses: A-Z, a-z, 0-9, -, _ (and optionally = for padding)
    # Legacy underscore encoding would have dots (.) for file extensions
    # and often readable words
    import re
    # If it contains a dot followed by common extension, it's likely legacy
    if re.search(r'\.(md|json|txt|py)$', s, re.IGNORECASE):
        return False
    # If it matches pure base64url pattern, it's likely encoded
    return bool(re.match(r'^[A-Za-z0-9_-]+$', s))


def get_project_state(project_root: Path) -> Optional[Dict[str, Any]]:
    """Load the current project state from various state files."""
    state = {
        "root": str(project_root),
        "phase": "unknown",
        "planning": None,
        "outline": None,
        "draft": None,
        "edit": None
    }

    state_dir = project_root / ".crucible" / "state"

    # Check for planning state
    planning_state = state_dir / "planning-state.json"
    if planning_state.exists():
        state["planning"] = safe_read_json(planning_state)
        state["phase"] = "planning"

    # Check for outline state
    outline_state = state_dir / "outline-state.json"
    if outline_state.exists():
        state["outline"] = safe_read_json(outline_state)
        state["phase"] = "outlining"

    # Check for draft state
    draft_state = state_dir / "draft-state.json"
    if draft_state.exists():
        state["draft"] = safe_read_json(draft_state)
        state["phase"] = "writing"

    # Check for edit state
    edit_state = state_dir / "edit-state.json"
    if edit_state.exists():
        state["edit"] = safe_read_json(edit_state)
        state["phase"] = "editing"

    return state


def format_output(data: Dict[str, Any], for_hook: bool = False) -> str:
    """Format output data as JSON, optionally for hook consumption."""
    if for_hook:
        return json.dumps(data, indent=2)
    return json.dumps(data, indent=2, ensure_ascii=False)


def get_python_command() -> str:
    """
    Get the appropriate Python command for the current platform.

    On Windows, 'python' is typically the correct command.
    On Unix-like systems, 'python3' is preferred to avoid Python 2.

    Returns:
        "python" on Windows, "python3" on Unix/Linux/macOS
    """
    if sys.platform == "win32":
        return "python"
    else:
        return "python3"


def get_backup_directory(project_root: Path, structure_type: str) -> Path:
    """
    Get the backup directory path based on project structure type.

    Args:
        project_root: Path to the project root
        structure_type: One of "dotcrucible", "rootlevel", or "legacy"

    Returns:
        Path to the backup directory (creates it if needed)

    For dotcrucible: .crucible/backups/
    For rootlevel/legacy: .crucible-backups/ (at project root to avoid conflicts)
    """
    if structure_type == "dotcrucible":
        backup_dir = project_root / ".crucible" / "backups"
    else:
        # For rootlevel and legacy, use a separate hidden directory
        # to avoid mixing with project files
        backup_dir = project_root / ".crucible-backups"

    ensure_directory(backup_dir)
    return backup_dir


def get_backup_paths_for_structure(project_root: Path, structure_type: str) -> Dict[str, List[Path]]:
    """
    Get the directories and files to backup based on project structure type.

    Args:
        project_root: Path to the project root
        structure_type: One of "dotcrucible", "rootlevel", or "legacy"

    Returns:
        Dict with "dirs" and "files" lists of paths to backup
    """
    dirs: List[Path] = []
    files: List[Path] = []

    # Common directories across all structures
    common_dirs = ["planning", "outline", "draft", "manuscript", "chapters"]
    for dirname in common_dirs:
        dir_path = project_root / dirname
        if dir_path.exists():
            dirs.append(dir_path)

    # Common files across all structures
    common_files = ["CLAUDE.md", "story-bible.json", "style-profile.json"]
    for filename in common_files:
        file_path = project_root / filename
        if file_path.exists():
            files.append(file_path)

    # Structure-specific paths
    if structure_type == "dotcrucible":
        crucible_dir = project_root / ".crucible"
        if crucible_dir.exists():
            # State directory
            state_dir = crucible_dir / "state"
            if state_dir.exists():
                dirs.append(state_dir)
            # Story bible markdown files
            story_bible_dir = crucible_dir / "story-bible"
            if story_bible_dir.exists():
                dirs.append(story_bible_dir)
            # Incremental backups (include in full backups)
            incremental_dir = crucible_dir / "backups" / "incremental"
            if incremental_dir.exists():
                dirs.append(incremental_dir)

    elif structure_type == "rootlevel":
        # state.json is at project root
        state_file = project_root / "state.json"
        if state_file.exists():
            files.append(state_file)

    elif structure_type == "legacy":
        # Legacy may have various configurations
        # Already handled by common dirs/files
        pass

    return {"dirs": dirs, "files": files}


# Module-level exports
__all__ = [
    "find_crucible_project",
    "find_crucible_project_with_type",
    "get_crucible_root",
    "get_plugin_root",
    "ensure_directory",
    "safe_read_json",
    "safe_write_json",
    "get_timestamp",
    "get_iso_timestamp",
    "count_words",
    "count_words_in_file",
    "extract_title_from_claude_md",
    "find_files",
    "backup_file",
    "encode_path_b64",
    "decode_path_b64",
    "is_base64_encoded_path",
    "get_project_state",
    "format_output",
    "get_python_command",
    "get_backup_directory",
    "get_backup_paths_for_structure"
]
