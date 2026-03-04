#!/usr/bin/env python3
"""
Script: draft_utils.py
Purpose: Centralized utilities for Crucible draft file handling
Crucible Suite Plugin

This module provides shared functionality used by anti-hallucination hooks:
- Draft file detection (is_draft_file)
- Debug logging configuration (setup_logging)
- Plugin root resolution with fallback (get_plugin_root)

Environment Variables:
- CRUCIBLE_DEBUG: Set to any truthy value to enable debug logging to stderr
- CLAUDE_PLUGIN_ROOT: Primary plugin root path (set by Claude Code)
"""

import logging
import os
import re
import sys
from typing import Optional

__all__ = [
    'DRAFT_PATH_PATTERNS',
    'is_draft_file',
    'setup_logging',
    'get_plugin_root',
]

# Ensure Python 3.8+
if sys.version_info < (3, 8):
    print("Error: Python 3.8+ required", file=sys.stderr)
    sys.exit(1)


# Patterns that identify draft/chapter files
DRAFT_PATH_PATTERNS = [
    r'[/\\]draft[/\\]',
    r'[/\\]chapters[/\\]',
    r'[/\\]manuscript[/\\]',
    r'chapter[_-]?\d+',
    r'ch[_-]?\d+',
]


def is_draft_file(file_path: str) -> bool:
    """
    Check if file path matches draft/chapter patterns.

    Args:
        file_path: Path to the file being checked

    Returns:
        True if this appears to be a draft/chapter file
    """
    if not file_path:
        return False

    path_lower = file_path.lower()

    for pattern in DRAFT_PATH_PATTERNS:
        if re.search(pattern, path_lower):
            return True

    return False


def setup_logging(name: str) -> logging.Logger:
    """
    Configure and return a logger for the specified module.

    Debug logging is enabled when the CRUCIBLE_DEBUG environment variable
    is set to any truthy value (e.g., "1", "true", "yes").

    When debug is disabled, uses NullHandler to suppress all output.
    When debug is enabled, logs to stderr at DEBUG level.

    Args:
        name: Logger name (typically __name__ or script name)

    Returns:
        Configured logging.Logger instance
    """
    logger = logging.getLogger(name)

    # Check if debug mode is enabled
    debug_enabled = os.environ.get('CRUCIBLE_DEBUG', '').lower() in ('1', 'true', 'yes', 'on')

    if debug_enabled:
        logger.setLevel(logging.DEBUG)
        # Add stderr handler if not already present
        if not logger.handlers:
            handler = logging.StreamHandler(sys.stderr)
            handler.setLevel(logging.DEBUG)
            formatter = logging.Formatter(
                '[%(name)s] %(levelname)s: %(message)s'
            )
            handler.setFormatter(formatter)
            logger.addHandler(handler)
    else:
        logger.setLevel(logging.WARNING)
        # Use NullHandler to suppress output
        if not logger.handlers:
            logger.addHandler(logging.NullHandler())

    return logger


def get_plugin_root() -> Optional[str]:
    """
    Get the plugin root directory path.

    Resolution order:
    1. CLAUDE_PLUGIN_ROOT environment variable (set by Claude Code)
    2. Relative path calculation from this script's location

    The relative fallback allows the scripts to work when run manually
    for testing, without requiring the environment variable.

    Returns:
        Absolute path to plugin root directory, or None if not found
    """
    # Primary: Check environment variable
    plugin_root = os.environ.get('CLAUDE_PLUGIN_ROOT', '')
    if plugin_root and os.path.isdir(plugin_root):
        return plugin_root

    # Fallback: Calculate relative to this script
    # This file is at: plugin/crucible-suite/scripts/draft_utils.py
    # Plugin root is at: plugin/crucible-suite/
    try:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        # Go up one directory from scripts/ to plugin root
        plugin_root = os.path.dirname(script_dir)

        # Validate: plugin root should contain .claude-plugin directory
        plugin_marker = os.path.join(plugin_root, '.claude-plugin')
        if os.path.isdir(plugin_marker):
            return plugin_root

        # Alternative validation: check for skills/ directory
        skills_dir = os.path.join(plugin_root, 'skills')
        if os.path.isdir(skills_dir):
            return plugin_root

    except (OSError, AttributeError):
        pass

    return None
