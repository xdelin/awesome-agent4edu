#!/usr/bin/env python3
"""
Script: validate_before_write.py
Purpose: PreToolUse hook for soft validation of anti-hallucination markers
Crucible Suite Plugin

This hook runs before Write operations and checks draft/chapter files for
potentially invented content that lacks proper [INVENTED: ...] or [VERIFY] markers.
It provides soft warnings (does NOT block writes).
"""

import json
import os
import re
import sys
from typing import Optional

from draft_utils import is_draft_file, setup_logging
from cross_platform import find_crucible_project
from pathlib import Path

# Ensure Python 3.8+
if sys.version_info < (3, 8):
    print("Error: Python 3.8+ required", file=sys.stderr)
    sys.exit(1)

# Setup logging (enable with CRUCIBLE_DEBUG=1)
logger = setup_logging('validate_before_write')


def load_established_names(project_root: Optional[Path]) -> set:
    """
    Load established character and location names from story bible.

    Args:
        project_root: Path to the project root containing story-bible.json

    Returns:
        Set of lowercase names that are established in the story
    """
    names = set()

    if project_root is None:
        return names

    story_bible_path = project_root / 'story-bible.json'
    if not story_bible_path.is_file():
        return names

    try:
        with open(story_bible_path, 'r', encoding='utf-8') as f:
            data = json.load(f)

        # Extract character names from character_states
        character_states = data.get('character_states', {})
        for name in character_states.keys():
            names.add(name.lower())
            # Also add individual words for multi-word names
            for word in name.split():
                if len(word) > 2:
                    names.add(word.lower())

        # Extract location names
        locations = data.get('locations', {})
        for name in locations.keys():
            names.add(name.lower())
            for word in name.split():
                if len(word) > 2:
                    names.add(word.lower())

    except (json.JSONDecodeError, IOError):
        pass

    return names


def find_unmarked_entities(content: str, established_names: set = None) -> list:
    """
    Find potentially invented names/places that lack markers.

    Scans for capitalized proper nouns that might be invented without
    being tagged with [INVENTED: ...] or [VERIFY].

    Args:
        content: The content being written

    Returns:
        List of warning strings (max 5)
    """
    warnings = []
    lines = content.split('\n')

    for i, line in enumerate(lines, 1):
        # Skip lines that already have markers
        if '[INVENTED' in line or '[VERIFY' in line:
            continue

        # Skip dialogue attribution and common phrases
        if line.strip().startswith('"') or line.strip().startswith("'"):
            continue

        # Look for patterns suggesting invented names
        # "named X", "called X", "the X Inn/Tavern/etc."
        invention_patterns = [
            (r'\bnamed\s+([A-Z][a-z]+(?:\s+[A-Z][a-z]+)?)\b', 'named'),
            (r'\bcalled\s+([A-Z][a-z]+(?:\s+[A-Z][a-z]+)?)\b', 'called'),
            (r'\bthe\s+([A-Z][a-z]+(?:\s+[A-Z][a-z]+)?)\s+(?:Inn|Tavern|Tower|Gate|Temple|Palace|Hall|District|Quarter|Street|Road|Bridge|River|Mountain|Forest|Valley)\b', 'location'),
        ]

        for pattern, pattern_type in invention_patterns:
            matches = re.findall(pattern, line)
            for match in matches:
                # Skip very common names/words
                if match.lower() in {'the', 'a', 'an', 'he', 'she', 'it', 'they'}:
                    continue
                # Skip names established in story bible
                if established_names and match.lower() in established_names:
                    continue
                warnings.append(
                    f"Line {i}: '{match}' ({pattern_type}) may need [INVENTED: ...] marker"
                )

    # Return first 5 warnings to avoid overwhelming output
    return warnings[:5]


def main():
    """
    Main entry point for PreToolUse hook.

    Reads JSON from stdin, checks if Write operation targets a draft file,
    scans for unmarked invented content, and outputs soft warnings.
    """
    # Read hook input from stdin
    try:
        input_data = json.load(sys.stdin)
    except json.JSONDecodeError:
        # Invalid input - allow write to proceed silently
        sys.exit(0)

    tool_name = input_data.get('tool_name', '')
    tool_input = input_data.get('tool_input', {})

    # Only check Write and Edit operations
    if tool_name not in ('Write', 'Edit'):
        sys.exit(0)

    file_path = tool_input.get('file_path', '')

    # Get content based on tool type:
    # - Edit tool: uses 'new_string' for the replacement text
    # - Write tool: uses 'content' for the full file content
    if tool_name == 'Edit':
        content = tool_input.get('new_string', '')
    else:
        content = tool_input.get('content', '')

    # Only check draft/chapter files
    if not is_draft_file(file_path):
        logger.debug(f'Skipping non-draft file: {file_path}')
        sys.exit(0)

    # Skip empty content
    if not content or len(content.strip()) < 100:
        sys.exit(0)

    # Load established names from story bible if project exists
    established_names = None
    start_path = Path(file_path).parent if file_path else None
    project_root = find_crucible_project(start_path)
    if project_root:
        established_names = load_established_names(project_root)

    # Check for unmarked entities
    warnings = find_unmarked_entities(content, established_names)

    logger.debug(f'Found {len(warnings)} potential unmarked entities in {file_path}')

    if warnings:
        # Output soft warning via systemMessage (does not block)
        # Per hooks.md: systemMessage is shown to user as context
        warning_text = "Anti-hallucination reminder: " + "; ".join(warnings[:3])
        if len(warnings) > 3:
            warning_text += f" (+{len(warnings) - 3} more)"

        output = {
            "systemMessage": warning_text
        }
        print(json.dumps(output))

    # Always exit 0 - this is a soft warning, never blocks
    sys.exit(0)


if __name__ == '__main__':
    main()
