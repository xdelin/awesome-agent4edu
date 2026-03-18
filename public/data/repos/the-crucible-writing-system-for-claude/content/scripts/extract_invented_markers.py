#!/usr/bin/env python3
"""
Script: extract_invented_markers.py
Purpose: PostToolUse hook for extracting [INVENTED: ...] and [VERIFY] markers from draft files
Crucible Suite Plugin

This hook runs after Write and Edit operations and extracts any [INVENTED: ...] markers
from draft/chapter files, recording them in the story bible for tracking
invented names, places, and other creative elements that require verification.
It also counts [VERIFY] markers to remind the author of uncertain details.

Only processes Write and Edit tool operations (defense-in-depth check).
Empty markers (e.g., [INVENTED: ]) are detected but skipped.
"""

import json
import os
import re
import subprocess
import sys
from typing import Optional, Tuple

from draft_utils import is_draft_file, setup_logging, get_plugin_root

# Ensure Python 3.8+
if sys.version_info < (3, 8):
    print("Error: Python 3.8+ required", file=sys.stderr)
    sys.exit(1)

# Setup logging (enable with CRUCIBLE_DEBUG=1)
logger = setup_logging('extract_invented_markers')


def extract_chapter_scene(file_path: str) -> Tuple[str, str]:
    """
    Extract chapter and scene numbers from file path.

    Handles patterns like:
    - chapter-5.md -> ('5', '0')
    - chapter_5.md -> ('5', '0')
    - ch_3_scene_2.md -> ('3', '2')
    - chapter5-scene2.md -> ('5', '2')
    - chapter10-scene3.md -> ('10', '3')

    Args:
        file_path: Path to the file

    Returns:
        Tuple of (chapter, scene) as strings. Defaults to ('0', '0') if no pattern matches.
    """
    if not file_path:
        return ('0', '0')

    basename = os.path.basename(file_path).lower()

    # Pattern for chapter with optional scene: chapter-5-scene-2, ch_3_scene_2, chapter10-scene3
    pattern_with_scene = r'(?:chapter|ch)[_-]?(\d+)[_-]?scene[_-]?(\d+)'
    match = re.search(pattern_with_scene, basename)
    if match:
        return (match.group(1), match.group(2))

    # Pattern for chapter only: chapter-5, chapter_5, ch_3, chapter10
    pattern_chapter_only = r'(?:chapter|ch)[_-]?(\d+)'
    match = re.search(pattern_chapter_only, basename)
    if match:
        return (match.group(1), '0')

    return ('0', '0')


def find_project_root(start_path: str) -> Optional[str]:
    """
    Find the project root by walking up directory tree.

    Looks for either:
    - story-bible.json (primary indicator)
    - .crucible/ directory (fallback indicator)

    Args:
        start_path: Starting file path to search from

    Returns:
        Directory path containing story-bible.json or .crucible/, or None if not found
    """
    if not start_path:
        return None

    # Start from directory containing the file
    current = os.path.dirname(os.path.abspath(start_path))

    # Walk up directory tree
    while current:
        # Primary: look for story-bible.json
        story_bible_path = os.path.join(current, 'story-bible.json')
        if os.path.isfile(story_bible_path):
            return current

        # Fallback: look for .crucible/ directory
        crucible_dir = os.path.join(current, '.crucible')
        if os.path.isdir(crucible_dir):
            return current

        parent = os.path.dirname(current)
        if parent == current:
            # Reached filesystem root
            break
        current = parent

    return None


def extract_invented_markers(content: str) -> list:
    """
    Extract [INVENTED: ...] markers from content.

    Empty markers (e.g., [INVENTED: ]) are matched but skipped and not recorded.

    Args:
        content: Text content to search for markers

    Returns:
        List of tuples: (description, category)
        Category is inferred from keywords or defaults to "general"
    """
    if not content:
        return []

    pattern = r'\[INVENTED:\s*([^\]]*)\]'
    matches = re.findall(pattern, content)

    results = []
    for match in matches:
        description = match.strip()

        # Skip empty markers
        if not description:
            continue

        description_lower = description.lower()

        # Infer category from keywords
        if any(kw in description_lower for kw in ['character', 'person', 'name']):
            category = 'character'
        elif any(kw in description_lower for kw in ['location', 'place', 'city', 'town', 'village']):
            category = 'location'
        elif any(kw in description_lower for kw in ['item', 'object', 'artifact', 'weapon']):
            category = 'item'
        else:
            category = 'general'

        results.append((description, category))

    return results


def extract_verify_markers(content: str) -> list:
    """
    Extract [VERIFY] markers from content.

    These markers indicate uncertain details that need verification against
    source documents (story bible, outline, world forge, etc.).

    Args:
        content: Text content to search for markers

    Returns:
        List of tuples: (line_number, context)
        Context is up to 50 chars before and after the marker
    """
    if not content:
        return []

    results = []
    lines = content.split('\n')

    for i, line in enumerate(lines, 1):
        # Find all [VERIFY] markers in this line
        for match in re.finditer(r'\[VERIFY\]', line):
            # Get context around the marker (50 chars before and after)
            start = max(0, match.start() - 50)
            end = min(len(line), match.end() + 50)
            context = line[start:end].strip()
            results.append((i, context))

    return results


def main():
    """
    Main entry point for PostToolUse hook.

    Reads hook input from stdin, extracts [INVENTED: ...] markers from draft files,
    and records them in the story bible via update_story_bible.py.
    """
    try:
        # Read JSON input from stdin (PostToolUse hook format)
        input_data = json.loads(sys.stdin.read())

        # Check tool_name for consistency (defense-in-depth)
        tool_name = input_data.get('tool_name', '')
        if tool_name not in ('Write', 'Edit'):
            sys.exit(0)

        tool_input = input_data.get('tool_input', {})
        file_path = tool_input.get('file_path', '')

        # Check if this is a draft file
        if not is_draft_file(file_path):
            logger.debug(f'Skipping non-draft file: {file_path}')
            sys.exit(0)

        # Read the file content if it exists
        if not os.path.isfile(file_path):
            sys.exit(0)

        with open(file_path, 'r', encoding='utf-8') as f:
            content = f.read()

        # Extract markers from content
        invented_markers = extract_invented_markers(content)
        verify_markers = extract_verify_markers(content)

        logger.debug(f'Found {len(invented_markers)} invented, {len(verify_markers)} verify markers in {file_path}')

        if not invented_markers and not verify_markers:
            # No markers found, exit silently
            sys.exit(0)

        # Find project root
        project_root = find_project_root(file_path)

        if not project_root:
            # No project root found, exit gracefully
            sys.exit(0)

        # Get path to update_story_bible.py (with fallback for testing)
        plugin_root = get_plugin_root()
        if not plugin_root:
            logger.debug('Plugin root not found (no CLAUDE_PLUGIN_ROOT or relative path)')
            sys.exit(0)

        update_script = os.path.join(
            plugin_root,
            'skills',
            'crucible-writer',
            'scripts',
            'update_story_bible.py'
        )

        if not os.path.isfile(update_script):
            # Script not found, exit gracefully
            sys.exit(0)

        # Extract chapter and scene from file path
        chapter, scene = extract_chapter_scene(file_path)

        # Record each invented marker in the story bible
        recorded_count = 0
        for description, category in invented_markers:
            try:
                subprocess.run(
                    [
                        sys.executable,
                        update_script,
                        project_root,
                        '--invented',
                        chapter,
                        scene,
                        description,
                        '--category',
                        category
                    ],
                    capture_output=True,
                    text=True,
                    timeout=30
                )
                recorded_count += 1
            except (subprocess.TimeoutExpired, subprocess.SubprocessError):
                # Continue with other markers if one fails
                continue

        # Output JSON with systemMessage if markers were found
        messages = []
        if recorded_count > 0:
            messages.append(f'Recorded {recorded_count} invented detail(s) in story bible')
        if verify_markers:
            messages.append(f'{len(verify_markers)} [VERIFY] marker(s) need verification')

        if messages:
            output = {
                'systemMessage': f'[Crucible] {"; ".join(messages)}.'
            }
            print(json.dumps(output))

        sys.exit(0)

    except json.JSONDecodeError:
        # Invalid JSON input, exit gracefully
        sys.exit(0)
    except Exception:
        # Any other error, exit gracefully (this is a soft process)
        sys.exit(0)


if __name__ == '__main__':
    main()
