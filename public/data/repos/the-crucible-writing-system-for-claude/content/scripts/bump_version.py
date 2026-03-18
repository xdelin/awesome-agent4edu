#!/usr/bin/env python3
"""
Script: bump_version.py
Purpose: Bump version in VERSION file and plugin.json
Usage:
    python bump_version.py patch   # 1.0.3 -> 1.0.4
    python bump_version.py minor   # 1.0.3 -> 1.1.0
    python bump_version.py major   # 1.0.3 -> 2.0.0
    python bump_version.py         # Show current version
"""

import json
import sys
from pathlib import Path


def get_script_dir() -> Path:
    return Path(__file__).parent


def get_plugin_root() -> Path:
    return get_script_dir().parent


def read_version() -> str:
    version_file = get_plugin_root() / "VERSION"
    if version_file.exists():
        return version_file.read_text().strip()
    return "0.0.0"


def write_version(version: str) -> None:
    plugin_root = get_plugin_root()

    # Update VERSION file
    version_file = plugin_root / "VERSION"
    version_file.write_text(f"{version}\n")

    # Update plugin.json
    plugin_json = plugin_root / ".claude-plugin" / "plugin.json"
    if plugin_json.exists():
        data = json.loads(plugin_json.read_text())
        data["version"] = version
        plugin_json.write_text(json.dumps(data, indent=2) + "\n")


def bump_version(current: str, bump_type: str) -> str:
    parts = current.split(".")
    if len(parts) != 3:
        parts = ["0", "0", "0"]

    major, minor, patch = int(parts[0]), int(parts[1]), int(parts[2])

    if bump_type == "major":
        major += 1
        minor = 0
        patch = 0
    elif bump_type == "minor":
        minor += 1
        patch = 0
    elif bump_type == "patch":
        patch += 1
    else:
        raise ValueError(f"Unknown bump type: {bump_type}")

    return f"{major}.{minor}.{patch}"


def main():
    current = read_version()

    if len(sys.argv) < 2:
        print(f"Current version: {current}")
        print("\nUsage:")
        print("  python bump_version.py patch   # 1.0.3 -> 1.0.4")
        print("  python bump_version.py minor   # 1.0.3 -> 1.1.0")
        print("  python bump_version.py major   # 1.0.3 -> 2.0.0")
        return

    bump_type = sys.argv[1].lower()

    if bump_type not in ("major", "minor", "patch"):
        print(f"Error: Unknown bump type '{bump_type}'")
        print("Use: major, minor, or patch")
        sys.exit(1)

    new_version = bump_version(current, bump_type)
    write_version(new_version)

    print(f"Version bumped: {current} -> {new_version}")
    print("\nUpdated files:")
    print("  - VERSION")
    print("  - .claude-plugin/plugin.json")
    print("\nDon't forget to update CHANGELOG.md!")


if __name__ == "__main__":
    main()
