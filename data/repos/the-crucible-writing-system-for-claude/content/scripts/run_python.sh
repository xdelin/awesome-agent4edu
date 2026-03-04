#!/bin/bash
# Script: run_python.sh
# Purpose: Cross-platform Python 3 executor for Crucible Suite hooks
# Crucible Suite Plugin
#
# This wrapper detects the appropriate Python 3 interpreter across platforms:
# - Unix/macOS: typically 'python3'
# - Windows (Git Bash/MSYS2): typically 'python' or 'py -3'
# - Windows (native): 'py -3' (Python Launcher)
#
# Usage: run_python.sh <script.py> [args...]

script="$1"
shift

# Check if script argument provided
if [ -z "$script" ]; then
    echo "Error: No script specified" >&2
    echo "Usage: run_python.sh <script.py> [args...]" >&2
    exit 1
fi

# Check if script exists
if [ ! -f "$script" ]; then
    echo "Error: Script not found: $script" >&2
    exit 1
fi

# Try Python interpreters in order of preference
# 1. python3 - standard on Unix/macOS, available on some Windows setups
# 2. py -3   - Windows Python Launcher (reliable on Windows)
# 3. python  - fallback, common on Windows, may be Python 3

if command -v python3 &>/dev/null; then
    exec python3 "$script" "$@"
elif command -v py &>/dev/null; then
    exec py -3 "$script" "$@"
elif command -v python &>/dev/null; then
    # Verify it's Python 3
    python_version=$(python --version 2>&1 | grep -oE '[0-9]+' | head -1)
    if [ "$python_version" = "3" ]; then
        exec python "$script" "$@"
    else
        echo "Error: 'python' is Python 2, need Python 3" >&2
        exit 1
    fi
else
    echo "Error: No Python 3 interpreter found" >&2
    echo "Please install Python 3 and ensure it's in your PATH" >&2
    exit 1
fi
