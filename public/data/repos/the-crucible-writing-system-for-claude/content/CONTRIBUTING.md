# Contributing to Crucible Suite

Thank you for your interest in contributing to Crucible Suite! This document provides guidelines and information for contributors.

## Table of Contents

- [Code of Conduct](#code-of-conduct)
- [Getting Started](#getting-started)
- [How to Contribute](#how-to-contribute)
  - [Reporting Bugs](#reporting-bugs)
  - [Suggesting Features](#suggesting-features)
  - [Submitting Changes](#submitting-changes)
- [Development Setup](#development-setup)
- [Pull Request Process](#pull-request-process)
- [Style Guide](#style-guide)
- [Questions](#questions)

## Code of Conduct

This project follows our [Code of Conduct](CODE_OF_CONDUCT.md). By participating, you are expected to uphold this code. Please report unacceptable behavior to the project maintainers.

## Getting Started

Before contributing, please:

1. Read this entire document
2. Check existing [issues](https://github.com/forsonny/The-Crucible-Writing-System-For-Claude/issues) to avoid duplicates
3. Familiarize yourself with the [Crucible Structure](README.md#the-crucible-structure)

## How to Contribute

### Reporting Bugs

Before submitting a bug report:

- Check the [troubleshooting section](README.md#troubleshooting) in the README
- Search existing issues to see if the bug has already been reported
- Collect information about your environment

When reporting a bug, include:

- **Clear title** describing the issue
- **Steps to reproduce** the behavior
- **Expected behavior** vs actual behavior
- **Environment details**:
  - Claude Code version
  - Python version (`python --version`)
  - Operating system
  - Plugin version (check `VERSION` file)
- **Error messages** or logs if available
- **Screenshots** if applicable

### Suggesting Features

Feature requests are welcome! When suggesting a feature:

- **Check existing requests** to avoid duplicates
- **Describe the problem** your feature would solve
- **Propose a solution** if you have one in mind
- **Consider the scope** - does it fit the Crucible Suite's purpose?

Good feature requests include:

- Clear use case explanation
- Proposed implementation approach (optional)
- Consideration of backward compatibility
- Impact on existing workflows

### Submitting Changes

We accept contributions in these areas:

- **Bug fixes** - Patches for existing functionality
- **Documentation** - Improvements to docs, comments, or examples
- **New features** - After discussion in an issue
- **Reference materials** - Writing guides, beat descriptions, etc.
- **Python scripts** - Automation improvements

## Development Setup

### Prerequisites

- Claude Code installed and configured
- Python 3.8 or higher
- Git

### Setup Steps

1. **Fork the repository**
   ```bash
   # Fork via GitHub UI, then clone your fork
   git clone https://github.com/YOUR-USERNAME/The-Crucible-Writing-System-For-Claude.git
   cd The-Crucible-Writing-System-For-Claude
   ```

2. **Create a branch**
   ```bash
   git checkout -b feature/your-feature-name
   # or
   git checkout -b fix/your-bug-fix
   ```

3. **Install as local plugin** (for testing)
   ```bash
   # In Claude Code:
   /plugin marketplace add ./plugin
   /plugin install crucible-suite@crucible-writing-system
   ```

4. **Test your changes**
   ```bash
   # Test Python scripts
   cd plugin/crucible-suite/scripts
   python detect_project.py
   python status_reporter.py
   ```

### Directory Structure

```
plugin/crucible-suite/
├── .claude-plugin/          # Plugin manifest
├── agents/                  # Review agent definitions
├── commands/                # Slash command definitions
├── skills/                  # Skill packages
│   ├── crucible-planner/
│   ├── crucible-outliner/
│   ├── crucible-writer/
│   └── crucible-editor/
├── rules/                   # Project rules
├── hooks/                   # Hook configuration
├── scripts/                 # Python automation
└── templates/               # Project templates
```

## Pull Request Process

### Before Submitting

1. **Update documentation** if you changed functionality
2. **Test your changes** thoroughly
3. **Follow the style guide** (see below)
4. **Keep changes focused** - one feature or fix per PR

### PR Checklist

- [ ] Code follows project style guidelines
- [ ] Self-reviewed the changes
- [ ] Added/updated documentation as needed
- [ ] Tested on both Windows and Unix (if applicable)
- [ ] No new warnings or errors
- [ ] Related issue linked (if applicable)

### PR Title Format

Use a clear, descriptive title:

```
feat: Add series support for multi-book projects
fix: Resolve backup path encoding on Windows
docs: Update installation instructions
refactor: Simplify status reporter logic
```

### Review Process

1. Submit your PR with a clear description
2. Maintainers will review within a few days
3. Address any requested changes
4. Once approved, maintainers will merge

## Style Guide

### Python Code

- **Python 3.8+ compatible** - No features requiring newer versions
- **PEP 8** - Follow standard Python style
- **Docstrings** - Document functions with purpose, args, and returns
- **Error handling** - Use try/except with specific exceptions
- **Cross-platform** - Test on Windows and Unix

Example:
```python
def find_project_root(start_path: Path) -> Optional[Path]:
    """
    Find the Crucible project root from a starting path.

    Args:
        start_path: Directory to start searching from

    Returns:
        Path to project root, or None if not found
    """
    current = start_path.resolve()
    while current != current.parent:
        if (current / ".crucible").exists():
            return current
        current = current.parent
    return None
```

### Markdown Files

- **Line length** - Keep under 120 characters
- **Headers** - Use ATX style (`#`, `##`, etc.)
- **Code blocks** - Always specify language
- **Links** - Use relative paths for internal links
- **Tables** - Align columns with pipes

### Commit Messages

Follow [Conventional Commits](https://www.conventionalcommits.org/):

```
<type>: <description>

[optional body]

[optional footer]
```

Types:
- `feat` - New feature
- `fix` - Bug fix
- `docs` - Documentation only
- `style` - Code style (formatting, etc.)
- `refactor` - Code change that neither fixes a bug nor adds a feature
- `test` - Adding or updating tests
- `chore` - Maintenance tasks

## Questions

- **General questions** - Open a [Discussion](https://github.com/forsonny/The-Crucible-Writing-System-For-Claude/discussions)
- **Bug reports** - Open an [Issue](https://github.com/forsonny/The-Crucible-Writing-System-For-Claude/issues)
- **Security issues** - See [SECURITY.md](SECURITY.md)

---

Thank you for contributing to Crucible Suite!
