# Apple Notes MCP Server

> Powerful tools for automating Apple Notes using Model Context Protocol (MCP)

**mcp-name: io.github.henilcalagiya/mcp-apple-notes**

[![PyPI version](https://badge.fury.io/py/mcp-apple-notes.svg)](https://pypi.org/project/mcp-apple-notes/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)

## Overview

Apple Notes MCP Server provides seamless integration of Apple Notes with any MCP-compatible client. It enables full note automation — including creating, reading, updating, and deleting notes — through a simple and secure AppleScript API layer.

## Features

* **Full CRUD support** for Apple Notes (Create, Read, Update, Delete)
* **Works with Continue.dev, Claude Desktop, Perplexity, and other MCP clients**
* **Native AppleScript integration** for reliable macOS automation
* **Comprehensive tools** for Apple Notes automation
* **Automatic installation** via `uvx` or manual setup
* **FastMCP implementation** with modern decorator-based API

## Requirements

* **macOS** - Required for AppleScript support
* **Python 3.10+** - Required for MCP SDK compatibility
* **Apple Notes application** - Must be installed and accessible
* **MCP-compatible client** (e.g., Continue.dev, Claude Desktop)

## Quick Start

### Step 1: Install uv (if not already installed)

**macOS/Linux:**
```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

**Windows (PowerShell):**
```powershell
powershell -c "irm https://astral.sh/uv/install.ps1 | iex"
```

**Alternative (pip):**
```bash
pip install uv
```

### Step 2: Add MCP Configuration

Add this configuration to your MCP client:

**For Perplexity, Cursor, Continue.dev, Claude Desktop:**
```json
{
  "mcpServers": {
    "apple-notes": {
      "command": "uvx",
      "args": ["mcp-apple-notes@latest"]
    }
  }
}
```

**🎉 That's it!** Your MCP client will automatically install and run the package when needed.

> **Note:** Ensure your MCP client has permission to access Apple Notes and execute AppleScript commands in System Preferences → Security & Privacy → Privacy → Accessibility.

## Available Tools

### **Note Management (6 tools)**
- `create_note` - Create notes with HTML content
- `read_note` - Read notes by ID with verification
- `update_note` - Update notes by ID with HTML content
- `delete_note` - Delete notes by ID with verification
- `move_note` - Move notes between folders
- `list_all_notes` - List all notes across all folders

### **Folder Management (5 tools)**
- `create_folder` - Create folders with path support
- `read_folder` - Read folder details by ID
- `rename_folder` - Rename folders by ID
- `delete_folder` - Delete folders by ID
- `move_folder` - Move folders between locations

### **Search & Structure (3 tools)**
- `search_notes` - Search notes by keywords
- `list_folder_with_structure` - Show folder hierarchy
- `list_notes_with_structure` - Show folders + notes hierarchy



## Content Support

**HTML Formatting**: `<h1-h6>`, `<b><i><u>`, `<p><div><br>`, `<ul><ol><li>`, `<table>`, `<a>`

**Special Features**:
- Unicode and emoji support (🚀, ✅, 📝)
- Nested folder paths (up to 5 levels)
- Automatic character escaping
- Rich content with headers, lists, tables

## Architecture

The server follows the MCP protocol specification and is built with a modular architecture:

- **AppleScript Layer** - Handles direct interaction with Apple Notes
- **Tools Layer** - Wraps AppleScript operations for MCP tools
- **FastMCP Server Layer** - Implements MCP protocol using decorators

## Troubleshooting

### Common Issues

1. **"No module named 'mcp'" Error**
   - Use `uvx` for automatic installation
   - Ensure `uv` is installed and in PATH

2. **AppleScript Permission Denied**
   - Grant permission to your terminal/MCP client in System Preferences → Security & Privacy → Privacy → Accessibility

3. **Notes Not Found**
   - Ensure Apple Notes app is installed and accessible
   - Check that notes exist in the default location

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Author

**Henil C Alagiya**

* **GitHub**: [@henilcalagiya](https://github.com/henilcalagiya)
* **LinkedIn**: [Henil C Alagiya](https://linkedin.com/in/henilcalagiya)

**Support & Contributions:**

* 🐛 **Report Issues**: [GitHub Issues](https://github.com/henilcalagiya/mcp-apple-notes/issues)
* 💬 **Questions**: Reach out on [LinkedIn](https://linkedin.com/in/henilcalagiya)
* 🤝 **Contributions**: Pull requests welcome!
