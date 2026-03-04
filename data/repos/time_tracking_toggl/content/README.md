# @louis030195/toggl-mcp

Dead simple MCP (Model Context Protocol) server for Toggl time tracking. Control your Toggl timer directly from Claude, ChatGPT, or any LLM that supports MCP.

[![Support Development](https://img.shields.io/badge/Support-Development-yellow?style=for-the-badge)](https://buy.stripe.com/14A14n0eZ1wyfix9a8gA802)

## Features

- ⏱️ Start/stop timers
- 📊 View current timer
- 📈 Get today's time entries
- 🗂️ List projects
- 🗑️ Delete time entries

## Installation

### Claude Desktop

Add to your Claude Desktop configuration (`~/Library/Application Support/Claude/claude_desktop_config.json`):

```json
{
  "mcpServers": {
    "toggl": {
      "command": "npx",
      "args": ["-y", "@louis030195/toggl-mcp"],
      "env": {
        "TOGGL_API_KEY": "your-toggl-api-key"
      }
    }
  }
}
```

### Claude Code

```bash
# Install globally in user scope with API key
claude mcp add -s user toggl npx -e TOGGL_API_KEY=your-toggl-api-key -- -y @louis030195/toggl-mcp
```

## Get Your Toggl API Key

1. Go to [Toggl Track Profile](https://track.toggl.com/profile)
2. Scroll down to "API Token"
3. Click "Click to reveal" and copy your token

## Usage

Once configured, you can use natural language to control Toggl:

- "Start tracking work on the MCP server project"
- "Stop the current timer"
- "What am I currently tracking?"
- "Show me today's time entries"
- "List all my projects"

## Tools

### `toggl_start`
Start a new timer with a description and optional project.

### `toggl_stop`
Stop the currently running timer.

### `toggl_current`
Get information about the currently running timer.

### `toggl_today`
Get all time entries for today with total duration.

### `toggl_projects`
List all projects in your workspace.

### `toggl_delete`
Delete a time entry by its ID.

### `toggl_weekly`
Get weekly time tracking summary with total hours and breakdowns.
- Parameters: `week_offset` (optional, number)
  - `0` = current week
  - `-1` = last week
  - `-2` = two weeks ago, etc.

Returns:
- Total hours for the week
- Daily breakdown (hours per day)
- Project breakdown (hours per project)
- Full list of entries

### `toggl_last_week`
Convenience function to get last week's time tracking summary (equivalent to `toggl_weekly` with `week_offset: -1`).

## Development

```bash
# Clone the repo
git clone https://github.com/louis030195/toggl-mcp.git
cd toggl-mcp

# Install dependencies
npm install

# Build
npm run build

# Run locally
TOGGL_API_KEY=your-api-key npm start
```

## License

MIT

## Author

[Louis Beaumont](https://twitter.com/louis030195)
