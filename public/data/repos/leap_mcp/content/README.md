# LEAP MCP

### AI Educational Video Generator

Transform any topic into a short explainer video with AI narration and 3Blue1Brown-style animations.

![Python](https://img.shields.io/badge/Python-3.8--3.11-blue)
![License](https://img.shields.io/badge/License-MIT-green)
![MCP](https://img.shields.io/badge/MCP-Compatible-purple)

![LEAP MCP Demo](assets/demo.gif)

## Quick Start

**Prerequisites:**

- Python 3.8-3.11
- [OpenAI API key](https://platform.openai.com/api-keys) - Create account and get API key
- [Claude Code](https://www.anthropic.com/claude-code) - Download and install

**Note:** Currently supports Claude Code only. Support for Cursor, Windsurf, and other MCP clients planned for future releases.

```bash
git clone https://github.com/sid-thephysicskid/leap-mcp
cd leap-mcp
python setup.py
```

The setup will ask for your OpenAI API key. If you skip it, manually add it to `.env`:

```bash
cp .env.example .env
# Edit .env and add: OPENAI_API_KEY=your_actual_key_here
```

**Next steps:**

1. Start Claude Code: `claude`
2. In Claude Code chat, type `/mcp` - you should see `leap-mcp` connected
3. If not connected, restart Claude Code and try `/mcp` again
4. Start creating videos:

```
Create an educational video about black holes
```

## Features

**Videos** - 2-minute structure: Hook → Phenomenon → Mechanism → Synthesis  
**AI Narration** - 6 natural voices powered by OpenAI  
**Animations** - Manim engine (same as 3Blue1Brown)  
**Any Topic** - From quantum physics to cooking recipes

## Examples

```
Create a video about the French Revolution for high schoolers
Make a video about Vincent van Gogh's painting techniques
Create a video explaining cryptocurrency to my grandmother
Make a video explaining MCP servers
```

**Voices:** `nova` (default) • `alloy` • `echo` • `fable` • `shimmer` • `onyx`

## Troubleshooting

**Python 3.13?** Manim's packages don't support Python 3.13+ supported yet - use 3.8-3.11  
**Connection issues?** Restart Claude Code, run `/mcp`  
**API errors?** Check OpenAI key in `.env` file

## Contributing

**Good First Issues:**

- [ ] Add Cursor/Windsurf/Qwen Agent support
- [ ] Create subject or domain specific templates (history, art, business, financial reports etc.)
- [ ] Improve error messages and debugging
- [ ] Add video quality/resolution options

## License

MIT - Built with [Manim](https://manim.community) • [FastMCP](https://github.com/jlowin/fastmcp) • [OpenAI](https://openai.com)
