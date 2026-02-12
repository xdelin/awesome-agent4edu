# Stella MCP Server

A [Model Context Protocol (MCP)](https://modelcontextprotocol.io/) server for creating and manipulating [Stella](https://www.iseesystems.com/store/products/stella-professional.aspx) system dynamics models. This enables AI assistants like Claude to programmatically build, read, validate, and save `.stmx` files in the XMILE format.

## What is this for?

**Stella** is a system dynamics modeling tool used for simulating complex systems in fields like ecology, biogeochemistry, economics, and engineering. This MCP server allows AI assistants to:

- **Create models from scratch** - Build stock-and-flow diagrams programmatically
- **Read existing models** - Parse and understand .stmx files
- **Validate models** - Check for errors like undefined variables or missing connections
- **Modify models** - Add stocks, flows, auxiliaries, and connectors
- **Save models** - Export valid XMILE files that open in Stella Professional

This is particularly useful for:
- Teaching system dynamics modeling
- Rapid prototyping of models through natural language
- Batch creation or modification of models
- Documenting and explaining existing models

## Installation

### From PyPI

```bash
pip install stella-mcp
```

### From source

```bash
git clone https://github.com/bradleylab/stella-mcp.git
cd stella-mcp
pip install -e .
```

### Requirements

- Python 3.10+
- `mcp>=1.0.0`

## Configuration

### Claude Desktop

Add to your `claude_desktop_config.json`:

```json
{
  "mcpServers": {
    "stella": {
      "command": "stella-mcp"
    }
  }
}
```

### Claude Code

Add to your `.claude/settings.json`:

```json
{
  "mcpServers": {
    "stella": {
      "command": "stella-mcp"
    }
  }
}
```

### Development mode

If running from source:

```json
{
  "mcpServers": {
    "stella": {
      "command": "python",
      "args": ["-m", "stella_mcp.server"],
      "cwd": "/path/to/stella-mcp"
    }
  }
}
```

## Available Tools

### Model Creation & I/O

| Tool | Description |
|------|-------------|
| `create_model` | Create a new model with name and time settings (start, stop, dt, method) |
| `read_model` | Load an existing .stmx file |
| `save_model` | Save model to a .stmx file |

### Model Building

| Tool | Description |
|------|-------------|
| `add_stock` | Add a stock (reservoir) with initial value and units |
| `add_flow` | Add a flow between stocks with an equation |
| `add_aux` | Add an auxiliary variable (parameter or calculation) |
| `add_connector` | Add a dependency connector between variables |

### Model Inspection

| Tool | Description |
|------|-------------|
| `list_variables` | List all stocks, flows, and auxiliaries |
| `validate_model` | Check for errors (undefined variables, missing connections, etc.) |
| `get_model_xml` | Preview the XMILE XML output |

## Example Usage

### Creating a simple population model

```
User: Create a simple exponential growth model with a population starting at 100
      and a growth rate of 0.1 per year

Claude: [Uses create_model, add_stock, add_aux, add_flow, add_connector, save_model]
        Creates population_growth.stmx with:
        - Stock: Population (initial=100)
        - Aux: growth_rate (0.1)
        - Flow: growth (Population * growth_rate) into Population
```

### Reading and analyzing an existing model

```
User: Read the carbon cycle model and explain what it does

Claude: [Uses read_model, list_variables]
        This model has 3 stocks (Atmosphere, Land Biota, Soil) and 6 flows
        representing carbon exchange through photosynthesis, respiration...
```

### Building a biogeochemical model

```
User: Create a two-box ocean model with surface and deep nutrients

Claude: [Uses create_model, add_stock (x4), add_aux (x8), add_flow (x6), save_model]
        Creates a model with nutrient cycling between surface and deep ocean
        including upwelling, downwelling, biological uptake, and remineralization
```

## Validation

The `validate_model` tool checks for:

- **Undefined variables** - References to variables that don't exist
- **Mass balance issues** - Stocks without flows, flows referencing non-existent stocks
- **Missing connections** - Equations using variables without connectors (warning)
- **Orphan flows** - Flows not connected to any stock
- **Circular dependencies** - Infinite loops in auxiliary calculations

## XMILE Compatibility

- Output files use the [XMILE standard](https://docs.oasis-open.org/xmile/xmile/v1.0/xmile-v1.0.html)
- Compatible with **Stella Professional 1.9+** and **Stella Architect**
- Auto-layout positions elements reasonably; adjust manually in Stella if needed
- Variable names with spaces are converted to underscores internally

## Project Structure

```
stella-mcp/
├── README.md
├── LICENSE
├── pyproject.toml
└── stella_mcp/
    ├── __init__.py
    ├── server.py      # MCP server implementation
    ├── xmile.py       # XMILE XML generation and parsing
    └── validator.py   # Model validation logic
```

## Contributing

Contributions are welcome! Please feel free to submit issues or pull requests.

## License

MIT License - see [LICENSE](LICENSE) for details.

## Acknowledgments

- [Model Context Protocol](https://modelcontextprotocol.io/) by Anthropic
- [ISEE Systems](https://www.iseesystems.com/) for Stella and the XMILE format
