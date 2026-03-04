# PyBullet MCP Server

A Model Context Protocol (MCP) server that provides physics simulation capabilities using PyBullet. This server allows external applications to interact with PyBullet physics simulations through standardized MCP protocol endpoints.

## Features

🎯 **Physics Simulation**: Create and manage 3D physics simulations
🧊 **Object Management**: Load URDF objects, set colors, apply forces
⚡ **Real-time Control**: Step simulation, set velocities, reset poses
📊 **State Monitoring**: Query simulation status and object information
🎲 **Scene Generation**: Create random test scenes with multiple objects

## Components

### Resources

The server provides access to simulation data through custom bullet:// URIs:
- `bullet://simulation/status` - Current simulation state and parameters
- `bullet://simulation/objects` - List of all objects in the simulation
- `bullet://simulation/physics_params` - Physics world configuration
- `bullet://objects/{id}` - Individual object details and current state

### Tools

The server implements comprehensive physics simulation tools:

- **create_simulation**: Initialize or reset the physics world
  - Configure gravity, time step, real-time mode
  - Choose between GUI or headless (DIRECT) mode

- **load_object**: Load URDF objects into the simulation
  - Support for standard PyBullet objects (cube, sphere, r2d2, etc.)
  - Configurable position, orientation, scaling, and name

- **set_object_color**: Change object appearance
  - Set RGBA color values for visual customization

- **apply_force**: Apply forces to objects
  - World frame or link frame force application
  - Configurable force vector and application point

- **set_velocity**: Control object motion
  - Set linear and angular velocities directly

- **step_simulation**: Advance the physics simulation
  - Single or multiple simulation steps
  - Respects real-time constraints

- **get_object_info**: Query detailed object state
  - Position, orientation, velocities, and metadata

- **create_random_scene**: Generate test scenarios
  - Configurable number of objects and area size
  - Random positions, orientations, and colors

- **reset_object_pose**: Reset object positions
  - Precise position and orientation control

## Installation

1. Clone or download this MCP server
2. Install dependencies using uv:
   ```bash
   uv sync --dev --all-extras
   ```

## Usage

### With Claude Desktop

Add to your Claude Desktop configuration:

**macOS**: `~/Library/Application\ Support/Claude/claude_desktop_config.json`
**Windows**: `%APPDATA%/Claude/claude_desktop_config.json`

```json
{
  "mcpServers": {
    "bullet-mcp-server": {
      "command": "uv",
      "args": [
        "--directory",
        "/path/to/bullet-mcp-server",
        "run",
        "bullet-mcp-server"
      ]
    }
  }
}
```

### Development

Run the server directly:
```bash
uv run bullet-mcp-server
```

Or in development mode:
```bash
uv run python -m bullet_mcp_server.server
```

## Example Usage

Once connected to Claude Desktop, you can:

1. **Create a physics world**:
   "Create a new physics simulation with standard gravity"

2. **Load objects**:
   "Load a cube at position [0, 0, 2] and a sphere at [1, 1, 3]"

3. **Apply physics**:
   "Apply a force of [10, 0, 0] to the cube and step the simulation 100 times"

4. **Create scenes**:
   "Create a random scene with 15 objects in a 10x10 area"

5. **Monitor state**:
   "Show me the current positions and velocities of all objects"

## Development

The server is built using:
- **PyBullet**: Physics simulation engine
- **MCP SDK**: Model Context Protocol implementation
- **uv**: Fast Python package manager
- **Python 3.8+**: Modern Python features

### Project Structure
```
bullet-mcp-server/
├── src/bullet_mcp_server/
│   ├── __init__.py
│   └── server.py          # Main MCP server implementation
├── .github/
│   └── copilot-instructions.md
├── .vscode/
│   └── mcp.json          # VS Code MCP configuration
├── pyproject.toml        # Project configuration
└── README.md
```

## Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests if applicable
5. Submit a pull request

## License

This project is open source. See the license file for details.

## Support

For issues and questions:
- Check the PyBullet documentation: https://pybullet.org/
- Review MCP protocol docs: https://modelcontextprotocol.io/
- Open an issue in this repository
      ]
    }
  }
  ```
</details>

<details>
  <summary>Published Servers Configuration</summary>
  ```
  "mcpServers": {
    "bullet-mcp-server": {
      "command": "uvx",
      "args": [
        "bullet-mcp-server"
      ]
    }
  }
  ```
</details>

## Development

### Building and Publishing

To prepare the package for distribution:

1. Sync dependencies and update lockfile:
```bash
uv sync
```

2. Build package distributions:
```bash
uv build
```

This will create source and wheel distributions in the `dist/` directory.

3. Publish to PyPI:
```bash
uv publish
```

Note: You'll need to set PyPI credentials via environment variables or command flags:
- Token: `--token` or `UV_PUBLISH_TOKEN`
- Or username/password: `--username`/`UV_PUBLISH_USERNAME` and `--password`/`UV_PUBLISH_PASSWORD`

### Debugging

Since MCP servers run over stdio, debugging can be challenging. For the best debugging
experience, we strongly recommend using the [MCP Inspector](https://github.com/modelcontextprotocol/inspector).


You can launch the MCP Inspector via [`npm`](https://docs.npmjs.com/downloading-and-installing-node-js-and-npm) with this command:

```bash
npx @modelcontextprotocol/inspector uv --directory /home/borges/dev/bullet-mcp-server run bullet-mcp-server
```


Upon launching, the Inspector will display a URL that you can access in your browser to begin debugging.