"""MCP server for Stella system dynamics models."""

import json
from pathlib import Path
from typing import Any

from mcp.server import Server
from mcp.server.stdio import stdio_server
from mcp.types import Tool, TextContent

from .xmile import GraphicalFunction, StellaModel, parse_stmx
from .validator import validate_model

# Global model state (one model at a time)
_current_model: StellaModel | None = None


def get_model() -> StellaModel:
    """Get the current model, raising an error if none exists."""
    if _current_model is None:
        raise ValueError("No model created. Use create_model first.")
    return _current_model


def build_graphical_function(data: dict | None) -> GraphicalFunction | None:
    """Build a GraphicalFunction from tool input."""
    if not data:
        return None
    ypts = data.get("ypts")
    if not ypts:
        raise ValueError("graphical_function requires non-empty ypts")
    xscale = data.get("xscale")
    xpts = data.get("xpts")
    if xscale is not None and xpts is not None:
        raise ValueError("graphical_function cannot define both xscale and xpts")
    yscale = data.get("yscale")
    return GraphicalFunction(
        ypts=[float(val) for val in ypts],
        xscale=(float(xscale["min"]), float(xscale["max"])) if xscale else None,
        xpts=[float(val) for val in xpts] if xpts is not None else None,
        yscale=(float(yscale["min"]), float(yscale["max"])) if yscale else None,
        gf_type=data.get("type"),
    )


# Create MCP server
server = Server("stella-mcp")


@server.list_tools()
async def list_tools() -> list[Tool]:
    """List available tools."""
    graphical_function_schema = {
        "type": "object",
        "description": "Graphical function (lookup table) definition",
        "properties": {
            "ypts": {
                "type": "array",
                "items": {"type": "number"},
                "description": "Y values for the lookup table",
            },
            "xscale": {
                "type": "object",
                "description": "X scale when x points are evenly spaced",
                "properties": {
                    "min": {"type": "number"},
                    "max": {"type": "number"},
                },
                "required": ["min", "max"],
            },
            "xpts": {
                "type": "array",
                "items": {"type": "number"},
                "description": "Explicit X values (same length as ypts)",
            },
            "yscale": {
                "type": "object",
                "description": "Optional Y scale for display",
                "properties": {
                    "min": {"type": "number"},
                    "max": {"type": "number"},
                },
                "required": ["min", "max"],
            },
            "type": {
                "type": "string",
                "description": "Graphical function type (e.g., continuous or discrete)",
            },
        },
        "required": ["ypts"],
    }
    return [
        Tool(
            name="create_model",
            description="Create a new Stella model with specified time settings",
            inputSchema={
                "type": "object",
                "properties": {
                    "name": {"type": "string", "description": "Model name"},
                    "start": {"type": "number", "description": "Simulation start time", "default": 0},
                    "stop": {"type": "number", "description": "Simulation stop time", "default": 100},
                    "dt": {"type": "number", "description": "Time step", "default": 0.25},
                    "method": {"type": "string", "description": "Integration method (Euler or RK4)", "default": "Euler"},
                    "time_units": {"type": "string", "description": "Time units", "default": "Years"},
                },
                "required": ["name"],
            },
        ),
        Tool(
            name="add_stock",
            description="Add a stock (reservoir) to the current model",
            inputSchema={
                "type": "object",
                "properties": {
                    "name": {"type": "string", "description": "Stock name"},
                    "initial_value": {"type": "string", "description": "Initial value (number or equation)"},
                    "units": {"type": "string", "description": "Units", "default": ""},
                    "non_negative": {"type": "boolean", "description": "Prevent negative values", "default": True},
                    "x": {"type": "number", "description": "X position (optional, auto-positioned if not specified)"},
                    "y": {"type": "number", "description": "Y position (optional, auto-positioned if not specified)"},
                },
                "required": ["name", "initial_value"],
            },
        ),
        Tool(
            name="add_flow",
            description="Add a flow between stocks in the current model",
            inputSchema={
                "type": "object",
                "properties": {
                    "name": {"type": "string", "description": "Flow name"},
                    "equation": {"type": "string", "description": "Flow rate equation"},
                    "units": {"type": "string", "description": "Units", "default": ""},
                    "from_stock": {"type": "string", "description": "Source stock (null for external source)"},
                    "to_stock": {"type": "string", "description": "Destination stock (null for external sink)"},
                    "non_negative": {"type": "boolean", "description": "Prevent negative values", "default": True},
                    "x": {"type": "number", "description": "X position (optional, auto-positioned if not specified)"},
                    "y": {"type": "number", "description": "Y position (optional, auto-positioned if not specified)"},
                    "graphical_function": graphical_function_schema,
                },
                "required": ["name", "equation"],
            },
        ),
        Tool(
            name="add_aux",
            description="Add an auxiliary variable (parameter or intermediate calculation) to the current model",
            inputSchema={
                "type": "object",
                "properties": {
                    "name": {"type": "string", "description": "Variable name"},
                    "equation": {"type": "string", "description": "Equation or constant value"},
                    "units": {"type": "string", "description": "Units", "default": ""},
                    "x": {"type": "number", "description": "X position (optional, auto-positioned if not specified)"},
                    "y": {"type": "number", "description": "Y position (optional, auto-positioned if not specified)"},
                    "graphical_function": graphical_function_schema,
                },
                "required": ["name", "equation"],
            },
        ),
        Tool(
            name="add_connector",
            description="Add a connector (dependency arrow) between variables",
            inputSchema={
                "type": "object",
                "properties": {
                    "from_var": {"type": "string", "description": "Source variable name"},
                    "to_var": {"type": "string", "description": "Target variable name (the one using from_var)"},
                },
                "required": ["from_var", "to_var"],
            },
        ),
        Tool(
            name="save_model",
            description="Save the current model to a .stmx file",
            inputSchema={
                "type": "object",
                "properties": {
                    "filepath": {"type": "string", "description": "Output file path (.stmx)"},
                },
                "required": ["filepath"],
            },
        ),
        Tool(
            name="read_model",
            description="Read an existing .stmx file and load it as the current model",
            inputSchema={
                "type": "object",
                "properties": {
                    "filepath": {"type": "string", "description": "Path to .stmx file"},
                },
                "required": ["filepath"],
            },
        ),
        Tool(
            name="validate_model",
            description="Validate the current model for errors and warnings",
            inputSchema={
                "type": "object",
                "properties": {},
            },
        ),
        Tool(
            name="list_variables",
            description="List all variables (stocks, flows, auxiliaries) in the current model",
            inputSchema={
                "type": "object",
                "properties": {},
            },
        ),
        Tool(
            name="get_model_xml",
            description="Get the XMILE XML representation of the current model (for preview)",
            inputSchema={
                "type": "object",
                "properties": {},
            },
        ),
    ]


@server.call_tool()
async def call_tool(name: str, arguments: dict[str, Any]) -> list[TextContent]:
    """Handle tool calls."""
    global _current_model

    try:
        if name == "create_model":
            _current_model = StellaModel(name=arguments["name"])
            _current_model.sim_specs.start = arguments.get("start", 0)
            _current_model.sim_specs.stop = arguments.get("stop", 100)
            _current_model.sim_specs.dt = arguments.get("dt", 0.25)
            _current_model.sim_specs.method = arguments.get("method", "Euler")
            _current_model.sim_specs.time_units = arguments.get("time_units", "Years")
            return [TextContent(
                type="text",
                text=f"Created model '{arguments['name']}' with time range {_current_model.sim_specs.start}-{_current_model.sim_specs.stop}, dt={_current_model.sim_specs.dt}"
            )]

        elif name == "add_stock":
            model = get_model()
            model.add_stock(
                name=arguments["name"],
                initial_value=arguments["initial_value"],
                units=arguments.get("units", ""),
                non_negative=arguments.get("non_negative", True),
                x=arguments.get("x"),
                y=arguments.get("y"),
            )
            pos_info = ""
            if arguments.get("x") is not None and arguments.get("y") is not None:
                pos_info = f" at position ({arguments['x']}, {arguments['y']})"
            return [TextContent(
                type="text",
                text=f"Added stock '{arguments['name']}' with initial value {arguments['initial_value']}{pos_info}"
            )]

        elif name == "add_flow":
            model = get_model()
            model.add_flow(
                name=arguments["name"],
                equation=arguments["equation"],
                units=arguments.get("units", ""),
                from_stock=arguments.get("from_stock"),
                to_stock=arguments.get("to_stock"),
                non_negative=arguments.get("non_negative", True),
                x=arguments.get("x"),
                y=arguments.get("y"),
                graphical_function=build_graphical_function(arguments.get("graphical_function")),
            )
            flow_desc = []
            if arguments.get("from_stock"):
                flow_desc.append(f"from {arguments['from_stock']}")
            if arguments.get("to_stock"):
                flow_desc.append(f"to {arguments['to_stock']}")
            flow_str = " ".join(flow_desc) if flow_desc else "(external)"
            pos_info = ""
            if arguments.get("x") is not None and arguments.get("y") is not None:
                pos_info = f" at position ({arguments['x']}, {arguments['y']})"
            return [TextContent(
                type="text",
                text=f"Added flow '{arguments['name']}' {flow_str}: {arguments['equation']}{pos_info}"
            )]

        elif name == "add_aux":
            model = get_model()
            model.add_aux(
                name=arguments["name"],
                equation=arguments["equation"],
                units=arguments.get("units", ""),
                x=arguments.get("x"),
                y=arguments.get("y"),
                graphical_function=build_graphical_function(arguments.get("graphical_function")),
            )
            pos_info = ""
            if arguments.get("x") is not None and arguments.get("y") is not None:
                pos_info = f" at position ({arguments['x']}, {arguments['y']})"
            return [TextContent(
                type="text",
                text=f"Added auxiliary '{arguments['name']}' = {arguments['equation']}{pos_info}"
            )]

        elif name == "add_connector":
            model = get_model()
            model.add_connector(
                from_var=arguments["from_var"],
                to_var=arguments["to_var"],
            )
            return [TextContent(
                type="text",
                text=f"Added connector from '{arguments['from_var']}' to '{arguments['to_var']}'"
            )]

        elif name == "save_model":
            model = get_model()
            filepath = Path(arguments["filepath"])
            if not filepath.suffix:
                filepath = filepath.with_suffix(".stmx")
            xml_content = model.to_xml()
            filepath.write_text(xml_content, encoding="utf-8")
            return [TextContent(
                type="text",
                text=f"Saved model to {filepath}"
            )]

        elif name == "read_model":
            filepath = Path(arguments["filepath"])
            _current_model = parse_stmx(str(filepath))
            n_stocks = len(_current_model.stocks)
            n_flows = len(_current_model.flows)
            n_aux = len(_current_model.auxs)
            return [TextContent(
                type="text",
                text=f"Loaded model '{_current_model.name}' with {n_stocks} stocks, {n_flows} flows, {n_aux} auxiliaries"
            )]

        elif name == "validate_model":
            model = get_model()
            errors = validate_model(model)
            if not errors:
                return [TextContent(type="text", text="Model validation passed with no errors or warnings.")]

            result_lines = ["Model validation results:"]
            for err in errors:
                prefix = "ERROR" if err.severity == "error" else "WARNING"
                result_lines.append(f"  [{prefix}] {err.category}: {err.message}")
            return [TextContent(type="text", text="\n".join(result_lines))]

        elif name == "list_variables":
            model = get_model()
            lines = [f"Model: {model.name}", ""]

            if model.stocks:
                lines.append("Stocks:")
                for name, stock in model.stocks.items():
                    lines.append(f"  - {stock.name} = {stock.initial_value} [{stock.units}]")
                lines.append("")

            if model.flows:
                lines.append("Flows:")
                for name, flow in model.flows.items():
                    from_str = flow.from_stock or "external"
                    to_str = flow.to_stock or "external"
                    lines.append(f"  - {flow.name}: {from_str} -> {to_str} = {flow.equation}")
                lines.append("")

            if model.auxs:
                lines.append("Auxiliaries:")
                for name, aux in model.auxs.items():
                    lines.append(f"  - {aux.name} = {aux.equation} [{aux.units}]")

            return [TextContent(type="text", text="\n".join(lines))]

        elif name == "get_model_xml":
            model = get_model()
            xml = model.to_xml()
            # Truncate if too long
            if len(xml) > 10000:
                xml = xml[:10000] + "\n... (truncated)"
            return [TextContent(type="text", text=xml)]

        else:
            return [TextContent(type="text", text=f"Unknown tool: {name}")]

    except Exception as e:
        return [TextContent(type="text", text=f"Error: {str(e)}")]


async def run_server():
    """Run the MCP server."""
    async with stdio_server() as (read_stream, write_stream):
        await server.run(read_stream, write_stream, server.create_initialization_options())


def main():
    """Entry point for the MCP server."""
    import asyncio
    asyncio.run(run_server())


if __name__ == "__main__":
    main()
