import base64
from mcp.server.models import InitializationOptions
from mcp.server.lowlevel import NotificationOptions, Server
import mcp.types as types
import mcp.server.stdio
from .wolfram_client import client
import httpx

server = Server("MCP-wolfram-alpha")


@server.list_prompts()
async def handle_list_prompts() -> list[types.Prompt]:
    """
    List available prompts.
    """
    return [
        types.Prompt(
            name="wa",
            description="Ask Wolfram Alpha a question",
            arguments=[
                types.PromptArgument(
                    name="query",
                    description="query to ask Wolfram Alpha",
                    required=True,
                )
            ],
        )
    ]


@server.get_prompt()
async def handle_get_prompt(
    name: str, arguments: dict[str, str] | None
) -> types.GetPromptResult:
    """
    Generate a prompt by combining arguments with server state.
    """

    # TODO: Implement checks for this
    assert arguments is not None, "Arguments are required"

    return types.GetPromptResult(
        description="Ask Wolfram Alpha a question",
        messages=[
            types.PromptMessage(
                role="user",
                content=types.TextContent(
                    type="text",
                    text=f"Use wolfram alpha to answer the following question: {arguments['query']}",
                ),
            )
        ],
    )


@server.list_tools()
async def handle_list_tools() -> list[types.Tool]:
    """
    List available tools.
    Each tool specifies its arguments using JSON Schema validation.
    """
    return [
        types.Tool(
            name="query-wolfram-alpha",
            description="Use Wolfram Alpha to answer a question. This tool should be used when you need complex math or symbolic intelligence.",
            inputSchema={
                "type": "object",
                "properties": {
                    "query": {"type": "string"}  # Correct property: `query` with type `string`
                },
                "required": ["query"]  # Marking `query` as required
            },
        )
    ]


@server.call_tool()
async def handle_call_tool(
    name: str, arguments: dict | None
) -> list[types.TextContent | types.ImageContent | types.EmbeddedResource]:
    """
    Handle tool execution requests.
    Tools can modify server state and notify clients of changes.
    """
    if not arguments:
        raise ValueError("Missing arguments")

    if name == "query-wolfram-alpha":
        results: list[types.TextContent | types.ImageContent | types.EmbeddedResource] = []
        query = arguments.get("query")
        if not query:
            raise ValueError("Missing 'query' parameter for Wolfram Alpha tool")

        try:
            response = await client.aquery(query)
        except Exception as e:
            raise Exception("Failed to query Wolfram Alpha") from e
        
        try:
            async with httpx.AsyncClient() as http_client:
                for pod in response.pods:
                    for subpod in pod.subpods:

                        if subpod.plaintext:  # Handle text content
                            results.append(types.TextContent(
                                type="text",
                                text=subpod.plaintext
                            ))
                            
                        elif subpod.img:  # Handle image content
                            img_url = subpod.img.get("src")
                            if img_url:
                                img_response = await http_client.get(img_url)
                                if img_response.status_code == 200:
                                    img_base64 = base64.b64encode(img_response.content).decode('utf-8')
                                    results.append(types.ImageContent(
                                        type="image",
                                        data=img_base64,
                                        mimeType="image/png"
                                    ))
        except Exception as e:
            raise Exception("Failed to parse response from Wolfram Alpha") from e

        return results

    raise ValueError(f"Unknown tool: {name}")


async def main():
    # Run the server using stdin/stdout streams
    async with mcp.server.stdio.stdio_server() as (read_stream, write_stream):
        await server.run(
            read_stream,
            write_stream,
            InitializationOptions(
                server_name="MCP-wolfram-alpha",
                server_version="0.2.2",
                capabilities=server.get_capabilities(
                    notification_options=NotificationOptions(),
                    experimental_capabilities={},
                ),
            ),
        )
