"""
Health check script for the MCP Simple arXiv HTTP endpoint.

Connects to the production MCP server, initializes a session,
lists available tools, and calls list_categories to verify functionality.
"""
import sys
import asyncio
from mcp.client.streamable_http import streamable_http_client
from mcp.client.session import ClientSession
from mcp.types import Implementation

MCP_URL = "https://mcp.andybrandt.net/arxiv"


async def main() -> None:
    """
    Connect to the MCP server, initialize, list tools, and call list_categories.
    """
    # Create the HTTP transport and session using async context managers
    async with streamable_http_client(MCP_URL) as (read_stream, write_stream, _):
        async with ClientSession(
            read_stream,
            write_stream,
            client_info=Implementation(name="healthcheck", version="0.1")
        ) as session:
            # 1) Initialize the session
            init_result = await session.initialize()

            print("=== initialize() ===")
            print(f"Protocol version: {init_result.protocolVersion}")
            print(f"Server: {init_result.serverInfo.name} v{init_result.serverInfo.version}")

            # 2) List available tools
            tools_result = await session.list_tools()

            print("\n=== tools/list ===")
            for tool in tools_result.tools:
                print(f"- {tool.name}: {tool.description or ''}")

            # 3) Call list_categories to verify a tool works
            categories_result = await session.call_tool("list_categories", {})

            print("\n=== list_categories() ===")
            # Tool results contain a list of content items
            # Only show first 500 chars since category list is long
            for content in categories_result.content:
                if hasattr(content, 'text'):
                    text = content.text
                    if len(text) > 500:
                        print(text[:500] + "\n... (truncated)")
                    else:
                        print(text)


if __name__ == "__main__":
    try:
        asyncio.run(main())
        print("\n✓ Health check passed")
    except Exception as e:
        import traceback
        print(f"ERROR: {e}")
        traceback.print_exc()
        sys.exit(1)

