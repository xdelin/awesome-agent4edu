import httpx
from config.settings import MCP_BASE

class MCPClient:
    """
    Handles communication with the MCP server over HTTP.
    Fetches schema and triggers tools like 'explain_concept', 'summarize_text', etc.
    """

    def __init__(self, base_url = MCP_BASE):
        self.base_url = base_url.rstrip("/")
        self.client = httpx.Client(timeout = 60.0)
    
    def call_tool(self, tool_name, payload):
        """Calls a specific tool endpoint on the MCP server."""
        url = f"{self.base_url}/{tool_name}"
        try:
            response = self.client.post(url, json = payload)
            response.raise_for_status()
            return response.json()
        except Exception as e:
            print(f"[MCPClient] Error calling tool '{tool_name}': {e}")
            return {"error": str(e)}