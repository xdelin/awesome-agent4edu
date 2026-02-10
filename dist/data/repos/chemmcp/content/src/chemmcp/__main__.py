import logging

logger = logging.getLogger(__name__)


# Import all the modules to register them with the MCP server
from .tools import *

from .utils.mcp_app import run_mcp_server

if __name__ == "__main__":
    run_mcp_server()
