"""
Handles environment variables and configuration constants for the project
"""
import os
from dotenv import load_dotenv

# load envt var from .env
load_dotenv()

# OpenAI API
OPENAI_API_KEY = os.getenv("OPENAI_API_KEY")
if not OPENAI_API_KEY:
    raise ValueError("OPENAI_API_KEY missing in .env file")

# default model name
MODEL_NAME = "gpt-4o-mini"

# MCP server and base url
MCP_PORT = 7860
MCP_BASE = f"http://localhost:{MCP_PORT}/api"