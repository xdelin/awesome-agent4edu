"""
main.py — Entry point for the TutorMCP project.

You can start either:
1. The MCP server (Gradio + FastAPI)
2. The command-line client (interactive)

Usage:
    python main.py server
    python main.py client
"""

import sys
import subprocess

def run_server():
    """Runs the TutorMCP server on the configured port."""
    print("Starting TutorMCP Server...")
    subprocess.run(["python", "-m", "server.tutor_mcp_server"])

def run_client():
    """Runs the TutorMCP client in terminal mode."""
    print("Launching TutorMCP Client...")
    subprocess.run(["python", "-m", "client.runner"])

def main():
    """Parses command-line args and launches the chosen mode."""
    if len(sys.argv) < 2:
        print("Usage: python main.py [server|client]")
        return

    mode = sys.argv[1].lower()

    if mode == "server":
        run_server()
    elif mode == "client":
        run_client()
    else:
        print("Invalid option. Use 'server' or 'client'.")

if __name__ == "__main__":
    main()
