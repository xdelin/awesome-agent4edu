# kaggle-mcp-server
This project provides an MCP (Model Context Protocol) server for interacting with Kaggle from Claude Desktop (or any MCP-compatible client)

# 🚀 Setup
# 1. Clone the repository and install dependencies
     pip install -r requirements.txt

# 2. Configure Kaggle API credentials

      Go to your Kaggle account -> Settings
      
      Scroll to API section → click Create New API Token.
      This will download a file called kaggle.json.

      Place kaggle.json in the following location on your system:
      
      Linux/Mac:
      
      ~/.kaggle/kaggle.json
      
      Windows:
      
      C:\Users\<your-username>\.kaggle\kaggle.json

# 3. Configure Claude Desktop

     To connect Claude with this MCP server, create or update the file:
     
     Claude/claude_desktop_config.json
     
     with the following contents:
     
     {
       "mcpServers": {
         "Kaggle": {
           "command": "<path-to-your-python-executable>",
           "args": ["<path-to-your-kaggle-server.py>"]
         }
       }
     }
     
     Replace <path-to-your-python-executable> with the path to your Python interpreter (e.g. python or C:/Python/python.exe).
     
     Replace <path-to-your-kaggle-server.py> with the path to this repo’s kaggle-server.py file.

# 📜 Usage

     Once configured:
     
     Start Claude Desktop.
     
     Start the Kaggle MCP server (Run kaggle-server.py as a Python File).
     
     You can now use MCP tools defined in kaggle-server.py, for example:
     
     get_competitions_list: Fetches the list of available Kaggle competitions and returns them in a JSON-friendly format.

# 🛠️ About the Code

     kaggle-server.py registers MCP tools using the @mcp.tool() decorator.
     
     The server communicates with Claude via MCP and exposes Kaggle API methods.
     
     Example tool:
     
     @mcp.tool()
     def get_competitions_list() -> list[dict]:
         """
         Fetch the list of available Kaggle competitions.
         """

# ✅ Requirements

      Python 3.10+
          
      Kaggle API credentials (kaggle.json)
          
      Dependencies listed in requirements.txt

# 🔒 Security Notes

     Never commit your kaggle.json file.
     
     Ensure the .kaggle folder has restricted permissions.
