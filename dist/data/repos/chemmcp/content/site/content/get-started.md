---
title: "Get Started"
draft: false
weight: 1
---

ChemMCP is a chemistry toolkit that empowers AI assistants with advanced chemistry capabilities.

This guide will help you quickly set up and start using ChemMCP.

## Quick Setup

1. **Install uv**

   We recommend using [uv](https://github.com/astral-sh/uv), a fast Python package and project manager. Install uv with the following commands:

   ```bash
   # On macOS and Linux
   curl -LsSf https://astral.sh/uv/install.sh | sh
   ```
   ```powershell
   # On Windows, use Powershell
   powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"
   ```

2. **Download ChemMCP**

   ```bash
   git clone https://github.com/OSU-NLP-Group/ChemMCP.git
   ```
   
3. **Set Up and Install**

   ```bash
   cd ChemMCP
   uv sync
   ```

   This above commands create a Python virtual environment in the .venv folder with all required packages installed. Then activate the environment:

   ```bash
   # On macOS and Linux
   source .venv/bin/activate
   ```

   ```cmd
   # On Windows (CMD)
   .venv\Scripts\activate.bat
   ```
   
   To install ChemMCP:
   
   ```bash
   uv pip install -e . --no-build-isolation
   ```

   **Optional: Use an Existing Environment**
   
   To install ChemMCP in an existing Python or Conda environment:
   
   ```bash
   # Install dependencies
   pip install -r requirements
   
   # Install ChemMCP
   pip install -e --no-build-isolation
   ```

4. **Install docker**

   If you are to use the [PythonExecutor](/tools/python_executor/) tool, which relies on docker, please follow [its official doc](/tools/python_executor/) to install it. Otherwise, please exclude it from the tools to run, using [QuickConfig](/quick-config), or the MCP server would not successfully start.

## Usage

ChemMCP supports two usage modes:

- **MCP**: Run ChemMCP as a [Model Context Protocol (MCP)](https://claude.ai/download) server to easily integrate with MCP-compatible clients like [Claude Desktop](https://github.com/punkpeye/awesome-mcp-clients) and [more](https://github.com/punkpeye/awesome-mcp-clients).

- **Python Calling**: Use ChemMCP as a standard Python package directly in your Python code.

The following sections will walk you through the basic usage of both modes, after setting up required API keys. 

### Required API Keys

Set the following environment variables to enable tool functionality:

- `CHEMSPACE_API_KEY`: Get from [ChemSpace](https://chem-space.com/)

- `RXN4CHEM_API_KEY`: Get from [IBM RXN4Chem](https://rxn.res.ibm.com)

- `TAVILY_API_KEY`: Get from [Tavily](https://tavily.com/)

- `LLM_MODEL_NAME`: Your preferred LLM model name (e.g., `openai/gpt-4o`) in [LiteLLM](https://docs.litellm.ai/docs/#basic-usage) format

- Other LLM credentials, such as `OPENAI_API_KEY`, or Azure-specific keys if using Azure OpenAI, detailed at [LiteLLM usage](https://docs.litellm.ai/docs/#basic-usage)


{{< alert "circle-info" >}}

To use all tools, you'll need all the above variables. If you only plan to use a subset, not all variables are required, and you can use our [**QuickConfig**](/quick-config) utility to auto-configure your tools.

{{< /alert >}}

### MCP Mode

You can run ChemMCP as an MCP server to integrate with different clients. Here we show two examples, and the configuration should be very similar for other senarios.

#### Claude Desktop Integration

Follow [this guide](https://modelcontextprotocol.io/quickstart/server#testing-your-server-with-claude-for-desktop) to set up the configuration file of Claude Desktop. The JSON config for ChemMCP:

```json
{
    "mcpServers": {
        "ChemMCP": {
            "command": "/ABSTRACT/PATH/TO/uv",  // Use `which uv` to get its path
            "args": ["run", "-m", "chemmcp"],
            "toolCallTimeoutMillis": 300000,  // Set this value because some tools may be slow in response of requests
            "env": {
                "CHEMSPACE_API_KEY": "API_KEY",
                "RXN4CHEM_API_KEY": "API_KEY",
                "TAVILY_API_KEY": "API_KEY",
                "LLM_MODEL_NAME": "openai/gpt-4o",  // Or any other LLM names supported by LiteLLM
                // Add required LLM credentials
                // ...
            }
        }
    }
}
```

#### Integration with LLM APIs (e.g., OpenAI Agents SDK)

LLM providers provide their SDKs to support MCP servers. Take [OpenAI Agents SDK](https://openai.github.io/openai-agents-python/mcp/) as an example, connect to ChemMCP with the following code:

```python
async with MCPServerStdio(
    params={
        "command": "/ABSTRACT/PATH/TO/uv",  # Use `which uv` to get its path
        "args": ["--directory", "/ABSTRACT/PATH/TO/ChemMCP", "run", "-m", "chemmcp"],
        "toolCallTimeoutMillis": 300000,  # Set this value because some tools may be slow in response of requests
        "env": {
        		"CHEMSPACE_API_KEY": "API_KEY",
            "RXN4CHEM_API_KEY": "API_KEY",
            "TAVILY_API_KEY": "API_KEY",
            "LLM_MODEL_NAME": "openai/gpt-4o-2024-08-06",  # Or any other LLM names supported by LiteLLM
          	# Aad required LLM credentials
          	# ...
        }
    }
) as server:
    tools = await server.list_tools()
```

### Python Calling Mode

ChemMCP tools can also be called directly in Python:

```python
import os
from chemmcp.tools import WebSearch  # or any of the tool names

# Of course, you can set only those variables required by the tools to use
envs = {
    "CHEMSPACE_API_KEY": "API_KEY",
    "RXN4CHEM_API_KEY": "API_KEY",
    "TAVILY_API_KEY": "API_KEY",
    "LLM_MODEL_NAME": "openai/gpt-4o-2024-08-06",  # Or any other LLM names supported by LiteLLM
    # Aad required LLM credentials
    # ...
}

for key, value in envs:
  os.environ[key] = value

web_search = WebSearch()
result = web_search.run_code('What is the boiling point of water?')
```

## Available Tools

Currently, ChemMCP tools fall into three categories:

1. **General Tools**: General and broad information retrieval and web searching.
   
2. **Molecule Tools**: Various analyses, predictions, and conversions related to chemical compounds and their properties.
   
3. **Reaction Tools**: Chemical reaction prediction and analysis.

See the full list and documentation on the [Tools](/tools) page.

## Next Steps

- Browse all tools on the [Tools](/tools) page.
- Use [QuickConfig](/quick-config) to auto-configure tools based on your needs.
- Contribute your own tools or improve existing ones (see [Dev Guide](/dev-guide)).

## Contact

Have questions or feedback?
- Open an issue for bug reports or feature requests on our [GitHub repository](https://github.com/OSU-NLP-Group/ChemMCP).

- Email us at `yu.3737 at osu.edu` -- we are eager to know your ideas and suggestions!