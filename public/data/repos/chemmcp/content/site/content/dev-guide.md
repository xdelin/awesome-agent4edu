---
title: "Development Guide"
draft: false
weight: 3
---



We sincerely welcome your contribution of new tools. This guide walks you through the steps to create, test, and submit a new MCP-compatible tool in the ChemMCP package.

## Step 0: Setup

To begin with, fork our [Github repo](https://github.com/OSU-NLP-Group/ChemMCP) into your own account and clone to your local machine. You could install uv and do `uv sync` as introduced [here](/get-started/#quick-setup), which creates a Python environment for you and smooths your development.

## Step 1: Choose Your Tool’s Names

You need to determine a tool name and a tool function name for your tool. Specifically:

- **Tool Name**: The name of the tool.

  - Use **PascalCase** (every word head capitalized), noun-phrase style. 
  - Examples: `WebSearch`, `BbbpPredictor`, `Smiles2Iupac`.

- **Function Name**: The Python function-like name of the tool, used in the MCP server.

  - Use **snake-case**, verb-phrase style.
  - Examples: `search_web`, `predict_bbbp`, `convert_smiles_to_iupac`.

- **File Name** (or module name): The Python module name, **matching your tool name**.
  - Lowercase, snake-case.
  - Examples: `web_search.py`, `bbbp_predictor.py`, `smiles2iupac.py`.

## Step 2: Create Your Tool File and Register

Create a Python file with your file name, under `src/chemmcp/tools/`. Then, register your tool to the module map in `src/chemmcp/tools/__init__.py`:

```python
_tool_module_map = {
    # … existing entries …
    "MyNewTool": "my_new_tool",  # Tool Name: File Name.
}
```

## Step 3: Implement Your Tool

Copy the following template to the file, and implement the mentioned variables and functions for your tool:

```python
import os
import logging
from typing import Optional

from openai import OpenAI

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPToolError  # your custom errors
from ..utils.mcp_app import ChemMCPManager, run_mcp_server

# Use `logger.info("msg")` or `logger.debug("msg")` to print messages if neceesary, instead of `print`.
logger = logging.getLogger(__name__)  


@ChemMCPManager.register_tool  # Use this decorator if you want it discoverable in MCP mode
class MyNewTool(BaseTool):  # Class name must match the tool name
    __version__      = "0.1.0"            # semver: MAJOR.MINOR.PATCH
    name             = "MyNewTool"        # must match class name
    func_name        = "do_something"     # snake-case verb phrase
    description      = "Brief description of what MyNewTool does."
    implementation_description = "Brief description of how the tool works, e.g., with what models/packages."
    oss_dependencies = [                  # leave empty if no direct open source software dependencies
        ("Dependency name", "Depenendency URL", "Dependency license (None if not found)")
    ]
    services_and_software = [             # leave empty if no external hosted services and software required
        ("Service or software name", "URL (None if no specific URL)"),
    ]
    categories       = ["General"]        # choose from Molecule | Reaction | General
    tags             = ["API", "Neural Networks", "Molecular Properties", "LLMs"]  # keywords for this tool
    required_envs    = [                  # if this tool needs users to set API keys or any environment variables; otherwise, leave it an empty list
        ("OPENAI_API_KEY", "API key for OpenAI"),
    ]

    # the function signature for the following _run_base function
    # four elements presented in strings:  input_domain_name,    type,   default value ("N/A" if no default),   the description of this input
    code_input_sig   = [
        ("smiles", "str", "N/A", "SMILES string of the molecule."),
        ("threshold",      "float", "0.5", "Cutoff threshold."),
    ]

    # the function signature for the following _run_text function
    # four elements presented in strings:  input_domain_name,    type,   default value ("N/A" if no default),   the description of this input
    text_input_sig   = [
        ("smiles_and_threshold", "str", "N/A", "The SMILES string of the moledule and the cutoff threshold, separated by a space."),
    ]

    # the output for both the _run_base and _run_code function
    # three elements presented in strings: output_domain_name,   type,   the description of the output
    output_sig       = [
        ("result", "float", "The score of the input molecule."),
    ]

    # the concrete example(s) of the input and output of this tool
    # at least one example
    examples         = [
        {
            "code_input": {  # the input to the _run_base function. must match your defined `code_input_sig`
                "molecule_smiles": "CCO",
                "threshold": 0.7
            },
            "text_input": {  # the input to the _run_text function. must match your defined `text_input_sig`
                "smiles_and_threshold": "CCO  0.7"
            },
            "output": {  # the output of both of the function. must match your defined `output_sig`
                "result": 0.13
            }
        },
    ]

    # You need to define your __init__ function, if you have any customized settings for this tool
    # Note: every argument must have a default value
    def __init__(
        self,
        openai_api_key: Optional[str] = None,  # if your tool needs api_keys or any envs
        init: bool = True,                     # whether to run _init_modules. must have, just put this line
        interface: str = "code"                # whether to use the code interface or text interface. must have, just put this line
    ):
        # Load API key or environment variable if any
        if openai_api_key is None:  # first from the environment variables
            openai_api_key = os.getenv("MY_TOOL_API_KEY")
        if openai_api_key is None:  # if no api_key found, raise an error
            # we recommend using the errors defined in `src/chemmcp/utils/errors.py`
            raise ChemMCPToolError("Please set `MY_TOOL_API_KEY`.")
        self.openai_api_key = openai_api_key
        
        super().__init__(init=init, interface=interface)

    # Optional. Recommeneded to load some checkpoints or initialize some external clients here.
    # This function will be called in super().__init__ if init==True.
    def _init_modules(self):
        """
        Load some checkpoints or initialize some external clients.
        """
        self.openai_client = OpenAI(self.openai_api_key)

    # Your must implement this function.
    # The input arguments are the input or your tool.
    # The output is the output of your tool.
    def _run_base(self, smiles: str, threshold: float = 0.5) -> float:
        """
        Core logic of your tool.
        """
        # … your implementation of the tool function …
        # You could call external APIs, load any models, or implement your own algorithms.
        score: float = self.openai_client.predict(smiles)  # this is just a simplified example

        return score

    # This function is for those applications that only support one string as the input.
    # Typically, you just need to parse the input string and get the inputs to the _run_base function, then call _run_base.
    # This is optional: if your _run_base only has one str input argument, then you don't need to manually implement this -- ChemMCP will do it for you.
    # Otherwise, please implement it.
    def _run_text(self, smiles_and_threshold: str) -> str:
        """
        Parse the only text str input, and then call _run_base. Return the result of _run_base.
        """

        # parse the text query into _run_base's input arguments
        smiles, threshold = smiles_and_threshold.strip().split()
        smiles = str(smiles)
        threshold = float(threshold)

        # call _run_base
        result = self._run_base(smiles, threshold)

        return result
```

At the bottom of your file, add:

```python
if __name__ == "__main__":
    run_mcp_server()
```

## Step 4: Test Your Tool

You can write a script, or use Jupyter Notebook, or any other ways you like, to test your tool. Basically, initialize an instance of your tool class, and check if the tool works as expected.

```python
import os
from chemmcp.tools import MyNewTool
os.environ['OPENAI_API_KEY'] = "YOUR_API_KEY"

my_new_tool = MyNewTool()

# it calls your _run_base function
result = my_new_tool.run_code(smiles="CCO", threshold=0.5)

# it calls your _run_text function
result = my_new_tool.run_text(smiles_and_threshold="CCO 0.5")
```

If your tool supports MCP (i.e., you used `@ChemMCPManager.register_tool`), then to ensure it works well under the MCP mode, use [MCP Inspector](https://github.com/modelcontextprotocol/inspector) to test your tool.
```bash
# Assume your file name is `my_new_tool.py`
npx @modelcontextprotocol/inspector uv run -m src.chemmcp.tools.my_new_tool
```

## Step 5: Submit Your Tool

After testing your tool, you can now use git to submit to your forked repo, and then submit a pull request to [our repo](https://github.com/OSU-NLP-Group/ChemMCP). We will check and merge your awesome work to ChemMCP.

Thank you very much 🥰


## Contact

Have questions or feedback?
- Open an issue for bug reports or feature requests on our [GitHub repository](https://github.com/OSU-NLP-Group/ChemMCP).

- Email us at `yu.3737 at osu.edu` -- we are eager to know your ideas and suggestions!
