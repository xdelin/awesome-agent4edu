import os
import logging
from typing import Optional, List, Union
import base64

from mcp.server.fastmcp import Image

from ...utils.base_tool import BaseTool
from ...utils.mcp_app import ChemMCPManager, run_mcp_server
from .jupyter_backbone import JupyterBackbone


logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class PythonExecutor(BaseTool):
    __version__ = "0.1.0"
    name = "PythonExecutor"
    func_name = 'run_code'
    description = "Execute Python code in a Jupyter notebook. New packages can be installed by running `!pip install <package_name>`."
    implementation_description = "Uses a Jupyter kernel to execute Python code and return the results. Supports text output, images (PNG/SVG), and error messages. The kernel is managed through a JupyterBackbone class that handles the communication with the Jupyter kernel gateway."
    categories = ["General"]
    tags = ["Code Execution"]
    required_envs = []
    text_input_sig = [("code", "str", "N/A", "The Python code to execute.")]
    code_input_sig = [("code", "str", "N/A", "The Python code to execute.")]
    output_sig = [("result", "str", "The result of the Python code execution.")]
    examples = [
        {'text_input': {'code': 'print("Hello, world!")'}, 'code_input': {'code': 'print("Hello, world!")'}, 'output': {'result': 'Hello, world!'}},
    ]
    oss_dependencies = [
        ("Jupyter Notebook", "https://github.com/jupyter/notebook", "BSD 3-Clause")
    ]
    services_and_software = [("docker", "https://www.docker.com/")]

    def __init__(self, kernel_url: Optional[str] = None, init=True, interface='code'):
        kernel_url = kernel_url or os.getenv("JUPYTER_KERNEL_GATEWAY_URL", None)
        self.jupyter = JupyterBackbone(kernel_url=kernel_url)
        super().__init__(init, interface=interface)

    def __del__(self):
        if hasattr(self, 'jupyter') and self.jupyter:
            self.jupyter.close()

    def _run_base(self, code: str) -> List[Union[str, Image]]:
        msgs = self.jupyter.execute_jupyter(code)

        parts = []
        for msg in msgs:
            mtype = msg['msg_type']
            content = msg['content']

            if mtype == 'stream':
                parts.append(content.get('text', ''))

            elif mtype in ('execute_result', 'display_data'):
                data = content.get('data', {})
                # PNG image
                if 'image/png' in data:
                    png_bytes = base64.b64decode(data['image/png'])
                    parts.append(Image(data=png_bytes, format='png'))
                # SVG fallback
                elif 'image/svg+xml' in data:
                    parts.append(data['image/svg+xml'])
                # Plain text fallback
                elif 'text/plain' in data:
                    parts.append(data['text/plain'])

            elif mtype == 'error':
                tb = '\n'.join(content.get('traceback', []))
                parts.append(tb)
        
        return parts


if __name__ == "__main__":
    run_mcp_server()
