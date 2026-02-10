import os
import logging
from typing import Optional

from tavily import TavilyClient

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPSearchFailError, ChemMCPApiNotFoundError
from ..utils.mcp_app import ChemMCPManager, run_mcp_server


logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class WebSearch(BaseTool):
    __version__ = "0.1.0"
    name = "WebSearch"
    func_name = 'search_web'
    description = "Search the web for any questions and knowledge and obtain a concise answer based on thesearch results."
    implementation_description = "Uses the [Tavily](https://tavily.com/) API to search the web for any questions and knowledge and obtain a concise answer based on the search results."
    oss_dependencies = []
    services_and_software = [("Tavily", "https://www.tavily.com/")]
    categories = ["General"]
    tags = ["Web Searching", "LLMs", "Neural Networks", "APIs"]
    required_envs = [("TAVILY_API_KEY", "The API key for [Tavily](https://tavily.com/).")]
    text_input_sig = [("query", "str", "N/A", "The search query.")]
    code_input_sig = [("query", "str", "N/A", "The search query.")]
    output_sig = [("result", "str", "The answer to the search query summarized by Tavily's LLM.")]
    examples = [
        {'text_input': {'query': 'What is the boiling point of water?'}, 'code_input': {'query': 'What is the boiling point of water?'}, 'output': {'result': 'The boiling point of water at sea level is 100°C (212°F).'}},
    ]

    def __init__(self, tavily_api_key: Optional[str] = None, init=True, interface='code'):
        if tavily_api_key is None:
            tavily_api_key = os.getenv("TAVILY_API_KEY", None)
        if tavily_api_key is None:
            raise ChemMCPApiNotFoundError("Cannot find the API key for Tavily. Please set the TAVILY_API_KEY environment variable.")
        self.client = TavilyClient(api_key=tavily_api_key)
        super().__init__(init, interface=interface)

    def _run_base(self, query: str) -> str:
        logger.info("Running WebSearch with query: %s", query)
        response = self.client.search(query, search_depth='advanced', include_answer=True)
        try:
            answer = response['answer']
        except KeyError as e:
            raise ChemMCPSearchFailError(f"Error searching the web: {e}") from e
        return answer


if __name__ == "__main__":
    run_mcp_server()
