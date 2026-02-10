import os

from litellm import completion

from ..utils.errors import ChemMCPToolInitError


def llm_completion(*args, llm_model=None, **kwargs):
    if llm_model is None:
        llm_model = os.getenv("LLM_MODEL_NAME", None)
    if llm_model is None:
        raise ChemMCPToolInitError("The LLM_MODEL_NAME environment variable is no set. To use an LLM in a tool, please set the variable to a model name supported by LiteLLM (https://docs.litellm.ai/docs/#basic-usage), and set the corresponding API credentials in the environment variables.")
    
    if len(args) == 0 and 'model' not in kwargs:
        kwargs['model'] = llm_model
    
    return completion(*args, **kwargs)

