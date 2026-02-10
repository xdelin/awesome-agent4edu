from abc import ABC, abstractmethod
import logging
import inspect
from typing import Dict, List, Tuple, Literal, Any, Optional
from pydantic import BaseModel, Field, field_validator, ConfigDict, ValidationError, ValidationInfo

from .errors import ChemMCPToolMetadataError


logger = logging.getLogger(__name__)


class ToolMeta(BaseModel):
    version: str = Field(..., description="The version of the tool.", pattern=r"^\d+\.\d+\.\d+$", alias='__version__')
    name: str = Field(..., description="The name of the tool.", min_length=1)
    func_name: str = Field(..., description="The function name of the tool.", min_length=1)
    description: str = Field(..., description="The description of the tool.", min_length=1)
    implementation_description: str = Field(..., description="The implementation description of the tool.", min_length=1)
    oss_dependencies: List[Tuple[str, str, Optional[str]]] = Field(..., description="The implementation sources of the tool. For example, if the code is borrowed from another project, then it should be listed here. Each element is a tuple of (source_name, source_url, source license (None if unknown)).")
    services_and_software: List[Tuple[str, Optional[str]]] = Field(..., description="The services and software that the tool depends on. Each element is a tuple of (service_name, service_url).")
    categories: List[Literal["Molecule", "Reaction", "General"]] = Field(..., description="The categories of the tool.", min_length=1)
    tags: List[str] = Field(..., description="The tags of the tool.", min_length=1)
    required_envs: List[Tuple[str, str]] = Field(..., description="The required environment variables for the tool.")
    text_input_sig: List[Tuple[str, str, str, str]] = Field(..., description="The text input signature of the tool. Each element is a tuple of (arg_name, arg_type, arg_default, arg_description).")
    code_input_sig: List[Tuple[str, str, str, str]] = Field(..., description="The code input signature of the tool. Each element is a tuple of (arg_name, arg_type, arg_default, arg_description).")
    output_sig: List[Tuple[str, str, str]] = Field(..., description="The output signature of the tool. Each element is a tuple of (output_name, output_type, output_description).")
    examples: List[Dict[Literal["text_input", "code_input", "output"], Dict[str, Any]]] = Field(..., description="The examples of the tool.")

    # Validate the entire 'examples' list once, *after* parsing.
    @field_validator("examples", mode="after")
    def check_example_keys(cls, examples_list, info: ValidationInfo):
        tool_class = info.context.get('tool_class')
        code_input_sig = tool_class.code_input_sig
        code_input_keys = set([k for k, _, _, _ in code_input_sig])
        text_input_sig = tool_class.text_input_sig
        text_input_keys = set([k for k, _, _, _ in text_input_sig])
        output_keys = set([k for k, _, _ in tool_class.output_sig])

        # Check if the keys match the tool's input signatures
        for ex in examples_list:
            if not {"text_input","code_input","output"}.issubset(ex.keys()):
                raise ValueError("each example must have text_input, code_input and output")
            
            # Check if the keys match the tool's input signatures
            if code_input_keys != set(ex['code_input'].keys()):
                raise ValueError("the keys of code_input must match the tool's input signatures: %s vs %s" % (code_input_keys, set(ex['code_input'].keys())))
            if text_input_keys != set(ex['text_input'].keys()):
                raise ValueError("the keys of text_input must match the tool's input signatures: %s vs %s" % (text_input_keys, set(ex['text_input'].keys())))
            if output_keys != set(ex['output'].keys()):
                raise ValueError("the keys of output must match the tool's output signatures: %s vs %s" % (output_keys, set(ex['output'].keys())))

        return examples_list
    
    @field_validator("name", mode="after")
    def check_name(cls, name, info: ValidationInfo):
        tool_class = info.context.get('tool_class')

        if not name.isidentifier():
            raise ValueError("the name of the tool must be a valid identifier")
        
        if name != tool_class.__name__:
            raise ValueError("the name of the tool must be the same as the class name")

        return name
    
    # allow populating 'version' via the class’s __version__ attribute
    model_config = ConfigDict(populate_by_name=True)


class BaseTool(ABC):
    _registered_tool = True
    _registered_mcp_tool = False
    
    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)

        # If this class is still abstract (i.e. it has abstract methods left unimplemented), skip the metadata check entirely.
        if inspect.isabstract(cls) or getattr(cls, '__abstract__', False):
            if hasattr(cls, '__abstract__'):
                # remove it so it doesn't propagate
                delattr(cls, "__abstract__")
            return

        # Pydantic v2: model_fields maps our declared names → FieldInfo
        meta_dict = {}
        for name, field_info in ToolMeta.model_fields.items():
            # field_info.alias is "__version__" for our version field; 
            # alias==None otherwise, so fallback to `name`.
            key = field_info.alias or name

            if not hasattr(cls, key):
                raise AttributeError(
                    f"{cls.__name__} is missing required class attribute `{key}`"
                )
            meta_dict[key] = getattr(cls, key)

        # one shot validation
        try:
            cls._meta = ToolMeta.model_validate(
                meta_dict, 
                context={'tool_class': cls}
            )
        except ValidationError as e:
            error_msg = str(e)
            error_msg = error_msg.split("\n")
            error_msg[0] = error_msg[0].replace(f"for {ToolMeta.__name__}", f"for {cls.__name__} ({cls.__module__})")
            error_msg = "\n".join(error_msg)
            raise ChemMCPToolMetadataError(error_msg) from e

    @classmethod
    def get_doc(cls, interface='code'):
        assert interface in ('text', 'code'), "Interface '%s' is not supported. Please use 'text' or 'code'." % interface
        inputs = ""
        for name, type_, default, description in (cls.code_input_sig if interface == 'code' else cls.text_input_sig):
            inputs += f"    {name} ({type_}) [default: {default}]: {description}\n"
        outputs = ""
        for name, type_, description in cls.output_sig:
            outputs += f"    {name} ({type_}): {description}\n"
        
        doc = f"""{cls.description}

Args:
{inputs}
Returns:
{outputs}
"""
        return doc

    def __init__(self, init=True, interface='code') -> None:
        assert interface in ('text', 'code'), "Interface '%s' is not supported. Please use 'text' or 'code'." % interface
        self.interface = interface
        super().__init__()
        if init:
            self._init_modules()

    def _init_modules(self):
        pass

    def __call__(self, *args, **kwargs):
        logger.debug("Calling `{}` in __call__".format(self.__class__.name))
        if self.interface == 'text':
            r = self.run_text(args[0])
        elif self.interface == 'code':
            r = self.run_code(*args, **kwargs)
        else:
            raise NotImplementedError("Interface '%s' is not supported. Please use 'text' or 'code'." % self.interface)
        logger.debug("Ending `{}` in __call__".format(self.__class__.name))
        return r

    def run_text(self, query, *args, **kwargs):
        return self._run_text(query, *args, **kwargs)
    
    def _run_text(self, query, *args, **kwargs):
        # Check the signature of self._run_base
        sig = inspect.signature(self._run_base)
        params = list(sig.parameters.values())
        
        # If the number of parameters is 1, and the type is str, then we assume the tool is text-compatible
        if len(params) == 1 and params[0].annotation == str:
            return self._run_base(query)
        else:
            raise NotImplementedError("Text interface is not implemented for this tool yet.")
    
    def run_code(self, *args, **kwargs):
        return self._run_code(*args, **kwargs)
    
    def _run_code(self, *args, **kwargs):
        return self._run_base(*args, **kwargs)
    
    @abstractmethod
    def _run_base(self, *args, **kwargs):
        raise NotImplementedError
