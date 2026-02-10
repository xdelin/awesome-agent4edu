import os
import sys
import json
import yaml
import argparse

# figure out “project_root/src” relative to this script:
HERE = os.path.dirname(__file__)
SRC  = os.path.abspath(os.path.join(HERE, "..", "src"))
sys.path.insert(0, SRC)

from chemmcp.tools import __all__ as all_tools
import chemmcp
from .modules_utils import get_class_file_path, git_last_commit_date_for_class


def get_tool_metainfo(tool_cls: type):
    frontmatter = {
        "title": f"{tool_cls.name} ({tool_cls.func_name})",
        "categories": list(tool_cls.categories),
        "tags": list(tool_cls.tags),
        "description": tool_cls.description,
        "weight": 2,
        "draft": False,
    }

    txt = f"""---
{yaml.dump(frontmatter)}
---
"""

    return txt


def get_tool_badges(tool_cls: type):
    # get the last updated date of the file
    last_updated_date = git_last_commit_date_for_class(tool_cls)

    # get the badges
    badges = []
    badges.append(f"Version: {tool_cls.__version__}")
    badges.append(f"Last Update: {last_updated_date}")
    badges.append("%sMCP Support" % ("" if tool_cls._registered_mcp_tool else "No "))
    badges.append("%sPython Calling Support" % ("" if tool_cls._registered_tool else "No "))
    
    badges_txt = """<div style="display: flex; flex-wrap: wrap; gap: 0.75rem; align-items: center;">
%s
</div>
""" % "\n".join(["  {{< badge >}}%s{{< /badge >}}" % badge for badge in badges])
    
    return badges_txt


def get_tool_description(tool_cls: type):
    description_txt = '{{< lead >}}\n**' + tool_cls.description + '**\n{{< /lead >}}\n\n'
    return description_txt


def get_tool_examples(tool_cls: type):
    tool_examples = tool_cls.examples
    if len(tool_examples) == 0:
        return ""
    
    def get_example_input_lines(tool_example_input: dict):
        example_input_lines = []
        for key, value in tool_example_input.items():
            example_input_lines.append(f"{key}: {repr(value)}")
        example_input_txt = "```yaml\n" + "\n".join(example_input_lines) + "\n```"
        return example_input_txt
    
    def get_example_output_lines(tool_example_output: dict):
        example_output_lines = []
        for key, value in tool_example_output.items():
            if isinstance(value, str):
                value_lines = value.strip().split('\n')
                if len(value_lines) > 1:
                    example_output_lines.append('"""' + value_lines[0])
                    for line in value_lines[1:-1]:
                        example_output_lines.append(line)
                    example_output_lines.append(value_lines[-1] + '"""')
                else:
                    example_output_lines.append(f"{repr(value)}")
            else:
                example_output_lines.append(f"{repr(value)}")
        example_output_txt = "```\n" + "\n".join(example_output_lines) + "\n```"
        return example_output_txt
    
    examples_txt = ""
    for idx, tool_example in enumerate(tool_examples, start=1):
        example_title = "Example" + ("" if len(tool_examples) == 1 else f" {idx}")
        example_title = "**" + example_title + "**"
        
        example_txt = """%s

Input:
%s

Text Input (used for the `run_text` function in the Python calling mode):
%s

Output:
%s

""" % (
        example_title,
        get_example_input_lines(tool_example['code_input']),
        get_example_input_lines(tool_example['text_input']),
        get_example_output_lines(tool_example['output'])
    )
    
    examples_txt += example_txt

    return examples_txt


def get_tool_mcp_usage(tool_cls: type, required_envs: dict, requires_llms: bool):
    txt = 'Configure your MCP client following its instructions with something like:\n'

    envs = {}
    for env in required_envs:
        envs[env] = "VALUE_TO_BE_SET"

    if requires_llms:
        envs["LLM_MODEL_NAME"] = "VALUE_TO_BE_SET"

    envs_txt = json.dumps(envs, indent=4)
    envs_txt_lines = envs_txt.strip().split("\n")
    for idx, line in enumerate(envs_txt_lines):
        if idx > 0:
            line = "    " + line
        envs_txt_lines[idx] = line

    if requires_llms:
        envs_txt_lines = envs_txt_lines[:-1] + ["        // Add required LLM credentials", "        // ..."] + envs_txt_lines[-1:]

    envs_txt = "\n".join(envs_txt_lines)

    txt += f"""```JSON
{{
    "command": "/ABSTRACT/PATH/TO/uv",  // Use `which uv` to get its path
    "args": ["--directory", "/ABSTRACT/PATH/TO/ChemMCP", "run", "--tools", "{tool_cls.name}"],
    "toolCallTimeoutMillis": 300000,
    "env": {envs_txt}
}}
```"""
        
    return txt


def get_tool_python_usage(tool_cls: type, required_envs: dict, requires_llms: bool):
    envs_setup = []
    for env in required_envs:
        envs_setup.append(f"os.environ['{env}'] = 'VALUE_TO_BE_SET'")

    if requires_llms:
        envs_setup.append("os.environ['LLM_MODEL_NAME'] = 'VALUE_TO_BE_SET'")
        envs_setup.append("# Also add LLM credentials required by the LLM model, such as `OPENAI_API_KEY`")

    if len(envs_setup) > 0:
        envs_setup = ["# Set the environment variables"] + envs_setup
        envs_setup = "\n".join(envs_setup) + '\n\n'
    else:
        envs_setup = ""

    code_input_signature = tool_cls.code_input_sig
    code_input_example = tool_cls.examples[0]['code_input']
    run_code_inputs = []
    for param_name, param_type, param_default, param_description in code_input_signature:
        run_code_inputs.append(f'    {param_name}={repr(code_input_example[param_name])}')
    run_code_inputs = '\n'.join(run_code_inputs)

    text_input_signature = tool_cls.text_input_sig
    text_input_example = tool_cls.examples[0]['text_input']
    run_text_inputs = []
    for param_name, param_type, param_default, param_description in text_input_signature:
        run_text_inputs.append(f'    {param_name}={repr(text_input_example[param_name])}')
    run_text_inputs = '\n'.join(run_text_inputs)

    txt = """```python
import os
from chemmcp.tools import {tool_name}

{envs_setup}# Initialize the tool
tool = {tool_name}()

# The tool has two alternative ways to run:
# 1. Run with separate input domains (recommended)
output = tool.run_code(
{run_code_inputs}
)
# 2. Run with text-only input
output = tool.run_text(
{run_text_inputs}
)
```
""".format(tool_name=tool_cls.name, envs_setup=envs_setup, run_code_inputs=run_code_inputs, run_text_inputs=run_text_inputs)
    
    txt += """

Each tool in ChemMCP has two ways to run:
- **`run_code`** (recommended): The inputs contain one or more domains, each of which can be a str, an int, a float, etc.
- **`run_text`**: The inputs are a single string in a specific format. The tool will parse the string to extract the input domains. This is useful in scenarios where an agent framework calls tools only with text input.
The output is the same in both cases.

For the input and output domains, please refer to the tool's [signature](#tool-signature).
"""
    return txt


def get_tool_usage(tool_cls: type):
    if tool_cls._registered_mcp_tool and tool_cls._registered_tool:
        support_description = "The tool supports both [MCP mode](#mcp-mode) and [Python calling mode](#python-calling-mode)."
    elif tool_cls._registered_mcp_tool:
        support_description = "The tool supports [MCP mode](#mcp-mode)."
    elif tool_cls._registered_tool:
        support_description = "The tool supports [Python calling mode](#python-calling-mode)."
    else:
        raise ValueError(f"Tool {tool_cls.name} is not registered as either MCP tool or Python tool.")
    
    # Get the instructions for the environment variables
    required_envs = dict(tool_cls.required_envs)
    try:
        requires_llms = required_envs.pop("__llms__")
    except KeyError:
        requires_llms = False
    else:
        requires_llms = True

    envs_instructions = ""
    if len(required_envs) > 0 or requires_llms:
        envs_instructions = "### Environment Variables\n"

        envs_instructions += "This tool requires the following environment variables:\n"
        for env, env_description in required_envs.items():
            envs_instructions += f"- **{env}**: {env_description}\n"

        if requires_llms:
            envs_instructions += f"- **LLM_MODEL_NAME**: The name of the LLM to use. See [LiteLLM](https://docs.litellm.ai/docs/#basic-usage) for more details.\n"
            envs_instructions += f"- Other LLM credentials are required to be set in the `env` field. See [LiteLLM](https://docs.litellm.ai/docs/#basic-usage) for more details.\n"
    
    if tool_cls._registered_mcp_tool:
        mcp_usage = get_tool_mcp_usage(tool_cls, required_envs, requires_llms)
    else:
        mcp_usage = "This tool does not support MCP mode."

    if tool_cls._registered_tool:
        python_usage = get_tool_python_usage(tool_cls, required_envs, requires_llms)
    else:
        python_usage = "This tool does not support Python calling mode."

    txt = """## Usage

{support_description}

{envs_instructions}

### MCP Mode

{mcp_usage}

### Python Calling Mode

{python_usage}
""".format(support_description=support_description, envs_instructions=envs_instructions, mcp_usage=mcp_usage, python_usage=python_usage)
    
    return txt

def get_tool_signature(tool_cls: type):
    input_table_lines = ["| Name | Type | Default | Description |", "| --- | --- | --- | --- |"]
    for param_name, param_type, param_default, param_description in tool_cls.code_input_sig:
        input_table_lines.append(f"| {param_name} | {param_type} | {param_default} | {param_description} |")
    input_table = "\n".join(input_table_lines)

    text_input_table_lines = ["| Name | Type | Default | Description |", "| --- | --- | --- | --- |"]
    for param_name, param_type, param_default, param_description in tool_cls.text_input_sig:
        text_input_table_lines.append(f"| {param_name} | {param_type} | {param_default} | {param_description} |")
    text_input_table = "\n".join(text_input_table_lines)

    output_table_lines = ["| Name | Type | Description |", "| --- | --- | --- |"]
    for param_name, param_type, param_description in tool_cls.output_sig:
        output_table_lines.append(f"| {param_name} | {param_type} | {param_description} |")
    output_table = "\n".join(output_table_lines)

    # Get the instructions for the environment variables
    required_envs = dict(tool_cls.required_envs)
    try:
        required_envs.pop("__llms__")
    except KeyError:
        pass
    else:
        required_envs["LLM_MODEL_NAME"] = "The name of the LLM to use. See [LiteLLM](https://docs.litellm.ai/docs/#basic-usage) for more details."
        required_envs["LLM credentials"] = "Other LLM credentials are required to be set in the `env` field. See [LiteLLM](https://docs.litellm.ai/docs/#basic-usage) for more details."

    if len(required_envs) > 0:
        envs_table_lines = ["| Name | Description |", "| --- | --- |"]
        for env, env_description in required_envs.items():
            envs_table_lines.append(f"| {env} | {env_description} |")
        envs_table = "\n".join(envs_table_lines)
    else:
        envs_table = "No required environment variables for this tool."

    signature_txt = """## Tool Signature\n\n

### Input
Used in the MCP mode, as well as the `run_code` function in the Python calling mode.
{input_table}

### Text Input
Used in the `run_text` function in the Python calling mode.
{text_input_table}

### Output
The output is the same in both input cases.
{output_table}

### Envs
{envs_table}
""".format(input_table=input_table, text_input_table=text_input_table, output_table=output_table, envs_table=envs_table)
    return signature_txt

def get_tool_implementation_details(tool_cls: type):
    txt = f"""## Implementation Details

- **Implementation Description**: {tool_cls.implementation_description}
"""
    oss_dependencies = tool_cls.oss_dependencies
    if len(oss_dependencies) > 0:
        oss_dependencies_txt = "- **Open-source dependencies** (code source or required libraries):\n"
        for dependency in oss_dependencies:
            oss_dependencies_txt += f"  - [{dependency[0]}]({dependency[1]}) ({'Unknown license' if dependency[2] is None else dependency[2]})\n"
        txt += oss_dependencies_txt
    else:
        txt += "- **Open-source dependencies** (code source or required libraries): None\n"

    services_and_software = tool_cls.services_and_software
    if len(services_and_software) > 0:
        services_and_software_txt = "- **Hosted services and software** (required for running the tool):\n"
        for service in services_and_software:
            if service[1] is None:
                services_and_software_txt += f"  - {service[0]}\n"
            else:
                services_and_software_txt += f"  - [{service[0]}]({service[1]})\n"
        txt += services_and_software_txt
    else:
        txt += "- **Hosted services and software** (required for running the tool): None\n"
    return txt

def generate_tool_doc(tool_name: str, save_dir: str='site/content/tools/'):
    tool_cls = getattr(chemmcp.tools, tool_name)
    tool_file_path = get_class_file_path(tool_cls)
    tool_file_base_name = os.path.basename(tool_file_path)
    tool_file_name = os.path.splitext(tool_file_base_name)[0]

    doc = ""
    doc += get_tool_metainfo(tool_cls)
    doc += get_tool_badges(tool_cls)
    doc += get_tool_description(tool_cls)
    doc += get_tool_examples(tool_cls)
    doc += get_tool_usage(tool_cls)
    doc += get_tool_signature(tool_cls)
    doc += get_tool_implementation_details(tool_cls)
    
    os.makedirs(save_dir, exist_ok=True)
    doc_file_path = os.path.join(save_dir, tool_file_name + '.md')
    with open(doc_file_path, 'w') as f:
        f.write(doc)

    return doc


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--tools', type=str, nargs='+', help='The tools to generate the documentation for.')
    parser.add_argument('--save-dir', type=str, default='site/content/tools/', help='The directory to save the tool documentation.')
    return parser.parse_args()


def main():
    args = parse_args()
    if args.tools is None:
        tools = all_tools
    else:
        tools = args.tools

    for tool_name in tools:
        generate_tool_doc(tool_name, save_dir=args.save_dir)


if __name__ == '__main__':
    main()
