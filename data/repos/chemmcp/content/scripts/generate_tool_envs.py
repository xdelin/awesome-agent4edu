import os
import sys
import json
import argparse
from typing import List

# figure out “project_root/src” relative to this script:
HERE = os.path.dirname(__file__)
SRC  = os.path.abspath(os.path.join(HERE, "..", "src"))
sys.path.insert(0, SRC)

from chemmcp.tools import __all__ as all_tools
import chemmcp


def generate_tool_envs(tool_names: List[str], save_path: str='site/static/data/tool_envs/all_tool_envs.json'):
    tools_envs_dict = {}
    for tool_name in tool_names:
        tool_cls = getattr(chemmcp.tools, tool_name)
        if not tool_cls._registered_mcp_tool:
            # Don't include tools that are not registered as MCP tools
            # This is for quick config for the MCP mode
            continue

        tool_envs = tool_cls.required_envs

        tool_envs_dict = {}
        for env_name, env_description in tool_envs:
            tool_envs_dict[env_name] = env_description
        
        tools_envs_dict[tool_name] = tool_envs_dict

    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    with open(save_path, 'w') as f:
        json.dump(tools_envs_dict, f, indent=4)
    
    return tools_envs_dict


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--tools', type=str, nargs='+', help='The tools to gather the envs for.')
    parser.add_argument('--save_path', type=str, default='site/static/data/tool_envs/all_tool_envs.json', help='The path to save the tool envs.')
    return parser.parse_args()


def main():
    args = parse_args()
    if args.tools is None:
        tools = all_tools
    else:
        tools = args.tools

    generate_tool_envs(tools, save_path=args.save_path)


if __name__ == '__main__':
    main()
