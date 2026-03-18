from mcp.server.fastmcp import FastMCP
from kaggle.api.kaggle_api_extended import KaggleApi
from typing import Any

mcp = FastMCP("Kaggle")

api = KaggleApi()
api.authenticate()

import json

def make_json_safe(data):
    return json.loads(str(data))


#   -------------------------
#   Competition Related Tools
#   -------------------------
    

@mcp.tool()
def get_competitions_list() -> Any:
    """
    Fetch the list of available Kaggle competitions.

    Returns:
        list: A list of competitions.
        str: An error message if the request fails.
    """
    try:
        competitions_raw = api.competitions_list()
        competitions_list = [make_json_safe(competition) for competition in competitions_raw]
        return competitions_list
    except Exception as e:
        return f"Exception occurred: {e}"


@mcp.tool()
def get_competition_leaderboard_view(competition_name: str) -> list | str:
    """
    Fetch the leaderboard for a specific competition.

    Args:
        competition_name (str): The name of the competition.

    Returns:
        list: Leaderboard data for the competition.
        str: An error message if the request fails.
    """
    try:
        competitons_leaderboard_raw = api.competition_leaderboard_view(competition=competition_name)
        competitions_leaderboard_list = [make_json_safe(competition) for competition in competitons_leaderboard_raw]
        return competitions_leaderboard_list
    except Exception as e:
        return f"Exception occurred: {e}"


@mcp.tool()
def get_valid_competition_categories() -> list | str:
    """
    Fetch the valid categories for competitions.

    Returns:
        list: A list of competition categories.
        str: An error message if the request fails.
    """
    try:
        return api.valid_competition_categories
    except Exception as e:
        return f"Exception occurred: {e}"


@mcp.tool()
def get_valid_competition_groups() -> list | str:
    """
    Fetch the valid groups for competitions.

    Returns:
        list: A list of competition groups.
        str: An error message if the request fails.
    """
    try:
        return api.valid_competition_groups
    except Exception as e:
        return f"Exception occurred: {e}"


#    ---------------------
#    Dataset Related Tools
#    ---------------------

@mcp.tool()
def get_dataset_list(dataset_category: str) -> list | str:
    """
    Fetch the list of datasets for a given category.

    Args:
        dataset_category (str): The category of datasets
            (e.g., Finance, Healthcare, Education).

    Returns:
        list: A list of datasets.
        str: An error message if the request fails.
    """
    try:
        return api.dataset_list(search=dataset_category)
    except Exception as e:
        return f"Exception occurred: {e}"


@mcp.tool()
def get_dataset_labels() -> list | str:
    """
    Fetch the list of dataset labels.

    Returns:
        list: A list of dataset labels.
        str: An error message if the request fails.
    """
    try:
        return api.dataset_labels
    except Exception as e:
        return f"Exception occurred: {e}"


#     -------------------
#     Model Related Tools
#     -------------------

@mcp.tool()
def get_model_list() -> list | str:
    """
    Fetch the list of available models.

    Returns:
        list: A list of models.
        str: An error message if the request fails.
    """
    try:
        return api.model_list()
    except Exception as e:
        return f"Exception occurred: {e}"


@mcp.tool()
def get_model_get(model_name: str):
    """
    Fetch details of a specific model.

    Args:
        model_name (str): The name of the model.

    Returns:
        dict: Model details.
        str: An error message if the request fails.
    """
    try:
        return api.model_get(model=model_name)
    except Exception as e:
        return f"Exception occurred: {e}"


#     ---------------------
#     Tools for Valid Types
#     ---------------------

@mcp.tool()
def get_valid_list_languages() -> list | str:
    """
    Fetch the list of valid languages.

    Returns:
        list: A list of languages.
        str: An error message if the request fails.
    """
    try:
        return api.valid_list_languages
    except Exception as e:
        return f"Exception occurred: {e}"


@mcp.tool()
def get_valid_list_kernel_types() -> list | str:
    """
    Fetch the list of valid kernel types.

    Returns:
        list: A list of kernel types.
        str: An error message if the request fails.
    """
    try:
        return api.valid_list_kernel_types
    except Exception as e:
        return f"Exception occurred: {e}"


@mcp.tool()
def get_valid_list_output_types() -> list | str:
    """
    Fetch the list of valid output types.

    Returns:
        list: A list of output types.
        str: An error message if the request fails.
    """
    try:
        return api.valid_list_output_types
    except Exception as e:
        return f"Exception occurred: {e}"


if __name__ == "__main__":
    print("Starting MCP server...")
    mcp.run(transport='stdio')
