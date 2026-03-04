import functools
import inspect


class ChemMCPError(Exception): ...
"""Errors related to the ChemMCP library."""

class ChemMCPToolMetadataError(ChemMCPError): ...
"""Errors related to the metadata of a ChemMCP tool."""

class ChemMCPApiNotFoundError(ChemMCPError): ...
"""The API key for an external resource is not found or correctly set."""

class ChemMCPToolInitError(ChemMCPError): ...
"""The tool cannot be initialized, e.g., checkpoint file not found or cannot be loaded."""

class ChemMCPRemoteServerDownError(ChemMCPError): ...
"""The remote server is down or unreachable."""

class ChemMCPToolProcessError(ChemMCPError): ...
"""The tool cannot correctly run or return a result due to a runtime error."""

class ChemMCPInputError(ChemMCPError): ...
"""The input to a tool is invalid."""

class ChemMCPSearchFailError(ChemMCPError): ...
"""The search on an external resource cannot get reasonable results."""


def catch_errors(func):
    """Wrap sync or async func so that MyError becomes a returned string."""
    if inspect.iscoroutinefunction(func):
        @functools.wraps(func)
        async def wrapper(*args, **kwargs):
            try:
                return await func(*args, **kwargs)
            except ChemMCPError as e:
                return f"Error ({e.__class__.__name__}): {e}"
        return wrapper
    else:
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            try:
                return func(*args, **kwargs)
            except ChemMCPError as e:
                return f"Error ({e.__class__.__name__}): {e}"
        return wrapper
