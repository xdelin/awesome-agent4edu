
from mcp.server.fastmcp import Context, FastMCP
from pydantic import Field

from mcp_apple_notes.applescript import ValidationUtils
from mcp_apple_notes.tools import NotesTools

# Initialize FastMCP server
mcp = FastMCP(name="mcp-apple-notes")

# Initialize tools
notes_tools = NotesTools()


@mcp.tool()
async def create_note(
    ctx: Context,
    name: str = Field(
        ...,
        description="Note title wrapped in <h1> tags (e.g., '<h1>My Note Title</h1>')",
    ),
    body: str = Field(
        ...,
        description="Note body content with appropriate HTML formatting (e.g., '<p>Content here</p>'). For proper spacing between two sections, use <br>.",
    ),
    folder_path: str = Field(
        default="Notes",
        description="Target folder path (e.g., 'Work' or 'Work/Projects/2024'). Defaults to 'Notes'",
    ),
) -> str:
    """Create a new note with specified name and content.

    Supported HTML: <h1-h6>, <b><i><u>, <p><div><br>, <ul><ol><li>, <table>, <a>
    Best Practices: Use semantic HTML, add <br> tags for spacing, avoid CSS styles
    Folders: Root level or nested paths (up to 5 levels deep)
    Limitations: No special chars in name, no complex CSS/JS

    Example:
    name: "<h1>Project Report</h1>"
    body: "<p>Status: <b>In Progress</b></p><br><ul><li>Task 1</li></ul>"
    """
    try:
        # Validate the title content in the name parameter
        ValidationUtils.validate_html_title_content(name)

        # Combine name (title) and body into complete HTML content with <br> spacing
        combined_content = name + "<br>" + body
        note = await notes_tools.create_note("note", combined_content, folder_path)
        return str(note)
    except ValueError as e:
        # Handle validation errors with clear messages
        error_msg = f"Invalid input: {str(e)}"
        await ctx.error(error_msg)
        raise ValueError(error_msg)
    except Exception as e:
        await ctx.error(f"Error creating note: {str(e)}")
        raise


@mcp.tool()
async def create_folder(
    ctx: Context,
    folder_name: str = Field(
        ..., description='Name of the folder to create (1-128 chars, no < > : " | ? *)'
    ),
    folder_path: str = Field(
        default="",
        description="Optional nested path (e.g., 'Work/Projects'). If empty, creates at root level. Max 5 levels deep.",
    ),
) -> str:
    """Create a folder in Apple Notes.

    Features:
    - Creates folders at root level or nested paths (up to 5 levels deep)
    - Unicode and emoji support for international characters
    - Duplicate name detection and comprehensive validation
    - iCloud account support
    - Only creates folder if parent path exists

    Validation:
    - Max 128 characters, no special chars: < > : " | ? *
    - Parent paths must exist, prevents duplicates
    - Will throw error if parent path doesn't exist
    """
    try:
        folder = await notes_tools.create_folder(folder_name, folder_path)

        # Return only name and primary key ID
        result = f"Folder Name: {folder.get('name', 'N/A')}\n"
        result += f"Folder ID: {folder.get('id', 'N/A')}\n"

        return result
    except ValueError as e:
        # Handle validation errors with clear messages
        error_msg = f"Invalid input: {str(e)}"
        await ctx.error(error_msg)
        raise ValueError(error_msg)
    except RuntimeError as e:
        # Handle AppleScript errors with helpful context
        error_msg = str(e)
        await ctx.error(error_msg)
        raise RuntimeError(error_msg)
    except Exception as e:
        # Handle unexpected errors
        error_msg = f"Unexpected error creating folder '{folder_name}' in path '{folder_path}': {str(e)}"
        await ctx.error(error_msg)
        raise


@mcp.tool()
async def read_note(
    ctx: Context,
    note_id: str = Field(
        ..., description="Primary key ID of the note to read (e.g., 'p1308')"
    ),
    note_name: str = Field(..., description="Name of the note to verify and read"),
) -> str:
    """Read a note by its primary key ID with AppleScript verification.

    Security Features:
    - Verifies note exists with the given ID and name
    - Uses primary key ID for precise identification
    - Returns full note content with metadata

    Output:
    - Note name, ID, folder, creation/modification dates
    - Full note content (title + body)
    - Status and read method information
    """
    try:
        note_data = await notes_tools.read_note(note_id, note_name)

        # Format the response
        result = "Note Content:\n"
        result += f"Note Name: {note_data.get('name', 'N/A')}\n"
        result += f"Note ID: {note_data.get('note_id', 'N/A')}\n"
        result += f"Folder: {note_data.get('folder', 'N/A')}\n\n"
        result += f"Full Content:\n{note_data.get('body', 'No content available')}\n"

        return result

    except ValueError as e:
        error_msg = f"Invalid input: {str(e)}"
        await ctx.error(error_msg)
        raise ValueError(error_msg)
    except RuntimeError as e:
        error_msg = f"Note not found or path error: {str(e)}"
        await ctx.error(error_msg)
        raise RuntimeError(error_msg)
    except Exception as e:
        await ctx.error(f"Error reading note: {str(e)}")
        raise


@mcp.tool()
async def read_folder(
    ctx: Context,
    folder_id: str = Field(
        ..., description="Primary key ID of the folder to read (e.g., 'p2330')"
    ),
    folder_name: str = Field(..., description="Name of the folder to verify and read"),
) -> str:
    """Read a folder by its primary key ID with AppleScript verification.

    Security Features:
    - Verifies folder exists with the given ID and name
    - Uses primary key ID for precise identification
    - Returns detailed folder information with contents

    Output:
    - Folder metadata (name, ID, creation/modification dates)
    - Direct child folders (names and IDs)
    - Notes in the folder (names and IDs)
    - Summary counts
    """
    try:
        folder_data = await notes_tools.read_folder(folder_id, folder_name)

        # Format the response
        result = "Folder Details:\n"
        result += f"Folder Name: {folder_data.get('name', 'N/A')}\n"
        result += f"Folder ID: {folder_data.get('folder_id', 'N/A')}\n"
        result += f"Creation Date: {folder_data.get('creation_date', 'N/A')}\n"
        result += (
            f"Modification Date: {folder_data.get('modification_date', 'N/A')}\n\n"
        )

        # Child folders section
        child_folders = folder_data.get("child_folders", [])
        child_folders_count = folder_data.get("child_folders_count", 0)
        result += f"Direct Child Folders ({child_folders_count}):\n"

        if child_folders:
            for i, child_folder in enumerate(child_folders, 1):
                result += f"  {i:2d}. {child_folder.get('name', 'N/A')}\n"
                result += f"      ID: {child_folder.get('id', 'N/A')}\n"
        else:
            result += "  No child folders\n"

        result += "\n"

        # Notes section
        notes = folder_data.get("notes", [])
        notes_count = folder_data.get("notes_count", 0)
        result += f"Notes ({notes_count}):\n"

        if notes:
            for i, note in enumerate(notes, 1):
                result += f"  {i:2d}. {note.get('name', 'N/A')}\n"
                result += f"      ID: {note.get('note_id', 'N/A')}\n"
        else:
            result += "  No notes\n"

        return result

    except ValueError as e:
        error_msg = f"Invalid input: {str(e)}"
        await ctx.error(error_msg)
        raise ValueError(error_msg)
    except RuntimeError as e:
        error_msg = f"Folder not found or verification error: {str(e)}"
        await ctx.error(error_msg)
        raise RuntimeError(error_msg)
    except Exception as e:
        await ctx.error(f"Error reading folder: {str(e)}")
        raise


@mcp.tool()
async def update_note(
    ctx: Context,
    note_id: str = Field(
        ..., description="Primary key ID of the note to update (e.g., 'p1234')"
    ),
    note_name: str = Field(
        ..., description="Current name of the note to verify and update"
    ),
    new_name: str = Field(
        ...,
        description="New note title wrapped in <h1> tags (e.g., '<h1>Updated Title</h1>')",
    ),
    new_body: str = Field(
        ...,
        description="New note body content with appropriate HTML formatting (e.g., '<p>Updated content</p>') For proper spacing between two sections, use <br>.",
    ),
) -> str:
    """Update an existing note by its primary key ID with AppleScript verification.

    Required: note_id, note_name, new_name, and new_body parameters
    Supported HTML: <h1-h6>, <b><i><u>, <p><div><br>, <ul><ol><li>, <table>, <a>
    Best Practices: Use semantic HTML, add <br> tags for spacing, avoid CSS styles
    Security: AppleScript verifies ID and name match before updating

    Example:
    note_id: "p1234"
    note_name: "Current Note Title"
    new_name: "<h1>Updated Report</h1>"
    new_body: "<p>Status: <b>Complete</b></p><br><ul><li>Done</li></ul>"
    """
    try:
        # Validate the title content in the new_name parameter
        ValidationUtils.validate_html_title_content(new_name)

        # Combine new_name (title) and new_body into complete HTML content with <br> spacing
        combined_content = new_name + "<br>" + new_body
        updated_note = await notes_tools.update_note(
            note_id, note_name, combined_content
        )

        # Format the response with primary key ID
        result = "Note Update Result:\n"
        result += f"Note Name: {updated_note.get('name', 'N/A')}\n"
        result += f"Note ID: {updated_note.get('note_id', 'N/A')}\n"
        result += f"Status: {updated_note.get('status', 'N/A')}\n"

        return result

    except ValueError as e:
        error_msg = f"Invalid input: {str(e)}"
        await ctx.error(error_msg)
        raise ValueError(error_msg)
    except RuntimeError as e:
        error_msg = f"Note not found or path error: {str(e)}"
        await ctx.error(error_msg)
        raise RuntimeError(error_msg)
    except Exception as e:
        await ctx.error(f"Error updating note: {str(e)}")
        raise


@mcp.tool()
async def rename_folder(
    ctx: Context,
    folder_id: str = Field(
        ..., description="Primary key ID of the folder to rename (e.g., 'p2330')"
    ),
    current_name: str = Field(
        ..., description="Current name of the folder to verify and rename"
    ),
    new_name: str = Field(
        ..., description='New name for the folder (1-128 chars, no < > : " | ? *)'
    ),
) -> str:
    """Rename a folder in Apple Notes by ID with enhanced name verification.

    Features:
    - Renames folders by ID with double name verification
    - Comprehensive validation and duplicate detection
    - Unicode and emoji support for folder names
    - Works with root level and nested paths

    Validation:
    - Max 128 characters, no special chars: < > : " | ? *
    - Prevents duplicate names in same location
    - New name cannot be same as current name
    - Verifies folder ID matches current name (double verification)
    - Gets actual folder name by ID for additional security
    """
    try:
        rename_result = await notes_tools.rename_folder(
            folder_id, current_name, new_name
        )

        # Format the response
        result = "Folder Rename Result:\n"
        result += f"Folder ID: {rename_result.get('folder_id', 'N/A')}\n"
        result += f"Old Name: {rename_result.get('current_name', 'N/A')}\n"
        result += f"New Name: {rename_result.get('new_name', 'N/A')}\n"
        result += f"Path: {rename_result.get('folder_path', 'N/A')}\n"
        result += f"Status: {rename_result.get('status', 'N/A')}\n"

        return result
    except ValueError as e:
        # Handle validation errors with clear messages
        error_msg = f"Invalid input: {str(e)}"
        await ctx.error(error_msg)
        raise ValueError(error_msg)
    except RuntimeError as e:
        # Handle AppleScript errors with helpful context
        error_msg = str(e)
        await ctx.error(error_msg)
        raise RuntimeError(error_msg)
    except Exception as e:
        # Handle unexpected errors
        error_msg = f"Unexpected error renaming folder '{current_name}' to '{new_name}': {str(e)}"
        await ctx.error(error_msg)
        raise


@mcp.tool()
async def move_folder(
    ctx: Context,
    folder_id: str = Field(
        ..., description="Primary key ID of the folder to move (e.g., 'p2330')"
    ),
    folder_name: str = Field(..., description="Name of the folder to verify and move"),
    target_path: str = Field(
        default="",
        description="Target path where to move the folder (e.g., 'Archive'). If empty, moves to root level.",
    ),
) -> str:
    """Move a folder from one location to another in Apple Notes by ID with AppleScript verification.

    Features:
    - Moves folders by ID with AppleScript verification
    - Comprehensive validation and duplicate detection
    - Unicode and emoji support for folder names
    - Works with root level and nested paths

    Validation:
    - Validates target path exists
    - Prevents duplicate names in target location
    - Enforces 5-level nesting depth limit
    - AppleScript verifies ID and name match the same folder
    """
    try:
        move_result = await notes_tools.move_folder(folder_id, folder_name, target_path)

        # Format the response
        result = "Folder Move Result:\n"
        result += f"Folder ID: {move_result.get('folder_id', 'N/A')}\n"
        result += f"Folder Name: {move_result.get('folder_name', 'N/A')}\n"
        result += f"Source Path: {move_result.get('source_path', 'N/A')}\n"
        result += f"Target Path: {move_result.get('target_path', 'N/A')}\n"
        result += f"Status: {move_result.get('status', 'N/A')}\n"

        return result
    except ValueError as e:
        # Handle validation errors with clear messages
        error_msg = f"Invalid input: {str(e)}"
        await ctx.error(error_msg)
        raise ValueError(error_msg)
    except RuntimeError as e:
        # Handle AppleScript errors with helpful context
        error_msg = str(e)
        await ctx.error(error_msg)
        raise RuntimeError(error_msg)
    except Exception as e:
        # Handle unexpected errors
        error_msg = f"Unexpected error moving folder '{folder_name}' to '{target_path}': {str(e)}"
        await ctx.error(error_msg)
        raise


@mcp.tool()
async def list_folder_with_structure(ctx: Context) -> str:
    """List the complete folder structure with hierarchical tree format.

    Features:
    - Shows all folders in hierarchical tree format
    - Displays folder nesting levels with visual indicators
    - Works with root level and nested folder structures

    Output Format:
    - Tree structure with ├── and └── indicators
    - Clear hierarchy visualization
    - Folder names with proper indentation

    Returns:
        Hierarchical tree structure of all folders in Apple Notes
    """
    try:
        folder_structure = await notes_tools.list_folder_with_structure()

        if not folder_structure:
            return "No folders found in Apple Notes"

        # Return filtered AppleScript result
        return f"Apple Notes Folder Structure:\n\n{folder_structure}"
    except Exception as e:
        await ctx.error(f"Error listing folder structure: {str(e)}")
        raise


@mcp.tool()
async def delete_note(
    ctx: Context,
    note_id: str = Field(
        ..., description="Primary key ID of the note to delete (e.g., 'p1308')"
    ),
    note_name: str = Field(..., description="Name of the note to verify and delete"),
) -> str:
    """Delete a note by its primary key ID with AppleScript verification.

    Security Features:
    - Verifies note exists with the given ID
    - Confirms note name matches before deletion
    - Uses primary key ID for precise identification
    - AppleScript handles ID and name verification automatically
    - Provides detailed error messages for troubleshooting

    Output:
    - Note name, ID, deletion status and method information
    """
    try:
        deleted_note = await notes_tools.delete_note(note_id, note_name)

        # Format the response
        result = "Note Deletion Result:\n"
        result += f"Note Name: {deleted_note.get('name', 'N/A')}\n"
        result += f"Note ID: {deleted_note.get('note_id', 'N/A')}\n"
        result += f"Status: {deleted_note.get('status', 'N/A')}\n"

        return result

    except ValueError as e:
        error_msg = f"Invalid input: {str(e)}"
        await ctx.error(error_msg)
        raise ValueError(error_msg)
    except RuntimeError as e:
        error_msg = f"Note not found or path error: {str(e)}"
        await ctx.error(error_msg)
        raise RuntimeError(error_msg)
    except Exception as e:
        await ctx.error(f"Error deleting note: {str(e)}")
        raise


@mcp.tool()
async def delete_folder(
    ctx: Context,
    folder_id: str = Field(
        ..., description="Primary key ID of the folder to delete (e.g., 'p2330')"
    ),
    folder_name: str = Field(
        ..., description="Name of the folder to verify and delete"
    ),
) -> str:
    """Delete a folder in Apple Notes by ID with AppleScript verification.

    Security Features:
    - Verifies folder exists with the given ID
    - Confirms folder name matches before deletion
    - Uses primary key ID for precise identification
    - AppleScript handles ID and name verification automatically
    - Provides detailed error messages for troubleshooting

    Output:
    - Folder name, ID, deletion status and method information
    """
    try:
        delete_result = await notes_tools.delete_folder(folder_id, folder_name)

        # Format the response
        result = "Folder Delete Result:\n"
        result += f"Folder Name: {delete_result.get('name', 'N/A')}\n"
        result += f"Folder ID: {delete_result.get('folder_id', 'N/A')}\n"
        result += f"Path: {delete_result.get('folder_path', 'N/A')}\n"
        result += f"Status: {delete_result.get('status', 'N/A')}\n"

        return result

    except ValueError as e:
        error_msg = f"Invalid input: {str(e)}"
        await ctx.error(error_msg)
        raise ValueError(error_msg)
    except RuntimeError as e:
        error_msg = f"Folder not found or name mismatch: {str(e)}"
        await ctx.error(error_msg)
        raise RuntimeError(error_msg)
    except Exception as e:
        await ctx.error(f"Error deleting folder: {str(e)}")
        raise


@mcp.tool()
async def list_notes_with_structure(ctx: Context) -> str:
    """List the complete folder structure with notes included in hierarchical tree format.

    Features:
    - Shows all folders and notes in hierarchical tree format
    - Displays folder nesting levels with visual indicators
    - Lists notes within each folder
    - Works with root level and nested folder structures

    Output Format:
    - Tree structure with ├── and └── indicators
    - Folder names with proper indentation
    - Notes listed under their respective folders

    Returns:
        Hierarchical tree structure of all folders and notes in Apple Notes
    """
    try:
        notes_structure = await notes_tools.list_notes_with_structure()

        if not notes_structure:
            return "No folders or notes found in Apple Notes"

        # Return filtered AppleScript result
        return f"Apple Notes Structure with Notes:\n\n{notes_structure}"
    except Exception as e:
        await ctx.error(f"Error listing notes structure: {str(e)}")
        raise


@mcp.tool()
async def list_all_notes(ctx: Context) -> str:
    """List all notes across all folders with their names and IDs.

    Features:
    - Lists ALL notes from ALL folders in Apple Notes
    - Shows note names, IDs, and folder locations
    - Includes notes from Recently Deleted folder
    - Provides comprehensive system overview

    Output Format:
    - Numbered list of all notes
    - Note names with emoji indicators
    - Note IDs for reference
    - Folder location for each note

    Returns:
        Complete list of all notes across all folders in Apple Notes
    """
    try:
        notes_list = await notes_tools.list_all_notes()

        if not notes_list:
            return "No notes found in Apple Notes"

        # Format the response
        result = f"All Notes ({len(notes_list)} total):\n\n"

        for i, note in enumerate(notes_list, 1):
            result += f"{i:3d}. {note.get('name', 'N/A')}\n"
            result += f"     ID: {note.get('note_id', 'N/A')}\n"
            result += f"     Folder: {note.get('folder', 'N/A')}\n"
            result += "\n"

        return result

    except Exception as e:
        await ctx.error(f"Error listing all notes: {str(e)}")
        raise


@mcp.tool()
async def move_note(
    ctx: Context,
    note_id: str = Field(
        ..., description="Primary key ID of the note to move (e.g., 'p1308')"
    ),
    note_name: str = Field(..., description="Name of the note to verify and move"),
    target_folder_path: str = Field(
        ...,
        description="Target folder path where to move the note (e.g., 'Archive' or 'Work/Completed')",
    ),
) -> str:
    """Move a note from one folder to another in Apple Notes by ID with AppleScript verification.

    Features:
    - Moves notes between root folders and nested paths (up to 5 levels deep)
    - Comprehensive validation and error handling
    - Supports all folder path types (root, simple, nested)
    - Maintains note content and metadata during move
    - AppleScript verifies ID and name match before moving

    Requirements:
    - Note must exist with the given ID and name
    - Target folder path must exist
    - AppleScript handles all verification automatically
    """
    try:
        move_result = await notes_tools.move_note(
            note_id, note_name, target_folder_path
        )

        # Format the response
        result = "Note Move Result:\n"
        result += f"Note Name: {move_result.get('name', 'N/A')}\n"
        result += f"Note ID: {move_result.get('note_id', 'N/A')}\n"
        result += f"Source Folder: {move_result.get('source_folder', 'N/A')}\n"
        result += f"Target Folder: {move_result.get('target_folder', 'N/A')}\n"
        result += f"Status: {move_result.get('status', 'N/A')}\n"
        result += f"Message: {move_result.get('message', 'N/A')}\n"

        return result
    except ValueError as e:
        # Handle validation errors with clear messages
        error_msg = f"Invalid input: {str(e)}"
        await ctx.error(error_msg)
        raise ValueError(error_msg)
    except RuntimeError as e:
        # Handle AppleScript errors with helpful context
        error_msg = str(e)
        await ctx.error(error_msg)
        raise RuntimeError(error_msg)
    except Exception as e:
        # Handle unexpected errors
        error_msg = f"Unexpected error moving note '{note_id}' to '{target_folder_path}': {str(e)}"
        await ctx.error(error_msg)
        raise


@mcp.tool()
async def search_notes(
    ctx: Context,
    keywords: str = Field(
        ..., description="Comma-separated keywords to search for in note content"
    ),
) -> str:
    """Search for notes containing the specified keywords.

    Features:
    - Searches through all notes in Apple Notes
    - Finds notes containing any of the specified keywords
    - Case-insensitive search in note content
    - Returns note details with matched keywords

    Output Format:
    - Numbered list of matching notes
    - Note names, IDs, and folder locations
    - Matched keywords for each note
    """
    try:
        # Parse keywords from comma-separated string
        keyword_list = [kw.strip() for kw in keywords.split(",") if kw.strip()]

        if not keyword_list:
            return (
                "No valid keywords provided. Please provide comma-separated keywords."
            )

        # Search for notes
        search_results = await notes_tools.search_notes(keyword_list)

        if not search_results:
            return f"No notes found containing the keywords: {', '.join(keyword_list)}"

        # Format the response
        result = f"Search Results ({len(search_results)} notes found):\n\n"
        result += f"Keywords searched: {', '.join(keyword_list)}\n\n"

        for i, note in enumerate(search_results, 1):
            result += f"{i:3d}. {note.get('name', 'N/A')}\n"
            result += f"     ID: {note.get('note_id', 'N/A')}\n"
            result += f"     Folder: {note.get('folder', 'N/A')}\n"
            result += f"     Matched: {note.get('matched_keyword', 'N/A')}\n"
            result += "\n"

        return result

    except ValueError as e:
        error_msg = f"Invalid search input: {str(e)}"
        await ctx.error(error_msg)
        raise ValueError(error_msg)
    except RuntimeError as e:
        error_msg = f"Search failed: {str(e)}"
        await ctx.error(error_msg)
        raise RuntimeError(error_msg)
    except Exception as e:
        await ctx.error(f"Error searching notes: {str(e)}")
        raise


def main():
    """Main entry point for the MCP server."""
    mcp.run()


# Run the server
if __name__ == "__main__":
    main()
