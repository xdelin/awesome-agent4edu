
from .base_operations import BaseAppleScriptOperations
from .note_id_utils import NoteIDUtils
from .validation_utils import ValidationUtils


class CreateNoteOperations(BaseAppleScriptOperations):
    """Operations for creating Apple Notes."""

    @staticmethod
    async def create_note(
        name: str, body: str, folder_path: str = "Notes"
    ) -> dict[str, str]:
        """Create a new note with specified name, body, and folder path.

        This unified method handles both simple folders and nested paths.
        The folder path must exist before creating the note.

        IMPORTANT: The body parameter should contain complete HTML-formatted content
        that combines both title and body content. This will be passed directly to Apple Notes.

        Args:
            name: Name of the note for Apple Notes internal reference (used for identification)
            body: Complete HTML content combining title and body (e.g., '<h1>Title</h1><p>Content</p>')
            folder_path: Folder path (e.g., "Work" or "Work/Projects/2024").
                        Must exist before creating note. Defaults to "Notes".
        """
        # Validate and clean inputs
        try:
            validated_name = ValidationUtils.validate_note_name(name)
        except ValueError as e:
            # If name is too long, truncate it intelligently
            if "exceeds Apple Notes limit" in str(e):
                validated_name = ValidationUtils.truncate_note_name(name)
                validated_name = ValidationUtils.validate_note_name(
                    validated_name
                )  # Validate the truncated name
            else:
                raise e

        body = ValidationUtils.validate_note_body(body)
        folder_path = folder_path.strip()

        # Use the body content directly - it contains combined title + body HTML
        html_content = body

        # Check if it's a simple folder (no slashes) or nested path
        if "/" not in folder_path:
            # Simple folder - use direct folder access
            return await CreateNoteOperations._create_note_in_simple_folder(
                validated_name, html_content, folder_path
            )
        else:
            # Nested path - validate path exists and create note
            return await CreateNoteOperations._create_note_in_nested_path(
                validated_name, html_content, folder_path
            )

    @staticmethod
    async def _create_note_in_simple_folder(
        name: str, html_content: str, folder_name: str
    ) -> dict[str, str]:
        """Create a new note in a simple folder (no nested paths)."""
        # Escape the HTML content for AppleScript
        escaped_html = ValidationUtils.create_applescript_quoted_string(html_content)
        escaped_folder = ValidationUtils.create_applescript_quoted_string(folder_name)

        script = f"""
        tell application "Notes"
            try
                set targetFolder to folder {escaped_folder}
                set newNote to make new note at targetFolder with properties {{body:{escaped_html}}}
                return {{name:(name of newNote), folder:{escaped_folder}, note_id:(id of newNote as string)}}
            on error errMsg
                return "error:" & errMsg
            end try
        end tell
        """
        result = await CreateNoteOperations.execute_applescript(script)

        # Check if there was an error
        if result.startswith("error:"):
            raise RuntimeError(f"Failed to create note: {result[6:]}")

        return CreateNoteOperations._parse_note_result(result, folder_name)

    @staticmethod
    async def _create_note_in_nested_path(
        name: str, html_content: str, folder_path: str
    ) -> dict[str, str]:
        """Create a new note in a nested folder path. Path must exist."""
        # Validate that the folder path exists
        path_exists = await ValidationUtils.check_path_exists(folder_path)
        if not path_exists:
            raise RuntimeError(
                f"Folder path '{folder_path}' does not exist. Please create the folder structure first."
            )

        # Parse the path components
        path_components = ValidationUtils.parse_folder_path(folder_path)

        # Escape the HTML content for AppleScript
        escaped_html = ValidationUtils.create_applescript_quoted_string(html_content)
        escaped_folder_path = ValidationUtils.create_applescript_quoted_string(
            folder_path
        )

        script = f"""
        tell application "Notes"
            try
                set primaryAccount to account "iCloud"
                set currentFolder to missing value
                set pathComponents to {{{", ".join([ValidationUtils.create_applescript_quoted_string(component) for component in path_components])}}}
                
                repeat with i from 1 to count of pathComponents
                    set componentName to item i of pathComponents
                    
                    if currentFolder is missing value then
                        -- Start from root folders in iCloud account
                        repeat with rootFolder in folders of primaryAccount
                            if name of rootFolder is componentName then
                                set currentFolder to rootFolder
                                exit repeat
                            end if
                        end repeat
                    else
                        -- Navigate into subfolders
                        repeat with subFolder in folders of currentFolder
                            if name of subFolder is componentName then
                                set currentFolder to subFolder
                                exit repeat
                            end if
                        end repeat
                    end if
                end repeat
                
                if currentFolder is missing value then
                    return "error:Target folder not found"
                end if
                
                set newNote to make new note at currentFolder with properties {{body:{escaped_html}}}
                return {{name:(name of newNote), folder:{escaped_folder_path}, note_id:(id of newNote as string)}}
            on error errMsg
                return "error:iCloud account not available. Please enable iCloud Notes sync - " & errMsg
            end try
        end tell
        """

        result = await CreateNoteOperations.execute_applescript(script)

        # Check if there was an error
        if result.startswith("error:"):
            raise RuntimeError(f"Failed to create note: {result[6:]}")

        return CreateNoteOperations._parse_note_result(result, folder_path)

    @staticmethod
    def _parse_note_result(result: str, folder_path: str) -> dict[str, str]:
        """Parse the AppleScript result and return note information."""
        try:
            # Extract the note information from the result
            name_start = result.find("name:") + 5
            name_end = result.find(", folder:", name_start)
            name = result[name_start:name_end].strip()

            # Remove quotes from the name if present
            if name.startswith('"') and name.endswith('"'):
                name = name[1:-1]

            folder_start = result.find("folder:") + 7
            folder_end = result.find(", note_id:", folder_start)
            folder = result[folder_start:folder_end].strip()

            note_id_start = result.find("note_id:") + 8
            full_note_id = result[note_id_start:].strip().rstrip(",")

            # Extract just the primary key from the full note ID
            primary_key = NoteIDUtils.extract_primary_key(full_note_id)

            return {"name": name, "folder": folder, "note_id": primary_key}
        except Exception as e:
            raise RuntimeError(f"Failed to parse created note result: {str(e)}")
