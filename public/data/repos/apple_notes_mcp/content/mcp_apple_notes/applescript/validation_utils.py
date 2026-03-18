
from .base_operations import BaseAppleScriptOperations


class ValidationUtils(BaseAppleScriptOperations):
    """Centralized validation utilities for Apple Notes operations."""

    # Constants
    MAX_NESTING_DEPTH = 5
    MAX_NOTE_NAME_LENGTH = 250
    INVALID_CHARS = ["<", ">", ":", '"', "|", "?", "*"]

    @staticmethod
    def validate_folder_name(folder_name: str) -> str:
        """Validate and clean folder name.

        Args:
            folder_name: The folder name to validate

        Returns:
            Cleaned folder name

        Raises:
            ValueError: If folder name is invalid
        """
        if not folder_name or not folder_name.strip():
            raise ValueError("Folder name cannot be empty or contain only whitespace")

        # Clean the name
        folder_name = folder_name.strip()

        # Check length limit (same as create_folder uses)
        if len(folder_name) > 128:
            raise ValueError(
                f"Folder name exceeds maximum length of 128 characters (current: {len(folder_name)} characters)"
            )

        # Check for invalid characters
        for char in ValidationUtils.INVALID_CHARS:
            if char in folder_name:
                raise ValueError(f"Folder name contains invalid character '{char}'")

        return folder_name

    @staticmethod
    def validate_folder_path(folder_path: str) -> str:
        """Validate and clean folder path.

        Args:
            folder_path: The folder path to validate

        Returns:
            Cleaned folder path

        Raises:
            ValueError: If path is invalid
        """
        if not folder_path:
            return ""

        # Clean the path
        folder_path = folder_path.strip()

        # Remove leading/trailing slashes
        folder_path = folder_path.strip("/")

        # Check for invalid patterns
        if "//" in folder_path:
            raise ValueError("Folder path contains invalid double slashes")

        # Check for invalid characters
        for char in ValidationUtils.INVALID_CHARS:
            if char in folder_path:
                raise ValueError(f"Folder path contains invalid character '{char}'")

        return folder_path

    @staticmethod
    def validate_note_name(name: str) -> str:
        """Validate and clean note name.

        Args:
            name: The note name to validate

        Returns:
            Cleaned note name

        Raises:
            ValueError: If note name is invalid
        """
        if not name or not name.strip():
            raise ValueError("Note name cannot be empty or contain only whitespace")

        # Clean the name
        name = name.strip()

        # Handle backtick-escaped names
        if name.startswith("`") and name.endswith("`"):
            # Extract the name from backticks and skip validation
            escaped_name = name[1:-1]  # Remove backticks
            if not escaped_name:
                raise ValueError(
                    "Note name cannot be empty when using backtick escaping"
                )

            # Check for Apple Notes title length limit
            if len(escaped_name) > ValidationUtils.MAX_NOTE_NAME_LENGTH:
                raise ValueError(
                    f"Note name exceeds Apple Notes limit of {ValidationUtils.MAX_NOTE_NAME_LENGTH} characters (current: {len(escaped_name)} characters)"
                )

            # Return the escaped name without backticks
            return escaped_name.strip()

        # Check for Apple Notes title length limit
        if len(name) > ValidationUtils.MAX_NOTE_NAME_LENGTH:
            raise ValueError(
                f"Note name exceeds Apple Notes limit of {ValidationUtils.MAX_NOTE_NAME_LENGTH} characters (current: {len(name)} characters)"
            )

        # Check for invalid characters (only for non-escaped names)
        for char in ValidationUtils.INVALID_CHARS:
            if char in name:
                raise ValueError(
                    f"Note name contains invalid character '{char}'. Use backticks (`name`) to escape special characters."
                )

        return name

    @staticmethod
    def validate_note_body(body: str) -> str:
        """Validate note body.

        Args:
            body: The note body to validate

        Returns:
            Validated note body
        """
        if body is None:
            body = ""
        return str(body)

    @staticmethod
    def extract_title_from_html(html_content: str) -> str:
        """Extract title text from HTML content that may contain <h1> tags.

        Args:
            html_content: HTML content that may contain <h1>title</h1>

        Returns:
            Extracted title text (empty string if no title found)
        """
        if not html_content:
            return ""

        # Simple regex to extract content from <h1> tags
        import re

        # Look for <h1>content</h1> pattern (case insensitive)
        h1_match = re.search(
            r"<h1[^>]*>(.*?)</h1>", html_content, re.IGNORECASE | re.DOTALL
        )
        if h1_match:
            title_content = h1_match.group(1).strip()
            # Remove any other HTML tags from title content
            title_content = re.sub(r"<[^>]+>", "", title_content).strip()
            return title_content

        return ""

    @staticmethod
    def validate_html_title_content(html_content: str) -> str:
        """Validate that HTML content contains a non-empty title.

        Args:
            html_content: HTML content that should contain a valid <h1> title

        Returns:
            Original HTML content if valid

        Raises:
            ValueError: If title is empty, contains only whitespace, or is missing
        """
        if not html_content:
            raise ValueError("Content cannot be empty")

        # Extract title from HTML
        title_text = ValidationUtils.extract_title_from_html(html_content)

        # Check if title is empty or whitespace only
        if not title_text or not title_text.strip():
            raise ValueError(
                "Title cannot be empty or contain only whitespace. Please provide a valid title in <h1> tags."
            )

        return html_content

    @staticmethod
    def parse_folder_path(folder_path: str) -> list[str]:
        """Parse a folder path into components.

        Args:
            folder_path: The folder path to parse

        Returns:
            List of folder path components
        """
        if not folder_path:
            return ["Notes"]  # Default folder
        return [part.strip() for part in folder_path.split("/") if part.strip()]

    @staticmethod
    def validate_nesting_depth(
        folder_path: str, folder_name: str = None, operation: str = "create"
    ) -> None:
        """Validate that the operation won't exceed maximum nesting depth.

        Args:
            folder_path: The path where the operation will occur
            folder_name: The name of the folder (optional, for better error messages)
            operation: The operation being performed (create, move, etc.)

        Raises:
            ValueError: If the nesting depth would exceed the maximum allowed
        """
        # Count the depth of the existing path
        path_depth = 0
        if folder_path:
            path_components = folder_path.split("/")
            path_depth = len([comp for comp in path_components if comp.strip()])

        # Total depth will be path_depth + 1 (for the new/moved folder)
        total_depth = path_depth + 1

        if total_depth > ValidationUtils.MAX_NESTING_DEPTH:
            folder_info = f"'{folder_name}'" if folder_name else "folder"
            raise ValueError(
                f"Cannot {operation} {folder_info} in path '{folder_path}'. "
                f"This would create a nesting depth of {total_depth} levels, "
                f"which exceeds the maximum allowed depth of {ValidationUtils.MAX_NESTING_DEPTH} levels. "
                f"Please {operation} the folder at a higher level in the hierarchy."
            )

    @staticmethod
    def validate_move_operation(
        source_path: str, target_path: str, folder_name: str
    ) -> None:
        """Validate that a move operation is valid.

        Args:
            source_path: The current path of the folder
            target_path: The target path where the folder will be moved
            folder_name: The name of the folder being moved

        Raises:
            ValueError: If the move operation is invalid
        """
        # Check if source and target are the same
        if source_path == target_path:
            raise ValueError(f"Cannot move folder '{folder_name}' to the same location")

        # Validate nesting depth for target path
        ValidationUtils.validate_nesting_depth(target_path, folder_name, "move")

    @staticmethod
    async def check_path_exists(folder_path: str) -> bool:
        """Check if a folder path exists.

        Args:
            folder_path: The folder path to check

        Returns:
            True if path exists, False otherwise
        """
        try:
            # Handle root level (empty path) - root level always exists
            if not folder_path or folder_path.strip() == "":
                return True

            path_components = ValidationUtils.parse_folder_path(folder_path)
            if not path_components or not path_components[0]:
                return False

            script = f"""
            tell application "Notes"
                try
                    set primaryAccount to account "iCloud"
                    set currentFolder to missing value
                    set pathComponents to {{{", ".join([f'"{component}"' for component in path_components])}}}
                    
                    repeat with i from 1 to count of pathComponents
                        set componentName to item i of pathComponents
                        
                        if currentFolder is missing value then
                            -- Check root folders in iCloud account
                            set found to false
                            repeat with rootFolder in folders of primaryAccount
                                if name of rootFolder is componentName then
                                    set currentFolder to rootFolder
                                    set found to true
                                    exit repeat
                                end if
                            end repeat
                            if not found then
                                return "error:Folder not found"
                            end if
                        else
                            -- Check subfolders
                            set found to false
                            repeat with subFolder in folders of currentFolder
                                if name of subFolder is componentName then
                                    set currentFolder to subFolder
                                    set found to true
                                    exit repeat
                                end if
                            end repeat
                            if not found then
                                return "error:Folder not found"
                            end if
                        end if
                    end repeat
                    
                    return "exists"
                on error errMsg
                    return "error:iCloud account not available. Please enable iCloud Notes sync - " & errMsg
                end try
            end tell
            """

            result = await ValidationUtils.execute_applescript(script)
            return result == "exists"
        except Exception:
            return False

    @staticmethod
    async def check_folder_exists_at_root(folder_name: str) -> bool:
        """Check if a folder exists at root level.

        Args:
            folder_name: The folder name to check

        Returns:
            True if folder exists at root, False otherwise
        """
        try:
            script = f"""
            tell application "Notes"
                try
                    set primaryAccount to account "iCloud"
                    repeat with rootFolder in folders of primaryAccount
                        if name of rootFolder is "{folder_name}" then
                            return "exists"
                        end if
                    end repeat
                    return "not_found"
                on error errMsg
                    return "error:iCloud account not available. Please enable iCloud Notes sync - " & errMsg
                end try
            end tell
            """

            result = await ValidationUtils.execute_applescript(script)
            return result == "exists"
        except Exception:
            return False

    @staticmethod
    def create_applescript_quoted_string(text: str) -> str:
        """Create a properly quoted string for AppleScript.

        Args:
            text: The text to quote

        Returns:
            Quoted string for AppleScript
        """
        if not text:
            return '""'

        # Escape the text for AppleScript string literals
        escaped_text = text.replace("\\", "\\\\").replace('"', '\\"')
        return f'"{escaped_text}"'

    @staticmethod
    def truncate_note_name(name: str, max_length: int = None) -> str:
        """Intelligently truncate note name to fit Apple Notes limit.

        Args:
            name: The note name to truncate
            max_length: Maximum length (defaults to MAX_NOTE_NAME_LENGTH)

        Returns:
            Truncated note name
        """
        if max_length is None:
            max_length = ValidationUtils.MAX_NOTE_NAME_LENGTH

        if len(name) <= max_length:
            return name

        # Try to truncate at word boundaries
        truncated = name[: max_length - 3]  # Leave room for "..."

        # Find the last space to avoid cutting words
        last_space = truncated.rfind(" ")
        if last_space > max_length * 0.7:  # If we can find a space in the last 30%
            truncated = truncated[:last_space]

        return truncated + "..."
