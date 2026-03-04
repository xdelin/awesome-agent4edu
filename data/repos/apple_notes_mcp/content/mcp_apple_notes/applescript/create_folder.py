from typing import Any

from .base_operations import BaseAppleScriptOperations
from .note_id_utils import NoteIDUtils
from .validation_utils import ValidationUtils


class CreateFolderOperations(BaseAppleScriptOperations):
    """Operations for creating Apple Notes folders."""

    # Maximum nesting depth allowed by Apple Notes
    MAX_NESTING_DEPTH = 5

    @staticmethod
    async def create_folder(folder_name: str, folder_path: str = "") -> dict[str, Any]:
        """Create a folder in Apple Notes.

        Args:
            folder_name: Name of the folder to create
            folder_path: Optional path where to create the folder (e.g., "Parent/Child")
        """
        # Validate and clean folder name
        folder_name = ValidationUtils.validate_folder_name(folder_name)

        if not folder_path.strip():
            # Create at root level - validate nesting depth
            ValidationUtils.validate_nesting_depth("", folder_name, "create")
            return await CreateFolderOperations._create_root_folder(folder_name)
        else:
            # Validate and clean folder path
            try:
                folder_path = ValidationUtils.validate_folder_path(folder_path)
                if not folder_path:
                    # Empty path after validation, create at root level
                    ValidationUtils.validate_nesting_depth("", folder_name, "create")
                    return await CreateFolderOperations._create_root_folder(folder_name)

                # Validate nesting depth
                ValidationUtils.validate_nesting_depth(
                    folder_path, folder_name, "create"
                )

                # Create in specified path - AppleScript will handle path existence and duplicate checks
                return await CreateFolderOperations._create_nested_folder(
                    folder_name, folder_path
                )
            except ValueError:
                # Re-raise validation errors
                raise
            except RuntimeError as e:
                # Provide more helpful error messages
                error_msg = str(e)
                if "not found" in error_msg.lower():
                    raise RuntimeError(
                        f"Invalid folder path '{folder_path}'. The specified path does not exist. Please check the path and try again."
                    )
                elif "permission" in error_msg.lower():
                    raise RuntimeError(
                        f"Permission denied when creating folder '{folder_name}' in path '{folder_path}'. Please check your Apple Notes permissions."
                    )
                else:
                    raise RuntimeError(
                        f"Failed to create folder '{folder_name}' in path '{folder_path}': {error_msg}"
                    )

    @staticmethod
    async def _create_root_folder(folder_name: str) -> dict[str, Any]:
        """Create folder at root level."""
        escaped_name = ValidationUtils.create_applescript_quoted_string(folder_name)

        script = f"""
        tell application "Notes"
            try
                set primaryAccount to account "iCloud"
                set newFolder to make new folder with properties {{name:{escaped_name}}} at primaryAccount
                set folderName to name of newFolder
                set folderId to id of newFolder as string
                return "name:" & folderName & ", id:" & folderId
            on error errMsg number errNum
                return "ERROR:iCloud account not available. Please enable iCloud Notes sync - " & errMsg
            end try
        end tell
        """

        result = await CreateFolderOperations.execute_applescript(script)

        if result.startswith("ERROR:"):
            raise RuntimeError(f"Failed to create folder: {result[6:]}")

        try:
            # Parse the result - new format: "name:FolderName, id:fullId, path:optionalPath"
            if result.startswith("name:"):
                # Extract name
                name_start = result.find("name:") + 5
                name_end = result.find(", id:", name_start)
                if name_end == -1:
                    name_end = result.find(", path:", name_start)
                if name_end == -1:
                    name_end = len(result)
                name = result[name_start:name_end].strip()

                # Extract ID and convert to primary key
                folder_id = ""
                id_start = result.find(", id:") + 5
                if id_start > 4:
                    id_end = result.find(", path:", id_start)
                    if id_end == -1:
                        id_end = len(result)
                    full_id = result[id_start:id_end].strip()
                    # Extract primary key from full Core Data ID
                    folder_id = NoteIDUtils.extract_folder_primary_key(full_id)

                # Extract path if present
                path = ""
                path_start = result.find(", path:") + 7
                if path_start > 6:
                    path = result[path_start:].strip()

                return {"name": name, "id": folder_id}
            else:
                # Fallback to old format parsing
                name_start = result.find("name:") + 5
                name_end = result.find("}", name_start)
                if name_end == -1:
                    name_end = len(result)
                name = result[name_start:name_end].strip()

                return {"name": name, "id": ""}  # No ID for old format fallback
        except Exception as e:
            raise RuntimeError(f"Failed to parse folder creation result: {str(e)}")

    @staticmethod
    async def _create_nested_folder(
        folder_name: str, folder_path: str
    ) -> dict[str, Any]:
        """Create folder in nested path."""
        path_components = [
            comp.strip() for comp in folder_path.split("/") if comp.strip()
        ]

        escaped_name = ValidationUtils.create_applescript_quoted_string(folder_name)

        # Build path navigation script
        path_navigation = []
        for i, component in enumerate(path_components):
            escaped_component = ValidationUtils.create_applescript_quoted_string(
                component
            )
            if i == 0:
                path_navigation.append(
                    f"set targetFolder to folder {escaped_component}"
                )
            else:
                path_navigation.append(
                    f"set targetFolder to folder {escaped_component} of targetFolder"
                )

        navigation_script = "\n                ".join(path_navigation)

        script = f"""
        tell application "Notes"
            try
                set primaryAccount to account "iCloud"
                {navigation_script}
                set newFolder to make new folder with properties {{name:{escaped_name}}} at targetFolder
                set folderName to name of newFolder
                set folderId to id of newFolder as string
                return "name:" & folderName & ", id:" & folderId & ", path:{folder_path}/{folder_name}"
            on error errMsg number errNum
                if errNum = -1728 then -- Folder not found
                    return "ERROR:Path not found: {folder_path}"
                else
                    return "ERROR:iCloud account not available. Please enable iCloud Notes sync - " & errMsg
                end if
            end try
        end tell
        """

        result = await CreateFolderOperations.execute_applescript(script)

        if result.startswith("ERROR:"):
            raise RuntimeError(f"Failed to create nested folder: {result[6:]}")

        try:
            # Parse the result - new format: "name:FolderName, id:fullId, path:optionalPath"
            if result.startswith("name:"):
                # Extract name
                name_start = result.find("name:") + 5
                name_end = result.find(", id:", name_start)
                if name_end == -1:
                    name_end = result.find(", path:", name_start)
                if name_end == -1:
                    name_end = len(result)
                name = result[name_start:name_end].strip()

                # Extract ID and convert to primary key
                folder_id = ""
                id_start = result.find(", id:") + 5
                if id_start > 4:
                    id_end = result.find(", path:", id_start)
                    if id_end == -1:
                        id_end = len(result)
                    full_id = result[id_start:id_end].strip()
                    # Extract primary key from full Core Data ID
                    folder_id = NoteIDUtils.extract_folder_primary_key(full_id)

                # Extract path if present
                path = ""
                path_start = result.find(", path:") + 7
                if path_start > 6:
                    path = result[path_start:].strip()

                return {"name": name, "id": folder_id}
            else:
                # Fallback to old format parsing
                name_start = result.find("name:") + 5
                name_end = result.find("}", name_start)
                if name_end == -1:
                    name_end = len(result)
                name = result[name_start:name_end].strip()

                return {"name": name, "id": ""}  # No ID for old format fallback
        except Exception as e:
            raise RuntimeError(f"Failed to parse nested folder result: {str(e)}")

    # Backward compatibility methods
    @staticmethod
    async def create_folder_in_existing_path(folder_path: str) -> dict[str, Any]:
        """Create a folder in an existing path (does not create parent folders)."""
        path_components = ValidationUtils.parse_folder_path(folder_path)
        if len(path_components) > 1:
            folder_name = path_components[-1]
            parent_path = "/".join(path_components[:-1])
            return await CreateFolderOperations._create_nested_folder(
                folder_name, parent_path
            )
        else:
            return await CreateFolderOperations._create_root_folder(path_components[0])
