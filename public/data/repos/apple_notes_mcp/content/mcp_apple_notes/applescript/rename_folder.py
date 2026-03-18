from typing import Any

from .base_operations import BaseAppleScriptOperations
from .validation_utils import ValidationUtils


class RenameFolderOperations(BaseAppleScriptOperations):
    """Operations for renaming folders in Apple Notes."""

    @staticmethod
    def _create_applescript_quoted_string(text: str) -> str:
        """Escape text for safe AppleScript usage."""
        return ValidationUtils.create_applescript_quoted_string(text)

    @staticmethod
    async def _check_duplicate_name(
        new_name: str, folder_path: str, current_name: str
    ) -> None:
        """Check if the new folder name would create a duplicate."""
        # Get all folder names in the target location
        if not folder_path:
            # Root level - check root folders
            script = """
            tell application "Notes"
                try
                    set primaryAccount to account "iCloud"
                    set folderNames to {}
                    repeat with rootFolder in folders of primaryAccount
                        set end of folderNames to name of rootFolder
                    end repeat
                    return folderNames
                on error errMsg
                    return "error:iCloud account not available. Please enable iCloud Notes sync - " & errMsg
                end try
            end tell
            """
        else:
            # Nested level - navigate to parent and check subfolders
            path_components = ValidationUtils.parse_folder_path(folder_path)
            script = f"""
            tell application "Notes"
                try
                    set primaryAccount to account "iCloud"
                    set currentFolder to missing value
                    set pathComponents to {{{", ".join([f'"{component}"' for component in path_components])}}}
                    
                    repeat with i from 1 to count of pathComponents
                        set componentName to item i of pathComponents
                        
                        if currentFolder is missing value then
                            repeat with rootFolder in folders of primaryAccount
                                if name of rootFolder is componentName then
                                    set currentFolder to rootFolder
                                    exit repeat
                                end if
                            end repeat
                        else
                            repeat with subFolder in folders of currentFolder
                                if name of subFolder is componentName then
                                    set currentFolder to subFolder
                                    exit repeat
                                end if
                            end repeat
                        end if
                    end repeat
                    
                    if currentFolder is missing value then
                        return "error:Path not found"
                    end if
                    
                    set folderNames to {{}}
                    repeat with subFolder in folders of currentFolder
                        set end of folderNames to name of subFolder
                    end repeat
                    return folderNames
                on error errMsg
                    return "error:iCloud account not available. Please enable iCloud Notes sync - " & errMsg
                end try
            end tell
            """

        result = await RenameFolderOperations.execute_applescript(script)

        if result.startswith("error:"):
            raise RuntimeError(f"Failed to check duplicates: {result[6:]}")

        # Parse folder names and check for duplicates
        if result and result != "{}":
            # Split the result by comma and check each name
            folder_names = [
                name.strip().strip('"') for name in result.split(",") if name.strip()
            ]

            for folder_name in folder_names:
                # Check if the new name already exists (excluding the current folder being renamed)
                if folder_name == new_name:
                    # Only raise error if it's not the same folder being renamed
                    # OR if the new name is different from current name
                    if folder_name != current_name:
                        location = (
                            "root level"
                            if not folder_path
                            else f"folder '{folder_path}'"
                        )
                        raise ValueError(
                            f"A folder named '{new_name}' already exists in {location}"
                        )

    @staticmethod
    async def rename_folder_by_id(
        folder_id: str, current_name: str, new_name: str
    ) -> dict[str, Any]:
        """Rename a folder by its primary key ID with AppleScript verification.

        This method relies on AppleScript's built-in verification:
        1. Validates input parameters (ID, current name, new name)
        2. AppleScript verifies ID and name match the same folder
        3. Performs the rename operation if verification passes

        Args:
            folder_id: Primary key ID of the folder (e.g., "p2330")
            current_name: Current name of the folder to verify and rename
            new_name: New name for the folder

        Returns:
            Rename result with status and details

        Raises:
            ValueError: If folder ID, current name, or new name is empty or invalid
            RuntimeError: If folder not found, name doesn't match, or duplicate name exists
        """
        # Validate inputs
        if not folder_id or not folder_id.strip():
            raise ValueError("Folder ID cannot be empty or contain only whitespace")

        if not current_name or not current_name.strip():
            raise ValueError(
                "Current folder name cannot be empty or contain only whitespace"
            )

        if not new_name or not new_name.strip():
            raise ValueError(
                "New folder name cannot be empty or contain only whitespace"
            )

        folder_id = folder_id.strip()
        current_name = current_name.strip()
        new_name = new_name.strip()

        # Validate folder names
        try:
            validated_current_name = ValidationUtils.validate_folder_name(current_name)
            validated_new_name = ValidationUtils.validate_folder_name(new_name)

            # Check if new name is same as current name
            if validated_current_name == validated_new_name:
                raise ValueError("New folder name cannot be the same as current name")

        except ValueError as e:
            raise ValueError(f"Invalid input: {str(e)}")

        # Check for duplicate folder names in the same location
        # Since we don't have folder_path, we'll check for duplicates at root level
        await RenameFolderOperations._check_duplicate_name(
            validated_new_name, "", validated_current_name
        )

        # Build full Core Data ID from primary key using dynamic store UUID
        # First get a sample folder to extract the store UUID
        script_get_uuid = """
        tell application "Notes"
            try
                set primaryAccount to account "iCloud"
                set sampleFolder to folder 1 of primaryAccount
                set sampleId to id of sampleFolder as string
                return sampleId
            on error errMsg
                return "error:iCloud account not available. Please enable iCloud Notes sync - " & errMsg
            end try
        end tell
        """

        sample_result = await RenameFolderOperations.execute_applescript(
            script_get_uuid
        )
        if sample_result.startswith("error:"):
            raise RuntimeError(
                f"Could not get store UUID for renaming: {sample_result[6:]}"
            )

        # Extract store UUID from sample ID
        store_uuid = sample_result.split("//")[1].split("/")[0]
        full_folder_id = f"x-coredata://{store_uuid}/ICFolder/{folder_id}"

        # Escape strings for safe AppleScript usage
        escaped_current_name = RenameFolderOperations._create_applescript_quoted_string(
            validated_current_name
        )
        escaped_new_name = RenameFolderOperations._create_applescript_quoted_string(
            validated_new_name
        )

        # Get folder info and verify name matches, then rename
        script = f"""
        tell application "Notes"
            try
                set primaryAccount to account "iCloud"
                set targetFolder to folder id "{full_folder_id}"
                
                set actualFolderName to name of targetFolder as string
                
                -- Verify the folder name matches
                if actualFolderName is not {escaped_current_name} then
                    return "error:Folder name mismatch. Expected: {validated_current_name}, Found: " & actualFolderName
                end if
                
                -- Rename the folder
                set name of targetFolder to {escaped_new_name}
                
                return "success:" & actualFolderName & ", " & {escaped_new_name} & ", " & "root level"
            on error errMsg
                return "error:iCloud account not available. Please enable iCloud Notes sync - " & errMsg
            end try
        end tell
        """

        result = await RenameFolderOperations.execute_applescript(script)

        if result.startswith("error:"):
            error_msg = result[6:]
            if "Folder name mismatch" in error_msg:
                raise RuntimeError(error_msg)
            elif "Folder not found" in error_msg or "not found" in error_msg:
                raise RuntimeError(f"Folder with ID '{folder_id}' not found")
            else:
                raise RuntimeError(f"Failed to rename folder: {error_msg}")

        return RenameFolderOperations._parse_rename_by_id_result(
            result, folder_id, validated_current_name, validated_new_name
        )

    @staticmethod
    def _parse_rename_by_id_result(
        result: str, folder_id: str, current_name: str, new_name: str
    ) -> dict[str, str]:
        """Parse the AppleScript result for rename by ID and return rename information."""
        try:
            # The result format is: success:current_name, new_name, folder_path
            if result.startswith("success:"):
                parts = result[8:].split(", ")  # Remove "success:" prefix

                if len(parts) >= 3:
                    return {
                        "folder_id": folder_id,
                        "current_name": parts[0],
                        "new_name": parts[1],
                        "folder_path": parts[2],
                        "status": "renamed",
                        "rename_method": "by_id_and_name",
                    }
                else:
                    return {
                        "folder_id": folder_id,
                        "current_name": current_name,
                        "new_name": new_name,
                        "folder_path": "root level",
                        "status": "renamed",
                        "rename_method": "by_id_and_name",
                    }
            else:
                raise RuntimeError(f"Unexpected result format: {result}")
        except Exception as e:
            raise RuntimeError(f"Failed to parse rename by ID result: {str(e)}")
