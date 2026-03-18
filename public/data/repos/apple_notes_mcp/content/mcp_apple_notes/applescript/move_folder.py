from typing import Any

from .base_operations import BaseAppleScriptOperations
from .validation_utils import ValidationUtils


class MoveFolderOperations(BaseAppleScriptOperations):
    """Operations for moving folders in Apple Notes."""

    @staticmethod
    def _validate_folder_name(folder_name: str) -> str:
        """Validate and clean folder name."""
        return ValidationUtils.validate_folder_name(folder_name)

    @staticmethod
    def _validate_folder_path(folder_path: str) -> str:
        """Validate and clean folder path."""
        return ValidationUtils.validate_folder_path(folder_path)

    @staticmethod
    def _validate_nesting_depth(
        source_path: str, target_path: str, folder_name: str
    ) -> None:
        """Validate that the move operation won't exceed maximum nesting depth.

        Args:
            source_path: The current path of the folder
            target_path: The target path where the folder will be moved
            folder_name: The name of the folder being moved

        Raises:
            ValueError: If the nesting depth would exceed the maximum allowed
        """
        ValidationUtils.validate_move_operation(source_path, target_path, folder_name)

    @staticmethod
    async def _check_folder_exists_at_root(folder_name: str) -> bool:
        """Check if a folder exists at root level."""
        return await ValidationUtils.check_folder_exists_at_root(folder_name)

    @staticmethod
    async def _check_path_exists(folder_path: str) -> bool:
        """Check if a folder path exists."""
        return await ValidationUtils.check_path_exists(folder_path)

    @staticmethod
    def _create_applescript_quoted_string(text: str) -> str:
        """Escape text for safe AppleScript usage."""
        return ValidationUtils.create_applescript_quoted_string(text)

    @staticmethod
    async def move_folder(
        source_path: str, folder_name: str, target_path: str = ""
    ) -> dict[str, Any]:
        """Move a folder from one location to another in Apple Notes.

        Args:
            source_path: The current path of the folder to move
            folder_name: The name of the folder to move
            target_path: The target path where to move the folder. If empty, moves to root level.
        """
        # Validate inputs
        folder_name = MoveFolderOperations._validate_folder_name(folder_name)
        source_path = MoveFolderOperations._validate_folder_path(source_path)
        target_path = MoveFolderOperations._validate_folder_path(target_path)

        # Validate nesting depth
        MoveFolderOperations._validate_nesting_depth(
            source_path, target_path, folder_name
        )

        # Check if source path exists
        source_exists = await MoveFolderOperations._check_path_exists(source_path)
        if not source_exists:
            raise RuntimeError(f"Source path '{source_path}' does not exist")

        # Check if target path exists (if provided)
        if target_path:
            target_exists = await MoveFolderOperations._check_path_exists(target_path)
            if not target_exists:
                raise RuntimeError(f"Target path '{target_path}' does not exist")

        # Check if folder exists in source path
        if source_path:
            # Source is a nested path
            source_folder_exists = await MoveFolderOperations._check_path_exists(
                f"{source_path}/{folder_name}"
            )
            if not source_folder_exists:
                raise RuntimeError(
                    f"Folder '{folder_name}' not found in path '{source_path}'"
                )
        else:
            # Source is root level
            source_folder_exists = (
                await MoveFolderOperations._check_folder_exists_at_root(folder_name)
            )
            if not source_folder_exists:
                raise RuntimeError(f"Folder '{folder_name}' not found at root level")

        # Check if target already has a folder with the same name
        if target_path:
            target_folder_exists = await MoveFolderOperations._check_path_exists(
                f"{target_path}/{folder_name}"
            )
            if target_folder_exists:
                raise RuntimeError(
                    f"Target path '{target_path}' already contains a folder named '{folder_name}'"
                )
        else:
            # Check root level for duplicate
            root_folder_exists = (
                await MoveFolderOperations._check_folder_exists_at_root(folder_name)
            )
            if root_folder_exists:
                raise RuntimeError(
                    f"Root level already contains a folder named '{folder_name}'"
                )

        # Perform the move operation
        if not target_path:
            # Move to root level
            return await MoveFolderOperations._move_to_root(source_path, folder_name)
        else:
            # Move to target path
            return await MoveFolderOperations._move_to_path(
                source_path, folder_name, target_path
            )

    @staticmethod
    async def _move_to_root(source_path: str, folder_name: str) -> dict[str, Any]:
        """Move a folder to root level."""
        # Escape strings for safe AppleScript usage
        escaped_folder_name = MoveFolderOperations._create_applescript_quoted_string(
            folder_name
        )
        escaped_source_path = MoveFolderOperations._create_applescript_quoted_string(
            source_path
        )

        script = f"""
        tell application "Notes"
            try
                set primaryAccount to account "iCloud"
                
                if {escaped_source_path} is "" then
                    return "error:Cannot move folder to same location"
                end if
                
                -- Find source folder and move to root
                set sourceFolder to missing value
                
                -- Navigate to source path
                set pathComponents to words of {escaped_source_path} delimited by "/"
                set currentFolder to missing value
                
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
                    return "error:Cannot navigate to source folder"
                end if
                
                -- Move the folder to root
                move folder {escaped_folder_name} of currentFolder to beginning of folders
                return "success:" & {escaped_folder_name} & ", " & {escaped_source_path} & ", root"
                
            on error errMsg
                return "error:iCloud account not available. Please enable iCloud Notes sync - " & errMsg
            end try
        end tell
        """

        result = await MoveFolderOperations.execute_applescript(script)

        if result.startswith("error:"):
            raise RuntimeError(f"Failed to move folder: {result[6:]}")

        return MoveFolderOperations._parse_move_result(result)

    @staticmethod
    async def _move_to_path(
        source_path: str, folder_name: str, target_path: str
    ) -> dict[str, Any]:
        """Move a folder to a specific target path."""
        # Escape strings for safe AppleScript usage
        escaped_folder_name = MoveFolderOperations._create_applescript_quoted_string(
            folder_name
        )
        escaped_source_path = MoveFolderOperations._create_applescript_quoted_string(
            source_path
        )
        escaped_target_path = MoveFolderOperations._create_applescript_quoted_string(
            target_path
        )

        script = f"""
        tell application "Notes"
            try
                set primaryAccount to account "iCloud"
                
                -- Find source folder
                set sourceFolder to missing value
                if {escaped_source_path} is "" then
                    -- Source is root level
                    repeat with rootFolder in folders of primaryAccount
                        if name of rootFolder is {escaped_folder_name} then
                            set sourceFolder to rootFolder
                            exit repeat
                        end if
                    end repeat
                else
                    -- Source is nested path
                    set pathComponents to words of {escaped_source_path} delimited by "/"
                    set currentFolder to missing value
                    
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
                    
                    if currentFolder is not missing value then
                        repeat with subFolder in folders of currentFolder
                            if name of subFolder is {escaped_folder_name} then
                                set sourceFolder to subFolder
                                exit repeat
                            end if
                        end repeat
                    end if
                end if
                
                if sourceFolder is missing value then
                    return "error:Source folder not found"
                end if
                
                -- Find target folder
                set targetFolder to missing value
                set targetPathComponents to words of {escaped_target_path} delimited by "/"
                set currentFolder to missing value
                
                repeat with i from 1 to count of targetPathComponents
                    set componentName to item i of targetPathComponents
                    
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
                    return "error:Cannot navigate to target folder"
                end if
                
                -- Move the folder
                move sourceFolder to beginning of folders of currentFolder
                return "success:" & {escaped_folder_name} & ", " & {escaped_source_path} & ", " & {escaped_target_path}
                
            on error errMsg
                return "error:iCloud account not available. Please enable iCloud Notes sync - " & errMsg
            end try
        end tell
        """

        result = await MoveFolderOperations.execute_applescript(script)

        if result.startswith("error:"):
            raise RuntimeError(f"Failed to move folder: {result[6:]}")

        return MoveFolderOperations._parse_move_result(result)

    @staticmethod
    def _parse_move_result(result: str) -> dict[str, str]:
        """Parse the AppleScript result and return move information."""
        try:
            # The result format is: success:folder_name, source_path, target_path
            if result.startswith("success:"):
                parts = result[8:].split(", ")  # Remove "success:" prefix

                if len(parts) >= 3:
                    return {
                        "folder_name": parts[0],
                        "source_path": parts[1],
                        "target_path": parts[2],
                        "status": "moved",
                        "message": f"Successfully moved folder '{parts[0]}' from '{parts[1]}' to '{parts[2]}'",
                    }
                else:
                    raise RuntimeError("Invalid move result format")
            else:
                raise RuntimeError(f"Unexpected result format: {result}")
        except Exception as e:
            raise RuntimeError(f"Failed to parse move result: {str(e)}")

    @staticmethod
    async def move_folder_by_id(
        folder_id: str, folder_name: str, target_path: str = ""
    ) -> dict[str, Any]:
        """Move a folder by its primary key ID with AppleScript verification.

        This method relies on AppleScript's built-in verification:
        1. Validates input parameters (ID, folder name, target path)
        2. AppleScript verifies ID and name match the same folder
        3. Performs the move operation if verification passes

        Args:
            folder_id: Primary key ID of the folder (e.g., "p2330")
            folder_name: Name of the folder to verify and move
            target_path: The target path where to move the folder. If empty, moves to root level.

        Returns:
            Move result with status and details

        Raises:
            ValueError: If folder ID or name is empty or invalid
            RuntimeError: If target path doesn't exist, folder not found, or name doesn't match
        """
        # Validate inputs
        if not folder_id or not folder_id.strip():
            raise ValueError("Folder ID cannot be empty or contain only whitespace")

        if not folder_name or not folder_name.strip():
            raise ValueError("Folder name cannot be empty or contain only whitespace")

        folder_id = folder_id.strip()
        folder_name = folder_name.strip()

        # Validate target path
        target_path = MoveFolderOperations._validate_folder_path(target_path)

        # Check if target path exists (if provided)
        if target_path:
            target_exists = await MoveFolderOperations._check_path_exists(target_path)
            if not target_exists:
                raise RuntimeError(f"Target path '{target_path}' does not exist")

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

        sample_result = await MoveFolderOperations.execute_applescript(script_get_uuid)
        if sample_result.startswith("error:"):
            raise RuntimeError(
                f"Could not get store UUID for moving: {sample_result[6:]}"
            )

        # Extract store UUID from sample ID
        store_uuid = sample_result.split("//")[1].split("/")[0]
        full_folder_id = f"x-coredata://{store_uuid}/ICFolder/{folder_id}"

        # Escape strings for safe AppleScript usage
        escaped_folder_name = MoveFolderOperations._create_applescript_quoted_string(
            folder_name
        )
        escaped_target_path = MoveFolderOperations._create_applescript_quoted_string(
            target_path
        )

        # Get folder info and verify name matches, then move
        script = f"""
        tell application "Notes"
            try
                set primaryAccount to account "iCloud"
                set targetFolder to folder id "{full_folder_id}"
                
                set actualFolderName to name of targetFolder as string
                
                -- Verify the folder name matches
                if actualFolderName is not {escaped_folder_name} then
                    return "error:Folder name mismatch. Expected: {folder_name}, Found: " & actualFolderName
                end if
                
                -- Move the folder
                if {escaped_target_path} is "" then
                    -- Move to root level
                    move targetFolder to beginning of folders
                    return "success:" & actualFolderName & ", root level"
                else
                    -- Move to target path (simplified approach)
                    set targetFolderPath to {escaped_target_path}
                    
                    -- Try to find target folder at root level first
                    set foundTarget to false
                    repeat with rootFolder in folders of primaryAccount
                        if name of rootFolder is targetFolderPath then
                            move targetFolder to beginning of folders of rootFolder
                            set foundTarget to true
                            exit repeat
                        end if
                    end repeat
                    
                    if not foundTarget then
                        return "error:Target folder not found: " & targetFolderPath
                    end if
                    
                    return "success:" & actualFolderName & ", " & {escaped_target_path}
                end if
                
            on error errMsg
                return "error:iCloud account not available. Please enable iCloud Notes sync - " & errMsg
            end try
        end tell
        """

        result = await MoveFolderOperations.execute_applescript(script)

        if result.startswith("error:"):
            error_msg = result[6:]
            if "Folder name mismatch" in error_msg:
                raise RuntimeError(error_msg)
            elif "Folder not found" in error_msg or "not found" in error_msg:
                raise RuntimeError(f"Folder with ID '{folder_id}' not found")
            else:
                raise RuntimeError(f"Failed to move folder: {error_msg}")

        return MoveFolderOperations._parse_move_by_id_result(
            result, folder_id, folder_name, target_path
        )

    @staticmethod
    def _parse_move_by_id_result(
        result: str, folder_id: str, folder_name: str, target_path: str
    ) -> dict[str, str]:
        """Parse the AppleScript result for move by ID and return move information."""
        try:
            # The result format is: success:folder_name, target_path
            if result.startswith("success:"):
                parts = result[8:].split(", ")  # Remove "success:" prefix

                if len(parts) >= 2:
                    return {
                        "folder_id": folder_id,
                        "folder_name": parts[0],
                        "target_path": parts[1],
                        "status": "moved",
                        "move_method": "by_id_and_name",
                    }
                else:
                    return {
                        "folder_id": folder_id,
                        "folder_name": folder_name,
                        "target_path": target_path if target_path else "root level",
                        "status": "moved",
                        "move_method": "by_id_and_name",
                    }
            else:
                raise RuntimeError(f"Unexpected result format: {result}")
        except Exception as e:
            raise RuntimeError(f"Failed to parse move by ID result: {str(e)}")
