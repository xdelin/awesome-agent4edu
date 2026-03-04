from typing import Any

from .base_operations import BaseAppleScriptOperations


class DeleteFolderOperations(BaseAppleScriptOperations):
    """Operations for deleting folders in Apple Notes by ID with name verification."""

    @staticmethod
    async def delete_folder_by_id_and_name(
        folder_id: str, folder_name: str
    ) -> dict[str, Any]:
        """Delete a folder by its primary key ID with AppleScript verification.

        This method relies on AppleScript's built-in verification:
        1. Validates input parameters (ID, folder name)
        2. AppleScript verifies ID and name match the same folder
        3. Performs the deletion if verification passes

        Args:
            folder_id: Primary key ID of the folder (e.g., "p2330")
            folder_name: Name of the folder to verify and delete

        Returns:
            Deletion result with status and details

        Raises:
            ValueError: If folder ID or name is empty or invalid
            RuntimeError: If folder not found, name doesn't match, or deletion fails
        """
        # Validate inputs
        if not folder_id or not folder_id.strip():
            raise ValueError("Folder ID cannot be empty or contain only whitespace")

        if not folder_name or not folder_name.strip():
            raise ValueError("Folder name cannot be empty or contain only whitespace")

        folder_id = folder_id.strip()
        folder_name = folder_name.strip()

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

        sample_result = await DeleteFolderOperations.execute_applescript(
            script_get_uuid
        )
        if sample_result.startswith("error:"):
            raise RuntimeError(
                f"Could not get store UUID for deletion: {sample_result[6:]}"
            )

        # Extract store UUID from sample ID
        store_uuid = sample_result.split("//")[1].split("/")[0]
        full_folder_id = f"x-coredata://{store_uuid}/ICFolder/{folder_id}"

        # Get folder info and verify name matches, then delete
        script = f"""
        tell application "Notes"
            try
                set primaryAccount to account "iCloud"
                set targetFolder to folder id "{full_folder_id}"
                
                set actualFolderName to name of targetFolder as string
                set folderId to id of targetFolder as string
                
                -- Verify the folder name matches
                if actualFolderName is not "{folder_name}" then
                    return "error:Folder name mismatch. Expected: {folder_name}, Found: " & actualFolderName
                end if
                
                delete targetFolder
                
                return "success:" & actualFolderName & ", " & folderId
            on error errMsg
                return "error:iCloud account not available. Please enable iCloud Notes sync - " & errMsg
            end try
        end tell
        """

        result = await DeleteFolderOperations.execute_applescript(script)

        if result.startswith("error:"):
            error_msg = result[6:]
            if "Folder name mismatch" in error_msg:
                raise RuntimeError(error_msg)
            elif "Folder not found" in error_msg or "not found" in error_msg:
                raise RuntimeError(f"Folder with ID '{folder_id}' not found")
            else:
                raise RuntimeError(f"Failed to delete folder: {error_msg}")

        return DeleteFolderOperations._parse_delete_by_id_result(
            result, folder_id, folder_name
        )

    @staticmethod
    def _parse_delete_by_id_result(
        result: str, folder_id: str, folder_name: str
    ) -> dict[str, str]:
        """Parse the AppleScript result for delete by ID and return deletion information."""
        try:
            # The result format is: success:name, full_id
            if result.startswith("success:"):
                parts = result[8:].split(", ")  # Remove "success:" prefix

                if len(parts) >= 2:
                    from .note_id_utils import NoteIDUtils

                    # Extract primary key from full ID
                    primary_key = NoteIDUtils.extract_folder_primary_key(parts[1])

                    return {
                        "name": parts[0],
                        "folder_id": primary_key,
                        "folder_path": "root level",
                        "status": "deleted",
                        "deletion_method": "by_id_and_name",
                    }
                else:
                    return {
                        "name": folder_name,
                        "folder_id": folder_id,
                        "folder_path": "root level",
                        "status": "deleted",
                        "deletion_method": "by_id_and_name",
                    }
            else:
                raise RuntimeError(f"Unexpected result format: {result}")
        except Exception as e:
            raise RuntimeError(f"Failed to parse delete by ID result: {str(e)}")
