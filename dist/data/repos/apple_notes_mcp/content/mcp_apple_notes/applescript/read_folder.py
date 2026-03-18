from typing import Any

from .base_operations import BaseAppleScriptOperations


class ReadFolderOperations(BaseAppleScriptOperations):
    """Operations for reading Apple Notes folders by ID with name verification."""

    @staticmethod
    async def read_folder_by_id_and_name(
        folder_id: str, folder_name: str
    ) -> dict[str, Any]:
        """Read a folder by its primary key ID with AppleScript verification.

        This method relies on AppleScript's built-in verification:
        1. Validates input parameters (ID, folder name)
        2. AppleScript verifies ID and name match the same folder
        3. Performs the read operation if verification passes

        Args:
            folder_id: Primary key ID of the folder (e.g., "p2330")
            folder_name: Name of the folder to verify and read

        Returns:
            Folder data with metadata, child folders, and notes

        Raises:
            ValueError: If folder ID or name is empty or invalid
            RuntimeError: If folder not found, name doesn't match, or read fails
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

        sample_result = await ReadFolderOperations.execute_applescript(script_get_uuid)
        if sample_result.startswith("error:"):
            raise RuntimeError(
                f"Could not get store UUID for reading: {sample_result[6:]}"
            )

        # Extract store UUID from sample ID
        store_uuid = sample_result.split("//")[1].split("/")[0]
        full_folder_id = f"x-coredata://{store_uuid}/ICFolder/{folder_id}"

        # Read the folder using the full folder ID with name verification
        script = f"""
        tell application "Notes"
            try
                set primaryAccount to account "iCloud"
                set targetFolder to folder id "{full_folder_id}"
                
                set actualFolderName to name of targetFolder
                
                -- Verify that the actual folder name matches the provided name
                if actualFolderName is not equal to "{folder_name}" then
                    return "error:Folder name mismatch. Expected: {folder_name}, Actual: " & actualFolderName
                end if
                
                set folderId to id of targetFolder
                
                -- Get direct child folders
                set childFolders to {{}}
                set childFolderCount to count of folders of targetFolder
                
                repeat with i from 1 to childFolderCount
                    set childFolder to folder i of targetFolder
                    set childFolderInfo to {{name of childFolder, id of childFolder}}
                    set end of childFolders to childFolderInfo
                end repeat
                
                -- Get notes in the folder
                set folderNotes to {{}}
                set noteCount to count of notes of targetFolder
                
                repeat with i from 1 to noteCount
                    set currentNote to note i of targetFolder
                    set noteInfo to {{name of currentNote, id of currentNote}}
                    set end of folderNotes to noteInfo
                end repeat
                
                -- Build result string
                set resultString to "success:" & actualFolderName & "|||" & folderId & "|||Unknown|||Unknown|||" & childFolderCount & "|||" & noteCount
                
                -- Add child folders info
                repeat with childFolderInfo in childFolders
                    set resultString to resultString & "|||" & item 1 of childFolderInfo & "|||" & item 2 of childFolderInfo
                end repeat
                
                -- Add notes info
                repeat with noteInfo in folderNotes
                    set resultString to resultString & "|||" & item 1 of noteInfo & "|||" & item 2 of noteInfo
                end repeat
                
                return resultString
            on error errMsg
                return "error:iCloud account not available. Please enable iCloud Notes sync - " & errMsg
            end try
        end tell
        """

        result = await ReadFolderOperations.execute_applescript(script)

        if result.startswith("error:"):
            error_msg = result[6:]
            if "Folder name mismatch" in error_msg:
                raise RuntimeError(f"Folder name mismatch: {error_msg}")
            elif "Folder not found" in error_msg or "not found" in error_msg:
                raise RuntimeError(f"Folder with ID '{folder_id}' not found")
            else:
                raise RuntimeError(f"Failed to read folder: {error_msg}")

        return ReadFolderOperations._parse_read_by_id_and_name_result(result, folder_id)

    @staticmethod
    def _parse_read_by_id_and_name_result(
        result: str, folder_id: str
    ) -> dict[str, Any]:
        """Parse the AppleScript result for read by ID and name verification."""
        try:
            # The result format is: success:name|||full_id|||creation_date|||modification_date|||child_folder_count|||note_count|||child_folder1_name|||child_folder1_id|||...|||note1_name|||note1_id|||...
            if result.startswith("success:"):
                content = result[8:]  # Remove "success:" prefix
                parts = content.split("|||")

                if len(parts) >= 6:
                    from .note_id_utils import NoteIDUtils

                    # Extract basic folder info
                    folder_name = parts[0]
                    full_folder_id = parts[1]
                    creation_date = parts[2]
                    modification_date = parts[3]
                    child_folder_count = int(parts[4])
                    note_count = int(parts[5])

                    # Extract primary key from full ID
                    primary_key = NoteIDUtils.extract_primary_key(full_folder_id)

                    # Parse child folders
                    child_folders = []
                    child_folder_start = 6
                    child_folder_end = child_folder_start + (child_folder_count * 2)

                    for i in range(child_folder_start, child_folder_end, 2):
                        if i + 1 < len(parts):
                            child_folder_name = parts[i]
                            child_folder_full_id = parts[i + 1]
                            child_folder_primary_key = NoteIDUtils.extract_primary_key(
                                child_folder_full_id
                            )
                            child_folders.append(
                                {
                                    "name": child_folder_name,
                                    "id": child_folder_primary_key,
                                }
                            )

                    # Parse notes
                    notes = []
                    note_start = child_folder_end

                    for i in range(note_start, len(parts), 2):
                        if i + 1 < len(parts):
                            note_name = parts[i]
                            note_full_id = parts[i + 1]
                            note_primary_key = NoteIDUtils.extract_primary_key(
                                note_full_id
                            )
                            notes.append(
                                {"name": note_name, "note_id": note_primary_key}
                            )

                    return {
                        "name": folder_name,
                        "folder_id": primary_key,
                        "creation_date": creation_date,
                        "modification_date": modification_date,
                        "child_folders_count": child_folder_count,
                        "notes_count": note_count,
                        "child_folders": child_folders,
                        "notes": notes,
                        "status": "found",
                        "read_method": "by_id_and_name",
                    }
                else:
                    return {
                        "name": "Unknown",
                        "folder_id": folder_id,
                        "creation_date": "Unknown",
                        "modification_date": "Unknown",
                        "child_folders_count": 0,
                        "notes_count": 0,
                        "child_folders": [],
                        "notes": [],
                        "status": "found",
                        "read_method": "by_id_and_name",
                    }
            else:
                raise RuntimeError(f"Unexpected result format: {result}")
        except Exception as e:
            raise RuntimeError(f"Failed to parse read by ID and name result: {str(e)}")
