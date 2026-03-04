
from .base_operations import BaseAppleScriptOperations
from .validation_utils import ValidationUtils


class NoteIDUtils(BaseAppleScriptOperations):
    """Utilities for note ID operations and duplicate handling."""

    @staticmethod
    async def get_all_notes_in_folder(
        folder_path: str = "Notes",
    ) -> list[dict[str, str]]:
        """Get all notes in a folder with their IDs and names.

        Args:
            folder_path: Folder path to search in

        Returns:
            List of dictionaries with note_id, name, and folder info

        Raises:
            RuntimeError: If folder path doesn't exist
        """
        folder_path = folder_path.strip()

        if "/" not in folder_path:
            return await NoteIDUtils._get_all_notes_in_simple_folder(folder_path)
        else:
            return await NoteIDUtils._get_all_notes_in_nested_path(folder_path)

    @staticmethod
    async def _get_all_notes_in_simple_folder(folder_name: str) -> list[dict[str, str]]:
        """Get all notes in a simple folder."""
        script = f"""
        tell application "Notes"
            try
                set targetFolder to folder "{folder_name}"
                set noteList to {{}}
                
                repeat with theNote in notes of targetFolder
                    set noteInfo to {{id:(id of theNote as string), name:(name of theNote as string), folder:"{folder_name}"}}
                    copy noteInfo to end of noteList
                end repeat
                
                return noteList
            on error errMsg
                return "error:" & errMsg
            end try
        end tell
        """

        result = await NoteIDUtils.execute_applescript(script)

        if result.startswith("error:"):
            raise RuntimeError(f"Failed to get notes: {result[6:]}")

        return NoteIDUtils._parse_notes_list(result, folder_name)

    @staticmethod
    async def _get_all_notes_in_nested_path(folder_path: str) -> list[dict[str, str]]:
        """Get all notes in a nested folder path."""
        # Validate that the folder path exists
        path_exists = await ValidationUtils.check_path_exists(folder_path)
        if not path_exists:
            raise RuntimeError(f"Folder path '{folder_path}' does not exist.")

        # Parse the path components
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
                
                set noteList to {{}}
                repeat with theNote in notes of currentFolder
                    set noteInfo to {{id:(id of theNote as string), name:(name of theNote as string), folder:"{folder_path}"}}
                    copy noteInfo to end of noteList
                end repeat
                
                return noteList
            on error errMsg
                return "error:iCloud account not available. Please enable iCloud Notes sync - " & errMsg
            end try
        end tell
        """

        result = await NoteIDUtils.execute_applescript(script)

        if result.startswith("error:"):
            raise RuntimeError(f"Failed to get notes: {result[6:]}")

        return NoteIDUtils._parse_notes_list(result, folder_path)

    @staticmethod
    def _parse_notes_list(result: str, folder_path: str) -> list[dict[str, str]]:
        """Parse the AppleScript result and return list of notes."""
        try:
            notes = []
            # Parse the AppleScript result format: {id:noteId, name:noteName, folder:folderPath}, {id:noteId2, name:noteName2, folder:folderPath}

            # Split by }, { to get individual note entries
            note_entries = result.split("}, {")

            for entry in note_entries:
                # Clean up the entry
                entry = entry.strip()
                if entry.startswith("{"):
                    entry = entry[1:]
                if entry.endswith("}"):
                    entry = entry[:-1]

                # Parse individual note data
                note_info = {}

                # Extract ID
                id_start = entry.find("id:") + 3
                id_end = entry.find(", name:", id_start)
                if id_end == -1:
                    id_end = entry.find(", folder:", id_start)
                note_info["id"] = entry[id_start:id_end].strip()

                # Extract name
                name_start = entry.find("name:") + 5
                name_end = entry.find(", folder:", name_start)
                if name_end == -1:
                    name_end = len(entry)
                name = entry[name_start:name_end].strip()
                # Remove quotes if present
                if name.startswith('"') and name.endswith('"'):
                    name = name[1:-1]
                note_info["name"] = name

                # Use the provided folder_path instead of parsing from result
                note_info["folder"] = folder_path

                notes.append(note_info)

            return notes

        except Exception as e:
            raise RuntimeError(f"Failed to parse notes list: {str(e)}")

    @staticmethod
    def extract_primary_key(full_note_id: str) -> str:
        """Extract just the primary key from a full Core Data ID.

        Args:
            full_note_id: Full Core Data ID like "x-coredata://UUID/ICNote/p123"

        Returns:
            Just the primary key part like "p123"
        """
        try:
            # Split by '/' and get the last part
            parts = full_note_id.split("/")
            if len(parts) > 0:
                return parts[-1]  # Gets "p123" part
            return full_note_id
        except:
            return full_note_id

    @staticmethod
    def extract_folder_primary_key(full_folder_id: str) -> str:
        """Extract just the primary key from a full folder Core Data ID.

        Args:
            full_folder_id: Full Core Data ID like "x-coredata://UUID/ICFolder/p2330"

        Returns:
            Just the primary key part like "p2330"
        """
        try:
            # Split by '/' and get the last part
            parts = full_folder_id.split("/")
            if len(parts) > 0:
                return parts[-1]  # Gets "p2330" part
            return full_folder_id
        except:
            return full_folder_id

    @staticmethod
    async def get_folder_name_by_id(folder_id: str) -> str:
        """Get folder name by its primary key ID.

        Args:
            folder_id: Primary key ID of the folder (e.g., "p2330")

        Returns:
            Folder name

        Raises:
            ValueError: If folder ID is empty or invalid
            RuntimeError: If folder not found by ID
        """
        # Validate folder ID
        if not folder_id or not folder_id.strip():
            raise ValueError("Folder ID cannot be empty or contain only whitespace")

        folder_id = folder_id.strip()

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

        sample_result = await NoteIDUtils.execute_applescript(script_get_uuid)
        if sample_result.startswith("error:"):
            raise RuntimeError(
                f"Could not get store UUID for folder lookup: {sample_result[6:]}"
            )

        # Extract store UUID from sample ID
        store_uuid = sample_result.split("//")[1].split("/")[0]
        full_folder_id = f"x-coredata://{store_uuid}/ICFolder/{folder_id}"

        # Get folder name by ID
        script = f"""
        tell application "Notes"
            try
                set primaryAccount to account "iCloud"
                set targetFolder to folder id "{full_folder_id}"
                
                set folderName to name of targetFolder as string
                return "success:" & folderName
            on error errMsg
                return "error:iCloud account not available. Please enable iCloud Notes sync - " & errMsg
            end try
        end tell
        """

        result = await NoteIDUtils.execute_applescript(script)

        if result.startswith("error:"):
            error_msg = result[6:]
            if "Folder not found" in error_msg or "not found" in error_msg:
                raise RuntimeError(f"Folder with ID '{folder_id}' not found")
            else:
                raise RuntimeError(f"Failed to get folder name: {error_msg}")

        if result.startswith("success:"):
            return result[8:]  # Remove "success:" prefix
        else:
            raise RuntimeError(f"Unexpected result format: {result}")
