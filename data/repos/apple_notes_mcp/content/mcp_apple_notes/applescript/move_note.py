
from .base_operations import BaseAppleScriptOperations
from .validation_utils import ValidationUtils


class MoveNoteOperations(BaseAppleScriptOperations):
    """Operations for moving Apple Notes between folders."""

    @staticmethod
    async def move_note_by_id_and_name(
        note_id: str, note_name: str, target_folder_path: str
    ) -> dict[str, str]:
        """Move a note by its primary key ID with AppleScript verification.

        This method relies on AppleScript's built-in verification:
        1. Validates input parameters (ID, note name, target path)
        2. AppleScript verifies ID and name match the same note
        3. Performs the move operation if verification passes

        Args:
            note_id: Primary key ID of the note to move (e.g., "p1308")
            note_name: Name of the note to verify and move
            target_folder_path: Target folder path where to move the note

        Returns:
            Move result with status and details

        Raises:
            ValueError: If inputs are invalid
            RuntimeError: If note not found, name doesn't match, or move fails
        """
        # Validate inputs
        if not note_id or not note_id.strip():
            raise ValueError("Note ID cannot be empty or contain only whitespace")

        if not note_name or not note_name.strip():
            raise ValueError("Note name cannot be empty or contain only whitespace")

        if not target_folder_path or not target_folder_path.strip():
            raise ValueError(
                "Target folder path cannot be empty or contain only whitespace"
            )

        note_id = note_id.strip()
        note_name = note_name.strip()
        target_folder_path = target_folder_path.strip()

        # Validate target folder path
        target_folder_path = ValidationUtils.validate_folder_path(target_folder_path)

        # Check if target folder path exists
        target_exists = await ValidationUtils.check_path_exists(target_folder_path)
        if not target_exists:
            raise RuntimeError(
                f"Target folder path '{target_folder_path}' does not exist"
            )

        # Perform move operation using AppleScript verification
        return await MoveNoteOperations._perform_move_operation_by_id_and_name(
            note_id, note_name, target_folder_path
        )

    @staticmethod
    async def _verify_note_in_folder(note_id: str, folder_path: str) -> bool:
        """Verify that a note exists in the specified folder."""
        # Escape parameters for AppleScript
        escaped_note_id = ValidationUtils.create_applescript_quoted_string(note_id)

        # Handle root level vs nested path
        if not folder_path or folder_path.strip() == "":
            # Root level - check in root folders
            script = f"""
            tell application "Notes"
                try
                    set primaryAccount to account "iCloud"
                    repeat with rootFolder in folders of primaryAccount
                        repeat with noteItem in notes of rootFolder
                            set noteId to id of noteItem as string
                            if noteId ends with {escaped_note_id} then
                                return "true"
                            end if
                        end repeat
                    end repeat
                    return "false"
                on error errMsg
                    return "error:iCloud account not available. Please enable iCloud Notes sync - " & errMsg
                end try
            end tell
            """
        else:
            # Nested path - navigate to folder
            path_components = ValidationUtils.parse_folder_path(folder_path)

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
                    
                    -- Get note ID from the target folder
                    if currentFolder is not missing value then
                        repeat with noteItem in notes of currentFolder
                            set noteId to id of noteItem as string
                            if noteId ends with {escaped_note_id} then
                                return noteId
                            end if
                        end repeat
                    end if
                    
                    return "error:Note not found"
                    
                on error errMsg
                    return "error:iCloud account not available. Please enable iCloud Notes sync - " & errMsg
                end try
            end tell
            """

        result = await MoveNoteOperations.execute_applescript(script)

        if result.startswith("error:"):
            raise RuntimeError(f"Error verifying note: {result[6:]}")

        return result.strip() == "true"

    @staticmethod
    async def _perform_move_operation_by_id_and_name(
        note_id: str, note_name: str, target_folder_path: str
    ) -> dict[str, str]:
        """Perform the move operation using AppleScript with ID and name verification."""
        # Build full Core Data ID from primary key using dynamic store UUID
        script_get_uuid = """
        tell application "Notes"
            try
                set primaryAccount to account "iCloud"
                set sampleNote to note 1 of primaryAccount
                set sampleId to id of sampleNote as string
                return sampleId
            on error errMsg
                return "error:iCloud account not available. Please enable iCloud Notes sync - " & errMsg
            end try
        end tell
        """

        sample_result = await MoveNoteOperations.execute_applescript(script_get_uuid)
        if sample_result.startswith("error:"):
            raise RuntimeError(
                f"Could not get store UUID for moving: {sample_result[6:]}"
            )

        # Extract store UUID from sample ID
        store_uuid = sample_result.split("//")[1].split("/")[0]
        full_note_id = f"x-coredata://{store_uuid}/ICNote/{note_id}"

        # Escape parameters for AppleScript
        escaped_full_note_id = ValidationUtils.create_applescript_quoted_string(
            full_note_id
        )
        escaped_target_folder = ValidationUtils.create_applescript_quoted_string(
            target_folder_path
        )
        escaped_note_name = ValidationUtils.create_applescript_quoted_string(note_name)

        # Move operation with name verification
        script = f"""
        tell application "Notes"
            try
                set primaryAccount to account "iCloud"
                set targetNote to note id {escaped_full_note_id}
                
                set actualNoteName to name of targetNote as string
                
                -- Verify the note name matches
                if actualNoteName is not {escaped_note_name} then
                    return "error:Note name mismatch. Expected: {note_name}, Found: " & actualNoteName
                end if
                
                -- Move the note to target folder
                move targetNote to folder {escaped_target_folder}
                
                return "success:" & actualNoteName & ", " & "{note_id}" & ", " & {escaped_target_folder}
            on error errMsg
                return "error:iCloud account not available. Please enable iCloud Notes sync - " & errMsg
            end try
        end tell
        """

        result = await MoveNoteOperations.execute_applescript(script)

        if result.startswith("error:"):
            error_msg = result[6:]
            if "Note name mismatch" in error_msg:
                raise RuntimeError(error_msg)
            elif "Note not found" in error_msg or "not found" in error_msg:
                raise RuntimeError(f"Note with ID '{note_id}' not found")
            else:
                raise RuntimeError(f"Failed to move note: {error_msg}")

        return MoveNoteOperations._parse_move_by_id_and_name_result(
            result, note_id, note_name, target_folder_path
        )

    @staticmethod
    async def _perform_move_operation(
        note_id: str, source_folder_path: str, target_folder_path: str
    ) -> dict[str, str]:
        """Perform the move operation using simple AppleScript."""
        # Get the full note ID (we need the complete x-coredata:// URL)
        full_note_id = await MoveNoteOperations._get_full_note_id(
            note_id, source_folder_path
        )

        # Escape parameters for AppleScript
        escaped_full_note_id = ValidationUtils.create_applescript_quoted_string(
            full_note_id
        )
        escaped_target_folder = ValidationUtils.create_applescript_quoted_string(
            target_folder_path
        )

        # Simple move operation
        script = f"""
        tell application "Notes"
            try
                move note id {escaped_full_note_id} to folder {escaped_target_folder}
                return "moved:success:{note_id}:{source_folder_path}:{target_folder_path}"
            on error errMsg
                return "error:" & errMsg
            end try
        end tell
        """

        result = await MoveNoteOperations.execute_applescript(script)

        if result.startswith("error:"):
            raise RuntimeError(f"Failed to move note: {result[6:]}")

        return MoveNoteOperations._parse_move_result(
            result, note_id, source_folder_path, target_folder_path
        )

    @staticmethod
    async def _get_full_note_id(note_id: str, folder_path: str) -> str:
        """Get the full note ID (x-coredata:// URL) from the short ID."""
        # Escape parameters for AppleScript
        escaped_note_id = ValidationUtils.create_applescript_quoted_string(note_id)

        # Handle root level vs nested path
        if not folder_path or folder_path.strip() == "":
            # Root level - check in root folders
            script = f"""
            tell application "Notes"
                try
                    set primaryAccount to account "iCloud"
                    repeat with rootFolder in folders of primaryAccount
                        repeat with noteItem in notes of rootFolder
                            set noteId to id of noteItem as string
                            if noteId ends with {escaped_note_id} then
                                return noteId
                            end if
                        end repeat
                    end repeat
                    return "error:Note not found"
                on error errMsg
                    return "error:" & errMsg
                end try
            end tell
            """
        else:
            # Nested path - navigate to folder
            path_components = ValidationUtils.parse_folder_path(folder_path)

            script = f"""
            tell application "Notes"
                try
                    set primaryAccount to account "iCloud"
                    set currentFolder to missing value
                    set pathComponents to {{{", ".join([ValidationUtils.create_applescript_quoted_string(component) for component in path_components])}}}
                    
                    repeat with i from 1 to count of pathComponents
                        set componentName to item i of pathComponents
                        
                        if currentFolder is missing value then
                            -- Start from root folders
                            repeat with rootFolder in folders
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
                    
                    -- Get note ID from the target folder
                    if currentFolder is not missing value then
                        repeat with noteItem in notes of currentFolder
                            set noteId to id of noteItem as string
                            if noteId ends with {escaped_note_id} then
                                return noteId
                            end if
                        end repeat
                    end if
                    
                    return "error:Note not found"
                    
                on error errMsg
                    return "error:" & errMsg
                end try
            end tell
            """

        result = await MoveNoteOperations.execute_applescript(script)

        if result.startswith("error:"):
            raise RuntimeError(f"Error getting full note ID: {result[6:]}")

        return result.strip()

    @staticmethod
    def _parse_move_by_id_and_name_result(
        result: str, note_id: str, note_name: str, target_folder_path: str
    ) -> dict[str, str]:
        """Parse the AppleScript result for move by ID and name and return structured data."""
        try:
            # The result format is: success:name, note_id, target_folder
            if result.startswith("success:"):
                parts = result[8:].split(", ")  # Remove "success:" prefix

                if len(parts) >= 3:
                    return {
                        "name": parts[0],
                        "note_id": parts[1],
                        "target_folder": parts[2],
                        "status": "moved",
                        "message": "Note moved successfully",
                        "move_method": "by_id_and_name",
                    }
                else:
                    return {
                        "name": note_name,
                        "note_id": note_id,
                        "target_folder": target_folder_path,
                        "status": "moved",
                        "message": "Note moved successfully",
                        "move_method": "by_id_and_name",
                    }
            else:
                raise RuntimeError(f"Unexpected result format: {result}")
        except Exception as e:
            raise RuntimeError(f"Failed to parse move by ID and name result: {str(e)}")

    @staticmethod
    def _parse_move_result(
        result: str, note_id: str, source_folder_path: str, target_folder_path: str
    ) -> dict[str, str]:
        """Parse the AppleScript result and return structured data."""
        try:
            if result.startswith("moved:success:"):
                # Parse the moved result format: "moved:success:noteId:sourceFolder:targetFolder"
                parts = result.split(":")
                if len(parts) >= 5:
                    note_id = parts[2]
                    source_folder = parts[3]
                    target_folder = parts[4]

                    return {
                        "name": f"Note {note_id}",
                        "note_id": note_id,
                        "source_folder": source_folder,
                        "target_folder": target_folder,
                        "status": "moved",
                        "message": "Note moved successfully",
                        "creation_date": "N/A",
                        "modification_date": "N/A",
                    }

            # Fallback to input-based response
            return {
                "name": f"Note {note_id}",
                "note_id": note_id,
                "source_folder": source_folder_path or "Root",
                "target_folder": target_folder_path or "Root",
                "status": "moved",
                "message": "Note moved successfully",
                "creation_date": "N/A",
                "modification_date": "N/A",
            }
        except Exception as e:
            raise RuntimeError(f"Error parsing move result: {str(e)}")
