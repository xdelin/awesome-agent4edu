
from .base_operations import BaseAppleScriptOperations


class DeleteNoteOperations(BaseAppleScriptOperations):
    """Operations for deleting Apple Notes by ID."""

    @staticmethod
    async def delete_note_by_id_and_name(
        note_id: str, note_name: str
    ) -> dict[str, str]:
        """Delete a note by its primary key ID with AppleScript verification.

        This method relies on AppleScript's built-in verification:
        1. Validates input parameters (ID, note name)
        2. AppleScript verifies ID and name match the same note
        3. Performs the deletion if verification passes

        Args:
            note_id: Primary key ID of the note (e.g., "p1308")
            note_name: Name of the note to verify and delete

        Returns:
            Deletion result with status and details

        Raises:
            ValueError: If note ID or name is empty or invalid
            RuntimeError: If note not found, name doesn't match, or deletion fails
        """
        # Validate inputs
        if not note_id or not note_id.strip():
            raise ValueError("Note ID cannot be empty or contain only whitespace")

        if not note_name or not note_name.strip():
            raise ValueError("Note name cannot be empty or contain only whitespace")

        note_id = note_id.strip()
        note_name = note_name.strip()

        # Build full Core Data ID from primary key using dynamic store UUID
        # First get a sample note to extract the store UUID
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

        sample_result = await DeleteNoteOperations.execute_applescript(script_get_uuid)
        if sample_result.startswith("error:"):
            raise RuntimeError(
                f"Could not get store UUID for deletion: {sample_result[6:]}"
            )

        # Extract store UUID from sample ID
        store_uuid = sample_result.split("//")[1].split("/")[0]
        full_note_id = f"x-coredata://{store_uuid}/ICNote/{note_id}"

        # Perform the deletion using the full note ID with name verification
        script = f"""
        tell application "Notes"
            try
                set primaryAccount to account "iCloud"
                set targetNote to note id "{full_note_id}"
                
                set actualNoteName to name of targetNote as string
                set noteId to id of targetNote as string
                
                -- Verify the note name matches
                if actualNoteName is not "{note_name}" then
                    return "error:Note name mismatch. Expected: {note_name}, Found: " & actualNoteName
                end if
                
                delete targetNote
                
                return "success:" & actualNoteName & ", " & noteId
            on error errMsg
                return "error:iCloud account not available. Please enable iCloud Notes sync - " & errMsg
            end try
        end tell
        """

        result = await DeleteNoteOperations.execute_applescript(script)

        if result.startswith("error:"):
            error_msg = result[6:]
            if "Note name mismatch" in error_msg:
                raise RuntimeError(error_msg)
            elif "Note not found" in error_msg or "not found" in error_msg:
                raise RuntimeError(f"Note with ID '{note_id}' not found")
            else:
                raise RuntimeError(f"Failed to delete note: {error_msg}")

        return DeleteNoteOperations._parse_delete_by_id_and_name_result(
            result, note_id, note_name
        )

    @staticmethod
    def _parse_delete_by_id_and_name_result(
        result: str, note_id: str, note_name: str
    ) -> dict[str, str]:
        """Parse the AppleScript result for delete by ID and name and return deletion information."""
        try:
            # The result format is: success:name, full_id
            if result.startswith("success:"):
                parts = result[8:].split(", ")  # Remove "success:" prefix

                if len(parts) >= 2:
                    from .note_id_utils import NoteIDUtils

                    # Extract primary key from full ID
                    primary_key = NoteIDUtils.extract_primary_key(parts[1])

                    return {
                        "name": parts[0],
                        "note_id": primary_key,
                        "status": "deleted",
                        "deletion_method": "by_id_and_name",
                    }
                else:
                    return {
                        "name": note_name,
                        "note_id": note_id,
                        "status": "deleted",
                        "deletion_method": "by_id_and_name",
                    }
            else:
                raise RuntimeError(f"Unexpected result format: {result}")
        except Exception as e:
            raise RuntimeError(
                f"Failed to parse delete by ID and name result: {str(e)}"
            )
