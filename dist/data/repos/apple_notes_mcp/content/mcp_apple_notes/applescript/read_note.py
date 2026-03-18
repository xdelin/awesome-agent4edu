
from .base_operations import BaseAppleScriptOperations


class ReadNoteOperations(BaseAppleScriptOperations):
    """Operations for reading Apple Notes by ID with name verification."""

    @staticmethod
    async def read_note_by_id_and_name(note_id: str, note_name: str) -> dict[str, str]:
        """Read a note by its primary key ID with AppleScript verification.

        This method relies on AppleScript's built-in verification:
        1. Validates input parameters (ID, note name)
        2. AppleScript verifies ID and name match the same note
        3. Performs the read operation if verification passes

        Args:
            note_id: Primary key ID of the note (e.g., "p1308")
            note_name: Name of the note to verify and read

        Returns:
            Note data with content, metadata, and verification info

        Raises:
            ValueError: If note ID or name is empty or invalid
            RuntimeError: If note not found, name doesn't match, or read fails
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

        sample_result = await ReadNoteOperations.execute_applescript(script_get_uuid)
        if sample_result.startswith("error:"):
            raise RuntimeError(
                f"Could not get store UUID for reading: {sample_result[6:]}"
            )

        # Extract store UUID from sample ID
        store_uuid = sample_result.split("//")[1].split("/")[0]
        full_note_id = f"x-coredata://{store_uuid}/ICNote/{note_id}"

        # Read the note using the full note ID with name verification
        script = f"""
        tell application "Notes"
            try
                set primaryAccount to account "iCloud"
                set targetNote to note id "{full_note_id}"
                
                set actualNoteName to name of targetNote
                
                -- Verify that the actual note name matches the provided name
                if actualNoteName is not equal to "{note_name}" then
                    return "error:Note name mismatch. Expected: {note_name}, Actual: " & actualNoteName
                end if
                
                set noteId to id of targetNote
                set noteBody to body of targetNote
                set creationDate to creation date of targetNote
                set modificationDate to modification date of targetNote
                
                return "success:" & actualNoteName & "|||" & noteId & "|||" & noteBody & "|||" & creationDate & "|||" & modificationDate
            on error errMsg
                return "error:iCloud account not available. Please enable iCloud Notes sync - " & errMsg
            end try
        end tell
        """

        result = await ReadNoteOperations.execute_applescript(script)

        if result.startswith("error:"):
            error_msg = result[6:]
            if "Note name mismatch" in error_msg:
                raise RuntimeError(f"Note name mismatch: {error_msg}")
            elif "Note not found" in error_msg or "not found" in error_msg:
                raise RuntimeError(f"Note with ID '{note_id}' not found")
            else:
                raise RuntimeError(f"Failed to read note: {error_msg}")

        return ReadNoteOperations._parse_read_by_id_and_name_result(result, note_id)

    @staticmethod
    def _parse_read_by_id_and_name_result(result: str, note_id: str) -> dict[str, str]:
        """Parse the AppleScript result for read by ID and name verification."""
        try:
            # The result format is: success:name|||full_id|||body|||creation_date|||modification_date
            if result.startswith("success:"):
                content = result[8:]  # Remove "success:" prefix
                parts = content.split("|||")

                if len(parts) >= 5:
                    from .note_id_utils import NoteIDUtils

                    # Extract primary key from full ID
                    primary_key = NoteIDUtils.extract_primary_key(parts[1])

                    return {
                        "name": parts[0],
                        "note_id": primary_key,
                        "body": parts[2],
                        "creation_date": parts[3],
                        "modification_date": parts[4],
                        "status": "found",
                        "read_method": "by_id_and_name",
                    }
                else:
                    return {
                        "name": "Unknown",
                        "note_id": note_id,
                        "body": "Could not retrieve content",
                        "creation_date": "Unknown",
                        "modification_date": "Unknown",
                        "status": "found",
                        "read_method": "by_id_and_name",
                    }
            else:
                raise RuntimeError(f"Unexpected result format: {result}")
        except Exception as e:
            raise RuntimeError(f"Failed to parse read by ID and name result: {str(e)}")
