
from .base_operations import BaseAppleScriptOperations
from .validation_utils import ValidationUtils


class UpdateNoteOperations(BaseAppleScriptOperations):
    """Operations for updating Apple Notes by ID."""

    @staticmethod
    async def update_note_by_id_and_name(
        note_id: str, note_name: str, combined_content: str
    ) -> dict[str, str]:
        """Update an existing note by its primary key ID with AppleScript verification.

        This method relies on AppleScript's built-in verification:
        1. Validates input parameters (ID, note name, content)
        2. AppleScript verifies ID and name match the same note
        3. Performs the update if verification passes

        Args:
            note_id: Primary key ID of the note (e.g., "p1234")
            note_name: Current name of the note to verify and update
            combined_content: Complete HTML content combining title and body (e.g., '<h1>Title</h1><p>Content</p>')

        Returns:
            Updated note metadata with primary key ID

        Raises:
            ValueError: If note ID or name is empty or invalid
            RuntimeError: If note not found, name doesn't match, or update fails
        """
        # Validate inputs
        if not note_id or not note_id.strip():
            raise ValueError("Note ID cannot be empty or contain only whitespace")

        if not note_name or not note_name.strip():
            raise ValueError("Note name cannot be empty or contain only whitespace")

        if not combined_content or not combined_content.strip():
            raise ValueError(
                "Combined content cannot be empty or contain only whitespace"
            )

        note_id = note_id.strip()
        note_name = note_name.strip()

        # Validate and prepare content
        html_content = ValidationUtils.validate_note_body(combined_content)

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

        sample_result = await UpdateNoteOperations.execute_applescript(script_get_uuid)
        if sample_result.startswith("error:"):
            raise RuntimeError(f"Could not get store UUID: {sample_result[6:]}")

        # Extract store UUID from sample ID
        store_uuid = sample_result.split("//")[1].split("/")[0]
        full_note_id = f"x-coredata://{store_uuid}/ICNote/{note_id}"

        # Escape HTML content for AppleScript
        escaped_html = ValidationUtils.create_applescript_quoted_string(html_content)

        # Use AppleScript to update note by ID with name verification
        script = f"""
        tell application "Notes"
            try
                set primaryAccount to account "iCloud"
                set targetNote to note id "{full_note_id}"
                
                set actualNoteName to name of targetNote as string
                
                -- Verify the note name matches
                if actualNoteName is not "{note_name}" then
                    return "error:Note name mismatch. Expected: {note_name}, Found: " & actualNoteName
                end if
                
                set body of targetNote to {escaped_html}
                
                return "success:" & actualNoteName & ", " & "{note_id}"
            on error errMsg
                return "error:iCloud account not available. Please enable iCloud Notes sync - " & errMsg
            end try
        end tell
        """

        result = await UpdateNoteOperations.execute_applescript(script)

        if result.startswith("error:"):
            error_msg = result[6:]
            if "Note name mismatch" in error_msg:
                raise RuntimeError(error_msg)
            elif "Note not found" in error_msg or "not found" in error_msg:
                raise RuntimeError(f"Note with ID '{note_id}' not found")
            else:
                raise RuntimeError(f"Failed to update note: {error_msg}")

        return UpdateNoteOperations._parse_update_by_id_and_name_result(
            result, note_id, note_name
        )

    @staticmethod
    def _parse_update_by_id_and_name_result(
        result: str, note_id: str, note_name: str
    ) -> dict[str, str]:
        """Parse AppleScript result for update by ID and name."""
        try:
            # The result format is: success:name, note_id
            if result.startswith("success:"):
                parts = result[8:].split(", ")  # Remove "success:" prefix

                if len(parts) >= 2:
                    return {
                        "name": parts[0],
                        "note_id": parts[1],
                        "status": "updated",
                        "update_method": "by_id_and_name",
                    }
                else:
                    return {
                        "name": note_name,
                        "note_id": note_id,
                        "status": "updated",
                        "update_method": "by_id_and_name",
                    }
            else:
                raise RuntimeError(f"Unexpected result format: {result}")
        except Exception as e:
            raise RuntimeError(
                f"Failed to parse update by ID and name result: {str(e)}"
            )

    @staticmethod
    def _parse_update_result(result: str) -> dict[str, str]:
        """Parse AppleScript result for update by ID."""
        try:
            # Extract information from result
            name_start = result.find("name:") + 5
            name_end = result.find(", note_id:", name_start)
            name = result[name_start:name_end].strip()

            note_id_start = result.find("note_id:") + 8
            note_id_end = result.find(", creation_date:", note_id_start)
            note_id = result[note_id_start:note_id_end].strip()

            creation_start = result.find("creation_date:") + 14
            creation_end = result.find(", modification_date:", creation_start)
            creation_date = result[creation_start:creation_end].strip()

            modification_start = result.find("modification_date:") + 18
            modification_date = result[modification_start:].strip().rstrip("}")

            # Clean up quotes
            if name.startswith('"') and name.endswith('"'):
                name = name[1:-1]
            if note_id.startswith('"') and note_id.endswith('"'):
                note_id = note_id[1:-1]

            return {
                "name": name,
                "note_id": note_id,
                "creation_date": creation_date,
                "modification_date": modification_date,
                "status": "updated",
            }
        except Exception as e:
            raise RuntimeError(f"Failed to parse update result: {str(e)}")
