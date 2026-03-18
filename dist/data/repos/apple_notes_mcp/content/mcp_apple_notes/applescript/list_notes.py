
from .base_operations import BaseAppleScriptOperations
from .note_id_utils import NoteIDUtils


class ListNotesOperations(BaseAppleScriptOperations):
    """Operations for listing notes across all folders."""

    @staticmethod
    async def list_all_notes() -> list[dict[str, str]]:
        """Get a list of all notes across all folders with their IDs and names.

        Returns:
            List of dictionaries with note_id, name, and folder info

        Raises:
            RuntimeError: If AppleScript execution fails
        """
        script = """
        tell application "Notes"
            set outputText to ""
            set iCloudAccount to account "iCloud"
            repeat with currentFolder in folders of iCloudAccount
                repeat with currentNote in notes of currentFolder
                    set noteName to name of currentNote
                    set noteID to id of currentNote
                    set noteInfo to ("Name: " & noteName & " | ID: " & noteID)
                    if outputText is "" then
                        set outputText to noteInfo
                    else
                        set outputText to outputText & return & noteInfo
                    end if
                end repeat
            end repeat
            return outputText
        end tell
        """

        result = await ListNotesOperations.execute_applescript(script)

        if result.startswith("error:"):
            raise RuntimeError(f"Failed to list all notes: {result[6:]}")

        return ListNotesOperations._parse_notes_list(result)

    @staticmethod
    def _parse_notes_list(result: str) -> list[dict[str, str]]:
        """Parse the AppleScript result into a list of note dictionaries.

        Args:
            result: Raw AppleScript result string in format
                    "Name: Note Name | ID: full_id"

        Returns:
            List of dictionaries with note_id and name
        """
        if not result.strip():
            return []

        notes: list[dict[str, str]] = []
        for line in result.strip().split("\r"):
            line = line.strip()
            if not line:
                continue
            try:
                # Split into parts like {"Name": "...", "ID": "..."}
                parts = dict(item.strip().split(": ", 1) for item in line.split(" | "))
                full_note_id = parts["ID"]
                notes.append(
                    {
                        "note_id": NoteIDUtils.extract_primary_key(full_note_id),
                        "name": parts["Name"],
                    }
                )
            except (ValueError, KeyError):
                # Skip malformed entries
                continue

        return notes
