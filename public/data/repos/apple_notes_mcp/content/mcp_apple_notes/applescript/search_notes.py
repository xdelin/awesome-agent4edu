import subprocess

from .base_operations import BaseAppleScriptOperations


class SearchNotesOperations(BaseAppleScriptOperations):
    """Operations for searching notes in Apple Notes."""

    @staticmethod
    async def search_notes(keywords: list[str]) -> list[dict[str, str]]:
        """Search for notes containing the specified keywords.

        Args:
            keywords: List of keywords to search for

        Returns:
            List of dictionaries with note_id, name, folder, and match info
        """
        try:
            # Convert keywords list to AppleScript format
            keywords_str = ", ".join([f'"{keyword}"' for keyword in keywords])

            # Build the AppleScript command
            script = f"""
            tell application "Notes"
                try
                    set primaryAccount to account "iCloud"
                    set keywords to {{{keywords_str}}}
                    set foundNotes to {{}}
                    
                    repeat with currentNote in every note of primaryAccount
                        set noteName to name of currentNote as string
                        set noteID to id of currentNote as string
                        set noteBody to body of currentNote as string
                        
                        -- Get the folder name for this note
                        set noteFolder to "Notes"
                        try
                            set noteFolder to name of container of currentNote as string
                        on error
                            set noteFolder to "Notes"
                        end try
                        
                        -- Check if any keyword matches
                        set foundKeyword to ""
                        repeat with keyword in keywords
                            if noteBody contains keyword then
                                set foundKeyword to keyword as string
                                exit repeat
                            end if
                        end repeat
                        
                        -- If keyword found, add to results
                        if foundKeyword is not "" then
                            copy noteID to end of foundNotes
                            copy noteName to end of foundNotes
                            copy noteFolder to end of foundNotes
                            copy foundKeyword to end of foundNotes
                        end if
                    end repeat
                    
                    return foundNotes
                on error errMsg
                    return "error:iCloud account not available. Please enable iCloud Notes sync - " & errMsg
                end try
            end tell
            """

            # Execute the AppleScript
            result = await SearchNotesOperations.execute_applescript(script)

            if result.startswith("error:"):
                raise RuntimeError(f"Failed to search notes: {result[6:]}")

            return SearchNotesOperations._parse_search_results(result)

        except subprocess.CalledProcessError as e:
            # Handle AppleScript errors
            error_message = e.stderr.decode("utf-8") if e.stderr else str(e)
            raise RuntimeError(f"Failed to search notes: {error_message}")
        except Exception as e:
            # Handle other errors
            raise RuntimeError(f"Unexpected error searching notes: {str(e)}")

    @staticmethod
    def _parse_search_results(result: str) -> list[dict[str, str]]:
        """Parse AppleScript result into list of search result dictionaries.

        Args:
            result: Raw AppleScript result

        Returns:
            List of dictionaries with note_id, name, folder, and matched_keyword
        """
        notes = []

        if not result or result == "{}":
            return notes

        # AppleScript returns lists as comma-separated values
        # Format: id1, name1, folder1, keyword1, id2, name2, folder2, keyword2, ...
        try:
            # Split by comma and process quadruplets
            parts = [part.strip() for part in result.split(",") if part.strip()]

            # Process quadruplets (id, name, folder, keyword)
            for i in range(0, len(parts), 4):
                if i + 3 < len(parts):
                    full_note_id = parts[i].strip('"')
                    note_name = parts[i + 1].strip('"')
                    folder_name = parts[i + 2].strip('"')
                    matched_keyword = parts[i + 3].strip('"')

                    # Extract just the primary key
                    from .note_id_utils import NoteIDUtils

                    short_id = NoteIDUtils.extract_primary_key(full_note_id)

                    # Skip empty entries and invalid IDs
                    if short_id and note_name and not short_id.startswith("item"):
                        notes.append(
                            {
                                "note_id": short_id,
                                "name": note_name,
                                "folder": folder_name,
                                "matched_keyword": matched_keyword,
                            }
                        )

        except Exception:
            # If parsing fails, return empty list
            return []

        return notes
