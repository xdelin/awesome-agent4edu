from .base_operations import BaseAppleScriptOperations
from .note_id_utils import NoteIDUtils


class NotesStructureOperations(BaseAppleScriptOperations):
    """Operations for getting complete Apple Notes structure with notes included."""

    @staticmethod
    async def get_notes_structure() -> str:
        """Get complete notes structure with folders and notes - return raw AppleScript data."""
        script = """
tell application "Notes"
	try
		set structureList to {}
		
		-- Direct iCloud account access
		set primaryAccount to account "iCloud"
		
		-- Get all root folders and their contents from iCloud account only
		repeat with rootFolder in folders of primaryAccount
			set rootName to name of rootFolder
			set rootId to id of rootFolder as string
			set structureList to structureList & {"Root Folder: " & rootName & " (ID: " & rootId & ")"}
			
			-- Get notes in root folder
			repeat with theNote in notes of rootFolder
				set noteName to name of theNote
				set structureList to structureList & {"  ├── Note: " & noteName}
			end repeat
			
			-- Get level 2 subfolders and their contents
			repeat with subFolder in folders of rootFolder
				set subName to name of subFolder
				set subId to id of subFolder as string
				set structureList to structureList & {"  ├── Subfolder: " & subName & " (ID: " & subId & ")"}
				
				-- Get notes in level 2 subfolder
				repeat with theNote in notes of subFolder
					set noteName to name of theNote
					set structureList to structureList & {"    ├── Note: " & noteName}
				end repeat
				
				-- Get level 3 subfolders and their contents
				repeat with subSubFolder in folders of subFolder
					set subSubName to name of subSubFolder
					set subSubId to id of subSubFolder as string
					set structureList to structureList & {"    ├── Sub-subfolder: " & subSubName & " (ID: " & subSubId & ")"}
					
					-- Get notes in level 3 subfolder
					repeat with theNote in notes of subSubFolder
						set noteName to name of theNote
						set structureList to structureList & {"      ├── Note: " & noteName}
					end repeat
					
					-- Get level 4 subfolders and their contents
					repeat with subSubSubFolder in folders of subSubFolder
						set subSubSubName to name of subSubSubFolder
						set subSubSubId to id of subSubSubFolder as string
						set structureList to structureList & {"      ├── Sub-sub-subfolder: " & subSubSubName & " (ID: " & subSubSubId & ")"}
						
						-- Get notes in level 4 subfolder
						repeat with theNote in notes of subSubSubFolder
							set noteName to name of theNote
							set structureList to structureList & {"        ├── Note: " & noteName}
						end repeat
						
						-- Get level 5 subfolders and their contents
						repeat with subSubSubSubFolder in folders of subSubSubFolder
							set subSubSubName to name of subSubSubSubFolder
							set subSubSubSubId to id of subSubSubSubFolder as string
							set structureList to structureList & {"        ├── Sub-sub-sub-subfolder: " & subSubSubName & " (ID: " & subSubSubSubId & ")"}
							
							-- Get notes in level 5 subfolder
							repeat with theNote in notes of subSubSubSubFolder
								set noteName to name of theNote
								set structureList to structureList & {"          ├── Note: " & noteName}
							end repeat
						end repeat
					end repeat
				end repeat
			end repeat
			
			set structureList to structureList & {""}
		end repeat
		
		-- Convert to string
		set AppleScript's text item delimiters to return
		set notesStructure to structureList as string
		set AppleScript's text item delimiters to ""
		
		return notesStructure
		
	on error errMsg
		return "error:iCloud account not available. Please enable iCloud Notes sync - " & errMsg
	end try
end tell
        """

        result = await NotesStructureOperations.execute_applescript(script)

        if result.startswith("error:"):
            raise RuntimeError(f"Failed to get notes structure: {result[6:]}")

        return result

    @staticmethod
    async def get_filtered_notes_structure() -> str:
        """Get filtered notes structure - remove root folders whose IDs appear in subfolders, but keep IDs visible."""
        # First get the complete notes structure
        complete_structure = await NotesStructureOperations.get_notes_structure()

        # Parse the structure to extract IDs
        lines = complete_structure.split("\r")
        root_folder_ids = set()
        subfolder_ids = set()

        for line in lines:
            line = line.strip()
            if not line:
                continue

            # Extract ID from the line
            if "(ID:" in line:
                # Extract the ID part
                start_idx = line.find("(ID:") + 4
                end_idx = line.find(")", start_idx)
                if start_idx > 3 and end_idx > start_idx:
                    folder_id = line[start_idx:end_idx].strip()

                    # Determine if it's a root folder or subfolder based on indentation
                    if line.startswith("Root Folder:"):
                        root_folder_ids.add(folder_id)
                    elif (
                        "├── Subfolder:" in line
                        or "├── Sub-subfolder:" in line
                        or "├── Sub-sub-subfolder:" in line
                        or "├── Sub-sub-sub-subfolder:" in line
                    ):
                        subfolder_ids.add(folder_id)

        # Filter out root folders whose IDs appear in subfolders
        filtered_root_ids = root_folder_ids - subfolder_ids

        # Rebuild the structure with only filtered root folders, but keep IDs visible
        filtered_lines = []
        current_root_id = None
        include_current_root = False

        for line in lines:
            line = line.strip()
            if not line:
                if include_current_root:
                    filtered_lines.append("")
                current_root_id = None
                include_current_root = False
                continue

            # Check if this is a root folder line
            if line.startswith("Root Folder:"):
                # Extract the ID
                if "(ID:" in line:
                    start_idx = line.find("(ID:") + 4
                    end_idx = line.find(")", start_idx)
                    if start_idx > 3 and end_idx > start_idx:
                        current_root_id = line[start_idx:end_idx].strip()
                        include_current_root = current_root_id in filtered_root_ids

                if include_current_root:
                    # Keep the ID in the display for root folders
                    clean_line = line.replace("Root Folder:", "").strip()
                    # Extract primary key from full ID
                    if "(ID:" in clean_line:
                        name_part = clean_line.split(" (ID:")[0]
                        full_id_part = clean_line.split(" (ID:")[1].split(")")[0]
                        primary_key = NoteIDUtils.extract_folder_primary_key(
                            full_id_part
                        )
                        clean_line = f"{name_part} (ID: {primary_key})"
                    filtered_lines.append(clean_line)
            elif include_current_root and line:
                # Include subfolder and note lines if we're including the current root
                # Keep ID in folder lines for display and format with tree symbols
                if "(ID:" in line:
                    # Extract the folder name and ID
                    name_part = line.split(" (ID:")[0]
                    full_id_part = line.split(" (ID:")[1].split(")")[0]
                    primary_key = NoteIDUtils.extract_folder_primary_key(full_id_part)

                    # Determine the level and format accordingly
                    if "├── Subfolder:" in name_part:
                        # Level 1: ├──
                        folder_name = name_part.replace("├── Subfolder:", "").strip()
                        formatted_line = f"├── {folder_name} (ID: {primary_key})"
                    elif "├── Sub-subfolder:" in name_part:
                        # Level 2: │   ├──
                        folder_name = name_part.replace(
                            "├── Sub-subfolder:", ""
                        ).strip()
                        formatted_line = f"│   ├── {folder_name} (ID: {primary_key})"
                    elif "├── Sub-sub-subfolder:" in name_part:
                        # Level 3: │   │   ├──
                        folder_name = name_part.replace(
                            "├── Sub-sub-subfolder:", ""
                        ).strip()
                        formatted_line = (
                            f"│   │   ├── {folder_name} (ID: {primary_key})"
                        )
                    elif "├── Sub-sub-sub-subfolder:" in name_part:
                        # Level 4: │   │   │   ├──
                        folder_name = name_part.replace(
                            "├── Sub-sub-sub-subfolder:", ""
                        ).strip()
                        formatted_line = (
                            f"│   │   │   ├── {folder_name} (ID: {primary_key})"
                        )
                    else:
                        formatted_line = f"{name_part} (ID: {primary_key})"

                    filtered_lines.append(formatted_line)
                else:
                    # Handle note lines
                    if "├── Note:" in line:
                        note_name = line.replace("├── Note:", "").strip()
                        # Determine indentation level based on the original line
                        if line.startswith("  ├── Note:"):
                            formatted_line = f"├── {note_name}"
                        elif line.startswith("    ├── Note:"):
                            formatted_line = f"│   ├── {note_name}"
                        elif line.startswith("      ├── Note:"):
                            formatted_line = f"│   │   ├── {note_name}"
                        elif line.startswith("        ├── Note:"):
                            formatted_line = f"│   │   │   ├── {note_name}"
                        elif line.startswith("          ├── Note:"):
                            formatted_line = f"│   │   │   │   ├── {note_name}"
                        else:
                            formatted_line = f"├── {note_name}"
                        filtered_lines.append(formatted_line)
                    else:
                        filtered_lines.append(line)

        # Convert back to string
        return "\r".join(filtered_lines)
