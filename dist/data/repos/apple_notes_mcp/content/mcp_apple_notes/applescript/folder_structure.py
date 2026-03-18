from .base_operations import BaseAppleScriptOperations
from .note_id_utils import NoteIDUtils


class FolderStructureOperations(BaseAppleScriptOperations):
    """Operations for getting complete Apple Notes folder structure."""

    @staticmethod
    async def get_folders_structure() -> str:
        """Get complete folder structure from primary account only - return raw AppleScript data."""
        script = """
tell application "Notes"
	try
		set folderList to {}
		
		-- Direct iCloud account access
		set primaryAccount to account "iCloud"
		
		-- Get all root folders and their subfolders from iCloud account only
		repeat with rootFolder in folders of primaryAccount
			set rootName to name of rootFolder
			set rootId to id of rootFolder as string
			set folderList to folderList & {"Root Folder: " & rootName & " (ID: " & rootId & ")"}
			
			-- Get level 2 subfolders
			repeat with subFolder in folders of rootFolder
				set subName to name of subFolder
				set subId to id of subFolder as string
				set folderList to folderList & {"  ├── Subfolder: " & subName & " (ID: " & subId & ")"}
				
				-- Get level 3 subfolders
				repeat with subSubFolder in folders of subFolder
					set subSubName to name of subSubFolder
					set subSubId to id of subSubFolder as string
					set folderList to folderList & {"    ├── Sub-subfolder: " & subSubName & " (ID: " & subSubId & ")"}
					
					-- Get level 4 subfolders
					repeat with subSubSubFolder in folders of subSubFolder
						set subSubSubName to name of subSubSubFolder
						set subSubSubId to id of subSubSubFolder as string
						set folderList to folderList & {"      ├── Sub-sub-subfolder: " & subSubSubName & " (ID: " & subSubSubId & ")"}
						
						-- Get level 5 subfolders
						repeat with subSubSubSubFolder in folders of subSubSubFolder
							set subSubSubName to name of subSubSubSubFolder
							set subSubSubSubId to id of subSubSubSubFolder as string
							set folderList to folderList & {"        ├── Sub-sub-sub-subfolder: " & subSubSubName & " (ID: " & subSubSubSubId & ")"}
						end repeat
					end repeat
				end repeat
			end repeat
			
			set folderList to folderList & {""}
		end repeat
		
		-- Convert to string
		set AppleScript's text item delimiters to return
		set folderStructure to folderList as string
		set AppleScript's text item delimiters to ""
		
		return folderStructure
		
	on error errMsg
		return "error:iCloud account not available. Please enable iCloud Notes sync - " & errMsg
	end try
end tell
        """

        result = await FolderStructureOperations.execute_applescript(script)

        if result.startswith("error:"):
            raise RuntimeError(f"Failed to get folders structure: {result[6:]}")

        return result

    @staticmethod
    async def get_filtered_folders_structure() -> str:
        """Get filtered folder structure - remove root folders whose IDs appear in subfolders, but keep IDs visible."""
        # First get the complete folder structure
        complete_structure = await FolderStructureOperations.get_folders_structure()

        # Parse the structure to build a proper hierarchy
        lines = complete_structure.split("\r")

        # Build a hierarchy map to track parent-child relationships
        hierarchy = {}
        all_folder_ids = set()
        root_folder_ids = set()
        subfolder_ids = set()

        current_root = None
        current_parent = None

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
                    all_folder_ids.add(folder_id)

                    # Determine the level and parent
                    if line.startswith("Root Folder:"):
                        # This is a root folder
                        current_root = folder_id
                        current_parent = None
                        root_folder_ids.add(folder_id)
                        # Extract clean name without ID
                        clean_name = (
                            line.replace("Root Folder:", "").split(" (ID:")[0].strip()
                        )
                        hierarchy[folder_id] = {
                            "parent": None,
                            "children": [],
                            "name": clean_name,
                        }
                    elif "├── Subfolder:" in line:
                        # This is a level 1 subfolder
                        subfolder_ids.add(folder_id)
                        if current_root:
                            hierarchy[current_root]["children"].append(folder_id)
                            # Extract clean name without ID
                            clean_name = (
                                line.replace("├── Subfolder:", "")
                                .split(" (ID:")[0]
                                .strip()
                            )
                            hierarchy[folder_id] = {
                                "parent": current_root,
                                "children": [],
                                "name": clean_name,
                            }
                            current_parent = folder_id
                    elif "├── Sub-subfolder:" in line:
                        # This is a level 2 subfolder
                        subfolder_ids.add(folder_id)
                        if current_parent:
                            hierarchy[current_parent]["children"].append(folder_id)
                            # Extract clean name without ID
                            clean_name = (
                                line.replace("├── Sub-subfolder:", "")
                                .split(" (ID:")[0]
                                .strip()
                            )
                            hierarchy[folder_id] = {
                                "parent": current_parent,
                                "children": [],
                                "name": clean_name,
                            }
                    elif "├── Sub-sub-subfolder:" in line:
                        # This is a level 3 subfolder
                        subfolder_ids.add(folder_id)
                        if current_parent:
                            # Extract clean name without ID
                            clean_name = (
                                line.replace("├── Sub-sub-subfolder:", "")
                                .split(" (ID:")[0]
                                .strip()
                            )
                            hierarchy[folder_id] = {
                                "parent": current_parent,
                                "children": [],
                                "name": clean_name,
                            }
                    elif "├── Sub-sub-sub-subfolder:" in line:
                        # This is a level 4 subfolder
                        subfolder_ids.add(folder_id)
                        if current_parent:
                            # Extract clean name without ID
                            clean_name = (
                                line.replace("├── Sub-sub-sub-subfolder:", "")
                                .split(" (ID:")[0]
                                .strip()
                            )
                            hierarchy[folder_id] = {
                                "parent": current_parent,
                                "children": [],
                                "name": clean_name,
                            }

        # Filter out root folders whose IDs appear as subfolders
        filtered_root_ids = root_folder_ids - subfolder_ids

        # Build the filtered structure
        filtered_lines = []

        for root_id in filtered_root_ids:
            if root_id in hierarchy:
                root_info = hierarchy[root_id]
                root_name = root_info["name"]

                # Extract primary key from full ID
                primary_key = NoteIDUtils.extract_folder_primary_key(root_id)
                filtered_lines.append(f"{root_name} (ID: {primary_key})")

                # Add children recursively
                FolderStructureOperations._add_children_to_structure(
                    filtered_lines, root_id, hierarchy, 1
                )

                # Add empty line between root folders
                filtered_lines.append("")

        # Convert back to string
        return "\r".join(filtered_lines)

    @staticmethod
    def _add_children_to_structure(
        lines: list, parent_id: str, hierarchy: dict, level: int
    ):
        """Recursively add children to the structure with proper indentation."""
        if parent_id not in hierarchy:
            return

        parent_info = hierarchy[parent_id]
        for child_id in parent_info["children"]:
            if child_id in hierarchy:
                child_info = hierarchy[child_id]
                child_name = child_info["name"]

                # Extract primary key from full ID
                primary_key = NoteIDUtils.extract_folder_primary_key(child_id)

                # Add proper indentation based on level
                if level == 1:
                    lines.append(f"├── {child_name} (ID: {primary_key})")
                elif level == 2:
                    lines.append(f"│   ├── {child_name} (ID: {primary_key})")
                elif level == 3:
                    lines.append(f"│   │   ├── {child_name} (ID: {primary_key})")
                elif level == 4:
                    lines.append(f"│   │   │   └── {child_name} (ID: {primary_key})")

                # Recursively add children
                FolderStructureOperations._add_children_to_structure(
                    lines, child_id, hierarchy, level + 1
                )
