from typing import Any

from mcp_apple_notes.applescript import (
    CreateFolderOperations,
    CreateNoteOperations,
    DeleteFolderOperations,
    DeleteNoteOperations,
    FolderStructureOperations,
    ListNotesOperations,
    MoveFolderOperations,
    MoveNoteOperations,
    NotesStructureOperations,
    ReadFolderOperations,
    ReadNoteOperations,
    RenameFolderOperations,
    SearchNotesOperations,
    UpdateNoteOperations,
)


class NotesTools:
    """Tools for Apple Notes operations."""

    async def create_note(
        self, name: str, body: str, folder_path: str = "Notes"
    ) -> dict[str, str]:
        """Create a new note with specified name, body, and folder path."""
        return await CreateNoteOperations.create_note(name, body, folder_path)

    async def create_folder(
        self, folder_name: str, folder_path: str = ""
    ) -> dict[str, Any]:
        """Create a folder in Apple Notes."""
        return await CreateFolderOperations.create_folder(folder_name, folder_path)

    async def read_note(self, note_id: str, note_name: str) -> dict[str, str]:
        """Read a note by its primary key ID with AppleScript verification."""
        return await ReadNoteOperations.read_note_by_id_and_name(note_id, note_name)

    async def read_folder(self, folder_id: str, folder_name: str) -> dict[str, Any]:
        """Read a folder by its primary key ID with AppleScript verification."""
        return await ReadFolderOperations.read_folder_by_id_and_name(
            folder_id, folder_name
        )

    async def rename_folder(
        self, folder_id: str, current_name: str, new_name: str
    ) -> dict[str, Any]:
        """Rename a folder in Apple Notes by ID with enhanced name verification."""
        return await RenameFolderOperations.rename_folder_by_id(
            folder_id, current_name, new_name
        )

    async def move_folder(
        self, folder_id: str, folder_name: str, target_path: str = ""
    ) -> dict[str, Any]:
        """Move a folder from one location to another in Apple Notes by ID with AppleScript verification."""
        return await MoveFolderOperations.move_folder_by_id(
            folder_id, folder_name, target_path
        )

    async def list_folder_with_structure(self) -> str:
        """List the complete folder structure with hierarchical tree format."""
        return await FolderStructureOperations.get_filtered_folders_structure()

    async def list_notes_with_structure(self) -> str:
        """List the complete folder structure with notes included in hierarchical tree format."""
        return await NotesStructureOperations.get_filtered_notes_structure()

    async def update_note(
        self, note_id: str, note_name: str, combined_content: str
    ) -> dict[str, str]:
        """Update an existing note by its primary key ID with AppleScript verification."""
        return await UpdateNoteOperations.update_note_by_id_and_name(
            note_id, note_name, combined_content
        )

    async def delete_note(self, note_id: str, note_name: str) -> dict[str, str]:
        """Delete a note by its primary key ID with AppleScript verification."""
        return await DeleteNoteOperations.delete_note_by_id_and_name(note_id, note_name)

    async def list_all_notes(self) -> list[dict[str, str]]:
        """Get a list of all notes across all folders with their IDs and names."""
        return await ListNotesOperations.list_all_notes()

    async def move_note(
        self, note_id: str, note_name: str, target_folder_path: str
    ) -> dict[str, str]:
        """Move a note from source folder to target folder by ID with AppleScript verification."""
        return await MoveNoteOperations.move_note_by_id_and_name(
            note_id, note_name, target_folder_path
        )

    async def delete_folder(self, folder_id: str, folder_name: str) -> dict[str, Any]:
        """Delete a folder in Apple Notes by ID with AppleScript verification."""
        return await DeleteFolderOperations.delete_folder_by_id_and_name(
            folder_id, folder_name
        )

    async def search_notes(self, keywords: list[str]) -> list[dict[str, str]]:
        """Search for notes containing the specified keywords."""
        return await SearchNotesOperations.search_notes(keywords)
