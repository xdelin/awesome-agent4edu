"""
MCP Apple Notes Server

A Model Context Protocol (MCP) server for Apple Notes integration using AppleScript.
"""

from .applescript import (
    BaseAppleScriptOperations,
    CreateFolderOperations,
    CreateNoteOperations,
    DeleteFolderOperations,
    DeleteNoteOperations,
    FolderStructureOperations,
    MoveFolderOperations,
    NoteIDUtils,
    NotesStructureOperations,
    ReadFolderOperations,
    ReadNoteOperations,
    RenameFolderOperations,
    UpdateNoteOperations,
    ValidationUtils,
)
from .tools.notes_tools import NotesTools

__version__ = "0.1.0"
__author__ = "Henil C Alagiya"
__email__ = "henilcalagiya@gmail.com"

__all__ = [
    "NotesTools",
    "BaseAppleScriptOperations",
    "CreateNoteOperations",
    "CreateFolderOperations",
    "ReadNoteOperations",
    "ReadFolderOperations",
    "DeleteFolderOperations",
    "RenameFolderOperations",
    "MoveFolderOperations",
    "FolderStructureOperations",
    "NotesStructureOperations",
    "UpdateNoteOperations",
    "DeleteNoteOperations",
    "ValidationUtils",
    "NoteIDUtils",
]
