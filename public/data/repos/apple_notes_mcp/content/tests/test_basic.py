"""Basic tests for mcp-apple-notes package."""

import pytest
from mcp_apple_notes import __version__


def test_version():
    """Test that version is defined."""
    assert __version__ is not None
    assert isinstance(__version__, str)


def test_import_server():
    """Test that server module can be imported."""
    from mcp_apple_notes.server import mcp
    assert mcp is not None


def test_import_tools():
    """Test that tools module can be imported."""
    from mcp_apple_notes.tools import notes_tools
    assert notes_tools is not None


def test_import_applescript():
    """Test that applescript module can be imported."""
    from mcp_apple_notes.applescript import create_note
    assert create_note is not None


def test_package_structure():
    """Test that all required modules exist."""
    import mcp_apple_notes
    
    # Check that main modules exist
    assert hasattr(mcp_apple_notes, 'server')
    assert hasattr(mcp_apple_notes, 'tools')
    assert hasattr(mcp_apple_notes, 'applescript')
