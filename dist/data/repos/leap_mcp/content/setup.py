#!/usr/bin/env python3
"""
LEAP MCP Setup for Claude Code
"""

import os
import json
import sys
from pathlib import Path
import subprocess

def check_python_version():
    """Check Python version compatibility"""
    version = sys.version_info
    if version.major != 3 or version.minor < 8:
        print("ERROR: Python 3.8 or higher is required")
        print(f"You have: Python {version.major}.{version.minor}.{version.micro}")
        sys.exit(1)
    elif version.minor >= 13:
        print("WARNING: Python 3.13+ detected")
        print("Manim has known compatibility issues with Python 3.13")
        print("We strongly recommend using Python 3.8-3.11")
        print()
        response = input("Continue anyway? (not recommended) (y/n): ").strip().lower()
        if response != 'y':
            print("\nTo use a compatible Python version:")
            print("  - Use pyenv: pyenv install 3.11.9 && pyenv local 3.11.9")
            print("  - Or conda: conda create -n leap python=3.11")
            sys.exit(1)
    else:
        print(f"[OK] Python {version.major}.{version.minor} detected - compatible")

def setup_env_file():
    """Setup .env file with OpenAI API key"""
    print("\n=== Setting up API Key ===")
    
    env_file = Path(".env")
    env_example = Path(".env.example")
    
    # Create .env.example
    with open(env_example, 'w') as f:
        f.write("# OpenAI API Configuration\n")
        f.write("OPENAI_API_KEY=your_openai_api_key_here\n")
    
    if env_file.exists():
        print("[OK] .env file already exists")
        # Check if API key is set
        with open(env_file, 'r') as f:
            content = f.read()
            if 'OPENAI_API_KEY=' in content and 'your_openai_api_key_here' not in content:
                print("[OK] OpenAI API key appears to be configured")
                return True
    
    print("\nYou need an OpenAI API key for text-to-speech narration")
    print("Get one at: https://platform.openai.com/api-keys")
    print()
    
    api_key = input("Enter your OpenAI API key (or press Enter to add it manually later): ").strip()
    
    with open(env_file, 'w') as f:
        f.write("# OpenAI API Configuration\n")
        if api_key and api_key.startswith('sk-'):
            f.write(f"OPENAI_API_KEY={api_key}\n")
            print("[OK] API key saved to .env")
            return True
        else:
            f.write("OPENAI_API_KEY=your_openai_api_key_here\n")
            print("[!] Add your API key to .env before generating videos")
            return False

def setup_mcp_config():
    """Create .mcp.json configuration for Claude Code"""
    print("\n=== Configuring Claude Code ===")
    
    project_dir = Path.cwd()
    python_path = sys.executable
    
    config = {
        "mcpServers": {
            "leap-mcp": {
                "command": str(python_path),
                "args": [str(project_dir / "leap_mcp.py")]
            }
        }
    }
    
    config_file = Path(".mcp.json")
    
    # Backup existing config if present
    if config_file.exists():
        backup_file = Path(".mcp.json.backup")
        import shutil
        shutil.copy(config_file, backup_file)
        print(f"[OK] Backed up existing config to {backup_file}")
    
    with open(config_file, 'w') as f:
        json.dump(config, f, indent=2)
    
    print(f"[OK] Created .mcp.json for Claude Code")
    print(f"    Using Python: {python_path}")

def install_dependencies():
    """Install Python dependencies"""
    print("\n=== Installing Dependencies ===")
    
    requirements_file = Path("requirements.txt")
    if not requirements_file.exists():
        print("ERROR: requirements.txt not found")
        return False
    
    print("Installing packages...")
    try:
        result = subprocess.run(
            [sys.executable, "-m", "pip", "install", "-r", "requirements.txt"],
            capture_output=True,
            text=True
        )
        if result.returncode == 0:
            print("[OK] All dependencies installed")
            return True
        else:
            print("ERROR: Failed to install dependencies")
            print(result.stderr)
            return False
    except Exception as e:
        print(f"ERROR: {e}")
        return False

def main():
    """Run setup for Claude Code"""
    print("LEAP MCP Setup for Claude Code")
    print("=" * 40)
    
    # Check Python version
    check_python_version()
    
    # Install dependencies
    if not install_dependencies():
        print("\n[!] Fix dependency issues and run setup again")
        sys.exit(1)
    
    # Setup environment file
    api_key_configured = setup_env_file()
    
    # Create MCP configuration
    setup_mcp_config()
    
    # Final instructions
    print("\n" + "=" * 40)
    print("Setup Complete!")
    print("\nNext steps:")
    
    if not api_key_configured:
        print("1. Add your OpenAI API key to .env file")
        print("2. Restart Claude Code")
        print("3. Run /mcp in chat to connect")
    else:
        print("1. Restart Claude Code") 
        print("2. Run /mcp in chat to connect")
    
    print("\nTest with: 'Create an educational video about gravity'")
    print("\nIf connection fails, run /doctor for diagnostics")

if __name__ == "__main__":
    main()