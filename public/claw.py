import argparse
import json
import os
import shutil
import sys
import urllib.request
import urllib.error

# Configuration
# Ideally, this would be a remote URL like "https://xdelin.github.io/awesome-agent4edu/registry.json"
# But for this local demo, we point to the file we just generated.
REGISTRY_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "skills-web-interface", "public", "registry.json")
LOCAL_REPO_BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "skills-web-interface", "public")
REMOTE_BASE_URL = "https://xdelin.github.io/awesome-agent4edu/"

def load_registry():
    # Try local first (dev mode)
    if os.path.exists(REGISTRY_PATH):
        try:
            with open(REGISTRY_PATH, 'r', encoding='utf-8') as f:
                return json.load(f), LOCAL_REPO_BASE
        except Exception as e:
            print(f"Error loading local registry: {e}")
    
    # Try remote
    registry_url = REMOTE_BASE_URL + "registry.json"
    print(f"Loading registry from {registry_url}...")
    try:
        with urllib.request.urlopen(registry_url) as response:
            return json.loads(response.read().decode('utf-8')), REMOTE_BASE_URL
    except Exception as e:
        print(f"Error loading remote registry: {e}")
        return {}, None

def install_skill(skill_name, target_dir="skills"):
    registry, base_url = load_registry()
    if not registry:
        return

    # Case-insensitive search
    skill_key = None
    for key in registry:
        if key.lower() == skill_name.lower():
            skill_key = key
            break
    
    if not skill_key:
        print(f"Skill '{skill_name}' not found in registry.")
        return

    skill_info = registry[skill_key]
    print(f"Found skill: {skill_key}")
    
    # Check if we have a direct zip path
    zip_rel_path = skill_info.get("path")
    if zip_rel_path and zip_rel_path.endswith('.zip'):
        # Determine install location
        install_path = os.path.abspath(os.path.join(target_dir, skill_name))
        
        # Determine source
        if base_url.startswith('http'):
            # Remote download
            zip_url = base_url + zip_rel_path
            print(f"Downloading {zip_url}...")
            try:
                temp_zip = f"{skill_name}.zip"
                urllib.request.urlretrieve(zip_url, temp_zip)
                
                os.makedirs(install_path, exist_ok=True)
                shutil.unpack_archive(temp_zip, install_path)
                os.remove(temp_zip)
                print(f"Successfully installed {skill_name}!")
            except Exception as e:
                 print(f"Error downloading/installing {skill_name}: {e}")
        else:
            # Local copy
            local_zip_path = os.path.join(base_url, zip_rel_path)
            os.makedirs(install_path, exist_ok=True)
            try:
                shutil.unpack_archive(local_zip_path, install_path)
                print(f"Successfully installed {skill_name}!")
            except Exception as e:
                print(f"Error installing {skill_name}: {e}")
        return


    if skill_info.get("type") == "hosted":
        # Install from manifest (Old logic)
        author = skill_info.get("author", "unknown")
        install_path = os.path.join(target_dir, author, skill_name)
        
        if os.path.exists(install_path):
            print(f"Skill already exists at {install_path}")
            # Ask to overwrite? For now, just warn.
            ans = input("Overwrite? (y/N): ")
            if ans.lower() != 'y':
                return

        # Fetch manifest
        manifest_rel_path = skill_info["manifest_url"]
        content_rel_base = skill_info["content_base_url"]
        
        manifest_path = os.path.join(LOCAL_REPO_BASE, manifest_rel_path)
        
        if not os.path.exists(manifest_path):
            print(f"Error: Manifest not found at {manifest_path}")
            return
            
        with open(manifest_path, 'r', encoding='utf-8') as f:
            manifest = json.load(f)
            
        # Recursive download
        file_count = 0
        
        def process_node(node, current_rel_path=""):
            nonlocal file_count
            if node["type"] == "file":
                src_file_path = os.path.join(LOCAL_REPO_BASE, content_rel_base, node["path"])
                dest_file_path = os.path.join(install_path, node["path"])
                
                os.makedirs(os.path.dirname(dest_file_path), exist_ok=True)
                
                try:
                    shutil.copy2(src_file_path, dest_file_path)
                    print(f"  + {node['path']}")
                    file_count += 1
                except Exception as e:
                    print(f"  ! Error copying {node['path']}: {e}")
                    
            elif node["type"] == "folder":
                for child in node["children"]:
                    process_node(child, os.path.join(current_rel_path, node["name"]))

        print(f"Installing {skill_key} to {install_path}...")
        for node in manifest["tree"]:
            process_node(node)
            
        print(f"\nSuccessfully installed {file_count} files.")
        print(f"Skill is ready at: {install_path}")

    else:
        print(f"This skill is hosted externally at: {skill_info.get('url')}")
        print("Automatic installation for external repos is not yet supported in this demo.")
        print("Please clone manually.")

def list_skills():
    registry, _ = load_registry()
    if not registry:
        return
    
    print(f"{'Skill Name':<35} | {'Type':<10} | {'Description'}")
    print("-" * 80)
    for name, info in registry.items():
        desc = info.get("description", "") or ""
        desc = desc[:50] + "..." if len(desc) > 50 else desc
        print(f"{name:<35} | {info.get('type'):<10} | {desc.replace(chr(10), ' ')}")

def main():
    parser = argparse.ArgumentParser(description="ClawHub CLI - Skill Manager for OpenClaw")
    subparsers = parser.add_subparsers(dest="command", help="Command to run")
    
    # Install command
    install_parser = subparsers.add_parser("install", help="Install a skill")
    install_parser.add_argument("skill", help="Name of the skill to install")
    install_parser.add_argument("--dir", default="skills", help="Target directory (default: ./skills)")

    # List command
    subparsers.add_parser("list", help="List available skills")

    args = parser.parse_args()
    
    if args.command == "install":
        install_skill(args.skill, args.dir)
    elif args.command == "list":
        list_skills()
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
