import os
import json
import shutil
import sys

# Configuration
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.path.join(BASE_DIR, "public", "data", "repos")
REPO_MAP_FILE = os.path.join(BASE_DIR, "src", "data", "repo_map.json")
MAX_FILE_SIZE = 1 * 1024 * 1024  # 1MB
SKIP_DIRS = {'.git', '.github', '.vscode', '.idea', 'node_modules', '__pycache__', 'dist', 'build'}

def get_dir_structure(root_dir):
    """
    Recursively scans the directory to build a tree structure for the manifest.
    """
    tree = []
    try:
        items = os.listdir(root_dir)
        # Sort: folders first, then files
        items.sort(key=lambda x: (not os.path.isdir(os.path.join(root_dir, x)), x.lower()))
    except PermissionError:
        return []
        
    for item in items:
        if item in SKIP_DIRS: 
            continue
            
        full_path = os.path.join(root_dir, item)
        rel_path = os.path.relpath(full_path, root_dir).replace("\\", "/")
        
        if os.path.isdir(full_path):
            children = get_dir_structure(full_path)
            if children: # Only add non-empty folders
                tree.append({
                    "name": item,
                    "type": "folder",
                    "path": rel_path,
                    "children": children
                })
        else:
            # File
            try:
                size = os.path.getsize(full_path)
                if size <= MAX_FILE_SIZE:
                    ext = os.path.splitext(item)[1].lower()
                    # Simple heuristic for language
                    tree.append({
                        "name": item,
                        "type": "file",
                        "path": rel_path,
                        "size": size,
                        "language": ext.replace('.', '') or 'text'
                    })
            except OSError:
                pass
                
    return tree

def import_repo(source_path, skill_name):
    if not os.path.exists(source_path):
        print(f"Error: Source path '{source_path}' does not exist.")
        return

    print(f"Importing '{skill_name}' from '{source_path}'...")

    # update repo map
    repo_id = skill_name.replace(" ", "_").replace("/", "_")
    target_dir = os.path.join(OUTPUT_DIR, repo_id)
    
    # Clean output
    if os.path.exists(target_dir):
        shutil.rmtree(target_dir)
    os.makedirs(target_dir)
    
    content_dir = os.path.join(target_dir, "content")
    os.makedirs(content_dir)

    # 1. Copy files
    # We use shutil.copytree but we want to skip ignored dirs and large files
    # It's easier to copy specific valid files or copy all then clean.
    # Let's copy all then clean for simplicity.
    
    print("  Copying files...")
    def ignore_patterns(path, names):
        return [n for n in names if n in SKIP_DIRS]

    try:
        shutil.copytree(source_path, content_dir, dirs_exist_ok=True, ignore=ignore_patterns)
    except Exception as e:
        print(f"  Warning during copy: {e}")

    # 2. Clean large files
    print("  Cleaning large files...")
    for root, dirs, files in os.walk(content_dir):
        for name in files:
            fp = os.path.join(root, name)
            try:
                if os.path.getsize(fp) > MAX_FILE_SIZE:
                    os.remove(fp)
            except:
                pass

    # 3. Generate Manifest
    print("  Generating manifest...")
    tree = get_dir_structure(content_dir)
    manifest_path = os.path.join(target_dir, "manifest.json")
    with open(manifest_path, 'w', encoding='utf-8') as f:
        json.dump({"name": skill_name, "tree": tree}, f, indent=2)

    # 4. Update Repo Map
    print("  Updating database...")
    if os.path.exists(REPO_MAP_FILE):
        with open(REPO_MAP_FILE, 'r', encoding='utf-8') as f:
            repo_map = json.load(f)
    else:
        repo_map = {}
    
    repo_map[skill_name] = repo_id
    
    with open(REPO_MAP_FILE, 'w', encoding='utf-8') as f:
        json.dump(repo_map, f, indent=2)

    print(f"Done! '{skill_name}' is now available in the viewer.")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python import_local_repo.py <path_to_local_folder> <perfect_skill_name_match>")
        print("Example: python import_local_repo.py C:/Projects/MySkill 'My Skill Name'")
    else:
        import_repo(sys.argv[1], sys.argv[2])
