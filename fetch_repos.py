import json
import os
import shutil
import subprocess
import concurrent.futures
import re
import stat
import time

# Configuration
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
SKILLS_FILE = os.path.join(BASE_DIR, "src", "data", "skills.json")
OUTPUT_DIR = os.path.join(BASE_DIR, "public", "data", "repos")
REPO_MAP_FILE = os.path.join(BASE_DIR, "src", "data", "repo_map.json")

# Limits
MAX_FILE_SIZE = 1 * 1024 * 1024  # 1MB
ALLOWED_EXTENSIONS = {
    # Code
    '.py', '.js', '.jsx', '.ts', '.tsx', '.html', '.css', '.scss', 
    '.java', '.c', '.cpp', '.h', '.cs', '.go', '.rs', '.php', '.rb', 
    '.swift', '.kt', '.scala', '.sh', '.bat', '.ps1', '.sql', '.r',
    # Config/Data
    '.json', '.yaml', '.yml', '.xml', '.toml', '.ini', '.conf', 
    '.dockerfile', 'dockerfile', '.gitignore', '.env.example',
    # Docs
    '.md', '.markdown', '.txt', '.rst',
    # Images (small ones)
    '.png', '.jpg', '.jpeg', '.gif', '.svg', '.ico'
}
SKIP_DIRS = {'.git', '.github', '.vscode', '.idea', 'node_modules', '__pycache__', 'dist', 'build', 'venv', 'env', '.docusaurus'}

def on_rm_error(func, path, exc_info):
    """
    Error handler for shutil.rmtree.
    If the error is due to an access error (read only file)
    it attempts to add write permission and then retries.
    """
    # Is the error an access error?
    if not os.access(path, os.W_OK):
        os.chmod(path, stat.S_IWRITE)
        func(path)
    else:
        raise

def load_skills():
    with open(SKILLS_FILE, 'r', encoding='utf-8') as f:
        return json.load(f)

def sanitize_name(name):
    # Sanitize more aggressively to ensure valid folder names
    return re.sub(r'[\\/*?:"<>|]', "", name).replace(" ", "_").lower().strip()

def get_dir_structure(root_dir, rel_path=""):
    """
    Recursively builds a JSON structure of the directory.
    """
    items = []
    try:
        # Sort directories first, then files
        entries = list(os.scandir(root_dir))
        entries.sort(key=lambda e: (not e.is_dir(), e.name.lower()))
        
        for entry in entries:
            if entry.name in SKIP_DIRS:
                continue
                
            full_path = os.path.join(root_dir, entry.name)
            item_rel_path = os.path.join(rel_path, entry.name).replace("\\", "/")
            
            if entry.is_dir():
                children = get_dir_structure(full_path, item_rel_path)
                if children: # Only add directories if they have content
                    items.append({
                        "name": entry.name,
                        "type": "folder",
                        "path": item_rel_path,
                        "children": children
                    })
            else:
                ext = os.path.splitext(entry.name)[1].lower()
                name_lower = entry.name.lower()
                
                is_allowed = (name_lower in ALLOWED_EXTENSIONS or ext in ALLOWED_EXTENSIONS)
                
                # Loose check for text files if extension is weird but small
                if not is_allowed and os.path.getsize(full_path) < 100 * 1024:
                     # Peek content? No, risky. 
                     # But some LICENSE files or makefiles don't have ext.
                     if "license" in name_lower or "makefile" in name_lower or "dockerfile" in name_lower:
                        is_allowed = True
                
                if is_allowed:
                    if os.path.getsize(full_path) <= MAX_FILE_SIZE:
                        items.append({
                            "name": entry.name,
                            "type": "file",
                            "path": item_rel_path,
                            "size": os.path.getsize(full_path),
                            "language": ext.replace('.', '') or 'text'
                        })
    except PermissionError:
        pass
        
    return items

def process_skill(skill):
    name = skill.get('name')
    url = skill.get('url')
    
    if not url or "github.com" not in url:
        return None

    sanitized = sanitize_name(name)
    target_dir = os.path.join(OUTPUT_DIR, sanitized)
    
    # Check if we already have a successful manifest
    manifest_path = os.path.join(target_dir, "manifest.json")
    if os.path.exists(manifest_path):
        return sanitized # Skip existing
        
    # Clean if partial
    if os.path.exists(target_dir):
        shutil.rmtree(target_dir, onerror=on_rm_error)
    
    os.makedirs(target_dir, exist_ok=True)
    
    try:
        temp_clone_dir = os.path.join(target_dir, "_temp_clone")
        
        # git clone --depth 1
        result = subprocess.run(
            ["git", "clone", "--depth", "1", url, temp_clone_dir], 
            stdout=subprocess.DEVNULL, 
            stderr=subprocess.DEVNULL,
            timeout=180 
        )
        
        if result.returncode != 0:
            shutil.rmtree(target_dir, onerror=on_rm_error) 
            return None

        # Clean .git using safe delete
        git_dir = os.path.join(temp_clone_dir, ".git")
        if os.path.exists(git_dir):
            shutil.rmtree(git_dir, onerror=on_rm_error)

        # Generate Tree
        tree = get_dir_structure(temp_clone_dir)
        
        # Save manifest
        with open(manifest_path, 'w', encoding='utf-8') as f:
            json.dump({"name": name, "tree": tree}, f, indent=2)
            
        # Move content
        content_dir = os.path.join(target_dir, "content")
        shutil.move(temp_clone_dir, content_dir)
        
        # Cleanup large files from content dir
        for root, dirs, files in os.walk(content_dir, topdown=False):
            for name in files:
                fp = os.path.join(root, name)
                try:
                    if os.path.getsize(fp) > MAX_FILE_SIZE:
                        os.remove(fp)
                except OSError:
                    pass
            
            for name in dirs:
                if name in SKIP_DIRS:
                    shutil.rmtree(os.path.join(root, name), onerror=on_rm_error)

        return sanitized
        
    except Exception as e:
        if os.path.exists(target_dir):
            try:
                shutil.rmtree(target_dir, onerror=on_rm_error) 
            except:
                pass
        return None

def main():
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
    
    skills_data = []
    with open(SKILLS_FILE, 'r', encoding='utf-8') as f:
        data = json.load(f)
        for cat in data:
            skills_data.extend(cat.get('skills', []))
            
    # Load existing map to update it
    repo_map = {}
    if os.path.exists(REPO_MAP_FILE):
        try:
            with open(REPO_MAP_FILE, 'r', encoding='utf-8') as f:
                repo_map = json.load(f)
        except:
            pass

    print(f"Starting download for {len(skills_data)} repositories...")
    
    # We use fewer workers to reduce IO/Network contention
    with concurrent.futures.ThreadPoolExecutor(max_workers=5) as executor:
        future_to_skill = {executor.submit(process_skill, skill): skill for skill in skills_data}
        
        completed_count = 0
        for future in concurrent.futures.as_completed(future_to_skill):
            skill = future_to_skill[future]
            try:
                result = future.result()
                if result:
                    repo_map[skill['name']] = result
                    completed_count += 1
                    print(f"[{completed_count}] Cached: {skill['name']}")
                else:
                    print(f"Failed: {skill['name']}")
            except Exception as exc:
                print(f"Exception for {skill['name']}: {exc}")

    # Save Mapping
    with open(REPO_MAP_FILE, 'w', encoding='utf-8') as f:
        json.dump(repo_map, f, indent=2)
    
    print(f"Done. Helper saved to {REPO_MAP_FILE}")

if __name__ == "__main__":
    main()
