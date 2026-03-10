"""
Trim repo cache to only keep README files.
Updates manifest.json to only show README entries,
and deletes all non-README files from content/.
"""
import json
import os
import shutil
import stat

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
REPOS_DIR = os.path.join(BASE_DIR, "public", "data", "repos")

def on_rm_error(func, path, exc_info):
    if not os.access(path, os.W_OK):
        os.chmod(path, stat.S_IWRITE)
        func(path)
    else:
        raise

def is_readme(name):
    """Check if a filename is a README variant."""
    lower = name.lower()
    return lower.startswith('readme') or lower == 'skill.md'

def filter_tree_for_readmes(tree):
    """Filter manifest tree to only include README files (at root level only)."""
    result = []
    for node in tree:
        if node['type'] == 'file' and is_readme(node['name']):
            result.append(node)
    return result

def process_repo(repo_dir):
    """Process a single repo directory: keep only README files."""
    manifest_path = os.path.join(repo_dir, "manifest.json")
    content_dir = os.path.join(repo_dir, "content")
    
    if not os.path.exists(manifest_path):
        return False
    
    # Load manifest
    with open(manifest_path, 'r', encoding='utf-8') as f:
        manifest = json.load(f)
    
    # Filter tree to only README files at root
    original_tree = manifest.get('tree', [])
    filtered_tree = filter_tree_for_readmes(original_tree)
    
    # Collect paths of README files to keep
    keep_files = set()
    for node in filtered_tree:
        keep_files.add(node['path'])
    
    # Delete all files in content/ except READMEs at root
    if os.path.exists(content_dir):
        for item in os.listdir(content_dir):
            item_path = os.path.join(content_dir, item)
            if item in keep_files:
                continue  # Keep this README
            try:
                if os.path.isdir(item_path):
                    shutil.rmtree(item_path, onerror=on_rm_error)
                else:
                    os.remove(item_path)
            except Exception as e:
                print(f"  Warning: Could not delete {item}: {e}")
    
    # Update manifest
    manifest['tree'] = filtered_tree
    with open(manifest_path, 'w', encoding='utf-8') as f:
        json.dump(manifest, f, indent=2, ensure_ascii=False)
    
    return len(filtered_tree) > 0

def main():
    if not os.path.exists(REPOS_DIR):
        print("No repos directory found!")
        return
    
    repos = [d for d in os.listdir(REPOS_DIR) 
             if os.path.isdir(os.path.join(REPOS_DIR, d))]
    
    print(f"Processing {len(repos)} repos...")
    
    has_readme = 0
    no_readme = 0
    
    for i, repo_name in enumerate(sorted(repos)):
        repo_path = os.path.join(REPOS_DIR, repo_name)
        result = process_repo(repo_path)
        if result:
            has_readme += 1
        else:
            no_readme += 1
        
        if (i + 1) % 20 == 0:
            print(f"  Processed {i+1}/{len(repos)}...")
    
    print(f"\nDone! {has_readme} repos have README, {no_readme} repos have no README.")
    
    # Calculate final size
    total_size = 0
    for root, dirs, files in os.walk(REPOS_DIR):
        for f in files:
            fp = os.path.join(root, f)
            total_size += os.path.getsize(fp)
    print(f"Total size: {total_size / 1024 / 1024:.2f} MB")

if __name__ == '__main__':
    main()
