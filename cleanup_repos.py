"""
cleanup_repos.py

Trim each repo's content/ directory to ONLY contain files listed in manifest.json.
This removes the bulk of cloned files that exceed our size budget.

Also:
- Removes leftover _clone_temp directories
- Uses Windows \\?\ extended path prefix for long-path support
- Deletes empty directories after file removal
"""

import os
import json
import shutil
import stat

BASE_DIR  = os.path.dirname(os.path.abspath(__file__))
REPOS_DIR = os.path.join(BASE_DIR, "public", "data", "repos")

# Repos that are simply too huge to serve — delete content/ entirely
# (their manifest is kept, so the UI can still show the file tree but files
# will gracefully fail to load and fall back to the GitHub link)
SKIP_FULL_REPOS = {'openfoam'}  # openfoam/openfoam-dev is millions of files


def on_rm_error(func, path, exc_info):
    if not os.access(path, os.W_OK):
        os.chmod(path, stat.S_IWRITE)
        func(path)
    else:
        raise exc_info[1]


def win_long(path: str) -> str:
    """Add \\?\ extended-length prefix for Windows paths > 260 chars."""
    if os.name == 'nt':
        path = os.path.abspath(path)
        if not path.startswith('\\\\?\\'):
            return '\\\\?\\' + path
    return path


def collect_manifest_paths(tree, result=None):
    """Recursively collect all file paths listed in the manifest tree."""
    if result is None:
        result = set()
    for node in tree:
        if node.get('type') == 'file':
            result.add(node['path'].replace('\\', '/'))
        elif node.get('type') == 'folder':
            collect_manifest_paths(node.get('children', []), result)
    return result


def safe_remove_file(path: str) -> bool:
    """Remove a file, with Windows long-path support."""
    try:
        os.remove(win_long(path))
        return True
    except Exception:
        return False


def safe_remove_dir(path: str):
    """Remove a directory tree, with Windows long-path support."""
    try:
        # Try \\?\ prefix first
        p = win_long(path)
        shutil.rmtree(p, onerror=on_rm_error)
    except Exception:
        pass


def trim_repo(repo_dir: str) -> tuple:
    """
    Trim content/ to only keep files listed in manifest.json.
    Returns (kept, deleted, errors) counts.
    """
    repo_name = os.path.basename(repo_dir)
    manifest_path = os.path.join(repo_dir, 'manifest.json')
    content_dir   = os.path.join(repo_dir, 'content')

    if not os.path.exists(manifest_path):
        return 0, 0, 0

    # Load manifest
    try:
        with open(manifest_path, 'r', encoding='utf-8') as f:
            manifest = json.load(f)
    except Exception as e:
        print(f"  [{repo_name}] ERROR reading manifest: {e}")
        return 0, 0, 1

    # For repos that are too huge: delete entire content/
    if repo_name.lower() in SKIP_FULL_REPOS:
        if os.path.exists(content_dir):
            print(f"  [{repo_name}] Deleting entire content/ (too large)...")
            safe_remove_dir(content_dir)
        return 0, -1, 0

    if not manifest.get('full'):
        return 0, 0, 0  # Not a full download, skip

    if not os.path.exists(content_dir):
        return 0, 0, 0

    # Collect file paths that should be kept
    kept_paths = collect_manifest_paths(manifest.get('tree', []))

    kept = deleted = errors = 0

    # Walk content/ using os.scandir for better long-path support
    for root, dirs, files in os.walk(content_dir, topdown=False, onerror=lambda e: None):
        # Skip known skip dirs (shouldn't be in content/ after clone, but just in case)
        skip_names = {'.git', 'node_modules', '__pycache__'}
        dirs[:] = [d for d in dirs if d not in skip_names]

        for fname in files:
            try:
                fpath = os.path.join(root, fname)
                rel   = os.path.relpath(fpath, content_dir).replace('\\', '/')
                if rel in kept_paths:
                    kept += 1
                else:
                    if safe_remove_file(fpath):
                        deleted += 1
                    else:
                        errors += 1
            except Exception:
                errors += 1

    # Remove empty directories (skip content_dir itself)
    for root, dirs, files in os.walk(content_dir, topdown=False, onerror=lambda e: None):
        if root == content_dir:
            continue
        try:
            p = win_long(root)
            if not os.listdir(p):
                os.rmdir(p)
        except Exception:
            pass

    return kept, deleted, errors


def cleanup_temp_dirs():
    """Remove leftover _clone_temp directories."""
    count = 0
    for repo_name in os.listdir(REPOS_DIR):
        temp_dir = os.path.join(REPOS_DIR, repo_name, '_clone_temp')
        if os.path.exists(temp_dir):
            print(f"Removing temp: {repo_name}/_clone_temp")
            safe_remove_dir(temp_dir)
            count += 1
    if count:
        print(f"Removed {count} temp dir(s).\n")
    else:
        print("No temp dirs found.\n")


def main():
    print(f"Repos dir: {REPOS_DIR}\n")

    # Step 1: clean up temp dirs
    print("=== Step 1: Clean up _clone_temp leftovers ===")
    cleanup_temp_dirs()

    # Step 2: trim content/ directories
    print("=== Step 2: Trim content/ to match manifest ===")
    repo_dirs = sorted([
        os.path.join(REPOS_DIR, d)
        for d in os.listdir(REPOS_DIR)
        if os.path.isdir(os.path.join(REPOS_DIR, d))
    ])

    total_kept = total_deleted = total_errors = 0
    for repo_dir in repo_dirs:
        repo_name = os.path.basename(repo_dir)
        kept, deleted, errors = trim_repo(repo_dir)
        if deleted > 0 or errors > 0 or deleted == -1:
            if deleted == -1:
                print(f"  [{repo_name}] Deleted entire content/ (too large)")
            else:
                print(f"  [{repo_name}] kept={kept}, deleted={deleted}, errors={errors}")
        total_kept    += max(kept, 0)
        total_deleted += max(deleted, 0)
        total_errors  += errors

    print(f"\n=== Done ===")
    print(f"  Total kept   : {total_kept}")
    print(f"  Total deleted: {total_deleted}")
    print(f"  Total errors : {total_errors}")


if __name__ == '__main__':
    main()
