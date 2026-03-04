"""
download_all_repos.py

Comprehensive script to:
1. Ensure ALL skills in skills.json have repo_map entries
2. Download content for all repos that are missing content
3. Uses ZIP download (fast) with proxy support
"""

import json
import os
import re
import io
import stat
import shutil
import zipfile
import time
import urllib.request
import concurrent.futures
from threading import Lock

BASE_DIR    = os.path.dirname(os.path.abspath(__file__))
SKILLS_FILE = os.path.join(BASE_DIR, "src", "data", "skills.json")
REPOS_DIR   = os.path.join(BASE_DIR, "public", "data", "repos")
REPO_MAP    = os.path.join(BASE_DIR, "src", "data", "repo_map.json")

MAX_FILE_SIZE  = 300 * 1024        # 300 KB per file
MAX_REPO_SIZE  = 5 * 1024 * 1024   # 5 MB total content per repo
DOWNLOAD_TIMEOUT = 120
WORKERS = 4

PROXY = "http://127.0.0.1:7890"

SKIP_DIRS = {
    '.git', '.github', '.vscode', '.idea',
    'node_modules', '__pycache__', 'dist', 'build', 'venv', 'env',
    '.next', '.nuxt', 'coverage', 'out', 'target', '.docusaurus',
    'vendor', 'site-packages', '.eggs', 'htmlcov', '.pytest_cache',
    '.mypy_cache', 'site', '_site', '.gradle', 'Pods', 'DerivedData',
    '.tox', 'buck-out', '.stack-work', 'elm-stuff',
}

SKIP_FILES = {
    'package-lock.json', 'yarn.lock', 'pnpm-lock.yaml', 'composer.lock',
    'poetry.lock', 'Pipfile.lock', 'Gemfile.lock', '.DS_Store', 'Thumbs.db',
}

SKIP_EXTENSIONS = {
    '.zip', '.tar', '.gz', '.bz2', '.rar', '.7z', '.tgz', '.xz',
    '.exe', '.msi', '.dmg', '.pkg', '.deb', '.rpm', '.apk',
    '.bin', '.dll', '.so', '.dylib', '.a', '.lib', '.o', '.obj',
    '.class', '.pyc', '.pyo', '.pyd', '.whl',
    '.pdf', '.doc', '.docx', '.xls', '.xlsx', '.ppt', '.pptx',
    '.mp4', '.mp3', '.wav', '.avi', '.mov', '.mkv', '.flac', '.ogg',
    '.webm', '.m4a', '.aac',
    '.ttf', '.otf', '.woff', '.woff2', '.eot',
    '.dat', '.npy', '.npz', '.h5', '.hdf5', '.pb', '.onnx',
    '.pth', '.pt', '.pkl', '.joblib', '.model', '.safetensors',
    '.lock',
}


def on_rm_error(func, path, exc_info):
    if not os.access(path, os.W_OK):
        os.chmod(path, stat.S_IWRITE)
        func(path)
    else:
        raise exc_info[1]


def remove_dir(path):
    if os.path.exists(path):
        try:
            shutil.rmtree(path, onerror=on_rm_error)
        except Exception:
            pass


def sanitize_name(name):
    return re.sub(r'[\\/*?:"<>|]', '', name).replace(' ', '_').replace('(', '').replace(')', '').lower().strip().rstrip('_')


def parse_github_url(url):
    """Returns (owner, repo, branch, subpath)."""
    url = url.rstrip('/')
    # Subdirectory: /tree/branch/path
    m = re.match(r'https://github\.com/([^/]+)/([^/]+)/tree/([^/]+)/(.+)', url)
    if m:
        return m.groups()
    # Top-level repo
    m = re.match(r'https://github\.com/([^/]+)/([^/]+)', url)
    if m:
        owner, repo = m.groups()
        return owner, repo, None, None
    return None, None, None, None


def download_zip(owner, repo, branch=None):
    """Download repo as zip, return bytes or None."""
    branches_to_try = [branch] if branch else ['main', 'master']
    
    proxy_handler = urllib.request.ProxyHandler({
        'http': PROXY, 'https': PROXY
    })
    opener = urllib.request.build_opener(proxy_handler)
    
    for b in branches_to_try:
        url = f"https://github.com/{owner}/{repo}/archive/refs/heads/{b}.zip"
        try:
            req = urllib.request.Request(url, headers={
                'User-Agent': 'Mozilla/5.0 (compatible; download_all/1.0)'
            })
            resp = opener.open(req, timeout=DOWNLOAD_TIMEOUT)
            return resp.read(), b
        except Exception as e:
            err = str(e)
            if '404' in err and b != branches_to_try[-1]:
                continue
            if '404' not in err:
                time.sleep(2)
                try:
                    resp = opener.open(req, timeout=DOWNLOAD_TIMEOUT)
                    return resp.read(), b
                except:
                    pass
    return None, None


def build_tree(root_dir):
    """Build file tree from directory, respecting size limits."""
    size_total = [0]
    
    def _walk(dir_path, rel_path):
        try:
            entries = sorted(os.scandir(dir_path), key=lambda e: (not e.is_dir(), e.name.lower()))
        except (PermissionError, OSError):
            return []
        
        result = []
        for entry in entries:
            name = entry.name
            if entry.is_dir():
                if name in SKIP_DIRS or (name.startswith('.') and name not in ('.devcontainer',)):
                    continue
                child_path = f'{rel_path}/{name}' if rel_path else name
                children = _walk(entry.path, child_path)
                if children:
                    result.append({
                        'name': name, 'type': 'folder',
                        'path': child_path, 'children': children,
                    })
            else:
                if name in SKIP_FILES:
                    continue
                if name.endswith('.min.js') or name.endswith('.min.css'):
                    continue
                ext = os.path.splitext(name)[1].lower()
                if ext in SKIP_EXTENSIONS:
                    continue
                try:
                    size = entry.stat().st_size
                except OSError:
                    continue
                if size > MAX_FILE_SIZE or size_total[0] + size > MAX_REPO_SIZE:
                    continue
                size_total[0] += size
                file_path = f'{rel_path}/{name}' if rel_path else name
                result.append({
                    'name': name, 'type': 'file',
                    'path': file_path, 'size': size,
                    'language': ext.lstrip('.') or 'text',
                })
        return result
    
    return _walk(root_dir, '')


def collect_tree_paths(tree):
    paths = set()
    for node in tree:
        if node['type'] == 'file':
            paths.add(node['path'])
        elif node['type'] == 'folder':
            paths.update(collect_tree_paths(node.get('children', [])))
    return paths


def trim_to_tree(content_dir, tree):
    kept = collect_tree_paths(tree)
    for root, dirs, files in os.walk(content_dir, topdown=False):
        for fname in files:
            fpath = os.path.join(root, fname)
            rel = os.path.relpath(fpath, content_dir).replace(os.sep, '/')
            if rel not in kept:
                try: os.remove(fpath)
                except: pass
        if root != content_dir:
            try:
                if not os.listdir(root): os.rmdir(root)
            except: pass


_lock = Lock()


def process_skill(skill, repo_id):
    """Download one skill's repo content. Returns (name, repo_id, status)."""
    name = skill.get('name', '')
    url = skill.get('url', '')
    
    if not url or 'github.com' not in url:
        return name, repo_id, 'skip_no_github'
    
    repo_dir = os.path.join(REPOS_DIR, repo_id)
    content_dir = os.path.join(repo_dir, 'content')
    manifest_path = os.path.join(repo_dir, 'manifest.json')
    
    # Check if already has content
    if os.path.exists(content_dir):
        count = sum(1 for _, _, fs in os.walk(content_dir) for f in fs)
        if count > 0:
            return name, repo_id, 'already_ok'
    
    # Parse URL
    owner, repo, branch, subpath = parse_github_url(url)
    if not owner or not repo:
        return name, repo_id, 'skip_bad_url'
    
    # Download zip
    zip_data, used_branch = download_zip(owner, repo, branch)
    if not zip_data:
        return name, repo_id, 'fail_download'
    
    # Extract
    remove_dir(content_dir)
    os.makedirs(content_dir, exist_ok=True)
    
    try:
        zf = zipfile.ZipFile(io.BytesIO(zip_data))
    except zipfile.BadZipFile:
        return name, repo_id, 'fail_bad_zip'
    
    members = zf.namelist()
    if not members:
        return name, repo_id, 'fail_empty_zip'
    
    # Find the common prefix (e.g., "repo-main/")
    prefix = members[0]
    if '/' in prefix:
        prefix = prefix.split('/')[0] + '/'
    
    # If subpath, adjust prefix
    if subpath:
        prefix = prefix + subpath.rstrip('/') + '/'
    
    extracted = 0
    total_size = 0
    for member in members:
        if not member.startswith(prefix):
            continue
        rel = member[len(prefix):]
        if not rel or rel.endswith('/'):
            continue
        
        parts = rel.split('/')
        if any(p in SKIP_DIRS for p in parts[:-1]):
            continue
        if any(p.startswith('.') and p not in ('.devcontainer',) for p in parts[:-1]):
            continue
        fname = parts[-1]
        if fname in SKIP_FILES:
            continue
        if fname.endswith('.min.js') or fname.endswith('.min.css'):
            continue
        ext = os.path.splitext(fname)[1].lower()
        if ext in SKIP_EXTENSIONS:
            continue
        
        info = zf.getinfo(member)
        if info.file_size > MAX_FILE_SIZE:
            continue
        if total_size + info.file_size > MAX_REPO_SIZE:
            continue
        
        target = os.path.join(content_dir, rel.replace('/', os.sep))
        os.makedirs(os.path.dirname(target), exist_ok=True)
        try:
            with zf.open(member) as src, open(target, 'wb') as dst:
                data = src.read()
                dst.write(data)
            total_size += len(data)
            extracted += 1
        except Exception:
            pass
    
    zf.close()
    
    if extracted == 0:
        # Try to at least get README
        return name, repo_id, 'fail_no_files'
    
    # Build tree & write manifest
    tree = build_tree(content_dir)
    trim_to_tree(content_dir, tree)
    
    manifest = {'name': name, 'full': True, 'tree': tree}
    os.makedirs(repo_dir, exist_ok=True)
    with open(manifest_path, 'w', encoding='utf-8') as f:
        json.dump(manifest, f, indent=2, ensure_ascii=False)
    
    return name, repo_id, f'ok ({extracted} files)'


def main():
    # Load skills
    with open(SKILLS_FILE, 'r', encoding='utf-8') as f:
        raw = json.load(f)
    all_skills = [s for cat in raw for s in cat.get('skills', [])]
    
    # Load or create repo map
    repo_map = {}
    if os.path.exists(REPO_MAP):
        with open(REPO_MAP, 'r', encoding='utf-8') as f:
            repo_map = json.load(f)
    
    # Ensure ALL skills have repo_map entries
    new_mappings = 0
    for skill in all_skills:
        name = skill['name']
        if name not in repo_map:
            repo_id = sanitize_name(name)
            # Avoid collisions
            if repo_id in repo_map.values():
                repo_id = repo_id + '_' + str(hash(name) % 10000)
            repo_map[name] = repo_id
            new_mappings += 1
    
    if new_mappings:
        print(f"Added {new_mappings} new repo_map entries")
        with open(REPO_MAP, 'w', encoding='utf-8') as f:
            json.dump(repo_map, f, indent=2, ensure_ascii=False)
    
    # Find skills that need downloading
    needs_download = []
    for skill in all_skills:
        name = skill['name']
        repo_id = repo_map[name]
        content_dir = os.path.join(REPOS_DIR, repo_id, 'content')
        
        if os.path.exists(content_dir):
            count = sum(1 for _, _, fs in os.walk(content_dir) for f in fs)
            if count > 0:
                continue
        
        needs_download.append((skill, repo_id))
    
    print(f"Total skills: {len(all_skills)}")
    print(f"Need download: {len(needs_download)}")
    print()
    
    if not needs_download:
        print("All repos have content! Nothing to do.")
        return
    
    ok = failed = skipped = 0
    total = len(needs_download)
    
    with concurrent.futures.ThreadPoolExecutor(max_workers=WORKERS) as executor:
        futures = {
            executor.submit(process_skill, skill, repo_id): (skill, repo_id)
            for skill, repo_id in needs_download
        }
        
        for i, future in enumerate(concurrent.futures.as_completed(futures), 1):
            try:
                name, rid, status = future.result()
            except Exception as e:
                skill, rid = futures[future]
                name = skill.get('name', '?')
                status = f'exception: {e}'
            
            tag = f'[{i}/{total}]'
            if status.startswith('ok'):
                ok += 1
                print(f'{tag} OK {name} -> {rid} - {status}')
            elif status.startswith('skip') or status == 'already_ok':
                skipped += 1
                print(f'{tag} SKIP {name} - {status}')
            else:
                failed += 1
                print(f'{tag} FAIL {name} - {status}')
    
    # Save updated repo_map
    with open(REPO_MAP, 'w', encoding='utf-8') as f:
        json.dump(repo_map, f, indent=2, ensure_ascii=False)
    
    print(f'\n{"="*50}')
    print(f'  Downloaded: {ok}')
    print(f'  Skipped:    {skipped}')
    print(f'  Failed:     {failed}')
    print(f'{"="*50}')


if __name__ == '__main__':
    main()
