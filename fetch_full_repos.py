"""
fetch_full_repos.py

Download complete repository files for all skills.
- Skips repos already fully downloaded (manifest has "full": true)
- Smart filtering: skip node_modules, dist, binaries, etc.
- Per-file cap: 300KB; per-repo cap: 10MB
- Uses git proxy via http://127.0.0.1:7890
- Supports subdirectory repos (github.com/owner/repo/tree/branch/path)
- Resumable: re-run safely, already-done repos are skipped
"""

import json
import os
import re
import shutil
import subprocess
import stat
import time
import concurrent.futures
from threading import Lock

# ─── Config ───────────────────────────────────────────────────────────────────
BASE_DIR    = os.path.dirname(os.path.abspath(__file__))
SKILLS_FILE = os.path.join(BASE_DIR, "src", "data", "skills.json")
OUTPUT_DIR  = os.path.join(BASE_DIR, "public", "data", "repos")
REPO_MAP_FILE = os.path.join(BASE_DIR, "src", "data", "repo_map.json")

MAX_FILE_SIZE  = 300 * 1024        # 300 KB per file
MAX_REPO_SIZE  = 5 * 1024 * 1024   # 5 MB total content per repo
MAX_FILES      = 400               # max files in manifest per repo
CLONE_TIMEOUT  = 180               # seconds
WORKERS        = 3                 # parallel git clone workers

# ─── Skip Lists ───────────────────────────────────────────────────────────────
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
    'poetry.lock', 'Pipfile.lock', 'Gemfile.lock', 'gemfile.lock',
    '.DS_Store', 'Thumbs.db', '.gitmodules', '.gitattributes',
    '.npmrc', '.yarnrc', '.nvmrc',
}

SKIP_EXTENSIONS = {
    # Archives
    '.zip', '.tar', '.gz', '.bz2', '.rar', '.7z', '.tgz', '.xz',
    # Executables / compiled
    '.exe', '.msi', '.dmg', '.pkg', '.deb', '.rpm', '.apk',
    '.bin', '.dll', '.so', '.dylib', '.a', '.lib', '.o', '.obj',
    # Bytecode
    '.class', '.pyc', '.pyo', '.pyd', '.whl',
    # Documents (serve as-is may not render well)
    '.pdf', '.doc', '.docx', '.xls', '.xlsx', '.ppt', '.pptx',
    # Media
    '.mp4', '.mp3', '.wav', '.avi', '.mov', '.mkv', '.flac', '.ogg',
    '.webm', '.m4a', '.aac',
    # Fonts
    '.ttf', '.otf', '.woff', '.woff2', '.eot',
    # ML / Data
    '.dat', '.npy', '.npz', '.h5', '.hdf5', '.pb', '.onnx',
    '.pth', '.pt', '.pkl', '.joblib', '.model', '.safetensors',
    # Lock / Generated
    '.lock',
}


# ─── Helpers ──────────────────────────────────────────────────────────────────

def on_rm_error(func, path, exc_info):
    """rmtree error handler: fix read-only files on Windows."""
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


def sanitize_name(name: str) -> str:
    return re.sub(r'[\\/*?:"<>|]', '', name).replace(' ', '_').lower().strip()


def parse_github_url(url: str):
    """
    Returns (clone_url, subpath).
    clone_url: e.g. https://github.com/owner/repo
    subpath:   None OR relative dir inside the repo, e.g. "skills/brand-guidelines"
    """
    url = url.rstrip('/')
    # Subdirectory: /tree/branch/path
    m = re.match(r'https://github\.com/([^/]+)/([^/]+)/tree/[^/]+/(.+)', url)
    if m:
        owner, repo, subpath = m.groups()
        return f'https://github.com/{owner}/{repo}', subpath.strip('/')
    # Top-level repo
    m = re.match(r'https://github\.com/([^/]+)/([^/]+)', url)
    if m:
        owner, repo = m.groups()
        return f'https://github.com/{owner}/{repo}', None
    return None, None


def needs_full_download(repo_dir: str) -> bool:
    """True if this repo directory doesn't have a complete download yet."""
    manifest_path = os.path.join(repo_dir, 'manifest.json')
    if not os.path.exists(manifest_path):
        return True
    try:
        with open(manifest_path, 'r', encoding='utf-8') as f:
            manifest = json.load(f)
        return not manifest.get('full', False)
    except Exception:
        return True


def get_dir_structure(root_dir: str, rel_path: str = '', size_tracker=None):
    """
    Recursively build JSON file-tree.
    size_tracker is a mutable list [current_total_bytes] for the per-repo cap.
    """
    if size_tracker is None:
        size_tracker = [0]

    items = []
    try:
        entries = sorted(
            os.scandir(root_dir),
            key=lambda e: (not e.is_dir(), e.name.lower())
        )
    except (PermissionError, OSError):
        return items

    for entry in entries:
        name = entry.name

        if entry.is_dir():
            if name in SKIP_DIRS or name.endswith('.egg-info'):
                continue
            child_path = f'{rel_path}/{name}' if rel_path else name
            children = get_dir_structure(entry.path, child_path, size_tracker)
            if children:
                items.append({
                    'name': name,
                    'type': 'folder',
                    'path': child_path,
                    'children': children,
                })
        else:
            if name in SKIP_FILES:
                continue
            # Skip minified JS/CSS
            if name.endswith('.min.js') or name.endswith('.min.css'):
                continue

            ext = os.path.splitext(name)[1].lower()
            if ext in SKIP_EXTENSIONS:
                continue

            try:
                size = entry.stat().st_size
            except OSError:
                continue

            if size > MAX_FILE_SIZE:
                continue

            # Per-repo size budget
            if size_tracker[0] + size > MAX_REPO_SIZE:
                continue

            # Per-repo file count budget
            if len(items) >= MAX_FILES:
                continue

            size_tracker[0] += size
            file_path = f'{rel_path}/{name}' if rel_path else name
            items.append({
                'name': name,
                'type': 'file',
                'path': file_path,
                'size': size,
                'language': ext.lstrip('.') or 'text',
            })

    return items


def _win_long(path: str) -> str:
    """Add \\?\ extended-length prefix for Windows long paths."""
    if os.name == 'nt':
        path = os.path.abspath(path)
        if not path.startswith('\\\\?\\'):
            return '\\\\?\\' + path
    return path


def _collect_manifest_paths(tree, result=None):
    if result is None:
        result = set()
    for node in tree:
        if node.get('type') == 'file':
            result.add(node['path'].replace('\\', '/'))
        elif node.get('type') == 'folder':
            _collect_manifest_paths(node.get('children', []), result)
    return result


def trim_content_to_manifest(content_dir: str, tree: list):
    """Delete files from content/ that are NOT listed in the manifest tree."""
    kept_paths = _collect_manifest_paths(tree)

    for root, dirs, files in os.walk(content_dir, topdown=False, onerror=lambda e: None):
        dirs[:] = [d for d in dirs if d not in SKIP_DIRS]
        for fname in files:
            try:
                fpath = os.path.join(root, fname)
                rel = os.path.relpath(fpath, content_dir).replace('\\', '/')
                if rel not in kept_paths:
                    try:
                        os.remove(_win_long(fpath))
                    except Exception:
                        pass
            except Exception:
                pass

    # Remove empty directories
    for root, dirs, files in os.walk(content_dir, topdown=False, onerror=lambda e: None):
        if root == content_dir:
            continue
        try:
            p = _win_long(root)
            if not os.listdir(p):
                os.rmdir(p)
        except Exception:
            pass




def clone_repo(clone_url: str, target_dir: str):
    """
    Clone repo to target_dir with --depth 1.
    Returns (success: bool, message: str).
    """
    env = os.environ.copy()
    # Inject proxy via git environment config
    env['GIT_CONFIG_COUNT'] = '2'
    env['GIT_CONFIG_KEY_0']   = 'http.proxy'
    env['GIT_CONFIG_VALUE_0'] = 'http://127.0.0.1:7890'
    env['GIT_CONFIG_KEY_1']   = 'https.proxy'
    env['GIT_CONFIG_VALUE_1'] = 'http://127.0.0.1:7890'

    try:
        result = subprocess.run(
            ['git', 'clone', '--depth', '1', clone_url, target_dir],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
            timeout=CLONE_TIMEOUT,
            env=env,
        )
        stderr = result.stderr.decode('utf-8', errors='replace').strip()
        return result.returncode == 0, stderr
    except subprocess.TimeoutExpired:
        return False, 'TimeoutExpired'
    except Exception as exc:
        return False, str(exc)


# ─── Main Processing ──────────────────────────────────────────────────────────

_repo_map_lock = Lock()


def process_skill(skill: dict, repo_map: dict):
    """
    Download a single skill's repo.
    Returns (skill_name, repo_id, status_string).
    """
    name = skill.get('name', '')
    url  = skill.get('url', '')

    if not url or 'github.com' not in url:
        return name, None, 'skip_no_url'

    repo_id  = repo_map.get(name) or sanitize_name(name)
    repo_dir = os.path.join(OUTPUT_DIR, repo_id)

    if not needs_full_download(repo_dir):
        return name, repo_id, 'already_done'

    clone_url, subpath = parse_github_url(url)
    if not clone_url:
        return name, repo_id, 'fail_invalid_url'

    os.makedirs(repo_dir, exist_ok=True)
    content_dir  = os.path.join(repo_dir, 'content')
    manifest_path = os.path.join(repo_dir, 'manifest.json')
    temp_dir     = os.path.join(repo_dir, '_clone_temp')

    # Cleanup leftover temp
    remove_dir(temp_dir)

    # ── Clone ──
    success, msg = clone_repo(clone_url, temp_dir)
    if not success:
        remove_dir(temp_dir)
        return name, repo_id, f'fail_clone: {msg[:150]}'

    # Remove .git directory
    remove_dir(os.path.join(temp_dir, '.git'))

    # Determine source: whole repo or a subpath
    if subpath:
        sub_dir = os.path.join(temp_dir, subpath.replace('/', os.sep))
        source_dir = sub_dir if os.path.isdir(sub_dir) else temp_dir
    else:
        source_dir = temp_dir

    # Build file tree manifest
    tree = get_dir_structure(source_dir)

    # Remove old content dir
    remove_dir(content_dir)

    # Move source → content/
    try:
        if source_dir == temp_dir:
            # Rename entire temp clone
            shutil.move(temp_dir, content_dir)
        else:
            # Copy subdirectory → content/, then delete temp
            shutil.copytree(source_dir, content_dir)
            remove_dir(temp_dir)
    except Exception as exc:
        remove_dir(temp_dir)
        return name, repo_id, f'fail_move: {exc}'

    # Write manifest with full=true marker
    manifest = {
        'name': name,
        'full': True,
        'tree': tree,
    }
    try:
        with open(manifest_path, 'w', encoding='utf-8') as f:
            json.dump(manifest, f, indent=2, ensure_ascii=False)
    except Exception as exc:
        return name, repo_id, f'fail_manifest: {exc}'

    # Trim content/ to only keep files listed in the manifest
    trim_content_to_manifest(content_dir, tree)

    return name, repo_id, 'done'


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Load skills
    with open(SKILLS_FILE, 'r', encoding='utf-8') as f:
        raw = json.load(f)
    all_skills = [s for cat in raw for s in cat.get('skills', [])]

    # Load existing repo map
    repo_map: dict = {}
    if os.path.exists(REPO_MAP_FILE):
        try:
            with open(REPO_MAP_FILE, 'r', encoding='utf-8') as f:
                repo_map = json.load(f)
        except Exception:
            pass

    total  = len(all_skills)
    done   = 0
    skip   = 0
    failed = 0
    counter = 0

    print(f"Processing {total} skills with {WORKERS} workers...\n")

    with concurrent.futures.ThreadPoolExecutor(max_workers=WORKERS) as executor:
        future_map = {
            executor.submit(process_skill, skill, repo_map): skill
            for skill in all_skills
        }

        for future in concurrent.futures.as_completed(future_map):
            counter += 1
            try:
                skill_name, repo_id, status = future.result()
            except Exception as exc:
                skill = future_map[future]
                skill_name = skill.get('name', '?')
                repo_id = None
                status  = f'exception: {exc}'

            # Update repo_map thread-safely
            if repo_id:
                with _repo_map_lock:
                    repo_map[skill_name] = repo_id

            tag = f'[{counter}/{total}]'
            if status == 'already_done':
                skip += 1
                print(f'{tag} ✓ (already full) {skill_name}')
            elif status == 'done':
                done += 1
                print(f'{tag} ✓ Downloaded: {skill_name}')
            elif status == 'skip_no_url':
                skip += 1
                print(f'{tag} – No URL: {skill_name}')
            else:
                failed += 1
                print(f'{tag} ✗ {status}: {skill_name}')

    # Persist updated repo map
    with open(REPO_MAP_FILE, 'w', encoding='utf-8') as f:
        json.dump(repo_map, f, indent=2, ensure_ascii=False)

    print(f'\n═══════════════════════════════════════')
    print(f'  Downloaded : {done}')
    print(f'  Skipped    : {skip}')
    print(f'  Failed     : {failed}')
    print(f'═══════════════════════════════════════')
    print(f'repo_map saved → {REPO_MAP_FILE}')


if __name__ == '__main__':
    main()
