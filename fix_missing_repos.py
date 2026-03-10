import os
import json
import shutil
import urllib.request
import zipfile
import io
import time

# Configurations
PROXY = "http://127.0.0.1:7890"

# Setup Proxy - REMOVED DEFAULT
# os.environ["http_proxy"] = PROXY
# os.environ["https_proxy"] = PROXY

def get_repo_id(name, repo_map):
    return repo_map.get(name)

def copy_local_skill(skill_name, local_subpath, repo_map):
    rid = get_repo_id(skill_name, repo_map)
    if not rid:
        print(f"Skipping {skill_name}: No ID found")
        return

    # Source: absolute path to local skills_demo
    src_base = os.path.abspath(os.path.join("..", "skills_demo", "anthropics-skills", "skills"))
    src_path = os.path.join(src_base, local_subpath)
    
    # Dest: public/data/repos/{id}/content
    dest_path = os.path.abspath(os.path.join("public", "data", "repos", rid, "content"))
    
    if os.path.exists(dest_path) and os.listdir(dest_path):
        print(f"Skipping copy for {skill_name}: Content already exists")
        return
    
    if os.path.exists(src_path):
        # Ignore node_modules and other junk
        shutil.copytree(src_path, dest_path, ignore=shutil.ignore_patterns('node_modules', '.git', '*.pyc', '__pycache__'))
        print(f"✅ Copied {skill_name} from local: {src_path}")
    else:
        print(f"❌ Local source not found for {skill_name}: {src_path}")

def create_placeholder(skill_name, url, repo_map):
    rid = get_repo_id(skill_name, repo_map)
    if not rid: return

    dest_path = os.path.join("public", "data", "repos", rid, "content")
    os.makedirs(dest_path, exist_ok=True)
    
    readme_path = os.path.join(dest_path, "README.md")
    content = f"# {skill_name}\n\nThis resource is a web service or proprietary platform. It does not have a public code repository.\n\nYou can access it here: [{url}]({url})\n"
    
    with open(readme_path, "w", encoding="utf-8") as f:
        f.write(content)
    print(f"📝 Created placeholder for {skill_name}")

def attempt_download_with_retry(url, dest_path, subdir, proxies_list):
    for proxy in proxies_list:
        try:
            if proxy:
                print(f"   Downloading zip (Proxy: {proxy})...")
                os.environ["http_proxy"] = proxy
                os.environ["https_proxy"] = proxy
                proxy_handler = urllib.request.ProxyHandler({'http': proxy, 'https': proxy})
                opener = urllib.request.build_opener(proxy_handler)
            else:
                print(f"   Downloading zip (Direct)...")
                # Remove proxy env vars for direct connection
                if "http_proxy" in os.environ: del os.environ["http_proxy"]
                if "https_proxy" in os.environ: del os.environ["https_proxy"]
                # Use empty proxy handler to be sure
                opener = urllib.request.build_opener(urllib.request.ProxyHandler({}))
            
            opener.addheaders = [('User-Agent', 'Mozilla/5.0')]
            
            with opener.open(url, timeout=30) as response:
                return response.read()
                
        except Exception as e:
            print(f"     Failed with {proxy or 'Direct'}: {e}")
    return None

def attempt_download(skill_name, repo_url, repo_map):
    rid = get_repo_id(skill_name, repo_map)
    if not rid: return

    dest_path = os.path.join("public", "data", "repos", rid, "content")
    if os.path.exists(dest_path) and len(os.listdir(dest_path)) > 0:
        print(f"Skipping {skill_name}: Already has content")
        return

    print(f"⬇️ Attempting download for {skill_name}: {repo_url}")
    
    # Basic logic: treat as github repo zip download
    if "github.com" not in repo_url:
        print(f"⚠️ {skill_name} is not a GitHub URL. Skipping.")
        return

    # Clean URL for zip download
    subdir = ""
    parts = repo_url.split('/')
    if len(parts) >= 5:
        owner = parts[3]
        repo = parts[4]
        # specific handling for tree/main
        zip_url = f"https://github.com/{owner}/{repo}/archive/refs/heads/main.zip"
        
        if "tree/main" in repo_url:
            try:
                main_idx = parts.index("main")
                subdir = "/".join(parts[main_idx+1:])
            except ValueError:
                pass
        
        # Try both 7890 and direct
        proxies = ["http://127.0.0.1:7890", None]
        
        zip_content = attempt_download_with_retry(zip_url, dest_path, subdir, proxies)
        
        if not zip_content:
            # Fallback to master
            print("   'main' branch failed, trying 'master'...")
            zip_url = zip_url.replace("main.zip", "master.zip")
            zip_content = attempt_download_with_retry(zip_url, dest_path, subdir, proxies)
            
        if not zip_content:
            print(f"❌ Failed to download {skill_name} after all attempts.")
            return

        try:
            with zipfile.ZipFile(io.BytesIO(zip_content)) as z:
                # Find root folder in zip
                root_in_zip = z.namelist()[0].split('/')[0]
                
                target_prefix = f"{root_in_zip}/{subdir}" if subdir else root_in_zip
                if not target_prefix.endswith('/'): target_prefix += '/'

                found_files = False
                os.makedirs(dest_path, exist_ok=True)
                
                for file_info in z.infolist():
                    if file_info.filename.startswith(target_prefix) and not file_info.is_dir():
                        # relative path inside content
                        rel_path = file_info.filename[len(target_prefix):]
                        if not rel_path: continue
                        
                        target_file = os.path.join(dest_path, rel_path)
                        os.makedirs(os.path.dirname(target_file), exist_ok=True)
                        with z.open(file_info) as source, open(target_file, "wb") as target:
                            shutil.copyfileobj(source, target)
                        found_files = True
                
                if not found_files:
                    # HEURISTIC: Try to find a matching folder if exact subdir failed
                    print(f"   ⚠️ Exact subdir '{subdir}' not found. Searching for alternatives...")
                    target_name = subdir.split('/')[-1] if subdir else ""
                    best_match = None
                    for name in z.namelist():
                        if file_info.is_dir() and target_name in name:
                             # simple heuristic
                             pass
                    
                    # Just dump top level directories to help debug or implement better search
                    top_dirs = set()
                    for name in z.namelist():
                        parts = name.split('/')
                        if len(parts) > 1: top_dirs.add(parts[1])
                    print(f"      Available top-level dirs in zip: {list(top_dirs)[:10]}")
                    
                    # Try 'src/{target_name}' or just 'src' or '{target_name}'
                    # For modelcontextprotocol/servers, it is usually src/
                    possible_subdirs = [f"src/{target_name}", target_name, f"servers/{target_name}"]
                    
                    for try_sub in possible_subdirs:
                        if try_sub == subdir: continue
                        print(f"      Trying alternative subdir: {try_sub}")
                        target_prefix = f"{root_in_zip}/{try_sub}"
                        if not target_prefix.endswith('/'): target_prefix += '/'
                        
                        found_alt = False
                        for file_info in z.infolist():
                            if file_info.filename.startswith(target_prefix) and not file_info.is_dir():
                                rel_path = file_info.filename[len(target_prefix):]
                                target_file = os.path.join(dest_path, rel_path)
                                os.makedirs(os.path.dirname(target_file), exist_ok=True)
                                with z.open(file_info) as source, open(target_file, "wb") as target:
                                    shutil.copyfileobj(source, target)
                                found_alt = True
                        
                        if found_alt:
                            print(f"✅ Found and extracted {skill_name} in {try_sub}")
                            found_files = True
                            break
                            
                if found_files:
                    print(f"✅ Downloaded and extracted {skill_name}")
                else:
                    print(f"❌ Content not found in zip for {skill_name} (subdir: {subdir})")

        except Exception as e:
            print(f"❌ Zip processing failed for {skill_name}: {e}")

def main():
    repo_map = json.load(open('src/data/repo_map.json', encoding='utf-8'))
    skills_data = json.load(open('src/data/skills.json', encoding='utf-8'))
    all_skills = [s for c in skills_data for s in c['skills']]
    
    # 1. Local Copies (Anthropic Skills)
    local_copies = {
        "Anthropic PPTX Skill": "pptx",
        "Anthropic XLSX Skill": "xlsx",
        "Anthropic DOCX Skill": "docx",
        "Anthropic PDF Skill": "pdf",
        "Anthropic Web App Testing": "webapp-testing",
        "Anthropic Algorithmic Art": "algorithmic-art",
        "Anthropic MCP Builder": "mcp-builder"
    }
    
    print("--- Phase 1: Local Copies ---")
    for name, subpath in local_copies.items():
        copy_local_skill(name, subpath, repo_map)
        
    # 2. Non-Repo Placeholders
    placeholders = [
        "LeetCode", "Khanmigo", "Claude for Education", 
        "Duolingo Max", "Speak", "Quizlet Q-Chat", "ChatGLM-Edu", "MathGPT", "BioContext"
    ]
    
    print("\n--- Phase 2: Placeholders ---")
    for s in all_skills:
        if s['name'] in placeholders:
            create_placeholder(s['name'], s.get('url', ''), repo_map)

    # 3. Retry Failures with likely corrections
    manual_corrections = {
        # Attempt to correct URLs if they are wrong or point to monorepos
        "WolframAlpha": "https://github.com/modelcontextprotocol/servers/tree/main/src/wolfram-alpha",
        "Sentry MCP": "https://github.com/modelcontextprotocol/servers/tree/main/src/sentry",
        "Kubernetes MCP": "https://github.com/modelcontextprotocol/servers/tree/main/src/kubernetes", # Try likely location
        "Neon MCP": "https://github.com/neondatabase/mcp-server-neon",
        "Todoist MCP": "https://github.com/actionsOnGoogle/todoist-mcp-server",
        "Discord MCP": "https://github.com/v-v-vishnu/discord-mcp-server",
        "UnrealZoo": "https://github.com/UnrealZoo/UnrealZoo",
        "claude-skill-data-cleaner": "https://github.com/brook-miller/claude-skill-data-cleaner"
    }
    
    print("\n--- Phase 3: Retry Downloads ---")
    for name, url in manual_corrections.items():
        attempt_download(name, url, repo_map)

if __name__ == "__main__":
    main()
