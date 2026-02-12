import json
import os
import time
import urllib.request
import urllib.error

# Configuration
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
SKILLS_FILE = os.path.join(BASE_DIR, "src", "data", "skills.json")
BACKUP_FILE = os.path.join(BASE_DIR, "src", "data", "skills.json.bak")

# Get Token from Env (or set it here)
GITHUB_TOKEN = os.environ.get("GITHUB_TOKEN")

def get_repo_info(owner, repo):
    url = f"https://api.github.com/repos/{owner}/{repo}"
    req = urllib.request.Request(url)
    if GITHUB_TOKEN:
        req.add_header("Authorization", f"Bearer {GITHUB_TOKEN}")
    req.add_header("User-Agent", "Skills-Web-Interface-Updater")
    
    try:
        with urllib.request.urlopen(req) as response:
            data = json.load(response)
            return {
                "stars": data.get("stargazers_count", 0),
                "forks": data.get("forks_count", 0),
                "author": data.get("owner", {}).get("login", ""),
                "authorAvatar": data.get("owner", {}).get("avatar_url", ""),
                "description": data.get("description", ""),
                "updatedAt": data.get("pushed_at", "")
            }
    except urllib.error.HTTPError as e:
        if e.code == 403:
            print(f"Rate limit exceeded for {owner}/{repo}. Waiting 60s...")
            time.sleep(60)
            return get_repo_info(owner, repo) # Retry once
        elif e.code == 404:
            print(f"Repo not found: {owner}/{repo}")
            return None
        else:
            print(f"Error fetching {owner}/{repo}: {e}")
            return None
    except Exception as e:
        print(f"Error: {e}")
        return None

def update_skills():
    if not os.path.exists(SKILLS_FILE):
        print("Skills file not found.")
        return

    # Backup
    with open(SKILLS_FILE, 'r', encoding='utf-8') as f:
        data = json.load(f)
    
    with open(BACKUP_FILE, 'w', encoding='utf-8') as f:
        json.dump(data, f, indent=2, ensure_ascii=False)

    total_updated = 0
    max_updates_for_demo = 15 # Limit for demo to avoid rate limits
    
    for category in data:

        print(f"Processing category: {category.get('category', {}).get('en', 'Unknown')}")
        for skill in category.get("skills", []):
            url = skill.get("url", "")
            if "github.com" in url:
                parts = url.rstrip("/").split("/")
                if len(parts) >= 5:
                    owner = parts[-2]
                    repo = parts[-1]
                    
                    # Skip if already has rich data and updated recently (optional, skipping for now to force update)
                    print(f"  Fetching: {owner}/{repo}...")
                    
                    meta = get_repo_info(owner, repo)
                    if meta:
                        skill["stars"] = meta["stars"]
                        skill["author"] = meta["author"]
                        # Optional: Use API description if manual one is missing?
                        # skill["description"]["en"] = meta["description"] 
                        
                        # Add a github specific object for more details
                        skill["github"] = meta
                        total_updated += 1
                        print(f"    -> Stars: {meta['stars']}, Author: {meta['author']}")

                        if total_updated >= max_updates_for_demo:
                            print(f"Reached demo limit of {max_updates_for_demo} updates. Stopping.")
                            break
        if total_updated >= max_updates_for_demo:
             break

    with open(SKILLS_FILE, 'w', encoding='utf-8') as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
    
    print(f"Done. Updated {total_updated} skills.")

if __name__ == "__main__":
    if not GITHUB_TOKEN:
        print("WARNING: No GITHUB_TOKEN set. Rate limits will be strict (60/hr).")
        print("Set it via: $env:GITHUB_TOKEN='your_token'")
        time.sleep(2)
    
    update_skills()
