"""
Update skills.json with GitHub stars/author for ALL skills.
Handles both standard repo URLs and subdirectory URLs.
Uses system proxy if available.
"""
import json
import os
import re
import time
import urllib.request
import urllib.error

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
SKILLS_FILE = os.path.join(BASE_DIR, "src", "data", "skills.json")

GITHUB_TOKEN = os.environ.get("GITHUB_TOKEN", "")

# Setup proxy
proxy_handler = urllib.request.ProxyHandler({
    'http': 'http://127.0.0.1:7890',
    'https': 'http://127.0.0.1:7890',
})
opener = urllib.request.build_opener(proxy_handler)
urllib.request.install_opener(opener)

def parse_github_url(url):
    """Extract owner/repo from various GitHub URL formats."""
    url = url.rstrip('/')
    # Remove https://github.com/
    path = re.sub(r'^https?://github\.com/', '', url)
    parts = path.split('/')
    if len(parts) >= 2:
        return parts[0], parts[1]
    return None, None

def get_repo_info(owner, repo):
    url = f"https://api.github.com/repos/{owner}/{repo}"
    req = urllib.request.Request(url)
    if GITHUB_TOKEN:
        req.add_header("Authorization", f"Bearer {GITHUB_TOKEN}")
    req.add_header("User-Agent", "Skills-Meta-Updater")
    
    try:
        with urllib.request.urlopen(req, timeout=15) as response:
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
            # Check rate limit
            try:
                remaining = e.headers.get('X-RateLimit-Remaining', '?')
                reset_time = int(e.headers.get('X-RateLimit-Reset', '0'))
                wait = max(reset_time - int(time.time()), 10)
                print(f"    Rate limited! Remaining: {remaining}. Waiting {wait}s...")
                time.sleep(min(wait, 120))
                return get_repo_info(owner, repo)  # retry
            except:
                print(f"    Rate limited, waiting 60s...")
                time.sleep(60)
                return get_repo_info(owner, repo)
        elif e.code == 404:
            print(f"    404 Not Found")
            return None
        else:
            print(f"    HTTP Error {e.code}")
            return None
    except Exception as e:
        print(f"    Error: {e}")
        return None

def main():
    with open(SKILLS_FILE, 'r', encoding='utf-8') as f:
        data = json.load(f)

    updated = 0
    skipped = 0
    failed = 0
    total_skills = sum(len(c.get('skills', [])) for c in data)
    
    print(f"Total skills: {total_skills}")
    if GITHUB_TOKEN:
        print("Using GITHUB_TOKEN (5000 req/hr)")
    else:
        print("WARNING: No GITHUB_TOKEN, limited to 60 req/hr")
    print()

    count = 0
    for category in data:
        cat_name = category.get('category', {}).get('en', 'Unknown')
        for skill in category.get("skills", []):
            count += 1
            url = skill.get("url", "")
            
            # Skip if already has valid data
            existing_stars = skill.get('stars') or (skill.get('github', {}) or {}).get('stars')
            existing_author = skill.get('author') or (skill.get('github', {}) or {}).get('author')
            if existing_stars and existing_author:
                skipped += 1
                continue
            
            if "github.com" not in url:
                # Try to extract author from URL for non-github
                skipped += 1
                continue

            owner, repo = parse_github_url(url)
            if not owner or not repo:
                print(f"[{count}/{total_skills}] {skill['name']}: Can't parse URL: {url}")
                failed += 1
                continue

            print(f"[{count}/{total_skills}] {skill['name']} -> {owner}/{repo}...", end=" ")
            
            meta = get_repo_info(owner, repo)
            if meta:
                skill["stars"] = meta["stars"]
                skill["author"] = meta["author"]
                skill["github"] = meta
                updated += 1
                print(f"⭐ {meta['stars']} by {meta['author']}")
            else:
                failed += 1
                # At least set author from URL
                if not existing_author:
                    skill["author"] = owner
                    print(f"(set author={owner} from URL)")
            
            # Small delay to avoid hammering
            time.sleep(0.5)

    # Save
    with open(SKILLS_FILE, 'w', encoding='utf-8') as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
    
    print(f"\nDone! Updated: {updated}, Skipped: {skipped}, Failed: {failed}")

if __name__ == '__main__':
    main()
