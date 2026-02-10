
import json
import os
import requests
import re
import time
from concurrent.futures import ThreadPoolExecutor

# Paths
base_path = r"C:\Users\ddsasd\OneDrive - 浮光浅夏\桌面\github\skills-web-interface"
skills_file = os.path.join(base_path, "src", "data", "skills.json")
readmes_dir = os.path.join(base_path, "public", "data", "readmes")
map_file = os.path.join(base_path, "src", "data", "readme_map.json")

def load_skills():
    with open(skills_file, 'r', encoding='utf-8') as f:
        return json.load(f)

def sanitize_filename(name):
    # simple sanitization
    return re.sub(r'[\\/*?:"<>|]', "", name).replace(" ", "_").lower() + ".md"

def fetch_readme(skill):
    name = skill.get('name')
    url = skill.get('url')
    
    if not url or "github.com" not in url:
        return name, None
    
    # Extract owner/repo
    match = re.search(r'github\.com/([^/]+)/([^/]+)', url)
    if not match:
        return name, None
    
    owner, repo = match.groups()
    repo = repo.rstrip('/') # remove trailing slash if any
    
    # Try different branches and filenames
    branches = ['main', 'master']
    filenames = ['README.md', 'readme.md', 'README.txt', 'readme.txt']
    
    content = None
    
    for branch in branches:
        for fname in filenames:
            raw_url = f"https://raw.githubusercontent.com/{owner}/{repo}/{branch}/{fname}"
            try:
                # Use a timeout to avoid hanging
                response = requests.get(raw_url, timeout=5)
                if response.status_code == 200:
                    content = response.text
                    break
            except Exception as e:
                print(f"Error fetching {name}: {e}")
        if content:
            break
            
    if content:
        # Save file
        filename = sanitize_filename(name)
        filepath = os.path.join(readmes_dir, filename)
        try:
            with open(filepath, 'w', encoding='utf-8') as f:
                f.write(content)
            print(f"✅ Fetched: {name}")
            return name, filename
        except Exception as e:
            print(f"❌ Write Error {name}: {e}")
            return name, None
    else:
        print(f"⚠️ Not Found: {name} (Tried {owner}/{repo})")
        return name, None

def main():
    if not os.path.exists(readmes_dir):
        os.makedirs(readmes_dir)
        
    data = load_skills()
    all_skills = []
    for cat in data:
        all_skills.extend(cat['skills'])
        
    readme_map = {}
    
    print(f"Fetching READMEs for {len(all_skills)} skills...")
    
    # Use threading to speed up
    with ThreadPoolExecutor(max_workers=10) as executor:
        results = executor.map(fetch_readme, all_skills)
        
    for name, filename in results:
        if filename:
            # We map the sanitized name relative to public/
            readme_map[name] = f"data/readmes/{filename}"
            
    # Save the map
    with open(map_file, 'w', encoding='utf-8') as f:
        json.dump(readme_map, f, indent=2)
        
    print(f"Done. Mapped {len(readme_map)} READMEs.")

if __name__ == "__main__":
    main()
