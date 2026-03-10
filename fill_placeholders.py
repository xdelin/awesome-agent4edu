import os
import json

def main():
    repo_map = json.load(open('src/data/repo_map.json', encoding='utf-8'))
    skills_data = json.load(open('src/data/skills.json', encoding='utf-8'))
    all_skills = [s for c in skills_data for s in c['skills']]

    count = 0
    for s in all_skills:
        name = s['name']
        rid = repo_map.get(name)
        if not rid: continue
        
        dest_path = os.path.join("public", "data", "repos", rid, "content")
        os.makedirs(dest_path, exist_ok=True)
        
        # Check if empty (or only has empty subdirs)
        is_empty = True
        if os.path.exists(dest_path):
            files = os.listdir(dest_path)
            if files:
                is_empty = False
        
        if is_empty:
            print(f"📝 Creating placeholder for {name} ({rid})...")
            readme_path = os.path.join(dest_path, "README.md")
            url = s.get('url', '#')
            content = f"""# {name}

This resource could not be automatically downloaded or does not have a public code repository viewable here.

**Project URL**: [{url}]({url})

Please visit the official website or repository for more details.
"""
            with open(readme_path, "w", encoding="utf-8") as f:
                f.write(content)
            count += 1
            
    print(f"✅ Filled {count} missing repositories with placeholders.")

if __name__ == "__main__":
    main()
