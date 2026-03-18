import json
import os
import shutil
import re

# Paths
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
SOURCE_SKILLS_DIR = os.path.normpath(os.path.join(BASE_DIR, "../1234/skills"))
CURATED_JSON_PATH = os.path.normpath(os.path.join(BASE_DIR, "../openclaw-skills-edu/curated_skills.json"))

DEST_REPOS_DIR = os.path.join(BASE_DIR, "public", "data", "repos")
SKILLS_JSON_PATH = os.path.join(BASE_DIR, "src", "data", "skills.json")
REPO_MAP_PATH = os.path.join(BASE_DIR, "src", "data", "repo_map.json")

def load_json(path):
    if os.path.exists(path):
        with open(path, 'r', encoding='utf-8') as f:
            return json.load(f)
    return {}

def save_json(path, data):
    with open(path, 'w', encoding='utf-8') as f:
        json.dump(data, f, indent=2, ensure_ascii=False)

def build_manifest(repo_path, prefix=""):
    tree = []
    for item in os.listdir(repo_path):
        item_path = os.path.join(repo_path, item)
        rel_path = os.path.join(prefix, item).replace("\\", "/")
        if os.path.isdir(item_path):
            tree.append({
                "name": item,
                "path": rel_path,
                "type": "folder",
                "children": build_manifest(item_path, rel_path)
            })
        else:
            tree.append({
                "name": item,
                "path": rel_path,
                "type": "file",
                "size": os.path.getsize(item_path),
                "url": "", # Local files don't need URL in manifest for this simple implementation
            })
    return tree

def normalize_category(cat_str):
    # Extract English part before (
    cat = cat_str.split('(')[0].strip()
    return cat

def main():
    print("Loading data...")
    curated_skills = load_json(CURATED_JSON_PATH)
    current_skills = load_json(SKILLS_JSON_PATH)
    repo_map = load_json(REPO_MAP_PATH)

    # Helper to find existing category in current_skills
    def find_category(cat_name):
        for entry in current_skills:
            if entry["category"]["en"] == cat_name:
                return entry
        return None

    # Counter
    total_imported = 0

    print(f"Processing skills from {CURATED_JSON_PATH}...")
    for cat_key, skills_list in curated_skills.items():
        cat_en = normalize_category(cat_key)
        
        # Map categories if needed
        # Just use name for now, create if not exists
        category_entry = find_category(cat_en)
        if not category_entry:
            # Create new category
            # Extract Chinese part if available
            cat_zh = cat_en
            match = re.search(r'\((.*?)\)', cat_key)
            if match:
                cat_zh = match.group(1)
            
            category_entry = {
                "category": { "en": cat_en, "zh": cat_zh },
                "description": { "en": "", "zh": "" },
                "skills": []
            }
            current_skills.append(category_entry)
        
        for skill in skills_list:
            title = skill.get("title")
            if not title: continue
            
            # Check if source exists
            source_path = os.path.join(SOURCE_SKILLS_DIR, title)
            if not os.path.exists(source_path):
                # Try finding if title is part of folder name? No, stick to exact match first
                # print(f"Skipping {title}: source folder not found in {SOURCE_SKILLS_DIR}")
                continue

            # Update repo map
            repo_id = title # Use title as repo_id
            repo_map[title] = repo_id

            # Prepare destination
            dest_path = os.path.join(DEST_REPOS_DIR, repo_id)
            dest_content_path = os.path.join(dest_path, "content")
            
            # Create content dir
            if os.path.exists(dest_content_path):
                try:
                    shutil.rmtree(dest_content_path)
                except PermissionError:
                    print(f"Warning: Could not delete existing content for {title}. Trying to overwrite.")
                except Exception as e:
                    print(f"Warning: Error deleting {dest_content_path}: {e}")
            
            if not os.path.exists(dest_content_path):
                 os.makedirs(dest_content_path)

            # Copy files
            try:
                # Copy contents of source_path to dest_content_path
                # shutil.copytree(source_path, dest_content_path, dirs_exist_ok=True) # py3.8+
                # For compatibility, iterate
                for item in os.listdir(source_path):
                    s = os.path.join(source_path, item)
                    d = os.path.join(dest_content_path, item)
                    if os.path.isdir(s):
                        shutil.copytree(s, d)
                    else:
                        shutil.copy2(s, d)
            except Exception as e:
                print(f"Error copying {title}: {e}")
                continue

            # Generate manifest
            manifest = {
                "tree": build_manifest(dest_content_path)
            }
            with open(os.path.join(dest_path, "manifest.json"), 'w', encoding='utf-8') as f:
                json.dump(manifest, f, indent=2)

            # Check if skill already in category
            existing_skill = next((s for s in category_entry["skills"] if s["name"] == title), None)
            
            # Determine import prompt
            install_prompt = f"claw install {title}"
            
            skill_entry = {
                "name": title,
                "url": "", # Local
                "description": {
                    "en": skill.get("description", ""),
                    "zh": skill.get("description", "")
                },
                "useCase": {
                    "en": "",
                    "zh": ""
                },
                "type": "Skill",
                "installPrompt": install_prompt,
                "stars": 0,
                "author": skill.get("author", "unknown"),
                "localPath": f"data/repos/{repo_id}"
            }

            if existing_skill:
                # update it
                existing_skill.update(skill_entry)
            else:
                category_entry["skills"].append(skill_entry)
            
            total_imported += 1

    print(f"Imported {total_imported} skills.")
    save_json(SKILLS_JSON_PATH, current_skills)
    save_json(REPO_MAP_PATH, repo_map)

    # Generate a public registry for the CLI
    registry = {}
    for cat in current_skills:
        for skill in cat["skills"]:
            name = skill["name"]
            # Determine source type
            if name in repo_map:
                # Hosted on this site
                registry[name] = {
                    "type": "hosted",
                    "manifest_url": f"data/repos/{repo_map[name]}/manifest.json",
                    "content_base_url": f"data/repos/{repo_map[name]}/content/",
                    "author": skill.get("author", "unknown"),
                    "description": skill.get("description", {}).get("en", "")
                }
            else:
                # External
                registry[name] = {
                    "type": "external",
                    "url": skill.get("url", ""),
                    "author": skill.get("author", "unknown"),
                    "description": skill.get("description", {}).get("en", "")
                }
    
    save_json(os.path.join(BASE_DIR, "public", "registry.json"), registry)
    print("Generated public/registry.json for CLI.")

if __name__ == "__main__":
    main()
