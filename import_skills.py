import json
import os
import shutil
import re

# Paths
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
SOURCE_SKILLS_DIR = os.path.normpath(os.path.join(BASE_DIR, "../1234/skills"))

DEST_REPOS_DIR = os.path.join(BASE_DIR, "public", "data", "repos")
SKILLS_JSON_PATH = os.path.join(BASE_DIR, "src", "data", "skills.json")
REGISTRY_JSON_PATH = os.path.join(BASE_DIR, "public", "registry.json")

def load_json(path):
    if os.path.exists(path):
        with open(path, 'r', encoding='utf-8') as f:
            return json.load(f)
    return {}

def save_json(path, data):
    with open(path, 'w', encoding='utf-8') as f:
        json.dump(data, f, indent=2, ensure_ascii=False)

def parse_frontmatter(file_path):
    metadata = {}
    content = ""
    if os.path.exists(file_path):
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                raw = f.read()
                # Simple regex for YAML frontmatter
                match = re.search(r'^---\s*\n(.*?)\n---\s*(.*)', raw, re.DOTALL)
                if match:
                    fm = match.group(1)
                    content = match.group(2)
                    for line in fm.split('\n'):
                        if ':' in line:
                            parts = line.split(':', 1)
                            key = parts[0].strip()
                            val = parts[1].strip().strip('"\'')
                            metadata[key] = val
                else:
                    content = raw
        except Exception as e:
            print(f"Error reading {file_path}: {e}")
    return metadata, content

import stat

def handle_remove_readonly(func, path, exc):
    excvalue = exc[1]
    if func in (os.rmdir, os.remove, os.unlink) and excvalue.errno == 13:
        os.chmod(path, stat.S_IWRITE)
        func(path)

def main():
    print(f"Source: {SOURCE_SKILLS_DIR}")
    print(f"Dest: {DEST_REPOS_DIR}")

    # Initialize data structures
    skills_categories = []
    registry = {}
    
    # "All Skills" category
    all_skills_cat = {
        "category": { "en": "All Skills", "zh": "全部技能" },
        "description": { "en": "Complete list of available skills", "zh": "所有可用技能列表" },
        "skills": []
    }

    if not os.path.exists(SOURCE_SKILLS_DIR):
        print("Source directory not found!")
        return

    # Clean destination
    if os.path.exists(DEST_REPOS_DIR):
        # Retry logic for windows permissions
        try:
             shutil.rmtree(DEST_REPOS_DIR, onerror=handle_remove_readonly)
        except Exception as e:
             print(f"Warning: Could not fully delete {DEST_REPOS_DIR}: {e}")
             
    os.makedirs(DEST_REPOS_DIR, exist_ok=True)

    skill_folders = [f for f in os.listdir(SOURCE_SKILLS_DIR) if os.path.isdir(os.path.join(SOURCE_SKILLS_DIR, f))]
    print(f"Found {len(skill_folders)} potential skills.")

    count = 0
    for skill_name in skill_folders:
        src_path = os.path.join(SOURCE_SKILLS_DIR, skill_name)
        skill_md_path = os.path.join(src_path, "SKILL.md")
        
        # Only process if SKILL.md exists
        if os.path.exists(skill_md_path):
            # Parse metadata first
            metadata, content = parse_frontmatter(skill_md_path)
            
            desc = metadata.get('description', 'No description provided.')
            author = metadata.get('author', 'Community')
            tags = [t.strip() for t in metadata.get('tags', '').split(',')] if metadata.get('tags') else []

            # Create destination zip
            # We want the zip to contain the contents at root level usually, 
            # Or inside a folder? standard is to have a folder inside. 
            # shutil.make_archive with root_dir=src_path creates zip with contents at root.
            # If we want folder inside, we use root_dir=parent and base_dir=folder.
            # Let's stick to root contents for simplicity, or check what claw expects.
            # Assuming claw unzips to a target folder.
            
            try:
                # Create zip file: public/data/repos/<skill_name>.zip
                zip_base_name = os.path.join(DEST_REPOS_DIR, skill_name)
                shutil.make_archive(zip_base_name, 'zip', src_path)
                
                # Also copy SKILL.md separately for web viewer (optional but good)
                os.makedirs(os.path.join(DEST_REPOS_DIR, skill_name), exist_ok=True)
                shutil.copy2(skill_md_path, os.path.join(DEST_REPOS_DIR, skill_name, "SKILL.md"))
                
            except Exception as e:
                print(f"Error packing {skill_name}: {e}")
                continue

            # Registry entry points to zip
            registry[skill_name] = {
                "type": "hosted",
                "path": f"data/repos/{skill_name}.zip",
                "description": desc,
                "homepage": ""
            }

            # Skills JSON entry (unchanged)
            skill_entry = {
                "name": skill_name,
                "description": desc,
                "tags": tags,
                "author": author,
                "avatar": f"https://ui-avatars.com/api/?name={skill_name}&background=random",
                "stars": 0
            }
            all_skills_cat["skills"].append(skill_entry)
            count += 1

    skills_categories.append(all_skills_cat)
    
    # Save files
    save_json(SKILLS_JSON_PATH, skills_categories)
    save_json(REGISTRY_JSON_PATH, registry)
    
    print(f"Successfully imported {count} skills.")

if __name__ == "__main__":
    main()
