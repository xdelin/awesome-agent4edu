import json
import os

skills_json_path = 'src/data/skills.json'

NEW_CATEGORIES = {
    # Keep existing ones mostly, but split Productivity
    "Education & Tutoring": {
        "en": "Education & Tutoring",
        "zh": "教育与辅导"
    },
    "Academic & Writing": {
        "en": "Academic & Writing",
        "zh": "学术与写作"
    },
    "Coding & Data": {
        "en": "Coding & Data",
        "zh": "编程与数据"
    },
    # New Split Category for PPT/Presentation
    "Visual & Presentation": {
        "en": "Visual & Presentation",
        "zh": "视觉与演示 (PPT)" 
    },
    "Productivity & Career": {
        "en": "Productivity & Career",
        "zh": "生产力与职业"
    },
    "MCP & Meta Skills": {
        "en": "MCP & Meta Skills",
        "zh": "MCP 与元技能"
    }
}

PPT_KEYWORDS = ["ppt", "powerpoint", "slide", "figma", "video", "youtube", "image", "visual"]

def reclassify_ppt():
    if not os.path.exists(skills_json_path):
        print(f"Error: {skills_json_path} not found.")
        return

    with open(skills_json_path, 'r', encoding='utf-8') as f:
        data = json.load(f)

    new_viz_cat = {
        "category": NEW_CATEGORIES["Visual & Presentation"],
        "description": {
            "en": "PowerPoint generation, visual design, and presentation tools.",
            "zh": "PowerPoint 生成、视觉设计及演示工具。"
        },
        "skills": []
    }

    # Helper to check if skill belongs to visual
    def is_visual(skill):
        text = (skill['name'] + " " + skill['description']['en']).lower()
        return any(k in text for k in PPT_KEYWORDS)

    # 1. Extract visual skills from existing categories
    for cat in data:
        remaining_skills = []
        for skill in cat['skills']:
            if is_visual(skill):
                new_viz_cat['skills'].append(skill)
                # print(f"Moved {skill['name']} to Visual & Presentation")
            else:
                remaining_skills.append(skill)
        cat['skills'] = remaining_skills

    # 2. Add the new category to the list (insert before Productivity)
    # Find index of Productivity
    prod_idx = -1
    for i, cat in enumerate(data):
        if cat['category']['en'] == "Productivity & Career":
            prod_idx = i
            break
    
    if prod_idx != -1:
        data.insert(prod_idx, new_viz_cat)
    else:
        # Append if not found
        data.append(new_viz_cat)

    # 3. Clean up empty categories if any (optional, but good practice)
    data = [c for c in data if len(c['skills']) > 0]

    with open(skills_json_path, 'w', encoding='utf-8') as f:
        json.dump(data, f, ensure_ascii=False, indent=2)
    print(f"Reclassified {len(new_viz_cat['skills'])} skills into Visual & Presentation.")

if __name__ == "__main__":
    reclassify_ppt()
