import json
import os

skills_json_path = 'src/data/skills.json'
backup_path = 'src/data/skills_backup.json'

# New Category Definitions
NEW_CATEGORIES = {
    "Education & Tutoring": {
        "zh": "教育与辅导",
        "sources": ["Educational AI & Tutoring", "Math & Science Assistant"]
    },
    "Academic & Writing": {
        "zh": "学术与写作",
        "sources": ["Research & Academic", "Writing & Content Creation"]
    },
    "Coding & Data": {
        "zh": "编程与数据",
        "sources": ["Coding & Technical Education", "Data Analysis & Visualization", "Accessibility & UX Design"]
    },
    "Productivity & Career": {
        "zh": "生产力与职业",
        "sources": ["Productivity & Second Brain", "Career & Resume"]
    },
    "MCP & Meta Skills": {
        "zh": "MCP 与元技能",
        "sources": ["MCP & Meta Skills"]
    }
}

def merge_skills():
    if not os.path.exists(skills_json_path):
        print(f"Error: {skills_json_path} not found.")
        return

    with open(skills_json_path, 'r', encoding='utf-8') as f:
        data = json.load(f)

    # Save backup
    with open(backup_path, 'w', encoding='utf-8') as f:
        json.dump(data, f, ensure_ascii=False, indent=2)
    print(f"Backed up to {backup_path}")

    new_data_map = {} # Key: New English Name, Value: Item dict

    for new_en, info in NEW_CATEGORIES.items():
        new_data_map[new_en] = {
            "category": {
                "en": new_en,
                "zh": info["zh"]
            },
            "description": {
                "en": "", # Will aggregate or set generic default
                "zh": ""
            },
            "skills": []
        }

    # Helper to find target category
    def get_target_cat(old_en):
        for new_en, info in NEW_CATEGORIES.items():
            if old_en in info["sources"]:
                return new_en
        return "MCP & Meta Skills" # Fallback

    # Process existing data
    for category_item in data:
        old_en = category_item['category']['en']
        target_cat = get_target_cat(old_en)
        
        # Append skills
        new_data_map[target_cat]['skills'].extend(category_item['skills'])
        
        # Combine descriptions (optional, or just pick one, or write new ones)
        # For now, let's keep the description of the "primary" source if empty
        current_desc = new_data_map[target_cat]['description']['en']
        if not current_desc:
            new_data_map[target_cat]['description'] = category_item['description']

    # Convert map back to list
    new_data_list = list(new_data_map.values())
    
    # Sort by the order defined in NEW_CATEGORIES
    ordered_list = []
    for key in NEW_CATEGORIES.keys():
        if key in new_data_map:
            ordered_list.append(new_data_map[key])

    # Update descriptions for merged categories to be more inclusive
    # (Simple hardcoded updates for better quality)
    descriptions = {
        "Education & Tutoring": {
            "en": "General teaching skills, math/science assistants, and educational tools.",
            "zh": "通用教学技能、数理科学助手及教育辅助工具。"
        },
        "Academic & Writing": {
            "en": "Academic research, paper writing, and creative content creation.",
            "zh": "学术研究、论文写作及创意内容创作。"
        },
        "Coding & Data": {
            "en": "Computer science education, data analysis, visualization, and technical design.",
            "zh": "计算机科学教育、数据分析、可视化及技术设计。"
        },
        "Productivity & Career": {
            "en": "Personal productivity, knowledge management, and career development.",
            "zh": "个人生产力、知识管理及职业发展。"
        },
        "MCP & Meta Skills": {
            "en": "Meta-skills for prompting and Model Context Protocol servers.",
            "zh": "提示词元技能及 Model Context Protocol (MCP) 服务器。"
        }
    }

    for item in ordered_list:
        cat_en = item['category']['en']
        if cat_en in descriptions:
            item['description'] = descriptions[cat_en]

    # Write back
    with open(skills_json_path, 'w', encoding='utf-8') as f:
        json.dump(ordered_list, f, ensure_ascii=False, indent=2)
    print(f"Updated {skills_json_path} with {len(ordered_list)} categories.")

if __name__ == "__main__":
    merge_skills()
