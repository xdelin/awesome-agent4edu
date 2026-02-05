import json
import os

skills_json_path = 'src/data/skills.json'

new_skills = [
    {
        "name": "Office-PowerPoint-MCP-Server",
        "url": "https://github.com/GongRzhe/Office-PowerPoint-MCP-Server",
        "description": {
            "en": "A powerful MCP server for PowerPoint manipulation using python-pptx. Create, edit, and manipulate presentations.",
            "zh": "一个强大的 PowerPoint 操作 MCP 服务器。使用 python-pptx 创建、编辑和处理演示文稿。"
        },
        "useCase": {
            "en": "Automating slide creation, updating charts, and bulk editing.",
            "zh": "自动化幻灯片创建、更新图表和批量编辑。"
        },
        "type": "MCP"
    },
    {
        "name": "mcp-server-okppt",
        "url": "https://github.com/NeekChaw/mcp-server-okppt",
        "description": {
            "en": "Generates high-quality PPTX slides from SVG images created by LLMs. Ensures vector quality.",
            "zh": "从 LLM 生成的 SVG 图像生成高质量 PPTX 幻灯片。确保矢量质量。"
        },
        "useCase": {
            "en": "Creating visually designing slides with complex diagrams directly from Claude.",
            "zh": "直接从 Claude 创建包含复杂图表的视觉设计幻灯片。"
        },
        "type": "MCP"
    },
    {
        "name": "pptx-mcp",
        "url": "https://github.com/samos123/pptx-mcp",
        "description": {
            "en": "Simple MCP server to create slides using Python PPTX library.",
            "zh": "使用 Python PPTX 库创建幻灯片的简单 MCP 服务器。"
        },
        "useCase": {
            "en": "Quickly generating basic slide decks from text outlines.",
            "zh": "从文本大纲快速生成基本幻灯片。"
        },
        "type": "MCP"
    }
]

def add_skills():
    if not os.path.exists(skills_json_path):
        print(f"Error: {skills_json_path} not found.")
        return

    with open(skills_json_path, 'r', encoding='utf-8') as f:
        data = json.load(f)

    target_category = "Productivity & Career"
    
    found = False
    for category in data:
        if category['category']['en'] == target_category:
            existing_names = {s['name'] for s in category['skills']}
            for skill in new_skills:
                if skill['name'] not in existing_names:
                    category['skills'].insert(0, skill) # Add to top
                    print(f"Added {skill['name']}")
                else:
                    print(f"Skipped {skill['name']} (already exists)")
            found = True
            break
    
    if not found:
        print(f"Category {target_category} not found!")

    with open(skills_json_path, 'w', encoding='utf-8') as f:
        json.dump(data, f, ensure_ascii=False, indent=2)
    print("Done.")

if __name__ == "__main__":
    add_skills()
