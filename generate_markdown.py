import json
import os
import re
import datetime

# Define paths
skills_json_path = 'src/data/skills.json'
output_dir = '../awesome-education-mcp'
readme_cn_path = os.path.join(output_dir, 'README_CN.md')
readme_en_path = os.path.join(output_dir, 'README.md')

# Emoji Map
emoji_map = {
    "Education & Tutoring": "🏫",
    "Academic & Writing": "📚",
    "Coding & Data": "💻",
    "Visual & Presentation": "🎨",
    "Productivity & Career": "🧠",
    "MCP & Meta Skills": "🛠️",
    "Commercial Education AI": "🏭",
    "Education LLMs": "🤖",
    "Embodied Intelligence & VR Education": "🥽",
    "Agent Frameworks & Educational Applications": "🤖"
}

def load_skills():
    with open(skills_json_path, 'r', encoding='utf-8') as f:
        return json.load(f)

def clean_anchor(text):
    # Remove emoji, special chars, keep only alphanumeric and hyphens
    # 1. Remove & 
    text = text.replace("&", "")
    # 2. Replace spaces with hyphens
    text = text.replace(" ", "-")
    # 3. Lowercase
    text = text.lower()
    # 4. Remove any other non-url friendly chars (simple regex)
    text = re.sub(r'[^a-z0-9\-]', '', text)
    # 5. Remove duplicate hyphens
    text = re.sub(r'-+', '-', text)
    return text.strip('-')

def generate_readme_cn(data):
    # Generates README_CN.md (All Content)
    content = "# Awesome Education AI: MCP, Skills & LLM Apps\n\n"
    content += "一份精心策划的教育 AI 资源列表，涵盖 **Model Context Protocol (MCP)** 服务器、**Claude Skills**、**教育类 LLM** 和 **Agent 框架**，专注于提升学术研究、教学效率和个性化学习体验。\n\n"
    content += "[English](./README.md) | [中文](./README_CN.md)\n\n"
    content += "## 目录\n\n"

    # TOC
    for category in data:
        # Include all items
        items = category['skills']
        if not items: continue
        cat_zh = category['category']['zh']
        cat_en = category['category']['en']
        emoji = emoji_map.get(cat_en, "📂")
        anchor = clean_anchor(cat_en)
        content += f"- [{emoji} {cat_zh}](#{anchor})\n"
    
    content += "\n---\n\n"

    for category in data:
        items = category['skills']
        if not items: continue

        cat_en = category['category']['en']
        cat_zh = category['category']['zh']
        emoji = emoji_map.get(cat_en, "📂")
        anchor = clean_anchor(cat_en)
        
        content += f"## <a id='{anchor}'></a>{emoji} {cat_en} ({cat_zh})\n\n"
        
        desc_zh = category['description']['zh']
        if desc_zh:
            content += f"> {desc_zh}\n\n"
        
        for skill in items:
            name = skill['name']
            url = skill['url']
            desc = skill['description']['zh']
            use_case = skill['useCase']['zh']
            # type = skill.get('type', 'Tool')
            
            content += f"- **[{name}]({url})**\n"
            content += f"  - **描述**: {desc}\n"
            content += f"  - **适用场景**: {use_case}\n\n"

    content += "---\n\n"
    content += f"*Last updated/最后更新: {datetime.date.today()}*\n"
    return content

def generate_readme_en(data):
    # Generates README.md (All Content)
    content = "# Awesome Education AI: MCP, Skills & LLM Apps\n\n"
    content += "A comprehensive curated list of AI resources for education, including **Model Context Protocol (MCP)** servers, **Claude Skills**, **Education LLMs**, and **Agent Frameworks**, focused on academic research, teaching efficiency, and personalized learning.\n\n"
    content += "[English](./README.md) | [中文](./README_CN.md)\n\n"
    content += "## Table of Contents\n\n"

    # TOC
    for category in data:
        items = category['skills']
        if not items: continue
        cat_en = category['category']['en']
        emoji = emoji_map.get(cat_en, "📂")
        anchor = clean_anchor(cat_en)
        content += f"- [{emoji} {cat_en}](#{anchor})\n"
    
    content += "\n---\n\n"

    for category in data:
        items = category['skills']
        if not items: continue

        cat_en = category['category']['en']
        emoji = emoji_map.get(cat_en, "📂")
        anchor = clean_anchor(cat_en)
        
        content += f"## <a id='{anchor}'></a>{emoji} {cat_en}\n\n"
        
        desc_en = category['description']['en']
        if desc_en:
            content += f"> {desc_en}\n\n"
        
        for skill in items:
            name = skill['name']
            url = skill['url']
            desc = skill['description']['en']
            use_case = skill['useCase']['en']
            
            content += f"- **[{name}]({url})**\n"
            content += f"  - **Description**: {desc}\n"
            content += f"  - **Use Case**: {use_case}\n\n"

    content += "---\n\n"
    content += f"*Last updated: {datetime.date.today()}*\n"
    return content

def main():
    data = load_skills()
    
    readme_cn = generate_readme_cn(data)
    with open(readme_cn_path, 'w', encoding='utf-8') as f:
        f.write(readme_cn)
    print(f"Generated {readme_cn_path}")

    readme_en = generate_readme_en(data)
    with open(readme_en_path, 'w', encoding='utf-8') as f:
        f.write(readme_en)
    print(f"Generated {readme_en_path}")

if __name__ == "__main__":
    main()
