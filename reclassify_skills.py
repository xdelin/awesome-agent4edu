import json

# Define the new categories and their descriptions
new_categories = {
    "Computer Science": {
        "en": "Computer Science",
        "zh": "计算机科学",
        "description_en": "Development tools, coding assistants, database management, and DevOps.",
        "description_zh": "开发工具、编程助手、数据库管理及 DevOps。",
        "skills": []
    },
    "Math & Science": {
        "en": "Math & Science",
        "zh": "数理科学",
        "description_en": "Mathematics, physics, STEM research, and scientific tools.",
        "description_zh": "数学、物理、STEM 研究及科学工具。",
        "skills": []
    },
    "Data & Analysis": {
        "en": "Data & Analysis",
        "zh": "数据与分析",
        "description_en": "Data science, web scraping, analytics, and big data tools.",
        "description_zh": "数据科学、网络抓取、分析及大数据工具。",
        "skills": []
    },
    "Academic & Writing": {
        "en": "Academic & Writing",
        "zh": "学术与写作",
        "description_en": "Academic research, paper writing, literature reviews, and creative writing.",
        "description_zh": "学术研究、论文写作、文献综述及创意写作。",
        "skills": []
    },
    "Visual & Presentation": {
        "en": "Visual & Presentation",
        "zh": "视觉与演示 (PPT)",
        "description_en": "PowerPoint generation, data visualization, UI design, and visual tools.",
        "description_zh": "PowerPoint 生成、数据可视化、UI 设计及视觉工具。",
        "skills": []
    },
    "Notes & Knowledge Base": {
        "en": "Notes & Knowledge Base",
        "zh": "笔记与知识库",
        "description_en": "Note-taking, personal knowledge management (PKM), and document processing.",
        "description_zh": "笔记整理、个人知识库管理 (PKM) 及文档处理。",
        "skills": []
    },
    "Career & Productivity": {
        "en": "Career & Productivity",
        "zh": "职业规划与生产力",
        "description_en": "Career development, project management, communication, and task tracking.",
        "description_zh": "职业发展、项目管理、沟通协作及任务追踪。",
        "skills": []
    },
    "Intelligent Tutoring": {
        "en": "Intelligent Tutoring",
        "zh": "智能导学",
        "description_en": "AI tutors, educational assessment, and personalized learning assistants.",
        "description_zh": "AI 导师、教育评估及个性化学习助手。",
        "skills": []
    },
    "Meta Skills": {
        "en": "Meta Skills",
        "zh": "元技能",
        "description_en": "Prompt engineering, MCP meta-lists, and skill development tools.",
        "description_zh": "提示词工程、MCP 元列表及技能开发工具。",
        "skills": []
    }
}

# Mapping of skill names to new categories
skill_mapping = {
    # Computer Science
    "GitHub MCP Server": "Computer Science",
    "Playwright MCP": "Computer Science",
    "Desktop Commander": "Computer Science",
    "Kubernetes MCP": "Computer Science",
    "Sentry MCP": "Computer Science",
    "Octocode MCP": "Computer Science",
    "MindsDB MCP": "Computer Science",
    "SQLite Explorer": "Computer Science",
    "DBHub": "Computer Science",
    "Neon MCP": "Computer Science",
    "awesome-vibe-coding": "Computer Science",
    "coding-agent-skills": "Computer Science",
    "test-generator": "Computer Science",
    "mlx-dev-skill": "Computer Science",
    "prompt-builder": "Computer Science",
    "paper2code-skill": "Computer Science",
    "Docker MCP": "Computer Science",
    "a11y-specialist-skills": "Computer Science",

    # Math & Science
    "Kvante": "Math & Science",
    "Arxiv MCP Server": "Math & Science",

    # Data & Analysis
    "Firecrawl MCP": "Data & Analysis",
    "BrowserBase MCP": "Data & Analysis",
    "GenAI Toolbox": "Data & Analysis",
    "claude-skill-data-cleaner": "Data & Analysis",
    "activitywatch-analysis-skill": "Data & Analysis",

    # Academic & Writing
    "writing-agent": "Academic & Writing",
    "The-Crucible-Writing-System-For-Claude": "Academic & Writing",
    "humanizer": "Academic & Writing",
    "ux-writing-skill": "Academic & Writing",
    "vibe-writing": "Academic & Writing",
    "academic-paper-skills": "Academic & Writing",
    "research-units-pipeline-skills": "Academic & Writing",
    "academic-writing-skills": "Academic & Writing",
    "GPT Researcher": "Academic & Writing",
    "ZotLink": "Academic & Writing",

    # Visual & Presentation
    "pptx-mcp": "Visual & Presentation",
    "mcp-server-okppt": "Visual & Presentation",
    "Office-PowerPoint-MCP-Server": "Visual & Presentation",
    "Figma Context MCP": "Visual & Presentation",
    "census-demographics-skill": "Visual & Presentation",
    "Observable-Plot-Claude-Skill": "Visual & Presentation",
    "chartjs-expert": "Visual & Presentation",
    "neurodivergent-visual-org": "Visual & Presentation",
    "Git MCP": "Visual & Presentation",
    "claude-dolphin": "Visual & Presentation",

    # Notes & Knowledge Base
    "claudian": "Notes & Knowledge Base",
    "claude-code-obsidian-starter": "Notes & Knowledge Base",
    "productivity-skills": "Notes & Knowledge Base",
    "NotebookLM MCP": "Notes & Knowledge Base",
    "Official Notion MCP Server": "Notes & Knowledge Base",
    "Anki MCP Server": "Notes & Knowledge Base",
    "PDF Reader MCP": "Notes & Knowledge Base",

    # Career & Productivity
    "Resume-Analysis-Assistant": "Career & Productivity",
    "JobCoach": "Career & Productivity",
    "linkedin-mcp": "Career & Productivity",
    "Linear MCP": "Career & Productivity",
    "Todoist MCP": "Career & Productivity",
    "Slack MCP": "Career & Productivity",
    "Discord MCP": "Career & Productivity",

    # Intelligent Tutoring
    "claude-educational-ai-skills": "Intelligent Tutoring",
    "EdTech-AI-Coach": "Intelligent Tutoring",
    "ELA-Tutor-Agent": "Intelligent Tutoring",
    "quiz-master": "Intelligent Tutoring",
    "apprenticemode": "Intelligent Tutoring",

    # Meta Skills
    "awesome-claude-skills (mejba13)": "Meta Skills",
    "awesome-claude-skills (ponderous)": "Meta Skills",
    "everything-claude-code": "Meta Skills",
    "skill-description-optimizer": "Meta Skills",
    "Skill Seekers": "Meta Skills",
    "claude-prompt-engineering-guide": "Meta Skills"
}

source_file = 'src/data/skills.json'
output_file = 'src/data/skills.json'
backup_file = 'src/data/skills_backup_v2.json'

with open(source_file, 'r', encoding='utf-8') as f:
    data = json.load(f)

# Backup
with open(backup_file, 'w', encoding='utf-8') as f:
    json.dump(data, f, ensure_ascii=False, indent=2)

# Reclassify
for group in data:
    for skill in group['skills']:
        skill_name = skill['name']
        if skill_name in skill_mapping:
            new_cat_key = skill_mapping[skill_name]
            new_categories[new_cat_key]['skills'].append(skill)
        else:
            print(f"Warning: Skill '{skill_name}' not mapped. Placing in Meta Skills.")
            new_categories["Meta Skills"]['skills'].append(skill)

# Construct final list
final_data = []
# Ensure order conforms to script definition order
ordered_keys = [
    "Intelligent Tutoring",
    "Math & Science",
    "Computer Science",
    "Data & Analysis",
    "Visual & Presentation",
    "Academic & Writing",
    "Notes & Knowledge Base",
    "Career & Productivity",
    "Meta Skills"
]

for key in ordered_keys:
    cat_data = new_categories[key]
    if cat_data['skills']:
        final_data.append({
            "category": {
                "en": cat_data['en'],
                "zh": cat_data['zh']
            },
            "description": {
                "en": cat_data['description_en'],
                "zh": cat_data['description_zh']
            },
            "skills": cat_data['skills']
        })

with open(output_file, 'w', encoding='utf-8') as f:
    json.dump(final_data, f, ensure_ascii=False, indent=2)

print(f"Successfully reclassified skills into {len(final_data)} categories.")
