# Awesome Education MCP / 教育领域 MCP 服务器精选

精选的 **Model Context Protocol (MCP)** 服务器列表，专注于教育、学术研究、生产力和知识管理。

如果您正在寻找 Claude Prompts & Skills，请查看 [CLAUDE_SKILLS.md](./CLAUDE_SKILLS.md)。

[English](./README.md) | [中文](./README_CN.md)

## 目录

- [📚 学术与写作](#academic-writing)
- [💻 编程与数据](#coding-data)
- [🎨 视觉与演示 (PPT)](#visual-presentation)
- [🧠 生产力与职业](#productivity-career)

---

## <a id='academic-writing'></a>📚 Academic & Writing (学术与写作)

> 学术研究、论文写作及创意内容创作。

- **[Arxiv MCP Server](https://github.com/blazickjp/arxiv-mcp-server)**
  - **描述**: 一个专门用于从 arXiv 数据库搜索、检索和分析学术论文的 MCP 服务器。
  - **适用场景**: 需要在 Claude 中直接快速访问最新预印本的 STEM 研究人员。

- **[ZotLink](https://github.com/TonybotNi/ZotLink)**
  - **描述**: 生产就绪的 Zotero MCP 服务器。可保存带有丰富元数据和智能 PDF 附件的开放预印本。
  - **适用场景**: 管理大量参考书目的学者，希望自动化“从论文到资料库”的工作流程。

- **[Firecrawl MCP](https://github.com/firecrawl/firecrawl-mcp-server)**
  - **描述**: 高级 Web 抓取/爬行 MCP，将网站转换为 LLM 可读的 Markdown。可处理动态内容。
  - **适用场景**: 从网络构建自定义数据集用于研究或分析。

- **[BrowserBase MCP](https://github.com/browserbase/mcp-server-browserbase)**
  - **描述**: 基于云的 AI 无头浏览器。允许在安全云环境中运行浏览器自动化脚本。
  - **适用场景**: 自动化需要完整浏览器上下文的复杂 Web 任务。

## <a id='coding-data'></a>💻 Coding & Data (编程与数据)

> 计算机科学教育、数据分析、可视化及技术设计。

- **[GitHub MCP Server](https://github.com/modelcontextprotocol/servers/tree/main/src/github)**
  - **描述**: GitHub 官方 MCP 协议集成。直接管理 Issues、PR 并分析仓库内容。
  - **适用场景**: 学习版本控制工作流，自动化代码审查。

- **[Playwright MCP](https://github.com/executeautomation/mcp-playwright)**
  - **描述**: 浏览器自动化服务器，让 Claude 能够浏览网页、测试站点并截图。
  - **适用场景**: 学习自动化测试，为项目抓取数据。

- **[Desktop Commander](https://github.com/wonderwhy-er/DesktopCommanderMCP)**
  - **描述**: 对桌面操作系统的深度控制。文件管理、终端执行和窗口控制。
  - **适用场景**: 学习命令行操作，自动化本地文件任务。

- **[Kubernetes MCP](https://github.com/fluxninja/mcp-kubernetes)**
  - **描述**: 与 Kubernetes 集群交互。查看日志、描述 Pod 和管理资源。
  - **适用场景**: 高级 DevOps 学习，集群监控。

- **[Sentry MCP](https://github.com/sentry-demos/mcp-server-sentry)**
  - **描述**: 检索 Sentry 问题和错误，以便结合上下文更快地调试应用程序。
  - **适用场景**: 学习生产环境调试和错误追踪。

- **[Octocode MCP](https://github.com/bgauryy/octocode-mcp)**
  - **描述**: 用于跨仓库进行语义代码研究的 MCP 服务器，允许你查找真实的代码实现。
  - **适用场景**: 计算机学生学习特定算法如何在开源项目中实现。

- **[MindsDB MCP](https://github.com/mindsdb/mindsdb)**
  - **描述**: 将 Claude 连接到 MindsDB 的数据库 AI 层。通过 SQL 训练模型并进行预测。
  - **适用场景**: 使用标准 SQL 语法学习应用机器学习。

- **[SQLite Explorer](https://github.com/hannesrudolph/sqlite-explorer-fastmcp-mcp-server)**
  - **描述**: 对 SQLite 数据库的安全只读访问。安全地探索表和架构。
  - **适用场景**: 初学者学习数据库架构和 SQL 查询。

- **[DBHub](https://github.com/bytebase/dbhub)**
  - **描述**: 数据库版本控制和协作工具。类似于数据库架构的 GitHub。
  - **适用场景**: 学习数据库迁移和架构生命周期管理。

- **[Neon MCP](https://github.com/neondatabase/mcp-server-neon)**
  - **描述**: 通过 MCP 进行 Serverless Postgres 操作。为测试扩展特定的数据分支。
  - **适用场景**: 现代云原生数据库开发。

## <a id='visual-presentation'></a>🎨 Visual & Presentation (视觉与演示 (PPT))

> PowerPoint 生成、视觉设计及演示工具。

- **[Git MCP](https://github.com/idosal/git-mcp)**
  - **描述**: 用于读取、搜索和操作本地 Git 仓库的工具。可视化历史记录和差异。
  - **适用场景**: 理解本地 Git 历史，调试复杂的合并问题。

- **[Docker MCP](https://github.com/ckreiling/mcp-server-docker)**
  - **描述**: 直接从 Claude 聊天中管理 Docker 容器、镜像和卷。
  - **适用场景**: 学习容器化概念，调试开发环境。

- **[Figma Context MCP](https://github.com/GLips/Figma-Context-MCP)**
  - **描述**: 从 Figma 设计中提取操作上下文，帮助 Claude 理解 UI 结构。
  - **适用场景**: 将设计模型转化为代码的学习。

- **[pptx-mcp](https://github.com/samos123/pptx-mcp)**
  - **描述**: 使用 Python PPTX 库创建幻灯片的简单 MCP 服务器。
  - **适用场景**: 从文本大纲快速生成基本幻灯片。

- **[mcp-server-okppt](https://github.com/NeekChaw/mcp-server-okppt)**
  - **描述**: 从 LLM 生成的 SVG 图像生成高质量 PPTX 幻灯片。确保矢量质量。
  - **适用场景**: 直接从 Claude 创建包含复杂图表的视觉设计幻灯片。

- **[Office-PowerPoint-MCP-Server](https://github.com/GongRzhe/Office-PowerPoint-MCP-Server)**
  - **描述**: 一个强大的 PowerPoint 操作 MCP 服务器。使用 python-pptx 创建、编辑和处理演示文稿。
  - **适用场景**: 自动化幻灯片创建、更新图表和批量编辑。

## <a id='productivity-career'></a>🧠 Productivity & Career (生产力与职业)

> 个人生产力、知识管理及职业发展。

- **[Slack MCP](https://github.com/modelcontextprotocol/servers/tree/main/src/slack)**
  - **描述**: 与 Slack 频道和消息交互。总结讨论串或起草回复。
  - **适用场景**: 学习团队沟通管理。

- **[Discord MCP](https://github.com/v-v-vishnu/discord-mcp-server)**
  - **描述**: 将 Claude 连接到 Discord 服务器。阅读消息并管理社区。
  - **适用场景**: 管理学习小组或学生社区。

- **[Linear MCP](https://github.com/jerhadf/linear-mcp-server)**
  - **描述**: 在 Linear 中管理任务和问题。从对话中创建工单。
  - **适用场景**: 学生学习敏捷/Scrum 项目管理。

- **[Todoist MCP](https://github.com/actionsOnGoogle/todoist-mcp-server)**
  - **描述**: 直接管理 Todoist 任务。添加作业或项目截止日期。
  - **适用场景**: 个人任务追踪和截止日期管理。

- **[NotebookLM MCP](https://github.com/PleasePrompto/notebooklm-mcp)**
  - **描述**: 将 Claude 与 Google NotebookLM 集成，允许你的 AI 代理通过有根据的回答来研究你自己的文档。
  - **适用场景**: 与你的课堂笔记或教科书“对话”，使用特定材料准备考试。

- **[Official Notion MCP Server](https://github.com/makenotion/notion-mcp-server)**
  - **描述**: Notion 官方 MCP 服务器。支持从 Claude 直接阅读页面、搜索数据库和创建内容。
  - **适用场景**: 整理课程表、更新研究 Wiki 或追踪申请进度。

- **[Anki MCP Server](https://github.com/nailuoGG/anki-mcp-server)**
  - **描述**: 通过 AnkiConnect 将 Claude 连接到 Anki。从对话生成抽认卡并自动添加到牌组。
  - **适用场景**: 将课堂笔记或词汇表即时转换为可复习的抽认卡。

- **[PDF Reader MCP](https://github.com/SylphxAI/pdf-reader-mcp)**
  - **描述**: 高性能 PDF 处理服务器，针对文本提取的速度（快 10 倍）和准确性进行了优化。
  - **适用场景**: 快速总结长篇教科书、提取数据表或转换作业。

- **[linkedin-mcp](https://github.com/pegasusheavy/linkedin-mcp)**
  - **描述**: Model Context Protocol (MCP) 服务器，用于管理 LinkedIn 个人资料、技能和教育背景。
  - **适用场景**: 个人品牌建设、自动更新职业履历。

---

*Last updated/最后更新: 2026-02-05*
