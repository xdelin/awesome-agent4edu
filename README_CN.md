# Awesome Education & Productivity MCP Servers

精选的 Model Context Protocol (MCP) 服务列表，专注于教育、校园生活、办公生产力、学术研究及个性化学习领域。

[English](./README.md) | [中文](./README_CN.md)

## 📖 目录

- [🏫 校园与学生管理 (Campus & Student Management)](#-校园与学生管理-campus--student-management)
  - [课程与 LMS 集成 (LMS Integration)](#课程与-lms-集成-lms-integration)
  - [学生事务与规划](#学生事务与规划)
- [📚 学习工具与助手 (Study Tools & Assistants)](#-学习工具与助手-study-tools--assistants)
  - [测验与备考 (Quiz & Exam Prep)](#测验与备考-quiz--exam-prep)
  - [语言学习 (Language Learning)](#语言学习-language-learning)
  - [学科专项 (Math & Science)](#学科专项-math--science)
- [🚀 生产力与知识库 (Productivity & Knowledge Base)](#-生产力与知识库-productivity--knowledge-base)
  - [Notion 集成](#notion-集成)
  - [Obsidian 笔记](#obsidian-笔记)
  - [总结与摘要工具](#总结与摘要工具)
- [🎨 教学设计 (Teaching Design)](#-教学设计-teaching-design)
  - [PPT 演示文稿](#ppt-演示文稿)
  - [动画制作](#动画制作)
  - [图表与可视化](#图表与可视化)
- [📝 学术写作 (Academic Writing)](#-学术写作-academic-writing)
  - [论文逻辑与结构](#论文逻辑与结构-1)
  - [排版与文献管理](#排版与文献管理-1)
- [💻 编程教育 (Coding Education)](#-编程教育-coding-education)

---

## 🏫 校园与学生管理 (Campus & Student Management)

连接学校系统，管理课程表和作业。

### 课程与 LMS 集成 (LMS Integration)
直接在 AI 中访问 Canvas, Brightspace 等学习管理系统。

- **[makenotion/notion-mcp-server](https://github.com/makenotion/notion-mcp-server)** (Official)
  通过官方 Notion MCP 管理学习笔记和任务看板。

- **[akshsgaur/CMUCanvasMCPSErver](https://github.com/akshsgaur/CMUCanvasMCPSErver)**
  Canvas LMS 集成，让学生通过 AI 对话查询课程、作业截止日期和待办事项（专为 CMU 设计但可适配）。

- **[pranav-vijayananth/brightspace-mcp-server](https://github.com/pranav-vijayananth/brightspace-mcp-server)**
  连接 Brightspace 账户，获取课程信息。

- **[MisterCommand/mcp-canvas-lms](https://github.com/MisterCommand/mcp-canvas-lms)**
  TypeScript 编写的 Canvas LMS 访问工具。

- **[markmusic27/stanford-mcp](https://github.com/markmusic27/stanford-mcp)**
  斯坦福课程目录查询工具，包含课程评分查询。

### 学生事务与规划
管理学业进度和校园生活。

- **[tejpalvirk/student](https://github.com/tejpalvirk/student)**
  基于知识图谱的学生上下文管理，关联课程、作业和考试信息。

- **[Nitishroy-7033/Student_MCP_Server](https://github.com/Nitishroy-7033/Student_MCP_Server)**
  模块化的学院管理系统接口，涵盖费用、学科和考试结构。

- **[abhijeetgiram/mcp-student-leave-copilot-ab](https://github.com/abhijeetgiram/mcp-student-leave-copilot-ab)**
  集成 GitHub Copilot 的学生请假申请管理工具。

---

## 📚 学习工具与助手 (Study Tools & Assistants)

### 测验与备考 (Quiz & Exam Prep)
自动生成测验，辅助备考。

- **[cardea-mcp/ExamPrepAgent](https://github.com/cardea-mcp/ExamPrepAgent)**
  智能备考伴侣，利用 LLM 协助学生准备各类认证考试。

- **[RidiculousBuffal/quizAgent](https://github.com/RidiculousBuffal/quizAgent)**
  AI 增强的测验代理，支持 SSO 登录。

- **[kkjzio/mcp-quiz-generator](https://github.com/kkjzio/mcp-quiz-generator)**
  根据需求生成 Markdown 格式的题目，并可转换为 HTML 和 Word 测验文件。

- **[sheikhcoders/interleaved-learning-mcp](https://github.com/sheikhcoders/interleaved-learning-mcp)**
  综合学习工具：包含学习日程安排、测验生成、闪卡洗牌和学习追踪。

- **[PratyayRajak/todopomo-mcp](https://github.com/PratyayRajak/todopomo-mcp)**
  结合待办事项与番茄工作法（Pomodoro）的学习辅助工具。

### 语言学习 (Language Learning)
辅助外语习得。

- **[w41ch0ng/MandarinMCP](https://github.com/w41ch0ng/MandarinMCP)**
  中文（普通话）学习工具，覆盖 HSK 1-6 级词汇，支持进度追踪和 Anki 导出。

- **[satish-kori-27/swedish-language-mcp-server](https://github.com/satish-kori-27/swedish-language-mcp-server)**
  瑞典语学习辅助。

- **[Jonathan-Racaud/Glossa](https://github.com/Jonathan-Racaud/Glossa)**
  通用的语言学习 MCP 服务器。

### 学科专项 (Math & Science)
特定领域的深度学习工具。

- **[clouatre-labs/math-mcp-learning-server](https://github.com/clouatre-labs/math-mcp-learning-server)**
  支持矩阵运算、数据可视化的数学教学服务器，带有持久化工作区。

- **[Augmented-Nature/GeneOntology-MCP-Server](https://github.com/Augmented-Nature/GeneOntology-MCP-Server)**
  生物信息学工具，访问基因本体论 (GO) 数据，进行本体分析和富集研究。

- **[cyanheads/clinicaltrialsgov-mcp-server](https://github.com/cyanheads/clinicaltrialsgov-mcp-server)**
  搜索和检索 ClinicalTrials.gov 的临床试验数据，适合医学研究。

- **[ebrandler/mcportfolio](https://github.com/ebrandler/mcportfolio)**
  金融学学生的利器，提供 9 种专业工具进行投资组合优化和量化分析。

---

## 🚀 生产力与知识库 (Productivity & Knowledge Base)

### Notion 集成
最流行的多合一笔记工具集成。

- **[makenotion/notion-mcp-server](https://github.com/makenotion/notion-mcp-server)** 🔥
  **官方** Notion MCP 服务器。

- **[suekou/mcp-notion-server](https://github.com/suekou/mcp-notion-server)**
  社区广受欢迎的实现，功能丰富。

- **[Badhansen/notion-mcp](https://github.com/Badhansen/notion-mcp)**
  专注于简单的个人待办事项管理。

### Obsidian 笔记
本地优先的知识管理与第二大脑。

- **[MarkusPfundstein/mcp-obsidian](https://github.com/MarkusPfundstein/mcp-obsidian)** 🔥
  (2.7k+ Stars) 通过 Local REST API 插件与 Obsidian 交互，极为流行。

- **[cyanheads/obsidian-mcp-server](https://github.com/cyanheads/obsidian-mcp-server)**
  功能全面的工具套件，支持读取、写入、搜索笔记以及管理 Frontmatter。

- **[entanglr/zettelkasten-mcp](https://github.com/entanglr/zettelkasten-mcp)**
  专门实现 **卢曼卡片盒笔记法 (Zettelkasten)** 的工具，帮助根据原子笔记建立知识连接。

### 总结与摘要工具
快速消化视频与文档内容。

- **[5hivanand/mcp-youtube-video-summary](https://github.com/5hivanand/mcp-youtube-video-summary)**
  获取 YouTube 视频字幕并在本地生成可定制的摘要。

- **[dEitY719/mcp-youtube-summary](https://github.com/dEitY719/mcp-youtube-summary)**
  基于 FastMCP 构建的 YouTube 视频智能摘要工具。

---

## 🎨 教学设计 (Teaching Design)

### PPT 演示文稿
- **[samos123/pptx-mcp](https://github.com/samos123/pptx-mcp)**
  自然语言生成 PPT，基于 python-pptx。
- **[Weichenleeeee123/ppt-mcp-server](https://github.com/Weichenleeeee123/ppt-mcp-server)**
  全面的 PPT 编辑功能。

### 动画制作
- **[wstcpyt/manim-mcp](https://github.com/wstcpyt/manim-mcp)**
  基于 Manim 引擎的数学动画制作。
- **[ampersante/spine2d-animation-mcp](https://github.com/ampersante/spine2d-animation-mcp)**
  Spine2D 骨骼动画工具。

### 图表与可视化
- **[skmprb/md-mermaid-chart-pdf-mcp](https://github.com/skmprb/md-mermaid-chart-pdf-mcp)**
  Mermaid 图表转 PDF。
- **[TakanariShimbo/quickchart-mcp-server](https://github.com/TakanariShimbo/quickchart-mcp-server)**
  调用 QuickChart 生成统计图。

---

## 📝 学术写作 (Academic Writing)

### 论文逻辑与结构
- **[shaneholloman/mcp-knowledge-graph](https://github.com/shaneholloman/mcp-knowledge-graph)**
  构建论文逻辑的本地知识图谱。

### 排版与文献管理
- **[cyanheads/docwriter-mcp-server](https://github.com/cyanheads/docwriter-mcp-server)**
  LaTeX 文档结构化创建与 Biber 引用管理。
- **[devroopsaha744/TexMCP](https://github.com/devroopsaha744/TexMCP)**
  LaTeX 代码片段渲染微服务。
- **[zongmin-yu/sqlite-literature-management](https://github.com/zongmin-yu/sqlite-literature-management-fastmcp-mcp-server)**
  基于 SQLite 的文献源管理。

---

## 💻 编程教育 (Coding Education)

- **[Shen-Ming-Hong/singular-blockly](https://github.com/Shen-Ming-Hong/singular-blockly)**
  可视化编程工具（Blockly），支持 Arduino 和 MicroPython，适合 STEM 教育。

- **[parth012001/Codlab-ai-platform](https://github.com/parth012001/Codlab-ai-platform)**
  AI 驱动的编程教育平台，包含挑战、自动评分和面试模拟。
