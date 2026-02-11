# Awesome Education AI: MCP, Skills & LLM Apps

一份精心策划的教育 AI 资源列表，涵盖 **Model Context Protocol (MCP)** 服务器、**Claude Skills**、**教育类 LLM** 和 **Agent 框架**，专注于提升学术研究、教学效率和个性化学习体验。

[English](./README.md) | [中文](./README_CN.md)

## 目录

- [智能导学](#intelligent-tutoring)
- [数理科学](#math-science)
- [生物与医学](#biology-medicine)
- [化学与材料](#chemistry-materials)
- [艺术、人文与社会科学](#arts-humanities-social-sciences)
- [计算机科学与工程](#computer-science-engineering)
- [数据与分析](#data-analysis)
- [视觉与演示 (PPT)](#visual-presentation)
- [学术与写作](#academic-writing)
- [笔记与知识库](#notes-knowledge-base)
- [职业规划与生产力](#career-productivity)
- [元技能](#meta-skills)
- [商业化教育 AI](#commercial-education-ai)
- [教育大模型](#education-llms)
- [Agent 框架与教育应用](#agent-frameworks-educational-applications)

---

## <a id='intelligent-tutoring'></a>Intelligent Tutoring (智能导学)

> AI 导师、教育评估及个性化学习助手。

- **[claude-educational-ai-skills](https://github.com/ChatGPT3a01/claude-educational-ai-skills)**
  - **描述**: 专为教育研究者、教师和数字学习平台开发者设计的 AI 技能集。包括课程设计、评估生成等功能。
  - **适用场景**: 备课、教学法研究、自动化出题。

- **[EdTech-AI-Coach](https://github.com/VaibhavKanojia3773/EdTech-AI-Coach)**
  - **描述**: 基于 Python 的工具，利用 Claude 分析学生测试数据 (JSON)，生成个性化、激励性的反馈报告 (PDF)。
  - **适用场景**: 教师不仅需要分数，还需要为每个学生提供具体的改进建议。

- **[ELA-Tutor-Agent](https://github.com/cozette-brown/ELA-Tutor-Agent-System-Prompting-Experiment)**
  - **描述**: 一个实验性的 System Prompt，将 Claude 转变为专业的英语语言艺术 (ELA) 导师。
  - **适用场景**: 针对 6-8 年级学生的阅读理解、写作指导和文学分析对话。

- **[quiz-master](https://github.com/tanmaypatil/quiz-master)**
  - **描述**: 利用 GenAI 和 Claude Code 生成通用知识测验的应用程序。
  - **适用场景**: 课堂快速测验、自我检测、Flashcard 制作。

- **[apprenticemode](https://github.com/yavzius/apprenticemode)**
  - **描述**: 一个独特的技能，用于检测用户的"知识缺口" (Knowledge Gaps)。通过你的 Prompt 反向推断你不知道什么。
  - **适用场景**: 元认知学习、查漏补缺、自我评估编程水平。

- **[EduBase](https://github.com/EduBase/MCP)**
  - **描述**: 与 EduBase 交互，这是一个具有测验和考试管理功能的电子学习平台。
  - **适用场景**: 参加测验、管理课程。

- **[Canvas LMS](https://github.com/ahnopologetic/canvas-lms-mcp)**
  - **描述**: 连接到 Canvas 学习管理系统 (LMS)，以访问课程数据、作业和公告。
  - **适用场景**: 检查成绩、管理作业。

- **[Educational Tutor MCP](https://github.com/iamisaackn/AI-TutorMCP)**
  - **描述**: 基于模型上下文协议 (MCP) 构建的 AI 教育工具包。
  - **适用场景**: 个性化辅导，教育辅助。

## <a id='math-science'></a>Math & Science (数理科学)

> 数学、物理、STEM 研究及科学工具。

- **[Kvante](https://github.com/MaxPowerBoy1980/Kvante)**
  - **描述**: AI 数学助手，帮助学生通过 OCR 识别手写数学题，并提供**分步解答**（而不是直接给出答案）。
  - **适用场景**: 数学辅导、家庭作业辅助、理解解题逻辑。

- **[MCPR (R Language)](https://github.com/phisanti/MCPR)**
  - **描述**: 允许 AI 代理参与交互式实时 R 会话以进行统计分析。
  - **适用场景**: 统计学课程、数据分析。

- **[HOCR (Math Capture)](https://github.com/Wooonster/hocr_mcp_server)**
  - **描述**: 从图像中提取手写数学公式并将其转换为 LaTeX。
  - **适用场景**: 数字化数学笔记、家庭作业帮助。

- **[Stella MCP](https://github.com/bradleylab/stella-mcp)**
  - **描述**: 创建并验证用于科学模拟的系统动力学模型 (.stmx 文件)。
  - **适用场景**: 系统思维、科学建模课程。

- **[Calculator](https://github.com/githejie/mcp-server-calculator)**
  - **描述**: 一个精确的数值计算器，防止 LLM 数学幻觉。
  - **适用场景**: 数学作业检查、基础计算。

- **[Manim Composer](https://github.com/openclaw/skills/tree/main/skills/inclinedadarsh/manim-composer)**
  - **描述**: 使用 Manim 引擎创建数学动画和可视化效果。
  - **适用场景**: 制作数学教程、可视化复杂概念。

- **[Quantum Lab](https://github.com/openclaw/skills/tree/main/skills/bramdo/quantum-lab)**
  - **描述**: 运行量子物理模拟和量子实验室实验的 Python 脚本。
  - **适用场景**: 物理学生、量子计算研究。

- **[NetKet MCP](https://github.com/pathintegral-institute/mcp.science)**
  - **描述**: 与 NetKet 交互，使用神经网络模拟量子多体系统。
  - **适用场景**: 量子物理研究、多体问题模拟。

- **[Fermat Math Engine](https://github.com/abhiphile/fermat-mcp)**
  - **描述**: 将 SymPy、NumPy 和 Matplotlib 统一到一个服务器中，用于符号代数、数值计算和可视化。
  - **适用场景**: 复杂数学问题、物理模拟。

- **[Phys-MCP](https://github.com/BlinkZer0/Phys-MCP)**
  - **描述**: 专注于物理的计算代数系统 (CAS)，用于解决具有单位感知计算的物理问题。
  - **适用场景**: 物理作业、单位换算、力学问题。

- **[Bullet MCP Server](https://github.com/devgabrielsborges/bullet-mcp-server)**
  - **描述**: 与 PyBullet 物理引擎交互，设置和控制物理模拟。
  - **适用场景**: 机器人模拟、物理演示、碰撞检测。

- **[NASA MCP Server](https://github.com/ProgramComputer/NASA-MCP-server)**
  - **描述**: 访问 NASA API 获取天文学数据，包括火星车照片和地球观测图像。
  - **适用场景**: 天文学教育、太空研究、探索宇宙。

- **[Wolfram Alpha MCP](https://github.com/SecretiveShell/MCP-wolfram-alpha)**
  - **描述**: 访问 Wolfram Alpha 的计算智能，用于数学、物理和数据查询。
  - **适用场景**: 检查复杂计算、科学事实检索。

- **[WolframAlpha](https://github.com/mcp-servers/wolfram-alpha)**
  - **描述**: 与 WolframAlpha 集成，提供计算知识和科学数据。
  - **适用场景**: 解决复杂的数学问题、物理计算、化学方程式。

- **[Manim](https://github.com/3b1b/manim)**
  - **描述**: 用于解释数学视频的动画引擎。通过编程创建精确的动画。
  - **适用场景**: 可视化数学概念，制作教育视频。

- **[OpenFOAM](https://github.com/OpenFOAM/OpenFOAM-dev)**
  - **描述**: 免费开源的计算流体动力学 (CFD) 软件。
  - **适用场景**: 工程仿真、流体动力学研究。

- **[Open-Meteo MCP](https://github.com/open-meteo/open-meteo)**
  - **描述**: 天气数据 API 集成。
  - **适用场景**: 地球科学数据收集、气候研究。

## <a id='biology-medicine'></a>Biology & Medicine (生物与医学)

> 基因组学、蛋白质分析、医学研究及生物信息学工具。

- **[Bio Agents](https://github.com/dogeplusplus/bio-agents-mcp)**
  - **描述**: 集成 Protein Data Bank (PDB) 和 ChEMBL，用于查询药物数据和蛋白质结构。
  - **适用场景**: 药物发现、蛋白质可视化、生物活性研究。

- **[Ensembl MCP](https://github.com/effieklimi/ensembl-mcp-server)**
  - **描述**: 访问 Ensembl 数据库，获取脊椎动物基因组、基因查找和序列比对数据。
  - **适用场景**: 基因组研究、基因序列分析。

- **[PubMed Search](https://github.com/andybrandt/mcp-simple-pubmed)**
  - **描述**: 从 PubMed 搜索和阅读医学及生命科学论文。
  - **适用场景**: 医学研究、生物学学习。

- **[ClinicalTrials.gov MCP](https://github.com/cyanheads/clinicaltrialsgov-mcp-server)**
  - **描述**: 查询 ClinicalTrials.gov API 获取临床研究记录和资格标准。
  - **适用场景**: 医学研究、查找临床试验、检查研究详情。

- **[AlphaFold MCP Server](https://github.com/Augmented-Nature/AlphaFold-MCP-Server)**
  - **描述**: 从 AlphaFold 蛋白质结构数据库检索预测的 3D 蛋白质结构。
  - **适用场景**: 分子生物学、药物设计、理解蛋白质功能。

- **[Reactome MCP Server](https://github.com/Augmented-Nature/Reactome-MCP-Server)**
  - **描述**: 与 Reactome 通路数据库交互，用于系统生物学分析。
  - **适用场景**: 研究生物通路、代谢研究。

- **[BioThings MCP](https://github.com/longevity-genie/biothings-mcp)**
  - **描述**: 与 BioThings API 交互以检索基因、变异、药物及化学数据。
  - **适用场景**: 生物信息学研究、遗传学课程。

- **[UCSC Genome Browser](https://github.com/hlydecker/ucsc-genome-mcp)**
  - **描述**: 与 UCSC 基因组浏览器 API 交互以查找基因组和染色体。
  - **适用场景**: 基因组学教育、染色体可视化。

- **[BioContext](https://github.com/biocontext/biocontext)**
  - **描述**: 用于在文本中查找和可视化生物学内容的工具/语境。
  - **适用场景**: 生物学研究辅助。

- **[Ensembl REST](https://github.com/Ensembl/ensembl-rest)**
  - **描述**: 通过 REST API 访问 Ensembl 基因组数据。
  - **适用场景**: 基因组数据检索和分析。

## <a id='chemistry-materials'></a>Chemistry & Materials (化学与材料)

> 分子建模、化学数据、材料科学及晶体结构。

- **[Materials Project MCP](https://github.com/pathintegral-institute/mcp.science)**
  - **描述**: 通过 Materials Project API 访问电池、分子和晶体结构数据。
  - **适用场景**: 材料发现、电池研究、晶体结构分析。

- **[MatMCP](https://github.com/ZuchGuillotine/MatMCP)**
  - **描述**: 用于材料科学查询和化学性质检索的专用服务器。
  - **适用场景**: 通用化学查询、材料性质查找。

- **[RDKit MCP Server](https://github.com/tandemai-inc/rdkit-mcp-server)**
  - **描述**: 封装 RDKit 库，用于 SMILES 生成和分子修改等化学信息学任务。
  - **适用场景**: 化学信息学、分子设计、化学性质计算。

- **[ChemMCP](https://github.com/OSU-NLP-Group/ChemMCP)**
  - **描述**: 一个用于分子分析和交互的综合化学工具包。
  - **适用场景**: 复杂化学问题、分子分析。

- **[XTB MCP Server](https://github.com/PhelanShao/xtb-mcp-server)**
  - **描述**: 用于生成扩展紧束缚 (xTB) 量子化学计算输入文件的工具。
  - **适用场景**: 量子化学、自动化计算工作流。

- **[Catalysis Hub](https://github.com/QuentinCody/catalysishub-mcp-server)**
  - **描述**: 从 Catalysis Hub 数据库搜索和检索科学数据。
  - **适用场景**: 化学研究、催化实验。

## <a id='arts-humanities-social-sciences'></a>Arts, Humanities & Social Sciences (艺术、人文与社会科学)

> 历史、地理、社会学及艺术。

- **[Google Maps MCP](https://github.com/cablate/mcp-google-map)**
  - **描述**: 通过 Google Maps 搜索地点、计算距离并获取位置详情。
  - **适用场景**: 地理课程、城市规划、旅行计划。

- **[DeepL API](https://github.com/DeepLcom/deepl-python)**
  - **描述**: 高质量机器翻译库。
  - **适用场景**: 语言学习、翻译研究。

- **[Music21](https://github.com/cuthbertLab/music21)**
  - **描述**: 计算机辅助音乐学工具包。
  - **适用场景**: 乐理分析、作曲。

- **[MuseScore](https://github.com/musescore/MuseScore)**
  - **描述**: 开源音乐制谱软件。
  - **适用场景**: 音乐教育、乐谱创作。

- **[Google Maps MCP](https://github.com/modelcontextprotocol/servers/tree/main/src/google-maps)**
  - **描述**: 访问地理数据和定位服务。
  - **适用场景**: 地理教育、文化研究。

## <a id='computer-science-engineering'></a>Computer Science & Engineering (计算机科学与工程)

> 开发工具、编程助手、数据库管理及 DevOps。 包含具身智能与 VR 环境。

- **[Anthropic Web App Testing](https://github.com/anthropics/anthropic-skills/tree/main/skills/webapp-testing)**
  - **描述**: 使用 Playwright/Puppeteer 概念的自动化 Web 应用程序测试技能。
  - **适用场景**: 学习 QA 自动化、回归测试。

- **[paper2code-skill](https://github.com/issol14/paper2code-skill)**
  - **描述**: 一个通过 Claude Code 将研究论文转换为可执行代码的技能。
  - **适用场景**: 计算机科学/AI 研究生复现论文算法。

- **[Git MCP Server](https://github.com/modelcontextprotocol/servers/tree/main/src/git)**
  - **描述**: 用于在本地读取、搜索和操作 Git 仓库的工具。
  - **适用场景**: 本地版本控制管理、提交更改和查看历史记录。

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

- **[awesome-vibe-coding](https://github.com/adriannoes/awesome-vibe-coding)**
  - **描述**: 为非技术背景的构建者提供的资源库。包含 Cursor/Claude 规则 (Rules)、AI Prompts 和自动化模板。
  - **适用场景**: 编程初学者、"Vibe Coding" (直觉式编程)、快速原型开发。

- **[coding-agent-skills](https://github.com/Miaoge-Ge/coding-agent-skills)**
  - **描述**: 通用的专家系统提示词库，包含中英双语的角色设定（如深度学习专家、架构师）。
  - **适用场景**: 计算机专业学生进行深度技术学习、代码审查。

- **[test-generator](https://github.com/mkdirrobert/test-generator)**
  - **描述**: Claude Code 的通用测试生成技能，支持自我学习并覆盖 9 种编程语言。
  - **适用场景**: 学习单元测试编写、提高代码质量。

- **[mlx-dev-skill](https://github.com/tkwn2080/mlx-dev-skill)**
  - **描述**: 专门用于在 Apple Silicon 上编写符合习惯的 Apple MLX 代码的技能。
  - **适用场景**: 学习 AI/ML 开发、Mac 用户。

- **[prompt-builder](https://github.com/Sorata-sensei/prompt-builder)**
  - **描述**: 一个交互式表单，帮助非技术学生生成高质量的 Claude/ChatGPT 提示词，用于构建 PWA 应用。
  - **适用场景**: 编程入门课程、通过自然语言构建软件的教学。

- **[Octocode MCP](https://github.com/bgauryy/octocode-mcp)**
  - **描述**: 用于跨仓库进行语义代码研究的 MCP 服务器，允许你查找真实的代码实现。
  - **适用场景**: 计算机学生学习特定算法如何在开源项目中实现。

- **[MCP JS](https://github.com/r33drichards/mcp-js)**
  - **描述**: 一个使用 V8 引擎的安全 JavaScript 执行沙箱，用于运行本地生成的代码。
  - **适用场景**: 学习 JS、安全测试代码片段。

- **[Neovim MCP](https://github.com/linw1995/nvim-mcp)**
  - **描述**: 允许 LLM 直接与 Neovim 实例交互，读取缓冲区和编辑文本。
  - **适用场景**: 学习 Vim 命令、终端 AI 辅助编程。

- **[VS Code MCP Server](https://github.com/juehang/vscode-mcp-server)**
  - **描述**: 允许外部代理控制 VS Code，读取文件树并查看 linter 错误。
  - **适用场景**: AI 结对编程、修复 Lint 错误。

- **[Smart Contracts Wizard](https://github.com/OpenZeppelin/contracts-wizard)**
  - **描述**: 基于 OpenZeppelin 模板生成安全的智能合约。
  - **适用场景**: 学习 Solidity 和 Web3 开发、合约审计。

- **[Go Language Server (gopls)](https://github.com/hloiseaufcms/mcp-gopls)**
  - **描述**: 与 Go 语言服务器协议交互，解释代码并查找定义。
  - **适用场景**: 学习 Go 语言、理解代码库结构。

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

- **[a11y-specialist-skills](https://github.com/masuP9/a11y-specialist-skills)**
  - **描述**: Claude Code 的无障碍 (Accessibility) 专家技能插件。
  - **适用场景**: 学习 WCAG 标准、进行网站无障碍审计。

- **[Docker MCP](https://github.com/ckreiling/mcp-server-docker)**
  - **描述**: 直接从 Claude 聊天中管理 Docker 容器、镜像和卷。
  - **适用场景**: 学习容器化概念，调试开发环境。

- **[LeetCode MCP](https://github.com/jinzcdev/leetcode-mcp-server)**
  - **描述**: 允许 AI 搜索和检索 LeetCode 编程题目、解答和提交数据。
  - **适用场景**: 算法学习、面试准备。

- **[OpenAPI Explorer](https://github.com/kadykov/mcp-openapi-schema-explorer)**
  - **描述**: 通过与 OpenAPI 规范交互来探索和理解 API 结构。
  - **适用场景**: 学习 API 设计、集成测试。

- **[Code Mentor](https://github.com/openclaw/skills/tree/main/skills/samuelkahessay/code-mentor)**
  - **描述**: 面向所有级别经验的综合型 AI 编程导师。
  - **适用场景**: 学习编程概念、代码审查、调试辅助。

- **[Backend Patterns](https://github.com/openclaw/skills/tree/main/skills/charmmm718/backend-patterns)**
  - **描述**: 后端架构模式、API 设计和数据库管理指南。
  - **适用场景**: 学习系统设计、后端开发最佳实践。

- **[Git Essentials](https://github.com/openclaw/skills/tree/main/skills/arnarsson/git-essentials)**
  - **描述**: 掌握版本控制所必备的 Git 命令和工作流。
  - **适用场景**: 学习 Git 版本控制、掌握分支和合并。

- **[Docker Essentials](https://github.com/openclaw/skills/tree/main/skills/arnarsson/docker-essentials)**
  - **描述**: 容器化所必备的 Docker 命令和工作流。
  - **适用场景**: 学习 Docker、容器管理、DevOps 基础。

- **[Postgres MCP](https://github.com/modelcontextprotocol/servers/tree/main/src/postgres)**
  - **描述**: 官方 Postgres MCP 服务器。检查模式并运行查询。
  - **适用场景**: 数据库管理教育、SQL 练习。

- **[LeetCode](https://leetcode.com)**
  - **描述**: 学习编程和准备技术面试的平台。
  - **适用场景**: 算法练习、数据结构学习。

- **[Kubernetes CLI](https://github.com/kubernetes/kubernetes)**
  - **描述**: 生产级容器编排。
  - **适用场景**: DevOps 教育、云工程。

- **[Chrome DevTools Protocol](https://github.com/ChromeDevTools/devtools-protocol)**
  - **描述**: 用于调试和分析 Web 应用程序的工具。
  - **适用场景**: Web 开发教育、前端调试。

- **[Git MCP](https://github.com/idosal/git-mcp)**
  - **描述**: 用于读取、搜索和操作本地 Git 仓库的工具。可视化历史记录和差异。
  - **适用场景**: 理解本地 Git 历史，调试复杂的合并问题。

- **[AI2-THOR](https://github.com/allenai/ai2thor)**
  - **描述**: 用于视觉 AI 的交互式 3D 环境，允许智能体导航并与物体交互。
  - **适用场景**: 机器人训练、视觉推理教育。

- **[ThreeDWorld (TDW)](https://github.com/threedworld-mit/tdw)**
  - **描述**: 一个用于高保真感官数据交互式物理模拟的平台。
  - **适用场景**: 物理模拟、空间智能研究。

- **[UnrealZoo](https://github.com/UnrealZoo/UnrealZoo)**
  - **描述**: 基于虚幻引擎 5 开发和评估具身智能体的平台 (ICCV 2025)。
  - **适用场景**: 照片级仿真、具身智能体训练。

## <a id='data-analysis'></a>Data & Analysis (数据与分析)

> 数据科学、网络抓取、分析及大数据工具。

- **[Anthropic XLSX Skill](https://github.com/anthropics/anthropic-skills/tree/main/skills/xlsx)**
  - **描述**: 用于分析和操作 Excel 电子表格的官方技能。计算数据并生成见解。
  - **适用场景**: 数据分析、财务建模、实验室结果处理。

- **[Firecrawl MCP](https://github.com/firecrawl/firecrawl-mcp-server)**
  - **描述**: 高级 Web 抓取/爬行 MCP，将网站转换为 LLM 可读的 Markdown。可处理动态内容。
  - **适用场景**: 从网络构建自定义数据集用于研究或分析。

- **[BrowserBase MCP](https://github.com/browserbase/mcp-server-browserbase)**
  - **描述**: 基于云的 AI 无头浏览器。允许在安全云环境中运行浏览器自动化脚本。
  - **适用场景**: 自动化需要完整浏览器上下文的复杂 Web 任务。

- **[GenAI Toolbox](https://github.com/googleapis/genai-toolbox)**
  - **描述**: Google 官方工具，用于 BigQuery、AlloyDB 等。使用自然语言与海量数据集交互。
  - **适用场景**: 学习云数据仓库和大数据分析。

- **[claude-skill-data-cleaner](https://github.com/brook-miller/claude-skill-data-cleaner)**
  - **描述**: 让 Claude 检查 CSV 文件并创建一个 Python 脚本来清理数据并准备分析。
  - **适用场景**: 数据科学入门、清洗脏数据。

- **[activitywatch-analysis-skill](https://github.com/BayramAnnakov/activitywatch-analysis-skill)**
  - **描述**: 使用 ActivityWatch 数据进行每周生产力分析。计算专注度分数，检测"死亡循环"并生成见解。
  - **适用场景**: 自我量化 (Quantified Self)、提高学习效率。

- **[Jupyter Notebook MCP](https://github.com/datalayer/jupyter-mcp-server)**
  - **描述**: 与 Jupyter Notebook 交互，主要用于数据分析任务。
  - **适用场景**: 数据科学工作流、交互式编程。

- **[NetworkX Graph Analysis](https://github.com/Bright-L01/networkx-mcp-server)**
  - **描述**: 使用 NetworkX 进行图分析和可视化（中心性、社区检测）。
  - **适用场景**: 网络科学、社交网络分析。

- **[Fetch MCP](https://github.com/zcaceres/fetch-mcp)**
  - **描述**: 灵活地从网络获取 JSON、文本和 HTML 数据。用于收集数据集素材。
  - **适用场景**: 网络抓取、数据收集。

- **[Vizro (Data Viz)](https://github.com/mckinsey/vizro/tree/main/vizro-mcp)**
  - **描述**: 用于创建经过验证且可维护的数据图表和仪表板的工具。
  - **适用场景**: 学习数据可视化、创建仪表板。

- **[Kaggle MCP](https://github.com/KrishnaPramodParupudi/kaggle-mcp-server)**
  - **描述**: 直接浏览 Kaggle 竞赛、排行榜、数据集和 Kernels。
  - **适用场景**: 数据科学练习、查找数据集。

## <a id='visual-presentation'></a>Visual & Presentation (视觉与演示 (PPT))

> PowerPoint 生成、数据可视化、UI 设计及视觉工具。

- **[Anthropic PPTX Skill](https://github.com/anthropics/anthropic-skills/tree/main/skills/pptx)**
  - **描述**: 用于以编程方式生成 PowerPoint 演示文稿的官方技能。
  - **适用场景**: 创建学术演示、从笔记生成讲义幻灯片。

- **[Anthropic Algorithmic Art](https://github.com/anthropics/anthropic-skills/tree/main/skills/algorithmic-art)**
  - **描述**: 使用数学函数生成基于 SVG 的算法艺术。
  - **适用场景**: 创意编程教育、生成抽象视觉效果。

- **[chartjs-expert](https://github.com/sjnims/chartjs-expert)**
  - **描述**: Chart.js 专家技能，支持交互式代码生成和多框架集成 (React/Vue 等)。
  - **适用场景**: 前端开发者学习图表库集成。

- **[census-demographics-skill](https://github.com/baofeng-dong/census-demographics-skill)**
  - **描述**: 为 Claude Code 设计的数据可视化技能，可创建展示教育、收入和人口统计数据的地图。
  - **适用场景**: 社会学研究、地理信息系统 (GIS) 课程、数据新闻学。

- **[Observable-Plot-Claude-Skill](https://github.com/jkoets/Observable-Plot-Claude-Skill)**
  - **描述**: 指导使用 Observable Plot (JavaScript 库) 创建数据可视化的技能。
  - **适用场景**: 学习现代 Web 数据可视化。

- **[Figma Context MCP](https://github.com/GLips/Figma-Context-MCP)**
  - **描述**: 从 Figma 设计中提取操作上下文，帮助 Claude 理解 UI 结构。
  - **适用场景**: 将设计模型转化为代码的学习。

- **[neurodivergent-visual-org](https://github.com/JackReis/neurodivergent-visual-org)**
  - **描述**: 为神经多样性思维 (Neurodivergent thinking) 设计的视觉组织工具 - 提供富有同理心的任务分解和决策树。
  - **适用场景**: ADHD 学生的时间管理、辅助执行功能障碍。

- **[pptx-mcp](https://github.com/samos123/pptx-mcp)**
  - **描述**: 使用 Python PPTX 库创建幻灯片的简单 MCP 服务器。
  - **适用场景**: 从文本大纲快速生成基本幻灯片。

- **[mcp-server-okppt](https://github.com/NeekChaw/mcp-server-okppt)**
  - **描述**: 从 LLM 生成的 SVG 图像生成高质量 PPTX 幻灯片。确保矢量质量。
  - **适用场景**: 直接从 Claude 创建包含复杂图表的视觉设计幻灯片。

- **[Office-PowerPoint-MCP-Server](https://github.com/GongRzhe/Office-PowerPoint-MCP-Server)**
  - **描述**: 一个强大的 PowerPoint 操作 MCP 服务器。使用 python-pptx 创建、编辑和处理演示文稿。
  - **适用场景**: 自动化幻灯片创建、更新图表和批量编辑。

- **[claude-dolphin](https://github.com/nyldn/claude-dolphin)**
  - **描述**: 综合 UX/UI 设计系统技能 - 包含 WCAG 2.2 AA 审计、设计一致性检查和 Refactoring UI 原则。
  - **适用场景**: UI/UX 设计学生、全栈开发者提升设计感。

- **[Markmap](https://github.com/jinzcdev/markmap-mcp-server)**
  - **描述**: 将 Markdown 内容转换为交互式思维导图，用于视觉学习。
  - **适用场景**: 创建学习大纲、头脑风暴会议。

- **[Excalidraw MCP App](https://github.com/antonpk1/excalidraw-mcp-app)**
  - **描述**: Excalidraw 官方 MCP 应用。以编程方式创建和编辑精美的手绘风格图表。
  - **适用场景**: 创建解释性图表、系统架构白板演示。

- **[MCP Mermaid](https://github.com/hustcc/mcp-mermaid)**
  - **描述**: 使用 AI 动态生成 Mermaid 图表。支持流程图、时序图等多种类型。
  - **适用场景**: 可视化流程、创建代码文档图表。

- **[Leap MCP](https://github.com/sid-thephysicskid/leap-mcp)**
  - **描述**: 将任何主题转化为带有 AI 旁白和 3Blue1Brown 风格动画的解说视频。
  - **适用场景**: 制作教育视频、可视化数学/物理概念。

- **[Manim MCP](https://github.com/wstcpyt/manim-mcp)**
  - **描述**: Manim 数学动画引擎的 MCP 服务器。生成精确的数学动画代码。
  - **适用场景**: 数学教育、动态几何演示。

- **[Kroki MCP](https://github.com/utain/kroki-mcp)**
  - **描述**: 通过 Kroki 将文本图表定义（PlantUML, Mermaid, BlockDiag 等）转换为图像。
  - **适用场景**: 教学 UML、可视化软件架构。

- **[Mindmap MCP Server](https://github.com/YuChenSSR/mindmap-mcp-server)**
  - **描述**: 生成和管理思维导图。有助于整理思路和头脑风暴。
  - **适用场景**: 头脑风暴、整理课堂笔记。

- **[FFmpeg MCP](https://github.com/egoist/ffmpeg-mcp)**
  - **描述**: FFmpeg 的 MCP 服务器。通过聊天即可处理视频和音频。
  - **适用场景**: 编辑讲座录音/录像、转换视频格式。

## <a id='academic-writing'></a>Academic & Writing (学术与写作)

> 学术研究、论文写作、文献综述及创意写作。

- **[Anthropic DOCX Skill](https://github.com/anthropics/anthropic-skills/tree/main/skills/docx)**
  - **描述**: 全面的文档创建、编辑和分析，支持修订模式和格式保留。
  - **适用场景**: 创建新文档、修改内容、处理修订记录。

- **[Anthropic PDF Skill](https://github.com/anthropics/anthropic-skills/tree/main/skills/pdf)**
  - **描述**: 全面的 PDF 操作工具包，用于提取文本/表格、创建 PDF 和处理表单。
  - **适用场景**: 从论文提取数据、填写表单、合并文档。

- **[writing-agent](https://github.com/dongbeixiaohuo/writing-agent)**
  - **描述**: 基于 Claude Code Skills 的"反 AI 味"写作系统。支持从选题生成、风格建模到发布评审的完整工作流。
  - **适用场景**: 让 AI 写出的文章像人一样自然，适合自媒体和博客写作。

- **[The-Crucible-Writing-System-For-Claude](https://github.com/forsonny/The-Crucible-Writing-System-For-Claude)**
  - **描述**: 一个集成的写作系统，包含三个 Claude Skills，引导作家从故事概念到完成初稿。基于 36 拍叙事框架。
  - **适用场景**: 小说创作、剧本写作、创意写作课程。

- **[humanizer](https://github.com/blader/humanizer)**
  - **描述**: 一个 Claude Code skill，专门用于移除文本中的"AI 生成痕迹"，使语言更加自然和人性化。
  - **适用场景**: 润色 AI 起草的草稿，使其更符合人类表达习惯。

- **[ux-writing-skill](https://github.com/content-designer/ux-writing-skill)**
  - **描述**: 系统化的 UX 写作技能，用于强制执行设计系统中的内容质量标准。
  - **适用场景**: 交互设计学生、产品文案规范学习。

- **[vibe-writing](https://github.com/pakholeung37/vibe-writing)**
  - **描述**: 一种与 AI 协作进行小说创作的 wiki 风格方法论。
  - **适用场景**: 探索人机协作创作的新模式。

- **[academic-paper-skills](https://github.com/lishix520/academic-paper-skills)**
  - **描述**: 使用 Claude Code 规划和撰写学术论文的系统框架。包含"策略家"(规划) 和 "作曲家"(写作) 技能。
  - **适用场景**: 硕博研究生论文写作、学术出版规划。

- **[research-units-pipeline-skills](https://github.com/WILLOSCAR/research-units-pipeline-skills)**
  - **描述**: 将研究管道视为语义执行单元。强调证据优先的方法，防止空洞的写作。
  - **适用场景**: 严谨的学术研究项目、文献综述。

- **[academic-writing-skills](https://github.com/bahayonghang/academic-writing-skills)**
  - **描述**: 专门针对学术写作场景优化的 Claude Code 和 Codex 辅助能力。
  - **适用场景**: 提升学术语言的规范性和逻辑性。

- **[GPT Researcher](https://github.com/assafelovic/gpt-researcher)**
  - **描述**: 一个自主代理，可以对任何主题进行深度研究，从网络来源抓取并聚合信息，并附带引用。
  - **适用场景**: 生成文献综述，收集公正信息，避免幻觉。

- **[ZotLink](https://github.com/TonybotNi/ZotLink)**
  - **描述**: 生产就绪的 Zotero MCP 服务器。可保存带有丰富元数据和智能 PDF 附件的开放预印本。
  - **适用场景**: 管理大量参考书目的学者，希望自动化“从论文到资料库”的工作流程。

- **[arXiv Search](https://github.com/andybrandt/mcp-simple-arxiv)**
  - **描述**: 直接从 arXiv 搜索和阅读研究论文。文献综述的必备工具。
  - **适用场景**: 文献综述、查找最新预印本。

- **[Open Library](https://github.com/8enSmith/mcp-open-library)**
  - **描述**: 通过 Internet Archive 的 Open Library API 搜索书籍和作者信息。
  - **适用场景**: 书目数据检索、查找书籍元数据。

- **[Mendeley MCP](https://github.com/pallaprolus/mendeley-mcp)**
  - **描述**: 在 Mendeley 中搜索资料库、浏览文件夹并管理参考文献。
  - **适用场景**: 为论文写作管理参考文献。

- **[TexMCP](https://github.com/devroopsaha744/TexMCP)**
  - **描述**: 将 LaTeX 转换为高质量的 PDF 文档或其他格式。
  - **适用场景**: 编译 LaTeX 论文、生成数学公式。

- **[Doc Co-authoring](https://github.com/openclaw/skills/tree/main/skills/seanphan/doc-coauthoring)**
  - **描述**: 与 AI 共同撰写文档的结构化工作流程，确保一致的语气和内容。
  - **适用场景**: 协作写作、学术论文、报告。

- **[Nano PDF](https://github.com/openclaw/skills/tree/main/skills/steipete/nano-pdf)**
  - **描述**: 使用 nano-pdf CLI 通过自然语言指令编辑 PDF。
  - **适用场景**: 通过聊天修改 PDF 内容、拆分/合并页面。

- **[Mineru PDF Parser](https://github.com/openclaw/skills/tree/main/skills/kesslerio/mineru-pdf-parser-clawdbot-skill)**
  - **描述**: 利用本地 CPU 快速、高质量地将 PDF 解析为 Markdown/JSON。
  - **适用场景**: 从复杂的学术论文中提取干净的文本以供 AI 分析。

- **[Nutrient Doc Processing](https://github.com/openclaw/skills/tree/main/skills/jdrhyne/nutrient-document-processing)**
  - **描述**: 企业级文档处理：格式转换、内容提取和 OCR。
  - **适用场景**: 数字化扫描的教科书、处理繁重的文档工作流。

- **[Markdown Converter](https://github.com/openclaw/skills/tree/main/skills/steipete/markdown-converter)**
  - **描述**: 将各种文档格式（HTML、文本等）转换为干净的 Markdown。
  - **适用场景**: 将学习资料标准化为统一的可读格式。

- **[Typst MCP](https://github.com/johannesbrandenburger/typst-mcp)**
  - **描述**: 编译并验证 Typst 代码（LaTeX 的现代替代品）。
  - **适用场景**: 撰写学术论文、学习排版。

- **[MarkItDown](https://github.com/microsoft/markitdown)**
  - **描述**: 将各种文件格式（PDF、PowerPoint、Word）转换为干净的 Markdown。
  - **适用场景**: 总结讲座、转换学习资料。

- **[Pandoc MCP](https://github.com/vivekVells/mcp-pandoc)**
  - **描述**: 使用 Pandoc 进行无缝文档格式转换（Docx ↔ Markdown ↔ HTML ↔ PDF）。
  - **适用场景**: 转换作业、格式化论文。

- **[Microsoft 365 MCP](https://github.com/softeria/ms-365-mcp-server)**
  - **描述**: 连接 Microsoft Graph API 以使用 Word、Excel 和 Outlook。
  - **适用场景**: 处理作业文档、管理学术邮件。

- **[Wikipedia MCP](https://github.com/Rudra-ravi/wikipedia-mcp)**
  - **描述**: 搜索并检索维基百科的摘要信息。
  - **适用场景**: 通用知识、事实核查、历史研究。

- **[PDFMathTranslate](https://github.com/PDFMathTranslate/PDFMathTranslate)**
  - **描述**: 在保留排版和公式的同时翻译科学论文。
  - **适用场景**: 阅读外文文献、学术翻译。

## <a id='notes-knowledge-base'></a>Notes & Knowledge Base (笔记与知识库)

> 笔记整理、个人知识库管理 (PKM) 及文档处理。

- **[claudian](https://github.com/YishenTu/claudian)**
  - **描述**: Obsidian 插件，将 Claude Code 嵌入到你的笔记库中，作为 AI 协作者。
  - **适用场景**: 学术笔记整理、Zettelkasten (卡片盒笔记法) 维护、论文灵感连接。

- **[Logseq Skill](https://github.com/openclaw/skills/tree/main/skills/juanirm/logseq)**
  - **描述**: 用于与本地 Logseq 图/知识库交互的命令。
  - **适用场景**: 知识管理、每日日记、连接想法。

- **[Rationality](https://github.com/openclaw/skills/tree/main/skills/xertrov/rationality)**
  - **描述**: 理性思维、偏见检查和决策的结构化框架。
  - **适用场景**: 批判性思维、论证分析、复杂决策支持。

- **[claude-code-obsidian-starter](https://github.com/ArtemXTech/claude-code-obsidian-starter)**
  - **描述**: Claude Code + Obsidian 的免费启动套件。预配置了项目管理、任务、每日例程的 Skills。
  - **适用场景**: 需要构建个人知识库 (PKM) 的大学生和研究员。

- **[productivity-skills](https://github.com/mcdow-webworks/productivity-skills)**
  - **描述**: 将 Claude 转化为 AI 生产力伙伴。包含笔记技能，支持对话式捕捉信息和智能搜索。
  - **适用场景**: 快速记录课堂重点、会议纪要整理。

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

- **[Obsidian REST MCP](https://github.com/MarkusPfundstein/mcp-obsidian)**
  - **描述**: 通过 REST API 与本地 Obsidian 库交互，进行知识管理。
  - **适用场景**: 以编程方式管理个人知识库。

- **[Google Drive MCP](https://github.com/isaacphi/mcp-gdrive)**
  - **描述**: 在 Google Drive 中读取和搜索文件。
  - **适用场景**: 访问云文档、整理课程作业。

- **[OpenZIM MCP](https://github.com/cameronrye/openzim-mcp)**
  - **描述**: 在没有网络的情况下离线访问 ZIM 格式的知识库（如维基百科）。
  - **适用场景**: 离线研究、百科知识查询。

- **[Rember (Flashcards)](https://github.com/rember/rember-mcp)**
  - **描述**: 在 Rember 中创建间隔重复抽认卡，以记住你学到的任何东西。
  - **适用场景**: 备考、词汇学习。

- **[Apple Notes MCP](https://github.com/henilcalagiya/mcp-apple-notes)**
  - **描述**: 直接读取和写入 Apple Notes。在 macOS 上整理课堂笔记。
  - **适用场景**: 整理课堂笔记、同步摘要。

## <a id='career-productivity'></a>Career & Productivity (职业规划与生产力)

> 职业发展、项目管理、沟通协作及任务追踪。

- **[Exa AI Search](https://github.com/exa-labs/exa-mcp-server)**
  - **描述**: 专为 AI 代理设计的搜索引擎，用于检索高质量、相关的内容。
  - **适用场景**: 深度研究、查找经过验证的来源。

- **[Brave Search](https://github.com/brave/brave-search-mcp-server)**
  - **描述**: 注重隐私的 Web 搜索，让代理在不被追踪的情况下查找最新信息。
  - **适用场景**: 一般研究、查找时事。

- **[Slack MCP](https://github.com/jtalk22/slack-mcp-server)**
  - **描述**: 与 Slack 频道和私信交互。非常适合学习小组或班级交流。
  - **适用场景**: 团队协作、学习小组管理。

- **[Time Tracking (Toggl)](https://github.com/louis030195/toggl-mcp)**
  - **描述**: 通过 Toggl Track 启动和停止计时器。适用于番茄工作法。
  - **适用场景**: 生产力追踪、学习时长记录。

- **[Filesystem (Go)](https://github.com/mark3labs/mcp-filesystem-server)**
  - **描述**: 强大的本地文件系统访问，允许代理读取你的项目文件（Go 实现）。
  - **适用场景**: 管理本地项目、整理下载内容。

- **[Discord MCP](https://github.com/v-v-vishnu/discord-mcp-server)**
  - **描述**: 将 Claude 连接到 Discord 服务器。阅读消息并管理社区。
  - **适用场景**: 管理学习小组或学生社区。

- **[Linear MCP](https://github.com/jerhadf/linear-mcp-server)**
  - **描述**: 在 Linear 中管理任务和问题。从对话中创建工单。
  - **适用场景**: 学生学习敏捷/Scrum 项目管理。

- **[Todoist MCP](https://github.com/actionsOnGoogle/todoist-mcp-server)**
  - **描述**: 直接管理 Todoist 任务。添加作业或项目截止日期。
  - **适用场景**: 个人任务追踪和截止日期管理。

- **[Resume-Analysis-Assistant](https://github.com/amogh-vardhan/Resume-Analysis-and-Recruitment-Assistant)**
  - **描述**: 使用 NLP 和 Claude 集成来分析简历与职位描述 (JD) 的匹配度。
  - **适用场景**: 毕业生求职、简历优化、模拟面试准备。

- **[JobCoach](https://github.com/sayaleepande/JobCoach)**
  - **描述**: 智能职业指导 Agent，用于简历审查、技能差距分析和学习路径规划。
  - **适用场景**: 职业规划、技能提升路线图。

- **[linkedin-mcp](https://github.com/pegasusheavy/linkedin-mcp)**
  - **描述**: Model Context Protocol (MCP) 服务器，用于管理 LinkedIn 个人资料、技能和教育背景。
  - **适用场景**: 个人品牌建设、自动更新职业履历。

## <a id='meta-skills'></a>Meta Skills (元技能)

> 提示词工程、MCP 元列表及技能开发工具。

- **[claude-prompt-engineering-guide](https://github.com/ThamJiaHe/claude-prompt-engineering-guide)**
  - **描述**: 编写专业 Claude 提示词的综合指南，包含 MCP、Skills 集成。符合 Anthropic 官方最佳实践。
  - **适用场景**: 学习如何更好地控制 AI 模型。

- **[awesome-claude-skills (mejba13)](https://github.com/mejba13/awesome-claude-skills)**
  - **描述**: 收集了大量经过验证的 Claude Skills。
  - **适用场景**: 好的技能集合源。

- **[awesome-claude-skills (ponderous)](https://github.com/ponderous-dustiness314/awesome-claude-skills)**
  - **描述**: 收集了大量经过验证的 Claude Skills。
  - **适用场景**: 好的技能集合源。

- **[everything-claude-code](https://github.com/affaan-m/everything-claude-code)**
  - **描述**: Claude Code 配置的终极集合 - Agents, Skills, Hooks, MCPs。
  - **适用场景**: 想要彻底定制自己 Claude 开发环境的高级计算机学生。

- **[skill-description-optimizer](https://github.com/LZMW/skill-description-optimizer)**
  - **描述**: 用于优化 Claude Skill 描述的"元技能" (Meta-skill)，使用数据驱动分析。
  - **适用场景**: 开发自己的 Claude 技能时使用。

- **[Skill Seekers](https://github.com/yusufkaraaslan/Skill_Seekers)**
  - **描述**: 将文档网站、Git 仓库和 PDF 转换为结构化的 Claude AI Skills，以教授 Claude 新领域的知识。
  - **适用场景**: 通过将文档转换为 AI 导师来掌握新框架。

- **[Anthropic MCP Builder](https://github.com/anthropics/anthropic-skills/tree/main/skills/mcp-builder)**
  - **描述**: 一个用于帮助你构建、测试和调试新 MCP 服务器的元技能。
  - **适用场景**: 为 Claude 开发自定义工具，学习 MCP 架构。

## <a id='commercial-education-ai'></a>Commercial Education AI (商业化教育 AI)

> 正在变革教育行业的领先商业 AI 产品。

- **[Khanmigo](https://khanmigo.ai/)**
  - **描述**: 可汗学院推出的 AI 导学工具，通过提问而非直接给出答案的方式引导学生学习。
  - **适用场景**: 个性化辅导、苏格拉底式教学。

- **[Claude for Education](https://www.anthropic.com)**
  - **描述**: Anthropic 为教育机构设计的负责任的 AI 助手。
  - **适用场景**: 行政任务、写作辅助、编程辅导。

- **[Duolingo Max](https://www.duolingo.com)**
  - **描述**: 多邻国基于 GPT-4 推出的订阅服务，提供"解释我的答案"和"角色扮演"功能。
  - **适用场景**: 语言学习、对话练习。

- **[Speak](https://www.speak.com)**
  - **描述**: 专注于口语对话练习的 AI 语言导师。
  - **适用场景**: 英语口语练习、流利度训练。

- **[Quizlet Q-Chat](https://quizlet.com)**
  - **描述**: Quizlet 内置的 AI 导师，用于交互式学习。
  - **适用场景**: 备考、主动回忆练习。

## <a id='education-llms'></a>Education LLMs (教育大模型)

> 专门针对教育任务、推理和辅导优化的大语言模型。

- **[EduChat](https://github.com/ECNU-ICALK/EduChat)**
  - **描述**: 华东师范大学开发的教育专用大模型，能识别心理和学习需求。
  - **适用场景**: 教育心理学、学生咨询、教学支持。

- **[ChatGLM-Edu (Zhipu)](https://chatglm.cn)**
  - **描述**: 智谱 AI 基于强大的 ChatGLM 基座推出的教育模型。
  - **适用场景**: 通用辅导、语文、数理化问答。

- **[MathGPT (TAL)](https://www.tal.com)**
  - **描述**: 好未来（学而思）推出的专精于数学解题的大模型。
  - **适用场景**: 数学解题、逻辑推理。

- **[DeepSeek-Math](https://github.com/deepseek-ai/DeepSeek-Math)**
  - **描述**: 一个专门的数学模型，推动了开源语言模型在数学推理方面的极限。
  - **适用场景**: 复杂数学问题、竞赛数学、自动推理。

- **[InternLM-Math](https://github.com/InternLM/InternLM-Math)**
  - **描述**: 最先进的双语（中英）开源数学推理大模型。
  - **适用场景**: 双语数学辅导、思维链推理。

- **[MetaMath](https://github.com/meta-math/MetaMath)**
  - **描述**: 一个用于为大语言模型自动生成（Bootstrap）数学问题的项目。
  - **适用场景**: 生成数学训练数据、微调教育模型。

## <a id='agent-frameworks-educational-applications'></a>Agent Frameworks & Educational Applications (Agent 框架与教育应用)

> 用于构建自主智能体的框架，以及基于 LangChain、CrewAI 等的特定教育实现。

- **[LangChain](https://github.com/langchain-ai/langchain)**
  - **描述**: 构建上下文感知推理应用程序和智能体的行业标准框架。
  - **适用场景**: 构建自定义 RAG 导师、教科书分析流程。

- **[CrewAI](https://github.com/crewAIInc/crewAI)**
  - **描述**: 用于编排角色扮演、协同工作的自主 AI 智能体的框架。
  - **适用场景**: 课堂模拟、多角度辩论、教案设计团队。

- **[AI-Tutor](https://github.com/Hawary00/AI-Tutor)**
  - **描述**: 利用 LangChain、LangGraph 和 RAG 进行科学学习的模块化教育助手。
  - **适用场景**: 科学辅导、个性化学习路径。

- **[Multi-Agent Study Assistant](https://github.com/A-R007/Multi-Agent-Study-Assistant)**
  - **描述**: 基于 Phidata 和 LangChain 构建的包含 6 个专业智能体（测验、路线图、问答）的平台。
  - **适用场景**: 综合学习辅助、备考。

- **[Evaluate Crew](https://github.com/poojithM/Multi-Agent-System-for-Educational-Assessment-and-Improvement)**
  - **描述**: 基于 CrewAI 的系统，用于测验生成、表现分析和质量保证。
  - **适用场景**: 自动评分、学生评估。

- **[Multi-Agent RAG (Medical)](https://github.com/SJ9VRF/Multi-Agent-RAG)**
  - **描述**: 基于 AutoGen 的多智能体 RAG 系统，专门用于医学教育文档。
  - **适用场景**: 医学生培训、复杂文档检索。

---

*Last updated/最后更新: 2026-02-12*
