# AI Expert Skills / AI 专家技能库

[English](#english) | [中文](#chinese)

<a name="english"></a>
## 🇬🇧 English

### Introduction
Welcome to the **AI Expert Skills** repository. This collection hosts a set of refined system prompts and skill definitions designed to turn your AI coding assistant into a domain-specific expert. Whether you use **Trae, Cursor, Claude Code, Open Code**, or other AI coding tools, you'll find the right persona here.

These skills are optimized to enhance your pair programming experience by providing structured, high-quality context for your AI.

### Directory Structure
- `skills/`: English skills. Each skill is a directory: `skills/<skill-slug>/SKILL.md`.
- `skills_cn/`: Chinese skills. Each skill is a directory: `skills_cn/<中文技能名>/SKILL.md` (Chinese directory names).

### Skills Included

#### General Coding Skills
| Skill | Description | File |
|-------|-------------|------|
| **Competitive Programming Expert** | Algorithmic problem solving, complexity analysis, and optimization strategies. | [View](skills/competitive-programming-expert/SKILL.md) |
| **C++ Programming Expert** | Modern C++ standards (17/20), memory safety, and performance tuning. | [View](skills/cpp-expert/SKILL.md) |
| **Python Programming Expert** | Pythonic code style, standard library usage, and performance optimization. | [View](skills/python-expert/SKILL.md) |
| **Deep Learning Expert** | Model architecture design, training strategies, and paper interpretation. | [View](skills/deep-learning-expert/SKILL.md) |
| **LLM Testing Expert** | Evaluation strategies, prompt engineering, and safety testing for LLMs. | [View](skills/llm-testing-expert/SKILL.md) |
| **Software Architect** | High-level system design, technology selection, and trade-off analysis. | [View](skills/software-architect/SKILL.md) |
| **GitHub Master** | Git workflows, repository management, and GitHub Actions CI/CD. | [View](skills/github-master/SKILL.md) |

#### Three.js Skills
| Skill | Description | File |
|------|-------------|------|
| **threejs-fundamentals** | Scene setup, cameras, renderer, Object3D hierarchy, coordinate systems. | [View](skills/threejs-fundamentals/SKILL.md) |
| **threejs-geometry** | Built-in shapes, BufferGeometry, custom geometry, instancing. | [View](skills/threejs-geometry/SKILL.md) |
| **threejs-materials** | PBR materials, basic/phong/standard materials, shader materials. | [View](skills/threejs-materials/SKILL.md) |
| **threejs-lighting** | Light types, shadows, environment lighting, light helpers. | [View](skills/threejs-lighting/SKILL.md) |
| **threejs-textures** | Texture types, UV mapping, environment maps, render targets. | [View](skills/threejs-textures/SKILL.md) |
| **threejs-animation** | Keyframe animation, skeletal animation, morph targets, animation mixing. | [View](skills/threejs-animation/SKILL.md) |
| **threejs-loaders** | GLTF/GLB loading, texture loading, async patterns, caching. | [View](skills/threejs-loaders/SKILL.md) |
| **threejs-shaders** | GLSL basics, ShaderMaterial, uniforms, custom effects. | [View](skills/threejs-shaders/SKILL.md) |
| **threejs-postprocessing** | EffectComposer, bloom, DOF, screen effects, custom passes. | [View](skills/threejs-postprocessing/SKILL.md) |
| **threejs-interaction** | Raycasting, camera controls, mouse/touch input, object selection. | [View](skills/threejs-interaction/SKILL.md) |

### How It Works
AI coding tools can load skill files from their skills directory to provide domain-specific context. Example triggers:
- Create a 3D scene → `threejs-fundamentals`
- Add lighting and shadows → `threejs-lighting`
- Load a GLTF model → `threejs-loaders` + `threejs-animation`
- Create custom visual effects → `threejs-shaders` + `threejs-postprocessing`

### How to Use

#### 1. Choose your language
Navigate to `skills/` (English) or `skills_cn/` (Chinese).

#### 2. Select a skill
Choose the skill file for the expert you need (e.g., `software-architect/SKILL.md`).

#### 3. Install the skill

**Project-level** (recommended for team collaboration):
- **Trae**: Place the file in `.trae/skills/` directory.
- **Cursor**: Place the file in `.cursor/skills/` directory (or use `.cursorrules` file in project root for legacy support).
- **Claude Code**: Place the file in `.claude/skills/` directory.
- **Open Code**: Place the file in `.opencode/skills/` directory (also supports `.claude/skills/` for cross-IDE compatibility).

**Global** (applies to all projects):
- **Trae**: `~/.trae/skills/`
- **Cursor**: `~/.cursor/skills/`
- **Claude Code**: `~/.claude/skills/`
- **Open Code**: `~/.config/opencode/skills/`

Tip: Open Code can read skills from `.claude/skills/`, so placing skills there enables sharing between Claude Code and Open Code.

### Skill File Structure

Each skill file uses YAML frontmatter + Markdown content:

```markdown
---
name: skill-name
description: When this skill should be activated
---

# Skill Title

## Quick Start

[Minimal working example]

## Core Concepts

[Key APIs and concepts]

## Common Patterns

[Real-world usage patterns]

## Performance Tips

[Optimization guidance]

## See Also

[Related skills]
```

---

<a name="chinese"></a>
## 🇨🇳 中文

### 简介
欢迎来到 **AI 专家技能库**。本项目收集了一系列精心打磨的系统提示词和技能定义，旨在将您的 AI 编程助手转化为特定领域的专家。无论您使用的是 **Trae、Cursor、Claude Code、Open Code** 还是其他 AI 编程工具，这里都能找到合适的角色。

这些技能旨在通过提供结构化、高质量的上下文，提升您的结对编程体验。

### 目录结构
- `skills/`：英文技能。每个技能是一个目录：`skills/<skill-slug>/SKILL.md`。
- `skills_cn/`：中文技能。每个技能是一个目录：`skills_cn/<中文技能名>/SKILL.md`（目录名使用中文）。

### 技能列表

#### 通用编程技能
| 技能 | 描述 | 文件 |
|------|------|------|
| **竞赛编程专家** | 算法解题、复杂度分析与代码优化策略。 | [查看](skills_cn/竞赛编程专家/SKILL.md) |
| **C++ 编程专家** | 现代 C++ 标准 (17/20)、内存安全与性能调优。 | [查看](skills_cn/C++编程专家/SKILL.md) |
| **Python 编程专家** | Pythonic 代码风格、标准库使用与性能优化。 | [查看](skills_cn/Python编程专家/SKILL.md) |
| **深度学习专家** | 模型架构设计、训练策略与论文解读。 | [查看](skills_cn/深度学习专家/SKILL.md) |
| **大语言模型测试专家** | LLM 评测策略、提示词工程与安全性测试。 | [查看](skills_cn/大语言模型测试专家/SKILL.md) |
| **软件架构师** | 高层系统设计、技术选型与架构权衡分析。 | [查看](skills_cn/软件架构师/SKILL.md) |
| **GitHub 大师** | Git 工作流、仓库管理与 GitHub Actions CI/CD。 | [查看](skills_cn/GitHub大师/SKILL.md) |

#### Three.js 技能
| 技能 | 描述 | 文件 |
|------|------|------|
| **Three.js基础** | 场景/相机/渲染器、坐标系与 Object3D 层级。 | [查看](skills_cn/Three.js基础/SKILL.md) |
| **Three.js几何体** | BufferGeometry、自定义几何与 Instancing。 | [查看](skills_cn/Three.js几何体/SKILL.md) |
| **Three.js材质** | 常用材质选型、PBR 参数、透明与渲染状态。 | [查看](skills_cn/Three.js材质/SKILL.md) |
| **Three.js灯光** | 灯光类型、阴影与调试辅助。 | [查看](skills_cn/Three.js灯光/SKILL.md) |
| **Three.js纹理** | 贴图加载、UV、色彩空间、环境贴图与渲染目标。 | [查看](skills_cn/Three.js纹理/SKILL.md) |
| **Three.js动画** | AnimationMixer、动作播放与混合。 | [查看](skills_cn/Three.js动画/SKILL.md) |
| **Three.js加载器** | GLTF/GLB 加载、资源处理与缓存。 | [查看](skills_cn/Three.js加载器/SKILL.md) |
| **Three.js着色器** | ShaderMaterial、uniforms 与自定义效果。 | [查看](skills_cn/Three.js着色器/SKILL.md) |
| **Three.js后处理** | EffectComposer 与常见屏幕特效链路。 | [查看](skills_cn/Three.js后处理/SKILL.md) |
| **Three.js交互** | Raycaster 拾取、输入事件与选中反馈。 | [查看](skills_cn/Three.js交互/SKILL.md) |

### 工作原理
AI 编程工具会从技能目录加载匹配上下文的技能文件，提供更准确的领域知识。常见触发示例：
- 搭建 3D 场景 → `threejs-fundamentals` / Three.js基础
- 添加灯光与阴影 → `threejs-lighting` / Three.js灯光
- 加载 GLTF 模型 → `threejs-loaders` + `threejs-animation` / Three.js加载器 + Three.js动画
- 自定义视觉效果 → `threejs-shaders` + `threejs-postprocessing` / Three.js着色器 + Three.js后处理

### 使用方法

#### 1. 选择语言
进入 `skills/` (英文) 或 `skills_cn/` (中文) 目录。

#### 2. 选择技能
选择您需要的专家技能文件（如 `软件架构师/SKILL.md`）。

#### 3. 安装技能

**项目级**（推荐用于团队协作）：
- **Trae**: 将文件放入 `.trae/skills/` 目录下。
- **Cursor**: 将文件放入 `.cursor/skills/` 目录下（传统方式支持在项目根目录放置 `.cursorrules` 文件）。
- **Claude Code**: 将文件放入 `.claude/skills/` 目录下。
- **Open Code**: 将文件放入 `.opencode/skills/` 目录下（也支持读取 `.claude/skills/` 实现跨 IDE 兼容）。

**全局**（适用于所有项目）：
- **Trae**: `~/.trae/skills/`
- **Cursor**: `~/.cursor/skills/`
- **Claude Code**: `~/.claude/skills/`
- **Open Code**: `~/.config/opencode/skills/`

提示：Open Code 可以读取 `.claude/skills/` 目录的技能文件，因此将文件放在该目录可实现 Claude Code 与 Open Code 共享。

### 技能文件结构

技能文件由 YAML frontmatter 与 Markdown 正文组成：

```markdown
---
name: skill-name
description: 何时触发该技能
---

# 技能标题

## 快速开始

[最小可运行示例]

## 核心概念

[关键 API 与概念]

## 常用模式

[工程化用法]

## 性能提示

[优化建议]

## 相关技能

[交叉引用]
```

---

### ⚠️ 重要说明 / Important Notes

1. **文件格式**: 技能文件通常包含 YAML frontmatter（元数据）和 Markdown 内容。确保文件名使用 `.md` 扩展名。

2. **Cursor 兼容性**: 
   - 新版 Cursor (0.46+) 支持 Skills 系统，使用 `.cursor/skills/` 目录
   - 旧版 Cursor 使用项目根目录的 `.cursorrules` 文件
   - 建议优先使用 `.cursor/skills/` 以支持多技能管理

3. **多 IDE 共享**: 如果团队使用多种 IDE，建议将技能文件放在 `.claude/skills/`，因为：
   - Claude Code 原生支持
   - Open Code 兼容读取
   - 其他工具（如 Trae、Cursor）可通过软链接或复制使用
