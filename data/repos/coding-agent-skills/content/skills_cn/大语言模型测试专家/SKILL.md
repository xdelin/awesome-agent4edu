---
name: llm-testing-expert
description: 当用户需要设计 LLM 评测方案、构建测试集、进行模型质量保障时使用此技能。触发关键词：LLM 测试、模型评测、提示词工程、回归测试、红队测试、幻觉检测、RAG 测试、Agent 测试、A/B 测试、LLM-as-a-Judge。适用于模型开发和应用落地的质量保障场景。
---

# 大语言模型测试与评测

## 描述
提供 LLM 测试方案设计、测试集构建、自动化评估和安全红队测试建议，确保模型在功能、性能、鲁棒性和安全性上满足生产要求。

## 何时使用
- 用户请求"设计 LLM 测试方案"或"如何评测模型质量"
- 用户需要构建测试集、设计测试用例或进行回归测试
- 用户寻求提示词优化、工程化和版本管理建议
- 用户需要评估模型性能指标（准确率、幻觉率、安全性、延迟、成本）
- 用户询问如何进行 A/B 测试或对比评测（多模型/多版本）
- 用户需要自动化测试流程或 CI/CD 集成方案
- 用户请求红队测试设计（越狱、提示注入、对抗攻击）
- 用户询问 RAG 系统或 Agent 应用的专项测试策略

## 何时不使用
- 用户只需要传统软件测试（单元测试、集成测试），不涉及 LLM 特性
- 用户的问题是模型训练或微调技术，而非测试评估
- 用户需要的是数据标注或数据清洗，而非测试方案设计
- 用户只是询问某个 LLM 框架或工具的使用方法，不涉及测试策略
- 用户的问题是纯粹的 Prompt Engineering（如何写更好的 Prompt），而非如何测试 Prompt 效果

## 输入
```typescript
{
  testTarget: {
    type: string                  // 测试对象类型(base-model/fine-tuned-model/rag-system/agent/prompt-template)
    modelInfo?: {
      name: string                // 模型名称(如 gpt-4, claude-3.5-sonnet)
      version?: string            // 版本号
      deployment?: string         // 部署方式(API/self-hosted)
    }
    applicationContext?: string   // 应用场景(如客服问答、代码生成、文档摘要)
  }
  testObjectives: {
    functional?: boolean          // 功能正确性测试
    performance?: boolean         // 性能测试(延迟、吞吐、成本)
    robustness?: boolean          // 鲁棒性测试(对抗样本、边界用例)
    safety?: boolean              // 安全性测试(越狱、注入、PII泄露)
    userExperience?: boolean      // 用户体验测试
  }
  constraints: {
    budget?: string               // 测试预算(API调用成本/人工标注成本)
    timeline?: string             // 测试周期
    existingTestAssets?: string[] // 现有测试资产(测试集、标注数据)
    complianceRequirements?: string // 合规要求(GDPR、行业监管)
  }
  riskAreas?: string[]            // 已知风险点(如幻觉、偏见、隐私泄露)
  existingMetrics?: string        // 现有评估指标和基线
}
```

## 输出
```typescript
{
  testStrategy: {
    scope: string                 // 测试范围定义
    approach: string              // 测试方法(黑盒/白盒/灰盒)
    testLevels: {
      unit?: string               // 单元测试(单个Prompt/单个工具调用)
      integration?: string        // 集成测试(多轮对话/工具链)
      system?: string             // 系统测试(端到端场景)
      acceptance?: string         // 验收测试(用户场景覆盖)
    }
  }
  testPlan: {
    functional: {
      testCases: {
        id: string
        scenario: string          // 测试场景
        input: string             // 输入样例
        expectedOutput: string    // 期望输出(可用模糊规则)
        passCriteria: string      // 通过标准
      }[]
      coverage: string[]          // 覆盖维度(指令遵循/格式输出/推理能力等)
    }
    performance: {
      metrics: {
        name: string              // 指标名称(latency/throughput/token-cost)
        target: string            // 目标值(如 p95 < 2s)
        measurement: string       // 测量方法
      }[]
      loadProfile: string         // 负载模型(并发用户数、请求模式)
    }
    robustness: {
      adversarialCases: string[]  // 对抗样本设计
      edgeCases: string[]         // 边界用例
      stressScenarios: string[]   // 压力场景(超长输入、极端参数)
    }
    safety: {
      redTeamScenarios: {
        type: string              // 攻击类型(jailbreak/injection/data-extraction)
        technique: string         // 攻击技术
        expectedDefense: string   // 期望防御措施
      }[]
      harmfulContentCategories: string[] // 有害内容类别(暴力/歧视/隐私)
    }
  }
  testDataset: {
    sources: string[]             // 数据来源(公开基准/领域数据/合成数据)
    composition: {
      positive: number            // 正向用例占比
      negative: number            // 负向/对抗用例占比
      edge: number                // 边界用例占比
    }
    sampleSize: string            // 样本量(按统计显著性计算)
    labelingStrategy: string      // 标注策略(人工/自动/混合)
  }
  evaluationMethod: {
    automated: {
      metrics: string[]           // 自动化指标(BLEU/ROUGE/exact-match/regex)
      tools: string[]             // 评估工具
    }
    humanEval: {
      criteria: string[]          // 人工评估标准
      raterGuidelines: string     // 标注员指南
      interRaterAgreement: string // 一致性要求(如 Kappa > 0.7)
    }
    llmAsJudge?: {
      judgeModel: string          // 评判模型
      rubric: string              // 评分规则
      calibration: string         // 校准方法(与人工标注对齐)
    }
  }
  regressionPlan: {
    triggerConditions: string[]   // 触发回归的条件(模型更新/Prompt变更)
    baselineVersion: string       // 基线版本
    comparisonMetrics: string[]   // 对比指标
    reportFormat: string          // 报告格式
  }
  cicdIntegration?: string        // CI/CD 集成方案
  specializedTests?: {
    rag?: {
      retrievalQuality: string    // 检索质量测试(召回率/排序)
      citationAccuracy: string    // 引用准确性测试
      faithfulness: string        // 忠实度测试(是否仅基于检索内容)
    }
    agent?: {
      toolCallCorrectness: string // 工具调用参数正确性
      planningRationality: string // 规划合理性
      errorRecovery: string       // 错误恢复能力
    }
  }
}
```

## 执行工作流

在开始执行前，复制以下清单，并在每一步完成后显式标记状态。

### Step 1: 测试目标与风险识别
- 明确测试对象类型(基础模型/微调模型/RAG/Agent/Prompt 模板)
- 确认核心测试目标(功能/性能/鲁棒性/安全性/用户体验)
- 识别已知风险区域(幻觉、偏见、隐私泄露、提示注入)
- 了解应用场景和约束(预算、时间、合规要求)

**反馈闭环**: 若测试目标不明确或冲突(如既要全面覆盖又要极低成本)，必须与用户对齐优先级。

### Step 2: 测试策略制定
- 选择测试方法(黑盒/白盒/灰盒)
- 定义测试层次(单元/集成/系统/验收)，采用测试金字塔模型
- 确定测试覆盖范围(功能维度、性能指标、安全场景)
- 选择评估方式(自动化/人工/LLM-as-a-Judge/混合)

**测试金字塔原则**:
- 底层(单元测试): 数量最多，成本最低，执行最快(如单 Prompt 功能测试)
- 中层(集成测试): 适中数量，测试组件交互(如多轮对话、工具链)
- 顶层(系统测试): 少量关键场景，端到端验证(如完整用户旅程)

**反馈闭环**: 若用户预算或时间极其有限，优先设计高优先级的冒烟测试集，而非追求全面覆盖。

### Step 3: 测试集构建
- 确定数据来源(公开基准如 MMLU/HumanEval、领域数据、脱敏生产日志、合成数据)
- 设计用例分布: 正向用例(60-70%) + 负向/对抗用例(20-30%) + 边界用例(10%)
- 为每个用例定义 ID、场景、输入、期望输出、通过标准
- 确定样本量(根据统计显著性要求，通常需要 100+ 样本)

**用例设计原则**:
- **正向用例**: 覆盖核心功能和常见场景
- **负向用例**: 测试拒答能力、错误处理、对抗样本
- **边界用例**: 超长输入、极端参数、多语言混合、特殊字符

**反馈闭环**: 若现有测试资产不足，优先从生产日志中采样真实用例，而非完全合成数据。

### Step 4: 评估方法设计
- 自动化评估: 选择合适指标(BLEU/ROUGE/exact-match/regex/结构化输出校验)
- 人工评估: 制定评估标准、标注员指南、确保一致性(Kappa > 0.7)
- LLM-as-a-Judge: 选择评判模型、设计评分规则、与人工标注校准

**评估方式选择**:
- **简单任务**(如分类、提取): 自动化评估(exact-match/F1)
- **复杂任务**(如摘要、创作): LLM-as-a-Judge + 人工抽样验证
- **安全性测试**: 人工评估(检测有害内容、越狱成功率)

**反馈闭环**: 若自动化指标与人工评估不一致，必须重新校准或调整指标权重。

### Step 5: 专项测试设计(如适用)

#### RAG 系统测试
- **检索质量**: 测试召回率(Recall@K)、排序质量(MRR/NDCG)
- **引用准确性**: 验证引用来源是否正确、是否存在幻觉引用
- **忠实度**: 检查回答是否仅基于检索内容，无额外捏造

#### Agent 系统测试
- **工具调用正确性**: 验证工具选择、参数传递是否正确
- **规划合理性**: 评估步骤规划是否高效、逻辑是否连贯
- **错误恢复**: 测试遇到工具失败或异常时的处理能力

**反馈闭环**: 若 RAG/Agent 测试暴露系统性问题(如检索质量差、工具调用失败率高)，应返回系统设计层面优化，而非仅调整测试。

### Step 6: 版本化与回归测试
- 记录模型版本、Prompt 版本、测试集版本
- 定义回归触发条件(模型更新/Prompt 变更/测试集扩展)
- 建立基线(Baseline)并对比新版本性能
- 生成可追溯的测试报告(包含通过率、性能对比、失败用例分析)

**回归测试原则**:
- 每次变更必须通过核心测试集(Golden Set)
- 新增功能必须补充对应测试用例
- 性能退化超过阈值(如准确率下降 >5%)应阻断发布

**反馈闭环**: 若回归测试发现性能退化，必须分析根因(模型问题/Prompt 问题/测试集问题)，而非直接回滚。

### Step 7: 安全红队测试(如适用)
- 设计越狱场景(Jailbreak): 测试绕过安全护栏的攻击
- 设计注入攻击(Prompt Injection): 测试恶意指令注入
- 设计数据提取攻击: 测试 PII 泄露、训练数据记忆
- 设计有害内容生成: 测试暴力、歧视、误导信息等

**红队测试方法**:
- **手工测试**: 由安全专家设计对抗样本
- **自动化攻击**: 使用工具(如 Garak、PromptInject)生成攻击样本
- **众包红队**: 邀请外部人员尝试攻击

**反馈闭环**: 若红队测试发现严重漏洞，必须优先修复并重新测试，而非掩盖或忽略。

## 失败处理

### 测试覆盖不足
- **现象**: 用户提供的测试目标过于宽泛，无法在有限资源下全面覆盖
- **处理**: 基于风险优先级(高风险场景优先)设计分层测试计划，明确哪些场景暂不覆盖

### 评估指标选择错误
- **现象**: 用户使用不适合的指标(如用 BLEU 评估对话质量)
- **处理**: 解释指标局限性，推荐更合适的指标组合(如对话质量用 LLM-as-a-Judge + 人工评估)

### 自动化评估与人工评估不一致
- **现象**: 自动化指标显示性能提升，但人工评估发现质量下降
- **处理**: 分析不一致原因(指标设计问题/标注偏差/模型过拟合)，重新校准评估体系

### 测试数据质量问题
- **现象**: 测试集包含大量噪声、重复或不现实的样本
- **处理**: 建立测试集质量检查流程(去重/异常检测/代表性验证)，优先使用脱敏的真实生产数据

### 红队测试发现严重漏洞
- **现象**: 模型可被轻易越狱或注入
- **处理**: 立即暂停部署，修复安全问题(加强输入过滤/输出审查/模型微调)，重新进行红队测试
