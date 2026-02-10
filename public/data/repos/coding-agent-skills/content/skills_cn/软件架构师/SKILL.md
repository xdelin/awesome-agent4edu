---
name: software-architect
description: 当用户需要系统架构设计、技术选型或架构评审时使用此技能。触发关键词：系统设计、架构方案、技术选型、微服务、分布式系统、高可用、扩展性、CAP定理、架构权衡、重构方案。适用于需要高层设计决策的场景，不适用于具体代码实现。
---

# 软件架构设计

## 描述
提供系统架构设计方案、技术选型建议和架构权衡分析，确保方案满足业务需求和非功能性要求。

## 何时使用
- 用户请求"设计一个XX系统的架构"或"帮我做技术选型"
- 用户询问系统扩展性、高可用性、性能优化策略
- 用户需要微服务拆分、模块划分、服务边界设计
- 用户提出架构评审或重构需求
- 用户询问 CAP 定理、一致性模型、分布式事务等架构理论在具体场景中的应用
- 用户对比不同架构模式或技术栈的优劣

## 何时不使用
- 用户只需要具体代码实现或编程技巧，不涉及架构层面决策
- 用户的问题是纯算法题或数据结构实现
- 用户需要的是 UI/UX 设计，而非系统架构
- 用户只是询问某个框架或库的使用方法，不涉及系统整体设计
- 用户的项目规模极小（如个人脚本、单文件工具），不需要架构设计

## 输入
```typescript
{
  requirements: {
    businessScenario: string      // 业务场景描述
    userScale?: string            // 用户规模(如 "100万DAU")
    qps?: string                  // 查询/请求量级(如 "1000 QPS峰值")
    dataVolume?: string           // 数据量级(如 "10TB 历史数据")
    latencyRequirement?: string   // 延迟要求(如 "p99 < 100ms")
    readWriteRatio?: string       // 读写比例(如 "读:写 = 9:1")
    consistencyRequirement?: string // 一致性要求(强一致/最终一致/因果一致)
  }
  constraints: {
    budget?: string               // 预算限制
    teamSize?: string             // 团队规模和技能栈
    timeline?: string             // 交付时间窗口
    existingTechStack?: string[]  // 现有技术栈
    complianceRequirements?: string // 合规要求(如 GDPR、等保)
  }
  goals: {
    priority: string              // 核心目标(性能优先/成本优先/快速上线)
    tradeoffs?: string            // 已知的权衡偏好
  }
  existingArchitecture?: string   // 现有架构(用于评审/重构场景)
}
```

## 输出
```typescript
{
  architectureProposal: {
    style: string                 // 架构风格(单体/分层/微服务/事件驱动/CQRS/Serverless)
    topology: string              // 架构拓扑描述(组件、数据流、调用关系)
    diagram?: string              // Mermaid图或C4模型描述
    componentBreakdown: {
      name: string
      responsibility: string
      technology: string
    }[]
  }
  techStackRecommendation: {
    category: string              // 类别(数据库/缓存/消息队列/负载均衡等)
    options: {
      name: string
      pros: string[]
      cons: string[]
      justification: string       // 选择理由
    }[]
    recommendation: string        // 推荐方案
    risks: string[]               // 潜在风险
    mitigationPlan: string        // 风险缓解措施
  }[]
  tradeoffAnalysis: {
    dimension: string             // 权衡维度(性能vs成本/一致性vs可用性等)
    options: {
      choice: string
      pros: string[]
      cons: string[]
      applicableScenarios: string[]
    }[]
    recommendation: string
  }[]
  nonFunctionalRequirements: {
    highAvailability?: string     // 高可用方案(故障转移/降级/熔断/限流)
    scalability?: string          // 扩展性方案(无状态设计/分片/缓存)
    observability?: string        // 可观测性(日志/监控/追踪/告警)
    security?: string             // 安全性(认证/授权/加密/审计)
  }
  evolutionPath?: string          // 架构演进路径(如何从单体迁移到微服务)
}
```

## 执行工作流

在开始执行前，复制以下清单，并在每一步完成后显式标记状态。

### Step 1: 需求澄清与约束识别
- 确认业务场景和核心用例
- 量化关键指标(用户规模、QPS、延迟、数据量)
- 识别读写特征和一致性要求
- 明确资源限制(预算、团队、时间、现有技术栈)
- 确定优先级目标(性能/成本/上线速度)

**反馈闭环**: 若关键约束信息不足(如缺少 QPS 或数据量级)，必须询问用户补充，避免基于假设设计。

### Step 2: 架构风格选择
- 列出候选架构风格(单体/分层/微服务/事件驱动/CQRS/Serverless)
- 对比各风格在当前场景下的适用性
- 说明推荐方案及理由

**反馈闭环**: 若业务场景存在多种合理架构选择，提供 2-3 个方案供用户权衡，而非单一答案。

### Step 3: 组件设计与技术选型
- 定义系统核心组件及其职责边界
- 描述组件间的数据流和调用关系
- 对每个技术选型类别(数据库/缓存/消息队列等)给出候选方案对比
- 说明推荐技术栈及选型理由(性能/成熟度/生态/团队熟悉度)

**反馈闭环**: 若推荐的技术栈与现有技术栈冲突，必须说明迁移成本和演进路径，或提供兼容方案。

### Step 4: 权衡分析
- 识别关键权衡维度(性能vs成本、一致性vs可用性、复杂度vs灵活性)
- 对每个维度给出不同选择的利弊分析
- 结合业务目标给出推荐决策

**示例权衡维度**:
- 性能 vs 成本: 垂直扩展 vs 水平扩展
- 一致性 vs 可用性: CAP 定理的应用(CP vs AP)
- 复杂度 vs 灵活性: 单体 vs 微服务
- 实时性 vs 资源消耗: 推模型 vs 拉模型

### Step 5: 非功能性需求设计
- 高可用方案: 故障转移、服务降级、熔断、限流
- 可扩展方案: 无状态设计、分片策略、缓存层次
- 可观测方案: 日志聚合、指标监控、分布式追踪、告警规则
- 安全性方案: 认证授权、数据加密、审计日志

**反馈闭环**: 若用户未明确非功能性需求优先级，提供默认合理方案，并标注可选增强项。

### Step 6: 风险识别与缓解
- 列出架构方案的潜在风险(技术风险、运维风险、成本风险)
- 为每个风险提供缓解措施
- 说明架构演进路径(如从单体到微服务的迁移策略)

**反馈闭环**: 若方案存在显著风险且无法完全缓解，必须明确告知用户，而非隐瞒或淡化。

## 失败处理

### 需求不明确
- **现象**: 用户只提供模糊描述(如"设计一个高并发系统")
- **处理**: 返回需求澄清问题清单，要求用户补充关键约束(用户规模、QPS、延迟、数据量)

### 约束冲突
- **现象**: 用户要求在极低预算下实现极高性能
- **处理**: 明确指出约束冲突，提供多个方案及其成本/性能权衡，由用户决策

### 技术栈限制
- **现象**: 用户要求使用不适合场景的技术栈(如用关系型数据库存储时序数据)
- **处理**: 说明技术栈的局限性，提供替代方案或折中方案(如 TimescaleDB)

### 架构评审失败
- **现象**: 现有架构存在严重设计缺陷
- **处理**: 明确指出问题点、影响范围和严重性，提供渐进式重构路径，而非要求全部推倒重来