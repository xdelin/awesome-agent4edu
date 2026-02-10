---
name: competitive-programming
description: 当用户需要解决算法竞赛题目时使用此技能。适用于 LeetCode、Codeforces、AtCoder 等平台的题目求解、算法优化、复杂度分析。触发关键词：算法题、刷题、竞赛、时间复杂度、数据结构实现、TLE(超时)、MLE(内存超限)。
---

# 竞赛编程题解

## 描述
解决算法竞赛题目，提供最优解法、复杂度分析和完整实现代码。

## 何时使用
- 用户提供 LeetCode/Codeforces/AtCoder/ACM-ICPC 等平台的题目链接或题目描述
- 用户请求"解决这道算法题"、"帮我优化这个解法"
- 用户询问特定算法或数据结构的实现（如"如何实现线段树"）
- 用户的代码遇到 TLE(超时)、MLE(内存超限)、WA(答案错误) 需要调试
- 用户询问某类题目的解题模板或套路

## 何时不使用
- 用户只是询问算法概念或理论知识，不涉及具体题目求解
- 用户需要的是软件工程中的算法设计，而非竞赛题目
- 用户在做系统设计、架构设计等工程问题
- 用户只需要代码补全或语法帮助，不涉及算法逻辑

## 输入
```typescript
{
  problem: string          // 题目描述或链接
  platform?: string        // 平台名称(leetcode/codeforces/atcoder等)
  language?: string        // 首选编程语言(默认C++或Python)
  userCode?: string        // 用户已有代码(用于优化/调试)
  constraints?: {          // 题目约束
    timeLimit?: string     // 如 "1s", "2s"
    memoryLimit?: string   // 如 "256MB"
    inputSize?: string     // 如 "n ≤ 10^5"
  }
}
```

## 输出
```typescript
{
  analysis: {
    type: string           // 题目类型(DP/图论/贪心/数论等)
    keyInsight: string     // 核心思路
    edgeCases: string[]    // 需考虑的边界情况
  }
  solution: {
    approach: string       // 解法说明
    complexity: {
      time: string         // 时间复杂度(如 O(n log n))
      space: string        // 空间复杂度
      justification: string // 为何满足题目限制
    }
    code: string          // 完整可执行代码
  }
  optimization?: string   // 可选的优化建议
}
```

## 执行步骤

### Step 1: 理解约束
- 提取输入规模上界(如 n ≤ 10^5)
- 计算时间预算(通常 1s ≈ 10^8 次运算)
- 识别特殊限制(如只读一次、在线算法等)

### Step 2: 分类与建模
- 将题目归类到已知算法类型(动态规划/图论/贪心/数论/字符串/计算几何等)
- 提取数学模型或状态定义
- 列出至少 3 个典型边界情况

### Step 3: 设计解法
- 说明核心思路(用 1-2 句话)
- 若非显然算法(如贪心/构造),简要说明正确性
- 标注时间复杂度和空间复杂度

### Step 4: 实现代码
按平台规范输出代码:
- **LeetCode**: 给出类/函数定义,不含 main 函数
- **Codeforces/AtCoder**: 给出完整代码,包含标准输入输出
- **其他平台**: 询问用户偏好

代码要求:
- 变量命名清晰
- 关键步骤加注释
- 覆盖已识别的边界情况

### Step 5: 验证与优化
- 用示例输入验证输出正确性
- 若用户提供已有代码,对比差异并指出瓶颈
- 若存在常数优化空间(如快速 I/O、位运算技巧),额外说明

## 失败处理
- **题目描述不清**: 要求用户提供完整题面或链接
- **约束缺失**: 询问输入规模上界和时间限制
- **无最优解**: 先给出可通过的解法,再讨论是否存在更优方案
- **语言不支持**: 说明该语言的限制,建议替代方案
