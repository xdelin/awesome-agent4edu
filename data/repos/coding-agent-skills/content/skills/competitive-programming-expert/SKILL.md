---
name: competitive-programming
description: Use this skill when user needs to solve competitive programming problems. Applicable to LeetCode, Codeforces, AtCoder and similar platforms. Triggers include: algorithm problem, coding challenge, time complexity analysis, data structure implementation, TLE (Time Limit Exceeded), MLE (Memory Limit Exceeded).
---

# Competitive Programming Problem Solver

## Description
Solve competitive programming problems with optimal solutions, complexity analysis, and complete code implementation.

## When to Use
- User provides problem link or description from LeetCode/Codeforces/AtCoder/ACM-ICPC platforms
- User requests "solve this algorithm problem" or "optimize this solution"
- User asks for specific algorithm/data structure implementation (e.g., "how to implement segment tree")
- User's code encounters TLE/MLE/WA and needs debugging
- User asks for solution templates or patterns for certain problem types

## When NOT to Use
- User only asks about algorithm concepts or theory without specific problem
- User needs algorithm design for software engineering, not competitive programming
- User is doing system design or architecture problems
- User only needs code completion or syntax help without algorithmic logic

## Input
```typescript
{
  problem: string          // Problem description or link
  platform?: string        // Platform name (leetcode/codeforces/atcoder, etc.)
  language?: string        // Preferred language (default: C++ or Python)
  userCode?: string        // User's existing code (for optimization/debugging)
  constraints?: {          // Problem constraints
    timeLimit?: string     // e.g., "1s", "2s"
    memoryLimit?: string   // e.g., "256MB"
    inputSize?: string     // e.g., "n ≤ 10^5"
  }
}
```

## Output
```typescript
{
  analysis: {
    type: string           // Problem type (DP/Graph/Greedy/Number Theory, etc.)
    keyInsight: string     // Core idea
    edgeCases: string[]    // Edge cases to consider
  }
  solution: {
    approach: string       // Solution explanation
    complexity: {
      time: string         // Time complexity (e.g., O(n log n))
      space: string        // Space complexity
      justification: string // Why it meets problem constraints
    }
    code: string          // Complete executable code
  }
  optimization?: string   // Optional optimization suggestions
}
```

## Execution Steps

### Step 1: Understand Constraints
- Extract input size upper bound (e.g., n ≤ 10^5)
- Calculate time budget (typically 1s ≈ 10^8 operations)
- Identify special restrictions (e.g., read-only once, online algorithm)

### Step 2: Classify and Model
- Categorize problem into known algorithm types (DP/Graph/Greedy/Number Theory/String/Computational Geometry)
- Extract mathematical model or state definition
- List at least 3 typical edge cases

### Step 3: Design Solution
- Explain core idea in 1-2 sentences
- For non-obvious algorithms (e.g., greedy/constructive), briefly justify correctness
- Specify time and space complexity

### Step 4: Implement Code
Output code according to platform conventions:
- **LeetCode**: Provide class/function definition without main function
- **Codeforces/AtCoder**: Provide complete code with standard I/O
- **Other platforms**: Ask user preference

Code requirements:
- Clear variable naming
- Comments on key steps
- Cover identified edge cases

### Step 5: Verify and Optimize
- Validate correctness with sample inputs
- If user provides existing code, compare differences and identify bottlenecks
- If constant-factor optimizations exist (e.g., fast I/O, bitwise tricks), mention separately

## Failure Handling
- **Unclear problem**: Request complete problem statement or link
- **Missing constraints**: Ask for input size bounds and time limits
- **No optimal solution exists**: Provide passable solution first, then discuss if better approach exists
- **Language not supported**: Explain language limitations and suggest alternatives