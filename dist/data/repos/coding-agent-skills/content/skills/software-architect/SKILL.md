---
name: software-architect
description: Use this skill when user needs system architecture design, technology stack selection, or architecture review. Trigger keywords: system design, architecture proposal, tech stack, microservices, distributed systems, high availability, scalability, CAP theorem, architecture tradeoffs, refactoring plan. Applicable to high-level design decisions, not specific code implementation.
---

# Software Architecture Design

## Description
Provide system architecture design proposals, technology stack recommendations, and architecture tradeoff analysis to ensure solutions meet business requirements and non-functional requirements.

## When to Use
- User requests "design architecture for XX system" or "help me choose tech stack"
- User asks about system scalability, high availability, or performance optimization strategies
- User needs microservices decomposition, module partitioning, or service boundary design
- User requests architecture review or refactoring proposal
- User asks how to apply architectural theories (CAP theorem, consistency models, distributed transactions) to specific scenarios
- User compares different architecture patterns or technology stacks

## When NOT to Use
- User only needs specific code implementation or programming tips without architectural decisions
- User's problem is a pure algorithm or data structure implementation
- User needs UI/UX design rather than system architecture
- User only asks how to use a framework or library without system-level design
- User's project is extremely small (e.g., personal script, single-file tool) and doesn't require architecture

## Input
```typescript
{
  requirements: {
    businessScenario: string      // Business scenario description
    userScale?: string            // User scale (e.g., "1M DAU")
    qps?: string                  // Query/Request volume (e.g., "1000 QPS peak")
    dataVolume?: string           // Data volume (e.g., "10TB historical data")
    latencyRequirement?: string   // Latency requirement (e.g., "p99 < 100ms")
    readWriteRatio?: string       // Read/Write ratio (e.g., "read:write = 9:1")
    consistencyRequirement?: string // Consistency (strong/eventual/causal)
  }
  constraints: {
    budget?: string               // Budget constraints
    teamSize?: string             // Team size and skill set
    timeline?: string             // Delivery timeline
    existingTechStack?: string[]  // Existing technology stack
    complianceRequirements?: string // Compliance (e.g., GDPR, SOC2)
  }
  goals: {
    priority: string              // Core objective (performance/cost/time-to-market)
    tradeoffs?: string            // Known tradeoff preferences
  }
  existingArchitecture?: string   // Existing architecture (for review/refactoring)
}
```

## Output
```typescript
{
  architectureProposal: {
    style: string                 // Architecture style (monolithic/layered/microservices/event-driven/CQRS/serverless)
    topology: string              // Topology description (components, data flow, call relationships)
    diagram?: string              // Mermaid diagram or C4 model description
    componentBreakdown: {
      name: string
      responsibility: string
      technology: string
    }[]
  }
  techStackRecommendation: {
    category: string              // Category (database/cache/message queue/load balancer, etc.)
    options: {
      name: string
      pros: string[]
      cons: string[]
      justification: string       // Selection rationale
    }[]
    recommendation: string        // Recommended solution
    risks: string[]               // Potential risks
    mitigationPlan: string        // Risk mitigation measures
  }[]
  tradeoffAnalysis: {
    dimension: string             // Tradeoff dimension (performance vs cost/consistency vs availability, etc.)
    options: {
      choice: string
      pros: string[]
      cons: string[]
      applicableScenarios: string[]
    }[]
    recommendation: string
  }[]
  nonFunctionalRequirements: {
    highAvailability?: string     // HA solution (failover/degradation/circuit breaking/rate limiting)
    scalability?: string          // Scalability solution (stateless design/sharding/caching)
    observability?: string        // Observability (logging/monitoring/tracing/alerting)
    security?: string             // Security (authentication/authorization/encryption/audit)
  }
  evolutionPath?: string          // Architecture evolution path (e.g., monolith to microservices migration)
}
```

## Execution Workflow

Copy the following checklist before starting, and explicitly mark status after completing each step.

### Step 1: Requirement Clarification and Constraint Identification
- Confirm business scenario and core use cases
- Quantify key metrics (user scale, QPS, latency, data volume)
- Identify read/write patterns and consistency requirements
- Clarify resource constraints (budget, team, timeline, existing tech stack)
- Determine priority objectives (performance/cost/time-to-market)

**Feedback Loop**: If critical constraint information is missing (e.g., QPS or data volume), MUST ask user to provide details. Avoid designing based on assumptions.

### Step 2: Architecture Style Selection
- List candidate architecture styles (monolithic/layered/microservices/event-driven/CQRS/serverless)
- Compare applicability of each style in current scenario
- Explain recommended approach and rationale

**Feedback Loop**: If business scenario allows multiple reasonable architectures, provide 2-3 options for user to weigh, rather than a single answer.

### Step 3: Component Design and Technology Selection
- Define core system components and their responsibility boundaries
- Describe data flow and call relationships between components
- For each technology category (database/cache/message queue, etc.), provide candidate comparisons
- Explain recommended tech stack and selection rationale (performance/maturity/ecosystem/team familiarity)

**Feedback Loop**: If recommended tech stack conflicts with existing stack, MUST explain migration cost and evolution path, or provide compatible alternatives.

### Step 4: Tradeoff Analysis
- Identify key tradeoff dimensions (performance vs cost, consistency vs availability, complexity vs flexibility)
- Provide pros/cons analysis for different choices in each dimension
- Give recommended decision based on business objectives

**Example Tradeoff Dimensions**:
- Performance vs Cost: Vertical scaling vs Horizontal scaling
- Consistency vs Availability: CAP theorem application (CP vs AP)
- Complexity vs Flexibility: Monolith vs Microservices
- Real-time vs Resource consumption: Push vs Pull model

### Step 5: Non-Functional Requirements Design
- High Availability: Failover, service degradation, circuit breaking, rate limiting
- Scalability: Stateless design, sharding strategy, caching layers
- Observability: Log aggregation, metrics monitoring, distributed tracing, alerting rules
- Security: Authentication/authorization, data encryption, audit logs

**Feedback Loop**: If user hasn't clarified non-functional requirement priorities, provide default reasonable solutions and mark optional enhancements.

### Step 6: Risk Identification and Mitigation
- List potential risks of architecture proposal (technical, operational, cost risks)
- Provide mitigation measures for each risk
- Explain architecture evolution path (e.g., monolith to microservices migration strategy)

**Feedback Loop**: If solution has significant unmitigable risks, MUST explicitly inform user rather than hide or downplay.

## Failure Handling

### Unclear Requirements
- **Symptom**: User only provides vague description (e.g., "design a high-concurrency system")
- **Action**: Return requirement clarification checklist, ask user to provide key constraints (user scale, QPS, latency, data volume)

### Conflicting Constraints
- **Symptom**: User demands extremely high performance with extremely low budget
- **Action**: Clearly point out constraint conflicts, provide multiple solutions with cost/performance tradeoffs, let user decide

### Tech Stack Limitations
- **Symptom**: User requires using unsuitable tech stack (e.g., relational database for time-series data)
- **Action**: Explain tech stack limitations, provide alternatives or compromises (e.g., TimescaleDB)

### Architecture Review Failure
- **Symptom**: Existing architecture has serious design flaws
- **Action**: Clearly point out issues, impact scope, and severity. Provide incremental refactoring path rather than complete rewrite.