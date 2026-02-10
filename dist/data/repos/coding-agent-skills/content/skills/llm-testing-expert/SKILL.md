---
name: llm-testing-expert
description: Use this skill when user needs to design LLM evaluation strategies, build test datasets, or ensure model quality. Trigger keywords: LLM testing, model evaluation, prompt engineering, regression testing, red teaming, hallucination detection, RAG testing, Agent testing, A/B testing, LLM-as-a-Judge. Applicable to quality assurance in model development and application deployment.
---

# LLM Testing and Evaluation

## Description
Provide LLM test strategy design, test dataset construction, automated evaluation, and security red teaming recommendations to ensure models meet production requirements for functionality, performance, robustness, and safety.

## When to Use
- User requests "design LLM testing strategy" or "how to evaluate model quality"
- User needs to build test datasets, design test cases, or perform regression testing
- User seeks prompt optimization, engineering, and version management recommendations
- User needs to evaluate model performance metrics (accuracy, hallucination rate, safety, latency, cost)
- User asks how to conduct A/B testing or comparative evaluation (multi-model/multi-version)
- User needs automated testing workflows or CI/CD integration
- User requests red team testing design (jailbreak, prompt injection, adversarial attacks)
- User asks about specialized testing strategies for RAG systems or Agent applications

## When NOT to Use
- User only needs traditional software testing (unit tests, integration tests) without LLM-specific concerns
- User's problem is about model training or fine-tuning techniques, not testing/evaluation
- User needs data annotation or data cleaning, not test strategy design
- User only asks how to use an LLM framework or tool, not testing strategy
- User's problem is pure Prompt Engineering (how to write better prompts), not how to test prompt effectiveness

## Input
```typescript
{
  testTarget: {
    type: string                  // Test object type (base-model/fine-tuned-model/rag-system/agent/prompt-template)
    modelInfo?: {
      name: string                // Model name (e.g., gpt-4, claude-3.5-sonnet)
      version?: string            // Version
      deployment?: string         // Deployment method (API/self-hosted)
    }
    applicationContext?: string   // Application scenario (e.g., customer service, code generation, document summarization)
  }
  testObjectives: {
    functional?: boolean          // Functional correctness testing
    performance?: boolean         // Performance testing (latency, throughput, cost)
    robustness?: boolean          // Robustness testing (adversarial samples, edge cases)
    safety?: boolean              // Safety testing (jailbreak, injection, PII leakage)
    userExperience?: boolean      // User experience testing
  }
  constraints: {
    budget?: string               // Testing budget (API call cost/manual annotation cost)
    timeline?: string             // Testing timeline
    existingTestAssets?: string[] // Existing test assets (test sets, annotated data)
    complianceRequirements?: string // Compliance requirements (GDPR, industry regulations)
  }
  riskAreas?: string[]            // Known risk areas (e.g., hallucination, bias, privacy leakage)
  existingMetrics?: string        // Existing evaluation metrics and baselines
}
```

## Output
```typescript
{
  testStrategy: {
    scope: string                 // Test scope definition
    approach: string              // Test approach (black-box/white-box/gray-box)
    testLevels: {
      unit?: string               // Unit testing (single prompt/single tool call)
      integration?: string        // Integration testing (multi-turn dialogue/tool chain)
      system?: string             // System testing (end-to-end scenarios)
      acceptance?: string         // Acceptance testing (user scenario coverage)
    }
  }
  testPlan: {
    functional: {
      testCases: {
        id: string
        scenario: string          // Test scenario
        input: string             // Input sample
        expectedOutput: string    // Expected output (can use fuzzy rules)
        passCriteria: string      // Pass criteria
      }[]
      coverage: string[]          // Coverage dimensions (instruction following/format output/reasoning ability, etc.)
    }
    performance: {
      metrics: {
        name: string              // Metric name (latency/throughput/token-cost)
        target: string            // Target value (e.g., p95 < 2s)
        measurement: string       // Measurement method
      }[]
      loadProfile: string         // Load model (concurrent users, request pattern)
    }
    robustness: {
      adversarialCases: string[]  // Adversarial sample design
      edgeCases: string[]         // Edge cases
      stressScenarios: string[]   // Stress scenarios (extra-long input, extreme parameters)
    }
    safety: {
      redTeamScenarios: {
        type: string              // Attack type (jailbreak/injection/data-extraction)
        technique: string         // Attack technique
        expectedDefense: string   // Expected defense measures
      }[]
      harmfulContentCategories: string[] // Harmful content categories (violence/discrimination/privacy)
    }
  }
  testDataset: {
    sources: string[]             // Data sources (public benchmarks/domain data/synthetic data)
    composition: {
      positive: number            // Positive case ratio
      negative: number            // Negative/adversarial case ratio
      edge: number                // Edge case ratio
    }
    sampleSize: string            // Sample size (calculated by statistical significance)
    labelingStrategy: string      // Labeling strategy (manual/automated/hybrid)
  }
  evaluationMethod: {
    automated: {
      metrics: string[]           // Automated metrics (BLEU/ROUGE/exact-match/regex)
      tools: string[]             // Evaluation tools
    }
    humanEval: {
      criteria: string[]          // Human evaluation criteria
      raterGuidelines: string     // Annotator guidelines
      interRaterAgreement: string // Consistency requirement (e.g., Kappa > 0.7)
    }
    llmAsJudge?: {
      judgeModel: string          // Judge model
      rubric: string              // Scoring rules
      calibration: string         // Calibration method (align with human annotation)
    }
  }
  regressionPlan: {
    triggerConditions: string[]   // Conditions triggering regression (model update/prompt change)
    baselineVersion: string       // Baseline version
    comparisonMetrics: string[]   // Comparison metrics
    reportFormat: string          // Report format
  }
  cicdIntegration?: string        // CI/CD integration solution
  specializedTests?: {
    rag?: {
      retrievalQuality: string    // Retrieval quality testing (recall/ranking)
      citationAccuracy: string    // Citation accuracy testing
      faithfulness: string        // Faithfulness testing (based only on retrieved content)
    }
    agent?: {
      toolCallCorrectness: string // Tool call parameter correctness
      planningRationality: string // Planning rationality
      errorRecovery: string       // Error recovery capability
    }
  }
}
```

## Execution Workflow

Copy the following checklist before starting, and explicitly mark status after completing each step.

### Step 1: Test Objective and Risk Identification
- Clarify test object type (base model/fine-tuned model/RAG/Agent/prompt template)
- Confirm core test objectives (functional/performance/robustness/safety/user experience)
- Identify known risk areas (hallucination, bias, privacy leakage, prompt injection)
- Understand application scenario and constraints (budget, timeline, compliance requirements)

**Feedback Loop**: If test objectives are unclear or conflicting (e.g., comprehensive coverage with extremely low cost), MUST align priorities with user.

### Step 2: Test Strategy Formulation
- Choose test approach (black-box/white-box/gray-box)
- Define test levels (unit/integration/system/acceptance), adopt test pyramid model
- Determine test coverage (functional dimensions, performance metrics, security scenarios)
- Select evaluation methods (automated/manual/LLM-as-a-Judge/hybrid)

**Test Pyramid Principle**:
- Bottom (unit tests): Most numerous, lowest cost, fastest execution (e.g., single prompt function test)
- Middle (integration tests): Moderate quantity, test component interaction (e.g., multi-turn dialogue, tool chain)
- Top (system tests): Few critical scenarios, end-to-end validation (e.g., complete user journey)

**Feedback Loop**: If user budget or timeline is extremely limited, prioritize designing high-priority smoke tests rather than pursuing comprehensive coverage.

### Step 3: Test Dataset Construction
- Determine data sources (public benchmarks like MMLU/HumanEval, domain data, anonymized production logs, synthetic data)
- Design case distribution: Positive cases (60-70%) + Negative/adversarial cases (20-30%) + Edge cases (10%)
- Define ID, scenario, input, expected output, pass criteria for each case
- Determine sample size (based on statistical significance, typically requires 100+ samples)

**Case Design Principles**:
- **Positive cases**: Cover core functionality and common scenarios
- **Negative cases**: Test rejection capability, error handling, adversarial samples
- **Edge cases**: Extra-long input, extreme parameters, multilingual mix, special characters

**Feedback Loop**: If existing test assets are insufficient, prioritize sampling real cases from production logs rather than fully synthetic data.

### Step 4: Evaluation Method Design
- Automated evaluation: Select appropriate metrics (BLEU/ROUGE/exact-match/regex/structured output validation)
- Manual evaluation: Develop evaluation criteria, annotator guidelines, ensure consistency (Kappa > 0.7)
- LLM-as-a-Judge: Select judge model, design scoring rules, calibrate with human annotation

**Evaluation Method Selection**:
- **Simple tasks** (e.g., classification, extraction): Automated evaluation (exact-match/F1)
- **Complex tasks** (e.g., summarization, creation): LLM-as-a-Judge + human sampling validation
- **Safety testing**: Manual evaluation (detect harmful content, jailbreak success rate)

**Feedback Loop**: If automated metrics are inconsistent with manual evaluation, MUST recalibrate or adjust metric weights.

### Step 5: Specialized Test Design (if applicable)

#### RAG System Testing
- **Retrieval quality**: Test recall (Recall@K), ranking quality (MRR/NDCG)
- **Citation accuracy**: Verify citation sources are correct, check for hallucinated citations
- **Faithfulness**: Check if answer is based only on retrieved content, no fabrication

#### Agent System Testing
- **Tool call correctness**: Verify tool selection and parameter passing are correct
- **Planning rationality**: Assess if step planning is efficient and logic is coherent
- **Error recovery**: Test handling capability when encountering tool failures or exceptions

**Feedback Loop**: If RAG/Agent testing reveals systematic issues (e.g., poor retrieval quality, high tool call failure rate), should return to system design level for optimization, not just adjust testing.

### Step 6: Versioning and Regression Testing
- Record model version, prompt version, test dataset version
- Define regression trigger conditions (model update/prompt change/test set expansion)
- Establish baseline and compare new version performance
- Generate traceable test reports (including pass rate, performance comparison, failure case analysis)

**Regression Testing Principles**:
- Every change must pass core test set (Golden Set)
- New features must add corresponding test cases
- Performance degradation exceeding threshold (e.g., accuracy drop >5%) should block release

**Feedback Loop**: If regression testing finds performance degradation, MUST analyze root cause (model issue/prompt issue/test set issue), not just rollback directly.

### Step 7: Security Red Team Testing (if applicable)
- Design jailbreak scenarios: Test bypassing safety guardrails
- Design injection attacks (Prompt Injection): Test malicious instruction injection
- Design data extraction attacks: Test PII leakage, training data memorization
- Design harmful content generation: Test violence, discrimination, misinformation, etc.

**Red Team Testing Methods**:
- **Manual testing**: Security experts design adversarial samples
- **Automated attacks**: Use tools (e.g., Garak, PromptInject) to generate attack samples
- **Crowdsourced red teaming**: Invite external parties to attempt attacks

**Feedback Loop**: If red team testing finds serious vulnerabilities, MUST prioritize fixing and retesting, not cover up or ignore.

## Failure Handling

### Insufficient Test Coverage
- **Symptom**: User-provided test objectives are too broad, cannot comprehensively cover with limited resources
- **Action**: Design tiered test plan based on risk priority (high-risk scenarios first), clarify which scenarios are not covered

### Wrong Evaluation Metrics
- **Symptom**: User uses inappropriate metrics (e.g., using BLEU to evaluate dialogue quality)
- **Action**: Explain metric limitations, recommend more suitable metric combinations (e.g., dialogue quality uses LLM-as-a-Judge + manual evaluation)

### Automated and Manual Evaluation Inconsistency
- **Symptom**: Automated metrics show performance improvement, but manual evaluation finds quality degradation
- **Action**: Analyze inconsistency reasons (metric design issue/annotation bias/model overfitting), recalibrate evaluation system

### Test Data Quality Issues
- **Symptom**: Test set contains large amounts of noise, duplicates, or unrealistic samples
- **Action**: Establish test set quality check process (deduplication/anomaly detection/representativeness verification), prioritize using anonymized real production data

### Red Team Testing Finds Serious Vulnerabilities
- **Symptom**: Model can be easily jailbroken or injected
- **Action**: Immediately pause deployment, fix security issues (strengthen input filtering/output review/model fine-tuning), re-conduct red team testing
