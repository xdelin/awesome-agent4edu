# Prompt Testing and Evaluation Framework

## Purpose

This document outlines the framework for systematically testing and evaluating different prompt strategies for the Kvante math tutoring system.

## Evaluation Dimensions

### 1. Educational Effectiveness
- **Learning Outcomes**: Does the response help students learn?
- **Pedagogical Soundness**: Are educational best practices followed?
- **Engagement**: Does the response maintain student interest?
- **Scaffolding**: Is appropriate support provided?

### 2. Mathematical Accuracy
- **Correctness**: Are all mathematical statements accurate?
- **Completeness**: Is important information included?
- **Precision**: Is mathematical language used correctly?
- **Consistency**: Are explanations internally consistent?

### 3. User Experience
- **Clarity**: Is the response easy to understand?
- **Helpfulness**: Does it address the student's need?
- **Appropriateness**: Is the difficulty level suitable?
- **Encouragement**: Does it maintain student motivation?

## Testing Methodology

### Test Problem Categories

#### Basic Algebra
- Linear equations
- Quadratic equations
- Systems of equations
- Inequalities

#### Geometry
- Area and perimeter
- Angle relationships  
- Similar triangles
- Coordinate geometry

#### Calculus
- Limits
- Derivatives
- Integrals
- Applications

#### Word Problems
- Rate problems
- Mixture problems
- Geometry applications
- Optimization

### Evaluation Rubric

Rate each dimension on a scale of 1-5:

**5 - Excellent**
- Exceeds expectations in all criteria
- Could serve as a model response
- Highly effective for learning

**4 - Good**
- Meets most criteria well
- Minor areas for improvement
- Generally effective

**3 - Satisfactory**
- Meets basic requirements
- Some weaknesses present
- Acceptable but not optimal

**2 - Needs Improvement**
- Significant gaps in quality
- Multiple issues to address
- Limited effectiveness

**1 - Poor**
- Fails to meet basic requirements
- Major problems present
- Not suitable for students

### Sample Evaluation Questions

#### For Step-by-Step Solutions:
1. Does the response break down the problem appropriately?
2. Are steps logically sequenced?
3. Is the level of detail appropriate for the student?
4. Does it avoid giving away the final answer?
5. Are guiding questions effective?

#### For Hints:
1. Is the hint helpful without being too revealing?
2. Does it address the specific point where the student is stuck?
3. Is it encouraging and supportive?
4. Does it guide toward the correct approach?

#### For Step Checking:
1. Is the evaluation of the student's work accurate?
2. Is feedback constructive and specific?
3. Does it help the student understand their error (if any)?
4. Is the guidance for next steps clear?

## Data Collection

### Quantitative Metrics
- Response time
- Student progression rate
- Accuracy of mathematical content
- Number of follow-up questions needed

### Qualitative Feedback
- Student satisfaction ratings
- Perceived helpfulness
- Clarity and comprehension
- Overall learning experience

### Test Scenarios

#### Scenario 1: Correct Student Work
- Student provides correct step
- Evaluate AI's recognition and encouragement
- Check for appropriate next-step guidance

#### Scenario 2: Incorrect Student Work  
- Student makes common error
- Evaluate AI's error identification
- Assess quality of corrective guidance

#### Scenario 3: Partial Understanding
- Student shows partial progress
- Evaluate AI's building on existing knowledge
- Check scaffolding effectiveness

#### Scenario 4: Complete Confusion
- Student is completely stuck
- Evaluate AI's ability to restart/redirect
- Check for appropriate simplification

## Implementation Process

### Phase 1: Baseline Testing
1. Test current prompts with standard problem set
2. Establish baseline metrics
3. Identify key areas for improvement

### Phase 2: Iterative Improvement
1. Develop alternative prompt variations
2. A/B test different approaches
3. Measure improvements in key metrics

### Phase 3: Validation
1. Test with broader range of problems
2. Validate improvements hold across categories
3. Conduct user acceptance testing

### Phase 4: Continuous Monitoring
1. Monitor performance in production
2. Collect ongoing feedback
3. Make incremental improvements

## Success Criteria

A prompt strategy is considered successful if it achieves:
- **Educational Effectiveness**: Average rating ≥ 4.0
- **Mathematical Accuracy**: 100% accuracy required
- **User Experience**: Average satisfaction ≥ 4.0
- **Consistency**: Performance stable across problem types