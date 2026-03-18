# 5-Step UX Audit Framework

Systematic methodology for evaluating user interface quality, accessibility, and design consistency.

## Overview

A UX audit is a systematic, data-driven assessment of a product's usability, accessibility, and design quality. This framework provides a repeatable process for identifying issues and prioritizing improvements.

## Step 1: Scope Definition

### Define Boundaries

**System-wide audit:**
- All user-facing components
- Complete user journeys
- Cross-cutting concerns (navigation, error handling)

**Section audit:**
- Single page or feature
- Specific component library
- Particular user flow

### Identify Stakeholders

- Who uses this interface?
- What are their goals?
- What are their constraints (devices, abilities)?

### Set Success Criteria

- WCAG 2.2 AA compliance
- Design system adherence
- Specific usability metrics

### Output: Scope Document

```markdown
## Audit Scope

**Target:** [Component/Page/System]
**User Personas:** [List primary users]
**Goals:** [What users try to accomplish]
**Criteria:**
- WCAG 2.2 AA compliance
- Design token usage > 95%
- Touch target compliance 100%
```

## Step 2: Heuristic Evaluation

Evaluate against established usability principles.

### Nielsen's 10 Heuristics

#### 1. Visibility of System Status
- [ ] User knows current state (page, form progress)
- [ ] Loading states clearly indicated
- [ ] Success/error feedback immediate
- [ ] Progress indicators for long operations

**Check:** Submit a form - is feedback immediate and clear?

#### 2. Match Between System and Real World
- [ ] Language familiar to users (not technical jargon)
- [ ] Icons and symbols widely recognized
- [ ] Information in natural, logical order
- [ ] Metaphors make sense

**Check:** Would a first-time user understand all labels?

#### 3. User Control and Freedom
- [ ] Easy to undo actions
- [ ] Easy to exit unwanted states
- [ ] Cancel buttons on forms/dialogs
- [ ] Back navigation always works

**Check:** Can user recover from any mistake easily?

#### 4. Consistency and Standards
- [ ] Same action = same result everywhere
- [ ] Terminology consistent throughout
- [ ] UI patterns match platform conventions
- [ ] Visual styling consistent

**Check:** Do similar elements behave the same way?

#### 5. Error Prevention
- [ ] Confirmation for destructive actions
- [ ] Validation before submission
- [ ] Good defaults reduce errors
- [ ] Constraints prevent invalid input

**Check:** What happens if user enters invalid data?

#### 6. Recognition Rather Than Recall
- [ ] Options visible, not hidden
- [ ] Help and instructions accessible
- [ ] Recent actions/items remembered
- [ ] Context maintained across sessions

**Check:** Does user need to remember anything from previous pages?

#### 7. Flexibility and Efficiency of Use
- [ ] Shortcuts for expert users
- [ ] Customization options available
- [ ] Frequent actions easily accessible
- [ ] Default settings sensible

**Check:** Can experienced users work faster?

#### 8. Aesthetic and Minimalist Design
- [ ] Only relevant information shown
- [ ] Visual hierarchy clear
- [ ] No clutter or noise
- [ ] White space used effectively

**Check:** Is every element necessary?

#### 9. Help Users Recognize, Diagnose, and Recover from Errors
- [ ] Error messages in plain language
- [ ] Error clearly indicates problem
- [ ] Solution suggested
- [ ] No error codes/technical details

**Check:** Are error messages helpful?

#### 10. Help and Documentation
- [ ] Help available when needed
- [ ] Task-focused documentation
- [ ] Search functionality works
- [ ] Common questions addressed

**Check:** Can users find help for common tasks?

### Output: Heuristic Findings

| Heuristic | Score (1-5) | Issues Found |
|-----------|-------------|--------------|
| Visibility | 4 | Loading states missing on 2 pages |
| Match | 5 | No issues |
| ... | ... | ... |

## Step 3: Accessibility Audit

Systematic check against WCAG 2.2 AA criteria.

### Automated Testing

Run these tools on every page:

**axe DevTools:**
```bash
# Browser extension or CLI
npx @axe-core/cli https://yoursite.com
```

**Lighthouse:**
```bash
# In Chrome DevTools > Lighthouse > Accessibility
# Or CLI:
npx lighthouse https://yoursite.com --only-categories=accessibility
```

**WAVE:**
- Browser extension for visual feedback
- Shows errors in context

### Manual Testing

**Keyboard Navigation:**
1. Tab through entire page
2. Verify focus visible at all times
3. Check logical tab order
4. Test all interactive elements with Enter/Space
5. Ensure no keyboard traps

**Screen Reader Testing:**
1. Navigate with headings (H key)
2. Navigate with landmarks (D key)
3. Check form labels announced
4. Verify images described
5. Test dynamic content announcements

**Visual Testing:**
1. Zoom to 200% - content still usable?
2. Resize to 320px width - no horizontal scroll?
3. High contrast mode - content visible?
4. Color blindness simulation - info clear?

### Output: Accessibility Report

| Criterion | Status | Issues |
|-----------|--------|--------|
| 1.1.1 Non-text Content | Pass | - |
| 1.4.3 Contrast (Minimum) | Fail | 3 buttons below 4.5:1 |
| 2.1.1 Keyboard | Partial | Dropdown not keyboard accessible |
| ... | ... | ... |

## Step 4: Consistency Analysis

Check design system adherence and pattern consistency.

### Design Token Audit

**Colors:**
```bash
# Search for hardcoded colors
grep -r "#[0-9a-fA-F]\{3,6\}" src/components/
grep -r "rgb\|hsl\|rgba" src/components/
```

Count: [X] hardcoded values found

**Spacing:**
```bash
# Search for hardcoded spacing
grep -r "px\|rem\|em" src/components/ | grep -v "var(--"
```

**Typography:**
- Font sizes from scale?
- Font weights from scale?
- Line heights consistent?

### Component Inventory

| Component | shadcn Available | Using shadcn | Custom |
|-----------|-----------------|--------------|--------|
| Button | Yes | Yes | - |
| Input | Yes | Yes | - |
| DatePicker | Yes | No | Custom |
| ... | ... | ... | ... |

**Opportunities:**
- [ ] Replace custom DatePicker with shadcn
- [ ] Consolidate 3 different card styles

### Pattern Consistency

**Buttons:**
- Same variants used consistently?
- Icon placement consistent?
- Loading states consistent?

**Forms:**
- Label placement consistent?
- Error display consistent?
- Required field indication consistent?

**Navigation:**
- Same patterns on all pages?
- Active states consistent?
- Mobile menu behavior consistent?

### Output: Consistency Report

```markdown
## Consistency Issues

### Critical (Must Fix)
- 12 hardcoded color values in components/

### High (Should Fix)
- Custom DatePicker should use shadcn calendar
- Button icon placement inconsistent (3 patterns)

### Medium (Nice to Have)
- Card border-radius varies (4px, 8px, 12px)
```

## Step 5: Prioritization & Recommendations

### Severity Classification

**Critical (Blocker):**
- Users cannot complete task
- Accessibility violation blocking users
- Example: Form submit button not keyboard accessible

**High (Major):**
- Significant difficulty or confusion
- Legal/compliance risk
- Example: Contrast ratio 3:1 instead of 4.5:1

**Medium (Moderate):**
- User can complete task but with friction
- Inconsistency that causes confusion
- Example: Different error message formats

**Low (Minor):**
- Polish opportunity
- Minor inconsistency
- Example: Slightly different border radius

### Effort Estimation

**Quick Win (< 1 hour):**
- Add alt text
- Fix single color value
- Add aria-label

**Small (1-4 hours):**
- Replace hardcoded colors with tokens
- Add focus states
- Fix form labels

**Medium (1-2 days):**
- Replace custom component with shadcn
- Implement loading states
- Add skip links

**Large (> 2 days):**
- Major refactor
- New component system
- Accessibility overhaul

### Priority Matrix

|              | Quick Win | Small | Medium | Large |
|--------------|-----------|-------|--------|-------|
| **Critical** | Do Now    | Do Now| Sprint | Plan  |
| **High**     | Do Now    | Sprint| Sprint | Backlog|
| **Medium**   | Sprint    | Sprint| Backlog| Backlog|
| **Low**      | Backlog   | Backlog| Later | Later |

### Final Report Template

```markdown
# UX Audit Report: [Scope]

**Date:** [Date]
**Auditor:** Claude Code + UX/UI Design System Skill

## Executive Summary
[2-3 sentence overview of findings]

## Critical Issues (Must Fix Immediately)
| Issue | Location | Recommendation | Effort |
|-------|----------|----------------|--------|

## High Priority (Fix This Sprint)
| Issue | Location | Recommendation | Effort |
|-------|----------|----------------|--------|

## Medium Priority (Backlog)
| Issue | Location | Recommendation | Effort |
|-------|----------|----------------|--------|

## Low Priority (Future Consideration)
| Issue | Location | Recommendation | Effort |
|-------|----------|----------------|--------|

## Passing Checks
- [List of areas that passed audit]

## Metrics
- WCAG violations: [X] critical, [Y] serious
- Design token compliance: [X]%
- Component coverage: [X]% using shadcn

## Recommendations
1. [Top priority action]
2. [Second priority action]
3. [Third priority action]
```

## Ongoing Maintenance

### Regular Audits
- Full audit: Quarterly
- Accessibility scan: After each release
- Consistency check: With each new component

### Automated Checks
- axe-core in CI/CD pipeline
- Lighthouse performance budget
- Visual regression testing

### Team Education
- Share findings with team
- Document patterns in style guide
- Review new components against checklist
