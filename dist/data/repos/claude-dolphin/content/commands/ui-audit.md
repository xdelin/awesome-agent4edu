---
description: Comprehensive UI audit for accessibility, consistency, and usability
argument-hint: "[scope: file|component|page|system]"
allowed-tools: Read, Glob, Grep, Bash(git:*)
---

# UI Audit Command

Perform a systematic UI audit on the specified scope.

**Scope options:**
- File path (e.g., `src/components/ui/button.tsx`)
- Component name (e.g., `Button`)
- Page path (e.g., `src/app/dashboard/page.tsx`)
- `system` for full codebase audit

## Audit Process

### Step 1: Scope Identification

Identify all UI-related files in the specified scope:

```
# Find component files
Glob: **/*.tsx, **/*.jsx in scope directory
# Find style files
Glob: **/*.css, **/*.scss in scope directory
# Find type definitions
Grep: interface.*Props in scope
```

### Step 2: Accessibility Audit (WCAG 2.2 AA)

For each component, check:

**Visual Accessibility:**
- [ ] Color contrast ratios meet minimums (4.5:1 body, 3:1 large text)
- [ ] Focus states visible (2px+ outline, high contrast)
- [ ] Touch targets >= 44×44px for interactive elements
- [ ] No color-only information conveyance

**Semantic Accessibility:**
- [ ] Proper heading hierarchy (no skipped levels)
- [ ] Form labels associated with inputs (htmlFor/id match)
- [ ] Button vs link distinction appropriate
- [ ] Alt text on images (empty for decorative)

**Interactive Accessibility:**
- [ ] Keyboard navigable (tabIndex, handlers)
- [ ] ARIA labels where semantic HTML insufficient
- [ ] Loading/error states announced (aria-live)
- [ ] No keyboard traps

### Step 3: Design Consistency Audit

Compare against project design system:

**Design Tokens:**
```
# Check for hardcoded colors
Grep: #[0-9a-fA-F]{3,6}|rgb\(|hsl\( in scope
# Check for hardcoded spacing
Grep: \d+px(?!.*var\() in scope
```

- [ ] Colors from design token palette
- [ ] Spacing from consistent scale (4, 8, 16, 24, 32, 48, 64px)
- [ ] Typography from defined scale
- [ ] Border radius consistent

**Component Patterns:**
- [ ] Using shadcn components where applicable
- [ ] Consistent button variants and sizes
- [ ] Uniform input/form styling
- [ ] Card/container patterns followed

### Step 4: Usability Review

Evaluate user experience quality:

- [ ] Visual hierarchy clear (primary action obvious)
- [ ] Information architecture logical
- [ ] Error states helpful (not just "Error occurred")
- [ ] Loading states present for async operations
- [ ] Empty states designed (not just blank)

### Step 5: Generate Report

Format findings as structured report:

```markdown
## UI Audit Report: [Scope]

**Audit Date:** [Current Date]
**Files Analyzed:** [Count]

### Critical Issues (Must Fix)
| Issue | Location | Line | Recommendation |
|-------|----------|------|----------------|
| [Issue description] | [file.tsx] | [42] | [How to fix] |

### Warnings (Should Fix)
| Issue | Location | Line | Recommendation |
|-------|----------|------|----------------|

### Suggestions (Nice to Have)
| Suggestion | Location | Recommendation |
|------------|----------|----------------|

### Passing Checks
- Semantic HTML structure ✓
- Form accessibility ✓
- [Other passing items]

### Metrics
- Accessibility issues: [X] critical, [Y] warnings
- Design token compliance: [X]%
- shadcn component usage: [X]%

### Priority Actions
1. [Highest priority fix]
2. [Second priority fix]
3. [Third priority fix]
```

## Special Checks by Component Type

### Buttons
- Has accessible name (text content or aria-label)
- Focus state visible
- Disabled state styled and has aria-disabled
- Loading state announced

### Forms
- All inputs have labels
- Required fields indicated (not color-only)
- Error messages linked to inputs (aria-describedby)
- Submit button present and accessible

### Modals/Dialogs
- Focus trapped within dialog
- Escape key closes dialog
- Focus returns to trigger on close
- aria-modal="true" present

### Navigation
- Current page indicated (aria-current)
- Skip link present for main content
- Consistent across pages
- Mobile-responsive

### Data Tables
- Headers use <th> with scope
- Caption or aria-label describes table
- Sortable columns have aria-sort
- Pagination accessible

## Using shadcn MCP

If shadcn MCP is available, also run:

```
search_items_in_registries: Find components that could replace custom implementations
view_items_in_registries: Compare custom vs shadcn implementations
get_audit_checklist: Run shadcn's built-in audit
```
