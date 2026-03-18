---
description: Analyze design token usage and component consistency
allowed-tools: Read, Glob, Grep
---

# Design System Review Command

Analyze the codebase for design token usage, component consistency, and design system adherence.

## Analysis Process

### Step 1: Identify Design System Configuration

Find and read configuration files:

```
# Tailwind config
Glob: tailwind.config.{js,ts,mjs}

# CSS variables
Glob: **/globals.css, **/variables.css

# shadcn config
Read: components.json
```

Document:
- Defined color tokens
- Spacing scale
- Typography scale
- Border radius values
- Shadow definitions

### Step 2: Audit Design Token Usage

#### Color Audit

Search for hardcoded colors:

```
# Hex colors
Grep: #[0-9a-fA-F]{3,6} in src/components/

# RGB/RGBA
Grep: rgb\(|rgba\( in src/components/

# HSL/HSLA
Grep: hsl\(|hsla\( in src/components/
```

For each hardcoded value, record:
- File location
- Line number
- Current value
- Suggested token replacement

#### Spacing Audit

Search for hardcoded spacing:

```
# Pixel values not in CSS variables
Grep: \d+px(?!.*var\() in src/components/

# Check for inconsistent patterns
Grep: (p|m|gap|space)-\d+ in src/components/
```

Compare against defined scale:
- 4px (p-1, m-1)
- 8px (p-2, m-2)
- 12px (p-3, m-3)
- 16px (p-4, m-4)
- 24px (p-6, m-6)
- 32px (p-8, m-8)
- 48px (p-12, m-12)
- 64px (p-16, m-16)

#### Typography Audit

Search for inconsistent typography:

```
# Font sizes
Grep: font-size:|text-\[|text-xs|text-sm|text-base|text-lg|text-xl in src/

# Font weights
Grep: font-weight:|font-(thin|light|normal|medium|semibold|bold) in src/

# Line heights
Grep: line-height:|leading- in src/
```

### Step 3: Component Inventory

#### shadcn Component Usage

If shadcn MCP available:
```
search_items_in_registries: List all available components
```

Manual check:
```
# List installed shadcn components
Glob: src/components/ui/*.tsx

# Count usages of each
Grep: from "@/components/ui/[component]" in src/
```

#### Custom Component Detection

```
# Find components not in ui/ directory
Glob: src/components/**/*.tsx
# Exclude ui/ directory
```

For each custom component, check:
- Does shadcn have equivalent?
- Why was custom chosen?
- Does it follow patterns?

### Step 4: Pattern Consistency Analysis

#### Button Patterns

```
Grep: <Button in src/
```

Check for consistency:
- [ ] Same variants used for same purposes
- [ ] Icon placement consistent (left vs right)
- [ ] Size usage consistent
- [ ] Loading state pattern consistent

#### Form Patterns

```
Grep: <Input|<Select|<Checkbox|<RadioGroup in src/
```

Check for consistency:
- [ ] Label placement (above vs inline)
- [ ] Error display pattern
- [ ] Required field indication
- [ ] Help text pattern

#### Card Patterns

```
Grep: <Card in src/
```

Check for consistency:
- [ ] Padding values
- [ ] Border radius
- [ ] Shadow usage
- [ ] Header/content/footer structure

#### Modal/Dialog Patterns

```
Grep: <Dialog|<Sheet|<AlertDialog in src/
```

Check for consistency:
- [ ] When to use Dialog vs Sheet
- [ ] Header structure
- [ ] Action button placement
- [ ] Close button pattern

### Step 5: Generate Report

```markdown
## Design System Review Report

**Date:** [Current Date]
**Files Analyzed:** [Count]

### Design Token Violations

#### Hardcoded Colors
| File | Line | Value | Suggested Token |
|------|------|-------|-----------------|
| button.tsx | 42 | #3b82f6 | bg-primary |
| card.tsx | 18 | #f3f4f6 | bg-muted |

**Total:** [X] hardcoded colors
**Token Compliance:** [Y]%

#### Hardcoded Spacing
| File | Line | Value | Suggested Token |
|------|------|-------|-----------------|
| header.tsx | 24 | 20px | p-5 (or 24px → p-6) |

**Total:** [X] off-scale values
**Scale Compliance:** [Y]%

### Component Coverage

| Component Type | shadcn Available | Using shadcn | Custom |
|---------------|------------------|--------------|--------|
| Button | ✓ | ✓ | - |
| Input | ✓ | ✓ | - |
| DatePicker | ✓ | ✗ | Custom |
| DataTable | ✓ | ✗ | Custom |

**shadcn Coverage:** [X]%

### Opportunities

**Replace with shadcn:**
- [ ] `CustomDatePicker` → `@shadcn/calendar`
- [ ] `CustomSelect` → `@shadcn/select`

**Consolidate patterns:**
- [ ] 3 different card padding patterns → standardize to p-6
- [ ] 2 button icon placements → standardize to left

### Pattern Inconsistencies

| Pattern | Variations Found | Recommendation |
|---------|-----------------|----------------|
| Button icons | Left (15), Right (8) | Standardize left |
| Card padding | p-4 (10), p-6 (5), p-8 (2) | Standardize p-6 |
| Form errors | Below (12), Inline (3) | Standardize below |

### Recommendations

1. **High Priority:** Replace [X] hardcoded colors with tokens
2. **Medium Priority:** Migrate CustomDatePicker to shadcn calendar
3. **Low Priority:** Consolidate card padding patterns

### Action Items

- [ ] Create PR to fix color token violations
- [ ] Evaluate shadcn calendar for date picking needs
- [ ] Document card padding standard in style guide
```

## Quick Checks

### Color Token Check
```bash
# Count hardcoded colors
grep -r "#[0-9a-fA-F]\{3,6\}" src/components/ | wc -l
```

### Spacing Scale Check
```bash
# Find non-standard spacing
grep -rE "(p|m|gap)-\[" src/components/
```

### Component Usage
```bash
# Count shadcn component imports
grep -r "from \"@/components/ui" src/ | cut -d'"' -f2 | sort | uniq -c
```
