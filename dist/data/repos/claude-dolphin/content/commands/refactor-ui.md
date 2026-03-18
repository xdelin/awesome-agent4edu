---
description: Apply tactical UI improvements using Refactoring UI principles
argument-hint: "[component-path]"
allowed-tools: Read, Edit, Write
---

# Refactor UI Command

Apply tactical UI improvements to the specified component using Refactoring UI principles.

**Required:** Specify a component file path.

## Refactoring Process

### Step 1: Current State Analysis

Read the target component and document:

**Visual Structure:**
- Current spacing values
- Typography (sizes, weights, colors)
- Color usage
- Border and shadow usage
- Layout approach

**Hierarchy:**
- What's the primary element?
- What's secondary?
- What's tertiary?
- Is hierarchy clear?

### Step 2: Apply Tactical Improvements

#### Spacing Improvements

**Before checking:** What spacing values are currently used?

**Principles:**
- Use consistent scale: 4, 8, 12, 16, 24, 32, 48, 64px
- Start with MORE whitespace than feels necessary
- Group related items closer together
- Separate unrelated items further apart

**Common fixes:**
```tsx
// Before: Random spacing
<div className="p-5 mb-3 gap-7">

// After: Scale-based spacing
<div className="p-6 mb-4 gap-8">
```

**Padding guidance:**
- Tight (dense UI): 8-12px
- Standard: 16-24px
- Spacious: 24-32px
- Hero/cards: 32-48px

#### Typography Improvements

**Principles:**
- Maximum 2-3 font sizes per component
- Create hierarchy through WEIGHT first, then size
- Don't use pure black for text (use #111 or similar)
- Body text minimum 16px for readability

**Font size scale:**
```
text-xs:   12px - Labels, captions
text-sm:   14px - Secondary text
text-base: 16px - Body text
text-lg:   18px - Large body
text-xl:   20px - Subheadings
text-2xl:  24px - Section headings
text-3xl:  30px - Page headings
```

**Weight hierarchy:**
```tsx
// Before: Size-only hierarchy
<h2 className="text-2xl">Title</h2>
<p className="text-base">Body</p>
<span className="text-sm">Meta</span>

// After: Weight + size hierarchy
<h2 className="text-xl font-semibold text-gray-900">Title</h2>
<p className="text-base font-normal text-gray-600">Body</p>
<span className="text-sm font-normal text-gray-500">Meta</span>
```

#### Color Improvements

**Principles:**
- Use fewer colors, more intentionally
- Saturated colors for interactive elements ONLY
- Gray/neutral for most content
- Create hierarchy with opacity, not hue

**Common fixes:**
```tsx
// Before: Too many colors
<div className="text-blue-500">
  <h2 className="text-purple-600">Title</h2>
  <p className="text-green-500">Status</p>
</div>

// After: Restrained color use
<div className="text-gray-900">
  <h2 className="text-gray-900 font-semibold">Title</h2>
  <p className="text-green-600">Status</p> {/* Only status uses color */}
</div>
```

**Text color hierarchy:**
```
text-gray-900: Primary text (headings)
text-gray-700: Secondary text (body)
text-gray-500: Tertiary text (meta, captions)
text-gray-400: Disabled/placeholder
```

#### Border & Shadow Improvements

**Principles:**
- Prefer spacing over borders for separation
- Use background contrast instead of borders when possible
- Shadows: consistent light direction (top-left)
- Layer shadows for depth

**Border alternatives:**
```tsx
// Before: Border for separation
<div className="border-b">Item 1</div>
<div className="border-b">Item 2</div>

// After: Spacing for separation
<div className="pb-4">Item 1</div>
<div className="pt-4">Item 2</div>

// Or: Background contrast
<div className="bg-gray-50 rounded-lg p-4">
  <div>Item 1</div>
  <div>Item 2</div>
</div>
```

**Shadow improvements:**
```tsx
// Before: Heavy single shadow
<div className="shadow-lg">

// After: Subtle layered shadow
<div className="shadow-sm ring-1 ring-gray-200/50">
```

#### Visual Hierarchy Improvements

**Principles:**
- ONE primary action per view
- Secondary actions should visually recede
- Use size + weight + color together
- Lead the eye with alignment and whitespace

**Button hierarchy:**
```tsx
// Primary action: Most prominent
<Button variant="default">Save Changes</Button>

// Secondary action: Less prominent
<Button variant="outline">Preview</Button>

// Tertiary action: Least prominent
<Button variant="ghost">Cancel</Button>
```

### Step 3: Accessibility Verification

After refactoring, verify:

- [ ] **Contrast maintained:** Check new colors meet 4.5:1 for text
- [ ] **Focus states preserved:** Ensure focus indicators still visible
- [ ] **Touch targets adequate:** Interactive elements still >= 44x44px
- [ ] **Semantic structure intact:** No accidental hierarchy breaks

### Step 4: Present Changes

Output format:

```markdown
## UI Refactoring: [Component Name]

### Summary
[Brief description of changes made]

### Before/After Comparison

#### Spacing
**Before:**
\`\`\`tsx
<div className="p-5 gap-7">
\`\`\`

**After:**
\`\`\`tsx
<div className="p-6 gap-8">  // Aligned to 8px scale
\`\`\`

**Rationale:** Aligned to 8px spacing scale for consistency.

#### Typography
**Before:**
\`\`\`tsx
<h2 className="text-2xl">Title</h2>
<p className="text-base text-black">Description</p>
\`\`\`

**After:**
\`\`\`tsx
<h2 className="text-xl font-semibold text-gray-900">Title</h2>
<p className="text-base text-gray-600">Description</p>
\`\`\`

**Rationale:** Added weight hierarchy, softened text colors for visual depth.

#### Colors
[If applicable]

#### Shadows/Borders
[If applicable]

### Accessibility Check
- [x] Contrast ratios verified
- [x] Focus states preserved
- [x] Touch targets adequate

### Full Diff
\`\`\`diff
- <div className="p-5 gap-7">
+ <div className="p-6 gap-8">
-   <h2 className="text-2xl">Title</h2>
+   <h2 className="text-xl font-semibold text-gray-900">Title</h2>
...
\`\`\`
```

## Quick Reference

### Spacing Scale
| Class | Pixels | Use Case |
|-------|--------|----------|
| p-1, gap-1 | 4px | Tight grouping |
| p-2, gap-2 | 8px | Small gaps |
| p-3, gap-3 | 12px | Internal padding |
| p-4, gap-4 | 16px | Standard |
| p-6, gap-6 | 24px | Section gaps |
| p-8, gap-8 | 32px | Major separation |

### Typography Hierarchy
| Role | Size | Weight | Color |
|------|------|--------|-------|
| Heading | xl-2xl | semibold-bold | gray-900 |
| Body | base | normal | gray-600-700 |
| Secondary | sm | normal | gray-500 |
| Caption | xs | normal | gray-400-500 |

### Color Roles
| Role | Usage |
|------|-------|
| gray-900 | Primary text, headings |
| gray-700 | Body text |
| gray-500 | Secondary text, metadata |
| primary | Interactive elements, links |
| destructive | Errors, delete actions |
| success | Success states |
