---
description: WCAG 2.2 AA accessibility compliance check
argument-hint: "[component-path]"
allowed-tools: Read, Grep
---

# Accessibility Check Command

Perform focused WCAG 2.2 AA compliance audit on the specified component or file.

## WCAG 2.2 AA Checklist

Read the target file and evaluate against each criterion below.

### 1. Perceivable

#### 1.1 Text Alternatives
- [ ] **1.1.1 Non-text Content** - All images have alt text
  - Check: `<img>` elements have `alt` attribute
  - Check: Icon buttons have aria-label or visible text
  - Check: Decorative images have `alt=""`

#### 1.3 Adaptable
- [ ] **1.3.1 Info and Relationships** - Semantic structure
  - Check: Headings use h1-h6 appropriately
  - Check: Lists use ul/ol/dl
  - Check: Tables have th headers
  - Check: Forms have proper label associations

- [ ] **1.3.2 Meaningful Sequence** - Reading order logical
  - Check: DOM order matches visual order

#### 1.4 Distinguishable
- [ ] **1.4.1 Use of Color** - Color not sole indicator
  - Check: Links have underline or other indicator
  - Check: Errors not indicated by color alone
  - Check: Required fields not indicated by color alone

- [ ] **1.4.3 Contrast (Minimum)** - Text contrast ratios
  - Check: Body text >= 4.5:1 against background
  - Check: Large text (18px+) >= 3:1
  - Check: Bold text (14px+) >= 3:1

- [ ] **1.4.11 Non-text Contrast** - UI element contrast
  - Check: Form inputs borders >= 3:1
  - Check: Icons >= 3:1
  - Check: Focus indicators >= 3:1

- [ ] **1.4.12 Text Spacing** - Adjustable spacing
  - Check: Content readable with 1.5x line height
  - Check: No overflow with 2x paragraph spacing

### 2. Operable

#### 2.1 Keyboard Accessible
- [ ] **2.1.1 Keyboard** - All functionality via keyboard
  - Check: Interactive elements have tabIndex or are native
  - Check: onClick handlers have onKeyDown equivalents
  - Check: Custom widgets support expected keys

- [ ] **2.1.2 No Keyboard Trap** - Can navigate away
  - Check: Focus can leave all components
  - Check: Modals have escape key handler

#### 2.4 Navigable
- [ ] **2.4.3 Focus Order** - Logical focus sequence
  - Check: tabIndex values are 0 or -1 (not positive)
  - Check: Focus order matches visual layout

- [ ] **2.4.6 Headings and Labels** - Descriptive
  - Check: Headings describe content
  - Check: Form labels describe purpose

- [ ] **2.4.7 Focus Visible** - Focus indicator visible
  - Check: :focus or :focus-visible styles defined
  - Check: outline not set to none without replacement

- [ ] **2.4.11 Focus Not Obscured** (WCAG 2.2)
  - Check: Focus indicator not hidden by other elements

#### 2.5 Input Modalities
- [ ] **2.5.3 Label in Name** - Accessible name includes visible text
  - Check: aria-label includes visible text if different

- [ ] **2.5.8 Target Size** (WCAG 2.2) - Touch targets
  - Check: Interactive elements >= 24x24px
  - Recommended: >= 44x44px for touch

### 3. Understandable

#### 3.1 Readable
- [ ] **3.1.1 Language of Page** - Language declared
  - Check: html element has lang attribute

#### 3.2 Predictable
- [ ] **3.2.1 On Focus** - No unexpected changes
  - Check: Focus doesn't trigger context change
  - Check: No auto-submit on focus

- [ ] **3.2.2 On Input** - Predictable behavior
  - Check: Input doesn't auto-submit without warning

#### 3.3 Input Assistance
- [ ] **3.3.1 Error Identification** - Errors described
  - Check: Error messages in text (not color only)
  - Check: Error linked to input (aria-describedby)

- [ ] **3.3.2 Labels or Instructions** - Guidance provided
  - Check: All inputs have visible labels
  - Check: Required fields indicated
  - Check: Format requirements explained

### 4. Robust

#### 4.1 Compatible
- [ ] **4.1.1 Parsing** - Valid markup
  - Check: No duplicate IDs
  - Check: Elements properly nested

- [ ] **4.1.2 Name, Role, Value** - ARIA correctness
  - Check: Custom widgets have appropriate roles
  - Check: State changes announced (aria-expanded, etc.)

- [ ] **4.1.3 Status Messages** - Updates announced
  - Check: role="status" or role="alert" for updates
  - Check: aria-live regions for dynamic content

## Output Format

Generate report with severity levels:

```markdown
## Accessibility Report: [Component]

### Critical (User Blocked)
| Criterion | Issue | Location | Fix |
|-----------|-------|----------|-----|
| 2.1.1 | Button not keyboard accessible | Line 42 | Add onKeyDown handler |

### Serious (Major Difficulty)
| Criterion | Issue | Location | Fix |
|-----------|-------|----------|-----|
| 1.4.3 | Text contrast 3.2:1, needs 4.5:1 | Line 18 | Darken text to #333 |

### Moderate (Inconvenient)
| Criterion | Issue | Location | Fix |
|-----------|-------|----------|-----|

### Minor (Best Practice)
| Criterion | Issue | Location | Fix |
|-----------|-------|----------|-----|

### Passing
- 1.1.1 Non-text Content ✓
- 2.4.7 Focus Visible ✓
- [etc.]

### Summary
- Critical: [X] issues
- Serious: [X] issues
- Moderate: [X] issues
- Minor: [X] issues
```

## Common Fixes

**Missing accessible name:**
```tsx
// Add aria-label
<button aria-label="Close dialog"><XIcon /></button>

// Or add visually hidden text
<button>
  <XIcon aria-hidden="true" />
  <span className="sr-only">Close dialog</span>
</button>
```

**Missing focus styles:**
```css
/* Add focus-visible styles */
.interactive:focus-visible {
  outline: 2px solid var(--ring);
  outline-offset: 2px;
}
```

**Color contrast fix:**
```css
/* Before: 3.2:1 */
color: #666;

/* After: 4.5:1+ */
color: #525252;
```

**Touch target size:**
```tsx
// Ensure minimum size
<button className="min-h-[44px] min-w-[44px] p-3">
  <Icon className="w-5 h-5" />
</button>
```
