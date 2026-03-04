# WCAG 2.2 AA Comprehensive Checklist

Complete accessibility checklist based on Web Content Accessibility Guidelines 2.2 Level AA.

## 1. Perceivable

### 1.1 Text Alternatives (Level A)

**Images:**
- [ ] All `<img>` elements have `alt` attribute
- [ ] Informative images have descriptive alt text
- [ ] Decorative images have empty `alt=""`
- [ ] Complex images (charts, graphs) have extended descriptions
- [ ] Image alt text doesn't start with "image of" or "picture of"

**Icons:**
- [ ] Icon buttons have accessible names (aria-label or visually hidden text)
- [ ] Decorative icons hidden from screen readers (`aria-hidden="true"`)

**Code Example:**
```tsx
// Informative image
<img src="chart.png" alt="Sales increased 25% in Q4 2024" />

// Decorative image
<img src="decorative-line.png" alt="" />

// Icon button
<button aria-label="Close dialog">
  <XIcon aria-hidden="true" />
</button>
```

### 1.3 Adaptable (Level A)

**Structure:**
- [ ] Content uses semantic HTML (`<header>`, `<main>`, `<nav>`, `<footer>`)
- [ ] Headings describe content and follow hierarchy
- [ ] Lists use `<ul>`, `<ol>`, or `<dl>` appropriately
- [ ] Tables have `<th>` headers with `scope` attribute

**Forms:**
- [ ] Form inputs have associated `<label>` elements
- [ ] Required fields indicated (not by color alone)
- [ ] Related inputs grouped with `<fieldset>` and `<legend>`

**Code Example:**
```tsx
// Proper form structure
<fieldset>
  <legend>Shipping Address</legend>
  <div>
    <label htmlFor="street">Street Address</label>
    <input id="street" type="text" required aria-required="true" />
  </div>
</fieldset>
```

### 1.4 Distinguishable (Level AA)

**Color Contrast:**
- [ ] Normal text: 4.5:1 contrast ratio minimum
- [ ] Large text (18px+ or 14px+ bold): 3:1 minimum
- [ ] UI components and icons: 3:1 minimum
- [ ] Focus indicators: 3:1 against adjacent colors

**Testing Tools:**
- WebAIM Contrast Checker: https://webaim.org/resources/contrastchecker/
- Chrome DevTools: Elements > Styles > contrast ratio
- axe DevTools browser extension

**Text Spacing (1.4.12):**
- [ ] Line height at least 1.5× font size
- [ ] Paragraph spacing at least 2× font size
- [ ] Letter spacing at least 0.12× font size
- [ ] Word spacing at least 0.16× font size

**Reflow (1.4.10):**
- [ ] Content viewable at 320px width without horizontal scrolling
- [ ] No loss of content or functionality at 400% zoom

**Non-text Contrast (1.4.11):**
- [ ] Form input borders: 3:1 against background
- [ ] Custom focus indicators: 3:1 against adjacent colors
- [ ] Icons conveying information: 3:1 minimum

## 2. Operable

### 2.1 Keyboard Accessible (Level A)

**Navigation:**
- [ ] All interactive elements focusable via Tab key
- [ ] Focus order matches visual layout (logical)
- [ ] No keyboard traps (can Tab away from all elements)
- [ ] Custom widgets support expected keyboard interactions

**Skip Links:**
- [ ] "Skip to main content" link as first focusable element
- [ ] Skip link visible on focus

**Code Example:**
```tsx
// Skip link
<a href="#main-content" className="sr-only focus:not-sr-only focus:absolute focus:top-4 focus:left-4 focus:z-50 focus:bg-white focus:p-4">
  Skip to main content
</a>

<main id="main-content">
  {/* Main content */}
</main>
```

### 2.4 Navigable (Level AA)

**Page Titles:**
- [ ] Each page has unique, descriptive `<title>`
- [ ] Title format: "Page Name | Site Name"

**Focus Visible (2.4.7):**
- [ ] Focus indicator visible on all interactive elements
- [ ] Focus indicator has 3:1 contrast
- [ ] Minimum 2px outline or equivalent

**Focus Appearance (2.4.11 - WCAG 2.2):**
- [ ] Focus indicator area at least as large as 2px perimeter
- [ ] Contrast ratio of 3:1 between focused and unfocused states

**Code Example:**
```css
/* Visible focus indicator */
:focus-visible {
  outline: 2px solid #005fcc;
  outline-offset: 2px;
}

/* Remove default, add custom */
button:focus {
  outline: none;
  box-shadow: 0 0 0 3px rgba(0, 95, 204, 0.5);
}
```

**Link Purpose (2.4.4):**
- [ ] Link text describes destination
- [ ] Avoid "click here" or "read more" alone
- [ ] If generic text needed, use aria-label

**Multiple Ways (2.4.5):**
- [ ] More than one way to find pages (navigation, search, sitemap)

### 2.5 Input Modalities (Level AA)

**Target Size (2.5.8 - WCAG 2.2):**
- [ ] Interactive elements at least 24×24px CSS pixels
- [ ] Recommended: 44×44px for touch interfaces
- [ ] Inline links exempt if sufficient spacing

**Pointer Gestures (2.5.1):**
- [ ] Multipoint gestures have single-pointer alternatives
- [ ] Path-based gestures have alternatives

**Dragging Movements (2.5.7 - WCAG 2.2):**
- [ ] Drag operations have non-dragging alternatives
- [ ] Example: drag-to-reorder also has buttons

**Code Example:**
```tsx
// Minimum touch target
<button className="min-h-[44px] min-w-[44px] p-3">
  <Icon className="w-5 h-5" />
</button>
```

## 3. Understandable

### 3.1 Readable (Level AA)

**Language:**
- [ ] Page language declared: `<html lang="en">`
- [ ] Language changes marked: `<span lang="es">Hola</span>`

### 3.2 Predictable (Level AA)

**Consistent Navigation (3.2.3):**
- [ ] Navigation appears in same location across pages
- [ ] Navigation order consistent

**Consistent Identification (3.2.4):**
- [ ] Same functionality has same labels across site
- [ ] Icons used consistently

**On Input (3.2.2):**
- [ ] No unexpected context changes on input
- [ ] Form submission requires explicit action (button click)

### 3.3 Input Assistance (Level AA)

**Error Identification (3.3.1):**
- [ ] Errors clearly identified in text (not color alone)
- [ ] Error messages describe the problem
- [ ] Error messages suggest correction

**Labels or Instructions (3.3.2):**
- [ ] All form fields have visible labels
- [ ] Required fields clearly indicated
- [ ] Format requirements explained (e.g., "MM/DD/YYYY")

**Error Suggestion (3.3.3):**
- [ ] When error detected, suggest correction if known
- [ ] Example: "Email must include @"

**Error Prevention (3.3.4):**
- [ ] Legal/financial submissions: confirm, review, or reversible
- [ ] User-entered data: confirm before submission

**Code Example:**
```tsx
// Accessible error handling
<div>
  <label htmlFor="email">
    Email <span className="text-red-600">*</span>
  </label>
  <input
    id="email"
    type="email"
    aria-invalid={hasError}
    aria-describedby={hasError ? "email-error" : undefined}
  />
  {hasError && (
    <p id="email-error" className="text-red-600" role="alert">
      Please enter a valid email address (e.g., name@example.com)
    </p>
  )}
</div>
```

## 4. Robust

### 4.1 Compatible (Level AA)

**Parsing (4.1.1):**
- [ ] Valid HTML (no duplicate IDs)
- [ ] Elements properly nested
- [ ] Tags properly closed

**Name, Role, Value (4.1.2):**
- [ ] Custom widgets have appropriate ARIA roles
- [ ] State changes announced (aria-expanded, aria-pressed)
- [ ] Name/label accessible to assistive technology

**Status Messages (4.1.3):**
- [ ] Status messages announced without focus change
- [ ] Use `role="status"` or `role="alert"` appropriately
- [ ] Toast notifications accessible

**Code Example:**
```tsx
// Accessible status message
<div role="status" aria-live="polite">
  {isLoading && "Loading..."}
  {isSuccess && "Changes saved successfully"}
</div>

// Alert for errors
<div role="alert">
  {error && `Error: ${error.message}`}
</div>
```

## Quick Testing Checklist

### Keyboard Testing
1. Tab through entire page - can you reach everything?
2. Is focus visible at all times?
3. Can you activate all buttons/links with Enter/Space?
4. Can you escape from all dialogs/menus?
5. Does Tab order match visual order?

### Screen Reader Testing
1. Do images have meaningful alt text?
2. Are headings properly structured?
3. Are form fields properly labeled?
4. Are errors announced?
5. Are dynamic changes announced?

### Visual Testing
1. Does text meet contrast requirements?
2. Is focus indicator visible?
3. Do touch targets meet size requirements?
4. Does content reflow at 320px width?
5. Is information conveyed without color alone?

## Automated Testing Tools

**Browser Extensions:**
- axe DevTools (Deque)
- WAVE (WebAIM)
- Lighthouse (Chrome DevTools)

**CI/CD Integration:**
- axe-core
- pa11y
- jest-axe

**Code Example:**
```tsx
// jest-axe test
import { axe, toHaveNoViolations } from 'jest-axe';

expect.extend(toHaveNoViolations);

test('component is accessible', async () => {
  const { container } = render(<MyComponent />);
  const results = await axe(container);
  expect(results).toHaveNoViolations();
});
```
