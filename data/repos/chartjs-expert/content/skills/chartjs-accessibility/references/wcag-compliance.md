# WCAG Compliance for Chart.js

This reference provides detailed WCAG 2.1 requirements and a compliance checklist for Chart.js visualizations.

## WCAG Levels Overview

| Level | Description | Target |
|-------|-------------|--------|
| A | Minimum accessibility | Required baseline |
| AA | Standard accessibility | Recommended for most sites |
| AAA | Enhanced accessibility | Specialized contexts |

Most organizations target **WCAG 2.1 Level AA** compliance.

## Relevant Success Criteria

### 1.1.1 Non-text Content (Level A)

Canvas charts are non-text content requiring text alternatives.

**Requirements:**

- Provide `aria-label` on canvas with chart summary
- Include fallback content inside `<canvas>` tags
- Use `aria-describedby` for detailed descriptions

**Implementation:**

```html
<canvas
  id="chart"
  role="img"
  aria-label="Sales increased 25% from Q1 to Q4 2024"
  aria-describedby="chart-desc"
>
  <!-- Fallback for no-JS or assistive tech -->
  <p>Chart showing quarterly sales growth in 2024.</p>
</canvas>
<div id="chart-desc" class="visually-hidden">
  Detailed breakdown: Q1 $12K, Q2 $19K, Q3 $15K, Q4 $25K.
</div>
```

### 1.4.1 Use of Color (Level A)

Color must not be the only visual means of conveying information.

**Requirements:**

- Use patterns, shapes, or labels in addition to color
- Ensure legend items are distinguishable without color

**Implementation:**

```javascript
// Use patterns with chartjs-plugin-patternomaly
backgroundColor: [
  pattern.draw('diagonal', '#E69F00'),
  pattern.draw('dot', '#56B4E9'),
  pattern.draw('cross', '#009E73')
]
```

### 1.4.3 Contrast (Minimum) (Level AA)

Text and graphical elements need sufficient contrast.

**Requirements:**

| Element Type | Minimum Ratio |
|--------------|---------------|
| Normal text | 4.5:1 |
| Large text (18px+ or 14px+ bold) | 3:1 |
| Graphical objects | 3:1 |
| User interface components | 3:1 |

**Implementation:**

```javascript
options: {
  plugins: {
    legend: {
      labels: { color: '#1a1a1a' }  // Dark text on light bg
    },
    title: {
      color: '#1a1a1a'
    }
  },
  scales: {
    x: {
      ticks: { color: '#1a1a1a' },
      grid: { color: '#666666' }  // 3:1 ratio minimum
    }
  }
}
```

### 1.4.11 Non-text Contrast (Level AA)

Graphical objects required to understand content need 3:1 contrast.

**Requirements:**

- Chart bars, lines, and points must have 3:1 contrast with background
- Adjacent data elements should be distinguishable

**Implementation:**

```javascript
// Ensure bar colors contrast with white background
backgroundColor: [
  '#0072B2',  // 7.3:1 on white
  '#E69F00',  // 3.0:1 on white (borderline - add border)
  '#009E73',  // 4.5:1 on white
  '#D55E00'   // 4.6:1 on white
]
```

### 2.1.1 Keyboard (Level A)

All functionality must be operable via keyboard.

**Requirements:**

- Interactive charts must support keyboard navigation
- Focus must be visible when navigating

**Implementation:**

```javascript
canvas.setAttribute('tabindex', '0');

canvas.addEventListener('keydown', (e) => {
  switch (e.key) {
    case 'ArrowRight':
      navigateToNext();
      break;
    case 'ArrowLeft':
      navigateToPrevious();
      break;
    case 'Enter':
    case ' ':
      activateCurrentElement();
      break;
    case 'Escape':
      clearSelection();
      break;
  }
});
```

### 2.3.1 Three Flashes or Below Threshold (Level A)

Content must not flash more than 3 times per second.

**Requirements:**

- Avoid rapid animations that could trigger seizures
- Animation durations should be reasonable

**Implementation:**

```javascript
options: {
  animation: {
    duration: 1000,  // 1 second - safe duration
    // Avoid repeating animations
  }
}
```

### 2.3.3 Animation from Interactions (Level AAA)

Motion animation can be disabled.

**Requirements:**

- Respect `prefers-reduced-motion` media query
- Provide option to disable animations

**Implementation:**

```javascript
const prefersReducedMotion = window.matchMedia(
  '(prefers-reduced-motion: reduce)'
).matches;

options: {
  animation: prefersReducedMotion ? false : { duration: 1000 }
}
```

### 4.1.2 Name, Role, Value (Level A)

User interface components must have accessible names and roles.

**Requirements:**

- Canvas must have `role="img"` for static charts
- Interactive charts need appropriate ARIA attributes

**Implementation:**

```html
<!-- Static chart -->
<canvas role="img" aria-label="Chart description"></canvas>

<!-- Interactive chart -->
<canvas
  role="application"
  aria-label="Interactive chart"
  aria-roledescription="data visualization"
  tabindex="0"
></canvas>
```

## Compliance Checklist

### Essential (Level A)

- [ ] Canvas has `role="img"` attribute
- [ ] Canvas has descriptive `aria-label`
- [ ] Fallback content inside `<canvas>` tags
- [ ] Color is not the only means of conveying data
- [ ] Interactive elements are keyboard accessible
- [ ] No content flashes more than 3 times/second
- [ ] Focus is visible when navigating

### Recommended (Level AA)

- [ ] Text contrast ratio is at least 4.5:1
- [ ] Graphical object contrast is at least 3:1
- [ ] `prefers-reduced-motion` is respected
- [ ] Data table alternative is available
- [ ] ARIA live regions announce changes
- [ ] Chart title is descriptive and visible

### Enhanced (Level AAA)

- [ ] Text contrast ratio is at least 7:1
- [ ] Extended audio descriptions available
- [ ] Sign language interpretation for video content
- [ ] All animations can be disabled

## Testing Tools

### Automated Testing

| Tool | Purpose |
|------|---------|
| axe DevTools | Browser extension for accessibility audits |
| WAVE | Web accessibility evaluation tool |
| Lighthouse | Chrome DevTools accessibility audit |
| Pa11y | Command-line accessibility testing |

### Manual Testing

| Method | Purpose |
|--------|---------|
| Screen reader | Test with NVDA, JAWS, or VoiceOver |
| Keyboard only | Navigate without mouse |
| Zoom to 200% | Verify content remains usable |
| Color blindness simulators | Test color accessibility |

### Color Contrast Tools

| Tool | URL |
|------|-----|
| WebAIM Contrast Checker | webaim.org/resources/contrastchecker |
| Colour Contrast Analyser | tpgi.com/color-contrast-checker |
| Stark (Figma/Sketch) | getstark.co |

## Common Pitfalls

### 1. Missing Text Alternatives

**Problem:** Canvas has no `aria-label` or fallback content.

**Solution:**

```html
<canvas
  id="chart"
  role="img"
  aria-label="Monthly revenue showing 25% growth"
>
  <p>Chart data: Jan $10K, Feb $12K, Mar $15K...</p>
</canvas>
```

### 2. Color-Only Differentiation

**Problem:** Multiple datasets distinguished only by color.

**Solution:**

```javascript
datasets: [
  {
    label: 'Product A',
    borderDash: [],           // Solid line
    pointStyle: 'circle'      // Circle markers
  },
  {
    label: 'Product B',
    borderDash: [5, 5],       // Dashed line
    pointStyle: 'rect'        // Square markers
  }
]
```

### 3. Insufficient Contrast

**Problem:** Light colors on white background.

**Solution:**

- Use contrast checker before selecting colors
- Add borders to low-contrast fills
- Prefer darker shades from accessible palettes

### 4. Ignoring Motion Preferences

**Problem:** Animations play regardless of user preferences.

**Solution:**

```javascript
const reduceMotion = window.matchMedia(
  '(prefers-reduced-motion: reduce)'
).matches;

new Chart(ctx, {
  options: {
    animation: reduceMotion ? false : { duration: 1000 }
  }
});
```

### 5. No Keyboard Access

**Problem:** Chart tooltips only accessible via mouse hover.

**Solution:** Implement keyboard navigation with arrow keys and programmatic tooltip display.

## Framework-Specific Considerations

### React (react-chartjs-2)

```jsx
<div role="img" aria-label={chartDescription}>
  <Bar data={data} options={options} />
</div>
```

### Vue (vue-chartjs)

```vue
<template>
  <div role="img" :aria-label="chartDescription">
    <Bar :data="data" :options="options" />
  </div>
</template>
```

### Angular (ng2-charts)

```html
<div role="img" [attr.aria-label]="chartDescription">
  <canvas baseChart [data]="data" [options]="options"></canvas>
</div>
```
