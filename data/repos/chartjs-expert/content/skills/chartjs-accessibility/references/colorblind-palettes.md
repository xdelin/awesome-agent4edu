# Colorblind-Safe Palettes for Chart.js

This reference provides ready-to-use color palettes optimized for users with color vision deficiencies (CVD), covering approximately 8% of men and 0.5% of women.

## Types of Color Blindness

| Type | Affected Colors | Prevalence |
|------|-----------------|------------|
| Deuteranopia | Red-Green (green weak) | 6% of men |
| Protanopia | Red-Green (red weak) | 2% of men |
| Tritanopia | Blue-Yellow | 0.01% of population |
| Achromatopsia | Complete color blindness | 0.003% of population |

## Okabe-Ito Palette

The most widely recommended colorblind-safe palette, designed by Masataka Okabe and Kei Ito.

```javascript
const okabeIto = {
  orange:      '#E69F00',
  skyBlue:     '#56B4E9',
  bluishGreen: '#009E73',
  yellow:      '#F0E442',
  blue:        '#0072B2',
  vermilion:   '#D55E00',
  reddishPurple: '#CC79A7',
  black:       '#000000'
};

// As array for Chart.js
const okabeItoPalette = [
  '#E69F00', '#56B4E9', '#009E73', '#F0E442',
  '#0072B2', '#D55E00', '#CC79A7', '#000000'
];
```

### Usage in Chart.js

```javascript
new Chart(ctx, {
  type: 'bar',
  data: {
    labels: ['A', 'B', 'C', 'D', 'E'],
    datasets: [{
      data: [12, 19, 15, 25, 22],
      backgroundColor: okabeItoPalette.slice(0, 5)
    }]
  }
});
```

## IBM Design Color Palette

IBM's accessibility-focused palette with carefully tested combinations.

```javascript
const ibmPalette = {
  ultramarine40: '#648FFF',
  magenta50:     '#DC267F',
  orange40:      '#FE6100',
  gold20:        '#FFB000',
  indigo50:      '#785EF0'
};

const ibmColors = [
  '#648FFF', '#DC267F', '#FE6100', '#FFB000', '#785EF0'
];
```

### Sequential Palette (Single Hue)

For heatmaps and gradients:

```javascript
const ibmBlueSequential = [
  '#EDF5FF', '#D0E2FF', '#A6C8FF', '#78A9FF',
  '#4589FF', '#0F62FE', '#0043CE', '#002D9C'
];
```

## Tableau Colorblind-Safe Palette

Tableau's 10-color palette optimized for data visualization.

```javascript
const tableauColorblind = [
  '#1f77b4', // Blue
  '#ff7f0e', // Orange
  '#2ca02c', // Green
  '#d62728', // Red
  '#9467bd', // Purple
  '#8c564b', // Brown
  '#e377c2', // Pink
  '#7f7f7f', // Gray
  '#bcbd22', // Olive
  '#17becf'  // Cyan
];
```

## Paul Tol's Palettes

Colorblind-safe palettes by Paul Tol for scientific visualization.

### Bright (7 colors)

```javascript
const tolBright = [
  '#4477AA', // Blue
  '#EE6677', // Red
  '#228833', // Green
  '#CCBB44', // Yellow
  '#66CCEE', // Cyan
  '#AA3377', // Purple
  '#BBBBBB'  // Grey
];
```

### Vibrant (7 colors)

```javascript
const tolVibrant = [
  '#EE7733', // Orange
  '#0077BB', // Blue
  '#33BBEE', // Cyan
  '#EE3377', // Magenta
  '#CC3311', // Red
  '#009988', // Teal
  '#BBBBBB'  // Grey
];
```

### Muted (10 colors)

```javascript
const tolMuted = [
  '#CC6677', // Rose
  '#332288', // Indigo
  '#DDCC77', // Sand
  '#117733', // Green
  '#88CCEE', // Cyan
  '#882255', // Wine
  '#44AA99', // Teal
  '#999933', // Olive
  '#AA4499', // Purple
  '#DDDDDD'  // Pale grey
];
```

## ColorBrewer Colorblind-Safe

Selected colorblind-safe palettes from ColorBrewer2.org.

### Qualitative (Categorical Data)

```javascript
// Set1 (8 colors) - Modified for colorblind safety
const colorBrewerSet = [
  '#e41a1c', '#377eb8', '#4daf4a', '#984ea3',
  '#ff7f00', '#ffff33', '#a65628', '#f781bf'
];

// Dark2 (8 colors) - Good for small areas
const colorBrewerDark2 = [
  '#1b9e77', '#d95f02', '#7570b3', '#e7298a',
  '#66a61e', '#e6ab02', '#a6761d', '#666666'
];
```

### Sequential (Ordered Data)

```javascript
// Blues (single hue)
const colorBrewerBlues = [
  '#f7fbff', '#deebf7', '#c6dbef', '#9ecae1',
  '#6baed6', '#4292c6', '#2171b5', '#084594'
];

// YlOrRd (yellow-orange-red)
const colorBrewerYlOrRd = [
  '#ffffcc', '#ffeda0', '#fed976', '#feb24c',
  '#fd8d3c', '#fc4e2a', '#e31a1c', '#b10026'
];
```

### Diverging (Data with Midpoint)

```javascript
// RdBu (red to blue, white center)
const colorBrewerRdBu = [
  '#b2182b', '#d6604d', '#f4a582', '#fddbc7',
  '#d1e5f0', '#92c5de', '#4393c3', '#2166ac'
];
```

## Pattern Fills with Patternomaly

When color alone is insufficient, add patterns.

### Installation

```html
<script src="https://cdn.jsdelivr.net/npm/patternomaly@1.3.2/dist/patternomaly.min.js"></script>
```

### Available Patterns

| Pattern | Best For |
|---------|----------|
| `diagonal` | Bar charts |
| `dot` | Highlighting specific bars |
| `cross` | Negative values |
| `zigzag` | Secondary data |
| `line` | Time series fills |
| `weave` | Complex patterns |
| `ring` | Circular emphasis |
| `square` | Grid-like data |

### Pattern Usage

```javascript
new Chart(ctx, {
  type: 'bar',
  data: {
    labels: ['Q1', 'Q2', 'Q3', 'Q4'],
    datasets: [{
      data: [12, 19, 15, 25],
      backgroundColor: [
        pattern.draw('diagonal', '#0072B2'),
        pattern.draw('dot', '#E69F00'),
        pattern.draw('zigzag', '#009E73'),
        pattern.draw('cross', '#D55E00')
      ]
    }]
  }
});
```

### Combining Color and Pattern

```javascript
// Create consistent styling
function createAccessibleBackground(color, patternType) {
  return pattern.draw(patternType, color);
}

const datasets = [{
  label: 'Product A',
  data: [10, 20, 30],
  backgroundColor: createAccessibleBackground('#0072B2', 'diagonal'),
  borderColor: '#0072B2',
  borderWidth: 2
}];
```

## Line Chart Differentiation

Use line styles and point shapes in addition to color.

### Line Styles

```javascript
datasets: [
  {
    label: 'Series A',
    borderColor: '#0072B2',
    borderDash: [],           // Solid
    pointStyle: 'circle'
  },
  {
    label: 'Series B',
    borderColor: '#E69F00',
    borderDash: [5, 5],       // Dashed
    pointStyle: 'rect'
  },
  {
    label: 'Series C',
    borderColor: '#009E73',
    borderDash: [10, 5],      // Long dash
    pointStyle: 'triangle'
  },
  {
    label: 'Series D',
    borderColor: '#D55E00',
    borderDash: [2, 2],       // Dotted
    pointStyle: 'star'
  }
]
```

### Available Point Styles

| Style | Shape |
|-------|-------|
| `circle` | Filled circle |
| `cross` | Cross (+) |
| `crossRot` | Rotated cross (x) |
| `dash` | Horizontal line |
| `line` | Vertical line |
| `rect` | Filled square |
| `rectRounded` | Rounded square |
| `rectRot` | Rotated square (diamond) |
| `star` | Star shape |
| `triangle` | Triangle |

## Testing Color Palettes

### Online Simulators

| Tool | URL |
|------|-----|
| Coblis | color-blindness.com/coblis-color-blindness-simulator |
| Pilestone | pilestone.com/pages/color-blindness-simulator-1 |
| Chromatic Vision Simulator | asada.website/cvsimulator/e |

### Browser Extensions

- **Colorblindly** (Chrome): Simulates 8 types of CVD
- **NoCoffee** (Chrome): Vision impairment simulator
- **Sim Daltonism** (macOS): System-wide color blindness simulation

### Testing Procedure

1. Select your color palette
2. Create a test chart with all colors visible
3. Run through each CVD type in simulator
4. Verify all data series remain distinguishable
5. Add patterns or shapes if colors are confusing

## Quick Reference: Safe Color Combinations

### 2-Color Charts

```javascript
['#0072B2', '#E69F00']  // Blue + Orange (highest contrast)
['#009E73', '#D55E00']  // Green + Vermilion
```

### 3-Color Charts

```javascript
['#0072B2', '#E69F00', '#009E73']  // Blue, Orange, Green
['#56B4E9', '#E69F00', '#CC79A7']  // Sky blue, Orange, Pink
```

### 4-Color Charts

```javascript
['#0072B2', '#E69F00', '#009E73', '#CC79A7']  // Okabe-Ito subset
['#648FFF', '#DC267F', '#FE6100', '#FFB000']  // IBM palette
```

### 5+ Color Charts

For 5 or more categories, strongly recommend adding patterns or shapes in addition to color differences.

## Complete Example

```javascript
// Accessible chart with colorblind-safe palette and patterns
const accessiblePalette = [
  { color: '#0072B2', pattern: 'diagonal' },
  { color: '#E69F00', pattern: 'dot' },
  { color: '#009E73', pattern: 'zigzag' },
  { color: '#D55E00', pattern: 'cross' },
  { color: '#CC79A7', pattern: 'line' }
];

new Chart(ctx, {
  type: 'bar',
  data: {
    labels: ['Jan', 'Feb', 'Mar', 'Apr', 'May'],
    datasets: accessiblePalette.map((item, i) => ({
      label: `Category ${i + 1}`,
      data: [Math.random() * 100],
      backgroundColor: pattern.draw(item.pattern, item.color),
      borderColor: item.color,
      borderWidth: 2
    }))
  }
});
```
