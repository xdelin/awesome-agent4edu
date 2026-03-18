# Chart.js Dataset Properties Reference

Complete reference for all dataset properties by chart type.

## Universal Dataset Properties

These properties work across all chart types:

| Property | Type | Description |
|----------|------|-------------|
| `label` | string | Dataset label for legend/tooltip |
| `data` | array | Data values (required) |
| `order` | number | Drawing order (lower = behind) |
| `hidden` | boolean | Initially hide dataset |
| `clip` | number/object/false | Clipping relative to chart area |

## Line Dataset Properties

```javascript
{
  // Basic
  label: 'My Dataset',
  data: [1, 2, 3],

  // Line styling
  borderColor: 'rgb(75, 192, 192)',
  borderWidth: 3,
  borderDash: [5, 5],          // Dashed line
  borderDashOffset: 0,
  borderCapStyle: 'butt',      // butt, round, square
  borderJoinStyle: 'miter',    // round, bevel, miter

  // Fill
  fill: false,                 // false, true, 'origin', 'start', 'end', index
  backgroundColor: 'rgba(75, 192, 192, 0.2)',

  // Line curve
  tension: 0.1,                // 0 = straight, 0.4 = curved
  cubicInterpolationMode: 'default',  // 'default', 'monotone'
  stepped: false,              // false, true, 'before', 'after', 'middle'

  // Points
  pointStyle: 'circle',        // circle, cross, crossRot, dash, line, rect, rectRounded, rectRot, star, triangle
  pointRadius: 3,
  pointBackgroundColor: 'rgb(75, 192, 192)',
  pointBorderColor: 'rgb(75, 192, 192)',
  pointBorderWidth: 1,
  pointHitRadius: 1,
  pointRotation: 0,

  // Hover state
  pointHoverRadius: 4,
  pointHoverBackgroundColor: 'rgb(75, 192, 192)',
  pointHoverBorderColor: 'rgb(75, 192, 192)',
  pointHoverBorderWidth: 1,

  // Segments (per-segment styling)
  segment: {
    borderColor: ctx => ctx.p0.skip || ctx.p1.skip ? 'gray' : undefined,
    borderDash: ctx => ctx.p0.skip || ctx.p1.skip ? [6, 6] : undefined
  },

  // Behavior
  showLine: true,
  spanGaps: false,             // Connect over null values

  // Axes
  xAxisID: 'x',
  yAxisID: 'y',
  indexAxis: 'x'               // 'x' or 'y' for vertical
}
```

## Bar Dataset Properties

```javascript
{
  // Basic
  label: 'My Dataset',
  data: [1, 2, 3],

  // Bar styling
  backgroundColor: 'rgba(54, 162, 235, 0.5)',
  borderColor: 'rgb(54, 162, 235)',
  borderWidth: 1,
  borderRadius: 0,             // number or { topLeft, topRight, bottomLeft, bottomRight }
  borderSkipped: 'start',      // start, end, middle, bottom, left, top, right, false

  // Size
  barThickness: undefined,     // number or 'flex'
  maxBarThickness: undefined,
  minBarLength: undefined,
  barPercentage: 0.9,
  categoryPercentage: 0.8,

  // Hover state
  hoverBackgroundColor: undefined,
  hoverBorderColor: undefined,
  hoverBorderWidth: 1,
  hoverBorderRadius: 0,

  // Behavior
  grouped: true,               // Group bars at same index
  skipNull: undefined,
  base: undefined,             // Base value for the bar

  // Axes
  xAxisID: 'x',
  yAxisID: 'y',
  indexAxis: 'x',              // 'y' for horizontal bars
  stack: 'bar'                 // Stack group ID
}
```

## Pie/Doughnut Dataset Properties

```javascript
{
  // Basic
  label: 'My Dataset',
  data: [300, 50, 100],

  // Arc styling
  backgroundColor: ['#ff6384', '#36a2eb', '#ffce56'],
  borderColor: '#fff',
  borderWidth: 2,
  borderAlign: 'center',       // center, inner
  borderRadius: 0,             // number or { outerStart, outerEnd, innerStart, innerEnd }
  borderJoinStyle: undefined,
  borderDash: [],
  borderDashOffset: 0,

  // Layout
  offset: 0,                   // Arc offset
  spacing: 0,                  // Space between arcs
  weight: 1,                   // Relative thickness

  // Per-dataset overrides
  rotation: undefined,         // Starting angle
  circumference: undefined,    // Arc sweep

  // Hover state
  hoverBackgroundColor: undefined,
  hoverBorderColor: undefined,
  hoverBorderWidth: undefined,
  hoverBorderDash: undefined,
  hoverBorderDashOffset: undefined,
  hoverOffset: 0
}
```

## Radar Dataset Properties

```javascript
{
  // Basic
  label: 'My Dataset',
  data: [65, 59, 90, 81, 56],

  // Line styling
  borderColor: 'rgb(255, 99, 132)',
  borderWidth: 3,
  borderDash: [],
  borderDashOffset: 0,
  borderCapStyle: 'butt',
  borderJoinStyle: 'miter',

  // Fill
  fill: true,
  backgroundColor: 'rgba(255, 99, 132, 0.2)',

  // Points
  pointStyle: 'circle',
  pointRadius: 3,
  pointBackgroundColor: 'rgb(255, 99, 132)',
  pointBorderColor: '#fff',
  pointBorderWidth: 1,
  pointHitRadius: 1,

  // Hover
  pointHoverRadius: 4,
  pointHoverBackgroundColor: 'rgb(255, 99, 132)',
  pointHoverBorderColor: '#fff',
  pointHoverBorderWidth: 1,

  // Behavior
  tension: 0,
  spanGaps: false
}
```

## Scatter/Bubble Dataset Properties

```javascript
{
  // Basic
  label: 'My Dataset',
  data: [
    { x: 10, y: 20 },           // Scatter
    { x: 10, y: 20, r: 5 }      // Bubble (r = radius)
  ],

  // Point styling
  backgroundColor: 'rgba(255, 99, 132, 0.5)',
  borderColor: 'rgb(255, 99, 132)',
  borderWidth: 1,

  pointStyle: 'circle',
  pointRadius: 3,               // Ignored for bubble (uses r from data)
  pointHoverRadius: 4,

  // Hover
  hoverBackgroundColor: undefined,
  hoverBorderColor: undefined,
  hoverBorderWidth: 1,

  // Axes
  xAxisID: 'x',
  yAxisID: 'y'
}
```

## Polar Area Dataset Properties

```javascript
{
  // Basic
  label: 'My Dataset',
  data: [11, 16, 7, 3, 14],

  // Arc styling
  backgroundColor: [
    'rgba(255, 99, 132, 0.5)',
    'rgba(75, 192, 192, 0.5)',
    'rgba(255, 205, 86, 0.5)'
  ],
  borderColor: '#fff',
  borderWidth: 2,
  borderAlign: 'center',

  // Hover
  hoverBackgroundColor: undefined,
  hoverBorderColor: undefined,
  hoverBorderWidth: undefined
}
```

## Scriptable Options

Many properties accept functions for dynamic values:

```javascript
{
  backgroundColor: (context) => {
    const value = context.dataset.data[context.dataIndex];
    return value > 50 ? 'green' : 'red';
  },
  borderWidth: (context) => {
    return context.active ? 3 : 1;
  }
}
```

### Context Object

```javascript
{
  active: boolean,          // Is element active (hovered)
  chart: Chart,             // Chart instance
  dataIndex: number,        // Index of current data
  dataset: object,          // Current dataset
  datasetIndex: number,     // Index of current dataset
  parsed: object,           // Parsed data values
  raw: any,                 // Raw data value
  element: Element          // Visual element
}
```

## Indexable Options

Array values apply to each data point:

```javascript
{
  backgroundColor: [
    'red',    // First bar
    'blue',   // Second bar
    'green'   // Third bar
  ],
  borderWidth: [1, 2, 3]
}
```
