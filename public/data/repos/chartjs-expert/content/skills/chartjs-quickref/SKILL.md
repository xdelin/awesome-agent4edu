---
name: chartjs-quickref
description: This skill should be used when the user asks "Chart.js cheat sheet", "Chart.js quick reference", "Chart.js snippets", "common Chart.js patterns", "Chart.js copy paste", "Chart.js quick tips", "Chart.js one-liners", "quick Chart.js examples", "Chart.js recipes", or needs copy-paste ready code snippets for common Chart.js v4.5.1 patterns without deep documentation.
---

# Chart.js Quick Reference (v4.5.1)

Copy-paste snippets for common Chart.js patterns.

## Data Formatting

### Currency Labels

```javascript
ticks: {
  callback: (value) => '$' + value.toLocaleString()
}
```

### Percentage Labels

```javascript
ticks: {
  callback: (value) => value + '%'
}
```

### Abbreviated Numbers (1K, 1M)

```javascript
ticks: {
  callback: (value) => {
    if (value >= 1000000) return (value / 1000000).toFixed(1) + 'M';
    if (value >= 1000) return (value / 1000).toFixed(1) + 'K';
    return value;
  }
}
```

### Date/Time Labels

```javascript
// Requires: import 'chartjs-adapter-date-fns'; (or luxon, moment)
scales: {
  x: {
    type: 'time',
    time: {
      unit: 'day',
      displayFormats: { day: 'MMM d' }
    }
  }
}
```

## Common Configurations

### Hide Legend

```javascript
plugins: { legend: { display: false } }
```

### Hide Tooltip

```javascript
plugins: { tooltip: { enabled: false } }
```

### Hide Title

```javascript
plugins: { title: { display: false } }
```

### Start Y-axis at Zero

```javascript
scales: { y: { beginAtZero: true } }
```

### Set Axis Range

```javascript
scales: { y: { min: 0, max: 100 } }
```

### Horizontal Bar Chart

```javascript
type: 'bar',
options: { indexAxis: 'y' }
```

### Disable Animation

```javascript
animation: false
```

### Stacked Bar/Line Chart

```javascript
scales: {
  x: { stacked: true },
  y: { stacked: true }
}
```

## Styling

### Chart.js Default Palette

```javascript
const colors = [
  'rgb(255, 99, 132)',   // red
  'rgb(54, 162, 235)',   // blue
  'rgb(255, 206, 86)',   // yellow
  'rgb(75, 192, 192)',   // teal
  'rgb(153, 102, 255)',  // purple
  'rgb(255, 159, 64)'    // orange
];
```

### Semi-Transparent Backgrounds

```javascript
backgroundColor: 'rgba(54, 162, 235, 0.5)',
borderColor: 'rgb(54, 162, 235)'
```

### Hide Grid Lines

```javascript
scales: {
  x: { grid: { display: false } },
  y: { grid: { display: false } }
}
```

### Custom Font

```javascript
Chart.defaults.font.family = "'Inter', sans-serif";
Chart.defaults.font.size = 14;
```

## Responsiveness

### Fill Container

```javascript
options: {
  responsive: true,
  maintainAspectRatio: false
}
```

```html
<div style="height: 400px;">
  <canvas id="myChart"></canvas>
</div>
```

### Fixed Aspect Ratio

```javascript
options: {
  responsive: true,
  aspectRatio: 2  // width:height = 2:1
}
```

### Square Chart

```javascript
options: { aspectRatio: 1 }
```

## Dynamic Updates

### Update Data

```javascript
chart.data.datasets[0].data = [1, 2, 3, 4, 5];
chart.update();
```

### Update Without Animation

```javascript
chart.update('none');
```

### Add Data Point

```javascript
chart.data.labels.push('New Label');
chart.data.datasets[0].data.push(42);
chart.update();
```

### Remove Data Point

```javascript
chart.data.labels.pop();
chart.data.datasets[0].data.pop();
chart.update();
```

### Destroy Chart

```javascript
chart.destroy();
```

## Event Handling

### Click Handler

```javascript
options: {
  onClick: (event, elements) => {
    if (elements.length > 0) {
      const index = elements[0].index;
      console.log('Clicked:', chart.data.labels[index]);
    }
  }
}
```

### Export as Image

```javascript
const imageUrl = chart.toBase64Image();
```

## Tree-Shaking Imports

### Bar Chart

```javascript
import { Chart, BarController, BarElement, CategoryScale, LinearScale } from 'chart.js';
Chart.register(BarController, BarElement, CategoryScale, LinearScale);
```

### Line Chart

```javascript
import { Chart, LineController, LineElement, PointElement, CategoryScale, LinearScale } from 'chart.js';
Chart.register(LineController, LineElement, PointElement, CategoryScale, LinearScale);
```

### Pie/Doughnut Chart

```javascript
import { Chart, PieController, ArcElement, Legend, Tooltip } from 'chart.js';
Chart.register(PieController, ArcElement, Legend, Tooltip);
// Use DoughnutController for doughnut charts
```

### With Legend and Tooltip

```javascript
import { Legend, Tooltip } from 'chart.js';
Chart.register(Legend, Tooltip);
```

## Quick Start Templates

### Minimal Bar Chart

```javascript
new Chart(ctx, {
  type: 'bar',
  data: {
    labels: ['A', 'B', 'C'],
    datasets: [{ data: [10, 20, 30] }]
  }
});
```

### Minimal Line Chart

```javascript
new Chart(ctx, {
  type: 'line',
  data: {
    labels: ['Jan', 'Feb', 'Mar'],
    datasets: [{ data: [10, 20, 15] }]
  }
});
```

### Minimal Pie Chart

```javascript
new Chart(ctx, {
  type: 'pie',
  data: {
    labels: ['Red', 'Blue', 'Yellow'],
    datasets: [{ data: [30, 50, 20] }]
  }
});
```
