---
name: chartjs-accessibility
description: This skill should be used when the user asks "Chart.js accessibility", "accessible charts", "WCAG Chart.js", "Chart.js ARIA labels", "Chart.js screen reader", "colorblind-safe charts", "Chart.js keyboard navigation", "Chart.js reduced motion", "prefers-reduced-motion Chart.js", "Chart.js alt text", "Chart.js color contrast", "accessible data visualization", "Chart.js focus management", "Chart.js fallback content", "aria-label canvas", or needs help making Chart.js visualizations accessible and WCAG-compliant in v4.5.1.
---

# Chart.js Accessibility (v4.5.1)

Making charts accessible ensures all users, including those with disabilities, can understand your data visualizations. This guide covers WCAG compliance, screen reader support, color accessibility, keyboard navigation, and motion sensitivity.

## Canvas Accessibility Basics

Canvas elements are inherently inaccessible to screen readers. Add proper ARIA attributes and fallback content.

### ARIA Labels and Roles

```html
<canvas
  id="myChart"
  role="img"
  aria-label="Bar chart showing quarterly sales: Q1 $12,000, Q2 $19,000, Q3 $15,000, Q4 $25,000"
></canvas>
```

### Fallback Content

Provide content inside the canvas tag for screen readers and browsers without canvas support:

```html
<canvas id="myChart" role="img" aria-label="Sales comparison chart">
  <p>Quarterly sales data: Q1 had $12,000, Q2 had $19,000, Q3 had $15,000, Q4 had $25,000. Q4 showed the highest performance.</p>
</canvas>
```

### Using aria-describedby

Link to a detailed data summary:

```html
<canvas
  id="myChart"
  role="img"
  aria-label="Monthly revenue chart"
  aria-describedby="chart-description"
></canvas>
<div id="chart-description" class="visually-hidden">
  Revenue grew steadily from January ($10K) to December ($45K),
  with a notable spike in November during the holiday season.
</div>
```

```css
.visually-hidden {
  position: absolute;
  width: 1px;
  height: 1px;
  padding: 0;
  margin: -1px;
  overflow: hidden;
  clip: rect(0, 0, 0, 0);
  white-space: nowrap;
  border: 0;
}
```

## Color Accessibility

Ensure charts are usable by people with color vision deficiencies.

### Colorblind-Safe Palettes

The Okabe-Ito palette is optimized for all types of color blindness:

```javascript
const okabe_ito = [
  '#E69F00', // Orange
  '#56B4E9', // Sky Blue
  '#009E73', // Bluish Green
  '#F0E442', // Yellow
  '#0072B2', // Blue
  '#D55E00', // Vermilion
  '#CC79A7', // Reddish Purple
  '#000000'  // Black
];

new Chart(ctx, {
  type: 'bar',
  data: {
    labels: ['A', 'B', 'C', 'D'],
    datasets: [{
      data: [12, 19, 15, 25],
      backgroundColor: okabe_ito.slice(0, 4)
    }]
  }
});
```

### Pattern Fills with chartjs-plugin-patternomaly

Add patterns to distinguish data without relying on color:

```html
<script src="https://cdn.jsdelivr.net/npm/patternomaly@1.3.2/dist/patternomaly.min.js"></script>
```

```javascript
new Chart(ctx, {
  type: 'bar',
  data: {
    labels: ['Q1', 'Q2', 'Q3', 'Q4'],
    datasets: [{
      data: [12, 19, 15, 25],
      backgroundColor: [
        pattern.draw('diagonal', '#E69F00'),
        pattern.draw('dot', '#56B4E9'),
        pattern.draw('zigzag', '#009E73'),
        pattern.draw('cross', '#0072B2')
      ]
    }]
  }
});
```

Available patterns: `plus`, `cross`, `dash`, `cross-dash`, `dot`, `dot-dash`, `disc`, `ring`, `line`, `line-vertical`, `weave`, `zigzag`, `zigzag-vertical`, `diagonal`, `diagonal-right-left`, `square`, `box`, `triangle`, `triangle-inverted`, `diamond`, `diamond-box`.

### Contrast Requirements

WCAG AA requires:

| Element | Minimum Ratio |
|---------|---------------|
| Text (normal) | 4.5:1 |
| Text (large 18px+) | 3:1 |
| Graphical objects | 3:1 |
| UI components | 3:1 |

Ensure data elements have sufficient contrast against background:

```javascript
options: {
  plugins: {
    legend: {
      labels: {
        color: '#1a1a1a' // High contrast text
      }
    }
  },
  scales: {
    x: {
      ticks: { color: '#1a1a1a' },
      grid: { color: '#666666' }
    },
    y: {
      ticks: { color: '#1a1a1a' },
      grid: { color: '#666666' }
    }
  }
}
```

## Keyboard Navigation

Enable keyboard access for interactive chart elements.

### Focus Management for Tooltips

```javascript
const chart = new Chart(ctx, {
  type: 'bar',
  data: chartData,
  options: {
    plugins: {
      tooltip: {
        enabled: true
      }
    }
  }
});

// Keyboard navigation
let currentIndex = 0;
const dataLength = chart.data.datasets[0].data.length;

document.addEventListener('keydown', (e) => {
  if (e.key === 'ArrowRight') {
    currentIndex = (currentIndex + 1) % dataLength;
    showTooltipAtIndex(currentIndex);
  } else if (e.key === 'ArrowLeft') {
    currentIndex = (currentIndex - 1 + dataLength) % dataLength;
    showTooltipAtIndex(currentIndex);
  } else if (e.key === 'Escape') {
    chart.setActiveElements([]);
    chart.tooltip.setActiveElements([], { x: 0, y: 0 });
    chart.update();
  }
});

function showTooltipAtIndex(index) {
  chart.setActiveElements([{ datasetIndex: 0, index }]);
  chart.tooltip.setActiveElements(
    [{ datasetIndex: 0, index }],
    { x: 0, y: 0 }
  );
  chart.update();

  // Announce to screen readers
  announceDataPoint(index);
}

function announceDataPoint(index) {
  const label = chart.data.labels[index];
  const value = chart.data.datasets[0].data[index];
  const announcement = `${label}: ${value}`;

  const liveRegion = document.getElementById('chart-live-region');
  liveRegion.textContent = announcement;
}
```

```html
<div id="chart-live-region" aria-live="polite" class="visually-hidden"></div>
```

### Focusable Canvas

Make the canvas focusable for keyboard users:

```html
<canvas id="myChart" tabindex="0" role="img" aria-label="Interactive sales chart. Use arrow keys to navigate data points."></canvas>
```

## Motion Sensitivity

Respect user preferences for reduced motion to prevent discomfort.

### Detecting prefers-reduced-motion

```javascript
const prefersReducedMotion = window.matchMedia(
  '(prefers-reduced-motion: reduce)'
).matches;

new Chart(ctx, {
  type: 'bar',
  data: chartData,
  options: {
    animation: prefersReducedMotion ? false : {
      duration: 1000,
      easing: 'easeOutQuart'
    },
    transitions: {
      active: {
        animation: {
          duration: prefersReducedMotion ? 0 : 300
        }
      }
    }
  }
});
```

### Listening for Preference Changes

```javascript
const motionQuery = window.matchMedia('(prefers-reduced-motion: reduce)');

motionQuery.addEventListener('change', (e) => {
  chart.options.animation = e.matches ? false : { duration: 1000 };
  chart.update();
});
```

### Graceful Animation Reduction

Instead of disabling completely, reduce animation:

```javascript
const animationConfig = prefersReducedMotion
  ? { duration: 0 }
  : {
      duration: 1000,
      easing: 'easeOutQuart',
      delay: (context) => context.dataIndex * 100
    };
```

## Data Table Alternatives

Provide tabular data alongside or as alternative to charts.

### Hidden Data Table

```html
<figure>
  <canvas id="myChart" role="img" aria-label="Quarterly sales chart"></canvas>
  <figcaption>
    <details>
      <summary>View data table</summary>
      <table>
        <caption>Quarterly Sales Data 2024</caption>
        <thead>
          <tr>
            <th scope="col">Quarter</th>
            <th scope="col">Sales ($)</th>
          </tr>
        </thead>
        <tbody>
          <tr><td>Q1</td><td>12,000</td></tr>
          <tr><td>Q2</td><td>19,000</td></tr>
          <tr><td>Q3</td><td>15,000</td></tr>
          <tr><td>Q4</td><td>25,000</td></tr>
        </tbody>
      </table>
    </details>
  </figcaption>
</figure>
```

### Generating Tables from Chart Data

```javascript
function generateDataTable(chart) {
  const labels = chart.data.labels;
  const datasets = chart.data.datasets;

  let html = '<table><thead><tr><th scope="col">Category</th>';
  datasets.forEach(ds => {
    html += `<th scope="col">${ds.label}</th>`;
  });
  html += '</tr></thead><tbody>';

  labels.forEach((label, i) => {
    html += `<tr><td>${label}</td>`;
    datasets.forEach(ds => {
      html += `<td>${ds.data[i]}</td>`;
    });
    html += '</tr>';
  });

  html += '</tbody></table>';
  return html;
}
```

### Screen Reader Summary Statistics

```html
<div id="chart-summary" class="visually-hidden">
  Chart summary: Sales ranged from $12,000 (Q1) to $25,000 (Q4).
  Average: $17,750. Total: $71,000. Trend: upward.
</div>
```

## Complete Accessible Chart Example

```javascript
// Check motion preference
const prefersReducedMotion = window.matchMedia(
  '(prefers-reduced-motion: reduce)'
).matches;

// Colorblind-safe palette
const accessibleColors = [
  '#0072B2', '#E69F00', '#009E73', '#D55E00'
];

const chart = new Chart(
  document.getElementById('accessibleChart'),
  {
    type: 'bar',
    data: {
      labels: ['Q1', 'Q2', 'Q3', 'Q4'],
      datasets: [{
        label: 'Sales ($)',
        data: [12000, 19000, 15000, 25000],
        backgroundColor: accessibleColors
      }]
    },
    options: {
      responsive: true,
      animation: prefersReducedMotion ? false : { duration: 1000 },
      plugins: {
        legend: {
          labels: { color: '#1a1a1a' }
        },
        title: {
          display: true,
          text: 'Quarterly Sales 2024',
          color: '#1a1a1a'
        }
      },
      scales: {
        x: { ticks: { color: '#1a1a1a' } },
        y: {
          ticks: {
            color: '#1a1a1a',
            callback: (v) => '$' + v.toLocaleString()
          },
          beginAtZero: true
        }
      }
    }
  }
);
```

```html
<canvas
  id="accessibleChart"
  role="img"
  tabindex="0"
  aria-label="Bar chart: Quarterly Sales 2024. Q1: $12,000, Q2: $19,000, Q3: $15,000, Q4: $25,000. Use arrow keys to navigate."
  aria-describedby="chart-details"
></canvas>
<div id="chart-details" class="visually-hidden">
  Sales increased overall through 2024, with Q4 showing the highest
  revenue at $25,000, representing a 67% increase from Q3.
</div>
```

## Additional Resources

- See `references/wcag-compliance.md` for detailed WCAG requirements and compliance checklist
- See `references/colorblind-palettes.md` for complete colorblind-safe palette options
- See `examples/accessible-bar-chart.html` for a complete working accessible chart
