# Chart.js Configuration Patterns Reference

Advanced patterns for configuring and customizing Chart.js through global defaults, helper functions, and common configuration strategies.

## Configuration Object Anatomy

Every Chart.js chart follows this structure:

```javascript
const config = {
  type: 'bar',           // Required: Chart type
  data: {                // Required: Data configuration
    labels: [],          // Data point labels
    datasets: [{         // Array of datasets
      label: '',         // Dataset label
      data: [],          // Data values
      // Dataset-specific options
    }]
  },
  options: {             // Optional: Chart options
    responsive: true,
    plugins: {},
    scales: {},
    // ... other options
  },
  plugins: []            // Optional: Custom plugins
};

const chart = new Chart(ctx, config);
```

### Configuration Hierarchy

Chart.js applies configuration in this order (last wins):

1. **Chart.js defaults** - Built-in defaults
2. **Global defaults** - Set via `Chart.defaults`
3. **Chart type defaults** - Type-specific defaults
4. **Dataset defaults** - Dataset-level configuration
5. **Chart options** - Options passed to `new Chart()`
6. **Dataset options** - Dataset-specific overrides

## Global Defaults

### Setting Global Defaults

Modify `Chart.defaults` to apply settings to all future charts:

```javascript
// Set default font family for all charts
Chart.defaults.font.family = "'Inter', 'Helvetica', 'Arial', sans-serif";
Chart.defaults.font.size = 14;

// Set default colors
Chart.defaults.color = '#666';
Chart.defaults.borderColor = 'rgba(0, 0, 0, 0.1)';

// Set default interaction mode
Chart.defaults.interaction.mode = 'nearest';
Chart.defaults.interaction.intersect = false;

// Disable animations globally
Chart.defaults.animation = false;

// Set default responsive behavior
Chart.defaults.responsive = true;
Chart.defaults.maintainAspectRatio = true;
```

**Important:** Changes to `Chart.defaults` only affect charts created **after** the change.

### Chart Type Defaults

Set defaults for specific chart types:

```javascript
// All line charts will have these settings
Chart.defaults.datasets.line.borderWidth = 2;
Chart.defaults.datasets.line.fill = false;
Chart.defaults.datasets.line.tension = 0.4;

// All bar charts will have these settings
Chart.defaults.datasets.bar.borderWidth = 0;
Chart.defaults.datasets.bar.borderRadius = 4;

// All bubble charts
Chart.defaults.datasets.bubble.borderWidth = 1;
Chart.defaults.datasets.bubble.radius = 5;
```

### Scale Defaults

Configure default scale options:

```javascript
// All linear scales
Chart.defaults.scales.linear.min = 0;
Chart.defaults.scales.linear.grid.display = true;
Chart.defaults.scales.linear.grid.color = 'rgba(0, 0, 0, 0.05)';

// All category scales
Chart.defaults.scales.category.grid.display = false;

// All radial linear scales (radar/polar)
Chart.defaults.scales.radialLinear.angleLines.display = true;
```

### Plugin Defaults

Set default plugin configuration:

```javascript
// Legend defaults
Chart.defaults.plugins.legend.display = true;
Chart.defaults.plugins.legend.position = 'top';
Chart.defaults.plugins.legend.align = 'center';

// Tooltip defaults
Chart.defaults.plugins.tooltip.enabled = true;
Chart.defaults.plugins.tooltip.mode = 'nearest';
Chart.defaults.plugins.tooltip.backgroundColor = 'rgba(0, 0, 0, 0.8)';
Chart.defaults.plugins.tooltip.padding = 12;
Chart.defaults.plugins.tooltip.cornerRadius = 4;

// Title defaults
Chart.defaults.plugins.title.display = false;
Chart.defaults.plugins.title.font.size = 18;
Chart.defaults.plugins.title.font.weight = 'bold';
```

## Helper Functions

Chart.js provides helper functions for common operations. Import them separately:

```javascript
import Chart from 'chart.js/auto';
import { getRelativePosition } from 'chart.js/helpers';
```

### Available Helpers

| Helper | Purpose | Import |
|--------|---------|--------|
| `getRelativePosition(event, chart)` | Get canvas coordinates from event | `import { getRelativePosition } from 'chart.js/helpers'` |
| `noop()` | Empty function | `import { noop } from 'chart.js/helpers'` |
| `uid()` | Generate unique ID | `import { uid } from 'chart.js/helpers'` |
| `isNullOrUndef(value)` | Check null/undefined | `import { isNullOrUndef } from 'chart.js/helpers'` |
| `isArray(value)` | Check if array | `import { isArray } from 'chart.js/helpers'` |
| `isObject(value)` | Check if object | `import { isObject } from 'chart.js/helpers'` |
| `valueOrDefault(value, defaultValue)` | Return value or default | `import { valueOrDefault } from 'chart.js/helpers'` |

### Converting Events to Data Values

Use `getRelativePosition` to convert mouse/touch events to chart data coordinates:

```javascript
import Chart from 'chart.js/auto';
import { getRelativePosition } from 'chart.js/helpers';

const chart = new Chart(ctx, {
  type: 'line',
  data: data,
  options: {
    onClick: (e) => {
      // Get click position relative to chart
      const canvasPosition = getRelativePosition(e, chart);

      // Convert to data values using scales
      const dataX = chart.scales.x.getValueForPixel(canvasPosition.x);
      const dataY = chart.scales.y.getValueForPixel(canvasPosition.y);

      console.log(`Clicked at data point: (${dataX}, ${dataY})`);
    }
  }
});
```

### Custom Tooltip Positioning

```javascript
import { getRelativePosition } from 'chart.js/helpers';

options: {
  plugins: {
    tooltip: {
      callbacks: {
        title: function(context) {
          const position = getRelativePosition(context[0].event, chart);
          return `Position: ${position.x}, ${position.y}`;
        }
      }
    }
  }
}
```

## Responsiveness Patterns

### Basic Responsive Configuration

```javascript
options: {
  responsive: true,              // Chart resizes with container (default: true)
  maintainAspectRatio: true,    // Maintain aspect ratio (default: true)
  aspectRatio: 2,               // Width/height ratio (default: 2, radial: 1)
  resizeDelay: 0,               // Debounce resize updates in ms (default: 0)

  // Callback when chart resizes
  onResize: (chart, size) => {
    console.log(`Chart resized to ${size.width}x${size.height}`);
  }
}
```

### Responsive with Custom Aspect Ratio

```javascript
// Square chart (common for radar, polar area)
options: {
  aspectRatio: 1
}

// Wide chart (good for timelines)
options: {
  aspectRatio: 3
}

// Tall chart (good for ranking/leaderboard)
options: {
  aspectRatio: 0.5
}
```

### Fixed Size (Non-Responsive)

```javascript
options: {
  responsive: false,
  maintainAspectRatio: false
}
```

Then set canvas dimensions directly:

```html
<canvas id="myChart" width="600" height="400"></canvas>
```

### Responsive Container Sizing

Chart.js gets dimensions from the canvas's **parent container**:

```html
<!-- Good: Chart takes 80% of container width -->
<div style="width: 80%; max-width: 1200px;">
  <canvas id="chart1"></canvas>
</div>

<!-- Bad: No container sizing -->
<canvas id="chart4"></canvas>
```

### Flexbox/Grid Layout

When using flexbox or grid, set `min-width: 0` on chart containers to prevent overflow:

```html
<!-- Grid layout with proper sizing -->
<div style="display: grid; grid-template-columns: 1fr 1fr;">
  <div style="min-width: 0;"><canvas id="chart1"></canvas></div>
  <div style="min-width: 0;"><canvas id="chart2"></canvas></div>
</div>

<!-- Flexbox layout with proper sizing -->
<div style="display: flex; gap: 20px;">
  <div style="flex: 1; min-width: 0;"><canvas id="chart3"></canvas></div>
  <div style="flex: 1; min-width: 0;"><canvas id="chart4"></canvas></div>
</div>
```

**Why `min-width: 0`?** Flex and grid items have an implicit `min-width: auto` which prevents them from shrinking below their content size. Charts can cause container overflow without this fix.

### Responsive Font Sizes

Use Chart.js responsive font sizing (RFS):

```javascript
Chart.defaults.font.size = 14;  // Base size

// Automatically scales down on smaller screens
options: {
  plugins: {
    title: {
      font: {
        size: 18  // Scales responsively
      }
    }
  }
}
```

### Media Query-Like Behavior

Adjust chart options based on container size using the `onResize` callback:

```javascript
options: {
  onResize: (chart, size) => {
    if (size.width < 600) {
      // Mobile: hide legend, reduce padding
      chart.options.plugins.legend.display = false;
      chart.options.layout.padding = 10;
    } else {
      // Desktop: show legend, increase padding
      chart.options.plugins.legend.display = true;
      chart.options.layout.padding = 20;
    }
    chart.update('none');
  }
}
```

## Common Configuration Patterns

### Dark Mode Support

```javascript
// Light theme (default)
const lightTheme = {
  color: '#666',
  borderColor: 'rgba(0, 0, 0, 0.1)',
  gridColor: 'rgba(0, 0, 0, 0.05)',
  backgroundColor: '#fff'
};

// Dark theme
const darkTheme = {
  color: '#ccc',
  borderColor: 'rgba(255, 255, 255, 0.1)',
  gridColor: 'rgba(255, 255, 255, 0.05)',
  backgroundColor: '#1a1a1a'
};

// Apply theme
function applyTheme(theme) {
  Chart.defaults.color = theme.color;
  Chart.defaults.borderColor = theme.borderColor;
  Chart.defaults.scales.linear.grid.color = theme.gridColor;
  Chart.defaults.plugins.tooltip.backgroundColor = theme.backgroundColor;
}

// Toggle theme
const isDark = window.matchMedia('(prefers-color-scheme: dark)').matches;
applyTheme(isDark ? darkTheme : lightTheme);
```

### Accessible Color Palettes

```javascript
// High contrast colors for accessibility
const accessibleColors = [
  '#1f77b4',  // Blue
  '#ff7f0e',  // Orange
  '#2ca02c',  // Green
  '#d62728',  // Red
  '#9467bd',  // Purple
  '#8c564b',  // Brown
  '#e377c2',  // Pink
  '#7f7f7f',  // Gray
];

// Apply to datasets
data: {
  datasets: datasets.map((ds, i) => ({
    ...ds,
    backgroundColor: accessibleColors[i % accessibleColors.length],
    borderColor: accessibleColors[i % accessibleColors.length]
  }))
}
```

### Performance Optimization

```javascript
options: {
  // Disable animations for better performance
  animation: false,

  // Reduce number of ticks
  scales: {
    x: {
      ticks: {
        maxTicksLimit: 10
      }
    },
    y: {
      ticks: {
        maxTicksLimit: 8
      }
    }
  },

  // Use decimation plugin for large datasets
  plugins: {
    decimation: {
      enabled: true,
      algorithm: 'lttb',  // Largest-Triangle-Three-Buckets
      samples: 500
    }
  }
}
```

### Print-Friendly Configuration

```javascript
// Detect print media
const printMediaQuery = window.matchMedia('print');

printMediaQuery.addEventListener('change', (e) => {
  if (e.matches) {
    // Printing: use high contrast, remove animations
    chart.options.animation = false;
    chart.options.plugins.legend.labels.color = '#000';
    chart.update('none');
  } else {
    // Screen: restore original settings
    chart.options.animation = true;
    chart.options.plugins.legend.labels.color = '#666';
    chart.update('none');
  }
});
```

### Locale Configuration

```javascript
// Set locale for number/date formatting
Chart.defaults.locale = 'de-DE';

// Custom number formatting
options: {
  scales: {
    y: {
      ticks: {
        callback: function(value) {
          return new Intl.NumberFormat('en-US', {
            style: 'currency',
            currency: 'USD'
          }).format(value);
        }
      }
    }
  }
}
```

## Configuration Reuse

### Creating Chart Templates

```javascript
// Base configuration for all company charts
const companyChartDefaults = {
  options: {
    responsive: true,
    maintainAspectRatio: true,
    plugins: {
      legend: {
        position: 'bottom',
        labels: {
          font: { family: "'Inter', sans-serif" },
          padding: 15
        }
      },
      tooltip: {
        backgroundColor: 'rgba(0, 0, 0, 0.85)',
        padding: 12,
        cornerRadius: 6
      }
    }
  }
};

// Merge with specific chart config
const myChart = new Chart(ctx, {
  type: 'line',
  data: myData,
  options: {
    ...companyChartDefaults.options,
    scales: {
      y: { beginAtZero: true }
    }
  }
});
```

### Configuration Factory Functions

```javascript
function createBarChartConfig(data, title) {
  return {
    type: 'bar',
    data: data,
    options: {
      responsive: true,
      plugins: {
        title: {
          display: true,
          text: title
        },
        legend: {
          display: false
        }
      },
      scales: {
        y: { beginAtZero: true }
      }
    }
  };
}

// Use factory
const chart1 = new Chart(ctx1, createBarChartConfig(salesData, 'Q1 Sales'));
const chart2 = new Chart(ctx2, createBarChartConfig(revenueData, 'Q1 Revenue'));
```
