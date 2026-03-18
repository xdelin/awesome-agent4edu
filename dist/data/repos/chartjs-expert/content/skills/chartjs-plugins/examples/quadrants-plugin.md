# Quadrants Plugin Example

Complete working example of a plugin that draws colored quadrant backgrounds on a scatter chart.

## Plugin Implementation

```javascript
const quadrants = {
  id: 'quadrants',
  beforeDraw(chart, args, options) {
    const {ctx, chartArea: {left, top, right, bottom}, scales: {x, y}} = chart;

    // Find the pixel coordinates for value 0 on each axis
    const midX = x.getPixelForValue(0);
    const midY = y.getPixelForValue(0);

    ctx.save();

    // Top-left quadrant (negative X, positive Y)
    ctx.fillStyle = options.topLeft || 'rgba(255, 0, 0, 0.1)';
    ctx.fillRect(left, top, midX - left, midY - top);

    // Top-right quadrant (positive X, positive Y)
    ctx.fillStyle = options.topRight || 'rgba(0, 255, 0, 0.1)';
    ctx.fillRect(midX, top, right - midX, midY - top);

    // Bottom-right quadrant (positive X, negative Y)
    ctx.fillStyle = options.bottomRight || 'rgba(0, 0, 255, 0.1)';
    ctx.fillRect(midX, midY, right - midX, bottom - midY);

    // Bottom-left quadrant (negative X, negative Y)
    ctx.fillStyle = options.bottomLeft || 'rgba(255, 255, 0, 0.1)';
    ctx.fillRect(left, midY, midX - left, bottom - midY);

    ctx.restore();
  },
  defaults: {
    topLeft: 'rgba(255, 0, 0, 0.1)',
    topRight: 'rgba(0, 255, 0, 0.1)',
    bottomRight: 'rgba(0, 0, 255, 0.1)',
    bottomLeft: 'rgba(255, 255, 0, 0.1)'
  }
};
```

## Chart Configuration

```javascript
import Chart from 'chart.js/auto';

// Register the plugin globally
Chart.register(quadrants);

// Create scatter chart with quadrants
const data = {
  datasets: [{
    label: 'Dataset 1',
    data: [
      {x: -50, y: 30},
      {x: 40, y: 60},
      {x: 70, y: -20},
      {x: -30, y: -40}
    ],
    borderColor: 'rgb(75, 192, 192)',
    backgroundColor: 'rgba(75, 192, 192, 0.5)'
  }]
};

const config = {
  type: 'scatter',
  data: data,
  options: {
    plugins: {
      quadrants: {
        topLeft: 'rgba(255, 99, 132, 0.15)',
        topRight: 'rgba(54, 162, 235, 0.15)',
        bottomRight: 'rgba(75, 192, 192, 0.15)',
        bottomLeft: 'rgba(255, 205, 86, 0.15)'
      }
    },
    scales: {
      x: {
        min: -100,
        max: 100
      },
      y: {
        min: -100,
        max: 100
      }
    }
  }
};

const chart = new Chart(
  document.getElementById('myChart'),
  config
);
```

## Usage with React (react-chartjs-2)

```jsx
import { Scatter } from 'react-chartjs-2';
import {
  Chart as ChartJS,
  LinearScale,
  PointElement,
  Tooltip,
  Legend
} from 'chart.js';

// Register plugin
const quadrants = {
  id: 'quadrants',
  beforeDraw(chart, args, options) {
    // ... plugin code above
  }
};

ChartJS.register(
  LinearScale,
  PointElement,
  Tooltip,
  Legend,
  quadrants  // Register custom plugin
);

function QuadrantChart() {
  const data = {
    datasets: [{
      label: 'My Data',
      data: [
        {x: -50, y: 30},
        {x: 40, y: 60},
        {x: 70, y: -20},
        {x: -30, y: -40}
      ]
    }]
  };

  const options = {
    plugins: {
      quadrants: {
        topLeft: 'rgba(255, 99, 132, 0.15)',
        topRight: 'rgba(54, 162, 235, 0.15)',
        bottomRight: 'rgba(75, 192, 192, 0.15)',
        bottomLeft: 'rgba(255, 205, 86, 0.15)'
      }
    },
    scales: {
      x: { min: -100, max: 100 },
      y: { min: -100, max: 100 }
    }
  };

  return <Scatter data={data} options={options} />;
}
```

## Usage with Vue (vue-chartjs)

```vue
<template>
  <Scatter :data="chartData" :options="chartOptions" />
</template>

<script setup>
import { Scatter } from 'vue-chartjs';
import {
  Chart as ChartJS,
  LinearScale,
  PointElement,
  Tooltip,
  Legend
} from 'chart.js';

const quadrants = {
  id: 'quadrants',
  beforeDraw(chart, args, options) {
    // ... plugin code
  }
};

ChartJS.register(
  LinearScale,
  PointElement,
  Tooltip,
  Legend,
  quadrants
);

const chartData = {
  datasets: [{
    label: 'My Data',
    data: [
      {x: -50, y: 30},
      {x: 40, y: 60},
      {x: 70, y: -20},
      {x: -30, y: -40}
    ]
  }]
};

const chartOptions = {
  plugins: {
    quadrants: {
      topLeft: 'rgba(255, 99, 132, 0.15)',
      topRight: 'rgba(54, 162, 235, 0.15)',
      bottomRight: 'rgba(75, 192, 192, 0.15)',
      bottomLeft: 'rgba(255, 205, 86, 0.15)'
    }
  },
  scales: {
    x: { min: -100, max: 100 },
    y: { min: -100, max: 100 }
  }
};
</script>
```

## Disabling the Plugin

```javascript
// Disable for specific chart
const chart = new Chart(ctx, {
  type: 'scatter',
  data: data,
  options: {
    plugins: {
      quadrants: false  // Disable quadrants plugin
    }
  }
});
```

## Dynamic Updates

```javascript
// Update quadrant colors dynamically
function updateQuadrantColors(chart, colors) {
  chart.options.plugins.quadrants = colors;
  chart.update();
}

// Usage
updateQuadrantColors(chart, {
  topLeft: 'rgba(100, 0, 0, 0.2)',
  topRight: 'rgba(0, 100, 0, 0.2)',
  bottomRight: 'rgba(0, 0, 100, 0.2)',
  bottomLeft: 'rgba(100, 100, 0, 0.2)'
});
```

## TypeScript Support

```typescript
import {Chart, ChartType, Plugin} from 'chart.js';

interface QuadrantsOptions {
  topLeft?: string;
  topRight?: string;
  bottomRight?: string;
  bottomLeft?: string;
}

declare module 'chart.js' {
  interface PluginOptionsByType<TType extends ChartType> {
    quadrants?: QuadrantsOptions;
  }
}

const quadrants: Plugin<'scatter'> = {
  id: 'quadrants',
  beforeDraw(chart, args, options) {
    // TypeScript knows about options.topLeft, etc.
    const {ctx, chartArea: {left, top, right, bottom}, scales: {x, y}} = chart;
    const midX = x.getPixelForValue(0);
    const midY = y.getPixelForValue(0);

    ctx.save();
    ctx.fillStyle = options.topLeft || 'rgba(255, 0, 0, 0.1)';
    ctx.fillRect(left, top, midX - left, midY - top);
    // ... rest of implementation
    ctx.restore();
  }
};
```

## Key Concepts

1. **beforeDraw hook** - Draws before chart elements, so quadrants appear as background
2. **getPixelForValue(0)** - Converts data value (0) to canvas pixel coordinate
3. **ctx.save() / ctx.restore()** - Preserves canvas state for other drawings
4. **chartArea** - Boundaries of the actual chart plotting area (excludes axes, labels)
5. **Plugin options** - Configured via `options.plugins.quadrants`

## Use Cases

- **Risk matrices** - Visualize risk quadrants (likelihood vs impact)
- **Portfolio analysis** - Growth/value vs large/small cap
- **Performance grids** - Sales vs profit margins
- **SWOT analysis** - Strengths, weaknesses, opportunities, threats
- **Eisenhower matrix** - Urgent/important task prioritization
