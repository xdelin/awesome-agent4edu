# Performance Optimization Reference

Comprehensive guide to optimizing Chart.js performance for large datasets and performance-sensitive applications.

## Data Structure Optimization

### Disable Parsing

Provide data in the internal format accepted by datasets and scales, then disable parsing:

```javascript
new Chart(ctx, {
  type: 'line',
  data: {
    datasets: [{
      data: [{x: 1, y: 10}, {x: 2, y: 20}]  // Pre-formatted
    }]
  },
  options: {
    parsing: false  // Skip parsing step
  }
});
```

### Normalize Data

Chart.js is fastest with unique, sorted, consistent indices across datasets:

```javascript
new Chart(ctx, {
  type: 'line',
  data: data,
  options: {
    normalized: true  // Tell Chart.js data is pre-normalized
  }
});
```

**Requirements for normalized data:**

- Unique x values (no duplicates)
- Sorted in ascending order
- Consistent across all datasets

## Data Decimation

The decimation plugin reduces points before rendering, improving both memory usage and render performance.

### Configuration

```javascript
options: {
  plugins: {
    decimation: {
      enabled: true,
      algorithm: 'lttb',  // or 'min-max'
      samples: 500,       // Target samples for LTTB
      threshold: 800      // Minimum points to trigger decimation
    }
  }
}
```

### Algorithms

| Algorithm | Description | Best For |
|-----------|-------------|----------|
| `lttb` | Largest Triangle Three Buckets - preserves visual shape | Trend visualization |
| `min-max` | Preserves peaks and valleys (up to 4 points/pixel) | Noisy signals, peak detection |

### Requirements

Decimation only works when ALL conditions are met:

1. Dataset `indexAxis` is `'x'`
2. Chart type is `line`
3. X axis is `linear` or `time` type
4. `parsing: false` is set
5. Dataset object is mutable
6. Point count exceeds threshold

**Note**: If any requirement is not met, decimation silently skips - the chart renders with all data points. Check the browser console for no errors but slow performance to diagnose missing requirements.

```javascript
// Complete decimation example
new Chart(ctx, {
  type: 'line',
  data: {
    datasets: [{
      data: largeDataset,  // Pre-parsed [{x, y}, ...]
      radius: 0            // Hide points for performance
    }]
  },
  options: {
    parsing: false,
    animation: false,
    plugins: {
      decimation: {
        enabled: true,
        algorithm: 'lttb',
        samples: 500
      }
    },
    scales: {
      x: { type: 'linear' }
    }
  }
});
```

## Animation Optimization

### Disable Animations

For large datasets, disable animations entirely:

```javascript
options: {
  animation: false
}
```

Line charts use Path2D caching when animations are disabled (where Path2D is available).

### Conditional Animation

```javascript
const dataPoints = data.datasets[0].data.length;

options: {
  animation: dataPoints > 1000 ? false : {
    duration: 1000,
    easing: 'easeOutQuart'
  }
}
```

## Line Chart Optimizations

### Disable Points

Hide point markers to reduce drawing operations:

```javascript
// Per dataset
datasets: [{
  pointRadius: 0,
  pointHoverRadius: 0
}]

// Or globally
options: {
  elements: {
    point: {
      radius: 0
    }
  }
}
```

### Disable Line Drawing

For scatter-like visualization with only points:

```javascript
datasets: [{
  showLine: false  // Only draw points
}]

// Or globally
options: {
  showLine: false
}
```

### Enable Span Gaps

Improves performance by disabling line segmentation:

```javascript
datasets: [{
  spanGaps: true  // Draw continuous line, ignore gaps
}]

// Or globally
options: {
  spanGaps: true
}
```

### Keep Default Line Settings

Automatic data decimation during draw works when these remain default:

- `tension: 0` (no Bezier curves)
- `stepped: false`
- `borderDash: []` (solid line)

Bezier curves are more expensive to render than straight lines.

## Scale Optimization

### Specify Min/Max

Avoid auto-calculation of scale limits:

```javascript
options: {
  scales: {
    x: {
      type: 'time',
      min: new Date('2024-01-01').valueOf(),
      max: new Date('2024-12-31').valueOf()
    },
    y: {
      type: 'linear',
      min: 0,
      max: 100
    }
  }
}
```

### Tick Sampling

Reduce label calculation overhead:

```javascript
options: {
  scales: {
    x: {
      ticks: {
        sampleSize: 10,     // Sample labels for sizing
        maxRotation: 0,     // Fixed rotation (avoid calculation)
        autoSkip: true
      }
    }
  }
}
```

## Web Workers

For heavy computation, offload to a web worker using OffscreenCanvas:

### Main Thread

```javascript
const canvas = document.getElementById('myChart');
const offscreen = canvas.transferControlToOffscreen();

const worker = new Worker('chart-worker.js');
worker.postMessage({canvas: offscreen, config}, [offscreen]);
```

### Worker (chart-worker.js)

```javascript
importScripts('https://cdn.jsdelivr.net/npm/chart.js');

onmessage = function(event) {
  const {canvas, config} = event.data;
  const chart = new Chart(canvas, config);

  // Manual resize (no DOM events in workers)
  canvas.width = 800;
  canvas.height = 400;
  chart.resize();
};
```

### Limitations

- No DOM access (no mouse interactions by default)
- Functions can't be transferred (strip from config)
- Manual resize handling required
- Check browser support for OffscreenCanvas

## Performance Checklist

| Optimization | Impact | When to Use |
|--------------|--------|-------------|
| `parsing: false` | High | Pre-formatted data |
| `normalized: true` | High | Sorted, unique x values |
| Decimation plugin | High | >1000 points |
| `animation: false` | Medium | Large datasets, frequent updates |
| `pointRadius: 0` | Medium | Line charts with many points |
| Specify `min`/`max` | Medium | Known data ranges |
| `spanGaps: true` | Low | Continuous data |
| Tick sampling | Low | Many axis labels |

## Streaming Data Pattern

Efficiently update charts with real-time data:

```javascript
function addDataPoint(chart, label, values) {
  chart.data.labels.push(label);
  chart.data.datasets.forEach((dataset, i) => {
    dataset.data.push(values[i]);
  });

  // Maintain fixed window (e.g., last 100 points)
  const maxPoints = 100;
  if (chart.data.labels.length > maxPoints) {
    chart.data.labels.shift();
    chart.data.datasets.forEach(dataset => {
      dataset.data.shift();
    });
  }

  // Use 'none' mode for no animation on updates
  chart.update('none');
}
```

## Babel Transpilation Note

If using Babel 7.9+, enable `loose` mode for better class construction performance. See [babel/babel#11356](https://github.com/babel/babel/issues/11356).
