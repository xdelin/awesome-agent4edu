# Advanced Axis Customization

Advanced patterns for dynamic axis styling, scriptable options, and precise positioning in Chart.js.

## Scriptable Options

Many axis options accept functions instead of static values, enabling dynamic styling based on data or context.

### Scriptable Context Object

When an option is scriptable, your function receives a context object:

```javascript
{
  tick: {           // Current tick
    value: 0,       // Numeric value
    label: '0',     // Formatted label string
    major: false    // Whether this is a major tick
  },
  index: 0,         // Tick index in array
  scale: {...},     // The scale instance
  chart: {...}      // The chart instance
}
```

### Scriptable Grid Options

| Option | Scriptable | Indexable | Use Case |
|--------|:----------:|:---------:|----------|
| `grid.color` | Yes | Yes | Highlight specific lines |
| `grid.lineWidth` | Yes | Yes | Emphasize zero line |
| `grid.tickColor` | Yes | Yes | Match tick to grid |
| `grid.tickBorderDash` | Yes | Yes | Dashed tick marks |

```javascript
// Highlight zero line
grid: {
  color: (context) => context.tick.value === 0
    ? 'rgba(0, 0, 0, 0.5)'
    : 'rgba(0, 0, 0, 0.1)',
  lineWidth: (context) => context.tick.value === 0 ? 2 : 1
}

// Alternating grid colors
grid: {
  color: (context) => context.index % 2 === 0
    ? 'rgba(0, 0, 0, 0.1)'
    : 'rgba(0, 0, 0, 0.05)'
}
```

### Scriptable Tick Options

```javascript
ticks: {
  // Conditional tick colors
  color: (context) => {
    if (context.tick.value < 0) return 'red';
    if (context.tick.value > 100) return 'green';
    return '#666';
  },

  // Dynamic font weight
  font: (context) => ({
    size: 12,
    weight: context.tick.major ? 'bold' : 'normal'
  })
}
```

## Tick Alignment

Control label positioning with `align` and `crossAlign`:

| Property | Values | Direction |
|----------|--------|-----------|
| `align` | `start`, `center`, `end` | Along axis direction |
| `crossAlign` | `near`, `center`, `far` | Perpendicular to axis |

```javascript
// Left-align Y-axis labels (useful for horizontal bar charts)
scales: {
  y: {
    ticks: {
      align: 'start',      // Label position along tick
      crossAlign: 'far'    // Push labels away from chart
    }
  }
}
```

**Note:** `crossAlign` only works when tick rotation is 0 and axis is at edge position.

## Scale Bounds

The `bounds` property controls how the scale range is determined:

| Value | Behavior |
|-------|----------|
| `'data'` | Scale fits to data; labels outside data range are hidden |
| `'ticks'` | Scale fits to tick marks; data outside may be truncated |

```javascript
scales: {
  x: {
    bounds: 'data',    // Default for time scale
    // vs
    bounds: 'ticks'    // Default for linear scale
  }
}
```

**Tip:** When switching from time scale to linear/logarithmic, add `bounds: 'ticks'` to maintain expected behavior.

## Advanced Axis Positioning

### Center Position

Place an axis in the center of the chart area:

```javascript
scales: {
  x: {
    position: 'center',  // Axis at vertical center
    // Must specify axis type when using center
    axis: 'x'
  },
  y: {
    position: 'center',  // Axis at horizontal center
    axis: 'y'
  }
}
```

### Dynamic Position (Relative to Data)

Position an axis at a specific data value:

```javascript
scales: {
  y: {
    position: { x: 0 }   // Y-axis crosses X at x=0
  },
  x: {
    position: { y: 50 }  // X-axis crosses Y at y=50
  }
}
```

## Major Ticks

Enable major ticks to distinguish significant values:

```javascript
scales: {
  y: {
    ticks: {
      major: {
        enabled: true    // Mark ticks at round numbers as major
      },
      // Style major ticks differently
      font: (context) => ({
        weight: context.tick.major ? 'bold' : 'normal'
      }),
      color: (context) => context.tick.major ? '#000' : '#666'
    }
  }
}
```

## Practical Patterns

### Financial Chart Axes

```javascript
scales: {
  y: {
    ticks: {
      // Red for negative, green for positive
      color: (ctx) => ctx.tick.value < 0 ? '#dc3545' : '#28a745',
      callback: (value) => '$' + value.toLocaleString()
    },
    grid: {
      // Bold zero line
      color: (ctx) => ctx.tick.value === 0
        ? 'rgba(0,0,0,0.3)'
        : 'rgba(0,0,0,0.1)',
      lineWidth: (ctx) => ctx.tick.value === 0 ? 2 : 1
    }
  }
}
```

### Dashboard with Threshold Lines

```javascript
const THRESHOLD = 75;

scales: {
  y: {
    grid: {
      color: (ctx) => {
        if (ctx.tick.value === THRESHOLD) return 'rgba(255,0,0,0.5)';
        return 'rgba(0,0,0,0.1)';
      },
      lineWidth: (ctx) => ctx.tick.value === THRESHOLD ? 2 : 1
    }
  }
}
```
