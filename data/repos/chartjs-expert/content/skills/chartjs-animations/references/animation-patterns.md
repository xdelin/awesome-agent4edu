# Animation Patterns Reference

Advanced animation patterns for Chart.js v4.5.1. These patterns go beyond basic configuration to create custom visual effects.

## Progressive Line Reveal

Draw a line chart progressively from left to right, creating a "drawing" effect.

### Basic Pattern

The key is using scriptable `delay` and `from` properties to stagger when each point appears:

```javascript
const totalDuration = 10000;
const delayBetweenPoints = totalDuration / data.length;

const animation = {
  x: {
    type: 'number',
    easing: 'linear',
    duration: delayBetweenPoints,
    from: NaN, // Point initially skipped (not drawn)
    delay(ctx) {
      if (ctx.type !== 'data' || ctx.xStarted) {
        return 0;
      }
      ctx.xStarted = true;
      return ctx.index * delayBetweenPoints;
    }
  },
  y: {
    type: 'number',
    easing: 'linear',
    duration: delayBetweenPoints,
    from: (ctx) => {
      // Start from previous point's Y position for smooth line
      if (ctx.index === 0) {
        return ctx.chart.scales.y.getPixelForValue(100);
      }
      return ctx.chart.getDatasetMeta(ctx.datasetIndex)
        .data[ctx.index - 1].getProps(['y'], true).y;
    },
    delay(ctx) {
      if (ctx.type !== 'data' || ctx.yStarted) {
        return 0;
      }
      ctx.yStarted = true;
      return ctx.index * delayBetweenPoints;
    }
  }
};
```

### Key Concepts

| Concept | Purpose |
|---------|---------|
| `from: NaN` | Makes point invisible until animation starts |
| `ctx.xStarted` flag | Prevents re-triggering delay on updates |
| `getDatasetMeta().data[].getProps()` | Gets previous point's animated position |
| Matching X and Y delays | Creates smooth line drawing effect |

### With Easing Variation

Apply easing to the delay calculation for non-linear reveal speed:

```javascript
const easing = helpers.easingEffects.easeOutQuad;
const totalDuration = 5000;

const delay = (ctx) => easing(ctx.index / data.length) * totalDuration;
const duration = (ctx) => easing(ctx.index / data.length) * totalDuration / data.length;
```

This makes the line slow down as it approaches the end (easeOut) or speed up (easeIn).

## Staggered Bar Animation

Animate bars appearing one at a time with a cascading effect.

### Basic Pattern

```javascript
let delayed;

const config = {
  type: 'bar',
  data: data,
  options: {
    animation: {
      onComplete: () => {
        delayed = true;
      },
      delay: (context) => {
        let delay = 0;
        if (context.type === 'data' && context.mode === 'default' && !delayed) {
          // Delay by data point position and dataset position
          delay = context.dataIndex * 300 + context.datasetIndex * 100;
        }
        return delay;
      }
    }
  }
};
```

### Delay Calculation Strategies

| Strategy | Formula | Effect |
|----------|---------|--------|
| By data point | `dataIndex * 100` | Left to right wave |
| By dataset | `datasetIndex * 200` | Front to back layers |
| Combined | `dataIndex * 100 + datasetIndex * 50` | Diagonal cascade |
| Reverse | `(dataCount - dataIndex) * 100` | Right to left |

### Preventing Re-animation

The `delayed` flag prevents staggered animation on subsequent updates:

```javascript
let delayed = false;

options: {
  animation: {
    onComplete: () => { delayed = true; },
    delay: (ctx) => {
      if (ctx.mode === 'default' && !delayed) {
        return ctx.dataIndex * 100;
      }
      return 0;
    }
  }
}
```

## Drop Animation

Make data points drop from the top of the chart and bounce into position.

### Basic Pattern

```javascript
options: {
  animations: {
    y: {
      easing: 'easeOutBounce',
      from: (ctx) => {
        if (ctx.type === 'data') {
          if (ctx.mode === 'default' && !ctx.dropped) {
            ctx.dropped = true;
            return 0; // Start from top (pixel 0)
          }
        }
      }
    }
  }
}
```

### Drop From Scale Maximum

Drop from the maximum value on the Y scale:

```javascript
animations: {
  y: {
    duration: 2000,
    from: (ctx) => {
      if (ctx.type === 'data' && ctx.mode === 'default') {
        return ctx.chart.scales.y.getPixelForValue(
          ctx.chart.scales.y.max
        );
      }
    },
    easing: 'easeOutBounce'
  }
}
```

### Easing Recommendations for Drop

| Easing | Effect | Use Case |
|--------|--------|----------|
| `easeOutBounce` | Bounces at bottom | Playful, attention-grabbing |
| `easeOutElastic` | Overshoots and springs back | Energetic, dynamic |
| `easeOutQuart` | Smooth deceleration | Professional, subtle |
| `easeOutBack` | Slight overshoot | Elegant entrance |

## Color Fade Transition

Animate colors from one value to another.

### Fade In From Transparent

```javascript
options: {
  animations: {
    backgroundColor: {
      type: 'color',
      duration: 1000,
      from: 'transparent',
      easing: 'easeInOutQuad'
    },
    borderColor: {
      type: 'color',
      duration: 1000,
      from: 'transparent',
      easing: 'easeInOutQuad'
    }
  }
}
```

### Color Pulse Effect

Continuously animate between two colors:

```javascript
options: {
  animations: {
    backgroundColor: {
      type: 'color',
      duration: 2000,
      from: 'rgba(255, 99, 132, 0.2)',
      to: 'rgba(255, 99, 132, 0.8)',
      loop: true,
      easing: 'easeInOutSine'
    }
  }
}
```

### Dynamic Color Based on Value

```javascript
datasets: [{
  backgroundColor: (ctx) => {
    const value = ctx.dataset.data[ctx.dataIndex];
    return value > 0 ? 'rgba(75, 192, 192, 0.5)' : 'rgba(255, 99, 132, 0.5)';
  },
  animation: {
    backgroundColor: {
      type: 'color',
      duration: 500
    }
  }
}]
```

## Custom Drawing During Animation

Use the `onProgress` callback with a custom plugin for advanced visual effects.

### Reveal Curtain Effect

Draw a white rectangle that shrinks as animation progresses:

```javascript
let animationProgress = 0;

const chart = new Chart(ctx, {
  type: 'line',
  data: data,
  options: {
    animation: {
      duration: 2000,
      onProgress: (animation) => {
        animationProgress = animation.currentStep / animation.numSteps;
      }
    }
  },
  plugins: [{
    id: 'revealCurtain',
    afterDraw: (chart) => {
      if (animationProgress < 1) {
        const { ctx, chartArea: { left, top, right, bottom } } = chart;
        const width = right - left;
        const revealX = left + (width * animationProgress);

        ctx.save();
        ctx.fillStyle = 'rgba(255, 255, 255, 0.9)';
        ctx.fillRect(revealX, top, right - revealX, bottom - top);
        ctx.restore();
      }
    }
  }]
});
```

### Progress Indicator Overlay

Show a circular progress indicator during animation:

```javascript
plugins: [{
  id: 'progressIndicator',
  afterDraw: (chart, args, options) => {
    if (animationProgress < 1) {
      const { ctx, width, height } = chart;
      const centerX = width / 2;
      const centerY = height / 2;
      const radius = 30;

      // Draw arc based on progress
      ctx.save();
      ctx.beginPath();
      ctx.arc(centerX, centerY, radius, -Math.PI / 2,
        -Math.PI / 2 + (2 * Math.PI * animationProgress));
      ctx.strokeStyle = 'rgba(75, 192, 192, 0.8)';
      ctx.lineWidth = 4;
      ctx.stroke();
      ctx.restore();
    }
  }
}]
```

## Combining Patterns

### Staggered Drop with Color Fade

```javascript
options: {
  animation: {
    delay: (ctx) => ctx.dataIndex * 100
  },
  animations: {
    y: {
      from: 0,
      duration: 1500,
      easing: 'easeOutBounce'
    },
    backgroundColor: {
      type: 'color',
      duration: 1000,
      delay: 500, // Start after drop begins
      from: 'transparent'
    }
  }
}
```

### Per-Dataset Staggered Animation

```javascript
datasets: [
  {
    label: 'Revenue',
    data: revenueData,
    animation: { delay: 0 }
  },
  {
    label: 'Expenses',
    data: expenseData,
    animation: { delay: 500 }
  },
  {
    label: 'Profit',
    data: profitData,
    animation: { delay: 1000 }
  }
]
```

## Performance Tips

| Pattern | Performance Impact | Optimization |
|---------|-------------------|--------------|
| Progressive line | High for large datasets | Reduce `totalDuration` |
| Staggered | Medium | Use smaller delays |
| Drop | Low | N/A |
| Color fade | Low | N/A |
| Custom drawing | Varies | Minimize canvas operations |

### Disable for Large Datasets

```javascript
const dataPoints = data.datasets[0].data.length;

options: {
  animation: dataPoints > 500 ? false : {
    duration: Math.max(500, 2000 - dataPoints * 3)
  }
}
```

## Accessibility Considerations

Always respect user preferences for reduced motion:

```javascript
const prefersReducedMotion = window.matchMedia(
  '(prefers-reduced-motion: reduce)'
).matches;

options: {
  animation: prefersReducedMotion ? false : {
    // Your animation config
  }
}
```
