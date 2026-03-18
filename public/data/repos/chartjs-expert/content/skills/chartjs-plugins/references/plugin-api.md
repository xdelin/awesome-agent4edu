# Chart.js Plugin API Reference

Complete reference for Chart.js plugin lifecycle hooks.

## Hook Execution Order

### Chart Initialization

```
beforeInit
  └─> afterInit
```

### Chart Update

```
beforeUpdate (cancelable)
  └─> beforeElementsUpdate
    └─> reset (resets all elements)
      └─> beforeDatasetsUpdate (cancelable)
        └─> beforeDatasetUpdate (per dataset, cancelable)
          └─> afterDatasetUpdate (per dataset)
        └─> afterDatasetsUpdate
      └─> beforeLayout (cancelable)
        └─> afterLayout
      └─> afterUpdate
```

### Chart Render

```
beforeRender (cancelable)
  └─> beforeDraw (cancelable)
    └─> [Draw chart elements]
      └─> afterDraw
    └─> afterRender
```

### Event Handling

```
beforeEvent (can modify event)
  └─> [Process event]
    └─> afterEvent (args.changed = true triggers re-render)
```

### Chart Destruction

```
beforeDestroy
  └─> stop (stops all animations)
    └─> afterDestroy
```

## Hook Parameters

All hooks receive three parameters:

```javascript
hookName: (chart, args, options) => {
  // chart: Chart instance
  // args: Hook-specific arguments
  // options: Plugin options from config
}
```

### Chart Parameter

The `chart` parameter provides access to:

```javascript
{
  ctx,              // CanvasRenderingContext2D
  canvas,           // HTMLCanvasElement
  width,            // number - chart width
  height,           // number - chart height
  aspectRatio,      // number
  currentDevicePixelRatio, // number
  chartArea,        // {left, top, right, bottom, width, height}
  scales,           // {[scaleId]: Scale}
  data,             // ChartData
  options,          // ChartOptions
  config,           // ChartConfiguration

  // Methods
  update: (mode?) => void,
  render: () => void,
  draw: () => void,
  stop: () => void,
  resize: (width?, height?) => void,
  clear: () => void,
  toBase64Image: (type?, quality?) => string,
  destroy: () => void,
  getElementsAtEventForMode: (e, mode, options, useFinalPosition) => Element[],
  setActiveElements: (elements) => void,
  getActiveElements: () => ActiveElement[]
}
```

### Args Parameter

Hook-specific arguments:

| Hook | Args Contents |
|------|---------------|
| `beforeEvent` / `afterEvent` | `{event, replay, changed, cancelable}` |
| `beforeDatasetUpdate` / `afterDatasetUpdate` | `{index, meta}` |
| Most render hooks | `{}` (empty object) |

### Options Parameter

Plugin-specific options from `options.plugins.{plugin-id}`:

```javascript
// Chart config
{
  options: {
    plugins: {
      myPlugin: {
        color: 'red',
        enabled: true
      }
    }
  }
}

// In plugin hook
beforeDraw: (chart, args, options) => {
  console.log(options.color);   // 'red'
  console.log(options.enabled); // true
}
```

## Cancelable Hooks

Some hooks can prevent subsequent operations by returning `false`:

### beforeUpdate

Return `false` to prevent chart update:

```javascript
beforeUpdate: (chart, args, options) => {
  if (someCondition) {
    return false;  // Cancel update
  }
}
```

### beforeRender

Return `false` to prevent rendering:

```javascript
beforeRender: (chart, args, options) => {
  if (!chart.data.datasets.length) {
    return false;  // Don't render empty chart
  }
}
```

### beforeDraw

Return `false` to prevent drawing:

```javascript
beforeDraw: (chart, args, options) => {
  if (options.enabled === false) {
    return false;  // Skip drawing
  }
}
```

### beforeLayout

Return `false` to prevent layout calculation:

```javascript
beforeLayout: (chart, args, options) => {
  if (chart.width < 100) {
    return false;  // Chart too small
  }
}
```

### beforeDatasetsUpdate

Return `false` to prevent dataset updates:

```javascript
beforeDatasetsUpdate: (chart, args, options) => {
  if (chart.data.datasets.length === 0) {
    return false;
  }
}
```

## Event Hook Details

### beforeEvent

Intercept and modify events before processing:

```javascript
beforeEvent: (chart, args, options) => {
  const {event} = args;

  // Modify event
  if (event.type === 'click') {
    event.x += 10;  // Offset click position
  }

  // Cancel event processing
  if (someCondition) {
    return false;
  }
}
```

### afterEvent

React to events and trigger re-renders:

```javascript
afterEvent: (chart, args, options) => {
  const {event, replay, changed} = args;

  if (event.type === 'mousemove') {
    // Do something
    args.changed = true;  // Trigger re-render
  }
}
```

## Canvas Context Methods

Common canvas operations in plugins:

### Drawing Shapes

```javascript
ctx.fillRect(x, y, width, height);
ctx.strokeRect(x, y, width, height);
ctx.clearRect(x, y, width, height);

ctx.beginPath();
ctx.arc(x, y, radius, startAngle, endAngle);
ctx.stroke();
ctx.fill();

ctx.moveTo(x, y);
ctx.lineTo(x, y);
ctx.stroke();
```

### Text Rendering

```javascript
ctx.font = '16px Arial';
ctx.fillStyle = 'black';
ctx.textAlign = 'center';      // 'left' | 'right' | 'center' | 'start' | 'end'
ctx.textBaseline = 'middle';   // 'top' | 'hanging' | 'middle' | 'alphabetic' | 'ideographic' | 'bottom'
ctx.fillText(text, x, y, maxWidth?);
ctx.strokeText(text, x, y, maxWidth?);

// Measure text
const metrics = ctx.measureText(text);
const width = metrics.width;
```

### Styling

```javascript
ctx.fillStyle = 'rgb(200, 0, 0)';
ctx.strokeStyle = 'rgba(0, 0, 200, 0.5)';
ctx.lineWidth = 2;
ctx.lineCap = 'round';     // 'butt' | 'round' | 'square'
ctx.lineJoin = 'round';    // 'round' | 'bevel' | 'miter'
ctx.setLineDash([5, 10]);  // [dash length, gap length]
ctx.lineDashOffset = 0;
```

### Transformations

```javascript
ctx.save();     // Save current state
ctx.restore();  // Restore saved state

ctx.translate(x, y);
ctx.rotate(angle);
ctx.scale(x, y);
ctx.transform(a, b, c, d, e, f);
```

### Compositing

```javascript
ctx.globalAlpha = 0.5;
ctx.globalCompositeOperation = 'source-over';
// Operations: source-over, destination-over, lighter, copy, etc.
```

## Plugin Defaults

Define default options for your plugin:

```javascript
const plugin = {
  id: 'myPlugin',
  beforeDraw: (chart, args, options) => {
    // options.color defaults to 'blue'
    // options.width defaults to 2
  },
  defaults: {
    color: 'blue',
    width: 2,
    nested: {
      enabled: true,
      style: 'solid'
    }
  }
};
```

Default merging follows deep merge strategy - user options override defaults recursively.

## Multiple Plugin Instances

Plugins execute in registration order. For global plugins:

```javascript
Chart.register(plugin1, plugin2, plugin3);

// Execution order for each hook:
// 1. plugin1.beforeDraw()
// 2. plugin2.beforeDraw()
// 3. plugin3.beforeDraw()
```

Inline plugins execute after global plugins:

```javascript
new Chart(ctx, {
  plugins: [inlinePlugin1, inlinePlugin2]
});

// Execution order:
// 1. All global plugins
// 2. inlinePlugin1
// 3. inlinePlugin2
```

## Performance Considerations

### Hook Frequency

| Hook Category | Frequency | Performance Impact |
|---------------|-----------|-------------------|
| Init hooks | Once per chart | Low |
| Update hooks | On data change | Medium |
| Render hooks | Every frame (animations) | **High** |
| Event hooks | On user interaction | Medium-High |

### Optimization Tips

1. **Cache in update hooks, use in render hooks:**

```javascript
let cachedData;

const plugin = {
  afterUpdate: (chart) => {
    // Expensive calculation once per update
    cachedData = expensiveCalculation(chart.data);
  },
  afterDraw: (chart) => {
    // Use cached result (called many times during animation)
    drawUsingCachedData(cachedData);
  }
};
```

2. **Early returns for disabled plugins:**

```javascript
beforeDraw: (chart, args, options) => {
  if (!options.enabled) return;
  // Expensive drawing logic
}
```

3. **Minimize canvas state changes:**

```javascript
// Good: Single save/restore
ctx.save();
ctx.fillStyle = 'red';
ctx.fillRect(0, 0, 10, 10);
ctx.fillRect(20, 20, 10, 10);
ctx.restore();

// Bad: Multiple save/restore cycles
ctx.save();
ctx.fillStyle = 'red';
ctx.fillRect(0, 0, 10, 10);
ctx.restore();
ctx.save();
ctx.fillStyle = 'red';
ctx.fillRect(20, 20, 10, 10);
ctx.restore();
```

## Debugging Plugins

```javascript
const debugPlugin = {
  id: 'debug',
  beforeInit: (chart) => console.log('beforeInit', chart.id),
  afterInit: (chart) => console.log('afterInit', chart.id),
  beforeUpdate: (chart) => console.log('beforeUpdate', chart.id),
  afterUpdate: (chart) => console.log('afterUpdate', chart.id),
  beforeDraw: (chart) => console.log('beforeDraw', chart.id),
  afterDraw: (chart) => console.log('afterDraw', chart.id),
  beforeEvent: (chart, args) => console.log('beforeEvent', args.event.type),
  afterEvent: (chart, args) => console.log('afterEvent', args.event.type)
};
```
