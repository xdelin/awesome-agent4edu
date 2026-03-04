# Custom Chart Types

Create custom chart types by extending `Chart.DatasetController`.

## Basic Structure

```javascript
class MyType extends Chart.DatasetController {
  // Required: unique identifier
  static id = 'myType';

  // Required: default configuration
  static defaults = {
    // Element type for the dataset (e.g., 'line', 'bar', null)
    datasetElementType: false,
    // Element type for each data point (e.g., 'point', 'bar', null)
    dataElementType: 'point'
  };

  // Required: update elements when data changes
  update(mode) {
    // mode: 'active', 'hide', 'reset', 'resize', 'show', or undefined
  }
}

Chart.register(MyType);

// Use the new chart type
new Chart(ctx, {
  type: 'myType',
  data: data,
  options: options
});
```

## Optional Methods

```javascript
class MyType extends Chart.DatasetController {
  // Draw dataset representation
  draw() {}

  // Initialize controller
  initialize() {}

  // Link to scales (override for single-scale charts like pie/doughnut)
  linkScales() {}

  // Parse data into controller metadata
  parse(start, count) {}
}
```

## Extending Built-in Types

Extend existing controllers for customization:

| Controller | Chart Type |
|------------|------------|
| `BarController` | bar |
| `BubbleController` | bubble |
| `DoughnutController` | doughnut |
| `LineController` | line |
| `PieController` | pie |
| `PolarAreaController` | polarArea |
| `RadarController` | radar |
| `ScatterController` | scatter |

### Example: Custom Bubble Chart

```javascript
import { BubbleController } from 'chart.js';

class CustomBubble extends BubbleController {
  draw() {
    // Call parent to draw all points
    super.draw(arguments);

    // Custom drawing: red box around first point
    const meta = this.getMeta();
    const pt0 = meta.data[0];
    const { x, y } = pt0.getProps(['x', 'y']);
    const { radius } = pt0.options;

    const ctx = this.chart.ctx;
    ctx.save();
    ctx.strokeStyle = 'red';
    ctx.lineWidth = 1;
    ctx.strokeRect(x - radius, y - radius, 2 * radius, 2 * radius);
    ctx.restore();
  }
}

CustomBubble.id = 'customBubble';
CustomBubble.defaults = BubbleController.defaults;

Chart.register(CustomBubble);

new Chart(ctx, {
  type: 'customBubble',
  data: data
});
```

## TypeScript Typings

Augment `ChartTypeRegistry` for new chart types:

```typescript
import { ChartTypeRegistry } from 'chart.js';

declare module 'chart.js' {
  interface ChartTypeRegistry {
    customBubble: ChartTypeRegistry['bubble']
  }
}
```

## UMD Access

In UMD builds, controllers are available directly on `Chart`:

```javascript
// Chart.BarController, Chart.LineController, etc.
class Custom extends Chart.BubbleController {
  // ...
}
```
