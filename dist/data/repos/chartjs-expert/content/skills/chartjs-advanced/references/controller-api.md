# Chart Controller API Reference

Complete API documentation for extending Chart.js controllers.

## DatasetController Interface

All custom chart types must extend `Chart.DatasetController` and implement this interface.

### Required Properties

| Property | Type | Description |
|----------|------|-------------|
| `id` | `string` | Unique identifier for the controller (used as `type` in chart config) |
| `defaults` | `object` | Default configuration options |

### Required Defaults

```javascript
{
  // Type of element to create for the dataset (e.g., 'line', 'bar')
  // Set to false or null if no dataset-level element needed
  datasetElementType: string | null | false,

  // Type of element to create for each data point (e.g., 'point', 'bar')
  // Set to false or null if no per-data elements needed
  dataElementType: string | null | false,
}
```

### Required Methods

#### update(mode)

Update elements in response to new data.

```javascript
update(mode) {
  const meta = this.getMeta();
  // mode can be: 'active', 'hide', 'reset', 'resize', 'show', or undefined
  this.updateElements(meta.data, 0, meta.data.length, mode);
}
```

### Optional Override Methods

#### draw()

Draw the dataset representation. Base implementation works for most cases.

```javascript
draw() {
  super.draw();

  // Custom drawing after base implementation
  const meta = this.getMeta();
  const ctx = this.chart.ctx;

  meta.data.forEach((element, index) => {
    // Custom drawing logic
  });
}
```

#### initialize()

Initialize the controller. Called once when chart is created.

```javascript
initialize() {
  super.initialize();
  // Custom initialization
}
```

#### linkScales()

Link the dataset to scales. Override for charts using single scale (polar, doughnut).

```javascript
linkScales() {
  // Default links to x and y scales
  super.linkScales();
}
```

#### parse(start, count)

Parse raw data into controller metadata. Override for non-Cartesian data formats.

```javascript
parse(start, count) {
  // Default Cartesian parsing
  super.parse(start, count);
}
```

## Built-in Controllers

Available controllers to extend:

| Controller | Chart Type | Elements |
|------------|-----------|----------|
| `BarController` | Bar/Column | BarElement |
| `BubbleController` | Bubble | PointElement |
| `DoughnutController` | Doughnut/Pie | ArcElement |
| `LineController` | Line | LineElement + PointElement |
| `PieController` | Pie | ArcElement |
| `PolarAreaController` | Polar Area | ArcElement |
| `RadarController` | Radar | LineElement + PointElement |
| `ScatterController` | Scatter | PointElement |

## Accessing Controller Data

### getMeta()

Returns metadata object for the dataset:

```javascript
const meta = this.getMeta();
// meta.data - array of visual elements
// meta.dataset - the dataset-level element (e.g., the line in a line chart)
// meta.hidden - whether dataset is hidden
// meta.index - dataset index
```

### Element Properties

Get element properties with `getProps()`:

```javascript
const element = meta.data[0];
const {x, y, base, width, height} = element.getProps(['x', 'y', 'base', 'width', 'height']);
```

### Chart Context

Access chart context through `this.chart`:

```javascript
const ctx = this.chart.ctx;        // Canvas 2D context
const chartArea = this.chart.chartArea;  // Drawing area bounds
const data = this.chart.data;      // Chart data configuration
```

## Registration

Always register custom controllers:

```javascript
// ES modules
import {Chart} from 'chart.js';
Chart.register(MyController);

// UMD
Chart.register(Chart.MyController);
```

## TypeScript Support

Add type declarations for custom chart types using declaration merging.

### Extending ChartTypeRegistry

Create a `.d.ts` file to register your custom chart type:

```typescript
// custom-charts.d.ts
import {ChartTypeRegistry} from 'chart.js';

declare module 'chart.js' {
  interface ChartTypeRegistry {
    // Extend existing type (inherits all config options)
    derivedBubble: ChartTypeRegistry['bubble'];

    // Or define custom type with specific options
    customBar: {
      chartOptions: ChartTypeRegistry['bar']['chartOptions'];
      datasetOptions: ChartTypeRegistry['bar']['datasetOptions'] & {
        customOption?: boolean;
      };
      defaultDataPoint: ChartTypeRegistry['bar']['defaultDataPoint'];
      metaExtensions: ChartTypeRegistry['bar']['metaExtensions'];
      parsedDataType: ChartTypeRegistry['bar']['parsedDataType'];
      scales: ChartTypeRegistry['bar']['scales'];
    };
  }
}
```

### Usage with TypeScript

```typescript
import {Chart, ChartConfiguration} from 'chart.js';
import './custom-charts.d.ts';

// TypeScript now recognizes 'derivedBubble' as valid type
const config: ChartConfiguration<'derivedBubble'> = {
  type: 'derivedBubble',
  data: {
    datasets: [{
      data: [{x: 10, y: 20, r: 5}]
    }]
  }
};

const chart = new Chart(ctx, config);
```

### Type-Safe Custom Options

For controllers with custom options:

```typescript
// Define custom options interface
interface CustomBarOptions {
  customOption?: boolean;
  boxStrokeStyle?: string;
}

// Access options in controller
class CustomBarController extends BarController {
  draw() {
    super.draw();

    // TypeScript knows about customOption
    const opts = this.options as CustomBarOptions;
    if (opts.customOption) {
      // Custom drawing logic
    }
  }
}
```

## Common Patterns

### Adding Custom Drawing

```javascript
draw() {
  super.draw();

  const ctx = this.chart.ctx;
  ctx.save();

  // Custom drawing
  ctx.fillStyle = 'rgba(255, 0, 0, 0.5)';
  // ... drawing operations

  ctx.restore();
}
```

### Custom Hover Effects

```javascript
// In defaults
defaults: {
  ...ParentController.defaults,
  hoverBorderWidth: 3,
  hoverBorderColor: 'blue'
}
```

### Conditional Element Styling

```javascript
updateElements(elements, start, count, mode) {
  for (let i = start; i < start + count; i++) {
    const element = elements[i];
    const properties = this.resolveDataElementOptions(i, mode);

    // Custom property logic
    if (this.chart.data.datasets[this.index].data[i] < 0) {
      properties.backgroundColor = 'red';
    }

    this.updateElement(element, i, properties, mode);
  }
}
```
