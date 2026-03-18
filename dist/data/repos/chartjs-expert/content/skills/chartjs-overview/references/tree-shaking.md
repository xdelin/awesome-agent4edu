# Chart.js Tree-Shaking Guide

Complete reference for registering only the components needed for each chart type.

## Component Categories

Chart.js components are organized into these categories:

- **Controllers**: Chart type controllers (BarController, LineController, etc.)
- **Elements**: Visual elements (BarElement, PointElement, ArcElement, LineElement)
- **Scales**: Axis types (CategoryScale, LinearScale, TimeScale, etc.)
- **Plugins**: Built-in plugins (Legend, Tooltip, Title, etc.)

## Chart Type Requirements

### Bar Chart

```javascript
import {
  Chart,
  BarController,
  BarElement,
  CategoryScale,
  LinearScale,
  Legend,
  Tooltip
} from 'chart.js';

Chart.register(
  BarController,
  BarElement,
  CategoryScale,    // x-axis
  LinearScale,      // y-axis
  Legend,
  Tooltip
);
```

### Line Chart

```javascript
import {
  Chart,
  LineController,
  LineElement,
  PointElement,
  CategoryScale,
  LinearScale,
  Legend,
  Tooltip,
  Filler          // Only if using fill: true
} from 'chart.js';

Chart.register(
  LineController,
  LineElement,
  PointElement,
  CategoryScale,
  LinearScale,
  Legend,
  Tooltip,
  Filler
);
```

### Pie Chart

```javascript
import {
  Chart,
  PieController,
  ArcElement,
  Legend,
  Tooltip
} from 'chart.js';

Chart.register(
  PieController,
  ArcElement,
  Legend,
  Tooltip
);
```

### Doughnut Chart

```javascript
import {
  Chart,
  DoughnutController,
  ArcElement,
  Legend,
  Tooltip
} from 'chart.js';

Chart.register(
  DoughnutController,
  ArcElement,
  Legend,
  Tooltip
);
```

### Radar Chart

```javascript
import {
  Chart,
  RadarController,
  RadialLinearScale,
  PointElement,
  LineElement,
  Legend,
  Tooltip,
  Filler
} from 'chart.js';

Chart.register(
  RadarController,
  RadialLinearScale,
  PointElement,
  LineElement,
  Legend,
  Tooltip,
  Filler
);
```

### Polar Area Chart

```javascript
import {
  Chart,
  PolarAreaController,
  RadialLinearScale,
  ArcElement,
  Legend,
  Tooltip
} from 'chart.js';

Chart.register(
  PolarAreaController,
  RadialLinearScale,
  ArcElement,
  Legend,
  Tooltip
);
```

### Scatter Chart

```javascript
import {
  Chart,
  ScatterController,
  PointElement,
  LinearScale,
  Legend,
  Tooltip
} from 'chart.js';

Chart.register(
  ScatterController,
  PointElement,
  LinearScale,    // Both x and y
  Legend,
  Tooltip
);
```

### Bubble Chart

```javascript
import {
  Chart,
  BubbleController,
  PointElement,
  LinearScale,
  Legend,
  Tooltip
} from 'chart.js';

Chart.register(
  BubbleController,
  PointElement,
  LinearScale,
  Legend,
  Tooltip
);
```

## Available Scales

### Cartesian Scales (x/y axes)

| Scale | Import | Use Case |
|-------|--------|----------|
| `CategoryScale` | `import { CategoryScale } from 'chart.js'` | String labels |
| `LinearScale` | `import { LinearScale } from 'chart.js'` | Numeric data |
| `LogarithmicScale` | `import { LogarithmicScale } from 'chart.js'` | Exponential data |
| `TimeScale` | `import { TimeScale } from 'chart.js'` | Date/time data |
| `TimeSeriesScale` | `import { TimeSeriesScale } from 'chart.js'` | Time series |

### Radial Scales

| Scale | Import | Use Case |
|-------|--------|----------|
| `RadialLinearScale` | `import { RadialLinearScale } from 'chart.js'` | Radar, polar area |

## Available Plugins

| Plugin | Import | Purpose |
|--------|--------|---------|
| `Legend` | `import { Legend } from 'chart.js'` | Chart legend |
| `Tooltip` | `import { Tooltip } from 'chart.js'` | Hover tooltips |
| `Title` | `import { Title } from 'chart.js'` | Chart title |
| `SubTitle` | `import { SubTitle } from 'chart.js'` | Chart subtitle |
| `Filler` | `import { Filler } from 'chart.js'` | Area fill under lines |
| `Decimation` | `import { Decimation } from 'chart.js'` | Data decimation |
| `Colors` | `import { Colors } from 'chart.js'` | Default color palette |

## Time Scale Requirements

When using TimeScale or TimeSeriesScale, install a date adapter:

```bash
npm install chartjs-adapter-date-fns date-fns
# or
npm install chartjs-adapter-moment moment
# or
npm install chartjs-adapter-luxon luxon
# or
npm install chartjs-adapter-dayjs-4 dayjs
```

Import the adapter after Chart.js:

```javascript
import {
  Chart,
  TimeScale,
  // ... other imports
} from 'chart.js';
import 'chartjs-adapter-date-fns';

Chart.register(TimeScale, /* ... */);
```

## Mixed Charts

For mixed charts, register all required controllers:

```javascript
import {
  Chart,
  BarController,
  LineController,
  BarElement,
  LineElement,
  PointElement,
  CategoryScale,
  LinearScale,
  Legend,
  Tooltip
} from 'chart.js';

Chart.register(
  BarController,
  LineController,
  BarElement,
  LineElement,
  PointElement,
  CategoryScale,
  LinearScale,
  Legend,
  Tooltip
);
```

## Checking Missing Components

If you forget to register a component, Chart.js shows a helpful error:

```
Unhandled Promise Rejection: Error: "bar" is not a registered controller.
```

This indicates you need to import and register `BarController`.

## Bundle Size Comparison

Approximate sizes (minified + gzipped):

| Import Method | Size |
|--------------|------|
| `chart.js/auto` | ~67 KB |
| Manual (bar chart only) | ~35 KB |
| Manual (line chart only) | ~38 KB |
| Manual (pie chart only) | ~28 KB |

Tree-shaking can reduce bundle size by 25-60% depending on chart types used.
