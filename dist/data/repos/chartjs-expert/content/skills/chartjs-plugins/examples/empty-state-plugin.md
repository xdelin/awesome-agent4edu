# Empty State Plugin Example

Plugin that displays a message or visual indicator when a chart has no data.

## Plugin Implementation

```javascript
const emptyStatePlugin = {
  id: 'emptyState',
  afterDraw(chart, args, options) {
    const {datasets} = chart.data;

    // Check if any dataset has data
    let hasData = false;
    for (let i = 0; i < datasets.length; i++) {
      if (datasets[i].data && datasets[i].data.length > 0) {
        hasData = true;
        break;
      }
    }

    // If no data, draw empty state
    if (!hasData) {
      const {ctx, chartArea: {left, top, right, bottom}} = chart;
      const centerX = (left + right) / 2;
      const centerY = (top + bottom) / 2;

      ctx.save();

      // Draw message
      ctx.font = options.font || '16px Arial';
      ctx.fillStyle = options.color || 'rgba(0, 0, 0, 0.5)';
      ctx.textAlign = 'center';
      ctx.textBaseline = 'middle';
      ctx.fillText(options.message || 'No data available', centerX, centerY);

      ctx.restore();
    }
  },
  defaults: {
    message: 'No data available',
    font: '16px Arial',
    color: 'rgba(0, 0, 0, 0.5)'
  }
};
```

## Empty Doughnut Variant

For doughnut charts, draw a ring to show the expected chart shape:

```javascript
const emptyDoughnutPlugin = {
  id: 'emptyDoughnut',
  afterDraw(chart, args, options) {
    const {datasets} = chart.data;
    let hasData = datasets.some(ds => ds.data && ds.data.length > 0);

    if (!hasData) {
      const {ctx, chartArea: {left, top, right, bottom}} = chart;
      const centerX = (left + right) / 2;
      const centerY = (top + bottom) / 2;
      const radius = Math.min(right - left, bottom - top) / 2;

      ctx.save();

      // Draw placeholder ring
      ctx.beginPath();
      ctx.lineWidth = options.width || 2;
      ctx.strokeStyle = options.color || 'rgba(200, 200, 200, 0.5)';
      ctx.arc(
        centerX,
        centerY,
        radius - (options.radiusDecrease || 0),
        0,
        2 * Math.PI
      );
      ctx.stroke();

      // Optional: Add text in center
      if (options.message) {
        ctx.font = options.font || '14px Arial';
        ctx.fillStyle = options.textColor || 'rgba(0, 0, 0, 0.5)';
        ctx.textAlign = 'center';
        ctx.textBaseline = 'middle';
        ctx.fillText(options.message, centerX, centerY);
      }

      ctx.restore();
    }
  },
  defaults: {
    color: 'rgba(200, 200, 200, 0.5)',
    width: 2,
    radiusDecrease: 20,
    message: '',
    font: '14px Arial',
    textColor: 'rgba(0, 0, 0, 0.5)'
  }
};
```

## Usage

```javascript
import Chart from 'chart.js/auto';

Chart.register(emptyStatePlugin);

const data = {
  labels: ['Jan', 'Feb', 'Mar'],
  datasets: [{
    label: 'Sales',
    data: []  // Empty dataset
  }]
};

const chart = new Chart(ctx, {
  type: 'bar',
  data: data,
  options: {
    plugins: {
      emptyState: {
        message: 'No sales data yet',
        font: '18px Arial',
        color: 'gray'
      }
    }
  }
});
```

## React Example

```jsx
import { Bar } from 'react-chartjs-2';
import {
  Chart as ChartJS,
  CategoryScale,
  LinearScale,
  BarElement
} from 'chart.js';

const emptyStatePlugin = {
  id: 'emptyState',
  afterDraw(chart, args, options) {
    // ... plugin implementation
  }
};

ChartJS.register(
  CategoryScale,
  LinearScale,
  BarElement,
  emptyStatePlugin
);

function SalesChart({ salesData }) {
  const data = {
    labels: ['Jan', 'Feb', 'Mar', 'Apr', 'May'],
    datasets: [{
      label: 'Sales',
      data: salesData  // Could be empty array
    }]
  };

  const options = {
    plugins: {
      emptyState: {
        message: salesData.length === 0
          ? 'No sales data yet'
          : 'Loading...',
        font: '18px Arial',
        color: '#999'
      }
    }
  };

  return <Bar data={data} options={options} />;
}
```

## Advanced: Loading State

Extend to show loading indicator:

```javascript
const loadingStatePlugin = {
  id: 'loadingState',
  afterDraw(chart, args, options) {
    // Check for data
    const {datasets} = chart.data;
    const hasData = datasets.some(ds => ds.data && ds.data.length > 0);

    if (!hasData) {
      const {ctx, chartArea: {left, top, right, bottom}} = chart;
      const centerX = (left + right) / 2;
      const centerY = (top + bottom) / 2;

      ctx.save();

      if (options.loading) {
        // Draw spinner
        const radius = 20;
        const lineWidth = 3;
        const rotation = (Date.now() / 10) % 360;

        ctx.strokeStyle = options.spinnerColor || '#3498db';
        ctx.lineWidth = lineWidth;
        ctx.lineCap = 'round';

        ctx.translate(centerX, centerY);
        ctx.rotate((rotation * Math.PI) / 180);

        ctx.beginPath();
        ctx.arc(0, 0, radius, 0, Math.PI * 1.5);
        ctx.stroke();

        ctx.setTransform(1, 0, 0, 1, 0, 0);  // Reset transform

        // Animate
        requestAnimationFrame(() => chart.update('none'));
      } else {
        // Draw empty state message
        ctx.font = options.font || '16px Arial';
        ctx.fillStyle = options.color || 'gray';
        ctx.textAlign = 'center';
        ctx.textBaseline = 'middle';
        ctx.fillText(options.message || 'No data available', centerX, centerY);
      }

      ctx.restore();
    }
  },
  defaults: {
    loading: false,
    message: 'No data available',
    font: '16px Arial',
    color: 'gray',
    spinnerColor: '#3498db'
  }
};
```

### Usage with Loading State

```javascript
// Show loading spinner
chart.options.plugins.loadingState.loading = true;
chart.update('none');

// Fetch data
fetch('/api/data')
  .then(res => res.json())
  .then(data => {
    chart.data.datasets[0].data = data;
    chart.options.plugins.loadingState.loading = false;
    chart.update();
  });
```

## Empty State with Icon

```javascript
const emptyStateWithIcon = {
  id: 'emptyStateIcon',
  afterDraw(chart, args, options) {
    const {datasets} = chart.data;
    const hasData = datasets.some(ds => ds.data?.length > 0);

    if (!hasData) {
      const {ctx, chartArea: {left, top, right, bottom}} = chart;
      const centerX = (left + right) / 2;
      const centerY = (top + bottom) / 2;

      ctx.save();

      // Draw icon (simple folder/chart icon)
      ctx.strokeStyle = options.iconColor || 'rgba(0, 0, 0, 0.3)';
      ctx.fillStyle = options.iconColor || 'rgba(0, 0, 0, 0.3)';
      ctx.lineWidth = 2;

      // Draw simple bar chart icon
      const iconSize = 40;
      const barWidth = 8;
      const spacing = 4;

      ctx.fillRect(centerX - iconSize/2, centerY - 10, barWidth, 20);
      ctx.fillRect(centerX - iconSize/2 + barWidth + spacing, centerY - 20, barWidth, 30);
      ctx.fillRect(centerX - iconSize/2 + (barWidth + spacing) * 2, centerY - 15, barWidth, 25);

      // Draw message below icon
      ctx.font = options.font || '14px Arial';
      ctx.fillStyle = options.color || 'rgba(0, 0, 0, 0.5)';
      ctx.textAlign = 'center';
      ctx.textBaseline = 'top';
      ctx.fillText(options.message || 'No data', centerX, centerY + 25);

      ctx.restore();
    }
  }
};
```

## TypeScript Support

```typescript
import {Plugin, ChartType} from 'chart.js';

interface EmptyStateOptions {
  message?: string;
  font?: string;
  color?: string;
  loading?: boolean;
  spinnerColor?: string;
}

declare module 'chart.js' {
  interface PluginOptionsByType<TType extends ChartType> {
    emptyState?: EmptyStateOptions;
  }
}

const emptyStatePlugin: Plugin = {
  id: 'emptyState',
  afterDraw(chart, args, options) {
    // ... implementation
  }
};
```

## Use Cases

- **Dashboard widgets** - Show friendly message when no data loaded
- **Real-time charts** - Indicate waiting for first data point
- **User-generated content** - Prompt user to add data
- **Filtered views** - Explain when filters produce no results
- **Loading states** - Animated spinner while fetching data
