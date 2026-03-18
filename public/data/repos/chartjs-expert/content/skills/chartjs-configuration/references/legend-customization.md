# Legend Customization

Detailed documentation for Chart.js legend customization including custom click handlers, the Legend Item Interface, and filtering/sorting.

## Complete Legend Options

Namespace: `options.plugins.legend`

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `display` | `boolean` | `true` | Show legend |
| `position` | `string` | `'top'` | Position (`top`, `bottom`, `left`, `right`, `chartArea`) |
| `align` | `string` | `'center'` | Alignment (`start`, `center`, `end`) |
| `maxHeight` | `number` | - | Maximum height in pixels |
| `maxWidth` | `number` | - | Maximum width in pixels |
| `fullSize` | `boolean` | `true` | Take full canvas width/height |
| `onClick` | `function` | | Click handler |
| `onHover` | `function` | | Hover handler |
| `onLeave` | `function` | | Mouse leave handler |
| `reverse` | `boolean` | `false` | Reverse dataset order |
| `labels` | `object` | | Label configuration |
| `rtl` | `boolean` | | Right-to-left rendering |
| `textDirection` | `string` | | Force text direction (`'rtl'` or `'ltr'`) |
| `title` | `object` | | Legend title configuration |

## Legend Labels Options

Namespace: `options.plugins.legend.labels`

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `boxWidth` | `number` | `40` | Color box width |
| `boxHeight` | `number` | `font.size` | Color box height |
| `color` | `Color` | `Chart.defaults.color` | Text color |
| `font` | `Font` | `Chart.defaults.font` | Font settings |
| `padding` | `number` | `10` | Padding between legend items |
| `generateLabels` | `function` | | Custom label generator |
| `filter` | `function` | | Filter legend items |
| `sort` | `function` | | Sort legend items |
| `pointStyle` | `PointStyle` | `'circle'` | Point style (when `usePointStyle: true`) |
| `textAlign` | `string` | `'center'` | Text alignment (`left`, `center`, `right`) |
| `usePointStyle` | `boolean` | `false` | Use point style instead of boxes |
| `pointStyleWidth` | `number` | | Point style width |
| `useBorderRadius` | `boolean` | `false` | Match dataset border radius |
| `borderRadius` | `number` | | Override border radius |

## Legend Title Options

Namespace: `options.plugins.legend.title`

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `display` | `boolean` | `false` | Show legend title |
| `text` | `string` | | Title text |
| `color` | `Color` | `Chart.defaults.color` | Text color |
| `font` | `Font` | `Chart.defaults.font` | Font settings |
| `padding` | `Padding` | `0` | Title padding |

## Legend Item Interface

Legend items returned by `generateLabels` and passed to click handlers:

```javascript
{
  text: string,                    // Display text
  borderRadius: number | object,   // Border radius
  datasetIndex: number,            // Associated dataset index
  fillStyle: Color,                // Box fill color
  fontColor: Color,                // Text color
  hidden: boolean,                 // Is dataset hidden (shows strikethrough)
  lineCap: string,                 // Box border cap style
  lineDash: number[],              // Box border dash pattern
  lineDashOffset: number,          // Box border dash offset
  lineJoin: string,                // Box border join style
  lineWidth: number,               // Box border width
  strokeStyle: Color,              // Box border color
  pointStyle: string | Image,      // Point style (if usePointStyle)
  rotation: number                 // Point rotation (degrees)
}
```

## Event Handlers

### onClick

Called when clicking a legend item:

```javascript
onClick: function(event, legendItem, legend) {
  // event: Native click event
  // legendItem: The clicked Legend Item
  // legend: The legend plugin instance

  const index = legendItem.datasetIndex;
  const chart = legend.chart;

  // Default behavior: toggle visibility
  if (chart.isDatasetVisible(index)) {
    chart.hide(index);
    legendItem.hidden = true;
  } else {
    chart.show(index);
    legendItem.hidden = false;
  }
}
```

### onHover

Called when hovering over a legend item:

```javascript
onHover: function(event, legendItem, legend) {
  const chart = legend.chart;
  chart.canvas.style.cursor = 'pointer';

  // Highlight the hovered dataset
  const datasetIndex = legendItem.datasetIndex;
  chart.data.datasets.forEach((dataset, i) => {
    if (i !== datasetIndex) {
      dataset.backgroundColor = 'rgba(128, 128, 128, 0.2)';
    }
  });
  chart.update();
}
```

### onLeave

Called when mouse leaves a legend item:

```javascript
onLeave: function(event, legendItem, legend) {
  const chart = legend.chart;
  chart.canvas.style.cursor = 'default';

  // Restore original colors
  chart.data.datasets.forEach((dataset, i) => {
    dataset.backgroundColor = originalColors[i];
  });
  chart.update();
}
```

## Custom Click Handlers

### Link Multiple Datasets

Toggle multiple datasets together:

```javascript
const defaultLegendClickHandler = Chart.defaults.plugins.legend.onClick;
const pieDoughnutLegendClickHandler = Chart.controllers.doughnut.overrides.plugins.legend.onClick;

plugins: {
  legend: {
    onClick: function(e, legendItem, legend) {
      const index = legendItem.datasetIndex;
      const type = legend.chart.config.type;

      // Link first two datasets
      if (index <= 1) {
        const chart = legend.chart;
        [chart.getDatasetMeta(0), chart.getDatasetMeta(1)].forEach(meta => {
          meta.hidden = meta.hidden === null
            ? !chart.data.datasets[index].hidden
            : null;
        });
        chart.update();
      } else {
        // Default behavior for other datasets
        if (type === 'pie' || type === 'doughnut') {
          pieDoughnutLegendClickHandler(e, legendItem, legend);
        } else {
          defaultLegendClickHandler(e, legendItem, legend);
        }
      }
    }
  }
}
```

### Single Selection Mode

Only one dataset visible at a time:

```javascript
plugins: {
  legend: {
    onClick: function(e, legendItem, legend) {
      const chart = legend.chart;
      const clickedIndex = legendItem.datasetIndex;

      chart.data.datasets.forEach((dataset, index) => {
        const meta = chart.getDatasetMeta(index);
        meta.hidden = index !== clickedIndex;
      });

      chart.update();
    }
  }
}
```

### Prevent Toggle

Make certain datasets always visible:

```javascript
plugins: {
  legend: {
    onClick: function(e, legendItem, legend) {
      const index = legendItem.datasetIndex;

      // Prevent hiding the first dataset
      if (index === 0 && legend.chart.isDatasetVisible(0)) {
        return;
      }

      // Default behavior
      Chart.defaults.plugins.legend.onClick.call(this, e, legendItem, legend);
    }
  }
}
```

## Filtering Legend Items

Remove items from legend:

```javascript
plugins: {
  legend: {
    labels: {
      filter: function(legendItem, data) {
        // Hide legend items for datasets with 'hidden' label
        return legendItem.text !== 'Hidden';
      }
    }
  }
}
```

### Filter by Dataset Property

```javascript
filter: function(legendItem, data) {
  const dataset = data.datasets[legendItem.datasetIndex];
  // Only show if dataset has showInLegend: true
  return dataset.showInLegend !== false;
}
```

## Sorting Legend Items

Sort legend item order:

```javascript
plugins: {
  legend: {
    labels: {
      sort: function(a, b, data) {
        // Alphabetical order
        return a.text.localeCompare(b.text);
      }
    }
  }
}
```

### Sort by Dataset Value

```javascript
sort: function(a, b, data) {
  // Sort by sum of dataset values (descending)
  const sumA = data.datasets[a.datasetIndex].data.reduce((s, v) => s + v, 0);
  const sumB = data.datasets[b.datasetIndex].data.reduce((s, v) => s + v, 0);
  return sumB - sumA;
}
```

## Custom Label Generation

Override default label generation:

```javascript
plugins: {
  legend: {
    labels: {
      generateLabels: function(chart) {
        const datasets = chart.data.datasets;
        return datasets.map((dataset, i) => ({
          text: `${dataset.label} (${dataset.data.length} points)`,
          fillStyle: dataset.backgroundColor,
          strokeStyle: dataset.borderColor,
          lineWidth: dataset.borderWidth,
          hidden: !chart.isDatasetVisible(i),
          datasetIndex: i
        }));
      }
    }
  }
}
```

## Point Style Legend

Use dataset point styles in legend:

```javascript
plugins: {
  legend: {
    labels: {
      usePointStyle: true,
      pointStyle: 'circle',        // Default point style
      pointStyleWidth: 20          // Width of point style
    }
  }
}
```

## Chart Type Overrides

Doughnut, pie, and polar charts have different legend defaults. Override them:

```javascript
// Access chart type overrides
Chart.overrides.doughnut.plugins.legend.position = 'right';
Chart.overrides.pie.plugins.legend.labels.usePointStyle = true;
```

## HTML Legend

For complex visual customizations, use an HTML legend:

```javascript
const htmlLegendPlugin = {
  id: 'htmlLegend',
  afterUpdate(chart, args, options) {
    const container = document.getElementById('legend-container');

    // Clear existing
    container.innerHTML = '';

    // Build HTML
    const ul = document.createElement('ul');
    chart.data.datasets.forEach((dataset, i) => {
      const li = document.createElement('li');
      li.innerHTML = `
        <span class="color" style="background: ${dataset.backgroundColor}"></span>
        <span class="text">${dataset.label}</span>
      `;
      li.onclick = () => {
        chart.setDatasetVisibility(i, !chart.isDatasetVisible(i));
        chart.update();
      };
      if (!chart.isDatasetVisible(i)) {
        li.classList.add('hidden');
      }
      ul.appendChild(li);
    });
    container.appendChild(ul);
  }
};

// Use the plugin
new Chart(ctx, {
  plugins: [htmlLegendPlugin],
  options: {
    plugins: {
      legend: { display: false }  // Hide default legend
    }
  }
});
```
