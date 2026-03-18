# Tooltip Customization

Detailed documentation for Chart.js tooltip customization including external tooltips, callbacks, and custom positioners.

## Complete Tooltip Options

Namespace: `options.plugins.tooltip`

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `enabled` | `boolean` | `true` | Enable on-canvas tooltips |
| `external` | `function` | `null` | External tooltip handler |
| `mode` | `string` | `interaction.mode` | Element selection mode |
| `intersect` | `boolean` | `interaction.intersect` | Require intersection |
| `position` | `string` | `'average'` | Positioning mode |
| `callbacks` | `object` | | Text callbacks |
| `itemSort` | `function` | | Sort tooltip items |
| `filter` | `function` | | Filter tooltip items |
| `backgroundColor` | `Color` | `'rgba(0,0,0,0.8)'` | Background color |
| `titleColor` | `Color` | `'#fff'` | Title text color |
| `titleFont` | `Font` | `{weight: 'bold'}` | Title font |
| `titleAlign` | `string` | `'left'` | Title alignment |
| `titleSpacing` | `number` | `2` | Title line spacing |
| `titleMarginBottom` | `number` | `6` | Margin below title |
| `bodyColor` | `Color` | `'#fff'` | Body text color |
| `bodyFont` | `Font` | `{}` | Body font |
| `bodyAlign` | `string` | `'left'` | Body alignment |
| `bodySpacing` | `number` | `2` | Body line spacing |
| `footerColor` | `Color` | `'#fff'` | Footer text color |
| `footerFont` | `Font` | `{weight: 'bold'}` | Footer font |
| `footerAlign` | `string` | `'left'` | Footer alignment |
| `footerSpacing` | `number` | `2` | Footer line spacing |
| `footerMarginTop` | `number` | `6` | Margin above footer |
| `padding` | `Padding` | `6` | Internal padding |
| `caretPadding` | `number` | `2` | Caret distance from tooltip |
| `caretSize` | `number` | `5` | Caret size (px) |
| `cornerRadius` | `number` | `6` | Border radius |
| `multiKeyBackground` | `Color` | `'#fff'` | Background behind color boxes |
| `displayColors` | `boolean` | `true` | Show color boxes |
| `boxWidth` | `number` | `bodyFont.size` | Color box width |
| `boxHeight` | `number` | `bodyFont.size` | Color box height |
| `boxPadding` | `number` | `1` | Padding around color box |
| `usePointStyle` | `boolean` | `false` | Use point style for boxes |
| `borderColor` | `Color` | `'rgba(0,0,0,0)'` | Border color |
| `borderWidth` | `number` | `0` | Border width |
| `rtl` | `boolean` | | Right-to-left |
| `textDirection` | `string` | | Force text direction |
| `xAlign` | `string` | | Caret X position |
| `yAlign` | `string` | | Caret Y position |

## Complete Callbacks List

All callbacks return `string | string[] | undefined`. Returning `undefined` uses the default.

| Callback | Arguments | Per-Dataset | Description |
|----------|-----------|-------------|-------------|
| `beforeTitle` | `TooltipItem[]` | No | Text before title |
| `title` | `TooltipItem[]` | No | Title text |
| `afterTitle` | `TooltipItem[]` | No | Text after title |
| `beforeBody` | `TooltipItem[]` | No | Text before body |
| `beforeLabel` | `TooltipItem` | Yes | Text before each label |
| `label` | `TooltipItem` | Yes | Label text |
| `labelColor` | `TooltipItem` | Yes | Color box styling |
| `labelTextColor` | `TooltipItem` | Yes | Label text color |
| `labelPointStyle` | `TooltipItem` | Yes | Point style for box |
| `afterLabel` | `TooltipItem` | Yes | Text after each label |
| `afterBody` | `TooltipItem[]` | No | Text after body |
| `beforeFooter` | `TooltipItem[]` | No | Text before footer |
| `footer` | `TooltipItem[]` | No | Footer text |
| `afterFooter` | `TooltipItem[]` | No | Text after footer |

### Tooltip Item Context

```javascript
{
  chart: Chart,              // Chart instance
  label: string,             // Label text
  parsed: object,            // Parsed data values
  raw: any,                  // Raw data value
  formattedValue: string,    // Formatted value
  dataset: object,           // Dataset object
  datasetIndex: number,      // Dataset index
  dataIndex: number,         // Data point index
  element: Element           // Chart element
}
```

### Callback Examples

```javascript
callbacks: {
  // Format as currency
  label: (context) => {
    const value = context.parsed.y;
    return new Intl.NumberFormat('en-US', {
      style: 'currency',
      currency: 'USD'
    }).format(value);
  },

  // Custom color box
  labelColor: (context) => ({
    borderColor: 'rgb(0, 0, 255)',
    backgroundColor: 'rgb(255, 0, 0)',
    borderWidth: 2,
    borderDash: [2, 2],
    borderRadius: 2
  }),

  // Custom point style
  labelPointStyle: (context) => ({
    pointStyle: 'triangle',
    rotation: 0
  }),

  // Sum in footer
  footer: (items) => {
    const sum = items.reduce((total, item) => total + item.parsed.y, 0);
    return 'Total: ' + sum;
  }
}
```

### Per-Dataset Callback Override

```javascript
data: {
  datasets: [{
    label: 'Sales',
    data: [10, 20, 30],
    tooltip: {
      callbacks: {
        label: (context) => 'Custom: ' + context.parsed.y
      }
    }
  }]
}
```

## Filter and Sort

```javascript
plugins: {
  tooltip: {
    // Filter out specific items
    filter: (tooltipItem, data) => {
      return tooltipItem.parsed.y > 0;  // Only positive values
    },

    // Sort items
    itemSort: (a, b, data) => {
      return b.parsed.y - a.parsed.y;  // Descending by value
    }
  }
}
```

## External (HTML) Tooltips

Create custom HTML tooltips instead of canvas-rendered:

```javascript
plugins: {
  tooltip: {
    enabled: false,  // Disable canvas tooltip
    external: function(context) {
      const { chart, tooltip } = context;

      // Get or create tooltip element
      let tooltipEl = document.getElementById('chartjs-tooltip');
      if (!tooltipEl) {
        tooltipEl = document.createElement('div');
        tooltipEl.id = 'chartjs-tooltip';
        document.body.appendChild(tooltipEl);
      }

      // Hide if no tooltip
      if (tooltip.opacity === 0) {
        tooltipEl.style.opacity = 0;
        return;
      }

      // Set content
      if (tooltip.body) {
        const titleLines = tooltip.title || [];
        const bodyLines = tooltip.body.map(b => b.lines);

        let html = '<div class="tooltip-title">' + titleLines.join('<br>') + '</div>';
        html += '<div class="tooltip-body">';
        bodyLines.forEach((lines, i) => {
          const colors = tooltip.labelColors[i];
          const style = `background:${colors.backgroundColor};border-color:${colors.borderColor}`;
          html += `<div><span class="color-box" style="${style}"></span>${lines}</div>`;
        });
        html += '</div>';
        tooltipEl.innerHTML = html;
      }

      // Position tooltip
      const position = chart.canvas.getBoundingClientRect();
      tooltipEl.style.opacity = 1;
      tooltipEl.style.position = 'absolute';
      tooltipEl.style.left = position.left + window.pageXOffset + tooltip.caretX + 'px';
      tooltipEl.style.top = position.top + window.pageYOffset + tooltip.caretY + 'px';
      tooltipEl.style.pointerEvents = 'none';
    }
  }
}
```

## Tooltip Model

The tooltip model passed to external handlers:

```javascript
{
  chart: Chart,
  dataPoints: TooltipItem[],   // Items in tooltip

  // Positioning
  xAlign: string,
  yAlign: string,
  x: number,                   // Top-left X
  y: number,                   // Top-left Y
  width: number,
  height: number,
  caretX: number,              // Caret X position
  caretY: number,              // Caret Y position

  // Content
  title: string[],
  beforeBody: string[],
  body: [{
    before: string[],
    lines: string[],
    after: string[]
  }],
  afterBody: string[],
  footer: string[],

  // Styling
  labelColors: [{
    borderColor: Color,
    backgroundColor: Color
  }],
  labelTextColors: Color[],
  labelPointStyles: [{
    pointStyle: PointStyle,
    rotation: number
  }],

  opacity: number,             // 0 = hidden
  options: object              // Tooltip options
}
```

## Custom Position Modes

Register custom tooltip positioners:

```javascript
import { Tooltip } from 'chart.js';

Tooltip.positioners.topLeft = function(elements, eventPosition) {
  const tooltip = this;  // Tooltip instance

  if (!elements.length) {
    return false;
  }

  return {
    x: 0,
    y: 0,
    xAlign: 'left',   // Optional: override xAlign
    yAlign: 'top'     // Optional: override yAlign
  };
};

// Use it
new Chart(ctx, {
  options: {
    plugins: {
      tooltip: {
        position: 'topLeft'
      }
    }
  }
});
```

### TypeScript Registration

```typescript
declare module 'chart.js' {
  interface TooltipPositionerMap {
    topLeft: TooltipPositionerFunction<ChartType>;
  }
}
```

## Tooltip Alignment

Control caret position:

```javascript
plugins: {
  tooltip: {
    xAlign: 'center',  // 'left', 'center', 'right'
    yAlign: 'bottom'   // 'top', 'center', 'bottom'
  }
}
```

## Per-Plugin Events

Limit which events trigger tooltips:

```javascript
options: {
  events: ['mousemove', 'mouseout', 'click', 'touchstart', 'touchmove'],
  plugins: {
    tooltip: {
      events: ['click']  // Tooltip only on click
    }
  }
}
```
