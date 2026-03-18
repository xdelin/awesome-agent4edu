# Tooltip Model Reference

The tooltip model contains all parameters available for rendering custom tooltips. This is
essential when building external (HTML) tooltips.

## Tooltip Model Interface

```javascript
{
  // The chart instance
  chart: Chart,

  // The tooltip items (see TooltipItem interface below)
  dataPoints: TooltipItem[],

  // Caret alignment
  xAlign: string,  // 'left' | 'center' | 'right'
  yAlign: string,  // 'top' | 'center' | 'bottom'

  // Position (top-left corner of tooltip)
  x: number,
  y: number,
  width: number,
  height: number,

  // Where the tooltip arrow points
  caretX: number,
  caretY: number,

  // Body content
  body: [{
    before: string[],  // Lines before main item
    lines: string[],   // Main item lines with color square
    after: string[]    // Lines after main item
  }],
  beforeBody: string[],  // Lines after title, before body
  afterBody: string[],   // Lines after body, before footer

  // Title and footer
  title: string[],
  footer: string[],

  // Styling for each body item
  labelColors: [{
    borderColor: Color,
    backgroundColor: Color,
    borderWidth: number,
    borderDash: number[],
    borderRadius: number
  }],
  labelTextColors: Color[],
  labelPointStyles: [{
    pointStyle: PointStyle,
    rotation: number
  }],

  // Visibility (0 = hidden)
  opacity: number,

  // Current tooltip options
  options: Object
}
```

## TooltipItem Interface

Each item in `dataPoints` implements this interface:

```javascript
{
  // The chart instance
  chart: Chart,

  // Label for the tooltip
  label: string,

  // Parsed data values for this point
  parsed: {
    x: number,  // For Cartesian charts
    y: number,
    r: number   // For bubble charts
  },

  // Raw data value (before parsing)
  raw: any,

  // Formatted value string
  formattedValue: string,

  // Dataset reference and index
  dataset: Object,
  datasetIndex: number,

  // Data point index within dataset
  dataIndex: number,

  // The chart element (point, arc, bar, etc.)
  element: Element
}
```

## Using the Model in External Tooltips

```javascript
const externalTooltipHandler = (context) => {
  const { chart, tooltip } = context;

  // Access the model directly
  const model = tooltip;

  // Check visibility
  if (model.opacity === 0) {
    hideTooltip();
    return;
  }

  // Get positioning
  const { caretX, caretY } = model;
  const position = chart.canvas.getBoundingClientRect();

  // Calculate absolute position
  const left = position.left + window.pageXOffset + caretX;
  const top = position.top + window.pageYOffset + caretY;

  // Access data points
  model.dataPoints.forEach((point) => {
    console.log(point.label, point.formattedValue);
    console.log('Dataset:', point.dataset.label);
    console.log('Value:', point.parsed.y);
  });

  // Access body lines (already formatted)
  model.body.forEach((bodyItem, i) => {
    const colors = model.labelColors[i];
    bodyItem.lines.forEach((line) => {
      // Render line with colors.backgroundColor
    });
  });

  // Access title and footer
  model.title.forEach((titleLine) => {
    // Render title
  });

  model.footer.forEach((footerLine) => {
    // Render footer
  });
};
```

## Common Patterns

### Get Current Font

```javascript
const bodyFont = Chart.helpers.toFont(tooltip.options.bodyFont);
tooltipEl.style.font = bodyFont.string;
```

### Apply Padding

```javascript
tooltipEl.style.padding = tooltip.options.padding + 'px';
```

### Handle Caret Alignment

```javascript
// Position based on alignment
if (tooltip.yAlign === 'top') {
  tooltipEl.style.transform = 'translateY(10px)';
} else if (tooltip.yAlign === 'bottom') {
  tooltipEl.style.transform = 'translateY(-10px)';
}
```
