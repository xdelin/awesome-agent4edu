# Chart.js Instance Methods Reference

Complete guide to methods available on Chart.js chart instances for dynamic updates and manipulation.

## Core Update Methods

### update()

Update the chart with new data or options, with optional animation control.

```javascript
const chart = new Chart(ctx, config);

// Standard update with default animation
chart.data.datasets[0].data = [10, 20, 30, 40];
chart.update();

// Update without animation
chart.update('none');

// Update with specific animation mode
chart.update('active');  // Only animate active elements
chart.update('show');    // Show animation
chart.update('hide');    // Hide animation
chart.update('resize');  // Resize animation
```

**Animation Modes:**

| Mode | Description |
|------|-------------|
| `undefined` (default) | Uses chart's configured animation |
| `'none'` | No animation |
| `'active'` | Only animates active elements |
| `'show'` | Show animation for new elements |
| `'hide'` | Hide animation for removed elements |
| `'resize'` | Animation for resize operations |

**Best Practice:** Use `update('none')` for rapid data changes or when performance is critical.

### reset()

Reset the chart to its initial state before any animations or interactions.

```javascript
chart.reset();
```

This is useful for:
- Returning chart to original appearance
- Resetting after user interactions
- Preparing chart for new animation sequence

### resize()

Manually trigger a resize of the chart canvas.

```javascript
// Auto-detect size from container
chart.resize();

// Set specific dimensions
chart.resize(800, 400);
```

**When to use:**
- Container size changed programmatically
- Chart inside collapsible/hidden element that becomes visible
- Manual override of responsive behavior

**Note:** Chart.js automatically handles window resize events when `responsive: true` (default).

### destroy()

Completely destroy the chart instance, removing all references and event listeners.

```javascript
chart.destroy();
```

**Important:** Always call `destroy()` before removing a chart from the DOM to prevent memory leaks.

**Use cases:**
- Removing chart from single-page applications
- Replacing one chart type with another
- Cleaning up before page navigation

```javascript
// Proper cleanup pattern
if (chart) {
  chart.destroy();
  chart = null;
}

// Create new chart
chart = new Chart(ctx, newConfig);
```

## Data Access Methods

### getElementsAtEventForMode()

Get chart elements at a specific event position using interaction modes.

```javascript
const elements = chart.getElementsAtEventForMode(
  event,           // Mouse/touch event
  'nearest',       // Interaction mode
  { intersect: true },  // Options
  true             // Use final position
);

// Returns array of elements
elements.forEach(element => {
  console.log('Dataset index:', element.datasetIndex);
  console.log('Data index:', element.index);
});
```

**Interaction Modes:**

| Mode | Description |
|------|-------------|
| `'point'` | Finds elements at exact point |
| `'nearest'` | Finds nearest element |
| `'index'` | Finds elements at same index across datasets |
| `'dataset'` | Finds elements in same dataset |
| `'x'` | Finds elements with matching x value |
| `'y'` | Finds elements with matching y value |

**Options:**

- `intersect: true` - Only return elements directly under the point
- `intersect: false` - Return elements in the same index/dataset even if not directly under point
- `axis: 'x'` or `axis: 'y'` - Restrict search to one axis

### getSortedVisibleDatasetMetas()

Get metadata for all visible datasets, sorted by drawing order.

```javascript
const metas = chart.getSortedVisibleDatasetMetas();

metas.forEach(meta => {
  console.log('Dataset label:', meta.label);
  console.log('Dataset type:', meta.type);
  console.log('Visible:', meta.visible);
});
```

## Export Methods

### toBase64Image()

Generate a base64-encoded PNG image of the current chart.

```javascript
const image = chart.toBase64Image();

// Use in img tag
document.getElementById('chartImage').src = image;

// Download as file
const link = document.createElement('a');
link.download = 'chart.png';
link.href = image;
link.click();
```

**Options:**

```javascript
// Specify image type and quality
const jpeg = chart.toBase64Image('image/jpeg', 0.8);  // 80% quality JPEG
const png = chart.toBase64Image('image/png');         // PNG (default)
```

**Browser Support:** Works in all modern browsers. Uses canvas `toDataURL()` internally.

## Configuration Access

### Chart Properties

Access chart configuration and state through instance properties:

```javascript
// Access chart data
console.log(chart.data);

// Access chart options
console.log(chart.options);

// Access chart configuration
console.log(chart.config);

// Access canvas element
console.log(chart.canvas);

// Access 2D context
console.log(chart.ctx);

// Access chart width/height
console.log(chart.width, chart.height);

// Access scales
console.log(chart.scales.x);
console.log(chart.scales.y);

// Check if chart is destroyed
console.log(chart.isDestroyed);
```

## Animation Control

### stop()

Stop all running animations immediately.

```javascript
chart.stop();
```

Useful for:
- Pausing animations on user request
- Stopping animations before destroying chart
- Performance optimization

### render()

Manually trigger a chart render without updating data.

```javascript
chart.render();
```

**Difference from update():**
- `render()` - Redraws chart without recalculating data/scales
- `update()` - Recalculates everything then renders

Use `render()` when only visual changes are needed without data changes.

## Common Patterns

### Periodic Data Updates

```javascript
setInterval(() => {
  // Add new data point
  chart.data.labels.push(new Date().toLocaleTimeString());
  chart.data.datasets[0].data.push(Math.random() * 100);

  // Keep only last 20 points
  if (chart.data.labels.length > 20) {
    chart.data.labels.shift();
    chart.data.datasets[0].data.shift();
  }

  chart.update('none');  // Update without animation for smooth real-time updates
}, 1000);
```

### Click Handling with Data Extraction

```javascript
chart.options.onClick = (event, activeElements) => {
  if (activeElements.length > 0) {
    const element = activeElements[0];
    const datasetIndex = element.datasetIndex;
    const index = element.index;
    const value = chart.data.datasets[datasetIndex].data[index];
    const label = chart.data.labels[index];

    console.log(`Clicked: ${label} = ${value}`);
  }
};
```

### Toggling Dataset Visibility

```javascript
// Hide a specific dataset
chart.hide(datasetIndex);

// Show a specific dataset
chart.show(datasetIndex);

// Toggle dataset visibility
chart.isDatasetVisible(datasetIndex)
  ? chart.hide(datasetIndex)
  : chart.show(datasetIndex);
```

### Responsive Container Handling

```javascript
// Chart inside a tab or collapsible panel
function showChart() {
  // Make container visible
  document.getElementById('chartContainer').style.display = 'block';

  // Trigger resize after container is visible
  setTimeout(() => {
    chart.resize();
  }, 100);
}
```

## Performance Tips

1. **Batch Updates:** Modify multiple datasets, then call `update()` once
2. **Disable Animation:** Use `update('none')` for frequent updates
3. **Destroy Unused Charts:** Always call `destroy()` to free memory
4. **Debounce Resize:** Don't call `resize()` on every pixel change

```javascript
// Good: Batch updates
chart.data.datasets[0].data = newData1;
chart.data.datasets[1].data = newData2;
chart.update();  // Single update call

// Bad: Multiple updates
chart.data.datasets[0].data = newData1;
chart.update();
chart.data.datasets[1].data = newData2;
chart.update();  // Unnecessary second update
```
