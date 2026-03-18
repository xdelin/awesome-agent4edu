# Updating Charts

Charts can be updated after creation. Chart.js animates to new values automatically.

## Adding and Removing Data

```javascript
function addData(chart, label, newData) {
  chart.data.labels.push(label);
  chart.data.datasets.forEach((dataset) => {
    dataset.data.push(newData);
  });
  chart.update();
}

function removeData(chart) {
  chart.data.labels.pop();
  chart.data.datasets.forEach((dataset) => {
    dataset.data.pop();
  });
  chart.update();
}
```

## Updating Options

### Mutating In Place

Preserves existing options and calculated values:

```javascript
function updateByMutating(chart) {
  chart.options.plugins.title.text = 'New Title';
  chart.update();
}
```

### Replacing Options Object

Like creating a new chart - old options are discarded:

```javascript
function updateAsNewObject(chart) {
  chart.options = {
    responsive: true,
    plugins: {
      title: {
        display: true,
        text: 'Chart.js'
      }
    },
    scales: {
      x: { display: true },
      y: { display: true }
    }
  };
  chart.update();
}
```

## Updating Scales

Pass complete scale configuration when updating:

```javascript
function updateScales(chart) {
  chart.options.scales = {
    newId: {
      display: true
    },
    y: {
      display: true,
      type: 'logarithmic'
    }
  };
  chart.update();

  // Update references - old ones are invalid
  const xScale = chart.scales.newId;
  const yScale = chart.scales.y;
}
```

Update single scale by ID:

```javascript
function updateSingleScale(chart) {
  chart.options.scales.y = {
    type: 'logarithmic'
  };
  chart.update();
}
```

## Update Modes

Control animation behavior:

```javascript
// With animation (default)
chart.update();

// Without animation
chart.update('none');

// Active elements only
chart.update('active');

// Per-dataset animation
chart.update(ctx => ctx.datasetIndex === 0 ? 'active' : 'none');
```

### Available Modes

| Mode | Behavior |
|------|----------|
| `undefined` | Default animation |
| `'none'` | Skip animation |
| `'active'` | Animate active elements |
| `'hide'` | Hide animation |
| `'show'` | Show animation |
| `'reset'` | Reset animation |
| `'resize'` | Resize animation |
