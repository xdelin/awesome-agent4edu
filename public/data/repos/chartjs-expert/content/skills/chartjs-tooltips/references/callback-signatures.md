# Tooltip Callback Signatures

Complete reference for all tooltip callbacks in Chart.js v4.5.1.

## Callback Overview

Namespace: `options.plugins.tooltip.callbacks`

All callbacks receive `this` bound to the tooltip object. Return `undefined` to use the default,
or an empty string to hide the element.

## Title Callbacks

### beforeTitle

```javascript
callbacks: {
  beforeTitle: function(tooltipItems) {
    // tooltipItems: TooltipItem[]
    return 'Before Title Text';  // string | string[] | undefined
  }
}
```

### title

```javascript
callbacks: {
  title: function(tooltipItems) {
    // Default: returns the label from the first item
    return tooltipItems[0].label;  // string | string[] | undefined
  }
}
```

### afterTitle

```javascript
callbacks: {
  afterTitle: function(tooltipItems) {
    return 'After Title Text';  // string | string[] | undefined
  }
}
```

## Body Callbacks

### beforeBody

```javascript
callbacks: {
  beforeBody: function(tooltipItems) {
    return 'Before all body items';  // string | string[] | undefined
  }
}
```

### beforeLabel

Supports per-dataset override.

```javascript
callbacks: {
  beforeLabel: function(tooltipItem) {
    // tooltipItem: single TooltipItem
    return `Dataset: ${tooltipItem.dataset.label}`;  // string | string[] | undefined
  }
}
```

### label

Supports per-dataset override.

```javascript
callbacks: {
  label: function(tooltipItem) {
    // Default format: "dataset.label: formattedValue"
    let label = tooltipItem.dataset.label || '';
    if (label) label += ': ';
    if (tooltipItem.parsed.y !== null) {
      label += tooltipItem.formattedValue;
    }
    return label;  // string | string[] | undefined
  }
}
```

### afterLabel

Supports per-dataset override.

```javascript
callbacks: {
  afterLabel: function(tooltipItem) {
    return 'After this item';  // string | string[] | undefined
  }
}
```

### afterBody

```javascript
callbacks: {
  afterBody: function(tooltipItems) {
    return 'After all body items';  // string | string[] | undefined
  }
}
```

## Footer Callbacks

### beforeFooter

```javascript
callbacks: {
  beforeFooter: function(tooltipItems) {
    return 'Before Footer';  // string | string[] | undefined
  }
}
```

### footer

```javascript
callbacks: {
  footer: function(tooltipItems) {
    // Calculate sum example
    let sum = 0;
    tooltipItems.forEach((item) => {
      sum += item.parsed.y;
    });
    return 'Sum: ' + sum;  // string | string[] | undefined
  }
}
```

### afterFooter

```javascript
callbacks: {
  afterFooter: function(tooltipItems) {
    return 'After Footer';  // string | string[] | undefined
  }
}
```

## Styling Callbacks

### labelColor

Supports per-dataset override.

```javascript
callbacks: {
  labelColor: function(tooltipItem) {
    return {
      borderColor: 'rgb(0, 0, 255)',
      backgroundColor: 'rgb(255, 0, 0)',
      borderWidth: 2,
      borderDash: [2, 2],
      borderRadius: 2
    };  // object | undefined
  }
}
```

### labelTextColor

Supports per-dataset override.

```javascript
callbacks: {
  labelTextColor: function(tooltipItem) {
    return '#543453';  // Color | undefined
  }
}
```

### labelPointStyle

Supports per-dataset override.

Requires `usePointStyle: true` in tooltip options.

```javascript
callbacks: {
  labelPointStyle: function(tooltipItem) {
    return {
      pointStyle: 'triangle',  // 'circle' | 'cross' | 'star' | 'triangle' | etc.
      rotation: 0
    };  // object | undefined
  }
}
```

## Per-Dataset Overrides

Callbacks marked with "Supports per-dataset override" can be defined per dataset:

```javascript
const data = {
  datasets: [{
    label: 'Sales',
    data: [10, 20, 30],
    tooltip: {
      callbacks: {
        label: function(context) {
          return `Sales: $${context.parsed.y}`;
        },
        labelColor: function(context) {
          return {
            backgroundColor: 'green',
            borderColor: 'darkgreen'
          };
        }
      }
    }
  }, {
    label: 'Expenses',
    data: [5, 15, 25],
    tooltip: {
      callbacks: {
        label: function(context) {
          return `Expenses: -$${context.parsed.y}`;
        },
        labelColor: function(context) {
          return {
            backgroundColor: 'red',
            borderColor: 'darkred'
          };
        }
      }
    }
  }]
};
```

## Filter and Sort Callbacks

### filter

Filter which items appear in the tooltip:

```javascript
options: {
  plugins: {
    tooltip: {
      filter: function(tooltipItem, data) {
        // Only show items with value > 10
        return tooltipItem.parsed.y > 10;
      }
    }
  }
}
```

### itemSort

Sort tooltip items:

```javascript
options: {
  plugins: {
    tooltip: {
      itemSort: function(a, b, data) {
        // Sort by value descending
        return b.parsed.y - a.parsed.y;
      }
    }
  }
}
```

## Multi-Line Returns

Any callback can return an array for multiple lines:

```javascript
callbacks: {
  label: function(tooltipItem) {
    return [
      `Value: ${tooltipItem.parsed.y}`,
      `Index: ${tooltipItem.dataIndex}`,
      `Dataset: ${tooltipItem.dataset.label}`
    ];
  },
  footer: function(tooltipItems) {
    const sum = tooltipItems.reduce((a, b) => a + b.parsed.y, 0);
    const avg = sum / tooltipItems.length;
    return [
      `Total: ${sum}`,
      `Average: ${avg.toFixed(2)}`
    ];
  }
}
```
