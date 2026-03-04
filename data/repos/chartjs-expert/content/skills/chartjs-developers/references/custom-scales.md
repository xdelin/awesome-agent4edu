# Custom Scales

Create custom axis types by extending `Chart.Scale`.

## Basic Structure

```javascript
class MyScale extends Chart.Scale {
  static id = 'myScale';
  static defaults = {
    // Default configuration
  };
}

Chart.register(MyScale);

new Chart(ctx, {
  type: 'line',
  options: {
    scales: {
      y: {
        type: 'myScale'  // Use custom scale
      }
    }
  }
});
```

## Required Interface

Custom scales must implement these methods:

```javascript
class MyScale extends Chart.Scale {
  // Determine data min/max - set this.min and this.max
  determineDataLimits() {}

  // Generate tick marks - create this.ticks array
  buildTicks() {}

  // Get label for a value
  getLabelForValue(value) {}

  // Get pixel position for tick index
  getPixelForTick(index) {}

  // Get pixel position for data value
  getPixelForValue(value, index) {}

  // Get data value for pixel position
  getValueForPixel(pixel) {}
}
```

## Optional Overrides

Base class provides default implementations:

```javascript
class MyScale extends Chart.Scale {
  // Add labels to tick objects
  // Default: calls this.options.ticks.callback(tick, index, ticks)
  generateTickLabels() {}

  // Calculate label rotation
  // Default: only rotates for horizontal scales
  calculateLabelRotation() {}

  // Fit scale into canvas
  // Use: this.maxWidth, this.maxHeight, this.margins
  // Set: this.minSize, this.width, this.height
  fit() {}

  // Draw scale on canvas
  // Use: this.left, this.right, this.top, this.bottom
  draw(chartArea) {}
}
```

## Scale Properties

Properties available during fitting:

```javascript
{
  // Bounding box edges
  left: number,
  right: number,
  top: number,
  bottom: number,
  width: number,   // right - left
  height: number,  // bottom - top

  // Margins (outside bounding box)
  margins: {
    left: number,
    right: number,
    top: number,
    bottom: number
  },

  // Padding (inside bounding box)
  paddingLeft: number,
  paddingRight: number,
  paddingTop: number,
  paddingBottom: number
}
```

## Utility Methods

Inherited from `Chart.Scale`:

```javascript
// Check if scale is horizontal
scale.isHorizontal();

// Get tick objects ({label, major})
scale.getTicks();
```

## Alternative Registration

For scales not extending `Chart.Scale`:

```javascript
// Explicit registration
Chart.registry.addScales(MyScale);
```
