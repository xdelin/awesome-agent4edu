# Scale API Reference

Complete API documentation for extending Chart.js scales.

## Scale Interface

All custom scale types must extend `Chart.Scale` and implement this interface.

### Required Properties

| Property | Type | Description |
|----------|------|-------------|
| `id` | `string` | Unique identifier for the scale (used as `type` in scale config) |
| `defaults` | `object` | Default configuration options |

### Required Methods

#### determineDataLimits()

Determines the data limits. Should set `this.min` and `this.max` to be the data max/min.

```javascript
determineDataLimits() {
  const {min, max} = this.getMinMax(true);
  this.min = min;
  this.max = max;
}
```

#### buildTicks()

Generate tick marks. Should create a ticks array on the axis instance.

```javascript
buildTicks() {
  const ticks = [];
  // Generate ticks based on min/max
  let value = this.min;
  while (value <= this.max) {
    ticks.push({value});
    value += this.stepSize;
  }
  return ticks;
}
```

#### getLabelForValue(value)

Get the label to show for the given value.

```javascript
getLabelForValue(value) {
  return value.toFixed(2);
}
```

#### getPixelForValue(value, index)

Get the pixel coordinate for a given value.

- For horizontal axis: returns x coordinate
- For vertical axis: returns y coordinate

```javascript
getPixelForValue(value, index) {
  const decimal = (value - this.min) / (this.max - this.min);
  return this.getPixelForDecimal(decimal);
}
```

#### getValueForPixel(pixel)

Get the value for a given pixel coordinate.

```javascript
getValueForPixel(pixel) {
  const decimal = this.getDecimalForPixel(pixel);
  return this.min + decimal * (this.max - this.min);
}
```

### Optional Override Methods

#### generateTickLabels()

Adds labels to objects in the ticks array. Default calls `this.options.ticks.callback(numericalTick, index, ticks)`.

#### calculateLabelRotation()

Determine how much the labels will rotate. Default only rotates labels if the scale is horizontal.

#### fit()

Fits the scale into the canvas.

- `this.maxWidth` and `this.maxHeight` provide maximum dimensions
- `this.margins` is the amount of space on either side for expansion
- Must set `this.minSize` to `{width, height}`
- Must set `this.width` and `this.height`

#### draw(chartArea)

Draws the scale onto the canvas.

- `this.(left|right|top|bottom)` indicate the drawing area
- `chartArea` contains `{left, right, top, bottom}` for drawing grid lines

## Scale Properties

Scale instances receive these properties during the fitting process:

```javascript
{
  left: number,    // left edge of the scale bounding box
  right: number,   // right edge of the bounding box
  top: number,
  bottom: number,
  width: number,   // same as right - left
  height: number,  // same as bottom - top

  // Margin on each side (outside the bounding box, like CSS)
  margins: {
    left: number,
    right: number,
    top: number,
    bottom: number
  },

  // Padding inside the bounding box (like CSS)
  paddingLeft: number,
  paddingRight: number,
  paddingTop: number,
  paddingBottom: number
}
```

## Utility Methods

The `Chart.Scale` base class provides these utility functions:

| Method | Description |
|--------|-------------|
| `isHorizontal()` | Returns true if the scale is horizontal |
| `getTicks()` | Returns the scale tick objects (`{label, major}`) |
| `getMinMax(canStack)` | Get min/max values from datasets (see Log2 example below) |
| `getPixelForDecimal(decimal)` | Convert 0-1 decimal to pixel position |
| `getDecimalForPixel(pixel)` | Convert pixel position to 0-1 decimal |

**Using `getMinMax()`**: Call in `determineDataLimits()` to get data bounds. Pass `true` to allow stacking consideration. See the Log2 scale example for typical usage pattern.

## Registration

Always register custom scales:

```javascript
// ES modules
import {Chart} from 'chart.js';
Chart.register(MyScale);

// Alternative explicit registration
Chart.registry.addScales(MyScale);
```

## Complete Example: Log2 Scale

```javascript
import {Chart, Scale} from 'chart.js';

class Log2Axis extends Scale {
  constructor(cfg) {
    super(cfg);
    this.start = undefined;
    this.end = undefined;
  }

  parse(raw, index) {
    const value = +raw;
    return isFinite(value) && value > 0 ? value : null;
  }

  determineDataLimits() {
    const {min, max} = this.getMinMax(true);
    this.min = min > 0 ? min : 1;
    this.max = max;
  }

  buildTicks() {
    const ticks = [];
    let power = Math.floor(Math.log2(this.min));
    let value = Math.pow(2, power);

    while (value <= this.max) {
      ticks.push({value});
      power++;
      value = Math.pow(2, power);
    }

    return ticks;
  }

  getPixelForValue(value) {
    if (value === undefined || value <= 0) {
      value = this.min;
    }

    const decimal = (Math.log2(value) - Math.log2(this.min)) /
                    (Math.log2(this.max) - Math.log2(this.min));
    return this.getPixelForDecimal(decimal);
  }

  getValueForPixel(pixel) {
    const decimal = this.getDecimalForPixel(pixel);
    return Math.pow(2, Math.log2(this.min) +
           decimal * (Math.log2(this.max) - Math.log2(this.min)));
  }
}

Log2Axis.id = 'log2';
Log2Axis.defaults = {};

Chart.register(Log2Axis);

// Usage
new Chart(ctx, {
  type: 'line',
  data: data,
  options: {
    scales: {
      y: {
        type: 'log2'
      }
    }
  }
});
```

## Built-in Scales

Available scales to extend:

| Scale | Type | Description |
|-------|------|-------------|
| `LinearScale` | `linear` | Linear numeric scale |
| `LogarithmicScale` | `logarithmic` | Logarithmic (base 10) scale |
| `CategoryScale` | `category` | Discrete category scale |
| `TimeScale` | `time` | Date/time scale |
| `TimeSeriesScale` | `timeseries` | Time series optimized scale |
| `RadialLinearScale` | `radialLinear` | Radial scale for radar/polar charts |
