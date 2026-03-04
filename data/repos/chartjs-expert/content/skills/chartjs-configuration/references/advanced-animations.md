# Advanced Animations

Detailed documentation for Chart.js animation configuration including transitions, callbacks, and the animation object.

## Animation Configuration Hierarchy

Chart.js animation is configured at three levels:

```javascript
options: {
  animation: { },     // Base settings for all animations
  animations: { },    // Per-property animation configuration
  transitions: { }    // Mode-specific animation overrides
}
```

## Animation Callbacks

The `onProgress` and `onComplete` callbacks receive an animation object:

```javascript
{
  chart: Chart,           // The chart instance
  currentStep: number,    // Current animation frame
  initial: boolean,       // True for initial chart render
  numSteps: number        // Total animation frames
}
```

### Progress Bar Example

```javascript
const chart = new Chart(ctx, {
  type: 'line',
  data: data,
  options: {
    animation: {
      onProgress: function(animation) {
        const percent = animation.currentStep / animation.numSteps;
        progressBar.style.width = (percent * 100) + '%';
      },
      onComplete: function(animation) {
        progressBar.style.width = '100%';
      }
    }
  }
});
```

## Transitions

Transitions control animations for specific chart modes. Core transitions:

| Mode | When Triggered |
|------|----------------|
| `'active'` | Element is hovered |
| `'hide'` | Dataset hidden via legend or API |
| `'show'` | Dataset shown via legend or API |
| `'reset'` | Chart reset |
| `'resize'` | Chart container resized |

### Default Transition Settings

| Mode | Setting | Value | Purpose |
|------|---------|-------|---------|
| `active` | `animation.duration` | `400` | Faster hover response |
| `resize` | `animation.duration` | `0` | Instant resize |
| `show` | `animations.colors.from` | `'transparent'` | Fade in colors |
| `show` | `animations.visible.duration` | `0` | Immediate visibility |
| `hide` | `animations.colors.to` | `'transparent'` | Fade out colors |
| `hide` | `animations.visible.easing` | `'easeInExpo'` | Late visibility change |

### Custom Transition Example

Animate from zero when showing datasets:

```javascript
options: {
  transitions: {
    show: {
      animations: {
        x: { from: 0 },
        y: { from: 0 }
      }
    },
    hide: {
      animations: {
        x: { to: 0 },
        y: { to: 0 }
      }
    }
  }
}
```

### Custom Transition Mode

Create custom transitions via the API:

```javascript
// Define a custom transition
chart.options.transitions.myCustomMode = {
  animation: {
    duration: 2000,
    easing: 'easeInOutBounce'
  }
};

// Use it when updating
chart.update('myCustomMode');
```

## Per-Property Animations

Configure animations for specific properties:

```javascript
animations: {
  // Animate numeric properties
  radius: {
    duration: 400,
    easing: 'linear',
    from: 0,
    to: 5
  },

  // Animate colors
  backgroundColor: {
    type: 'color',
    duration: 1000,
    from: 'rgba(255, 0, 0, 0)',
    to: 'rgba(255, 0, 0, 1)'
  },

  // Animate multiple properties together
  myAnimation: {
    properties: ['x', 'y', 'radius'],
    duration: 500
  }
}
```

### Animation Property Options

| Property | Type | Description |
|----------|------|-------------|
| `properties` | `string[]` | Properties this animation applies to |
| `type` | `string` | `'number'`, `'color'`, or `'boolean'` |
| `from` | `any` | Start value |
| `to` | `any` | End value |
| `fn` | `function` | Custom interpolation function |
| `duration` | `number` | Animation duration (ms) |
| `easing` | `string` | Easing function |
| `delay` | `number` | Delay before starting |
| `loop` | `boolean` | Loop endlessly |

### Default Animations

Chart.js defines these by default (overridden by most controllers):

```javascript
animations: {
  numbers: {
    properties: ['x', 'y', 'borderWidth', 'radius', 'tension'],
    type: 'number'
  },
  colors: {
    properties: ['color', 'borderColor', 'backgroundColor'],
    type: 'color'
  }
}
```

### Looping Animation Example

Create a continuously animating tension:

```javascript
options: {
  animations: {
    tension: {
      duration: 1000,
      easing: 'linear',
      from: 1,
      to: 0,
      loop: true
    }
  },
  scales: {
    y: { min: 0, max: 100 }  // Fix scale to prevent jumps
  }
}
```

### Custom Interpolation Function

```javascript
animations: {
  customProp: {
    fn: (from, to, factor) => {
      // factor goes from 0 to 1
      return from + (to - from) * Math.pow(factor, 2);
    }
  }
}
```

## Disabling Animations

```javascript
// Disable all animations
chart.options.animation = false;

// Disable specific property animations
chart.options.animations.colors = false;
chart.options.animations.x = false;

// Disable animation for specific mode
chart.options.transitions.active.animation.duration = 0;

// Disable per-update
chart.update('none');
```

## Scriptable Animation Options

Animation options can be scriptable (functions):

```javascript
animation: {
  delay: (context) => context.dataIndex * 100,  // Stagger by index
  duration: (context) => context.active ? 400 : 1000
}
```

### Context Object

The animation context includes:

```javascript
{
  chart: Chart,
  dataIndex: number,
  datasetIndex: number,
  active: boolean,
  type: string  // 'data', 'dataset', etc.
}
```

## Programmatic Animation Control

```javascript
// Update with specific transition
chart.update('active');
chart.update('resize');
chart.update('none');  // No animation

// Stop all animations
chart.stop();

// Check if animating
const isAnimating = chart.isAnimating();
```
