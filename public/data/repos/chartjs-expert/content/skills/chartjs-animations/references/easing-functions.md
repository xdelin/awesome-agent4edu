# Easing Functions Visual Guide

Complete reference for Chart.js easing functions with descriptions and use cases.

## Understanding Easing

Easing functions control the rate of change during an animation:
- **easeIn**: Starts slow, accelerates toward the end
- **easeOut**: Starts fast, decelerates toward the end (most natural)
- **easeInOut**: Starts slow, accelerates, then decelerates

## Linear

**Function**: `'linear'`

- Constant speed throughout animation
- No acceleration or deceleration
- Mathematical: `y = x`

**Use Cases**:
- Loading indicators
- Progress bars
- Mechanical/robotic motion
- Looping animations (tension, rotation)

## Quadratic (Quad)

Mathematical formula: `y = x²`

| Function | Pattern | Use Case |
|----------|---------|----------|
| `easeInQuad` | Slow → Fast | Accelerating into view |
| `easeOutQuad` | Fast → Slow | Natural deceleration |
| `easeInOutQuad` | Slow → Fast → Slow | Smooth both ends |

**Characteristics**: Gentle curve, subtle acceleration.

## Cubic

Mathematical formula: `y = x³`

| Function | Pattern | Use Case |
|----------|---------|----------|
| `easeInCubic` | Slow → Fast | Stronger acceleration |
| `easeOutCubic` | Fast → Slow | Smooth landing |
| `easeInOutCubic` | Slow → Fast → Slow | Balanced, professional |

**Characteristics**: More pronounced than Quad, still feels natural.

## Quartic (Quart)

Mathematical formula: `y = x⁴`

| Function | Pattern | Use Case |
|----------|---------|----------|
| `easeInQuart` | Very Slow → Very Fast | Dramatic entrance |
| `easeOutQuart` | Very Fast → Very Slow | **Default** - most natural |
| `easeInOutQuart` | Very Slow → Very Fast → Very Slow | Emphasized motion |

**Characteristics**: Strong acceleration, feels powerful yet smooth.

## Quintic (Quint)

Mathematical formula: `y = x⁵`

| Function | Pattern | Use Case |
|----------|---------|----------|
| `easeInQuint` | Extremely Slow → Extremely Fast | Maximum drama |
| `easeOutQuint` | Extremely Fast → Extremely Slow | Heavy object settling |
| `easeInOutQuint` | Extreme curve | Exaggerated motion |

**Characteristics**: Very strong acceleration, can feel exaggerated.

## Sinusoidal (Sine)

Mathematical formula: Based on sine wave

| Function | Pattern | Use Case |
|----------|---------|----------|
| `easeInSine` | Gentle start | Soft entrance |
| `easeOutSine` | Gentle end | Soft exit |
| `easeInOutSine` | Gentle both | Smooth, wave-like |

**Characteristics**: Smooth, gentle, wave-like motion. Very natural.

## Exponential (Expo)

Mathematical formula: Exponential function

| Function | Pattern | Use Case |
|----------|---------|----------|
| `easeInExpo` | Very slow start, explosive end | Rocket launch |
| `easeOutExpo` | Explosive start, quick settle | Sudden stop |
| `easeInOutExpo` | Explosive middle | Dramatic transitions |

**Characteristics**: Extreme acceleration, feels abrupt.

## Circular (Circ)

Mathematical formula: Based on circular arc

| Function | Pattern | Use Case |
|----------|---------|----------|
| `easeInCirc` | Slow → Aggressive end | Sharp acceleration |
| `easeOutCirc` | Aggressive start → Slow | Sharp deceleration |
| `easeInOutCirc` | Sharp both ends | Angular motion |

**Characteristics**: Sharp changes, geometric feel.

## Elastic

Mathematical formula: Oscillating decay

| Function | Pattern | Use Case |
|----------|---------|----------|
| `easeInElastic` | Pulls back, then springs forward | Spring release |
| `easeOutElastic` | Springs forward, oscillates to rest | Playful bounce |
| `easeInOutElastic` | Oscillates both ends | Springy, attention-grabbing |

**Characteristics**: Overshoots target and oscillates back. Playful, energetic.

**Warning**: Can be distracting or unprofessional in serious contexts.

## Back

Mathematical formula: Overshoots then returns

| Function | Pattern | Use Case |
|----------|---------|----------|
| `easeInBack` | Pulls back before moving forward | Wind-up motion |
| `easeOutBack` | Overshoots, then settles | Elastic snap |
| `easeInOutBack` | Pulls back, overshoots, settles | Emphasized motion |

**Characteristics**: Slight overshoot, feels anticipatory or playful.

## Bounce

Mathematical formula: Bouncing ball physics

| Function | Pattern | Use Case |
|----------|---------|----------|
| `easeInBounce` | Bounces before arriving | Bouncing into place |
| `easeOutBounce` | Bounces after arriving | **Popular** - bouncing ball |
| `easeInOutBounce` | Bounces at both ends | Double bounce |

**Characteristics**: Simulates bouncing physics. Fun, playful, attention-grabbing.

## Easing Function Selection Guide

### Professional / Corporate

```javascript
// Recommended for dashboards, analytics, business charts
easing: 'easeOutQuart'    // Default, most natural
easing: 'easeInOutCubic'  // Smooth, balanced
easing: 'easeOutCubic'    // Clean deceleration
```

### Playful / Consumer Apps

```javascript
// Recommended for games, consumer apps, creative projects
easing: 'easeOutBounce'   // Fun bounce effect
easing: 'easeOutElastic'  // Spring effect
easing: 'easeOutBack'     // Slight overshoot
```

### Data Visualization

```javascript
// Recommended for data updates, real-time charts
easing: 'easeOutQuart'    // Natural feel
easing: 'linear'          // Constant speed
easing: 'easeInOutQuad'   // Gentle acceleration
```

### Emphasis / Attention

```javascript
// Recommended for important updates, alerts
easing: 'easeOutBounce'   // Draws attention
easing: 'easeInOutQuart'  // Emphasized motion
easing: 'easeInOutElastic' // Maximum attention
```

### Loading / Progress

```javascript
// Recommended for loaders, progress indicators
easing: 'linear'          // Constant progress
easing: 'easeOutSine'     // Gentle completion
```

## Combining Easing with Properties

Different properties can have different easing:

```javascript
animations: {
  x: {
    duration: 1000,
    easing: 'easeOutQuart'  // Natural horizontal motion
  },
  y: {
    duration: 1000,
    easing: 'easeOutBounce'  // Bouncy vertical motion
  },
  colors: {
    duration: 500,
    easing: 'linear'  // Smooth color transition
  }
}
```

## Duration Recommendations by Easing

| Easing Type | Recommended Duration | Reason |
|-------------|---------------------|--------|
| Linear | 500-1000ms | Constant speed doesn't need to be long |
| Quad/Cubic/Quart | 1000-1500ms | Natural timing for smooth curves |
| Elastic/Bounce | 1500-2500ms | Need time for oscillations |
| Expo | 800-1200ms | Fast enough to feel snappy |
| Sine | 1000-1500ms | Gentle waves benefit from time |

## Testing Easing Functions

Quick test to visualize all easing functions:

```javascript
const easings = [
  'linear', 'easeInQuad', 'easeOutQuad', 'easeInOutQuad',
  'easeInCubic', 'easeOutCubic', 'easeInOutCubic',
  'easeInQuart', 'easeOutQuart', 'easeInOutQuart',
  'easeInElastic', 'easeOutElastic', 'easeInOutElastic',
  'easeOutBounce'
];

easings.forEach(easing => {
  console.log(`Testing: ${easing}`);
  chart.options.animation.easing = easing;
  chart.update();
});
```

## Common Mistakes

### 1. Using Elastic/Bounce for Professional Dashboards

```javascript
// ❌ Too playful for business context
easing: 'easeOutBounce'

// ✅ Professional and natural
easing: 'easeOutQuart'
```

### 2. Long Duration with Fast Easing

```javascript
// ❌ Feels sluggish
duration: 3000,
easing: 'easeOutExpo'  // Settles quickly, but waits 3s

// ✅ Matched timing
duration: 1000,
easing: 'easeOutExpo'
```

### 3. Elastic on Every Update

```javascript
// ❌ Distracting on frequent updates
easing: 'easeOutElastic'  // With live data updating every second

// ✅ Subtle on frequent updates
easing: 'easeOutQuad'
```

## Accessibility Considerations

Users with vestibular disorders may be sensitive to motion:

```javascript
// Check user preferences
const prefersReducedMotion = window.matchMedia('(prefers-reduced-motion: reduce)').matches;

const chart = new Chart(ctx, {
  options: {
    animation: prefersReducedMotion ? {
      duration: 0  // Disable animations
    } : {
      duration: 1000,
      easing: 'easeOutQuart'
    }
  }
});
```

Avoid these for accessibility-sensitive contexts:
- ❌ `easeInOutElastic` - oscillates too much
- ❌ `easeOutBounce` - too much movement
- ✅ `easeOutQuart` - smooth and gentle
- ✅ `linear` - predictable motion
