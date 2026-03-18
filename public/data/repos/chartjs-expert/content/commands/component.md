---
name: component
description: Generate Chart.js chart code interactively
argument-hint: [chart-type]
allowed-tools: Read, Write, Edit, AskUserQuestion
---

# Chart.js Component Generator

Generate production-ready Chart.js v4.5.1 chart code with proper imports and configuration.

## Instructions

Generate a complete Chart.js chart component based on user requirements.

If a chart type is provided as `$ARGUMENTS`, use that type. Otherwise, ask the user what type of chart they want.

### Chart Type Selection

If no chart type specified, ask the user to select from:
- **Line** - Trends over time, continuous data
- **Bar** - Comparing categories, discrete data
- **Pie** - Parts of a whole (no center)
- **Doughnut** - Parts of a whole (with center)
- **Radar** - Multivariate data comparison
- **Polar Area** - Cyclical data patterns
- **Scatter** - Correlation between variables
- **Bubble** - Three-dimensional data

### Framework Selection

Ask which framework/environment to generate code for:
- **Vanilla JS** - Plain JavaScript with CDN or ES modules
- **React** - Using react-chartjs-2
- **Vue 3** - Using vue-chartjs
- **Angular** - Using ng2-charts
- **Rails 8** - Stimulus controller with importmaps

### Code Generation Requirements

Generate complete, working code that includes:

1. **Proper imports** - Tree-shaken imports for production or chart.js/auto for prototyping
2. **Component registration** - Register required Chart.js components
3. **Sample data** - Realistic placeholder data that can be easily replaced
4. **Configuration** - Common options like responsive, legend, title
5. **Type safety** - TypeScript types for Angular/React if appropriate
6. **Comments** - Brief comments explaining customization points

### Framework-Specific Patterns

**Vanilla JS:**
```javascript
import Chart from 'chart.js/auto';
// Or tree-shaken imports for production
```

**React (react-chartjs-2):**
```jsx
import { Bar } from 'react-chartjs-2';
import { Chart as ChartJS, ... } from 'chart.js';
ChartJS.register(...);
```

**Vue 3 (vue-chartjs):**
```vue
<script setup>
import { Bar } from 'vue-chartjs';
import { Chart as ChartJS, ... } from 'chart.js';
ChartJS.register(...);
</script>
```

**Angular (ng2-charts):**
```typescript
import { BaseChartDirective } from 'ng2-charts';
import { ChartConfiguration } from 'chart.js';
```

**Rails 8 (Stimulus):**
```javascript
import { Controller } from "@hotwired/stimulus";
import Chart from "chart.js/auto";
```

### Output Format

Present the generated code in a code block with the appropriate language identifier. Include:

1. The complete component/controller code
2. Example usage in a template/view (if applicable)
3. Brief customization notes

### Follow-up

After generating code, ask if the user wants to:
- Customize the data structure
- Add additional chart options (animations, interactions)
- Add multiple datasets
- Configure axes or scales
- Add a custom plugin
