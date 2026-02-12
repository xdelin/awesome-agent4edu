# Chart.js Expert Plugin for Claude Code

[![Markdown Lint](https://github.com/sjnims/chartjs-expert/actions/workflows/markdownlint.yml/badge.svg)](https://github.com/sjnims/chartjs-expert/actions/workflows/markdownlint.yml)
[![HTML Lint](https://github.com/sjnims/chartjs-expert/actions/workflows/html-lint.yml/badge.svg)](https://github.com/sjnims/chartjs-expert/actions/workflows/html-lint.yml)
[![YAML Lint](https://github.com/sjnims/chartjs-expert/actions/workflows/yaml-lint.yml/badge.svg)](https://github.com/sjnims/chartjs-expert/actions/workflows/yaml-lint.yml)
[![Links](https://github.com/sjnims/chartjs-expert/actions/workflows/links.yml/badge.svg)](https://github.com/sjnims/chartjs-expert/actions/workflows/links.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Chart.js](https://img.shields.io/badge/Chart.js-v4.5.1-ff6384.svg)](https://www.chartjs.org/)

A comprehensive Claude Code plugin providing expert-level assistance for Chart.js v4.5.1.

## Features

- **11 Specialized Skills**: Deep knowledge across all Chart.js domains
- **Interactive Command**: `/chartjs:component` - Generate chart code for any framework
- **Proactive Agent**: Automatically assists when working with charts or data visualization

## Skills

| Skill | Coverage |
|-------|----------|
| `chartjs-overview` | Installation, setup, tree-shaking, global configuration |
| `chartjs-chart-types` | Line, bar, pie, doughnut, radar, polar, bubble, scatter, mixed charts |
| `chartjs-configuration` | Options, responsive design, interactions, layout, fonts, colors |
| `chartjs-axes` | Cartesian/radial axes, time scales, logarithmic scales, ticks, grid lines |
| `chartjs-animations` | Animation options, easing functions, callbacks, transitions |
| `chartjs-tooltips` | Tooltip configuration, callbacks, positioning, external HTML tooltips |
| `chartjs-plugins` | Plugin architecture, lifecycle hooks, creating custom plugins |
| `chartjs-integrations` | React, Vue, Angular, Rails 8 (Hotwire/Stimulus), vanilla JS |
| `chartjs-developers` | Chart instance API, extending Chart.js, TypeScript support |
| `chartjs-advanced` | Gradients, custom chart types, custom scales, programmatic control |
| `chartjs-accessibility` | WCAG compliance, ARIA labels, colorblind palettes, keyboard navigation |

## Command

### `/chartjs:component`

Interactive chart code generator. Prompts for:

- Chart type (line, bar, pie, etc.)
- Target framework (React, Vue, Angular, Rails 8, vanilla JS)
- Generates production-ready code with proper imports and tree-shaking

## Agent

### `chartjs-expert`

Triggers proactively when:

- Working with Chart.js code or configuration
- Implementing data visualization features
- Troubleshooting chart rendering issues
- Asking about charting best practices

## Installation

Clone and reference the plugin directory:

```bash
git clone https://github.com/sjnims/chartjs-expert.git
claude --plugin-dir /path/to/chartjs-expert
```

Or add to your Claude Code settings to load automatically.

## Contributing

Contributions are welcome! See [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines
on adding skills, commands, examples, and documentation.

## Chart.js Version

This plugin targets **Chart.js v4.5.1**. Examples use tree-shaking patterns
for production code (manual component registration rather than `chart.js/auto`).

## License

[MIT](LICENSE)
