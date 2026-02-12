# Observable Plot Skill for Claude

A Claude skill that provides expert guidance for creating data visualizations with Observable Plot, a JavaScript library for exploratory data analysis.

## What It Does

This skill helps Claude assist you with:
- Creating charts and graphs (line charts, bar charts, scatter plots, area charts, histograms, heatmaps, maps)
- Working with `Plot.plot()`, marks, scales, and transforms
- Debugging visualization code
- Styling and customizing charts
- Converting concepts from other tools (like Tableau) to Observable Plot

## Installation

1. Go to [claude.ai](https://claude.ai) and sign in
2. Click your profile icon (bottom-left)
3. Select **Settings** → **Skills**
4. Click **Add Skill** and upload `observable-plot.skill`

The skill activates automatically when working with Observable Plot.

## Usage

Ask Claude natural questions about your visualization needs:

```
"I have sales data with date, revenue, and region columns.
Create a line chart showing revenue over time with a separate line per region."
```

```
"How do I add a trend line to my scatter plot?"
```

```
"My bars are in the wrong order - how do I sort by value?"
```

For more detailed examples and tips, see [observable-plot-getting-started.md](observable-plot-getting-started.md).

## Files

- `observable-plot.skill` - The skill file to install in Claude
- `observable-plot-getting-started.md` - Comprehensive guide with examples and tips

## Resources

- [Observable Plot Documentation](https://observablehq.com/plot/)
- [Observable Plot Gallery](https://observablehq.com/@observablehq/plot-gallery)
