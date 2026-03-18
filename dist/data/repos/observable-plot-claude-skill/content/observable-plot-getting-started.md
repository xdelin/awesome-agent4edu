# Getting Started with the Observable Plot Skill

Welcome! This guide will help you get the most out of Claude and the Observable Plot skill for creating data visualizations.

## What is a Skill?

Think of a skill as giving Claude specialized training in a particular area. Without the Observable Plot skill, Claude knows *about* the library but might occasionally get syntax wrong or suggest outdated approaches. With the skill installed, Claude has a reliable reference guide loaded and ready whenever you're working on visualizations.

## Installing the Skill

1. Go to [claude.ai](https://claude.ai) and sign in
2. Click on your profile icon (bottom-left corner)
3. Select **Settings**
4. Find **Skills** in the menu
5. Click **Add Skill** and upload the `observable-plot.skill` file

That's it! The skill activates automatically whenever your conversation involves Observable Plot.

---

## What Can You Ask Claude?

Here are practical examples organized by the kind of work you probably do regularly. Feel free to copy and modify these prompts!

### Converting Your Ideas to Code

**Start simple and describe what you want:**

> "I have sales data with columns for date, revenue, and region. I want a line chart showing revenue over time, with a separate line for each region."

> "Create a bar chart comparing Q1 vs Q2 performance across our 5 departments."

> "I need a scatter plot showing the relationship between marketing spend and conversions."

**Claude will generate working code you can use directly.**

---

### Working with Data You Have

**Share your data structure and ask for help:**

> "Here's what my data looks like:
> ```
> [
>   {month: "Jan", product: "Widget A", sales: 1200},
>   {month: "Jan", product: "Widget B", sales: 800},
>   {month: "Feb", product: "Widget A", sales: 1400},
>   ...
> ]
> ```
> How do I create a stacked bar chart showing total sales by month, colored by product?"

**Tip:** You can paste sample rows of your actual data (removing sensitive info) and Claude will write code that matches your column names.

---

### Tableau-to-Plot Translation

**If you know how to do something in Tableau, ask Claude to translate:**

> "In Tableau, I'd create a dual-axis chart with bars and a line. How do I do the same thing in Observable Plot?"

> "I want something like Tableau's small multiples / trellis chart - breaking out the same chart by category."

> "How do I get a Tableau-style tooltip that shows details when I hover?"

---

### Debugging and Fixing Issues

**When something isn't working:**

> "This code isn't showing my chart - what's wrong?
> ```javascript
> Plot.plot({
>   marks: [
>     Plot.barY(data, {x: "category", y: "value"})
>   ]
> })
> ```
> My data is: [paste a few rows]"

> "My line chart has gaps where there shouldn't be any. The data has some null values - how do I handle that?"

> "The bars in my chart are in the wrong order. How do I sort them by value instead of alphabetically?"

---

### Styling and Polish

**Making charts look professional:**

> "How do I change the colors in my chart to match our brand? Our primary color is #1B4D7A."

> "The axis labels are getting cut off - how do I add more margin on the left?"

> "Can you add a title and subtitle to this chart?"

> "How do I format the y-axis to show currency (like $1.2M instead of 1200000)?"

---

### Learning and Understanding

**When you want to understand the concepts:**

> "What's the difference between barY and rectY? When should I use each one?"

> "Explain how scales work in Observable Plot - I'm confused about domains and ranges."

> "What are 'transforms' and when would I use them?"

> "Show me a simple example of faceting with an explanation of each part."

---

## Tips for Better Results

### Be Specific About Your Data
Instead of: "Make me a chart"

Try: "I have monthly revenue data with columns `date`, `amount`, and `category`. Create a line chart with the x-axis showing months and y-axis showing revenue, with different colored lines for each category."

### Share Sample Data
Even just 3-5 rows helps Claude write code that matches your actual column names:

```
My data looks like:
{date: "2024-01", region: "North", revenue: 45000}
{date: "2024-01", region: "South", revenue: 38000}
...
```

### Ask Follow-Up Questions
After Claude generates a chart, you can ask things like:
- "Now add a trend line"
- "Change the color scheme to use blues"
- "Can you make the font bigger?"
- "Add gridlines to make it easier to read"

### Iterate on Designs
It's totally fine to go back and forth:
- "Actually, can we try this as a horizontal bar chart instead?"
- "The legend is taking up too much space - can we move it or remove it?"
- "This is close, but can you make it look more minimal?"

---

## Quick Reference: Common Chart Types

| You Want... | Ask For... |
|-------------|------------|
| Trend over time | Line chart (`Plot.line` or `Plot.lineY`) |
| Compare categories | Bar chart (`Plot.barY` or `Plot.barX`) |
| Show distribution | Histogram (`Plot.rectY` with `binX`) |
| Correlation between two measures | Scatter plot (`Plot.dot`) |
| Part-of-whole over time | Stacked area chart (`Plot.areaY` with `stackY`) |
| Side-by-side comparison by group | Faceted chart (small multiples) |
| Geographic data | Map (`Plot.geo`) |
| Matrix/heatmap | Cell chart (`Plot.cell`) |

---

## Getting Help

If Claude's response doesn't quite work or you're stuck:

1. **Share the error message** - If you see an error in your browser console, paste it
2. **Show what you expected** - "I wanted X but got Y instead"
3. **Ask for alternatives** - "Is there another way to do this?"

---

## Helpful Links

- [Observable Plot Documentation](https://observablehq.com/plot/) - Official docs with lots of examples
- [Observable Plot Gallery](https://observablehq.com/@observablehq/plot-gallery) - Visual examples you can browse for inspiration
