---
name: ux-ui-design-system
description: |
  Use this skill when the user mentions ANY of: "design", "redesign", "UI", "UX", "styling", "component", "accessibility", "a11y", "WCAG", "contrast", "spacing", "typography", "color", "layout", "responsive", "touch target", "focus", "keyboard", "screen reader", "design system", "shadcn", "refactor UI", "audit", "review design", "improve UI", "fix styling", "make accessible", "design tokens", "visual hierarchy", "button", "input", "form", "card", "modal", "dialog", "navigation", "header", "footer", "sidebar", "table", "list", "grid", "flex", "animation", "transition", "hover", "shadow", "border", "radius", "padding", "margin", "theme", "dark mode", "light mode", "icon", "image", "avatar", "badge", "tooltip", "dropdown", "select", "checkbox", "radio", "switch", "slider", "progress", "skeleton", "loading", "error state", "empty state", "radix", "tailwind", "CSS", "styles", "visual", "interface", "usability", "user experience", "user interface".

  Also use when the user asks to "review", "audit", "improve", "fix", "refactor", "check", "analyze", or "make accessible" any user interface element.

  NOT for: Backend logic, API design, database queries, or non-visual functionality.
---

# Claude Dolphin - UX/UI Design System

Systematic approach to UI/UX quality, accessibility compliance, and design consistency for Claude Code.

## When to Apply

Activate for ANY design-related work:
- Creating or modifying UI components
- Reviewing existing interfaces
- Accessibility compliance checks
- Design system consistency audits
- Visual refinement and polish
- Component styling changes
- Layout modifications

## Core Methodology

### Hierarchy of Concerns (Priority Order)

1. **Accessibility** - WCAG 2.2 AA compliance (non-negotiable)
2. **Consistency** - Design system alignment
3. **Usability** - Intuitive interaction patterns
4. **Aesthetics** - Visual polish and refinement

### Before Any UI Work

1. **Assess the current state** - Read existing code before suggesting changes
2. **Check for design tokens** - Use existing colors, spacing, typography
3. **Verify component library** - Use shadcn/existing components before creating custom
4. **Consider accessibility first** - Every change must maintain or improve a11y

## Systematic Checks

### Accessibility (WCAG 2.2 AA)

**Color Contrast:**
- Body text: 4.5:1 minimum ratio
- Large text (18px+ or 14px+ bold): 3:1 minimum ratio
- UI components and icons: 3:1 minimum ratio
- Focus indicators: 3:1 against adjacent colors

**Interactive Elements:**
- Touch targets: minimum 44×44px (WCAG 2.2)
- Focus states: visible, 2px+ outline, high contrast
- Keyboard navigation: logical tab order, no traps
- Skip links for repetitive content

**Semantic Structure:**
- Proper heading hierarchy (h1 > h2 > h3, no skips)
- Form labels associated with inputs
- Alt text for images (empty for decorative)
- ARIA labels where semantic HTML insufficient

**Motion & Time:**
- Respect `prefers-reduced-motion`
- No auto-playing media
- Provide pause controls for animations

### Design Consistency

**Design Tokens:**
- Colors from defined palette (CSS variables or Tailwind config)
- Spacing from consistent scale (4, 8, 12, 16, 24, 32, 48, 64px)
- Typography from defined scale
- Border radius consistent across components

**Component Patterns:**
- Reuse existing components (shadcn/ui preferred)
- Consistent button styles and sizes
- Uniform form input styling
- Card/container patterns followed

### Refactoring UI Principles

**Spacing:**
- Start with MORE whitespace than feels necessary
- Use consistent scale (base 4px or 8px)
- Group related elements, separate unrelated
- Padding inside, margin outside

**Typography:**
- Limit to 2-3 font sizes per component
- Create hierarchy through weight FIRST, then size
- Line height: 1.5 for body, 1.2-1.3 for headings
- Avoid pure black (#000); use dark gray (#111 or similar)

**Color:**
- Use fewer colors, more intentionally
- Saturated colors for interactive elements only
- Neutral colors for most content
- Create hierarchy with opacity/saturation, not hue

**Borders & Shadows:**
- Prefer spacing over borders for separation
- Use background contrast instead of borders when possible
- Shadows: offset toward light source (slight up-left)
- Layer shadows for depth (ambient + direct)

**Visual Hierarchy:**
- One primary action per view (most prominent)
- Secondary actions visually recede
- Use size + weight + color together
- Lead the eye with alignment and whitespace

## Available Commands

- `/ui-audit [scope]` - Comprehensive UI audit (accessibility, consistency, usability)
- `/accessibility-check [component]` - WCAG 2.2 AA compliance check
- `/design-system-review` - Design token and component consistency analysis
- `/refactor-ui [component]` - Apply tactical UI improvements

## shadcn MCP Integration

When shadcn MCP is available, use these tools:

- `search_items_in_registries` - Find existing components before creating custom
- `view_items_in_registries` - Inspect component implementation details
- `get_item_examples_from_registries` - Find usage patterns and demos
- `get_add_command_for_items` - Get install commands for components

**Before creating any UI component:**
1. Search shadcn registry for existing solution
2. Check if installed components cover the use case
3. Only create custom if no suitable component exists

## Reference Materials

Detailed guidance available in references directory:
- `references/wcag-checklist.md` - Complete WCAG 2.2 AA checklist
- `references/refactoring-ui.md` - Tactical design principles
- `references/shadcn-patterns.md` - Component best practices
- `references/audit-framework.md` - 5-step audit methodology

## Quick Reference

### Contrast Ratios
| Element | Minimum Ratio |
|---------|---------------|
| Body text | 4.5:1 |
| Large text (18px+) | 3:1 |
| Bold text (14px+) | 3:1 |
| UI components | 3:1 |
| Focus indicators | 3:1 |

### Spacing Scale
```
4px  - tight spacing, icons
8px  - small gaps, inline elements
12px - component internal padding
16px - standard gap, card padding
24px - section separation
32px - major section gaps
48px - page section separation
64px - hero/major divisions
```

### Typography Scale
```
12px - captions, labels
14px - secondary text, metadata
16px - body text (minimum for readability)
18px - large body, emphasis
20px - subheadings
24px - section headings
30px - page headings
36px - hero headings
48px - display text
```

## Anti-Patterns to Avoid

- Hardcoded color values (use tokens)
- Inconsistent spacing (use scale)
- Missing focus states
- Color-only information conveyance
- Touch targets smaller than 44×44px
- Skipped heading levels
- Missing alt text on meaningful images
- Auto-playing animations without pause
- Forms without proper labels
- Creating custom components when shadcn exists
