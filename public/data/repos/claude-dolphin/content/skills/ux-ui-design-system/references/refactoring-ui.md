# Refactoring UI - Tactical Design Principles

Practical design principles for building better user interfaces, based on the Refactoring UI methodology.

## Visual Hierarchy

### The Problem
Without clear hierarchy, everything competes for attention. Users don't know where to look or what's important.

### Primary, Secondary, Tertiary

Every interface needs clear levels of emphasis:

**Primary:** The main action or information
- Largest, boldest, most colorful
- One primary element per view
- Example: Main CTA button, page title

**Secondary:** Supporting information
- Smaller, lighter, less saturated
- Can have multiple secondary elements
- Example: Secondary buttons, subtitles

**Tertiary:** De-emphasized details
- Smallest, lightest, most muted
- Metadata, timestamps, labels
- Example: "Posted 2 hours ago", field hints

### Hierarchy Through Typography

**Don't rely on size alone.** Use:

1. **Font weight** - Bold for emphasis, light for de-emphasis
2. **Color** - Dark for important, lighter gray for secondary
3. **Size** - Vary, but not as primary differentiator

```css
/* Good hierarchy */
.title { font-size: 24px; font-weight: 700; color: #111; }
.subtitle { font-size: 16px; font-weight: 500; color: #444; }
.meta { font-size: 14px; font-weight: 400; color: #666; }

/* Bad - size only */
.title { font-size: 32px; }
.subtitle { font-size: 20px; }
.meta { font-size: 14px; }
```

### Semantic vs Visual Hierarchy

**Semantic:** What HTML tag (h1, h2, p) - for accessibility
**Visual:** How it looks - for user experience

These don't have to match! A card title might be visually prominent but semantically an h3.

## Spacing & Layout

### The Spacing System

**Pick a base unit** (4px or 8px recommended) and derive all spacing from it:

```
4px   - Tiny gaps, icon padding
8px   - Small gaps, tight grouping
12px  - Component internal spacing
16px  - Default gap, standard padding
24px  - Section gaps, card padding
32px  - Major section separation
48px  - Page section dividers
64px  - Hero sections, major breaks
96px  - Full-page section separation
```

### Spacing Principles

**Start with too much whitespace**, then reduce if needed. It's easier to identify "too much space" than "not enough."

**Use space to create relationships:**
- Close together = related
- Far apart = unrelated
- Consistent spacing = professional

**The "squint test":** Squint at your UI. Related items should blur together; sections should remain distinct.

### Dense vs Spacious

**Dense interfaces** (dashboards, data tables):
- Tighter spacing
- More information visible
- For expert users who need efficiency

**Spacious interfaces** (marketing pages, onboarding):
- Generous whitespace
- Easier scanning
- For casual users or first impressions

## Typography

### The Type Scale

Don't calculate - pick a hand-crafted scale:

```
12px  - Captions, badges, labels
14px  - Secondary text, table data, metadata
16px  - Body text (never smaller for readability)
18px  - Large body, emphasis
20px  - Subheadings, card titles
24px  - Section headings
30px  - Page headings
36px  - Hero subheadings
48px  - Hero headings, display
60px  - Marketing headlines
```

**Limit sizes:** Use 2-3 sizes per component maximum.

### Font Weight

Weight creates hierarchy more subtly than size:

- **400 (Regular):** Body text, most content
- **500 (Medium):** Slight emphasis, labels
- **600 (Semibold):** Subheadings, buttons
- **700 (Bold):** Headings, primary emphasis

**Avoid extremes:** 300 (light) and 800+ (black) are hard to read.

### Line Height

- **Body text:** 1.5-1.75 (spacious, readable)
- **Headings:** 1.2-1.3 (tighter, looks better)
- **UI elements:** 1.0-1.2 (buttons, badges)

### Letter Spacing

- **Body text:** Normal (0)
- **All caps:** +0.05em to +0.1em (improves readability)
- **Large headings:** -0.02em to -0.01em (looks tighter)

## Color

### Building a Palette

**Grays (9+ shades):**
```
gray-50:  #fafafa  - Subtle backgrounds
gray-100: #f5f5f5  - Card backgrounds
gray-200: #e5e5e5  - Borders, dividers
gray-300: #d4d4d4  - Disabled states
gray-400: #a3a3a3  - Placeholder text
gray-500: #737373  - Secondary text
gray-600: #525252  - Body text
gray-700: #404040  - Emphasis text
gray-800: #262626  - Headings
gray-900: #171717  - Primary text
```

**Primary color (5-7 shades):**
```
primary-50:  Subtle backgrounds
primary-100: Hover backgrounds
primary-200: Borders, outlines
primary-300: Secondary icons
primary-400: -
primary-500: Primary buttons, links
primary-600: Hover states
primary-700: Active states
```

### Color Usage Principles

**Reserve saturated colors for actions:**
- Primary buttons
- Links
- Active states
- Error/success indicators

**Use neutrals for content:**
- Text
- Backgrounds
- Borders
- Most UI chrome

**Don't use pure gray for text.** Add a slight tint:
- Warm: Mix with orange/yellow
- Cool: Mix with blue

### Accessible Color Combinations

Always verify contrast ratios:

| Background | Text | Minimum Ratio |
|------------|------|---------------|
| White (#fff) | gray-600+ | 4.5:1 for body |
| gray-100 | gray-700+ | 4.5:1 for body |
| Primary-500 | White | 4.5:1 for body |

## Borders & Shadows

### Borders

**Use borders sparingly.** Alternatives:
- Background color contrast
- Box shadows
- Increased spacing

**When borders are needed:**
- 1px for subtle
- Light colors (gray-200 or lighter)
- Don't compete with content

### Shadows

**Shadow creates elevation.** Higher elevation = larger, more diffuse shadow.

**Two-layer shadows** look more realistic:
```css
/* Ambient light (soft, large) */
box-shadow: 0 10px 25px -5px rgba(0,0,0,0.05);

/* Direct light (sharp, small) */
box-shadow: 0 4px 6px -2px rgba(0,0,0,0.1);

/* Combined */
box-shadow:
  0 10px 25px -5px rgba(0,0,0,0.05),
  0 4px 6px -2px rgba(0,0,0,0.1);
```

**Shadow direction:** Consistent top-left light source (shadows go down-right).

### Elevation Scale

```css
/* Level 0: Flat */
box-shadow: none;

/* Level 1: Raised (cards) */
box-shadow: 0 1px 3px rgba(0,0,0,0.12), 0 1px 2px rgba(0,0,0,0.08);

/* Level 2: Floating (dropdowns) */
box-shadow: 0 4px 6px rgba(0,0,0,0.1), 0 2px 4px rgba(0,0,0,0.06);

/* Level 3: Overlay (modals) */
box-shadow: 0 10px 20px rgba(0,0,0,0.15), 0 3px 6px rgba(0,0,0,0.1);

/* Level 4: Dialog */
box-shadow: 0 25px 50px rgba(0,0,0,0.25);
```

## Images & Icons

### Images

**Overlay text on images:**
- Add dark overlay (rgba(0,0,0,0.4))
- Or add gradient overlay
- Ensure text contrast

**Constrain proportions:** Set max-width/height, use object-fit.

### Icons

**Consistency is key:**
- Same style (outline/filled/duotone)
- Same stroke width
- Same visual weight
- From same icon set

**Icon sizing:**
- 16px: Inline with text
- 20px: In buttons
- 24px: Standalone
- 32px+: Feature icons

**Icon colors:**
- Match text color hierarchy
- Primary color for actions
- Gray for decorative

## Common Patterns

### Buttons

**Primary button:**
- Filled background
- High contrast text
- One per view

**Secondary button:**
- Outline or subtle background
- Lower contrast
- Multiple allowed

**Tertiary button:**
- Text only or very subtle
- Lowest emphasis
- For less important actions

### Cards

**Card spacing:**
- 16-24px internal padding
- 16-32px gap between cards
- Consistent border-radius (4px, 8px, or 12px)

**Card hierarchy:**
- Title (20px, semibold)
- Description (14-16px, regular, gray-600)
- Meta (12-14px, regular, gray-500)
- Action (button or link)

### Forms

**Form spacing:**
- 24px between form groups
- 8px between label and input
- 4px between input and help text

**Form labels:**
- Always visible (not just placeholders)
- Above input (not inline)
- Semibold or medium weight

## Quick Checklist

- [ ] Clear visual hierarchy (primary/secondary/tertiary)
- [ ] Consistent spacing from defined scale
- [ ] Typography limited to 2-3 sizes per component
- [ ] Colors from defined palette (no hardcoded values)
- [ ] Shadows consistent with elevation scale
- [ ] One primary action per view
- [ ] Text contrast meets accessibility requirements
- [ ] Interactive elements have visible hover/focus states
