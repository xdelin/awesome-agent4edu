# Claude Dolphin

<p align="center">
  <img src="https://img.shields.io/badge/Claude_Code-Plugin-blueviolet" alt="Claude Code Plugin">
  <img src="https://img.shields.io/badge/License-MIT-green" alt="MIT License">
  <img src="https://img.shields.io/badge/WCAG-2.2_AA-blue" alt="WCAG 2.2 AA">
</p>

**Comprehensive UX/UI design system skill for Claude Code** - WCAG 2.2 AA accessibility auditing, design consistency checks, Refactoring UI principles, and shadcn/ui integration.

```
   ___  ___  _    ____  _   _ ___ _   _
  |   \/ _ \| |  |  _ \| | | |_ _| \ | |
  | |) | (_) | |__| |_) | |_| || ||  \| |
  |___/ \___/|____|____/ \___/|___|_|\_|
```

## Features

- **WCAG 2.2 AA Compliance** - Complete accessibility auditing with actionable fixes
- **Design System Consistency** - Audit design token usage and component patterns
- **Refactoring UI Principles** - Apply tactical visual improvements
- **shadcn/ui Integration** - MCP-powered component discovery and best practices
- **5-Step Audit Framework** - Systematic UX evaluation methodology

## Table of Contents

- [Installation](#installation)
- [Quick Start](#quick-start)
- [Available Commands](#available-commands)
- [Reference Materials](#reference-materials)
- [When to Use](#when-to-use)
- [Core Methodology](#core-methodology)
- [Quick Reference](#quick-reference)
- [Contributing](#contributing)
- [License](#license)

## Installation

### Automated Installation (via Claude Code)

In a Claude Code session, simply ask:

```
Install the claude-dolphin plugin from https://github.com/nyldn/claude-dolphin
```

### Manual Installation

```bash
# Clone to Claude plugins directory
git clone https://github.com/nyldn/claude-dolphin.git ~/.claude/plugins/claude-dolphin
```

### Clone Anywhere and Symlink

```bash
# Clone to your preferred location
git clone https://github.com/nyldn/claude-dolphin.git ~/git/claude-dolphin

# Create symlink in Claude plugins
ln -s ~/git/claude-dolphin ~/.claude/plugins/claude-dolphin
```

## Quick Start

Once installed, the skill activates automatically for design-related work. Use the slash commands for specific auditing tasks:

```bash
# Comprehensive UI audit
/ui-audit src/components/

# WCAG 2.2 AA accessibility check
/accessibility-check src/components/ui/button.tsx

# Design token and consistency review
/design-system-review

# Apply tactical UI improvements
/refactor-ui src/components/landing/hero.tsx
```

## Available Commands

| Command | Description |
|---------|-------------|
| `/ui-audit [scope]` | Comprehensive UI audit (accessibility, consistency, usability) |
| `/accessibility-check [component]` | WCAG 2.2 AA compliance check |
| `/design-system-review` | Design token and component consistency analysis |
| `/refactor-ui [component]` | Apply tactical UI improvements |

### Command Scope Options

- **File path**: `src/components/ui/button.tsx`
- **Component name**: `Button`
- **Page path**: `src/app/dashboard/page.tsx`
- **System-wide**: `system` for full codebase audit

## Reference Materials

The plugin includes comprehensive reference guides:

| Reference | Description |
|-----------|-------------|
| `wcag-checklist.md` | Complete WCAG 2.2 AA checklist |
| `refactoring-ui.md` | Tactical design principles |
| `shadcn-patterns.md` | Component best practices |
| `audit-framework.md` | 5-step audit methodology |

## When to Use

The skill activates automatically for ANY design-related work:

- Creating or modifying UI components
- Reviewing existing interfaces
- Accessibility compliance checks
- Design system consistency audits
- Visual refinement and polish
- Component styling changes
- Layout modifications

**Trigger keywords**: design, redesign, UI, UX, styling, component, accessibility, a11y, WCAG, contrast, spacing, typography, color, layout, responsive, touch target, focus, keyboard, screen reader, design system, shadcn, refactor UI, audit, review design, improve UI, fix styling, make accessible

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

## Quick Reference

### Contrast Ratios (WCAG 2.2 AA)

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

### Touch Target Sizes

| Context | Minimum Size |
|---------|-------------|
| WCAG 2.2 requirement | 24×24px |
| Recommended for touch | 44×44px |
| Apple HIG | 44×44pt |
| Material Design | 48×48dp |

## Anti-Patterns to Avoid

- ❌ Hardcoded color values (use tokens)
- ❌ Inconsistent spacing (use scale)
- ❌ Missing focus states
- ❌ Color-only information conveyance
- ❌ Touch targets smaller than 44×44px
- ❌ Skipped heading levels
- ❌ Missing alt text on meaningful images
- ❌ Auto-playing animations without pause
- ❌ Forms without proper labels
- ❌ Creating custom components when shadcn exists

## shadcn MCP Integration

When shadcn MCP is available, use these tools:

| Tool | Purpose |
|------|---------|
| `search_items_in_registries` | Find existing components before creating custom |
| `view_items_in_registries` | Inspect component implementation details |
| `get_item_examples_from_registries` | Find usage patterns and demos |
| `get_add_command_for_items` | Get install commands for components |

**Before creating any UI component:**
1. Search shadcn registry for existing solution
2. Check if installed components cover the use case
3. Only create custom if no suitable component exists

## Project Structure

```
claude-dolphin/
├── plugin.json                 # Plugin manifest
├── LICENSE                     # MIT License
├── README.md                   # This file
├── commands/                   # Slash commands
│   ├── ui-audit.md
│   ├── accessibility-check.md
│   ├── design-system-review.md
│   └── refactor-ui.md
└── skills/
    └── ux-ui-design-system/
        ├── SKILL.md            # Main skill file
        └── references/         # Reference materials
            ├── wcag-checklist.md
            ├── refactoring-ui.md
            ├── shadcn-patterns.md
            └── audit-framework.md
```

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add some amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

### Areas for Contribution

- Additional WCAG criteria documentation
- Framework-specific patterns (Vue, Svelte, etc.)
- Additional design system integrations
- Automated testing helpers

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- [WCAG 2.2](https://www.w3.org/TR/WCAG22/) - Web Content Accessibility Guidelines
- [Refactoring UI](https://refactoringui.com/) - Tactical design advice
- [shadcn/ui](https://ui.shadcn.com/) - Re-usable components
- [Claude Code](https://claude.ai/code) - Anthropic's CLI for Claude

---

<p align="center">
  Made with care by <a href="https://github.com/nyldn">nyldn</a>
</p>
