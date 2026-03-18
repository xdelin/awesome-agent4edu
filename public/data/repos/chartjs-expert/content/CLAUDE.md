# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a Claude Code plugin providing Chart.js v4.5.1 expertise with 12 skills, 1 command, and 1 agent.

## Repository Structure

```text
.claude-plugin/
└── plugin.json                   # Plugin manifest

agents/
└── chartjs-expert.md             # Proactive agent for Chart.js work

commands/
└── component.md                  # /chartjs:component command

skills/                           # 11 specialized knowledge skills
└── chartjs-*/                    # Each skill directory contains:
    ├── SKILL.md                  # Core knowledge (required)
    ├── references/               # Deep-dive documentation (optional)
    └── examples/                 # Working code examples (optional)
```

Skill domains: overview, chart-types, configuration, axes, animations, tooltips, plugins, integrations, developers, advanced, accessibility, quickref.

## Linting Commands

```bash
# Markdown
markdownlint '**/*.md'

# HTML (example files)
npx htmlhint 'skills/**/examples/*.html'

# YAML
yamllint -c .yamllint.yml .github/ .claude-plugin/

# Broken links
lychee --cache '**/*.md'

# GitHub Actions
actionlint
```

## Local Development

Test the plugin locally:

```bash
claude --plugin-dir .
```

This loads the plugin for testing skills, commands, and agent behavior.

## Plugin Component Rules

When editing plugin components, follow these conventions:

### Skills (`skills/*/SKILL.md`)

- YAML frontmatter requires `name` and `description`
- `description` must start with "This skill should be used when..." and include trigger phrases
- Body contains core knowledge, not exhaustive documentation
- Use `references/` subdirectory for deep-dive documentation on specific topics
- Use `examples/` subdirectory for working code samples (`.html` or `.md`)

### Commands (`commands/*.md`)

- YAML frontmatter requires `name`, `description`, `allowed-tools`
- Optional: `argument-hint` for showing expected arguments
- `description` must be 60 characters or fewer
- Use imperative voice in body ("Do X" not "You should X")
- `allowed-tools` should be minimal and appropriate for the command's purpose

### Agents (`agents/*.md`)

- YAML frontmatter requires `name`, `description`, `model`, `color`, `tools`
- `description` must include 2-4 `<example>` blocks with Context, user/assistant dialogue, and `<commentary>`
- `model: inherit` uses the parent model
- `tools` is a comma-separated list (e.g., `Read, Write, Edit, Grep, Glob, Bash`)

## Version Management

Version is tracked in `/.claude-plugin/plugin.json`. The project uses semantic versioning and maintains a CHANGELOG.md using Keep a Changelog format.

## CI Workflows

PRs run linting automatically based on changed file types:

- `markdownlint.yml` - Markdown files
- `html-lint.yml` - HTML example files
- `yaml-lint.yml` - YAML config files
- `links.yml` - Broken link checking (weekly + on MD changes)
- `component-validation.yml` - Claude-powered validation of plugin components
- `claude-pr-review.yml` - Claude-powered PR review
- `validate-workflows.yml` - GitHub Actions syntax validation
- `semantic-labeler.yml` - Auto-labeling PRs by path
- `ci-failure-analysis.yml` - Claude analysis of CI failures
- `stale.yml` - Stale issue/PR management
- `greet.yml` - New contributor greeting
- `sync-labels.yml` - Sync labels from .github/labels.yml

## Linting Configuration

| Linter | Config File |
|--------|-------------|
| markdownlint | `.markdownlint.json` |
| htmlhint | `.htmlhintrc` |
| yamllint | `.yamllint.yml` |

## Conventions

**Commits:** Conventional Commits format (`feat:`, `fix:`, `docs:`, `chore:`, etc.)

**Branches:** `docs/description`, `feat/description`, `fix/description`, `chore/description`

## Chart.js Reference

This plugin targets **Chart.js v4.5.1**. When adding examples or updating skills:

- Use tree-shaking patterns for production code (manual component registration)
- `chart.js/auto` is acceptable for prototyping examples only
- Include required components: Controllers, Elements, Scales, Plugins
