# Contributing to Chart.js Expert Plugin

Thank you for your interest in contributing! This document provides
guidelines for contributing to this Claude Code plugin.

## Code of Conduct

This project adheres to our [Code of Conduct](CODE_OF_CONDUCT.md).
By participating, you agree to uphold this code.

## Types of Contributions

| Type | Description |
|------|-------------|
| **Skills** | Add new Chart.js expertise (e.g., new chart types, techniques) |
| **Commands** | Enhance `/chartjs:component` or add new commands |
| **Agent** | Improve the proactive agent's triggers or responses |
| **Examples** | Add working HTML/JS examples demonstrating Chart.js patterns |
| **Documentation** | Fix errors, add examples, improve clarity |

## Plugin Structure

```text
.claude-plugin/
└── plugin.json                   # Plugin manifest

agents/
└── chartjs-expert.md             # Proactive agent

commands/
└── component.md                  # /chartjs:component command

skills/                           # Knowledge skills
└── chartjs-*/                    # Each skill has:
    ├── SKILL.md                  # Core knowledge (required)
    ├── references/               # Deep-dive documentation (optional)
    └── examples/                 # Working code samples (optional)
```

See [CLAUDE.md](CLAUDE.md) for detailed component conventions.

## Development Setup

1. Clone the repository
2. Review [CLAUDE.md](CLAUDE.md) for component rules
3. Test locally: `claude --plugin-dir .`

## Component Guidelines

### Skills

- YAML frontmatter requires `name` and `description`
- `description` must start with "This skill should be used when..."
- Body contains core knowledge, not exhaustive documentation
- Put detailed docs in `references/`, examples in `examples/`

### Commands

- YAML frontmatter requires `name`, `description`, `allowed-tools`
- `description` must be 60 characters or fewer
- Use imperative voice ("Do X" not "You should X")

### Agents

- YAML frontmatter requires `name`, `description`, `model`, `color`, `tools`
- Include 2-4 `<example>` blocks with user/assistant dialogue

## Linting

Run these before submitting:

```bash
# Markdown (required)
markdownlint '**/*.md'

# HTML examples
npx htmlhint 'skills/**/examples/*.html'

# YAML files
yamllint -c .yamllint.yml .github/ .claude-plugin/

# Broken links
lychee --cache '**/*.md'
```

## Pull Request Process

1. Fork and create a feature branch (`docs/description`, `feat/description`, etc.)
2. Make changes following the component guidelines above
3. Ensure all linters pass
4. Submit PR using the provided template
5. Address review feedback

## Chart.js Version

This plugin targets **Chart.js v4.5.1**. When adding examples:

- Use tree-shaking patterns for production code
- `chart.js/auto` is acceptable for prototyping only
- Include required components (Controllers, Elements, Scales, Plugins)

## Questions?

Open a [Discussion](https://github.com/sjnims/chartjs-expert/discussions)
for questions or ideas.
