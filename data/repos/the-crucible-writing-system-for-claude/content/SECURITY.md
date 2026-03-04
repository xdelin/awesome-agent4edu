# Security Policy

## Supported Versions

| Version | Supported          |
| ------- | ------------------ |
| 1.0.x   | :white_check_mark: |

## Reporting a Vulnerability

We take security seriously. If you discover a security vulnerability in Crucible Suite, please report it responsibly.

### How to Report

1. **Do NOT open a public GitHub issue** for security vulnerabilities
2. Instead, contact the maintainers directly through GitHub's private vulnerability reporting feature, or reach out via the repository's Discussions with a private message request

### What to Include

When reporting a vulnerability, please include:

- Description of the vulnerability
- Steps to reproduce
- Potential impact
- Suggested fix (if you have one)

### Response Timeline

- **Initial response**: Within 48 hours
- **Status update**: Within 7 days
- **Resolution target**: Within 30 days for critical issues

### Scope

This security policy covers:

- The Crucible Suite plugin code
- Python automation scripts
- Hook configurations
- Any code that executes on user systems

### Out of Scope

- Claude Code itself (report to Anthropic)
- Third-party dependencies (report to respective maintainers)
- Issues in user-generated content (planning documents, manuscripts, etc.)

## Security Best Practices for Users

When using Crucible Suite:

1. **Review hooks before enabling** - Hooks execute shell commands; ensure you trust the source
2. **Keep Python updated** - Use Python 3.8+ with latest security patches
3. **Backup your work** - The plugin creates automatic backups, but maintain your own as well
4. **Don't share state files** - `.crucible/state/` may contain project-specific data

## Acknowledgments

We appreciate responsible disclosure and will acknowledge security researchers who report valid vulnerabilities (with their permission).
