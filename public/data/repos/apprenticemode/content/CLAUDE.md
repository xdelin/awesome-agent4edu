# apprenticemode

Claude Code plugin marketplace for learning and growth.

## Structure

```
apprenticemode/
├── .claude-plugin/
│   └── marketplace.json    # Marketplace manifest
├── plugins/
│   └── claudebrief/        # Knowledge gap detection plugin
│       ├── .claude-plugin/
│       │   └── plugin.json
│       ├── commands/
│       │   └── reflect.md
│       └── CLAUDE.md
└── README.md
```

## Plugins

### claudebrief v1.0.0

Surface knowledge gaps from your Claude Code sessions. Detects "dangerous verbs" that hide required decisions.

**Commands:**
- `/claudebrief:reflect` — Analyze recent sessions
- `/claudebrief:reflect current` — Analyze this session

**Dangerous Verbs:**
- retry, backoff → What's the retry policy?
- cache, memoize → What's the invalidation policy?
- optimize, faster → What's the performance target?
- prod-ready, ship → What's the rollback plan?
- refactor, clean up → What must stay the same?
- async, parallel → What's the ordering guarantee?

## Development

```bash
# Install from local marketplace
claude plugin marketplace add /path/to/apprenticemode/.claude-plugin/marketplace.json
claude plugin install claudebrief@apprenticemode
```
