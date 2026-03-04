# claudebrief

**Surface knowledge gaps from your Claude Code sessions.**

Detects "dangerous verbs" that hide required decisions, tracks patterns, and compounds learnings into your CLAUDE.md.

---

## The Problem

You say "add retry logic" and Claude builds it. But later you realize:
- You never specified the retry policy
- Claude had to assume exponential backoff
- You didn't think about idempotency

These **dangerous verbs** sound complete but hide required decisions:

| Verb | What's Usually Missing |
|------|------------------------|
| retry, backoff | Retry policy (idempotency, max attempts, failure handling) |
| cache, memoize | Invalidation policy (TTL, consistency, eviction) |
| optimize, faster | Performance target (baseline, target, measurement) |
| prod-ready, ship | Rollback plan (monitoring, alerting, SLOs) |
| refactor, clean up | Invariants (what must NOT change, tests) |
| async, parallel | Ordering guarantee (race conditions, cancellation) |

---

## How It Works

1. **Run** `/claudebrief:reflect` after a session
2. **Scans** for dangerous verbs + evidence (Claude asked/assumed something)
3. **Shows** max 1 gap with evidence and quick rule
4. **Confirms** — only saves if you say yes
5. **Compounds** — on 3rd occurrence, offers CLAUDE.md patch

---

## Install

```bash
# Add the marketplace
claude plugin marketplace add yavzius/apprenticemode

# Install the plugin
claude plugin install claudebrief@apprenticemode
```

---

## Usage

```bash
# Analyze recent sessions
/claudebrief:reflect

# Analyze current session
/claudebrief:reflect current
```

### Example Output

```
───────────────────────────────────────────────────────────────
Gap: RETRIES (HIGH confidence)
───────────────────────────────────────────────────────────────

Your prompt:
"add retry logic to the stripe webhook handler"

Evidence:
Claude asked: "Should I use exponential backoff? Do you want idempotency checks?"

Missing decision:
What's the retry policy?

Quick rule:
When asking for retries, specify: backoff strategy + max attempts + failure handling

───────────────────────────────────────────────────────────────
```

---

## What Makes It Different

- **Dangerous verbs focus** — Not generic "vague prompt" detection. Targets specific high-risk patterns.
- **Tier A evidence only** — Only surfaces gaps where Claude actually asked or assumed something.
- **Precision over recall** — Max 1 gap per run. No overwhelming reports.
- **Compounds into CLAUDE.md** — Patterns become permanent rules (on 3rd occurrence).

---

## Evidence Tiers

| Tier | What Counts | Shown? |
|------|-------------|--------|
| A | Claude asked clarifying Q, stated assumption, user corrected | Yes |
| B | Surrounding context | On-demand |
| C | Keywords only | Never |

**No evidence = no gap shown.**

---

## Repository Structure

```
apprenticemode/
├── .claude-plugin/
│   └── marketplace.json
├── plugins/
│   └── claudebrief/
│       ├── .claude-plugin/
│       │   └── plugin.json
│       ├── commands/
│       │   └── reflect.md
│       └── CLAUDE.md
└── README.md
```

---

## License

MIT
