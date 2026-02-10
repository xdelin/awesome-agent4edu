# claudebrief

Claude Code plugin that surfaces knowledge gaps from your sessions. Detects "dangerous verbs" that hide required decisions, tracks patterns, and compounds learnings into CLAUDE.md.

## How It Works

1. Run `/claudebrief:reflect` after a session (or `current` for this session)
2. Plugin scans for dangerous verbs + Tier A evidence
3. Shows max 1 gap with evidence and quick rule
4. User confirms → gap saved
5. On 3rd occurrence → offers CLAUDE.md patch

## Dangerous Verbs

| Verb | Missing Decision |
|------|------------------|
| retry, backoff | What's the retry policy? |
| cache, memoize | What's the invalidation policy? |
| optimize, faster | What's the performance target? |
| prod-ready, ship | What's the rollback plan? |
| refactor, clean up | What must stay the same? |
| async, parallel | What's the ordering guarantee? |

## Evidence Tiers

| Tier | What Counts | Action |
|------|-------------|--------|
| A | Claude asked clarifying Q, stated assumption, user corrected | Show |
| B | Surrounding context | On-demand only |
| C | Keywords only | Never show |

**Rule:** Only surface gaps with Tier A evidence.

## Commands

- `/claudebrief:reflect` — Analyze recent sessions
- `/claudebrief:reflect current` — Analyze this session

## Data

- `~/.claude/claudebrief/confirmed.json` — Saved gaps with recurrence tracking

## Key Constraints

- **Tier A evidence required** — No evidence = no surface
- **Max 1 gap per run** — Precision over recall
- **User confirms** — Only save if user says yes
- **3+ to suggest CLAUDE.md** — Patterns, not one-offs
