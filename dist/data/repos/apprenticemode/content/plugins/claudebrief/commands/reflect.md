---
description: Surface knowledge gaps from your Claude Code sessions
---

# claudebrief reflect

Analyzes your sessions to find knowledge gaps — specifically the "dangerous verbs" that hide required decisions.

---

## What It Detects

**Dangerous verbs** — words that sound complete but hide required constraints:

| Verb | What's Usually Missing |
|------|------------------------|
| retry, backoff | Retry policy (idempotency, max attempts, failure handling) |
| cache, memoize | Invalidation policy (TTL, consistency, eviction) |
| optimize, faster | Performance target (baseline, target metric, how to measure) |
| prod-ready, ship | Rollback plan (monitoring, alerting, SLOs) |
| refactor, clean up | Invariants (what must NOT change, test requirements) |
| async, parallel | Ordering guarantee (race conditions, cancellation, backpressure) |

---

## Mode Selection

Check `$ARGUMENTS`:

- `current` → Analyze THIS conversation (no file reads needed)
- Empty → Analyze recent session files

---

## Current Session Mode

If argument is `current`:

1. Look at the conversation history in this session
2. Scan for dangerous verbs in your messages
3. Check: Did Claude ask a clarifying question or state an assumption?
4. If yes → that's a gap (you didn't specify something Claude needed)

Present findings:

```
Quick take on this session...

{findings or "Looks clean - no gaps detected."}
```

Skip to "Present Gap" if gap found.

---

## Recent Sessions Mode

### Step 1: Find Sessions

```bash
find ~/.claude/projects -name "*.jsonl" -type f 2>/dev/null | grep -v "agent-" | xargs ls -t 2>/dev/null | head -10
```

### Step 2: Analyze (First 3 Sessions)

For each session, extract ordered messages:

```bash
cat "{path}" | jq -r 'select(.type == "user" or .type == "assistant") | "\(.type): \(.message.content)"' 2>/dev/null | head -150
```

Apply detection logic:
1. Find dangerous verbs in USER messages
2. Check if user specified the required decision
3. Check if Claude asked about it or made an assumption

**Stop at first HIGH confidence gap.**

---

## Detection Logic

### Dangerous Verb Patterns

| Category | Triggers | Required Decision | User Should Specify |
|----------|----------|-------------------|---------------------|
| RETRIES | retry, backoff, exponential | retry_policy | idempotency, dedup, max attempts, circuit breaker |
| CACHING | cache, memoize, redis | invalidation_policy | TTL, invalidation trigger, consistency |
| PERFORMANCE | faster, slow, optimize | performance_target | baseline metric, target metric, measurement |
| PRODUCTION | production, prod-ready, ship | rollback_plan | rollback, monitoring, alerting |
| REFACTORING | refactor, clean up, restructure | invariants | what stays same, test requirements |
| CONCURRENCY | async, parallel, concurrent | ordering_guarantee | ordering, error handling, cancellation |

### Evidence Tiers (Only Show Tier A)

| Tier | Evidence Type | Examples |
|------|---------------|----------|
| **A** | Claude asked clarifying question | "Should I use exponential backoff?", "Do you want TTL-based invalidation?" |
| **A** | Claude stated assumption | "I'll assume we want...", "By default I'll..." |
| **A** | User corrected Claude | "No, that's not...", "Wrong, it should be..." |
| **A** | User frustration | "ugh", "that's not what I meant", "again?" |
| B | Surrounding context | Other messages near the gap |
| C | Keyword heuristics | Pattern matching only |

**Rule:** Only surface gaps with Tier A evidence. No evidence = no gap shown.

### Confidence Levels

- **HIGH**: Claude explicitly asked about the missing decision OR Claude made an assumption
- **MEDIUM**: Dangerous verb present + user didn't specify + relevant context
- **LOW**: Pattern match only → don't surface

---

## Present Gap

When a gap is found:

```
───────────────────────────────────────────────────────────────
Gap: {CATEGORY} ({confidence} confidence)
───────────────────────────────────────────────────────────────

Your prompt:
"{user's original message}"

Evidence:
{Claude asked: "..." OR Claude assumed: "..." OR You corrected: "..."}

Missing decision:
{What's the retry policy? / What's the invalidation policy? / etc.}

Quick rule:
{When asking for X, specify: A, B, and C}

───────────────────────────────────────────────────────────────
```

**AskUserQuestion:**
- question: "Does this resonate?"
- header: "Gap"
- options:
  1. label: "Yes, save it" / description: "I should include this next time"
  2. label: "No, dismiss" / description: "Not relevant here"

---

## Handle Response

**If "Yes, save it":**

1. Ensure directory exists:
```bash
mkdir -p ~/.claude/claudebrief
```

2. Load existing gaps:
```bash
cat ~/.claude/claudebrief/confirmed.json 2>/dev/null || echo '{"gaps":[]}'
```

3. Add new gap to the array:
```json
{
  "id": "gap_{timestamp}",
  "category": "{RETRIES|CACHING|PERFORMANCE|PRODUCTION|REFACTORING|CONCURRENCY}",
  "evidence": "{exact quote}",
  "missing": "{the decision that was missing}",
  "one_liner": "{quick rule}",
  "timestamp": "{ISO timestamp}",
  "project": "{project name from path}"
}
```

4. Write updated gaps file

5. Count occurrences of this category. If >= 3:
   → Go to "CLAUDE.md Suggestion"

**If "No, dismiss":**
- Say "Got it, dismissed." and end.

---

## CLAUDE.md Suggestion

If this gap category has 3+ confirmed occurrences:

```
This is the {n}th time this has come up.

Suggested addition to your CLAUDE.md:

┌─────────────────────────────────────────────────────┐
│ ## {Category}                                       │
│ When asking for {trigger}, always specify:          │
│ - {first consideration}                             │
│ - {second consideration}                            │
│ - {third consideration}                             │
└─────────────────────────────────────────────────────┘
```

**Templates by category:**

**RETRIES:**
```markdown
## Retry Logic
When asking for retry logic, always specify:
- Idempotency requirements (can this operation be safely repeated?)
- Deduplication strategy (how to detect/prevent duplicates)
- Failure handling (what happens after max retries?)
```

**CACHING:**
```markdown
## Caching
When asking for caching, always specify:
- Invalidation strategy (when/how cache is cleared)
- TTL requirements (how long data stays fresh)
- Consistency requirements (eventual vs strong)
```

**PERFORMANCE:**
```markdown
## Performance
When asking for optimization, always provide:
- Current baseline metrics (what's the current latency/throughput?)
- Target metrics (what's the goal?)
- How to measure (what benchmark/profile to run)
```

**PRODUCTION:**
```markdown
## Production Readiness
When preparing for production, specify:
- Rollback strategy (how to undo if broken)
- Monitoring/alerting requirements (what to watch)
- SLO targets (what uptime/latency is required)
```

**REFACTORING:**
```markdown
## Refactoring
When asking for refactoring, specify:
- What behavior must NOT change (invariants)
- Test coverage expectations (what tests must pass)
- Breaking change tolerance (is breaking external API ok?)
```

**CONCURRENCY:**
```markdown
## Concurrency
When adding async/parallel code, specify:
- Ordering requirements (does sequence matter?)
- Error handling per task (fail-fast vs collect all?)
- Cancellation behavior (how to abort in-flight work)
```

**AskUserQuestion:**
- question: "Add this to your CLAUDE.md?"
- header: "CLAUDE.md"
- options:
  1. label: "Yes, add it" / description: "Append to ~/.claude/CLAUDE.md"
  2. label: "No, skip" / description: "Maybe later"

**If "Yes, add it":**
1. Read `~/.claude/CLAUDE.md` (create if doesn't exist)
2. Append the template
3. Confirm: "Added to ~/.claude/CLAUDE.md"

---

## Wrap Up

After processing:

```
Done.

{If gap saved: "Gap saved. Run /claudebrief:reflect again in a few sessions to track patterns."}
{If CLAUDE.md updated: "CLAUDE.md updated. Claude will see this rule in future sessions."}
{If no gaps: "Clean session - no obvious gaps found."}
```

---

## Key Constraints

- **Tier A evidence required** — No evidence = no surface
- **Max 1 gap per run** — Stop at first HIGH confidence gap
- **User confirms** — Only save if user says yes
- **3+ to suggest CLAUDE.md** — Don't suggest too early
- **No subagents** — All analysis in main context
- **Dangerous verbs focus** — Not generic "vague prompt" detection
