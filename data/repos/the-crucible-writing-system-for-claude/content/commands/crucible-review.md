---
allowed-tools: Read, Glob, Grep, Bash, Task
argument-hint: [chapter range] | "latest"
description: Manually trigger a bi-chapter review using 5 specialized agents. Use between automatic reviews or for specific chapters.
---

# /crucible-review

Trigger a manual review of your chapters using specialized review agents.

## Usage

- `/crucible-suite:crucible-review latest` - Review most recent 2 chapters
- `/crucible-suite:crucible-review 5-6` - Review chapters 5 and 6
- `/crucible-suite:crucible-review 1-4` - Review chapters 1 through 4

## Execution Instructions

**IMPORTANT:** When this command is invoked, you MUST:

0. **ALWAYS use the AskUserQuestion tool** for presenting options to the user (NOT plain text A/B/C options)

1. **Detect project root** - Find the Crucible project directory:
   - Use Glob to search for `.crucible/` directory: `Glob(pattern="**/.crucible")`
   - The project root is the PARENT of the `.crucible/` directory found
   - Store this as `PROJECT_ROOT` (absolute path) for all subsequent steps
   - If no `.crucible/` found, inform user this must be run from within a Crucible project

2. **Determine chapter range** from the argument:
   - `latest` → Read `{PROJECT_ROOT}/.crucible/state/draft-state.json` to find the last 2 completed chapters
   - `X-Y` → Use chapters X through Y
   - Single number → Review that chapter only

3. **Resolve all file paths** - Build absolute paths for agents:
   - Chapter files: `{PROJECT_ROOT}/draft/chapters/ch{NN}.md`
   - Story bible: `{PROJECT_ROOT}/story-bible.json`
   - Style profile: `{PROJECT_ROOT}/style-profile.json`
   - Master outline: `{PROJECT_ROOT}/outline/master-outline.md`
   - Scene breakdown: `{PROJECT_ROOT}/outline/scene-breakdown.md`
   - World forge: `{PROJECT_ROOT}/planning/world-forge.md`

4. **Launch ALL 5 review agents in parallel** using the Task tool with ABSOLUTE PATHS:

```
Task(subagent_type="crucible-suite:voice-checker", prompt="Review chapters X-Y for voice consistency.

PROJECT ROOT: {PROJECT_ROOT}

CHAPTERS TO REVIEW:
- {PROJECT_ROOT}/draft/chapters/ch05.md
- {PROJECT_ROOT}/draft/chapters/ch06.md

REFERENCE FILES:
- Style profile: {PROJECT_ROOT}/style-profile.json

Analyze for voice consistency against the style profile. Return a structured report following your output format.")

Task(subagent_type="crucible-suite:continuity-checker", prompt="Review chapters X-Y for continuity errors.

PROJECT ROOT: {PROJECT_ROOT}

CHAPTERS TO REVIEW:
- {PROJECT_ROOT}/draft/chapters/ch05.md
- {PROJECT_ROOT}/draft/chapters/ch06.md

REFERENCE FILES:
- Story bible: {PROJECT_ROOT}/story-bible.json
- Planning docs: {PROJECT_ROOT}/planning/

Check for continuity errors against the story bible. Return a structured report following your output format.")

Task(subagent_type="crucible-suite:outline-checker", prompt="Review chapters X-Y for outline adherence.

PROJECT ROOT: {PROJECT_ROOT}

CHAPTERS TO REVIEW:
- {PROJECT_ROOT}/draft/chapters/ch05.md
- {PROJECT_ROOT}/draft/chapters/ch06.md

REFERENCE FILES:
- Master outline: {PROJECT_ROOT}/outline/master-outline.md
- Scene breakdown: {PROJECT_ROOT}/outline/scene-breakdown.md
- Forge points: {PROJECT_ROOT}/planning/forge-points/

Verify adherence to the outline. Return a structured report following your output format.")

Task(subagent_type="crucible-suite:timeline-checker", prompt="Review chapters X-Y for timeline consistency.

PROJECT ROOT: {PROJECT_ROOT}

CHAPTERS TO REVIEW:
- {PROJECT_ROOT}/draft/chapters/ch05.md
- {PROJECT_ROOT}/draft/chapters/ch06.md

REFERENCE FILES:
- Story bible (timeline data): {PROJECT_ROOT}/story-bible.json
- World forge (travel times): {PROJECT_ROOT}/planning/world-forge.md

Analyze chronological consistency and temporal logic. Return a structured report following your output format.")

Task(subagent_type="crucible-suite:prose-checker", prompt="Review chapters X-Y for prose quality.

PROJECT ROOT: {PROJECT_ROOT}

CHAPTERS TO REVIEW:
- {PROJECT_ROOT}/draft/chapters/ch05.md
- {PROJECT_ROOT}/draft/chapters/ch06.md

REFERENCE FILES:
- Style profile: {PROJECT_ROOT}/style-profile.json

Analyze prose craft including pacing, show vs tell, and dialogue. Return a structured report following your output format.")
```

**NOTE:** Replace `{PROJECT_ROOT}` with the actual absolute path discovered in step 1, and adjust chapter numbers/filenames based on step 2.

5. **Wait for all agents to complete** and collect their reports

6. **Compile a consolidated report** showing:
   - All critical issues (from any agent)
   - All warnings (grouped by type)
   - All suggestions
   - Overall scores from each agent

7. **Update the draft state** to record the review using the update script:
```bash
python "${CLAUDE_PLUGIN_ROOT}/skills/crucible-writer/scripts/update_story_bible.py" "{PROJECT_ROOT}" --mark-review-complete <last_chapter_reviewed>
```
This sets `review_pending = false` and `last_review_at_chapter` to the last chapter number reviewed.

## What Each Agent Does

### 1. Voice Checker (`voice-checker`)
- Analyzes voice consistency against your style sample
- Identifies tone shifts
- Checks POV adherence
- Evaluates character voice distinction

### 2. Continuity Checker (`continuity-checker`)
- Verifies character details match story bible
- Checks physical descriptions
- Validates relationship states
- Flags potential plot holes

### 3. Outline Checker (`outline-checker`)
- Verifies all planned scenes are present
- Checks beat coverage
- Identifies deviations from outline
- Notes acceptable creative changes

### 4. Timeline Checker (`timeline-checker`)
- Verifies chronological consistency
- Checks time references
- Validates travel times
- Tracks character age/state progression

### 5. Prose Checker (`prose-checker`)
- Evaluates show vs tell balance
- Identifies pacing issues
- Assesses dialogue effectiveness
- Highlights strongest passages

## Review Output

Present a consolidated report:

```
═══════════════════════════════════════
BI-CHAPTER REVIEW COMPLETE
Chapters: [X-Y]
═══════════════════════════════════════

CRITICAL ISSUES (Must Fix)
─────────────────────────────────────
[Collected from all agents]

WARNINGS (Should Fix)
─────────────────────────────────────
Voice: [issues]
Continuity: [issues]
Outline: [issues]
Timeline: [issues]
Prose: [issues]

SUGGESTIONS (Consider)
─────────────────────────────────────
[Collected from all agents]

OVERALL SCORES
─────────────────────────────────────
Voice Consistency: X/10
Continuity: X/10
Outline Fidelity: X/10
Timeline: X/10
Prose Quality: X/10
AVERAGE: X/10

NEXT STEPS
─────────────────────────────────────
```

Then use AskUserQuestion:
```json
{
  "questions": [
    {
      "header": "Next Steps",
      "question": "How would you like to proceed after this review?",
      "options": [
        {"label": "Address issues", "description": "Fix critical issues before continuing"},
        {"label": "Continue writing", "description": "Note issues and continue (return later)"},
        {"label": "Show details", "description": "View detailed report from specific agent"}
      ],
      "multiSelect": false
    }
  ]
}
```

## When to Use

- When you want feedback before the automatic bi-chapter trigger
- After making significant revisions
- When something feels "off" but you can't identify it
- Before sending chapters to beta readers
- To check specific problem areas

## Prerequisites

Requires written chapters to review. Works best with:
- Style profile on file (`style-profile.json` at project root)
- Story bible populated (`story-bible.json` at project root)
- Chapter outline available (`outline/master-outline.md` at project root)
