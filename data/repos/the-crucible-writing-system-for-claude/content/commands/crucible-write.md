---
allowed-tools: Read, Write, Edit, Glob, Grep, Bash, Task
argument-hint: [chapter number] | "continue"
description: Draft novel prose scene-by-scene from Crucible outlines. Includes automatic bi-chapter reviews every 2 chapters.
---

# /crucible-write

Write prose drafts scene-by-scene from your Crucible outlines.

## Execution Instructions

**IMPORTANT:** When this command is invoked, you MUST:

1. **Invoke the `crucible-suite:crucible-writer` skill** using the Skill tool
2. The skill will guide you through the complete writing workflow
3. Verify outline is complete before writing
4. Follow the skill's scene-by-scene protocol
5. **Track chapter completion** and trigger bi-chapter reviews per the skill instructions
6. When triggering reviews, use the Task tool with the 5 review agent subagent_types

## Usage

- `/crucible-suite:crucible-write 1` - Start writing chapter 1
- `/crucible-suite:crucible-write 5` - Jump to chapter 5
- `/crucible-suite:crucible-write continue` - Resume from last position

## What This Does

1. Activates the crucible-writer skill (you must follow its instructions)
2. Loads your chapter outline and style profile
3. Writes scene-by-scene, matching your voice
4. Verifies each scene against the outline
5. Updates the story bible after each chapter
6. Triggers bi-chapter reviews every 2 chapters

## Bi-Chapter Reviews

Every 2 chapters, the system automatically triggers a comprehensive review using 5 specialized agents:

| Agent | Focus |
|-------|-------|
| voice-checker | Style/voice consistency |
| continuity-checker | Plot/character continuity |
| outline-checker | Outline fidelity |
| timeline-checker | Chronological consistency |
| prose-checker | Craft-level feedback |

Reviews help catch issues before they compound.

## Prerequisites

Requires:
1. Completed chapter outlines from `/crucible-suite:crucible-outline`
2. Style sample (2,000+ words of your writing) OR style preferences
3. Crucible Summary Card (for quick reference)

## Word Count Targets

Default targets per scene type:
- Action scenes: 1,500-2,500 words
- Dialogue scenes: 1,000-2,000 words
- Reflection scenes: 800-1,500 words
- Transitional scenes: 500-1,000 words

## Anti-Hallucination

The writer follows strict rules:
- Never invents plot points not in the outline
- Never creates characters not in the Constellation Bible
- Flags any [INVENTED] minor details for your review
- Asks you if anything seems missing rather than guessing

## What Happens Next

After writing completes, you can:
- `/crucible-suite:crucible-edit` - Begin revision and editing
- `/crucible-suite:crucible-review` - Trigger a manual review
- Continue writing more chapters
