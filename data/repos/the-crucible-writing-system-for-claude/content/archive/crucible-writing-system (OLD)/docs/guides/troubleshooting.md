# Troubleshooting Guide

Solutions to common issues when using the Crucible Writing System.

## Planning Issues

### Problem: Claude doesn't recognize I want to use the Crucible System

**Symptoms:**
- Claude offers generic writing help
- No multi-choice questions appear
- Planning documents aren't generated

**Solutions:**
1. Be explicit: "Use the Crucible Planner skill to plan my novel"
2. Mention the Crucible Structure directly
3. Check that the `crucible-planner` skill is enabled

---

### Problem: The multi-choice options don't fit my story

**Symptoms:**
- None of the A-E options work
- Your vision differs from suggestions

**Solution:**
Always use option "Other (describe)" and explain what you want. The options are starting points, not limitations.

---

### Problem: I want to change an earlier answer

**Symptoms:**
- Realized a better choice
- Earlier answer conflicts with later decisions

**Solution:**
Tell Claude: "I want to change my answer for [specific question]. Instead of [old answer], I want [new answer]."

Claude will update and may ask follow-up questions to maintain coherence.

---

### Problem: Planning session broke mid-process

**Symptoms:**
- Connection lost
- Had to close the session
- Returning after a break

**Solution:**
1. Start a new session
2. Say: "I want to continue my Crucible planning"
3. Upload your project state if you have it
4. Claude will resume from the last saved point

If you don't have saved state:
- Describe where you were
- Upload any documents already generated
- Claude can reconstruct and continue

---

### Problem: The generated documents have gaps

**Symptoms:**
- Missing information in planning docs
- Sections that say "TBD" or are incomplete

**Solution:**
Tell Claude: "The [document name] is missing [specific information]. Can we fill that in?"

Claude will ask the necessary questions to complete it.

---

## Outlining Issues

### Problem: The beat-to-chapter mapping doesn't work for my story

**Symptoms:**
- Standard 25-chapter structure feels wrong
- Certain beats need more/less space

**Solutions:**
1. Request a different mapping (18, 25, or 35 chapters)
2. Ask Claude to adjust specific sections
3. Say: "I need [Movement X] to be longer because [reason]"

The mappings are flexible within the overall proportions.

---

### Problem: Can't figure out how to handle a Forge Point chapter

**Symptoms:**
- Unsure how to make all three strands converge
- The sacrifice isn't clear
- Crises feel forced

**Solutions:**
1. Review your Forge Point Blueprint from planning
2. Ask Claude: "Help me design the convergence for Forge Point [X]"
3. Ensure each strand has genuine stakes before convergence

Common issue: Crises not truly simultaneous. Fix by identifying a single event that triggers all three.

---

### Problem: Foreshadowing tracker is overwhelming

**Symptoms:**
- Too many threads to track
- Unsure what needs planting
- Payoffs feel disconnected

**Solution:**
Focus on **required** threads first:
- Each mercy (4 threads)
- Dark Mirror reveal
- Fire mastery path
- Core bond test
- Sacrifice foreshadowing

Optional threads can be added but aren't essential to structure.

---

### Problem: Outline feels too detailed/not detailed enough

**Symptoms:**
- Outline reads like prose (too detailed)
- Scenes are vague (not detailed enough)

**Solutions:**

Too detailed:
- Outlines should describe WHAT happens, not HOW it's written
- Remove prose-level description
- Keep goals, conflicts, turns

Not detailed enough:
- Each scene needs: goal, conflict, turn
- Key moments should be listed
- Plants/payoffs should be identified

---

## Writing Issues

### Problem: Claude's prose doesn't match my voice

**Symptoms:**
- Writing feels generic
- Word choices are wrong
- Rhythm is off

**Solutions:**
1. Provide a longer style sample (3,000+ words)
2. Ask Claude to analyze your sample and show the profile
3. Point out specific mismatches: "My style uses shorter sentences" or "I never use the word 'orbs' for eyes"
4. Request a recalibration sample before continuing

---

### Problem: Claude invented something not in my outline

**Symptoms:**
- New character appeared
- Plot event wasn't planned
- World detail contradicts planning docs

**Solutions:**
1. All inventions should be flagged `[INVENTED]`
2. If you like it: approve it for the story bible
3. If you don't: request revision to match outline
4. Remind Claude: "The outline is law—don't invent plot elements"

---

### Problem: Session broke mid-chapter

**Symptoms:**
- Lost connection during writing
- Scene incomplete
- Unsure what was saved

**Solution:**
1. Resume session
2. Say: "Continue my writing session for [title]"
3. Claude will load the story bible and resume
4. If uncertain, ask: "Where did we leave off?"

**Prevention:**
- Save after every scene
- Keep the story bible updated
- Note exact stopping points

---

### Problem: Word counts are wrong

**Symptoms:**
- Chapter running too long
- Chapter too short
- Pacing feels off

**Solutions:**

**Too long:**
- Request tighter prose
- Ask Claude to trust the reader more
- Identify padding that can be cut
- Say: "This scene should be ~1,500 words, not 2,500"

**Too short:**
- Request expanded sensory details
- Add interiority at key moments
- Lengthen dialogue exchanges
- Say: "This chapter needs more breathing room"

---

### Problem: Continuity error detected

**Symptoms:**
- Character in wrong place
- Event contradicts earlier scene
- Information character shouldn't have

**Solutions:**
1. Check story bible for character states
2. Revise the problematic section
3. Update story bible after fix
4. If error is in story bible, correct it first

---

### Problem: Claude keeps asking me to verify instead of writing

**Symptoms:**
- Too many clarification requests
- Progress feels slow
- Questions about obvious things

**Solutions:**
1. Provide more detailed outlines
2. Say: "For minor details, make reasonable choices and flag with [INVENTED]"
3. Front-load answers in the scene brief
4. Confirm: "Proceed without asking unless there's a genuine gap"

---

## Technical Issues

### Problem: Can't find my saved documents

**Symptoms:**
- Documents not appearing
- Lost previous work
- File structure unclear

**Solution:**
Ask Claude: "Show me the file structure for my project"

Standard locations:
```
planning/           <- Planning documents
outline/            <- Outline documents  
manuscript/         <- Written prose
story-bible.json    <- Continuity tracker
```

---

### Problem: Documents won't upload

**Symptoms:**
- Upload fails
- Claude doesn't see attached files
- Format issues

**Solutions:**
1. Use standard formats (.md, .txt, .docx)
2. Try pasting content directly if upload fails
3. Break large documents into smaller parts
4. Ensure file isn't corrupted

---

### Problem: Claude seems confused about my project

**Symptoms:**
- Mixed up details from different stories
- Forgot what we discussed
- Asking questions already answered

**Solutions:**
1. Start with clear context: "We're working on [title], a [description]"
2. Upload the Crucible Summary Card at session start
3. Reference specific documents when relevant
4. Keep sessions focused on one project

---

## General Tips

### If Something Isn't Working

1. **Be explicit** — Tell Claude exactly what's wrong
2. **Reference documents** — "According to my Crucible Thesis..."
3. **Request specific changes** — Not "fix this" but "change X to Y"
4. **Save frequently** — Prevents major loss
5. **Take breaks** — Fresh eyes help

### Prevention

- Complete each phase before moving on
- Review documents before proceeding
- Keep backups of important documents
- Note where you stopped in each session
- Use clear, consistent naming

### When to Start Over

Consider restarting a phase if:
- Fundamental premise has changed
- Multiple documents have major errors
- Coherence is lost
- You'd spend more time fixing than redoing

Usually better to adjust than restart, but sometimes a clean slate helps.

---

## Getting Help

If you encounter issues not covered here:

1. **Describe the problem clearly** — What happened? What did you expect?
2. **Share relevant documents** — Context helps
3. **Specify the phase** — Planning, outlining, or writing?
4. **Note what you've tried** — Prevents repeated suggestions

Claude can usually diagnose and fix issues with enough information.

---

## Related Documentation

- [Getting Started](getting-started.md)
- [Crucible Planner](../skills/crucible-planner.md)
- [Crucible Outliner](../skills/crucible-outliner.md)
- [Crucible Writer](../skills/crucible-writer.md)
