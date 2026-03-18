# Vibe Writing Guide

A complete tutorial on the wiki-style fiction writing methodology.

## Philosophy

Traditional fiction writing often follows a linear process: outline, draft, revise. Vibe Writing takes a different approach inspired by wiki systems and world-building games.

**Core Insight**: A story is a web of interconnected elements. Characters know other characters. Locations have histories. Events cause other events. By making these connections explicit and navigable, you create a richer, more consistent narrative.

**The "Vibe"**: Instead of forcing every detail upfront, you capture the essence—the vibe—of story elements as they emerge. Red links serve as creative promises: "something goes here." This reduces pressure and encourages organic development.

## The Method

### Step 1: Start with Structure

Before any prose, establish your foundation:

1. **OUTLINE.md**: High-level story structure (acts, major beats)
2. **TIMELINE.md**: Chronological order of events
3. **TODO.md**: Tracking file for red links and tasks

This outline-first approach ensures you know where the story is going before building the wiki.

### Step 2: Create Core Entries

Identify the essential elements:

- **Protagonist**: Who drives the story?
- **Antagonist/Conflict**: What opposes them?
- **Key Locations**: Where does the action happen?
- **Central Concept**: What's the story really about?

Create wiki entries for each. Don't worry about completeness—capture what you know.

### Step 3: Write with Links

As you write entries or chapters, use `[[Double Brackets]]` for:

- Character names: `[[Detective Chen]]`
- Locations: `[[Gray Manor]]`
- Concepts: `[[The Gray Inheritance]]`
- Events: `[[The Night of the Murder]]`

If the linked file doesn't exist, that's a **red link**. Add it to `TODO.md`.

### Step 4: Fill Red Links

Red links are opportunities, not obligations. Fill them when:

- You need that information for current writing
- Inspiration strikes
- The story naturally develops that direction

Each filled red link often creates new red links—this is healthy growth.

### Step 5: Iterate

Writing is rewriting. As the story develops:

- Update entries with new information
- Add cross-references you missed
- Resolve contradictions
- Deepen shallow entries

## Working with AI

### Setting Up

The `CLAUDE.md` file provides context for AI assistants. It contains:

- Project overview and genre
- Linking and formatting rules
- Document templates
- Creative guidelines
- Current project status

When starting a session, ensure the AI has read `CLAUDE.md`.

### Effective Prompts

**World-Building Expansion**:
```
Read [[Character Name]] and expand their background section.
Consider their relationship with [[Other Character]] and
their connection to [[Location]].
```

**Chapter Writing**:
```
Write Chapter 3 following the outline. The POV character is
[[Detective Chen]]. The scene takes place at [[Gray Manor]].
Reference the established character voices.
```

**Consistency Check**:
```
Review all character files and check for contradictions in
relationships. Cross-reference with TIMELINE.md.
```

**Brainstorming**:
```
I need a reason for [[Marcus Wells]] to be at [[Gray Manor]]
the night of the murder that isn't the obvious one. Give me
three possibilities with different implications.
```

**Deep Discussion**:
```
Let's discuss [[Elena Ross]]'s motivation. Is her loyalty to
Victoria believable? What would make her character more complex?
```

### AI Best Practices

1. **Provide Context**: Reference specific files rather than expecting memory
2. **Be Specific**: "Expand the personality section" beats "make it better"
3. **Iterate**: Build on AI suggestions rather than accepting first drafts
4. **Verify**: Check AI additions against established facts
5. **Challenge**: Ask "why" and "what if" to deepen ideas

## File Organization

### Naming Convention

- **PascalCase** for all content files: `DetectiveChen.md`
- **Chapters** use numbered prefix: `01.TheCall.md`
- **System files** use UPPERCASE: `README.md`, `CLAUDE.md`

### Link Syntax

Links use the display name, which maps to PascalCase filename:

| Link Syntax | File |
|-------------|------|
| `[[Detective Chen]]` | `DetectiveChen.md` |
| `[[Gray Manor]]` | `GrayManor.md` |
| `[[The Gray Inheritance]]` | `TheGrayInheritance.md` |

### YAML Frontmatter

Every content file starts with:

```yaml
---
type: character | location | concept | chapter
tags: [relevant, tags]
---
```

Types help with organization. Tags enable filtering and grouping.

## Workflow Example

### Day 1: Foundation
1. Write `OUTLINE.md` with three-act structure
2. Create protagonist entry with basic info
3. Create primary location entry
4. Note red links in `TODO.md`

### Day 2: Expansion
1. Fill the most critical red link
2. Write first chapter draft
3. Add new red links discovered while writing
4. Update `TIMELINE.md` with chapter events

### Day 3: Deepening
1. Use AI to expand a character's background
2. Cross-reference for consistency
3. Add sensory details to location
4. Connect new details to existing entries

### Ongoing
- Pick a red link that interests you
- Write when inspired, organize when not
- Let the wiki grow with the story
- Revisit and refine as you learn more about your world

## Tips and Tricks

### Dealing with Red Links

- **Too many**: Focus on plot-critical ones first
- **Stuck on one**: Skip it, write around it, return later
- **Proliferating**: That's usually good—rich world

### Maintaining Consistency

- Re-read related entries before writing new content
- Keep `TIMELINE.md` updated
- Use AI for consistency checks
- Trust your wiki more than your memory

### Avoiding Over-Planning

- You don't need every entry complete before writing
- Minimum viable entries are fine
- Let prose writing reveal what needs expansion
- Red links are permission to figure it out later

### When to Break Rules

This method serves the writing, not vice versa. If something isn't working:

- Reorganize files as needed
- Add new entry types
- Modify templates
- Skip formalities when inspired

## Common Patterns

### Character Relationship Web
Create entries for all major characters first, even if sparse. Add relationship sections that link to each other. This reveals gaps and connections.

### Location as Character
Treat significant locations like characters with history, personality (atmosphere), and relationships (who inhabits/visits them).

### Concept as Anchor
Abstract ideas (themes, mysteries, rules) deserve their own entries. They help maintain focus and consistency.

### Chapter as Proof
Each chapter should demonstrate that the wiki entries work together. If you can't write the chapter, something is missing or wrong in the wiki.

## Getting Started

1. Fork or copy this demo repository
2. Read through all files to understand the format
3. Replace the demo content with your story
4. Start with `OUTLINE.md` for your plot
5. Create entries for your core elements
6. Write your first chapter
7. Let the red links guide your expansion

Happy writing!
