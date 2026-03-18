# Claude Code Instructions

This file provides context and guidelines for AI collaboration on this fiction project.

## Project Overview

**Project**: The Gray Inheritance
**Genre**: Modern Mystery / Detective Fiction
**Method**: Wiki-style world-building with interconnected entries

This repository uses the "Vibe Writing" methodology where every story element (character, location, concept) is a wiki entry linked via `[[Name]]` syntax.

## Global Rules

### 1. Double-Bracket Linking
- Use `[[Name]]` to reference other entries
- Link names must match file names in PascalCase
- Example: `[[Detective Chen]]` → `DetectiveChen.md`
- Create red links freely; track them in `TODO.md`

### 2. YAML Frontmatter
Every content file must begin with:
```yaml
---
type: character | location | concept | chapter
tags: [relevant, tags, here]
---
```

### 3. File Naming
- Use PascalCase: `DetectiveChen.md`, `GrayManor.md`
- Chapters use numbered prefix: `01.TheCall.md`, `02.TheScene.md`
- No spaces or special characters in filenames

### 4. Consistency
- Check existing entries before adding new details
- Respect established timeline in `TIMELINE.md`
- Follow character traits and relationships defined in their files
- Flag contradictions rather than silently overwriting

## Document Templates

### Character Template
```markdown
---
type: character
tags: [role-tags]
---

# Character Name

## Overview
Brief introduction and role in the story.

## Background
History, profession, key life events.

## Personality
Traits, motivations, flaws.

## Relationships
- [[Other Character]]: Nature of relationship

## Appearance
Physical description.

## Notes
Writer notes, development ideas.
```

### Location Template
```markdown
---
type: location
tags: [location-type]
---

# Location Name

## Overview
What this place is and its significance.

## Description
Physical details, atmosphere, sensory elements.

## History
Relevant background.

## Key Areas
Notable sub-locations or rooms.

## Associated Characters
- [[Character]]: Their connection to this place

## Notes
Writer notes.
```

### Concept Template
```markdown
---
type: concept
tags: [concept-type]
---

# Concept Name

## Overview
What this concept represents in the story.

## Details
Explanation and mechanics.

## Relevance
How it connects to the plot and characters.

## Notes
Writer notes.
```

### Chapter Template
```markdown
---
type: chapter
tags: [act-number, themes]
---

# Chapter Title

## Summary
Brief chapter overview (spoilers OK here).

## Scene Breakdown
- Scene 1: Description
- Scene 2: Description

## Characters Present
- [[Character 1]]
- [[Character 2]]

## Locations
- [[Location]]

---

## Chapter Content

[The actual prose goes here]
```

## Creative Guidelines

### Voice and Tone
- Third person limited, following Detective Chen
- Atmospheric, noir-influenced but modern
- Dialogue should reveal character
- Show don't tell for emotions

### Mystery Writing
- Plant clues fairly; readers should be able to solve it
- Red herrings must have alternate explanations
- Character motivations must be believable
- The solution should be surprising but inevitable in hindsight

### World Consistency
- Modern day, unnamed American city
- Realistic procedural elements (police work, forensics)
- The Gray family represents old money and secrets
- Weather and atmosphere should reflect mood

## AI Collaboration Modes

### World-Building Expansion
When asked to expand a wiki entry:
- Read related entries first for consistency
- Add details that create hooks for other entries
- Suggest new red links for TODO.md
- Maintain established facts

### Chapter Writing
When asked to write chapter content:
- Follow the outline in OUTLINE.md
- Reference character files for voice and personality
- Include sensory details from location files
- Weave in clues and foreshadowing naturally

### Consistency Checking
When asked to verify consistency:
- Cross-reference timeline
- Check character relationship symmetry
- Verify location details match across files
- Flag any contradictions found

### Creative Brainstorming
When asked for ideas:
- Propose multiple options with trade-offs
- Consider impact on existing story elements
- Think about thematic resonance
- Suggest concrete next steps

### Deep Discussion
When discussing characters or themes:
- Explore psychological depth
- Consider real-world parallels
- Discuss narrative purpose
- Challenge assumptions constructively

## Current Status

Check `TODO.md` for:
- Current writing stage
- Red links needing content
- Plot tasks in progress

Check `OUTLINE.md` for:
- Overall story structure
- What has been written
- What comes next
