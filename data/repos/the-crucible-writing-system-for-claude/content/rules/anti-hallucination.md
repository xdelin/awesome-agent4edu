---
paths:
  # Primary Crucible draft locations (where prose is written)
  - draft/**/*.md
  - draft/chapters/**/*.md
  - manuscript/**/*.md
  # Chapter files by naming pattern
  - "**/chapter-*.md"
  - "**/chapter_*.md"
  - "**/ch-*.md"
  - "**/ch_*.md"
  # Scene files
  - "**/scene-*.md"
  - "**/scene_*.md"
---

# Anti-Hallucination Rules

## Core Principle

**Never invent what should exist in the source documents.**

The Crucible system maintains planning documents, outlines, and story bibles specifically to prevent hallucination. These documents are the single source of truth.

## Verification Requirements

### Before Writing Any Scene

1. **Load the scene outline** - Know what must happen
2. **Check character details** - Verify against Constellation Bible
3. **Check location details** - Verify against World Forge
4. **Check timeline** - Verify against chapter sequence
5. **Note any gaps** - Flag missing information

### During Writing

1. **Use only established names** - No inventing character names
2. **Use only established places** - No inventing location names
3. **Follow magic rules** - No new abilities without establishment
4. **Track what you write** - Note new details for story bible

### After Writing

1. **Verify against outline** - All required elements present?
2. **Check for inventions** - Did you add unplanned details?
3. **Update story bible** - Record any new minor details
4. **Flag uncertainties** - Mark [INVENTED] or [VERIFY] tags

## The Three Laws of Anti-Hallucination

### Law 1: If It's Not in the Outline, Don't Write It

The outline defines what happens. If a scene isn't outlined:
- ASK the author before writing it
- Don't assume it should exist
- Don't improvise plot points

### Law 2: If You're Unsure, Ask

When uncertain about any detail:
- Character appearance? Check story bible. Still unsure? ASK.
- Location description? Check world forge. Still unsure? ASK.
- Plot point? Check outline. Still unsure? ASK.

**Asking is always better than guessing.**

### Law 3: If You Must Invent, Flag It

Minor details sometimes need invention (a servant's name, a meal's description). When this happens:

```
The servant—[INVENTED: name "Mira"]—brought the tea.
```

Or in a note:
```
[INVENTED DETAILS THIS SCENE:
- Servant name: Mira
- Tea type: jasmine
- Room has blue curtains]
```

Author can then approve, reject, or modify.

## What Can Be Invented (Minor)

- Incidental character names (servants, shopkeepers, crowds)
- Specific food/drink descriptions
- Minor environmental details (weather, time of day when not specified)
- Transitional actions (walking, sitting, etc.)

## What Cannot Be Invented (Major)

- Named character traits or history
- Plot events or revelations
- Magic system capabilities
- Relationship states or changes
- World-building rules
- Character motivations or secrets

## Verification Checkpoints

### Per Scene
- [ ] Scene outline loaded
- [ ] Required elements listed
- [ ] Character details verified
- [ ] Location details verified
- [ ] No major inventions made
- [ ] Minor inventions flagged

### Per Chapter
- [ ] All scenes match outline
- [ ] Character states consistent
- [ ] Timeline logical
- [ ] Story bible updated
- [ ] Inventions reviewed with author

### Per Session
- [ ] Working from correct outline version
- [ ] Story bible is current
- [ ] Previous chapter summary loaded
- [ ] No assumptions made without verification

## When in Doubt

The golden rule: **ASK THE AUTHOR**

Better to pause and verify than to write something that contradicts established canon. The author knows their story better than any AI can infer.
