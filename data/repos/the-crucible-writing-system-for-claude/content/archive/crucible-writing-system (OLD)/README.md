# The Crucible Writing System

A comprehensive AI-assisted novel planning, outlining, and writing system for epic fantasy fiction.

## What is the Crucible Writing System?

The Crucible Writing System is an integrated suite of three Claude skills designed to guide writers from initial story concept to completed first draft. Built around the **Crucible Structure**—a 36-beat narrative framework with three interwoven story strands—the system ensures thematic coherence, structural integrity, and narrative craft throughout the entire writing process.

## The Three Skills

| Skill | Purpose | When to Use |
|-------|---------|-------------|
| **[Crucible Planner](docs/skills/crucible-planner.md)** | Transform a premise into 7 comprehensive planning documents | Starting a new novel; "plan my fantasy book" |
| **[Crucible Outliner](docs/skills/crucible-outliner.md)** | Convert planning documents into chapter-by-chapter outlines | After planning is complete; "outline my book" |
| **[Crucible Writer](docs/skills/crucible-writer.md)** | Draft prose scene-by-scene from your outline | After outlining; "write my novel" |

## The Workflow

```
┌─────────────────────────────────────────────────────────────────────────┐
│                        CRUCIBLE WRITING SYSTEM                          │
├─────────────────────────────────────────────────────────────────────────┤
│                                                                         │
│   PHASE 1: PLANNING                    PHASE 2: OUTLINING              │
│   ┌─────────────────────┐              ┌─────────────────────┐         │
│   │  Crucible Planner   │              │  Crucible Outliner  │         │
│   │                     │              │                     │         │
│   │  Premise → 7 Docs   │    ────►     │  Docs → Chapters    │         │
│   │                     │              │                     │         │
│   │  • Crucible Thesis  │              │  • Beat mapping     │         │
│   │  • Strand Maps (3)  │              │  • Scene breakdown  │         │
│   │  • Forge Points     │              │  • Foreshadowing    │         │
│   │  • Dark Mirror      │              │  • Pacing           │         │
│   │  • Constellation    │              │                     │         │
│   │  • Mercy Ledger     │              │                     │         │
│   │  • World Forge      │              │                     │         │
│   └─────────────────────┘              └─────────────────────┘         │
│              │                                    │                     │
│              │                                    │                     │
│              └─────────────────┬──────────────────┘                     │
│                                │                                        │
│                                ▼                                        │
│                       PHASE 3: WRITING                                  │
│                       ┌─────────────────────┐                          │
│                       │   Crucible Writer   │                          │
│                       │                     │                          │
│                       │  Outline → Prose    │                          │
│                       │                     │                          │
│                       │  • Scene-by-scene   │                          │
│                       │  • Style matching   │                          │
│                       │  • Continuity       │                          │
│                       │  • Story bible      │                          │
│                       └─────────────────────┘                          │
│                                │                                        │
│                                ▼                                        │
│                      📚 COMPLETED MANUSCRIPT                            │
│                                                                         │
└─────────────────────────────────────────────────────────────────────────┘
```

## The Crucible Structure

At the heart of the system is the **Crucible Structure**, a 36-beat narrative architecture with:

### Three Story Strands

| Strand | Focus | What It Tracks |
|--------|-------|----------------|
| **QUEST** | External | The mission, burden, or objective |
| **FIRE** | Internal | The power, curse, or transformation |
| **CONSTELLATION** | Relational | The bonds that anchor or break |

### Five Movements + Coda

| Movement | Name | % of Book | Function |
|----------|------|-----------|----------|
| I | Ignition | 10% | Light the forge |
| II | First Tempering | 20% | Shape through failure |
| III | Scattering | 25% | Expand, harden, fragment |
| IV | Brightest Burning | 25% | Master, gather, choose |
| V | Final Forging | 15% | Converge, fail, transcend |
| Coda | Tempered Blade | 5% | Reveal what was made |

### Four Forge Points + Apex

Critical convergence moments where all three strands must be in simultaneous crisis. The protagonist cannot resolve all three—they must sacrifice one to save the others.

➡️ **[Read the complete Crucible Structure](docs/framework/crucible-structure.md)**

## Quick Start

### 1. Start Planning

Tell Claude:
> "I want to plan an epic fantasy novel using the Crucible Structure. My premise is: [your story idea]"

### 2. Answer Questions

Claude will guide you through multi-choice questions to build your planning documents. Expect ~60 questions across all documents.

### 3. Generate Outline

Once planning is complete:
> "Create a chapter outline from my Crucible planning documents"

### 4. Write Your Novel

With your outline ready:
> "Help me write my novel starting from Chapter 1"

## Documentation

### Core Framework
- **[The Crucible Structure](docs/framework/crucible-structure.md)** — The 36-beat narrative architecture
- **[Forge Points](docs/framework/forge-points.md)** — Strand convergence mechanics
- **[Dark Mirror](docs/framework/dark-mirror.md)** — Antagonist design principles
- **[Mercy Engine](docs/framework/mercy-engine.md)** — How mercy enables victory

### Skill Guides
- **[Crucible Planner Guide](docs/skills/crucible-planner.md)** — Planning skill documentation
- **[Crucible Outliner Guide](docs/skills/crucible-outliner.md)** — Outlining skill documentation
- **[Crucible Writer Guide](docs/skills/crucible-writer.md)** — Writing skill documentation

### User Guides
- **[Getting Started](docs/guides/getting-started.md)** — First-time user guide
- **[Multi-Book Series](docs/guides/series-planning.md)** — Planning trilogies and beyond
- **[Troubleshooting](docs/guides/troubleshooting.md)** — Common issues and solutions

### Templates
- **[Planning Document Templates](docs/templates/planning-templates.md)**
- **[Outline Templates](docs/templates/outline-templates.md)**
- **[Story Bible Templates](docs/templates/story-bible-templates.md)**

## Key Features

### ✅ Structural Integrity
Every beat, every chapter, every scene connects to the overall narrative architecture.

### ✅ Theme Integration
Your story's theme weaves through every element—never stated, always demonstrated.

### ✅ Foreshadowing Tracking
Plants and payoffs are systematically tracked to ensure satisfying setups and resolutions.

### ✅ Character Consistency
The story bible maintains character states, preventing continuity errors.

### ✅ Style Preservation
Your voice is captured and maintained throughout the entire manuscript.

### ✅ Anti-Hallucination
Rigorous verification protocols ensure Claude writes your story, not its own.

## Requirements

- Claude with the three Crucible skills enabled:
  - `crucible-planner`
  - `crucible-outliner`  
  - `crucible-writer`
- Your story premise or concept
- Time for the interactive planning process (~2-3 hours for complete planning)

## Trigger Phrases

Claude will activate the appropriate skill when you say things like:

**Planning:**
- "Plan my fantasy novel"
- "I have a story premise..."
- "Use the Crucible Structure to plan..."
- "Help me plan my book"

**Outlining:**
- "Outline my book"
- "Create chapter outlines"
- "Turn my plan into chapters"
- "Outline Book 1"

**Writing:**
- "Write my novel"
- "Draft chapter 1"
- "Start writing from my outline"
- "Help me write my book"

## Best Practices

1. **Complete each phase before moving on** — Don't skip planning or outlining
2. **Save frequently** — Sessions can break; saved work persists
3. **Provide style samples** — Better samples = better prose
4. **Answer "Other" thoughtfully** — The multi-choice questions are starting points
5. **Review generated documents** — Adjust before proceeding
6. **Trust the structure** — The 36 beats work; resist the urge to "improve" them mid-process

## Version

Crucible Writing System v1.0

---

*The Crucible Writing System was designed to help writers forge their stories through fire—to create narratives where transformation is earned, mercy matters, and wounds become gifts.*
