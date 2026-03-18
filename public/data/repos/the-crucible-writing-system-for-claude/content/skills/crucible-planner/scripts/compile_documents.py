#!/usr/bin/env python3
"""
Compile all planning documents from collected answers.

Usage:
    python compile_documents.py <project_path>

Generates all 9 planning documents plus a summary card.

Answer format support:
    Handles both old format (string values) and new format (dict with question/answer/description).
"""

import json
import os
import sys
from datetime import datetime


def load_state(project_path: str) -> dict:
    """Load state from project directory, checking new location first."""
    # New location: .crucible/state/planning-state.json
    new_path = os.path.join(project_path, ".crucible", "state", "planning-state.json")
    if os.path.exists(new_path):
        with open(new_path, 'r') as f:
            return json.load(f)

    # Legacy location: state.json at project root
    legacy_path = os.path.join(project_path, "state.json")
    if os.path.exists(legacy_path):
        with open(legacy_path, 'r') as f:
            return json.load(f)

    raise FileNotFoundError(
        f"No state file found. Checked:\n"
        f"  - {new_path}\n"
        f"  - {legacy_path}"
    )


def get_answer(answers: dict, key: str, default: str = '') -> str:
    """
    Extract answer value from answers dict.

    Handles both formats:
    - Old format: answers[key] is a string
    - New format: answers[key] is {"question": "...", "answer": "...", "description": "..."}

    Args:
        answers: The answers dict for a document
        key: The answer key to look up
        default: Default value if key not found

    Returns:
        The answer string value
    """
    value = answers.get(key)
    if value is None:
        return default
    if isinstance(value, dict):
        return value.get('answer', default)
    return value


def generate_crucible_thesis(state: dict, output_dir: str) -> str:
    """Generate the Crucible Thesis document."""
    answers = state["answers"]["doc1_crucible_thesis"]
    project = state["project"]

    content = f"""# Crucible Thesis: {project['title']}

## The Forging Question

> Can the protagonist become {get_answer(answers, 'forging_become', '[what they must become]')} without becoming {get_answer(answers, 'dark_mirror_represents', '[what the dark mirror represents]')}?

---

## The Three Strands

### The Burden (Quest Strand)
**Nature:** {get_answer(answers, 'burden_type', '[burden type]')}

The external weight that must be carried toward resolution.

### The Fire (Fire Strand)
**Nature:** {get_answer(answers, 'fire_type', '[fire type]')}

The internal power that threatens to consume.

### The Core Bond (Constellation Strand)
**Bond Type:** {get_answer(answers, 'core_bond_type', '[bond type]')}

The anchor that cannot break.

---

## The Dark Mirror

**Connection:** {get_answer(answers, 'dark_mirror_connection', '[connection]')}
**Their Truth:** "{get_answer(answers, 'antagonist_truth', '[antagonist philosophy]')}"

What the protagonist could become if the forging fails.

---

## The Theme

> {get_answer(answers, 'theme', '[theme statement]')}

---

## The Surrender

At the Apex, the protagonist must willingly surrender: **{get_answer(answers, 'surrender_type', '[surrender]')}**

---

## The Blade's Purpose

After the forging, the protagonist becomes: **{get_answer(answers, 'blade_purpose', '[purpose]')}**

---

## Original Premise

{project['premise']}
"""
    
    filepath = os.path.join(output_dir, "crucible-thesis.md")
    with open(filepath, 'w') as f:
        f.write(content)
    
    return filepath


def generate_quest_strand(state: dict, output_dir: str) -> str:
    """Generate the Quest Strand Map."""
    answers = state["answers"]["doc2_quest_strand"]
    project = state["project"]

    content = f"""# Quest Strand Map: {project['title']}

## The Burden

**Origin:** {get_answer(answers, 'burden_origin', '[origin]')}
**Why This Protagonist:** {get_answer(answers, 'why_protagonist', '[reason]')}

---

## Quest Arc by Movement

### Movement I — Discovery
The burden is discovered/received.

### Movement II — First Attempt
First significant obstacle; first failure.

### Movement III — Expansion
Quest scope expands: {get_answer(answers, 'quest_escalation', '[escalation]')}

### Movement IV — Convergence
All threads converging toward resolution.

### Movement V — Resolution
The burden ends: {get_answer(answers, 'resolution_method', '[resolution]')}

---

## The Antagonist's Stake

**Why They Oppose:** {get_answer(answers, 'antagonist_stake', '[stake]')}

---

## The Impossible Element

**What Makes It Hopeless:** {get_answer(answers, 'impossible_requirement', '[requirement]')}
"""
    
    filepath = os.path.join(output_dir, "strand-maps", "quest-strand.md")
    with open(filepath, 'w') as f:
        f.write(content)
    
    return filepath


def generate_fire_strand(state: dict, output_dir: str) -> str:
    """Generate the Fire Strand Map."""
    answers = state["answers"]["doc3_fire_strand"]
    project = state["project"]

    content = f"""# Fire Strand Map: {project['title']}

## The Fire

**Manifestation:** {get_answer(answers, 'fire_manifestation', '[manifestation]')}

---

## The Danger

**If Unmastered:** {get_answer(answers, 'fire_danger', '[danger]')}

---

## The Cost

**Each Use Costs:** {get_answer(answers, 'cost_of_use', '[cost]')}

---

## Fire Arc by Movement

### Movement I — Awakening
The Fire first manifests.

### Movement II — Instability
Control is partial/unreliable.

### Movement III — Hardening
**The Line Crossed:** {get_answer(answers, 'hardening_line', '[line]')}

### Movement IV — Mastery
**How Achieved:** {get_answer(answers, 'mastery_method', '[method]')}

### Movement V — Surrender/Transformation
The Fire is transformed or surrendered at the Apex.

---

## Mastery Requirements

**What Mastery Requires:** {get_answer(answers, 'mastery_path', '[requirements]')}
"""
    
    filepath = os.path.join(output_dir, "strand-maps", "fire-strand.md")
    with open(filepath, 'w') as f:
        f.write(content)
    
    return filepath


def generate_constellation_strand(state: dict, output_dir: str) -> str:
    """Generate the Constellation Strand Map."""
    answers = state["answers"]["doc4_constellation_strand"]
    project = state["project"]

    content = f"""# Constellation Strand Map: {project['title']}

## Core Constellation

### The Core Bond
From Crucible Thesis.

### The Faithful
**Who:** {get_answer(answers, 'faithful_companion', '[faithful]')}

### The Sacrifice
**Who Dies:** {get_answer(answers, 'sacrifice_character', '[sacrifice]')}

---

## Constellation Arc by Movement

### Movement I — Formation
Core companions commit.

### Movement II — Testing
First strain on relationships.

### Movement III — Expansion & Fracture
**New Allies:** {get_answer(answers, 'constellation_expansion', '[allies]')}
**The Strain:** {get_answer(answers, 'betrayal_source', '[strain]')}

### Movement IV — Sacrifice & Anchor
The sacrifice occurs; core bond proves unbreakable.

### Movement V — Resolution
Constellation holds through final trial.

### Coda — New Constellation
**After the Forging:** {get_answer(answers, 'constellation_fate', '[fate]')}

---

## The Bond That Saves

**Which Bond Anchors:** {get_answer(answers, 'bond_that_saves', '[bond]')}
"""
    
    filepath = os.path.join(output_dir, "strand-maps", "constellation-strand.md")
    with open(filepath, 'w') as f:
        f.write(content)
    
    return filepath


def generate_forge_points(state: dict, output_dir: str) -> list:
    """Generate all Forge Point documents."""
    fp_answers = state["answers"]["doc5_forge_points"]
    project = state["project"]
    filepaths = []

    forge_points = [
        ("fp0_ignition", "0", "Ignition", "I", "10%"),
        ("fp1_first", "1", "First Crucible", "II", "30%"),
        ("fp2_second", "2", "Second Crucible", "III", "55%"),
        ("fp3_third", "3", "Third Crucible", "IV", "75%"),
        ("apex", "Apex", "Willed Surrender", "V", "85%")
    ]

    for key, num, name, movement, percentage in forge_points:
        answers = fp_answers.get(key, {})

        content = f"""# Forge Point {num}: {name}

**Placement:** End of Movement {movement} (~{percentage})

---

## The Three Crises

### Quest Crisis
{get_answer(answers, 'quest_crisis', '[quest crisis]')}

### Fire Crisis
{get_answer(answers, 'fire_crisis', '[fire crisis]')}

### Constellation Crisis
{get_answer(answers, 'constellation_crisis', '[constellation crisis]')}

---

## The Impossible Choice

The protagonist cannot resolve all three.

---

## The Sacrifice

**What Is Sacrificed:** {get_answer(answers, 'sacrifice', '[sacrifice]')}

---

## The Aftermath

This shapes the protagonist and propels the story into the next movement.
"""
        
        filename = f"fp{num.lower().replace(' ', '-')}-{name.lower().replace(' ', '-')}.md"
        filepath = os.path.join(output_dir, "forge-points", filename)
        with open(filepath, 'w') as f:
            f.write(content)
        filepaths.append(filepath)
    
    return filepaths


def generate_dark_mirror(state: dict, output_dir: str) -> str:
    """Generate the Dark Mirror Profile."""
    answers = state["answers"]["doc6_dark_mirror"]
    project = state["project"]

    content = f"""# Dark Mirror Profile: {project['title']}

## The Parallel Path

**Origin Parallel:** {get_answer(answers, 'origin_parallel', '[parallel]')}

---

## The Divergence Point

**The Moment:** {get_answer(answers, 'divergence', '[divergence]')}

---

## The Antagonist's Truth

> "{get_answer(answers, 'antagonist_want', '[philosophy]')}"

**Why It's Compelling:** The truth within their worldview.

---

## The Dark Mirror's Offer (Beat 26)

**The Offer:** {get_answer(answers, 'compelling_offer', '[offer]')}
**Why It's Tempting:** {get_answer(answers, 'why_tempting', '[temptation]')}
**The Hidden Cost:** {get_answer(answers, 'hidden_cost', '[cost]')}

---

## The Defeat

**How They're Defeated:** {get_answer(answers, 'defeat_method', '[method]')}
**The Antagonist's End:** {get_answer(answers, 'antagonist_end', '[end]')}
"""
    
    filepath = os.path.join(output_dir, "dark-mirror-profile.md")
    with open(filepath, 'w') as f:
        f.write(content)
    
    return filepath


def generate_constellation_bible(state: dict, output_dir: str) -> str:
    """Generate the Constellation Bible."""
    answers = state["answers"]["doc7_constellation_bible"]
    project = state["project"]
    protagonist = answers.get("protagonist", {})
    characters = answers.get("characters", [])

    content = f"""# Constellation Bible: {project['title']}

## The Protagonist

### Before the Forging
**Unlit State:** {get_answer(protagonist, 'unlit_state', '[state]')}
**The Wound:** {get_answer(protagonist, 'wound', '[wound]')}
**The Lie:** {get_answer(protagonist, 'lie', '[lie]')}

### After the Forging
**Forged State:** To be determined through the story.
**The Truth:** What they learn.
**The Blade's Purpose:** From Crucible Thesis.

---

## Core Constellation

Characters and their roles in the constellation.

"""

    for char in characters:
        content += f"""### {get_answer(char, 'name', 'Character')}
**Role:** {get_answer(char, 'role', '[role]')}
**Relationship:** {get_answer(char, 'relationship', '[relationship]')}

"""
    
    filepath = os.path.join(output_dir, "constellation-bible.md")
    with open(filepath, 'w') as f:
        f.write(content)
    
    return filepath


def generate_mercy_ledger(state: dict, output_dir: str) -> str:
    """Generate the Mercy Ledger."""
    answers = state["answers"]["doc8_mercy_ledger"]
    project = state["project"]

    content = f"""# Mercy Ledger: {project['title']}

## Mercy Engine Overview

The protagonist shows mercy four times. Each costs something. Each pays off unexpectedly.

"""

    mercies = [
        ("mercy_1", "1", "The Seed", "II", "IV"),
        ("mercy_2", "2", "The Investment", "III", "V"),
        ("mercy_3", "3", "The Risk", "IV", "V"),
        ("mercy_4", "4", "The Impossible Gift", "V", "V (Beat 31)")
    ]

    for key, num, name, movement, payoff_movement in mercies:
        mercy = answers.get(key, {})
        content += f"""---

## Mercy {num}: {name}

**Movement:** {movement}
**Recipient:** {get_answer(mercy, 'recipient', '[recipient]')}

**The Merciful Act:** {get_answer(mercy, 'act', '[act]')}

**Immediate Cost:** {get_answer(mercy, 'cost', '[cost]')}

**Later Payoff (Movement {payoff_movement}):** {get_answer(mercy, 'payoff', '[payoff]')}

"""
    
    content += """---

## The Unexpected Agents (Beat 31)

At the moment of three failures, those who received mercy act to enable victory.
"""
    
    filepath = os.path.join(output_dir, "mercy-ledger.md")
    with open(filepath, 'w') as f:
        f.write(content)
    
    return filepath


def generate_world_forge(state: dict, output_dir: str) -> str:
    """Generate the World Forge document."""
    answers = state["answers"]["doc9_world_forge"]
    project = state["project"]

    content = f"""# World Forge: {project['title']}

## The World's Wound

**The Cosmic Problem:** {get_answer(answers, 'world_wound', '[wound]')}

---

## Power System

### Source
**Type:** {get_answer(answers, 'power_source', '[source]')}

### Limitations
{get_answer(answers, 'power_limitations', '[limitations]')}

### Previous Wielders
{get_answer(answers, 'previous_wielders', '[previous wielders]')}

---

## World-Protagonist Mirror

How the world reflects the protagonist's state:
{get_answer(answers, 'world_mirror', '[mirror description]')}

---

## Key Locations

{get_answer(answers, 'key_locations', '[locations]')}

---

## History

**Relevant Past:** {get_answer(answers, 'history', '[history]')}
"""
    
    filepath = os.path.join(output_dir, "world-forge.md")
    with open(filepath, 'w') as f:
        f.write(content)
    
    return filepath


def generate_summary(state: dict, output_dir: str) -> str:
    """Generate the Crucible Summary Card."""
    project = state["project"]
    thesis = state["answers"]["doc1_crucible_thesis"]

    content = f"""# Crucible Summary: {project['title']}

## The Forging Question
> Can the protagonist become {get_answer(thesis, 'forging_become', '[X]')} without becoming {get_answer(thesis, 'dark_mirror_represents', '[Y]')}?

## Three Strands
- **Quest:** {get_answer(thesis, 'burden_type', '[burden]')}
- **Fire:** {get_answer(thesis, 'fire_type', '[fire]')}
- **Constellation:** {get_answer(thesis, 'core_bond_type', '[bond]')}

## Four Forge Points + Apex
1. **Ignition:** Ordinary world destroyed
2. **First Crucible:** First impossible choice
3. **Second Crucible:** Trust broken
4. **Third Crucible:** Someone dies
5. **Apex:** Willed surrender

## The Mercy Engine
Four costly mercies enable victory through unexpected agents.

## Theme
> {get_answer(thesis, 'theme', '[theme]')}

## The Blade's Purpose
{get_answer(thesis, 'blade_purpose', '[purpose]')}

---

*Generated by Crucible Planner*
*{datetime.now().strftime('%Y-%m-%d')}*
"""
    
    filepath = os.path.join(output_dir, "crucible-summary.md")
    with open(filepath, 'w') as f:
        f.write(content)
    
    return filepath


def compile_all(project_path: str) -> list:
    """Compile all documents."""
    state = load_state(project_path)
    output_dir = os.path.join(project_path, "planning")
    
    # Ensure directories exist
    os.makedirs(os.path.join(output_dir, "strand-maps"), exist_ok=True)
    os.makedirs(os.path.join(output_dir, "forge-points"), exist_ok=True)
    
    files_created = []
    
    print("Generating planning documents...")

    files_created.append(generate_crucible_thesis(state, output_dir))
    print("   [OK] Crucible Thesis")

    files_created.append(generate_quest_strand(state, output_dir))
    print("   [OK] Quest Strand Map")

    files_created.append(generate_fire_strand(state, output_dir))
    print("   [OK] Fire Strand Map")

    files_created.append(generate_constellation_strand(state, output_dir))
    print("   [OK] Constellation Strand Map")

    files_created.extend(generate_forge_points(state, output_dir))
    print("   [OK] Forge Point Blueprints (5)")

    files_created.append(generate_dark_mirror(state, output_dir))
    print("   [OK] Dark Mirror Profile")

    files_created.append(generate_constellation_bible(state, output_dir))
    print("   [OK] Constellation Bible")

    files_created.append(generate_mercy_ledger(state, output_dir))
    print("   [OK] Mercy Ledger")

    files_created.append(generate_world_forge(state, output_dir))
    print("   [OK] World Forge")

    files_created.append(generate_summary(state, output_dir))
    print("   [OK] Summary Card")

    print(f"\n[OK] Generated {len(files_created)} documents in {output_dir}")
    
    return files_created


def main():
    if len(sys.argv) < 2:
        print("Usage: python compile_documents.py <project_path>")
        sys.exit(1)
    
    project_path = sys.argv[1]
    
    try:
        files = compile_all(project_path)
        print("\nFiles created:")
        for f in files:
            print(f"   {f}")
    except FileNotFoundError as e:
        print(f"[ERROR] {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
