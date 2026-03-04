# Question to State Key Mapping

Use this mapping to save answers to the correct state keys after each AskUserQuestion response.

## How to Use

After receiving an answer from AskUserQuestion:
1. Find the question header in the table below
2. Use the corresponding Document and State Key
3. Run the save command with **5 parameters** (full context):

```bash
python scripts/save_state.py <project_path> --answer <document> <state_key> "<question_text>" "<answer>" "<description>"
```

**Parameters:**
- `document` - Document key from table below
- `state_key` - State key from table below
- `question_text` - The full question you asked
- `answer` - The option label the user selected
- `description` - The description of the selected option

---

## Document 1: Crucible Thesis

**Document key:** `doc1_crucible_thesis`

| Question # | Header | State Key | Save Command |
|------------|--------|-----------|--------------|
| Q1.1 | Burden | `burden_type` | `--answer doc1_crucible_thesis burden_type` |
| Q1.2 | Fire | `fire_type` | `--answer doc1_crucible_thesis fire_type` |
| Q1.3 | Core Bond | `core_bond_type` | `--answer doc1_crucible_thesis core_bond_type` |
| Q1.4 | Dark Mirror | `dark_mirror_connection` | `--answer doc1_crucible_thesis dark_mirror_connection` |
| Q1.5 | Forging | `forging_become` | `--answer doc1_crucible_thesis forging_become` |
| Q1.5b | Dark Mirror Represents | `dark_mirror_represents` | `--answer doc1_crucible_thesis dark_mirror_represents` |
| Q1.6 | Antagonist Truth | `antagonist_truth` | `--answer doc1_crucible_thesis antagonist_truth` |
| Q1.7 | Surrender | `surrender_type` | `--answer doc1_crucible_thesis surrender_type` |
| Q1.8 | Theme | `theme` | `--answer doc1_crucible_thesis theme` |
| Q1.9 | Blade Purpose | `blade_purpose` | `--answer doc1_crucible_thesis blade_purpose` |
| Q1.10 | Verification | *(no save - confirmation only)* | `--complete 1` after confirmed |

---

## Document 2: Quest Strand Map

**Document key:** `doc2_quest_strand`

| Question # | Header | State Key | Save Command |
|------------|--------|-----------|--------------|
| Q2.1 | Burden Origin | `burden_origin` | `--answer doc2_quest_strand burden_origin` |
| Q2.2 | Why Protagonist | `why_protagonist` | `--answer doc2_quest_strand why_protagonist` |
| Q2.3 | Antagonist Stake | `antagonist_stake` | `--answer doc2_quest_strand antagonist_stake` |
| Q2.4 | Quest Escalation | `quest_escalation` | `--answer doc2_quest_strand quest_escalation` |
| Q2.5 | Impossible | `impossible_requirement` | `--answer doc2_quest_strand impossible_requirement` |
| Q2.6 | Resolution | `resolution_method` | `--answer doc2_quest_strand resolution_method` |
| Q2.7 | Verification | *(no save - confirmation only)* | `--complete 2` after confirmed |

---

## Document 3: Fire Strand Map

**Document key:** `doc3_fire_strand`

| Question # | Header | State Key | Save Command |
|------------|--------|-----------|--------------|
| Q3.1 | Manifestation | `fire_manifestation` | `--answer doc3_fire_strand fire_manifestation` |
| Q3.2 | Danger | `fire_danger` | `--answer doc3_fire_strand fire_danger` |
| Q3.3 | Cost | `cost_of_use` | `--answer doc3_fire_strand cost_of_use` |
| Q3.4 | Mastery Path | `mastery_path` | `--answer doc3_fire_strand mastery_path` |
| Q3.5 | Hardening | `hardening_line` | `--answer doc3_fire_strand hardening_line` |
| Q3.6 | Mastery Moment | `mastery_method` | `--answer doc3_fire_strand mastery_method` |
| Q3.7 | Verification | *(no save - confirmation only)* | `--complete 3` after confirmed |

---

## Document 4: Constellation Strand Map

**Document key:** `doc4_constellation_strand`

| Question # | Header | State Key | Save Command |
|------------|--------|-----------|--------------|
| Q4.1 | Faithful | `faithful_companion` | `--answer doc4_constellation_strand faithful_companion` |
| Q4.2 | Sacrifice | `sacrifice_character` | `--answer doc4_constellation_strand sacrifice_character` |
| Q4.3 | Betrayal | `betrayal_source` | `--answer doc4_constellation_strand betrayal_source` |
| Q4.4 | Expansion | `constellation_expansion` | `--answer doc4_constellation_strand constellation_expansion` |
| Q4.5 | Bond Saves | `bond_that_saves` | `--answer doc4_constellation_strand bond_that_saves` |
| Q4.6 | Fate | `constellation_fate` | `--answer doc4_constellation_strand constellation_fate` |
| Q4.7 | Verification | *(no save - confirmation only)* | `--complete 4` after confirmed |

---

## Document 5: Forge Point Blueprints

### FP0: Ignition

**Document key:** `doc5_forge_points.fp0_ignition`

| Header | State Key | Save Command |
|--------|-----------|--------------|
| Quest Crisis | `quest_crisis` | `--answer doc5_forge_points.fp0_ignition quest_crisis` |
| Fire Crisis | `fire_crisis` | `--answer doc5_forge_points.fp0_ignition fire_crisis` |
| Constellation Crisis | `constellation_crisis` | `--answer doc5_forge_points.fp0_ignition constellation_crisis` |
| Sacrifice | `sacrifice` | `--answer doc5_forge_points.fp0_ignition sacrifice` |

### FP1: First Crucible

**Document key:** `doc5_forge_points.fp1_first`

| Header | State Key | Save Command |
|--------|-----------|--------------|
| Quest Crisis | `quest_crisis` | `--answer doc5_forge_points.fp1_first quest_crisis` |
| Fire Crisis | `fire_crisis` | `--answer doc5_forge_points.fp1_first fire_crisis` |
| Constellation Crisis | `constellation_crisis` | `--answer doc5_forge_points.fp1_first constellation_crisis` |
| Sacrifice | `sacrifice` | `--answer doc5_forge_points.fp1_first sacrifice` |

### FP2: Second Crucible

**Document key:** `doc5_forge_points.fp2_second`

| Header | State Key | Save Command |
|--------|-----------|--------------|
| Quest Crisis | `quest_crisis` | `--answer doc5_forge_points.fp2_second quest_crisis` |
| Fire Crisis | `fire_crisis` | `--answer doc5_forge_points.fp2_second fire_crisis` |
| Constellation Crisis | `constellation_crisis` | `--answer doc5_forge_points.fp2_second constellation_crisis` |
| Sacrifice | `sacrifice` | `--answer doc5_forge_points.fp2_second sacrifice` |

### FP3: Third Crucible

**Document key:** `doc5_forge_points.fp3_third`

| Header | State Key | Save Command |
|--------|-----------|--------------|
| Quest Crisis | `quest_crisis` | `--answer doc5_forge_points.fp3_third quest_crisis` |
| Fire Crisis | `fire_crisis` | `--answer doc5_forge_points.fp3_third fire_crisis` |
| Constellation Crisis | `constellation_crisis` | `--answer doc5_forge_points.fp3_third constellation_crisis` |
| Sacrifice | `sacrifice` | `--answer doc5_forge_points.fp3_third sacrifice` |

### Apex: Willed Surrender

**Document key:** `doc5_forge_points.apex`

| Header | State Key | Save Command |
|--------|-----------|--------------|
| Quest Crisis | `quest_crisis` | `--answer doc5_forge_points.apex quest_crisis` |
| Fire Crisis | `fire_crisis` | `--answer doc5_forge_points.apex fire_crisis` |
| Constellation Crisis | `constellation_crisis` | `--answer doc5_forge_points.apex constellation_crisis` |
| Sacrifice | `sacrifice` | `--answer doc5_forge_points.apex sacrifice` |

After all 5 Forge Points: `--complete 5`

---

## Document 6: Dark Mirror Profile

**Document key:** `doc6_dark_mirror`

| Question # | Header | State Key | Save Command |
|------------|--------|-----------|--------------|
| Q6.1 | Origin | `origin_parallel` | `--answer doc6_dark_mirror origin_parallel` |
| Q6.2 | Divergence | `divergence` | `--answer doc6_dark_mirror divergence` |
| Q6.3 | Want | `antagonist_want` | `--answer doc6_dark_mirror antagonist_want` |
| Q6.4 | Offer | `compelling_offer` | `--answer doc6_dark_mirror compelling_offer` |
| Q6.5 | Tempting | `why_tempting` | `--answer doc6_dark_mirror why_tempting` |
| Q6.6 | Cost | `hidden_cost` | `--answer doc6_dark_mirror hidden_cost` |
| Q6.7 | Defeat | `defeat_method` | `--answer doc6_dark_mirror defeat_method` |
| Q6.8 | End | `antagonist_end` | `--answer doc6_dark_mirror antagonist_end` |
| Q6.9 | Verification | *(no save - confirmation only)* | `--complete 6` after confirmed |

---

## Document 7: Constellation Bible

### Protagonist Profile

**Document key:** `doc7_constellation_bible.protagonist`

| Header | State Key | Save Command |
|--------|-----------|--------------|
| Unlit State | `unlit_state` | `--answer doc7_constellation_bible.protagonist unlit_state` |
| Wound | `wound` | `--answer doc7_constellation_bible.protagonist wound` |
| Lie | `lie` | `--answer doc7_constellation_bible.protagonist lie` |

### Character Entries

For each character, use `doc7_constellation_bible` with character data:
- `name` - Character name
- `role` - Their role in the constellation
- `relationship` - Their relationship to protagonist

**Note:** Characters are stored in an array. The save_state.py handles character additions automatically.

After all characters: `--complete 7`

---

## Document 8: Mercy Ledger

### Mercy 1

**Document key:** `doc8_mercy_ledger.mercy_1`

| Header | State Key | Save Command |
|--------|-----------|--------------|
| Recipient | `recipient` | `--answer doc8_mercy_ledger.mercy_1 recipient` |
| Act | `act` | `--answer doc8_mercy_ledger.mercy_1 act` |
| Cost | `cost` | `--answer doc8_mercy_ledger.mercy_1 cost` |
| Payoff | `payoff` | `--answer doc8_mercy_ledger.mercy_1 payoff` |

### Mercy 2

**Document key:** `doc8_mercy_ledger.mercy_2`

| Header | State Key | Save Command |
|--------|-----------|--------------|
| Recipient | `recipient` | `--answer doc8_mercy_ledger.mercy_2 recipient` |
| Act | `act` | `--answer doc8_mercy_ledger.mercy_2 act` |
| Cost | `cost` | `--answer doc8_mercy_ledger.mercy_2 cost` |
| Payoff | `payoff` | `--answer doc8_mercy_ledger.mercy_2 payoff` |

### Mercy 3

**Document key:** `doc8_mercy_ledger.mercy_3`

| Header | State Key | Save Command |
|--------|-----------|--------------|
| Recipient | `recipient` | `--answer doc8_mercy_ledger.mercy_3 recipient` |
| Act | `act` | `--answer doc8_mercy_ledger.mercy_3 act` |
| Cost | `cost` | `--answer doc8_mercy_ledger.mercy_3 cost` |
| Payoff | `payoff` | `--answer doc8_mercy_ledger.mercy_3 payoff` |

### Mercy 4

**Document key:** `doc8_mercy_ledger.mercy_4`

| Header | State Key | Save Command |
|--------|-----------|--------------|
| Recipient | `recipient` | `--answer doc8_mercy_ledger.mercy_4 recipient` |
| Act | `act` | `--answer doc8_mercy_ledger.mercy_4 act` |
| Cost | `cost` | `--answer doc8_mercy_ledger.mercy_4 cost` |
| Payoff | `payoff` | `--answer doc8_mercy_ledger.mercy_4 payoff` |

After all 4 mercies: `--complete 8`

---

## Document 9: World Forge

**Document key:** `doc9_world_forge`

| Question # | Header | State Key | Save Command |
|------------|--------|-----------|--------------|
| Q9.1 | World Wound | `world_wound` | `--answer doc9_world_forge world_wound` |
| Q9.2 | Power Source | `power_source` | `--answer doc9_world_forge power_source` |
| Q9.3 | Limitations | `power_limitations` | `--answer doc9_world_forge power_limitations` |
| Q9.4 | Previous | `previous_wielders` | `--answer doc9_world_forge previous_wielders` |
| Q9.5 | Locations | `key_locations` | `--answer doc9_world_forge key_locations` |
| Q9.6 | Mirror | `world_mirror` | `--answer doc9_world_forge world_mirror` |
| Q9.7 | Cultures | `cultures` | `--answer doc9_world_forge cultures` |
| Q9.8 | History | `history` | `--answer doc9_world_forge history` |
| Q9.9 | Verification | *(no save - confirmation only)* | `--complete 9` after confirmed |

---

## Quick Reference: Full Command Examples

**New 5-parameter format** - saves question, answer, AND description:

```bash
# Save an answer to Document 1 (Burden question)
python scripts/save_state.py "./crucible-project" --answer doc1_crucible_thesis burden_type "What form does the external burden take in your story?" "Physical object" "An artifact that must be destroyed or protected"

# Save Fire type
python scripts/save_state.py "./crucible-project" --answer doc1_crucible_thesis fire_type "What is the nature of the protagonist's Fire?" "Magical ability" "A power that corrupts with use"

# Save nested answer (Forge Point)
python scripts/save_state.py "./crucible-project" --answer doc5_forge_points.fp0_ignition quest_crisis "What Quest crisis emerges at Ignition?" "Artifact stolen" "The burden is taken by enemy forces"

# Save nested answer (Mercy)
python scripts/save_state.py "./crucible-project" --answer doc8_mercy_ledger.mercy_1 recipient "Who receives this mercy?" "Enemy soldier" "A combatant who surrendered"

# Save protagonist info
python scripts/save_state.py "./crucible-project" --answer doc7_constellation_bible.protagonist unlit_state "Who is the protagonist before Ignition?" "Restless but trapped" "Has potential but circumstances prevent growth"

# Mark document complete
python scripts/save_state.py "./crucible-project" --complete 1

# Update progress position
python scripts/save_state.py "./crucible-project" --progress 2 3
```

**Stored format:**
```json
{
  "question": "What form does the external burden take in your story?",
  "answer": "Physical object",
  "description": "An artifact that must be destroyed or protected"
}
```
