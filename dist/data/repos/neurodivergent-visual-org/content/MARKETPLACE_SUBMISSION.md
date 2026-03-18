# Community Marketplace Submission Guide

## Status: Ready for Submission

Your plugin is now marketplace-ready and can be submitted to the `jeremylongshore/claude-code-plugins-plus` community marketplace.

## Pre-Submission Checklist

✅ **Valid plugin.json** - Complete with v3.1.1 metadata
✅ **marketplace.json** - Added for discoverability
✅ **Comprehensive README.md** - Installation, usage, features documented
✅ **LICENSE file** - MIT license included
✅ **Tested functionality** - Skill loads correctly in Claude Code
✅ **Security compliance** - No hardcoded credentials
✅ **Version control** - Pushed to GitHub repository

## Plugin Details

- **Name**: neurodivergent-visual-org
- **Category**: productivity
- **Version**: 3.1.1
- **Repository**: https://github.com/JackReis/neurodivergent-visual-org
- **License**: MIT
- **Components**: 1 skill

## Installation Methods

Users can install via:

1. **Direct marketplace** (recommended):
   ```
   /plugin marketplace add JackReis/neurodivergent-visual-org
   ```

2. **Community marketplace** (after submission approval):
   ```
   /plugin marketplace add jeremylongshore/claude-code-plugins-plus
   ```
   Then browse and install via `/plugin` menu

3. **Manual installation**:
   ```bash
   cd ~/.claude/plugins/
   git clone https://github.com/JackReis/neurodivergent-visual-org.git
   ```

## Submission Process

### Option A: Manual PR Submission

1. Fork https://github.com/jeremylongshore/claude-code-plugins-plus

2. Clone your fork:
   ```bash
   git clone https://github.com/YOUR_USERNAME/claude-code-plugins-plus.git
   cd claude-code-plugins-plus
   ```

3. Create feature branch:
   ```bash
   git checkout -b add-neurodivergent-visual-org
   ```

4. Copy plugin to productivity category:
   ```bash
   cp -r /path/to/neurodivergent-visual-org-plugin \
         plugins/productivity/neurodivergent-visual-org
   ```

5. Add entry to `.claude-plugin/marketplace.extended.json`:
   ```json
   {
     "name": "neurodivergent-visual-org",
     "source": "./plugins/productivity/neurodivergent-visual-org",
     "description": "Create ADHD-friendly visual organizational tools (Mermaid diagrams) optimized for neurodivergent thinking patterns with accessibility modes",
     "version": "3.1.1",
     "category": "productivity",
     "keywords": [
       "adhd",
       "neurodivergent",
       "visual-organization",
       "mermaid",
       "accessibility",
       "task-management",
       "decision-making",
       "executive-function",
       "productivity",
       "cognitive-accessibility"
     ],
     "author": {
       "name": "Jack Reis",
       "email": "hello@jack.digital"
     },
     "components": {
       "skills": 1
     }
   }
   ```

6. Run sync script:
   ```bash
   pnpm run sync-marketplace
   ```

7. Commit and push:
   ```bash
   git add .
   git commit -m "feat(productivity): add neurodivergent-visual-org plugin v3.1.1"
   git push origin add-neurodivergent-visual-org
   ```

8. Create pull request on GitHub with description:
   ```
   ## Plugin: neurodivergent-visual-org

   **Category**: Productivity
   **Version**: 3.1.1
   **License**: MIT

   ### Description
   ADHD-friendly visual organizational tools using Mermaid diagrams with
   adaptive cognitive modes and accessibility features.

   ### Features
   - 22 Mermaid diagram types for different cognitive needs
   - Neurodivergent/Neurotypical mode system with auto-detection
   - Colorblind-safe and monochrome accessibility modes
   - Research-backed design principles (ADHD neuroscience, cognitive load theory)
   - Task breakdown, decision-making, time management support

   ### Components
   - 1 skill: neurodivergent-visual-org

   ### Testing
   - ✅ Skill loads correctly in Claude Code
   - ✅ All documentation complete
   - ✅ MIT licensed
   - ✅ No security issues
   ```

### Option B: Automated Submission (Future)

The repository may add automated submission workflows in the future. Check their CONTRIBUTING.md for updates.

## Post-Submission

Once the PR is approved and merged:

1. Users can install via the community marketplace
2. Your plugin will be listed in the catalog of 254+ plugins
3. Plugin will be discoverable via `/plugin` menu in Claude Code

## Direct Usage (No Submission Required)

Your plugin is **already usable** without community marketplace approval via:

```
/plugin marketplace add JackReis/neurodivergent-visual-org
```

This command points directly to your GitHub repository, allowing immediate installation by anyone who knows the repository name.

## Repository Maintenance

Keep your standalone repository updated:
- Tag releases with semantic versioning (v3.1.1, v3.2.0, etc.)
- Update marketplace.json version numbers
- Document changes in README.md
- Maintain changelog for user visibility

---

**Created**: 2025-11-23
**Plugin Repository**: https://github.com/JackReis/neurodivergent-visual-org
**Marketplace**: https://github.com/jeremylongshore/claude-code-plugins-plus
