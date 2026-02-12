# Awesome Claude Skills

A comprehensive, well-documented collection of Claude skills covering development, business automation, content creation, and data analysis. Each skill includes usage examples, configuration guides, and best practices.

## Repository Details

**Name:** awesome-claude-skills

**Description:** A comprehensive, well-documented collection of Claude skills covering development, business automation, content creation, and data analysis. Each skill includes usage examples, configuration guides, and best practices.

**Topics/Tags:**
- claude
- claude-ai
- ai-skills
- mcp-servers
- automation
- developer-tools
- awesome-list
- productivity

---

## Installed Skills

### 1. Content Trend Researcher

**Location:** `.claude/skills/content-trend-researcher/`

**Description:** Advanced content and topic research skill that analyzes trends across Google Analytics, Google Trends, Substack, Medium, Reddit, LinkedIn, X, blogs, podcasts, and YouTube to generate data-driven article outlines based on user intent analysis.

**Capabilities:**
- Multi-platform trend analysis across 10+ platforms
- User intent analysis (informational, commercial, transactional, navigational)
- Content gap identification
- SEO-optimized article outline generation
- Platform-specific publishing strategies

**Source:** [alirezarezvani/claude-code-skill-factory](https://github.com/alirezarezvani/claude-code-skill-factory/tree/main/generated-skills/content-trend-researcher)

---

### 2. Frontend Dev Guidelines

**Location:** `.claude/skills/frontend-dev-guidelines/`

**Description:** Frontend development guidelines for React/TypeScript applications. Modern patterns including Suspense, lazy loading, useSuspenseQuery, file organization with features directory, MUI v7 styling, TanStack Router, performance optimization, and TypeScript best practices.

**Capabilities:**
- React component patterns with TypeScript
- Suspense-based data fetching with TanStack Query
- Lazy loading and code splitting strategies
- MUI v7 styling patterns and best practices
- TanStack Router implementation
- Performance optimization techniques
- File organization and project structure
- TypeScript standards and type safety
- Loading and error state management
- Complete working examples

**When to Use:**
- Creating new React components or pages
- Building new features with proper structure
- Setting up data fetching with TanStack Query
- Implementing routing with TanStack Router
- Styling components with MUI v7
- Performance optimization
- Organizing frontend code architecture

**Source:** [diet103/claude-code-infrastructure-showcase](https://github.com/diet103/claude-code-infrastructure-showcase/tree/main/.claude/skills/frontend-dev-guidelines)

**GitHub Stats:** ⭐ 3,532 stars | 🍴 480 forks

---

### 3. Blog Post Writer

**Location:** `.claude/skills/blog-post-writer/`

**Description:** Transform brain dumps into polished blog posts in Nick Nisi's voice. Use when the user wants to write a blog post with scattered ideas, talking points, and conclusions that need organization into a cohesive narrative with Nick's conversational, authentic, and thoughtful tone.

**Capabilities:**
- Transform unstructured brain dumps into polished blog posts
- Apply Nick Nisi's conversational and authentic writing style
- Organize scattered ideas into coherent narratives
- Follow Story Circle narrative framework when appropriate
- Balance vulnerability and technical depth
- Create engaging openings and conclusions
- Use specific examples and real details
- Vary sentence and paragraph length for rhythm
- Include technical content naturally (code, commands, tools)

**Writing Style Features:**
- Conversational yet substantive tone
- Vulnerable and authentic voice
- Journey-based narratives (problem → experience → resolution)
- Mix of short punchy sentences and longer explanations
- Self-aware humor
- Specific technical details without over-explaining
- Honest about limitations and uncertainty

**When to Use:**
- Writing blog posts from scattered thoughts
- Organizing technical content into readable narratives
- Creating developer-focused content
- Transforming notes into polished articles
- Applying a consistent voice and tone
- Structuring personal or technical stories

**Source:** [nicknisi/dotfiles](https://github.com/nicknisi/dotfiles/tree/main/home/.claude/skills/blog-post-writer)

**GitHub Stats:** ⭐ 2,890 stars | 🍴 371 forks

---

### 4. Writing Effective Prompts

**Location:** `.claude/skills/writing-effective-prompts/`

**Description:** Structure Claude prompts for clarity and better results using roles, explicit instructions, context, positive framing, and strategic organization. Use when crafting prompts for complex tasks, long documents, tool workflows, or code generation.

**Capabilities:**
- Role-based prompt structuring for behavioral context
- Explicit instruction formatting for specific requirements
- Context and motivation explanations for better generalization
- Positive framing (what TO do vs what NOT to do)
- Example alignment with `<example>` tags
- XML tag structuring for complex multi-section outputs
- Task chaining for multi-step processes
- Long context optimization (document placement strategies)
- Tool usage strategy and parallelization guidance
- Code generation quality enhancement techniques

**Core Principles:**
1. **Start with a Role** - Set behavioral context upfront
2. **Be Explicit** - Replace vague requests with specific requirements
3. **Add Context** - Explain why to help Claude generalize
4. **Use Positive Framing** - Tell what TO do, not what NOT to do
5. **Provide Examples** - Use `<example>` tags to show desired output
6. **Structure Output** - Use XML tags for complex responses
7. **Chain Tasks** - Break complex processes into phases
8. **Optimize Context** - Put large documents first, queries last

**When to Use:**
- Crafting prompts for complex tasks
- Working with long documents (20K+ tokens)
- Defining tool workflows
- Code generation projects
- Multi-step processes
- Multi-document analysis
- Improving prompt clarity and results

**Source:** [CaptainCrouton89/.claude](https://github.com/CaptainCrouton89/.claude/tree/main/skills/prompting-effectively)

**GitHub Stats:** ⭐ 479 stars | 🍴 68 forks

---

## Installation

### Prerequisites
- Claude Code CLI installed
- Python 3.8+ (for Python-based skills)
- Git (optional, for cloning updates)

### Setup

The skills are already installed in this repository. To use them with Claude Code:

1. **Navigate to this directory:**
   ```bash
   cd /Volumes/AI\ Engineer/awesome-claude-skills
   ```

2. **Verify skill installation:**
   ```bash
   ls -la .claude/skills/
   ```

3. **Start using skills with Claude Code:**
   ```bash
   claude
   ```

---

## How to Use Skills

### Using Content Trend Researcher

**Basic Usage:**
```
@content-trend-researcher

Topic: AI automation for small businesses
Platforms: Google Trends, Reddit, LinkedIn, YouTube
```

**Advanced Usage with Detailed Parameters:**
```
@content-trend-researcher

Topic: "Sustainable fashion for millennials"
Platforms: ALL
Intent focus: informational
Target audience: Environmentally conscious millennials aged 25-35
Content type: blog
Analysis depth: deep
Number of outlines: 3
```

**Quick Topic Research:**
```
@content-trend-researcher

Topic: "Remote work productivity tools"
Platforms: Google Trends, YouTube, LinkedIn
Analysis depth: quick
Number of outlines: 1
```

### Output Format

The skill returns comprehensive research including:
- **Topic Overview:** Search volume, trend direction, competition level
- **Platform Insights:** Platform-specific best practices and metrics
- **User Intent Analysis:** Intent breakdown and top questions
- **Content Gaps:** Underserved topics with opportunity scores
- **Article Outlines:** Complete SEO-optimized outlines with structure
- **Recommendations:** Publishing schedule and promotion strategy

For detailed usage examples, see: `.claude/skills/content-trend-researcher/HOW_TO_USE.md`

---

### Using Frontend Dev Guidelines

**Basic Usage - Create a Component:**
```
@frontend-dev-guidelines

I need to create a user profile component that fetches user data and displays it with proper loading states.
```

**Create a New Feature:**
```
@frontend-dev-guidelines

I'm building a new "posts" feature with:
- API service layer for fetching posts
- Components for listing and viewing posts
- Proper TypeScript types
- Routing with TanStack Router
- Lazy loading and Suspense
```

**Data Fetching Pattern:**
```
@frontend-dev-guidelines

Show me the correct pattern for fetching data using useSuspenseQuery with TanStack Query.
Include error handling and cache management.
```

**Styling Guidance:**
```
@frontend-dev-guidelines

I need to style a complex form layout with MUI v7.
Should I use inline styles or a separate file? Show me the pattern.
```

**Performance Optimization:**
```
@frontend-dev-guidelines

My component is re-rendering too often.
Show me how to optimize it with useMemo, useCallback, and React.memo.
```

### Quick Reference

The skill includes comprehensive guides for:
- **Component Patterns** - Modern React component structure
- **Data Fetching** - useSuspenseQuery with cache strategies
- **File Organization** - Features vs components structure
- **Styling** - MUI v7 patterns and best practices
- **Routing** - TanStack Router setup
- **Loading States** - Suspense boundaries and error handling
- **Performance** - Optimization techniques
- **TypeScript** - Type safety standards
- **Common Patterns** - Forms, auth, DataGrid wrappers
- **Complete Examples** - Full working code samples

For detailed documentation, see: `.claude/skills/frontend-dev-guidelines/SKILL.md`

---

### Using Blog Post Writer

**Basic Usage - Transform Brain Dump:**
```
@blog-post-writer

Here are my scattered thoughts on using AI coding tools:
- cursor is in IDE, feels familiar
- but claude code is in terminal, my natural environment
- tried cursor first, felt weird leaving vim
- claude code met me where I was
- not about which is better, about workflow fit
- some devs love IDE integration
- I need terminal access
- conclusion: use what fits YOUR workflow
```

**Technical Blog Post:**
```
@blog-post-writer

Brain dump on my experience with TanStack Query:
- Started using it last month
- Coming from Redux, felt weird at first
- Cache invalidation patterns are different
- useSuspenseQuery changed my approach to loading states
- No more isLoading checks everywhere
- Performance improved dramatically
- Would recommend to anyone doing data fetching
```

**Career/Personal Story:**
```
@blog-post-writer

Notes from my first month at the new company:
- Impostor syndrome hit hard
- Everyone seemed so experienced
- Pair programming helped break the ice
- Made first PR contribution after week 2
- Code review feedback was constructive
- Starting to feel like I belong
- Key learning: everyone feels uncertain sometimes
```

**Journey-Based Narrative:**
```
@blog-post-writer

My path from skeptic to AI tool advocate:
- Thought AI would replace developers
- Tried Claude Code reluctantly
- Watched it debug my code using the same tools I use
- Realized it's augmentation, not replacement
- Now using it daily for routine tasks
- Freed up time for creative problem-solving
- Perspective: we're not being replaced, we're being amplified
```

### What You'll Get

The skill will:
1. **Read voice-tone.md** - Understand Nick's writing style
2. **Check story-circle.md** - Identify narrative opportunities
3. **Organize content** - Structure ideas into coherent sections
4. **Apply voice characteristics** - Write with conversational authenticity
5. **Review and refine** - Ensure quality and narrative flow

**Output includes:**
- Polished blog post with engaging opening
- Clear narrative arc (problem → journey → resolution)
- Technical details with specific examples
- Vulnerability and authenticity
- Forward-looking conclusion
- Varied paragraph lengths for rhythm

For detailed documentation, see: `.claude/skills/blog-post-writer/SKILL.md`

---

### Using Writing Effective Prompts

**Role-Based Prompting:**
```
@writing-effective-prompts

I need to create a prompt for generating React components. How should I structure it?

Task: Generate a user profile card component
Requirements:
- TypeScript
- MUI v7 styling
- Responsive design
- Loading states
```

**Complex Task Structuring:**
```
@writing-effective-prompts

Help me structure a prompt for a multi-step data migration task:
- Read from legacy PostgreSQL database
- Transform data structures
- Validate and clean data
- Write to new MongoDB schema
- Generate migration report
```

**Long Document Analysis:**
```
@writing-effective-prompts

I need to analyze 3 large documents (annual report, competitor analysis, market research). How should I structure the prompt for best results?
```

**Code Generation Quality:**
```
@writing-effective-prompts

Improve this prompt for better code generation:
"Create a dashboard"

I want:
- Real-time data
- Interactive filters
- Charts and graphs
- Export functionality
```

**Tool Workflow Optimization:**
```
@writing-effective-prompts

Help me write a prompt for Claude Code to:
- Search multiple files for a pattern
- Update matching code
- Run tests
- Create a PR

How should I structure this for parallel execution?
```

### Key Techniques

**1. Start with a Role:**
```
You are an expert software test engineer.
Help me write comprehensive unit tests.
```

**2. Be Explicit with Instructions:**
```
Create a dashboard with:
- Real-time data visualization
- Interactive filtering
- Export to CSV/PDF
Include as many relevant features as possible.
```

**3. Add Context:**
```
Response will be read via text-to-speech,
so avoid ellipses and use complete sentences.
```

**4. Use XML Tags:**
```
<code_quality>
Assess overall quality
</code_quality>

<security_review>
Review security concerns
</security_review>
```

**5. Provide Examples:**
```
<example>
Input: "Added JWT auth"
Output: feat(auth): implement JWT authentication
</example>
```

For detailed documentation, see: `.claude/skills/writing-effective-prompts/SKILL.md`

---

## Skill Documentation

Each skill includes:
- **SKILL.md** - Complete skill documentation with capabilities and examples
- **HOW_TO_USE.md** - Step-by-step usage guide with invocation patterns
- **sample_input.json** - Example input format
- **expected_output.json** - Example output format
- Python modules for skill functionality

---

## Directory Structure

```
awesome-claude-skills/
├── README.md                                    # This file
├── QUICK_START.md                               # Quick start guide
├── .gitignore                                   # Git ignore file
├── .claude/
│   └── skills/
│       ├── content-trend-researcher/
│       │   ├── SKILL.md                        # Skill documentation
│       │   ├── HOW_TO_USE.md                   # Usage guide
│       │   ├── sample_input.json               # Sample inputs
│       │   ├── expected_output.json            # Sample outputs
│       │   ├── intent_analyzer.py              # Intent analysis module
│       │   ├── outline_generator.py            # Outline generation module
│       │   ├── platform_insights.py            # Platform analysis module
│       │   └── trend_analyzer.py               # Trend analysis module
│       ├── frontend-dev-guidelines/
│       │   ├── SKILL.md                        # Main skill documentation
│       │   └── resources/                      # Modular resource files
│       │       ├── common-patterns.md          # Forms, auth, DataGrid patterns
│       │       ├── complete-examples.md        # Full working examples
│       │       ├── component-patterns.md       # React component structure
│       │       ├── data-fetching.md            # TanStack Query patterns
│       │       ├── file-organization.md        # Project structure guide
│       │       ├── loading-and-error-states.md # Suspense and error handling
│       │       ├── performance.md              # Optimization techniques
│       │       ├── routing-guide.md            # TanStack Router setup
│       │       ├── styling-guide.md            # MUI v7 styling patterns
│       │       └── typescript-standards.md     # TypeScript best practices
│       ├── blog-post-writer/
│       │   ├── SKILL.md                        # Main skill documentation
│       │   └── references/                     # Writing style references
│       │       ├── voice-tone.md               # Nick Nisi's voice and tone guide
│       │       └── story-circle.md             # Story Circle narrative framework
│       └── writing-effective-prompts/
│           └── SKILL.md                        # Prompt engineering guide
└── .idea/                                       # IDE configuration
```

---

## Adding More Skills

To add more skills from the [claude-code-skill-factory](https://github.com/alirezarezvani/claude-code-skill-factory) repository:

1. **Browse available skills:**
   - Visit: https://github.com/alirezarezvani/claude-code-skill-factory/tree/main/generated-skills

2. **Download a skill:**
   ```bash
   # Example: Adding the 'technical-writing' skill
   cd .claude/skills
   curl -L https://github.com/alirezarezvani/claude-code-skill-factory/archive/main.tar.gz | tar xz --strip=3 "claude-code-skill-factory-main/generated-skills/SKILL_NAME"
   ```

3. **Or use the installation script:**
   ```bash
   # Coming soon: install-skill.sh
   ./install-skill.sh SKILL_NAME
   ```

---

## Tips for Best Results

### For Content Trend Researcher:

1. **Be Specific with Topics**
   - Bad: "Marketing"
   - Good: "Email marketing automation for e-commerce businesses"

2. **Choose Relevant Platforms**
   - Only select platforms where your target audience consumes content

3. **Define Your Audience**
   - Include demographic details
   - Specify industry or niche
   - Mention experience level

4. **Specify Analysis Depth**
   - **quick** (15-20 min): High-level overview, 1 outline
   - **standard** (30-40 min): Comprehensive analysis, 1-2 outlines
   - **deep** (45-60 min): Complete competitive analysis, 3+ outlines

### For Frontend Dev Guidelines:

1. **Use the Skill When:**
   - Creating new React components or pages
   - Building features with proper structure
   - Need data fetching patterns
   - Setting up routing
   - Styling with MUI v7
   - Optimizing performance

2. **Follow the Core Principles:**
   - Lazy load heavy components (DataGrid, charts, editors)
   - Use Suspense for loading states (no early returns)
   - Use useSuspenseQuery for data fetching
   - Organize features with api/, components/, hooks/, helpers/ subdirectories
   - Use import aliases: @/, ~types, ~components, ~features

3. **Quick Checklist for Components:**
   - Use `React.FC<Props>` pattern
   - Lazy load if heavy component
   - Wrap in `<SuspenseLoader>`
   - Use `useSuspenseQuery` for data
   - No early returns with loading spinners
   - Use `useCallback` for event handlers

4. **Navigation Guide:**
   - Creating component → Read `resources/component-patterns.md`
   - Fetching data → Read `resources/data-fetching.md`
   - File organization → Read `resources/file-organization.md`
   - Styling → Read `resources/styling-guide.md`
   - Performance → Read `resources/performance.md`
   - Full examples → Read `resources/complete-examples.md`

### For Blog Post Writer:

1. **Provide Brain Dumps:**
   - Don't worry about organization
   - Include scattered thoughts, bullet points, random ideas
   - Add code examples, commands, or technical details
   - Include conclusions or takeaways

2. **Writing Style Applied:**
   - Conversational and authentic voice
   - Journey-based narratives (problem → experience → resolution)
   - Vulnerability and honesty about uncertainty
   - Specific examples with real details
   - Mix of short and long sentences

3. **Content Types:**
   - Technical blog posts (tools, frameworks, coding patterns)
   - Career/personal stories (onboarding, learning, growth)
   - Tool comparisons and reviews
   - Journey from skepticism to understanding
   - Problem-solving narratives

4. **Voice Characteristics:**
   - Writes like talking to a peer over coffee
   - Admits uncertainty or being wrong
   - Uses specific examples (tool names, commands, numbers)
   - Shows vulnerability appropriately
   - Self-aware humor, not forced
   - Ends with forward momentum

### For Writing Effective Prompts:

1. **Always Start with a Role:**
   - Define Claude's expertise upfront
   - Example: "You are an expert software test engineer"
   - Sets behavioral context for better responses

2. **Be Explicit, Not Vague:**
   - Bad: "Create a dashboard"
   - Good: "Create a dashboard with real-time data visualization, interactive filtering, responsive design, and export functionality"
   - List specific requirements

3. **Use Positive Framing:**
   - Tell what TO do, not what NOT to do
   - Bad: "Don't use complex language"
   - Good: "Use clear, simple language"

4. **Add Context and Motivation:**
   - Explain WHY requirements matter
   - Example: "Response will be read via TTS, so use complete sentences"
   - Helps Claude generalize better

5. **Leverage XML Tags:**
   - For complex multi-section outputs
   - Use `<code_quality>`, `<security_review>`, etc.
   - Structures responses clearly

6. **Optimize Long Context:**
   - Put large documents (20K+ tokens) FIRST
   - Place queries LAST
   - Can improve quality by up to 30%

7. **Provide Aligned Examples:**
   - Use `<example>` tags
   - Show exact desired output format
   - Examples powerfully shape responses

---

## Platform Combinations

### For Thought Leadership
```
Platforms: LinkedIn, Medium, Substack, Podcasts
```

### For Technical Content
```
Platforms: Blogs, YouTube, Reddit, Medium
```

### For Viral/Engagement
```
Platforms: X, LinkedIn, YouTube, TikTok
```

### For Community Building
```
Platforms: Reddit, Substack, Podcasts, YouTube
```

---

## Troubleshooting

### Skill Not Found
- Ensure the skill is in `.claude/skills/` directory
- Check that SKILL.md file exists
- Restart Claude Code if needed

### Analysis Too Generic
- Be more specific with your topic and target audience
- Provide detailed audience demographics
- Specify intent focus

### Need More Outlines
- Explicitly request the number: `Number of outlines: 5`

### Platform Recommendations Irrelevant
- Only select 3-4 most relevant platforms for your audience
- Or ask skill to help identify best platforms first

---

## Common Workflows

### Weekly Newsletter Research
```
Monday: Research topic with @content-trend-researcher
Tuesday: Review outlines and select best angle
Wednesday: Write newsletter draft
Thursday: Create social promotion snippets
Friday: Schedule and publish
```

### Monthly Content Calendar
```
Week 1: Research 4 topics (one per week)
Week 2: Draft all articles
Week 3: Optimize and add multimedia
Week 4: Schedule publication and promotion
```

### React Feature Development
```
Day 1: Use @frontend-dev-guidelines to set up feature structure
       Create api/, components/, hooks/, helpers/, types/ directories
Day 2: Implement API service layer with TypeScript types
Day 3: Build components with Suspense and lazy loading
Day 4: Add routing with TanStack Router
Day 5: Performance optimization and testing
```

### Component Refactoring
```
Step 1: @frontend-dev-guidelines - Review current component
Step 2: Identify performance issues (unnecessary re-renders)
Step 3: Apply optimization patterns (useMemo, useCallback, React.memo)
Step 4: Refactor data fetching to useSuspenseQuery
Step 5: Update loading states with Suspense boundaries
Step 6: Test and measure performance improvements
```

### Technical Blog Writing
```
Step 1: Collect scattered notes and thoughts throughout the week
Step 2: @blog-post-writer [paste brain dump with technical details]
Step 3: Review generated post for accuracy and voice
Step 4: Request revisions if needed
Step 5: Publish to blog platform
Step 6: Share on LinkedIn, X, or dev communities
```

### Content Creation Pipeline
```
Monday: Research topic with @content-trend-researcher
Tuesday: Collect notes and personal experiences
Wednesday: @blog-post-writer - Transform notes into blog post
Thursday: Build accompanying demo with @frontend-dev-guidelines
Friday: Publish blog post with demo link
```

### Prompt Engineering Workflow
```
Step 1: Draft initial prompt for complex task
Step 2: @writing-effective-prompts - Improve prompt structure
Step 3: Add role definition and context
Step 4: Include specific examples with <example> tags
Step 5: Test refined prompt
Step 6: Iterate based on results
```

### Code Generation with Quality
```
Step 1: @writing-effective-prompts - Structure code generation prompt
Step 2: Add role: "You are an expert React developer"
Step 3: Define explicit requirements with XML tags
Step 4: @frontend-dev-guidelines - Generate code
Step 5: Review and refine
```

---

## Contributing

To contribute skills to this collection:

1. Fork this repository
2. Add your skill to `.claude/skills/`
3. Include proper documentation (SKILL.md, HOW_TO_USE.md)
4. Test the skill thoroughly
5. Submit a pull request

---

## Resources

- [Claude Code Documentation](https://docs.anthropic.com/claude/docs/claude-code)
- [Claude Code Skill Factory](https://github.com/alirezarezvani/claude-code-skill-factory)
- [Skills Marketplace](https://skillsmp.com)
- [Claude API Documentation](https://docs.anthropic.com/claude/reference/getting-started-with-the-api)

---

## License

This repository follows the licensing of individual skills. Please check each skill's source repository for specific license information.

---

## Support

For issues or questions:
1. Review the skill's SKILL.md documentation
2. Check HOW_TO_USE.md for usage examples
3. Review sample_input.json and expected_output.json
4. Open an issue on GitHub

---

## Acknowledgments

- Skills sourced from:
  - [claude-code-skill-factory](https://github.com/alirezarezvani/claude-code-skill-factory) by alirezarezvani
  - [claude-code-infrastructure-showcase](https://github.com/diet103/claude-code-infrastructure-showcase) by diet103
  - [nicknisi/dotfiles](https://github.com/nicknisi/dotfiles) by nicknisi
  - [CaptainCrouton89/.claude](https://github.com/CaptainCrouton89/.claude) by CaptainCrouton89
- [Skills Marketplace](https://skillsmp.com) for skill discovery
- Anthropic for Claude AI and Claude Code

---

**Last Updated:** November 11, 2024

**Version:** 1.0.0
