# Quick Start Guide

Get started with your Claude Skills in 5 minutes!

## Step 1: Verify Installation

Check that the skills are installed correctly:

```bash
ls -la .claude/skills/
```

You should see four skill directories:

**content-trend-researcher/**
- SKILL.md, HOW_TO_USE.md
- Python modules (intent_analyzer.py, outline_generator.py, etc.)
- sample_input.json, expected_output.json

**frontend-dev-guidelines/**
- SKILL.md
- resources/ directory with 10 guide files

**blog-post-writer/**
- SKILL.md
- references/ directory with voice-tone.md and story-circle.md

**writing-effective-prompts/**
- SKILL.md
- Prompt engineering best practices

## Step 2: Start Claude Code

Navigate to this directory and start Claude:

```bash
cd /Volumes/AI\ Engineer/awesome-claude-skills
claude
```

## Step 3: Use Your First Skill

### Content Trend Researcher

Try this simple example:

```
@content-trend-researcher

Topic: Python programming tutorials
Platforms: YouTube, Medium, Reddit
Analysis depth: quick
Number of outlines: 1
```

Claude will analyze trends and provide:
- Search volume and trend data
- User intent breakdown
- Platform-specific insights
- One complete article outline

### Frontend Dev Guidelines

Try creating a component:

```
@frontend-dev-guidelines

I need to create a user profile component that:
- Fetches user data with useSuspenseQuery
- Displays user info with MUI v7 components
- Shows proper loading states with Suspense
- Uses TypeScript for type safety
```

Claude will provide:
- Complete component code with best practices
- Proper TypeScript types
- Data fetching pattern with TanStack Query
- MUI v7 styling examples
- Loading state management

### Blog Post Writer

Transform scattered notes into a polished blog post:

```
@blog-post-writer

Brain dump on my first week using Claude Code:
- Terminal-based, feels natural
- Integrates with my vim workflow
- Can see it using same commands I use (rg, git, npm)
- Made me rethink how I work with AI
- Not replacing me, augmenting my workflow
- Freed up time for creative problem-solving
```

Claude will provide:
- Polished blog post with engaging opening
- Conversational and authentic tone
- Journey-based narrative structure
- Technical details with specific examples
- Vulnerability and honesty
- Forward-looking conclusion

### Writing Effective Prompts

Improve your prompt quality:

```
@writing-effective-prompts

I need to create a prompt for generating a React dashboard. How should I structure it?

Current prompt: "Create a dashboard"

I want it to include:
- Real-time data visualization
- Interactive filtering
- Responsive design
- Export functionality
```

Claude will provide:
- Role-based prompt structure
- Explicit instruction formatting
- Context and motivation guidance
- XML tag examples for complex outputs
- Quality enhancement techniques
- Before/after prompt comparison

## Step 4: Try a More Advanced Request

```
@content-trend-researcher

Topic: "AI tools for content creators in 2025"
Platforms: ALL
Intent focus: commercial
Target audience: Content creators and marketers looking to improve productivity
Content type: blog
Analysis depth: deep
Number of outlines: 3
Include: SEO keyword analysis, content gap identification, platform-specific strategies
```

This will give you comprehensive research with multiple article outlines!

## Common Use Cases

### Content Trend Researcher

#### 1. Blog Post Research
```
@content-trend-researcher

Topic: [Your blog topic]
Platforms: Google Trends, Medium, Blogs, Reddit
Content type: blog
Number of outlines: 2
```

#### 2. Newsletter Topic Ideas
```
@content-trend-researcher

Topic: [Your niche]
Platforms: Substack, Medium, LinkedIn
Content type: newsletter
Number of outlines: 4
```

#### 3. Video Content Planning
```
@content-trend-researcher

Topic: [Your video topic]
Platforms: YouTube, TikTok, X
Content type: video script
Number of outlines: 1
```

#### 4. Competitive Analysis
```
@content-trend-researcher

Topic: [Your topic]
Platforms: ALL
Analysis depth: deep
Focus: Identify content gaps and underserved audiences
```

### Frontend Dev Guidelines

#### 1. Create a New Component
```
@frontend-dev-guidelines

Create a product card component with:
- Props for product data
- Lazy loading for images
- Add to cart button with loading state
- Responsive design with MUI Grid
```

#### 2. Set Up Data Fetching
```
@frontend-dev-guidelines

Show me the pattern for fetching a list of posts using useSuspenseQuery.
Include TypeScript types, error handling, and cache management.
```

#### 3. Build a Complete Feature
```
@frontend-dev-guidelines

I'm building a "posts" feature. Set up the directory structure with:
- API service layer (postsApi.ts)
- Components (PostList, PostDetail)
- Custom hooks (useSuspensePost)
- TypeScript types
- Routing with lazy loading
```

#### 4. Performance Optimization
```
@frontend-dev-guidelines

My DataGrid component is slow. Show me how to:
- Lazy load the DataGrid component
- Memoize expensive computations
- Optimize re-renders with useCallback
- Add proper Suspense boundaries
```

#### 5. Styling with MUI v7
```
@frontend-dev-guidelines

Create a responsive dashboard layout using MUI v7:
- Grid layout with sidebar and main content
- Theme-aware styling
- Mobile-first responsive design
```

### Blog Post Writer

#### 1. Technical Blog Post
```
@blog-post-writer

Random thoughts on switching from Redux to TanStack Query:
- Redux felt like overkill for my use case
- Too much boilerplate for simple data fetching
- TanStack Query's cache-first approach clicked
- useSuspenseQuery eliminated loading state complexity
- Performance improved noticeably
- Reduced code by ~40%
- Wish I'd switched sooner
```

#### 2. Career Story
```
@blog-post-writer

Notes from my first month at new company:
- Impostor syndrome was real
- Everyone seemed 10x more experienced
- First PR was terrifying
- Code review feedback was actually helpful, not harsh
- Pair programming broke the ice
- Made real contributions by week 3
- Learned: everyone feels uncertain sometimes
```

#### 3. Tool Comparison
```
@blog-post-writer

Cursor vs Claude Code - my experience:
- Tried Cursor first, everyone talking about it
- IDE integration felt disconnected from my workflow
- I live in the terminal (vim, tmux, zsh)
- Claude Code met me where I already work
- Not about which is "better"
- About what fits YOUR workflow
- Use tools that match your environment
```

#### 4. Learning Journey
```
@blog-post-writer

My path from "AI will replace us" to daily AI user:
- Started as skeptic
- Thought AI would eliminate developer jobs
- Tried Claude Code out of curiosity
- Watched it use rg, npm test - same tools I use
- Realized it's augmentation, not replacement
- Now use it for repetitive tasks
- Frees me up for creative problem-solving
- Changed perspective entirely
```

### Writing Effective Prompts

#### 1. Improve Code Generation Prompt
```
@writing-effective-prompts

Current: "Create a login form"

Improve this prompt to generate production-quality code with:
- Form validation
- Error handling
- Loading states
- Accessibility features
```

#### 2. Structure Multi-Step Task
```
@writing-effective-prompts

Help me structure a prompt for a data migration task:
- Extract data from MongoDB
- Transform to new schema
- Validate data integrity
- Load into PostgreSQL
- Generate migration report

How should I organize this for Claude Code?
```

#### 3. Long Document Analysis
```
@writing-effective-prompts

I need to analyze 3 large documents (50K tokens total):
- Annual report
- Competitor analysis
- Market research

What's the optimal prompt structure for best results?
```

#### 4. Tool Workflow Optimization
```
@writing-effective-prompts

Optimize this prompt for parallel execution:
- Search 10 files for deprecated API calls
- Update to new API syntax
- Run tests on each file
- Generate change summary
```

#### 5. Role-Based Prompting
```
@writing-effective-prompts

I need help creating a prompt for security code review.
What role should I assign and what structure works best?
```

## Tips for Success

### For Content Trend Researcher:
1. **Be Specific:** "AI automation for small businesses" beats "AI"
2. **Know Your Audience:** Define who you're writing for
3. **Choose Relevant Platforms:** Only select where your audience is active
4. **Start Simple:** Begin with "quick" analysis, then try "deep"

### For Frontend Dev Guidelines:
1. **Be Clear:** Describe exactly what component or feature you need
2. **Mention Tech Stack:** Specify React, TypeScript, MUI v7, TanStack Query
3. **Reference Resources:** Ask for specific guide files when needed
4. **Follow Patterns:** Stick to the established patterns (Suspense, lazy loading, etc.)

### For Blog Post Writer:
1. **Provide Brain Dumps:** Don't organize - just dump scattered thoughts
2. **Include Details:** Add code examples, tool names, specific numbers
3. **Show Journey:** Include before/after, problem/solution elements
4. **Be Honest:** Include uncertainties, mistakes, learning moments
5. **Add Context:** Who you were writing for, what you learned

### For Writing Effective Prompts:
1. **Start with Role:** Always define Claude's expertise first
2. **Be Explicit:** List specific requirements, not vague goals
3. **Use Positive Framing:** Say what TO do, not what NOT to do
4. **Add Context:** Explain WHY requirements matter
5. **Provide Examples:** Use `<example>` tags to show desired output
6. **Long Context:** Put large documents FIRST, queries LAST
7. **XML Tags:** Structure complex responses with tags

## What You'll Get

### Content Trend Researcher:
- **Trend Analysis:** What's hot and what's not
- **Intent Analysis:** What users really want
- **Content Gaps:** Opportunities competitors missed
- **Article Outlines:** Complete structure with SEO optimization
- **Platform Strategy:** Where and when to publish
- **Recommendations:** Next steps and follow-up topics

### Frontend Dev Guidelines:
- **Code Examples:** Production-ready React/TypeScript code
- **Best Practices:** Industry-standard patterns
- **Type Safety:** Complete TypeScript types and interfaces
- **Performance Tips:** Optimization techniques
- **File Structure:** Proper organization patterns
- **Resource Links:** References to detailed guide files

### Blog Post Writer:
- **Polished Blog Post:** Well-structured narrative with engaging flow
- **Conversational Tone:** Authentic voice that sounds like talking to a peer
- **Journey Structure:** Problem → experience → resolution framework
- **Technical Details:** Specific examples with code, commands, tool names
- **Vulnerability:** Honest about uncertainty and learning
- **Varied Rhythm:** Mix of short punchy sentences and longer explanations
- **Strong Ending:** Forward-looking conclusion with momentum

### Writing Effective Prompts:
- **Structured Prompt:** Clear role definition and explicit instructions
- **Quality Enhancement:** Techniques to improve output quality
- **Context Optimization:** Better organization for long documents
- **XML Tag Examples:** Templates for complex multi-section responses
- **Example Formatting:** `<example>` tag patterns
- **Tool Workflow:** Parallelization strategies
- **Before/After Comparison:** See the improvement in prompt quality
- **Best Practices:** Checklist of dos and don'ts

## Next Steps

### For Content Trend Researcher:
1. Read the full documentation: `.claude/skills/content-trend-researcher/SKILL.md`
2. Review usage examples: `.claude/skills/content-trend-researcher/HOW_TO_USE.md`
3. Check sample files: `sample_input.json` and `expected_output.json`
4. Experiment with different topics and platforms
5. Create your own invocation templates for recurring needs

### For Frontend Dev Guidelines:
1. Read the main guide: `.claude/skills/frontend-dev-guidelines/SKILL.md`
2. Explore resource guides in `.claude/skills/frontend-dev-guidelines/resources/`:
   - `component-patterns.md` - Component structure
   - `data-fetching.md` - TanStack Query patterns
   - `file-organization.md` - Project structure
   - `styling-guide.md` - MUI v7 styling
   - `performance.md` - Optimization techniques
   - `complete-examples.md` - Full working examples
3. Apply patterns to your React projects
4. Reference the quick checklist in SKILL.md

### For Blog Post Writer:
1. Read the main guide: `.claude/skills/blog-post-writer/SKILL.md`
2. Study the voice guide: `.claude/skills/blog-post-writer/references/voice-tone.md`
3. Understand narrative structure: `.claude/skills/blog-post-writer/references/story-circle.md`
4. Practice with brain dumps throughout the week
5. Build a library of your polished posts

### For Writing Effective Prompts:
1. Read the guide: `.claude/skills/writing-effective-prompts/SKILL.md`
2. Study the core principles (role, explicit instructions, context)
3. Learn XML tag structuring for complex outputs
4. Practice with simple prompts first, then complex tasks
5. Apply to your daily Claude interactions
6. Build a library of effective prompt templates

## Skill Combinations

Use skills together for powerful workflows:

```
# Content Creation Pipeline
1. Research topic with @content-trend-researcher
2. Collect notes and experiences throughout the week
3. Transform notes into blog post with @blog-post-writer
4. Build accompanying demo with @frontend-dev-guidelines
5. Publish and promote

# Full-Stack Developer Blog
1. Research trending React patterns with @content-trend-researcher
2. Build demo app with @frontend-dev-guidelines
3. Write tutorial blog post with @blog-post-writer
4. Share on developer communities

# Tech Thought Leadership
1. Document learnings and experiences daily
2. Use @blog-post-writer weekly to create posts
3. Build code examples with @frontend-dev-guidelines
4. Research related topics with @content-trend-researcher
5. Maintain consistent publishing schedule

# Prompt Engineering Master Flow
1. Start with basic prompt idea
2. @writing-effective-prompts - Structure and enhance prompt
3. Add role definition and context
4. Test prompt with actual task
5. Iterate based on results
6. Save to personal prompt library

# Quality Code Generation
1. @writing-effective-prompts - Create structured prompt
2. Define role: "You are an expert [language] developer"
3. Add explicit requirements with XML tags
4. @frontend-dev-guidelines - Generate code
5. Review and refine output
6. Apply learnings to next prompt
```

## Need Help?

- Review README.md for complete documentation
- Check the Troubleshooting section in README.md
- Review each skill's documentation files
- Browse the resources/ directory for detailed guides

Happy creating!
