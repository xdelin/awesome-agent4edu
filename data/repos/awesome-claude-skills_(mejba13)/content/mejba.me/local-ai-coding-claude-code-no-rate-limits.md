---
title: "How to Run Claude Code Locally Without Rate Limits: Build AI Apps for Free"
slug: local-ai-coding-claude-code-no-rate-limits
tags: AI Coding, Local AI, Claude Code, LM Studio, Self-Hosted AI
---

# How to Run Claude Code Locally Without Rate Limits: Build AI Apps for Free

You're in the zone, building a complex feature with your AI coding assistant. The code is flowing, the architecture is coming together, and then—**"Rate limit exceeded. Please try again in 60 seconds."**

Your momentum crashes. Your API credits are burning through your budget. And you're stuck waiting for permission to continue coding on your own machine.

I experienced this frustration firsthand, which led me to discover something game-changing: **running Claude Code entirely on local hardware with unlimited sessions, zero API costs, and complete privacy.** Even better, I built a fully functional AI-powered PDF chat application where the same local AI model that wrote the code also powers the app itself.

In this guide, I'll show you exactly how I set up a local AI coding environment using Claude Code Router, LM Studio, and local language models, then used it to build a production-ready Next.js application. You'll learn the technical setup, real challenges I faced, performance comparisons, and when local AI makes sense versus cloud alternatives.

## Why Local AI Coding Matters: Beyond Just Saving Money

### The Hidden Costs of Cloud AI Development

API rate limits aren't just annoying—they're expensive and disruptive. Here's what developers face with cloud-based AI coding tools:

- **Rate throttling during peak productivity**: Hit your limit mid-session when you're most productive
- **Escalating costs**: $20-100+/month for API access, scaling with usage
- **Privacy concerns**: Your proprietary code and documents sent to external servers
- **Internet dependency**: No coding assistance without stable connectivity
- **Token anxiety**: Constantly monitoring context limits and optimizing prompts to save money

### The Local AI Advantage

Running AI models locally on your own hardware eliminates these barriers:

**Unlimited Sessions**: Code for hours without interruption or throttling. No rate limits, no waiting periods, no usage caps.

**Zero Recurring Costs**: After initial setup, no monthly API fees. Your GPU becomes your AI infrastructure.

**Complete Privacy**: Code, documents, and proprietary information never leave your machine. Perfect for sensitive projects or regulated industries.

**Customization Freedom**: Choose models optimized for your needs, adjust context lengths, experiment without financial penalties.

**Learning Opportunity**: Understand how AI coding agents work under the hood, gaining insights into prompt engineering and model capabilities.

### When to Choose Local vs Cloud AI

**Local AI is ideal for:**
- Long coding sessions and exploratory development
- Sensitive codebases requiring privacy
- Learning AI coding workflows without budget constraints
- Teams wanting to control AI infrastructure
- Projects where internet connectivity is unreliable

**Cloud AI excels at:**
- Complex architectural decisions requiring cutting-edge models
- Rapid prototyping when setup time matters
- Occasional use where hardware investment isn't justified
- Collaborative features and cloud integrations

## Technical Setup: Building Your Local AI Coding Environment

### Prerequisites and Hardware Requirements

Before diving in, here's what you'll need:

**Hardware:**
- Modern GPU with 8GB+ VRAM (NVIDIA preferred)
  - Quen 3 (7B parameters): ~8GB VRAM, 50,000 token context
  - Quen 7B: ~16GB VRAM, 250,000 token context
- 16GB+ system RAM
- 50GB+ free storage for models

**Software:**
- Windows 10/11 (or Linux/macOS)
- Windows Subsystem for Linux (WSL) for Windows users
- LM Studio
- Claude Code Router
- Node.js 18+ (for the PDF app project)

### Step 1: Setting Up Windows Subsystem for Linux (WSL)

Why WSL? **AI coding agents are optimized for bash commands, not Windows PowerShell.** Running Claude Code in a Linux environment provides better compatibility and fewer errors.

**Install WSL on Windows:**

```bash
# Open PowerShell as Administrator
wsl --install

# Install Ubuntu distribution
wsl --install -d Ubuntu

# Launch Ubuntu and create your user account
wsl
```

**Configure WSL for AI development:**

```bash
# Update packages
sudo apt update && sudo apt upgrade -y

# Install essential development tools
sudo apt install build-essential git curl -y

# Install Node.js (for our PDF app)
curl -fsSL https://deb.nodesource.com/setup_18.x | sudo -E bash -
sudo apt install -y nodejs
```

**Why this matters:** When the AI agent runs commands like `mkdir`, `npm install`, or `git init`, it executes them as bash commands. WSL provides the native Linux environment the AI expects, reducing errors and improving compatibility.

### Step 2: Installing and Configuring LM Studio

LM Studio is your local model server, providing an OpenAI-compatible API for running models on your hardware.

**Installation:**

1. Download LM Studio from [lmstudio.ai](https://lmstudio.ai)
2. Install and launch the application
3. Navigate to the model search/download section

**Download AI Models:**

Search for and download these models:

- **Quen 3 (7B parameters)**: Great starting point, ~50k token context, runs on 8GB VRAM
- **Quen 7B**: Larger context (250k tokens), better for document processing, requires 16GB VRAM

**Configure the Local Server:**

```
1. Click "Local Server" tab in LM Studio
2. Select your downloaded model (start with Quen 3)
3. Configure settings:
   - Context Length: 50000 (Quen 3) or 250000 (Quen 7B)
   - GPU Layers: Max (for best performance)
   - Temperature: 0.7 (balanced creativity)
4. Start the server (default port: 1234)
5. Note the endpoint: http://localhost:1234/v1
```

**Verify the server:**

```bash
# Test the local API
curl http://localhost:1234/v1/models

# Should return your loaded model details
```

### Step 3: Routing Claude Code to Local Models

Now connect Claude Code to your local LM Studio server instead of cloud APIs.

**Install Claude Code Router:**

```bash
# In your WSL Ubuntu environment
npm install -g @anthropic/claude-code-router
```

**Configure the router:**

Create a configuration file pointing to your local LM Studio:

```bash
# Create config directory
mkdir -p ~/.claude-code

# Create router config
cat > ~/.claude-code/router-config.json << 'EOF'
{
  "provider": "openai-compatible",
  "baseURL": "http://localhost:1234/v1",
  "model": "quen-3-7b",
  "apiKey": "not-needed-for-local"
}
EOF
```

**Launch Claude Code with local routing:**

```bash
claude-code --router-config ~/.claude-code/router-config.json
```

**Success!** You now have Claude Code running entirely on local infrastructure, powered by your own GPU.

### Step 4: Understanding Permission Modes

Claude Code has two operating modes for command execution:

**Normal Mode (Recommended for production):**
- Requires approval before executing commands
- Safer for untrusted environments
- Better for learning what the AI is doing

**Unleashed Mode (Use with caution):**
- Skips command permission prompts
- Dramatically faster development cycles
- **Only use in isolated/safe environments**

```bash
# Start in unleashed mode (isolated WSL environment)
claude-code --dangerously-skip-permissions --router-config ~/.claude-code/router-config.json
```

I used unleashed mode in my isolated WSL environment to speed up the PDF app development, but always review generated code afterward.

## Building the PDF Chat Application: Real-World Local AI Coding

Now for the practical demonstration: building an AI-powered PDF chat application where the local model both writes the code and powers the app.

### Initial Prompt Engineering

The key to effective AI coding is providing a clear, detailed initial prompt. Here's what I used:

```
Build a Next.js 13 application that allows users to upload PDF files
and ask questions about the content. Requirements:

1. Tech stack: Next.js 13, React, Node.js
2. PDF processing: Extract text from PDF pages
3. AI integration: Use local LM Studio API for question answering
4. Features:
   - Upload PDF files to a local directory
   - Display PDF page by page
   - Ask questions about the PDF content
   - Inject relevant pages into AI context for answers
   - Show page numbers with responses
5. UI: Clean, minimal interface
6. API: Create Next.js API route to communicate with LM Studio

Start by creating the project structure and specification file.
```

### AI-Generated Project Scaffolding

The AI immediately created a comprehensive specification file and folder structure:

```
my-pdf-chat-app/
├── pages/
│   ├── api/
│   │   └── chat.js          # AI integration endpoint
│   ├── index.js             # Main PDF viewer interface
│   └── _app.js              # Next.js app wrapper
├── components/
│   ├── PDFViewer.js         # PDF display component
│   ├── ChatInterface.js     # Q&A interface
│   └── PageNavigator.js     # Page controls
├── lib/
│   ├── pdfProcessor.js      # PDF text extraction
│   └── aiClient.js          # LM Studio API client
├── public/
│   └── pdfs/                # PDF storage directory
├── package.json
└── README.md
```

The AI also generated a to-do list:

1. ✅ Create Next.js project structure
2. ✅ Install dependencies (pdf-parse, axios)
3. ✅ Build PDF processing utilities
4. ✅ Create API route for AI chat
5. ✅ Build React components
6. ⚠️ Test and debug routing
7. ✅ Implement PDF question answering

### Challenge 1: Next.js 13 Routing Confusion

When I first ran the application, I hit an immediate **404 error**. The AI had generated code using older Next.js routing conventions that weren't compatible with Next.js 13's new app directory structure.

**The Problem:**
```javascript
// AI generated this (old pattern):
pages/api/chat.js

// Next.js 13 expects:
app/api/chat/route.js
```

**The local AI struggled** to understand and fix the routing issue, repeatedly suggesting the same outdated pattern.

**The Solution:**
I switched to a cloud AI model temporarily to fix the routing structure. The cloud AI:
- Quickly identified the Next.js 13 routing issue
- Restructured the file organization correctly
- Updated import paths and API endpoint conventions

**Key Lesson:** Local AI models (especially smaller ones) can struggle with complex framework-specific issues. The hybrid approach—local AI for scaffolding and iteration, cloud AI for complex debugging—proved most effective.

### Challenge 2: Token Limits with Large Documents

With the app running, I tested it with a programming book PDF. When asking about content, I hit this error:

```
Error: Context length exceeded
Current context: 52,847 tokens
Maximum supported: 50,000 tokens
```

The Quen 3 model's 50,000 token limit couldn't handle injecting entire PDFs into context.

**First Solution: Selective Page Injection**

Instead of loading entire PDFs, inject only relevant pages:

```javascript
// lib/pdfProcessor.js
async function getRelevantPages(pdfPath, query, maxPages = 5) {
  const pdfData = await pdfParse(pdfPath);
  const pages = pdfData.text.split('\n\n'); // Simple page split

  // Simple relevance scoring (in production, use embeddings)
  const scoredPages = pages.map((page, index) => ({
    content: page,
    pageNum: index + 1,
    score: calculateRelevance(page, query) // Keyword matching
  }));

  return scoredPages
    .sort((a, b) => b.score - a.score)
    .slice(0, maxPages);
}
```

This worked for targeted questions but limited the AI's ability to answer questions requiring broader context.

**Better Solution: Upgrade to Quen 7B (250k tokens)**

I switched to the larger Quen 7B model with 250,000 token context:

```bash
# In LM Studio:
# 1. Download Quen 7B model
# 2. Load it in the local server
# 3. Update context length to 250000
# 4. Restart the server
```

**Results:**
- ✅ Entire 200-page technical book loaded into context
- ✅ Questions answered with full document awareness
- ⚠️ GPU usage at 100%, VRAM maxed at 16GB
- ⚠️ Response time increased from 2-3s to 8-12s

**Trade-off Analysis:**

| Model | Context Tokens | VRAM | Response Time | Best For |
|-------|---------------|------|---------------|----------|
| Quen 3 | 50,000 | 8GB | 2-3s | Targeted questions, code snippets |
| Quen 7B | 250,000 | 16GB | 8-12s | Full document analysis, books |

### Challenge 3: Page Numbering Discrepancies

The app displayed page numbers based on PDF footer text, which didn't always match the internal PDF page index.

```javascript
// Problem: User asks "What's on page 28?"
// PDF footer says "Page 28"
// But it's actually PDF internal page 30

// Solution: Display both numbers
<PageIndicator>
  Display Page: {displayPage} | PDF Page: {internalPage}
</PageIndicator>
```

This required manual UI tweaks the AI couldn't anticipate.

## Real-World Performance: Local vs Cloud AI

After building the complete application, here's my honest performance comparison:

### Development Speed

**Local AI (Quen 3/7B):**
- Project scaffolding: ✅ Excellent (5 minutes)
- Writing boilerplate code: ✅ Very good
- Understanding framework conventions: ⚠️ Struggles with latest versions
- Complex debugging: ❌ Often gets stuck in loops
- Overall: Great for 70% of coding tasks

**Cloud AI (Claude Sonnet):**
- Project scaffolding: ✅ Excellent (3 minutes)
- Writing boilerplate code: ✅ Very good
- Understanding framework conventions: ✅ Excellent, up-to-date
- Complex debugging: ✅ Excellent, understands context deeply
- Overall: Better for complex problem-solving

### Response Time

| Task | Local (Quen 3) | Local (Quen 7B) | Cloud (Sonnet) |
|------|----------------|-----------------|----------------|
| Generate component | 2-3s | 8-12s | 1-2s |
| Debug routing error | 5-8s | 15-20s | 2-3s |
| Explain code block | 1-2s | 5-8s | 1s |
| Full file refactor | 10-15s | 30-45s | 3-5s |

### Cost Comparison (100 hours of coding)

**Local AI:**
- Hardware: $800 GPU (one-time, if not already owned)
- Electricity: ~$15-30/month @ $0.12/kWh
- Total year 1: $1,000-1,200

**Cloud AI:**
- API costs: $50-150/month (varies by usage)
- Total year 1: $600-1,800

**Break-even:** 6-12 months, then local becomes pure savings

### Privacy & Security

**Local AI:** ✅ Complete privacy, no data leaves your machine
**Cloud AI:** ⚠️ Data sent to external servers, check privacy policies

## Best Practices for Local AI Coding

Through building this PDF chat app locally, I learned these essential practices:

### 1. Effective Prompt Engineering

**Bad prompt:**
```
Build a PDF reader app
```

**Good prompt:**
```
Build a Next.js 13 application with these requirements:
- Tech stack: Next.js 13 app directory, React 18, TypeScript
- PDF processing: Use pdf-parse library to extract text
- AI integration: LM Studio local API at localhost:1234
- Features: [detailed list]
- File structure: [specify organization]
- Error handling: [specify approach]
```

Specificity dramatically improves AI-generated code quality.

### 2. Managing Dependencies and Code Quality

AI-generated package.json files often include:
- ❌ Outdated package versions
- ❌ Unnecessary dependencies
- ❌ Missing peer dependencies

**Always review and update:**

```bash
# Check for outdated packages
npm outdated

# Update to latest compatible versions
npm update

# Audit for security vulnerabilities
npm audit fix
```

### 3. Document Chunking and Vector Embeddings

For production PDF chat applications, loading entire documents into context is impractical.

**Better approach:**

```javascript
// 1. Split document into chunks
const chunks = splitDocument(pdfText, { chunkSize: 1000, overlap: 200 });

// 2. Generate embeddings for each chunk
const embeddings = await generateEmbeddings(chunks);

// 3. Store in vector database
await vectorDB.insert(embeddings);

// 4. Query with user question
const relevantChunks = await vectorDB.search(userQuery, { topK: 5 });

// 5. Inject only relevant chunks into AI context
const answer = await ai.complete({
  context: relevantChunks.join('\n'),
  question: userQuery
});
```

This approach:
- ✅ Scales to documents of any size
- ✅ Reduces token usage dramatically
- ✅ Improves answer relevance
- ✅ Works with smaller context models

Libraries to explore: LangChain, llamaindex, Pinecone, ChromaDB

### 4. Balancing Context Length vs Response Speed

Larger context windows enable richer responses but slow generation:

**Optimization strategies:**
- Use smaller models (Quen 3) for rapid iteration
- Switch to larger models (Quen 7B) when full context needed
- Implement smart context pruning (remove irrelevant sections)
- Cache frequently accessed document sections

### 5. Security Considerations for Unleashed Mode

**Never use unleashed/permission-skipping mode if:**
- Working with production codebases
- Connected to production databases
- On a shared system
- Handling sensitive credentials

**Safe for unleashed mode:**
- Isolated WSL/container environments
- Personal development projects
- Sandboxed testing environments
- When you audit generated commands after execution

## Limitations and When to Switch to Cloud

Local AI coding isn't perfect. Here's when I switched back to cloud models:

### When Cloud AI Is Superior

1. **Complex architectural decisions**: Cloud models better understand system design patterns
2. **Latest framework features**: Cloud models trained on more recent data
3. **Multi-file refactoring**: Better project-wide awareness
4. **Natural language understanding**: More nuanced interpretation of requirements
5. **Debugging complex errors**: Deeper reasoning about edge cases

### Hybrid Workflow (Best of Both Worlds)

My optimal workflow:

```
1. Initial scaffolding → Local AI (fast, unlimited iteration)
2. Complex routing/architecture → Cloud AI (better understanding)
3. Component implementation → Local AI (cost-effective, private)
4. Debugging framework issues → Cloud AI (more up-to-date)
5. Refactoring and optimization → Local AI (iterative, unlimited)
6. Production deployment review → Cloud AI (final quality check)
```

This hybrid approach saved ~60% on API costs while maintaining high code quality.

## The Future of Local AI Development

Based on this hands-on experience, here's where local AI coding is heading:

**Near-term improvements (6-12 months):**
- Smaller models with comparable performance (better efficiency)
- Longer context windows standard (500k+ tokens)
- Better framework/library awareness in open models
- Specialized coding models optimized for local deployment

**Exciting possibilities:**
- Team-hosted AI coding servers (shared local infrastructure)
- Domain-specific fine-tuned models (your codebase's style)
- Real-time code review agents running locally
- Offline-first AI development workflows

**Current limitations closing:**
- Model quality gap narrowing between local and cloud
- Hardware becoming more affordable (GPU prices stabilizing)
- Tooling maturing (LM Studio, Claude Code Router, etc.)

## Conclusion: Taking Control of Your AI Coding Future

Running Claude Code locally with LM Studio and local models isn't just about saving money on API costs—though that's a significant benefit. It's about **taking control of your development workflow**, eliminating external dependencies, protecting your privacy, and coding without artificial limits.

Through building this PDF chat application, I demonstrated that local AI coding is:
- ✅ **Viable for real projects**: Built a production-ready Next.js app entirely locally
- ✅ **Cost-effective**: Zero API costs after initial hardware investment
- ✅ **Private**: Code and documents never leave your machine
- ⚠️ **Not perfect**: Cloud AI still superior for complex reasoning
- ✅ **Rapidly improving**: Model quality advancing quickly

### Your Next Steps

Ready to try local AI coding?

1. **Start small**: Install LM Studio, download Quen 3, test basic prompts
2. **Build something real**: Pick a small project (Todo app, calculator, API wrapper)
3. **Learn the limits**: Push the model until you find where it struggles
4. **Develop hybrid workflow**: Identify which tasks work locally vs cloud
5. **Iterate and optimize**: Adjust models, context lengths, and prompts

The future of AI-assisted development is hybrid—local for privacy, cost, and iteration; cloud for complexity and cutting-edge capabilities. By mastering both, you'll code faster, smarter, and with complete control.

The era of unlimited AI coding is here. Your GPU is waiting.

## 🤝 Hire / Work with me:

* 🔗 **Fiverr** (custom builds, integrations, performance): [fiverr.com/s/EgxYmWD](https://www.fiverr.com/s/EgxYmWD)
* 🌐 **Mejba Personal Portfolio**: [mejba.me](https://www.mejba.me)
* 🏢 **Ramlit Limited**: [ramlit.com](https://www.ramlit.com)
* 🎨 **ColorPark Creative Agency**: [colorpark.io](https://www.colorpark.io)
* 🛡 **xCyberSecurity Global Services**: [xcybersecurity.io](https://www.xcybersecurity.io)