# Claude Project Guide for Kvante

Kvante is a Danish-developed AI-powered math assistant that helps students solve math problems step by step — without giving away the final answer. The project emphasizes solving tasks by hand: students write their solutions on paper, take a photo of their work, and receive stepwise feedback based on their progress over time.

## Project Goals
- Enable students to scan handwritten math problems using an intuitive UI
- Use Claude/GPT to generate step-by-step guidance — never just the final answer
- Collect user feedback for evaluation and improvement
- Build a test dashboard for evaluating prompt strategies

## Key Technologies
- Frontend: React, Tailwind, Vite
- Backend: Node.js + Express
- Claude via API (Anthropic), OpenAI API (optional fallback)
- OCR via Tesseract.js

## Code Style
- Use ES modules (import/export)
- Use async/await for asynchronous code
- Format with Prettier, lint with ESLint
- Keep components modular and reusable
- Use `.ts` and `.tsx` where possible

## Prompts and LLM usage
- Prompts are stored in `/prompts/`
- Always include system-level instructions to avoid hallucination
- Do not generate answers to math problems — only step-by-step reasoning

## Testing & Evaluation
- Use `testing/` directory for prompt experiments
- Write logs to local JSON for now
- Use `apps/dashboard` for visualizing prompt results and feedback

## Claude Code Behavior
- Claude may generate, refactor, or document code when prompted
- Claude should follow the structure defined above
- Claude should NOT refactor code that is explicitly marked as final or legacy
