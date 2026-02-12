# Kvante

Kvante is a Danish-developed AI-powered math assistant that helps students solve math problems step by step — without giving away the final answer. Students write their solutions on paper, take photos of their work, and receive stepwise feedback based on their progress over time.

## Project Structure

```
kvante/
├── backend/          # Node.js Express server
├── frontend/         # React Vite application
├── dashboard/        # Prompt testing dashboard
├── ocr/             # OCR utilities using Tesseract.js
├── prompts/         # Claude/OpenAI prompt collection
├── tests/           # Automated tests
├── data/            # Local JSON data storage
└── CLAUDE.md        # Project instructions for Claude
```

## Getting Started

1. Install dependencies:
   ```bash
   npm install
   ```

2. Start development servers:
   ```bash
   npm run dev
   ```

3. Run tests:
   ```bash
   npm test
   ```

## Key Features

- Handwritten math problem scanning via intuitive UI
- Step-by-step guidance using Claude/GPT APIs
- User feedback collection for evaluation and improvement
- Test dashboard for evaluating prompt strategies

## Technologies

- **Frontend**: React, Tailwind CSS, Vite
- **Backend**: Node.js, Express, TypeScript
- **AI**: Claude API (Anthropic), OpenAI API
- **OCR**: Tesseract.js
