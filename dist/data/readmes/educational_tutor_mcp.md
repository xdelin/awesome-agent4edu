# TutorMCP – AI Tutor Toolkit

TutorMCP is an **AI-powered educational toolkit** built on the **Model Context Protocol (MCP)**.
It exposes four interactive learning tools via a **Gradio interface** and an **MCP-compatible API**, enabling seamless interaction between AI agents and human learners.

## Overview

TutorMCP provides four core AI-driven tools:

| Tool                     | Description                                                      | Arguments                                               |
|------------------------ | ---------------------------------------------------------------- | ------------------------------------------------------- |
|  `explain_concept`     | Streams an explanation of a concept at a chosen depth (1–5).     | `question: str`, `level: int`                           |
|  `summarize_text`      | Summarizes text to a target length ratio.                        | `text: str`, `compression_ratio: float (0.1–0.8)`       |
|  `generate_flashcards` | Generates JSON-lines flashcards for any topic.                   | `topic: str`, `num_cards: int (1–20)`                   |
|  `quiz_me`             | Creates a quiz with multiple-choice questions and an answer key. | `topic: str`, `level: int`, `num_questions: int (1–15)` |

Each function **streams partial responses in real time** for a smooth, low-latency learning experience.

## Project Structure

```
tutormcp/
│
├── main.py                  # Entry point (run as client or server)
│
├── requirements.txt
│
├── config/
│   └── settings.py          # Environment variables and constants
│
├── server/
│   ├── tutor_mcp_server.py  # Main MCP Gradio server
│   ├── explain.py           # Explain tool logic
│   ├── summarize.py         # Summarization logic
│   ├── flashcards.py        # Flashcard generation
│   ├── quiz.py              # Quiz generation
│   └── ui.py                # Gradio app UI definition
│
├── client/
│   ├── agent.py             # Smart MCP AI agent
│   ├── mcp_client.py        # Connects to MCP server
│   ├── runner.py            # Interactive chat loop
│   └── utils.py             # Helpers like print_markdown()
│
└── tests/
    ├── test_server.py
    ├── test_client.py
    └── test_tools.py
```

##  Installation

### 1. Clone the repository

```bash
git clone https://github.com/iamisaackn/tutormcp.git
cd tutormcp
```

### 2. Set up a virtual environment

```bash
python -m venv venv
source venv/bin/activate   # (Linux/Mac)
venv\Scripts\activate      # (Windows)
```

### 3. Install dependencies

```bash
pip install -r requirements.txt
```

### 4. Create a `.env` file

```bash
OPENAI_API_KEY=your_openai_api_key_here
```

---

## Usage

### Start the MCP server (Gradio interface)

```bash
python main.py --mode server
```

This launches a local Gradio interface at `http://localhost:7860`
and exposes the MCP endpoints for your AI tools.

### Start the client agent

```bash
python main.py --mode client
```

This starts an **interactive chat session** where your MCP Agent connects to the server and routes user queries to the correct tool.

---

### Example Commands

**Explain a concept:**

```
User: Explain quantum tunneling like I’m 10.
```

**Summarize text:**

```
User: Summarize this article to 20%.
```

**Generate flashcards:**

```
User: Create 10 flashcards on Python data structures.
```

**Take a quiz:**

```
User: Quiz me on neural networks, level 4, 5 questions.
```

---

## Tech Stack

* **Python 3.10+**
* **Gradio** – User interface and MCP server hosting
* **OpenAI API** – Core LLM for responses
* **httpx / requests** – For async API calls
* **dotenv** – Secure environment management

---

## Key Design Choices

* **MCP Integration** – Enables multi-tool orchestration and real-time communication via Server-Sent Events (SSE).
* **Streaming Responses** – Partial tokens are streamed for faster user feedback.
* **Extensible Architecture** – Easily add new tools (e.g., code explainer, diagram generator) under `server/`.

---

## Development Notes

Run code formatters and lint checks:

```bash
black .
flake8
```

Run tests:

```bash
pytest tests/
```

---

## License

MIT License © 2025 Isaac K. Ngugi

---

## Author

**Isaac K. Ngugi**
AI Engineer & Data Scientist – Nairobi, Kenya
📧 [[your.email@example.com](mailto:your.email@example.com)]
💼 [LinkedIn/GitHub link here]
