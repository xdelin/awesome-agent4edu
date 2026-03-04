"""
Main entry point for running the AI Tutor MCP server.
Exposes all learning tools (explain, summarize, flashcards, quiz) via Gradio MCP.
"""
import gradio as gr
from fastapi import FastAPI, Request
from fastapi.middleware.cors import CORSMiddleware
import uvicorn

from server.ui import build_demo
from config.settings import MCP_PORT

from server.explain import explain_concept
from server.summarize import summarize_text
from server.flashcards import generate_flashcards
from server.quiz import quiz_me

app = FastAPI()

# allow client to connect safely
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

@app.post("/api/explain_concept")
async def api_explain(request: Request):
    data = await request.json()
    question = data.get("question", "")
    level = int(data.get("level", 3))
    result = explain_concept(question, level)
    return {"output": result}

@app.post("/api/summarize_text")
async def api_summarize(request: Request):
    data = await request.json()
    text = data.get("text", "")
    ratio = float(data.get("ratio", 0.3))
    result = summarize_text(text, ratio)
    return {"output": result}

@app.post("/api/generate_flashcards")
async def api_flashcards(request: Request):
    data = await request.json()
    topic = data.get("topic", "")
    n = int(data.get("n", 5))
    result = generate_flashcards(topic, n)
    return {"output": result}

@app.post("/api/quiz_me")
async def api_quiz(request: Request):
    data = await request.json()
    topic = data.get("topic", "")
    level = int(data.get("level", 3))
    n = int(data.get("n", 5))
    result = quiz_me(topic, level, n)
    return {"output": result}

# Mount Gradio interface at root
gradio_app = build_demo()
app = gr.mount_gradio_app(app, gradio_app, path="/")

if __name__ == "__main__":
    print(f"Starting AI Tutor REST + Gradio server on port {MCP_PORT}")
    uvicorn.run(app, host="0.0.0.0", port=MCP_PORT)