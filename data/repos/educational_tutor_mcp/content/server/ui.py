"""
Builds the Gradio interface for all tools (Explain, Summarize, Flashcards, Quiz).
"""

import gradio as gr
from server.explain import explain_concept
from server.summarize import summarize_text
from server.flashcards import generate_flashcards
from server.quiz import quiz_me

def build_demo():
    """
    Builds and returns a Gradio Blocks UI with tabs for each tool.
    """

    with gr.Blocks() as demo:
        gr.Markdown("# AI TutorMCP")

        # Tab 1: Explain Concept
        with gr.Tab("Explain Concept"):
            q = gr.Textbox(label="Concept / Question", placeholder="e.g. What is quantum computing?")
            lvl = gr.Slider(1, 5, value=3, step=1, label="Explanation Level (1 = simple, 5 = advanced)")
            out1 = gr.Markdown()
            gr.Button("Explain").click(explain_concept, inputs=[q, lvl], outputs=out1, api_name="explain_concept")
        
        # Tab 2: Summarize Text
        with gr.Tab("Summarize Text"):
            txt = gr.Textbox(lines=8, label="Long Text", placeholder="Paste your text here...")
            ratio = gr.Slider(0.1, 0.8, value=0.3, step=0.05, label="Compression Ratio")
            out2 = gr.Markdown()
            gr.Button("Summarize").click(summarize_text, inputs=[txt, ratio], outputs=out2, api_name="summarize_text")

        # Tab 3: Flashcards
        with gr.Tab("Generate Flashcards"):
            topic_fc = gr.Textbox(label="Topic", placeholder="e.g. Machine Learning Basics")
            n_fc = gr.Slider(1, 20, value=5, step=1, label="# of Flashcards")
            out3 = gr.Markdown()
            gr.Button("Generate").click(generate_flashcards, inputs=[topic_fc, n_fc], outputs=out3, api_name="generate_flashcards")
        
        # Tab 4: Quiz
        with gr.Tab("Quiz Me"):
            topic_q = gr.Textbox(label="Topic", placeholder="e.g. Neural Networks")
            lvl_q = gr.Slider(1, 5, value=3, step=1, label="Difficulty Level")
            n_q = gr.Slider(1, 15, value=5, step=1, label="# of Questions")
            out4 = gr.Markdown()
            gr.Button("Start Quiz").click(quiz_me, inputs=[topic_q, lvl_q, n_q], outputs=out4, api_name="quiz_me")
    
    return demo