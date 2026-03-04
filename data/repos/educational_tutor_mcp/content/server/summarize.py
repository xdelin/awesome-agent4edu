"""
Defines the 'summarize_text' function which generates concise summaries of text.
"""

from openai import OpenAI
from config.settings import OPENAI_API_KEY, MODEL_NAME

client = OpenAI(api_key = OPENAI_API_KEY)

def summarize_text(text:str, compression_ratio:float = 0.3) -> str:
    """
    Generates a summary of *text* compressed to roughly compression_ratio length.
    compression_ratio must be between 0.1 and 0.8.
    """
    if not text.strip():
        return "Error: text cannot be blank"
    
    # ensure the ratio is within a valid bounds
    ratio = max(0.1, min(compression_ratio, 0.8))
    system_prompt = (
        "You are a world-class summarizer. Reduce the following text to about "
        f"{int(ratio * 100)}% of its original length while preserving key ideas."
    )

    stream = client.chat.completions.create(
        model = MODEL_NAME,
        messages = [
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": text},
        ],
        stream=True,
        temperature=0.5,
    )

    summary = ""
    for chunk in stream:
        delta = getattr(chunk.choices[0].delta, "content", None)
        if delta:
            summary += delta
    
    return summary
