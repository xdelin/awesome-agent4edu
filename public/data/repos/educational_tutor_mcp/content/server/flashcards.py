"""
Defines the 'generate_flashcards' function which streams study flashcards in JSON format.
"""

from openai import OpenAI
from typing import Generator
from config.settings import OPENAI_API_KEY, MODEL_NAME

client = OpenAI(api_key = OPENAI_API_KEY)

def generate_flashcards(topic: str, num_cards: int = 5) -> Generator[str, None, None]:
    """
    Stream *num_cards* Q/A flashcards for *topic* in JSON lines format.
    Example: {"q": "What is...", "a": "It is..."}
    """
    if not topic.strip():
        yield "Error: topic canot be blank."
        return
    if not (1 <= num_cards <= 20):
        yield "Error: num_cards must be between 1 and 20."
        return
    
    system_prompt = (
        "You are an AI that generates study flashcards. "
        'Return each flashcard on its own line as JSON: {"q": <question>, "a": <answer>}'
    )
    user_prompt = f"Create {num_cards} flashcards about {topic}."

    stream = client.chat.completions.create(
        model=MODEL_NAME,
        messages=[
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": user_prompt},
        ],
        stream=True,
        temperature=0.8,
    )

    output = ""
    for chunk in stream:
        delta = getattr(chunk.choices[0].delta, "content", None)
        if delta:
            output += delta
    yield output

