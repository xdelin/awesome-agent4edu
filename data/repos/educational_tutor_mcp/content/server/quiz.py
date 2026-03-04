"""
Defines the 'quiz_me' function which creates interactive multiple-choice quizzes.
"""

from openai import OpenAI
from typing import Generator
from config.settings import OPENAI_API_KEY, MODEL_NAME, EXPLANATION_LEVELS

client = OpenAI(api_key=OPENAI_API_KEY)

def quiz_me(topic: str, level: int = 3, num_questions: int = 5) -> Generator[str, None, None]:
    """
    Stream a quiz with multiple-choice questions, followed by an ANSWER KEY.
    """

    if not topic.strip():
        yield "ErrorL topic cannot be blank."
        return
    if not (1 <= num_questions <= 15):
        yield "Error: num_questions must be between 1 and 15."
        return
    
    level_desc = EXPLANATION_LEVELS.get(level, "at an intermediate level")
    system_prompt = (
        f"You are an AI quiz master. Create a {num_questions}-question multiple-choice quiz "
        f"about '{topic}', explained {level_desc}. Each question should have 4 options (A–D). "
        "After listing all questions, provide a section titled '### ANSWER KEY' showing the correct answers clearly."
    )

    stream = client.chat.completions.create(
        model=MODEL_NAME,
        messages=[
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": f"Create a quiz about: {topic}"},
        ],
        stream=True,
        temperature=0.7,
    )

    partial = ""
    for chunk in stream:
        delta = getattr(chunk.choices[0].delta, "content", None)
        if delta:
            partial += delta
            yield partial