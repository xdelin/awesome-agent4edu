"""
Defines the 'explain_concept' function which streams explanations at different complexity levels.
"""

from openai import OpenAI
from typing import Generator
from config.settings import OPENAI_API_KEY, MODEL_NAME

# initialize OpenAI client
client = OpenAI(api_key = OPENAI_API_KEY)

# map explanation levels to audience understanding
EXPLANATION_LEVELS = {
    1: "like I'm 5 years old",
    2: "like I'm 10 years old",
    3: "like a high school student",
    4: "like a college student",
    5: "like an expert in the field"
}

def explain_concept(question:str, level:int) -> str:
    """
    explain the *question* at the requested *level* (1–5).
    Level 1 = simple, Level 5 = advanced.
    """

    # guard against empty input
    if not question.strip():
        return "Error: question cannot be blank"
    
    # pick explanation level
    level_desc = EXPLANATION_LEVELS.get(level, "clearly and concisely")
    system_prompt = f"You are a helpful AI Tutor. Explain the following concept {level_desc}."

    try:
        response = client.chat.completions.create(
            model=MODEL_NAME,
            messages=[
                {"role": "system", "content": system_prompt},
                {"role": "user", "content": question},
            ],
            temperature=0.7,
            timeout = 60,
        )
        return response.choices[0].message.content.strip()
    
    except Exception as e:
        return f"Server error: {e}"