import sys
import os
import json
from load_reports import load_all_reports
from process_student_data import process_reports_to_summary
from prompt_templates import generate_feedback_prompt

from google import genai

# Set your Gemini API key here or ensure it's in your environment variables securely
GEMINI_API_KEY = "API KEY"

client = genai.Client(api_key=GEMINI_API_KEY)

def call_gemini_api(prompt: str) -> str:
    response = client.models.generate_content(
        model="gemini-2.0-flash",
        contents=prompt
    )
    return response.text.strip()

def load_single_report(filepath: str) -> dict:
    """Load a single JSON report from the given filepath, 
    unwrap if it's a list with one dict."""
    with open(filepath, "r", encoding="utf-8") as f:
        data = json.load(f)
    if isinstance(data, list) and len(data) == 1 and isinstance(data[0], dict):
        data = data[0]
    return data

def main(data_folder="data", single_file=None):
    if single_file:
        # Load only one report from the given file path
        filepath = os.path.join(data_folder, single_file)
        report = load_single_report(filepath)
        reports = [report]
    else:
        # Load all reports from folder
        reports = load_all_reports(data_folder)

    # Process reports to structured summary
    summary = process_reports_to_summary(reports)

    # Generate prompt text for LLM
    prompt = generate_feedback_prompt(summary)

    print("\n=== Prompt sent to LLM ===\n")
    print(prompt)

    # Call Gemini API
    print("\nCalling Gemini API...\n")
    llm_response = call_gemini_api(prompt)

    print("\n=== Gemini Response ===\n")
    print(llm_response)

    return summary, llm_response

if __name__ == "__main__":
    folder = sys.argv[1] if len(sys.argv) > 1 else "data"
    file_arg = sys.argv[2] if len(sys.argv) > 2 else None
    main(folder, single_file=file_arg)
