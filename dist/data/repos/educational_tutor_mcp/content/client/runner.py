"""
- Lets you talk to your tutor in the terminal.
- Uses menu-based interaction to select which tool to run.
- for testing before connecting a full UI.
"""

from client.agent import AITutorAgent
from client.utils import print_markdown

def main():
    """
    Simple interactive loop for testing the MCP tutor client.
    Run: python -m client.runner
    """

    print_markdown("Welcome to the AI Tutor MCP Client!")
    agent = AITutorAgent()

    while True:
        user_input = input("\nAsk me something (or type 'quit'): ").strip()
        if user_input.lower() in ["quit", "exit"]:
            print_markdown("Goodbye, future genius!")
            break

        print("\nChoose mode:")
        print("1. Explain concept")
        print("2. Summarize text")
        print("3. Generate flashcards")
        print("4. Quiz me")
        choice = input("> ").strip()

        if choice == "1":
            level = int(input("Explanation level (1–5): ") or "3")
            response = agent.explain(user_input, level)
        elif choice == "2":
            response = agent.summarize(user_input)
        elif choice == "3":
            response = agent.flashcards(user_input)
        elif choice == "4":
            response = agent.quiz(user_input)
        else:
            response = "Invalid choice."

        print_markdown(f"\n**Response:**\n{response}")

if __name__ == "__main__":
    main()
