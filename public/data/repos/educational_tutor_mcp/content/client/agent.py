from client.mcp_client import MCPClient
from client.utils import print_markdown

class AITutorAgent:
    """
    Smart assistant that interacts with the MCP server to provide tutoring help.
    Handles logic for explaining, summarizing, flashcards and quizzing.
    """

    def __init__(self):
        self.client = MCPClient()
    
    def explain(self, concept, level = 3):
        """Uses the MCP tool 'explain_concept'."""
        print_markdown(f"Explaining: **{concept}** at level {level}")
        result = self.client.call_tool("explain_concept", {"question": concept, "level": level})
        return result.get("output", "No response from server")
    
    def summarize(self, text):
        """Uses the MCP tool 'summarize_text'."""
        print_markdown("Summarizing text...")
        result = self.client.call_tool("summarize_text", {"text": text})
        return result.get("output", "No response from server")
    
    def flashcards(self, topic):
        """Generates flashcards."""
        print_markdown(f"Generating flashcards for: **{topic}**")
        result = self.client.call_tool("generate_flashcards", {"topic": topic})
        return result.get("output", [])
    
    def quiz(self, topic):
        """Generates a quiz for a topic."""
        print_markdown(f"Creating quiz for: **{topic}**")
        result = self.client.call_tool("quiz_me", {"topic": topic})
        return result.get("output", [])
