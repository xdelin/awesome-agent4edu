from rich.console import Console

console = Console()

def print_markdown(text):
    """Nicely formats text in Markdown style for terminal output."""
    console.print(text, style = "bold cyan")