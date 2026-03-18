"""
MCP server for accessing arXiv papers.
"""

import sys
sys.stdout.reconfigure(encoding='utf-8')
sys.stdin.reconfigure(encoding='utf-8')
from importlib.metadata import version
import asyncio
import logging
from datetime import datetime, date
from typing import Optional, Tuple

from fastmcp import FastMCP

from .arxiv_client import ArxivClient, SearchResult, SortBy, SortOrder
from .update_taxonomy import load_taxonomy, update_taxonomy_file

_version = version("mcp-simple-arxiv")
logger = logging.getLogger(__name__)

def get_first_sentence(text: str, max_len: int = 200) -> str:
    """
    Extract the first sentence from text, limiting length.

    Looks for common sentence endings (period, exclamation, question mark).
    If no sentence ending is found within max_len characters, truncates
    the text and appends ellipsis.

    Args:
        text: The input text to extract from.
        max_len: Maximum length of the returned string.

    Returns:
        The first sentence or a truncated version of the text.
    """
    # Look for common sentence endings
    for end in ['. ', '! ', '? ']:
        pos = text.find(end)
        if pos != -1 and pos < max_len:
            return text[:pos + 1]
    # If no sentence ending found, just take first max_len chars
    if len(text) > max_len:
        return text[:max_len].rstrip() + '...'
    return text


def parse_date_filter(
    date_from: Optional[str] = None,
    date_to: Optional[str] = None
) -> Tuple[str, Optional[str]]:
    """
    Parse and validate date range parameters for arXiv API.

    Converts user-friendly YYYY-MM-DD dates to arXiv's required format
    (YYYYMMDDTTTT where TTTT is 24-hour time in GMT).

    Args:
        date_from: Start date in YYYY-MM-DD format (uses 00:00 GMT).
        date_to: End date in YYYY-MM-DD format (uses 23:59 GMT).

    Returns:
        Tuple of (filter_string, error_message).
        If successful: (filter_string, None)
        If error: ("", error_message)
        If no dates provided: ("", None)
    """
    if not date_from and not date_to:
        return ("", None)

    # Validate and parse dates
    arxiv_date_from = None
    arxiv_date_to = None

    if date_from:
        try:
            parsed = datetime.strptime(date_from, "%Y-%m-%d")
            # arXiv format: YYYYMMDD + time (0000 for start of day)
            arxiv_date_from = parsed.strftime("%Y%m%d") + "0000"
        except ValueError:
            return ("", f"Invalid date_from format: '{date_from}'. Expected YYYY-MM-DD.")

    if date_to:
        try:
            parsed = datetime.strptime(date_to, "%Y-%m-%d")
            # arXiv format: YYYYMMDD + time (2359 for end of day)
            arxiv_date_to = parsed.strftime("%Y%m%d") + "2359"
        except ValueError:
            return ("", f"Invalid date_to format: '{date_to}'. Expected YYYY-MM-DD.")

    # Set defaults for open-ended ranges
    if arxiv_date_from and not arxiv_date_to:
        # From date_from to today
        arxiv_date_to = date.today().strftime("%Y%m%d") + "2359"
    elif arxiv_date_to and not arxiv_date_from:
        # From arXiv founding (August 1991) to date_to
        arxiv_date_from = "199108010000"

    # Validate date order
    if arxiv_date_from and arxiv_date_to:
        if arxiv_date_from > arxiv_date_to:
            return ("", "date_from cannot be after date_to.")

    # Build arXiv filter string
    filter_str = f"submittedDate:[{arxiv_date_from} TO {arxiv_date_to}]"
    return (filter_str, None)


def create_app() -> FastMCP:
    """
    Create and configure the FastMCP application instance.

    This factory function creates the MCP server and registers all tools
    for interacting with arXiv: search_papers, get_paper_data,
    get_full_paper_text, list_categories, and update_categories.

    Returns:
        A configured FastMCP application instance ready to run.
    """
    app = FastMCP("arxiv-server", version=_version)
    arxiv_client = ArxivClient()

    @app.tool(
        annotations={
            "title": "Search arXiv Papers",
            "readOnlyHint": True,
            "openWorldHint": True
        }
    )
    async def search_papers(
        query: str,
        max_results: int = 10,
        sort_by: str = "submitted_date",
        sort_order: str = "descending",
        date_from: str = None,
        date_to: str = None
    ) -> str:
        """
Search for papers on arXiv.

IMPORTANT - DEFAULT BEHAVIOR WARNING:
ArXiv treats space-separated words as OR by default, returning papers matching
ANY word. This often returns thousands of irrelevant results. Use field
prefixes (especially ti:) for precise searches.

SEARCH STRATEGY (in order of precision):
1. Start with ti: (title) searches - fastest and most relevant results
2. Add cat: (category) to filter by field - use list_categories tool first!
3. Use au: (author) when you know specific researchers
4. Combine multiple terms with AND for best results
5. Avoid plain keyword searches without field prefixes

QUERY OPERATORS:
- ti:"text"   - Search in title only (RECOMMENDED FOR PRECISION)
- abs:"text"  - Search in abstract
- au:"name"   - Search by author
- cat:CODE    - Filter by category (e.g., cat:cs.AI, cat:quant-ph)
- Combine with: AND, OR, ANDNOT

EXAMPLES (from most to least precise):
- ti:"neural networks" AND cat:cs.AI     - Title phrase + category (BEST)
- ti:"deep learning" AND au:bengio       - Title + author
- cat:cs.AI AND ti:transformer           - Category + title keyword
- ti:"machine learning"                  - Title phrase only
- "machine learning"                     - All fields (broad, use sparingly)

TROUBLESHOOTING - Too many irrelevant results?
1. Use ti:"exact phrase" instead of bare keywords
2. Add cat:CATEGORY to filter by field (run list_categories first)
3. Use AND to combine multiple specific terms
4. Avoid generic terms without ti: or cat: prefixes

DATE FILTERING:
Filter papers by submission date using date_from and/or date_to parameters.
- Papers from 2024: date_from="2024-01-01", date_to="2024-12-31"
- Recent papers (2025 onwards): date_from="2025-01-01"
- Historical papers (before 2020): date_to="2019-12-31"

Args:
    query: Search query string (use field prefixes for precision).
    max_results: Maximum results to return (1-100, default 10).
    sort_by: Sort field - "submitted_date", "updated_date", or "relevance".
    sort_order: Sort direction - "descending" or "ascending".
    date_from: Filter papers submitted on or after this date (YYYY-MM-DD format).
    date_to: Filter papers submitted on or before this date (YYYY-MM-DD format).
        """
        max_results = min(max_results, 10)

        # Validate sort_by
        sort_by_mapping = {
            "submitted_date": SortBy.SUBMITTED_DATE,
            "updated_date": SortBy.UPDATED_DATE,
            "relevance": SortBy.RELEVANCE,
        }
        if sort_by not in sort_by_mapping:
            valid_options = ", ".join(sort_by_mapping.keys())
            return f"Invalid sort_by value: '{sort_by}'. Valid options: {valid_options}"
        sort_by_enum = sort_by_mapping[sort_by]

        # Validate sort_order
        sort_order_mapping = {
            "descending": SortOrder.DESCENDING,
            "ascending": SortOrder.ASCENDING,
        }
        if sort_order not in sort_order_mapping:
            valid_options = ", ".join(sort_order_mapping.keys())
            return f"Invalid sort_order value: '{sort_order}'. Valid options: {valid_options}"
        sort_order_enum = sort_order_mapping[sort_order]

        # Build date filter if provided
        date_filter, date_error = parse_date_filter(date_from, date_to)
        if date_error:
            return date_error

        # Combine query with date filter
        if date_filter:
            final_query = f"({query}) AND {date_filter}"
        else:
            final_query = query

        search_result: SearchResult = await arxiv_client.search(
            final_query,
            max_results,
            sort_by=sort_by_enum,
            sort_order=sort_order_enum
        )

        if search_result.total_results == 0:
            return "No papers found matching your query."
        
        # Header with total count
        result = f"Found {search_result.total_results} total results"
        if search_result.results_returned < search_result.total_results:
            result += f", showing first {search_result.results_returned}"
        result += ".\n\n"

        # Format results in a readable way
        for i, paper in enumerate(search_result.papers, 1):
            result += f"{i}. {paper['title']}\n"
            result += f"   Authors: {', '.join(paper['authors'])}\n"
            result += f"   ID: {paper['id']}\n"
            result += f"   Categories: "
            if paper['primary_category']:
                result += f"Primary: {paper['primary_category']}"
            if paper['categories']:
                result += f", Additional: {', '.join(paper['categories'])}"
            result += f"\n   Published: {paper['published']}\n"
            
            # Add first sentence of abstract
            abstract_preview = get_first_sentence(paper['summary'])
            result += f"   Preview: {abstract_preview}\n"
            result += "\n"
        
        return result

    @app.tool(
        annotations={
            "title": "Get arXiv Paper Data",
            "readOnlyHint": True,
            "openWorldHint": True
        }
    )
    async def get_paper_data(paper_id: str) -> str:
        """Get detailed information about a specific paper including abstract and available formats."""
        paper = await arxiv_client.get_paper(paper_id)
        
        # Format paper details in a readable way with clear sections
        result = f"Title: {paper['title']}\n\n"
        
        # Metadata section
        result += "Metadata:\n"
        result += f"- Authors: {', '.join(paper['authors'])}\n"
        result += f"- Published: {paper['published']}\n"
        result += f"- Last Updated: {paper['updated']}\n"
        result += "- Categories: "
        if paper['primary_category']:
            result += f"Primary: {paper['primary_category']}"
        if paper['categories']:
            result += f", Additional: {', '.join(paper['categories'])}"
        result += "\n"
        
        if paper['doi']:
            result += f"- DOI: {paper['doi']}\n"
        if paper["journal_ref"]:
            result += f"- Journal Reference: {paper['journal_ref']}\n"
        
        # Abstract section
        result += "\nAbstract:\n"
        result += paper["summary"]
        result += "\n"
        
        # Access options section
        result += "\nAccess Options:\n"
        result += "- Abstract page: " + paper["abstract_url"] + "\n"
        if paper["html_url"]:  # Add HTML version if available
            result += "- Full text HTML version: " + paper["html_url"] + "\n"
        if paper["pdf_url"]:
            result += "- PDF version: " + paper["pdf_url"] + "\n"
        
        # Additional information section
        if paper["comment"]:
            result += "\nAdditional Information:\n"
            if paper["comment"]:
                result += "- Comment: " + paper["comment"] + "\n"
                
        return result

    @app.tool(
        task=True,
        annotations={
            "title": "Get full paper text as Markdown",
            "readOnlyHint": True,
            "openWorldHint": True
        }
    )
    async def get_full_paper_text(paper_id: str) -> str:
        """Get the full paper text as Markdown
        
        Downloads and converts the paper PDF to Markdown format using Docling.
        This operation takes 30-90 seconds depending on paper length.
        
        Important considerations:
        - Papers can be very large (even 10k-50k+ tokens) and may overwhelm your context window
        - Complex equations and figures will most likely not convert correctly to Markdown
        - Use get_paper_data first to review abstract before fetching full text
        """
        paper = await arxiv_client.get_paper_text_from_pdf(paper_id)
        return paper

    @app.tool(
        annotations={
            "title": "List arXiv Categories",
            "readOnlyHint": True,
            "openWorldHint": False
        }
    )
    def list_categories(primary_category: str = None) -> str:
        """List all available arXiv categories for use with cat: filter in search_papers.

        CALL THIS FIRST before using cat: in search queries to find valid category codes.

        Common categories:
        - cs.AI (Artificial Intelligence)
        - cs.LG (Machine Learning)
        - cs.CL (Computation and Language / NLP)
        - stat.ML (Statistics - Machine Learning)
        - quant-ph (Quantum Physics)
        - q-bio.NC (Quantitative Biology - Neurons and Cognition)

        Args:
            primary_category: Optional filter to show only subcategories of a specific
                primary category (e.g., "cs", "physics", "q-bio").
        """
        try:
            taxonomy = load_taxonomy()
        except Exception as e:
            logger.error(f"Error loading taxonomy: {e}")
            return f"Error loading category taxonomy. Try using update_categories tool to refresh it."

        result = "arXiv Categories:\n\n"
        
        for primary, data in taxonomy.items():
            if primary_category and primary != primary_category:
                continue
                
            result += f"{primary}: {data['name']}\n"
            for code, desc in data['subcategories'].items():
                result += f"  {primary}.{code}: {desc}\n"
            result += "\n"
            
        result += "\nUsage in search:\n"
        result += '- Search in specific category: cat:cs.AI\n'
        result += '- Combine with other terms: "neural networks" AND cat:cs.AI\n'
        result += '- Multiple categories: (cat:cs.AI OR cat:cs.LG)\n'
        result += '\nNote: If categories seem outdated, use the update_categories tool to refresh them.\n'
        
        return result

    @app.tool(
        annotations={
            "title": "Update arXiv Categories",
            "readOnlyHint": False,
            "openWorldHint": True
        }
    )
    def update_categories() -> str:
        """Update the stored category taxonomy by fetching the latest version from arxiv.org"""
        try:
            taxonomy = update_taxonomy_file()
            result = "Successfully updated category taxonomy.\n\n"
            result += f"Found {len(taxonomy)} primary categories:\n"
            for primary, data in taxonomy.items():
                result += f"- {primary}: {data['name']} ({len(data['subcategories'])} subcategories)\n"
            return result
        except Exception as e:
            logger.error(f"Error updating taxonomy: {e}")
            # FastMCP will handle raising this as a proper JSON-RPC error
            raise e

    return app

app = create_app()

def main():
    """Run the MCP server."""
    app.run()

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
