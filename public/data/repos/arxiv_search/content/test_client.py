import asyncio
import logging
import sys
from dataclasses import dataclass, field
from typing import Optional
from fastmcp.client import Client, StdioTransport

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


@dataclass
class TestResult:
    """Stores the result of a single test."""
    name: str
    passed: bool
    error_message: Optional[str] = None


@dataclass
class TestSummary:
    """Collects and summarizes test results."""
    results: list = field(default_factory=list)

    def add_result(self, name: str, passed: bool, error_message: Optional[str] = None) -> None:
        """Add a test result to the summary."""
        self.results.append(TestResult(name=name, passed=passed, error_message=error_message))

    def print_summary(self) -> bool:
        """
        Print a summary of all test results.

        Returns:
            True if all tests passed, False otherwise.
        """
        total = len(self.results)
        passed = sum(1 for r in self.results if r.passed)
        failed = total - passed

        logging.info("\n" + "=" * 60)
        logging.info("TEST SUMMARY")
        logging.info("=" * 60)
        logging.info(f"Total tests: {total}")
        logging.info(f"Passed:      {passed}")
        logging.info(f"Failed:      {failed}")

        if failed > 0:
            logging.info("\nFailed tests:")
            for result in self.results:
                if not result.passed:
                    logging.info(f"  ❌ {result.name}")
                    if result.error_message:
                        logging.info(f"     Error: {result.error_message}")

        logging.info("-" * 60)
        if failed == 0:
            logging.info("✅ ALL TESTS PASSED")
        else:
            logging.info(f"❌ {failed} TEST(S) FAILED")
        logging.info("=" * 60)

        return failed == 0

async def main():
    """
    Test client for the mcp-simple-arxiv server.
    Connects to the stdio server, lists tools, and calls each tool to verify functionality.

    Returns:
        True if all tests passed, False otherwise.
    """
    logging.info("Starting test client for mcp-simple-arxiv...")
    summary = TestSummary()

    # Configure the stdio transport to run the server as a module
    # Use sys.executable to ensure we use the same Python interpreter as the test
    transport = StdioTransport(
        command=sys.executable,
        args=["-m", "mcp_simple_arxiv"]
    )

    # Create a client with the transport
    client = Client(transport)

    async with client:
        # 1. List available tools
        test_name = "tools/list"
        try:
            logging.info(f"--- Testing {test_name} ---")
            tools = await client.list_tools()
            logging.info(f"Found {len(tools)} tools:")
            for tool in tools:
                logging.info(f"- {tool.name}: {tool.description.splitlines()[0]}")
            assert len(tools) == 5, "Expected 5 tools"
            logging.info(f"✅ {test_name} test PASSED")
            summary.add_result(test_name, passed=True)
        except Exception as e:
            logging.error(f"❌ {test_name} test FAILED: {e}")
            summary.add_result(test_name, passed=False, error_message=str(e))

        # 2. Test search_papers
        test_name = "search_papers"
        try:
            logging.info(f"\n--- Testing {test_name} ---")
            query = "electron"
            logging.info(f"Calling search_papers with query: '{query}'")
            result = await client.call_tool("search_papers", {"query": query, "max_results": 2})
            logging.info(f"Result:\n{result.data}")
            assert "Found" in result.data and "total results" in result.data
            logging.info(f"✅ {test_name} test PASSED")
            summary.add_result(test_name, passed=True)
        except Exception as e:
            logging.error(f"❌ {test_name} test FAILED: {e}")
            summary.add_result(test_name, passed=False, error_message=str(e))

        # Test search_papers sorting options
        test_name = "search_papers_sorting"
        try:
            logging.info(f"\n--- Testing {test_name} ---")
            result = await client.call_tool("search_papers", {
                "query": "neural networks",
                "max_results": 2,
                "sort_by": "relevance"
            })
            logging.info(f"Result:\n{result.data}")
            assert "Found" in result.data and "total results" in result.data

            result = await client.call_tool("search_papers", {
                "query": "neural networks",
                "max_results": 2,
                "sort_order": "ascending"
            })
            logging.info(f"Result:\n{result.data}")
            assert "Found" in result.data and "total results" in result.data

            # Test invalid sort_by
            result = await client.call_tool("search_papers", {
                "query": "test",
                "sort_by": "invalid"
            })
            logging.info(f"Invalid sort_by result:\n{result.data}")
            assert "Invalid sort_by value" in result.data

            # Test invalid sort_order
            result = await client.call_tool("search_papers", {
                "query": "test",
                "sort_order": "invalid"
            })
            logging.info(f"Invalid sort_order result:\n{result.data}")
            assert "Invalid sort_order value" in result.data
            logging.info(f"✅ {test_name} test PASSED")
            summary.add_result(test_name, passed=True)
        except Exception as e:
            logging.error(f"❌ {test_name} test FAILED: {e}")
            summary.add_result(test_name, passed=False, error_message=str(e))

        # Test search_papers date filtering
        test_name = "search_papers_date_filtering"
        try:
            logging.info(f"\n--- Testing {test_name} ---")

            # Test with both dates
            result = await client.call_tool("search_papers", {
                "query": "ti:transformer",
                "max_results": 2,
                "date_from": "2024-01-01",
                "date_to": "2024-12-31"
            })
            logging.info(f"Date range result:\n{result.data}")
            assert "Found" in result.data and "total results" in result.data

            # Test with only date_from
            result = await client.call_tool("search_papers", {
                "query": "ti:quantum",
                "max_results": 2,
                "date_from": "2024-06-01"
            })
            logging.info(f"Date from only result:\n{result.data}")
            assert "Found" in result.data and "total results" in result.data

            # Test invalid date format
            result = await client.call_tool("search_papers", {
                "query": "test",
                "date_from": "01-01-2024"
            })
            logging.info(f"Invalid date format result:\n{result.data}")
            assert "Invalid date_from format" in result.data

            # Test date_from after date_to
            result = await client.call_tool("search_papers", {
                "query": "test",
                "date_from": "2024-12-31",
                "date_to": "2024-01-01"
            })
            logging.info(f"Invalid date range result:\n{result.data}")
            assert "date_from cannot be after date_to" in result.data

            logging.info(f"✅ {test_name} test PASSED")
            summary.add_result(test_name, passed=True)
        except Exception as e:
            logging.error(f"❌ {test_name} test FAILED: {e}")
            summary.add_result(test_name, passed=False, error_message=str(e))

        # 3. Test get_paper_data
        test_name = "get_paper_data"
        try:
            logging.info(f"\n--- Testing {test_name} ---")
            paper_id = "0808.3772"  # A known paper
            logging.info(f"Calling get_paper_data with paper_id: '{paper_id}'")
            result = await client.call_tool("get_paper_data", {"paper_id": paper_id})
            logging.info(f"Result:\n{result.data}")
            assert "A common mass scale for satellite galaxies of the Milky Way" in result.data
            logging.info(f"✅ {test_name} test PASSED")
            summary.add_result(test_name, passed=True)
        except Exception as e:
            logging.error(f"❌ {test_name} test FAILED: {e}")
            summary.add_result(test_name, passed=False, error_message=str(e))

        # 4. Test get_full_paper_text (this takes 30-90 seconds)
        test_name = "get_full_paper_text"
        try:
            logging.info(f"\n--- Testing {test_name} ---")
            paper_id = "0808.3772"  # Same paper, relatively short
            logging.info(f"Calling get_full_paper_text with paper_id: '{paper_id}'")
            logging.info("(This may take 30-90 seconds as it downloads and converts the PDF...)")
            result = await client.call_tool("get_full_paper_text", {"paper_id": paper_id})
            # Check that we got markdown content back (should contain the title)
            assert "common mass scale" in result.data.lower() or "satellite galaxies" in result.data.lower()
            logging.info(f"Result length: {len(result.data)} characters")
            logging.info(f"✅ {test_name} test PASSED")
            summary.add_result(test_name, passed=True)
        except Exception as e:
            logging.error(f"❌ {test_name} test FAILED: {e}")
            summary.add_result(test_name, passed=False, error_message=str(e))

        # 5. Test list_categories (no filter)
        test_name = "list_categories"
        try:
            logging.info(f"\n--- Testing {test_name} ---")
            logging.info("Calling list_categories without a filter...")
            result = await client.call_tool("list_categories")
            logging.info(f"Result snippet:\n{result.data[:300]}...")
            assert "arXiv Categories" in result.data
            logging.info(f"✅ {test_name} (no filter) test PASSED")
            summary.add_result(test_name, passed=True)
        except Exception as e:
            logging.error(f"❌ {test_name} test FAILED: {e}")
            summary.add_result(test_name, passed=False, error_message=str(e))

        # 5b. Test list_categories (with filter)
        test_name = "list_categories_filtered"
        try:
            logging.info(f"\n--- Testing {test_name} ---")
            logging.info("Calling list_categories with filter 'cs'...")
            result = await client.call_tool("list_categories", {"primary_category": "cs"})
            logging.info(f"Result snippet:\n{result.data[:300]}...")
            assert "cs: Computer Science" in result.data
            assert "math: Mathematics" not in result.data
            logging.info(f"✅ {test_name} test PASSED")
            summary.add_result(test_name, passed=True)
        except Exception as e:
            logging.error(f"❌ {test_name} test FAILED: {e}")
            summary.add_result(test_name, passed=False, error_message=str(e))

        # 6. Test update_categories
        test_name = "update_categories"
        try:
            logging.info(f"\n--- Testing {test_name} ---")
            logging.info("Calling update_categories...")
            result = await client.call_tool("update_categories")
            logging.info(f"Result:\n{result.data}")
            assert "Successfully updated category taxonomy" in result.data
            logging.info(f"✅ {test_name} test PASSED")
            summary.add_result(test_name, passed=True)
        except Exception as e:
            logging.error(f"❌ {test_name} test FAILED: {e}")
            summary.add_result(test_name, passed=False, error_message=str(e))

    # Print final summary
    all_passed = summary.print_summary()
    return all_passed

if __name__ == "__main__":
    success = asyncio.run(main())
    sys.exit(0 if success else 1) 