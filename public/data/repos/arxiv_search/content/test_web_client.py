import asyncio
import httpx
import logging
import subprocess
import sys
import signal
from dataclasses import dataclass, field
from typing import Optional

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

SERVER_URL = "http://127.0.0.1:8000/mcp"
HEADERS = {
    "Accept": "application/json, text/event-stream",
    "Content-Type": "application/json"
}

async def check_server_ready(client: httpx.AsyncClient):
    """Polls the server until it is ready to accept connections."""
    for _ in range(20):  # Poll for up to 10 seconds
        try:
            response = await client.post(SERVER_URL, json={"jsonrpc": "2.0", "id": 0, "method": "tools/list"}, headers=HEADERS)
            if response.status_code == 200:
                logging.info("Web server is up and running.")
                return True
        except httpx.ConnectError:
            pass
        await asyncio.sleep(0.5)
    logging.error("Web server did not start in time.")
    return False

def parse_sse_response(response_text: str) -> dict:
    """
    Parse Server-Sent Events response and extract JSON data.

    Args:
        response_text: Raw SSE response text.

    Returns:
        Parsed JSON data from the response.

    Raises:
        ValueError: If no valid data event is found.
    """
    import json
    for line in response_text.strip().split('\n'):
        if line.startswith('data:'):
            return json.loads(line[len('data:'):].strip())
    raise ValueError("Did not receive a valid data event from the server.")


async def call_tool(client: httpx.AsyncClient, tool_name: str, params: dict = None) -> dict:
    """
    Helper function to call a tool via JSON-RPC.

    Args:
        client: HTTP client instance.
        tool_name: Name of the tool to call.
        params: Optional parameters for the tool.

    Returns:
        Parsed JSON-RPC response.
    """
    method = "tools/call" if tool_name != "tools/list" else "tools/list"

    payload = {
        "jsonrpc": "2.0",
        "id": 1,
        "method": method,
    }

    if method == "tools/call":
        payload["params"] = {"name": tool_name, "arguments": params or {}}

    response = await client.post(SERVER_URL, json=payload, headers=HEADERS)
    response.raise_for_status()
    return parse_sse_response(response.text)


async def call_background_tool(
    client: httpx.AsyncClient,
    tool_name: str,
    params: dict = None,
    poll_interval: float = 2.0,
    timeout: float = 180.0
) -> str:
    """
    Call a tool that runs as a background task, polling until complete.

    For tools with task=True, the server returns a task ID immediately.
    This function polls for task completion and returns the final result.

    Args:
        client: HTTP client instance.
        tool_name: Name of the tool to call.
        params: Optional parameters for the tool.
        poll_interval: Seconds between status checks.
        timeout: Maximum seconds to wait for completion.

    Returns:
        The tool's result string.

    Raises:
        TimeoutError: If task doesn't complete within timeout.
        ValueError: If task fails or returns unexpected response.
    """
    import json

    # Initial call - may return task info or immediate result
    response_json = await call_tool(client, tool_name, params)

    # Check if this is a background task response
    # Background tasks return a notification with task info
    result = response_json.get('result', {})

    # If we got a direct result (non-background), return it
    if 'structuredContent' in result:
        return result['structuredContent']['result']

    # For background tasks, we need to poll for completion
    # The response format for background tasks includes task metadata
    if 'task' not in result:
        # Try to extract task ID from the response structure
        # FastMCP may use different response formats
        raise ValueError(f"Unexpected response format: {response_json}")

    task_id = result['task']['id']
    logging.info(f"Background task started with ID: {task_id}")

    # Poll for task completion
    elapsed = 0.0
    while elapsed < timeout:
        await asyncio.sleep(poll_interval)
        elapsed += poll_interval

        # Check task status
        status_payload = {
            "jsonrpc": "2.0",
            "id": 1,
            "method": "tasks/get",
            "params": {"id": task_id}
        }
        status_response = await client.post(SERVER_URL, json=status_payload, headers=HEADERS)
        status_response.raise_for_status()
        status_json = parse_sse_response(status_response.text)

        task_result = status_json.get('result', {})
        task_status = task_result.get('status', '')

        logging.info(f"Task status: {task_status} (elapsed: {elapsed:.0f}s)")

        if task_status == 'completed':
            # Extract the result from completed task
            task_output = task_result.get('output', {})
            if 'structuredContent' in task_output:
                return task_output['structuredContent']['result']
            elif 'result' in task_output:
                return task_output['result']
            else:
                return str(task_output)

        elif task_status == 'failed':
            error_msg = task_result.get('error', 'Unknown error')
            raise ValueError(f"Background task failed: {error_msg}")

    raise TimeoutError(f"Background task {task_id} did not complete within {timeout} seconds")


async def main():
    """
    Test client for the mcp-simple-arxiv web server.
    Starts the server, runs tests, and then stops it.

    Returns:
        True if all tests passed, False otherwise.
    """
    server_process = None
    summary = TestSummary()

    try:
        logging.info("Starting web server process...")
        # Start the server as a subprocess
        server_process = subprocess.Popen(
            [sys.executable, "-m", "mcp_simple_arxiv.web_server"],
            stdout=sys.stdout,
            stderr=sys.stderr
        )

        async with httpx.AsyncClient(timeout=120.0) as client:  # Increased for PDF conversion
            if not await check_server_ready(client):
                raise RuntimeError("Could not connect to the web server.")

            # 1. List available tools
            test_name = "tools/list"
            try:
                logging.info(f"\n--- Testing {test_name} ---")
                response_json = await call_tool(client, "tools/list")
                tools = response_json['result']['tools']
                logging.info(f"Found {len(tools)} tools.")
                assert len(tools) == 5
                logging.info(f"✅ {test_name} test PASSED")
                summary.add_result(test_name, passed=True)
            except Exception as e:
                logging.error(f"❌ {test_name} test FAILED: {e}")
                summary.add_result(test_name, passed=False, error_message=str(e))

            # 2. Test search_papers
            test_name = "search_papers"
            try:
                logging.info(f"\n--- Testing {test_name} ---")
                query = "dark matter"
                response_json = await call_tool(client, "search_papers", {"query": query, "max_results": 1})
                result = response_json['result']['structuredContent']['result']
                logging.info(f"Result for '{query}':\n{result}")
                assert "Found" in result and "total results" in result
                logging.info(f"✅ {test_name} test PASSED")
                summary.add_result(test_name, passed=True)
            except Exception as e:
                logging.error(f"❌ {test_name} test FAILED: {e}")
                summary.add_result(test_name, passed=False, error_message=str(e))

            # Test search_papers sorting options
            test_name = "search_papers_sorting"
            try:
                logging.info(f"\n--- Testing {test_name} ---")
                response_json = await call_tool(client, "search_papers", {
                    "query": "neural networks",
                    "max_results": 2,
                    "sort_by": "relevance"
                })
                result = response_json['result']['structuredContent']['result']
                logging.info(f"Result:\n{result}")
                assert "Found" in result and "total results" in result

                response_json = await call_tool(client, "search_papers", {
                    "query": "neural networks",
                    "max_results": 2,
                    "sort_order": "ascending"
                })
                result = response_json['result']['structuredContent']['result']
                logging.info(f"Result:\n{result}")
                assert "Found" in result and "total results" in result

                # Test invalid sort_by
                response_json = await call_tool(client, "search_papers", {
                    "query": "test",
                    "sort_by": "invalid"
                })
                result = response_json['result']['structuredContent']['result']
                logging.info(f"Invalid sort_by result:\n{result}")
                assert "Invalid sort_by value" in result

                # Test invalid sort_order
                response_json = await call_tool(client, "search_papers", {
                    "query": "test",
                    "sort_order": "invalid"
                })
                result = response_json['result']['structuredContent']['result']
                logging.info(f"Invalid sort_order result:\n{result}")
                assert "Invalid sort_order value" in result
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
                response_json = await call_tool(client, "search_papers", {
                    "query": "ti:transformer",
                    "max_results": 2,
                    "date_from": "2024-01-01",
                    "date_to": "2024-12-31"
                })
                result = response_json['result']['structuredContent']['result']
                logging.info(f"Date range result:\n{result}")
                assert "Found" in result and "total results" in result

                # Test with only date_from
                response_json = await call_tool(client, "search_papers", {
                    "query": "ti:quantum",
                    "max_results": 2,
                    "date_from": "2024-06-01"
                })
                result = response_json['result']['structuredContent']['result']
                logging.info(f"Date from only result:\n{result}")
                assert "Found" in result and "total results" in result

                # Test invalid date format
                response_json = await call_tool(client, "search_papers", {
                    "query": "test",
                    "date_from": "01-01-2024"
                })
                result = response_json['result']['structuredContent']['result']
                logging.info(f"Invalid date format result:\n{result}")
                assert "Invalid date_from format" in result

                # Test date_from after date_to
                response_json = await call_tool(client, "search_papers", {
                    "query": "test",
                    "date_from": "2024-12-31",
                    "date_to": "2024-01-01"
                })
                result = response_json['result']['structuredContent']['result']
                logging.info(f"Invalid date range result:\n{result}")
                assert "date_from cannot be after date_to" in result

                logging.info(f"✅ {test_name} test PASSED")
                summary.add_result(test_name, passed=True)
            except Exception as e:
                logging.error(f"❌ {test_name} test FAILED: {e}")
                summary.add_result(test_name, passed=False, error_message=str(e))

            # 3. Test get_paper_data
            test_name = "get_paper_data"
            try:
                logging.info(f"\n--- Testing {test_name} ---")
                paper_id = "0808.3772"  # Using the same ID as the stdio test for consistency
                response_json = await call_tool(client, "get_paper_data", {"paper_id": paper_id})
                result = response_json['result']['structuredContent']['result']
                logging.info(f"Result for paper '{paper_id}':\n{result}")
                assert "A common mass scale for satellite galaxies of the Milky Way" in result
                logging.info(f"✅ {test_name} test PASSED")
                summary.add_result(test_name, passed=True)
            except Exception as e:
                logging.error(f"❌ {test_name} test FAILED: {e}")
                summary.add_result(test_name, passed=False, error_message=str(e))

            # 4. Test get_full_paper_text (background task - takes 30-90 seconds)
            test_name = "get_full_paper_text"
            try:
                logging.info(f"\n--- Testing {test_name} ---")
                paper_id = "0808.3772"  # Same paper, relatively short
                logging.info(f"Calling get_full_paper_text with paper_id: '{paper_id}'")
                logging.info("(This runs as a background task - polling for completion...)")
                # Use background task helper which handles polling
                result = await call_background_tool(
                    client,
                    "get_full_paper_text",
                    {"paper_id": paper_id},
                    poll_interval=3.0,
                    timeout=180.0
                )
                # Check that we got markdown content back (should contain the title)
                assert "common mass scale" in result.lower() or "satellite galaxies" in result.lower()
                logging.info(f"Result length: {len(result)} characters")
                logging.info(f"✅ {test_name} test PASSED")
                summary.add_result(test_name, passed=True)
            except Exception as e:
                logging.error(f"❌ {test_name} test FAILED: {e}")
                summary.add_result(test_name, passed=False, error_message=str(e))

            # 5. Test list_categories
            test_name = "list_categories"
            try:
                logging.info(f"\n--- Testing {test_name} ---")
                response_json = await call_tool(client, "list_categories")
                result = response_json['result']['structuredContent']['result']
                logging.info("Result snippet:\n" + result[:200] + "...")
                assert "arXiv Categories" in result
                logging.info(f"✅ {test_name} test PASSED")
                summary.add_result(test_name, passed=True)
            except Exception as e:
                logging.error(f"❌ {test_name} test FAILED: {e}")
                summary.add_result(test_name, passed=False, error_message=str(e))

            # 6. Test update_categories
            test_name = "update_categories"
            try:
                logging.info(f"\n--- Testing {test_name} ---")
                response_json = await call_tool(client, "update_categories")
                result = response_json['result']['structuredContent']['result']
                logging.info(f"Result:\n{result}")
                assert "Successfully updated category taxonomy" in result
                logging.info(f"✅ {test_name} test PASSED")
                summary.add_result(test_name, passed=True)
            except Exception as e:
                logging.error(f"❌ {test_name} test FAILED: {e}")
                summary.add_result(test_name, passed=False, error_message=str(e))

    except Exception as e:
        logging.error(f"An error occurred during testing: {e}", exc_info=True)
    finally:
        if server_process:
            logging.info("\nStopping web server process...")
            server_process.send_signal(signal.SIGINT)  # Send Ctrl+C
            try:
                server_process.wait(timeout=10)
                logging.info("Web server stopped gracefully.")
            except subprocess.TimeoutExpired:
                logging.warning("Web server did not stop gracefully, killing.")
                server_process.kill()

    # Print final summary
    all_passed = summary.print_summary()
    return all_passed

if __name__ == "__main__":
    success = asyncio.run(main())
    sys.exit(0 if success else 1) 