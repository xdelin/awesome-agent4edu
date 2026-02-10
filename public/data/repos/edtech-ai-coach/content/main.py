from load_reports import load_all_reports
from analytics import prepare_llm_context

# Load all submissions
reports = load_all_reports()

# Convert all to context objects
contexts = [prepare_llm_context(report) for report in reports]

