# load_reports.py
"""
Utility to load all student submission JSON files from the 'data/' folder.
"""

import json
from pathlib import Path
from typing import List, Dict
# load_reports.py

def load_single_report(filepath: str) -> dict:
    path = Path(filepath)
    if not path.exists():
        raise FileNotFoundError(f"File {filepath} not found.")
    with open(path, "r", encoding="utf-8") as f:
        data = json.load(f)
        if isinstance(data, list) and len(data) == 1 and isinstance(data[0], dict):
            data = data[0]
    return data


def load_all_reports(folder: str = "data") -> List[Dict]:
    """
    Load all JSON reports from a folder (assumes 4 files).
    
    Parameters:
    -----------
    folder : str
        Folder path where the JSON files are stored.
    
    Returns:
    --------
    List[Dict]
        List of parsed JSON dictionaries.
    """
    path = Path(folder)
    if not path.exists():
        raise FileNotFoundError(f"Folder {folder} not found.")

    files = sorted(path.glob("*.json"))
    if len(files) != 4:
        raise ValueError(f"Expected 4 JSON files, found {len(files)} in {folder}.")

    reports = []
    for file in files:
        with open(file, "r", encoding="utf-8") as f:
            data = json.load(f)
            # ✅ Fix if it's a list of 1 dict
            if isinstance(data, list) and len(data) == 1 and isinstance(data[0], dict):
                data = data[0]
            reports.append(data)

    return reports


# CLI usage (optional for dev testing)
if __name__ == "__main__":
    reports = load_all_reports()
    print(f"✅ Loaded {len(reports)} reports:")
    for i, r in enumerate(reports, 1):
        print(f"  Report {i} - Score: {r['totalMarkScored']} / {r['test']['totalMarks']}")
