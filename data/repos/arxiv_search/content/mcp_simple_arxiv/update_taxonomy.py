"""
Script to verify and optionally update arXiv category taxonomy.
"""

import json
import logging
from pathlib import Path
from .categories import CATEGORIES

logger = logging.getLogger(__name__)

TAXONOMY_FILE = Path(__file__).parent / "taxonomy.json"


def update_taxonomy_file() -> dict:
    """
    Create taxonomy.json from the built-in categories.

    Returns:
        dict: The taxonomy dictionary containing all arXiv categories.
    """
    logger.info(f"Creating taxonomy file at {TAXONOMY_FILE}...")
    with open(TAXONOMY_FILE, 'w', encoding='utf-8') as f:
        json.dump(CATEGORIES, f, indent=2, ensure_ascii=False)
    logger.info("Done!")
    return CATEGORIES


def load_taxonomy() -> dict:
    """
    Load taxonomy from the JSON file.

    If the file doesn't exist, it will be created from built-in categories.

    Returns:
        dict: The taxonomy dictionary containing all arXiv categories.
    """
    if not TAXONOMY_FILE.exists():
        logger.info(f"Taxonomy file not found at {TAXONOMY_FILE}, creating it...")
        return update_taxonomy_file()

    logger.debug(f"Loading taxonomy from {TAXONOMY_FILE}")
    with open(TAXONOMY_FILE, 'r', encoding='utf-8') as f:
        return json.load(f)


if __name__ == "__main__":
    # When run directly, create/update the taxonomy file
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
    logger.info("Creating taxonomy file from built-in categories...")
    taxonomy = update_taxonomy_file()
    logger.info(f"Created taxonomy with {len(taxonomy)} primary categories:")
    for primary, data in taxonomy.items():
        logger.info(f"- {primary}: {data['name']} ({len(data['subcategories'])} subcategories)")