#!/usr/bin/env python3

import requests
import json
from typing import Dict, List, Optional
from urllib.parse import urljoin
import time

class GNoMEExplorer:
    BASE_URL = "https://optimade-gnome.odbx.science/v1"
    
    def __init__(self):
        self.session = requests.Session()
        self.session.headers.update({
            "Accept": "application/vnd.api+json",
            "User-Agent": "MaterialsMCP/0.1.0"
        })
    
    def get_info(self) -> Dict:
        """Fetch provider information from the /info endpoint."""
        response = self.session.get(urljoin(self.BASE_URL, "info"))
        response.raise_for_status()
        return response.json()
    
    def query_structures(self, filter_str: str, page_size: int = 10) -> Dict:
        """Query structures with a given filter string."""
        params = {
            "filter": filter_str,
            "page_limit": page_size,
            "response_format": "json"
        }
        response = self.session.get(
            urljoin(self.BASE_URL, "structures"),
            params=params
        )
        response.raise_for_status()
        return response.json()
    
    def get_structure_details(self, structure_id: str) -> Dict:
        """Get full details for a specific structure."""
        response = self.session.get(
            urljoin(self.BASE_URL, f"structures/{structure_id}")
        )
        response.raise_for_status()
        return response.json()

def print_json(data: Dict, indent: int = 2) -> None:
    """Pretty print JSON data."""
    print(json.dumps(data, indent=indent))

def explore_api():
    explorer = GNoMEExplorer()
    
    print("=== Fetching Provider Information ===")
    try:
        info = explorer.get_info()
        print("\nProvider Info:")
        print(f"Provider: {info.get('meta', {}).get('provider', {}).get('name', 'Unknown')}")
        print(f"API Version: {info.get('meta', {}).get('api_version', 'Unknown')}")
        print("\nAvailable Properties:")
        properties = info.get('meta', {}).get('properties', {})
        for prop, details in properties.items():
            print(f"- {prop}: {details.get('description', 'No description')}")
    except Exception as e:
        print(f"Error fetching info: {e}")
        return

    print("\n=== Testing Various Queries ===")
    
    # List of queries to test
    queries = [
        {
            "name": "Stable Binary Compounds (Low Formation Energy)",
            "filter": 'nelements = 2 AND _gnome_formation_energy_per_atom < -1.0'
        },
        {
            "name": "Stable Compounds with Band Gap",
            "filter": '_gnome_formation_energy_per_atom < -2.0 AND _gnome_bandgap > 0.0'
        },
        {
            "name": "Very Stable Compounds",
            "filter": '_gnome_formation_energy_per_atom < -3.0 AND _gnome_decomposition_energy_per_atom < -0.1'
        },
        {
            "name": "Wide Band Gap Materials",
            "filter": '_gnome_bandgap > 3.0 AND _gnome_formation_energy_per_atom < -2.0'
        },
        {
            "name": "GaN Compounds (Alternative Query)",
            "filter": 'chemical_formula_reduced CONTAINS "Ga" AND chemical_formula_reduced CONTAINS "N"'
        }
    ]
    
    for query in queries:
        print(f"\nTesting Query: {query['name']}")
        print(f"Filter: {query['filter']}")
        try:
            start_time = time.time()
            results = explorer.query_structures(query['filter'])
            end_time = time.time()
            
            data = results.get('data', [])
            print(f"Found {len(data)} structures")
            print(f"Query took {end_time - start_time:.2f} seconds")
            
            if data:
                # Print first structure's properties
                print("\nAvailable properties in first structure:")
                attrs = data[0].get('attributes', {})
                for key in attrs.keys():
                    print(f"- {key}")
                
                # Get full details of first structure
                structure_id = data[0]['id']
                print(f"\nFetching full details for structure {structure_id}")
                details = explorer.get_structure_details(structure_id)
                print("\nFull structure details:")
                print_json(details)
                
        except Exception as e:
            print(f"Error executing query: {e}")

if __name__ == "__main__":
    explore_api() 