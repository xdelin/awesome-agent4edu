#!/usr/bin/env python3
"""
Census Data Fetching Module

This module provides functions to fetch demographic data from the US Census API.
"""

import os
import sys
from typing import Dict, List, Optional, Tuple
import pandas as pd
from census import Census
from us import states


class CensusDataFetcher:
    """Fetches demographic data from US Census API."""

    # Common variable mappings
    VARIABLES = {
        'population': {
            'total': 'B01003_001E',
            'density': 'B01003_001E',  # Will calculate with area
        },
        'age': {
            'median': 'B01002_001E',
            'total': 'B01001_001E',
            'under_18': 'B01001_003E',
            'over_65': 'B01001_020E',
        },
        'race': {
            'white': 'B02001_002E',
            'black': 'B02001_003E',
            'native': 'B02001_004E',
            'asian': 'B02001_005E',
            'pacific': 'B02001_006E',
            'other': 'B02001_007E',
            'two_or_more': 'B02001_008E',
            'hispanic': 'B03003_003E',
        },
        'education': {
            'total_25_over': 'B15003_001E',
            'high_school': 'B15003_017E',
            'some_college': 'B15003_019E',
            'associates': 'B15003_021E',
            'bachelors': 'B15003_022E',
            'masters': 'B15003_023E',
            'professional': 'B15003_024E',
            'doctorate': 'B15003_025E',
        },
        'income': {
            'median_household': 'B19013_001E',
            'per_capita': 'B19301_001E',
            'mean_household': 'B19025_001E',
        },
        'housing': {
            'total_units': 'B25001_001E',
            'median_value': 'B25077_001E',
            'median_rent': 'B25064_001E',
            'owner_occupied': 'B25003_002E',
            'renter_occupied': 'B25003_003E',
        },
        'employment': {
            'in_labor_force': 'B23025_002E',
            'employed': 'B23025_004E',
            'unemployed': 'B23025_005E',
        }
    }

    def __init__(self, api_key: Optional[str] = None):
        """
        Initialize the Census Data Fetcher.

        Args:
            api_key: Census API key. If None, will try to read from
                    CENSUS_API_KEY environment variable.
        """
        self.api_key = api_key or os.environ.get('CENSUS_API_KEY')
        if not self.api_key:
            raise ValueError(
                "Census API key required. Get one at "
                "https://api.census.gov/data/key_signup.html\n"
                "Set it via CENSUS_API_KEY environment variable or pass to constructor."
            )
        self.census = Census(self.api_key)

    def get_variables_for_category(self, category: str) -> List[Tuple[str, str]]:
        """
        Get variable codes for a demographic category.

        Args:
            category: Category name (e.g., 'population', 'age', 'race')

        Returns:
            List of (variable_name, variable_code) tuples
        """
        if category not in self.VARIABLES:
            available = ', '.join(self.VARIABLES.keys())
            raise ValueError(
                f"Unknown category '{category}'. Available: {available}"
            )

        return [
            (name, code)
            for name, code in self.VARIABLES[category].items()
        ]

    def fetch_state_data(
        self,
        variables: List[str],
        year: int = 2021,
        state_fips: str = '*'
    ) -> pd.DataFrame:
        """
        Fetch state-level data.

        Args:
            variables: List of Census variable codes
            year: Data year (default: 2021)
            state_fips: State FIPS code or '*' for all states

        Returns:
            DataFrame with state data
        """
        fields = ('NAME',) + tuple(variables)
        data = self.census.acs5.state(
            fields=fields,
            state_fips=state_fips,
            year=year
        )

        df = pd.DataFrame(data)
        df['GEOID'] = df['state']
        return df

    def fetch_county_data(
        self,
        variables: List[str],
        year: int = 2021,
        state_fips: str = '*',
        county_fips: str = '*'
    ) -> pd.DataFrame:
        """
        Fetch county-level data.

        Args:
            variables: List of Census variable codes
            year: Data year (default: 2021)
            state_fips: State FIPS code or '*' for all states
            county_fips: County FIPS code or '*' for all counties

        Returns:
            DataFrame with county data
        """
        fields = ('NAME',) + tuple(variables)
        data = self.census.acs5.state_county(
            fields=fields,
            state_fips=state_fips,
            county_fips=county_fips,
            year=year
        )

        df = pd.DataFrame(data)
        df['GEOID'] = df['state'] + df['county']
        return df

    def fetch_tract_data(
        self,
        variables: List[str],
        state_fips: str,
        county_fips: str = '*',
        year: int = 2021
    ) -> pd.DataFrame:
        """
        Fetch census tract-level data.

        Args:
            variables: List of Census variable codes
            state_fips: State FIPS code (required)
            county_fips: County FIPS code or '*' for all counties in state
            year: Data year (default: 2021)

        Returns:
            DataFrame with tract data
        """
        fields = ('NAME',) + tuple(variables)
        data = self.census.acs5.state_county_tract(
            fields=fields,
            state_fips=state_fips,
            county_fips=county_fips,
            tract='*',
            year=year
        )

        df = pd.DataFrame(data)
        df['GEOID'] = df['state'] + df['county'] + df['tract']
        return df

    def calculate_percentage(
        self,
        df: pd.DataFrame,
        numerator_col: str,
        denominator_col: str,
        result_col: str
    ) -> pd.DataFrame:
        """
        Calculate percentage from two columns.

        Args:
            df: Input DataFrame
            numerator_col: Column to use as numerator
            denominator_col: Column to use as denominator
            result_col: Name for result column

        Returns:
            DataFrame with new percentage column
        """
        df[result_col] = (
            pd.to_numeric(df[numerator_col], errors='coerce') /
            pd.to_numeric(df[denominator_col], errors='coerce') * 100
        )
        return df

    def get_state_fips(self, state_name: str) -> str:
        """
        Get FIPS code for a state.

        Args:
            state_name: State name or abbreviation

        Returns:
            State FIPS code
        """
        state = states.lookup(state_name)
        if state is None:
            raise ValueError(f"Unknown state: {state_name}")
        return state.fips


def main():
    """Example usage of CensusDataFetcher."""

    # Check for API key
    api_key = os.environ.get('CENSUS_API_KEY')
    if not api_key:
        print("Error: CENSUS_API_KEY environment variable not set")
        print("Get an API key at: https://api.census.gov/data/key_signup.html")
        print("Then set it: export CENSUS_API_KEY='your_key_here'")
        sys.exit(1)

    # Initialize fetcher
    fetcher = CensusDataFetcher(api_key)

    # Example: Fetch population data for all states
    print("Fetching state population data...")
    pop_vars = ['B01003_001E']  # Total population
    state_data = fetcher.fetch_state_data(pop_vars, year=2021)
    print(f"Fetched data for {len(state_data)} states")
    print(state_data.head())

    # Example: Fetch county data for California
    print("\nFetching California county data...")
    ca_fips = fetcher.get_state_fips('CA')
    county_data = fetcher.fetch_county_data(
        pop_vars,
        year=2021,
        state_fips=ca_fips
    )
    print(f"Fetched data for {len(county_data)} counties")
    print(county_data.head())


if __name__ == '__main__':
    main()
