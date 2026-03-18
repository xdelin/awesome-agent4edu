#!/usr/bin/env python3
"""
Create a choropleth map showing median household income by census tract in
Island County, Washington using US Census ACS 5-Year data.
"""

import requests
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from datetime import datetime
import os
import json

def fetch_income_data():
    """
    Fetch median household income data for Island County census tracts from Census API.
    Uses ACS 5-Year estimates.
    Island County FIPS: 53029
    """
    print("Fetching median household income data for Island County census tracts...")

    # ACS 5-Year Detailed Tables
    # B19013_001E: Median household income in the past 12 months (in inflation-adjusted dollars)
    base_url = "https://api.census.gov/data/2022/acs/acs5"

    params = {
        'get': 'NAME,B19013_001E',  # Tract name and median household income
        'for': 'tract:*',
        'in': 'state:53 county:029',  # Washington (53), Island County (029)
    }

    try:
        response = requests.get(base_url, params=params, timeout=30)
        response.raise_for_status()
        data = response.json()

        # Convert to DataFrame
        df = pd.DataFrame(data[1:], columns=data[0])

        # Clean and process data
        df['tract_name'] = df['NAME']
        df['median_income'] = pd.to_numeric(df['B19013_001E'], errors='coerce')
        df['tract_fips'] = df['tract']
        df['county_fips'] = df['county']
        df['state_fips'] = df['state']
        df['full_fips'] = df['state'] + df['county'] + df['tract']

        # Remove rows with missing data or invalid values (Census uses negative values for suppressed data)
        df = df.dropna(subset=['median_income'])
        df = df[df['median_income'] > 0]  # Filter out negative/invalid values

        print(f"Successfully fetched data for {len(df)} census tracts in Island County")
        return df[['tract_name', 'median_income', 'full_fips', 'tract_fips', 'county_fips', 'state_fips']]

    except requests.exceptions.RequestException as e:
        print(f"Error fetching data: {e}")
        return None

def format_currency(value):
    """Format value as currency."""
    return f"${value:,.0f}"

def create_map(df, output_dir):
    """Create a choropleth map of median household income by census tract."""
    print("Creating map...")

    try:
        # Load US census tracts shapefile from Census Bureau
        print("Loading census tract boundaries...")
        tracts_url = "https://www2.census.gov/geo/tiger/GENZ2022/shp/cb_2022_53_tract_500k.zip"
        tracts = gpd.read_file(tracts_url)

        # Filter for Island County (county FIPS 029)
        island_tracts = tracts[tracts['COUNTYFP'] == '029'].copy()
        island_tracts['full_fips'] = island_tracts['STATEFP'] + island_tracts['COUNTYFP'] + island_tracts['TRACTCE']

        # Merge with income data
        island_map = island_tracts.merge(df, on='full_fips', how='left')

        # Create figure with larger size for detail
        fig, ax = plt.subplots(1, 1, figsize=(16, 12))

        # Plot the choropleth
        island_map.plot(
            column='median_income',
            cmap='RdYlGn',
            linewidth=1.0,
            ax=ax,
            edgecolor='0.3',
            legend=True,
            legend_kwds={
                'label': "Median Household Income (dollars)",
                'orientation': "horizontal",
                'pad': 0.05,
                'shrink': 0.8,
                'format': '${x:,.0f}'
            }
        )

        # Add tract labels with income
        island_map['coords'] = island_map['geometry'].apply(
            lambda x: x.representative_point().coords[:][0]
        )

        for idx, row in island_map.iterrows():
            if pd.notna(row['median_income']):
                # Create shorter label - just tract number and income
                tract_num = row['TRACTCE']
                ax.annotate(
                    text=f"Tract {tract_num}\n{format_currency(row['median_income'])}",
                    xy=row['coords'],
                    horizontalalignment='center',
                    fontsize=8,
                    weight='bold',
                    bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8, edgecolor='gray', linewidth=0.5)
                )

        ax.axis('off')
        ax.set_title(
            'Median Household Income by Census Tract\nIsland County, Washington',
            fontsize=18,
            weight='bold',
            pad=20
        )

        # Add source and timestamp
        plt.figtext(
            0.5, 0.02,
            'Data Source: U.S. Census Bureau, American Community Survey 5-Year Estimates (2022)\n'
            f'Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}',
            ha='center',
            fontsize=9,
            style='italic'
        )

        plt.tight_layout()

        # Save the map
        map_file = os.path.join(output_dir, 'island_county_income_tracts_map.png')
        plt.savefig(map_file, dpi=300, bbox_inches='tight')
        print(f"Map saved to: {map_file}")

        plt.close()

        return island_map

    except Exception as e:
        print(f"Error creating map: {e}")
        import traceback
        traceback.print_exc()
        return None

def save_data(df, output_dir):
    """Save the data to CSV and JSON files."""
    print("Saving data files...")

    # Sort by income descending
    df_sorted = df.sort_values('median_income', ascending=False)

    # Save to CSV
    csv_file = os.path.join(output_dir, 'island_county_income_tracts_data.csv')
    df_sorted.to_csv(csv_file, index=False)
    print(f"CSV saved to: {csv_file}")

    # Save to JSON
    json_file = os.path.join(output_dir, 'island_county_income_tracts_data.json')
    df_sorted.to_json(json_file, orient='records', indent=2)
    print(f"JSON saved to: {json_file}")

    # Print summary statistics
    print("\n" + "="*70)
    print("ISLAND COUNTY CENSUS TRACT INCOME SUMMARY")
    print("="*70)
    print(f"Total Census Tracts: {len(df_sorted)}")
    print(f"Mean: {format_currency(df_sorted['median_income'].mean())}")
    print(f"Median: {format_currency(df_sorted['median_income'].median())}")
    print(f"Min: {format_currency(df_sorted['median_income'].min())}")
    print(f"Max: {format_currency(df_sorted['median_income'].max())}")
    print(f"Standard Deviation: {format_currency(df_sorted['median_income'].std())}")

    print("\n" + "="*70)
    print("TOP 5 CENSUS TRACTS BY MEDIAN INCOME")
    print("="*70)
    for idx, row in df_sorted.head(5).iterrows():
        print(f"Tract {row['tract_fips']:6s}: {format_currency(row['median_income']):>12s}")

    print("\n" + "="*70)
    print("BOTTOM 5 CENSUS TRACTS BY MEDIAN INCOME")
    print("="*70)
    for idx, row in df_sorted.tail(5).iterrows():
        print(f"Tract {row['tract_fips']:6s}: {format_currency(row['median_income']):>12s}")

    print("\n" + "="*70)
    print("ALL CENSUS TRACTS (sorted by income)")
    print("="*70)
    for idx, row in df_sorted.iterrows():
        print(f"Tract {row['tract_fips']:6s}: {format_currency(row['median_income']):>12s}")

def main():
    """Main execution function."""
    # Create output directory
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = f"island_county_income_tracts_{timestamp}"
    os.makedirs(output_dir, exist_ok=True)
    print(f"Output directory: {output_dir}")

    # Fetch data
    df = fetch_income_data()
    if df is None:
        print("Failed to fetch data. Exiting.")
        return

    # Create map
    map_data = create_map(df, output_dir)
    if map_data is None:
        print("Failed to create map. Exiting.")
        return

    # Save data files
    save_data(df, output_dir)

    print(f"\n✓ All outputs saved to: {output_dir}/")

if __name__ == "__main__":
    main()
