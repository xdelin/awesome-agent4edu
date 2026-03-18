#!/usr/bin/env python3
"""
Create a choropleth map showing median household income by county in Washington
using US Census ACS 5-Year data.
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
    Fetch median household income data for Washington counties from Census API.
    Uses ACS 5-Year estimates (Table S1903 - Median Income).
    """
    print("Fetching median household income data from US Census API...")

    # ACS 5-Year Subject Table S1903
    # S1903_C03_001E: Median income (dollars) -- Households
    base_url = "https://api.census.gov/data/2022/acs/acs5/subject"

    params = {
        'get': 'NAME,S1903_C03_001E',  # County name and median household income
        'for': 'county:*',
        'in': 'state:53',  # Washington FIPS code
    }

    try:
        response = requests.get(base_url, params=params, timeout=30)
        response.raise_for_status()
        data = response.json()

        # Convert to DataFrame
        df = pd.DataFrame(data[1:], columns=data[0])

        # Clean and process data
        df['county_name'] = df['NAME'].str.replace(' County, Washington', '')
        df['median_income'] = pd.to_numeric(df['S1903_C03_001E'], errors='coerce')
        df['county_fips'] = df['county']
        df['state_fips'] = df['state']
        df['full_fips'] = df['state'] + df['county']

        # Remove rows with missing data
        df = df.dropna(subset=['median_income'])

        print(f"Successfully fetched data for {len(df)} Washington counties")
        return df[['county_name', 'median_income', 'full_fips', 'county_fips', 'state_fips']]

    except requests.exceptions.RequestException as e:
        print(f"Error fetching data: {e}")
        return None

def format_currency(value):
    """Format value as currency."""
    return f"${value:,.0f}"

def create_map(df, output_dir):
    """Create a choropleth map of median household income by county."""
    print("Creating map...")

    try:
        # Load US counties shapefile from Census Bureau
        print("Loading county boundaries...")
        counties_url = "https://www2.census.gov/geo/tiger/GENZ2022/shp/cb_2022_us_county_500k.zip"
        counties = gpd.read_file(counties_url)

        # Filter for Washington (FIPS code 53)
        washington_counties = counties[counties['STATEFP'] == '53'].copy()
        washington_counties['full_fips'] = washington_counties['STATEFP'] + washington_counties['COUNTYFP']

        # Merge with income data
        washington_map = washington_counties.merge(df, on='full_fips', how='left')

        # Create figure
        fig, ax = plt.subplots(1, 1, figsize=(16, 10))

        # Plot the choropleth
        washington_map.plot(
            column='median_income',
            cmap='RdYlGn',
            linewidth=0.8,
            ax=ax,
            edgecolor='0.8',
            legend=True,
            legend_kwds={
                'label': "Median Household Income (dollars)",
                'orientation': "horizontal",
                'pad': 0.05,
                'shrink': 0.8,
                'format': '${x:,.0f}'
            }
        )

        # Add county labels
        washington_map['coords'] = washington_map['geometry'].apply(
            lambda x: x.representative_point().coords[:][0]
        )

        for idx, row in washington_map.iterrows():
            if pd.notna(row['median_income']):
                ax.annotate(
                    text=f"{row['county_name']}\n{format_currency(row['median_income'])}",
                    xy=row['coords'],
                    horizontalalignment='center',
                    fontsize=7,
                    weight='bold',
                    bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.7, edgecolor='none')
                )

        ax.axis('off')
        ax.set_title(
            'Median Household Income in Washington Counties',
            fontsize=16,
            weight='bold',
            pad=20
        )

        # Add source and timestamp
        plt.figtext(
            0.5, 0.02,
            'Data Source: U.S. Census Bureau, American Community Survey 5-Year Estimates (2022)\n'
            f'Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}',
            ha='center',
            fontsize=8,
            style='italic'
        )

        plt.tight_layout()

        # Save the map
        map_file = os.path.join(output_dir, 'washington_income_map.png')
        plt.savefig(map_file, dpi=300, bbox_inches='tight')
        print(f"Map saved to: {map_file}")

        plt.close()

        return washington_map

    except Exception as e:
        print(f"Error creating map: {e}")
        return None

def save_data(df, output_dir):
    """Save the data to CSV and JSON files."""
    print("Saving data files...")

    # Sort by income descending
    df_sorted = df.sort_values('median_income', ascending=False)

    # Save to CSV
    csv_file = os.path.join(output_dir, 'washington_income_data.csv')
    df_sorted.to_csv(csv_file, index=False)
    print(f"CSV saved to: {csv_file}")

    # Save to JSON
    json_file = os.path.join(output_dir, 'washington_income_data.json')
    df_sorted.to_json(json_file, orient='records', indent=2)
    print(f"JSON saved to: {json_file}")

    # Print summary statistics
    print("\n" + "="*60)
    print("SUMMARY STATISTICS")
    print("="*60)
    print(f"Mean: {format_currency(df_sorted['median_income'].mean())}")
    print(f"Median: {format_currency(df_sorted['median_income'].median())}")
    print(f"Min: {format_currency(df_sorted['median_income'].min())} ({df_sorted.iloc[-1]['county_name']})")
    print(f"Max: {format_currency(df_sorted['median_income'].max())} ({df_sorted.iloc[0]['county_name']})")
    print("\n" + "="*60)
    print("TOP 5 COUNTIES BY MEDIAN INCOME")
    print("="*60)
    for idx, row in df_sorted.head(5).iterrows():
        print(f"{row['county_name']:20s}: {format_currency(row['median_income'])}")
    print("\n" + "="*60)
    print("BOTTOM 5 COUNTIES BY MEDIAN INCOME")
    print("="*60)
    for idx, row in df_sorted.tail(5).iterrows():
        print(f"{row['county_name']:20s}: {format_currency(row['median_income'])}")

def main():
    """Main execution function."""
    # Create output directory
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = f"washington_income_{timestamp}"
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
