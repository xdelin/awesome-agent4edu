#!/usr/bin/env python3
"""
Create a population density map for all US states.
"""

import os
import sys
from datetime import datetime
from census_fetch import CensusDataFetcher
from map_visualize import CensusMapVisualizer

def main():
    print("=" * 60)
    print("US STATE POPULATION DENSITY MAP")
    print("=" * 60)
    print()

    # Get API key
    api_key = os.environ.get('CENSUS_API_KEY')
    if not api_key:
        print("Error: CENSUS_API_KEY environment variable not set")
        sys.exit(1)

    # Initialize fetcher and visualizer
    print("Initializing Census data fetcher...")
    fetcher = CensusDataFetcher(api_key)

    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    output_dir = f'./census_output/state_density_{timestamp}'
    viz = CensusMapVisualizer(output_dir=output_dir)
    print(f"Output directory: {output_dir}")
    print()

    # Fetch state population data
    print("Fetching population data for all US states...")
    pop_vars = ['B01003_001E']  # Total population
    state_data = fetcher.fetch_state_data(pop_vars, year=2021)
    print(f"✓ Fetched data for {len(state_data)} states/territories")
    print()

    # Load state geometries
    print("Loading state geographic boundaries...")
    states_geo = viz.load_geometries(level='state', year=2021)
    print(f"✓ Loaded geometries for {len(states_geo)} states")
    print()

    # Merge data with geometries
    print("Merging data with geometries...")
    merged = viz.merge_data(states_geo, state_data)

    # Calculate area in square miles
    print("Calculating population density...")
    # Project to Albers Equal Area for accurate area calculation
    merged_albers = merged.to_crs('EPSG:5070')  # NAD83 / Conus Albers
    merged['area_sq_mi'] = merged_albers.geometry.area / 2589988.110336  # Convert sq meters to sq miles

    # Calculate density
    merged['population'] = merged['B01003_001E'].astype(float)
    merged['density'] = merged['population'] / merged['area_sq_mi']

    # Filter to 50 states + DC (exclude territories)
    state_fips = [
        '01', '04', '05', '06', '08', '09', '10', '11', '12', '13',  # AL-GA
        '15', '16', '17', '18', '19', '20', '21', '22', '23', '24',  # HI-MD
        '25', '26', '27', '28', '29', '30', '31', '32', '33', '34',  # MA-NJ
        '35', '36', '37', '38', '39', '40', '41', '42', '44', '45',  # NM-SC
        '46', '47', '48', '49', '50', '51', '53', '54', '55', '56'   # SD-WY
    ]
    merged = merged[merged['GEOID'].isin(state_fips)]

    print(f"✓ Calculated density for {len(merged)} states")
    print()

    # Print some statistics
    print("Top 10 Most Dense States:")
    top_dense = merged.nlargest(10, 'density')
    for idx, row in top_dense.iterrows():
        name = row['NAME'] if 'NAME' in row and row['NAME'] else 'Unknown'
        print(f"  {name:<20} {row['density']:>10,.1f} people/sq mi")
    print()

    print("Top 10 Least Dense States:")
    low_dense = merged.nsmallest(10, 'density')
    for idx, row in low_dense.iterrows():
        name = row['NAME'] if 'NAME' in row and row['NAME'] else 'Unknown'
        print(f"  {name:<20} {row['density']:>10,.1f} people/sq mi")
    print()

    # Create the choropleth map
    print("Creating population density map...")
    fig, ax = viz.create_choropleth(
        merged,
        column='density',
        title='US State Population Density (2021)\nPeople per Square Mile',
        cmap='YlOrRd',
        figsize=(20, 12),
        legend_label='Population Density (people/sq mi)',
        scheme='quantiles',
        k=7
    )

    # Save the map
    print("Saving map...")
    map_file = viz.save_map(fig, 'us_state_population_density')
    print()

    # Create summary statistics
    print("Generating summary statistics...")
    stats_file = viz.create_summary_stats(
        merged,
        ['population', 'area_sq_mi', 'density'],
        filename='density_statistics.txt'
    )
    print()

    # Export data
    print("Exporting data to CSV...")
    # Select only the columns that exist
    cols_to_export = ['GEOID', 'population', 'area_sq_mi', 'density']
    if 'NAME' in merged.columns:
        cols_to_export.insert(0, 'NAME')
    elif 'STUSPS' in merged.columns:
        cols_to_export.insert(0, 'STUSPS')

    export_df = merged[cols_to_export].copy()
    export_df = export_df.sort_values('density', ascending=False)
    if 'NAME' in export_df.columns:
        export_df.rename(columns={'NAME': 'State'}, inplace=True)
    elif 'STUSPS' in export_df.columns:
        export_df.rename(columns={'STUSPS': 'State'}, inplace=True)
    csv_file = viz.export_data(export_df, filename='state_density_data.csv')
    print()

    # Create README
    readme_path = os.path.join(output_dir, 'README.txt')
    with open(readme_path, 'w') as f:
        f.write("US STATE POPULATION DENSITY MAP\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        f.write("DATA SOURCE\n")
        f.write("-" * 60 + "\n")
        f.write("Population: US Census Bureau, American Community Survey\n")
        f.write("            2021 5-Year Estimates (Table B01003)\n")
        f.write("Boundaries: US Census Bureau TIGER/Line Shapefiles, 2021\n\n")
        f.write("METHODOLOGY\n")
        f.write("-" * 60 + "\n")
        f.write("Population density calculated as: Total Population / Land Area\n")
        f.write("Area calculated using NAD83 Albers Equal Area projection\n")
        f.write("Map includes 50 US states + District of Columbia\n")
        f.write("Territories (Puerto Rico, Guam, etc.) excluded\n\n")
        f.write("FILES\n")
        f.write("-" * 60 + "\n")
        f.write("- us_state_population_density.png : Main map visualization\n")
        f.write("- state_density_data.csv          : Complete data table\n")
        f.write("- density_statistics.txt          : Statistical summary\n")
        f.write("- README.txt                      : This file\n\n")
        f.write("INTERPRETATION\n")
        f.write("-" * 60 + "\n")
        f.write("Darker colors indicate higher population density.\n")
        f.write("The map uses a quantile classification with 7 classes,\n")
        f.write("ensuring roughly equal numbers of states in each category.\n\n")

    print(f"✓ README saved to: {readme_path}")
    print()

    print("=" * 60)
    print("COMPLETE!")
    print("=" * 60)
    print()
    print(f"All outputs saved to: {output_dir}")
    print()
    print("Files created:")
    print(f"  • {map_file}")
    print(f"  • {csv_file}")
    print(f"  • {stats_file}")
    print(f"  • {readme_path}")
    print()

if __name__ == '__main__':
    main()
