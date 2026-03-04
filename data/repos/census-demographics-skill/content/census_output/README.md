# Census Output - Example Scripts

This directory contains example scripts demonstrating various census data visualizations.

## Example Scripts

### 1. Oregon Education Map (`oregon_education_map.py`)
Creates a map showing educational attainment (% with bachelor's degree or higher) by county in Oregon.

**Usage:**
```bash
python oregon_education_map.py
```

**Output:**
- High-resolution choropleth map
- CSV data file with all counties
- JSON data file
- Summary statistics

### 2. Washington Education Map (`washington_education_map.py`)
Similar to Oregon but for Washington state counties.

**Usage:**
```bash
python washington_education_map.py
```

### 3. Washington Income Map (`washington_income_map.py`)
Shows median household income by county in Washington state.

**Usage:**
```bash
python washington_income_map.py
```

### 4. Island County Tract-Level Income Map (`island_county_income_tracts_map.py`)
Detailed census tract-level analysis of median household income for Island County, WA.

**Usage:**
```bash
python island_county_income_tracts_map.py
```

**Features:**
- Census tract-level granularity
- 23 tracts in Island County (Whidbey and Camano Islands)
- Filters out invalid/suppressed data

## Output Files

Each script creates a timestamped directory containing:

1. **Map (PNG)** - High-resolution visualization (300 DPI)
2. **Data (CSV)** - All data in tabular format
3. **Data (JSON)** - Machine-readable data format

## Generated Directories

When you run these scripts, they create output directories with timestamps like:
- `oregon_education_YYYYMMDD_HHMMSS/`
- `washington_income_YYYYMMDD_HHMMSS/`
- `island_county_income_tracts_YYYYMMDD_HHMMSS/`

**Note:** These directories are gitignored to keep the repository clean. Only the source scripts are tracked.

## Data Sources

All scripts use:
- **US Census Bureau API** - For demographic data
- **TIGER/Line Shapefiles** - For geographic boundaries
- **ACS 5-Year Estimates (2022)** - Most recent available data

## Requirements

These scripts require the dependencies listed in the main `requirements.txt`:
- requests
- pandas
- geopandas
- matplotlib

## Customization

You can easily modify these scripts to:
- Change the state or county
- Use different census variables
- Adjust color schemes
- Modify map styling
- Export different file formats

## Notes

- No Census API key required for these specific examples (they use public endpoints)
- Internet connection required to download geographic boundaries
- Processing time varies by geographic level (tracts are slower than counties)
