# Census Demographics Visualization Skill

A comprehensive skill for downloading US Census data and creating demographic maps showing population density, educational attainment, age distribution, race/ethnicity, income levels, and more.

## Features

- **Comprehensive Demographic Data**: Access population, age, race, education, income, housing, and employment statistics
- **Multiple Geographic Levels**: State, county, and census tract level data
- **Beautiful Visualizations**: Create choropleth maps with customizable color schemes
- **Automatic Data Integration**: Seamlessly combines Census API data with geographic boundaries
- **Summary Statistics**: Generate detailed statistical summaries of your data
- **Export Capabilities**: Save maps as high-resolution images and data as CSV files

## Setup

### 1. Install Dependencies

```bash
pip install -r requirements.txt
```

### 2. Get a Census API Key

1. Visit https://api.census.gov/data/key_signup.html
2. Sign up for a free API key
3. Set your API key as an environment variable:

```bash
export CENSUS_API_KEY='your_api_key_here'
```

Or add it to your `.bashrc` or `.zshrc` for persistence:

```bash
echo 'export CENSUS_API_KEY="your_api_key_here"' >> ~/.bashrc
source ~/.bashrc
```

### 3. Invoke the Skill

In Claude Code, use the skill command:

```
/skill census-demographics
```

Or simply type `census-demographics` in your conversation.

## Usage Examples

### Example 1: Population Density Map

```
Show me a population density map of California counties
```

The skill will:
1. Fetch population data for California counties
2. Download county boundary geometries
3. Calculate population density
4. Create and save a choropleth map

### Example 2: Educational Attainment Comparison

```
Compare bachelor's degree attainment rates across all US states
```

The skill will:
1. Fetch education data from ACS 5-year estimates
2. Calculate percentage with bachelor's degrees
3. Create state-level choropleth map
4. Generate summary statistics

### Example 3: Age Distribution

```
Map median age by county in Texas
```

The skill will:
1. Fetch age data for Texas counties
2. Create visualization showing median age
3. Provide statistical summary

### Example 4: Multiple Demographics

```
Show me maps of population, median income, and education levels for counties in the Northeast
```

The skill will:
1. Fetch multiple demographic variables
2. Create a multi-panel map showing all three metrics
3. Export data to CSV for further analysis

### Example 5: Detailed Tract-Level Analysis

```
Create a detailed map of racial diversity in Cook County, Illinois at the tract level
```

The skill will:
1. Fetch race/ethnicity data at census tract level
2. Create high-resolution tract-level map
3. Calculate diversity indices

## Available Demographic Categories

### Population
- Total population
- Population density

### Age Distribution
- Median age
- Under 18 population
- Over 65 population
- Age cohorts

### Race and Ethnicity
- White alone
- Black or African American
- American Indian and Alaska Native
- Asian
- Native Hawaiian and Pacific Islander
- Hispanic or Latino
- Two or more races

### Educational Attainment (25+ years old)
- High school graduate
- Some college
- Associate's degree
- Bachelor's degree
- Master's degree
- Professional degree
- Doctorate degree

### Income
- Median household income
- Per capita income
- Mean household income

### Housing
- Total housing units
- Median home value
- Median rent
- Owner-occupied vs. renter-occupied

### Employment
- Labor force participation
- Employment rate
- Unemployment rate

## Census Variable Reference

Common Census variable codes used by this skill:

| Category | Variable | Code |
|----------|----------|------|
| Population | Total | B01003_001E |
| Age | Median | B01002_001E |
| Race | White | B02001_002E |
| Race | Black | B02001_003E |
| Race | Asian | B02001_005E |
| Ethnicity | Hispanic | B03003_003E |
| Education | Bachelor's | B15003_022E |
| Income | Median HH | B19013_001E |
| Housing | Median Value | B25077_001E |

For a complete list of variables, visit: https://api.census.gov/data/2021/acs/acs5/variables.html

## Geographic Levels

### State Level
- 50 states + DC + Puerto Rico
- Best for national comparisons
- Fast to process

### County Level
- ~3,200 counties nationwide
- Good balance of detail and performance
- Most commonly used level

### Census Tract Level
- ~84,000 tracts nationwide
- Neighborhood-level detail
- Requires state/county specification
- Slower to process for large areas

## Output Files

When you run the skill, it creates a timestamped output directory with:

1. **demographic_map.png** - High-resolution map (300 DPI)
2. **summary_statistics.txt** - Descriptive statistics
3. **census_data.csv** - Raw data table
4. **README.txt** - Explanation of the visualization

## Tips for Best Results

1. **Start Broad, Then Narrow**: Begin with state or county level, then drill down to tracts if needed

2. **Specify Time Period**: ACS 5-year estimates are most detailed but lag by 2-3 years

3. **Use Appropriate Color Schemes**:
   - Sequential (YlOrRd, Blues) for single-variable continuous data
   - Diverging (RdBu) for data with meaningful midpoint
   - Categorical (Set3) for discrete categories

4. **Consider Data Privacy**: Block-level data may be suppressed for small populations

5. **Combine Multiple Indicators**: Create multi-panel maps to show relationships

## Troubleshooting

### "API key not found"
- Ensure `CENSUS_API_KEY` environment variable is set
- Verify the key is valid at https://api.census.gov/data.html

### "Could not load geometries"
- Check internet connection (downloads shapefiles from Census Bureau)
- Verify the year is valid (typically 2010-present)

### "No data returned"
- Check FIPS codes are correct
- Verify variable code exists for the specified year/dataset
- Some variables only available at certain geographic levels

### Maps look cluttered
- Try a higher-level geography (counties instead of tracts)
- Use classification schemes (quantiles, natural breaks)
- Filter to a specific region instead of entire US

## Advanced Usage

### Custom Variable Codes

You can specify custom Census variable codes:

```
Fetch variable B25064_001E (median rent) for all counties and create a map
```

### Multiple Years Comparison

```
Compare median household income between 2015 and 2021 for California counties
```

### Calculated Metrics

```
Calculate and map the percentage of population with graduate degrees for each state
```

### Custom Color Schemes

```
Create a population map using the 'viridis' color scheme
```

## Python API Usage

You can also use the modules directly in Python scripts:

```python
from census_fetch import CensusDataFetcher
from map_visualize import CensusMapVisualizer

# Fetch data
fetcher = CensusDataFetcher(api_key='your_key')
data = fetcher.fetch_county_data(
    variables=['B01003_001E'],  # Total population
    state_fips='06',  # California
    year=2021
)

# Create visualization
viz = CensusMapVisualizer(output_dir='./my_maps')
counties = viz.load_geometries(level='county')
merged = viz.merge_data(counties, data)

fig, ax = viz.create_choropleth(
    merged,
    column='B01003_001E',
    title='California County Population',
    cmap='YlOrRd'
)

viz.save_map(fig, 'ca_population')
```

## Data Sources

All data comes from the US Census Bureau:

- **American Community Survey (ACS)**: Most detailed demographic data
  - 5-Year Estimates: Most reliable, covers all geographies
  - 1-Year Estimates: Most current, only for areas with 65k+ population

- **Decennial Census**: Complete population count every 10 years (2010, 2020)

- **TIGER/Line Shapefiles**: Geographic boundaries updated annually

## Resources

- Census API Documentation: https://www.census.gov/data/developers/data-sets.html
- Variable Search: https://api.census.gov/data.html
- Geography Reference: https://www.census.gov/programs-surveys/geography.html
- ACS Handbook: https://www.census.gov/programs-surveys/acs/guidance.html

## License

This skill uses publicly available US Census data. Census Bureau data are free from copyright restrictions.

## Support

For issues with this skill:
1. Check the troubleshooting section above
2. Verify your Census API key is valid
3. Consult the Census Bureau's data documentation
4. Ask Claude Code for help with specific error messages

## Version

Version: 1.0.0
Last Updated: 2025
Compatible with: Claude Code
