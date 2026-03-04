#!/usr/bin/env python3
"""
Census Data Visualization Module

This module provides functions to create choropleth maps from Census data.
"""

import os
from datetime import datetime
from typing import Optional, Dict, List, Tuple
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Patch
import contextily as ctx
import numpy as np


class CensusMapVisualizer:
    """Creates demographic visualizations from Census data."""

    # Color schemes for different data types
    COLOR_SCHEMES = {
        'sequential': ['YlOrRd', 'YlGnBu', 'PuBu', 'Purples', 'Blues', 'Greens'],
        'diverging': ['RdBu', 'RdYlGn', 'PiYG', 'BrBG'],
        'categorical': ['Set3', 'Pastel1', 'Pastel2', 'Accent']
    }

    def __init__(self, output_dir: str = './census_output'):
        """
        Initialize the visualizer.

        Args:
            output_dir: Directory to save output files
        """
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)

    def load_geometries(
        self,
        level: str = 'county',
        state: Optional[str] = None,
        year: int = 2021
    ) -> gpd.GeoDataFrame:
        """
        Load geographic boundaries from Census TIGER/Line.

        Args:
            level: Geographic level ('state', 'county', 'tract')
            state: State FIPS code (required for tract level)
            year: Vintage year

        Returns:
            GeoDataFrame with geometries
        """
        base_url = f'https://www2.census.gov/geo/tiger/TIGER{year}'

        if level == 'state':
            url = f'{base_url}/STATE/tl_{year}_us_state.zip'
        elif level == 'county':
            url = f'{base_url}/COUNTY/tl_{year}_us_county.zip'
        elif level == 'tract':
            if not state:
                raise ValueError("State FIPS required for tract-level geometries")
            url = f'{base_url}/TRACT/tl_{year}_{state}_tract.zip'
        else:
            raise ValueError(f"Unknown level: {level}")

        print(f"Loading {level}-level geometries from Census TIGER/Line...")
        gdf = gpd.read_file(url)

        # Ensure consistent GEOID column
        if level == 'state':
            gdf['GEOID'] = gdf['STATEFP']
        elif level == 'county':
            gdf['GEOID'] = gdf['STATEFP'] + gdf['COUNTYFP']
        elif level == 'tract':
            gdf['GEOID'] = gdf['STATEFP'] + gdf['COUNTYFP'] + gdf['TRACTCE']

        print(f"Loaded {len(gdf)} {level} geometries")
        return gdf

    def merge_data(
        self,
        geometries: gpd.GeoDataFrame,
        data: pd.DataFrame,
        on: str = 'GEOID'
    ) -> gpd.GeoDataFrame:
        """
        Merge Census data with geometries.

        Args:
            geometries: GeoDataFrame with geometries
            data: DataFrame with Census data
            on: Column to join on

        Returns:
            Merged GeoDataFrame
        """
        merged = geometries.merge(data, on=on, how='left')
        print(f"Merged data: {len(merged)} features")
        return merged

    def create_choropleth(
        self,
        gdf: gpd.GeoDataFrame,
        column: str,
        title: str,
        cmap: str = 'YlOrRd',
        figsize: Tuple[int, int] = (15, 10),
        legend_label: Optional[str] = None,
        scheme: Optional[str] = None,
        k: int = 5,
        add_basemap: bool = False,
        vmin: Optional[float] = None,
        vmax: Optional[float] = None
    ) -> Tuple[plt.Figure, plt.Axes]:
        """
        Create a choropleth map.

        Args:
            gdf: GeoDataFrame with data and geometries
            column: Column to visualize
            title: Map title
            cmap: Matplotlib colormap name
            figsize: Figure size (width, height)
            legend_label: Label for legend
            scheme: Classification scheme ('quantiles', 'equal_interval', 'natural_breaks')
            k: Number of classes for classification
            add_basemap: Whether to add a basemap
            vmin: Minimum value for color scale
            vmax: Maximum value for color scale

        Returns:
            Tuple of (figure, axes)
        """
        fig, ax = plt.subplots(1, 1, figsize=figsize)

        # Ensure column is numeric
        gdf[column] = pd.to_numeric(gdf[column], errors='coerce')

        # Remove null geometries
        gdf = gdf[gdf.geometry.notna()]

        # Project to Web Mercator for basemap
        if add_basemap:
            gdf = gdf.to_crs(epsg=3857)

        # Create the choropleth
        if scheme:
            gdf.plot(
                column=column,
                ax=ax,
                legend=True,
                cmap=cmap,
                edgecolor='0.8',
                linewidth=0.5,
                scheme=scheme,
                k=k,
                legend_kwds={'loc': 'lower left', 'frameon': False}
            )
        else:
            gdf.plot(
                column=column,
                ax=ax,
                legend=True,
                cmap=cmap,
                edgecolor='0.8',
                linewidth=0.5,
                vmin=vmin,
                vmax=vmax,
                legend_kwds={'label': legend_label or column, 'shrink': 0.5}
            )

        # Add basemap if requested
        if add_basemap:
            try:
                ctx.add_basemap(
                    ax,
                    source=ctx.providers.CartoDB.Positron,
                    zoom='auto'
                )
            except Exception as e:
                print(f"Warning: Could not add basemap: {e}")

        ax.set_title(title, fontsize=16, fontweight='bold', pad=20)
        ax.axis('off')

        plt.tight_layout()
        return fig, ax

    def create_multi_map(
        self,
        gdf: gpd.GeoDataFrame,
        columns: List[str],
        titles: List[str],
        cmap: str = 'YlOrRd',
        figsize: Tuple[int, int] = (20, 12),
        suptitle: Optional[str] = None
    ) -> Tuple[plt.Figure, List[plt.Axes]]:
        """
        Create multiple maps in a grid.

        Args:
            gdf: GeoDataFrame with data and geometries
            columns: List of columns to visualize
            titles: List of subplot titles
            cmap: Matplotlib colormap name
            figsize: Figure size (width, height)
            suptitle: Overall title

        Returns:
            Tuple of (figure, list of axes)
        """
        n_maps = len(columns)
        n_cols = min(2, n_maps)
        n_rows = (n_maps + n_cols - 1) // n_cols

        fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize)
        if n_maps == 1:
            axes = [axes]
        else:
            axes = axes.flatten()

        for idx, (col, title) in enumerate(zip(columns, titles)):
            ax = axes[idx]

            # Ensure column is numeric
            gdf[col] = pd.to_numeric(gdf[col], errors='coerce')

            gdf.plot(
                column=col,
                ax=ax,
                legend=True,
                cmap=cmap,
                edgecolor='0.8',
                linewidth=0.5,
                missing_kwds={'color': 'lightgrey'},
                legend_kwds={'shrink': 0.5}
            )

            ax.set_title(title, fontsize=12, fontweight='bold')
            ax.axis('off')

        # Hide empty subplots
        for idx in range(n_maps, len(axes)):
            axes[idx].axis('off')

        if suptitle:
            fig.suptitle(suptitle, fontsize=16, fontweight='bold', y=0.98)

        plt.tight_layout()
        return fig, axes

    def save_map(
        self,
        fig: plt.Figure,
        filename: str,
        dpi: int = 300
    ) -> str:
        """
        Save map to file.

        Args:
            fig: Matplotlib figure
            filename: Output filename (without extension)
            dpi: Resolution in dots per inch

        Returns:
            Full path to saved file
        """
        filepath = os.path.join(self.output_dir, f'{filename}.png')
        fig.savefig(filepath, dpi=dpi, bbox_inches='tight')
        print(f"Map saved to: {filepath}")
        return filepath

    def create_summary_stats(
        self,
        df: pd.DataFrame,
        columns: List[str],
        filename: str = 'summary_statistics.txt'
    ) -> str:
        """
        Create summary statistics file.

        Args:
            df: DataFrame with data
            columns: Columns to summarize
            filename: Output filename

        Returns:
            Full path to saved file
        """
        filepath = os.path.join(self.output_dir, filename)

        with open(filepath, 'w') as f:
            f.write("CENSUS DATA SUMMARY STATISTICS\n")
            f.write("=" * 60 + "\n")
            f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Total records: {len(df)}\n\n")

            for col in columns:
                if col not in df.columns:
                    continue

                # Convert to numeric
                data = pd.to_numeric(df[col], errors='coerce')
                data = data.dropna()

                if len(data) == 0:
                    continue

                f.write(f"\n{col}\n")
                f.write("-" * 60 + "\n")
                f.write(f"  Count:   {len(data):>12,}\n")
                f.write(f"  Mean:    {data.mean():>12,.2f}\n")
                f.write(f"  Median:  {data.median():>12,.2f}\n")
                f.write(f"  Std Dev: {data.std():>12,.2f}\n")
                f.write(f"  Min:     {data.min():>12,.2f}\n")
                f.write(f"  Max:     {data.max():>12,.2f}\n")

                # Percentiles
                percentiles = data.quantile([0.25, 0.5, 0.75])
                f.write(f"  25th %:  {percentiles[0.25]:>12,.2f}\n")
                f.write(f"  50th %:  {percentiles[0.50]:>12,.2f}\n")
                f.write(f"  75th %:  {percentiles[0.75]:>12,.2f}\n")

        print(f"Summary statistics saved to: {filepath}")
        return filepath

    def export_data(
        self,
        df: pd.DataFrame,
        filename: str = 'census_data.csv'
    ) -> str:
        """
        Export data to CSV.

        Args:
            df: DataFrame to export
            filename: Output filename

        Returns:
            Full path to saved file
        """
        filepath = os.path.join(self.output_dir, filename)

        # Drop geometry column if present (can't export to CSV)
        if 'geometry' in df.columns:
            df = df.drop(columns=['geometry'])

        df.to_csv(filepath, index=False)
        print(f"Data exported to: {filepath}")
        return filepath


def main():
    """Example usage of CensusMapVisualizer."""

    # Create visualizer
    viz = CensusMapVisualizer(output_dir='./example_output')

    # Load county geometries
    counties = viz.load_geometries(level='county', year=2021)

    # Create sample data (in real use, this would come from census_fetch.py)
    sample_data = pd.DataFrame({
        'GEOID': counties['GEOID'].head(100),
        'population': np.random.randint(10000, 500000, 100),
        'median_income': np.random.randint(30000, 100000, 100)
    })

    # Merge data
    merged = viz.merge_data(counties, sample_data)

    # Create choropleth map
    fig, ax = viz.create_choropleth(
        merged,
        column='population',
        title='Sample Population Map',
        cmap='YlOrRd',
        legend_label='Population'
    )

    # Save map
    viz.save_map(fig, 'sample_population_map')

    # Create summary statistics
    viz.create_summary_stats(merged, ['population', 'median_income'])

    # Export data
    viz.export_data(sample_data)

    print("\nExample visualizations complete!")


if __name__ == '__main__':
    main()
