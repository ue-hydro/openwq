#!/usr/bin/env python3
"""
Prepare Calibration Observations from GRQA Data

This script converts extracted GRQA water quality data to the format required
by the OpenWQ calibration workflow. It performs spatial matching of observation
stations to model river network reaches.

Workflow:
1. Load GRQA extracted stations and observations (from grqa_extract_stations.py)
2. Load model river network shapefile with reach IDs
3. Match stations to nearest reaches (within tolerance)
4. Convert species names using the mapping
5. Convert units if needed
6. Output calibration-ready observation file

Required Input Files:
- GRQA stations shapefile (from grqa_extract_stations.py)
- GRQA observations CSV (from grqa_extract_stations.py)
- Model river network shapefile with reach_id attribute
- Species mapping JSON (from grqa_extract_stations.py or custom)

Output Format:
- datetime: YYYY-MM-DD HH:MM:SS
- reach_id: integer (model reach ID)
- species: string (model species name)
- value: float (concentration)
- units: string (e.g., mg/l)
- source: string (station ID)
- uncertainty: float (optional, estimated or from data)
- quality_flag: string (GOOD, SUSPECT, BAD)

Author: OpenWQ Development Team
Date: 2024
"""

import os
import sys
import json
import argparse
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union
from datetime import datetime
import warnings

warnings.filterwarnings('ignore')

try:
    import pandas as pd
    import geopandas as gpd
    import numpy as np
    from shapely.geometry import Point
    from shapely.ops import nearest_points
except ImportError as e:
    print(f"Error: Missing required package: {e.name}")
    print("Install with: pip install pandas geopandas numpy shapely")
    sys.exit(1)


# Unit conversion factors to mg/L (target unit for OpenWQ)
UNIT_CONVERSIONS = {
    # Already in mg/L
    'mg/l': 1.0,
    'mg/L': 1.0,
    'MG/L': 1.0,
    'ppm': 1.0,

    # From µg/L
    'ug/l': 0.001,
    'ug/L': 0.001,
    'µg/l': 0.001,
    'µg/L': 0.001,
    'ppb': 0.001,

    # From g/L
    'g/l': 1000.0,
    'g/L': 1000.0,

    # Temperature (no conversion, just passthrough)
    '°C': 1.0,
    'C': 1.0,
    'degC': 1.0,

    # pH (no conversion)
    'pH': 1.0,
    '-': 1.0,

    # Percentage (for DO saturation)
    '%': 1.0,
    'percent': 1.0,
}

# Default uncertainty as fraction of value (if not provided)
DEFAULT_UNCERTAINTY_FRACTION = 0.10  # 10%


class StationReachMatcher:
    """
    Matches observation stations to model river reaches.
    """

    def __init__(self,
                 river_network_path: str,
                 reach_id_column: str = 'seg_id',
                 max_distance_m: float = 1000):
        """
        Initialize the matcher.

        Parameters
        ----------
        river_network_path : str
            Path to river network shapefile
        reach_id_column : str
            Column name containing reach IDs
        max_distance_m : float
            Maximum distance for matching (meters)
        """
        self.max_distance_m = max_distance_m
        self.reach_id_column = reach_id_column

        # Load river network
        print(f"Loading river network: {river_network_path}")
        self.river_gdf = gpd.read_file(river_network_path)

        if reach_id_column not in self.river_gdf.columns:
            available = list(self.river_gdf.columns)
            raise ValueError(
                f"Reach ID column '{reach_id_column}' not found. "
                f"Available columns: {available}"
            )

        # Ensure WGS84
        if self.river_gdf.crs is None:
            self.river_gdf = self.river_gdf.set_crs("EPSG:4326")
        elif self.river_gdf.crs.to_epsg() != 4326:
            self.river_gdf = self.river_gdf.to_crs("EPSG:4326")

        # Determine UTM zone for distance calculations
        centroid = self.river_gdf.dissolve().centroid.iloc[0]
        utm_zone = int((centroid.x + 180) / 6) + 1
        hemisphere = 'north' if centroid.y >= 0 else 'south'
        self.utm_epsg = 32600 + utm_zone if hemisphere == 'north' else 32700 + utm_zone

        # Project for accurate distance calculations
        self.river_gdf_utm = self.river_gdf.to_crs(f"EPSG:{self.utm_epsg}")

        print(f"  Loaded {len(self.river_gdf)} reaches")
        print(f"  Reach ID column: {reach_id_column}")
        print(f"  Using UTM zone {utm_zone} for distance calculations")

    def match_stations(self,
                       stations_gdf: gpd.GeoDataFrame,
                       station_id_column: str = 'site_id') -> pd.DataFrame:
        """
        Match stations to nearest reaches.

        Parameters
        ----------
        stations_gdf : gpd.GeoDataFrame
            Station locations
        station_id_column : str
            Column containing station IDs

        Returns
        -------
        pd.DataFrame
            Mapping of station_id to reach_id with distance
        """
        print(f"\nMatching {len(stations_gdf)} stations to reaches...")
        print(f"  Max matching distance: {self.max_distance_m}m")

        # Ensure same CRS
        if stations_gdf.crs is None:
            stations_gdf = stations_gdf.set_crs("EPSG:4326")
        elif stations_gdf.crs.to_epsg() != 4326:
            stations_gdf = stations_gdf.to_crs("EPSG:4326")

        # Project to UTM
        stations_utm = stations_gdf.to_crs(f"EPSG:{self.utm_epsg}")

        matches = []

        for idx, station in stations_utm.iterrows():
            station_id = station[station_id_column]
            station_point = station.geometry

            # Find nearest reach
            min_dist = float('inf')
            nearest_reach_id = None

            for _, reach in self.river_gdf_utm.iterrows():
                dist = station_point.distance(reach.geometry)
                if dist < min_dist:
                    min_dist = dist
                    nearest_reach_id = reach[self.reach_id_column]

            if min_dist <= self.max_distance_m:
                matches.append({
                    'station_id': station_id,
                    'reach_id': int(nearest_reach_id),
                    'distance_m': round(min_dist, 1),
                    'matched': True
                })
            else:
                matches.append({
                    'station_id': station_id,
                    'reach_id': None,
                    'distance_m': round(min_dist, 1),
                    'matched': False
                })

        matches_df = pd.DataFrame(matches)

        n_matched = matches_df['matched'].sum()
        print(f"  Matched: {n_matched}/{len(matches_df)} stations")

        if n_matched < len(matches_df):
            unmatched = matches_df[~matches_df['matched']]
            print(f"  Unmatched stations (distance > {self.max_distance_m}m):")
            for _, row in unmatched.iterrows():
                print(f"    {row['station_id']}: {row['distance_m']}m to nearest reach")

        return matches_df


class CalibrationObservationConverter:
    """
    Converts GRQA observations to OpenWQ calibration format.
    """

    def __init__(self,
                 species_mapping: Dict[str, str],
                 target_units: str = 'mg/l',
                 default_quality_flag: str = 'GOOD'):
        """
        Initialize converter.

        Parameters
        ----------
        species_mapping : dict
            Mapping from GRQA parameters to model species names
        target_units : str
            Target units for output (default: mg/l)
        default_quality_flag : str
            Default quality flag if not provided
        """
        self.species_mapping = species_mapping
        self.target_units = target_units
        self.default_quality_flag = default_quality_flag

    def convert(self,
                observations_df: pd.DataFrame,
                station_reach_mapping: pd.DataFrame,
                grqa_columns: Dict[str, str] = None) -> pd.DataFrame:
        """
        Convert observations to calibration format.

        Parameters
        ----------
        observations_df : pd.DataFrame
            GRQA observations
        station_reach_mapping : pd.DataFrame
            Station to reach mapping
        grqa_columns : dict, optional
            Column name mapping for GRQA data

        Returns
        -------
        pd.DataFrame
            Calibration-ready observations
        """
        # Default column mapping
        if grqa_columns is None:
            grqa_columns = {
                'site_id': 'site_id',
                'obs_date': 'obs_date',
                'obs_time': 'obs_time',
                'obs_value': 'obs_value',
                'unit': 'unit',
                'model_species': 'model_species',
                'grqa_parameter': 'grqa_parameter',
                'outlier': 'obs_iqr_outlier',
            }

        # Auto-detect columns
        grqa_columns = self._detect_columns(observations_df, grqa_columns)

        print(f"\nConverting {len(observations_df)} observations...")

        # Merge with reach mapping
        site_col = grqa_columns['site_id']
        obs_with_reach = observations_df.merge(
            station_reach_mapping[['station_id', 'reach_id', 'matched']],
            left_on=site_col,
            right_on='station_id',
            how='left'
        )

        # Filter to matched only
        obs_matched = obs_with_reach[obs_with_reach['matched'] == True].copy()
        print(f"  Observations with matched reaches: {len(obs_matched)}")

        if len(obs_matched) == 0:
            print("  WARNING: No observations could be matched to reaches!")
            return pd.DataFrame()

        # Build output dataframe
        calibration_obs = []

        date_col = grqa_columns.get('obs_date', 'obs_date')
        time_col = grqa_columns.get('obs_time', 'obs_time')
        value_col = grqa_columns.get('obs_value', 'obs_value')
        unit_col = grqa_columns.get('unit', 'unit')
        species_col = grqa_columns.get('model_species', 'model_species')
        outlier_col = grqa_columns.get('outlier', 'obs_iqr_outlier')

        for _, row in obs_matched.iterrows():
            # Parse datetime
            try:
                if time_col in row and pd.notna(row[time_col]):
                    dt_str = f"{row[date_col]} {row[time_col]}"
                else:
                    dt_str = f"{row[date_col]} 00:00:00"
                dt = pd.to_datetime(dt_str)
            except:
                dt = pd.to_datetime(row[date_col])

            # Get species name
            if species_col in row:
                species = row[species_col]
            else:
                # Try to map from grqa_parameter
                grqa_param = row.get(grqa_columns.get('grqa_parameter', 'grqa_parameter'), '')
                species = self.species_mapping.get(grqa_param, grqa_param)

            # Get and convert value
            value = float(row[value_col])
            unit = str(row.get(unit_col, 'mg/l'))

            # Convert units
            conversion_factor = UNIT_CONVERSIONS.get(unit, 1.0)
            value_converted = value * conversion_factor

            # Estimate uncertainty (10% of value if not provided)
            uncertainty = value_converted * DEFAULT_UNCERTAINTY_FRACTION

            # Quality flag
            quality_flag = self.default_quality_flag
            if outlier_col in row:
                if row[outlier_col] == True or row[outlier_col] == 1:
                    quality_flag = 'SUSPECT'

            calibration_obs.append({
                'datetime': dt.strftime('%Y-%m-%d %H:%M:%S'),
                'reach_id': int(row['reach_id']),
                'species': species,
                'value': round(value_converted, 6),
                'units': self.target_units,
                'source': str(row[site_col]),
                'uncertainty': round(uncertainty, 6),
                'quality_flag': quality_flag
            })

        result_df = pd.DataFrame(calibration_obs)

        # Sort by datetime and reach
        result_df = result_df.sort_values(['datetime', 'reach_id', 'species'])

        print(f"  Output observations: {len(result_df)}")
        print(f"  Species: {result_df['species'].unique().tolist()}")
        print(f"  Reaches: {result_df['reach_id'].nunique()}")
        print(f"  Date range: {result_df['datetime'].min()} to {result_df['datetime'].max()}")

        return result_df

    def _detect_columns(self, df: pd.DataFrame, col_map: Dict) -> Dict:
        """Auto-detect column names."""
        detected = {}
        columns_lower = {c.lower(): c for c in df.columns}

        for key, default in col_map.items():
            if default in df.columns:
                detected[key] = default
            elif default.lower() in columns_lower:
                detected[key] = columns_lower[default.lower()]
            else:
                # Try common alternatives
                alternatives = {
                    'site_id': ['site_id', 'station_id', 'site_no'],
                    'obs_date': ['obs_date', 'date', 'sample_date'],
                    'obs_value': ['obs_value', 'value', 'result', 'concentration'],
                    'unit': ['unit', 'units'],
                }
                for alt in alternatives.get(key, []):
                    if alt in df.columns:
                        detected[key] = alt
                        break
                    elif alt.lower() in columns_lower:
                        detected[key] = columns_lower[alt.lower()]
                        break

        return detected


def prepare_calibration_observations(
    # Input files
    grqa_stations_path: str,
    grqa_observations_path: str,
    river_network_path: str,
    species_mapping_path: str,

    # Output
    output_path: str,

    # Options
    reach_id_column: str = 'seg_id',
    max_match_distance_m: float = 1000,
    target_units: str = 'mg/l',

    # Filtering
    years: Tuple[int, int] = None,
    species_filter: List[str] = None,
    quality_filter: bool = True,
) -> pd.DataFrame:
    """
    Complete pipeline to prepare calibration observations.

    Parameters
    ----------
    grqa_stations_path : str
        Path to GRQA stations shapefile or CSV
    grqa_observations_path : str
        Path to GRQA observations CSV
    river_network_path : str
        Path to model river network shapefile
    species_mapping_path : str
        Path to species mapping JSON
    output_path : str
        Output CSV path
    reach_id_column : str
        Column name for reach IDs in river network
    max_match_distance_m : float
        Maximum distance for station-reach matching
    target_units : str
        Target units for output
    years : tuple, optional
        (start_year, end_year) filter
    species_filter : list, optional
        List of species to include
    quality_filter : bool
        Filter out outliers

    Returns
    -------
    pd.DataFrame
        Calibration-ready observations
    """
    print("=" * 60)
    print("PREPARE CALIBRATION OBSERVATIONS")
    print("=" * 60)

    # Load species mapping
    print(f"\nLoading species mapping: {species_mapping_path}")
    with open(species_mapping_path, 'r') as f:
        mapping_data = json.load(f)

    if 'species_mapping' in mapping_data:
        species_mapping = mapping_data['species_mapping']
    else:
        species_mapping = mapping_data

    print(f"  {len(species_mapping)} species mappings loaded")

    # Load stations
    print(f"\nLoading GRQA stations: {grqa_stations_path}")
    if grqa_stations_path.endswith('.shp'):
        stations_gdf = gpd.read_file(grqa_stations_path)
    else:
        stations_df = pd.read_csv(grqa_stations_path)
        # Detect lat/lon columns
        lat_col = next((c for c in stations_df.columns if c.lower() in ['lat_wgs84', 'lat', 'latitude']), None)
        lon_col = next((c for c in stations_df.columns if c.lower() in ['lon_wgs84', 'lon', 'longitude']), None)
        if lat_col and lon_col:
            geometry = [Point(xy) for xy in zip(stations_df[lon_col], stations_df[lat_col])]
            stations_gdf = gpd.GeoDataFrame(stations_df, geometry=geometry, crs="EPSG:4326")
        else:
            raise ValueError("Could not find lat/lon columns in stations CSV")
    print(f"  {len(stations_gdf)} stations loaded")

    # Load observations
    print(f"\nLoading GRQA observations: {grqa_observations_path}")
    observations_df = pd.read_csv(grqa_observations_path, low_memory=False)
    print(f"  {len(observations_df)} observations loaded")

    # Apply filters
    if years:
        print(f"\nFiltering by years: {years[0]}-{years[1]}")
        date_col = next((c for c in observations_df.columns if c.lower() in ['obs_date', 'date']), None)
        if date_col:
            observations_df[date_col] = pd.to_datetime(observations_df[date_col], errors='coerce')
            observations_df = observations_df[
                (observations_df[date_col].dt.year >= years[0]) &
                (observations_df[date_col].dt.year <= years[1])
            ]
            print(f"  After year filter: {len(observations_df)} observations")

    if species_filter:
        print(f"\nFiltering by species: {species_filter}")
        species_col = next((c for c in observations_df.columns if c.lower() in ['model_species', 'species']), None)
        if species_col:
            observations_df = observations_df[observations_df[species_col].isin(species_filter)]
            print(f"  After species filter: {len(observations_df)} observations")

    if quality_filter:
        outlier_col = next((c for c in observations_df.columns if 'outlier' in c.lower()), None)
        if outlier_col:
            print(f"\nFiltering outliers (column: {outlier_col})")
            before = len(observations_df)
            observations_df = observations_df[
                (observations_df[outlier_col] != True) &
                (observations_df[outlier_col] != 1)
            ]
            print(f"  Removed {before - len(observations_df)} outliers")

    # Match stations to reaches
    matcher = StationReachMatcher(
        river_network_path=river_network_path,
        reach_id_column=reach_id_column,
        max_distance_m=max_match_distance_m
    )

    # Detect station ID column
    site_col = next((c for c in stations_gdf.columns if c.lower() in ['site_id', 'station_id']), 'site_id')
    station_reach_mapping = matcher.match_stations(stations_gdf, station_id_column=site_col)

    # Convert to calibration format
    converter = CalibrationObservationConverter(
        species_mapping=species_mapping,
        target_units=target_units
    )

    calibration_df = converter.convert(
        observations_df=observations_df,
        station_reach_mapping=station_reach_mapping
    )

    if calibration_df.empty:
        print("\nWARNING: No observations could be converted!")
        return calibration_df

    # Save output
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    calibration_df.to_csv(output_path, index=False)
    print(f"\nSaved calibration observations to: {output_path}")

    # Save station-reach mapping for reference
    mapping_output = output_path.parent / f"{output_path.stem}_station_reach_mapping.csv"
    station_reach_mapping.to_csv(mapping_output, index=False)
    print(f"Saved station-reach mapping to: {mapping_output}")

    # Generate summary
    _generate_summary(calibration_df, output_path)

    # Generate calibration stations shapefile
    generate_calibration_stations_shapefile(
        calibration_df=calibration_df,
        stations_gdf=stations_gdf,
        station_reach_mapping=station_reach_mapping,
        output_path=output_path,
        station_id_column=site_col
    )

    print("\n" + "=" * 60)
    print("COMPLETE")
    print("=" * 60)

    return calibration_df


def generate_calibration_stations_shapefile(
    calibration_df: pd.DataFrame,
    stations_gdf: gpd.GeoDataFrame,
    station_reach_mapping: pd.DataFrame,
    output_path: Path,
    station_id_column: str = 'site_id'
) -> gpd.GeoDataFrame:
    """
    Generate a shapefile of calibration stations with metadata.

    Creates a shapefile containing all stations used in calibration with:
    - Station name/ID
    - Corresponding reach_id (for OpenWQ/host model mapping)
    - Parameters available at each station
    - Number of observations per parameter
    - Total observations
    - Distance to matched reach

    Parameters
    ----------
    calibration_df : pd.DataFrame
        Calibration observations dataframe
    stations_gdf : gpd.GeoDataFrame
        Original stations GeoDataFrame with geometries
    station_reach_mapping : pd.DataFrame
        Mapping of stations to reaches
    output_path : Path
        Output path for the shapefile
    station_id_column : str
        Column name for station IDs in stations_gdf

    Returns
    -------
    gpd.GeoDataFrame
        Stations shapefile with metadata
    """
    print("\nGenerating calibration stations shapefile...")

    # Get unique stations from calibration observations
    station_stats = calibration_df.groupby('source').agg({
        'reach_id': 'first',
        'species': lambda x: ','.join(sorted(x.unique())),
        'value': 'count'
    }).reset_index()

    station_stats.columns = ['station_id', 'reach_id', 'parameters', 'total_obs']

    # Count observations per parameter for each station
    param_counts = calibration_df.groupby(['source', 'species']).size().unstack(fill_value=0)
    param_counts = param_counts.reset_index()
    param_counts.columns = ['station_id'] + [f'n_{col}' for col in param_counts.columns[1:]]

    # Merge parameter counts with station stats
    station_stats = station_stats.merge(param_counts, on='station_id', how='left')

    # Get distance from station-reach mapping
    distance_df = station_reach_mapping[['station_id', 'distance_m']].copy()
    station_stats = station_stats.merge(distance_df, on='station_id', how='left')

    # Count number of parameters
    station_stats['n_params'] = station_stats['parameters'].apply(lambda x: len(x.split(',')))

    # Merge with station geometries
    # First, find the matching column in stations_gdf
    if station_id_column in stations_gdf.columns:
        id_col = station_id_column
    else:
        # Try to find a matching column
        id_col = next((c for c in stations_gdf.columns
                      if c.lower() in ['site_id', 'station_id', 'site_no']), None)
        if id_col is None:
            print("  WARNING: Could not find station ID column in stations shapefile")
            return None

    # Convert station_id to same type for merge
    stations_gdf[id_col] = stations_gdf[id_col].astype(str)
    station_stats['station_id'] = station_stats['station_id'].astype(str)

    # Merge geometry
    calibration_stations = stations_gdf[[id_col, 'geometry']].merge(
        station_stats,
        left_on=id_col,
        right_on='station_id',
        how='inner'
    )

    if len(calibration_stations) == 0:
        print("  WARNING: No stations matched between observations and shapefile")
        return None

    # Create GeoDataFrame
    calibration_stations_gdf = gpd.GeoDataFrame(
        calibration_stations,
        geometry='geometry',
        crs=stations_gdf.crs
    )

    # Rename columns for shapefile compatibility (max 10 chars)
    # Use 'spatial_id' as generic name (reach_id for mizuRoute, hruId for SUMMA)
    column_renames = {
        'station_id': 'stn_id',
        'reach_id': 'spatial_id',
        'parameters': 'params',
        'total_obs': 'total_obs',
        'distance_m': 'dist_m',
        'n_params': 'n_params'
    }

    # Rename observation count columns (shorten species names)
    for col in calibration_stations_gdf.columns:
        if col.startswith('n_') and col not in column_renames:
            # Shorten to max 10 chars for shapefile
            short_name = col[:10]
            column_renames[col] = short_name

    calibration_stations_gdf = calibration_stations_gdf.rename(columns=column_renames)

    # Drop the duplicate id column if present
    if id_col in calibration_stations_gdf.columns and 'stn_id' in calibration_stations_gdf.columns:
        if id_col != 'stn_id':
            calibration_stations_gdf = calibration_stations_gdf.drop(columns=[id_col])

    # Save shapefile
    shp_output = output_path.parent / f"{output_path.stem}_stations.shp"
    calibration_stations_gdf.to_file(shp_output)
    print(f"  Saved calibration stations shapefile: {shp_output}")
    print(f"  Total stations: {len(calibration_stations_gdf)}")
    print(f"  Attributes: {list(calibration_stations_gdf.columns)}")

    # Also save as GeoJSON for easier inspection
    geojson_output = output_path.parent / f"{output_path.stem}_stations.geojson"
    calibration_stations_gdf.to_file(geojson_output, driver='GeoJSON')
    print(f"  Saved GeoJSON: {geojson_output}")

    return calibration_stations_gdf


def generate_stations_from_csv(
    calibration_csv_path: str,
    river_network_path: str,
    output_path: str,
    reach_id_column: str = 'seg_id'
) -> gpd.GeoDataFrame:
    """
    Generate calibration stations shapefile from a user-provided CSV.

    When users provide their own CSV observations (not GRQA), this function
    creates a stations shapefile by extracting unique stations and matching
    them to reach centroids.

    Parameters
    ----------
    calibration_csv_path : str
        Path to user-provided calibration CSV
    river_network_path : str
        Path to river network shapefile
    output_path : str
        Output path for the stations shapefile
    reach_id_column : str
        Column in river network containing reach IDs

    Returns
    -------
    gpd.GeoDataFrame
        Stations shapefile with metadata
    """
    print("\nGenerating stations shapefile from CSV observations...")

    # Load calibration CSV
    calibration_df = pd.read_csv(calibration_csv_path)
    print(f"  Loaded {len(calibration_df)} observations from CSV")

    # Load river network
    river_gdf = gpd.read_file(river_network_path)
    print(f"  Loaded {len(river_gdf)} reaches from river network")

    # Get station statistics
    station_stats = calibration_df.groupby('source').agg({
        'reach_id': 'first',
        'species': lambda x: ','.join(sorted(x.unique())),
        'value': 'count'
    }).reset_index()

    station_stats.columns = ['station_id', 'reach_id', 'parameters', 'total_obs']

    # Count observations per parameter
    param_counts = calibration_df.groupby(['source', 'species']).size().unstack(fill_value=0)
    param_counts = param_counts.reset_index()
    param_counts.columns = ['station_id'] + [f'n_{col}' for col in param_counts.columns[1:]]

    station_stats = station_stats.merge(param_counts, on='station_id', how='left')
    station_stats['n_params'] = station_stats['parameters'].apply(lambda x: len(x.split(',')))

    # Get geometry from reach centroids (since we don't have original station locations)
    # Match reach_id to get centroid
    reach_centroids = river_gdf.copy()
    reach_centroids['centroid'] = reach_centroids.geometry.centroid
    reach_centroids = reach_centroids[[reach_id_column, 'centroid']].rename(
        columns={reach_id_column: 'reach_id', 'centroid': 'geometry'}
    )

    # Merge with station stats
    stations_with_geom = station_stats.merge(
        reach_centroids,
        on='reach_id',
        how='left'
    )

    # Filter out stations without geometry (unmatched reaches)
    stations_with_geom = stations_with_geom.dropna(subset=['geometry'])

    if len(stations_with_geom) == 0:
        print("  WARNING: No stations could be matched to reach geometries")
        return None

    # Create GeoDataFrame
    stations_gdf = gpd.GeoDataFrame(
        stations_with_geom,
        geometry='geometry',
        crs=river_gdf.crs
    )

    # Rename for shapefile compatibility
    # Use 'spatial_id' as generic name (reach_id for mizuRoute, hruId for SUMMA)
    column_renames = {
        'station_id': 'stn_id',
        'reach_id': 'spatial_id',
        'parameters': 'params',
        'total_obs': 'total_obs',
        'n_params': 'n_params'
    }

    for col in stations_gdf.columns:
        if col.startswith('n_') and col not in column_renames:
            column_renames[col] = col[:10]

    stations_gdf = stations_gdf.rename(columns=column_renames)

    # Save shapefile
    output_path = Path(output_path)
    shp_output = output_path.parent / f"{output_path.stem}_stations.shp"
    stations_gdf.to_file(shp_output)
    print(f"  Saved stations shapefile: {shp_output}")
    print(f"  Total stations: {len(stations_gdf)}")

    # Also save as GeoJSON
    geojson_output = output_path.parent / f"{output_path.stem}_stations.geojson"
    stations_gdf.to_file(geojson_output, driver='GeoJSON')
    print(f"  Saved GeoJSON: {geojson_output}")

    # Note about geometry source
    print("\n  NOTE: Station locations are based on reach centroids (not original")
    print("        station coordinates) since original locations were not provided.")

    return stations_gdf


def _generate_summary(calibration_df: pd.DataFrame, output_path: Path):
    """Generate summary statistics."""
    summary_path = output_path.parent / f"{output_path.stem}_summary.txt"

    with open(summary_path, 'w') as f:
        f.write("=" * 60 + "\n")
        f.write("CALIBRATION OBSERVATIONS SUMMARY\n")
        f.write("=" * 60 + "\n\n")

        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Output file: {output_path}\n\n")

        f.write(f"Total observations: {len(calibration_df)}\n")
        f.write(f"Number of reaches: {calibration_df['reach_id'].nunique()}\n")
        f.write(f"Number of stations (sources): {calibration_df['source'].nunique()}\n")

        dates = pd.to_datetime(calibration_df['datetime'])
        f.write(f"Date range: {dates.min()} to {dates.max()}\n\n")

        f.write("-" * 60 + "\n")
        f.write("OBSERVATIONS BY SPECIES\n")
        f.write("-" * 60 + "\n")

        for species in sorted(calibration_df['species'].unique()):
            species_data = calibration_df[calibration_df['species'] == species]
            n_obs = len(species_data)
            n_reach = species_data['reach_id'].nunique()
            val_min = species_data['value'].min()
            val_max = species_data['value'].max()
            val_mean = species_data['value'].mean()
            f.write(f"\n{species}:\n")
            f.write(f"  Observations: {n_obs}\n")
            f.write(f"  Reaches: {n_reach}\n")
            f.write(f"  Value range: {val_min:.4f} - {val_max:.4f} (mean: {val_mean:.4f})\n")

        f.write("\n" + "-" * 60 + "\n")
        f.write("OBSERVATIONS BY REACH (top 10)\n")
        f.write("-" * 60 + "\n\n")

        reach_counts = calibration_df.groupby('reach_id').size().sort_values(ascending=False)
        for reach_id, count in reach_counts.head(10).items():
            f.write(f"  Reach {reach_id}: {count} observations\n")

        f.write("\n" + "-" * 60 + "\n")
        f.write("QUALITY FLAGS\n")
        f.write("-" * 60 + "\n\n")

        for flag, count in calibration_df['quality_flag'].value_counts().items():
            pct = count / len(calibration_df) * 100
            f.write(f"  {flag}: {count} ({pct:.1f}%)\n")

    print(f"Saved summary to: {summary_path}")


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Prepare GRQA observations for OpenWQ calibration.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage
  python prepare_calibration_observations.py \\
      --stations grqa_stations.shp \\
      --observations grqa_observations.csv \\
      --river-network river_network.shp \\
      --species-mapping species_mapping.json \\
      --output calibration_observations.csv

  # With reach ID column and distance threshold
  python prepare_calibration_observations.py \\
      --stations grqa_stations.shp \\
      --observations grqa_observations.csv \\
      --river-network river_network.shp \\
      --species-mapping species_mapping.json \\
      --output calibration_observations.csv \\
      --reach-id-column seg_id \\
      --max-distance 500

  # Filter by years and species
  python prepare_calibration_observations.py \\
      --stations grqa_stations.shp \\
      --observations grqa_observations.csv \\
      --river-network river_network.shp \\
      --species-mapping species_mapping.json \\
      --output calibration_observations.csv \\
      --years 2010 2020 \\
      --species NO3 NH4 PO4

Output Format (CSV):
  datetime,reach_id,species,value,units,source,uncertainty,quality_flag
  2010-06-01 00:00:00,1200014181,NO3,2.50,mg/l,USGS_station_A,0.25,GOOD
        """
    )

    # Required arguments
    parser.add_argument(
        '--stations', '-s',
        required=True,
        help='Path to GRQA stations shapefile or CSV'
    )

    parser.add_argument(
        '--observations', '-obs',
        required=True,
        help='Path to GRQA observations CSV'
    )

    parser.add_argument(
        '--river-network', '-r',
        required=True,
        help='Path to model river network shapefile'
    )

    parser.add_argument(
        '--species-mapping', '-m',
        required=True,
        help='Path to species mapping JSON'
    )

    parser.add_argument(
        '--output', '-o',
        required=True,
        help='Output CSV path'
    )

    # Optional arguments
    parser.add_argument(
        '--reach-id-column',
        default='seg_id',
        help='Column name for reach IDs (default: seg_id)'
    )

    parser.add_argument(
        '--max-distance',
        type=float,
        default=1000,
        help='Maximum station-reach matching distance in meters (default: 1000)'
    )

    parser.add_argument(
        '--target-units',
        default='mg/l',
        help='Target units for output (default: mg/l)'
    )

    parser.add_argument(
        '--years',
        type=int,
        nargs=2,
        metavar=('START', 'END'),
        help='Filter by year range'
    )

    parser.add_argument(
        '--species',
        nargs='+',
        help='Filter by species (model names)'
    )

    parser.add_argument(
        '--include-outliers',
        action='store_true',
        help='Include observations flagged as outliers'
    )

    args = parser.parse_args()

    # Validate inputs
    for path_arg in ['stations', 'observations', 'river_network', 'species_mapping']:
        path = getattr(args, path_arg.replace('-', '_'))
        if not os.path.exists(path):
            print(f"Error: File not found: {path}")
            sys.exit(1)

    # Run conversion
    prepare_calibration_observations(
        grqa_stations_path=args.stations,
        grqa_observations_path=args.observations,
        river_network_path=args.river_network,
        species_mapping_path=args.species_mapping,
        output_path=args.output,
        reach_id_column=args.reach_id_column,
        max_match_distance_m=args.max_distance,
        target_units=args.target_units,
        years=tuple(args.years) if args.years else None,
        species_filter=args.species,
        quality_filter=not args.include_outliers
    )


if __name__ == "__main__":
    main()
