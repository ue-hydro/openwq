#!/usr/bin/env python3
"""
GRQA Station Extractor for River Networks - Calibration Workflow Integration

This script extracts water quality monitoring stations and observations from the
Global River Water Quality Archive (GRQA) dataset based on a river network shapefile,
with support for species mapping to model variables for calibration workflows.

GRQA Dataset: https://zenodo.org/records/15335450
Reference: Virro et al. (2021) - GRQA: Global River Water Quality Archive
           Earth System Science Data, 13, 5483-5507

Workflow Integration:
1. User provides river network shapefile
2. User provides species mapping (GRQA parameters -> Model species names)
3. Script downloads ONLY the required GRQA parameter files
4. Extracts stations within buffer of river network
5. Outputs calibration-ready files with model species names

Author: OpenWQ Development Team
Date: 2024
"""

import os
import sys
import json
import argparse
import zipfile
import requests
from pathlib import Path
from typing import List, Dict, Optional, Tuple, Union
from datetime import datetime
import warnings

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore')

try:
    import pandas as pd
    import geopandas as gpd
    from shapely.geometry import Point, box
    from shapely.ops import unary_union
    import numpy as np
except ImportError as e:
    print(f"Error: Missing required package. Please install: {e.name}")
    print("Install with: pip install pandas geopandas shapely numpy requests")
    sys.exit(1)


# GRQA Dataset Information
GRQA_ZENODO_RECORD = "15335450"
GRQA_ZENODO_BASE = f"https://zenodo.org/records/{GRQA_ZENODO_RECORD}"
GRQA_DATA_URL = f"{GRQA_ZENODO_BASE}/files/GRQA_data_v1.4.zip"

# GRQA Parameters (43 total as of v1.4) with units
GRQA_PARAMETERS = {
    # Nutrients - Nitrogen
    'TN': {'name': 'Total Nitrogen', 'unit': 'mg/L', 'group': 'nitrogen'},
    'DN': {'name': 'Dissolved Nitrogen', 'unit': 'mg/L', 'group': 'nitrogen'},
    'NO3': {'name': 'Nitrate', 'unit': 'mg/L', 'group': 'nitrogen'},
    'NO2': {'name': 'Nitrite', 'unit': 'mg/L', 'group': 'nitrogen'},
    'NH4': {'name': 'Ammonium', 'unit': 'mg/L', 'group': 'nitrogen'},
    'NO3_NO2': {'name': 'Nitrate + Nitrite', 'unit': 'mg/L', 'group': 'nitrogen'},
    'DIN': {'name': 'Dissolved Inorganic Nitrogen', 'unit': 'mg/L', 'group': 'nitrogen'},
    'DON': {'name': 'Dissolved Organic Nitrogen', 'unit': 'mg/L', 'group': 'nitrogen'},
    'TKN': {'name': 'Total Kjeldahl Nitrogen', 'unit': 'mg/L', 'group': 'nitrogen'},
    'PN': {'name': 'Particulate Nitrogen', 'unit': 'mg/L', 'group': 'nitrogen'},

    # Nutrients - Phosphorus
    'TP': {'name': 'Total Phosphorus', 'unit': 'mg/L', 'group': 'phosphorus'},
    'DP': {'name': 'Dissolved Phosphorus', 'unit': 'mg/L', 'group': 'phosphorus'},
    'DIP': {'name': 'Dissolved Inorganic Phosphorus', 'unit': 'mg/L', 'group': 'phosphorus'},
    'PO4': {'name': 'Phosphate', 'unit': 'mg/L', 'group': 'phosphorus'},
    'DOP': {'name': 'Dissolved Organic Phosphorus', 'unit': 'mg/L', 'group': 'phosphorus'},
    'PP': {'name': 'Particulate Phosphorus', 'unit': 'mg/L', 'group': 'phosphorus'},

    # Oxygen
    'DO': {'name': 'Dissolved Oxygen', 'unit': 'mg/L', 'group': 'oxygen'},
    'DO_sat': {'name': 'Dissolved Oxygen Saturation', 'unit': '%', 'group': 'oxygen'},
    'BOD': {'name': 'Biochemical Oxygen Demand', 'unit': 'mg/L', 'group': 'oxygen'},
    'BOD5': {'name': 'Biochemical Oxygen Demand (5-day)', 'unit': 'mg/L', 'group': 'oxygen'},
    'COD': {'name': 'Chemical Oxygen Demand', 'unit': 'mg/L', 'group': 'oxygen'},
    'CODMn': {'name': 'Chemical Oxygen Demand (Permanganate)', 'unit': 'mg/L', 'group': 'oxygen'},

    # Carbon
    'DOC': {'name': 'Dissolved Organic Carbon', 'unit': 'mg/L', 'group': 'carbon'},
    'TOC': {'name': 'Total Organic Carbon', 'unit': 'mg/L', 'group': 'carbon'},
    'TC': {'name': 'Total Carbon', 'unit': 'mg/L', 'group': 'carbon'},
    'TIC': {'name': 'Total Inorganic Carbon', 'unit': 'mg/L', 'group': 'carbon'},
    'DIC': {'name': 'Dissolved Inorganic Carbon', 'unit': 'mg/L', 'group': 'carbon'},
    'POC': {'name': 'Particulate Organic Carbon', 'unit': 'mg/L', 'group': 'carbon'},
    'PC': {'name': 'Particulate Carbon', 'unit': 'mg/L', 'group': 'carbon'},

    # Sediments
    'TSS': {'name': 'Total Suspended Solids', 'unit': 'mg/L', 'group': 'sediment'},
    'TDS': {'name': 'Total Dissolved Solids', 'unit': 'mg/L', 'group': 'sediment'},
    'SS': {'name': 'Suspended Solids', 'unit': 'mg/L', 'group': 'sediment'},
    'Turbidity': {'name': 'Turbidity', 'unit': 'NTU', 'group': 'sediment'},

    # Physical
    'Temp': {'name': 'Water Temperature', 'unit': '°C', 'group': 'physical'},
    'pH': {'name': 'pH', 'unit': 'pH', 'group': 'physical'},
    'EC': {'name': 'Electrical Conductivity', 'unit': 'µS/cm', 'group': 'physical'},
    'SC': {'name': 'Specific Conductance', 'unit': 'µS/cm', 'group': 'physical'},

    # Other
    'Chl_a': {'name': 'Chlorophyll-a', 'unit': 'µg/L', 'group': 'biological'},
    'Secchi': {'name': 'Secchi Depth', 'unit': 'm', 'group': 'physical'},
    'Si': {'name': 'Silica', 'unit': 'mg/L', 'group': 'other'},
    'SO4': {'name': 'Sulfate', 'unit': 'mg/L', 'group': 'other'},
    'Cl': {'name': 'Chloride', 'unit': 'mg/L', 'group': 'other'},
    'Alk': {'name': 'Alkalinity', 'unit': 'mg/L CaCO3', 'group': 'other'},
}

# Default species mapping templates for common model configurations
DEFAULT_SPECIES_MAPPINGS = {
    'nitrogen_cycle': {
        'NO3': 'NO3',
        'NH4': 'NH4',
        'NO2': 'NO2',
        'TN': 'TN',
        'DON': 'DON',
    },
    'phosphorus_cycle': {
        'PO4': 'PO4',
        'TP': 'TP',
        'DP': 'DP',
    },
    'carbon_cycle': {
        'DOC': 'DOC',
        'DIC': 'DIC',
        'TOC': 'TOC',
    },
    'oxygen': {
        'DO': 'DO',
        'BOD5': 'BOD',
    },
    'full_nutrient': {
        'NO3': 'NO3',
        'NH4': 'NH4',
        'PO4': 'PO4',
        'DO': 'DO',
        'DOC': 'DOC',
        'Temp': 'TEMP',
    }
}


class SpeciesMapper:
    """
    Handles mapping between GRQA parameter names and model species names.
    """

    def __init__(self, mapping: Dict[str, str] = None, mapping_file: str = None):
        """
        Initialize species mapper.

        Parameters
        ----------
        mapping : dict, optional
            Dictionary mapping GRQA parameters to model species names
            e.g., {'NO3': 'nitrate', 'NH4': 'ammonium'}
        mapping_file : str, optional
            Path to JSON file containing the mapping
        """
        self.mapping = {}
        self.reverse_mapping = {}

        if mapping_file:
            self.load_mapping(mapping_file)
        elif mapping:
            self.set_mapping(mapping)

    def set_mapping(self, mapping: Dict[str, str]):
        """Set the species mapping."""
        self.mapping = mapping
        self.reverse_mapping = {v: k for k, v in mapping.items()}

    def load_mapping(self, filepath: str):
        """Load mapping from JSON file."""
        with open(filepath, 'r') as f:
            data = json.load(f)

        # Support both flat mapping and structured format
        if 'species_mapping' in data:
            self.mapping = data['species_mapping']
        else:
            self.mapping = data

        self.reverse_mapping = {v: k for k, v in self.mapping.items()}
        print(f"Loaded species mapping from: {filepath}")
        print(f"  {len(self.mapping)} species mapped")

    def save_mapping(self, filepath: str, include_metadata: bool = True):
        """Save mapping to JSON file."""
        if include_metadata:
            data = {
                'description': 'GRQA to Model Species Mapping for Calibration',
                'created': datetime.now().isoformat(),
                'grqa_source': GRQA_ZENODO_BASE,
                'species_mapping': self.mapping,
                'grqa_parameters_info': {
                    k: GRQA_PARAMETERS[k] for k in self.mapping.keys()
                    if k in GRQA_PARAMETERS
                }
            }
        else:
            data = self.mapping

        with open(filepath, 'w') as f:
            json.dump(data, f, indent=2)

        print(f"Saved species mapping to: {filepath}")

    def get_grqa_parameters(self) -> List[str]:
        """Get list of GRQA parameters needed based on mapping."""
        return list(self.mapping.keys())

    def get_model_name(self, grqa_param: str) -> str:
        """Get model species name for a GRQA parameter."""
        return self.mapping.get(grqa_param, grqa_param)

    def get_grqa_name(self, model_species: str) -> str:
        """Get GRQA parameter name for a model species."""
        return self.reverse_mapping.get(model_species, model_species)

    def validate(self) -> Tuple[bool, List[str]]:
        """Validate that all mapped GRQA parameters exist."""
        invalid = [p for p in self.mapping.keys() if p not in GRQA_PARAMETERS]
        return len(invalid) == 0, invalid

    @staticmethod
    def create_template(template_name: str = None) -> Dict[str, str]:
        """Create a mapping template."""
        if template_name and template_name in DEFAULT_SPECIES_MAPPINGS:
            return DEFAULT_SPECIES_MAPPINGS[template_name].copy()
        else:
            # Return all parameters with same name
            return {k: k for k in GRQA_PARAMETERS.keys()}

    @staticmethod
    def list_templates():
        """List available mapping templates."""
        print("\nAvailable Species Mapping Templates:")
        print("=" * 50)
        for name, mapping in DEFAULT_SPECIES_MAPPINGS.items():
            print(f"\n{name}:")
            for grqa, model in mapping.items():
                info = GRQA_PARAMETERS.get(grqa, {})
                print(f"  {grqa:12s} -> {model:12s} ({info.get('name', 'N/A')})")


class GRQADownloader:
    """
    Handles selective downloading of GRQA data files, or loading from local storage.
    """

    def __init__(self, cache_dir: str, local_data_path: str = None):
        """
        Initialize downloader.

        Parameters
        ----------
        cache_dir : str
            Directory to cache downloaded files
        local_data_path : str, optional
            Path to already-downloaded GRQA data folder.
            If provided, data will be loaded from here instead of downloading.
            Expected structure:
              local_data_path/
                ├── GRQA_data_v1.3/ or GRQA_data_v1.4/
                │   ├── GRQA_data_v1.3_NO3.csv (or similar pattern)
                │   └── ...
                └── GRQA_meta_v1.3.csv (optional)
        """
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)

        self.local_data_path = Path(local_data_path) if local_data_path else None
        self.use_local = self.local_data_path is not None

        if self.use_local:
            if not self.local_data_path.exists():
                raise FileNotFoundError(
                    f"Local GRQA data path not found: {self.local_data_path}"
                )
            print(f"Using local GRQA data from: {self.local_data_path}")
            self._detect_local_structure()

    def _detect_local_structure(self):
        """Detect the structure of local GRQA data and find parameter files."""
        self.local_file_patterns = []
        self.local_version = None

        # Look for GRQA data directories or files
        for pattern in ['GRQA_data_v*', 'grqa_data_v*']:
            for path in self.local_data_path.glob(pattern):
                if path.is_dir():
                    # Found a versioned directory
                    version = path.name.split('_')[-1]  # e.g., "v1.4"
                    self.local_version = version
                    print(f"  Detected GRQA version: {version}")

                    # Look for CSV files inside
                    csv_files = list(path.glob('*.csv'))
                    if csv_files:
                        # Detect file naming pattern
                        sample = csv_files[0].name
                        if '_GRQA.csv' in sample:
                            self.local_file_patterns.append(
                                (path, '{param}_GRQA.csv')
                            )
                        elif 'GRQA_data_' in sample:
                            # Pattern: GRQA_data_v1.3_NO3.csv
                            self.local_file_patterns.append(
                                (path, f'GRQA_data_{version}_{{param}}.csv')
                            )
                        print(f"  Found {len(csv_files)} parameter files")
                    break

        # Also check for flat structure (files directly in local_data_path)
        direct_csvs = list(self.local_data_path.glob('*_GRQA.csv'))
        if direct_csvs:
            self.local_file_patterns.append(
                (self.local_data_path, '{param}_GRQA.csv')
            )
            print(f"  Found {len(direct_csvs)} parameter files in root")

        if not self.local_file_patterns:
            print("  Warning: No GRQA CSV files found. Will attempt download.")
            self.use_local = False

    def _find_local_file(self, parameter: str) -> Optional[Path]:
        """Find a parameter file in local data."""
        if not self.use_local:
            return None

        for base_path, pattern in self.local_file_patterns:
            filename = pattern.format(param=parameter)
            filepath = base_path / filename
            if filepath.exists():
                return filepath

            # Try case-insensitive search
            for f in base_path.glob(f'*{parameter}*'):
                if f.suffix.lower() == '.csv':
                    return f

        return None

    def get_file_url(self, parameter: str) -> str:
        """Get direct download URL for a parameter file."""
        # Individual files within the zip follow this pattern
        return f"{GRQA_DATA_URL}"

    def download_parameter_file(self, parameter: str, force: bool = False) -> Optional[Path]:
        """
        Get a specific parameter file from GRQA (local or download).

        If local_data_path was provided, uses local files.
        Otherwise, downloads from Zenodo (one-time download of full archive).
        """
        # First, check for local data
        if self.use_local:
            local_file = self._find_local_file(parameter)
            if local_file:
                return local_file
            else:
                print(f"  Warning: {parameter} not found in local data, will try download")

        # Check if already extracted to cache
        param_file = self.cache_dir / f"{parameter}_GRQA.csv"
        if param_file.exists() and not force:
            return param_file

        # Check if zip exists
        zip_path = self.cache_dir / "GRQA_data_v1.4.zip"

        if not zip_path.exists():
            print(f"Downloading GRQA data archive...")
            print(f"  This is a one-time download (~1.2 GB)")
            self._download_with_progress(GRQA_DATA_URL, zip_path)

        # Extract specific file
        return self._extract_parameter(zip_path, parameter)

    def _download_with_progress(self, url: str, dest: Path):
        """Download file with progress indicator."""
        response = requests.get(url, stream=True)
        total_size = int(response.headers.get('content-length', 0))

        with open(dest, 'wb') as f:
            downloaded = 0
            for chunk in response.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)
                    downloaded += len(chunk)
                    if total_size > 0:
                        pct = (downloaded / total_size) * 100
                        mb = downloaded / 1e6
                        print(f"\r  Downloaded: {mb:.1f} MB ({pct:.1f}%)", end="", flush=True)
        print()

    def _extract_parameter(self, zip_path: Path, parameter: str) -> Optional[Path]:
        """Extract a specific parameter file from the zip."""
        target_pattern = f"{parameter}_GRQA.csv"

        with zipfile.ZipFile(zip_path, 'r') as zf:
            for name in zf.namelist():
                if name.endswith(target_pattern):
                    # Extract to cache directory
                    zf.extract(name, self.cache_dir)
                    # Move to flat structure
                    extracted = self.cache_dir / name
                    dest = self.cache_dir / target_pattern
                    if extracted != dest:
                        extracted.rename(dest)
                    return dest

        print(f"  Warning: Parameter {parameter} not found in archive")
        return None

    def download_required_parameters(self,
                                      parameters: List[str],
                                      force: bool = False) -> Dict[str, Path]:
        """
        Get all required parameter files (from local storage or download).

        Parameters
        ----------
        parameters : list
            List of GRQA parameter codes
        force : bool
            Force re-download (ignored for local files)

        Returns
        -------
        dict
            Mapping of parameter codes to file paths
        """
        source_type = "local" if self.use_local else "remote"
        print(f"\nPreparing GRQA data for {len(parameters)} parameters ({source_type})...")

        available = {}
        missing = []

        for param in parameters:
            # First check local data
            if self.use_local:
                local_file = self._find_local_file(param)
                if local_file:
                    available[param] = local_file
                    print(f"  {param}: local")
                    continue

            # Then check cache
            param_file = self.cache_dir / f"{param}_GRQA.csv"
            if param_file.exists() and not force:
                available[param] = param_file
                print(f"  {param}: cached")
            else:
                missing.append(param)

        # Download missing (only if not all local)
        if missing:
            print(f"\nDownloading {len(missing)} parameter files...")
            for param in missing:
                filepath = self.download_parameter_file(param, force)
                if filepath:
                    available[param] = filepath
                    print(f"  {param}: downloaded")

        return available


class GRQACalibrationExtractor:
    """
    Extract GRQA water quality data for model calibration workflows.
    """

    def __init__(self,
                 output_dir: str,
                 species_mapper: SpeciesMapper,
                 cache_dir: str = None,
                 local_data_path: str = None,
                 buffer_distance_m: float = 1000,
                 years: Optional[Tuple[int, int]] = None):
        """
        Initialize the calibration extractor.

        Parameters
        ----------
        output_dir : str
            Directory to save extracted calibration data
        species_mapper : SpeciesMapper
            Species mapping configuration
        cache_dir : str, optional
            Directory to cache GRQA downloads (default: output_dir/grqa_cache)
        local_data_path : str, optional
            Path to already-downloaded GRQA data folder.
            If provided, data will be loaded from here instead of downloading.
        buffer_distance_m : float
            Buffer distance around river network in meters
        years : tuple, optional
            (start_year, end_year) to filter observations
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        self.species_mapper = species_mapper
        self.buffer_distance_m = buffer_distance_m
        self.years = years
        self.local_data_path = local_data_path

        # Set up cache and downloader
        self.cache_dir = Path(cache_dir) if cache_dir else self.output_dir / "grqa_cache"
        self.downloader = GRQADownloader(self.cache_dir, local_data_path=local_data_path)

        # Column mapping (will be detected from data)
        self.column_map = None

        # Results storage
        self.stations_gdf = None
        self.observations_df = None

    def load_river_network(self, shapefile_path: str) -> gpd.GeoDataFrame:
        """Load and validate river network shapefile."""
        print(f"\nLoading river network: {shapefile_path}")

        gdf = gpd.read_file(shapefile_path)

        if gdf.crs is None:
            print("  Warning: No CRS detected, assuming EPSG:4326 (WGS84)")
            gdf = gdf.set_crs("EPSG:4326")
        elif gdf.crs.to_epsg() != 4326:
            print(f"  Reprojecting from {gdf.crs} to EPSG:4326")
            gdf = gdf.to_crs("EPSG:4326")

        print(f"  Features: {len(gdf)}")
        print(f"  Bounds: {gdf.total_bounds}")

        return gdf

    def create_buffer(self, river_gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
        """Create buffer around river network for station selection."""
        print(f"\nCreating {self.buffer_distance_m}m buffer...")

        # Determine appropriate UTM zone
        centroid = river_gdf.dissolve().centroid.iloc[0]
        utm_zone = int((centroid.x + 180) / 6) + 1
        hemisphere = 'north' if centroid.y >= 0 else 'south'
        epsg_code = 32600 + utm_zone if hemisphere == 'north' else 32700 + utm_zone

        # Project, buffer, dissolve, reproject
        river_projected = river_gdf.to_crs(f"EPSG:{epsg_code}")
        buffered = river_projected.buffer(self.buffer_distance_m)
        buffer_gdf = gpd.GeoDataFrame(
            geometry=[unary_union(buffered)],
            crs=f"EPSG:{epsg_code}"
        ).to_crs("EPSG:4326")

        print(f"  Using UTM zone {utm_zone}{hemisphere[0].upper()}")

        return buffer_gdf

    def _detect_columns(self, df: pd.DataFrame) -> Dict[str, str]:
        """Auto-detect column names in GRQA CSV."""
        columns = [c.lower() for c in df.columns]
        original = df.columns.tolist()

        patterns = {
            'site_id': ['site_id', 'site_no', 'station_id'],
            'site_name': ['site_name', 'station_name'],
            'lat': ['lat_wgs84', 'lat', 'latitude'],
            'lon': ['lon_wgs84', 'lon', 'longitude'],
            'obs_date': ['obs_date', 'date', 'sample_date'],
            'obs_time': ['obs_time', 'time'],
            'obs_value': ['obs_value', 'value', 'result'],
            'unit': ['unit', 'units'],
            'source': ['source', 'data_source'],
            'outlier': ['obs_iqr_outlier', 'outlier', 'is_outlier'],
        }

        col_map = {}
        for field, candidates in patterns.items():
            for i, col in enumerate(columns):
                if col in candidates:
                    col_map[field] = original[i]
                    break

        return col_map

    def extract_stations_and_observations(self,
                                           buffer_gdf: gpd.GeoDataFrame) -> Tuple[gpd.GeoDataFrame, pd.DataFrame]:
        """
        Extract stations and observations for mapped species.

        Returns
        -------
        tuple
            (stations_gdf, observations_df)
        """
        parameters = self.species_mapper.get_grqa_parameters()
        print(f"\nExtracting data for {len(parameters)} parameters:")
        for p in parameters:
            model_name = self.species_mapper.get_model_name(p)
            info = GRQA_PARAMETERS.get(p, {})
            print(f"  {p} -> {model_name} ({info.get('name', 'N/A')})")

        # Download required files
        param_files = self.downloader.download_required_parameters(parameters)

        if not param_files:
            print("Error: No parameter files available!")
            return gpd.GeoDataFrame(), pd.DataFrame()

        # Extract stations and observations
        all_stations = []
        all_observations = []
        buffer_union = buffer_gdf.geometry.unary_union

        for param, filepath in param_files.items():
            print(f"\nProcessing {param}...")

            try:
                df = pd.read_csv(filepath, low_memory=False)

                # Detect columns on first file
                if self.column_map is None:
                    self.column_map = self._detect_columns(df)
                    print(f"  Detected columns: {list(self.column_map.keys())}")

                lat_col = self.column_map.get('lat', 'lat_wgs84')
                lon_col = self.column_map.get('lon', 'lon_wgs84')
                site_col = self.column_map.get('site_id', 'site_id')
                date_col = self.column_map.get('obs_date', 'obs_date')
                value_col = self.column_map.get('obs_value', 'obs_value')

                # Create geometry and filter by buffer
                geometry = [Point(xy) for xy in zip(df[lon_col], df[lat_col])]
                gdf = gpd.GeoDataFrame(df, geometry=geometry, crs="EPSG:4326")
                gdf_filtered = gdf[gdf.geometry.within(buffer_union)].copy()

                if len(gdf_filtered) == 0:
                    print(f"  No observations within buffer")
                    continue

                # Filter by years if specified
                if self.years:
                    gdf_filtered[date_col] = pd.to_datetime(gdf_filtered[date_col], errors='coerce')
                    gdf_filtered = gdf_filtered[
                        (gdf_filtered[date_col].dt.year >= self.years[0]) &
                        (gdf_filtered[date_col].dt.year <= self.years[1])
                    ]

                # Add model species name
                model_name = self.species_mapper.get_model_name(param)
                gdf_filtered['grqa_parameter'] = param
                gdf_filtered['model_species'] = model_name

                # Collect unique stations
                station_cols = [site_col, lat_col, lon_col, 'geometry']
                if self.column_map.get('site_name') in gdf_filtered.columns:
                    station_cols.insert(1, self.column_map['site_name'])

                stations = gdf_filtered[station_cols].drop_duplicates(subset=[site_col])
                stations['has_' + model_name] = True
                all_stations.append(stations)

                # Collect observations
                all_observations.append(gdf_filtered)

                print(f"  Found {len(gdf_filtered)} observations at {len(stations)} stations")

            except Exception as e:
                print(f"  Error: {e}")
                continue

        if not all_stations:
            print("\nNo data found within the buffer area!")
            return gpd.GeoDataFrame(), pd.DataFrame()

        # Combine stations (merge to track which species each station has)
        site_col = self.column_map.get('site_id', 'site_id')
        stations_combined = all_stations[0]
        for s in all_stations[1:]:
            stations_combined = stations_combined.merge(
                s.drop(columns='geometry'),
                on=site_col,
                how='outer',
                suffixes=('', '_dup')
            )
            # Clean up duplicate columns
            dup_cols = [c for c in stations_combined.columns if c.endswith('_dup')]
            stations_combined.drop(columns=dup_cols, inplace=True)

        # Fill NaN in has_* columns with False
        has_cols = [c for c in stations_combined.columns if c.startswith('has_')]
        stations_combined[has_cols] = stations_combined[has_cols].fillna(False)

        # Recreate geometry for combined stations
        lat_col = self.column_map.get('lat', 'lat_wgs84')
        lon_col = self.column_map.get('lon', 'lon_wgs84')
        stations_combined = gpd.GeoDataFrame(
            stations_combined.drop(columns='geometry'),
            geometry=[Point(xy) for xy in zip(
                stations_combined[lon_col], stations_combined[lat_col]
            )],
            crs="EPSG:4326"
        )

        # Combine observations
        observations_combined = pd.concat(all_observations, ignore_index=True)
        observations_df = pd.DataFrame(observations_combined.drop(columns='geometry'))

        self.stations_gdf = stations_combined
        self.observations_df = observations_df

        print(f"\n  Total stations: {len(stations_combined)}")
        print(f"  Total observations: {len(observations_df)}")

        return stations_combined, observations_df

    def generate_calibration_files(self, prefix: str = "calibration"):
        """
        Generate calibration-ready output files.

        Creates:
        - stations.shp/csv: Station locations with available species
        - observations.csv: All observations with model species names
        - {species}_timeseries.csv: Per-species time series at each station
        - calibration_config.json: Configuration for calibration workflow
        """
        if self.stations_gdf is None or self.observations_df is None:
            print("Error: No data to save. Run extraction first.")
            return

        print(f"\nGenerating calibration files in: {self.output_dir}")

        # 1. Save stations
        stations_shp = self.output_dir / f"{prefix}_stations.shp"
        stations_csv = self.output_dir / f"{prefix}_stations.csv"
        self.stations_gdf.to_file(stations_shp)
        self.stations_gdf.drop(columns='geometry').to_csv(stations_csv, index=False)
        print(f"  Stations: {stations_shp}")

        # 2. Save all observations
        obs_csv = self.output_dir / f"{prefix}_observations.csv"
        self.observations_df.to_csv(obs_csv, index=False)
        print(f"  Observations: {obs_csv}")

        # 3. Create per-species time series files
        timeseries_dir = self.output_dir / f"{prefix}_timeseries"
        timeseries_dir.mkdir(exist_ok=True)

        site_col = self.column_map.get('site_id', 'site_id')
        date_col = self.column_map.get('obs_date', 'obs_date')
        value_col = self.column_map.get('obs_value', 'obs_value')

        for model_species in self.observations_df['model_species'].unique():
            species_df = self.observations_df[
                self.observations_df['model_species'] == model_species
            ].copy()

            # Pivot to wide format: rows=dates, columns=stations
            if date_col in species_df.columns and site_col in species_df.columns:
                pivot_df = species_df.pivot_table(
                    index=date_col,
                    columns=site_col,
                    values=value_col,
                    aggfunc='mean'  # Average if multiple obs per day
                )
                ts_file = timeseries_dir / f"{model_species}_timeseries.csv"
                pivot_df.to_csv(ts_file)
                print(f"  Time series: {ts_file}")

        # 4. Generate calibration configuration
        config = self._generate_calibration_config(prefix)
        config_file = self.output_dir / f"{prefix}_config.json"
        with open(config_file, 'w') as f:
            json.dump(config, f, indent=2)
        print(f"  Config: {config_file}")

        # 5. Generate summary report
        self._generate_summary_report(prefix)

    def _generate_calibration_config(self, prefix: str) -> Dict:
        """Generate calibration workflow configuration."""
        site_col = self.column_map.get('site_id', 'site_id')
        date_col = self.column_map.get('obs_date', 'obs_date')

        # Get date range
        dates = pd.to_datetime(self.observations_df[date_col], errors='coerce')

        config = {
            'metadata': {
                'created': datetime.now().isoformat(),
                'grqa_version': '1.4',
                'grqa_source': GRQA_ZENODO_BASE,
                'buffer_distance_m': self.buffer_distance_m,
                'year_filter': self.years,
            },
            'species_mapping': self.species_mapper.mapping,
            'stations': {
                'count': len(self.stations_gdf),
                'file': f"{prefix}_stations.csv",
                'shapefile': f"{prefix}_stations.shp",
                'id_column': site_col,
            },
            'observations': {
                'count': len(self.observations_df),
                'file': f"{prefix}_observations.csv",
                'date_range': {
                    'start': str(dates.min().date()) if not dates.isna().all() else None,
                    'end': str(dates.max().date()) if not dates.isna().all() else None,
                },
                'species': list(self.observations_df['model_species'].unique()),
            },
            'timeseries_dir': f"{prefix}_timeseries",
            'calibration_targets': []
        }

        # Add per-species statistics
        for species in self.observations_df['model_species'].unique():
            species_data = self.observations_df[
                self.observations_df['model_species'] == species
            ]
            value_col = self.column_map.get('obs_value', 'obs_value')

            config['calibration_targets'].append({
                'model_species': species,
                'grqa_parameter': self.species_mapper.get_grqa_name(species),
                'observation_count': len(species_data),
                'station_count': species_data[site_col].nunique(),
                'value_stats': {
                    'min': float(species_data[value_col].min()),
                    'max': float(species_data[value_col].max()),
                    'mean': float(species_data[value_col].mean()),
                    'median': float(species_data[value_col].median()),
                },
                'timeseries_file': f"{prefix}_timeseries/{species}_timeseries.csv"
            })

        return config

    def _generate_summary_report(self, prefix: str):
        """Generate human-readable summary report."""
        report_file = self.output_dir / f"{prefix}_summary.txt"

        with open(report_file, 'w') as f:
            f.write("=" * 60 + "\n")
            f.write("GRQA Calibration Data Extraction Summary\n")
            f.write("=" * 60 + "\n\n")

            f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"GRQA Source: {GRQA_ZENODO_BASE}\n")
            f.write(f"Buffer Distance: {self.buffer_distance_m} m\n")
            if self.years:
                f.write(f"Year Range: {self.years[0]} - {self.years[1]}\n")
            f.write("\n")

            f.write("-" * 60 + "\n")
            f.write("SPECIES MAPPING\n")
            f.write("-" * 60 + "\n")
            for grqa, model in self.species_mapper.mapping.items():
                info = GRQA_PARAMETERS.get(grqa, {})
                f.write(f"  {grqa:12s} -> {model:12s} ({info.get('name', 'N/A')})\n")
            f.write("\n")

            f.write("-" * 60 + "\n")
            f.write("DATA SUMMARY\n")
            f.write("-" * 60 + "\n")
            f.write(f"  Total Stations: {len(self.stations_gdf)}\n")
            f.write(f"  Total Observations: {len(self.observations_df)}\n\n")

            f.write("  Per-Species Summary:\n")
            site_col = self.column_map.get('site_id', 'site_id')
            for species in sorted(self.observations_df['model_species'].unique()):
                species_data = self.observations_df[
                    self.observations_df['model_species'] == species
                ]
                n_obs = len(species_data)
                n_sta = species_data[site_col].nunique()
                f.write(f"    {species:12s}: {n_obs:8,} obs at {n_sta:4} stations\n")

            f.write("\n")
            f.write("-" * 60 + "\n")
            f.write("OUTPUT FILES\n")
            f.write("-" * 60 + "\n")
            f.write(f"  {prefix}_stations.shp      - Station locations (shapefile)\n")
            f.write(f"  {prefix}_stations.csv      - Station locations (CSV)\n")
            f.write(f"  {prefix}_observations.csv  - All observations\n")
            f.write(f"  {prefix}_timeseries/       - Per-species time series\n")
            f.write(f"  {prefix}_config.json       - Calibration configuration\n")
            f.write("\n")

        print(f"  Summary: {report_file}")

    def run(self, shapefile_path: str, prefix: str = "calibration"):
        """
        Run the complete calibration data extraction pipeline.

        Parameters
        ----------
        shapefile_path : str
            Path to river network shapefile
        prefix : str
            Output file prefix
        """
        print("\n" + "=" * 60)
        print("GRQA Calibration Data Extractor")
        print("=" * 60)

        # Validate species mapping
        valid, invalid = self.species_mapper.validate()
        if not valid:
            print(f"\nWarning: Unknown GRQA parameters: {invalid}")
            print("These will be skipped.")

        # Load river network
        river_gdf = self.load_river_network(shapefile_path)

        # Create buffer
        buffer_gdf = self.create_buffer(river_gdf)

        # Extract data
        stations, observations = self.extract_stations_and_observations(buffer_gdf)

        if stations.empty:
            print("\nNo data found. Check your inputs.")
            return None, None

        # Generate calibration files
        self.generate_calibration_files(prefix)

        print("\n" + "=" * 60)
        print("Extraction Complete!")
        print("=" * 60 + "\n")

        return stations, observations


def create_example_mapping_file(output_path: str):
    """Create an example species mapping JSON file."""
    example = {
        "description": "Example GRQA to Model Species Mapping",
        "instructions": [
            "Edit the 'species_mapping' section below",
            "Keys are GRQA parameter codes (see 'available_parameters')",
            "Values are your model's species names",
            "Remove any parameters you don't need"
        ],
        "species_mapping": {
            "NO3": "nitrate",
            "NH4": "ammonium",
            "PO4": "phosphate",
            "DO": "dissolved_oxygen",
            "DOC": "dissolved_organic_carbon",
            "Temp": "water_temperature"
        },
        "available_parameters": {
            code: info['name'] for code, info in GRQA_PARAMETERS.items()
        }
    }

    with open(output_path, 'w') as f:
        json.dump(example, f, indent=2)

    print(f"Created example mapping file: {output_path}")
    print("Edit this file to configure your species mapping.")


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Extract GRQA water quality data for model calibration.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Create example species mapping file
  python grqa_extract_stations.py --create-mapping my_mapping.json

  # Extract with custom mapping
  python grqa_extract_stations.py river.shp -m my_mapping.json -o ./calibration_data

  # Use built-in template
  python grqa_extract_stations.py river.shp --template nitrogen_cycle -o ./output

  # Specify years and buffer
  python grqa_extract_stations.py river.shp -m mapping.json -o ./output \\
      --years 2015 2023 --buffer 500

  # List available templates
  python grqa_extract_stations.py --list-templates

  # List all GRQA parameters
  python grqa_extract_stations.py --list-parameters

Species Mapping Templates:
  nitrogen_cycle   : NO3, NH4, NO2, TN, DON
  phosphorus_cycle : PO4, TP, DP
  carbon_cycle     : DOC, DIC, TOC
  oxygen           : DO, BOD5
  full_nutrient    : NO3, NH4, PO4, DO, DOC, Temp

GRQA Dataset: https://zenodo.org/records/15335450
        """
    )

    parser.add_argument(
        'shapefile',
        type=str,
        nargs='?',
        help='Path to river network shapefile'
    )

    parser.add_argument(
        '-m', '--mapping',
        type=str,
        help='Path to species mapping JSON file'
    )

    parser.add_argument(
        '-t', '--template',
        type=str,
        choices=list(DEFAULT_SPECIES_MAPPINGS.keys()),
        help='Use a built-in species mapping template'
    )

    parser.add_argument(
        '-o', '--output',
        type=str,
        default='./grqa_calibration',
        help='Output directory (default: ./grqa_calibration)'
    )

    parser.add_argument(
        '-b', '--buffer',
        type=float,
        default=1000,
        help='Buffer distance in meters (default: 1000)'
    )

    parser.add_argument(
        '--years',
        type=int,
        nargs=2,
        metavar=('START', 'END'),
        help='Year range filter (e.g., --years 2010 2023)'
    )

    parser.add_argument(
        '--cache',
        type=str,
        help='Directory to cache GRQA downloads'
    )

    parser.add_argument(
        '--prefix',
        type=str,
        default='calibration',
        help='Output file prefix (default: calibration)'
    )

    parser.add_argument(
        '--create-mapping',
        type=str,
        metavar='FILE',
        help='Create example species mapping file and exit'
    )

    parser.add_argument(
        '--list-templates',
        action='store_true',
        help='List available mapping templates and exit'
    )

    parser.add_argument(
        '--list-parameters',
        action='store_true',
        help='List all GRQA parameters and exit'
    )

    args = parser.parse_args()

    # Handle info commands
    if args.list_templates:
        SpeciesMapper.list_templates()
        return

    if args.list_parameters:
        print("\nGRQA Parameters (43 total):")
        print("=" * 70)
        for group in ['nitrogen', 'phosphorus', 'oxygen', 'carbon', 'sediment', 'physical', 'biological', 'other']:
            params = [k for k, v in GRQA_PARAMETERS.items() if v.get('group') == group]
            if params:
                print(f"\n{group.upper()}:")
                for p in params:
                    info = GRQA_PARAMETERS[p]
                    print(f"  {p:12s} - {info['name']:40s} [{info['unit']}]")
        return

    if args.create_mapping:
        create_example_mapping_file(args.create_mapping)
        return

    # Validate required arguments
    if not args.shapefile:
        parser.error("shapefile is required (unless using --create-mapping, --list-templates, or --list-parameters)")

    if not os.path.exists(args.shapefile):
        print(f"Error: Shapefile not found: {args.shapefile}")
        sys.exit(1)

    # Set up species mapping
    if args.mapping:
        if not os.path.exists(args.mapping):
            print(f"Error: Mapping file not found: {args.mapping}")
            sys.exit(1)
        mapper = SpeciesMapper(mapping_file=args.mapping)
    elif args.template:
        mapper = SpeciesMapper(mapping=DEFAULT_SPECIES_MAPPINGS[args.template])
        print(f"Using template: {args.template}")
    else:
        print("Error: Either --mapping or --template is required")
        print("Use --create-mapping to create an example mapping file")
        print("Use --list-templates to see available templates")
        sys.exit(1)

    # Run extraction
    extractor = GRQACalibrationExtractor(
        output_dir=args.output,
        species_mapper=mapper,
        cache_dir=args.cache,
        buffer_distance_m=args.buffer,
        years=tuple(args.years) if args.years else None
    )

    extractor.run(args.shapefile, prefix=args.prefix)


if __name__ == "__main__":
    main()
