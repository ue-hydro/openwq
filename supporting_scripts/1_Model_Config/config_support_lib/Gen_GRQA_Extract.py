#!/usr/bin/env python3
"""
Gen_GRQA_Extract.py — Standalone GRQA Station Extractor for OpenWQ Config Reports

Extracts water quality monitoring stations and observations from the
Global River Water Quality Archive (GRQA) dataset based on a river network
or basin shapefile.  Produces clipped CSV files that are used by the
config report (Gen_Report.py) for the observation data overlay and by the
calibration workflow.

This module is self-contained — it does NOT depend on calibration_lib or
any other external OpenWQ package.  Required third-party packages:
    pandas, geopandas, shapely, numpy, requests, pyproj

GRQA Dataset: https://zenodo.org/records/15335450
Reference: Virro et al. (2021) - GRQA: Global River Water Quality Archive
           Earth System Science Data, 13, 5483-5507

Author: OpenWQ Development Team
"""

import os
import sys
import json
import zipfile
import requests
from pathlib import Path
from typing import List, Dict, Optional, Tuple
from datetime import datetime
import warnings

warnings.filterwarnings('ignore')

try:
    import pandas as pd
    import geopandas as gpd
    from shapely.geometry import Point, box
    from shapely.ops import unary_union
    import numpy as np
    _DEPS_AVAILABLE = True
except ImportError:
    _DEPS_AVAILABLE = False

# ---------------------------------------------------------------------------
# GRQA Dataset Constants
# ---------------------------------------------------------------------------

GRQA_ZENODO_RECORD = "15335450"
GRQA_ZENODO_BASE = f"https://zenodo.org/records/{GRQA_ZENODO_RECORD}"
GRQA_DATA_URL = f"{GRQA_ZENODO_BASE}/files/GRQA_data_v1.4.zip"

# GRQA v1.4 filename mapping: parameter code -> filename stem in the zip
GRQA_V14_FILENAME_MAP = {
    'NO3': 'NO3N',
    'NH4': 'NH4N',
    'NO2': 'NO2N',
    'NH3': 'NH3N',
    'DO_sat': 'DOSAT',
    'Temp': 'TEMP',
    'PO4': 'DIP',
    'DP': 'TDP',
    'PP': 'TPP',
    'NO3_NO2': 'TIN',
    'SS': 'TSS',
    'PN': 'PN',
    'TKN': 'TKN',
    'DIN': 'DIN',
}

# 43 GRQA parameters with human-readable names, units, and group
GRQA_PARAMETERS = {
    'TN':  {'name': 'Total Nitrogen',               'unit': 'mg/L',       'group': 'nitrogen'},
    'DN':  {'name': 'Dissolved Nitrogen',            'unit': 'mg/L',       'group': 'nitrogen'},
    'NO3': {'name': 'Nitrate',                       'unit': 'mg/L',       'group': 'nitrogen'},
    'NO2': {'name': 'Nitrite',                       'unit': 'mg/L',       'group': 'nitrogen'},
    'NH4': {'name': 'Ammonium',                      'unit': 'mg/L',       'group': 'nitrogen'},
    'NO3_NO2': {'name': 'Nitrate + Nitrite',         'unit': 'mg/L',       'group': 'nitrogen'},
    'DIN': {'name': 'Dissolved Inorganic Nitrogen',  'unit': 'mg/L',       'group': 'nitrogen'},
    'DON': {'name': 'Dissolved Organic Nitrogen',    'unit': 'mg/L',       'group': 'nitrogen'},
    'TKN': {'name': 'Total Kjeldahl Nitrogen',       'unit': 'mg/L',       'group': 'nitrogen'},
    'PN':  {'name': 'Particulate Nitrogen',          'unit': 'mg/L',       'group': 'nitrogen'},
    'TP':  {'name': 'Total Phosphorus',              'unit': 'mg/L',       'group': 'phosphorus'},
    'DP':  {'name': 'Dissolved Phosphorus',          'unit': 'mg/L',       'group': 'phosphorus'},
    'DIP': {'name': 'Dissolved Inorganic Phosphorus','unit': 'mg/L',       'group': 'phosphorus'},
    'PO4': {'name': 'Phosphate',                     'unit': 'mg/L',       'group': 'phosphorus'},
    'DOP': {'name': 'Dissolved Organic Phosphorus',  'unit': 'mg/L',       'group': 'phosphorus'},
    'PP':  {'name': 'Particulate Phosphorus',        'unit': 'mg/L',       'group': 'phosphorus'},
    'DO':  {'name': 'Dissolved Oxygen',              'unit': 'mg/L',       'group': 'oxygen'},
    'DO_sat': {'name': 'Dissolved Oxygen Saturation','unit': '%',          'group': 'oxygen'},
    'BOD': {'name': 'Biochemical Oxygen Demand',     'unit': 'mg/L',       'group': 'oxygen'},
    'BOD5':{'name': 'Biochemical Oxygen Demand (5-day)','unit': 'mg/L',    'group': 'oxygen'},
    'COD': {'name': 'Chemical Oxygen Demand',        'unit': 'mg/L',       'group': 'oxygen'},
    'CODMn':{'name': 'Chemical Oxygen Demand (Permanganate)','unit':'mg/L','group': 'oxygen'},
    'DOC': {'name': 'Dissolved Organic Carbon',      'unit': 'mg/L',       'group': 'carbon'},
    'TOC': {'name': 'Total Organic Carbon',          'unit': 'mg/L',       'group': 'carbon'},
    'TC':  {'name': 'Total Carbon',                  'unit': 'mg/L',       'group': 'carbon'},
    'TIC': {'name': 'Total Inorganic Carbon',        'unit': 'mg/L',       'group': 'carbon'},
    'DIC': {'name': 'Dissolved Inorganic Carbon',    'unit': 'mg/L',       'group': 'carbon'},
    'POC': {'name': 'Particulate Organic Carbon',    'unit': 'mg/L',       'group': 'carbon'},
    'PC':  {'name': 'Particulate Carbon',            'unit': 'mg/L',       'group': 'carbon'},
    'TSS': {'name': 'Total Suspended Solids',        'unit': 'mg/L',       'group': 'sediment'},
    'TDS': {'name': 'Total Dissolved Solids',        'unit': 'mg/L',       'group': 'sediment'},
    'SS':  {'name': 'Suspended Solids',              'unit': 'mg/L',       'group': 'sediment'},
    'Turbidity': {'name': 'Turbidity',               'unit': 'NTU',        'group': 'sediment'},
    'Temp':{'name': 'Water Temperature',             'unit': '\u00b0C',    'group': 'physical'},
    'pH':  {'name': 'pH',                            'unit': 'pH',         'group': 'physical'},
    'EC':  {'name': 'Electrical Conductivity',       'unit': '\u00b5S/cm', 'group': 'physical'},
    'SC':  {'name': 'Specific Conductance',          'unit': '\u00b5S/cm', 'group': 'physical'},
    'Chl_a':{'name': 'Chlorophyll-a',                'unit': '\u00b5g/L',  'group': 'biological'},
    'Secchi':{'name': 'Secchi Depth',                'unit': 'm',          'group': 'physical'},
    'Si':  {'name': 'Silica',                        'unit': 'mg/L',       'group': 'other'},
    'SO4': {'name': 'Sulfate',                       'unit': 'mg/L',       'group': 'other'},
    'Cl':  {'name': 'Chloride',                      'unit': 'mg/L',       'group': 'other'},
    'Alk': {'name': 'Alkalinity',                    'unit': 'mg/L CaCO3', 'group': 'other'},
}

# Default species mapping templates
DEFAULT_SPECIES_MAPPINGS = {
    'nitrogen_cycle': {
        'NO3': 'NO3', 'NH4': 'NH4', 'NO2': 'NO2',
        'TN': 'TN', 'DON': 'DON',
    },
    'phosphorus_cycle': {
        'PO4': 'PO4', 'TP': 'TP', 'DP': 'DP',
    },
    'carbon_cycle': {
        'DOC': 'DOC', 'DIC': 'DIC', 'TOC': 'TOC',
    },
    'oxygen': {
        'DO': 'DO', 'BOD5': 'BOD',
    },
    'full_nutrient': {
        'NO3': 'NO3', 'NH4': 'NH4', 'PO4': 'PO4',
        'DO': 'DO', 'DOC': 'DOC', 'Temp': 'TEMP',
    },
}


# ---------------------------------------------------------------------------
# SpeciesMapper
# ---------------------------------------------------------------------------

class SpeciesMapper:
    """Maps GRQA parameter codes to model species names and vice versa."""

    def __init__(self, mapping: Dict[str, str] = None, mapping_file: str = None):
        self.mapping: Dict[str, str] = {}
        self.reverse_mapping: Dict[str, str] = {}
        if mapping_file:
            self.load_mapping(mapping_file)
        elif mapping:
            self.set_mapping(mapping)

    def set_mapping(self, mapping: Dict[str, str]):
        self.mapping = mapping
        self.reverse_mapping = {v: k for k, v in mapping.items()}

    def load_mapping(self, filepath: str):
        with open(filepath, 'r') as f:
            data = json.load(f)
        if 'species_mapping' in data:
            self.mapping = data['species_mapping']
        else:
            self.mapping = data
        self.reverse_mapping = {v: k for k, v in self.mapping.items()}

    def save_mapping(self, filepath: str, include_metadata: bool = True):
        if include_metadata:
            data = {
                'description': 'GRQA to Model Species Mapping',
                'created': datetime.now().isoformat(),
                'grqa_source': GRQA_ZENODO_BASE,
                'species_mapping': self.mapping,
                'grqa_parameters_info': {
                    k: GRQA_PARAMETERS[k] for k in self.mapping
                    if k in GRQA_PARAMETERS
                },
            }
        else:
            data = self.mapping
        with open(filepath, 'w') as f:
            json.dump(data, f, indent=2)

    def get_grqa_parameters(self) -> List[str]:
        return list(self.mapping.keys())

    def get_model_name(self, grqa_param: str) -> str:
        return self.mapping.get(grqa_param, grqa_param)

    def get_grqa_name(self, model_species: str) -> str:
        return self.reverse_mapping.get(model_species, model_species)

    def validate(self) -> Tuple[bool, List[str]]:
        invalid = [p for p in self.mapping if p not in GRQA_PARAMETERS]
        return len(invalid) == 0, invalid

    @staticmethod
    def create_template(template_name: str = None) -> Dict[str, str]:
        if template_name and template_name in DEFAULT_SPECIES_MAPPINGS:
            return DEFAULT_SPECIES_MAPPINGS[template_name].copy()
        return {k: k for k in GRQA_PARAMETERS}

    @staticmethod
    def list_templates():
        print("\nAvailable Species Mapping Templates:")
        print("=" * 50)
        for name, mapping in DEFAULT_SPECIES_MAPPINGS.items():
            print(f"\n{name}:")
            for grqa, model in mapping.items():
                info = GRQA_PARAMETERS.get(grqa, {})
                print(f"  {grqa:12s} -> {model:12s} ({info.get('name', 'N/A')})")


# ---------------------------------------------------------------------------
# GRQADownloader
# ---------------------------------------------------------------------------

class GRQADownloader:
    """Handles selective downloading / locating of GRQA CSV files."""

    def __init__(self, cache_dir: str, local_data_path: str = None):
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.local_data_path = Path(local_data_path) if local_data_path else None
        self.use_local = self.local_data_path is not None
        self.local_file_patterns: list = []
        self.local_version: Optional[str] = None

        if self.use_local:
            if not self.local_data_path.exists():
                raise FileNotFoundError(
                    f"Local GRQA data path not found: {self.local_data_path}")
            print(f"Using local GRQA data from: {self.local_data_path}")
            self._detect_local_structure()

    # ── internal ──

    def _detect_local_structure(self):
        self.local_file_patterns = []
        self.local_version = None
        for pattern in ['GRQA_data_v*', 'grqa_data_v*']:
            for path in self.local_data_path.glob(pattern):
                if path.is_dir():
                    version = path.name.split('_')[-1]
                    self.local_version = version
                    print(f"  Detected GRQA version: {version}")
                    csv_files = list(path.glob('*.csv'))
                    if csv_files:
                        sample = csv_files[0].name
                        if '_GRQA.csv' in sample:
                            self.local_file_patterns.append((path, '{param}_GRQA.csv'))
                        elif 'GRQA_data_' in sample:
                            self.local_file_patterns.append(
                                (path, f'GRQA_data_{version}_{{param}}.csv'))
                        print(f"  Found {len(csv_files)} parameter files")
                    break
        direct_csvs = list(self.local_data_path.glob('*_GRQA.csv'))
        if direct_csvs:
            self.local_file_patterns.append(
                (self.local_data_path, '{param}_GRQA.csv'))
            print(f"  Found {len(direct_csvs)} parameter files in root")
        if not self.local_file_patterns:
            print("  Warning: No GRQA CSV files found. Will attempt download.")
            self.use_local = False

    def _find_local_file(self, parameter: str) -> Optional[Path]:
        if not self.use_local:
            return None
        candidates = [parameter]
        if parameter in GRQA_V14_FILENAME_MAP:
            candidates.append(GRQA_V14_FILENAME_MAP[parameter])
        for base_path, pat in self.local_file_patterns:
            for cand in candidates:
                fp = base_path / pat.format(param=cand)
                if fp.exists():
                    return fp
            for cand in candidates:
                for f in base_path.glob(f'*{cand}*'):
                    if f.suffix.lower() == '.csv':
                        return f
        return None

    def _download_with_progress(self, url: str, dest: Path):
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
                        print(f"\r  Downloaded: {mb:.1f} MB ({pct:.1f}%)",
                              end="", flush=True)
        print()

    def _extract_parameter(self, zip_path: Path, parameter: str) -> Optional[Path]:
        candidates = [parameter]
        if parameter in GRQA_V14_FILENAME_MAP:
            candidates.append(GRQA_V14_FILENAME_MAP[parameter])
        with zipfile.ZipFile(zip_path, 'r') as zf:
            for cand in candidates:
                target = f"{cand}_GRQA.csv"
                for name in zf.namelist():
                    if name.endswith(target):
                        zf.extract(name, self.cache_dir)
                        extracted = self.cache_dir / name
                        dest = self.cache_dir / f"{parameter}_GRQA.csv"
                        if extracted != dest:
                            if dest.exists():
                                dest.unlink()
                            extracted.rename(dest)
                        if cand != parameter:
                            print(f"  Mapped {parameter} -> {cand} (GRQA v1.4 naming)")
                        return dest
        print(f"  Warning: Parameter {parameter} not found in archive "
              f"(tried: {', '.join(c + '_GRQA.csv' for c in candidates)})")
        return None

    # ── public ──

    def download_parameter_file(self, parameter: str,
                                force: bool = False) -> Optional[Path]:
        if self.use_local:
            local_file = self._find_local_file(parameter)
            if local_file:
                return local_file
            else:
                print(f"  Warning: {parameter} not found in local data, will try download")
        param_file = self.cache_dir / f"{parameter}_GRQA.csv"
        if param_file.exists() and not force:
            return param_file
        zip_path = self.cache_dir / "GRQA_data_v1.4.zip"
        need_download = not zip_path.exists()
        if zip_path.exists():
            try:
                with zipfile.ZipFile(zip_path, 'r') as zf:
                    zf.namelist()
            except (zipfile.BadZipFile, Exception):
                print("  Cached zip is corrupted — re-downloading...")
                zip_path.unlink()
                need_download = True
        if need_download:
            print("Downloading GRQA data archive...")
            print("  This is a one-time download (~1.2 GB)")
            self._download_with_progress(GRQA_DATA_URL, zip_path)
        return self._extract_parameter(zip_path, parameter)

    def download_required_parameters(self, parameters: List[str],
                                      force: bool = False) -> Dict[str, Path]:
        source_type = "local" if self.use_local else "remote"
        print(f"\nPreparing GRQA data for {len(parameters)} parameters ({source_type})...")
        available: Dict[str, Path] = {}
        missing: List[str] = []
        for param in parameters:
            if self.use_local:
                local_file = self._find_local_file(param)
                if local_file:
                    available[param] = local_file
                    print(f"  {param}: local")
                    continue
            param_file = self.cache_dir / f"{param}_GRQA.csv"
            if param_file.exists() and not force:
                available[param] = param_file
                print(f"  {param}: cached")
            else:
                missing.append(param)
        if missing:
            print(f"\nDownloading {len(missing)} parameter files...")
            for param in missing:
                filepath = self.download_parameter_file(param, force)
                if filepath:
                    available[param] = filepath
                    print(f"  {param}: downloaded")
        return available


# ---------------------------------------------------------------------------
# GRQACalibrationExtractor
# ---------------------------------------------------------------------------

class GRQACalibrationExtractor:
    """Extract GRQA water quality data for config reports and calibration."""

    def __init__(self,
                 output_dir: str,
                 species_mapper: SpeciesMapper,
                 cache_dir: str = None,
                 local_data_path: str = None,
                 buffer_distance_m: float = 1000,
                 years: Optional[Tuple[int, int]] = None):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.species_mapper = species_mapper
        self.buffer_distance_m = buffer_distance_m
        self.years = years
        self.local_data_path = local_data_path
        self.cache_dir = Path(cache_dir) if cache_dir else self.output_dir / "grqa_cache"
        self.downloader = GRQADownloader(self.cache_dir, local_data_path=local_data_path)
        self.column_map: Optional[Dict[str, str]] = None
        self.stations_gdf = None
        self.observations_df = None

    # ── CRS helpers ──

    @staticmethod
    def _get_wgs84_crs():
        WGS84_WKT = (
            'GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",'
            'SPHEROID["WGS_1984",6378137.0,298.257223563]],'
            'PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]]'
        )
        try:
            from pyproj import CRS
            return CRS.from_wkt(WGS84_WKT)
        except Exception:
            try:
                return {'init': 'epsg:4326'}
            except Exception:
                return None

    @staticmethod
    def _get_utm_crs(epsg_code: int):
        if epsg_code >= 32700:
            zone = epsg_code - 32700
            south = ' +south'
        else:
            zone = epsg_code - 32600
            south = ''
        proj4 = f'+proj=utm +zone={zone}{south} +datum=WGS84 +units=m +no_defs'
        try:
            from pyproj import CRS
            return CRS.from_proj4(proj4)
        except Exception:
            return proj4

    # ── shapefile loading ──

    def load_river_network(self, shapefile_path: str):
        """Load and validate river network shapefile, return GeoDataFrame."""
        print(f"\nLoading river network: {shapefile_path}")
        wgs84 = self._get_wgs84_crs()
        try:
            gdf = gpd.read_file(shapefile_path)
        except Exception:
            import fiona
            from shapely.geometry import shape as shp_shape
            with fiona.open(shapefile_path) as src:
                records = []
                for feat in src:
                    props = dict(feat['properties'])
                    props['geometry'] = shp_shape(feat['geometry'])
                    records.append(props)
            gdf = gpd.GeoDataFrame(records)

        if gdf.crs is None and wgs84 is not None:
            bounds = gdf.total_bounds
            if abs(bounds[0]) > 360 or abs(bounds[1]) > 360:
                print("  Warning: No CRS metadata but coordinates appear projected.")
            gdf.crs = wgs84
        elif gdf.crs is not None:
            try:
                epsg = gdf.crs.to_epsg()
                if epsg is not None and epsg != 4326:
                    print(f"  Reprojecting from EPSG:{epsg} to WGS84")
                    gdf = gdf.to_crs("EPSG:4326")
                elif epsg is None:
                    bounds = gdf.total_bounds
                    if abs(bounds[0]) > 360 or abs(bounds[1]) > 360:
                        gdf = gdf.to_crs("EPSG:4326")
            except Exception:
                bounds = gdf.total_bounds
                if abs(bounds[0]) > 360 or abs(bounds[1]) > 360:
                    try:
                        gdf = gdf.to_crs("EPSG:4326")
                    except Exception:
                        pass

        print(f"  Features: {len(gdf)}")
        print(f"  Bounds: {gdf.total_bounds}")
        return gdf

    # ── buffer creation ──

    def create_buffer(self, river_gdf):
        """Create buffer around river network for station selection."""
        from shapely.ops import transform as shp_transform
        from pyproj import Transformer

        print(f"\nCreating {self.buffer_distance_m}m buffer...")
        wgs84 = self._get_wgs84_crs()
        union_geom = unary_union(river_gdf.geometry)
        centroid = union_geom.centroid
        utm_zone = int((centroid.x + 180) / 6) + 1
        hemisphere = 'north' if centroid.y >= 0 else 'south'
        epsg_code = 32600 + utm_zone if hemisphere == 'north' else 32700 + utm_zone
        utm_crs = self._get_utm_crs(epsg_code)

        fwd = Transformer.from_crs(wgs84, utm_crs, always_xy=True)
        rev = Transformer.from_crs(utm_crs, wgs84, always_xy=True)

        projected = shp_transform(fwd.transform, union_geom)
        buffered = projected.buffer(self.buffer_distance_m)
        buffer_wgs84 = shp_transform(rev.transform, buffered)

        buffer_gdf = gpd.GeoDataFrame(geometry=[buffer_wgs84], crs=wgs84)
        print(f"  Using UTM zone {utm_zone}{hemisphere[0].upper()}")
        return buffer_gdf

    # ── column detection ──

    def _detect_columns(self, df) -> Dict[str, str]:
        columns = [c.lower() for c in df.columns]
        original = df.columns.tolist()
        patterns = {
            'site_id':   ['site_id', 'site_no', 'station_id'],
            'site_name': ['site_name', 'station_name'],
            'lat':       ['lat_wgs84', 'lat', 'latitude'],
            'lon':       ['lon_wgs84', 'lon', 'longitude'],
            'obs_date':  ['obs_date', 'date', 'sample_date'],
            'obs_time':  ['obs_time', 'time'],
            'obs_value': ['obs_value', 'value', 'result'],
            'unit':      ['unit', 'units'],
            'source':    ['source', 'data_source'],
            'outlier':   ['obs_iqr_outlier', 'outlier', 'is_outlier'],
        }
        col_map: Dict[str, str] = {}
        for field, cands in patterns.items():
            for i, col in enumerate(columns):
                if col in cands:
                    col_map[field] = original[i]
                    break
        return col_map

    # ── core extraction ──

    def extract_stations_and_observations(self, buffer_gdf):
        """Extract stations and observations for mapped species.

        Returns (stations_gdf, observations_df).
        """
        parameters = self.species_mapper.get_grqa_parameters()
        print(f"\nExtracting data for {len(parameters)} parameters:")
        for p in parameters:
            model_name = self.species_mapper.get_model_name(p)
            info = GRQA_PARAMETERS.get(p, {})
            print(f"  {p} -> {model_name} ({info.get('name', 'N/A')})")

        param_files = self.downloader.download_required_parameters(parameters)
        if not param_files:
            print("Error: No parameter files available!")
            return gpd.GeoDataFrame(), pd.DataFrame()

        all_stations = []
        all_observations = []
        buffer_union = buffer_gdf.geometry.unary_union

        for param, filepath in param_files.items():
            print(f"\nProcessing {param}...")
            try:
                try:
                    df = pd.read_csv(filepath, sep=';', low_memory=False)
                    if len(df.columns) < 3:
                        df = pd.read_csv(filepath, sep=',', low_memory=False)
                except Exception:
                    df = pd.read_csv(filepath, sep=',', low_memory=False)

                if self.column_map is None:
                    self.column_map = self._detect_columns(df)
                    print(f"  Detected columns: {list(self.column_map.keys())}")

                lat_col = self.column_map.get('lat', 'lat_wgs84')
                lon_col = self.column_map.get('lon', 'lon_wgs84')
                site_col = self.column_map.get('site_id', 'site_id')
                date_col = self.column_map.get('obs_date', 'obs_date')

                geometry = [Point(xy) for xy in zip(df[lon_col], df[lat_col])]
                gdf = gpd.GeoDataFrame(df, geometry=geometry,
                                       crs=self._get_wgs84_crs())
                gdf_filtered = gdf[gdf.geometry.within(buffer_union)].copy()

                if len(gdf_filtered) == 0:
                    print(f"  No observations within buffer")
                    continue

                if self.years:
                    gdf_filtered[date_col] = pd.to_datetime(
                        gdf_filtered[date_col], errors='coerce')
                    gdf_filtered = gdf_filtered[
                        (gdf_filtered[date_col].dt.year >= self.years[0]) &
                        (gdf_filtered[date_col].dt.year <= self.years[1])
                    ]

                model_name = self.species_mapper.get_model_name(param)
                gdf_filtered['grqa_parameter'] = param
                gdf_filtered['model_species'] = model_name

                station_cols = [site_col, lat_col, lon_col, 'geometry']
                if self.column_map.get('site_name') in gdf_filtered.columns:
                    station_cols.insert(1, self.column_map['site_name'])

                stations = gdf_filtered[station_cols].drop_duplicates(
                    subset=[site_col])
                stations['has_' + model_name] = True
                all_stations.append(stations)
                all_observations.append(gdf_filtered)

                print(f"  Found {len(gdf_filtered)} observations "
                      f"at {len(stations)} stations")
            except Exception as e:
                print(f"  Error: {e}")
                continue

        if not all_stations:
            print("\nNo data found within the buffer area!")
            return gpd.GeoDataFrame(), pd.DataFrame()

        # Combine stations
        site_col = self.column_map.get('site_id', 'site_id')
        for i, s in enumerate(all_stations):
            if site_col in s.columns:
                all_stations[i][site_col] = s[site_col].astype(str)

        stations_combined = all_stations[0]
        for s in all_stations[1:]:
            stations_combined = stations_combined.merge(
                s.drop(columns='geometry'),
                on=site_col, how='outer', suffixes=('', '_dup'))
            dup_cols = [c for c in stations_combined.columns
                        if c.endswith('_dup')]
            for col in dup_cols:
                base = col.removesuffix('_dup')
                if base in stations_combined.columns:
                    stations_combined[base] = stations_combined[base].fillna(
                        stations_combined[col])
            stations_combined.drop(columns=dup_cols, inplace=True)

        has_cols = [c for c in stations_combined.columns if c.startswith('has_')]
        stations_combined[has_cols] = stations_combined[has_cols].fillna(False)

        lat_col = self.column_map.get('lat', 'lat_wgs84')
        lon_col = self.column_map.get('lon', 'lon_wgs84')
        stations_combined = gpd.GeoDataFrame(
            stations_combined.drop(columns='geometry'),
            geometry=[Point(xy) for xy in zip(
                stations_combined[lon_col], stations_combined[lat_col])],
            crs=self._get_wgs84_crs())

        observations_combined = pd.concat(all_observations, ignore_index=True)
        observations_df = pd.DataFrame(
            observations_combined.drop(columns='geometry'))

        self.stations_gdf = stations_combined
        self.observations_df = observations_df

        print(f"\n  Total stations: {len(stations_combined)}")
        print(f"  Total observations: {len(observations_df)}")
        return stations_combined, observations_df

    # ── file generation ──

    def generate_calibration_files(self, prefix: str = "calibration"):
        """Generate calibration-ready output files (stations, obs, timeseries)."""
        if self.stations_gdf is None or self.observations_df is None:
            print("Error: No data to save. Run extraction first.")
            return

        print(f"\nGenerating calibration files in: {self.output_dir}")

        stations_shp = self.output_dir / f"{prefix}_stations.shp"
        stations_csv = self.output_dir / f"{prefix}_stations.csv"
        self.stations_gdf.to_file(stations_shp)
        self.stations_gdf.drop(columns='geometry').to_csv(stations_csv,
                                                           index=False)
        print(f"  Stations: {stations_shp}")

        obs_csv = self.output_dir / f"{prefix}_observations.csv"
        self.observations_df.to_csv(obs_csv, index=False)
        print(f"  Observations: {obs_csv}")

        timeseries_dir = self.output_dir / f"{prefix}_timeseries"
        timeseries_dir.mkdir(exist_ok=True)
        site_col = self.column_map.get('site_id', 'site_id')
        date_col = self.column_map.get('obs_date', 'obs_date')
        value_col = self.column_map.get('obs_value', 'obs_value')

        for model_species in self.observations_df['model_species'].unique():
            species_df = self.observations_df[
                self.observations_df['model_species'] == model_species].copy()
            if date_col in species_df.columns and site_col in species_df.columns:
                pivot_df = species_df.pivot_table(
                    index=date_col, columns=site_col,
                    values=value_col, aggfunc='mean')
                ts_file = timeseries_dir / f"{model_species}_timeseries.csv"
                pivot_df.to_csv(ts_file)
                print(f"  Time series: {ts_file}")

        config = self._generate_calibration_config(prefix)
        config_file = self.output_dir / f"{prefix}_config.json"
        with open(config_file, 'w') as f:
            json.dump(config, f, indent=2)
        print(f"  Config: {config_file}")

        self._generate_summary_report(prefix)

    def _generate_calibration_config(self, prefix: str) -> Dict:
        site_col = self.column_map.get('site_id', 'site_id')
        date_col = self.column_map.get('obs_date', 'obs_date')
        dates = pd.to_datetime(self.observations_df[date_col], errors='coerce')
        config: Dict = {
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
            'calibration_targets': [],
        }
        for species in self.observations_df['model_species'].unique():
            sp_data = self.observations_df[
                self.observations_df['model_species'] == species]
            value_col = self.column_map.get('obs_value', 'obs_value')
            config['calibration_targets'].append({
                'model_species': species,
                'grqa_parameter': self.species_mapper.get_grqa_name(species),
                'observation_count': len(sp_data),
                'station_count': sp_data[site_col].nunique(),
                'value_stats': {
                    'min': float(sp_data[value_col].min()),
                    'max': float(sp_data[value_col].max()),
                    'mean': float(sp_data[value_col].mean()),
                    'median': float(sp_data[value_col].median()),
                },
                'timeseries_file': f"{prefix}_timeseries/{species}_timeseries.csv",
            })
        return config

    def _generate_summary_report(self, prefix: str):
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
            f.write("\n" + "-" * 60 + "\nSPECIES MAPPING\n" + "-" * 60 + "\n")
            for grqa, model in self.species_mapper.mapping.items():
                info = GRQA_PARAMETERS.get(grqa, {})
                f.write(f"  {grqa:12s} -> {model:12s} "
                        f"({info.get('name', 'N/A')})\n")
            f.write("\n" + "-" * 60 + "\nDATA SUMMARY\n" + "-" * 60 + "\n")
            f.write(f"  Total Stations: {len(self.stations_gdf)}\n")
            f.write(f"  Total Observations: {len(self.observations_df)}\n\n")
            f.write("  Per-Species Summary:\n")
            site_col = self.column_map.get('site_id', 'site_id')
            for species in sorted(self.observations_df['model_species'].unique()):
                sp = self.observations_df[
                    self.observations_df['model_species'] == species]
                f.write(f"    {species:12s}: {len(sp):8,} obs at "
                        f"{sp[site_col].nunique():4} stations\n")
            f.write("\n" + "-" * 60 + "\nOUTPUT FILES\n" + "-" * 60 + "\n")
            f.write(f"  {prefix}_stations.shp      - Station locations\n")
            f.write(f"  {prefix}_stations.csv      - Station locations (CSV)\n")
            f.write(f"  {prefix}_observations.csv  - All observations\n")
            f.write(f"  {prefix}_timeseries/       - Per-species time series\n")
            f.write(f"  {prefix}_config.json       - Calibration configuration\n\n")
        print(f"  Summary: {report_file}")

    def run(self, shapefile_path: str, prefix: str = "calibration"):
        """Run the complete extraction pipeline."""
        print("\n" + "=" * 60)
        print("GRQA Calibration Data Extractor")
        print("=" * 60)
        valid, invalid = self.species_mapper.validate()
        if not valid:
            print(f"\nWarning: Unknown GRQA parameters: {invalid}")
        river_gdf = self.load_river_network(shapefile_path)
        buffer_gdf = self.create_buffer(river_gdf)
        stations, observations = self.extract_stations_and_observations(buffer_gdf)
        if stations.empty:
            print("\nNo data found. Check your inputs.")
            return None, None
        self.generate_calibration_files(prefix)
        print("\n" + "=" * 60)
        print("Extraction Complete!")
        print("=" * 60 + "\n")
        return stations, observations
