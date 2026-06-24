# Copyright 2026, Diogo Costa, diogo.costa@uevora.pt
# This file is part of OpenWQ model.

# This program, openWQ, is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# !/usr/bin/env python3
"""
OpenWQ Source/Sink Configuration Generator
Generate JSON configuration for chemical sources and sinks

PATCHED VERSION: Includes comprehensive workaround for Python 3.9 tarfile compatibility issues
"""

# Comprehensive workaround for tarfile compatibility issues in Python 3.9
import os
import sys
import tarfile

# Add data_filter for Python < 3.12 compatibility
if not hasattr(tarfile, 'data_filter'):
    # Create a dummy data_filter that just returns the tarinfo unchanged
    def data_filter(member, dest_path):
        return member


    tarfile.data_filter = data_filter

# Monkey patch backports.tarfile if it tries to import
if sys.version_info >= (3, 9):
    try:
        sys.modules['backports.tarfile'] = tarfile
    except:
        pass

import json
import calendar
from pathlib import Path
from typing import List, Dict, Union, Optional
import geopandas as gpd
import rasterio
from rasterio.transform import from_bounds
from rasterio.mask import mask as rasterio_mask
import xarray as xr
import numpy as np
import pandas as pd
from collections import defaultdict
import re
import warnings

warnings.filterwarnings('ignore')


###################
# Set SS from user prescribed data in csv
##################
def set_ss_from_csv(

        ss_config_filepath: str,

        json_header_comment: List[str],

        # Metadata
        ss_metadata_source: str,
        ss_metadata_comment: str,

        # List of source/sink configurations
        ss_method_csv_config: List[Dict[str, Union[str, int]]],

        # BGC engine info (for METADATA traceability)
        bgc_engine_label: str = "",
        chemical_species_list: Optional[List[str]] = None

) -> None:
    """
    Create OpenWQ Source/Sink configuration JSON file.

    Parameters:
        ss_config_filepath: Full path where the JSON file will be saved
        json_header_comment: List of comment lines to add at the top of the file
        ss_metadata_comment: Comment for metadata section
        ss_metadata_source: Source identifier for metadata section
        ss_method_csv_config: List of dictionaries, each containing:
            - Chemical_name: Name of the chemical species
            - Compartment_name: Name of the compartment
            - Type: "source" or "sink"
            - Units: Units (e.g., "kg")
            - Filepath: Path to CSV data file (header row auto-detected by 'YYYY' column)
            - Delimiter: CSV delimiter (e.g., ",")
    """

    # Create the destination directory path if it doesn't exist
    output_path = Path(ss_config_filepath)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # ── Validate CSV configurations ──────────────────────────────────
    # Required columns the C++ ASCII parser looks up by header name.
    # If any are missing the parser accesses index -1 → crash (SIGSEGV).
    _REQUIRED_COLS = {"YYYY", "MM", "DD", "HH", "MIN", "SEC",
                      "IX", "IY", "IZ", "LOAD", "LOAD_TYPE"}

    for idx, ss_config in enumerate(ss_method_csv_config, start=1):
        filepath = ss_config.get("Filepath", "")
        delimiter = ss_config.get("Delimiter", ",")
        label = ss_config.get("Chemical_name", "?")

        # Check file exists
        if not os.path.isfile(filepath):
            print(f"  ⚠  WARNING: SS entry {idx} ({label}) "
                  f"CSV file not found: {filepath}")
            print(f"     The model will skip this entry at runtime.")
            continue

        # Read file and validate header + data
        try:
            with open(filepath, 'r') as _f:
                all_lines = _f.readlines()

            # Find the header row (first line containing "YYYY")
            header_line = None
            header_lineno = None
            data_lines = []
            for li, line in enumerate(all_lines, start=1):
                if header_line is None and "YYYY" in line.upper():
                    header_line = line.strip().upper()
                    header_lineno = li
                elif header_line is not None:
                    stripped = line.strip()
                    if stripped:
                        data_lines.append((li, stripped))

            if header_line is None:
                print(f"  ⚠  WARNING: SS entry {idx} ({label}) CSV has no "
                      f"header row containing 'YYYY': {filepath}")
                print(f"     The C++ parser requires a header row with columns:")
                print(f"     YYYY,MM,DD,HH,MIN,SEC,IX,IY,IZ,LOAD,LOAD_TYPE")
                continue

            # Parse header columns
            header_cols = [c.strip() for c in header_line.split(delimiter)]
            header_set = set(header_cols)

            # Check for missing required columns
            missing = _REQUIRED_COLS - header_set
            if missing:
                print(f"  ⚠  WARNING: SS entry {idx} ({label}) CSV is missing "
                      f"required columns: {sorted(missing)}")
                print(f"     File: {filepath}")
                print(f"     Header found (line {header_lineno}): {header_line}")
                print(f"     Required columns: "
                      f"YYYY,MM,DD,HH,MIN,SEC,IX,IY,IZ,LOAD,LOAD_TYPE")
                print(f"     ─── This WILL crash the model (SIGSEGV). "
                      f"Please fix the CSV before running. ───")

            # Check data rows have the same number of columns as the header
            n_hdr = len(header_cols)
            for lineno, dline in data_lines[:3]:  # check first 3 data rows
                n_data = len(dline.split(delimiter))
                if n_data != n_hdr:
                    print(f"  ⚠  WARNING: SS entry {idx} ({label}) CSV data "
                          f"row {lineno} has {n_data} columns but header "
                          f"has {n_hdr}: {filepath}")

        except Exception as _e:
            print(f"  ⚠  WARNING: SS entry {idx} ({label}) could not "
                  f"validate CSV: {_e}")

        # Check species name matches BGC config
        chem_name = ss_config.get("Chemical_name", "")
        if chemical_species_list and chem_name not in chemical_species_list:
            print(f"  ⚠  WARNING: SS entry {idx} species '{chem_name}' "
                  f"not in BGC species list: {chemical_species_list}")

    # Build the complete JSON structure (no METADATA key — C++ parser counts
    # all top-level keys via .size(), so METADATA would be treated as an entry)
    config = {}

    # Metadata goes into // comment lines instead (C++ skips with skip_comments=true)
    metadata_comments = [
        f"// METADATA - Comment: {ss_metadata_comment}",
        f"// METADATA - Source: {ss_metadata_source}",
        f"// METADATA - BGC_Engine: {bgc_engine_label}",
        f"// METADATA - Chemical_Species: {', '.join(chemical_species_list) if chemical_species_list else ''}",
    ]

    # Add each source/sink configuration with numbered keys
    for idx, ss_config in enumerate(ss_method_csv_config, start=1):
        config[str(idx)] = {
            "CHEMICAL_NAME": ss_config["Chemical_name"],
            "COMPARTMENT_NAME": ss_config.get("Compartment_name",
                                               ss_config.get("ss_method_copernicus_compartment_name_for_load", "")),
            "TYPE": ss_config["Type"],
            "UNITS": ss_config["Units"],
            "DATA_FORMAT": "ASCII",
            "DATA": {
                "FILEPATH": ss_config["Filepath"],
                "DELIMITER": ss_config["Delimiter"]
            }
        }

    # Convert to JSON string with standard formatting
    json_string = json.dumps(config, indent=4)

    # Custom formatting function
    def compress_array(match):
        array_content = match.group(0)
        compressed = re.sub(r'\s+', ' ', array_content)
        compressed = re.sub(r'\[\s+', '[', compressed)
        compressed = re.sub(r'\s+\]', ']', compressed)
        compressed = re.sub(r'\s*,\s*', ', ', compressed)
        return compressed

    # Apply compression to all arrays
    json_string = re.sub(r'\[[^\[\]]*\]', compress_array, json_string)

    # Write to file
    with open(ss_config_filepath, 'w') as f:
        # Write comment lines first
        for comment in json_header_comment:
            f.write(comment + "\n")
        # Write metadata as comment lines (C++ parser ignores these)
        for mc in metadata_comments:
            f.write(mc + "\n")
        # Write the JSON content
        f.write(json_string)
        # Add newline at end of file
        f.write("\n")

    print(f"✓ Source/Sink config file saved to: {ss_config_filepath}")


###################
# Optimized Land Use/Land Cover Analysis for Hydrological Response Units (HRUs)
##################
# This script uses rasterstats for efficient zonal statistics calculation.
# Calculates the area of different land use classes from Copernicus LULC rasters
# (NetCDF format) for each polygon in a shapefile representing HRUs.


class OptimizedLULCAnalyzer:
    """Optimized analyzer for land use/land cover in hydrological response units."""

    def __init__(self, shapefile_info, output_dir='output'):
        """
        Initialize the LULC analyzer.

        Parameters:
        -----------
        shapefile_info : dict
            Dictionary containing shapefile information (REQUIRED):
            - 'path_to_shp': Path to the shapefile
            - 'mapping_key': Column name for HRU/basin ID (REQUIRED - must exist in shapefile)
        output_dir : str
            Directory to save output files
        """
        self.shapefile_info = shapefile_info
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)

        # Load shapefile
        # Use fiona-based loading to avoid pyproj CRS resolution issues
        # on older PROJ databases (e.g., PROJ 6.2.1 can't resolve epsg:4326)
        print(f"Loading shapefile: {shapefile_info['path_to_shp']}")
        try:
            self.gdf = gpd.read_file(shapefile_info['path_to_shp'])
        except Exception as crs_err:
            if 'Invalid projection' in str(crs_err) or 'CRSError' in type(crs_err).__name__:
                print(f"  CRS resolution failed ({crs_err}), loading without CRS...")
                import fiona
                from shapely.geometry import shape
                with fiona.open(shapefile_info['path_to_shp']) as src:
                    crs_wkt = src.crs_wkt
                    records = []
                    for feat in src:
                        props = dict(feat['properties'])
                        props['geometry'] = shape(feat['geometry'])
                        records.append(props)
                self.gdf = gpd.GeoDataFrame(records, geometry='geometry')
                # Set CRS from WKT string (bypasses EPSG lookup)
                if crs_wkt:
                    try:
                        self.gdf.crs = crs_wkt
                    except Exception:
                        print(f"  WARNING: Could not set CRS from WKT, proceeding without CRS")
            else:
                raise
        print(f"Loaded {len(self.gdf)} HRU polygons")
        print(f"CRS: {self.gdf.crs}")

        # Require mapping_key to be provided
        if 'mapping_key' not in shapefile_info:
            raise ValueError(
                "Missing 'mapping_key' in shapefile_info. "
                "Please specify which column in the shapefile contains the HRU IDs.\n"
                f"Available columns: {list(self.gdf.columns)}"
            )

        mapping_key = shapefile_info['mapping_key']

        # Check if the mapping_key column exists in the shapefile
        if mapping_key not in self.gdf.columns:
            raise ValueError(
                f"Column '{mapping_key}' specified in 'mapping_key' not found in shapefile.\n"
                f"Available columns: {list(self.gdf.columns)}\n"
                f"Please check the column name and try again."
            )

        # Use the specified mapping key column as HRU_ID
        self.shp_hru_id_column = mapping_key

        print(f"Using '{self.shp_hru_id_column}' column for HRU identification")

    def netcdf_to_geotiff(self, netcdf_path, output_tiff=None, use_temp=False):
        """
        Convert NetCDF to GeoTIFF for easier processing with rasterio.

        Parameters:
        -----------
        netcdf_path : str
            Path to the NetCDF file
        output_tiff : str, optional
            Path for output GeoTIFF. If None, creates temp file.
        use_temp : bool
            If True, save to temporary directory instead of output_dir

        Returns:
        --------
        str : Path to created GeoTIFF
        """
        print(f"  → Opening NetCDF file...", end='', flush=True)
        # Open NetCDF
        ds = xr.open_dataset(netcdf_path)
        print(" ✓")

        # Find land cover variable
        print(f"  → Finding land cover variable...", end='', flush=True)
        lc_var = None
        for var in ['lccs_class', 'land_cover', 'LC', 'lc', 'Band1']:
            if var in ds.variables:
                lc_var = var
                break

        if lc_var is None:
            print(" ✗")
            raise ValueError(f"Could not find LC variable. Available: {list(ds.variables.keys())}")
        print(f" ✓ (found '{lc_var}')")

        # Get data and squeeze extra dimensions
        print(f"  → Loading data...", end='', flush=True)
        lc_data = ds[lc_var]

        # Handle extra dimensions (time, band, etc.)
        # Keep only lat/lon dimensions
        while lc_data.ndim > 2:
            # Squeeze out dimensions of size 1, or take first slice for others
            for dim in lc_data.dims:
                if dim not in ['lat', 'lon', 'latitude', 'longitude']:
                    if lc_data.sizes[dim] == 1:
                        lc_data = lc_data.squeeze(dim)
                    else:
                        # Take first slice if dimension has multiple values
                        lc_data = lc_data.isel({dim: 0})
                    break  # Re-evaluate after each change

        # Now lc_data should be 2D (lat, lon)
        print(f" ✓ (shape: {lc_data.shape})")

        # Get coordinates
        print(f"  → Extracting coordinates...", end='', flush=True)
        if 'lat' in ds.dims and 'lon' in ds.dims:
            _latd, _lond = 'lat', 'lon'
            lats = ds.lat.values
            lons = ds.lon.values
        elif 'latitude' in ds.dims and 'longitude' in ds.dims:
            _latd, _lond = 'latitude', 'longitude'
            lats = ds.latitude.values
            lons = ds.longitude.values
        else:
            print(" ✗")
            raise ValueError(f"Could not find lat/lon. Available dims: {list(ds.dims.keys())}")
        print(" ✓")

        # ── Window to the basin bounding box (+ margin) ──────────────────────
        # ESA-CCI global rasters are ~64800×129600 ≈ 8.4e9 cells; materializing
        # one (the lc_data.values write below) needs ~8 GB RAM and a multi-GB
        # temp write — slow and failure-prone (e.g. the lulc_<year>_temp.tif
        # errors).  Read ONLY the cells around the basin instead.  On ANY problem
        # fall back to the full extent so behaviour is never worse than before.
        try:
            import numpy as _np
            _full_lat, _full_lon = len(lats), len(lons)
            _gdf_ll = self.gdf if getattr(self.gdf, 'crs', None) is None \
                else self.gdf.to_crs(4326)
            _minx, _miny, _maxx, _maxy = [float(v) for v in _gdf_ll.total_bounds]
            _m = 0.05  # ~5 km margin so edge HRUs keep full coverage
            _lat_in = _np.where((lats >= _miny - _m) & (lats <= _maxy + _m))[0]
            _lon_in = _np.where((lons >= _minx - _m) & (lons <= _maxx + _m))[0]
            if _lat_in.size and _lon_in.size:
                _l0, _l1 = int(_lat_in.min()), int(_lat_in.max()) + 1
                _o0, _o1 = int(_lon_in.min()), int(_lon_in.max()) + 1
                lc_data = lc_data.isel({_latd: slice(_l0, _l1),
                                        _lond: slice(_o0, _o1)})
                lats = lats[_l0:_l1]
                lons = lons[_o0:_o1]
                print(f"  → Windowed to basin bbox: {tuple(lc_data.shape)} "
                      f"(from full {_full_lat}×{_full_lon} globe)")
        except Exception as _we:
            print(f"  → (bbox windowing skipped: {_we}; using full extent)")

        # Calculate transform
        print(f"  → Calculating geotransform...", end='', flush=True)
        lon_res = abs(lons[1] - lons[0])
        lat_res = abs(lats[1] - lats[0])

        transform = from_bounds(
            lons.min() - lon_res / 2,
            lats.min() - lat_res / 2,
            lons.max() + lon_res / 2,
            lats.max() + lat_res / 2,
            len(lons),
            len(lats)
        )
        print(" ✓")

        # Create output path
        if output_tiff is None:
            year = self._extract_year(netcdf_path)
            if use_temp:
                import tempfile
                temp_dir = Path(tempfile.gettempdir()) / 'lulc_temp'
                temp_dir.mkdir(exist_ok=True)
                output_tiff = temp_dir / f'lulc_{year}_temp.tif'
            else:
                output_tiff = self.output_dir / f'lulc_{year}.tif'

        # Remove any corrupted leftover from a previous run
        output_tiff = Path(output_tiff)
        if output_tiff.exists():
            output_tiff.unlink()

        # Write to GeoTIFF
        # Use WKT format to avoid PROJ database issues
        wgs84_wkt = 'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]]'

        # Show spinner during long write operation
        import threading
        import time

        stop_spinner = threading.Event()

        def spinner():
            """Display animated spinner during long operations"""
            spinner_chars = ['⠋', '⠙', '⠹', '⠸', '⠼', '⠴', '⠦', '⠧', '⠇', '⠏']
            idx = 0
            while not stop_spinner.is_set():
                print(f'\r  → Writing GeoTIFF... {spinner_chars[idx % len(spinner_chars)]}', end='', flush=True)
                idx += 1
                time.sleep(0.1)

        spinner_thread = threading.Thread(target=spinner, daemon=True)
        spinner_thread.start()

        try:
            with rasterio.open(
                    output_tiff,
                    'w',
                    driver='GTiff',
                    height=len(lats),
                    width=len(lons),
                    count=1,
                    dtype=lc_data.dtype,
                    crs=wgs84_wkt,  # Use WKT instead of 'EPSG:4326'
                    transform=transform,
                    compress='lzw'
            ) as dst:
                # Ensure we write a 2D array
                data_to_write = lc_data.values
                if data_to_write.ndim > 2:
                    # Flatten extra dimensions
                    data_to_write = data_to_write.squeeze()
                dst.write(data_to_write, 1)
        finally:
            stop_spinner.set()
            spinner_thread.join(timeout=1)
            print(f'\r  → Writing GeoTIFF... ✓ ', flush=True)

        ds.close()
        print(f"  ✓ Created GeoTIFF: {Path(output_tiff).name}")
        return str(output_tiff)

    def calculate_lulc_areas(self, raster_path):
        """
        Calculate LULC areas for each HRU using zonal statistics.

        Parameters:
        -----------
        raster_path : str
            Path to the raster file (GeoTIFF)

        Returns:
        --------
        pd.DataFrame : Results with areas per HRU and LC class
        """
        year = self._extract_year(raster_path)
        print(f"  → Calculating LULC areas for {len(self.gdf)} HRUs (year {year})...")

        # Open raster
        with rasterio.open(raster_path) as src:
            raster_crs = src.crs
            transform = src.transform

            # Reproject shapefile to match raster
            print(f"    • Reprojecting HRUs to raster CRS...", end='', flush=True)
            if self.gdf.crs is not None and raster_crs is not None:
                gdf_reprojected = self.gdf.to_crs(raster_crs)
            else:
                # CRS not available (e.g., old PROJ database can't resolve EPSG codes)
                # Assume both are in WGS84/EPSG:4326 and skip reprojection
                print(" (skipping — CRS not set, assuming WGS84)", end='', flush=True)
                gdf_reprojected = self.gdf
            print(" ✓")

            # Get pixel size in meters
            pixel_area_m2 = self._get_pixel_area_m2(src)
            print(f"    • Pixel area: {pixel_area_m2:.2f} m²")

            results = []

            # Process each HRU
            print(f"    • Processing HRUs: ", end='', flush=True)
            for idx, row in gdf_reprojected.iterrows():
                # Get HRU ID using the specified column
                hru_id = row.get(self.shp_hru_id_column, idx)

                # Show progress every 10 HRUs
                if (idx + 1) % 10 == 0:
                    print(f"{idx + 1}", end='', flush=True)
                    if (idx + 1) % 50 == 0:
                        print(f"/{len(gdf_reprojected)}", end='', flush=True)
                    print("...", end='', flush=True)

                # Mask raster with polygon
                try:
                    out_image, out_transform = rasterio_mask(
                        src,
                        [row.geometry],
                        crop=True,
                        nodata=src.nodata
                    )

                    # Get land cover values
                    lc_values = out_image[0]

                    # Remove nodata values
                    if src.nodata is not None:
                        lc_values = lc_values[lc_values != src.nodata]

                    # Remove NaN
                    lc_values = lc_values[~np.isnan(lc_values)]

                    if len(lc_values) == 0:
                        continue

                    # Count pixels per class
                    unique, counts = np.unique(lc_values, return_counts=True)

                    # Calculate areas
                    for lc_class, count in zip(unique, counts):
                        area_m2 = count * pixel_area_m2
                        area_ha = area_m2 / 10_000
                        area_km2 = area_m2 / 1_000_000

                        results.append({
                            self.shp_hru_id_column: hru_id,
                            'Year': year,
                            'LC_Class': int(lc_class),
                            'Pixel_Count': int(count),
                            'Area_m2': area_m2,
                            'Area_ha': area_ha,
                            'Area_km2': area_km2
                        })

                except Exception as e:
                    print(f"\n    ⚠ Warning: Error processing HRU {hru_id}: {e}")
                    continue

            print(f" ✓ ({len(gdf_reprojected)} HRUs processed)")

        return pd.DataFrame(results)

    def _get_pixel_area_m2(self, raster_src):
        """
        Calculate pixel area in square meters.

        Parameters:
        -----------
        raster_src : rasterio.DatasetReader
            Opened raster dataset

        Returns:
        --------
        float : Pixel area in m²
        """
        if raster_src.crs is None or raster_src.crs.is_geographic:
            # For geographic CRS (lat/lon) or unknown CRS, use approximate area
            # Copernicus nominal resolution is 300m
            return 300 * 300
        else:
            # For projected CRS, calculate from transform
            transform = raster_src.transform
            pixel_width = abs(transform.a)
            pixel_height = abs(transform.e)
            return pixel_width * pixel_height

    def process_netcdf_files(self, netcdf_files, convert_to_tiff=True):
        """
        Process multiple NetCDF files.

        Parameters:
        -----------
        netcdf_files : list
            List of paths to NetCDF files
        convert_to_tiff : bool
            Whether to convert NetCDF to GeoTIFF first (recommended)

        Returns:
        --------
        pd.DataFrame : Combined results
        """
        all_results = []

        for nc_file in netcdf_files:
            print(f"\n{'=' * 60}")
            print(f"Processing: {Path(nc_file).name}")
            print(f"{'=' * 60}")

            is_already_tif = str(nc_file).lower().endswith(('.tif', '.tiff'))
            if convert_to_tiff and not is_already_tif:
                # Convert to GeoTIFF
                tiff_path = self.netcdf_to_geotiff(nc_file)
                raster_path = tiff_path
            else:
                raster_path = nc_file

            # Calculate areas
            df_result = self.calculate_lulc_areas(raster_path)
            all_results.append(df_result)

        # Combine all results
        df_combined = pd.concat(all_results, ignore_index=True)

        return df_combined

    def process_netcdf_files_and_clip(self, netcdf_files, clipped_output_dir, clipped_paths_list, convert_to_tiff=True):
        """
        Process multiple NetCDF files and create clipped rasters for entire basin.

        Parameters:
        -----------
        netcdf_files : list
            List of paths to NetCDF files
        clipped_output_dir : Path
            Directory to save clipped rasters
        clipped_paths_list : list
            List to append clipped raster paths to
        convert_to_tiff : bool
            Whether to convert NetCDF to GeoTIFF first (recommended)

        Returns:
        --------
        pd.DataFrame : Combined results
        """
        all_results = []
        temp_files = []  # Track temporary files for cleanup

        for nc_file in netcdf_files:
            print(f"\n{'=' * 60}")
            print(f"Processing: {Path(nc_file).name}")
            print(f"{'=' * 60}")

            is_already_tif = str(nc_file).lower().endswith(('.tif', '.tiff'))
            if convert_to_tiff and not is_already_tif:
                # Convert to GeoTIFF in temporary location
                tiff_path = self.netcdf_to_geotiff(nc_file, use_temp=True)
                raster_path = tiff_path
                temp_files.append(tiff_path)  # Track for cleanup
            else:
                raster_path = nc_file

            # Create clipped raster for entire basin
            year = self._extract_year(raster_path)
            clipped_raster_path = self.clip_raster_to_basin(
                raster_path,
                clipped_output_dir,
                year
            )
            clipped_paths_list.append(clipped_raster_path)

            # Calculate areas
            df_result = self.calculate_lulc_areas(raster_path)
            all_results.append(df_result)

        # Clean up temporary files
        print(f"\nCleaning up temporary files...")
        for temp_file in temp_files:
            try:
                Path(temp_file).unlink()
                print(f"  Deleted: {Path(temp_file).name}")
            except Exception as e:
                print(f"  Warning: Could not delete {temp_file}: {e}")

        # Combine all results
        df_combined = pd.concat(all_results, ignore_index=True)

        return df_combined

    def clip_raster_to_basin(self, raster_path, output_dir, year):
        """
        Clip raster to the extent of the entire shapefile (basin boundary).

        Parameters:
        -----------
        raster_path : str
            Path to the raster file to clip
        output_dir : Path
            Directory to save clipped raster
        year : int
            Year for naming the output file

        Returns:
        --------
        str : Path to the clipped raster
        """
        from rasterio.mask import mask as rasterio_mask
        from shapely.ops import unary_union

        # Create output path
        output_path = output_dir / f'lulc_basin_clipped_{year}.tif'

        print(f"  → Clipping raster to basin boundary...")

        print(f"    • Opening raster...", end='', flush=True)
        with rasterio.open(raster_path) as src:
            print(" ✓")

            # Reproject shapefile to match raster CRS
            print(f"    • Reprojecting shapefile to raster CRS...", end='', flush=True)
            if self.gdf.crs is not None and src.crs is not None:
                gdf_reprojected = self.gdf.to_crs(src.crs)
            else:
                # CRS not available (e.g., old PROJ database can't resolve EPSG codes)
                # Assume both are in WGS84/EPSG:4326 and skip reprojection
                print(" (skipping — CRS not set, assuming WGS84)", end='', flush=True)
                gdf_reprojected = self.gdf
            print(" ✓")

            # Create a single polygon representing the entire basin
            print(f"    • Creating basin boundary polygon...", end='', flush=True)
            basin_boundary = unary_union(gdf_reprojected.geometry)
            print(" ✓")

            # Clip raster to basin boundary
            import threading
            import time

            stop_spinner = threading.Event()

            def spinner():
                """Display animated spinner during long operations"""
                spinner_chars = ['⠋', '⠙', '⠹', '⠸', '⠼', '⠴', '⠦', '⠧', '⠇', '⠏']
                idx = 0
                while not stop_spinner.is_set():
                    print(
                        f'\r    • Clipping raster (this may take a moment)... {spinner_chars[idx % len(spinner_chars)]}',
                        end='', flush=True)
                    idx += 1
                    time.sleep(0.1)

            spinner_thread = threading.Thread(target=spinner, daemon=True)
            spinner_thread.start()

            try:
                # Choose a nodata value compatible with the raster dtype
                if src.nodata is not None:
                    _nodata = src.nodata
                else:
                    _nodata = 0 if np.issubdtype(src.dtypes[0], np.unsignedinteger) else -9999
                out_image, out_transform = rasterio_mask(
                    src,
                    [basin_boundary],
                    crop=True,
                    nodata=_nodata
                )
            finally:
                stop_spinner.set()
                spinner_thread.join(timeout=1)
                print(f'\r    • Clipping raster (this may take a moment)... ✓ ', flush=True)

            # Update metadata
            print(f"    • Preparing metadata...", end='', flush=True)
            out_meta = src.meta.copy()
            out_meta.update({
                "driver": "GTiff",
                "height": out_image.shape[1],
                "width": out_image.shape[2],
                "transform": out_transform,
                "compress": "lzw"
            })
            print(" ✓")

            # Write clipped raster
            print(f"    • Writing clipped raster...", end='', flush=True)
            with rasterio.open(output_path, "w", **out_meta) as dest:
                dest.write(out_image)
            print(" ✓")

        print(f"  ✓ Clipped raster saved: {output_path.name}")
        return str(output_path)

    def _extract_year(self, filename):
        """Extract year from filename."""
        match = re.search(r'(\d{4})', Path(filename).name)
        if match:
            return int(match.group(1))
        return None

    def export_results(self, df, filename='lulc_analysis.csv'):
        """Export results to CSV."""
        output_path = self.output_dir / filename
        df.to_csv(output_path, index=False)
        print(f"\n✓ Results exported to: {output_path}")
        return output_path

    def create_summary_statistics(self, df):
        """Create comprehensive summary statistics."""
        summaries = {}

        # Use the actual HRU ID column name from the shapefile
        hru_col = self.shp_hru_id_column

        # Summary by HRU and Year
        summaries['by_hru_year'] = df.groupby([hru_col, 'Year']).agg({
            'Area_km2': 'sum',
            'LC_Class': 'count'
        }).rename(columns={'LC_Class': 'Num_LC_Classes'}).reset_index()

        # Summary by LC Class and Year
        summaries['by_lc_year'] = df.groupby(['LC_Class', 'Year']).agg({
            'Area_km2': 'sum',
            hru_col: 'count'
        }).rename(columns={hru_col: 'Num_HRUs'}).reset_index()

        # Pivot table: HRUs x LC Classes
        summaries['pivot_hru_lc'] = df.pivot_table(
            values='Area_km2',
            index=[hru_col, 'Year'],
            columns='LC_Class',
            fill_value=0,
            aggfunc='sum'
        )

        return summaries

    def generate_report(self, df):
        """Generate a summary report."""
        print("\n" + "=" * 60)
        print("LULC ANALYSIS REPORT")
        print("=" * 60)

        # Use the actual HRU ID column name
        hru_col = self.shp_hru_id_column

        print(f"\nTotal records: {len(df)}")
        print(f"Number of HRUs: {df[hru_col].nunique()}")
        print(f"Years analyzed: {sorted(df['Year'].unique())}")
        print(f"Land cover classes found: {sorted(df['LC_Class'].unique())}")

        print("\n--- Total Area by Land Cover Class (km²) ---")
        lc_totals = df.groupby('LC_Class')['Area_km2'].sum().sort_values(ascending=False)
        for lc_class, area in lc_totals.items():
            print(f"  Class {lc_class:3d}: {area:10.2f} km²")

        print("\n--- Total Area by Year (km²) ---")
        year_totals = df.groupby('Year')['Area_km2'].sum()
        for year, area in year_totals.items():
            print(f"  {year}: {area:10.2f} km²")


def _clip_netcdf_to_bbox(src_path: str, dst_path: str,
                         bbox: List[float]) -> None:
    """Clip a global NetCDF/GeoTIFF raster to a bounding box using rasterio.

    Uses a windowed read so only the basin extent is loaded into memory,
    even for multi-GB global files.  For ESA CCI LC NetCDF files the land-
    cover data lives in the ``lccs_class`` subdataset; this function detects
    that automatically.

    Parameters
    ----------
    src_path : str
        Path to the source raster (NetCDF or GeoTIFF).
    dst_path : str
        Path where the clipped raster will be saved (GeoTIFF).
    bbox : list of float
        Bounding box as ``[west, south, east, north]`` in the raster's CRS
        (assumed WGS84 / EPSG:4326).
    """
    from rasterio.windows import from_bounds
    from rasterio.crs import CRS

    west, south, east, north = bbox

    # For NetCDF files with subdatasets, open the land-cover subdataset
    open_path = src_path
    if src_path.lower().endswith('.nc'):
        with rasterio.open(src_path) as probe:
            if hasattr(probe, 'subdatasets') and probe.subdatasets:
                # Prefer lccs_class, fall back to first subdataset
                for sd in probe.subdatasets:
                    if 'lccs_class' in sd:
                        open_path = sd
                        break
                else:
                    open_path = probe.subdatasets[0]

    with rasterio.open(open_path) as src:
        # Compute the pixel window that covers the bounding box
        window = from_bounds(west, south, east, north, src.transform)
        # Round to integer pixel indices
        window = window.round_offsets().round_lengths()

        # Read only the bands we need (windowed — memory-efficient)
        data = src.read(window=window)
        transform = src.window_transform(window)

        # Determine CRS — ESA CCI LC uses WGS84 but may not embed it
        crs = src.crs if src.crs is not None else CRS.from_epsg(4326)

        profile = src.profile.copy()
        profile.update(
            driver='GTiff',
            height=data.shape[1],
            width=data.shape[2],
            transform=transform,
            crs=crs,
            compress='lzw',
        )
        # Remove NetCDF-specific keys if present
        for key in ['crs_wkt']:
            profile.pop(key, None)

        with rasterio.open(dst_path, 'w', **profile) as dst:
            dst.write(data)


def _download_copernicus_lc_from_cds(
        basin_shapefile_path: str,
        years: List[int],
        cache_dir: str,
) -> str:
    """Download ESA CCI Land Cover data from the Copernicus Climate Data Store.

    Downloads the global file for each needed year, then clips locally to the
    basin bounding box (with a 0.5° buffer).  The clipped file (a few MB) is
    cached; the large global download (~2.2 GB) is deleted immediately after
    clipping.

    Requires:
        - ``pip install cdsapi``
        - CDS API key configured in ``~/.cdsapirc``
          (see https://cds.climate.copernicus.eu/how-to-api)

    Parameters
    ----------
    basin_shapefile_path : str
        Path to the basin/HRU shapefile (used to compute bounding box).
    years : list of int
        Years to download (e.g. [1993, 1994]).
    cache_dir : str
        Directory to store clipped NetCDF/GeoTIFF files (reused on subsequent runs).

    Returns
    -------
    str
        Path to the cache directory containing the clipped files.
    """
    try:
        import cdsapi
    except ImportError:
        raise ImportError(
            "The 'cdsapi' package is required for auto-downloading Copernicus data.\n"
            "Install with:  pip install cdsapi\n"
            "Then configure your API key: https://cds.climate.copernicus.eu/how-to-api\n"
            "Alternatively, set ss_method_copernicus_nc_lc_dir to a local directory."
        )

    import os, zipfile, shutil

    os.makedirs(cache_dir, exist_ok=True)

    # ── Compute bounding box from basin shapefile ──────────────────────
    basin_gdf = gpd.read_file(basin_shapefile_path)
    if basin_gdf.crs is not None:
        try:
            epsg = basin_gdf.crs.to_epsg()
            if epsg is not None and epsg != 4326:
                basin_gdf = basin_gdf.to_crs("EPSG:4326")
            elif epsg is None:
                bounds = basin_gdf.total_bounds
                if abs(bounds[0]) > 360 or abs(bounds[1]) > 360:
                    basin_gdf = basin_gdf.to_crs("EPSG:4326")
        except Exception:
            pass

    bounds = basin_gdf.total_bounds  # [west, south, east, north]
    # Add 0.5° buffer and clamp to valid range
    bbox = [
        max(float(bounds[0]) - 0.5, -180),   # west
        max(float(bounds[1]) - 0.5, -90),    # south
        min(float(bounds[2]) + 0.5, 180),    # east
        min(float(bounds[3]) + 0.5, 90),     # north
    ]
    print(f"  Basin clip bbox: W={bbox[0]:.2f} S={bbox[1]:.2f} E={bbox[2]:.2f} N={bbox[3]:.2f}")

    # ── Download each year (global) then clip locally ──────────────────
    client = cdsapi.Client()

    for year in sorted(years):
        # Clipped file uses .tif extension (GeoTIFF with LZW compression)
        clipped_filename = f"ESACCI-LC-L4-LCCS-Map-300m-P1Y-{year}-v2.0.7cds.tif"
        clipped_path = os.path.join(cache_dir, clipped_filename)

        # Also check for legacy .nc cached files
        nc_legacy = os.path.join(cache_dir,
                                 f"ESACCI-LC-L4-LCCS-Map-300m-P1Y-{year}-v2.0.7cds.nc")

        if os.path.exists(clipped_path):
            sz = os.path.getsize(clipped_path) / 1e6
            print(f"  Year {year}: cached ({sz:.1f} MB)")
            continue
        if os.path.exists(nc_legacy):
            sz = os.path.getsize(nc_legacy) / 1e6
            print(f"  Year {year}: cached .nc ({sz:.1f} MB)")
            continue

        # CDS versions: v2_0_7cds for 1992-2015, v2_1_1 for 2016+
        version = 'v2_0_7cds' if year <= 2015 else 'v2_1_1'
        print(f"  Downloading year {year} (version {version}) — global file from CDS...")

        zip_path = os.path.join(cache_dir, f"lc_{year}_global.zip")
        try:
            client.retrieve(
                'satellite-land-cover',
                {
                    'variable': 'all',
                    'year': str(year),
                    'version': version,
                    'format': 'zip',
                },
                zip_path,
            )
        except Exception as e:
            print(f"  WARNING: CDS download failed for year {year}: {e}")
            if os.path.exists(zip_path):
                os.remove(zip_path)
            continue

        # Extract .nc file from zip
        global_nc_path = None
        try:
            with zipfile.ZipFile(zip_path, 'r') as zf:
                nc_names = [n for n in zf.namelist() if n.endswith('.nc')]
                if nc_names:
                    zf.extract(nc_names[0], cache_dir)
                    global_nc_path = os.path.join(cache_dir, nc_names[0])
                else:
                    print(f"  WARNING: No .nc file found in zip for year {year}")
                    continue
        except zipfile.BadZipFile:
            print(f"  WARNING: Downloaded zip for year {year} is corrupted — skipping")
            continue
        finally:
            if os.path.exists(zip_path):
                os.remove(zip_path)

        # Clip global file to basin bounding box
        if global_nc_path and os.path.exists(global_nc_path):
            try:
                print(f"  Clipping year {year} to basin extent...", end='', flush=True)
                _clip_netcdf_to_bbox(global_nc_path, clipped_path, bbox)
                sz = os.path.getsize(clipped_path) / 1e6
                print(f" done ({sz:.1f} MB)")
            except Exception as e:
                print(f"\n  WARNING: Clipping failed for year {year}: {e}")
                # Fall back: keep the full global file under the legacy .nc name
                if not os.path.exists(nc_legacy):
                    shutil.move(global_nc_path, nc_legacy)
                    sz = os.path.getsize(nc_legacy) / 1e6
                    print(f"  Kept global file as fallback ({sz:.1f} MB)")
                global_nc_path = None  # prevent deletion below
            finally:
                # Delete global file after clipping (saves ~2 GB per year)
                if global_nc_path and os.path.exists(global_nc_path):
                    os.remove(global_nc_path)

    return cache_dir


def calc_copernicus_lulc(
        ss_config_filepath: str,
        ss_method_copernicus_basins_hrus: Dict[str, str],
        ss_method_copernicus_nc_lc_dir: Optional[str],
        ss_method_copernicus_period: List[Union[int, float]],
        recursive: bool = False,
        file_pattern: str = 'ESACCI-LC-*.nc'
):
    """
    Process Copernicus LULC NetCDF files from a local directory **or** auto-download
    from the Copernicus Climate Data Store (CDS).

    Parameters:
    -----------
    ss_config_filepath : str
        Path where outputs will be saved (directory will be extracted from this)
    ss_method_copernicus_basins_hrus : dict
        Dictionary containing:
        - 'path_to_shp': Path to the shapefile containing HRU/basin polygons
        - 'mapping_key': Column name in shapefile for HRU ID (e.g., 'HRU_ID', 'ID', 'FID')
    ss_method_copernicus_nc_lc_dir : str or None
        Directory containing Copernicus LULC NetCDF files.
        If **None**, data is auto-downloaded from the CDS for the needed years
        (clipped to the basin bounding box — typically a few MB per year
        instead of the full 2.2 GB global files).
        Requires ``cdsapi`` and a CDS API key configured in ``~/.cdsapirc``.
    ss_method_copernicus_period : list
        [year_start, year_end] - Filter files to include only this period
    recursive : bool, default False
        Whether to search subdirectories recursively
    file_pattern : str, default 'ESACCI-LC-*.nc'
        Glob pattern to match NetCDF files

    Returns:
    --------
    tuple : (results_df, summaries, clipped_rasters) - processed results, summary statistics, and paths to clipped rasters
    """

    year_start = ss_method_copernicus_period[0] if ss_method_copernicus_period else None
    year_end = ss_method_copernicus_period[1] if len(ss_method_copernicus_period) > 1 else None

    # Extract output directory from ss_config_filepath
    output_base_dir = Path(ss_config_filepath).parent
    output_base_dir.mkdir(parents=True, exist_ok=True)

    # Create ss_copernicus_files subdirectory for all outputs
    ss_output_dir = output_base_dir / 'ss_copernicus_files'
    ss_output_dir.mkdir(exist_ok=True)

    # Create subdirectory for clipped rasters
    clipped_rasters_dir = ss_output_dir / 'lulc_clipped_rasters'
    clipped_rasters_dir.mkdir(exist_ok=True)

    # ── Reuse precomputed LULC areas if present ─────────────────────────────
    # The LULC areas are climate-INDEPENDENT (they depend only on the basin +
    # ESA-CCI maps).  When lulc_areas_all.csv already exists (produced by the
    # model-config run and seeded into each calibration eval), reuse it instead
    # of re-reading the multi-GB ESA-CCI rasters.  Only the downstream load /
    # climate-adjustment steps re-run — so dynamic coefficients still affect the
    # loading, just without the slow re-clip.  Delete the CSV to force a re-clip.
    _areas_csv = ss_output_dir / 'lulc_areas_all.csv'

    # Offer a clean rebuild when LULC outputs already exist: ask whether to
    # delete everything and start from scratch (otherwise the cache is
    # reused/validated below).  Interactive only — a non-interactive run
    # (calibration / HPC) keeps the existing reuse behaviour.
    _has_lulc = _areas_csv.is_file() or any(clipped_rasters_dir.glob('*.tif'))
    # Skip the prompt during calibration / non-interactive reruns (Gen_Input_Driver
    # sets this when force_regenerate=True) so it is asked at most once per run,
    # not once per evaluation.
    _suppress_prompt = os.environ.get('OPENWQ_SUPPRESS_PROMPTS') == '1'
    if _has_lulc and not _suppress_prompt:
        import sys as _sys
        _wipe = False
        if _sys.stdin and _sys.stdin.isatty():
            try:
                _ans = input(
                    f"\n  LULC outputs already exist in {ss_output_dir}.\n"
                    f"  Delete ALL of them and rebuild from scratch? [y/N]: "
                ).strip().lower()
                _wipe = _ans in ('y', 'yes')
            except (EOFError, KeyboardInterrupt):
                _wipe = False
        if _wipe:
            import shutil as _shutil
            for _p in ss_output_dir.iterdir():
                try:
                    _shutil.rmtree(_p) if _p.is_dir() else _p.unlink()
                except Exception:
                    pass
            clipped_rasters_dir.mkdir(exist_ok=True)
            print("  ↳ Cleared LULC outputs — rebuilding from scratch.")

    if _areas_csv.is_file():
        # Before reusing the cache, VERIFY it was built for the CURRENT basin:
        # compare the cached GRU_IDs against the basin shapefile's mapping_key
        # values.  Acknowledge the outcome in the console either way — reuse when
        # they match, reprocess when they don't (e.g. the basin shapefile was
        # swapped, lumped GRU → delineated per-reach; reusing would silently
        # mis-key every load and produce all-zero output).
        _reuse_cache = True
        _msg = (f"  Reusing cached LULC areas ({_areas_csv.name}) — "
                "skipping the slow ESA-CCI re-clip.")
        try:
            _bi = ss_method_copernicus_basins_hrus or {}
            _shp = _bi.get('path_to_shp')
            _mk = _bi.get('mapping_key', 'GRU_ID')

            def _norm_id(v):
                try:
                    _fv = float(v)
                    return str(int(_fv)) if _fv.is_integer() else str(v)
                except (TypeError, ValueError):
                    return str(v)

            _cur_ids = set()
            _shp_name = os.path.basename(_shp) if _shp else '?'
            if _shp and os.path.isfile(_shp):
                import geopandas as _gpd
                _bg = _gpd.read_file(_shp)
                if _mk in _bg.columns:
                    _cur_ids = {_norm_id(v) for v in _bg[_mk].tolist()}
            # The cache's id column is named after the mapping_key used when it
            # was BUILT.  If the CURRENT mapping_key column isn't present, the
            # cache was built for a different basin/key (e.g. lumped 'SubId' →
            # delineated 'GRU_ID') — reprocess instead of reusing, which would
            # otherwise KeyError downstream on the missing column.
            _cache_cols = list(pd.read_csv(_areas_csv, nrows=0).columns)
            if _mk not in _cache_cols:
                _reuse_cache = False
                _msg = (f"  ↻ Cached LULC areas lack the current mapping_key "
                        f"column '{_mk}' (cache has {_cache_cols}) — built for a "
                        f"different basin/key; reprocessing for {_shp_name}.")
            else:
                _cached = pd.read_csv(_areas_csv, usecols=[_mk])
                _cache_ids = {_norm_id(v) for v in _cached[_mk].unique()}
                if _cur_ids and _cache_ids:
                    if _cur_ids == _cache_ids:
                        _msg = (f"  ✓ Cached LULC areas match the current basin "
                                f"({len(_cache_ids)} unit(s) in {_shp_name}) — "
                                "reusing them, skipping the ESA-CCI re-clip.")
                    else:
                        _reuse_cache = False
                        _msg = (f"  ↻ Cached LULC areas were built for a DIFFERENT "
                                f"basin ({len(_cache_ids)} cached unit(s) vs "
                                f"{len(_cur_ids)} in {_shp_name}) — reprocessing the "
                                "LULC for the current basin.")
                else:
                    _msg = (f"  Reusing cached LULC areas ({_areas_csv.name}) — could "
                            "not verify against the basin shapefile; assuming a match.")
        except Exception:
            _reuse_cache = True  # verification failed → keep the safe reuse default

        # Year-coverage check (independent of the basin check above): a cache
        # for the CORRECT basin is still stale if it covers fewer years than
        # the clipped rasters available for this basin.  Reusing such a cache
        # silently collapses every simulation year onto the single cached year
        # (the downstream proxy-year mapping then maps e.g. 1981..2009 → 1993).
        # Drop & recompute whenever the clipped rasters span years the cache
        # lacks.
        if _reuse_cache:
            try:
                _raster_years = {
                    _extract_year_from_filename(p.name)
                    for p in clipped_rasters_dir.glob('*.tif')
                }
                _raster_years.discard(None)
                _cache_years = {
                    int(y) for y in
                    pd.read_csv(_areas_csv, usecols=['Year'])['Year']
                      .dropna().unique()
                }
                _missing_years = _raster_years - _cache_years
                if _raster_years and _missing_years:
                    _reuse_cache = False
                    _msg = (
                        f"  ↻ Cached LULC areas cover only {sorted(_cache_years)} "
                        f"but clipped rasters span {min(_raster_years)}–"
                        f"{max(_raster_years)} ({len(_raster_years)} years) — "
                        "recomputing so the dynamic loads use every available "
                        "land-cover year.")
            except Exception:
                pass  # if the check can't run, keep the basin-based decision

        print("\n" + _msg)
        if _reuse_cache:
            results_df = pd.read_csv(_areas_csv)
            summaries = {}
            try:
                _an = OptimizedLULCAnalyzer(
                    ss_method_copernicus_basins_hrus, output_dir=str(ss_output_dir))
                summaries = _an.create_summary_statistics(results_df)
            except Exception:
                summaries = {}
            existing_rasters = sorted(str(p) for p in clipped_rasters_dir.glob('*.tif'))
            return results_df, summaries, existing_rasters
        else:
            try:
                _areas_csv.unlink()   # drop the stale cache, then re-clip below
            except Exception:
                pass

    # ── Auto-download from CDS if no local directory provided ──────────
    if ss_method_copernicus_nc_lc_dir is None:
        cache_dir = str(ss_output_dir / 'copernicus_cache')
        years_needed = list(range(
            int(year_start) if year_start else 1992,
            (int(year_end) if year_end else 2022) + 1,
        ))
        print(f"\nNo local Copernicus LULC directory provided — downloading from CDS.")
        print(f"  Years needed: {years_needed}")
        print(f"  Cache directory: {cache_dir}")
        ss_method_copernicus_nc_lc_dir = _download_copernicus_lc_from_cds(
            basin_shapefile_path=ss_method_copernicus_basins_hrus['path_to_shp'],
            years=years_needed,
            cache_dir=cache_dir,
        )

    print(f"Base output directory: {output_base_dir}")
    print(f"SS Copernicus files directory: {ss_output_dir}")
    print(f"Clipped rasters will be saved to: {clipped_rasters_dir}")

    # Convert directory path to Path object
    lc_dir = Path(ss_method_copernicus_nc_lc_dir)

    # Check if directory exists.  This point is reached ONLY when the cached
    # LULC areas (ss_copernicus_files/lulc_areas_all.csv) were NOT reused — a
    # fresh basin or a missing/mismatched cache — so the (multi-GB) ESA-CCI
    # rasters are genuinely needed.  On HPC that source dir is usually a local
    # path that was never copied, so fail with actionable guidance.
    if not lc_dir.exists():
        raise FileNotFoundError(
            "Copernicus LULC source directory not found:\n"
            f"    {ss_method_copernicus_nc_lc_dir}\n"
            "It is needed only because the cached LULC areas "
            "(ss_copernicus_files/lulc_areas_all.csv) were NOT reused for this "
            "basin.  On HPC the ESA-CCI rasters are typically NOT copied (tens of "
            "GB).  Fix by EITHER:\n"
            "  (a) provide a VALID lulc_areas_all.csv built locally for THIS basin "
            "shapefile (its GRU_IDs must match the basin's mapping_key values) so "
            "the re-clip is skipped; or\n"
            "  (b) copy the ESA-CCI raster(s) to a path that exists here and set "
            "ss_method_copernicus_nc_lc_dir to it.")

    if not lc_dir.is_dir():
        raise NotADirectoryError(f"Path is not a directory: {ss_method_copernicus_nc_lc_dir}")

    # Find all raster files (NetCDF or GeoTIFF from auto-download clipping)
    print(f"Searching for raster files in: {lc_dir}")
    print(f"Pattern: {file_pattern}")
    print(f"Recursive: {recursive}")

    if recursive:
        # Search recursively
        netcdf_files = sorted(lc_dir.rglob(file_pattern))
    else:
        # Search only in the specified directory
        netcdf_files = sorted(lc_dir.glob(file_pattern))

    # Also pick up .tif files (produced by auto-download local clipping)
    if not netcdf_files or file_pattern.endswith('.nc'):
        tif_pattern = file_pattern.replace('.nc', '.tif')
        if recursive:
            tif_files = sorted(lc_dir.rglob(tif_pattern))
        else:
            tif_files = sorted(lc_dir.glob(tif_pattern))
        # Merge, avoiding duplicates (same stem)
        existing_stems = {f.stem for f in netcdf_files}
        for tf in tif_files:
            if tf.stem not in existing_stems:
                netcdf_files.append(tf)
        netcdf_files = sorted(netcdf_files, key=lambda p: p.name)

    if not netcdf_files:
        # Try alternative patterns
        alternative_patterns = ['*.nc', '*.tif', 'ESACCI-*.nc', 'ESACCI-*.tif', 'LC-*.nc']
        print(f"\nNo files found with pattern '{file_pattern}'")
        print("Trying alternative patterns...")

        for alt_pattern in alternative_patterns:
            if recursive:
                netcdf_files = sorted(lc_dir.rglob(alt_pattern))
            else:
                netcdf_files = sorted(lc_dir.glob(alt_pattern))

            if netcdf_files:
                print(f"Found {len(netcdf_files)} files with pattern '{alt_pattern}'")
                break

    if not netcdf_files:
        raise FileNotFoundError(
            f"No NetCDF files found in {ss_method_copernicus_nc_lc_dir}\n"
            f"Searched pattern: {file_pattern}\n"
            f"Recursive: {recursive}\n"
            f"Please check the directory path and ensure it contains .nc files"
        )

    print(f"\nFound {len(netcdf_files)} NetCDF files")

    # Filter by year if specified
    if year_start is not None or year_end is not None:
        all_netcdf_files = list(netcdf_files)   # keep unfiltered copy for fallback
        filtered_files = []

        for nc_file in netcdf_files:
            # Extract year from filename
            year = _extract_year_from_filename(nc_file.name)

            if year is not None:
                include = True

                if year_start is not None and year < year_start:
                    include = False
                if year_end is not None and year > year_end:
                    include = False

                if include:
                    filtered_files.append(nc_file)
            else:
                # Include files where year cannot be determined
                print(f"Warning: Could not extract year from {nc_file.name}, including anyway")
                filtered_files.append(nc_file)

        netcdf_files = filtered_files
        print(f"After year filtering ({year_start}-{year_end}): {len(netcdf_files)} files")

        # ── Nearest-year fallback (proxy) ──────────────────────────────
        if not netcdf_files and all_netcdf_files:
            available_years = sorted(set(
                _extract_year_from_filename(f.name)
                for f in all_netcdf_files
                if _extract_year_from_filename(f.name) is not None
            ))
            if available_years:
                target = year_start if year_start is not None else year_end
                nearest = min(available_years, key=lambda y: abs(y - target))
                print(f"     ⚠  No files found for requested year range "
                      f"{year_start}–{year_end}.")
                print(f"        Available years in directory: {available_years}")
                print(f"        → Using nearest available year ({nearest}) as proxy.")
                netcdf_files = [
                    f for f in all_netcdf_files
                    if _extract_year_from_filename(f.name) == nearest
                ]
                print(f"        Files selected: {len(netcdf_files)}")

    if not netcdf_files:
        raise ValueError(
            f"No files remaining after filtering for years {year_start}-{year_end}"
        )

    # Convert to string paths
    netcdf_files = [str(f) for f in netcdf_files]

    # Display found files
    print(f"\nNetCDF files to process:")
    for i, nc_file in enumerate(netcdf_files, 1):
        year = _extract_year_from_filename(Path(nc_file).name)
        year_str = f" (Year: {year})" if year else ""
        print(f"  {i:2d}. {Path(nc_file).name}{year_str}")

    # Initialize analyzer
    print(f"\n{'=' * 60}")
    print(f"Initializing analyzer")
    print(f"{'=' * 60}")
    print(f"Shapefile: {ss_method_copernicus_basins_hrus['path_to_shp']}")
    print(f"HRU ID column (mapping_key): {ss_method_copernicus_basins_hrus['mapping_key']}")
    print(f"Output directory: {ss_output_dir}")

    analyzer = OptimizedLULCAnalyzer(
        ss_method_copernicus_basins_hrus,
        output_dir=str(ss_output_dir)
    )

    # Process files and create clipped rasters
    print(f"\n{'=' * 60}")
    print(f"Processing {len(netcdf_files)} NetCDF files")
    print(f"{'=' * 60}")

    clipped_raster_paths = []
    results_df = analyzer.process_netcdf_files_and_clip(
        netcdf_files,
        clipped_rasters_dir,
        clipped_raster_paths
    )

    # Export main results
    print("\n" + "=" * 60)
    print("EXPORTING RESULTS")
    print("=" * 60)
    analyzer.export_results(results_df, 'lulc_areas_all.csv')

    # Create and export summaries
    print("\nCreating summary statistics...")
    summaries = analyzer.create_summary_statistics(results_df)

    analyzer.export_results(summaries['by_hru_year'], 'summary_by_hru_year.csv')
    analyzer.export_results(summaries['by_lc_year'], 'summary_by_lc_year.csv')
    analyzer.export_results(summaries['pivot_hru_lc'], 'pivot_hru_lc.csv')

    # Generate report
    analyzer.generate_report(results_df)

    print("\n" + "=" * 60)
    print("✓ PROCESSING COMPLETE!")
    print("=" * 60)
    print(f"All results saved to: {ss_output_dir}/")
    print(f"  - lulc_areas_all.csv (detailed results)")
    print(f"  - summary_by_hru_year.csv")
    print(f"  - summary_by_lc_year.csv")
    print(f"  - pivot_hru_lc.csv")
    print(f"\nClipped rasters saved to: {clipped_rasters_dir}/")
    for raster_path in clipped_raster_paths:
        print(f"  - {Path(raster_path).name}")

    # Create nutrient loads reference file
    _create_nutrient_loads_reference(ss_output_dir)

    return results_df, summaries, clipped_raster_paths


def _create_nutrient_loads_reference(output_dir: Path):
    """
    Create a reference file with nutrient load coefficients by land use class.

    Parameters:
    -----------
    output_dir : Path
        Directory where the reference file will be saved
    """
    reference_content = """### Nutrient loads per hectare by land use class

Values below are typical **annual export coefficients** (runoff loads) normalized per hectare. They vary with climate, soils, and management, but give a usable starting point.

### Approximate export ranges (kg/ha/year)

| Land use class               | Total N (TN) | NO3–N | NO2–N | NH4–N | Total P (TP) | PO4–P | Citations |
|-----------------------------|-------------|-------|-------|-------|--------------|-------|-----------|
| Forest / natural catchment  | ~1–5        | ~0.5–3  | ≲0.01 | ~0.05–0.3 | ~0.05–0.3   | ~0.005–0.02 | Yoon et al. 2010; Clesceri et al. 1986; Wit et al. 2020 |
| Mixed rural (forest + ag)   | ~15–25      | ~5–15  | ≲0.05 | ~0.2–1.0   | ~0.4–1.0    | ~0.05–0.2  | Salvia‑Castellví et al. 2005; Yoon et al. 2010; Clesceri et al. 1986 |
| Grassland / pasture         | ~1–11       | ~0.5–5 | ≲0.02 | ~0.05–0.4 | ~0.1–0.3    | ~0.01–0.05 | Hopkins et al. 2025; Harmel et al. 2008; Clesceri et al. 1986 |
| Row‑crop / cultivated ag    | ~10–30 (often ~14) | ~5–25 | ≲0.1  | ~0.2–2.0 | ~1–5 (often ~2) | ~0.1–0.8 | Harmel et al. 2008; Hopkins et al. 2025; Salvia‑Castellví et al. 2005; Wang et al. 2015 |
| Managed boreal forest (ditched, some arable) | ~2.3 | ~0.8–1.5 | ≲0.01 | ~0.1–0.3 | ~0.10 | ~0.01–0.03 | Aaltonen et al. 2021; Yli‑Halla et al. 2021 |

**Figure 1:** Indicative TN, TP and inorganic species exports by land use

### Key basis for NO3, NO2, NH4 and PO4 values

- A mixed forest watershed in Korea exported **36.5 kg TN/ha/yr**, partitioned into **26.9 kg NO3–N**, **1.6 kg NH4–N**, and **0.1 kg PO4–P** per ha per year, with NO2–N negligible; this underpins the forest/natural and mixed‑rural orders of magnitude   (Yoon et al., 2010).  
- Wisconsin export‑coefficient work shows that **TIN (largely NO3–N)** and **ortho‑P (PO4–P)** exports rise systematically from forest to mixed to agricultural basins, with OP and TIN typically making up a substantial share of TP and TN, respectively   (Clesceri et al., 1986).  
- Field‑scale datasets compiled in the MANAGE database indicate that **>90% of inorganic N exports from tile‑drained row‑crop systems are NO3–N**, with NH4–N and NO2–N small fractions, supporting the higher NO3–N ranges assigned to cultivated agriculture relative to grassland   (Mercado & Engel, 2021; Hopkins et al., 2025).  
- In managed boreal organic soils, about **41% of TN was NO3–N and 23% NH4–N**, giving NH4–N loads on the order of **0.3–4 kg/ha/yr** where TN is 10–20 kg/ha/yr   (Yli‐Halla et al., 2021).  
- Export‑coefficient modeling for dissolved nutrients in the Three Gorges Reservoir Region suggests that **dissolved inorganic N (primarily NO3–N) and PO4–P from dryland and paddy agriculture dominate dissolved loads**, consistent with the upper ranges used here for row‑crop NO3–N and PO4–P   (Wang et al., 2015).  

These species‑level coefficients are intended as **order‑of‑magnitude planning values**. For detailed studies, replace them with region‑specific ratios of NO3–N, NO2–N, NH4–N and PO4–P to TN and TP derived from local monitoring or calibrated models.

## References
 
Hopkins, A., Harmel, R., Kleinman, P., Sahoo, D., & Ippolito, J. (2025). Nutrient Runoff From Agricultural Lands in North American Ecoregions. *JAWRA Journal of the American Water Resources Association*, 61. https://doi.org/10.1111/1752-1688.70004
 
Mercado, J., & Engel, B. (2021). Multi-Scale Analysis of the Dependence of Water Quality on Land Use Using Linear and Mixed Models. *Water*. https://doi.org/10.3390/w13192618
 
Yoon, S., Chung, S., Oh, D., & Lee, J. (2010). Monitoring of non-point source pollutants load from a mixed forest land use.. *Journal of environmental sciences*, 22 6, 801-5. https://doi.org/10.1016/s1001-0742(09)60180-7
 
Yli‐Halla, M., Lötjönen, T., Kekkonen, J., Virtanen, S., Marttila, H., Liimatainen, M., Saari, M., Mikkola, J., Suomela, R., & Joki-Tokola, E. (2021). Thickness of peat influences the leaching of substances and greenhouse gas emissions from a cultivated organic soil.. *The Science of the total environment*, 806 Pt 1, 150499. https://doi.org/10.1016/j.scitotenv.2021.150499
 
Clesceri, N., Curran, S., & Sedlák, R. (1986). NUTRIENT LOADS TO WISCONSIN LAKES: PART I. NITROGEN AND PHOSPHORUS EXPORT COEFFICIENTS1. *Journal of The American Water Resources Association*, 22, 983-990. https://doi.org/10.1111/j.1752-1688.1986.tb00769.x
 
Wang, J., Shao, J., Wang, D., Ni, J., & Xie, D. (2015). Simulation of the dissolved nitrogen and phosphorus loads in different land uses in the Three Gorges Reservoir Region--based on the improved export coefficient model.. *Environmental science. Processes & impacts*, 17 11, 1976-89. https://doi.org/10.1039/c5em00380f
"""

    reference_file = output_dir / 'SS_lulc_loads_reference.md'

    with open(reference_file, 'w', encoding='utf-8') as f:
        f.write(reference_content)

    print(f"\n✓ Nutrient loads reference file created: SS_lulc_loads_reference.md")


def _extract_year_from_filename(filename: str) -> Optional[int]:
    """
    Extract year from filename.

    Parameters:
    -----------
    filename : str
        Filename to extract year from

    Returns:
    --------
    int or None : Extracted year or None if not found
    """
    # Look for 4-digit year
    match = re.search(r'(\d{4})', filename)
    if match:
        year = int(match.group(1))
        # Validate it's a reasonable year (1900-2100)
        if 1900 <= year <= 2100:
            return year

    return None


def list_available_files(
        netcdf_copernicus_lc_dir: str,
        file_pattern: str = 'ESACCI-LC-*.nc',
        recursive: bool = False
) -> List[str]:
    """
    List all available NetCDF files in directory without processing.

    Parameters:
    -----------
    netcdf_copernicus_lc_dir : str
        Directory containing NetCDF files
    file_pattern : str
        Glob pattern to match files
    recursive : bool
        Whether to search recursively

    Returns:
    --------
    list : List of file paths
    """
    lc_dir = Path(netcdf_copernicus_lc_dir)

    if not lc_dir.exists():
        raise FileNotFoundError(f"Directory not found: {netcdf_copernicus_lc_dir}")

    if recursive:
        files = sorted(lc_dir.rglob(file_pattern))
    else:
        files = sorted(lc_dir.glob(file_pattern))

    print(f"Found {len(files)} files:")
    for i, f in enumerate(files, 1):
        year = _extract_year_from_filename(f.name)
        year_str = f" ({year})" if year else ""
        print(f"  {i:2d}. {f.name}{year_str}")

    return [str(f) for f in files]


# Main wrapper function
def set_ss_from_copernicus_lulc(
        ss_config_filepath: str,
        ss_method_copernicus_basins_hrus: Dict[str, str],
        ss_method_copernicus_nc_lc_dir: Optional[str],
        ss_method_copernicus_period: List[Union[int, float]],
        recursive: bool = False,
        file_pattern: str = 'ESACCI-LC-*.nc'
):
    """
    Main function to set source/sink from Copernicus LULC data.

    Parameters:
    -----------
    ss_config_filepath : str
        Path where configuration and outputs will be saved
    ss_method_copernicus_basins_hrus : dict
        Dictionary containing:
        - 'path_to_shp': Path to the shapefile containing HRU/basin polygons
        - 'mapping_key': Column name in shapefile for HRU ID (e.g., 'HRU_ID', 'ID', 'FID')
    ss_method_copernicus_nc_lc_dir : str
        Directory containing Copernicus LULC NetCDF files
    ss_method_copernicus_period : list
        [year_start, year_end] - Filter files to include only this period
    recursive : bool, default False
        Whether to search subdirectories recursively
    file_pattern : str, default 'ESACCI-LC-*.nc'
        Glob pattern to match NetCDF files

    Returns:
    --------
    tuple : (results_df, summaries, clipped_rasters)
    """
    # Calculate lulc areas using copernicus global data
    results, summaries, rasters = calc_copernicus_lulc(
        ss_config_filepath=ss_config_filepath,
        ss_method_copernicus_basins_hrus=ss_method_copernicus_basins_hrus,
        ss_method_copernicus_nc_lc_dir=ss_method_copernicus_nc_lc_dir,
        ss_method_copernicus_period=ss_method_copernicus_period,
        recursive=recursive,
        file_pattern=file_pattern
    )

    return results, summaries, rasters


def calculate_nutrient_loads(
        results_df: pd.DataFrame,
        load_coefficients: Dict[int, Dict[str, float]],
        output_dir: Path,
        shp_hru_id_column: str = 'HRU_ID'
) -> pd.DataFrame:
    """
    Calculate nutrient loads based on land cover areas and export coefficients.

    Parameters:
    -----------
    results_df : pd.DataFrame
        DataFrame with LULC areas per HRU (output from calc_copernicus_lulc)
    load_coefficients : dict
        Dictionary mapping LC_Class to nutrient coefficients (kg/ha/yr).
        Example:
        {
            10: {'TN': 14.0, 'TP': 2.0},  # Cropland
            50: {'TN': 3.0, 'TP': 0.15},  # Forest
            130: {'TN': 6.0, 'TP': 0.2}   # Grassland
        }
        Classes not in this dict will have loads set to 0.
    output_dir : Path
        Directory to save the loads CSV file
    shp_hru_id_column : str
        Name of the HRU ID column

    Returns:
    --------
    pd.DataFrame : Nutrient loads per HRU, year, and nutrient type
    """
    print("\n" + "=" * 60)
    print("CALCULATING NUTRIENT LOADS")
    print("=" * 60)

    # Get unique LC classes from results
    unique_lc_classes = sorted(results_df['LC_Class'].unique())
    print(f"\nLand cover classes found in data: {unique_lc_classes}")

    # Check which classes have coefficients
    classes_with_coeffs = [lc for lc in unique_lc_classes if lc in load_coefficients]
    classes_without_coeffs = [lc for lc in unique_lc_classes if lc not in load_coefficients]

    print(f"Classes with load coefficients: {classes_with_coeffs}")
    if classes_without_coeffs:
        print(f"⚠ Warning: Classes without coefficients (will be set to 0): {classes_without_coeffs}")

    # Get all nutrient types from coefficients
    nutrient_types = set()
    for coeffs in load_coefficients.values():
        nutrient_types.update(coeffs.keys())
    nutrient_types = sorted(nutrient_types)
    print(f"Nutrient types to calculate: {nutrient_types}")

    # Calculate loads
    loads_data = []

    print(f"\nCalculating loads for {len(results_df)} records...")
    for idx, row in results_df.iterrows():
        lc_class = int(row['LC_Class'])
        area_ha = row['Area_ha']
        hru_id = row[shp_hru_id_column]
        year = row['Year']

        # Get coefficients for this LC class
        if lc_class in load_coefficients:
            coeffs = load_coefficients[lc_class]

            # Calculate load for each nutrient type
            for nutrient in nutrient_types:
                coeff = coeffs.get(nutrient, 0.0)  # Default to 0 if nutrient not specified
                load_kg_yr = area_ha * coeff  # kg/yr = ha * (kg/ha/yr)

                loads_data.append({
                    shp_hru_id_column: hru_id,
                    'Year': year,
                    'LC_Class': lc_class,
                    'Area_ha': area_ha,
                    'Nutrient': nutrient,
                    'Coefficient_kg_ha_yr': coeff,
                    'Load_kg_yr': load_kg_yr
                })
        else:
            # No coefficients for this class - set load to 0
            for nutrient in nutrient_types:
                loads_data.append({
                    shp_hru_id_column: hru_id,
                    'Year': year,
                    'LC_Class': lc_class,
                    'Area_ha': area_ha,
                    'Nutrient': nutrient,
                    'Coefficient_kg_ha_yr': 0.0,
                    'Load_kg_yr': 0.0
                })

        if (idx + 1) % 1000 == 0:
            print(f"  Processed {idx + 1}/{len(results_df)} records...")

    loads_df = pd.DataFrame(loads_data)

    # Export detailed loads
    loads_file = output_dir / 'nutrient_loads_detailed.csv'
    loads_df.to_csv(loads_file, index=False)
    print(f"\n✓ Detailed loads exported to: {loads_file.name}")

    # Create summary by HRU and Year
    summary_by_hru = loads_df.groupby([shp_hru_id_column, 'Year', 'Nutrient'])['Load_kg_yr'].sum().reset_index()
    summary_file = output_dir / 'nutrient_loads_summary_by_hru.csv'
    summary_by_hru.to_csv(summary_file, index=False)
    print(f"✓ HRU summary exported to: {summary_file.name}")

    # Create pivot table: HRUs x Nutrients
    pivot_loads = loads_df.pivot_table(
        values='Load_kg_yr',
        index=[shp_hru_id_column, 'Year'],
        columns='Nutrient',
        aggfunc='sum',
        fill_value=0
    )
    pivot_file = output_dir / 'nutrient_loads_pivot.csv'
    pivot_loads.to_csv(pivot_file)
    print(f"✓ Pivot table exported to: {pivot_file.name}")

    # Print summary statistics
    print("\n--- Nutrient Load Summary (kg/yr) ---")
    for nutrient in nutrient_types:
        nutrient_loads = loads_df[loads_df['Nutrient'] == nutrient]
        total_load = nutrient_loads['Load_kg_yr'].sum()
        print(f"  {nutrient}: {total_load:,.2f} kg/yr")

    return loads_df, summary_by_hru


def _calculate_temporal_load_distribution(
        annual_load_kg: float,
        year: int,
        method: str = "uniform",
        climate_data: Optional[Dict] = None,
        precip_scaling_power: float = 1.0,
        temp_q10: float = 2.0,
        temp_reference_c: float = 15.0
) -> List[tuple]:
    """
    Calculate temporal distribution of annual nutrient load.

    Parameters:
    -----------
    annual_load_kg : float
        Total annual load in kg
    year : int
        Year for the load
    method : str
        Distribution method:
        - "uniform": Spread load uniformly throughout the year (equal 1/12 each month)
        - "seasonal": Climate-weighted distribution using ss_climate_data
          (precipitation + temperature). Falls back to a sine-curve approximation
          if no climate data is available for the requested year.
    climate_data : dict, optional
        Dict keyed by year: {year: {'precip_mm': [12 floats], 'temp_c': [12 floats]}}.
        Used when method="seasonal" to weight monthly loads by precipitation
        and temperature.
    precip_scaling_power : float
        Exponent for precipitation scaling (default 1.0 = linear)
    temp_q10 : float
        Q10 temperature sensitivity factor (default 2.0)
    temp_reference_c : float
        Reference temperature [°C] at which temperature factor = 1.0

    Returns:
    --------
    list of tuples : [(month, day, hour, load_kg), ...]
        List of temporal load entries.  The sum of all monthly loads
        always equals *annual_load_kg*.
    """
    import math

    # "harmonic"/"monthly" use a FLAT base (equal 1/12 each month); the seasonal
    # shape is imposed downstream by the calibration's per-month multipliers.
    if method in ("uniform", "harmonic", "monthly"):
        # Distribute load uniformly across 12 months (1st of each month at 1:00 AM)
        monthly_load = annual_load_kg / 12.0

        temporal_entries = []
        for month in range(1, 13):
            # Using day=1, hour=1, allows us to use "all" for other time fields
            temporal_entries.append((month, 1, 1, monthly_load))

        return temporal_entries

    elif method == "seasonal":
        # Preferred: use climate data (precip + temperature) when available.
        # Resolve this year's monthly cycle from an explicit per-year entry, a
        # string-keyed year, OR a single convenience entry under "all" (same
        # monthly cycle applied to every simulation year).
        year_data = None
        if isinstance(climate_data, dict) and climate_data:
            year_data = (climate_data.get(year)
                         or climate_data.get(str(year))
                         or climate_data.get("all")
                         or climate_data.get("ALL"))
        if year_data:
            if ('precip_mm' in year_data and 'temp_c' in year_data
                    and len(year_data['precip_mm']) == 12
                    and len(year_data['temp_c']) == 12):
                return _apply_climate_adjustment(
                    annual_load_kg=annual_load_kg,
                    year=year,
                    monthly_precip_mm=year_data['precip_mm'],
                    monthly_temp_c=year_data['temp_c'],
                    precip_scaling_power=precip_scaling_power,
                    temp_q10=temp_q10,
                    temp_reference_c=temp_reference_c
                )

        # Fallback: sine-curve approximation (no climate data)
        monthly_loads = []
        for month in range(1, 13):
            # Sine function centred on July (month 7)
            angle = 2 * math.pi * (month - 4) / 12
            # Amplitude: min is 0.2 of mean, max is 1.8 of mean
            relative_load = 1.0 + 0.8 * math.sin(angle)
            monthly_loads.append(relative_load)

        # Normalize so sum equals annual load
        total_relative = sum(monthly_loads)
        monthly_loads_kg = [(annual_load_kg * rel / total_relative) for rel in monthly_loads]

        # Create temporal entries (1st of each month at 1:00 AM)
        temporal_entries = []
        for month, load_kg in enumerate(monthly_loads_kg, start=1):
            temporal_entries.append((month, 1, 1, load_kg))

        return temporal_entries

    else:
        raise ValueError(
            f"Unknown temporal distribution method: '{method}'. "
            f"Valid options are 'uniform', 'seasonal', 'harmonic', or 'monthly'."
        )


def _apply_climate_adjustment(
        annual_load_kg: float,
        year: int,
        monthly_precip_mm: List[float],
        monthly_temp_c: List[float],
        precip_scaling_power: float = 1.0,
        temp_q10: float = 2.0,
        temp_reference_c: float = 15.0
) -> List[tuple]:
    """
    Distribute annual nutrient load across months using climate adjustments.

    The approach combines precipitation-scaling and temperature-scaling:
      - Precipitation scaling: months with more precipitation export proportionally more
        nutrients (runoff-driven), raised to a power to control nonlinearity.
      - Temperature scaling: biological processes (mineralization, denitrification)
        increase with temperature following a Q10 relationship.

    The combined monthly weight is: w_m = P_m^alpha * Q10^((T_m - T_ref) / 10)
    Monthly loads are then: L_m = annual_load * w_m / sum(w_m)

    Parameters:
    -----------
    annual_load_kg : float
        Total annual load in kg
    year : int
        Year for the load
    monthly_precip_mm : list of 12 floats
        Monthly precipitation totals [mm] for the year (Jan-Dec)
    monthly_temp_c : list of 12 floats
        Monthly mean temperature [deg C] for the year (Jan-Dec)
    precip_scaling_power : float, default 1.0
        Exponent for precipitation scaling (1.0 = linear, >1 = more concentrated
        in wet months, <1 = more uniform). Literature range: 0.5-2.0
    temp_q10 : float, default 2.0
        Q10 factor for temperature scaling. Q10=2.0 means biological reaction
        rates double per 10 deg C increase. Literature range: 1.5-3.0
    temp_reference_c : float, default 15.0
        Reference temperature [deg C] at which the temperature factor equals 1.0

    Returns:
    --------
    list of tuples : [(month, day, hour, load_kg), ...]
        List of 12 monthly temporal load entries
    """
    import math

    if len(monthly_precip_mm) != 12:
        raise ValueError(f"monthly_precip_mm must have 12 values, got {len(monthly_precip_mm)}")
    if len(monthly_temp_c) != 12:
        raise ValueError(f"monthly_temp_c must have 12 values, got {len(monthly_temp_c)}")

    # Compute monthly weights
    monthly_weights = []
    for month_idx in range(12):
        # Precipitation factor: P^alpha (clamp P to minimum 0.1 mm to avoid zero)
        precip = max(monthly_precip_mm[month_idx], 0.1)
        precip_factor = precip ** precip_scaling_power

        # Temperature factor: Q10^((T - Tref) / 10)
        temp = monthly_temp_c[month_idx]
        temp_factor = temp_q10 ** ((temp - temp_reference_c) / 10.0)

        # Combined weight
        monthly_weights.append(precip_factor * temp_factor)

    # Normalize weights so they sum to the annual load
    total_weight = sum(monthly_weights)
    if total_weight <= 0:
        # Fallback to uniform if all weights are zero
        monthly_loads_kg = [annual_load_kg / 12.0] * 12
    else:
        monthly_loads_kg = [(annual_load_kg * w / total_weight) for w in monthly_weights]

    # Create temporal entries (1st of each month at 1:00 AM)
    temporal_entries = []
    for month, load_kg in enumerate(monthly_loads_kg, start=1):
        temporal_entries.append((month, 1, 1, load_kg))

    return temporal_entries


def _read_climate_series(spec: Dict[str, str], kind: str):
    """Read ONE climate variable (precip or temp) into a pandas Series indexed
    by datetime, from a NetCDF or CSV file.

    ``spec`` keys:
      file_type          : 'nc' or 'csv'
      path               : path to the file
      nc_key_or_column   : NetCDF variable name, or CSV column name (the values)
      time_key_or_column : NetCDF time variable/dim, or CSV date column.  This is
                           what aligns the climate to the SIMULATION period, so
                           the right months/years are used.  Optional: defaults
                           to 'time' (NetCDF) or the first date-like CSV column.

    For ``kind='precip'`` the series is returned as a DEPTH per record [mm]:
    rate inputs (e.g. ``kg m-2 s-1`` / ``mm/s`` / ``m/s``) are converted to
    depth using the record's own time step.  For ``kind='temp'`` the series is
    in °C (Kelvin auto-converted).  This is host-model agnostic — the climate is
    supplied as data, not pulled from SUMMA/mizuRoute.
    """
    import os
    import pandas as pd

    ftype = str(spec.get("file_type", "")).strip().lower()
    path = spec.get("path")
    key = spec.get("nc_key_or_column")
    tkey = spec.get("time_key_or_column")      # explicit time axis (optional)
    if not path or not os.path.exists(str(path)):
        raise FileNotFoundError(
            f"ss_climate_data_source[{kind}]: file not found: {path}")
    if not key:
        raise ValueError(
            f"ss_climate_data_source[{kind}]: 'nc_key_or_column' is required.")

    if ftype == "nc" or str(path).lower().endswith(".nc"):
        import xarray as xr
        ds = xr.open_dataset(path)
        try:
            if key not in ds.variables:
                raise KeyError(
                    f"Variable '{key}' not in {path}. "
                    f"Available: {list(ds.data_vars)}")
            da = ds[key]
            # Time dimension: explicit time_key_or_column wins, else 'time',
            # else the variable's first dim.
            if tkey and tkey in da.dims:
                tdim = tkey
            elif "time" in da.dims:
                tdim = "time"
            else:
                tdim = da.dims[0] if da.dims else None
            reduce_dims = [d for d in da.dims if d != tdim]   # basin-average
            if reduce_dims:
                da = da.mean(dim=reduce_dims)
            ser = da.to_series()
            # Index by the named time coordinate when available, so the
            # alignment is explicit rather than positional.
            tcoord = (tkey if (tkey and tkey in ds.variables)
                      else (tdim if (tdim and tdim in ds.variables) else None))
            if tcoord is not None:
                ser.index = pd.to_datetime(ds[tcoord].values)
            units = str(da.attrs.get("units", "")).lower()
        finally:
            ds.close()
    elif ftype == "csv" or str(path).lower().endswith(".csv"):
        df = pd.read_csv(path)
        if tkey and tkey in df.columns:
            dtcol = tkey
        else:
            dtcol = next((c for c in df.columns
                          if str(c).lower() in ("datetime", "date", "time", "timestamp")),
                         df.columns[0])
        if key not in df.columns:
            raise KeyError(
                f"Column '{key}' not in {path}. Available: {list(df.columns)}")
        idx = pd.to_datetime(df[dtcol], errors="coerce")
        ser = pd.Series(pd.to_numeric(df[key], errors="coerce").values,
                        index=idx).dropna()
        units = ""
    else:
        raise ValueError(
            f"ss_climate_data_source[{kind}]: file_type must be 'nc' or 'csv' "
            f"(got '{ftype}').")

    # A usable datetime index is REQUIRED so the monthly aggregation aligns with
    # the simulation period; otherwise tell the user to set time_key_or_column.
    if not isinstance(ser.index, pd.DatetimeIndex):
        try:
            ser.index = pd.to_datetime(ser.index)
        except Exception:
            raise ValueError(
                f"ss_climate_data_source[{kind}]: could not interpret the time "
                f"axis of '{path}'. Set 'time_key_or_column' to the NetCDF time "
                f"variable / CSV date column.")

    ser = ser.sort_index()
    ser = ser[ser.index.notna()]
    # Record time step in seconds (for rate→depth conversion).
    try:
        dt_s = float(pd.Series(ser.index).diff().dropna()
                     .dt.total_seconds().median())
    except Exception:
        dt_s = 3600.0
    if not dt_s or dt_s != dt_s:
        dt_s = 3600.0

    if kind == "precip":
        rate_tokens = ("/s", "s-1", "s^-1", "mm/s", "m/s", "kg m-2 s-1",
                       "kg/m2/s", "kg m**-2 s**-1")
        is_rate = any(tok in units for tok in rate_tokens)
        if not units:                       # no units attr → magnitude heuristic
            try:
                is_rate = float(ser.mean()) < 1e-3   # kg/m2/s ~ 1e-5
            except Exception:
                is_rate = False
        if is_rate:
            ser = ser * dt_s                # rate [.. /s] × step [s] = depth [mm]
        return ser
    else:  # temp
        try:
            if float(ser.mean()) > 60.0:    # clearly Kelvin
                ser = ser - 273.15
        except Exception:
            pass
        return ser


def build_climate_data_from_timeseries(
        climate_source: Dict[str, Dict[str, str]],
        years: Optional[List[int]] = None,
        cache: bool = True) -> Dict[int, Dict[str, List[float]]]:
    """Build ``{year: {'precip_mm':[12], 'temp_c':[12]}}`` from user-supplied
    precipitation + temperature TIME SERIES (NetCDF or CSV).

    This is the data-driven alternative to the hand-entered ``ss_climate_data``
    dict: the dynamic climate-adjusted source/sink loads get their monthly
    precip/temperature straight from a file, keeping openWQ INDEPENDENT of any
    host model (the climate is provided as data, not read from SUMMA/mizuRoute).

    Parameters
    ----------
    climate_source : dict
        ``{'precip': {file_type, path, nc_key_or_column},
            'temp':   {file_type, path, nc_key_or_column}}``.
    years : list[int], optional
        Restrict to these years.  Default: all years present in BOTH series.
    cache : bool
        Cache the computed monthly dict (keyed by source spec) in the temp dir
        and reuse it while the source files are unchanged — the result does not
        depend on calibration parameters, so this avoids re-reading large
        NetCDFs on every evaluation.

    Returns
    -------
    dict : ``{year: {'precip_mm': [12 floats], 'temp_c': [12 floats]}}``
    """
    import os
    import json
    import hashlib
    import tempfile

    if (not isinstance(climate_source, dict)
            or "precip" not in climate_source or "temp" not in climate_source):
        raise ValueError(
            "ss_climate_data_source must be a dict with 'precip' and 'temp' "
            "entries, each {'file_type': 'nc'|'csv', 'path': ..., "
            "'nc_key_or_column': ...}.")
    p_spec = climate_source["precip"]
    t_spec = climate_source["temp"]

    cache_path = None
    if cache:
        try:
            key = json.dumps({"p": p_spec, "t": t_spec,
                              "y": sorted(int(y) for y in years) if years else None},
                             sort_keys=True)
            h = hashlib.md5(key.encode()).hexdigest()[:10]
            cache_path = os.path.join(tempfile.gettempdir(),
                                      f"ss_climate_monthly_{h}.json")
            src_mtime = max(os.path.getmtime(p_spec["path"]),
                            os.path.getmtime(t_spec["path"]))
            if (os.path.exists(cache_path)
                    and os.path.getmtime(cache_path) >= src_mtime):
                with open(cache_path) as f:
                    return {int(k): v for k, v in json.load(f).items()}
        except Exception:
            cache_path = None

    import pandas as pd
    p_ser = _read_climate_series(p_spec, "precip")   # depth per record [mm]
    t_ser = _read_climate_series(t_spec, "temp")     # °C

    p_mon = p_ser.groupby([p_ser.index.year, p_ser.index.month]).sum()
    t_mon = t_ser.groupby([t_ser.index.year, t_ser.index.month]).mean()

    common = sorted({y for (y, _m) in p_mon.index}
                    & {y for (y, _m) in t_mon.index})
    if years:
        want = {int(y) for y in years}
        common = [y for y in common if int(y) in want]

    out: Dict[int, Dict[str, List[float]]] = {}
    for y in common:
        precip = [float(p_mon.get((y, m), 0.0)) for m in range(1, 13)]
        temp = [float(t_mon.get((y, m), float("nan"))) for m in range(1, 13)]
        valid = [v for v in temp if v == v]
        fill = (sum(valid) / len(valid)) if valid else 15.0
        temp = [v if v == v else fill for v in temp]
        out[int(y)] = {"precip_mm": [round(v, 3) for v in precip],
                       "temp_c": [round(v, 3) for v in temp]}

    if not out:
        raise ValueError(
            "build_climate_data_from_timeseries: no overlapping years between the "
            "precip/temp series and the requested simulation period. Check the "
            "ss_climate_data_source paths/keys and the period.")

    if cache_path:
        try:
            with open(cache_path, "w") as f:
                json.dump(out, f)
        except Exception:
            pass
    return out


def build_daily_climate_from_timeseries(
        climate_source: Dict[str, Dict[str, str]],
        years: Optional[List[int]] = None,
        cache: bool = True) -> Dict[int, Dict[str, List]]:
    """Like :func:`build_climate_data_from_timeseries` but at DAILY resolution.

    Used for ``ss_climate_data_type = "time_series"`` so the source/sink load is
    adjusted at EVERY time step (each day) — tracking the actual day-to-day
    precip/temperature — instead of being held constant within a calendar month.

    Returns ``{year: {'month': [...], 'day': [...], 'precip_mm': [...],
    'temp_c': [...]}}`` with one aligned record per day present in the data for
    that year (precip summed per day, temperature averaged per day).
    """
    import os
    import json
    import hashlib
    import tempfile

    if (not isinstance(climate_source, dict)
            or "precip" not in climate_source or "temp" not in climate_source):
        raise ValueError(
            "ss_climate_data_source must be a dict with 'precip' and 'temp' "
            "entries, each {'file_type': 'nc'|'csv', 'path': ..., "
            "'nc_key_or_column': ...}.")
    p_spec = climate_source["precip"]
    t_spec = climate_source["temp"]

    cache_path = None
    if cache:
        try:
            key = json.dumps({"p": p_spec, "t": t_spec,
                              "y": sorted(int(y) for y in years) if years else None,
                              "res": "daily"}, sort_keys=True)
            h = hashlib.md5(key.encode()).hexdigest()[:10]
            cache_path = os.path.join(tempfile.gettempdir(),
                                      f"ss_climate_daily_{h}.json")
            src_mtime = max(os.path.getmtime(p_spec["path"]),
                            os.path.getmtime(t_spec["path"]))
            if (os.path.exists(cache_path)
                    and os.path.getmtime(cache_path) >= src_mtime):
                with open(cache_path) as f:
                    return {int(k): v for k, v in json.load(f).items()}
        except Exception:
            cache_path = None

    import pandas as pd
    p_ser = _read_climate_series(p_spec, "precip")   # depth per record [mm]
    t_ser = _read_climate_series(t_spec, "temp")     # °C

    p_day = p_ser.resample("1D").sum()
    t_day = t_ser.resample("1D").mean()
    df = pd.DataFrame({"precip_mm": p_day, "temp_c": t_day})
    # Keep temperature finite (a NaN day would poison its Q10 weight).
    if df["temp_c"].isna().any():
        _fill = df["temp_c"].mean() if df["temp_c"].notna().any() else 15.0
        df["temp_c"] = df["temp_c"].fillna(_fill)
    df["precip_mm"] = df["precip_mm"].fillna(0.0)

    want = {int(y) for y in years} if years else None
    out: Dict[int, Dict[str, List]] = {}
    for ts, r in df.iterrows():
        y = int(ts.year)
        if want is not None and y not in want:
            continue
        rec = out.setdefault(y, {"month": [], "day": [],
                                 "precip_mm": [], "temp_c": []})
        rec["month"].append(int(ts.month))
        rec["day"].append(int(ts.day))
        rec["precip_mm"].append(round(float(r["precip_mm"]), 4))
        rec["temp_c"].append(round(float(r["temp_c"]), 3))

    if not out:
        raise ValueError(
            "build_daily_climate_from_timeseries: no overlapping days between the "
            "precip/temp series and the requested period. Check the "
            "ss_climate_data_source paths/keys and the period.")

    if cache_path:
        try:
            with open(cache_path, "w") as f:
                json.dump(out, f)
        except Exception:
            pass
    return out


def _apply_climate_adjustment_daily(
        annual_load_kg: float,
        months: List[int],
        days: List[int],
        daily_precip_mm: List[float],
        daily_temp_c: List[float],
        precip_scaling_power: float = 1.0,
        temp_q10: float = 2.0,
        temp_reference_c: float = 15.0) -> List[tuple]:
    """Distribute an annual load across the DAYS of a year by daily climate
    weight: ``w_d = P_d^alpha * Q10^((T_d - T_ref)/10)``; ``L_d = annual *
    w_d / sum(w_d)``.  Returns ``[(month, day, load_kg), ...]`` (one per day),
    summing exactly to ``annual_load_kg``.  If a year has no precip at all the
    load is spread uniformly so the annual total is preserved.
    """
    n = len(daily_precip_mm)
    if n == 0:
        return []
    weights = []
    for i in range(n):
        p = max(float(daily_precip_mm[i]), 0.0)
        t = float(daily_temp_c[i])
        weights.append((p ** precip_scaling_power)
                       * (temp_q10 ** ((t - temp_reference_c) / 10.0)))
    tot = sum(weights)
    if tot <= 0:
        share = annual_load_kg / n
        return [(int(months[i]), int(days[i]), share) for i in range(n)]
    return [(int(months[i]), int(days[i]), annual_load_kg * weights[i] / tot)
            for i in range(n)]


def build_subannual_climate_from_timeseries(
        climate_source: Dict[str, Dict[str, str]],
        years: Optional[List[int]] = None,
        resolution: str = "hourly",
        cache: bool = True) -> Dict[int, Dict[str, List]]:
    """Like :func:`build_daily_climate_from_timeseries` but at HOURLY resolution,
    or at the file's NATIVE step (``resolution='native'`` → no resampling).

    Used by ``ss_climate_data_source_adjusting_resolution in ('hourly','native')``
    so the source/sink load is adjusted at every (sub-daily) time step.

    Returns ``{year: {'month':[], 'day':[], 'hour':[], 'precip_mm':[],
    'temp_c':[]}}`` with one aligned record per period present for that year
    (precip summed, temperature averaged).
    """
    import os
    import json
    import hashlib
    import tempfile

    if (not isinstance(climate_source, dict)
            or "precip" not in climate_source or "temp" not in climate_source):
        raise ValueError(
            "ss_climate_data_source must be a dict with 'precip' and 'temp' "
            "entries, each {'file_type': 'nc'|'csv', 'path': ..., "
            "'nc_key_or_column': ...}.")
    p_spec = climate_source["precip"]
    t_spec = climate_source["temp"]
    res = str(resolution).lower()

    cache_path = None
    if cache:
        try:
            key = json.dumps({"p": p_spec, "t": t_spec,
                              "y": sorted(int(y) for y in years) if years else None,
                              "res": res}, sort_keys=True)
            h = hashlib.md5(key.encode()).hexdigest()[:10]
            cache_path = os.path.join(tempfile.gettempdir(),
                                      f"ss_climate_{res}_{h}.json")
            src_mtime = max(os.path.getmtime(p_spec["path"]),
                            os.path.getmtime(t_spec["path"]))
            if (os.path.exists(cache_path)
                    and os.path.getmtime(cache_path) >= src_mtime):
                with open(cache_path) as f:
                    return {int(k): v for k, v in json.load(f).items()}
        except Exception:
            cache_path = None

    import pandas as pd
    p_ser = _read_climate_series(p_spec, "precip")   # depth per record [mm]
    t_ser = _read_climate_series(t_spec, "temp")     # °C

    if res == "native":
        # No resampling — use the file's own timestamps as-is.
        df = pd.DataFrame({"precip_mm": p_ser, "temp_c": t_ser})
    else:  # hourly
        df = pd.DataFrame({"precip_mm": p_ser.resample("1h").sum(),
                           "temp_c": t_ser.resample("1h").mean()})
    if df["temp_c"].isna().any():
        _fill = df["temp_c"].mean() if df["temp_c"].notna().any() else 15.0
        df["temp_c"] = df["temp_c"].fillna(_fill)
    df["precip_mm"] = df["precip_mm"].fillna(0.0)

    want = {int(y) for y in years} if years else None
    out: Dict[int, Dict[str, List]] = {}
    for ts, r in df.iterrows():
        y = int(ts.year)
        if want is not None and y not in want:
            continue
        rec = out.setdefault(y, {"month": [], "day": [], "hour": [],
                                 "precip_mm": [], "temp_c": []})
        rec["month"].append(int(ts.month))
        rec["day"].append(int(ts.day))
        rec["hour"].append(int(ts.hour))
        rec["precip_mm"].append(round(float(r["precip_mm"]), 4))
        rec["temp_c"].append(round(float(r["temp_c"]), 3))

    if not out:
        raise ValueError(
            "build_subannual_climate_from_timeseries: no overlapping periods "
            "between the precip/temp series and the requested period. Check the "
            "ss_climate_data_source paths/keys and the period.")

    if cache_path:
        try:
            with open(cache_path, "w") as f:
                json.dump(out, f)
        except Exception:
            pass
    return out


def _apply_climate_adjustment_periodic(
        annual_load_kg: float,
        months: List[int],
        days: List[int],
        hours: List[int],
        precip_mm: List[float],
        temp_c: List[float],
        precip_scaling_power: float = 1.0,
        temp_q10: float = 2.0,
        temp_reference_c: float = 15.0) -> List[tuple]:
    """Distribute an annual load across sub-daily periods by climate weight
    ``w = P^alpha * Q10^((T - T_ref)/10)``.  Returns
    ``[(month, day, hour, load_kg), ...]`` summing exactly to
    ``annual_load_kg`` (uniform fallback when the year has no precip weight).
    """
    n = len(precip_mm)
    if n == 0:
        return []
    weights = []
    for i in range(n):
        p = max(float(precip_mm[i]), 0.0)
        t = float(temp_c[i])
        weights.append((p ** precip_scaling_power)
                       * (temp_q10 ** ((t - temp_reference_c) / 10.0)))
    tot = sum(weights)
    if tot <= 0:
        share = annual_load_kg / n
        return [(int(months[i]), int(days[i]), int(hours[i]), share)
                for i in range(n)]
    return [(int(months[i]), int(days[i]), int(hours[i]),
             annual_load_kg * weights[i] / tot) for i in range(n)]


def _expand_climate_data_all(climate_data, period):
    """Expand a climate-data dict that uses the convenience key ``"all"`` into
    explicit per-year entries spanning the simulation ``period``
    (``[year_start, year_end]``).

    Lets the user write::

        ss_climate_data = {"all": {"precip_mm": [...12...],
                                   "temp_c":    [...12...]}}

    instead of repeating the same monthly cycle for every year &mdash; the same
    coefficients are then applied to ALL simulation years.  Any explicit
    per-year entries already present win; every remaining year in the period
    falls back to the ``"all"`` cycle.  Returns a dict keyed by ``int`` year
    (the ``"all"`` key is removed) so downstream code that sorts / casts the
    keys keeps working.  A dict without an ``"all"`` key is returned unchanged.
    """
    if not isinstance(climate_data, dict) or not climate_data:
        return climate_data
    all_key = next((k for k in climate_data
                    if isinstance(k, str) and k.strip().lower() == "all"), None)
    if all_key is None:
        return climate_data
    all_cycle = climate_data[all_key]
    expanded = {}
    for k, v in climate_data.items():
        if k == all_key:
            continue
        try:
            expanded[int(k)] = v
        except (ValueError, TypeError):
            expanded[k] = v
    try:
        y0, y1 = int(period[0]), int(period[1])
    except (TypeError, ValueError, IndexError):
        y0 = y1 = None
    if y0 is not None and y1 is not None:
        for yr in range(min(y0, y1), max(y0, y1) + 1):
            expanded.setdefault(yr, all_cycle)
    print(f"  ss_climate_data: 'all' cycle applied to "
          f"{sorted(k for k in expanded if isinstance(k, int))}")
    return expanded


def set_ss_climate_adjusted_export_coefficients(
        ss_config_filepath: str,
        json_header_comment: List[str],
        ss_method_copernicus_basins_hrus: Dict[str, str],
        ss_method_copernicus_nc_lc_dir: Optional[str],
        ss_method_copernicus_period: List[Union[int, float]],
        ss_method_copernicus_default_loads_bool: bool,
        ss_method_copernicus_compartment_name_for_load: str,
        climate_data: Dict[int, Dict[str, List[float]]],
        precip_scaling_power: float = 1.0,
        temp_q10: float = 2.0,
        temp_reference_c: float = 15.0,
        ss_metadata_comment: str = "Climate-adjusted nutrient loading from LULC",
        ss_metadata_source: str = "Copernicus LULC + climate adjustment",
        optional_load_coefficients: Optional[Dict[int, Dict[str, float]]] = None,
        recursive: bool = False,
        file_pattern: str = 'ESACCI-LC-*.nc',
        simulation_period: Optional[List[int]] = None,
        bgc_engine_label: str = "",
        chemical_species_list: Optional[List[str]] = None,
        # Where the monthly climate came from (for the report metadata):
        #   climate_data_type   = "fixed_parameters" | "time_series"
        #   climate_data_source = the ss_climate_data_source dict (time_series only)
        climate_data_type: str = "fixed_parameters",
        climate_data_source: Optional[Dict[str, Dict[str, str]]] = None,
        # When provided (time_series mode), the annual load is distributed at
        # DAILY resolution from this {year:{month,day,precip_mm,temp_c}} dict —
        # one SS entry per day — so the load is adjusted at every time step
        # instead of being constant within a month.
        daily_climate_data: Optional[Dict[int, Dict[str, List]]] = None,
        # Sub-daily (hourly / native) climate: {year:{month,day,hour,precip_mm,
        # temp_c}} — one SS entry per sub-daily period, unit = subannual_unit.
        subannual_climate_data: Optional[Dict[int, Dict[str, List]]] = None,
        subannual_unit: str = "hour"
) -> None:
    """
    Generate climate-adjusted source/sink JSON from Copernicus LULC data.

    This extends the Copernicus LULC method by modulating the static export
    coefficients with monthly precipitation and temperature data. This approach
    is inspired by SWAT's temporal nutrient dynamics and produces more realistic
    seasonal loading patterns than simple uniform or sinusoidal distributions.

    The monthly load for each HRU and nutrient is computed as:
        L_m = L_annual * w_m / sum(w_m)
    where the monthly weight w_m combines precipitation and temperature effects:
        w_m = P_m^alpha * Q10^((T_m - T_ref) / 10)

    Parameters:
    -----------
    ss_config_filepath : str
        Path where JSON configuration will be saved
    json_header_comment : List[str]
        Comment lines for JSON header
    ss_method_copernicus_basins_hrus : dict
        Shapefile info: {'path_to_shp': ..., 'mapping_key': ...}
    ss_method_copernicus_nc_lc_dir : str
        Directory with Copernicus LULC NetCDF files
    ss_method_copernicus_period : list
        [year_start, year_end]
    ss_method_copernicus_default_loads_bool : bool
        If True, use default export coefficients; if False, use optional_load_coefficients
    ss_method_copernicus_compartment_name_for_load : str
        OpenWQ compartment name (e.g., 'RIVER_NETWORK_REACHES')
    climate_data : dict
        Monthly climate data indexed by year. Format:
        {
            2020: {
                'precip_mm': [p1, p2, ..., p12],  # monthly precipitation [mm]
                'temp_c': [t1, t2, ..., t12]       # monthly mean temperature [deg C]
            },
            2021: { ... }
        }
        If a year from the LULC period is missing, the nearest available year is used.
    precip_scaling_power : float, default 1.0
        Exponent for precipitation scaling (1.0=linear, >1=more concentrated in wet months)
    temp_q10 : float, default 2.0
        Q10 temperature scaling factor (rate doubling per 10 deg C increase)
    temp_reference_c : float, default 15.0
        Reference temperature [deg C] at which temperature factor = 1.0
    ss_metadata_comment : str
        Comment for JSON metadata section
    ss_metadata_source : str
        Source identifier for JSON metadata section
    optional_load_coefficients : dict, optional
        Custom load coefficients (used if ss_method_copernicus_default_loads_bool=False)
    recursive : bool
        Whether to search subdirectories for NetCDF files
    file_pattern : str
        Glob pattern for NetCDF files

    """

    print("\n" + "=" * 60)
    print("CLIMATE-ADJUSTED EXPORT COEFFICIENT METHOD")
    print("=" * 60)
    print(f"Precipitation scaling power (alpha): {precip_scaling_power}")
    print(f"Temperature Q10: {temp_q10}")
    print(f"Reference temperature: {temp_reference_c} deg C")
    # Allow ss_climate_data = {"all": {...}} as shorthand for "apply this one
    # monthly cycle to every simulation year"; expand it to explicit per-year
    # entries so the rest of the pipeline (and the per-year lookup) is unchanged.
    climate_data = _expand_climate_data_all(
        climate_data, ss_method_copernicus_period)
    print(f"Climate data years available: "
          f"{sorted((climate_data or {}).keys(), key=str)}")

    # Step 1: Calculate LULC areas (reuses existing Copernicus pipeline)
    results, summaries, rasters = set_ss_from_copernicus_lulc(
        ss_config_filepath=ss_config_filepath,
        ss_method_copernicus_basins_hrus=ss_method_copernicus_basins_hrus,
        ss_method_copernicus_nc_lc_dir=ss_method_copernicus_nc_lc_dir,
        ss_method_copernicus_period=ss_method_copernicus_period,
        recursive=recursive,
        file_pattern=file_pattern
    )

    # Step 2: Determine export coefficients
    if ss_method_copernicus_default_loads_bool:
        print("\nUsing DEFAULT Copernicus land cover nutrient coefficients...")
        load_coefficients = get_default_copernicus_load_coefficients()
    else:
        if optional_load_coefficients is None:
            raise ValueError(
                "ss_method_copernicus_default_loads_bool=False but no optional_load_coefficients provided."
            )
        print("\nUsing CUSTOM load coefficients...")
        load_coefficients = optional_load_coefficients

    # Step 3: Calculate annual loads per HRU
    output_dir = Path(ss_config_filepath).parent / 'ss_copernicus_files'
    shp_hru_id_column = ss_method_copernicus_basins_hrus.get('mapping_key', 'HRU_ID')

    loads_df, pivot_loads_df = calculate_nutrient_loads(
        results_df=results,
        load_coefficients=load_coefficients,
        output_dir=output_dir,
        shp_hru_id_column=shp_hru_id_column
    )

    # Step 4: Build climate-adjusted JSON
    pivot_csv_path = output_dir / 'nutrient_loads_pivot.csv'
    pivot_df = pd.read_csv(pivot_csv_path)

    exclude_cols = [shp_hru_id_column, 'Year']
    nutrient_cols = [col for col in pivot_df.columns if col not in exclude_cols]

    # Filter to only species that exist in the BGC template
    if chemical_species_list:
        _skipped = [c for c in nutrient_cols if c not in chemical_species_list]
        nutrient_cols = [c for c in nutrient_cols if c in chemical_species_list]
        if _skipped:
            print(f"\n  ⚠ Skipping SS species not in BGC template: {', '.join(_skipped)}")

    cop_years_in_csv = sorted(int(y) for y in pivot_df['Year'].unique())

    # No METADATA key in JSON body — C++ parser counts all top-level keys
    # via .size(), so METADATA would be treated as a numbered entry.
    # Metadata goes into // comment lines instead (C++ skips with skip_comments=true).
    config = {}

    metadata_comments = [
        f"// METADATA - Comment: {ss_metadata_comment}",
        f"// METADATA - Source: {ss_metadata_source}",
        f"// METADATA - BGC_Engine: {bgc_engine_label}",
        f"// METADATA - Chemical_Species: {', '.join(chemical_species_list) if chemical_species_list else ''}",
        f"// METADATA - Climate_Adjustment: precip_scaling_power={precip_scaling_power}, "
        f"temp_Q10={temp_q10}, temp_reference_C={temp_reference_c}",
        f"// METADATA - Climate_Data_Type: {climate_data_type}",
        f"// METADATA - Climate_Resolution: "
        f"{(subannual_unit + 'ly (' + subannual_unit + '-step)') if subannual_climate_data else ('daily (per-time-step)' if daily_climate_data else 'monthly')}",
    ]
    # When the climate comes from a time series, record which file(s)/variables
    # were read, so the report shows it (kept colon-free so the report's
    # "// METADATA - Key: Value" parser splits cleanly).
    if str(climate_data_type).lower() == "time_series" and isinstance(climate_data_source, dict):
        _src_parts = []
        for _kind in ("precip", "temp"):
            _s = climate_data_source.get(_kind) or {}
            if _s:
                _src_parts.append(
                    f"{_kind} = {_s.get('path', '?')} "
                    f"[{_s.get('file_type', '?')}, var={_s.get('nc_key_or_column', '?')}, "
                    f"time={_s.get('time_key_or_column', 'time')}]")
        if _src_parts:
            metadata_comments.append(
                "// METADATA - Climate_Data_Source: " + "; ".join(_src_parts))

    # Helper: find nearest climate data year
    climate_years = sorted(climate_data.keys())

    def _get_climate_for_year(yr):
        if yr in climate_data:
            return climate_data[yr]
        # Find nearest year
        nearest = min(climate_years, key=lambda y: abs(y - yr))
        print(f"    (using climate data from {nearest} for year {yr})")
        return climate_data[nearest]

    # Daily-resolution lookup (time_series mode): one record-set per day.
    _use_daily = bool(daily_climate_data)
    _daily_years = sorted(int(k) for k in daily_climate_data.keys()) if _use_daily else []

    def _get_daily_climate_for_year(yr):
        if yr in daily_climate_data:
            return daily_climate_data[yr]
        nearest = min(_daily_years, key=lambda y: abs(y - yr))
        return daily_climate_data[nearest]

    # Sub-daily (hourly / native) lookup: one record-set per sub-daily period.
    _use_subannual = bool(subannual_climate_data)
    _sub_years = (sorted(int(k) for k in subannual_climate_data.keys())
                  if _use_subannual else [])

    def _get_subannual_climate_for_year(yr):
        if yr in subannual_climate_data:
            return subannual_climate_data[yr]
        nearest = min(_sub_years, key=lambda y: abs(y - yr))
        return subannual_climate_data[nearest]

    # Build year mapping: simulation years → Copernicus data years
    if simulation_period is not None:
        sim_start, sim_end = int(simulation_period[0]), int(simulation_period[1])
        sim_years = list(range(sim_start, sim_end + 1))
    else:
        sim_years = cop_years_in_csv

    year_mapping = {}
    for sy in sim_years:
        year_mapping[sy] = min(cop_years_in_csv, key=lambda cy: abs(cy - sy))

    proxy_years = {sy: cy for sy, cy in year_mapping.items() if sy != cy}
    if proxy_years:
        print(f"\n  ⚠  Proxy year mapping (sim year → Copernicus data year):")
        for sy in sorted(proxy_years):
            print(f"      {sy} → {proxy_years[sy]}")
        # Add proxy year warning as a comment line
        proxy_str = ', '.join(f'{sy}→{cy}' for sy, cy in sorted(proxy_years.items()))
        metadata_comments.append(
            f"// METADATA - Proxy_Years: {proxy_str}"
        )

    entry_idx = 1

    for sim_year in sim_years:
        cop_year = year_mapping[sim_year]
        is_proxy = (sim_year != cop_year)
        year_data = pivot_df[pivot_df['Year'] == cop_year]
        climate = _get_climate_for_year(sim_year)

        for nutrient in nutrient_cols:
            proxy_tag = f" (proxy→{cop_year})" if is_proxy else ""
            print(f"\n  Processing: Year {sim_year}{proxy_tag}, Nutrient {nutrient}")

            nutrient_loads = year_data[[shp_hru_id_column, nutrient]].copy()
            nutrient_loads = nutrient_loads[nutrient_loads[nutrient] > 0]

            if len(nutrient_loads) == 0:
                print(f"    No non-zero loads, skipping...")
                continue

            data_entries = {}
            sub_idx = 1

            for _, row in nutrient_loads.iterrows():
                hru_id = str(int(row[shp_hru_id_column]))
                annual_load_kg = row[nutrient]

                # Use cell_id directly - OpenWQ C++ will do the lookup
                spatial_id = hru_id

                if _use_subannual:
                    # HOURLY / NATIVE (time_series) — one SS entry per sub-daily
                    # period, scaled by that period's precip + temperature.  The
                    # entry pins month+day+hour and carries unit=subannual_unit.
                    sc = _get_subannual_climate_for_year(sim_year)
                    for month, day, hour, prd_load_kg in _apply_climate_adjustment_periodic(
                            annual_load_kg=annual_load_kg,
                            months=sc['month'], days=sc['day'], hours=sc['hour'],
                            precip_mm=sc['precip_mm'], temp_c=sc['temp_c'],
                            precip_scaling_power=precip_scaling_power,
                            temp_q10=temp_q10,
                            temp_reference_c=temp_reference_c):
                        data_entries[str(sub_idx)] = [
                            sim_year, int(month), int(day), int(hour), "all", "all",
                            spatial_id, "all", "all",
                            float(prd_load_kg),
                            "continuous", subannual_unit
                        ]
                        sub_idx += 1
                    continue

                if _use_daily:
                    # DAILY (time_series) — one SS entry per day, scaled by that
                    # day's precip + temperature, so the load is adjusted at
                    # every time step rather than constant within a month.
                    dc = _get_daily_climate_for_year(sim_year)
                    for month, day, day_load_kg in _apply_climate_adjustment_daily(
                            annual_load_kg=annual_load_kg,
                            months=dc['month'], days=dc['day'],
                            daily_precip_mm=dc['precip_mm'],
                            daily_temp_c=dc['temp_c'],
                            precip_scaling_power=precip_scaling_power,
                            temp_q10=temp_q10,
                            temp_reference_c=temp_reference_c):
                        data_entries[str(sub_idx)] = [
                            sim_year, int(month), int(day), "all", "all", "all",
                            spatial_id, "all", "all",
                            float(day_load_kg),
                            "continuous", "day"
                        ]
                        sub_idx += 1
                    continue

                # MONTHLY (fixed_parameters) — constant daily rate within month.
                temporal_entries = _apply_climate_adjustment(
                    annual_load_kg=annual_load_kg,
                    year=sim_year,
                    monthly_precip_mm=climate['precip_mm'],
                    monthly_temp_c=climate['temp_c'],
                    precip_scaling_power=precip_scaling_power,
                    temp_q10=temp_q10,
                    temp_reference_c=temp_reference_c
                )

                for month, day, hour, load_kg in temporal_entries:
                    days = calendar.monthrange(sim_year, int(month))[1]
                    daily_rate = load_kg / days
                    data_entries[str(sub_idx)] = [
                            sim_year, int(month), "all", "all", "all", "all",
                            spatial_id, "all", "all",
                            float(daily_rate),
                            "continuous", "day"
                        ]
                    sub_idx += 1

            if len(data_entries) == 0:
                continue

            # Build COMMENT — note proxy source when year is outside Copernicus range
            if is_proxy:
                comment = (f"Climate-adjusted loading for {nutrient} in {sim_year} "
                           f"(proxy from Copernicus year {cop_year})")
            else:
                comment = f"Climate-adjusted loading for {nutrient} in {sim_year}"

            config[str(entry_idx)] = {
                "CHEMICAL_NAME": nutrient,
                "COMPARTMENT_NAME": ss_method_copernicus_compartment_name_for_load,
                "COMMENT": comment,
                "TYPE": "source",
                "UNITS": "kg",
                "DATA_FORMAT": "JSON",
                "DATA": data_entries
            }

            print(f"    Added {len(data_entries)} climate-adjusted entries")
            entry_idx += 1

    # Write JSON
    json_string = json.dumps(config, indent=4)

    def compress_array(match):
        array_content = match.group(0)
        compressed = re.sub(r'\s+', ' ', array_content)
        compressed = re.sub(r'\[\s+', '[', compressed)
        compressed = re.sub(r'\s+\]', ']', compressed)
        compressed = re.sub(r'\s*,\s*', ', ', compressed)
        return compressed

    json_string = re.sub(r'\[[^\[\]]*\]', compress_array, json_string)

    output_path = Path(ss_config_filepath)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(ss_config_filepath, 'w') as f:
        for comment in json_header_comment:
            f.write(comment + "\n")
        # Write metadata as comment lines (C++ parser ignores these)
        for mc in metadata_comments:
            f.write(mc + "\n")
        f.write(json_string)
        f.write("\n")

    print("\n" + "=" * 60)
    print("CLIMATE-ADJUSTED SS JSON CREATED!")
    print("=" * 60)
    print(f"File: {ss_config_filepath}")
    print(f"Total source entries: {entry_idx - 1}")
    print(f"Climate parameters: alpha={precip_scaling_power}, Q10={temp_q10}, Tref={temp_reference_c}")

def create_openwq_ss_json_from_loads(
        pivot_loads_csv_path: Path,
        ss_config_filepath: str,
        json_header_comment: List[str],
        ss_method_copernicus_compartment_name_for_load: str,
        shp_hru_id_column: str,
        ss_method_copernicus_annual_to_seasonal_loads_method: str = "uniform",
        ss_metadata_comment: str = "Nutrient loading from land use/land cover",
        ss_metadata_source: str = "Copernicus LULC analysis",
        simulation_period: Optional[List[int]] = None,
        bgc_engine_label: str = "",
        chemical_species_list: Optional[List[str]] = None,
        climate_data: Optional[Dict] = None,
        precip_scaling_power: float = 1.0,
        temp_q10: float = 2.0,
        temp_reference_c: float = 15.0
) -> None:
    """
    Generate OpenWQ source/sink JSON configuration from nutrient loads.

    Parameters:
    -----------
    pivot_loads_csv_path : Path
        Path to the nutrient_loads_pivot.csv file
    ss_config_filepath : str
        Full path where the JSON file will be saved
    json_header_comment : List[str]
        List of comment lines to add at the top of the file
    ss_method_copernicus_compartment_name_for_load : str
        Name of the compartment (e.g., 'RIVER_NETWORK_REACHES')
    ss_metadata_comment : str
        Comment for metadata section
    ss_metadata_source : str
        Source identifier for metadata section

    """

    # Expand ss_climate_data = {"all": {...}} (and mixed int/"all" keys) into
    # explicit per-year entries up front, so the prints + per-year lookups below
    # never sort/compare a str key ("all") against int years.
    climate_data = _expand_climate_data_all(climate_data, simulation_period)

    print("\n" + "=" * 60)
    print("GENERATING OPENWQ SOURCE/SINK JSON")
    print("=" * 60)

    # Validate temporal distribution method. "harmonic"/"monthly" generate a
    # FLAT (uniform) monthly base here; their seasonal shape is imposed later by
    # the calibration (per-month load multipliers), so generation treats them as
    # uniform.
    valid_methods = ["uniform", "seasonal", "harmonic", "monthly"]
    if ss_method_copernicus_annual_to_seasonal_loads_method not in valid_methods:
        raise ValueError(
            f"Invalid temporal distribution method: '{ss_method_copernicus_annual_to_seasonal_loads_method}'. "
            f"Valid options are: {valid_methods}"
        )

    print(f"\nTemporal distribution method: {ss_method_copernicus_annual_to_seasonal_loads_method}")
    print("  → Load type: continuous (kg/day rate per month)")
    if ss_method_copernicus_annual_to_seasonal_loads_method == "uniform":
        print("  → Annual loads distributed uniformly across 12 months")
    elif ss_method_copernicus_annual_to_seasonal_loads_method in ("harmonic", "monthly"):
        print("  → Flat monthly base; seasonal shape set by calibration "
              f"({ss_method_copernicus_annual_to_seasonal_loads_method})")
    elif ss_method_copernicus_annual_to_seasonal_loads_method == "seasonal":
        if climate_data:
            print("  → Annual loads climate-weighted (precipitation + temperature)")
            print(f"  → Climate data provided for years: {sorted(climate_data.keys(), key=str)}")
            print(f"  → Parameters: alpha={precip_scaling_power}, Q10={temp_q10}, Tref={temp_reference_c}°C")
        else:
            print("  → Annual loads follow seasonal pattern (sine-curve fallback)")
            print("  → TIP: provide ss_climate_data for data-driven distribution")

    # Load pivot loads data
    print(f"\nLoading pivot loads from: {pivot_loads_csv_path.name}")
    pivot_df = pd.read_csv(pivot_loads_csv_path)

    # Get nutrient columns (exclude HRU_ID and Year columns)
    exclude_cols = [shp_hru_id_column, 'Year']
    nutrient_cols = [col for col in pivot_df.columns if col not in exclude_cols]

    # Filter to only species that exist in the BGC template
    if chemical_species_list:
        _skipped = [c for c in nutrient_cols if c not in chemical_species_list]
        nutrient_cols = [c for c in nutrient_cols if c in chemical_species_list]
        if _skipped:
            print(f"\n  ⚠ Skipping SS species not in BGC template: {', '.join(_skipped)}")

    print(f"\nNutrient species found: {nutrient_cols}")
    cop_years_in_csv = sorted(int(y) for y in pivot_df['Year'].unique())
    print(f"Copernicus years in CSV: {cop_years_in_csv}")

    # Build year mapping: simulation years → Copernicus data years
    # Years inside the Copernicus range use actual data; years outside use
    # the nearest available Copernicus year as proxy.
    if simulation_period is not None:
        sim_start, sim_end = int(simulation_period[0]), int(simulation_period[1])
        sim_years = list(range(sim_start, sim_end + 1))
    else:
        sim_years = cop_years_in_csv  # no remapping needed

    year_mapping = {}  # {sim_year: cop_year}
    for sy in sim_years:
        year_mapping[sy] = min(cop_years_in_csv, key=lambda cy: abs(cy - sy))

    # Log proxy years
    proxy_years = {sy: cy for sy, cy in year_mapping.items() if sy != cy}
    if proxy_years:
        print(f"\n  ⚠  Proxy year mapping (sim year → Copernicus data year):")
        for sy in sorted(proxy_years):
            print(f"      {sy} → {proxy_years[sy]}")
    print(f"  Simulation years to write: {sim_years[0]}–{sim_years[-1]}")

    # No METADATA key in JSON body — C++ parser counts all top-level keys
    # via .size(), so METADATA would be treated as a numbered entry.
    # Metadata goes into // comment lines instead (C++ skips with skip_comments=true).
    config = {}

    metadata_comments = [
        f"// METADATA - Comment: {ss_metadata_comment}",
        f"// METADATA - Source: {ss_metadata_source}",
        f"// METADATA - BGC_Engine: {bgc_engine_label}",
        f"// METADATA - Chemical_Species: {', '.join(chemical_species_list) if chemical_species_list else ''}",
    ]

    # Add proxy year warning as a comment line
    if proxy_years:
        proxy_str = ', '.join(f'{sy}→{cy}' for sy, cy in sorted(proxy_years.items()))
        metadata_comments.append(
            f"// METADATA - Proxy_Years: {proxy_str}"
        )

    entry_idx = 1

    # Group by simulation year and process each nutrient
    for sim_year in sim_years:
        cop_year = year_mapping[sim_year]
        is_proxy = (sim_year != cop_year)
        year_data = pivot_df[pivot_df['Year'] == cop_year]

        for nutrient in nutrient_cols:
            proxy_tag = f" (proxy→{cop_year})" if is_proxy else ""
            print(f"\n  Processing: Year {sim_year}{proxy_tag}, Nutrient {nutrient}")

            # Get non-zero loads for this Copernicus year and nutrient
            nutrient_loads = year_data[[shp_hru_id_column, nutrient]].copy()
            nutrient_loads = nutrient_loads[nutrient_loads[nutrient] > 0]  # Only non-zero loads

            if len(nutrient_loads) == 0:
                print(f"    ⚠ No non-zero loads for {nutrient} in {cop_year}, skipping...")
                continue

            # Create data entries
            data_entries = {}
            sub_idx = 1

            for _, row in nutrient_loads.iterrows():
                hru_id = row[shp_hru_id_column]
                hru_id = str(int(hru_id))
                annual_load_kg = row[nutrient]

                # Use cell_id directly - OpenWQ C++ will do the lookup
                spatial_id = hru_id

                # Calculate temporal distribution of annual load
                temporal_entries = _calculate_temporal_load_distribution(
                    annual_load_kg=annual_load_kg,
                    year=sim_year,
                    method=ss_method_copernicus_annual_to_seasonal_loads_method,
                    climate_data=climate_data,
                    precip_scaling_power=precip_scaling_power,
                    temp_q10=temp_q10,
                    temp_reference_c=temp_reference_c
                )

                # Create data entries for each time point
                # Format: [YYYY, MM, "all","all","all","all", cell_id, iy, iz,
                #          daily_rate, "continuous", "day"]
                # YYYY is the simulation year (not the Copernicus year)
                for month, day, hour, load_kg in temporal_entries:
                    days = calendar.monthrange(sim_year, int(month))[1]
                    daily_rate = load_kg / days
                    data_entries[str(sub_idx)] = [
                            sim_year, int(month), "all", "all", "all", "all",
                            spatial_id, "all", "all",
                            float(daily_rate),
                            "continuous", "day"
                        ]
                    sub_idx += 1

            if len(data_entries) == 0:
                print(f"    ⚠ No valid entries for {nutrient} in {sim_year}, skipping...")
                continue

            # Build COMMENT — note proxy source when year is outside Copernicus range
            if is_proxy:
                comment = (f"Copernicus LULC-based loading for {nutrient} in {sim_year} "
                           f"(proxy from Copernicus year {cop_year})")
            else:
                comment = f"Copernicus LULC-based loading for {nutrient} in {sim_year}"

            # Add to config
            config[str(entry_idx)] = {
                "CHEMICAL_NAME": nutrient,
                "COMPARTMENT_NAME": ss_method_copernicus_compartment_name_for_load,
                "COMMENT": comment,
                "TYPE": "source",
                "UNITS": "kg",
                "DATA_FORMAT": "JSON",
                "DATA": data_entries
            }

            print(f"    Added {len(data_entries)} load entries for {nutrient} in {sim_year}")
            entry_idx += 1

    # Convert to JSON string with standard formatting
    json_string = json.dumps(config, indent=4)

    # Custom formatting function
    def compress_array(match):
        array_content = match.group(0)
        compressed = re.sub(r'\s+', ' ', array_content)
        compressed = re.sub(r'\[\s+', '[', compressed)
        compressed = re.sub(r'\s+\]', ']', compressed)
        compressed = re.sub(r'\s*,\s*', ', ', compressed)
        return compressed

    # Apply compression to all arrays
    json_string = re.sub(r'\[[^\[\]]*\]', compress_array, json_string)

    # Create output directory if needed
    output_path = Path(ss_config_filepath)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Write to file
    with open(ss_config_filepath, 'w') as f:
        # Write comment lines first
        for comment in json_header_comment:
            f.write(comment + "\n")
        # Write metadata as comment lines (C++ parser ignores these)
        for mc in metadata_comments:
            f.write(mc + "\n")
        # Write the JSON content
        f.write(json_string)
        # Add newline at end of file
        f.write("\n")

    print("\n" + "=" * 60)
    print("✓ OpenWQ SOURCE/SINK JSON CREATED!")
    print("=" * 60)
    print(f"File saved to: {ss_config_filepath}")
    print(f"Total source entries: {entry_idx - 1}")
    print(f"Nutrients: {nutrient_cols}")
    print(f"Compartment: {ss_method_copernicus_compartment_name_for_load}")


def get_default_copernicus_load_coefficients() -> Dict[int, Dict[str, float]]:
    """
    Get default nutrient load coefficients for Copernicus land cover classes.

    Based on ESA CCI Land Cover classification and literature values.
    Coefficients are in kg/ha/yr for Total Nitrogen (TN) and Total Phosphorus (TP).

    Returns:
    --------
    dict : Load coefficients by LC class

    References:
    -----------
    See SS_lulc_loads_reference.md for sources and details.
    """
    return  {
        # Cropland classes (Row-crop / cultivated agriculture)
        10: {'TN': 14.0, 'TP': 2.0, 'NO3-N': 10.0, 'NO2-N': 0.05, 'NH4-N': 0.6, 'PO4-P': 0.3},   # Cropland, rainfed
        11: {'TN': 14.0, 'TP': 2.0, 'NO3-N': 10.0, 'NO2-N': 0.05, 'NH4-N': 0.6, 'PO4-P': 0.3},   # Herbaceous cover
        12: {'TN': 12.0, 'TP': 1.5, 'NO3-N': 8.0, 'NO2-N': 0.04, 'NH4-N': 0.5, 'PO4-P': 0.25},   # Tree or shrub cover
        20: {'TN': 18.0, 'TP': 2.5, 'NO3-N': 13.0, 'NO2-N': 0.06, 'NH4-N': 0.8, 'PO4-P': 0.4},   # Cropland, irrigated or post-flooding
        30: {'TN': 12.0, 'TP': 1.2, 'NO3-N': 8.0, 'NO2-N': 0.04, 'NH4-N': 0.5, 'PO4-P': 0.2},    # Mosaic cropland (>50%) / natural vegetation
        40: {'TN': 8.0, 'TP': 0.6, 'NO3-N': 5.0, 'NO2-N': 0.03, 'NH4-N': 0.3, 'PO4-P': 0.1},     # Mosaic natural vegetation (>50%) / cropland (mixed rural)

        # Forest classes (Forest / natural catchment)
        50: {'TN': 3.0, 'TP': 0.15, 'NO3-N': 1.5, 'NO2-N': 0.005, 'NH4-N': 0.15, 'PO4-P': 0.01},   # Tree cover, broadleaved, evergreen
        60: {'TN': 3.0, 'TP': 0.15, 'NO3-N': 1.5, 'NO2-N': 0.005, 'NH4-N': 0.15, 'PO4-P': 0.01},   # Tree cover, broadleaved, deciduous
        61: {'TN': 3.0, 'TP': 0.15, 'NO3-N': 1.5, 'NO2-N': 0.005, 'NH4-N': 0.15, 'PO4-P': 0.01},   # Tree cover, broadleaved, deciduous, closed
        62: {'TN': 3.0, 'TP': 0.15, 'NO3-N': 1.5, 'NO2-N': 0.005, 'NH4-N': 0.15, 'PO4-P': 0.01},   # Tree cover, broadleaved, deciduous, open
        70: {'TN': 2.5, 'TP': 0.12, 'NO3-N': 1.2, 'NO2-N': 0.005, 'NH4-N': 0.12, 'PO4-P': 0.008},  # Tree cover, needleleaved, evergreen
        71: {'TN': 2.5, 'TP': 0.12, 'NO3-N': 1.2, 'NO2-N': 0.005, 'NH4-N': 0.12, 'PO4-P': 0.008},  # Tree cover, needleleaved, evergreen, closed
        72: {'TN': 2.5, 'TP': 0.12, 'NO3-N': 1.2, 'NO2-N': 0.005, 'NH4-N': 0.12, 'PO4-P': 0.008},  # Tree cover, needleleaved, evergreen, open
        80: {'TN': 2.5, 'TP': 0.12, 'NO3-N': 1.2, 'NO2-N': 0.005, 'NH4-N': 0.12, 'PO4-P': 0.008},  # Tree cover, needleleaved, deciduous
        81: {'TN': 2.5, 'TP': 0.12, 'NO3-N': 1.2, 'NO2-N': 0.005, 'NH4-N': 0.12, 'PO4-P': 0.008},  # Tree cover, needleleaved, deciduous, closed
        82: {'TN': 2.5, 'TP': 0.12, 'NO3-N': 1.2, 'NO2-N': 0.005, 'NH4-N': 0.12, 'PO4-P': 0.008},  # Tree cover, needleleaved, deciduous, open
        90: {'TN': 3.0, 'TP': 0.15, 'NO3-N': 1.5, 'NO2-N': 0.005, 'NH4-N': 0.15, 'PO4-P': 0.01},   # Tree cover, mixed leaf type
        100: {'TN': 4.0, 'TP': 0.2, 'NO3-N': 2.0, 'NO2-N': 0.008, 'NH4-N': 0.2, 'PO4-P': 0.015},   # Mosaic tree and shrub (>50%) / herbaceous

        # Grassland and herbaceous classes (Grassland / pasture)
        110: {'TN': 5.0, 'TP': 0.2, 'NO3-N': 2.5, 'NO2-N': 0.01, 'NH4-N': 0.2, 'PO4-P': 0.025},   # Mosaic herbaceous cover (>50%) / tree and shrub
        120: {'TN': 4.0, 'TP': 0.15, 'NO3-N': 2.0, 'NO2-N': 0.01, 'NH4-N': 0.15, 'PO4-P': 0.02},  # Shrubland
        121: {'TN': 4.0, 'TP': 0.15, 'NO3-N': 2.0, 'NO2-N': 0.01, 'NH4-N': 0.15, 'PO4-P': 0.02},  # Shrubland evergreen
        122: {'TN': 4.0, 'TP': 0.15, 'NO3-N': 2.0, 'NO2-N': 0.01, 'NH4-N': 0.15, 'PO4-P': 0.02},  # Shrubland deciduous
        130: {'TN': 6.0, 'TP': 0.2, 'NO3-N': 3.0, 'NO2-N': 0.01, 'NH4-N': 0.25, 'PO4-P': 0.03},   # Grassland
        140: {'TN': 2.0, 'TP': 0.1, 'NO3-N': 1.0, 'NO2-N': 0.005, 'NH4-N': 0.1, 'PO4-P': 0.01},   # Lichens and mosses

        # Sparse vegetation and bare (Forest / natural catchment - low end)
        150: {'TN': 1.0, 'TP': 0.05, 'NO3-N': 0.5, 'NO2-N': 0.005, 'NH4-N': 0.05, 'PO4-P': 0.005},  # Sparse vegetation
        151: {'TN': 1.0, 'TP': 0.05, 'NO3-N': 0.5, 'NO2-N': 0.005, 'NH4-N': 0.05, 'PO4-P': 0.005},  # Sparse tree
        152: {'TN': 1.0, 'TP': 0.05, 'NO3-N': 0.5, 'NO2-N': 0.005, 'NH4-N': 0.05, 'PO4-P': 0.005},  # Sparse shrub
        153: {'TN': 1.0, 'TP': 0.05, 'NO3-N': 0.5, 'NO2-N': 0.005, 'NH4-N': 0.05, 'PO4-P': 0.005},  # Sparse herbaceous

        # Flooded/wetland classes (Mixed rural - moderate values)
        160: {'TN': 5.0, 'TP': 0.3, 'NO3-N': 2.5, 'NO2-N': 0.02, 'NH4-N': 0.3, 'PO4-P': 0.05},   # Tree cover, flooded, fresh or brackish water
        170: {'TN': 5.0, 'TP': 0.3, 'NO3-N': 2.5, 'NO2-N': 0.02, 'NH4-N': 0.3, 'PO4-P': 0.05},   # Tree cover, flooded, saline water
        180: {'TN': 5.0, 'TP': 0.3, 'NO3-N': 2.5, 'NO2-N': 0.02, 'NH4-N': 0.3, 'PO4-P': 0.05},   # Shrub or herbaceous cover, flooded

        # Urban and bare (Mixed rural to urban - moderate to high)
        190: {'TN': 8.0, 'TP': 1.0, 'NO3-N': 5.0, 'NO2-N': 0.03, 'NH4-N': 0.4, 'PO4-P': 0.15},   # Urban areas
        200: {'TN': 0.5, 'TP': 0.05, 'NO3-N': 0.25, 'NO2-N': 0.002, 'NH4-N': 0.025, 'PO4-P': 0.005},  # Bare areas
        201: {'TN': 0.5, 'TP': 0.05, 'NO3-N': 0.25, 'NO2-N': 0.002, 'NH4-N': 0.025, 'PO4-P': 0.005},  # Consolidated bare areas
        202: {'TN': 0.5, 'TP': 0.05, 'NO3-N': 0.25, 'NO2-N': 0.002, 'NH4-N': 0.025, 'PO4-P': 0.005},  # Unconsolidated bare areas

        # Water and snow/ice
        210: {'TN': 0.0, 'TP': 0.0, 'NO3-N': 0.0, 'NO2-N': 0.0, 'NH4-N': 0.0, 'PO4-P': 0.0},   # Water bodies
        220: {'TN': 0.0, 'TP': 0.0, 'NO3-N': 0.0, 'NO2-N': 0.0, 'NH4-N': 0.0, 'PO4-P': 0.0},   # Permanent snow and ice
    }


def set_ss_from_copernicus_lulc_with_loads(
        ss_config_filepath: str,
        json_header_comment: List[str],
        ss_method_copernicus_basins_hrus: Dict[str, str],
        ss_method_copernicus_nc_lc_dir: Optional[str],
        ss_method_copernicus_period: List[Union[int, float]],
        ss_method_copernicus_default_loads_bool: bool,

        ss_method_copernicus_compartment_name_for_load: str,
        ss_method_copernicus_annual_to_seasonal_loads_method: str = "uniform",
        ss_metadata_comment: str = "Nutrient loading from land use/land cover",
        ss_metadata_source: str = "Copernicus LULC analysis",
        optional_load_coefficients: Optional[Dict[int, Dict[str, float]]] = None,
        recursive: bool = False,
        file_pattern: str = 'ESACCI-LC-*.nc',
        simulation_period: Optional[List[int]] = None,
        bgc_engine_label: str = "",
        chemical_species_list: Optional[List[str]] = None,
        climate_data: Optional[Dict] = None,
        precip_scaling_power: float = 1.0,
        temp_q10: float = 2.0,
        temp_reference_c: float = 15.0
):
    """
    Process Copernicus LULC data, calculate nutrient loads, and generate OpenWQ JSON.

    This is a comprehensive wrapper that:
    1. Calculates LULC areas
    2. Calculates nutrient loads based on export coefficients
    3. Generates OpenWQ source/sink JSON configuration
    4. Exports all results

    Parameters:
    -----------
    ss_config_filepath : str
        Path where JSON configuration will be saved
    ss_method_copernicus_basins_hrus : dict
        Dictionary containing:
        - 'path_to_shp': Path to the shapefile
        - 'mapping_key': Column name for HRU ID (must match HDF5 variable)
    ss_method_copernicus_nc_lc_dir : str
        Directory containing Copernicus LULC NetCDF files
    ss_method_copernicus_period : list
        [year_start, year_end]
    ss_method_copernicus_default_loads_bool : bool
        If True, use default pre-coded load coefficients.
        If False, use optional_load_coefficients.
    json_header_comment : List[str]
        Comment lines to add at the top of the JSON file
    ss_method_copernicus_compartment_name_for_load : str
        OpenWQ compartment name (e.g., 'RIVER_NETWORK_REACHES')
    ss_metadata_comment : str
        Comment for JSON metadata section
    ss_metadata_source : str
        Source identifier for JSON metadata section
    optional_load_coefficients : dict, optional
        Custom load coefficients (only used if ss_method_copernicus_default_loads_bool=False)
    recursive : bool
        Whether to search subdirectories
    file_pattern : str
        Glob pattern for NetCDF files
    s_method_copernicus_annual_to_seasonal_loads_method : str, default "uniform"
        Temporal distribution method for annual loads:
        - 'uniform': Distribute loads evenly across 12 months (equal 1/12)
        - 'seasonal': Climate-weighted distribution using precipitation + temperature
          data from climate_data.  Falls back to a sine-curve if climate_data is
          not provided for a given year.
        In both cases the sum of 12 monthly loads equals the annual load.
    climate_data : dict, optional
        Dict keyed by year: {year: {'precip_mm': [12 floats], 'temp_c': [12 floats]}}.
        Used when method='seasonal' to weight monthly loads.
    precip_scaling_power : float
        Exponent for precipitation scaling (default 1.0)
    temp_q10 : float
        Q10 temperature sensitivity (default 2.0)
    temp_reference_c : float
        Reference temperature [°C] (default 15.0)

    Returns:
    --------
    tuple : (results_df, summaries, clipped_rasters, loads_df)
    """
    # Calculate LULC areas
    results, summaries, rasters = set_ss_from_copernicus_lulc(
        ss_config_filepath=ss_config_filepath,
        ss_method_copernicus_basins_hrus=ss_method_copernicus_basins_hrus,
        ss_method_copernicus_nc_lc_dir=ss_method_copernicus_nc_lc_dir,
        ss_method_copernicus_period=ss_method_copernicus_period,
        recursive=recursive,
        file_pattern=file_pattern
    )

    # Determine which coefficients to use
    if ss_method_copernicus_default_loads_bool:
        print("\nUsing DEFAULT pre-coded Copernicus land cover nutrient coefficients...")
        load_coefficients = get_default_copernicus_load_coefficients()
    else:
        if optional_load_coefficients is None:
            raise ValueError(
                "ss_method_copernicus_default_loads_bool=False but no optional_load_coefficients provided. "
                "Please provide custom load coefficients or set ss_method_copernicus_default_loads_bool=True."
            )
        print("\nUsing CUSTOM load coefficients from optional_load_coefficients...")
        print("Classes not specified in optional_load_coefficients will have loads set to 0.")
        load_coefficients = optional_load_coefficients

    # Get output directory
    output_dir = Path(ss_config_filepath).parent / 'ss_copernicus_files'

    # Get HRU ID column name
    shp_hru_id_column = ss_method_copernicus_basins_hrus.get('mapping_key', 'HRU_ID')

    # Calculate nutrient loads
    loads_df, pivot_loads_df = calculate_nutrient_loads(
        results_df=results,
        load_coefficients=load_coefficients,
        output_dir=output_dir,
        shp_hru_id_column=shp_hru_id_column
    )

    # Generate OpenWQ JSON configuration
    pivot_csv_path = output_dir / 'nutrient_loads_pivot.csv'

    create_openwq_ss_json_from_loads(
        pivot_loads_csv_path=pivot_csv_path,
        ss_config_filepath=ss_config_filepath,
        json_header_comment=json_header_comment,
        ss_method_copernicus_compartment_name_for_load=ss_method_copernicus_compartment_name_for_load,
        shp_hru_id_column=shp_hru_id_column,
        ss_method_copernicus_annual_to_seasonal_loads_method=ss_method_copernicus_annual_to_seasonal_loads_method,
        ss_metadata_comment=ss_metadata_comment,
        ss_metadata_source=ss_metadata_source,
        simulation_period=simulation_period,
        bgc_engine_label=bgc_engine_label,
        chemical_species_list=chemical_species_list,
        climate_data=climate_data,
        precip_scaling_power=precip_scaling_power,
        temp_q10=temp_q10,
        temp_reference_c=temp_reference_c,
    )

    return results, summaries, rasters, loads_df