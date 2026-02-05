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
        ss_method_csv_config: List[Dict[str, Union[str, int]]]

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

    # Build the complete JSON structure
    config = {
        "METADATA": {
            "Comment": ss_metadata_comment,
            "Source": ss_metadata_source
        }
    }

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
        print(f"Loading shapefile: {shapefile_info['path_to_shp']}")
        self.gdf = gpd.read_file(shapefile_info['path_to_shp'])
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
            lats = ds.lat.values
            lons = ds.lon.values
        elif 'latitude' in ds.dims and 'longitude' in ds.dims:
            lats = ds.latitude.values
            lons = ds.longitude.values
        else:
            print(" ✗")
            raise ValueError(f"Could not find lat/lon. Available dims: {list(ds.dims.keys())}")
        print(" ✓")

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
            gdf_reprojected = self.gdf.to_crs(raster_crs)
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
        if raster_src.crs.is_geographic:
            # For geographic CRS (lat/lon), use approximate area
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

            if convert_to_tiff:
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

            if convert_to_tiff:
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
            gdf_reprojected = self.gdf.to_crs(src.crs)
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
                out_image, out_transform = rasterio_mask(
                    src,
                    [basin_boundary],
                    crop=True,
                    nodata=src.nodata if src.nodata is not None else -9999
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


def calc_copernicus_lulc(
        ss_config_filepath: str,
        ss_method_copernicus_basin_info: Dict[str, str],
        ss_method_copernicus_nc_lc_dir: str,
        ss_method_copernicus_period: List[Union[int, float]],
        recursive: bool = False,
        file_pattern: str = 'ESACCI-LC-*.nc'
):
    """
    Process all Copernicus LULC NetCDF files from a directory.

    Parameters:
    -----------
    ss_config_filepath : str
        Path where outputs will be saved (directory will be extracted from this)
    ss_method_copernicus_basin_info : dict
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

    print(f"Base output directory: {output_base_dir}")
    print(f"SS Copernicus files directory: {ss_output_dir}")
    print(f"Clipped rasters will be saved to: {clipped_rasters_dir}")

    # Convert directory path to Path object
    lc_dir = Path(ss_method_copernicus_nc_lc_dir)

    # Check if directory exists
    if not lc_dir.exists():
        raise FileNotFoundError(f"Directory not found: {ss_method_copernicus_nc_lc_dir}")

    if not lc_dir.is_dir():
        raise NotADirectoryError(f"Path is not a directory: {ss_method_copernicus_nc_lc_dir}")

    # Find all NetCDF files
    print(f"Searching for NetCDF files in: {lc_dir}")
    print(f"Pattern: {file_pattern}")
    print(f"Recursive: {recursive}")

    if recursive:
        # Search recursively
        netcdf_files = sorted(lc_dir.rglob(file_pattern))
    else:
        # Search only in the specified directory
        netcdf_files = sorted(lc_dir.glob(file_pattern))

    if not netcdf_files:
        # Try alternative patterns
        alternative_patterns = ['*.nc', 'ESACCI-*.nc', 'LC-*.nc']
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
    print(f"Shapefile: {ss_method_copernicus_basin_info['path_to_shp']}")
    print(f"HRU ID column (mapping_key): {ss_method_copernicus_basin_info['mapping_key']}")
    print(f"Output directory: {ss_output_dir}")

    analyzer = OptimizedLULCAnalyzer(
        ss_method_copernicus_basin_info,
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
 
_These search results were found and analyzed using Consensus, an AI-powered search engine for research. Try it at https://consensus.app. © 2026 Consensus NLP, Inc. Personal, non-commercial use only; redistribution requires copyright holders’ consent._
 
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
        ss_method_copernicus_basin_info: Dict[str, str],
        ss_method_copernicus_nc_lc_dir: str,
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
    ss_method_copernicus_basin_info : dict
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
        ss_method_copernicus_basin_info=ss_method_copernicus_basin_info,
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
        method: str = "uniform"
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
        - "uniform": Spread load uniformly throughout the year (monthly)
        - "seasonal": Sine curve with peak in growing season (summer), low in winter

    Returns:
    --------
    list of tuples : [(month, day, hour, load_kg), ...]
        List of temporal load entries
    """
    import math

    if method == "uniform":
        # Distribute load uniformly across 12 months (1st of each month at 1:00 AM)
        monthly_load = annual_load_kg / 12.0

        temporal_entries = []
        for month in range(1, 13):
            # Using day=1, hour=1, allows us to use "all" for other time fields
            temporal_entries.append((month, 1, 1, monthly_load))

        return temporal_entries

    elif method == "seasonal":
        # Sine curve: peak in July (growing season), minimum in January (winter)
        monthly_loads = []
        for month in range(1, 13):
            # Sine function centered on July (month 7)
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
            f"Valid options are 'uniform', 'seasonal', or 'climate_adjusted'."
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


def set_ss_climate_adjusted_export_coefficients(
        ss_config_filepath: str,
        json_header_comment: List[str],
        ss_method_copernicus_basin_info: Dict[str, str],
        ss_method_copernicus_nc_lc_dir: str,
        ss_method_copernicus_period: List[Union[int, float]],
        ss_method_copernicus_default_loads_bool: bool,
        ss_method_copernicus_compartment_name_for_load: str,
        ss_method_copernicus_openwq_h5_results_file_example_for_mapping_key: Dict[str, str],
        climate_data: Dict[int, Dict[str, List[float]]],
        precip_scaling_power: float = 1.0,
        temp_q10: float = 2.0,
        temp_reference_c: float = 15.0,
        ss_metadata_comment: str = "Climate-adjusted nutrient loading from LULC",
        ss_metadata_source: str = "Copernicus LULC + climate adjustment",
        optional_load_coefficients: Optional[Dict[int, Dict[str, float]]] = None,
        recursive: bool = False,
        file_pattern: str = 'ESACCI-LC-*.nc'
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
    ss_method_copernicus_basin_info : dict
        Shapefile info: {'path_to_shp': ..., 'mapping_key': ...}
    ss_method_copernicus_nc_lc_dir : str
        Directory with Copernicus LULC NetCDF files
    ss_method_copernicus_period : list
        [year_start, year_end]
    ss_method_copernicus_default_loads_bool : bool
        If True, use default export coefficients; if False, use optional_load_coefficients
    ss_method_copernicus_compartment_name_for_load : str
        OpenWQ compartment name (e.g., 'RIVER_NETWORK_REACHES')
    ss_method_copernicus_openwq_h5_results_file_example_for_mapping_key : dict
        HDF5 file info: {'path_to_file': ..., 'mapping_key': ...}
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
    import h5py

    print("\n" + "=" * 60)
    print("CLIMATE-ADJUSTED EXPORT COEFFICIENT METHOD")
    print("=" * 60)
    print(f"Precipitation scaling power (alpha): {precip_scaling_power}")
    print(f"Temperature Q10: {temp_q10}")
    print(f"Reference temperature: {temp_reference_c} deg C")
    print(f"Climate data years available: {sorted(climate_data.keys())}")

    # Step 1: Calculate LULC areas (reuses existing Copernicus pipeline)
    results, summaries, rasters = set_ss_from_copernicus_lulc(
        ss_config_filepath=ss_config_filepath,
        ss_method_copernicus_basin_info=ss_method_copernicus_basin_info,
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
    shp_hru_id_column = ss_method_copernicus_basin_info.get('mapping_key', 'HRU_ID')

    loads_df, pivot_loads_df = calculate_nutrient_loads(
        results_df=results,
        load_coefficients=load_coefficients,
        output_dir=output_dir,
        shp_hru_id_column=shp_hru_id_column
    )

    # Step 4: Load HDF5 spatial mapping
    openwq_hru_mapping_key = ss_method_copernicus_openwq_h5_results_file_example_for_mapping_key['mapping_key']
    openwq_h5_result_path = ss_method_copernicus_openwq_h5_results_file_example_for_mapping_key['path_to_file']

    print(f"\nLoading spatial mapping from: {openwq_h5_result_path}")
    with h5py.File(openwq_h5_result_path, 'r') as h5f:
        xyz_elements = h5f['xyz_elements'][:]
        mapping_ids = h5f[openwq_hru_mapping_key][:]

    hru_to_xyz = {}
    for idx in range(len(mapping_ids)):
        hru_id = mapping_ids[idx].decode('utf-8')
        ix, iy, iz = xyz_elements[:, idx]
        hru_to_xyz[hru_id] = (int(ix), int(iy), int(iz))

    print(f"  Mapped {len(hru_to_xyz)} HRU/segments")

    # Step 5: Build climate-adjusted JSON
    pivot_csv_path = output_dir / 'nutrient_loads_pivot.csv'
    pivot_df = pd.read_csv(pivot_csv_path)

    exclude_cols = [shp_hru_id_column, 'Year']
    nutrient_cols = [col for col in pivot_df.columns if col not in exclude_cols]

    config = {
        "METADATA": {
            "Comment": ss_metadata_comment,
            "Source": ss_metadata_source,
            "Climate_Adjustment": {
                "precip_scaling_power": precip_scaling_power,
                "temp_Q10": temp_q10,
                "temp_reference_C": temp_reference_c
            }
        }
    }

    # Helper: find nearest climate data year
    climate_years = sorted(climate_data.keys())

    def _get_climate_for_year(yr):
        if yr in climate_data:
            return climate_data[yr]
        # Find nearest year
        nearest = min(climate_years, key=lambda y: abs(y - yr))
        print(f"    (using climate data from {nearest} for year {yr})")
        return climate_data[nearest]

    entry_idx = 1

    for year in sorted(pivot_df['Year'].unique()):
        year_data = pivot_df[pivot_df['Year'] == year]
        climate = _get_climate_for_year(int(year))

        for nutrient in nutrient_cols:
            print(f"\n  Processing: Year {year}, Nutrient {nutrient}")

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

                if hru_id not in hru_to_xyz:
                    continue

                ix, iy, iz = hru_to_xyz[hru_id]

                # Climate-adjusted temporal distribution
                temporal_entries = _apply_climate_adjustment(
                    annual_load_kg=annual_load_kg,
                    year=int(year),
                    monthly_precip_mm=climate['precip_mm'],
                    monthly_temp_c=climate['temp_c'],
                    precip_scaling_power=precip_scaling_power,
                    temp_q10=temp_q10,
                    temp_reference_c=temp_reference_c
                )

                for month, day, hour, load_kg in temporal_entries:
                    data_entries[str(sub_idx)] = [
                        int(year), int(month), int(day), int(hour), "all", "all",
                        ix, iy, iz,
                        float(load_kg),
                        "discrete"
                    ]
                    sub_idx += 1

            if len(data_entries) == 0:
                continue

            config[str(entry_idx)] = {
                "CHEMICAL_NAME": nutrient,
                "COMPARTMENT_NAME": ss_method_copernicus_compartment_name_for_load,
                "COMMENT": f"Climate-adjusted loading for {nutrient}",
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
        ss_method_copernicus_openwq_h5_results_file_example_for_mapping_key: Dict[str, str],
        shp_hru_id_column: str,
        ss_method_copernicus_annual_to_seasonal_loads_method: str = "uniform",
        ss_metadata_comment: str = "Nutrient loading from land use/land cover",
        ss_metadata_source: str = "Copernicus LULC analysis"
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
    openwq_hru_mapping_key : str
        Name of the HRU/segment ID column (e.g., 'Seg_ID', 'HRU_ID')
    ss_method_copernicus_openwq_h5_results_file_example_for_mapping_key : Dist[str, str]
        Path to OpenWQ HDF5 output file containing xyz_elements and mapping_key data
    ss_metadata_comment : str
        Comment for metadata section
    ss_metadata_source : str
        Source identifier for metadata section
    """
    import h5py

    openwq_hru_mapping_key = ss_method_copernicus_openwq_h5_results_file_example_for_mapping_key['mapping_key']
    openwq_h5_result_path = ss_method_copernicus_openwq_h5_results_file_example_for_mapping_key['path_to_file']

    print("\n" + "=" * 60)
    print("GENERATING OPENWQ SOURCE/SINK JSON")
    print("=" * 60)

    # Validate temporal distribution method
    valid_methods = ["uniform", "seasonal"]
    if ss_method_copernicus_annual_to_seasonal_loads_method not in valid_methods:
        raise ValueError(
            f"Invalid temporal distribution method: '{ss_method_copernicus_annual_to_seasonal_loads_method}'. "
            f"Valid options are: {valid_methods}"
        )

    print(f"\nTemporal distribution method: {ss_method_copernicus_annual_to_seasonal_loads_method}")
    if ss_method_copernicus_annual_to_seasonal_loads_method == "uniform":
        print("  → Annual loads will be distributed uniformly across 12 months")
        print("  → Using OpenWQ 'all' keyword for compact JSON (min=0, sec=0)")
    elif ss_method_copernicus_annual_to_seasonal_loads_method == "seasonal":
        print("  → Annual loads will follow seasonal pattern (sine curve)")
        print("  → Peak loads in summer (July), lowest in winter (January)")
        print("  → Using OpenWQ 'all' keyword for compact JSON (min=0, sec=0)")

    # Load pivot loads data
    print(f"\nLoading pivot loads from: {pivot_loads_csv_path.name}")
    pivot_df = pd.read_csv(pivot_loads_csv_path)

    # Load HDF5 mapping
    print(f"Loading spatial mapping from: {openwq_h5_result_path}")
    with h5py.File(openwq_h5_result_path, 'r') as h5f:
        # Load xyz_elements
        if 'xyz_elements' not in h5f:
            raise ValueError(f"'xyz_elements' not found in HDF5 file. Available: {list(h5f.keys())}")

        xyz_elements = h5f['xyz_elements'][:]

        # Load mapping key (e.g., Seg_ID)
        if openwq_hru_mapping_key not in h5f:
            raise ValueError(
                f"'{openwq_hru_mapping_key}' not found in HDF5 file. Available: {list(h5f.keys())}\n"
                f"Make sure the mapping_key matches a variable in the HDF5 file."
            )

        mapping_ids = h5f[openwq_hru_mapping_key][:]

    print(f"  Loaded {len(xyz_elements)} spatial elements")
    print(f"  xyz_elements shape: {xyz_elements.shape}")
    print(f"  {openwq_hru_mapping_key} shape: {mapping_ids.shape}")

    # Create mapping dictionary: HRU_ID -> (ix, iy, iz)
    hru_to_xyz = {}
    for idx in range(len(mapping_ids)):
        hru_id = mapping_ids[idx].decode('utf-8')
        ix, iy, iz = xyz_elements[:, idx]
        hru_to_xyz[hru_id] = (int(ix), int(iy), int(iz))

    print(f"  Created mapping for {len(hru_to_xyz)} HRU/segments")

    # Get nutrient columns (exclude HRU_ID and Year columns)
    exclude_cols = [shp_hru_id_column, 'Year']
    nutrient_cols = [col for col in pivot_df.columns if col not in exclude_cols]

    print(f"\nNutrient species found: {nutrient_cols}")
    print(f"Years found: {sorted(pivot_df['Year'].unique())}")

    # Build JSON structure
    config = {
        "METADATA": {
            "Comment": ss_metadata_comment,
            "Source": ss_metadata_source
        }
    }

    entry_idx = 1

    # Group by year and process each nutrient
    for year in sorted(pivot_df['Year'].unique()):
        year_data = pivot_df[pivot_df['Year'] == year]

        for nutrient in nutrient_cols:
            print(f"\n  Processing: Year {year}, Nutrient {nutrient}")

            # Get non-zero loads for this year and nutrient
            nutrient_loads = year_data[[shp_hru_id_column, nutrient]].copy()
            nutrient_loads = nutrient_loads[nutrient_loads[nutrient] > 0]  # Only non-zero loads

            if len(nutrient_loads) == 0:
                print(f"    ⚠ No non-zero loads for {nutrient} in {year}, skipping...")
                continue

            # Create data entries
            data_entries = {}
            sub_idx = 1

            for _, row in nutrient_loads.iterrows():
                hru_id = row[shp_hru_id_column]
                hru_id = str(int(hru_id))
                annual_load_kg = row[nutrient]

                # Get spatial coordinates
                if hru_id not in hru_to_xyz.keys():
                    print(f"    ⚠ Warning: HRU {hru_id} not found in HDF5 mapping, skipping...")
                    continue

                ix, iy, iz = hru_to_xyz[hru_id]

                # Calculate temporal distribution of annual load
                temporal_entries = _calculate_temporal_load_distribution(
                    annual_load_kg=annual_load_kg,
                    year=int(year),
                    method=ss_method_copernicus_annual_to_seasonal_loads_method
                )

                # Create data entries for each time point
                # Format: [YYYY, MM, DD, HH, min, sec, ix, iy, iz, load, "discrete"]
                # Use "all" for minutes and seconds (always 0) to compress JSON
                for month, day, hour, load_kg in temporal_entries:
                    data_entries[str(sub_idx)] = [
                        int(year), int(month), int(day), int(hour), "all", "all",  # Temporal (use "all" for min, sec)
                        ix, iy, iz,  # Spatial coordinates
                        float(load_kg),  # Load in kg
                        "discrete"
                    ]
                    sub_idx += 1

            if len(data_entries) == 0:
                print(f"    ⚠ No valid entries for {nutrient} in {year}, skipping...")
                continue

            # Add to config
            config[str(entry_idx)] = {
                "CHEMICAL_NAME": nutrient,
                "COMPARTMENT_NAME": ss_method_copernicus_compartment_name_for_load,
                "COMMENT": f"Copernicus LULC-based loading for {nutrient} in {year}",
                "TYPE": "source",
                "UNITS": "kg",
                "DATA_FORMAT": "JSON",
                "DATA": data_entries
            }

            print(f"    Added {len(data_entries)} load entries for {nutrient} in {year}")
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
        ss_method_copernicus_basin_info: Dict[str, str],
        ss_method_copernicus_nc_lc_dir: str,
        ss_method_copernicus_period: List[Union[int, float]],
        ss_method_copernicus_default_loads_bool: bool,

        ss_method_copernicus_compartment_name_for_load: str,
        ss_method_copernicus_openwq_h5_results_file_example_for_mapping_key: Dict[str, str],
        ss_method_copernicus_annual_to_seasonal_loads_method: str = "uniform",
        ss_metadata_comment: str = "Nutrient loading from land use/land cover",
        ss_metadata_source: str = "Copernicus LULC analysis",
        optional_load_coefficients: Optional[Dict[int, Dict[str, float]]] = None,
        recursive: bool = False,
        file_pattern: str = 'ESACCI-LC-*.nc'
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
    ss_method_copernicus_basin_info : dict
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
    ss_method_copernicus_openwq_h5_results_file_example_for_mapping_key : List[str, str]
        Path to OpenWQ HDF5 output file with xyz_elements and mapping_key
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
        - 'uniform': Distribute loads evenly across 12 months
        - 'seasonal': Use sine curve with peak in summer (July), minimum in winter (January)
        Uses OpenWQ 'all' keyword for minutes/seconds to compress JSON file size

    Returns:
    --------
    tuple : (results_df, summaries, clipped_rasters, loads_df)
    """
    # Calculate LULC areas
    results, summaries, rasters = set_ss_from_copernicus_lulc(
        ss_config_filepath=ss_config_filepath,
        ss_method_copernicus_basin_info=ss_method_copernicus_basin_info,
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
    shp_hru_id_column = ss_method_copernicus_basin_info.get('mapping_key', 'HRU_ID')

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
        ss_method_copernicus_openwq_h5_results_file_example_for_mapping_key=ss_method_copernicus_openwq_h5_results_file_example_for_mapping_key,
        shp_hru_id_column=shp_hru_id_column,
        ss_method_copernicus_annual_to_seasonal_loads_method=ss_method_copernicus_annual_to_seasonal_loads_method,
        ss_metadata_comment=ss_metadata_comment,
        ss_metadata_source=ss_metadata_source
    )

    return results, summaries, rasters, loads_df