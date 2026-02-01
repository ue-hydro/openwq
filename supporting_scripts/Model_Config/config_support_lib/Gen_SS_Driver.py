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
        ss_method_csv_metadata_source: str,
        ss_method_csv_metadata_comment: str,

        # List of source/sink configurations
        ss_method_csv_config: List[Dict[str, Union[str, int]]]

) -> None:
    """
    Create OpenWQ Source/Sink configuration JSON file.

    Parameters:
        ss_config_filepath: Full path where the JSON file will be saved
        json_header_comment: List of comment lines to add at the top of the file
        metadata_comment: Comment for metadata section
        metadata_source: Source identifier for metadata section
        ss_method_csv_config: List of dictionaries, each containing:
            - Chemical_name: Name of the chemical species
            - Compartment_name: Name of the compartment
            - Type: "source" or "sink"
            - Units: Units (e.g., "kg")
            - Data_Format: Format type (e.g., "ASCII")
            - Filepath: Path to data file
            - Delimiter: CSV delimiter (e.g., ",")
            - Number_of_header_rows: Number of header rows
            - Header_key_row: Row number containing headers
    """

    # Create the destination directory path if it doesn't exist
    output_path = Path(ss_config_filepath)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Build the complete JSON structure
    config = {
        "METADATA": {
            "Comment": ss_method_csv_metadata_comment,
            "Source": ss_method_csv_metadata_source
        }
    }

    # Add each source/sink configuration with numbered keys
    for idx, ss_config in enumerate(ss_method_csv_config, start=1):
        config[str(idx)] = {
            "Chemical_name": ss_config["Chemical_name"],
            "Compartment_name": ss_config["Compartment_name"],
            "Type": ss_config["Type"],
            "Units": ss_config["Units"],
            "Data_Format": "ASCII",
            "Data": {
                "Filepath": ss_config["Filepath"],
                "Delimiter": ss_config["Delimiter"],
                "Number_of_header_rows": ss_config["Number_of_header_rows"],
                "Header_key_row": ss_config["Header_key_row"]
            }
        }

    # Convert to JSON string with standard formatting
    json_string = json.dumps(config, indent=4)

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

    def __init__(self, shapefile_path, output_dir='output'):
        """
        Initialize the LULC analyzer.

        Parameters:
        -----------
        shapefile_path : str
            Path to the shapefile containing HRU polygons
        output_dir : str
            Directory to save output files
        """
        self.shapefile_path = shapefile_path
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)

        # Load shapefile
        print(f"Loading shapefile: {shapefile_path}")
        self.gdf = gpd.read_file(shapefile_path)
        print(f"Loaded {len(self.gdf)} HRU polygons")
        print(f"CRS: {self.gdf.crs}")

        # Ensure there's an ID column
        if 'HRU_ID' not in self.gdf.columns and 'ID' not in self.gdf.columns:
            self.gdf['HRU_ID'] = range(len(self.gdf))

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
                # Get HRU ID
                hru_id = row.get('HRU_ID', row.get('ID', idx))

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
                            'HRU_ID': hru_id,
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

        # Summary by HRU and Year
        summaries['by_hru_year'] = df.groupby(['HRU_ID', 'Year']).agg({
            'Area_km2': 'sum',
            'LC_Class': 'count'
        }).rename(columns={'LC_Class': 'Num_LC_Classes'})

        # Summary by LC Class and Year
        summaries['by_lc_year'] = df.groupby(['LC_Class', 'Year']).agg({
            'Area_km2': 'sum',
            'HRU_ID': 'count'
        }).rename(columns={'HRU_ID': 'Num_HRUs'})

        # Pivot table: HRUs x LC Classes
        summaries['pivot_hru_lc'] = df.pivot_table(
            values='Area_km2',
            index=['HRU_ID', 'Year'],
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

        print(f"\nTotal records: {len(df)}")
        print(f"Number of HRUs: {df['HRU_ID'].nunique()}")
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


def set_ss_from_copericus_lulc(
        basin_shapefile_path: str,
        netcdf_copernicus_lc_dir: str,
        ss_config_filepath: str,
        year_start: Optional[int] = None,
        year_end: Optional[int] = None,
        recursive: bool = False,
        file_pattern: str = 'ESACCI-LC-*.nc'
):
    """
    Process all Copernicus LULC NetCDF files from a directory.

    Parameters:
    -----------
    basin_shapefile_path : str
        Path to the shapefile containing HRU/basin polygons
    netcdf_copernicus_lc_dir : str
        Directory containing Copernicus LULC NetCDF files
    ss_config_filepath : str
        Path where outputs will be saved (directory will be extracted from this)
    year_start : int, optional
        Filter files to only include years >= this value
    year_end : int, optional
        Filter files to only include years <= this value
    recursive : bool, default False
        Whether to search subdirectories recursively
    file_pattern : str, default 'ESACCI-LC-*.nc'
        Glob pattern to match NetCDF files

    Returns:
    --------
    tuple : (results_df, summaries, clipped_rasters) - processed results, summary statistics, and paths to clipped rasters

    Examples:
    ---------
    # Process all files in directory
    results, summaries, rasters = set_ss_from_copericus_lulc(
        'basin.shp',
        'data/copernicus/',
        'output/config.json'
    )

    # Process only years 2000-2020
    results, summaries, rasters = set_ss_from_copericus_lulc(
        'basin.shp',
        'data/copernicus/',
        'output/config.json',
        year_start=2000,
        year_end=2020
    )
    """

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
    lc_dir = Path(netcdf_copernicus_lc_dir)

    # Check if directory exists
    if not lc_dir.exists():
        raise FileNotFoundError(f"Directory not found: {netcdf_copernicus_lc_dir}")

    if not lc_dir.is_dir():
        raise NotADirectoryError(f"Path is not a directory: {netcdf_copernicus_lc_dir}")

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
            f"No NetCDF files found in {netcdf_copernicus_lc_dir}\n"
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
    print(f"Shapefile: {basin_shapefile_path}")
    print(f"Output directory: {ss_output_dir}")

    analyzer = OptimizedLULCAnalyzer(
        basin_shapefile_path,
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

    return results_df, summaries, clipped_raster_paths


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