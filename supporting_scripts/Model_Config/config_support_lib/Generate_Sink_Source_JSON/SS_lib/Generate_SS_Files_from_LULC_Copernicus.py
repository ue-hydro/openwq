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

################################
# Generating SS json file from Copernicus LULC dataset
################################

import rasterio
import rasterio.mask
import geopandas as gpd
import numpy as np
import pandas as pd
from pathlib import Path
from tqdm import tqdm

# ESA CCI Land Cover Classes
LULC_CLASSES = {
    0: "No data",
    10: "Cropland, rainfed",
    11: "Cropland, rainfed, herbaceous cover",
    12: "Cropland, rainfed, tree or shrub cover",
    20: "Cropland, irrigated or post-flooding",
    30: "Mosaic cropland (>50%) / natural vegetation (<50%)",
    40: "Mosaic natural vegetation (>50%) / cropland (<50%)",
    50: "Tree cover, broadleaved, evergreen, closed to open (>15%)",
    60: "Tree cover, broadleaved, deciduous, closed to open (>15%)",
    61: "Tree cover, broadleaved, deciduous, closed (>40%)",
    62: "Tree cover, broadleaved, deciduous, open (15-40%)",
    70: "Tree cover, needleleaved, evergreen, closed to open (>15%)",
    71: "Tree cover, needleleaved, evergreen, closed (>40%)",
    72: "Tree cover, needleleaved, evergreen, open (15-40%)",
    80: "Tree cover, needleleaved, deciduous, closed to open (>15%)",
    81: "Tree cover, needleleaved, deciduous, closed (>40%)",
    82: "Tree cover, needleleaved, deciduous, open (15-40%)",
    90: "Tree cover, mixed leaf type",
    100: "Mosaic tree and shrub (>50%) / herbaceous cover (<50%)",
    110: "Mosaic herbaceous cover (>50%) / tree and shrub (<50%)",
    120: "Shrubland",
    121: "Shrubland evergreen",
    122: "Shrubland deciduous",
    130: "Grassland",
    140: "Lichens and mosses",
    150: "Sparse vegetation (tree, shrub, herbaceous cover) (<15%)",
    152: "Sparse shrub (<15%)",
    153: "Sparse herbaceous cover (<15%)",
    160: "Tree cover, flooded, fresh or brackish water",
    170: "Tree cover, flooded, saline water",
    180: "Shrub or herbaceous cover, flooded, fresh/saline/brackish water",
    190: "Urban areas",
    200: "Bare areas",
    201: "Consolidated bare areas",
    202: "Unconsolidated bare areas",
    210: "Water bodies",
    220: "Permanent snow and ice"
}


def calculate_pixel_area_km2(lat, pixel_size_deg=0.002777777777778):
    """Calculate pixel area in km² at given latitude"""
    R = 6371.0  # Earth's radius in km
    lat_rad = np.radians(lat)
    height_km = pixel_size_deg * (np.pi * R / 180.0)
    width_km = pixel_size_deg * (np.pi * R / 180.0) * np.cos(lat_rad)
    return height_km * width_km


def calculate_areas_single_basin(raster_path, polygon_gdf, basin_id):
    """
    Calculate LULC areas for a single basin

    Parameters:
    -----------
    raster_path : str
        Path to the GeoTIFF file
    polygon_gdf : GeoDataFrame
        GeoDataFrame containing the basin polygon
    basin_id : str
        Basin identifier

    Returns:
    --------
    dict : Dictionary with basin_id and area for each LULC class
    """

    result = {'basin_id': basin_id}

    try:
        with rasterio.open(raster_path) as src:
            # Reproject polygon to match raster CRS if needed
            if polygon_gdf.crs != src.crs:
                polygon_gdf = polygon_gdf.to_crs(src.crs)

            # Clip raster by polygon
            clipped_data, clipped_transform = rasterio.mask.mask(
                src, polygon_gdf.geometry, crop=True, all_touched=True
            )
            clipped_data = clipped_data[0]  # Get first band

            # Get latitude for each row
            height, width = clipped_data.shape
            rows = np.arange(height)
            lats = clipped_transform[5] + (rows + 0.5) * clipped_transform[4]

            # Get pixel resolution
            pixel_size_deg = abs(clipped_transform[4])

            # Calculate area for each row (latitude)
            areas_per_row = np.array([calculate_pixel_area_km2(lat, pixel_size_deg) for lat in lats])

            # Get unique classes and calculate areas
            unique_classes, counts = np.unique(clipped_data[clipped_data != 0], return_counts=True)

            # Calculate area for each class
            for lulc_class, count in zip(unique_classes, counts):
                class_mask = (clipped_data == lulc_class)
                class_area = np.sum(class_mask * areas_per_row[:, np.newaxis])

                # Use class number as column name
                result[f'class_{int(lulc_class)}'] = float(class_area)

        return result

    except Exception as e:
        print(f"  Error processing basin {basin_id}: {e}")
        return {'basin_id': basin_id, 'error': str(e)}


def batch_calculate_lulc_areas(raster_path, polygons_folder, output_csv, id_column=None):
    """
    Calculate LULC areas for all polygons in a folder

    Parameters:
    -----------
    raster_path : str
        Path to the GeoTIFF file
    polygons_folder : str
        Path to folder containing polygon files (GPKG, shapefiles, etc.)
    output_csv : str
        Path to output CSV file
    id_column : str, optional
        Column name in polygon files to use as basin ID. If None, uses filename

    Returns:
    --------
    DataFrame : Results with basins as rows and LULC classes as columns
    """

    polygons_folder = Path(polygons_folder)

    # Find all polygon files
    polygon_files = []
    for ext in ['*.gpkg', '*.shp', '*.geojson', '*.json']:
        polygon_files.extend(list(polygons_folder.glob(ext)))

    if not polygon_files:
        raise ValueError(f"No polygon files found in {polygons_folder}")

    print(f"Found {len(polygon_files)} polygon files")
    print(f"Raster: {raster_path}")
    print(f"Output: {output_csv}\n")

    # Process each polygon file
    all_results = []

    for polygon_file in tqdm(polygon_files, desc="Processing basins"):
        try:
            # Read polygon - try pyogrio first (faster), fallback to fiona
            import warnings
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore')
                try:
                    # Try pyogrio engine (faster and more stable)
                    gdf = gpd.read_file(str(polygon_file), engine='pyogrio')
                except:
                    # Fallback to default fiona engine
                    gdf = gpd.read_file(str(polygon_file))

            # Determine basin ID
            if id_column and id_column in gdf.columns:
                # Use specified column
                for idx, row in gdf.iterrows():
                    basin_id = str(row[id_column])
                    basin_gdf = gdf.iloc[[idx]]
                    result = calculate_areas_single_basin(raster_path, basin_gdf, basin_id)
                    all_results.append(result)
            else:
                # Use filename as basin ID
                basin_id = polygon_file.stem
                result = calculate_areas_single_basin(raster_path, gdf, basin_id)
                all_results.append(result)

        except Exception as e:
            print(f"\nError reading {polygon_file.name}: {e}")
            # Add failed basin to results
            all_results.append({
                'basin_id': polygon_file.stem,
                'error': str(e)
            })
            continue

    # Create DataFrame
    df = pd.DataFrame(all_results)

    # Check if we have any results
    if len(df) == 0 or 'basin_id' not in df.columns:
        print("\nERROR: No basins were successfully processed!")
        print("This could be due to:")
        print("  - Incompatible geopandas/fiona versions")
        print("  - Corrupted polygon files")
        print("  - CRS mismatch")
        print("\nTry updating packages:")
        print("  pip install --upgrade geopandas fiona")
        return pd.DataFrame()

    # Move basin_id to first column
    cols = ['basin_id'] + [col for col in df.columns if col != 'basin_id' and col != 'error']
    if 'error' in df.columns:
        cols.append('error')
    df = df[cols]

    # Fill NaN with 0 (classes not present in basin)
    for col in df.columns:
        if col.startswith('class_'):
            df[col] = df[col].fillna(0)

    # Add class names as a separate mapping
    class_name_mapping = {f'class_{k}': v for k, v in LULC_CLASSES.items()}

    # Sort columns by class number
    class_cols = sorted([col for col in df.columns if col.startswith('class_')],
                        key=lambda x: int(x.split('_')[1]))
    other_cols = [col for col in df.columns if not col.startswith('class_')]
    df = df[other_cols + class_cols]

    # Save to CSV
    df.to_csv(output_csv, index=False)

    print(f"\n{'=' * 70}")
    print(f"Results saved to: {output_csv}")
    print(f"Total basins processed: {len(df)}")
    print(f"Total LULC classes found: {len(class_cols)}")
    print(f"{'=' * 70}")

    # Save class name mapping to separate file
    mapping_file = Path(output_csv).with_suffix('.class_mapping.csv')
    mapping_df = pd.DataFrame([
        {'class_column': k, 'class_number': int(k.split('_')[1]), 'class_name': v}
        for k, v in class_name_mapping.items()
        if k in df.columns
    ])
    mapping_df.to_csv(mapping_file, index=False)
    print(f"Class name mapping saved to: {mapping_file}")

    return df

def GenSSjson_fromCopernicus_driver(
    hru_shp=None,
    copernicus_lulc_dir=None):

    # Find all tif files that correspond to all years
    rasterYear_folder = Path(copernicus_lulc_dir)

    # Find all polygon files
    rasterYear_files = []
    for ext in ['*.tif']:
        rasterYear_files.extend(list(rasterYear_folder.glob(ext)))

    for rastYear_f in rasterYear_files:
        print("******************************")
        print(f"Processing: {rastYear_f}")
        print("******************************")

        df = batch_calculate_lulc_areas(
            raster_path=rastYear_f,
            polygons_folder=hru_shp,
            output_csv=Path(str(rastYear_f)[0:-4] + ".csv")
        )

    # Show summary
    print("\nFirst 5 rows:")
    print(df.head())

    # Show total area per class across all basins
    class_cols = [col for col in df.columns if col.startswith('class_')]
    print("\nTotal area per class (km²):")
    for col in class_cols:
        class_num = int(col.split('_')[1])
        class_name = LULC_CLASSES.get(class_num, 'Unknown')
        total_area = df[col].sum()
        if total_area > 0:
            print(f"  {class_name}: {total_area:,.2f}")