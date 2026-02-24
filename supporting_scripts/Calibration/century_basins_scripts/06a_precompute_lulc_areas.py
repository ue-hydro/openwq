#!/usr/bin/env python3
# Copyright 2026, Diogo Costa, diogo.costa@uevora.pt
# This file is part of OpenWQ model.
#
# This program, openWQ, is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) later version.

"""
06a_precompute_lulc_areas.py
=============================

Phase 1 of the optimized SS JSON generation pipeline.

Precomputes LULC areas per HRU for all basins from Copernicus ESA-CCI
Land Cover NetCDFs, storing results as one CSV per basin.

KEY OPTIMIZATION: Year-outer loop.
  - Opens each global NetCDF ONCE per year (not 111 or 333 times)
  - Extracts a regional subset covering all Canadian basins (~385 MB in memory)
  - Clips to all 111 basins in one pass per year
  - No intermediate GeoTIFF files written to disk

PERFORMANCE:
  - Old approach: ~55 hours (333 variants × 17 years × 2.2GB reads each)
  - This script:  ~15-25 minutes (17 years × 1 read + 111 basin clips per year)

OUTPUT:
  {basin_dir}/shapefiles/{basin_id}_lulc_areas.csv
  Columns: GRU_ID, Year, LC_Class, Pixel_Count, Area_m2, Area_ha, Area_km2

USAGE:
  # Test on one basin
  python 06a_precompute_lulc_areas.py \\
      --prepared-dir /path/to/century_basins_prepared \\
      --copernicus-lc-dir /path/to/ESACCI-LC \\
      --single-basin CAN_01AM001_meso

  # Run all 111 basins
  python 06a_precompute_lulc_areas.py \\
      --prepared-dir /path/to/century_basins_prepared \\
      --copernicus-lc-dir /path/to/ESACCI-LC
"""

import sys
import os
import re
import glob
import time
import argparse
import numpy as np
import pandas as pd
from pathlib import Path

try:
    import xarray as xr
except ImportError:
    print("ERROR: xarray is required. Install: pip install xarray netCDF4")
    sys.exit(1)

try:
    import rasterio
    from rasterio.io import MemoryFile
    from rasterio.mask import mask as rasterio_mask
    from rasterio.transform import from_bounds
except ImportError:
    print("ERROR: rasterio is required. Install: pip install rasterio")
    sys.exit(1)

try:
    import geopandas as gpd
except ImportError:
    print("ERROR: geopandas is required. Install: pip install geopandas")
    sys.exit(1)


# WGS84 CRS in WKT format (avoids pyproj EPSG lookup issues on old PROJ)
WGS84_WKT = (
    'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,'
    'AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],'
    'PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],'
    'UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],'
    'AXIS["Latitude",NORTH],AXIS["Longitude",EAST],AUTHORITY["EPSG","4326"]]'
)

# Copernicus ESA-CCI nominal pixel resolution
PIXEL_AREA_M2 = 300 * 300  # 300m × 300m = 90,000 m²


def extract_year_from_filename(filename):
    """Extract year from Copernicus LULC filename.

    Example: ESACCI-LC-L4-LCCS-Map-300m-P1Y-1993-v2.0.7cds.nc → 1993
    """
    match = re.search(r'P1Y-(\d{4})', filename)
    if match:
        return int(match.group(1))
    # Fallback: find any 4-digit year between 1990-2025
    match = re.search(r'(199\d|200\d|201\d|202\d)', filename)
    if match:
        return int(match.group(1))
    return None


def load_gpkg_with_crs_fallback(gpkg_path):
    """Load a GeoPackage file with fallback for CRS resolution issues.

    Mirrors the fiona fallback in Gen_SS_Driver.py lines 178-204.
    """
    try:
        gdf = gpd.read_file(gpkg_path)
        return gdf
    except Exception as crs_err:
        if 'Invalid projection' in str(crs_err) or 'CRSError' in type(crs_err).__name__:
            import fiona
            from shapely.geometry import shape

            with fiona.open(gpkg_path) as src:
                crs_wkt = src.crs_wkt
                records = []
                for feat in src:
                    props = dict(feat['properties'])
                    props['geometry'] = shape(feat['geometry'])
                    records.append(props)

            gdf = gpd.GeoDataFrame(records, geometry='geometry')
            if crs_wkt:
                try:
                    gdf.crs = crs_wkt
                except Exception:
                    pass  # CRS couldn't be set, proceed without it
            return gdf
        else:
            raise


def find_basin_gpkg(basin_dir, basin_id):
    """Find the basin HRU GeoPackage file."""
    gpkg_path = os.path.join(basin_dir, "shapefiles", f"{basin_id}_basinHru.gpkg")
    if os.path.exists(gpkg_path):
        return gpkg_path
    # Fallback: search for any .gpkg in shapefiles/
    gpkg_files = glob.glob(os.path.join(basin_dir, "shapefiles", "*basinHru*.gpkg"))
    if gpkg_files:
        return gpkg_files[0]
    return None


def compute_super_bbox(basin_data, buffer_deg=0.5):
    """Compute the bounding box of all basins with a buffer.

    Returns dict with 'min_lat', 'max_lat', 'min_lon', 'max_lon'.
    """
    all_bounds = []
    for bdata in basin_data.values():
        bounds = bdata['gdf'].total_bounds  # (minx, miny, maxx, maxy) = (min_lon, min_lat, max_lon, max_lat)
        all_bounds.append(bounds)

    all_bounds = np.array(all_bounds)
    return {
        'min_lon': all_bounds[:, 0].min() - buffer_deg,
        'min_lat': all_bounds[:, 1].min() - buffer_deg,
        'max_lon': all_bounds[:, 2].max() + buffer_deg,
        'max_lat': all_bounds[:, 3].max() + buffer_deg,
    }


def calculate_lulc_areas_from_src(raster_src, gdf, hru_col, year, pixel_area_m2):
    """Calculate LULC areas per HRU from an already-open rasterio dataset.

    This is a streamlined version of OptimizedLULCAnalyzer.calculate_lulc_areas()
    (Gen_SS_Driver.py lines 375-473) that operates on an in-memory rasterio source.

    Parameters
    ----------
    raster_src : rasterio.DatasetReader
        Already-open rasterio dataset (can be MemoryFile)
    gdf : GeoDataFrame
        Basin HRU polygons
    hru_col : str
        Name of the HRU ID column (e.g., 'GRU_ID')
    year : int
        Year for this LULC data
    pixel_area_m2 : float
        Area per pixel in m²

    Returns
    -------
    pd.DataFrame with columns: [hru_col, Year, LC_Class, Pixel_Count, Area_m2, Area_ha, Area_km2]
    """
    results = []
    nodata = raster_src.nodata

    for idx, row in gdf.iterrows():
        hru_id = row[hru_col]
        geom = row.geometry

        if geom is None or geom.is_empty:
            continue

        try:
            out_image, _ = rasterio_mask(raster_src, [geom], crop=True, nodata=nodata)
            lc_values = out_image[0].flatten()

            # Remove nodata and NaN
            if nodata is not None:
                lc_values = lc_values[lc_values != nodata]
            lc_values = lc_values[~np.isnan(lc_values.astype(float))]
            lc_values = lc_values[lc_values > 0]  # Remove 0 (typically nodata for LC)

            if len(lc_values) == 0:
                continue

            unique, counts = np.unique(lc_values, return_counts=True)
            for lc_class, count in zip(unique, counts):
                area_m2 = int(count) * pixel_area_m2
                results.append({
                    hru_col: hru_id,
                    'Year': year,
                    'LC_Class': int(lc_class),
                    'Pixel_Count': int(count),
                    'Area_m2': area_m2,
                    'Area_ha': area_m2 / 10_000,
                    'Area_km2': area_m2 / 1_000_000,
                })

        except Exception as e:
            # Geometry might not intersect with raster extent
            if 'Input shapes do not overlap raster' in str(e):
                continue
            print(f"    Warning: Error processing HRU {hru_id}: {e}")

    return pd.DataFrame(results)


def main():
    parser = argparse.ArgumentParser(
        description="Phase 1: Precompute LULC areas per HRU from Copernicus ESA-CCI NetCDFs"
    )

    parser.add_argument('--prepared-dir', required=True,
                       help='Path to century_basins_prepared directory')
    parser.add_argument('--copernicus-lc-dir', required=True,
                       help='Path to Copernicus LULC NetCDF directory')
    parser.add_argument('--period', nargs=2, type=int, default=[1993, 2019],
                       metavar=('START', 'END'),
                       help='Year range (default: 1993 2009)')
    parser.add_argument('--single-basin', default=None,
                       help='Process a single basin (for testing)')
    parser.add_argument('--resume', action='store_true',
                       help='Skip basins that already have complete CSV')
    parser.add_argument('--force', action='store_true',
                       help='Overwrite existing CSVs')
    parser.add_argument('--dry-run', action='store_true',
                       help='Preview what would be processed')

    args = parser.parse_args()

    prepared_dir = os.path.abspath(args.prepared_dir)
    copernicus_lc_dir = os.path.abspath(args.copernicus_lc_dir)
    year_start, year_end = args.period

    # ──────────────────────────────────────────────────────────────
    # 1. Find NetCDF files
    # ──────────────────────────────────────────────────────────────
    nc_files = sorted(glob.glob(os.path.join(copernicus_lc_dir, "ESACCI-LC-*.nc")))
    if not nc_files:
        nc_files = sorted(glob.glob(os.path.join(copernicus_lc_dir, "*.nc")))
    if not nc_files:
        print(f"ERROR: No NetCDF files found in {copernicus_lc_dir}")
        sys.exit(1)

    # Filter by period
    nc_year_map = {}
    for nc in nc_files:
        yr = extract_year_from_filename(os.path.basename(nc))
        if yr and year_start <= yr <= year_end:
            nc_year_map[yr] = nc

    if not nc_year_map:
        print(f"ERROR: No NetCDFs found for period {year_start}-{year_end}")
        sys.exit(1)

    sorted_years = sorted(nc_year_map.keys())

    # ──────────────────────────────────────────────────────────────
    # 2. Find basin directories
    # ──────────────────────────────────────────────────────────────
    if args.single_basin:
        basin_dirs = [os.path.join(prepared_dir, f"basin_{args.single_basin}")]
        if not os.path.isdir(basin_dirs[0]):
            print(f"ERROR: Basin directory not found: {basin_dirs[0]}")
            sys.exit(1)
    else:
        basin_dirs = sorted(glob.glob(os.path.join(prepared_dir, "basin_*")))

    if not basin_dirs:
        print(f"ERROR: No basin directories found in {prepared_dir}")
        sys.exit(1)

    # ──────────────────────────────────────────────────────────────
    # 3. Preload all basin GeoPackages
    # ──────────────────────────────────────────────────────────────
    print("=" * 70)
    print("PHASE 1: PRECOMPUTE LULC AREAS PER BASIN")
    print("=" * 70)
    print(f"\nPrepared dir:   {prepared_dir}")
    print(f"LULC dir:       {copernicus_lc_dir}")
    print(f"Period:         {year_start}-{year_end}")
    print(f"NetCDF files:   {len(nc_year_map)} (years: {sorted_years[0]}-{sorted_years[-1]})")
    print(f"Basins:         {len(basin_dirs)}")
    print()

    print("Loading basin shapefiles...")
    basin_data = {}
    load_errors = []

    for basin_dir in basin_dirs:
        basin_name = os.path.basename(basin_dir)
        basin_id = basin_name.replace("basin_", "")

        gpkg_path = find_basin_gpkg(basin_dir, basin_id)
        if not gpkg_path:
            load_errors.append(f"  {basin_id}: No .gpkg found")
            continue

        csv_path = os.path.join(basin_dir, "shapefiles", f"{basin_id}_lulc_areas.csv")

        # Check resume/force
        if os.path.exists(csv_path) and not args.force:
            if args.resume:
                # Check if CSV has all expected years
                try:
                    existing = pd.read_csv(csv_path)
                    existing_years = set(existing['Year'].unique())
                    if existing_years >= set(sorted_years):
                        print(f"  {basin_id}: CSV already complete ({len(existing_years)} years), skipping")
                        continue
                    else:
                        missing = set(sorted_years) - existing_years
                        print(f"  {basin_id}: CSV has {len(existing_years)} years, missing {sorted(missing)}")
                except Exception:
                    pass  # CSV is corrupted/unreadable, re-process
            else:
                # Default: skip if CSV exists at all
                size_kb = os.path.getsize(csv_path) / 1024
                if size_kb > 1:
                    print(f"  {basin_id}: CSV already exists ({size_kb:.0f} KB), skipping (use --force to overwrite)")
                    continue

        try:
            gdf = load_gpkg_with_crs_fallback(gpkg_path)
            n_hrus = len(gdf)
            print(f"  {basin_id}: {n_hrus} HRUs loaded")

            basin_data[basin_id] = {
                'gdf': gdf,
                'basin_dir': basin_dir,
                'gpkg_path': gpkg_path,
                'csv_path': csv_path,
                'results': [],  # accumulate DataFrames across years
            }
        except Exception as e:
            load_errors.append(f"  {basin_id}: {e}")

    if load_errors:
        print(f"\nWarnings ({len(load_errors)} basins with issues):")
        for err in load_errors:
            print(err)

    if not basin_data:
        print("\nNo basins to process (all skipped or errored).")
        return

    print(f"\nBasins to process: {len(basin_data)}")
    total_hrus = sum(len(bd['gdf']) for bd in basin_data.values())
    print(f"Total HRUs:        {total_hrus}")

    if args.dry_run:
        print(f"\n*** DRY RUN — would process {len(nc_year_map)} years × {len(basin_data)} basins ***")
        for yr in sorted_years:
            print(f"  Year {yr}: {os.path.basename(nc_year_map[yr])}")
        return

    # ──────────────────────────────────────────────────────────────
    # 4. Compute super-bounding-box
    # ──────────────────────────────────────────────────────────────
    bbox = compute_super_bbox(basin_data, buffer_deg=0.5)
    print(f"\nSuper bounding box (all basins + 0.5° buffer):")
    print(f"  Lat: {bbox['min_lat']:.2f} to {bbox['max_lat']:.2f}")
    print(f"  Lon: {bbox['min_lon']:.2f} to {bbox['max_lon']:.2f}")

    lat_range = bbox['max_lat'] - bbox['min_lat']
    lon_range = bbox['max_lon'] - bbox['min_lon']
    # Approx pixels: lat_range / 0.002778° per pixel, lon_range / 0.002778°
    est_rows = int(lat_range / 0.002778)
    est_cols = int(lon_range / 0.002778)
    est_mb = (est_rows * est_cols) / 1e6
    print(f"  Estimated regional subset: {est_rows} × {est_cols} pixels ({est_mb:.0f} MB)")

    # ──────────────────────────────────────────────────────────────
    # 5. Year-outer loop: process each NetCDF once for all basins
    # ──────────────────────────────────────────────────────────────
    t_total_start = time.time()

    for year_idx, year in enumerate(sorted_years):
        nc_path = nc_year_map[year]
        nc_name = os.path.basename(nc_path)

        print(f"\n{'=' * 60}")
        print(f"[{year_idx+1}/{len(sorted_years)}] Year {year}: {nc_name}")
        print(f"{'=' * 60}")

        t_year_start = time.time()

        # ── 5a. Open NetCDF and extract regional subset ──
        print(f"  Opening NetCDF...", end='', flush=True)
        ds = xr.open_dataset(nc_path)
        print(" ✓")

        # Find the land cover variable
        lc_var = None
        for name in ['lccs_class', 'land_cover', 'LC', 'lc', 'Band1']:
            if name in ds.data_vars:
                lc_var = name
                break
        if lc_var is None:
            print(f"  ERROR: No land cover variable found. Vars: {list(ds.data_vars)}")
            ds.close()
            continue

        # Get coordinate names
        lat_name = None
        lon_name = None
        for name in ds.coords:
            if name.lower() in ('lat', 'latitude', 'y'):
                lat_name = name
            elif name.lower() in ('lon', 'longitude', 'x'):
                lon_name = name

        if not lat_name or not lon_name:
            print(f"  ERROR: Could not identify lat/lon coordinates. Coords: {list(ds.coords)}")
            ds.close()
            continue

        # Check if lat is descending (typical for global NetCDFs)
        lats = ds[lat_name].values
        lat_descending = lats[0] > lats[-1]

        # Extract regional subset
        print(f"  Extracting regional subset...", end='', flush=True)
        if lat_descending:
            lat_slice = slice(bbox['max_lat'], bbox['min_lat'])
        else:
            lat_slice = slice(bbox['min_lat'], bbox['max_lat'])
        lon_slice = slice(bbox['min_lon'], bbox['max_lon'])

        # Select the variable (handle time dimension if present)
        da = ds[lc_var]
        if 'time' in da.dims:
            da = da.isel(time=0)

        regional = da.sel({lat_name: lat_slice, lon_name: lon_slice})
        regional_data = regional.values.astype(np.uint8)  # LC classes are 0-255
        regional_lats = regional[lat_name].values
        regional_lons = regional[lon_name].values
        ds.close()

        print(f" ✓ ({regional_data.shape[0]} × {regional_data.shape[1]} pixels, "
              f"{regional_data.nbytes / 1e6:.0f} MB)")

        # ── 5b. Compute geotransform for regional subset ──
        # Resolution from coordinate spacing
        if len(regional_lats) > 1:
            lat_res = abs(regional_lats[1] - regional_lats[0])
        else:
            lat_res = 0.002778  # ~300m default
        if len(regional_lons) > 1:
            lon_res = abs(regional_lons[1] - regional_lons[0])
        else:
            lon_res = 0.002778

        # Upper-left corner (pixel center → edge)
        if lat_descending:
            ul_lat = regional_lats[0] + lat_res / 2
        else:
            ul_lat = regional_lats[-1] + lat_res / 2
        ul_lon = regional_lons[0] - lon_res / 2

        regional_transform = rasterio.transform.from_origin(
            west=ul_lon,
            north=ul_lat,
            xsize=lon_res,
            ysize=lat_res
        )

        # ── 5c. Create in-memory rasterio dataset ──
        print(f"  Creating in-memory raster...", end='', flush=True)
        memfile = MemoryFile()
        mem_dataset = memfile.open(
            driver='GTiff',
            height=regional_data.shape[0],
            width=regional_data.shape[1],
            count=1,
            dtype=regional_data.dtype,
            crs=WGS84_WKT,
            transform=regional_transform,
            nodata=0,
        )
        mem_dataset.write(regional_data, 1)
        print(" ✓")

        # ── 5d. Process each basin ──
        print(f"  Processing {len(basin_data)} basins: ", end='', flush=True)
        basins_ok = 0
        basins_err = 0

        for basin_idx, (basin_id, bdata) in enumerate(basin_data.items()):
            try:
                df = calculate_lulc_areas_from_src(
                    raster_src=mem_dataset,
                    gdf=bdata['gdf'],
                    hru_col='GRU_ID',
                    year=year,
                    pixel_area_m2=PIXEL_AREA_M2,
                )
                if len(df) > 0:
                    bdata['results'].append(df)
                basins_ok += 1
            except Exception as e:
                basins_err += 1
                if basins_err <= 3:
                    print(f"\n    Warning: {basin_id}: {e}")

            # Progress indicator
            if (basin_idx + 1) % 10 == 0:
                print(f"{basin_idx+1}", end='', flush=True)
            elif (basin_idx + 1) % 5 == 0:
                print(".", end='', flush=True)

        # ── 5e. Cleanup ──
        mem_dataset.close()
        memfile.close()
        del regional_data  # Free memory

        t_year_elapsed = time.time() - t_year_start
        print(f" ✓ ({basins_ok} ok, {basins_err} err, {t_year_elapsed:.1f}s)")

    # ──────────────────────────────────────────────────────────────
    # 6. Write output CSVs
    # ──────────────────────────────────────────────────────────────
    print(f"\n{'=' * 60}")
    print("Writing CSV files")
    print(f"{'=' * 60}")

    csvs_written = 0
    csvs_empty = 0

    for basin_id, bdata in basin_data.items():
        if not bdata['results']:
            print(f"  {basin_id}: No LULC data found (empty results)")
            csvs_empty += 1
            continue

        combined_df = pd.concat(bdata['results'], ignore_index=True)

        # Sort by GRU_ID, Year, LC_Class
        combined_df = combined_df.sort_values(['GRU_ID', 'Year', 'LC_Class']).reset_index(drop=True)

        csv_path = bdata['csv_path']
        os.makedirs(os.path.dirname(csv_path), exist_ok=True)
        combined_df.to_csv(csv_path, index=False)

        n_years = combined_df['Year'].nunique()
        n_hrus = combined_df['GRU_ID'].nunique()
        size_kb = os.path.getsize(csv_path) / 1024
        print(f"  {basin_id}: {len(combined_df)} rows, {n_years} years, "
              f"{n_hrus} HRUs ({size_kb:.0f} KB)")
        csvs_written += 1

    # ──────────────────────────────────────────────────────────────
    # 7. Summary
    # ──────────────────────────────────────────────────────────────
    t_total_elapsed = time.time() - t_total_start

    print(f"\n{'=' * 70}")
    print("PHASE 1 SUMMARY")
    print(f"{'=' * 70}")
    print(f"\nYears processed:  {len(sorted_years)} ({sorted_years[0]}-{sorted_years[-1]})")
    print(f"Basins processed: {len(basin_data)}")
    print(f"CSVs written:     {csvs_written}")
    if csvs_empty > 0:
        print(f"CSVs empty:       {csvs_empty}")
    print(f"Total time:       {t_total_elapsed:.1f}s ({t_total_elapsed/60:.1f} min)")
    print(f"\nNext step: Run 06b_generate_ss_from_csv.py to generate SS JSONs from these CSVs")


if __name__ == '__main__':
    main()
