#!/usr/bin/env python3
# Copyright 2026, Diogo Costa, diogo.costa@uevora.pt
# This file is part of OpenWQ model.
#
# This program, openWQ, is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) later version.

"""
05_extract_shapefiles.py
========================

Extract SUMMA HRU polygons and mizuRoute river network polylines from
century basin tar.gz archives and save them as GeoPackage (.gpkg) files.

For each basin, creates:
  century_basins_prepared/basin_{ID}/
    shapefiles/
      hru_polygons.gpkg          # SUMMA HRU/GRU catchment polygons
      river_network.gpkg         # mizuRoute river segment polylines

USAGE:
  python 05_extract_shapefiles.py \\
      --basins-dir /path/to/century_basins \\
      --prepared-dir /path/to/century_basins_prepared \\
      [--basin-filter CAN]  # optional: process only CAN or USA basins
"""

import sys
import os
import argparse
import glob
import tarfile
import tempfile
from pathlib import Path
from typing import Optional, Tuple

try:
    import fiona
    HAS_FIONA = True
except ImportError:
    HAS_FIONA = False

try:
    import geopandas as gpd
    HAS_GEOPANDAS = True
except ImportError:
    HAS_GEOPANDAS = False


def parse_basin_id(tar_path: str) -> str:
    """Extract basin ID from tar.gz filename.

    e.g. 'domain_CAN_01AM001_meso.tar.gz' -> 'CAN_01AM001_meso'
    """
    fname = Path(tar_path).name
    basename = fname.replace('.tar.gz', '')
    parts = basename.split('_')
    if len(parts) >= 4 and parts[0] == 'domain':
        return f"{parts[1]}_{parts[2]}_{parts[3]}"
    return basename


def find_shapefile_members(tar, pattern_subdir: str):
    """Find shapefile component members (shp, shx, dbf, prj, cpg) in a tar archive
    that match a given subdirectory pattern.

    Parameters
    ----------
    tar : tarfile.TarFile
        Open tar archive
    pattern_subdir : str
        Subdirectory to search for (e.g., 'shapefiles/catchment' or 'shapefiles/river_network')

    Returns
    -------
    list of tarfile.TarInfo
        Matching members
    """
    shp_extensions = {'.shp', '.shx', '.dbf', '.prj', '.cpg'}
    members = []
    for m in tar.getmembers():
        if pattern_subdir in m.name:
            ext = os.path.splitext(m.name)[1].lower()
            if ext in shp_extensions:
                members.append(m)
    return members


def find_shp_file(directory: str, pattern_subdir: str) -> Optional[str]:
    """Find a .shp file in extracted directory matching a subdirectory pattern.

    Parameters
    ----------
    directory : str
        Root of extracted archive
    pattern_subdir : str
        Subdirectory name to look in (e.g., 'catchment', 'river_network')

    Returns
    -------
    str or None
        Path to .shp file, or None if not found
    """
    matches = glob.glob(
        os.path.join(directory, '**', pattern_subdir, '*.shp'),
        recursive=True
    )
    if matches:
        return matches[0]
    return None


def _get_wkt_from_prj(shp_path: str) -> Optional[str]:
    """Read WKT CRS string from .prj sidecar file.

    Parameters
    ----------
    shp_path : str
        Path to .shp file

    Returns
    -------
    str or None
        WKT CRS string, or None if .prj not found
    """
    prj_path = os.path.splitext(shp_path)[0] + '.prj'
    if os.path.exists(prj_path):
        with open(prj_path, 'r') as f:
            return f.read().strip()
    return None


# WGS84 WKT fallback (all century basins are in EPSG:4326)
WGS84_WKT = (
    'GEOGCS["GCS_WGS_1984",'
    'DATUM["D_WGS_1984",'
    'SPHEROID["WGS_1984",6378137.0,298.257223563]],'
    'PRIMEM["Greenwich",0.0],'
    'UNIT["Degree",0.0174532925199433]]'
)


def convert_shp_to_gpkg(shp_path: str, gpkg_path: str) -> bool:
    """Convert a shapefile to GeoPackage format.

    Uses fiona as primary backend to avoid pyproj EPSG database issues
    on systems with old PROJ installations. Falls back to geopandas with
    WKT-based CRS workaround.

    Parameters
    ----------
    shp_path : str
        Path to input .shp file
    gpkg_path : str
        Path to output .gpkg file

    Returns
    -------
    bool
        True if conversion succeeded
    """
    os.makedirs(os.path.dirname(gpkg_path), exist_ok=True)

    # Read WKT CRS from .prj file (avoids EPSG lookup)
    wkt_crs = _get_wkt_from_prj(shp_path) or WGS84_WKT

    if HAS_FIONA:
        # Use fiona directly with crs_wkt (WKT string) to bypass EPSG database
        with fiona.open(shp_path) as src:
            schema = src.schema.copy()
            crs_wkt_str = src.crs_wkt  # WKT string avoids PROJ EPSG lookup

            # Pre-scan features to detect mixed geometry types
            # (e.g., Polygon + MultiPolygon → promote to Multi)
            features = list(src)

        geom_types = set()
        for f in features:
            if f.get('geometry') and f['geometry'].get('type'):
                geom_types.add(f['geometry']['type'])

        # Promote schema to Multi geometry if mixed types found
        declared_geom = schema.get('geometry', '')
        if 'MultiPolygon' in geom_types and declared_geom == 'Polygon':
            schema['geometry'] = 'MultiPolygon'
        elif 'MultiLineString' in geom_types and declared_geom == 'LineString':
            schema['geometry'] = 'MultiLineString'
        elif 'MultiPoint' in geom_types and declared_geom == 'Point':
            schema['geometry'] = 'MultiPoint'

        target_geom = schema['geometry']

        with fiona.open(gpkg_path, 'w', driver='GPKG',
                       schema=schema, crs=crs_wkt_str) as dst:
            for feature in features:
                geom = feature.get('geometry')
                if geom and geom.get('type'):
                    # Promote single geometries to Multi if schema is Multi
                    if target_geom == 'MultiPolygon' and geom['type'] == 'Polygon':
                        geom = {
                            'type': 'MultiPolygon',
                            'coordinates': [geom['coordinates']]
                        }
                        feature['geometry'] = geom
                    elif target_geom == 'MultiLineString' and geom['type'] == 'LineString':
                        geom = {
                            'type': 'MultiLineString',
                            'coordinates': [geom['coordinates']]
                        }
                        feature['geometry'] = geom
                    elif target_geom == 'MultiPoint' and geom['type'] == 'Point':
                        geom = {
                            'type': 'MultiPoint',
                            'coordinates': [geom['coordinates']]
                        }
                        feature['geometry'] = geom
                dst.write(feature)
        return True
    elif HAS_GEOPANDAS:
        # Geopandas with WKT CRS workaround
        import pyproj
        gdf = gpd.read_file(shp_path)
        # Override CRS with WKT to avoid EPSG lookup
        gdf.crs = pyproj.CRS.from_wkt(wkt_crs)
        gdf.to_file(gpkg_path, driver="GPKG")
        return True
    else:
        raise RuntimeError(
            "Neither fiona nor geopandas is available. "
            "Install one of them: pip install geopandas"
        )


def process_basin(
    tar_path: str,
    prepared_dir: str,
) -> dict:
    """Extract shapefiles from a basin tar.gz and convert to GeoPackage.

    Parameters
    ----------
    tar_path : str
        Path to the basin tar.gz archive
    prepared_dir : str
        Path to century_basins_prepared output directory

    Returns
    -------
    dict
        Result with basin_id, status, and any errors
    """
    basin_id = parse_basin_id(tar_path)
    basin_output_dir = os.path.join(prepared_dir, f"basin_{basin_id}")
    shapefiles_dir = os.path.join(basin_output_dir, "shapefiles")

    result = {
        'basin_id': basin_id,
        'hru_polygons': False,
        'river_network': False,
        'errors': [],
    }

    with tempfile.TemporaryDirectory(prefix=f"shp_{basin_id}_") as tmp_dir:
        try:
            with tarfile.open(tar_path, 'r:gz') as tar:
                # Extract catchment (HRU) shapefiles
                catchment_members = find_shapefile_members(tar, 'shapefiles/catchment')
                # Extract river network shapefiles
                river_members = find_shapefile_members(tar, 'shapefiles/river_network')

                all_members = catchment_members + river_members
                if not all_members:
                    result['errors'].append("No shapefile members found in archive")
                    return result

                tar.extractall(tmp_dir, members=all_members)

        except Exception as e:
            result['errors'].append(f"Extraction error: {e}")
            return result

        # Convert HRU polygons
        hru_shp = find_shp_file(tmp_dir, 'catchment')
        if hru_shp:
            try:
                gpkg_path = os.path.join(shapefiles_dir, f"{basin_id}_basinHru.gpkg")
                convert_shp_to_gpkg(hru_shp, gpkg_path)
                result['hru_polygons'] = True
            except Exception as e:
                result['errors'].append(f"HRU conversion error: {e}")
        else:
            result['errors'].append("No HRU shapefile found in catchment/")

        # Convert river network
        river_shp = find_shp_file(tmp_dir, 'river_network')
        if river_shp:
            try:
                gpkg_path = os.path.join(shapefiles_dir, f"{basin_id}_riverNetwork.gpkg")
                convert_shp_to_gpkg(river_shp, gpkg_path)
                result['river_network'] = True
            except Exception as e:
                result['errors'].append(f"River network conversion error: {e}")
        else:
            result['errors'].append("No river network shapefile found in river_network/")

    return result


def main():
    parser = argparse.ArgumentParser(
        description="Extract HRU polygons and river network polylines as GeoPackage",
    )

    parser.add_argument('--basins-dir', required=True,
                       help='Directory containing basin tar.gz archives')
    parser.add_argument('--prepared-dir', required=True,
                       help='Output directory (century_basins_prepared)')
    parser.add_argument('--basin-filter', default=None,
                       help='Filter basins by country code (CAN, USA)')

    args = parser.parse_args()

    # Find all basin archives
    tar_files = sorted(glob.glob(os.path.join(args.basins_dir, 'domain_*.tar.gz')))

    if args.basin_filter:
        tar_files = [t for t in tar_files
                    if f"domain_{args.basin_filter}_" in os.path.basename(t)]

    if not tar_files:
        print(f"ERROR: No basin archives found in {args.basins_dir}")
        sys.exit(1)

    print("=" * 70)
    print("SHAPEFILE EXTRACTION (to GeoPackage)")
    print("=" * 70)
    print(f"\nBasins:       {len(tar_files)}")
    print(f"Source:       {args.basins_dir}")
    print(f"Output:       {args.prepared_dir}")
    if HAS_FIONA:
        print(f"Backend:      fiona")
    elif HAS_GEOPANDAS:
        print(f"Backend:      geopandas (WKT CRS workaround)")
    else:
        print("ERROR: Neither fiona nor geopandas available!")
        sys.exit(1)
    print()

    # Process all basins
    success_hru = 0
    success_river = 0
    errors = []

    for i, tar_path in enumerate(tar_files, 1):
        basin_id = parse_basin_id(tar_path)
        print(f"[{i}/{len(tar_files)}] {basin_id}")

        result = process_basin(tar_path, args.prepared_dir)

        if result['hru_polygons']:
            print(f"  [✓] {basin_id}_basinHru.gpkg")
            success_hru += 1
        if result['river_network']:
            print(f"  [✓] {basin_id}_riverNetwork.gpkg")
            success_river += 1
        if result['errors']:
            for err in result['errors']:
                print(f"  [✗] {err}")
            errors.append((basin_id, result['errors']))

    # Summary
    print()
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Total basins:         {len(tar_files)}")
    print(f"HRU polygons:         {success_hru}/{len(tar_files)}")
    print(f"River network:        {success_river}/{len(tar_files)}")
    print(f"Basins with errors:   {len(errors)}")

    if errors:
        print("\nErrors:")
        for basin_id, errs in errors:
            for err in errs:
                print(f"  {basin_id}: {err}")

    print(f"\nOutput: {args.prepared_dir}/basin_*/shapefiles/")


if __name__ == '__main__':
    main()
