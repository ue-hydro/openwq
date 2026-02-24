#!/usr/bin/env python3
# Copyright 2026, Diogo Costa, diogo.costa@uevora.pt
# This file is part of OpenWQ model.
#
# This program, openWQ, is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) later version.

"""
06_generate_ss_jsons.py
========================

Generate the missing source/sink JSON files for all 333 variant configs
using Copernicus ESA-CCI Land Cover NetCDFs.

This is a lightweight script that:
  1. Iterates over existing variant directories in century_basins_prepared/
  2. Reads config_raw.json for each variant to extract SS parameters
  3. Uses the permanent .gpkg shapefiles (with GRU_ID column)
  4. Calls Gen_SS_Driver.set_ss_from_copernicus_lulc_with_loads() directly
  5. Generates only the missing openWQ_SS_*.json files (+ clipped rasters)

This avoids re-running the full 02_generate_openwq_configs.py which would
re-extract tar.gz archives and overwrite all existing BGC/transport configs.

USAGE:
  # Test on one basin first
  python 06_generate_ss_jsons.py \\
      --prepared-dir /path/to/century_basins_prepared \\
      --copernicus-lc-dir /path/to/ESACCI-LC \\
      --single-basin CAN_01AM001_meso

  # Run all basins
  python 06_generate_ss_jsons.py \\
      --prepared-dir /path/to/century_basins_prepared \\
      --copernicus-lc-dir /path/to/ESACCI-LC

  # Force overwrite existing SS JSONs
  python 06_generate_ss_jsons.py --prepared-dir ... --copernicus-lc-dir ... --force

  # Dry run (check what would be generated)
  python 06_generate_ss_jsons.py --prepared-dir ... --copernicus-lc-dir ... --dry-run
"""

import sys
import os
import json
import glob
import argparse
import traceback
from pathlib import Path

# Add Model_Config to path for Gen_SS_Driver import
SCRIPT_DIR = Path(__file__).resolve().parent
CALIBRATION_DIR = SCRIPT_DIR.parent
OPENWQ_SCRIPTS_DIR = CALIBRATION_DIR.parent  # supporting_scripts/
MODEL_CONFIG_DIR = OPENWQ_SCRIPTS_DIR / "Model_Config"
CONFIG_LIB_DIR = MODEL_CONFIG_DIR / "config_support_lib"

sys.path.insert(0, str(MODEL_CONFIG_DIR))
sys.path.insert(0, str(CONFIG_LIB_DIR))

try:
    import Gen_SS_Driver as ssJSON_lib
    SS_LIB_AVAILABLE = True
except ImportError as e:
    print(f"WARNING: Could not import Gen_SS_Driver: {e}")
    print("  The script requires the OpenWQ Model_Config library.")
    print(f"  Expected at: {CONFIG_LIB_DIR}")
    SS_LIB_AVAILABLE = False


# Default SS JSON filename (must match what openWQ_master.json references)
SS_JSON_FILENAME = "openWQ_SS_using_copernicus_lulc_with_static_coeff.json"


def find_basin_gpkg(basin_dir, basin_id):
    """Find the basin HRU GeoPackage file.

    Returns the path to {basin_id}_basinHru.gpkg in the shapefiles/ directory.
    """
    gpkg_path = os.path.join(basin_dir, "shapefiles", f"{basin_id}_basinHru.gpkg")
    if os.path.exists(gpkg_path):
        return gpkg_path

    # Fallback: search for any .gpkg in shapefiles/
    gpkg_files = glob.glob(os.path.join(basin_dir, "shapefiles", "*basinHru*.gpkg"))
    if gpkg_files:
        return gpkg_files[0]

    return None


def generate_ss_for_variant(variant_dir, basin_id, basin_dir,
                            copernicus_lc_dir, force=False, dry_run=False):
    """Generate the SS JSON for a single variant directory.

    Parameters
    ----------
    variant_dir : str
        Path to variant directory (e.g., .../variant_A_nitrogen_full/)
    basin_id : str
        Basin identifier (e.g., CAN_01AM001_meso)
    basin_dir : str
        Path to basin directory (e.g., .../basin_CAN_01AM001_meso/)
    copernicus_lc_dir : str
        Path to Copernicus LULC NetCDF directory
    force : bool
        If True, overwrite existing SS JSON
    dry_run : bool
        If True, don't write files

    Returns
    -------
    str or None
        Path to generated SS JSON, or None on failure
    """
    variant_name = os.path.basename(variant_dir)

    # Target SS JSON path
    openwq_in_dir = os.path.join(variant_dir, "openwq_in")
    ss_json_path = os.path.join(openwq_in_dir, SS_JSON_FILENAME)

    # Check if already exists
    if os.path.exists(ss_json_path) and not force:
        size_kb = os.path.getsize(ss_json_path) / 1024
        if size_kb > 1:  # >1 KB means it's likely a real file
            print(f"      [{variant_name}] SS JSON already exists ({size_kb:.0f} KB), skipping")
            return ss_json_path

    if dry_run:
        print(f"      [{variant_name}] Would generate: {SS_JSON_FILENAME}")
        return ss_json_path

    # Find basin HRU GeoPackage
    gpkg_path = find_basin_gpkg(basin_dir, basin_id)
    if not gpkg_path:
        print(f"      [{variant_name}] ERROR: No basinHru .gpkg found in {basin_dir}/shapefiles/")
        return None

    # Read config_raw.json for variant-specific parameters
    config_raw_path = os.path.join(variant_dir, "config_raw.json")
    if os.path.exists(config_raw_path):
        with open(config_raw_path, 'r') as f:
            config = json.load(f)
    else:
        config = {}

    # Extract parameters (override old values with correct ones)
    ss_period = config.get('ss_method_copernicus_period', [1993, 2019])
    compartment_name = config.get('ss_method_copernicus_compartment_name_for_load',
                                   'RIVER_NETWORK_REACHES')
    default_loads = config.get('ss_method_copernicus_default_loads_bool', True)
    seasonal_method = config.get('ss_method_copernicus_annual_to_seasonal_loads_method',
                                  'uniform')
    use_cellid = config.get('ss_use_cellid_mapping', True)

    # Metadata from config
    project_name = config.get('project_name', f'Century_{basin_id}')
    comment = config.get('comment', f'Auto-generated for {basin_id}')
    ss_source = config.get('ss_metadata_source', f'Auto-generated for {basin_id}')
    ss_comment = config.get('ss_metadata_comment', variant_name)

    # Ensure openwq_in directory exists
    os.makedirs(openwq_in_dir, exist_ok=True)

    # Build header comment
    json_header_comment = [
        f"// OpenWQ Source/Sink Configuration",
        f"// Project: {project_name}",
        f"// Basin: {basin_id}",
        f"// Variant: {variant_name}",
        f"// Generated by: 06_generate_ss_jsons.py",
        f"// LULC source: Copernicus ESA-CCI Land Cover ({ss_period[0]}-{ss_period[1]})",
    ]

    # Basin info (use permanent .gpkg with correct column name)
    basin_info = {
        'path_to_shp': gpkg_path,
        'mapping_key': 'GRU_ID',
    }

    print(f"      [{variant_name}] Generating SS JSON...")
    print(f"        LULC dir: {copernicus_lc_dir}")
    print(f"        Shapefile: {os.path.basename(gpkg_path)}")
    print(f"        Period: {ss_period[0]}-{ss_period[1]}")
    print(f"        Output: {ss_json_path}")

    try:
        result = ssJSON_lib.set_ss_from_copernicus_lulc_with_loads(
            ss_config_filepath=ss_json_path,
            json_header_comment=json_header_comment,
            ss_metadata_source=ss_source,
            ss_metadata_comment=ss_comment,
            ss_method_copernicus_basin_info=basin_info,
            ss_method_copernicus_nc_lc_dir=copernicus_lc_dir,
            ss_method_copernicus_period=ss_period,
            ss_method_copernicus_default_loads_bool=default_loads,
            ss_method_copernicus_compartment_name_for_load=compartment_name,
            ss_method_copernicus_openwq_h5_results_file_example_for_mapping_key=None,
            ss_method_copernicus_annual_to_seasonal_loads_method=seasonal_method,
            use_cellid_mapping=use_cellid,
        )

        if os.path.exists(ss_json_path):
            size_kb = os.path.getsize(ss_json_path) / 1024
            print(f"      [{variant_name}] SUCCESS ({size_kb:.0f} KB)")
            return ss_json_path
        else:
            print(f"      [{variant_name}] WARNING: Function returned but no file created")
            return None

    except Exception as e:
        print(f"      [{variant_name}] ERROR: {e}")
        traceback.print_exc()
        return None


def process_basin(basin_dir, copernicus_lc_dir, force=False, dry_run=False):
    """Process all variants for a single basin.

    Returns
    -------
    dict : {'basin_id': str, 'success': int, 'skipped': int, 'failed': int, 'errors': list}
    """
    basin_name = os.path.basename(basin_dir)
    basin_id = basin_name.replace("basin_", "")

    result = {
        'basin_id': basin_id,
        'success': 0,
        'skipped': 0,
        'failed': 0,
        'errors': [],
    }

    # Find all variant directories
    variant_dirs = sorted(glob.glob(os.path.join(basin_dir, "variant_*")))
    if not variant_dirs:
        print(f"  [{basin_id}] No variant directories found, skipping")
        return result

    print(f"  [{basin_id}] ({len(variant_dirs)} variants)")

    for variant_dir in variant_dirs:
        ss_path = generate_ss_for_variant(
            variant_dir=variant_dir,
            basin_id=basin_id,
            basin_dir=basin_dir,
            copernicus_lc_dir=copernicus_lc_dir,
            force=force,
            dry_run=dry_run,
        )

        if ss_path and os.path.exists(ss_path) and os.path.getsize(ss_path) > 100:
            result['success'] += 1
        elif ss_path and dry_run:
            result['success'] += 1
        elif ss_path and not os.path.exists(ss_path):
            result['failed'] += 1
            result['errors'].append(f"{os.path.basename(variant_dir)}: no file created")
        elif ss_path is None:
            result['failed'] += 1
        else:
            result['skipped'] += 1

    return result


def main():
    parser = argparse.ArgumentParser(
        description="Generate missing source/sink JSONs from Copernicus LULC data"
    )

    parser.add_argument('--prepared-dir', required=True,
                       help='Path to century_basins_prepared directory')
    parser.add_argument('--copernicus-lc-dir', required=True,
                       help='Path to Copernicus LULC NetCDF directory')
    parser.add_argument('--single-basin', default=None,
                       help='Process a single basin (for testing)')
    parser.add_argument('--force', action='store_true',
                       help='Overwrite existing SS JSON files')
    parser.add_argument('--dry-run', action='store_true',
                       help='Preview what would be generated')

    args = parser.parse_args()

    if not SS_LIB_AVAILABLE:
        print("ERROR: Gen_SS_Driver library is not available. Cannot proceed.")
        sys.exit(1)

    prepared_dir = os.path.abspath(args.prepared_dir)
    copernicus_lc_dir = os.path.abspath(args.copernicus_lc_dir)

    # Verify LULC directory has NetCDF files
    nc_files = sorted(glob.glob(os.path.join(copernicus_lc_dir, "ESACCI-LC-*.nc")))
    if not nc_files:
        # Try broader pattern
        nc_files = sorted(glob.glob(os.path.join(copernicus_lc_dir, "*.nc")))
    if not nc_files:
        print(f"ERROR: No NetCDF files found in {copernicus_lc_dir}")
        print("  Expected pattern: ESACCI-LC-*.nc")
        sys.exit(1)

    # Find basin directories
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

    # Count total variants
    total_variants = 0
    for bd in basin_dirs:
        total_variants += len(glob.glob(os.path.join(bd, "variant_*")))

    print("=" * 70)
    print("SOURCE/SINK JSON GENERATION FROM COPERNICUS LULC")
    print("=" * 70)
    print(f"\nPrepared dir:     {prepared_dir}")
    print(f"LULC dir:         {copernicus_lc_dir}")
    print(f"LULC files:       {len(nc_files)} NetCDFs")
    print(f"  First: {os.path.basename(nc_files[0])}")
    print(f"  Last:  {os.path.basename(nc_files[-1])}")
    print(f"Basins:           {len(basin_dirs)}")
    print(f"Total variants:   {total_variants}")
    if args.force:
        print(f"Mode:             FORCE (overwrite existing)")
    elif args.dry_run:
        print(f"Mode:             DRY RUN (no files written)")
    else:
        print(f"Mode:             Incremental (skip existing)")
    print()

    # Process each basin
    results = []
    for i, basin_dir in enumerate(basin_dirs):
        basin_name = os.path.basename(basin_dir)
        print(f"\n[{i+1}/{len(basin_dirs)}] {basin_name}")

        result = process_basin(
            basin_dir=basin_dir,
            copernicus_lc_dir=copernicus_lc_dir,
            force=args.force,
            dry_run=args.dry_run,
        )
        results.append(result)

    # Summary
    total_success = sum(r['success'] for r in results)
    total_failed = sum(r['failed'] for r in results)
    total_skipped = sum(r['skipped'] for r in results)
    total_errors = sum(len(r['errors']) for r in results)

    print("\n" + "=" * 70)
    print("SS JSON GENERATION SUMMARY")
    print("=" * 70)
    print(f"\nBasins processed:  {len(results)}")
    print(f"Variants total:    {total_variants}")
    print(f"  Success:         {total_success}")
    print(f"  Skipped:         {total_skipped}")
    print(f"  Failed:          {total_failed}")

    if total_errors > 0:
        print(f"\nErrors ({total_errors}):")
        for r in results:
            for err in r['errors']:
                print(f"  [{r['basin_id']}] {err}")

    # Verify output
    ss_json_count = len(glob.glob(
        os.path.join(prepared_dir, "basin_*", "variant_*", "openwq_in", "openWQ_SS_*.json")
    ))
    print(f"\nTotal SS JSON files on disk: {ss_json_count}")
    print(f"Expected:                    {total_variants}")

    if args.dry_run:
        print("\n*** DRY RUN — no files were written ***")


if __name__ == '__main__':
    main()
