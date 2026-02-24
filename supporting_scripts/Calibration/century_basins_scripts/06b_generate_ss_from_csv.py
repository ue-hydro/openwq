#!/usr/bin/env python3
# Copyright 2026, Diogo Costa, diogo.costa@uevora.pt
# This file is part of OpenWQ model.
#
# This program, openWQ, is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) later version.

"""
06b_generate_ss_from_csv.py
============================

Phase 2 of the optimized SS JSON generation pipeline.

Reads precomputed LULC area CSVs (from 06a_precompute_lulc_areas.py)
and generates the OpenWQ source/sink JSON files for all 333 variants.

NO RASTER PROCESSING — pure CSV → JSON transformation using existing
Gen_SS_Driver functions:
  - get_default_copernicus_load_coefficients()
  - calculate_nutrient_loads()
  - create_openwq_ss_json_from_loads()

PERFORMANCE:
  ~1 second per variant → ~6 minutes for all 333 variants.

PREREQUISITE:
  Run 06a_precompute_lulc_areas.py first to generate the LULC area CSVs.

USAGE:
  # Test on one basin
  python 06b_generate_ss_from_csv.py \\
      --prepared-dir /path/to/century_basins_prepared \\
      --single-basin CAN_01AM001_meso --force

  # Run all variants
  python 06b_generate_ss_from_csv.py \\
      --prepared-dir /path/to/century_basins_prepared

  # Dry run
  python 06b_generate_ss_from_csv.py \\
      --prepared-dir /path/to/century_basins_prepared --dry-run
"""

import sys
import os
import json
import glob
import time
import argparse
import traceback
from pathlib import Path

import pandas as pd

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
    print(f"  Expected at: {CONFIG_LIB_DIR}")
    SS_LIB_AVAILABLE = False

# Default SS JSON filename (must match what openWQ_master.json references)
SS_JSON_FILENAME = "openWQ_SS_using_copernicus_lulc_with_static_coeff.json"


def generate_ss_for_variant(variant_dir, basin_id, basin_dir, lulc_csv_path,
                            force=False, dry_run=False):
    """Generate the SS JSON for a single variant from precomputed LULC areas CSV.

    Parameters
    ----------
    variant_dir : str
        Path to variant directory (e.g., .../variant_A_nitrogen_full/)
    basin_id : str
        Basin identifier (e.g., CAN_01AM001_meso)
    basin_dir : str
        Path to basin directory
    lulc_csv_path : str
        Path to precomputed LULC areas CSV
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
    openwq_in_dir = os.path.join(variant_dir, "openwq_in")
    ss_json_path = os.path.join(openwq_in_dir, SS_JSON_FILENAME)

    # Check if already exists
    if os.path.exists(ss_json_path) and not force:
        size_kb = os.path.getsize(ss_json_path) / 1024
        if size_kb > 1:
            return ss_json_path  # Already exists, skip silently

    if dry_run:
        print(f"      [{variant_name}] Would generate: {SS_JSON_FILENAME}")
        return ss_json_path

    # Read config_raw.json for variant parameters
    config_raw_path = os.path.join(variant_dir, "config_raw.json")
    if os.path.exists(config_raw_path):
        with open(config_raw_path, 'r') as f:
            config = json.load(f)
    else:
        config = {}

    # Extract parameters
    compartment_name = config.get('ss_method_copernicus_compartment_name_for_load',
                                   'RIVER_NETWORK_REACHES')
    default_loads = config.get('ss_method_copernicus_default_loads_bool', True)
    seasonal_method = config.get('ss_method_copernicus_annual_to_seasonal_loads_method',
                                  'uniform')
    use_cellid = config.get('ss_use_cellid_mapping', True)

    # Metadata
    project_name = config.get('project_name', f'Century_{basin_id}')
    ss_source = config.get('ss_metadata_source', f'Auto-generated for {basin_id}')
    ss_comment = config.get('ss_metadata_comment', variant_name)
    ss_period = config.get('ss_method_copernicus_period', [1993, 2019])

    # Ensure output directories exist
    os.makedirs(openwq_in_dir, exist_ok=True)
    ss_files_dir = Path(openwq_in_dir) / 'ss_copernicus_files'
    ss_files_dir.mkdir(parents=True, exist_ok=True)

    # Build header comment
    json_header_comment = [
        f"// OpenWQ Source/Sink Configuration",
        f"// Project: {project_name}",
        f"// Basin: {basin_id}",
        f"// Variant: {variant_name}",
        f"// Generated by: 06b_generate_ss_from_csv.py",
        f"// LULC source: Copernicus ESA-CCI Land Cover ({ss_period[0]}-{ss_period[1]})",
    ]

    try:
        # 1. Read precomputed LULC areas CSV
        lulc_df = pd.read_csv(lulc_csv_path)

        # 2. Get load coefficients
        if default_loads:
            load_coefficients = ssJSON_lib.get_default_copernicus_load_coefficients()
        else:
            load_coefficients = ssJSON_lib.get_default_copernicus_load_coefficients()

        # 3. Calculate nutrient loads → writes pivot CSV to ss_files_dir
        loads_df, summary_df = ssJSON_lib.calculate_nutrient_loads(
            results_df=lulc_df,
            load_coefficients=load_coefficients,
            output_dir=ss_files_dir,
            shp_hru_id_column='GRU_ID',
        )

        # 4. Generate SS JSON from pivot CSV
        pivot_csv_path = ss_files_dir / 'nutrient_loads_pivot.csv'
        ssJSON_lib.create_openwq_ss_json_from_loads(
            pivot_loads_csv_path=pivot_csv_path,
            ss_config_filepath=ss_json_path,
            json_header_comment=json_header_comment,
            ss_method_copernicus_compartment_name_for_load=compartment_name,
            ss_method_copernicus_openwq_h5_results_file_example_for_mapping_key=None,
            shp_hru_id_column='GRU_ID',
            ss_method_copernicus_annual_to_seasonal_loads_method=seasonal_method,
            ss_metadata_comment=ss_comment,
            ss_metadata_source=ss_source,
            use_cellid_mapping=use_cellid,
        )

        if os.path.exists(ss_json_path):
            return ss_json_path
        else:
            return None

    except Exception as e:
        print(f"      [{variant_name}] ERROR: {e}")
        traceback.print_exc()
        return None


def main():
    parser = argparse.ArgumentParser(
        description="Phase 2: Generate SS JSONs from precomputed LULC area CSVs"
    )

    parser.add_argument('--prepared-dir', required=True,
                       help='Path to century_basins_prepared directory')
    parser.add_argument('--single-basin', default=None,
                       help='Process a single basin (for testing)')
    parser.add_argument('--single-variant', default=None,
                       help='Process a single variant name (e.g., variant_A_nitrogen_full)')
    parser.add_argument('--force', action='store_true',
                       help='Overwrite existing SS JSON files')
    parser.add_argument('--dry-run', action='store_true',
                       help='Preview what would be generated')

    args = parser.parse_args()

    if not SS_LIB_AVAILABLE:
        print("ERROR: Gen_SS_Driver library is not available. Cannot proceed.")
        sys.exit(1)

    prepared_dir = os.path.abspath(args.prepared_dir)

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

    # Count totals
    total_variants = 0
    for bd in basin_dirs:
        if args.single_variant:
            v_dir = os.path.join(bd, args.single_variant)
            if os.path.isdir(v_dir):
                total_variants += 1
        else:
            total_variants += len(glob.glob(os.path.join(bd, "variant_*")))

    print("=" * 70)
    print("PHASE 2: GENERATE SS JSONS FROM PRECOMPUTED LULC AREAS")
    print("=" * 70)
    print(f"\nPrepared dir:     {prepared_dir}")
    print(f"Basins:           {len(basin_dirs)}")
    print(f"Total variants:   {total_variants}")
    if args.force:
        print(f"Mode:             FORCE (overwrite existing)")
    elif args.dry_run:
        print(f"Mode:             DRY RUN")
    else:
        print(f"Mode:             Incremental (skip existing)")
    print()

    # Process each basin
    t_start = time.time()
    total_success = 0
    total_skipped = 0
    total_failed = 0
    total_no_csv = 0

    for i, basin_dir in enumerate(basin_dirs):
        basin_name = os.path.basename(basin_dir)
        basin_id = basin_name.replace("basin_", "")

        # Find precomputed LULC CSV
        lulc_csv_path = os.path.join(basin_dir, "shapefiles", f"{basin_id}_lulc_areas.csv")
        if not os.path.exists(lulc_csv_path):
            print(f"[{i+1}/{len(basin_dirs)}] {basin_id}: NO LULC CSV — run 06a first!")
            total_no_csv += 1
            continue

        # Find variant directories
        if args.single_variant:
            variant_dirs = [os.path.join(basin_dir, args.single_variant)]
            variant_dirs = [v for v in variant_dirs if os.path.isdir(v)]
        else:
            variant_dirs = sorted(glob.glob(os.path.join(basin_dir, "variant_*")))

        if not variant_dirs:
            continue

        basin_success = 0
        basin_skipped = 0
        basin_failed = 0

        for variant_dir in variant_dirs:
            variant_name = os.path.basename(variant_dir)
            ss_json_path = os.path.join(variant_dir, "openwq_in", SS_JSON_FILENAME)

            # Quick skip check (before calling generate)
            if os.path.exists(ss_json_path) and not args.force and not args.dry_run:
                size_kb = os.path.getsize(ss_json_path) / 1024
                if size_kb > 1:
                    basin_skipped += 1
                    continue

            result = generate_ss_for_variant(
                variant_dir=variant_dir,
                basin_id=basin_id,
                basin_dir=basin_dir,
                lulc_csv_path=lulc_csv_path,
                force=args.force,
                dry_run=args.dry_run,
            )

            if result and (os.path.exists(result) or args.dry_run):
                basin_success += 1
            else:
                basin_failed += 1

        # Summary line per basin
        parts = []
        if basin_success:
            parts.append(f"{basin_success} generated")
        if basin_skipped:
            parts.append(f"{basin_skipped} skipped")
        if basin_failed:
            parts.append(f"{basin_failed} failed")
        status = ", ".join(parts)
        print(f"[{i+1}/{len(basin_dirs)}] {basin_id}: {status}")

        total_success += basin_success
        total_skipped += basin_skipped
        total_failed += basin_failed

    # Summary
    t_elapsed = time.time() - t_start

    print(f"\n{'=' * 70}")
    print("PHASE 2 SUMMARY")
    print(f"{'=' * 70}")
    print(f"\nBasins processed:  {len(basin_dirs)}")
    print(f"Variants total:    {total_variants}")
    print(f"  Generated:       {total_success}")
    print(f"  Skipped:         {total_skipped}")
    print(f"  Failed:          {total_failed}")
    if total_no_csv > 0:
        print(f"  No LULC CSV:     {total_no_csv}")
    print(f"Total time:        {t_elapsed:.1f}s ({t_elapsed/60:.1f} min)")

    # Verify output
    ss_json_count = len(glob.glob(
        os.path.join(prepared_dir, "basin_*", "variant_*", "openwq_in", "openWQ_SS_*.json")
    ))
    print(f"\nTotal SS JSON files on disk: {ss_json_count}")
    print(f"Expected:                    {total_variants + total_skipped}")

    if args.dry_run:
        print("\n*** DRY RUN — no files were written ***")


if __name__ == '__main__':
    main()
