#!/usr/bin/env python3
# Copyright 2026, Diogo Costa, diogo.costa@uevora.pt
# This file is part of OpenWQ model.
#
# This program, openWQ, is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) later version.

"""
07_fixup_hpc_paths.py
======================

Rewrite local macOS paths in all generated JSON and Python configs to
HPC-compatible absolute paths before transferring to the cluster.

This rewrites:
  - openWQ_master.json  (CONFIG_FILEPATH, SS FILEPATH, MODULE_CONFIG_FILEPATH,
                         RESULTS_FOLDERPATH)
  - config_raw.json     (dir2save_input_files, path2selected_...,
                         path_to_shp, ss_method_copernicus_nc_lc_dir)
  - calibration_config_{A,B,C,D}.py  (CALIBRATION_LIB, HPC_BASE,
                                       apptainer_sif_path, apptainer_bind_path)

USAGE:
  python 07_fixup_hpc_paths.py \\
      --prepared-dir /path/to/century_basins_prepared \\
      --hpc-base /projects/F202407877CPCAA1/diogo/openwq/century_calibration \\
      --hpc-lulc-dir /projects/F202407877CPCAA1/diogo/global_wq/ESACCI_landuse/esa_cci_landcover \\
      --hpc-openwq-root /projects/F202407877CPCAA1/diogo/openwq/openwq \\
      --hpc-sif-path /projects/F202407877CPCAA1/diogo/openwq/openwq.sif \\
      [--dry-run]   # Preview changes without writing
"""

import sys
import os
import re
import json
import glob
import argparse
from pathlib import Path


# ─────────────────────────────────────────────────────────────────────
# JSON path rewriting
# ─────────────────────────────────────────────────────────────────────

def rewrite_json_paths(json_path, local_prepared_dir, hpc_base,
                       hpc_lulc_dir, hpc_openwq_root, dry_run=False):
    """Rewrite absolute paths in a JSON config file.

    Parameters
    ----------
    json_path : str
        Path to JSON file
    local_prepared_dir : str
        Local path to century_basins_prepared (e.g. /Users/.../century_basins_prepared)
    hpc_base : str
        HPC path base (e.g. /projects/.../century_calibration)
    hpc_lulc_dir : str
        HPC path to ESACCI land cover data
    hpc_openwq_root : str
        HPC path to OpenWQ source code root
    dry_run : bool
        If True, report changes but don't write

    Returns
    -------
    int
        Number of paths rewritten
    """
    with open(json_path, 'r') as f:
        content = f.read()

    original = content

    # Strip comment lines (JSON with comments from Gen_Input_Driver)
    lines = content.split('\n')
    clean_lines = [l for l in lines if not l.strip().startswith('//')]
    try:
        data = json.loads('\n'.join(clean_lines))
    except json.JSONDecodeError:
        # Try original content
        try:
            data = json.loads(content)
        except json.JSONDecodeError:
            print(f"  WARNING: Could not parse JSON: {json_path}")
            return 0

    count = 0

    def rewrite_value(val):
        nonlocal count
        if not isinstance(val, str):
            return val

        new_val = val

        # Replace local prepared dir → HPC base + basins/
        if local_prepared_dir in new_val:
            new_val = new_val.replace(local_prepared_dir, os.path.join(hpc_base, "basins"))
            count += 1

        # Replace OpenWQ source path → HPC openwq root
        # Pattern: /Users/.../openwq/supporting_scripts/...
        openwq_pattern = re.search(
            r'/[^"]*?/openwq/openwq/supporting_scripts/',
            new_val
        )
        if openwq_pattern:
            prefix = openwq_pattern.group(0)
            new_val = new_val.replace(
                prefix[:prefix.rindex('/openwq/openwq/supporting_scripts/')],
                hpc_openwq_root.rstrip('/').rsplit('/supporting_scripts', 1)[0]
                    if '/supporting_scripts' in hpc_openwq_root
                    else hpc_openwq_root
            )
            # Reconstruct: hpc_openwq_root + /supporting_scripts/...
            count += 1

        # Replace ESACCI-LC directory
        if '/data/ESACCI-LC/' in new_val:
            new_val = new_val.replace('/data/ESACCI-LC/', hpc_lulc_dir.rstrip('/') + '/')
            count += 1
        # Also catch the local ESACCI-LC dir
        esacci_local_pattern = re.search(r'/Users/[^"]*?/ESACCI-LC/', new_val)
        if esacci_local_pattern:
            new_val = new_val.replace(
                esacci_local_pattern.group(0),
                hpc_lulc_dir.rstrip('/') + '/'
            )
            count += 1

        # Replace temp dir paths (from shapefile extraction)
        temp_pattern = re.search(r'/var/folders/[^"]+', new_val)
        if temp_pattern:
            # These temp paths should point to the basin's shapefiles/ dir
            # Try to extract basin info
            basin_match = re.search(r'(CAN|USA)_\w+_\w+', new_val)
            if basin_match:
                basin_id = basin_match.group(0)
                new_val = os.path.join(
                    hpc_base, "basins", f"basin_{basin_id}",
                    "shapefiles", f"{basin_id}_basinHru.gpkg"
                )
            count += 1

        return new_val

    def walk_and_rewrite(obj):
        if isinstance(obj, dict):
            for k, v in obj.items():
                if isinstance(v, str):
                    obj[k] = rewrite_value(v)
                elif isinstance(v, (dict, list)):
                    walk_and_rewrite(v)
        elif isinstance(obj, list):
            for i, v in enumerate(obj):
                if isinstance(v, str):
                    obj[i] = rewrite_value(v)
                elif isinstance(v, (dict, list)):
                    walk_and_rewrite(v)

    walk_and_rewrite(data)

    if count > 0:
        if dry_run:
            print(f"  [DRY-RUN] {os.path.basename(json_path)}: {count} paths would be rewritten")
        else:
            # Write back with comment header if original had comments
            header_lines = []
            for line in lines:
                if line.strip().startswith('//'):
                    header_lines.append(line)
                else:
                    break

            with open(json_path, 'w') as f:
                if header_lines:
                    f.write('\n'.join(header_lines) + '\n\n')
                json.dump(data, f, indent=4)
                f.write('\n')
            print(f"  [✓] {os.path.basename(json_path)}: {count} paths rewritten")

    return count


# ─────────────────────────────────────────────────────────────────────
# Python calibration config rewriting
# ─────────────────────────────────────────────────────────────────────

def rewrite_python_paths(py_path, hpc_base, hpc_openwq_root, hpc_sif_path,
                         dry_run=False):
    """Rewrite HPC paths in calibration_config_*.py files.

    Parameters
    ----------
    py_path : str
        Path to Python calibration config
    hpc_base : str
        HPC base path for century calibration
    hpc_openwq_root : str
        HPC path to OpenWQ source
    hpc_sif_path : str
        HPC path to Apptainer SIF container
    dry_run : bool
        If True, report changes but don't write

    Returns
    -------
    int
        Number of lines modified
    """
    with open(py_path, 'r') as f:
        content = f.read()

    original = content
    count = 0

    # Replace CALIBRATION_LIB path
    content, n = re.subn(
        r'CALIBRATION_LIB\s*=\s*"[^"]+"',
        f'CALIBRATION_LIB = "{hpc_openwq_root}/supporting_scripts/Calibration"',
        content
    )
    count += n

    # Replace HPC_BASE
    content, n = re.subn(
        r'HPC_BASE\s*=\s*"[^"]+"',
        f'HPC_BASE = "{hpc_base}"',
        content
    )
    count += n

    # Replace apptainer_sif_path
    content, n = re.subn(
        r'apptainer_sif_path\s*=\s*"[^"]+"',
        f'apptainer_sif_path = "{hpc_sif_path}"',
        content
    )
    count += n

    # Replace apptainer_bind_path
    content, n = re.subn(
        r'apptainer_bind_path\s*=\s*[^\n]+',
        f'apptainer_bind_path = HPC_BASE + ":/code"',
        content
    )
    count += n

    if content != original:
        if dry_run:
            print(f"  [DRY-RUN] {os.path.basename(py_path)}: {count} replacements would be made")
        else:
            with open(py_path, 'w') as f:
                f.write(content)
            print(f"  [✓] {os.path.basename(py_path)}: {count} replacements")

    return count


# ─────────────────────────────────────────────────────────────────────
# Plain-text config rewriting (mizuroute.control, fileManager.txt, etc.)
# ─────────────────────────────────────────────────────────────────────

def rewrite_text_file_paths(text_path, local_prepared_dir, hpc_base,
                            dry_run=False):
    """Rewrite local absolute paths in plain-text config files.

    Performs a simple string replacement of the local prepared-dir prefix
    with the HPC basins path.  Works for mizuroute.control, SUMMA
    fileManager.txt, and similar files.

    Parameters
    ----------
    text_path : str
        Path to the text file
    local_prepared_dir : str
        Local path to century_basins_prepared
    hpc_base : str
        HPC base path
    dry_run : bool
        Preview only

    Returns
    -------
    int
        Number of replacements made
    """
    with open(text_path, 'r') as f:
        content = f.read()

    original = content
    hpc_basins = os.path.join(hpc_base, "basins")

    count = content.count(local_prepared_dir)
    if count > 0:
        content = content.replace(local_prepared_dir, hpc_basins)

    if content != original:
        if dry_run:
            print(f"  [DRY-RUN] {os.path.basename(text_path)}: "
                  f"{count} path(s) would be rewritten")
        else:
            with open(text_path, 'w') as f:
                f.write(content)
            print(f"  [✓] {os.path.basename(text_path)}: "
                  f"{count} path(s) rewritten")

    return count


# ─────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Rewrite local paths to HPC paths in all generated configs",
    )

    parser.add_argument('--prepared-dir', required=True,
                       help='Path to century_basins_prepared directory')
    parser.add_argument('--hpc-base', required=True,
                       help='HPC base path (e.g. /projects/.../century_calibration)')
    parser.add_argument('--hpc-lulc-dir', required=True,
                       help='HPC path to ESACCI land cover data directory')
    parser.add_argument('--hpc-openwq-root', required=True,
                       help='HPC path to OpenWQ source root (contains supporting_scripts/)')
    parser.add_argument('--hpc-sif-path', default=None,
                       help='HPC path to Apptainer .sif container')
    parser.add_argument('--dry-run', action='store_true',
                       help='Preview changes without writing files')

    args = parser.parse_args()

    prepared_dir = os.path.abspath(args.prepared_dir)

    print("=" * 70)
    print("HPC PATH FIXUP")
    print("=" * 70)
    print(f"\nLocal prepared dir:  {prepared_dir}")
    print(f"HPC base:            {args.hpc_base}")
    print(f"HPC LULC dir:        {args.hpc_lulc_dir}")
    print(f"HPC OpenWQ root:     {args.hpc_openwq_root}")
    if args.hpc_sif_path:
        print(f"HPC SIF path:        {args.hpc_sif_path}")
    if args.dry_run:
        print(f"\n*** DRY RUN MODE ***")
    print()

    # Find all basin directories
    basin_dirs = sorted(glob.glob(os.path.join(prepared_dir, "basin_*")))
    if not basin_dirs:
        print(f"ERROR: No basin directories found in {prepared_dir}")
        sys.exit(1)

    print(f"Found {len(basin_dirs)} basin directories\n")

    total_json = 0
    total_py = 0
    total_text = 0
    total_paths = 0

    for basin_dir in basin_dirs:
        basin_name = os.path.basename(basin_dir)
        basin_id = basin_name.replace("basin_", "")

        # Find all variant directories
        variant_dirs = sorted(glob.glob(os.path.join(basin_dir, "variant_*")))

        if variant_dirs:
            print(f"[{basin_id}]")

        for variant_dir in variant_dirs:
            variant_name = os.path.basename(variant_dir)

            # 1. Fix openWQ_master.json
            master_json = os.path.join(variant_dir, "openWQ_master.json")
            if os.path.exists(master_json):
                n = rewrite_json_paths(
                    master_json, prepared_dir, args.hpc_base,
                    args.hpc_lulc_dir, args.hpc_openwq_root,
                    dry_run=args.dry_run
                )
                total_paths += n
                total_json += 1

            # 2. Fix config_raw.json
            config_raw = os.path.join(variant_dir, "config_raw.json")
            if os.path.exists(config_raw):
                n = rewrite_json_paths(
                    config_raw, prepared_dir, args.hpc_base,
                    args.hpc_lulc_dir, args.hpc_openwq_root,
                    dry_run=args.dry_run
                )
                total_paths += n
                total_json += 1

            # 3. Fix any SS JSON files
            ss_jsons = glob.glob(os.path.join(variant_dir, "openwq_in", "openWQ_SS_*.json"))
            for ss_json in ss_jsons:
                n = rewrite_json_paths(
                    ss_json, prepared_dir, args.hpc_base,
                    args.hpc_lulc_dir, args.hpc_openwq_root,
                    dry_run=args.dry_run
                )
                total_paths += n
                total_json += 1

            # 4. Fix mizuRoute control file
            ctrl_file = os.path.join(variant_dir, "mizuroute_in",
                                     "settings", "mizuroute.control")
            if os.path.exists(ctrl_file):
                n = rewrite_text_file_paths(
                    ctrl_file, prepared_dir, args.hpc_base,
                    dry_run=args.dry_run
                )
                total_paths += n
                total_text += 1

            # 5. Fix fileManager.txt in mizuroute_in/ and summa_in/
            for fm_subpath in ["mizuroute_in/fileManager.txt",
                               "summa_in/fileManager.txt"]:
                fm_file = os.path.join(variant_dir, fm_subpath)
                if os.path.exists(fm_file):
                    n = rewrite_text_file_paths(
                        fm_file, prepared_dir, args.hpc_base,
                        dry_run=args.dry_run
                    )
                    total_paths += n
                    total_text += 1

        # 6. Fix calibration config Python files
        calib_dir = os.path.join(basin_dir, "calibration")
        calib_configs = sorted(glob.glob(os.path.join(calib_dir, "calibration_config_*.py")))
        for calib_py in calib_configs:
            if args.hpc_sif_path:
                n = rewrite_python_paths(
                    calib_py, args.hpc_base, args.hpc_openwq_root,
                    args.hpc_sif_path, dry_run=args.dry_run
                )
            else:
                n = rewrite_python_paths(
                    calib_py, args.hpc_base, args.hpc_openwq_root,
                    f"{args.hpc_base}/openwq.sif",
                    dry_run=args.dry_run
                )
            total_paths += n
            total_py += 1

    # Summary
    print()
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"JSON files processed:    {total_json}")
    print(f"Text files processed:    {total_text}")
    print(f"Python files processed:  {total_py}")
    print(f"Total path rewrites:     {total_paths}")
    if args.dry_run:
        print(f"\n*** DRY RUN — no files were modified ***")
    print()


if __name__ == '__main__':
    main()
