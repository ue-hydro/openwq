#!/usr/bin/env python3
# Copyright 2026, Diogo Costa, diogo.costa@uevora.pt
# This file is part of OpenWQ model.
#
# This program, openWQ, is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) later version.

"""
08_extract_mizuroute_summa_inputs.py
====================================

Extract mizuRoute and SUMMA input files from domain archives (domain_*.tar.gz)
into the century_basins_prepared directory structure, and set the simulation
period to 1993-2019.

For each basin archive this script:
  1. Extracts ERA5 forcing data to basin_{ID}/forcing/  (shared across variants)
  2. For each variant (A, B, C):
     - Creates mizuroute_in/  with control file, topology, remapping, params
     - Creates summa_in/      with all SUMMA settings (NC files, text configs, TBLs)
     - Generates updated mizuroute.control  (new sim period + local paths)
     - Generates updated SUMMA fileManager.txt  (new sim period + local paths)
     - Places a copy of fileManager.txt in mizuroute_in/ (entry point for
       the coupled executable, as expected by calibration_config_*.py)
  3. Creates empty output directories (simulations/run_1/SUMMA/ and .../mizuRoute/)

Local absolute paths are written initially; run 07_fixup_hpc_paths.py afterwards
to convert them to HPC paths.

USAGE:
  python 08_extract_mizuroute_summa_inputs.py \\
      --prepared-dir /path/to/century_basins_prepared \\
      --archives-dir /path/to/century_basins/century_basins \\
      [--sim-start "1993-01-01 01:00"] \\
      [--sim-end "2019-12-31 23:00"] \\
      [--force]   # Overwrite existing files
"""

import tarfile
import os
import sys
import glob
import re
import argparse
from pathlib import Path


# ─────────────────────────────────────────────────────────────────────
# Constants
# ─────────────────────────────────────────────────────────────────────

DEFAULT_SIM_START = "1993-01-01 01:00"
DEFAULT_SIM_END   = "2019-12-31 23:00"

VARIANT_DIRS = [
    "variant_A_nitrogen_full",
    "variant_B_swat",
    "variant_C_thermodynamic",
    "variant_D_phreeqc",
]

# Files to extract from archive → target locations
MIZUROUTE_ANCILLARY = ["topology.nc", "routing_remap.nc"]
MIZUROUTE_INPUT     = ["param.nml.default"]

SUMMA_NC_FILES = [
    "attributes.nc", "coldState.nc", "trialParams.nc",
]
SUMMA_TEXT_FILES = [
    "modelDecisions.txt", "outputControl.txt",
    "localParamInfo.txt", "basinParamInfo.txt",
]
SUMMA_TABLE_FILES = [
    "TBL_VEGPARM.TBL", "TBL_SOILPARM.TBL",
    "TBL_GENPARM.TBL", "TBL_MPTABLE.TBL",
]

ALL_SUMMA_EXTRACT = SUMMA_NC_FILES + SUMMA_TEXT_FILES + SUMMA_TABLE_FILES


# ─────────────────────────────────────────────────────────────────────
# Tar helpers
# ─────────────────────────────────────────────────────────────────────

def extract_member_bytes(tar, member_path):
    """Read a tar member into memory (bytes). Returns None on failure."""
    try:
        member = tar.getmember(member_path)
        f = tar.extractfile(member)
        if f:
            return f.read()
    except (KeyError, AttributeError):
        pass
    return None


def extract_member_to_disk(tar, member_path, dest_path):
    """Write a tar member to disk. Returns True on success."""
    data = extract_member_bytes(tar, member_path)
    if data is not None:
        os.makedirs(os.path.dirname(dest_path), exist_ok=True)
        with open(dest_path, 'wb') as f:
            f.write(data)
        return True
    return False


# ─────────────────────────────────────────────────────────────────────
# Control-file rewriters
# ─────────────────────────────────────────────────────────────────────

def rewrite_mizuroute_control(raw_bytes, prepared_dir, basin_id,
                               variant_dir, sim_start, sim_end):
    """Return a rewritten mizuroute.control string.

    Updates <sim_start>, <sim_end>, <ancil_dir>, <input_dir>, <output_dir>
    while preserving every other line verbatim.
    """
    basin_path   = os.path.join(prepared_dir, f"basin_{basin_id}")
    variant_path = os.path.join(basin_path, variant_dir)

    lines = raw_bytes.decode('utf-8').split('\n')
    out = []

    for line in lines:
        stripped = line.strip()

        # Keep blanks & pure comments
        if not stripped or stripped.startswith('!'):
            out.append(line)
            continue

        # Try to match <key>  value  ! comment
        m = re.match(r'(\s*<(\w+)>\s+)(.*?)(\s+!.*)?$', line)
        if not m:
            out.append(line)
            continue

        prefix  = m.group(1)
        key     = m.group(2)
        value   = m.group(3).strip()
        comment = m.group(4) or ''

        if key == 'ancil_dir':
            value = os.path.join(variant_path, "mizuroute_in", "ancillary_data") + "/"
        elif key == 'input_dir':
            value = os.path.join(variant_path, "simulations", "run_1", "SUMMA") + "/"
        elif key == 'output_dir':
            value = os.path.join(variant_path, "simulations", "run_1", "mizuRoute") + "/"
        elif key == 'sim_start':
            value = sim_start
        elif key == 'sim_end':
            value = sim_end

        out.append(f"{prefix}{value}    {comment}".rstrip())

    return '\n'.join(out)


def rewrite_summa_filemanager(raw_bytes, prepared_dir, basin_id,
                               variant_dir, sim_start, sim_end):
    """Return a rewritten SUMMA fileManager.txt string.

    Updates simStartTime, simEndTime, settingsPath, forcingPath, outputPath.
    """
    basin_path   = os.path.join(prepared_dir, f"basin_{basin_id}")
    variant_path = os.path.join(basin_path, variant_dir)

    lines = raw_bytes.decode('utf-8').split('\n')
    out = []

    for line in lines:
        stripped = line.strip()
        if not stripped:
            out.append(line)
            continue

        parts = stripped.split(None, 1)
        if len(parts) < 2:
            out.append(line)
            continue

        key = parts[0]

        if key == 'simStartTime':
            out.append(f"simStartTime         '{sim_start}'")
        elif key == 'simEndTime':
            out.append(f"simEndTime           '{sim_end}'")
        elif key == 'settingsPath':
            out.append(f"settingsPath         '{variant_path}/summa_in/'")
        elif key == 'forcingPath':
            out.append(f"forcingPath          '{basin_path}/forcing/'")
        elif key == 'outputPath':
            out.append(f"outputPath           '{variant_path}/simulations/run_1/SUMMA/'")
        else:
            out.append(line)

    return '\n'.join(out)


def rewrite_forcing_file_list(raw_bytes):
    """Return a rewritten forcingFileList.txt string.

    Keeps only the basename(s) — forcingPath in the fileManager already
    points to the correct directory.
    """
    lines = raw_bytes.decode('utf-8').strip().split('\n')
    out = []
    for line in lines:
        fname = line.strip().strip("'\"")
        if fname and not fname.startswith('!'):
            out.append(os.path.basename(fname))
    return '\n'.join(out) + '\n'


# ─────────────────────────────────────────────────────────────────────
# Per-basin processor
# ─────────────────────────────────────────────────────────────────────

def process_basin(tar, archive_prefix, basin_id, prepared_dir,
                  sim_start, sim_end, force=False):
    """Extract and set up mizuroute_in/, summa_in/, forcing/ for one basin.

    Returns a stats dict or None if the basin directory is missing.
    """
    basin_dir = os.path.join(prepared_dir, f"basin_{basin_id}")

    if not os.path.isdir(basin_dir):
        return None

    # Discover which variants exist
    variants = [v for v in VARIANT_DIRS
                if os.path.isdir(os.path.join(basin_dir, v))]
    if not variants:
        print(f"  WARNING: no variant dirs in {basin_dir}")
        return None

    stats = {"forcing": 0, "mizuroute": 0, "summa": 0,
             "control": 0, "skipped": 0, "missing": 0}

    # ── 1. Forcing (shared across variants) ──────────────────────────
    forcing_dir = os.path.join(basin_dir, "forcing")
    os.makedirs(forcing_dir, exist_ok=True)

    # Find .nc files under forcing/SUMMA_input/ in the archive
    forcing_prefix = f"{archive_prefix}/forcing/SUMMA_input/"
    forcing_members = [n for n in tar.getnames()
                       if n.startswith(forcing_prefix) and n.endswith('.nc')]

    for fm in forcing_members:
        fname = os.path.basename(fm)
        dest  = os.path.join(forcing_dir, fname)
        if os.path.exists(dest) and not force:
            stats["skipped"] += 1
        elif extract_member_to_disk(tar, fm, dest):
            sz = os.path.getsize(dest) / (1024 * 1024)
            stats["forcing"] += 1

    # ── 2. Read original control files from the archive ──────────────
    ctrl_bytes = extract_member_bytes(
        tar, f"{archive_prefix}/settings/mizuRoute/mizuroute.control")
    fm_bytes   = extract_member_bytes(
        tar, f"{archive_prefix}/settings/SUMMA/fileManager.txt")
    fl_bytes   = extract_member_bytes(
        tar, f"{archive_prefix}/settings/SUMMA/forcingFileList.txt")

    if ctrl_bytes is None:
        print(f"  WARNING: no mizuroute.control in archive")
        stats["missing"] += 1
    if fm_bytes is None:
        print(f"  WARNING: no fileManager.txt in archive")
        stats["missing"] += 1

    # ── 3. Per-variant setup ─────────────────────────────────────────
    for variant in variants:
        variant_path = os.path.join(basin_dir, variant)
        miz_in       = os.path.join(variant_path, "mizuroute_in")
        summa_in     = os.path.join(variant_path, "summa_in")

        # --- mizuRoute ancillary data (topology, remapping) ---
        for fname in MIZUROUTE_ANCILLARY:
            src  = f"{archive_prefix}/settings/mizuRoute/{fname}"
            dest = os.path.join(miz_in, "ancillary_data", fname)
            if os.path.exists(dest) and not force:
                stats["skipped"] += 1
            elif extract_member_to_disk(tar, src, dest):
                stats["mizuroute"] += 1
            else:
                stats["missing"] += 1

        # --- mizuRoute input (param namelist) ---
        for fname in MIZUROUTE_INPUT:
            src  = f"{archive_prefix}/settings/mizuRoute/{fname}"
            dest = os.path.join(miz_in, "input", fname)
            if os.path.exists(dest) and not force:
                stats["skipped"] += 1
            elif extract_member_to_disk(tar, src, dest):
                stats["mizuroute"] += 1
            else:
                stats["missing"] += 1

        # --- Rewritten mizuroute.control ---
        ctrl_dest = os.path.join(miz_in, "settings", "mizuroute.control")
        if ctrl_bytes and (not os.path.exists(ctrl_dest) or force):
            new_ctrl = rewrite_mizuroute_control(
                ctrl_bytes, prepared_dir, basin_id, variant,
                sim_start, sim_end)
            os.makedirs(os.path.dirname(ctrl_dest), exist_ok=True)
            with open(ctrl_dest, 'w') as f:
                f.write(new_ctrl)
            stats["control"] += 1

        # --- SUMMA binary & text files ---
        for fname in ALL_SUMMA_EXTRACT:
            src  = f"{archive_prefix}/settings/SUMMA/{fname}"
            dest = os.path.join(summa_in, fname)
            if os.path.exists(dest) and not force:
                stats["skipped"] += 1
            elif extract_member_to_disk(tar, src, dest):
                stats["summa"] += 1
            else:
                stats["missing"] += 1

        # --- Rewritten forcingFileList.txt ---
        fl_dest = os.path.join(summa_in, "forcingFileList.txt")
        if fl_bytes and (not os.path.exists(fl_dest) or force):
            new_fl = rewrite_forcing_file_list(fl_bytes)
            os.makedirs(summa_in, exist_ok=True)
            with open(fl_dest, 'w') as f:
                f.write(new_fl)
            stats["summa"] += 1

        # --- Rewritten SUMMA fileManager.txt in summa_in/ ---
        summa_fm_dest = os.path.join(summa_in, "fileManager.txt")
        if fm_bytes and (not os.path.exists(summa_fm_dest) or force):
            new_fm = rewrite_summa_filemanager(
                fm_bytes, prepared_dir, basin_id, variant,
                sim_start, sim_end)
            os.makedirs(summa_in, exist_ok=True)
            with open(summa_fm_dest, 'w') as f:
                f.write(new_fm)
            stats["control"] += 1

        # --- Copy of fileManager.txt in mizuroute_in/ ---
        # (This is the entry point for the coupled executable;
        #  calibration_config_*.py sets file_manager_path here.)
        miz_fm_dest = os.path.join(miz_in, "fileManager.txt")
        if fm_bytes and (not os.path.exists(miz_fm_dest) or force):
            new_fm = rewrite_summa_filemanager(
                fm_bytes, prepared_dir, basin_id, variant,
                sim_start, sim_end)
            os.makedirs(miz_in, exist_ok=True)
            with open(miz_fm_dest, 'w') as f:
                f.write(new_fm)
            stats["control"] += 1

        # --- Create empty output dirs ---
        os.makedirs(os.path.join(variant_path, "simulations", "run_1", "SUMMA"),
                    exist_ok=True)
        os.makedirs(os.path.join(variant_path, "simulations", "run_1", "mizuRoute"),
                    exist_ok=True)

    return stats


# ─────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Extract mizuRoute + SUMMA inputs from domain archives "
                    "into century_basins_prepared/")

    parser.add_argument('--prepared-dir', required=True,
                        help='Path to century_basins_prepared/')
    parser.add_argument('--archives-dir', required=True,
                        help='Path to directory containing domain_*.tar.gz')
    parser.add_argument('--sim-start', default=DEFAULT_SIM_START,
                        help=f'Simulation start (default: {DEFAULT_SIM_START})')
    parser.add_argument('--sim-end', default=DEFAULT_SIM_END,
                        help=f'Simulation end (default: {DEFAULT_SIM_END})')
    parser.add_argument('--force', action='store_true',
                        help='Overwrite existing files')

    args = parser.parse_args()

    prepared_dir = os.path.abspath(args.prepared_dir)
    archives_dir = os.path.abspath(args.archives_dir)

    # Discover archives
    archives = sorted(glob.glob(os.path.join(archives_dir, "domain_*.tar.gz")))
    if not archives:
        print(f"ERROR: No domain_*.tar.gz found in {archives_dir}")
        sys.exit(1)

    print("=" * 70)
    print("EXTRACT MIZUROUTE + SUMMA INPUTS FROM DOMAIN ARCHIVES")
    print("=" * 70)
    print(f"\nArchives dir:   {archives_dir}")
    print(f"Prepared dir:   {prepared_dir}")
    print(f"Sim period:     {args.sim_start}  →  {args.sim_end}")
    print(f"Archives found: {len(archives)}")
    print(f"Force:          {args.force}")
    print()

    # Totals
    total = {"basins": 0, "forcing": 0, "mizuroute": 0, "summa": 0,
             "control": 0, "skipped": 0, "missing": 0, "no_match": 0}

    for i, archive_path in enumerate(archives, 1):
        archive_name = os.path.basename(archive_path)

        # Extract basin ID: domain_{BASIN_ID}.tar.gz → BASIN_ID
        basin_id = archive_name.replace("domain_", "").replace(".tar.gz", "")
        archive_prefix = f"domain_{basin_id}"

        print(f"[{i:3d}/{len(archives)}] {basin_id} ... ", end="", flush=True)

        # Check prepared dir exists
        basin_dir = os.path.join(prepared_dir, f"basin_{basin_id}")
        if not os.path.isdir(basin_dir):
            print("SKIP (no prepared dir)")
            total["no_match"] += 1
            continue

        try:
            with tarfile.open(archive_path, 'r:gz') as tar:
                stats = process_basin(
                    tar, archive_prefix, basin_id, prepared_dir,
                    args.sim_start, args.sim_end, args.force)
        except Exception as e:
            print(f"ERROR: {e}")
            continue

        if stats is None:
            print("SKIP (no variants)")
            total["no_match"] += 1
            continue

        total["basins"]    += 1
        total["forcing"]   += stats["forcing"]
        total["mizuroute"] += stats["mizuroute"]
        total["summa"]     += stats["summa"]
        total["control"]   += stats["control"]
        total["skipped"]   += stats["skipped"]
        total["missing"]   += stats["missing"]

        parts = []
        if stats["forcing"]:
            parts.append(f"forcing={stats['forcing']}")
        if stats["mizuroute"]:
            parts.append(f"miz={stats['mizuroute']}")
        if stats["summa"]:
            parts.append(f"summa={stats['summa']}")
        if stats["control"]:
            parts.append(f"ctrl={stats['control']}")
        if stats["skipped"]:
            parts.append(f"skip={stats['skipped']}")
        if stats["missing"]:
            parts.append(f"miss={stats['missing']}")

        print(", ".join(parts) if parts else "nothing to do")

    # Summary
    print()
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Basins processed:    {total['basins']}")
    print(f"Forcing files:       {total['forcing']}")
    print(f"mizuRoute files:     {total['mizuroute']}")
    print(f"SUMMA files:         {total['summa']}")
    print(f"Control files:       {total['control']}")
    print(f"Skipped (exist):     {total['skipped']}")
    print(f"Missing in archive:  {total['missing']}")
    if total["no_match"]:
        print(f"No prepared dir:     {total['no_match']}")
    print()


if __name__ == '__main__':
    main()
