#!/usr/bin/env python3
# Copyright 2026, Diogo Costa, diogo.costa@uevora.pt
# This file is part of OpenWQ model.
#
# This program, openWQ, is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) later version.

"""
01_extract_grqa_all.py
======================

Extract GRQA water quality observations for all 111 century basins.

For each basin tar.gz:
  1. Extract the river network shapefile
  2. Match GRQA stations within 500m buffer
  3. Extract NO3, NH4, PO4, DO, TN, TP observations (1982-2019)
  4. Output observations.csv + feasibility report per basin

Output structure:
  century_basins_prepared/
    basin_CAN_01AM001_meso/
      observations/
        observations.csv
        calibration_stations.shp
        calibration_summary.txt
        calibration_config.json
    basin_CAN_01BJ010_meso/
      ...
    feasibility_report.csv   (summary across all basins)

USAGE:
  python 01_extract_grqa_all.py --basins-dir /path/to/century_basins
                                --output-dir /path/to/century_basins_prepared
                                [--grqa-data /path/to/local/GRQA]
                                [--buffer 500]
                                [--min-obs 10]
                                [--workers 4]
                                [--basin-filter CAN]
"""

import sys
import os
import argparse
import glob
import tarfile
import tempfile
import shutil
import json
import time
from pathlib import Path
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Dict, List, Optional, Tuple

# Add parent directories to path for imports
SCRIPT_DIR = Path(__file__).resolve().parent
CALIBRATION_DIR = SCRIPT_DIR.parent
sys.path.insert(0, str(CALIBRATION_DIR))

try:
    import pandas as pd
    import numpy as np
except ImportError:
    print("ERROR: pandas and numpy required. Install with: pip install pandas numpy")
    sys.exit(1)

try:
    from calibration_lib.observation_data.grqa_extract_stations import (
        GRQACalibrationExtractor,
        SpeciesMapper,
        GRQA_PARAMETERS,
    )
    GRQA_AVAILABLE = True
except ImportError:
    GRQA_AVAILABLE = False
    print("WARNING: GRQA extraction library not available. "
          "Make sure calibration_lib is accessible.")


# =============================================================================
# CONFIGURATION
# =============================================================================

# Species mapping for all three variants:
# We extract a superset of species needed across all variants
SPECIES_MAPPING = {
    # Nitrogen (needed by all 3 variants)
    "NO3": "NO3-N",
    "NH4": "NH4-N",
    "NO2": "NO2-N",
    "TN": "TN",
    "DON": "DON",

    # Phosphorus (needed by Variant B SWAT)
    "PO4": "PO4-P",
    "TP": "TP",

    # Oxygen (needed by Variants A and C for O2 limitation)
    "DO": "DO",

    # Carbon (needed by Variant C thermodynamic for DOC half-saturation)
    "DOC": "DOC",
}

# Time period matching the SUMMA+mizuRoute simulation range
DEFAULT_TIME_PERIOD = (1982, 2019)

# Minimum observations per station
DEFAULT_MIN_OBS = 10

# Buffer distance (m) for matching stations to river network
DEFAULT_BUFFER_M = 500


# =============================================================================
# BASIN EXTRACTION
# =============================================================================

def parse_basin_name(tar_path: str) -> Dict[str, str]:
    """
    Parse basin metadata from tar.gz filename.

    Example: domain_CAN_01AM001_meso.tar.gz
    Returns: {
        'filename': 'domain_CAN_01AM001_meso.tar.gz',
        'basename': 'domain_CAN_01AM001_meso',
        'country': 'CAN',
        'station_id': '01AM001',
        'scale': 'meso',
        'basin_id': 'CAN_01AM001_meso'
    }
    """
    fname = Path(tar_path).name
    basename = fname.replace('.tar.gz', '')
    parts = basename.split('_')

    # Expected format: domain_COUNTRY_STATIONID_SCALE
    if len(parts) >= 4 and parts[0] == 'domain':
        return {
            'filename': fname,
            'basename': basename,
            'country': parts[1],
            'station_id': parts[2],
            'scale': parts[3],
            'basin_id': f"{parts[1]}_{parts[2]}_{parts[3]}",
        }
    else:
        return {
            'filename': fname,
            'basename': basename,
            'country': 'UNKNOWN',
            'station_id': basename,
            'scale': 'unknown',
            'basin_id': basename,
        }


def find_river_network_shapefile(extract_dir: str) -> Optional[str]:
    """
    Find the river network shapefile inside an extracted basin directory.

    Expected location: {basin}/shapefiles/river_network/*_riverNetwork_lumped.shp
    """
    patterns = [
        os.path.join(extract_dir, '**', 'river_network', '*.shp'),
        os.path.join(extract_dir, '**', 'shapefiles', 'river_network', '*.shp'),
    ]

    for pattern in patterns:
        matches = glob.glob(pattern, recursive=True)
        if matches:
            return matches[0]

    return None


def find_topology_nc(extract_dir: str) -> Optional[str]:
    """
    Find the mizuRoute topology.nc file for time range extraction.
    """
    patterns = [
        os.path.join(extract_dir, '**', 'mizuRoute', 'topology.nc'),
        os.path.join(extract_dir, '**', 'settings', 'mizuRoute', 'topology.nc'),
    ]

    for pattern in patterns:
        matches = glob.glob(pattern, recursive=True)
        if matches:
            return matches[0]

    return None


def extract_basin_grqa(
    tar_path: str,
    output_dir: str,
    grqa_local_path: Optional[str] = None,
    grqa_cache_dir: Optional[str] = None,
    buffer_m: float = DEFAULT_BUFFER_M,
    min_obs: int = DEFAULT_MIN_OBS,
    time_period: Tuple[int, int] = DEFAULT_TIME_PERIOD,
) -> Dict:
    """
    Extract GRQA observations for a single basin.

    Parameters
    ----------
    tar_path : str
        Path to basin tar.gz archive
    output_dir : str
        Base output directory
    grqa_local_path : str, optional
        Path to local GRQA data
    grqa_cache_dir : str, optional
        Shared cache directory for GRQA downloads (avoids re-downloading per basin)
    buffer_m : float
        Buffer distance in meters
    min_obs : int
        Minimum observations per station
    time_period : tuple
        (start_year, end_year)

    Returns
    -------
    dict
        Extraction result with basin metadata and observation statistics
    """
    basin_info = parse_basin_name(tar_path)
    basin_id = basin_info['basin_id']
    result = {
        'basin_id': basin_id,
        'country': basin_info['country'],
        'station_id': basin_info['station_id'],
        'scale': basin_info['scale'],
        'status': 'pending',
        'n_stations': 0,
        'n_observations': 0,
        'species_found': [],
        'error': None,
    }

    # Output paths
    basin_output_dir = os.path.join(output_dir, f"basin_{basin_id}")
    obs_output_dir = os.path.join(basin_output_dir, "observations")
    os.makedirs(obs_output_dir, exist_ok=True)

    # Create temporary directory for extraction
    with tempfile.TemporaryDirectory(prefix=f"grqa_{basin_id}_") as tmp_dir:
        try:
            # 1. Extract tar.gz (only shapefiles needed)
            print(f"  [{basin_id}] Extracting river network...")
            with tarfile.open(tar_path, 'r:gz') as tar:
                # Extract only shapefile members to save time/space
                members = [m for m in tar.getmembers()
                          if 'river_network' in m.name or 'topology.nc' in m.name]
                if not members:
                    # Fall back to extracting all if river_network not found
                    members = tar.getmembers()
                tar.extractall(tmp_dir, members=members)

            # 2. Find river network shapefile
            river_shp = find_river_network_shapefile(tmp_dir)
            if river_shp is None:
                result['status'] = 'no_shapefile'
                result['error'] = 'River network shapefile not found'
                print(f"  [{basin_id}] WARNING: No river network shapefile found")
                return result

            # 3. Run GRQA extraction
            print(f"  [{basin_id}] Running GRQA extraction...")
            mapper = SpeciesMapper(mapping=SPECIES_MAPPING)

            extractor = GRQACalibrationExtractor(
                output_dir=obs_output_dir,
                species_mapper=mapper,
                cache_dir=grqa_cache_dir,
                local_data_path=grqa_local_path,
                buffer_distance_m=buffer_m,
                years=time_period,
            )

            stations, observations = extractor.run(
                shapefile_path=river_shp,
                prefix="calibration"
            )

            # 4. Check results
            if stations is None or (hasattr(stations, 'empty') and stations.empty):
                result['status'] = 'no_data'
                result['error'] = 'No GRQA data within buffer'
                print(f"  [{basin_id}] No GRQA data found within {buffer_m}m buffer")
                return result

            # Count results
            n_stations = len(stations) if stations is not None else 0
            n_obs = len(observations) if observations is not None else 0

            # Filter by minimum observations
            if n_obs < min_obs:
                result['status'] = 'insufficient_data'
                result['n_stations'] = n_stations
                result['n_observations'] = n_obs
                result['error'] = f'Only {n_obs} observations (min: {min_obs})'
                print(f"  [{basin_id}] Insufficient data: {n_obs} obs (min: {min_obs})")
                return result

            # Species available
            species_found = list(observations['model_species'].unique()) if 'model_species' in observations.columns else []

            result['status'] = 'success'
            result['n_stations'] = n_stations
            result['n_observations'] = n_obs
            result['species_found'] = species_found
            print(f"  [{basin_id}] Success: {n_stations} stations, {n_obs} observations, "
                  f"species: {species_found}")

        except Exception as e:
            result['status'] = 'error'
            result['error'] = str(e)
            print(f"  [{basin_id}] ERROR: {e}")

    return result


# =============================================================================
# MAIN WORKFLOW
# =============================================================================

def generate_feasibility_report(results: List[Dict], output_dir: str):
    """
    Generate a CSV feasibility report summarizing all basins.
    """
    df = pd.DataFrame(results)

    # Sort by status then basin_id
    status_order = {'success': 0, 'insufficient_data': 1, 'no_data': 2,
                    'no_shapefile': 3, 'error': 4, 'pending': 5}
    df['status_order'] = df['status'].map(status_order).fillna(99)
    df = df.sort_values(['status_order', 'basin_id']).drop(columns='status_order')

    # Save full report
    report_path = os.path.join(output_dir, "feasibility_report.csv")
    df.to_csv(report_path, index=False)

    # Print summary
    print("\n" + "=" * 70)
    print("FEASIBILITY REPORT")
    print("=" * 70)
    print(f"\nTotal basins: {len(df)}")
    print(f"\nStatus breakdown:")
    for status, count in df['status'].value_counts().items():
        print(f"  {status:20s}: {count:3d}")

    success = df[df['status'] == 'success']
    if not success.empty:
        print(f"\nSuccessful basins: {len(success)}")
        print(f"  Total stations:     {success['n_stations'].sum()}")
        print(f"  Total observations: {success['n_observations'].sum()}")
        print(f"  By country:")
        for country, group in success.groupby('country'):
            print(f"    {country}: {len(group)} basins")
        print(f"  By scale:")
        for scale, group in success.groupby('scale'):
            print(f"    {scale}: {len(group)} basins")

        # Species coverage
        all_species = set()
        for sp_list in success['species_found']:
            if isinstance(sp_list, list):
                all_species.update(sp_list)
        print(f"\n  Species found: {sorted(all_species)}")

    print(f"\nReport saved to: {report_path}")
    return df


def main():
    parser = argparse.ArgumentParser(
        description="Extract GRQA observations for all century basins",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage with local GRQA data
  python 01_extract_grqa_all.py \\
      --basins-dir /path/to/century_basins \\
      --output-dir /path/to/century_basins_prepared \\
      --grqa-data /path/to/GRQA

  # Filter to Canadian basins only
  python 01_extract_grqa_all.py \\
      --basins-dir /path/to/century_basins \\
      --output-dir /path/to/century_basins_prepared \\
      --basin-filter CAN

  # Custom buffer and parallel workers
  python 01_extract_grqa_all.py \\
      --basins-dir /path/to/century_basins \\
      --output-dir /path/to/century_basins_prepared \\
      --buffer 1000 --workers 8

  # Single basin test
  python 01_extract_grqa_all.py \\
      --basins-dir /path/to/century_basins \\
      --output-dir /path/to/century_basins_prepared \\
      --single-basin domain_CAN_01DG003_headwater.tar.gz
        """
    )

    parser.add_argument('--basins-dir', required=True,
                       help='Directory containing basin tar.gz archives')
    parser.add_argument('--output-dir', required=True,
                       help='Output directory for prepared basins')
    parser.add_argument('--grqa-data', default=None,
                       help='Path to local GRQA data (or downloads from Zenodo)')
    parser.add_argument('--buffer', type=float, default=DEFAULT_BUFFER_M,
                       help=f'Buffer distance in meters (default: {DEFAULT_BUFFER_M})')
    parser.add_argument('--min-obs', type=int, default=DEFAULT_MIN_OBS,
                       help=f'Minimum observations per basin (default: {DEFAULT_MIN_OBS})')
    parser.add_argument('--start-year', type=int, default=DEFAULT_TIME_PERIOD[0],
                       help=f'Start year for observations (default: {DEFAULT_TIME_PERIOD[0]})')
    parser.add_argument('--end-year', type=int, default=DEFAULT_TIME_PERIOD[1],
                       help=f'End year for observations (default: {DEFAULT_TIME_PERIOD[1]})')
    parser.add_argument('--workers', type=int, default=1,
                       help='Number of parallel workers (default: 1, sequential)')
    parser.add_argument('--basin-filter', default=None,
                       help='Filter basins by country code (e.g., CAN, USA)')
    parser.add_argument('--single-basin', default=None,
                       help='Process a single basin file (for testing)')
    parser.add_argument('--resume', action='store_true',
                       help='Skip basins that already have observations.csv')

    args = parser.parse_args()

    if not GRQA_AVAILABLE:
        print("ERROR: GRQA extraction library not available.")
        print("Make sure you're running from the correct directory and dependencies are installed.")
        sys.exit(1)

    # Validate inputs
    basins_dir = Path(args.basins_dir)
    if not basins_dir.exists():
        print(f"ERROR: Basins directory not found: {basins_dir}")
        sys.exit(1)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    time_period = (args.start_year, args.end_year)

    # Find all basin archives
    if args.single_basin:
        tar_files = [str(basins_dir / args.single_basin)]
        if not os.path.exists(tar_files[0]):
            print(f"ERROR: Basin file not found: {tar_files[0]}")
            sys.exit(1)
    else:
        tar_files = sorted(glob.glob(str(basins_dir / "domain_*.tar.gz")))

    if not tar_files:
        print(f"ERROR: No basin archives found in {basins_dir}")
        sys.exit(1)

    # Apply country filter
    if args.basin_filter:
        tar_files = [f for f in tar_files
                    if f"_{args.basin_filter}_" in Path(f).name]
        print(f"Filtered to {len(tar_files)} basins matching '{args.basin_filter}'")

    # Resume: skip basins with existing observations
    if args.resume:
        remaining = []
        for tf in tar_files:
            basin_info = parse_basin_name(tf)
            obs_file = output_dir / f"basin_{basin_info['basin_id']}" / "observations" / "calibration_observations.csv"
            if not obs_file.exists():
                remaining.append(tf)
            else:
                print(f"  Skipping {basin_info['basin_id']} (already extracted)")
        tar_files = remaining
        print(f"Resuming: {len(remaining)} basins remaining")

    print("\n" + "=" * 70)
    print("GRQA EXTRACTION FOR CENTURY BASINS")
    print("=" * 70)
    print(f"\nBasins directory:  {basins_dir}")
    print(f"Output directory:  {output_dir}")
    print(f"GRQA data source:  {args.grqa_data or 'Zenodo (download)'}")
    print(f"Buffer distance:   {args.buffer}m")
    print(f"Min observations:  {args.min_obs}")
    print(f"Time period:       {time_period[0]}-{time_period[1]}")
    print(f"Total basins:      {len(tar_files)}")
    print(f"Workers:           {args.workers}")
    print(f"\nSpecies mapping:")
    for grqa_name, model_name in SPECIES_MAPPING.items():
        print(f"  {grqa_name:6s} -> {model_name}")

    print(f"\n{'=' * 70}")
    start_time = time.time()

    # Shared GRQA cache directory (avoids re-downloading 1.2GB per basin)
    shared_cache_dir = str(output_dir / "grqa_shared_cache")
    os.makedirs(shared_cache_dir, exist_ok=True)

    # If a basin already has the downloaded zip, copy to shared cache
    existing_zip = output_dir / "basin_CAN_01DG003_headwater" / "observations" / "grqa_cache" / "GRQA_data_v1.4.zip"
    shared_zip = Path(shared_cache_dir) / "GRQA_data_v1.4.zip"
    if existing_zip.exists() and not shared_zip.exists():
        print(f"Copying existing GRQA zip to shared cache...")
        shutil.copy2(str(existing_zip), str(shared_zip))

    # Process basins
    results = []

    if args.workers <= 1:
        # Sequential processing
        for i, tar_path in enumerate(tar_files):
            basin_info = parse_basin_name(tar_path)
            print(f"\n[{i+1}/{len(tar_files)}] Processing {basin_info['basin_id']}...")
            result = extract_basin_grqa(
                tar_path=tar_path,
                output_dir=str(output_dir),
                grqa_local_path=args.grqa_data,
                grqa_cache_dir=shared_cache_dir,
                buffer_m=args.buffer,
                min_obs=args.min_obs,
                time_period=time_period,
            )
            results.append(result)
    else:
        # Parallel processing
        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            futures = {}
            for tar_path in tar_files:
                future = executor.submit(
                    extract_basin_grqa,
                    tar_path=tar_path,
                    output_dir=str(output_dir),
                    grqa_local_path=args.grqa_data,
                    grqa_cache_dir=shared_cache_dir,
                    buffer_m=args.buffer,
                    min_obs=args.min_obs,
                    time_period=time_period,
                )
                futures[future] = tar_path

            for i, future in enumerate(as_completed(futures)):
                tar_path = futures[future]
                basin_info = parse_basin_name(tar_path)
                try:
                    result = future.result()
                    results.append(result)
                    print(f"[{i+1}/{len(tar_files)}] Completed {basin_info['basin_id']}: "
                          f"{result['status']}")
                except Exception as e:
                    print(f"[{i+1}/{len(tar_files)}] FAILED {basin_info['basin_id']}: {e}")
                    results.append({
                        'basin_id': basin_info['basin_id'],
                        'country': basin_info['country'],
                        'station_id': basin_info['station_id'],
                        'scale': basin_info['scale'],
                        'status': 'error',
                        'n_stations': 0,
                        'n_observations': 0,
                        'species_found': [],
                        'error': str(e),
                    })

    elapsed = time.time() - start_time

    # Generate feasibility report
    report_df = generate_feasibility_report(results, str(output_dir))

    print(f"\nTotal time: {elapsed/60:.1f} minutes")
    print(f"Output directory: {output_dir}")


if __name__ == "__main__":
    main()
