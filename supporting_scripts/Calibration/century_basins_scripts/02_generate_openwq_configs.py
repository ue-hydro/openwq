#!/usr/bin/env python3
# Copyright 2026, Diogo Costa, diogo.costa@uevora.pt
# This file is part of OpenWQ model.
#
# This program, openWQ, is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) later version.

"""
02_generate_openwq_configs.py
=============================

Generate OpenWQ configuration files for all century basins, creating 4 variant
config sets per basin (one per N-cycle formulation).

For each basin:
  variant_A_nitrogen_full/     -> nitrogen_full.json BGC template (NATIVE_BGC_FLEX)
  variant_B_swat/              -> SWAT_full_nutrients.json BGC template (NATIVE_BGC_FLEX)
  variant_C_thermodynamic/     -> Combined nitrification + denitrification + anammox (NATIVE_BGC_FLEX)
  variant_D_phreeqc/           -> nitrogen_speciation.pqi PHREEQC template

Output structure:
  century_basins_prepared/
    basin_CAN_01AM001_meso/
      variant_A_nitrogen_full/
        openwq_in/
          openWQ_MODULE_NATIVE_BGC_FLEX.json
          openWQ_SS_copernicus.json (or CSV-based)
          openWQ_MODULE_OPENWQ_NATIVE_TD_ADVDISP.json
        openWQ_master.json
      variant_B_swat/
        openwq_in/
          ... (SWAT-based BGC)
        openWQ_master.json
      variant_C_thermodynamic/
        openwq_in/
          ... (thermodynamic BGC)
        openWQ_master.json
      variant_D_phreeqc/
        openwq_in/
          openWQ_MODULE_PHREEQC.json
          nitrogen_speciation.pqi
          phreeqc.dat
        openWQ_master.json
      observations/
        ... (from step 01)

USAGE:
  python 02_generate_openwq_configs.py \\
      --basins-dir /path/to/century_basins \\
      --prepared-dir /path/to/century_basins_prepared \\
      [--openwq-dir /path/to/openwq/supporting_scripts/Model_Config] \\
      [--variants A,B,C,D] \\
      [--ss-method copernicus] \\
      [--basin-filter CAN]
"""

import sys
import os
import argparse
import glob
import json
import tarfile
import tempfile
import shutil
import copy
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# Add parent paths for imports
SCRIPT_DIR = Path(__file__).resolve().parent
CALIBRATION_DIR = SCRIPT_DIR.parent
OPENWQ_DIR = CALIBRATION_DIR.parent
MODEL_CONFIG_DIR = OPENWQ_DIR / "Model_Config"

sys.path.insert(0, str(MODEL_CONFIG_DIR))
sys.path.insert(0, str(MODEL_CONFIG_DIR / "config_support_lib"))
sys.path.insert(0, str(CALIBRATION_DIR))  # For calibration_lib imports

try:
    import Gen_Input_Driver as gJSON_lib
    from Gen_Input_Driver import uniform_param
    CONFIG_LIB_AVAILABLE = True
except ImportError as e:
    CONFIG_LIB_AVAILABLE = False
    print(f"WARNING: Config generation library not available: {e}")


# =============================================================================
# BGC TEMPLATE PATHS (relative to config_support_lib)
# =============================================================================

CONFIG_SUPPORT_LIB = MODEL_CONFIG_DIR / "config_support_lib"
BGC_TEMPLATES = CONFIG_SUPPORT_LIB / "BGC_templates" / "NATIVE_BGC_FLEX"
PHREEQC_TEMPLATES = CONFIG_SUPPORT_LIB / "BGC_templates" / "PHREEQC"
# phreeqc.dat thermodynamic database (used by all PHREEQC variants)
# parents[7] navigates from century_basins_scripts → 6_mizuroute_cslm_openwq (repo root)
PHREEQC_DATABASE = Path(__file__).resolve().parents[7] / "test_bgc_templates" / "phreeqc_01_nitrogen" / "openwq_in" / "phreeqc.dat"

VARIANT_TEMPLATES = {
    'A': {
        'name': 'nitrogen_full',
        'label': 'Variant A: Full N Cycle',
        'dir_name': 'variant_A_nitrogen_full',
        'template': str(BGC_TEMPLATES / "individual_chemicals" / "nitrogen_full.json"),
        'species': ["NH4-N", "NO3-N", "NO2-N", "DON", "PON", "N2O", "DO"],
    },
    'B': {
        'name': 'swat_full_nutrients',
        'label': 'Variant B: SWAT Full Nutrients',
        'dir_name': 'variant_B_swat',
        'template': str(BGC_TEMPLATES / "popular_models" / "SWAT_full_nutrients.json"),
        'species': [
            "NO3-N", "NH4-N", "ORG_N_active", "ORG_N_stable", "ORG_N_fresh",
            "PO4_P_sol", "PO4_P_active", "PO4_P_stable", "ORG_P_humus", "ORG_P_fresh",
            "ORG_C_active", "ORG_C_stable", "ORG_C_fresh",
            "Sediment", "DO",
        ],
    },
    'C': {
        'name': 'thermodynamic',
        'label': 'Variant C: Thermodynamic N Cycle',
        'dir_name': 'variant_C_thermodynamic',
        'templates': [
            str(BGC_TEMPLATES / "thermodynamic" / "templates" / "nitrification.json"),
            str(BGC_TEMPLATES / "thermodynamic" / "templates" / "denitrification_thermodynamic.json"),
            str(BGC_TEMPLATES / "thermodynamic" / "templates" / "anammox.json"),
        ],
        'species': [
            "NH4-N", "NO2-N", "NO3-N", "N2", "DO", "DOC",
            "AOB_biomass", "NOB_biomass", "DENIT_biomass",
        ],
    },
    'D': {
        'name': 'phreeqc_nitrogen',
        'label': 'Variant D: PHREEQC N Cycle',
        'dir_name': 'variant_D_phreeqc',
        'bgc_module': 'PHREEQC',
        'pqi_template': str(PHREEQC_TEMPLATES / "nitrogen_speciation.pqi"),
        'database': str(PHREEQC_DATABASE),
        'mobile_species': ["N", "O", "H", "C"],
        'species': ["NO3-N", "NH4-N", "NO2-N", "N2", "DO"],
    },
}


# =============================================================================
# THERMODYNAMIC TEMPLATE MERGING
# =============================================================================

def merge_thermodynamic_templates(template_paths: List[str]) -> dict:
    """
    Merge multiple thermodynamic BGC templates into a single config.

    Each template has its own CYCLING_FRAMEWORKS. We merge them by:
    1. Union of all CHEMICAL_SPECIES
    2. Concatenation of all CYCLING_FRAMEWORKS (renumbering reactions)
    3. Merge compartment assignments

    Parameters
    ----------
    template_paths : list
        Paths to individual template JSON files

    Returns
    -------
    dict
        Merged BGC configuration
    """
    merged = {
        "CHEMICAL_SPECIES": {},
        "CYCLING_FRAMEWORKS": {},
    }

    all_species_list = {}      # Numeric-keyed species dict
    all_mobile_species = []    # BGC_GENERAL_MOBILE_SPECIES list
    all_species_extra = {}     # Any other CHEMICAL_SPECIES sub-keys
    all_frameworks = {}

    for tpath in template_paths:
        with open(tpath, 'r') as f:
            template = json.load(f)

        # Merge species (union of all species across templates)
        if "CHEMICAL_SPECIES" in template:
            cs = template["CHEMICAL_SPECIES"]
            # Handle "list" dict (numeric keys -> species names)
            if "list" in cs:
                for _k, species_name in cs["list"].items():
                    if species_name not in all_species_list.values():
                        next_key = str(len(all_species_list) + 1)
                        all_species_list[next_key] = species_name
            # Handle BGC_GENERAL_MOBILE_SPECIES list
            if "BGC_GENERAL_MOBILE_SPECIES" in cs:
                for s in cs["BGC_GENERAL_MOBILE_SPECIES"]:
                    if s not in all_mobile_species:
                        all_mobile_species.append(s)
            # Handle any other sub-keys (e.g., named species dicts)
            for name, spec in cs.items():
                if name not in ("list", "BGC_GENERAL_MOBILE_SPECIES"):
                    if name not in all_species_extra:
                        all_species_extra[name] = spec

        # Merge cycling frameworks (each template has its own named framework)
        if "CYCLING_FRAMEWORKS" in template:
            for fw_name, fw_data in template["CYCLING_FRAMEWORKS"].items():
                if fw_name not in all_frameworks:
                    all_frameworks[fw_name] = copy.deepcopy(fw_data)
                else:
                    # Same framework name from different templates - merge reactions
                    existing = all_frameworks[fw_name]
                    numeric_keys = [int(k) for k in existing.keys() if k.isdigit()]
                    max_key = max(numeric_keys) if numeric_keys else 0
                    for rxn_key, rxn_data in fw_data.items():
                        if rxn_key.isdigit():
                            new_key = str(max_key + int(rxn_key))
                            existing[new_key] = rxn_data
                        elif rxn_key not in existing:
                            existing[rxn_key] = rxn_data

    # Build merged species
    merged_species = {}
    if all_species_list:
        merged_species["list"] = all_species_list
    if all_mobile_species:
        merged_species["BGC_GENERAL_MOBILE_SPECIES"] = all_mobile_species
    merged_species.update(all_species_extra)

    merged["CHEMICAL_SPECIES"] = merged_species
    merged["CYCLING_FRAMEWORKS"] = all_frameworks

    return merged


# =============================================================================
# CONFIG GENERATION FOR A SINGLE BASIN + VARIANT
# =============================================================================

def generate_variant_config(
    basin_id: str,
    variant_key: str,
    variant_info: dict,
    output_dir: str,
    basin_extract_dir: str,
    ss_method: str = "using_copernicus_lulc_with_static_coeff",
    copernicus_lc_dir: Optional[str] = None,
) -> str:
    """
    Generate OpenWQ config files for one basin + one variant.

    Parameters
    ----------
    basin_id : str
        Basin identifier (e.g., 'CAN_01AM001_meso')
    variant_key : str
        Variant letter ('A', 'B', 'C', or 'D')
    variant_info : dict
        Variant configuration from VARIANT_TEMPLATES
    output_dir : str
        Output directory for this variant's configs
    basin_extract_dir : str
        Temporary directory with extracted basin files
    ss_method : str
        Source/sink method
    copernicus_lc_dir : str, optional
        Path to Copernicus land cover data

    Returns
    -------
    str
        Path to generated config directory
    """
    variant_dir = os.path.join(output_dir, variant_info['dir_name'])
    os.makedirs(variant_dir, exist_ok=True)

    # Handle BGC template
    bgc_json_path = None
    if variant_key == 'C':
        # Merge thermodynamic templates
        bgc_config = merge_thermodynamic_templates(variant_info['templates'])
        bgc_json_path = os.path.join(variant_dir, "merged_thermodynamic_bgc.json")
        with open(bgc_json_path, 'w') as f:
            json.dump(bgc_config, f, indent=2)
    elif variant_key == 'D':
        # PHREEQC variant: no BGC JSON template, uses .pqi file instead
        pass
    else:
        bgc_json_path = variant_info['template']

    # Find basin HRU shapefile for Copernicus SS method
    # Prefer permanent .gpkg from shapefiles/ dir over temp-extracted shapefiles
    basin_hru_gpkg = os.path.join(output_dir, "shapefiles", f"{basin_id}_basinHru.gpkg")
    if os.path.exists(basin_hru_gpkg):
        river_shp = basin_hru_gpkg
    else:
        # Fallback: search in extracted archive
        river_shp = None
        for pattern in ['**/river_basins/*.shp', '**/river_network/*.shp']:
            matches = glob.glob(os.path.join(basin_extract_dir, pattern), recursive=True)
            if matches:
                river_shp = matches[0]
                break

    # Find topology.nc for reach IDs
    topology_nc = None
    for pattern in ['**/mizuRoute/topology.nc']:
        matches = glob.glob(os.path.join(basin_extract_dir, pattern), recursive=True)
        if matches:
            topology_nc = matches[0]
            break

    # Determine species list and output format based on variant
    species = variant_info['species']

    # Build Gen_Input_Driver config
    config = {
        'project_name': f"Century_{basin_id}_{variant_info['name']}",
        'geographical_location': basin_id.split('_')[0],  # CAN or USA
        'authors': "OpenWQ Century Calibration",
        'date': "2026",
        'comment': f"Auto-generated for {variant_info['label']} - Basin {basin_id}",
        'hostmodel': "mizuroute",

        'dir2save_input_files': variant_dir,
        'running_on_docker': False,

        # Computational
        'solver': "FORWARD_EULER",
        'run_mode_debug': False,
        'use_num_threads': 4,

        # Initial conditions
        'ic_all_value': 0.1,
        'ic_all_units': "mg/l",

        # BGC module (set below based on variant)
        'bgc_module_name': variant_info.get('bgc_module', "NATIVE_BGC_FLEX"),
        'path2selected_NATIVE_BGC_FLEX_framework': bgc_json_path or "",

        # Transport
        'td_module_name': "OPENWQ_NATIVE_TD_ADVDISP",
        'td_module_dispersion_xyz': [0.5, 0.5, 0.5],
        'td_module_characteristic_length_m': 100.0,

        # Lateral exchange (not used for mizuroute-only runs)
        'le_module_name': "NONE",
        'le_module_config': [],

        # Sediment transport (enable for SWAT variant only)
        'ts_module_name': "NONE",

        # Sorption (not used in initial calibration)
        'si_module_name': "NONE",

        # Source/sink
        'ss_method': ss_method,
        'ss_metadata_source': f"Auto-generated for {basin_id}",
        'ss_metadata_comment': f"{variant_info['label']}",
        'ss_use_cellid_mapping': True,

        # External water fluxes
        'ewf_method': "none",

        # Output
        'output_format': "HDF5",
        'chemical_species': species,
        'units': "MG/L",
        'no_water_conc_flag': -9999,
        'export_sediment': (variant_key == 'B'),
        'timestep': [1, "day"],
        'compartments_and_cells': {
            "RIVER_NETWORK_REACHES": {
                "1": ["all", "all", "all"]
            }
        },
    }

    # Add PHREEQC-specific parameters for Variant D
    if variant_key == 'D':
        config['phreeqc_input_filepath'] = variant_info['pqi_template']
        config['phreeqc_database_filepath'] = variant_info['database']
        config['phreeqc_mobile_species'] = variant_info['mobile_species']
        config['phreeqc_component_h2o'] = True
        config['phreeqc_temperature_mapping'] = {"RIVER_NETWORK_REACHES": "Treach_K"}
        config['phreeqc_chemical_species_names'] = species
        # PHREEQC cannot use sorption isotherm module (handled natively)
        config['si_module_name'] = "NONE"

    # Add Copernicus SS config if applicable
    if ss_method.startswith("using_copernicus_lulc"):
        if river_shp:
            config['ss_method_copernicus_basin_info'] = {
                'path_to_shp': river_shp,
                'mapping_key': 'GRU_ID',
            }
        config['ss_method_copernicus_compartment_name_for_load'] = "RIVER_NETWORK_REACHES"
        config['ss_method_copernicus_nc_lc_dir'] = copernicus_lc_dir or '/data/ESACCI-LC/'
        config['ss_method_copernicus_period'] = [1993, 2019]
        config['ss_method_copernicus_default_loads_bool'] = True
        config['ss_method_copernicus_annual_to_seasonal_loads_method'] = 'uniform'
    elif ss_method == "load_from_csv":
        config['ss_method_csv_config'] = []  # Will be populated per basin

    # Generate config files using Gen_Input_Driver
    if CONFIG_LIB_AVAILABLE:
        try:
            gJSON_lib.Gen_Input_Driver(**config)
            print(f"    [{variant_key}] Config generated: {variant_dir}")
        except Exception as e:
            print(f"    [{variant_key}] Config generation error: {e}")
            # Save raw config as fallback
            config_path = os.path.join(variant_dir, "config_raw.json")
            # Convert non-serializable types
            serializable = {}
            for k, v in config.items():
                try:
                    json.dumps(v)
                    serializable[k] = v
                except (TypeError, ValueError):
                    serializable[k] = str(v)
            with open(config_path, 'w') as f:
                json.dump(serializable, f, indent=2)
            print(f"    [{variant_key}] Raw config saved: {config_path}")
    else:
        # Save raw config when library not available
        config_path = os.path.join(variant_dir, "config_raw.json")
        serializable = {}
        for k, v in config.items():
            try:
                json.dumps(v)
                serializable[k] = v
            except (TypeError, ValueError):
                serializable[k] = str(v)
        with open(config_path, 'w') as f:
            json.dump(serializable, f, indent=2)
        print(f"    [{variant_key}] Raw config saved (library not available): {config_path}")

    return variant_dir


# =============================================================================
# MAIN WORKFLOW
# =============================================================================

def process_basin(
    tar_path: str,
    prepared_dir: str,
    variants: List[str],
    ss_method: str,
    copernicus_lc_dir: Optional[str],
) -> Dict:
    """Process a single basin: extract and generate configs for all variants."""
    from calibration_lib.observation_data.grqa_extract_stations import (
        GRQACalibrationExtractor,
    )

    # Parse basin name
    fname = Path(tar_path).name
    basename = fname.replace('.tar.gz', '')
    parts = basename.split('_')
    if len(parts) >= 4 and parts[0] == 'domain':
        basin_id = f"{parts[1]}_{parts[2]}_{parts[3]}"
    else:
        basin_id = basename

    basin_output_dir = os.path.join(prepared_dir, f"basin_{basin_id}")
    result = {
        'basin_id': basin_id,
        'variants_generated': [],
        'errors': [],
    }

    # Extract basin archive to temp dir
    with tempfile.TemporaryDirectory(prefix=f"cfg_{basin_id}_") as tmp_dir:
        try:
            print(f"  Extracting {basin_id}...")
            with tarfile.open(tar_path, 'r:gz') as tar:
                # Extract shapefiles and settings
                members = [m for m in tar.getmembers()
                          if any(x in m.name for x in
                                 ['shapefile', 'river_network', 'river_basins',
                                  'mizuRoute', 'topology'])]
                if not members:
                    members = tar.getmembers()
                tar.extractall(tmp_dir, members=members)

            # Generate configs for each variant
            for v_key in variants:
                v_info = VARIANT_TEMPLATES[v_key]
                try:
                    generate_variant_config(
                        basin_id=basin_id,
                        variant_key=v_key,
                        variant_info=v_info,
                        output_dir=basin_output_dir,
                        basin_extract_dir=tmp_dir,
                        ss_method=ss_method,
                        copernicus_lc_dir=copernicus_lc_dir,
                    )
                    result['variants_generated'].append(v_key)
                except Exception as e:
                    result['errors'].append(f"Variant {v_key}: {e}")
                    print(f"    [{v_key}] ERROR: {e}")

        except Exception as e:
            result['errors'].append(f"Extraction: {e}")
            print(f"  ERROR extracting {basin_id}: {e}")

    return result


def main():
    parser = argparse.ArgumentParser(
        description="Generate OpenWQ configs for 3 variants per century basin",
    )

    parser.add_argument('--basins-dir', required=True,
                       help='Directory containing basin tar.gz archives')
    parser.add_argument('--prepared-dir', required=True,
                       help='Output directory (century_basins_prepared)')
    parser.add_argument('--openwq-dir', default=str(MODEL_CONFIG_DIR),
                       help='Path to Model_Config directory')
    parser.add_argument('--variants', default='A,B,C,D',
                       help='Comma-separated variant letters (default: A,B,C,D)')
    parser.add_argument('--ss-method', default='using_copernicus_lulc_with_static_coeff',
                       choices=['load_from_csv', 'using_copernicus_lulc_with_static_coeff',
                               'using_copernicus_lulc_with_dynamic_coeff'],
                       help='Source/sink method')
    parser.add_argument('--copernicus-lc-dir', default=None,
                       help='Path to Copernicus land cover NetCDF directory')
    parser.add_argument('--basin-filter', default=None,
                       help='Filter basins by country code (CAN, USA)')
    parser.add_argument('--single-basin', default=None,
                       help='Process a single basin (for testing)')
    parser.add_argument('--resume', action='store_true',
                       help='Skip basins that already have all variant configs')

    args = parser.parse_args()

    # Parse variants
    variants = [v.strip().upper() for v in args.variants.split(',')]
    for v in variants:
        if v not in VARIANT_TEMPLATES:
            print(f"ERROR: Unknown variant '{v}'. Must be A, B, or C.")
            sys.exit(1)

    basins_dir = Path(args.basins_dir)
    prepared_dir = Path(args.prepared_dir)
    prepared_dir.mkdir(parents=True, exist_ok=True)

    # Find basin archives
    if args.single_basin:
        tar_files = [str(basins_dir / args.single_basin)]
    else:
        tar_files = sorted(glob.glob(str(basins_dir / "domain_*.tar.gz")))

    if args.basin_filter:
        tar_files = [f for f in tar_files if f"_{args.basin_filter}_" in Path(f).name]

    # Resume support
    if args.resume:
        remaining = []
        for tf in tar_files:
            fname = Path(tf).name.replace('.tar.gz', '')
            parts = fname.split('_')
            if len(parts) >= 4:
                basin_id = f"{parts[1]}_{parts[2]}_{parts[3]}"
            else:
                basin_id = fname
            # Check if all variant dirs exist
            all_exist = all(
                (prepared_dir / f"basin_{basin_id}" / VARIANT_TEMPLATES[v]['dir_name']).exists()
                for v in variants
            )
            if not all_exist:
                remaining.append(tf)
        tar_files = remaining
        print(f"Resuming: {len(remaining)} basins remaining")

    if not tar_files:
        print("No basins to process.")
        sys.exit(0)

    print("\n" + "=" * 70)
    print("OPENWQ CONFIG GENERATION FOR CENTURY BASINS")
    print("=" * 70)
    print(f"\nBasins:     {len(tar_files)}")
    print(f"Variants:   {', '.join(VARIANT_TEMPLATES[v]['label'] for v in variants)}")
    print(f"SS method:  {args.ss_method}")
    print(f"Output:     {prepared_dir}")

    print("\nVariant BGC templates:")
    for v in variants:
        info = VARIANT_TEMPLATES[v]
        print(f"  [{v}] {info['label']}")
        if 'template' in info:
            print(f"      Template: {info['template']}")
        elif 'templates' in info:
            for t in info['templates']:
                print(f"      Template: {t}")
        print(f"      Species:  {info['species'][:5]}{'...' if len(info['species']) > 5 else ''}")

    # Process each basin
    results = []
    for i, tar_path in enumerate(tar_files):
        fname = Path(tar_path).name.replace('.tar.gz', '')
        parts = fname.split('_')
        basin_id = f"{parts[1]}_{parts[2]}_{parts[3]}" if len(parts) >= 4 else fname
        print(f"\n[{i+1}/{len(tar_files)}] {basin_id}")

        result = process_basin(
            tar_path=tar_path,
            prepared_dir=str(prepared_dir),
            variants=variants,
            ss_method=args.ss_method,
            copernicus_lc_dir=args.copernicus_lc_dir,
        )
        results.append(result)

    # Summary
    import pandas as pd
    df = pd.DataFrame(results)
    total = len(df)
    successful = df[df['variants_generated'].apply(lambda x: len(x) == len(variants))]

    print("\n" + "=" * 70)
    print("CONFIG GENERATION SUMMARY")
    print("=" * 70)
    print(f"\nTotal basins processed: {total}")
    print(f"All variants generated: {len(successful)}")
    print(f"Partial/failed:         {total - len(successful)}")

    # Save summary
    summary_path = prepared_dir / "config_generation_summary.json"
    summary = {
        'total_basins': total,
        'variants': variants,
        'ss_method': args.ss_method,
        'results': results,
    }
    # Serialize results
    for r in summary['results']:
        r['errors'] = [str(e) for e in r.get('errors', [])]
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"\nSummary saved to: {summary_path}")


if __name__ == "__main__":
    main()
