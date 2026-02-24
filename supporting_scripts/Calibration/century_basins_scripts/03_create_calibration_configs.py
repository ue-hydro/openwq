#!/usr/bin/env python3
# Copyright 2026, Diogo Costa, diogo.costa@uevora.pt
# This file is part of OpenWQ model.
#
# This program, openWQ, is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) later version.

"""
03_create_calibration_configs.py
================================

Generate calibration configuration Python scripts for each basin+variant.

Creates 4 calibration configs per basin, each referencing:
  - Its variant-specific OpenWQ config directory
  - Shared observation data
  - Variant-specific parameter lists and evaluation counts

Output structure:
  century_basins_prepared/
    basin_CAN_01AM001_meso/
      calibration/
        calibration_config_A.py   -> 18 params, 500 DDS evals (NATIVE_BGC_FLEX)
        calibration_config_B.py   -> 21 params, 600 DDS evals (NATIVE_BGC_FLEX)
        calibration_config_C.py   -> 16 params, 450 DDS evals (NATIVE_BGC_FLEX)
        calibration_config_D.py   -> 16 params, 450 DDS evals (PHREEQC)

USAGE:
  python 03_create_calibration_configs.py \\
      --prepared-dir /path/to/century_basins_prepared \\
      [--hpc-base /scratch/user/century_calibration] \\
      [--sif-path /path/to/openwq.sif] \\
      [--variants A,B,C,D] \\
      [--basin-filter CAN]
"""

import sys
import os
import argparse
import json
import glob
from pathlib import Path
from typing import Dict, List, Optional

SCRIPT_DIR = Path(__file__).resolve().parent


# =============================================================================
# VARIANT-SPECIFIC CALIBRATION PARAMETER DEFINITIONS
# =============================================================================

# Common source/sink parameters used by all 3 variants
COMMON_SS_PARAMS = [
    {
        "name": "cropland_TN_coeff",
        "file_type": "ss_copernicus_static",
        "path": {"lulc_class": 10, "species": "TN"},
        "initial": 25.0,
        "bounds": (5, 80),
        "transform": "linear",
    },
    {
        "name": "cropland_TP_coeff",
        "file_type": "ss_copernicus_static",
        "path": {"lulc_class": 10, "species": "TP"},
        "initial": 2.0,
        "bounds": (0.5, 10),
        "transform": "linear",
    },
    {
        "name": "forest_TN_coeff",
        "file_type": "ss_copernicus_static",
        "path": {"lulc_class": 50, "species": "TN"},
        "initial": 2.0,
        "bounds": (0.5, 10),
        "transform": "linear",
    },
    {
        "name": "forest_TP_coeff",
        "file_type": "ss_copernicus_static",
        "path": {"lulc_class": 50, "species": "TP"},
        "initial": 0.2,
        "bounds": (0.05, 1),
        "transform": "linear",
    },
    {
        "name": "urban_TN_coeff",
        "file_type": "ss_copernicus_static",
        "path": {"lulc_class": 190, "species": "TN"},
        "initial": 15.0,
        "bounds": (5, 50),
        "transform": "linear",
    },
    {
        "name": "urban_TP_coeff",
        "file_type": "ss_copernicus_static",
        "path": {"lulc_class": 190, "species": "TP"},
        "initial": 1.5,
        "bounds": (0.5, 5),
        "transform": "linear",
    },
    {
        "name": "grassland_TN_coeff",
        "file_type": "ss_copernicus_static",
        "path": {"lulc_class": 30, "species": "TN"},
        "initial": 10.0,
        "bounds": (2, 30),
        "transform": "linear",
    },
]

# Variant A: Full N Cycle (nitrogen_full.json)
# Uses named reaction keys and PARAMETERS->param_name->VALUE structure
# 11 BGC params + 7 SS params = 18 total
VARIANT_A_PARAMS = [
    {
        "name": "k_hydrol",
        "file_type": "bgc_json",
        "path": ["BIOGEOCHEMISTRY_CONFIGURATION", "CYCLING_FRAMEWORKS", "FULL_N_CYCLE",
                 "PON_HYDROLYSIS", "PARAMETERS", "k_hydrol", "VALUE"],
        "initial": 0.1,
        "bounds": (0.02, 0.5),
        "transform": "log",
    },
    {
        "name": "k_min",
        "file_type": "bgc_json",
        "path": ["BIOGEOCHEMISTRY_CONFIGURATION", "CYCLING_FRAMEWORKS", "FULL_N_CYCLE",
                 "DON_MINERALIZATION", "PARAMETERS", "k_min", "VALUE"],
        "initial": 0.08,
        "bounds": (0.02, 0.3),
        "transform": "log",
    },
    {
        "name": "k_nitrif1",
        "file_type": "bgc_json",
        "path": ["BIOGEOCHEMISTRY_CONFIGURATION", "CYCLING_FRAMEWORKS", "FULL_N_CYCLE",
                 "NITRIFICATION_STEP1", "PARAMETERS", "k_nitrif1", "VALUE"],
        "initial": 0.15,
        "bounds": (0.05, 0.5),
        "transform": "log",
    },
    {
        "name": "k_nitrif2",
        "file_type": "bgc_json",
        "path": ["BIOGEOCHEMISTRY_CONFIGURATION", "CYCLING_FRAMEWORKS", "FULL_N_CYCLE",
                 "NITRIFICATION_STEP2", "PARAMETERS", "k_nitrif2", "VALUE"],
        "initial": 0.5,
        "bounds": (0.2, 2.0),
        "transform": "log",
    },
    {
        "name": "K_O2_nitrif",
        "file_type": "bgc_json",
        "path": ["BIOGEOCHEMISTRY_CONFIGURATION", "CYCLING_FRAMEWORKS", "FULL_N_CYCLE",
                 "NITRIFICATION_STEP1", "PARAMETERS", "K_O2_nitrif", "VALUE"],
        "initial": 0.5,
        "bounds": (0.2, 2.0),
        "transform": "log",
    },
    {
        "name": "k_denitrif",
        "file_type": "bgc_json",
        "path": ["BIOGEOCHEMISTRY_CONFIGURATION", "CYCLING_FRAMEWORKS", "FULL_N_CYCLE",
                 "DENITRIFICATION", "PARAMETERS", "k_denitrif", "VALUE"],
        "initial": 0.1,
        "bounds": (0.01, 0.5),
        "transform": "log",
    },
    {
        "name": "K_O2_inhib",
        "file_type": "bgc_json",
        "path": ["BIOGEOCHEMISTRY_CONFIGURATION", "CYCLING_FRAMEWORKS", "FULL_N_CYCLE",
                 "DENITRIFICATION", "PARAMETERS", "K_O2_inhib", "VALUE"],
        "initial": 0.3,
        "bounds": (0.1, 1.0),
        "transform": "log",
    },
    {
        "name": "k_vol",
        "file_type": "bgc_json",
        "path": ["BIOGEOCHEMISTRY_CONFIGURATION", "CYCLING_FRAMEWORKS", "FULL_N_CYCLE",
                 "VOLATILIZATION", "PARAMETERS", "k_vol", "VALUE"],
        "initial": 0.01,
        "bounds": (0.0, 0.1),
        "transform": "linear",
    },
    {
        "name": "k_uptake",
        "file_type": "bgc_json",
        "path": ["BIOGEOCHEMISTRY_CONFIGURATION", "CYCLING_FRAMEWORKS", "FULL_N_CYCLE",
                 "PLANT_UPTAKE_NH4", "PARAMETERS", "k_uptake", "VALUE"],
        "initial": 0.1,
        "bounds": (0.0, 0.5),
        "transform": "linear",
    },
    {
        "name": "k_immob",
        "file_type": "bgc_json",
        "path": ["BIOGEOCHEMISTRY_CONFIGURATION", "CYCLING_FRAMEWORKS", "FULL_N_CYCLE",
                 "IMMOBILIZATION", "PARAMETERS", "k_immob", "VALUE"],
        "initial": 0.02,
        "bounds": (0.001, 0.1),
        "transform": "log",
    },
    {
        "name": "k_fixation",
        "file_type": "bgc_json",
        "path": ["BIOGEOCHEMISTRY_CONFIGURATION", "CYCLING_FRAMEWORKS", "FULL_N_CYCLE",
                 "FIXATION", "PARAMETERS", "k_fix", "VALUE"],
        "initial": 0.01,
        "bounds": (0.001, 0.05),
        "transform": "log",
    },
]

# Variant B: SWAT Full Nutrients
# Uses named reaction keys and PARAMETERS->param_name->VALUE structure
# 14 BGC params + 7 SS params = 21 total
VARIANT_B_PARAMS = [
    {
        "name": "beta_min",
        "file_type": "bgc_json",
        "path": ["BIOGEOCHEMISTRY_CONFIGURATION", "CYCLING_FRAMEWORKS", "SWAT_NITROGEN_CYCLE",
                 "MINERALIZATION_ACTIVE", "PARAMETERS", "beta_min", "VALUE"],
        "initial": 0.0003,
        "bounds": (0.0001, 0.001),
        "transform": "log",
    },
    {
        "name": "beta_rsd",
        "file_type": "bgc_json",
        "path": ["BIOGEOCHEMISTRY_CONFIGURATION", "CYCLING_FRAMEWORKS", "SWAT_NITROGEN_CYCLE",
                 "MINERALIZATION_FRESH", "PARAMETERS", "beta_rsd", "VALUE"],
        "initial": 0.05,
        "bounds": (0.02, 0.1),
        "transform": "log",
    },
    {
        "name": "f_hp",
        "file_type": "bgc_json",
        "path": ["BIOGEOCHEMISTRY_CONFIGURATION", "CYCLING_FRAMEWORKS", "SWAT_NITROGEN_CYCLE",
                 "HUMIFICATION", "PARAMETERS", "f_hp", "VALUE"],
        "initial": 0.2,
        "bounds": (0.1, 0.4),
        "transform": "linear",
    },
    {
        "name": "k_nit",
        "file_type": "bgc_json",
        "path": ["BIOGEOCHEMISTRY_CONFIGURATION", "CYCLING_FRAMEWORKS", "SWAT_NITROGEN_CYCLE",
                 "NITRIFICATION", "PARAMETERS", "k_nit", "VALUE"],
        "initial": 0.1,
        "bounds": (0.01, 0.5),
        "transform": "log",
    },
    {
        "name": "k_den",
        "file_type": "bgc_json",
        "path": ["BIOGEOCHEMISTRY_CONFIGURATION", "CYCLING_FRAMEWORKS", "SWAT_NITROGEN_CYCLE",
                 "DENITRIFICATION", "PARAMETERS", "k_den", "VALUE"],
        "initial": 0.05,
        "bounds": (0.01, 0.2),
        "transform": "log",
    },
    {
        "name": "k_vol_swat",
        "file_type": "bgc_json",
        "path": ["BIOGEOCHEMISTRY_CONFIGURATION", "CYCLING_FRAMEWORKS", "SWAT_NITROGEN_CYCLE",
                 "VOLATILIZATION", "PARAMETERS", "k_vol", "VALUE"],
        "initial": 0.02,
        "bounds": (0.0, 0.1),
        "transform": "linear",
    },
    {
        "name": "U_N",
        "file_type": "bgc_json",
        "path": ["BIOGEOCHEMISTRY_CONFIGURATION", "CYCLING_FRAMEWORKS", "SWAT_NITROGEN_CYCLE",
                 "PLANT_UPTAKE_NO3", "PARAMETERS", "U_N", "VALUE"],
        "initial": 0.1,
        "bounds": (0.01, 0.5),
        "transform": "log",
    },
    {
        "name": "f_NO3",
        "file_type": "bgc_json",
        "path": ["BIOGEOCHEMISTRY_CONFIGURATION", "CYCLING_FRAMEWORKS", "SWAT_NITROGEN_CYCLE",
                 "PLANT_UPTAKE_NO3", "PARAMETERS", "f_NO3", "VALUE"],
        "initial": 0.7,
        "bounds": (0.5, 1.0),
        "transform": "linear",
    },
    {
        "name": "P_min_rate",
        "file_type": "bgc_json",
        "path": ["BIOGEOCHEMISTRY_CONFIGURATION", "CYCLING_FRAMEWORKS", "SWAT_PHOSPHORUS_CYCLE",
                 "P_MINERALIZATION_ACTIVE", "PARAMETERS", "beta_min", "VALUE"],
        "initial": 0.001,
        "bounds": (0.0001, 0.01),
        "transform": "log",
    },
    {
        "name": "P_sorption_rate",
        "file_type": "bgc_json",
        "path": ["BIOGEOCHEMISTRY_CONFIGURATION", "CYCLING_FRAMEWORKS", "SWAT_PHOSPHORUS_CYCLE",
                 "P_SORPTION", "PARAMETERS", "k_ads_P", "VALUE"],
        "initial": 0.5,
        "bounds": (0.01, 5.0),
        "transform": "log",
    },
    {
        "name": "C_min_rate",
        "file_type": "bgc_json",
        "path": ["BIOGEOCHEMISTRY_CONFIGURATION", "CYCLING_FRAMEWORKS", "SWAT_CARBON_CYCLE",
                 "C_MINERALIZATION_ACTIVE", "PARAMETERS", "beta_min", "VALUE"],
        "initial": 0.0003,
        "bounds": (0.0001, 0.001),
        "transform": "log",
    },
    {
        "name": "algae_growth",
        "file_type": "bgc_json",
        "path": ["BIOGEOCHEMISTRY_CONFIGURATION", "CYCLING_FRAMEWORKS", "SWAT_INSTREAM",
                 "ALGAE_GROWTH", "PARAMETERS", "mu_max", "VALUE"],
        "initial": 1.5,
        "bounds": (0.5, 3.0),
        "transform": "log",
    },
    {
        "name": "algae_respiration",
        "file_type": "bgc_json",
        "path": ["BIOGEOCHEMISTRY_CONFIGURATION", "CYCLING_FRAMEWORKS", "SWAT_INSTREAM",
                 "ALGAE_RESPIRATION", "PARAMETERS", "rho_a", "VALUE"],
        "initial": 0.05,
        "bounds": (0.01, 0.2),
        "transform": "log",
    },
    {
        "name": "erosion_index",
        "file_type": "bgc_json",
        "path": ["BIOGEOCHEMISTRY_CONFIGURATION", "CYCLING_FRAMEWORKS", "SWAT_SEDIMENT_TRANSPORT",
                 "EROSION", "PARAMETERS", "erosion_index", "VALUE"],
        "initial": 0.5,
        "bounds": (0.1, 2.0),
        "transform": "linear",
    },
]

# Variant C: Thermodynamic N Cycle
# Uses numeric reaction keys and parameter_values (not PARAMETERS->VALUE)
# 9 BGC params + 7 SS params = 16 total
VARIANT_C_PARAMS = [
    {
        "name": "k_max_AOB",
        "file_type": "bgc_json",
        "path": ["CYCLING_FRAMEWORKS", "NITRIFICATION_THERMO", "1", "parameter_values", "k_aob"],
        "initial": 15.0,
        "bounds": (0.5, 50.0),
        "transform": "log",
    },
    {
        "name": "Km_NH4_AOB",
        "file_type": "bgc_json",
        "path": ["CYCLING_FRAMEWORKS", "NITRIFICATION_THERMO", "1", "parameter_values", "Km_NH4"],
        "initial": 0.5,
        "bounds": (0.1, 5.0),
        "transform": "log",
    },
    {
        "name": "k_max_NOB",
        "file_type": "bgc_json",
        "path": ["CYCLING_FRAMEWORKS", "NITRIFICATION_THERMO", "2", "parameter_values", "k_nob"],
        "initial": 20.0,
        "bounds": (0.5, 50.0),
        "transform": "log",
    },
    {
        "name": "k_max_denit",
        "file_type": "bgc_json",
        "path": ["CYCLING_FRAMEWORKS", "DENITRIFICATION_THERMODYNAMIC", "1", "parameter_values", "k_max"],
        "initial": 1.0,
        "bounds": (0.1, 5.0),
        "transform": "log",
    },
    {
        "name": "Km_NO3_denit",
        "file_type": "bgc_json",
        "path": ["CYCLING_FRAMEWORKS", "DENITRIFICATION_THERMODYNAMIC", "1", "parameter_values", "Km_NO3"],
        "initial": 0.5,
        "bounds": (0.1, 2.0),
        "transform": "log",
    },
    {
        "name": "Km_DOC_denit",
        "file_type": "bgc_json",
        "path": ["CYCLING_FRAMEWORKS", "DENITRIFICATION_THERMODYNAMIC", "1", "parameter_values", "Km_DOC"],
        "initial": 2.0,
        "bounds": (0.5, 5.0),
        "transform": "log",
    },
    {
        "name": "k_max_anammox",
        "file_type": "bgc_json",
        "path": ["CYCLING_FRAMEWORKS", "ANAMMOX_THERMO", "1", "parameter_values", "k_max"],
        "initial": 0.1,
        "bounds": (0.01, 0.5),
        "transform": "log",
    },
    {
        "name": "Km_O2_inhib",
        "file_type": "bgc_json",
        "path": ["CYCLING_FRAMEWORKS", "DENITRIFICATION_THERMODYNAMIC", "1", "parameter_values", "Km_O2_inhib"],
        "initial": 0.3,
        "bounds": (0.1, 1.0),
        "transform": "log",
    },
    {
        "name": "Y_biomass",
        "file_type": "bgc_json",
        "path": ["CYCLING_FRAMEWORKS", "NITRIFICATION_THERMO", "3", "parameter_values", "Y_aob"],
        "initial": 0.15,
        "bounds": (0.05, 0.3),
        "transform": "linear",
    },
]

# Variant D: PHREEQC Nitrogen Cycle
# Uses phreeqc_pqi file_type to modify .pqi KINETICS and SOLUTION blocks
# 9 PHREEQC params + 7 SS params = 16 total
VARIANT_D_PARAMS = [
    # --- KINETICS block: nitrification (NH4 -> NO3) ---
    {
        "name": "k_nitrif_phreeqc",
        "file_type": "phreeqc_pqi",
        "path": {"block": "KINETICS", "reaction": "nitrification", "param": "-parms", "index": 0},
        "initial": 0.1,
        "bounds": (0.01, 1.0),
        "transform": "log",
    },
    {
        "name": "Km_O2_nitrif",
        "file_type": "phreeqc_pqi",
        "path": {"block": "KINETICS", "reaction": "nitrification", "param": "-parms", "index": 1},
        "initial": 0.5,
        "bounds": (0.1, 2.0),
        "transform": "log",
    },
    # --- KINETICS block: denitrification (NO3 -> N2) ---
    {
        "name": "k_denitrif_phreeqc",
        "file_type": "phreeqc_pqi",
        "path": {"block": "KINETICS", "reaction": "denitrification", "param": "-parms", "index": 0},
        "initial": 0.05,
        "bounds": (0.005, 0.5),
        "transform": "log",
    },
    {
        "name": "Km_O2_inhib_denitrif",
        "file_type": "phreeqc_pqi",
        "path": {"block": "KINETICS", "reaction": "denitrification", "param": "-parms", "index": 1},
        "initial": 0.2,
        "bounds": (0.05, 1.0),
        "transform": "log",
    },
    # --- KINETICS block: DNRA (NO3 -> NH4) ---
    {
        "name": "k_dnra",
        "file_type": "phreeqc_pqi",
        "path": {"block": "KINETICS", "reaction": "dnra", "param": "-parms", "index": 0},
        "initial": 0.01,
        "bounds": (0.001, 0.1),
        "transform": "log",
    },
    # --- SOLUTION block: initial concentrations ---
    {
        "name": "init_NO3",
        "file_type": "phreeqc_pqi",
        "path": {"block": "SOLUTION", "species": "N(5)"},
        "initial": 5.0,
        "bounds": (0.1, 20.0),
        "transform": "log",
    },
    {
        "name": "init_NH4",
        "file_type": "phreeqc_pqi",
        "path": {"block": "SOLUTION", "species": "N(-3)"},
        "initial": 0.2,
        "bounds": (0.01, 5.0),
        "transform": "log",
    },
    {
        "name": "init_DO",
        "file_type": "phreeqc_pqi",
        "path": {"block": "SOLUTION", "species": "O(0)"},
        "initial": 8.0,
        "bounds": (2.0, 14.0),
        "transform": "linear",
    },
    {
        "name": "pe_initial",
        "file_type": "phreeqc_pqi",
        "path": {"block": "SOLUTION", "species": "pe"},
        "initial": 12.0,
        "bounds": (-2.0, 14.0),
        "transform": "linear",
    },
]

# Calibration settings per variant
VARIANT_SETTINGS = {
    'A': {
        'params': VARIANT_A_PARAMS + COMMON_SS_PARAMS,
        'max_evaluations': 500,
        'label': 'nitrogen_full',
        'dir_name': 'variant_A_nitrogen_full',
        'targets_species': ["NO3-N", "NH4-N"],
        'objective_weights': {"NO3-N": 1.0, "NH4-N": 0.5},
    },
    'B': {
        'params': VARIANT_B_PARAMS + COMMON_SS_PARAMS,
        'max_evaluations': 600,
        'label': 'swat_full_nutrients',
        'dir_name': 'variant_B_swat',
        'targets_species': ["NO3-N", "NH4-N", "PO4-P"],
        'objective_weights': {"NO3-N": 1.0, "NH4-N": 0.5, "PO4-P": 0.5},
    },
    'C': {
        'params': VARIANT_C_PARAMS + COMMON_SS_PARAMS,
        'max_evaluations': 450,
        'label': 'thermodynamic',
        'dir_name': 'variant_C_thermodynamic',
        'targets_species': ["NO3-N", "NH4-N"],
        'objective_weights': {"NO3-N": 1.0, "NH4-N": 0.5},
    },
    'D': {
        'params': VARIANT_D_PARAMS + COMMON_SS_PARAMS,
        'max_evaluations': 450,
        'label': 'phreeqc_nitrogen',
        'dir_name': 'variant_D_phreeqc',
        'targets_species': ["NO3-N", "NH4-N"],
        'objective_weights': {"NO3-N": 1.0, "NH4-N": 0.5},
    },
}


# =============================================================================
# CALIBRATION CONFIG FILE GENERATOR
# =============================================================================

def generate_calibration_config_py(
    basin_id: str,
    variant_key: str,
    variant_settings: dict,
    basin_dir: str,
    hpc_base: str,
    sif_path: str,
    executable_name: str = "mizuroute_lakes_openwq_fast",
) -> str:
    """
    Generate a calibration_config_{variant}.py file for one basin+variant.

    Returns the generated file path.
    """
    calib_dir = os.path.join(basin_dir, "calibration")
    os.makedirs(calib_dir, exist_ok=True)

    config_path = os.path.join(calib_dir, f"calibration_config_{variant_key}.py")
    v = variant_settings

    # Build parameter list as Python code string
    params_str = "calibration_parameters = [\n"
    for p in v['params']:
        bounds_str = f"({p['bounds'][0]}, {p['bounds'][1]})"
        path_str = repr(p['path'])
        params_str += f"""    {{
        "name": "{p['name']}",
        "file_type": "{p['file_type']}",
        "path": {path_str},
        "initial": {p['initial']},
        "bounds": {bounds_str},
        "transform": "{p['transform']}",
    }},
"""
    params_str += "]\n"

    # Objective weights
    weights_str = repr(v['objective_weights'])

    # Target species
    targets_str = repr(v['targets_species'])

    # Resolve absolute path to calibration library
    # This path will work both locally and on HPC
    calibration_lib_dir = str(SCRIPT_DIR.parent)  # .../supporting_scripts/Calibration/

    # Write the calibration config Python file
    content = f'''#!/usr/bin/env python3
# Auto-generated calibration config for {basin_id} - Variant {variant_key} ({v['label']})
# Generated by 03_create_calibration_configs.py
#
# USAGE:
#   python calibration_config_{variant_key}.py                    # Run calibration
#   python calibration_config_{variant_key}.py --dry-run          # Validate config
#   python calibration_config_{variant_key}.py --sensitivity-only # SA only
#   python calibration_config_{variant_key}.py --resume           # Resume from checkpoint

import sys
import os
import argparse

# Add calibration library to path
# NOTE: On HPC, update CALIBRATION_LIB to point to the scripts location
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CALIBRATION_LIB = "{calibration_lib_dir}"
# Fallback: try relative path from HPC basin directory
if not os.path.isdir(CALIBRATION_LIB):
    CALIBRATION_LIB = os.environ.get("OPENWQ_CALIBRATION_LIB", CALIBRATION_LIB)
sys.path.insert(0, CALIBRATION_LIB)

from calibration_lib.calibration_driver import run_calibration, run_sensitivity_analysis
from calibration_lib.parameter_defaults import validate_parameter_bounds

# =============================================================================
# PATH CONFIGURATION
# =============================================================================

# Basin and variant identifiers
BASIN_ID = "{basin_id}"
VARIANT = "{variant_key}"
VARIANT_LABEL = "{v['label']}"

# Local paths (for dry-run validation)
LOCAL_BASIN_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))
LOCAL_VARIANT_DIR = os.path.join(LOCAL_BASIN_DIR, "{v['dir_name']}")

# HPC paths (used during actual calibration)
HPC_BASE = "{hpc_base}"
HPC_BASIN_DIR = os.path.join(HPC_BASE, "basins", "basin_" + BASIN_ID)

# Model config directory (variant-specific)
base_model_config_dir = os.path.join(HPC_BASIN_DIR, "{v['dir_name']}")

# Calibration workspace
calibration_work_dir = os.path.join(HPC_BASIN_DIR, "calibration_workspace_{variant_key}")

# Observation data (shared across variants)
observation_data_path = os.path.join(HPC_BASIN_DIR, "observations", "calibration_observations.csv")

# Test case directory (variant-specific OpenWQ configs)
test_case_dir = os.path.join(HPC_BASIN_DIR, "{v['dir_name']}")

# =============================================================================
# CONTAINER RUNTIME (Apptainer for HPC)
# =============================================================================

container_runtime = "apptainer"
apptainer_sif_path = "{sif_path}"
apptainer_bind_path = HPC_BASE + ":/code"
executable_name = "{executable_name}"
executable_args = "-g 1 1"
file_manager_path = os.path.join("/code", "basins", "basin_" + BASIN_ID,
                                  "{v['dir_name']}", "mizuroute_in", "fileManager.txt")

# Docker settings (for local testing)
docker_container_name = "docker_openwq"
docker_compose_path = ""

# =============================================================================
# CALIBRATION PARAMETERS - Variant {variant_key}: {v['label']}
# Total: {len(v['params'])} parameters -> {v['max_evaluations']} DDS evaluations
# =============================================================================

{params_str}

# =============================================================================
# CALIBRATION SETTINGS
# =============================================================================

algorithm = "DDS"
max_evaluations = {v['max_evaluations']}
n_parallel = 1
objective_function = "KGE"
temporal_resolution = "monthly"
aggregation_method = "mean"

objective_weights = {weights_str}

calibration_targets = {{
    "species": {targets_str},
    "reach_ids": "all",
    "compartments": ["RIVER_NETWORK_REACHES"],
}}

random_seed = 42

# =============================================================================
# SENSITIVITY ANALYSIS (Morris screening before DDS)
# =============================================================================

run_sensitivity_first = True
sensitivity_method = "morris"
sensitivity_morris_trajectories = 10
sensitivity_morris_levels = 4
sensitivity_sobol_samples = 1024
sensitivity_threshold = 0.1

# =============================================================================
# HPC SETTINGS
# =============================================================================

hpc_enabled = True
hpc_scheduler = "slurm"
hpc_partition = "standard"
hpc_walltime = "24:00:00"
hpc_nodes = 1
hpc_tasks_per_node = 4
hpc_memory = "8G"
hpc_max_concurrent_jobs = 50

# =============================================================================
# OBSERVATION DATA
# =============================================================================

observation_data_source = "csv"  # Pre-extracted by 01_extract_grqa_all.py


# =============================================================================
# RUN CALIBRATION
# =============================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=f"OpenWQ Calibration: {{BASIN_ID}} - Variant {{VARIANT}} ({{VARIANT_LABEL}})"
    )
    parser.add_argument("--resume", action="store_true", help="Resume from checkpoint")
    parser.add_argument("--sensitivity-only", action="store_true", help="SA only")
    parser.add_argument("--dry-run", action="store_true", help="Validate config")
    args = parser.parse_args()

    config = {{
        "base_model_config_dir": base_model_config_dir,
        "calibration_work_dir": calibration_work_dir,
        "observation_data_path": observation_data_path,
        "test_case_dir": test_case_dir,
        "container_runtime": container_runtime,
        "docker_container_name": docker_container_name,
        "docker_compose_path": docker_compose_path,
        "apptainer_sif_path": apptainer_sif_path,
        "apptainer_bind_path": apptainer_bind_path,
        "executable_name": executable_name,
        "executable_args": executable_args,
        "file_manager_path": file_manager_path,
        "calibration_parameters": calibration_parameters,
        "algorithm": algorithm,
        "max_evaluations": max_evaluations,
        "n_parallel": n_parallel,
        "objective_function": objective_function,
        "objective_weights": objective_weights,
        "calibration_targets": calibration_targets,
        "random_seed": random_seed,
        "temporal_resolution": temporal_resolution,
        "aggregation_method": aggregation_method,
        "run_sensitivity_first": run_sensitivity_first,
        "sensitivity_method": sensitivity_method,
        "sensitivity_morris_trajectories": sensitivity_morris_trajectories,
        "sensitivity_morris_levels": sensitivity_morris_levels,
        "sensitivity_sobol_samples": sensitivity_sobol_samples,
        "sensitivity_threshold": sensitivity_threshold,
        "hpc_enabled": hpc_enabled,
        "hpc_scheduler": hpc_scheduler,
        "hpc_partition": hpc_partition,
        "hpc_walltime": hpc_walltime,
        "hpc_nodes": hpc_nodes,
        "hpc_tasks_per_node": hpc_tasks_per_node,
        "hpc_memory": hpc_memory,
        "hpc_max_concurrent_jobs": hpc_max_concurrent_jobs,
        "observation_data_source": observation_data_source,
    }}

    if args.dry_run:
        print("=" * 60)
        print(f"DRY RUN: {{BASIN_ID}} - Variant {{VARIANT}} ({{VARIANT_LABEL}})")
        print("=" * 60)
        print(f"Parameters: {{len(calibration_parameters)}}")
        for p in calibration_parameters:
            print(f"  - {{p['name']}}: {{p['bounds']}} ({{p['transform']}})")
            warnings = validate_parameter_bounds(p)
            for w in warnings:
                print(f"    ! {{w}}")
        print(f"Max evaluations: {{max_evaluations}}")
        print(f"Objective: {{objective_function}}")
        print(f"Temporal resolution: {{temporal_resolution}}")
        sys.exit(0)

    if args.sensitivity_only:
        print(f"Running sensitivity analysis: {{BASIN_ID}} - Variant {{VARIANT}}")
        results = run_sensitivity_analysis(**config)
    else:
        print(f"Running calibration: {{BASIN_ID}} - Variant {{VARIANT}}")
        results = run_calibration(resume=args.resume, **config)

    print("\\n" + "=" * 60)
    print(f"COMPLETE: {{BASIN_ID}} - Variant {{VARIANT}}")
    print("=" * 60)
'''

    with open(config_path, 'w') as f:
        f.write(content)

    os.chmod(config_path, 0o755)
    return config_path


# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Generate calibration configs for century basins",
    )
    parser.add_argument('--prepared-dir', required=True,
                       help='Prepared basins directory')
    parser.add_argument('--hpc-base', default='/scratch/user/century_calibration',
                       help='HPC base path for calibration workspace')
    parser.add_argument('--sif-path', default='/scratch/user/century_calibration/openwq.sif',
                       help='Path to Apptainer .sif file on HPC')
    parser.add_argument('--executable', default='mizuroute_lakes_openwq_fast',
                       help='Executable name')
    parser.add_argument('--variants', default='A,B,C,D',
                       help='Comma-separated variants')
    parser.add_argument('--basin-filter', default=None,
                       help='Filter by country code')

    args = parser.parse_args()

    variants = [v.strip().upper() for v in args.variants.split(',')]
    prepared_dir = Path(args.prepared_dir)

    # Find all basin directories
    basin_dirs = sorted(glob.glob(str(prepared_dir / "basin_*")))
    if not basin_dirs:
        print(f"No basin directories found in {prepared_dir}")
        sys.exit(1)

    if args.basin_filter:
        basin_dirs = [d for d in basin_dirs if args.basin_filter in Path(d).name]

    print("\n" + "=" * 70)
    print("CALIBRATION CONFIG GENERATION")
    print("=" * 70)
    print(f"\nBasins:       {len(basin_dirs)}")
    print(f"Variants:     {variants}")
    print(f"HPC base:     {args.hpc_base}")
    print(f"SIF path:     {args.sif_path}")

    for v in variants:
        vs = VARIANT_SETTINGS[v]
        print(f"\n  Variant {v} ({vs['label']}):")
        print(f"    Parameters:   {len(vs['params'])}")
        print(f"    Evaluations:  {vs['max_evaluations']}")
        print(f"    Targets:      {vs['targets_species']}")

    # Generate configs
    total_generated = 0
    total_errors = 0

    for i, basin_dir in enumerate(basin_dirs):
        basin_id = Path(basin_dir).name.replace('basin_', '')
        print(f"\n[{i+1}/{len(basin_dirs)}] {basin_id}")

        for v_key in variants:
            v_settings = VARIANT_SETTINGS[v_key]

            # Check if variant config directory exists
            variant_dir = os.path.join(basin_dir, v_settings['dir_name'])
            if not os.path.exists(variant_dir):
                print(f"  [{v_key}] WARNING: Variant dir missing: {variant_dir}")
                # Still generate the calibration config for HPC
                # (variant configs will be generated separately)

            # Check if observations exist
            obs_file = os.path.join(basin_dir, "observations", "calibration_observations.csv")
            if not os.path.exists(obs_file):
                print(f"  [{v_key}] WARNING: No observations file")

            try:
                config_path = generate_calibration_config_py(
                    basin_id=basin_id,
                    variant_key=v_key,
                    variant_settings=v_settings,
                    basin_dir=basin_dir,
                    hpc_base=args.hpc_base,
                    sif_path=args.sif_path,
                    executable_name=args.executable,
                )
                print(f"  [{v_key}] Generated: {os.path.basename(config_path)}")
                total_generated += 1
            except Exception as e:
                print(f"  [{v_key}] ERROR: {e}")
                total_errors += 1

    print("\n" + "=" * 70)
    print("CALIBRATION CONFIG SUMMARY")
    print("=" * 70)
    print(f"\nTotal configs generated: {total_generated}")
    print(f"Errors: {total_errors}")
    print(f"\nOutput: {prepared_dir}/basin_*/calibration/calibration_config_[A|B|C].py")


if __name__ == "__main__":
    main()
