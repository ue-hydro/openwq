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

"""
OpenWQ Calibration Configuration Template
=========================================

Copy this file to your working directory and modify for your calibration setup.
Run: python your_calibration_config.py [--resume]

SECTIONS:
  1. Path Configuration
  2. Container Runtime Configuration
  3. Calibration Parameters
  4. Calibration Settings
  5. Sensitivity Analysis Settings
  6. HPC Settings (Apptainer)
  7. Run Calibration

Full documentation: https://openwq.readthedocs.io
"""

import sys
import os
import argparse

# Add calibration library to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from calibration_driver import run_calibration, run_sensitivity_analysis


# =============================================================================
# 1. PATH CONFIGURATION
# =============================================================================

# Path to the directory containing your model configuration
# This should be the directory containing template_model_config.py
base_model_config_dir = "/Users/diogocosta/Documents/openwq_code/6_mizuroute_cslm_openwq/route/build/openwq/openwq/supporting_scripts/Model_Config"

# Working directory for calibration runs
# Each evaluation will create a subdirectory here
calibration_work_dir = "/Users/diogocosta/Documents/openwq_code/6_mizuroute_cslm_openwq/calibration_workspace"

# Observation data file (CSV format)
# Format: datetime,reach_id,species,value,units,source
observation_data_path = "/path/to/your/observations.csv"

# Path to the test case directory (containing mizuroute_in, etc.)
test_case_dir = "/Users/diogocosta/Documents/openwq_code/6_mizuroute_cslm_openwq/test_case_2"


# =============================================================================
# 2. CONTAINER RUNTIME CONFIGURATION
# =============================================================================

# Container runtime: "docker" (local) or "apptainer" (HPC)
container_runtime = "docker"

# --- Docker settings ---
docker_container_name = "docker_openwq"
docker_compose_path = "/Users/diogocosta/Documents/openwq_code/6_mizuroute_cslm_openwq/route/build/openwq/openwq/docker-compose.yml"

# --- Apptainer settings (for HPC) ---
apptainer_sif_path = "/path/to/openwq.sif"
apptainer_bind_path = "/scratch/user/openwq_code:/code"  # host:container

# --- Executable settings ---
executable_name = "mizuroute_lakes_cslm_openwq_fast"  # or mizuroute_lakes_cslm_openwq_debug
executable_args = "-g 1 1"  # Grid arguments
file_manager_path = "/code/openwq_code/6_mizuroute_cslm_openwq/test_case_2/mizuroute_in/fileManager.txt"


# =============================================================================
# 3. CALIBRATION PARAMETERS
# =============================================================================
#
# Each parameter is defined with:
#   - name: Unique identifier (used in results)
#   - file_type: Type of configuration file to modify
#   - path: Location of parameter within the file
#   - initial: Starting value
#   - bounds: (min, max) tuple
#   - transform: "linear" (default) or "log" for log-transformed optimization
#   - dtype: "float" (default) or "int" for integer parameters
#
# Supported file_types:
#   - "bgc_json": NATIVE_BGC_FLEX rate constants
#   - "phreeqc_pqi": PHREEQC input file parameters
#   - "sorption_json": Freundlich/Langmuir isotherm parameters
#   - "sediment_json": HYPE_HBVSED or HYPE_MMF parameters
#   - "transport_json": Dispersion coefficients
#   - "lateral_exchange_json": Exchange coefficients
#   - "ss_csv_scale": Source/sink CSV load scaling
#   - "ss_copernicus_static": Copernicus LULC export coefficients
#   - "ss_copernicus_dynamic": Climate-adjusted export coefficients
#   - "ss_climate_param": Climate response parameters
#   - "ss_ml_param": ML model hyperparameters
# =============================================================================

calibration_parameters = [
    # -------------------------------------------------------------------------
    # BGC Rate Constants (NATIVE_BGC_FLEX)
    # Path: ["CYCLING_FRAMEWORKS", "<framework>", "<rxn_num>", "parameter_values", "<param>"]
    # -------------------------------------------------------------------------
    {
        "name": "k_nitrification",
        "file_type": "bgc_json",
        "path": ["CYCLING_FRAMEWORKS", "N_cycle", "3", "parameter_values", "k"],
        "initial": 0.03,
        "bounds": (0.001, 0.5),
        "transform": "log"
    },
    {
        "name": "k_denitrification",
        "file_type": "bgc_json",
        "path": ["CYCLING_FRAMEWORKS", "N_cycle", "9", "parameter_values", "k"],
        "initial": 0.000001,
        "bounds": (1e-8, 0.01),
        "transform": "log"
    },
    {
        "name": "k_mineralization_active",
        "file_type": "bgc_json",
        "path": ["CYCLING_FRAMEWORKS", "N_cycle", "4", "parameter_values", "k"],
        "initial": 0.1,
        "bounds": (0.01, 1.0),
        "transform": "log"
    },

    # -------------------------------------------------------------------------
    # PHREEQC Parameters (if using bgc_module_name = "PHREEQC")
    # Uncomment and modify as needed
    # -------------------------------------------------------------------------
    # {
    #     "name": "initial_NO3",
    #     "file_type": "phreeqc_pqi",
    #     "path": {"block": "SOLUTION", "species": "N(5)"},
    #     "initial": 2.0,
    #     "bounds": (0.1, 10.0)
    # },
    # {
    #     "name": "pCO2_log",
    #     "file_type": "phreeqc_pqi",
    #     "path": {"block": "EQUILIBRIUM_PHASES", "phase": "CO2(g)", "field": "si"},
    #     "initial": -3.5,
    #     "bounds": (-4.5, -2.5)
    # },

    # -------------------------------------------------------------------------
    # Sorption Isotherm Parameters (Freundlich)
    # -------------------------------------------------------------------------
    # {
    #     "name": "NH4_Kfr",
    #     "file_type": "sorption_json",
    #     "path": {"module": "FREUNDLICH", "species": "NH4-N", "param": "Kfr"},
    #     "initial": 0.5,
    #     "bounds": (0.01, 10.0),
    #     "transform": "log"
    # },

    # -------------------------------------------------------------------------
    # Sediment Transport Parameters (HYPE_HBVSED)
    # -------------------------------------------------------------------------
    # {
    #     "name": "erosion_index",
    #     "file_type": "sediment_json",
    #     "path": {"module": "HYPE_HBVSED", "param": "EROSION_INDEX"},
    #     "initial": 0.4,
    #     "bounds": (0.1, 1.0)
    # },

    # -------------------------------------------------------------------------
    # Transport Parameters
    # -------------------------------------------------------------------------
    # {
    #     "name": "dispersion_x",
    #     "file_type": "transport_json",
    #     "path": {"param": "DISPERSION_X_M2_S"},
    #     "initial": 0.3,
    #     "bounds": (0.01, 10.0),
    #     "transform": "log"
    # },

    # -------------------------------------------------------------------------
    # Source/Sink: CSV Load Scaling
    # -------------------------------------------------------------------------
    {
        "name": "fertilizer_N_scale",
        "file_type": "ss_csv_scale",
        "path": {"species": "all"},  # Scale all species, or specify e.g., "N_ORG_active"
        "initial": 1.0,
        "bounds": (0.1, 5.0)
    },

    # -------------------------------------------------------------------------
    # Source/Sink: Copernicus LULC Export Coefficients (static)
    # -------------------------------------------------------------------------
    # {
    #     "name": "cropland_TN_coeff",
    #     "file_type": "ss_copernicus_static",
    #     "path": {"lulc_class": 10, "species": "TN"},
    #     "initial": 20.0,
    #     "bounds": (5.0, 50.0)
    # },

    # -------------------------------------------------------------------------
    # Source/Sink: Climate Response Parameters (dynamic)
    # -------------------------------------------------------------------------
    # {
    #     "name": "precip_scaling_power",
    #     "file_type": "ss_climate_param",
    #     "path": "ss_climate_precip_scaling_power",
    #     "initial": 1.0,
    #     "bounds": (0.5, 2.0)
    # },
    # {
    #     "name": "Q10_biological",
    #     "file_type": "ss_climate_param",
    #     "path": "ss_climate_temp_q10",
    #     "initial": 2.0,
    #     "bounds": (1.5, 3.0)
    # },
]


# =============================================================================
# 4. CALIBRATION SETTINGS
# =============================================================================

# Optimization algorithm: "DDS" (recommended), "RANDOM"
algorithm = "DDS"

# Maximum number of model evaluations
max_evaluations = 500

# Number of parallel evaluations
# For Docker: typically 1 (sequential)
# For Apptainer/HPC: can be 10-100 depending on cluster resources
n_parallel = 1

# Objective function: "RMSE", "NSE", "KGE"
# KGE (Kling-Gupta Efficiency) is recommended for hydrological/WQ models
objective_function = "KGE"

# Species weights for multi-species calibration
# Higher weight = more influence on objective function
objective_weights = {
    "NO3-N": 1.0,
    "NH4-N": 0.5,
}

# Calibration targets: which species and reaches to match
calibration_targets = {
    "species": ["NO3-N", "NH4-N"],
    "reach_ids": "all",  # or list of specific IDs: [1200014181, 200014181]
    "compartments": ["RIVER_NETWORK_REACHES"]
}

# Random seed for reproducibility (None for random)
random_seed = 42


# =============================================================================
# 5. SENSITIVITY ANALYSIS SETTINGS (OPTIONAL)
# =============================================================================
#
# Run sensitivity analysis before calibration to identify influential parameters.
# This is HIGHLY RECOMMENDED for large parameter sets (10+ parameters).
# =============================================================================

# Run sensitivity analysis before optimization
run_sensitivity_first = False

# Sensitivity analysis method: "morris" (screening) or "sobol" (quantitative)
sensitivity_method = "morris"

# Morris method settings
sensitivity_morris_trajectories = 10  # Number of trajectories (r)
sensitivity_morris_levels = 4  # Number of grid levels (p)

# Sobol method settings
sensitivity_sobol_samples = 1024  # Number of samples (N)

# Threshold for influential parameters (relative to max mu_star)
sensitivity_threshold = 0.1


# =============================================================================
# 6. HPC SETTINGS (for Apptainer)
# =============================================================================

# Enable HPC mode (generates SLURM/PBS job scripts)
hpc_enabled = False

# Job scheduler: "slurm" or "pbs"
hpc_scheduler = "slurm"

# SLURM/PBS settings
hpc_partition = "standard"
hpc_walltime = "24:00:00"
hpc_nodes = 1
hpc_tasks_per_node = 4
hpc_memory = "8G"

# Job array settings (for sensitivity analysis)
hpc_max_concurrent_jobs = 50


# =============================================================================
# 7. RUN CALIBRATION
# =============================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="OpenWQ Calibration Framework")
    parser.add_argument("--resume", action="store_true",
                       help="Resume from checkpoint")
    parser.add_argument("--sensitivity-only", action="store_true",
                       help="Run only sensitivity analysis")
    parser.add_argument("--dry-run", action="store_true",
                       help="Validate configuration without running")
    args = parser.parse_args()

    # Collect all configuration into a dictionary
    config = {
        # Paths
        "base_model_config_dir": base_model_config_dir,
        "calibration_work_dir": calibration_work_dir,
        "observation_data_path": observation_data_path,
        "test_case_dir": test_case_dir,

        # Container runtime
        "container_runtime": container_runtime,
        "docker_container_name": docker_container_name,
        "docker_compose_path": docker_compose_path,
        "apptainer_sif_path": apptainer_sif_path,
        "apptainer_bind_path": apptainer_bind_path,
        "executable_name": executable_name,
        "executable_args": executable_args,
        "file_manager_path": file_manager_path,

        # Calibration parameters
        "calibration_parameters": calibration_parameters,

        # Calibration settings
        "algorithm": algorithm,
        "max_evaluations": max_evaluations,
        "n_parallel": n_parallel,
        "objective_function": objective_function,
        "objective_weights": objective_weights,
        "calibration_targets": calibration_targets,
        "random_seed": random_seed,

        # Sensitivity analysis
        "run_sensitivity_first": run_sensitivity_first,
        "sensitivity_method": sensitivity_method,
        "sensitivity_morris_trajectories": sensitivity_morris_trajectories,
        "sensitivity_morris_levels": sensitivity_morris_levels,
        "sensitivity_sobol_samples": sensitivity_sobol_samples,
        "sensitivity_threshold": sensitivity_threshold,

        # HPC settings
        "hpc_enabled": hpc_enabled,
        "hpc_scheduler": hpc_scheduler,
        "hpc_partition": hpc_partition,
        "hpc_walltime": hpc_walltime,
        "hpc_nodes": hpc_nodes,
        "hpc_tasks_per_node": hpc_tasks_per_node,
        "hpc_memory": hpc_memory,
        "hpc_max_concurrent_jobs": hpc_max_concurrent_jobs,
    }

    if args.dry_run:
        print("=" * 60)
        print("DRY RUN - Configuration Validation")
        print("=" * 60)
        print(f"Parameters to calibrate: {len(calibration_parameters)}")
        for p in calibration_parameters:
            print(f"  - {p['name']}: {p['bounds']} (transform={p.get('transform', 'linear')})")
        print(f"Algorithm: {algorithm}")
        print(f"Max evaluations: {max_evaluations}")
        print(f"Container runtime: {container_runtime}")
        print(f"Objective function: {objective_function}")
        print("Configuration is valid.")
        sys.exit(0)

    if args.sensitivity_only:
        print("Running sensitivity analysis only...")
        results = run_sensitivity_analysis(**config)
    else:
        print("Running calibration...")
        results = run_calibration(resume=args.resume, **config)

    print("\n" + "=" * 60)
    print("CALIBRATION COMPLETE")
    print("=" * 60)
    print(f"Results saved to: {calibration_work_dir}/results/")
