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
  3. Calibration Parameters (with scientifically-grounded min/max bounds)
  4. Calibration Settings
  5. Sensitivity Analysis Settings
  6. HPC Settings (Apptainer)
  7. Run Calibration

Parameter Reference:
  See calibration_lib/parameter_defaults.py for comprehensive list of all
  calibratable parameters with their scientific ranges and documentation.

Full documentation: https://openwq.readthedocs.io
"""

import sys
import os
import argparse

# Add calibration library to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from calibration_driver import run_calibration, run_sensitivity_analysis
from calibration_lib.parameter_defaults import (
    BGC_PARAMETERS, PHREEQC_PARAMETERS, SORPTION_PARAMETERS,
    SEDIMENT_PARAMETERS, TRANSPORT_PARAMETERS, LATERAL_EXCHANGE_PARAMETERS,
    SOURCE_SINK_PARAMETERS, create_parameter_entry, validate_parameter_bounds
)

# GRQA extraction tools (optional - only needed if observation_data_source = "grqa")
try:
    from observation_data.grqa_extract_stations import GRQACalibrationExtractor
    from observation_data.prepare_calibration_observations import CalibrationObservationConverter
    GRQA_AVAILABLE = True
except ImportError:
    GRQA_AVAILABLE = False

# Copernicus workflow tools (optional - only needed if observation_data_source = "copernicus")
try:
    from observation_data.copernicus_observations import CopernicusObservationGenerator
    COPERNICUS_AVAILABLE = True
except ImportError:
    COPERNICUS_AVAILABLE = False


# =============================================================================
# 1. PATH CONFIGURATION
# =============================================================================

# Path to the directory containing your model configuration
# This should be the directory containing template_model_config.py
base_model_config_dir = "/Users/diogocosta/Documents/openwq_code/6_mizuroute_cslm_openwq/route/build/openwq/openwq/supporting_scripts/Model_Config"

# Working directory for calibration runs
# Each evaluation will create a subdirectory here
calibration_work_dir = "/Users/diogocosta/Documents/openwq_code/6_mizuroute_cslm_openwq/calibration_workspace"

# Path to the test case directory (containing mizuroute_in, etc.)
test_case_dir = "/Users/diogocosta/Documents/openwq_code/6_mizuroute_cslm_openwq/test_case_2"


# =============================================================================
# 1B. OBSERVATION DATA CONFIGURATION
# =============================================================================
#
# Three options for observation data:
#   Option A: "csv" - Provide pre-formatted observation CSV file
#   Option B: "grqa" - Extract from GRQA database (Global River Water Quality Archive)
#   Option C: "copernicus" - Generate synthetic observations from Copernicus data
#
# The calibration framework requires observations in CSV format with columns:
#   datetime, reach_id, species, value, units, source, uncertainty, quality_flag

# --- Select observation data source ---
# Options: "csv", "grqa", "copernicus"
observation_data_source = "csv"

# --- Option A: Pre-formatted observation file (used when observation_data_source = "csv") ---
# Provide path to existing observation CSV file
observation_data_path = "/path/to/your/observations.csv"

# --- Option B: Extract from GRQA database (used when observation_data_source = "grqa") ---
# Automatically downloads and processes data from GRQA (https://zenodo.org/records/15335450)
use_grqa = False  # Deprecated: use observation_data_source = "grqa" instead

# GRQA Configuration (only used if observation_data_source = "grqa")
grqa_config = {
    # -------------------------------------------------------------------------
    # DATA SOURCE: Choose between downloading or using local files
    # -------------------------------------------------------------------------
    # Option 1: Download from Zenodo (default)
    #   Set "local_data_path" to None - data will be downloaded automatically
    #
    # Option 2: Use already-downloaded GRQA data
    #   Set "local_data_path" to the folder containing GRQA CSV files
    #   Expected structure:
    #     local_data_path/
    #       ├── GRQA_data_v1.3/
    #       │   ├── GRQA_data_v1.3_NO3.csv
    #       │   ├── GRQA_data_v1.3_NH4.csv
    #       │   ├── GRQA_data_v1.3_TN.csv
    #       │   └── ... (other parameter files)
    #       └── GRQA_meta_v1.3.csv (station metadata)
    # -------------------------------------------------------------------------
    "local_data_path": None,  # e.g., "/data/GRQA" or None to download

    # River network shapefile for spatial matching
    "river_network_shapefile": "/path/to/river_network.shp",

    # Column in shapefile containing reach IDs
    "reach_id_column": "seg_id",

    # Bounding box to filter stations [min_lon, min_lat, max_lon, max_lat]
    # Set to None to use shapefile extent
    "bounding_box": None,  # e.g., [-125, 24, -66, 50] for CONUS

    # Time period for observations (ISO format)
    "start_date": "2000-01-01",
    "end_date": "2020-12-31",

    # Maximum distance (meters) to match stations to river reaches
    "max_station_distance_m": 500,

    # Minimum number of observations required per station
    "min_observations": 10,

    # Output directory for GRQA extraction
    "output_dir": None,  # Default: calibration_work_dir/observation_data

    # Species mapping: GRQA parameter names -> Model species names
    # Only parameters listed here will be downloaded and processed
    # See GRQA documentation for full parameter list:
    # https://zenodo.org/records/15335450
    "species_mapping": {
        # Nitrogen species
        "NO3": "NO3-N",      # Nitrate
        "NH4": "NH4-N",      # Ammonium
        "TN": "TN",          # Total Nitrogen
        "NO2": "NO2-N",      # Nitrite
        "DON": "DON",        # Dissolved Organic Nitrogen

        # Phosphorus species
        "PO4": "PO4-P",      # Orthophosphate
        "TP": "TP",          # Total Phosphorus
        "DP": "DP",          # Dissolved Phosphorus

        # General water quality
        "DO": "DO",          # Dissolved Oxygen
        "BOD": "BOD",        # Biochemical Oxygen Demand
        "COD": "COD",        # Chemical Oxygen Demand
        "DOC": "DOC",        # Dissolved Organic Carbon
        "TOC": "TOC",        # Total Organic Carbon
        "TSS": "TSS",        # Total Suspended Solids

        # Uncomment additional parameters as needed:
        # "Chl-a": "Chla",   # Chlorophyll-a
        # "EC": "EC",        # Electrical Conductivity
        # "pH": "pH",        # pH
        # "TDS": "TDS",      # Total Dissolved Solids
        # "Temp": "Temp",    # Temperature
    },

    # Unit conversions applied during extraction
    # Format: {grqa_unit: {target_unit: conversion_factor}}
    "unit_conversions": {
        "ug/l": {"mg/l": 0.001},
        "µg/l": {"mg/l": 0.001},
        "mg/l": {"mg/l": 1.0},
        "ppm": {"mg/l": 1.0},
    },

    # Target units for each species (optional)
    # If not specified, original units are preserved
    "target_units": {
        "NO3-N": "mg/l",
        "NH4-N": "mg/l",
        "TN": "mg/l",
        "TP": "mg/l",
        "DO": "mg/l",
    },

    # Quality control flags to include (GRQA quality flags)
    # Good/Acceptable values only
    "accepted_quality_flags": ["GOOD", "ACCEPTABLE", "A", "B", None],
}

# --- Option C: Copernicus workflow (used when observation_data_source = "copernicus") ---
# Generate synthetic observations using Copernicus land cover and climate data
# This uses export coefficients to estimate nutrient loads based on land use
copernicus_config = {
    # -------------------------------------------------------------------------
    # DATA PATHS: Specify paths to Copernicus data files
    # -------------------------------------------------------------------------
    # Land Cover Data:
    #   Download from: https://land.copernicus.eu/global/products/lc
    #   Supported formats: GeoTIFF (.tif), NetCDF (.nc)
    #   Resolution: 100m recommended for catchment-scale analysis
    #
    # Climate Data (optional):
    #   Download from: https://cds.climate.copernicus.eu/
    #   Supported: ERA5 NetCDF files with precipitation (tp) and temperature (t2m)
    # -------------------------------------------------------------------------

    # River network shapefile for spatial matching
    "river_network_shapefile": "/path/to/river_network.shp",

    # Column in shapefile containing reach IDs
    "reach_id_column": "seg_id",

    # Copernicus land cover data path (GeoTIFF or NetCDF)
    # Can be a single file or directory containing yearly files
    # Examples:
    #   Single file: "/data/copernicus/CGLS-LC100_2019.tif"
    #   Directory:   "/data/copernicus/land_cover/" (will use files matching time period)
    "land_cover_path": "/path/to/copernicus_land_cover.tif",

    # Climate data for dynamic export coefficients (optional)
    # ERA5 or other climate reanalysis data in NetCDF format
    # Set to None to use static export coefficients (no climate adjustment)
    # Examples:
    #   Single file: "/data/era5/era5_daily_2010-2020.nc"
    #   Directory:   "/data/era5/" (will load files matching time period)
    "climate_data_path": None,  # Set to enable climate-adjusted exports

    # Time period for synthetic observations
    "start_date": "2000-01-01",
    "end_date": "2020-12-31",

    # Temporal resolution for synthetic observations
    # Options: "daily", "weekly", "monthly", "annual"
    "temporal_resolution": "monthly",

    # Output directory for Copernicus workflow
    "output_dir": None,  # Default: calibration_work_dir/observation_data

    # Species to generate (must match model species names)
    "target_species": ["TN", "TP", "TSS"],

    # Export coefficients by land cover class (kg/ha/yr)
    # Copernicus Global Land Cover classes:
    #   10: Cropland, 20: Forest, 30: Grassland, 40: Shrubland,
    #   50: Wetland, 60: Water, 70: Tundra, 80: Artificial,
    #   90: Bare/sparse vegetation, 100: Snow/Ice
    "export_coefficients": {
        # Cropland - highest nutrient exports
        10: {"TN": 25.0, "TP": 2.0, "TSS": 500.0},
        # Forest - low background exports
        20: {"TN": 2.0, "TP": 0.2, "TSS": 50.0},
        # Grassland - moderate exports
        30: {"TN": 10.0, "TP": 0.8, "TSS": 150.0},
        # Shrubland
        40: {"TN": 5.0, "TP": 0.4, "TSS": 100.0},
        # Wetland - nutrient sink/source
        50: {"TN": 3.0, "TP": 0.3, "TSS": 30.0},
        # Water bodies
        60: {"TN": 0.5, "TP": 0.05, "TSS": 10.0},
        # Artificial/Urban - high runoff
        80: {"TN": 15.0, "TP": 1.5, "TSS": 300.0},
        # Bare/sparse vegetation
        90: {"TN": 1.0, "TP": 0.1, "TSS": 200.0},
    },

    # Climate response parameters (only used if climate_data_path is set)
    "climate_params": {
        # Precipitation scaling power (1.0 = linear, >1.0 = nonlinear erosion response)
        "precip_scaling_power": 1.5,
        # Temperature Q10 factor for biological processes
        "temp_q10": 2.0,
        # Reference temperature (°C)
        "temp_reference_c": 15.0,
    },

    # Uncertainty estimation
    "uncertainty_fraction": 0.3,  # 30% uncertainty on synthetic values

    # Add noise to synthetic observations for realism
    "add_noise": True,
    "noise_cv": 0.2,  # Coefficient of variation for noise
}


# =============================================================================
# 2. CONTAINER RUNTIME CONFIGURATION
# =============================================================================

# Container runtime: "docker" (local) or "apptainer" (HPC)
container_runtime = "docker"

# --- Docker settings ---
docker_container_name = "docker_openwq"
docker_compose_path = "/Users/diogocosta/Documents/openwq_code/6_mizuroute_cslm_openwq/route/build/openwq/openwq/containers/docker-compose.yml"

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
#   - bounds: (min, max) tuple - see parameter_defaults.py for scientific ranges
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
#
# SCIENTIFIC PARAMETER RANGES (from calibration_lib/parameter_defaults.py):
# -------------------------------------------------------------------------
# BGC (NATIVE_BGC_FLEX) - 1/day rates:
#   k_nitrification:        (0.001, 0.5)    log - Chapra: 0.01-0.2/day typical
#   k_denitrification:      (1e-8, 0.01)    log - Highly variable, anoxic
#   k_mineralization_active:(0.01, 1.0)     log - Fast pool decomposition
#   k_mineralization_stable:(0.0001, 0.01)  log - Slow pool, 10-100x slower
#   k_immobilization:       (0.001, 0.5)    log - Microbial uptake
#   k_plant_uptake_N:       (0.0, 0.5)      linear - Seasonal
#   k_volatilization:       (0.0, 0.1)      linear - pH/temp dependent
#   k_P_mineralization:     (0.001, 0.2)    log - Slower than N
#   k_P_adsorption:         (0.01, 5.0)     log - Fast process
#   k_BOD_decay:            (0.01, 1.0)     log - 0.1-0.5/day typical
#   k_reaeration:           (0.1, 20.0)     log - Depth/velocity dependent
#   theta (temp coeff):     (1.02, 1.12)    linear - Arrhenius factors
#
# PHREEQC - concentrations mg/L, SI dimensionless:
#   Ca:                     (5.0, 500.0)    log - 10-200 mg/L typical
#   Mg:                     (1.0, 200.0)    log - 5-50 mg/L typical
#   Na:                     (1.0, 10000.0)  log - Freshwater to brackish
#   N(5) as NO3:            (0.01, 50.0)    log - Pristine to agricultural
#   N(-3) as NH4:           (0.001, 10.0)   log - Surface water <0.5 mg/L
#   P as PO4:               (0.001, 5.0)    log - Oligo to eutrophic
#   pCO2_log (atm):         (-4.5, -1.5)    linear - Atmospheric to organic-rich
#   Calcite_SI:             (-3.0, 1.0)     linear - Under to supersaturated
#
# Sorption - Freundlich/Langmuir:
#   Kfr (L/kg):             (0.01, 100.0)   log - Sandy 0.1-1, clay 5-50
#   Nfr (exponent):         (0.3, 1.0)      linear - 1.0 = linear
#   qmax (mg/kg):           (10.0, 2000.0)  log - Sandy to clay
#   KL (L/mg):              (0.001, 1.0)    log - Affinity coefficient
#   Kadsdes (1/s):          (1e-6, 0.1)     log - Kinetic rate
#   bulk_density (kg/m3):   (1000, 2000)    linear
#
# Sediment (HYPE_HBVSED/MMF):
#   erosion_index:          (0.1, 2.0)      linear - Calibration coefficient
#   slope (m/m):            (0.001, 0.5)    linear - Land slope
#   slope_erosion_exponent: (0.5, 3.0)      linear - USLE ~1.4
#   precip_erosion_exponent:(0.5, 3.0)      linear - Erosivity ~1.5
#   cohesion (kPa):         (1.0, 30.0)     linear - Sandy to clay
#   erodibility (g/J):      (0.5, 10.0)     linear - MMF model
#   monthly factors:        (-0.5, 2.0)     linear - Seasonal adjustment
#
# Transport Dissolved:
#   dispersion_x (m2/s):    (0.01, 100.0)   log - Rivers 1-100 m2/s
#   dispersion_y (m2/s):    (0.001, 10.0)   log - 10x smaller than x
#   dispersion_z (m2/s):    (1e-5, 0.1)     log - Stratified to mixed
#   characteristic_length:  (10.0, 10000.0) log - Reach scale
#
# Lateral Exchange:
#   K_val (1/s):            (1e-14, 1e-6)   log - Very slow diffusive
#
# Source/Sink:
#   csv_scale (all/species):(0.1, 5.0)      linear - Load multiplier
#   cropland_TN (kg/ha/yr): (5.0, 80.0)     linear - Literature 10-50
#   cropland_TP (kg/ha/yr): (0.5, 10.0)     linear - Literature 0.5-5
#   forest_TN (kg/ha/yr):   (0.5, 10.0)     linear - Background 1-5
#   forest_TP (kg/ha/yr):   (0.05, 1.0)     linear - Background 0.1-0.5
#   urban_TN (kg/ha/yr):    (5.0, 50.0)     linear - Stormwater 5-30
#   urban_TP (kg/ha/yr):    (0.5, 5.0)      linear - Stormwater 0.5-3
#   precip_scaling_power:   (0.5, 2.5)      linear - Nonlinear erosion
#   Q10_biological:         (1.5, 4.0)      linear - Temp response
#   T_reference (degC):     (10.0, 25.0)    linear - Baseline temp
#
# ML Model Hyperparameters:
#   n_estimators:           (50, 500)       linear, int - Random Forest
#   max_depth:              (3, 20)         linear, int - Tree depth
#   learning_rate:          (0.01, 0.3)     log - Gradient boosting
# =============================================================================

calibration_parameters = [
    # =========================================================================
    # BGC Rate Constants (NATIVE_BGC_FLEX)
    # Path: ["CYCLING_FRAMEWORKS", "<framework>", "<rxn_num>", "parameter_values", "<param>"]
    # Scientific ranges from Chapra (2008) and Thomann & Mueller (1987)
    # =========================================================================
    {
        "name": "k_nitrification",
        "file_type": "bgc_json",
        "path": ["CYCLING_FRAMEWORKS", "N_cycle", "3", "parameter_values", "k"],
        "initial": 0.03,                    # Typical: 0.01-0.1/day
        "bounds": (0.001, 0.5),             # Chapra: up to 0.5/day in warm water
        "transform": "log"
    },
    {
        "name": "k_denitrification",
        "file_type": "bgc_json",
        "path": ["CYCLING_FRAMEWORKS", "N_cycle", "9", "parameter_values", "k"],
        "initial": 1e-6,                    # Highly variable
        "bounds": (1e-8, 0.01),             # SWAT: 0.001-0.01/day anoxic sediments
        "transform": "log"
    },
    {
        "name": "k_mineralization_active",
        "file_type": "bgc_json",
        "path": ["CYCLING_FRAMEWORKS", "N_cycle", "4", "parameter_values", "k"],
        "initial": 0.1,                     # Fast pool
        "bounds": (0.01, 1.0),              # Thomann & Mueller: 0.05-0.5/day
        "transform": "log"
    },

    # =========================================================================
    # PHREEQC Parameters (if using bgc_module_name = "PHREEQC")
    # Units: mg/L internally, converted to mol/kgw for PHREEQC
    # Uncomment and modify as needed for your geochemistry model
    # =========================================================================
    # {
    #     "name": "initial_Ca",
    #     "file_type": "phreeqc_pqi",
    #     "path": {"block": "SOLUTION", "species": "Ca"},
    #     "initial": 40.0,                  # mg/L - typical freshwater
    #     "bounds": (5.0, 500.0),           # Natural range
    #     "transform": "log"
    # },
    # {
    #     "name": "initial_NO3",
    #     "file_type": "phreeqc_pqi",
    #     "path": {"block": "SOLUTION", "species": "N(5)"},
    #     "initial": 2.0,                   # mg/L as N
    #     "bounds": (0.01, 50.0),           # Pristine to agricultural
    #     "transform": "log"
    # },
    # {
    #     "name": "initial_NH4",
    #     "file_type": "phreeqc_pqi",
    #     "path": {"block": "SOLUTION", "species": "N(-3)"},
    #     "initial": 0.1,                   # mg/L as N
    #     "bounds": (0.001, 10.0),          # Surface water typically <0.5
    #     "transform": "log"
    # },
    # {
    #     "name": "initial_PO4",
    #     "file_type": "phreeqc_pqi",
    #     "path": {"block": "SOLUTION", "species": "P"},
    #     "initial": 0.05,                  # mg/L as P
    #     "bounds": (0.001, 5.0),           # Oligotrophic to eutrophic
    #     "transform": "log"
    # },
    # {
    #     "name": "pCO2_log",
    #     "file_type": "phreeqc_pqi",
    #     "path": {"block": "EQUILIBRIUM_PHASES", "phase": "CO2(g)", "field": "si"},
    #     "initial": -3.5,                  # log10(atm) - atmospheric
    #     "bounds": (-4.5, -1.5),           # -4.5 (low) to -1.5 (organic-rich)
    #     "transform": "linear"
    # },
    # {
    #     "name": "pO2_log",
    #     "file_type": "phreeqc_pqi",
    #     "path": {"block": "EQUILIBRIUM_PHASES", "phase": "O2(g)", "field": "si"},
    #     "initial": -0.7,                  # log10(atm) - atmospheric O2
    #     "bounds": (-6.0, 0.0),            # Anoxic (-6) to supersaturated (0)
    #     "transform": "linear"
    # },
    # {
    #     "name": "kinetic_denitrif_rate",
    #     "file_type": "phreeqc_pqi",
    #     "path": {"block": "KINETICS", "reaction": "Denitrification", "param": "-parms", "index": 0},
    #     "initial": 1e-9,                  # mol/L/s
    #     "bounds": (1e-12, 1e-6),          # Microbially-mediated
    #     "transform": "log"
    # },

    # =========================================================================
    # Sorption Isotherm Parameters (Freundlich/Langmuir)
    # Sandy soils: lower Kfr; clay soils: higher Kfr
    # =========================================================================
    # {
    #     "name": "NH4_Kfr",
    #     "file_type": "sorption_json",
    #     "path": {"module": "FREUNDLICH", "species": "NH4-N", "param": "Kfr"},
    #     "initial": 1.0,                   # L/kg - moderate soil
    #     "bounds": (0.01, 100.0),          # Sandy 0.1-1, clay 5-50 L/kg
    #     "transform": "log"
    # },
    # {
    #     "name": "NH4_Nfr",
    #     "file_type": "sorption_json",
    #     "path": {"module": "FREUNDLICH", "species": "NH4-N", "param": "Nfr"},
    #     "initial": 0.7,                   # dimensionless
    #     "bounds": (0.3, 1.0),             # 1.0 = linear isotherm
    #     "transform": "linear"
    # },
    # {
    #     "name": "PO4_Kfr",
    #     "file_type": "sorption_json",
    #     "path": {"module": "FREUNDLICH", "species": "PO4-P", "param": "Kfr"},
    #     "initial": 10.0,                  # L/kg - P strongly sorbed
    #     "bounds": (0.1, 500.0),           # Depends on Fe/Al oxide content
    #     "transform": "log"
    # },
    # {
    #     "name": "NH4_qmax",
    #     "file_type": "sorption_json",
    #     "path": {"module": "LANGMUIR", "species": "NH4-N", "param": "qmax_mg_per_kg"},
    #     "initial": 200.0,                 # mg/kg
    #     "bounds": (10.0, 2000.0),         # Sandy 50-200, clay 500-2000
    #     "transform": "log"
    # },
    # {
    #     "name": "bulk_density",
    #     "file_type": "sorption_json",
    #     "path": {"module": "FREUNDLICH", "param": "BULK_DENSITY_KG_M3"},
    #     "initial": 1500.0,                # kg/m3
    #     "bounds": (1000.0, 2000.0),       # Sandy 1600-1800, clay 1200-1400
    #     "transform": "linear"
    # },

    # =========================================================================
    # Sediment Transport Parameters (HYPE_HBVSED / HYPE_MMF)
    # Based on HYPE model documentation and USLE
    # =========================================================================
    # {
    #     "name": "erosion_index",
    #     "file_type": "sediment_json",
    #     "path": {"module": "HYPE_HBVSED", "param": "EROSION_INDEX"},
    #     "initial": 0.5,                   # dimensionless
    #     "bounds": (0.1, 2.0),             # HYPE calibration range
    #     "transform": "linear"
    # },
    # {
    #     "name": "slope_erosion_exponent",
    #     "file_type": "sediment_json",
    #     "path": {"module": "HYPE_HBVSED", "param": "SLOPE_EROSION_FACTOR_EXPONENT"},
    #     "initial": 1.5,                   # dimensionless
    #     "bounds": (0.5, 3.0),             # USLE uses ~1.4
    #     "transform": "linear"
    # },
    # {
    #     "name": "precip_erosion_exponent",
    #     "file_type": "sediment_json",
    #     "path": {"module": "HYPE_HBVSED", "param": "PRECIP_EROSION_FACTOR_EXPONENT"},
    #     "initial": 1.5,                   # dimensionless
    #     "bounds": (0.5, 3.0),             # Erosivity exponent
    #     "transform": "linear"
    # },
    # {
    #     "name": "monthly_erosion_winter",
    #     "file_type": "sediment_json",
    #     "path": {"module": "HYPE_HBVSED", "param": "MONTHLY_EROSION_FACTOR", "index": 0},
    #     "initial": 0.0,                   # January factor
    #     "bounds": (-0.5, 2.0),            # Frozen ground = reduced erosion
    #     "transform": "linear"
    # },
    # # MMF module parameters
    # {
    #     "name": "cohesion",
    #     "file_type": "sediment_json",
    #     "path": {"module": "HYPE_MMF", "param": "COHESION"},
    #     "initial": 7.5,                   # kPa
    #     "bounds": (1.0, 30.0),            # Sandy 2-5, loam 5-15, clay 15-30
    #     "transform": "linear"
    # },
    # {
    #     "name": "erodibility",
    #     "file_type": "sediment_json",
    #     "path": {"module": "HYPE_MMF", "param": "ERODIBILITY"},
    #     "initial": 2.0,                   # g/J
    #     "bounds": (0.5, 10.0),            # MMF model
    #     "transform": "linear"
    # },

    # =========================================================================
    # Transport Parameters (Dispersion)
    # Scale-dependent; rivers typically 1-100 m2/s longitudinal
    # =========================================================================
    # {
    #     "name": "dispersion_x",
    #     "file_type": "transport_json",
    #     "path": {"param": "DISPERSION_X_M2_S"},
    #     "initial": 1.0,                   # m2/s
    #     "bounds": (0.01, 100.0),          # Rivers 1-100 m2/s
    #     "transform": "log"
    # },
    # {
    #     "name": "dispersion_y",
    #     "file_type": "transport_json",
    #     "path": {"param": "DISPERSION_Y_M2_S"},
    #     "initial": 0.1,                   # m2/s - typically 10x smaller
    #     "bounds": (0.001, 10.0),
    #     "transform": "log"
    # },
    # {
    #     "name": "dispersion_z",
    #     "file_type": "transport_json",
    #     "path": {"param": "DISPERSION_Z_M2_S"},
    #     "initial": 0.001,                 # m2/s - vertical mixing
    #     "bounds": (1e-5, 0.1),            # Stratified 1e-5, mixed 0.01-0.1
    #     "transform": "log"
    # },
    # {
    #     "name": "characteristic_length",
    #     "file_type": "transport_json",
    #     "path": {"param": "CHARACTERISTIC_LENGTH_M"},
    #     "initial": 100.0,                 # m
    #     "bounds": (10.0, 10000.0),        # Reach scale
    #     "transform": "log"
    # },

    # =========================================================================
    # Lateral Exchange Parameters (NATIVE_LE_BOUNDMIX)
    # Very slow diffusive exchange between compartments
    # =========================================================================
    # {
    #     "name": "K_val_surface_soil",
    #     "file_type": "lateral_exchange_json",
    #     "path": {"exchange_id": 1, "param": "K_val"},
    #     "initial": 1e-9,                  # 1/s
    #     "bounds": (1e-12, 1e-6),          # Slow diffusive exchange
    #     "transform": "log"
    # },
    # {
    #     "name": "K_val_soil_groundwater",
    #     "file_type": "lateral_exchange_json",
    #     "path": {"exchange_id": 2, "param": "K_val"},
    #     "initial": 1e-10,                 # 1/s - even slower
    #     "bounds": (1e-14, 1e-8),
    #     "transform": "log"
    # },

    # =========================================================================
    # Source/Sink: CSV Load Scaling
    # Global multiplier for loads from CSV files
    # =========================================================================
    {
        "name": "fertilizer_N_scale",
        "file_type": "ss_csv_scale",
        "path": {"species": "all"},         # Scale all species
        "initial": 1.0,                     # No scaling
        "bounds": (0.1, 5.0),               # Calibration range
        "transform": "linear"
    },
    # {
    #     "name": "TN_load_scale",
    #     "file_type": "ss_csv_scale",
    #     "path": {"species": "TN"},          # Scale only TN
    #     "initial": 1.0,
    #     "bounds": (0.1, 5.0),
    #     "transform": "linear"
    # },

    # =========================================================================
    # Source/Sink: Copernicus LULC Export Coefficients (static)
    # Literature-based export coefficients in kg/ha/yr
    # =========================================================================
    # {
    #     "name": "cropland_TN_coeff",
    #     "file_type": "ss_copernicus_static",
    #     "path": {"lulc_class": 10, "species": "TN"},
    #     "initial": 25.0,                  # kg/ha/yr
    #     "bounds": (5.0, 80.0),            # Literature: 10-50, fertilized up to 80
    #     "transform": "linear"
    # },
    # {
    #     "name": "cropland_TP_coeff",
    #     "file_type": "ss_copernicus_static",
    #     "path": {"lulc_class": 10, "species": "TP"},
    #     "initial": 2.0,                   # kg/ha/yr
    #     "bounds": (0.5, 10.0),            # Literature: 0.5-5, erosion-prone higher
    #     "transform": "linear"
    # },
    # {
    #     "name": "forest_TN_coeff",
    #     "file_type": "ss_copernicus_static",
    #     "path": {"lulc_class": 50, "species": "TN"},
    #     "initial": 2.0,                   # kg/ha/yr
    #     "bounds": (0.5, 10.0),            # Natural background: 1-5
    #     "transform": "linear"
    # },
    # {
    #     "name": "grassland_TN_coeff",
    #     "file_type": "ss_copernicus_static",
    #     "path": {"lulc_class": 30, "species": "TN"},
    #     "initial": 10.0,                  # kg/ha/yr
    #     "bounds": (2.0, 30.0),            # 5-20, grazed higher
    #     "transform": "linear"
    # },
    # {
    #     "name": "urban_TN_coeff",
    #     "file_type": "ss_copernicus_static",
    #     "path": {"lulc_class": 190, "species": "TN"},
    #     "initial": 15.0,                  # kg/ha/yr
    #     "bounds": (5.0, 50.0),            # Stormwater: 5-30
    #     "transform": "linear"
    # },

    # =========================================================================
    # Source/Sink: Climate Response Parameters (dynamic Copernicus)
    # Modulate loads based on precipitation and temperature
    # =========================================================================
    # {
    #     "name": "precip_scaling_power",
    #     "file_type": "ss_climate_param",
    #     "path": "ss_climate_precip_scaling_power",
    #     "initial": 1.0,                   # Linear relationship
    #     "bounds": (0.5, 2.5),             # Nonlinear erosion: 1.5-2.0
    #     "transform": "linear"
    # },
    # {
    #     "name": "Q10_biological",
    #     "file_type": "ss_climate_param",
    #     "path": "ss_climate_temp_q10",
    #     "initial": 2.0,                   # Standard Q10
    #     "bounds": (1.5, 4.0),             # Microbial 1.5-3.0
    #     "transform": "linear"
    # },
    # {
    #     "name": "T_reference",
    #     "file_type": "ss_climate_param",
    #     "path": "ss_climate_temp_reference_c",
    #     "initial": 15.0,                  # degC
    #     "bounds": (10.0, 25.0),           # Regional baseline
    #     "transform": "linear"
    # },

    # =========================================================================
    # Source/Sink: ML Model Hyperparameters
    # For machine learning-based load estimation
    # =========================================================================
    # {
    #     "name": "ml_n_estimators",
    #     "file_type": "ss_ml_param",
    #     "path": "ss_ml_n_estimators",
    #     "initial": 200,
    #     "bounds": (50, 500),
    #     "transform": "linear",
    #     "dtype": "int"
    # },
    # {
    #     "name": "ml_max_depth",
    #     "file_type": "ss_ml_param",
    #     "path": "ss_ml_max_depth",
    #     "initial": 8,
    #     "bounds": (3, 20),
    #     "transform": "linear",
    #     "dtype": "int"
    # },
    # {
    #     "name": "ml_learning_rate",
    #     "file_type": "ss_ml_param",
    #     "path": "ss_ml_learning_rate",
    #     "initial": 0.1,
    #     "bounds": (0.01, 0.3),
    #     "transform": "log"
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

def prepare_grqa_observations(config: dict, grqa_config: dict) -> str:
    """
    Download and process GRQA data for calibration.

    Returns path to the prepared observation CSV file.
    """
    if not GRQA_AVAILABLE:
        raise ImportError(
            "GRQA extraction tools not available. Install required dependencies:\n"
            "  pip install geopandas requests\n"
            "Or provide pre-formatted observation_data_path."
        )

    print("=" * 60)
    print("GRQA DATA EXTRACTION")
    print("=" * 60)

    # Set output directory
    output_dir = grqa_config.get("output_dir")
    if output_dir is None:
        output_dir = os.path.join(config["calibration_work_dir"], "observation_data")
    os.makedirs(output_dir, exist_ok=True)

    # Step 1: Extract GRQA stations and observations
    print("\n[1/2] Extracting GRQA stations and observations...")

    # Check for local data path
    local_data_path = grqa_config.get("local_data_path")
    if local_data_path:
        print(f"  Using local GRQA data: {local_data_path}")
    else:
        print("  Will download from Zenodo if needed")

    extractor = GRQACalibrationExtractor(
        river_network_shapefile=grqa_config["river_network_shapefile"],
        species_mapping=grqa_config["species_mapping"],
        output_dir=output_dir,
        local_data_path=local_data_path,
        bounding_box=grqa_config.get("bounding_box"),
        time_period=(grqa_config.get("start_date"), grqa_config.get("end_date"))
    )

    grqa_output = extractor.run_extraction()
    print(f"  Extracted {grqa_output.get('total_observations', 0)} observations "
          f"from {grqa_output.get('total_stations', 0)} stations")

    # Step 2: Convert to calibration format
    print("\n[2/2] Converting to calibration format...")
    converter = CalibrationObservationConverter(
        grqa_stations_path=os.path.join(output_dir, "grqa_stations.csv"),
        grqa_observations_path=os.path.join(output_dir, "grqa_observations.csv"),
        river_network_shapefile=grqa_config["river_network_shapefile"],
        reach_id_column=grqa_config.get("reach_id_column", "seg_id"),
        output_dir=output_dir
    )

    result = converter.convert(
        max_distance_m=grqa_config.get("max_station_distance_m", 500),
        min_observations=grqa_config.get("min_observations", 10),
        target_units=grqa_config.get("target_units"),
        accepted_quality_flags=grqa_config.get("accepted_quality_flags")
    )

    calibration_obs_path = result.get("calibration_observations_path")
    print(f"  Generated calibration observations: {calibration_obs_path}")
    print(f"  Total observations: {result.get('total_observations', 0)}")
    print(f"  Matched reaches: {result.get('matched_reaches', 0)}")
    print(f"  Species: {result.get('species_list', [])}")

    return calibration_obs_path


def prepare_copernicus_observations(config: dict, copernicus_config: dict) -> str:
    """
    Generate synthetic observations from Copernicus land cover and climate data.

    Returns path to the generated observation CSV file.
    """
    if not COPERNICUS_AVAILABLE:
        raise ImportError(
            "Copernicus observation tools not available. Install required dependencies:\n"
            "  pip install geopandas rasterio xarray\n"
            "Or provide pre-formatted observation_data_path or use GRQA extraction."
        )

    print("=" * 60)
    print("COPERNICUS OBSERVATION GENERATION")
    print("=" * 60)

    # Set output directory
    output_dir = copernicus_config.get("output_dir")
    if output_dir is None:
        output_dir = os.path.join(config["calibration_work_dir"], "observation_data")
    os.makedirs(output_dir, exist_ok=True)

    # Generate observations using Copernicus workflow
    print("\nGenerating synthetic observations from land cover data...")
    generator = CopernicusObservationGenerator(
        river_network_shapefile=copernicus_config["river_network_shapefile"],
        reach_id_column=copernicus_config.get("reach_id_column", "seg_id"),
        land_cover_path=copernicus_config["land_cover_path"],
        climate_data_path=copernicus_config.get("climate_data_path"),
        output_dir=output_dir
    )

    result = generator.generate(
        target_species=copernicus_config["target_species"],
        export_coefficients=copernicus_config["export_coefficients"],
        start_date=copernicus_config.get("start_date", "2000-01-01"),
        end_date=copernicus_config.get("end_date", "2020-12-31"),
        temporal_resolution=copernicus_config.get("temporal_resolution", "monthly"),
        climate_params=copernicus_config.get("climate_params"),
        uncertainty_fraction=copernicus_config.get("uncertainty_fraction", 0.3),
        add_noise=copernicus_config.get("add_noise", True),
        noise_cv=copernicus_config.get("noise_cv", 0.2)
    )

    calibration_obs_path = result.get("calibration_observations_path")
    print(f"  Generated calibration observations: {calibration_obs_path}")
    print(f"  Total observations: {result.get('total_observations', 0)}")
    print(f"  Reaches covered: {result.get('reaches_covered', 0)}")
    print(f"  Species: {result.get('species_list', [])}")
    print(f"  Time period: {result.get('start_date')} to {result.get('end_date')}")

    return calibration_obs_path


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="OpenWQ Calibration Framework")
    parser.add_argument("--resume", action="store_true",
                       help="Resume from checkpoint")
    parser.add_argument("--sensitivity-only", action="store_true",
                       help="Run only sensitivity analysis")
    parser.add_argument("--dry-run", action="store_true",
                       help="Validate configuration without running")
    parser.add_argument("--extract-grqa-only", action="store_true",
                       help="Only extract GRQA data without running calibration")
    parser.add_argument("--generate-copernicus-only", action="store_true",
                       help="Only generate Copernicus observations without running calibration")
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

        # Observation data settings
        "observation_data_source": observation_data_source,
        "use_grqa": use_grqa,  # Deprecated, kept for backwards compatibility
        "grqa_config": grqa_config,
        "copernicus_config": copernicus_config,
    }

    # ==========================================================================
    # Observation Data Preparation (based on observation_data_source)
    # ==========================================================================

    # Determine effective source (handle backwards compatibility with use_grqa)
    effective_source = observation_data_source
    if use_grqa and observation_data_source == "csv":
        effective_source = "grqa"  # Backwards compatibility

    if effective_source == "grqa" or args.extract_grqa_only:
        print("Preparing observation data from GRQA database...")
        try:
            obs_path = prepare_grqa_observations(config, grqa_config)
            config["observation_data_path"] = obs_path
            print(f"\nObservation data ready: {obs_path}")
        except Exception as e:
            print(f"\nERROR: GRQA extraction failed: {e}")
            if args.extract_grqa_only:
                sys.exit(1)
            else:
                print("Falling back to configured observation_data_path...")

        if args.extract_grqa_only:
            print("\n" + "=" * 60)
            print("GRQA EXTRACTION COMPLETE")
            print("=" * 60)
            print(f"Observation data saved to: {config['observation_data_path']}")
            print("\nTo run calibration, use: python your_config.py")
            sys.exit(0)

    elif effective_source == "copernicus" or args.generate_copernicus_only:
        print("Generating synthetic observations from Copernicus data...")
        try:
            obs_path = prepare_copernicus_observations(config, copernicus_config)
            config["observation_data_path"] = obs_path
            print(f"\nObservation data ready: {obs_path}")
        except Exception as e:
            print(f"\nERROR: Copernicus observation generation failed: {e}")
            if args.generate_copernicus_only:
                sys.exit(1)
            else:
                print("Falling back to configured observation_data_path...")

        if args.generate_copernicus_only:
            print("\n" + "=" * 60)
            print("COPERNICUS OBSERVATION GENERATION COMPLETE")
            print("=" * 60)
            print(f"Observation data saved to: {config['observation_data_path']}")
            print("\nTo run calibration, use: python your_config.py")
            sys.exit(0)

    elif effective_source == "csv":
        print(f"Using pre-formatted observation file: {observation_data_path}")

    else:
        print(f"WARNING: Unknown observation_data_source '{effective_source}'. "
              f"Valid options: 'csv', 'grqa', 'copernicus'")
        print(f"Falling back to: {observation_data_path}")

    if args.dry_run:
        print("=" * 60)
        print("DRY RUN - Configuration Validation")
        print("=" * 60)
        print(f"Parameters to calibrate: {len(calibration_parameters)}")

        # Validate each parameter
        all_warnings = []
        for p in calibration_parameters:
            print(f"  - {p['name']}: {p['bounds']} (transform={p.get('transform', 'linear')})")
            warnings = validate_parameter_bounds(p)
            all_warnings.extend(warnings)

        print(f"\nAlgorithm: {algorithm}")
        print(f"Max evaluations: {max_evaluations}")
        print(f"Container runtime: {container_runtime}")
        print(f"Objective function: {objective_function}")

        if all_warnings:
            print("\n" + "=" * 60)
            print("WARNINGS:")
            print("=" * 60)
            for w in all_warnings:
                print(f"  ! {w}")
            print("\nConsider fixing these issues before running calibration.")
        else:
            print("\nConfiguration is valid. All parameters have reasonable bounds.")

        # Show available parameter reference
        print("\n" + "=" * 60)
        print("PARAMETER REFERENCE")
        print("=" * 60)
        print("See calibration_lib/parameter_defaults.py for:")
        print("  - BGC_PARAMETERS: 20+ biogeochemical rate constants")
        print("  - PHREEQC_PARAMETERS: 20+ geochemistry parameters")
        print("  - SORPTION_PARAMETERS: Freundlich/Langmuir isotherm params")
        print("  - SEDIMENT_PARAMETERS: HYPE_HBVSED/MMF erosion params")
        print("  - TRANSPORT_PARAMETERS: Dispersion coefficients")
        print("  - LATERAL_EXCHANGE_PARAMETERS: Compartment exchange rates")
        print("  - SOURCE_SINK_PARAMETERS: LULC export coefficients, climate params")
        print("\nUse create_parameter_entry() helper for automatic bounds lookup.")

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
