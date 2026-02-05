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

import sys
import os
sys.path.insert(0, 'config_support_lib')

import Gen_Input_Driver as gJSON_lib
from Gen_Input_Driver import uniform_param

# ============================================================================
# OpenWQ Configuration File
# ============================================================================
#
# SECTIONS:
#   1. General Information
#   2. Computational Settings
#   3. Initial Conditions
#   4. Modules:
#      a. Biogeochemistry  (BGC)   - NATIVE_BGC_FLEX or PHREEQC
#      b. Transport         (TD)   - NATIVE_TD_ADV, OPENWQ_NATIVE_TD_ADVDISP, or NONE
#      c. Lateral Exchange  (LE)   - NATIVE_LE_BOUNDMIX or NONE
#      d. Sediment Transport (TS)  - HYPE_MMF, HYPE_HBVSED, or NONE
#      e. Sorption Isotherm (SI)   - FREUNDLICH, LANGMUIR, or NONE
#   5. Source/Sink Configuration
#   6. External Water Fluxes
#   7. Output Settings
#
# Full documentation: https://openwq.readthedocs.io
# ============================================================================


# ======================== 1. GENERAL INFORMATION ============================

project_name = "Demonstration"
geographical_location = "NA"
authors = "Diogo Costa"
date = "January, 2026"
comment = "This a demonstration input file setup"
hostmodel = "mizuroute"   # "mizuroute" or "summa"

# Directory to save generated input files
dir2save_input_files = "/Users/diogocosta/Documents/openwq_code/6_mizuroute_cslm_openwq/test_case_2/"

# Docker: set True when running inside the Docker container.
# Paths are auto-corrected from host to container using docker-compose.yml volume mount.
running_on_docker = True


# ====================== 2. COMPUTATIONAL SETTINGS ==========================

solver = "FORWARD_EULER"   # "FORWARD_EULER" or "SUNDIALS"
run_mode_debug = True      # True = print solver derivatives
use_num_threads = 4        # number of threads (or "all")


# ====================== 3. INITIAL CONDITIONS ==============================

ic_all_value = 2           # initial concentration/mass value
ic_all_units = "mg/l"      # units of ic_all_value


# ====================== 4a. BIOGEOCHEMISTRY MODULE ==========================
# Options: "NATIVE_BGC_FLEX" (user-defined kinetic expressions)
#          "PHREEQC"         (full geochemical equilibrium engine)

bgc_module_name = "NATIVE_BGC_FLEX"

# -- NATIVE_BGC_FLEX --
path2selected_NATIVE_BGC_FLEX_framework = "config_support_lib/examples_BGC_frameworks/openWQ_Ncycling_example.json"

# -- PHREEQC --
# See wiki: https://openwq.readthedocs.io (Modules > PHREEQC)
phreeqc_input_filepath = "config_support_lib/examples_PHREEQC/phreeqc_river.pqi"
phreeqc_database_filepath = "config_support_lib/examples_PHREEQC/phreeqc.dat"
phreeqc_mobile_species = ["H", "O", "Charge", "Ca"]   # must match PHREEQC components
phreeqc_component_h2o = True
phreeqc_temperature_mapping = {"RIVER_NETWORK_REACHES": "air_temperature"}  # or None
phreeqc_pressure_mapping = None
phreeqc_chemical_species_names = []  # empty = use .pqi defaults; run once with [] to see component names in log


# ====================== 4b. TRANSPORT DISSOLVED MODULE ======================
# Options: "NATIVE_TD_ADV", "OPENWQ_NATIVE_TD_ADVDISP", "NONE"

td_module_name = "OPENWQ_NATIVE_TD_ADVDISP"
td_module_dispersion_xyz = [0.3, 0.3, 0.3]      # [Dx, Dy, Dz] in m2/s
td_module_characteristic_length_m = 100.0         # distance between cell centers [m]


# ====================== 4c. LATERAL EXCHANGE MODULE =========================
# Options: "NATIVE_LE_BOUNDMIX", "NONE"

le_module_name = "NATIVE_LE_BOUNDMIX"
le_module_config = [
    {
        "direction": "z",
        "upper_compartment": "RUNOFF",
        "lower_compartment": "ILAYERVOLFRACWAT_SOIL",
        "K_val": 0.000000001
    }
]


# ====================== 4d. SEDIMENT TRANSPORT MODULE =======================
# Options: "HYPE_MMF", "HYPE_HBVSED", "NONE"
# See wiki: https://openwq.readthedocs.io (Modules > Transport Sediments)

ts_module_name = "HYPE_HBVSED"
ts_sediment_compartment = "RIVER_NETWORK_REACHES"
ts_transport_compartment = "RIVER_NETWORK_REACHES"
ts_direction = "z"
ts_erosion_inhibit_compartment = "RIVER_NETWORK_REACHES"
ts_data_format = "JSON"

# -- HYPE_MMF parameters (used only if ts_module_name = "HYPE_MMF") --
# Defaults applied to all cells; override per-cell in ts_mmf_parameters below
ts_mmf_defaults = {
    "COHESION": 7.5,  "ERODIBILITY": 2.0, "SREROEXP": 1.2,
    "CROPCOVER": 0.5, "GROUNDCOVER": 0.5, "SLOPE": 0.4,     # CROPCOVER + GROUNDCOVER must = 1.0
    "TRANSPORT_FACTOR_1": 0.2, "TRANSPORT_FACTOR_2": 0.5
}

# Spatially-varying overrides (leave {} to use defaults everywhere)
ts_mmf_parameters = {k: uniform_param(v) for k, v in ts_mmf_defaults.items()}

# -- HYPE_HBVSED parameters (used only if ts_module_name = "HYPE_HBVSED") --
ts_hbvsed_defaults = {
    "SLOPE": 0.4, "EROSION_INDEX": 0.4,
    "SOIL_EROSION_FACTOR_LAND_DEPENDENCE": 0.4, "SOIL_EROSION_FACTOR_SOIL_DEPENDENCE": 0.4,
    "SLOPE_EROSION_FACTOR_EXPONENT": 1.5, "PRECIP_EROSION_FACTOR_EXPONENT": 1.5,
    "PARAM_SCALING_EROSION_INDEX": 0.5
}
# Spatially-varying overrides (leave {} to use defaults everywhere)
ts_hbvsed_parameters = {k: uniform_param(v) for k, v in ts_hbvsed_defaults.items()}

# Monthly erosion factor (12 values, Jan-Dec). erodmonth = 1.0 + factor[month]
ts_hbvsed_monthly_erosion_factor = [0.0] * 12


# ====================== 4e. SORPTION ISOTHERM MODULE =======================
# Options: "FREUNDLICH", "LANGMUIR", "NONE"
# WARNING: Only compatible with NATIVE_BGC_FLEX. Set "NONE" when using PHREEQC
#          (PHREEQC handles sorption via SURFACE/EXCHANGE blocks in the .pqi file).

si_module_name = "NONE"
si_sediment_compartment = "RIVER_NETWORK_REACHES"
si_bulk_density_kg_m3 = 1500.0    # [kg/m3]
si_layer_thickness_m = 1.0        # [m]

# Per-species params (set None when si_module_name = "NONE")
# FREUNDLICH: {"species": {"Kfr": float, "Nfr": float, "Kadsdes_1_per_s": float}}
# LANGMUIR:   {"species": {"qmax_mg_per_kg": float, "KL_L_per_mg": float, "Kadsdes_1_per_s": float}}
si_species_params = None


# ====================== 5. SOURCE/SINK CONFIGURATION =======================
# Options:
#   "load_from_csv"                             - CSV time series
#   "using_copernicus_lulc_with_static_coeff"   - Copernicus LULC + static export coefficients
#   "using_copernicus_lulc_with_dynamic_coeff"  - Copernicus LULC + climate-modulated coefficients
#   "ml_model"                                  - Train ML model from monitoring data
#   "none"                                      - No source/sink
# Note: When using PHREEQC, set "none" unless CSV references valid PHREEQC components.

ss_method = "load_from_csv"

# -- Common metadata --
ss_metadata_source = "Just for demonstration"
ss_metadata_comment = "Leave any comments needed for future reference"

# -- CSV method settings (ss_method = "load_from_csv") --
ss_method_csv_config = [
    {
        "Chemical_name": "N_ORG_active",
        "Compartment_name": "RIVER_NETWORK_REACHES",
        "Type": "source",
        "Units": "kg",
        "Filepath": "/Users/diogocosta/Documents/openwq_code/6_mizuroute_cslm_openwq/route/build/openwq/openwq/supporting_scripts/Model_Config/config_support_lib/examples_SS_in_cvs/N_ORG_active_fertilizer.csv",
        "Delimiter": ","
    },
    {
        "Chemical_name": "NH4-N",
        "Compartment_name": "RIVER_NETWORK_REACHES",
        "Type": "source",
        "Units": "kg",
        "Filepath": "/Users/diogocosta/Documents/openwq_code/6_mizuroute_cslm_openwq/route/build/openwq/openwq/supporting_scripts/Model_Config/config_support_lib/examples_SS_in_cvs/NH4_fertilizer.csv",
        "Delimiter": ","
    }
]

# -- Copernicus LULC settings (ss_method = "using_copernicus_lulc_with_*") --
# See wiki: https://openwq.readthedocs.io (Source/Sink > Copernicus LULC)
ss_method_copernicus_basin_info = {
    'path_to_shp': '/Users/diogocosta/Documents/openwq_code/6_mizuroute_cslm_openwq/test_case/mizuroute_in/shapefiles/finalcat_info_v1-0.shp',
    'mapping_key': 'SubId'
}
ss_method_copernicus_openwq_h5_results_file_example_for_mapping_key = {
    'path_to_file': '/Users/diogocosta/Documents/openwq_code/6_mizuroute_cslm_openwq/test_case/openwq_out/HDF5/RIVER_NETWORK_REACHES@N_ORG_ACTIVE#MG|L-main.h5',
    'mapping_key': 'reachID'
}
ss_method_copernicus_compartment_name_for_load = "RIVER_NETWORK_REACHES"
ss_method_copernicus_nc_lc_dir = '/Users/diogocosta/Documents/ESACCI-LC/'   # Copernicus LULC NetCDF directory
ss_method_copernicus_period = [1993, 1994]                                   # [year_start, year_end]
ss_method_copernicus_default_loads_bool = True                               # True = use built-in export coefficients
# Custom load coefficients (only when ss_method_copernicus_default_loads_bool = False)
# Format: {lulc_class_code: {"species": load_kg_ha_yr, ...}}
ss_method_copernicus_optional_custom_annual_load_coeffs_per_lulc_class = {
    10: {'TN': 20.0, 'TP': 3.0, 'NH4': 2.5},    # cropland
    50: {'TN': 2.0,  'TP': 0.1, 'NH4': 0.2},     # forest
    130: {'TN': 5.0, 'TP': 0.15, 'NH4': 0.5}     # grassland
}
ss_method_copernicus_annual_to_seasonal_loads_method = 'uniform'   # 'uniform' or 'seasonal'

# -- Climate-adjusted settings (ss_method = "using_copernicus_lulc_with_dynamic_coeff") --
# Monthly weight: w_m = P_m^alpha * Q10^((T_m - T_ref) / 10)
ss_climate_data = {
    1993: {
        'precip_mm': [80, 70, 90, 100, 110, 120, 80, 60, 70, 90, 100, 85],
        'temp_c':    [2,   3,  7,  12,  17,  22, 25, 24, 19, 13,   7,  3]
    },
    1994: {
        'precip_mm': [85, 75, 95, 105, 115, 125, 85, 65, 75, 95, 105, 90],
        'temp_c':    [1,   4,  8,  13,  18,  23, 26, 25, 20, 14,   8,  2]
    }
}
ss_climate_precip_scaling_power = 1.0    # precipitation exponent (1.0 = linear)
ss_climate_temp_q10 = 2.0               # biological rate doubling per 10 C
ss_climate_temp_reference_c = 15.0       # reference temperature [C]

# -- ML model settings (ss_method = "ml_model") --
# Requires: pip install scikit-learn xgboost
ss_ml_training_data_csv = None         # path to CSV (columns: date, discharge_m3s, precip_mm, temp_c, <species>)
ss_ml_model_type = "xgboost"           # "xgboost" or "random_forest"
ss_ml_target_species = None            # e.g. ["NO3-N", "TP"] (None = auto-detect)
ss_ml_feature_columns = None           # e.g. ["discharge_m3s", "precip_mm"] (None = auto-detect)
ss_ml_n_estimators = 200               # number of trees
ss_ml_max_depth = 6                    # max tree depth


# ====================== 6. EXTERNAL WATER FLUXES ===========================
# Options: "fixed_value" or "none"
# Note: When using PHREEQC, set "none" unless species name is a valid PHREEQC component.

ewf_method = "none"
ewf_method_fixedval_comment = "External model: summa"
ewf_method_fixedval_source = "Just for demonstration purposes"
ewf_method_fixedval_chem_name = "NO3-N"
ewf_method_fixedval_value = 2.5
ewf_method_fixedval_units = "mg/l"
ewf_method_fixedval_external_inputflux_name = "SUMMA_RUNOFF"


# ====================== 7. OUTPUT SETTINGS ==================================

output_format = "HDF5"              # "HDF5" or "CSV" (HDF5 recommended)
chemical_species = ["NO3-N", "NH4-N", "N_ORG_fresh", "N_ORG_stable"]
units = "MG/L"
no_water_conc_flag = -9999          # flag for cells without water
export_sediment = True
timestep = [1, "hour"]              # printing timestep [value, unit]
compartments_and_cells = {
    "RIVER_NETWORK_REACHES": {
        "1": ["all", "all", "all"]
    }
}


# ============================================================================
# DO NOT EDIT BELOW THIS LINE -- Collects settings and generates input files
# ============================================================================

# Gather all configuration variables and call the generator
_config = {k: v for k, v in locals().items()
           if not k.startswith('_') and k not in ('sys', 'os', 'gJSON_lib', 'uniform_param')}
gJSON_lib.Gen_Input_Driver(**_config)
