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

"""
Example: Generate all OpenWQ input files
"""

################
# General information
################
project_name = "Demonstration"
geographical_location = "NA"
authors = "Diogo Costa"
date = "January, 2026"
comment = "This a demonstration input file setup"

hostmodel = "mizuroute" # mizuroute or summa

# Directory to save the master file and the folder with the other input files,
dir2save_input_files = "/Users/diogocosta/Documents/openwq_code/6_mizuroute_cslm_openwq/test_case_2/"

# Docker support: set to True when running inside the Docker container.
# When True, all file paths (dir2save_input_files, CSV paths, shapefiles, etc.)
# will be automatically corrected from host paths to container paths using
# the volume mount defined in docker-compose.yml.
# The docker-compose.yml volume mount '../../../../../../:/code:Z' maps:
#   Host: /Users/diogocosta/Documents -> Container: /code
running_on_docker = True

################
# Computational settings
################
solver = "FORWARD_EULER"           # "FORWARD_EULER" or "SUNDIALS"
run_mode_debug = True   # "True" will print the solver derivatives that lead to the concentrations
use_num_threads = 4     # number of threads, accepts "all"

################
# General configuration
################
ic_all_value = 2          # initial conditions value, concentration or mass
ic_all_units = "mg/l"     # set the units of ic_all_value

################
### MODULES
################

"""
--------------------------------------
Biogeochemistry module
Options:
    NATIVE_BGC_FLEX   - User-defined kinetic expressions (simple, fast)
    PHREEQC           - Full geochemical equilibrium engine (thermodynamic)
# in this case, "NONE" is not an option
--------------------------------------
"""
bgc_module_name = "NATIVE_BGC_FLEX"    # Change to "PHREEQC" to use the PHREEQC geochemical engine

# ---- NATIVE_BGC_FLEX settings ----
# Path to the BGC cycling framework JSON file (only used when bgc_module_name = "NATIVE_BGC_FLEX")
path2selected_NATIVE_BGC_FLEX_framework = "config_support_lib/examples_BGC_frameworks/openWQ_Ncycling_example.json"

# ---- PHREEQC settings ----
# Only used when bgc_module_name = "PHREEQC"
# Path to PHREEQC input file (.pqi) defining solutions, reactions, kinetics
phreeqc_input_filepath = "config_support_lib/examples_PHREEQC/phreeqc_river.pqi"
# Path to PHREEQC thermodynamic database (.dat)
# Download full database from: https://www.usgs.gov/software/phreeqc-version-3
phreeqc_database_filepath = "config_support_lib/examples_PHREEQC/phreeqc.dat"
# Mobile species indices (1-indexed). These species will be transported.
# The species order is determined by PhreeqcRM::GetComponents() after parsing
# the database + input file. Run once to see the component list in the log.
# IMPORTANT: Indices must not exceed the number of components found by PHREEQC.
phreeqc_mobile_species = [1, 2, 3, 4]
# Whether to include H2O as a component in PHREEQC (default True)
phreeqc_component_h2o = True
# Temperature dependency mapping: compartment name -> dependency variable name
# Set to None if not using temperature-dependent reactions
phreeqc_temperature_mapping = {
    "RIVER_NETWORK_REACHES": "air_temperature"
}
# Pressure dependency mapping: compartment name -> dependency variable name
# Set to None if not using pressure-dependent reactions
phreeqc_pressure_mapping = None
# Chemical species names for initial conditions in the config file.
# These should match the component names discovered by PHREEQC from the database.
# IMPORTANT: Run once with an empty list [] to see the actual component names in the log,
# then set this list accordingly. The number of species must match what PHREEQC finds.
# Common PHREEQC components (with full database): H, O, Charge, Ca, Mg, Na, K, N, C, Cl, S
# With the minimal example database, fewer components will be found.
# Set to empty list [] to skip initial conditions (PHREEQC will use .pqi defaults).
phreeqc_chemical_species_names = []

"""
--------------------------------------
Transport Dissolved module
Options:
    NATIVE_TD_ADV,
    OPENWQ_NATIVE_TD_ADVDISP,
    "NONE"
--------------------------------------
"""
td_module_name = "NATIVE_TD_ADV"
# if td_module_name=OPENWQ_NATIVE_TD_ADVDISP
td_module_dispersion_xyz = [0.3, 0.3, 0.3]


"""
--------------------------------------
Lateral Exchange module
Options:
    NATIVE_LE_BOUNDMIX
    NONE
--------------------------------------
"""
le_module_name = "NATIVE_LE_BOUNDMIX"
le_module_config=[
        {
            "direction": "z",
            "upper_compartment": "RUNOFF",
            "lower_compartment": "ILAYERVOLFRACWAT_SOIL",
            "K_val": 0.000000001
        }]

"""
--------------------------------------
Transport Sediments module
Options:
    HYPE_MMF
    HYPE_HBVSED
    NONE
--------------------------------------
"""
ts_module_name = "HYPE_HBVSED"
ts_sediment_compartment = "RIVER_NETWORK_REACHES"
ts_transport_compartment = "RIVER_NETWORK_REACHES"

# ---- Common TS configuration ----
ts_direction = "z"                                          # Direction of exchange: "x", "y", or "z"
ts_erosion_inhibit_compartment = "RIVER_NETWORK_REACHES"    # Compartment that inhibits erosion (e.g., snow layer)
ts_data_format = "JSON"                                     # Data format for parameters: "JSON" or "ASCII"

# ---- HYPE_MMF settings (Morgan-Morgan-Finney erosion model) ----
# Only used when ts_module_name = "HYPE_MMF"
# Default parameter values (applied to all cells unless overridden in ts_mmf_parameters)
# IMPORTANT: CROPCOVER + GROUNDCOVER must sum to 1.0
ts_mmf_defaults = {
    "COHESION": 7.5,              # Soil cohesion (kPa)
    "ERODIBILITY": 2.0,           # Soil erodibility (g/J)
    "SREROEXP": 1.2,              # Surface runoff erosion exponent (-)
    "CROPCOVER": 0.5,             # Crop cover fraction (0-1)
    "GROUNDCOVER": 0.5,           # Ground cover fraction (0-1)
    "SLOPE": 0.4,                 # Basin slope
    "TRANSPORT_FACTOR_1": 0.2,    # Transport capacity factor 1 (-)
    "TRANSPORT_FACTOR_2": 0.5     # Transport capacity factor 2 (-)
}
# Spatially-varying parameters (cell-specific overrides)
# Format: { "PARAM_NAME": { "0": ["IX","IY","IZ","VALUE"], "1": [ix, iy, iz, value], ... } }
# Use "ALL" for ix/iy/iz to apply to all cells in that dimension
# Row "0" is the header; rows "1","2",... are data entries
# Leave empty {} to use defaults everywhere (warnings shown, but defaults still apply)
# Example with spatially-varying overrides:
ts_mmf_parameters = {
    "COHESION": {
        "0": ["IX", "IY", "IZ", "VALUE"],
        "1": ["ALL", 1, 1, 7.5]
    },
    "ERODIBILITY": {
        "0": ["IX", "IY", "IZ", "VALUE"],
        "1": ["ALL", 1, 1, 2.0]
    },
    "SREROEXP": {
        "0": ["IX", "IY", "IZ", "VALUE"],
        "1": ["ALL", 1, 1, 1.2]
    },
    "CROPCOVER": {
        "0": ["IX", "IY", "IZ", "VALUE"],
        "1": ["ALL", 1, 1, 0.5]
    },
    "GROUNDCOVER": {
        "0": ["IX", "IY", "IZ", "VALUE"],
        "1": ["ALL", 1, 1, 0.5]
    },
    "SLOPE": {
        "0": ["IX", "IY", "IZ", "VALUE"],
        "1": ["ALL", 1, 1, 0.4]
    },
    "TRANSPORT_FACTOR_1": {
        "0": ["IX", "IY", "IZ", "VALUE"],
        "1": ["ALL", 1, 1, 0.2]
    },
    "TRANSPORT_FACTOR_2": {
        "0": ["IX", "IY", "IZ", "VALUE"],
        "1": ["ALL", 1, 1, 0.5]
    }
}

# ---- HYPE_HBVSED settings (HBV-sed erosion model) ----
# Only used when ts_module_name = "HYPE_HBVSED"
# Default parameter values (applied to all cells unless overridden in ts_hbvsed_parameters)
ts_hbvsed_defaults = {
    "SLOPE": 0.4,                                   # Basin slope
    "EROSION_INDEX": 0.4,                           # Erosion index
    "SOIL_EROSION_FACTOR_LAND_DEPENDENCE": 0.4,     # Land use dependent factor (lusepar)
    "SOIL_EROSION_FACTOR_SOIL_DEPENDENCE": 0.4,     # Soil dependent factor (soilpar)
    "SLOPE_EROSION_FACTOR_EXPONENT": 1.5,           # Slope factor exponent (slopepar)
    "PRECIP_EROSION_FACTOR_EXPONENT": 1.5,          # Precipitation exponent (precexppar)
    "PARAM_SCALING_EROSION_INDEX": 0.5              # Erosion index scaling parameter (eroindexpar)
}
# Spatially-varying parameters (same format as HYPE_MMF)
# Example with spatially-varying overrides:
ts_hbvsed_parameters = {
    "SLOPE": {
        "0": ["IX", "IY", "IZ", "VALUE"],
        "1": ["ALL", 1, 1, 0.4]
    },
    "EROSION_INDEX": {
        "0": ["IX", "IY", "IZ", "VALUE"],
        "1": ["ALL", 1, 1, 0.4]
    },
    "SOIL_EROSION_FACTOR_LAND_DEPENDENCE": {
        "0": ["IX", "IY", "IZ", "VALUE"],
        "1": ["ALL", 1, 1, 0.4]
    },
    "SOIL_EROSION_FACTOR_SOIL_DEPENDENCE": {
        "0": ["IX", "IY", "IZ", "VALUE"],
        "1": ["ALL", 1, 1, 0.4]
    },
    "SLOPE_EROSION_FACTOR_EXPONENT": {
        "0": ["IX", "IY", "IZ", "VALUE"],
        "1": ["ALL", 1, 1, 1.5]
    },
    "PRECIP_EROSION_FACTOR_EXPONENT": {
        "0": ["IX", "IY", "IZ", "VALUE"],
        "1": ["ALL", 1, 1, 1.5]
    },
    "PARAM_SCALING_EROSION_INDEX": {
        "0": ["IX", "IY", "IZ", "VALUE"],
        "1": ["ALL", 1, 1, 0.5]
    }
}
# Monthly erosion factor (12 values, Jan-Dec)
# Used as: erodmonth = 1.0 + MONTHLY_EROSION_FACTOR[month]
# Default: all zeros (erodmonth = 1.0 year-round)
ts_hbvsed_monthly_erosion_factor = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

"""
Sorption Isotherm module
Options:
    FREUNDLICH
    LANGMUIR
    NONE
"""
si_module_name = "NONE"
si_sediment_compartment = "SUMMA_RUNOFF"
# TODO - SETUP FOR FREUNDLICH AND LANGMUIR OPTIONS

################
# Sink sources
################
ss_method = "load_from_csv" # "load_from_csv", "using_copernicus_lulc", or "none"
                    # Note: When using PHREEQC, set to "none" unless your CSV references
                    # valid PHREEQC component names (not NATIVE_BGC_FLEX species like NO3-N)
################
# if ss_method = "load_from_csv"
ss_metadata_source = "Just for demonstration"
ss_metadata_comment = "Leave any comments needed for future reference"
ss_method_csv_config = source_sink_configs=[
        {
            "Chemical_name": "N_ORG_active",
            "Compartment_name": "RIVER_NETWORK_REACHES",
            "Type": "source",
            "Units": "kg",
            "Filepath": "/Users/diogocosta/Documents/openwq_code/6_mizuroute_cslm_openwq/route/build/openwq/openwq/supporting_scripts/Model_Config/config_support_lib/examples_SS_in_cvs/N_ORG_active_fertilizer.csv",
            "Delimiter": ",",
            "Number_of_header_rows": 3,
            "Header_key_row": 3
        },
        {
            "Chemical_name": "NH4-N",
            "Compartment_name": "RIVER_NETWORK_REACHES",
            "Type": "source",
            "Units": "kg",
            "Filepath": "/Users/diogocosta/Documents/openwq_code/6_mizuroute_cslm_openwq/route/build/openwq/openwq/supporting_scripts/Model_Config/config_support_lib/examples_SS_in_cvs/NH4_fertilizer.csv",
            "Delimiter": ",",
            "Number_of_header_rows": 3,
            "Header_key_row": 3
        }
    ]
################
# if ss_method = load_from_copernicus
ss_method_copernicus_basin_info = {
    'path_to_shp': '/Users/diogocosta/Documents/openwq_code/6_mizuroute_cslm_openwq/test_case/mizuroute_in/shapefiles/finalcat_info_v1-0.shp',
    'mapping_key': 'SubId'
}
ss_method_copernicus_openwq_h5_results_file_example_for_mapping_key = {
    'path_to_file': '/Users/diogocosta/Documents/openwq_code/6_mizuroute_cslm_openwq/test_case/openwq_out/HDF5/RIVER_NETWORK_REACHES@N_ORG_ACTIVE#MG|L-main.h5',
    'mapping_key': 'reachID'
}
ss_method_copernicus_compartment_name_for_load = "RIVER_NETWORK_REACHES"
# download nc with lulc from ESA CCI LC (copernicus: https://cds.climate.copernicus.eu/datasets/satellite-land-cover?tab=overview)
# and point the location to that folder below in 'ss_method_copernicus_nc_lc_dir'
ss_method_copernicus_nc_lc_dir = '/Users/diogocosta/Documents/ESACCI-LC/'
ss_method_copernicus_period = [1993, 1994]
ss_method_copernicus_default_loads_bool = True
# if ss_method_copernicus_default_loads_bool=False, then need to set ss_method_copernicus_optional_custom_annual_load_coeffs_per_lulc_class
# key is the class code
# values are: chemical_species: load (kg/ha/year)
ss_method_copernicus_optional_custom_annual_load_coeffs_per_lulc_class = {
    10: {'TN': 20.0, 'TP': 3.0, 'NH4': 2.5},  # Custom cropland values
    50: {'TN': 2.0, 'TP': 0.1, 'NH4': 0.2},   # Custom forest values
    130: {'TN': 5.0, 'TP': 0.15, 'NH4': 0.5}  # Custom grassland values
}
ss_method_copernicus_annual_to_seasonal_loads_method = 'uniform' # options: 'uniform' or 'seasonal'

################
# External water fluxes
################
ewf_method = "none" # "fixed_value" or "none"
                     # Note: When using PHREEQC, set to "none" unless the chemical name
                     # matches a valid PHREEQC component (not NATIVE_BGC_FLEX species like NO3-N)
ewf_method_fixedval_comment = "External model: summa"
ewf_method_fixedval_source = "Just for demonstration purposes"
ewf_method_fixedval_chem_name = "NO3-N"
ewf_method_fixedval_value = 2.5
ewf_method_fixedval_units = "mg/l"
ewf_method_fixedval_external_inputflux_name = "SUMMA_RUNOFF"

################
# Output settings
################
output_format = "HDF5"              # HDF5 or CSV (but HDF5 is recommended as can be used with other scripts
chemical_species = [1, 2, 3, 4]     # Chemical species indices to export (must not exceed num components)
units = "MG/L"                      # Units (concentrations or mass)
no_water_conc_flag = -9999          # Flag for cells without water, so no concentration values
export_sediment = True             # Set True to export sediments results
timestep = [1, "hour"]              # Printing timestep
compartments_and_cells = {          # Compartments and cells to export
    "RIVER_NETWORK_REACHES": {
        "1": ["all", "all", "all"]
    }
}

########################################################
# DON'T CHANGE BELOW THIS POINT-------
########################################################

# Call the function with all individual arguments
############################

gJSON_lib.Gen_Input_Driver(
    dir2save_input_files=dir2save_input_files,
    running_on_docker=running_on_docker,
    project_name=project_name,
    geographical_location=geographical_location,
    authors=authors,
    date=date,
    comment=comment,
    hostmodel=hostmodel,
    solver=solver,
    run_mode_debug=run_mode_debug,
    use_num_threads=use_num_threads,
    ic_all_value=ic_all_value,
    ic_all_units=ic_all_units,
    bgc_module_name=bgc_module_name,
    path2selected_NATIVE_BGC_FLEX_framework=path2selected_NATIVE_BGC_FLEX_framework,
    # PHREEQC-specific settings (ignored if bgc_module_name != "PHREEQC")
    phreeqc_input_filepath=phreeqc_input_filepath,
    phreeqc_database_filepath=phreeqc_database_filepath,
    phreeqc_mobile_species=phreeqc_mobile_species,
    phreeqc_component_h2o=phreeqc_component_h2o,
    phreeqc_temperature_mapping=phreeqc_temperature_mapping,
    phreeqc_pressure_mapping=phreeqc_pressure_mapping,
    phreeqc_chemical_species_names=phreeqc_chemical_species_names,
    # Transport and other modules
    td_module_name=td_module_name,
    td_module_dispersion_xyz=td_module_dispersion_xyz,
    le_module_name=le_module_name,
    le_module_config=le_module_config,
    ts_module_name=ts_module_name,
    ts_sediment_compartment=ts_sediment_compartment,
    ts_transport_compartment=ts_transport_compartment,
    ts_direction=ts_direction,
    ts_erosion_inhibit_compartment=ts_erosion_inhibit_compartment,
    ts_data_format=ts_data_format,
    ts_mmf_defaults=ts_mmf_defaults,
    ts_mmf_parameters=ts_mmf_parameters,
    ts_hbvsed_defaults=ts_hbvsed_defaults,
    ts_hbvsed_parameters=ts_hbvsed_parameters,
    ts_hbvsed_monthly_erosion_factor=ts_hbvsed_monthly_erosion_factor,
    si_module_name=si_module_name,
    si_sediment_compartment=si_sediment_compartment,
    ss_method=ss_method,
    ss_metadata_source=ss_metadata_source,
    ss_metadata_comment=ss_metadata_comment,
    ss_method_csv_config=ss_method_csv_config,
    ss_method_copernicus_basin_info=ss_method_copernicus_basin_info,
    ss_method_copernicus_openwq_h5_results_file_example_for_mapping_key=ss_method_copernicus_openwq_h5_results_file_example_for_mapping_key,
    ss_method_copernicus_compartment_name_for_load=ss_method_copernicus_compartment_name_for_load,
    ss_method_copernicus_nc_lc_dir=ss_method_copernicus_nc_lc_dir,
    ss_method_copernicus_period=ss_method_copernicus_period,
    ss_method_copernicus_default_loads_bool=ss_method_copernicus_default_loads_bool,
    ss_method_copernicus_optional_custom_annual_load_coeffs_per_lulc_class=ss_method_copernicus_optional_custom_annual_load_coeffs_per_lulc_class,
    ss_method_copernicus_annual_to_seasonal_loads_method=ss_method_copernicus_annual_to_seasonal_loads_method,
    ewf_method=ewf_method,
    ewf_method_fixedval_comment=ewf_method_fixedval_comment,
    ewf_method_fixedval_source=ewf_method_fixedval_source,
    ewf_method_fixedval_chem_name=ewf_method_fixedval_chem_name,
    ewf_method_fixedval_value=ewf_method_fixedval_value,
    ewf_method_fixedval_units=ewf_method_fixedval_units,
    ewf_method_fixedval_external_inputflux_name=ewf_method_fixedval_external_inputflux_name,
    output_format=output_format,
    chemical_species=chemical_species,
    units=units,
    no_water_conc_flag=no_water_conc_flag,
    export_sediment=export_sediment,
    compartments_and_cells=compartments_and_cells,
    timestep=timestep
)
