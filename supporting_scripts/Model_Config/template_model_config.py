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
# 🐷 General information
################
project_name = "Demonstration"
geographical_location = "NA"
authors = "Diogo Costa"
date = "January, 2026"
comment = "This a demonstration input file setup"

hostmodel = "mizuroute" # mizuroute or summa

# Directory to save the master file and the folder with the other input files,
dir2save_input_files = "/Users/diogocosta/Documents/openwq_code/test_input_file_generator/"

################
# 🐨 Computational settings
################
solver = "BE"           # "BE" or "SUNDIALS"
run_mode_debug = True   # "True" will print the solver derivatives that lead to the concentrations
use_num_threads = 4     # number of threads, accepts "all"

################
# 🐹 General configuration
################
ic_all_value=2          # initial conditions value, concentration or mass
ic_all_units="mg/l"     # set the units of ic_all_value

################
### MODULES
################

"""
--------------------------------------
🐮 Biogeochemistry module
Options: 
    NATIVE_BGC_FLEX
    PHREEQC"
# in this case, "NONE" is not an option
--------------------------------------
"""
bgc_module_name = "NATIVE_BGC_FLEX"
path2selected_NATIVE_BGC_FLEX_framework = "config_support_lib/examples_BGC_frameworks/openWQ_Ncycling_example.json"
# TODO - set up for PHREQC

"""
--------------------------------------
🐸 Transport Dissolved module
Options: 
    NATIVE_TD_ADV, 
    OPENWQ_NATIVE_TD_ADVDISP, 
    "NONE"
--------------------------------------
"""
td_module_name = "NATIVE_TD_ADV"
# if td_module_name=OPENWQ_NATIVE_TD_ADVDISP
td_module_dispersion_xyz = [0.3,0.3,0.3]


"""
--------------------------------------
🦁 Lateral Exchange module
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
🐥 Transport Sediments module
Options:
    HYPE_MMF
    HYPE_HBVSED
    NONE
--------------------------------------
"""
ts_module_name = "NONE"
ts_sediment_compartment = "RIVER_NETWORK_REACHES"
# TODO - SETUP FOR HYPE_MMF AND HYPE_HBVSED OPTIONS

"""
🐞 Sorption Isotherm module
Options:
    FREUNDLICH
    LANGMUIR
    NONE
"""
si_module_name = "NONE"
si_sediment_compartment = "SUMMA_RUNOFF"
# TODO - SETUP FOR FREUNDLICH AND LANGMUIR OPTIONS

################
# 🦞 Output settings
################
output_format = "HDF5"              # HDF5 or CSV (but HDF5 is recommended as can be used with other scripts
chemical_species = [1, 2, 3, 4, 5]  # Chemical
units = "MG/L"                      # Units (concentrations or mass)
no_water_conc_flag = -9999          # Flag for cells without water, so no concentration values
export_sediment = False             # Set True to export sediments results
timestep = [1, "hour"]              # Printing timestep
compartments_and_cells = {          # Compartments and cells to export
    "RIVER_NETWORK_REACHES": {
        "1": ["all", "all", "all"]
    }
}

################
# 🦎 Sink sources
################
ss_method = "using_copernicus_lulc" # "load_from_csv" or "using_copernicus_lulc"
################
# 👉🏼 if ss_method = "load_from_csv"
ss_metadata_source = "Just for demonstration"
ss_metadata_comment = "Leave any comments needed for future reference"
ss_method_csv_config = source_sink_configs=[
        {
            "Chemical_name": "NO3-N",
            "Compartment_name": "RIVER_NETWORK_REACHES",
            "Type": "source",
            "Units": "kg",
            "Filepath": "/Users/diogocosta/Documents/openwq_code/6_mizuroute_cslm_openwq/route/build/openwq/openwq/supporting_scripts/Model_Config/config_support_lib/examples_SS_in_cvs/NO3_fertilizer.csv",
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
# 👉🏼 if ss_method = load_from_copernicus
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
ss_method_copernicus_nc_lc_dir ='/Users/diogocosta/Documents/ESACCI-LC/'
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
################
# 🐠 External water fluxes
################
ewf_method = "fixed_value" # for now only "fixed value" option is available
ewf_method_fixedval_comment = "External model: summa"
ewf_method_fixedval_source = "Just for demonstration purposes"
ewf_method_fixedval_chem_name = "NO3-N"
ewf_method_fixedval_value = 2.5
ewf_method_fixedval_units = "mg/l"
ewf_method_fixedval_external_inputflux_name = "SUMMA_RUNOFF"


# ⚠️⚠️⚠️ ⚠️⚠️⚠️ ⚠️⚠️⚠️ ⚠️⚠️⚠️ ⚠️⚠️⚠️ ⚠️⚠️⚠️
########################################################
# DON'T CHANGE BELOW THIS POINT-------
########################################################
# ⚠️⚠️⚠️ ⚠️⚠️⚠️ ⚠️⚠️⚠️ ⚠️⚠️⚠️ ⚠️⚠️⚠️ ⚠️⚠️⚠️

# Call the function with all individual arguments
############################

gJSON_lib.Gen_Input_Driver(
    dir2save_input_files=dir2save_input_files,
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
    td_module_name=td_module_name,
    td_module_dispersion_xyz=td_module_dispersion_xyz,
    le_module_name=le_module_name,
    le_module_config=le_module_config,
    ts_module_name=ts_module_name,
    ts_sediment_compartment=ts_sediment_compartment,
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
