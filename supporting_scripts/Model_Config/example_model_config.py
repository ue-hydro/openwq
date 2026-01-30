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
dir2save_input_files = "/Users/diogocosta/Documents/openwq_code/test_input_file_generator/"

################
# Computational settings
################
solver = "BE"           # "BE" or "SUNDIALS"
run_mode_debug = True   # "True" will print the solver derivatives that lead to the concentrations
use_num_threads = 4     # number of threads, accepts "all"

################
### General configuration
################
ic_all_value=2          # initial conditions value, concentration or mass
ic_all_units="mg/l"     # set the units of ic_all_value

################
### MODULES
################

"""
--------------------------------------
Biogeochemistry module
Options: 
    NATIVE_BGC_FLEX
    PHREEQC"
# in this case, "NONE" is not an option
--------------------------------------
"""
bgc_module_name = "NATIVE_BGC_FLEX"
path2selected_NATIVE_BGC_FLEX_framework = "config_support_lib/BGC_framework_examples/openWQ_Ncycling_example.json"
# TODO - set up for PHREQC

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
td_module_dispersion_xyz = [0.3,0.3,0.3]


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
ts_module_name = "NONE"
ts_sediment_compartment = "RIVER_NETWORK_REACHES"
# TODO - SETUP FOR HYPE_MMF AND HYPE_HBVSED OPTIONS

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
### Output settings
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
# External water fluxes
external_water_fluxes = {
    "1": {
        "LABEL": "SUMMA",
        "FILEPATH": "openwq_in/openwq_EWF_mizu.json"
    }
}
################
# Sink sources
sink_sources = {
    "1": {
        "LABEL": "fertilizer_N_test",
        "FILEPATH": "openwq_in/openWQ_fertilizer_N.json"
    }
}


############################
# ⚠️⚠️⚠️
# DON'T CHANGE BELOW THIS POINT-------
# ⚠️⚠️⚠️
############################
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
    external_water_fluxes=external_water_fluxes,
    sink_sources=sink_sources,
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
    output_format=output_format,
    chemical_species=chemical_species,
    units=units,
    no_water_conc_flag=no_water_conc_flag,
    export_sediment=export_sediment,
    compartments_and_cells=compartments_and_cells,
    timestep=timestep
)
