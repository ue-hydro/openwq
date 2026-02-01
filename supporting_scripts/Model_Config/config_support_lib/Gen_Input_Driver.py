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
# !/usr/bin/env python3

from typing import Dict, List, Union, Any
import datetime
import os

import Gen_Master_file as mJSON_lib
import Gen_Config_file as cJSON_lib
import Gen_TDmodule_file as tdmJSON_lib
import Gen_LEmodule_file as lemJSON_lib
import Load_BGQmodule_file as bgqmJSON_lib
import Gen_SS_Driver as ssJSON_lib
import Gen_EWF_Driver as ewfJSON_lib

def Gen_Input_Driver(
        # where to save input files
        dir2save_input_files: str,

        # Top-level settings
        project_name: str,
        geographical_location: str,
        authors: str,
        date: str,
        comment: str,
        hostmodel:str,

        # Computational settings
        solver: str,
        run_mode_debug: bool,
        use_num_threads: Union[int, str],

        # Initial conditions: set equal to all
        ic_all_value: float,
        ic_all_units: str,

        # Biogeochemistry module
        bgc_module_name: str,
        path2selected_NATIVE_BGC_FLEX_framework: str,

        # Transport Dissolved module
        td_module_name: str,
        td_module_dispersion_xyz: List[Union[int, float]],

        # Lateral Exchange module
        le_module_name: str,
        le_module_config: List[Dict[str, Union[str, int, float]]],

        # Transport Sediments module
        ts_module_name: str,
        ts_sediment_compartment: str,

        # Sorption Isotherm module
        si_module_name: str,
        si_sediment_compartment: str,

        # Sink and Source settings
        ss_method: str,
        ss_method_csv_metadata_source: str,
        ss_method_csv_metadata_comment: str,
        ss_method_csv_config:List[Dict[str, Union[str, int]]],

        # External water fluxes
        ewf_method: str,
        ewf_method_fixedval_comment: str,
        ewf_method_fixedval_source: str,
        ewf_method_fixedval_chem_name: str,
        ewf_method_fixedval_value: str,
        ewf_method_fixedval_units: str,
        ewf_method_fixedval_external_inputflux_name: str,

        # Output settings
        output_format: str,
        chemical_species: List[int],
        units: str,
        no_water_conc_flag: int,
        export_sediment: bool,
        compartments_and_cells: Dict[str, Dict[str, List]],
        timestep: List[Union[int, str]]

) -> None:

    print(f"Project: {project_name}")
    print(f"Authors: {authors}")

    # Add comment lines at the top
    json_header_comment = [
        "",
        "// #####################",
        "// JSON file automatically generated using OpenWQ's supporting scripts",
        f"// {datetime.datetime.now().replace(microsecond=0)}\n"
        "// #####################",
        ""
    ]

    # Generate paths
    master_file_fullpath = os.path.join(f'{dir2save_input_files}', "openWQ_master.json")
    general_json_input_dir = os.path.join(f'{dir2save_input_files}', "openwq_in")
    config_file_fullpath = os.path.join(f'{general_json_input_dir}', "openWQ_config.json")
    bgc_config_filepath = os.path.join(f'{general_json_input_dir}', f"openWQ_MODULE_{bgc_module_name}.json")
    td_config_filepath = os.path.join(f'{general_json_input_dir}', f"openWQ_MODULE_{td_module_name}.json")
    le_config_filepath = os.path.join(f'{general_json_input_dir}', f"openWQ_MODULE_{le_module_name}.json")
    ts_config_filepath = os.path.join(f'{general_json_input_dir}', f"openWQ_MODULE_{ts_module_name}.json")
    si_config_filepath = os.path.join(f'{general_json_input_dir}', f"openWQ_MODULE_{si_module_name}.json")
    ss_config_filepath = os.path.join(f'{general_json_input_dir}', f"openWQ_SS_csv.json")
    ewf_config_filepath = os.path.join(f'{general_json_input_dir}', f"openWQ_EWF_fixed_value.json")
    output_file_fullpath = os.path.join(f'{dir2save_input_files}', "openwq_out/")

    ###############
    # Call create_master_json
    ###############

    mJSON_lib.create_master_json(
        json_header_comment=json_header_comment,
        master_file_fullpath=master_file_fullpath,
        config_file_fullpath=config_file_fullpath,
        bgc_config_filepath=bgc_config_filepath,
        td_config_filepath=td_config_filepath,
        le_config_filepath=le_config_filepath,
        ts_config_filepath=ts_config_filepath,
        si_config_filepath=si_config_filepath,
        ss_config_filepath=ss_config_filepath,
        ewf_config_filepath=ewf_config_filepath,
        output_file_fullpath=output_file_fullpath,
        project_name=project_name,
        geographical_location=geographical_location,
        authors=authors,
        date=date,
        comment=comment,
        solver=solver,
        run_mode_debug=run_mode_debug,
        use_num_threads=use_num_threads,
        bgc_module_name=bgc_module_name,
        td_module_name=td_module_name,
        le_module_name=le_module_name,
        ts_module_name=ts_module_name,
        ts_sediment_compartment=ts_sediment_compartment,
        si_module_name=si_module_name,
        si_sediment_compartment=si_sediment_compartment,
        ss_method_csv_metadata_source=ss_method_csv_metadata_source,
        ewf_method_fixedval_source=ewf_method_fixedval_source,
        output_format=output_format,
        chemical_species=chemical_species,
        units=units,
        no_water_conc_flag=no_water_conc_flag,
        export_sediment=export_sediment,
        compartments_and_cells=compartments_and_cells,
        timestep=timestep
    )

    ###############
    # Call create_config_json
    ###############

    if (hostmodel=="mizuroute"):
        compartment_names = ["RIVER_NETWORK_REACHES"]
    elif (hostmodel=="summa"):
        compartment_names = ["SCALARCANOPYWAT",
                             "ILAYERVOLFRACWAT_SNOW",
                             "RUNOFF",
                             "ILAYERVOLFRACWAT_SOIL",
                             "SCALARAQUIFER"]

    cJSON_lib.create_config_json(
        config_file_fullpath=config_file_fullpath,
        json_header_comment=json_header_comment,
        compartment_names=compartment_names,
        cycling_framework=["N_cycle"],
        chemical_species_names=["NO3-N", "NH4-N", "N_ORG_fresh", "N_ORG_stable", "N_ORG_active", "OUT"],
        ic_all_value=ic_all_value,
        ic_all_units=ic_all_units,
    )

    ###############
    # Call create_ts_module_json
    ###############

    tdmJSON_lib.create_td_module_json(
        td_config_filepath=td_config_filepath,
        json_header_comment =json_header_comment,
        td_module_name=td_module_name,
        td_module_dispersion_xyz=td_module_dispersion_xyz
    )

    ###############
    # Call create_le_module_json
    ###############

    lemJSON_lib.create_le_module_json(
        le_config_filepath=le_config_filepath,
        json_header_comment=json_header_comment,
        le_module_name=le_module_name,
        le_module_config=le_module_config
    )

    ###############
    # Call load_bgq_module_json
    ###############

    bgqmJSON_lib.load_bgq_module_json(
        json_header_comment=json_header_comment,
        path2selected_NATIVE_BGC_FLEX_framework=path2selected_NATIVE_BGC_FLEX_framework,
        bgc_config_filepath=bgc_config_filepath
    )


    ###############
    # Call gen_ss_driver
    ###############
    if (ss_method=="load_from_csv"):

        ssJSON_lib.set_ss_from_csv(
            ss_config_filepath=ss_config_filepath,
            json_header_comment=json_header_comment,
            ss_method_csv_metadata_source=ss_method_csv_metadata_source,
            ss_method_csv_metadata_comment=ss_method_csv_metadata_comment,
            ss_method_csv_config=ss_method_csv_config,
        )
    elif (ss_method=="using_copernicus_lulc"):

        # CONTINUE HERE
        results, summaries, rasters = ssJSON_lib.set_ss_from_copericus_lulc(
            ss_config_filepath=ss_config_filepath,
            basin_shapefile_path='/Users/diogocosta/Documents/openwq_code/6_mizuroute_cslm_openwq/test_case/mizuroute_in/shapefiles/finalcat_info_v1-0.shp',
            netcdf_copernicus_lc_dir='/Users/diogocosta/Documents/ESACCI-LC/',
            year_start=1993,
            year_end=1994
        )

        results = results

    else:
        print(f"WARNING: The SS method '{ss_method} is unkown or not available for automatic generation. Only method 'load_from_csv' is available")

    ###############
    # Call gen_ewf_driver
    ###############
    if (ewf_method=="fixed_value"):
        ewfJSON_lib.set_ewf_fixed_value(
            ewf_config_filepath=ewf_config_filepath,
            json_header_comment=json_header_comment,
            ewf_method_fixedval_comment=ewf_method_fixedval_comment,
            ewf_method_fixedval_source=ewf_method_fixedval_source,
            ewf_method_fixedval_chem_name=ewf_method_fixedval_chem_name,
            ewf_method_fixedval_value=ewf_method_fixedval_value,
            ewf_method_fixedval_units=ewf_method_fixedval_units,
            ewf_method_fixedval_external_inputflux_name=ewf_method_fixedval_external_inputflux_name
        )
    else:
        print(
            f"WARNING: The EWF method '{ewf_method} is unkown or not available for automatic generation. Only method 'fixed_value' is available")


