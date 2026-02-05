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

from typing import Dict, List, Union, Any, Optional
import datetime
import os
import re

import Gen_Master_file as mJSON_lib
import Gen_Config_file as cJSON_lib
import Gen_TDmodule_file as tdmJSON_lib
import Gen_LEmodule_file as lemJSON_lib
import Load_BGQmodule_file as bgqmJSON_lib
import Gen_PHREEQCmodule_file as phreeqcJSON_lib
import Gen_TSmodule_file as tsmJSON_lib
import Gen_SS_Driver as ssJSON_lib
import Gen_EWF_Driver as ewfJSON_lib

def _parse_docker_volume_mount(docker_compose_path: str) -> tuple:
    """
    Parse docker-compose.yml to extract the volume mount mapping.
    Returns (host_root, container_root) or (None, None) if parsing fails.

    The docker-compose.yml volume format is: 'host_path:container_path[:options]'
    where host_path can be relative (resolved from docker-compose.yml location).
    """
    try:
        with open(docker_compose_path, 'r') as f:
            content = f.read()

        # Find volume mount lines matching pattern like '- host_path:container_path:options'
        volume_pattern = re.compile(r'^\s*-\s+([^:]+):([^:]+)(?::.*)?$', re.MULTILINE)
        matches = volume_pattern.findall(content)

        if not matches:
            print("WARNING: No volume mounts found in docker-compose.yml")
            return None, None

        # Use the first volume mount
        host_path_raw, container_path = matches[0]
        host_path_raw = host_path_raw.strip()
        container_path = container_path.strip()

        # Resolve relative host path from docker-compose.yml directory
        compose_dir = os.path.dirname(os.path.abspath(docker_compose_path))
        host_root = os.path.normpath(os.path.join(compose_dir, host_path_raw))

        # Ensure paths end with '/'
        if not host_root.endswith('/'):
            host_root += '/'
        if not container_path.endswith('/'):
            container_path += '/'

        return host_root, container_path

    except FileNotFoundError:
        print(f"WARNING: docker-compose.yml not found at: {docker_compose_path}")
        return None, None
    except Exception as e:
        print(f"WARNING: Failed to parse docker-compose.yml: {e}")
        return None, None


def _correct_path_for_docker(path: str, host_root: str, container_root: str) -> str:
    """
    Replace host root prefix with container root in a file path.
    If the path doesn't start with host_root, return it unchanged.
    """
    if path and path.startswith(host_root):
        return container_root + path[len(host_root):]
    return path


def _correct_dict_paths_for_docker(d: dict, host_root: str, container_root: str) -> dict:
    """
    Recursively correct file paths in a dictionary.
    Looks for string values that contain the host root path and replaces them.
    """
    if d is None:
        return d
    corrected = {}
    for key, value in d.items():
        if isinstance(value, str) and host_root in value:
            corrected[key] = _correct_path_for_docker(value, host_root, container_root)
        elif isinstance(value, dict):
            corrected[key] = _correct_dict_paths_for_docker(value, host_root, container_root)
        elif isinstance(value, list):
            corrected[key] = _correct_list_paths_for_docker(value, host_root, container_root)
        else:
            corrected[key] = value
    return corrected


def _correct_list_paths_for_docker(lst: list, host_root: str, container_root: str) -> list:
    """
    Recursively correct file paths in a list.
    """
    if lst is None:
        return lst
    corrected = []
    for item in lst:
        if isinstance(item, str) and host_root in item:
            corrected.append(_correct_path_for_docker(item, host_root, container_root))
        elif isinstance(item, dict):
            corrected.append(_correct_dict_paths_for_docker(item, host_root, container_root))
        elif isinstance(item, list):
            corrected.append(_correct_list_paths_for_docker(item, host_root, container_root))
        else:
            corrected.append(item)
    return corrected


def Gen_Input_Driver(
        # where to save input files
        dir2save_input_files: str,

        # Docker support: correct paths from host to container
        running_on_docker: bool,

        # Top-level settings
        project_name: str,
        geographical_location: str,
        authors: str,
        date: str,
        comment: str,
        hostmodel: str,

        # Computational settings
        solver: str,
        run_mode_debug: bool,
        use_num_threads: Union[int, str],

        # Initial conditions: set equal to all
        ic_all_value: float,
        ic_all_units: str,

        # Biogeochemistry module
        bgc_module_name: str,
        path2selected_NATIVE_BGC_FLEX_framework: str = "",

        # PHREEQC-specific settings (used when bgc_module_name = "PHREEQC")
        phreeqc_input_filepath: str = "",
        phreeqc_database_filepath: str = "",
        phreeqc_mobile_species: Optional[List[str]] = None,
        phreeqc_component_h2o: bool = True,
        phreeqc_temperature_mapping: Optional[Dict[str, str]] = None,
        phreeqc_pressure_mapping: Optional[Dict[str, str]] = None,
        phreeqc_chemical_species_names: Optional[List[str]] = None,

        # Transport Dissolved module
        td_module_name: str = "NATIVE_TD_ADV",
        td_module_dispersion_xyz: List[Union[int, float]] = None,
        td_module_characteristic_length_m: Union[int, float] = 1.0,

        # Lateral Exchange module
        le_module_name: str = "NONE",
        le_module_config: List[Dict[str, Union[str, int, float]]] = None,

        # Transport Sediments module
        ts_module_name: str = "NONE",
        ts_sediment_compartment: str = "RIVER_NETWORK_REACHES",
        ts_transport_compartment: str = "RUNOFF",
        ts_direction: str = "z",
        ts_erosion_inhibit_compartment: str = "ILAYERVOLFRACWAT_SNOW",
        ts_data_format: str = "JSON",
        ts_mmf_defaults: Optional[Dict[str, float]] = None,
        ts_mmf_parameters: Optional[Dict[str, Dict[str, list]]] = None,
        ts_hbvsed_defaults: Optional[Dict[str, float]] = None,
        ts_hbvsed_parameters: Optional[Dict[str, Dict[str, list]]] = None,
        ts_hbvsed_monthly_erosion_factor: Optional[List[float]] = None,

        # Sorption Isotherm module
        si_module_name: str = "NONE",
        si_sediment_compartment: str = "RIVER_NETWORK_REACHES",

        # Sink and Source settings
        ss_method: str = "load_from_csv",
        ss_metadata_source: str = "",
        ss_metadata_comment: str = "",
        ss_method_csv_config: List[Dict[str, Union[str, int]]] = None,
        ss_method_copernicus_basin_info: Dict[str, str] = None,
        ss_method_copernicus_openwq_h5_results_file_example_for_mapping_key: Dict[str, str] = None,
        ss_method_copernicus_compartment_name_for_load: str = "",
        ss_method_copernicus_nc_lc_dir: str = "",
        ss_method_copernicus_period: List[Union[int, float]] = None,
        ss_method_copernicus_default_loads_bool: bool = True,
        ss_method_copernicus_annual_to_seasonal_loads_method: str = "uniform",

        # External water fluxes
        ewf_method: str = "fixed_value",
        ewf_method_fixedval_comment: str = "",
        ewf_method_fixedval_source: str = "",
        ewf_method_fixedval_chem_name: str = "",
        ewf_method_fixedval_value: str = "",
        ewf_method_fixedval_units: str = "",
        ewf_method_fixedval_external_inputflux_name: str = "",

        # Output settings
        output_format: str = "HDF5",
        chemical_species: List[str] = None,
        units: str = "MG/L",
        no_water_conc_flag: int = -9999,
        export_sediment: bool = False,
        compartments_and_cells: Dict[str, Dict[str, List]] = None,
        timestep: List[Union[int, str]] = None,

        # Optional parameters (MUST be at the end)
        ss_method_copernicus_optional_custom_annual_load_coeffs_per_lulc_class: Optional[Dict[int, Dict[str, float]]] = None

) -> None:

    print(f"Project: {project_name}")
    print(f"Authors: {authors}")

    # Docker path correction setup
    # When running_on_docker=True, the script runs on the HOST but generates JSON files
    # with paths that will be used INSIDE the Docker container. So we need two sets of paths:
    #   - Host paths: for writing files to disk (file I/O)
    #   - Docker paths: for embedding inside JSON content (read by OpenWQ at runtime in container)
    _docker_host_root = None
    _docker_container_root = None

    if running_on_docker:
        # Locate docker-compose.yml relative to this script
        # This script is in: .../supporting_scripts/Model_Config/config_support_lib/
        # docker-compose.yml is in: .../  (openwq root)
        script_dir = os.path.dirname(os.path.abspath(__file__))
        docker_compose_path = os.path.normpath(
            os.path.join(script_dir, '..', '..', '..', 'docker-compose.yml'))

        _docker_host_root, _docker_container_root = _parse_docker_volume_mount(docker_compose_path)

        if _docker_host_root and _docker_container_root:
            print(f"  Docker mode: mapping host '{_docker_host_root}' -> container '{_docker_container_root}'")

            # Correct paths that are ONLY embedded in JSON content (NOT used for file I/O)
            # These are paths that OpenWQ reads at runtime inside the container.
            #
            # NOTE: The following paths are NOT corrected here because they are used
            # for file I/O during config generation (the script reads/copies these files):
            #   - path2selected_NATIVE_BGC_FLEX_framework (read by Load_BGQmodule_file.py)
            #   - phreeqc_input_filepath / phreeqc_database_filepath (copied by Gen_PHREEQCmodule_file.py)
            #   - ss_method_copernicus_basin_info (shapefile read by Gen_SS_Driver.py)
            #   - ss_method_copernicus_nc_lc_dir (NetCDF files read by Gen_SS_Driver.py)
            #   - ss_method_copernicus_openwq_h5_results_file_example_for_mapping_key (HDF5 read)
            #
            # Only CSV source/sink Filepath entries are purely embedded in JSON (the CSV files
            # themselves are read by OpenWQ at runtime, not by the config generator).

            # Correct CSV source/sink file paths (embedded in SS JSON, read at runtime by OpenWQ)
            if ss_method_csv_config is not None:
                ss_method_csv_config = _correct_list_paths_for_docker(
                    ss_method_csv_config, _docker_host_root, _docker_container_root)
        else:
            print("WARNING: running_on_docker=True but could not determine Docker path mapping. "
                  "Paths will NOT be corrected.")

    # Set defaults for mutable args
    if td_module_dispersion_xyz is None:
        td_module_dispersion_xyz = [0.0, 0.0, 0.0]
    if le_module_config is None:
        le_module_config = []
    if ss_method_csv_config is None:
        ss_method_csv_config = []
    if chemical_species is None:
        chemical_species = ["NO3-N"]
    if timestep is None:
        timestep = [1, "hour"]

    # Add comment lines at the top
    json_header_comment = [
        "",
        "// #####################",
        "// JSON file automatically generated using OpenWQ's supporting scripts",
        f"// {datetime.datetime.now().replace(microsecond=0)}\n"
        "// #####################",
        ""
    ]

    # Generate paths for file I/O (host filesystem — where files are actually written)
    master_file_fullpath = os.path.join(f'{dir2save_input_files}', "openWQ_master.json")
    general_json_input_dir = os.path.join(f'{dir2save_input_files}', "openwq_in")
    config_file_fullpath = os.path.join(f'{general_json_input_dir}', "openWQ_config.json")
    bgc_config_filepath = os.path.join(f'{general_json_input_dir}', f"openWQ_MODULE_{bgc_module_name}.json")
    td_config_filepath = os.path.join(f'{general_json_input_dir}', f"openWQ_MODULE_{td_module_name}.json")
    le_config_filepath = os.path.join(f'{general_json_input_dir}', f"openWQ_MODULE_{le_module_name}.json")
    ts_config_filepath = os.path.join(f'{general_json_input_dir}', f"openWQ_MODULE_{ts_module_name}.json")
    si_config_filepath = os.path.join(f'{general_json_input_dir}', f"openWQ_MODULE_{si_module_name}.json")
    ss_config_filepath = os.path.join(f'{general_json_input_dir}', f"openWQ_SS_{ss_method}.json")
    ewf_config_filepath = os.path.join(f'{general_json_input_dir}', f"openWQ_EWF_fixed_value.json")
    output_file_fullpath = os.path.join(f'{dir2save_input_files}', "openwq_out/")

    # When running_on_docker, create Docker-corrected paths for JSON content
    # These paths are what OpenWQ sees at runtime inside the container
    if running_on_docker and _docker_host_root and _docker_container_root:
        docker_config_file_fullpath = _correct_path_for_docker(config_file_fullpath, _docker_host_root, _docker_container_root)
        docker_bgc_config_filepath = _correct_path_for_docker(bgc_config_filepath, _docker_host_root, _docker_container_root)
        docker_td_config_filepath = _correct_path_for_docker(td_config_filepath, _docker_host_root, _docker_container_root)
        docker_le_config_filepath = _correct_path_for_docker(le_config_filepath, _docker_host_root, _docker_container_root)
        docker_ts_config_filepath = _correct_path_for_docker(ts_config_filepath, _docker_host_root, _docker_container_root)
        docker_si_config_filepath = _correct_path_for_docker(si_config_filepath, _docker_host_root, _docker_container_root)
        docker_ss_config_filepath = _correct_path_for_docker(ss_config_filepath, _docker_host_root, _docker_container_root)
        docker_ewf_config_filepath = _correct_path_for_docker(ewf_config_filepath, _docker_host_root, _docker_container_root)
        docker_output_file_fullpath = _correct_path_for_docker(output_file_fullpath, _docker_host_root, _docker_container_root)
    else:
        # No Docker correction — JSON content paths are the same as host paths
        docker_config_file_fullpath = config_file_fullpath
        docker_bgc_config_filepath = bgc_config_filepath
        docker_td_config_filepath = td_config_filepath
        docker_le_config_filepath = le_config_filepath
        docker_ts_config_filepath = ts_config_filepath
        docker_si_config_filepath = si_config_filepath
        docker_ss_config_filepath = ss_config_filepath
        docker_ewf_config_filepath = ewf_config_filepath
        docker_output_file_fullpath = output_file_fullpath

    ###############
    # Call create_master_json
    ###############

    mJSON_lib.create_master_json(
        json_header_comment=json_header_comment,
        master_file_fullpath=master_file_fullpath,           # Host path: where to write the file
        config_file_fullpath=docker_config_file_fullpath,    # Docker path: embedded in JSON content
        bgc_config_filepath=docker_bgc_config_filepath,
        td_config_filepath=docker_td_config_filepath,
        le_config_filepath=docker_le_config_filepath,
        ts_config_filepath=docker_ts_config_filepath,
        si_config_filepath=docker_si_config_filepath,
        ss_config_filepath=docker_ss_config_filepath,
        ewf_config_filepath=docker_ewf_config_filepath,
        output_file_fullpath=docker_output_file_fullpath,
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
        ts_transport_compartment=ts_transport_compartment,
        si_module_name=si_module_name,
        si_sediment_compartment=si_sediment_compartment,
        ss_method=ss_method,
        ss_metadata_source=ss_metadata_source,
        ewf_method=ewf_method,
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
    # Determine compartment names from host model
    ###############

    if (hostmodel == "mizuroute"):
        compartment_names = ["RIVER_NETWORK_REACHES"]
    elif (hostmodel == "summa"):
        compartment_names = ["SCALARCANOPYWAT",
                             "ILAYERVOLFRACWAT_SNOW",
                             "RUNOFF",
                             "ILAYERVOLFRACWAT_SOIL",
                             "SCALARAQUIFER"]

    ###############
    # Call create_config_json
    # For PHREEQC: species come from the user-provided list (auto-discovered at runtime);
    #              no CYCLING_FRAMEWORK is needed.
    # For NATIVE_BGC_FLEX: species and cycling frameworks are hardcoded from the framework file.
    ###############

    if bgc_module_name == "PHREEQC":

        # PHREEQC species are discovered at runtime from the database/input file.
        # The user provides expected species names for initial conditions.
        if phreeqc_chemical_species_names is None or len(phreeqc_chemical_species_names) == 0:
            print("WARNING: phreeqc_chemical_species_names is empty. "
                  "Config file will have no initial conditions. "
                  "PHREEQC will use default initial conditions from the .pqi file.")
            chem_names = []
        else:
            chem_names = phreeqc_chemical_species_names

        cJSON_lib.create_config_json(
            config_file_fullpath=config_file_fullpath,
            json_header_comment=json_header_comment,
            compartment_names=compartment_names,
            cycling_framework=None,  # PHREEQC does not use cycling frameworks
            chemical_species_names=chem_names,
            ic_all_value=ic_all_value,
            ic_all_units=ic_all_units,
        )
    else:
        # NATIVE_BGC_FLEX path (original)
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
    # Call create_td_module_json
    ###############

    tdmJSON_lib.create_td_module_json(
        td_config_filepath=td_config_filepath,
        json_header_comment=json_header_comment,
        td_module_name=td_module_name,
        td_module_dispersion_xyz=td_module_dispersion_xyz,
        td_module_characteristic_length_m=td_module_characteristic_length_m
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
    # Call create_ts_module_json
    ###############

    tsmJSON_lib.create_ts_module_json(
        ts_config_filepath=ts_config_filepath,
        json_header_comment=json_header_comment,
        ts_module_name=ts_module_name,
        ts_direction=ts_direction,
        ts_erosion_inhibit_compartment=ts_erosion_inhibit_compartment,
        ts_data_format=ts_data_format,
        ts_mmf_defaults=ts_mmf_defaults,
        ts_mmf_parameters=ts_mmf_parameters,
        ts_hbvsed_defaults=ts_hbvsed_defaults,
        ts_hbvsed_parameters=ts_hbvsed_parameters,
        ts_hbvsed_monthly_erosion_factor=ts_hbvsed_monthly_erosion_factor
    )

    ###############
    # Call BGC module setup (NATIVE_BGC_FLEX or PHREEQC)
    ###############

    if bgc_module_name == "PHREEQC":
        # Generate PHREEQC module JSON + copy database and input files
        phreeqcJSON_lib.create_phreeqc_module_json(
            bgc_config_filepath=bgc_config_filepath,
            json_header_comment=json_header_comment,
            phreeqc_input_filepath=phreeqc_input_filepath,
            phreeqc_database_filepath=phreeqc_database_filepath,
            phreeqc_mobile_species=phreeqc_mobile_species if phreeqc_mobile_species else ["H"],
            phreeqc_component_h2o=phreeqc_component_h2o,
            phreeqc_temperature_mapping=phreeqc_temperature_mapping,
            phreeqc_pressure_mapping=phreeqc_pressure_mapping,
        )
    elif bgc_module_name == "NATIVE_BGC_FLEX":
        # Copy the selected BGC framework file
        bgqmJSON_lib.load_bgq_module_json(
            json_header_comment=json_header_comment,
            path2selected_NATIVE_BGC_FLEX_framework=path2selected_NATIVE_BGC_FLEX_framework,
            bgc_config_filepath=bgc_config_filepath
        )
    else:
        print(f"WARNING: Unknown bgc_module_name '{bgc_module_name}'. "
              f"Valid options are: NATIVE_BGC_FLEX, PHREEQC")

    ###############
    # Call gen_ss_driver
    ###############
    if (ss_method == "load_from_csv"):

        ssJSON_lib.set_ss_from_csv(
            ss_config_filepath=ss_config_filepath,
            json_header_comment=json_header_comment,
            ss_metadata_source=ss_metadata_source,
            ss_metadata_comment=ss_metadata_comment,
            ss_method_csv_config=ss_method_csv_config,
        )
    elif (ss_method == "using_copernicus_lulc"):

        ssJSON_lib.set_ss_from_copernicus_lulc_with_loads(
            ss_config_filepath=ss_config_filepath,
            json_header_comment=json_header_comment,
            ss_metadata_source=ss_metadata_source,
            ss_metadata_comment=ss_metadata_comment,
            ss_method_copernicus_basin_info=ss_method_copernicus_basin_info,
            ss_method_copernicus_nc_lc_dir=ss_method_copernicus_nc_lc_dir,
            ss_method_copernicus_openwq_h5_results_file_example_for_mapping_key=ss_method_copernicus_openwq_h5_results_file_example_for_mapping_key,
            ss_method_copernicus_compartment_name_for_load=ss_method_copernicus_compartment_name_for_load,
            ss_method_copernicus_period=ss_method_copernicus_period,
            ss_method_copernicus_default_loads_bool=ss_method_copernicus_default_loads_bool,
            optional_load_coefficients=ss_method_copernicus_optional_custom_annual_load_coeffs_per_lulc_class,
            ss_method_copernicus_annual_to_seasonal_loads_method=ss_method_copernicus_annual_to_seasonal_loads_method
        )
    elif (ss_method == "none"):
        print("  Skipping sink/source generation (ss_method='none')")

    else:
        print(
            f"WARNING: The SS method '{ss_method}' is unknown or not available for automatic generation. "
            f"Available methods: 'load_from_csv', 'using_copernicus_lulc', 'none'")

    ###############
    # Call gen_ewf_driver
    ###############
    if (ewf_method == "fixed_value"):
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
    elif (ewf_method == "none"):
        print("  Skipping external water flux generation (ewf_method='none')")
    else:
        print(
            f"WARNING: The EWF method '{ewf_method}' is unknown or not available for automatic generation. "
            f"Available methods: 'fixed_value', 'none'")
