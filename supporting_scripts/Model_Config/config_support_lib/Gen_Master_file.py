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
"""
OpenWQ Configuration Generator
Simple function to create OpenWQ JSON configuration from individual arguments
"""

import json
import re
from typing import Dict, List, Union, Any
import os

def create_master_json(

        json_header_comment: [],

        # paths
        master_file_fullpath: str,
        config_file_fullpath: str,
        bgc_config_filepath: str,
        td_config_filepath: str,
        le_config_filepath: str,
        ts_config_filepath: str,
        si_config_filepath: str,
        ss_config_filepath: str,
        ewf_config_filepath: str,
        output_file_fullpath: str,

        # Top-level settings
        project_name: str,
        geographical_location: str,
        authors: str,
        date: str,
        comment: str,
        solver: str,

        # Computational settings
        run_mode_debug: bool,
        use_num_threads: Union[int, str],

        # Biogeochemistry module
        bgc_module_name: str,

        # Transport Dissolved module
        td_module_name: str,

        # Lateral Exchange module
        le_module_name: str,

        # Transport Sediments module
        ts_module_name: str,
        ts_sediment_compartment: str,

        # Sorption Isotherm module
        si_module_name: str,
        si_sediment_compartment: str,

        # Sink and Source
        ss_method_csv_metadata_source: str,

        # External water fluxes
        ewf_method_fixedval_source: str,

        # Output settings
        output_format: str,
        chemical_species: List[int],
        units: str,
        no_water_conc_flag: int,
        export_sediment: bool,
        compartments_and_cells: Dict[str, Dict[str, List]],
        timestep: List[Union[int, str]]

) -> None:
    """
    Create OpenWQ configuration JSON file from individual arguments.

    All parameters are required - no defaults are provided.
    """

    # Build the JSON structure
    config = {
        "PROJECT_NAME": project_name,
        "GEOGRAPHICAL_LOCATION": geographical_location,
        "AUTHORS": authors,
        "DATE": date,
        "COMMENT": comment,
        "COMPUTATIONAL_SETTINGS": {
            "RUN_MODE_DEBUG": run_mode_debug,
            "USE_NUM_THREADS": use_num_threads
        },
        "SOLVER": solver,
        "OPENWQ_INPUT": {
            "CONFIG_FILEPATH": config_file_fullpath,
            "EXTERNAL_WATER_FLUXES": {
                "1": {
                    "LABEL": ewf_method_fixedval_source,
                    "FILEPATH": ewf_config_filepath,
                }
            },
            "SINK_SOURCE": {
                "1": {
                    "LABEL": ss_method_csv_metadata_source,
                    "FILEPATH": ss_config_filepath,
                }
            }
        },
        "MODULES": {
            "BIOGEOCHEMISTRY": {
                "MODULE_NAME": bgc_module_name,
                "MODULE_CONFIG_FILEPATH": bgc_config_filepath
            },
            "TRANSPORT_DISSOLVED": {
                "MODULE_NAME": td_module_name,
                "MODULE_CONFIG_FILEPATH": td_config_filepath
            },
            "LATERAL_EXCHANGE": {
                "MODULE_NAME": le_module_name,
                "MODULE_CONFIG_FILEPATH": le_config_filepath
            },
            "TRANSPORT_SEDIMENTS": {
                "MODULE_NAME": ts_module_name,
                "SEDIMENT_COMPARTMENT": ts_sediment_compartment,
                "MODULE_CONFIG_FILEPATH": ts_config_filepath
            },
            "SORPTION_ISOTHERM": {
                "MODULE_NAME": si_module_name,
                "SEDIMENT_COMPARTMENT": si_sediment_compartment,
                "MODULE_CONFIG_FILEPATH": si_config_filepath
            }
        },
        "OPENWQ_OUTPUT": {
            "RESULTS_FOLDERPATH": output_file_fullpath,
            "FORMAT": output_format,
            "CHEMICAL_SPECIES": chemical_species,
            "UNITS": units,
            "NO_WATER_CONC_FLAG": no_water_conc_flag,
            "EXPORT_SEDIMENT": export_sediment,
            "COMPARTMENTS_AND_CELLS": compartments_and_cells,
            "TIMESTEP": timestep
        }
    }

    # Save to file with lists formatted horizontally

    # First convert to JSON string with standard formatting
    json_string = json.dumps(config, indent=4)

    # Regex pattern to find arrays/lists and format them on single lines
    # This finds patterns like: [\n        item,\n        item\n    ]
    def compact_list(match):
        # Extract the list content and format it on one line
        list_content = match.group(0)
        # Remove newlines and extra spaces within the list
        compact = re.sub(r'\[\s+', '[', list_content)
        compact = re.sub(r'\s+\]', ']', compact)
        compact = re.sub(r',\s+', ', ', compact)
        return compact

    # Pattern matches lists (arrays) with their contents
    pattern = r'\[[\s\n]*(?:[^\[\]]*?)[\s\n]*\]'

    # Keep applying the pattern until no more changes (for nested structures)
    prev_string = ""
    while prev_string != json_string:
        prev_string = json_string
        json_string = re.sub(pattern, compact_list, json_string)

    # Write to file with comment lines at the top
    with open(master_file_fullpath, 'w') as f:
        # Write comment lines first
        for comment in json_header_comment:
            f.write(comment + "\n")
        # Write the JSON content
        f.write(json_string)

    print(f"✓ Master file saved to: {master_file_fullpath}")