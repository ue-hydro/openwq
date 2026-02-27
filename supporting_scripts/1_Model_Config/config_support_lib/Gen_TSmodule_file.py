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
OpenWQ Transport Sediments Module Configuration Generator
Generates JSON configuration files for HYPE_MMF and HYPE_HBVSED erosion models.
"""

import json
import re
from pathlib import Path
from typing import List, Dict, Union, Optional


def create_ts_module_json(

        ts_config_filepath: str,

        json_header_comment: List[str],

        # Module name: "HYPE_MMF", "HYPE_HBVSED", or "NONE"
        ts_module_name: str,

        # Common configuration
        ts_direction: str = "z",
        ts_erosion_inhibit_compartment: str = "ILAYERVOLFRACWAT_SNOW",
        ts_data_format: str = "JSON",

        # HYPE_MMF parameters (defaults)
        ts_mmf_defaults: Optional[Dict[str, float]] = None,

        # HYPE_MMF spatially-varying parameters
        ts_mmf_parameters: Optional[Dict[str, Dict[str, list]]] = None,

        # HYPE_HBVSED parameters (defaults)
        ts_hbvsed_defaults: Optional[Dict[str, float]] = None,

        # HYPE_HBVSED spatially-varying parameters
        ts_hbvsed_parameters: Optional[Dict[str, Dict[str, list]]] = None,

        # HYPE_HBVSED monthly erosion factor (12 values, Jan-Dec)
        ts_hbvsed_monthly_erosion_factor: Optional[List[float]] = None

) -> None:
    """
    Create OpenWQ Transport Sediments module configuration JSON file.

    Parameters:
        ts_config_filepath: Full path where the JSON file will be saved
        json_header_comment: List of comment lines to add at the top of the file
        ts_module_name: Name of the TS module ("HYPE_MMF", "HYPE_HBVSED", or "NONE")
        ts_direction: Direction of exchange ("x", "y", or "z")
        ts_erosion_inhibit_compartment: Compartment that inhibits erosion (e.g., snow)
        ts_data_format: Data format for parameters ("JSON" or "ASCII")
        ts_mmf_defaults: Default parameter values for HYPE_MMF
        ts_mmf_parameters: Spatially-varying parameters for HYPE_MMF
        ts_hbvsed_defaults: Default parameter values for HYPE_HBVSED
        ts_hbvsed_parameters: Spatially-varying parameters for HYPE_HBVSED
        ts_hbvsed_monthly_erosion_factor: Monthly erosion factors (12 values, Jan-Dec)
    """

    if ts_module_name == "NONE":
        print("  Skipping TS module generation (ts_module_name='NONE')")
        return

    # Create the directory path if it doesn't exist
    output_path = Path(ts_config_filepath)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Build the configuration section
    configuration = {
        "direction": ts_direction,
        "EROSION_INHIBIT_COMPARTMENT": ts_erosion_inhibit_compartment,
        "Data_Format": ts_data_format
    }

    if ts_module_name == "HYPE_MMF":

        # Set defaults
        if ts_mmf_defaults is None:
            ts_mmf_defaults = {
                "COHESION": 7.5,
                "ERODIBILITY": 2.0,
                "SREROEXP": 1.2,
                "CROPCOVER": 0.5,
                "GROUNDCOVER": 0.5,
                "SLOPE": 0.4,
                "TRANSPORT_FACTOR_1": 0.2,
                "TRANSPORT_FACTOR_2": 0.5
            }

        if ts_mmf_parameters is None:
            ts_mmf_parameters = {}

        # Build the complete JSON structure
        config = {
            "MODULE_NAME": ts_module_name,
            "CONFIGURATION": configuration,
            "PARAMETER_DEFAULTS": ts_mmf_defaults,
            "PARAMETERS": ts_mmf_parameters
        }

    elif ts_module_name == "HYPE_HBVSED":

        # Set defaults
        if ts_hbvsed_defaults is None:
            ts_hbvsed_defaults = {
                "SLOPE": 0.4,
                "EROSION_INDEX": 0.4,
                "SOIL_EROSION_FACTOR_LAND_DEPENDENCE": 0.4,
                "SOIL_EROSION_FACTOR_SOIL_DEPENDENCE": 0.4,
                "SLOPE_EROSION_FACTOR_EXPONENT": 1.5,
                "PRECIP_EROSION_FACTOR_EXPONENT": 1.5,
                "PARAM_SCALING_EROSION_INDEX": 0.5
            }

        if ts_hbvsed_parameters is None:
            ts_hbvsed_parameters = {}

        if ts_hbvsed_monthly_erosion_factor is None:
            ts_hbvsed_monthly_erosion_factor = [0.0] * 12

        # Build the complete JSON structure
        config = {
            "MODULE_NAME": ts_module_name,
            "CONFIGURATION": configuration,
            "PARAMETER_DEFAULTS": ts_hbvsed_defaults,
            "MONTHLY_EROSION_FACTOR": ts_hbvsed_monthly_erosion_factor,
            "PARAMETERS": ts_hbvsed_parameters
        }

    else:
        print(f"WARNING: Unknown ts_module_name '{ts_module_name}'. "
              f"Valid options are: HYPE_MMF, HYPE_HBVSED, NONE")
        return

    # Convert to JSON string with standard formatting
    json_string = json.dumps(config, indent=4)

    # Compact arrays onto single lines (match format of other generated JSON files)
    def compact_list(match):
        list_content = match.group(0)
        compact = re.sub(r'\[\s+', '[', list_content)
        compact = re.sub(r'\s+\]', ']', compact)
        compact = re.sub(r',\s+', ', ', compact)
        return compact

    pattern = r'\[[\s\n]*(?:[^\[\]]*?)[\s\n]*\]'
    prev_string = ""
    while prev_string != json_string:
        prev_string = json_string
        json_string = re.sub(pattern, compact_list, json_string)

    # Write to file
    with open(ts_config_filepath, 'w') as f:
        # Write comment lines first
        for comment in json_header_comment:
            f.write(comment + "\n")
        # Write the JSON content
        f.write(json_string)

    print(f"  Transport Sediments config file saved to: {ts_config_filepath}")
