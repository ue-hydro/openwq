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
OpenWQ Lateral Exchange Configuration Generator
Simple function to create lateral exchange JSON configuration from individual arguments
"""

import json
from pathlib import Path
from typing import List, Union, Dict


def create_le_module_json(

        le_config_filepath: str,

        json_header_comment: List[str],

        # Module name
        le_module_name: str,

        # Lateral Exchange configuration parameters (list of le_module_config)
        le_module_config: List[Dict[str, Union[str, int, float]]]

) -> None:
    """
    Create OpenWQ Lateral Exchange configuration JSON file from individual arguments.

    Parameters:
        le_config_filepath: Full path where the JSON file will be saved (directories will be created if needed)
        json_header_comment: List of comment lines to add at the top of the file
        le_module_name: Name of the lateral exchange module (e.g., "NATIVE_LE_BOUNDMIX")
        le_module_config: List of configuration dictionaries, each containing:
            - direction: Direction of exchange (e.g., "z")
            - upper_compartment: Name of upper compartment
            - lower_compartment: Name of lower compartment
            - K_val: Exchange coefficient value
    """

    # Create the directory path if it doesn't exist
    output_path = Path(le_config_filepath)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Build the CONFIGURATION section with numbered entries
    configuration_dict = {}
    for idx, config in enumerate(le_module_config, start=1):
        configuration_dict[str(idx)] = {
            "direction": config["direction"],
            "upper_compartment": config["upper_compartment"],
            "lower_compartment": config["lower_compartment"],
            "K_val": config["K_val"]
        }

    # Build the complete JSON structure
    config = {
        "MODULE_NAME": le_module_name,
        "CONFIGURATION": configuration_dict
    }

    # Convert to JSON string with standard formatting
    json_string = json.dumps(config, indent=4)

    # Write to file
    with open(le_config_filepath, 'w') as f:
        # Write comment lines first
        for comment in json_header_comment:
            f.write(comment + "\n")
        # Write the JSON content
        f.write(json_string)

    print(f"✓ Lateral Exchange config file saved to: {le_config_filepath}")