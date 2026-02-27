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
OpenWQ Biogeochemistry Configuration Generator
Simple function to create biogeochemistry JSON configuration from individual arguments
"""

import json
import re
from pathlib import Path
from typing import Dict, List, Union


def create_config_json(

        config_file_fullpath: str,

        json_header_comment: List[str],

        # Compartment settings
        compartment_names: List[str],  # Changed to list of compartment names

        # Cycling framework (None for PHREEQC, list for NATIVE_BGC_FLEX)
        cycling_framework: Union[List[str], None],

        # Chemical species
        chemical_species_names: List[str],

        # Initial conditions (same for all species)
        ic_all_value: Union[int, float],
        ic_all_units: str,

) -> None:
    """
    Create OpenWQ Biogeochemistry configuration JSON file from individual arguments.

    Parameters:
        config_file_fullpath: Full path where the JSON file will be saved (directories will be created if needed)
        json_header_comment: List of comment lines to add at the top of the file
        compartment_names: List of compartment names (e.g., ["RIVER_NETWORK_REACHES", "SNOW_LAYER", "SOIL_LAYER"])
        cycling_framework: List of cycling frameworks (e.g., ["N_cycle"]) or None for PHREEQC
        chemical_species_names: List of chemical species names (e.g., ["NO3-N", "NH4-N", ...])
        ic_all_value: Initial value for all species (same value for all)
        ic_all_units: Unit for initial conditions (e.g., "mg/l")
    """

    # Create the directory path if it doesn't exist
    output_path = Path(config_file_fullpath)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Build the initial conditions data for each chemical species
    initial_conditions_data = {}

    for species in chemical_species_names:
        initial_conditions_data[species] = {
            "1": ["all", "all", "all", ic_all_value, ic_all_units]
        }

    # Build compartment configuration (same for all compartments)
    # PHREEQC: no CYCLING_FRAMEWORK needed (reactions defined in .pqi file)
    # NATIVE_BGC_FLEX: CYCLING_FRAMEWORK references the framework defined in the BGC module file
    compartment_config = {}
    if cycling_framework is not None:
        compartment_config["CYCLING_FRAMEWORK"] = cycling_framework

    compartment_config["INITIAL_CONDITIONS"] = {
        "Data_Format": "JSON",
        "Data": initial_conditions_data
    }

    # Build the complete JSON structure with all compartments
    biogeochemistry_config = {}
    for compartment_name in compartment_names:
        biogeochemistry_config[compartment_name] = compartment_config.copy()

    config = {
        "BIOGEOCHEMISTRY_CONFIGURATION": biogeochemistry_config
    }

    # Save to file with lists formatted horizontally

    # First convert to JSON string with standard formatting
    json_string = json.dumps(config, indent=4)

    # Regex pattern to find arrays/lists and format them on single lines
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

    # Write to file
    with open(config_file_fullpath, 'w') as f:
        # Write comment lines first
        for comment in json_header_comment:
            f.write(comment + "\n")
        # Write the JSON content
        f.write(json_string)

    print(f"✓ Config file saved to: {config_file_fullpath}")