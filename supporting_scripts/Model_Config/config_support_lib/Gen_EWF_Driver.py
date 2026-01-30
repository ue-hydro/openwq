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
OpenWQ External Input Flux Configuration Generator
Generate JSON configuration for external input fluxes
"""

import json
import re
from pathlib import Path
from typing import List, Union


def set_ewf_fixed_value(

        ewf_config_filepath: str,

        json_header_comment: List[str],

        # Metadata
        ewf_method_fixedval_comment: str,
        ewf_method_fixedval_source: str,

        # Configuration parameters
        ewf_method_fixedval_units: str,
        ewf_method_fixedval_chem_name: str,
        ewf_method_fixedval_value: Union[int, float],
        ewf_method_fixedval_external_inputflux_name: str

) -> None:
    """
    Create OpenWQ External Input Flux configuration JSON file.

    Parameters:
        ewf_config_filepath: Full path where the JSON file will be saved
        json_header_comment: List of comment lines to add at the top of the file
        ewf_method_fixedval_comment: Comment for metadata section
        ewf_method_fixedval_source: ewf_method_fixedval_source identifier for metadata section
        ewf_method_fixedval_units: Units (e.g., "mg/l")
        external_inputflux_name: Name of external input flux (e.g., "SUMMA_RUNOFF")
        ewf_method_fixedval_chem_name: Name of the chemical species (e.g., "NO3-N")
        ewf_method_fixedval_value: Numeric ewf_method_fixedval_value for the flux
    """

    # Create the destination directory path if it doesn't exist
    output_path = Path(ewf_config_filepath)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Build the complete JSON structure
    config = {
        "METADATA": {
            "Comment": ewf_method_fixedval_comment,
            "ewf_method_fixedval_source": ewf_method_fixedval_source
        },
        "1": {
            "Data_Format": "JSON",
            "Units": ewf_method_fixedval_units,
            "External_InputFlux_name": ewf_method_fixedval_external_inputflux_name,
            "ewf_method_fixedval_chem_name": ewf_method_fixedval_chem_name,
            "Data": {
                "1": ["all", "all", "all", "all", "all", "all", "all", "all", "all", ewf_method_fixedval_value]
            }
        }
    }

    # Convert to JSON string with standard formatting
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
    with open(ewf_config_filepath, 'w') as f:
        # Write comment lines first
        for comment_line in json_header_comment:
            f.write(comment_line + "\n")
        # Write the JSON content
        f.write(json_string)
        # Add newline at end of file
        f.write("\n")

    print(f"✓ External flux config file saved to: {ewf_config_filepath}")