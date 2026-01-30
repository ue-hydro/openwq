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
OpenWQ Source/Sink Configuration Generator
Generate JSON configuration for chemical sources and sinks
"""

import json
from pathlib import Path
from typing import List, Dict, Union


def set_ss_from_csv(

        ss_config_filepath: str,

        json_header_comment: List[str],

        # Metadata
        ss_method_csv_metadata_source: str,
        ss_method_csv_metadata_comment: str,

        # List of source/sink configurations
        ss_method_csv_config: List[Dict[str, Union[str, int]]]

) -> None:
    """
    Create OpenWQ Source/Sink configuration JSON file.

    Parameters:
        ss_config_filepath: Full path where the JSON file will be saved
        json_header_comment: List of comment lines to add at the top of the file
        metadata_comment: Comment for metadata section
        metadata_source: Source identifier for metadata section
        ss_method_csv_config: List of dictionaries, each containing:
            - Chemical_name: Name of the chemical species
            - Compartment_name: Name of the compartment
            - Type: "source" or "sink"
            - Units: Units (e.g., "kg")
            - Data_Format: Format type (e.g., "ASCII")
            - Filepath: Path to data file
            - Delimiter: CSV delimiter (e.g., ",")
            - Number_of_header_rows: Number of header rows
            - Header_key_row: Row number containing headers
    """

    # Create the destination directory path if it doesn't exist
    output_path = Path(ss_config_filepath)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Build the complete JSON structure
    config = {
        "METADATA": {
            "Comment": ss_method_csv_metadata_comment,
            "Source": ss_method_csv_metadata_source
        }
    }

    # Add each source/sink configuration with numbered keys
    for idx, ss_config in enumerate(ss_method_csv_config, start=1):
        config[str(idx)] = {
            "Chemical_name": ss_config["Chemical_name"],
            "Compartment_name": ss_config["Compartment_name"],
            "Type": ss_config["Type"],
            "Units": ss_config["Units"],
            "Data_Format": "ASCII",
            "Data": {
                "Filepath": ss_config["Filepath"],
                "Delimiter": ss_config["Delimiter"],
                "Number_of_header_rows": ss_config["Number_of_header_rows"],
                "Header_key_row": ss_config["Header_key_row"]
            }
        }

    # Convert to JSON string with standard formatting
    json_string = json.dumps(config, indent=4)

    # Write to file
    with open(ss_config_filepath, 'w') as f:
        # Write comment lines first
        for comment in json_header_comment:
            f.write(comment + "\n")
        # Write the JSON content
        f.write(json_string)
        # Add newline at end of file
        f.write("\n")

    print(f"✓ Source/Sink config file saved to: {ss_config_filepath}")