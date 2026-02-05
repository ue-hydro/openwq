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
OpenWQ PHREEQC Module Configuration Generator

Creates the PHREEQC BGC module JSON file and copies the required
PHREEQC input file (.pqi) and database file (.dat) to the output directory.
"""

import json
import re
import shutil
from pathlib import Path
from typing import Dict, List, Optional, Union


def create_phreeqc_module_json(

        bgc_config_filepath: str,
        json_header_comment: List[str],

        # PHREEQC-specific settings
        phreeqc_input_filepath: str,
        phreeqc_database_filepath: str,
        phreeqc_mobile_species: List[str],
        phreeqc_component_h2o: bool = True,
        phreeqc_temperature_mapping: Optional[Dict[str, str]] = None,
        phreeqc_pressure_mapping: Optional[Dict[str, str]] = None,

) -> None:
    """
    Create OpenWQ PHREEQC module configuration JSON file and copy
    the required PHREEQC input (.pqi) and database (.dat) files
    to the output directory.

    Parameters:
        bgc_config_filepath: Full path where the PHREEQC module JSON will be saved
        json_header_comment: List of comment lines to add at the top of the file
        phreeqc_input_filepath: Path to the PHREEQC input file (.pqi)
        phreeqc_database_filepath: Path to the PHREEQC thermodynamic database (.dat)
        phreeqc_mobile_species: List of mobile species names (must match PHREEQC components)
        phreeqc_component_h2o: Whether to include H2O as a component (default True)
        phreeqc_temperature_mapping: Optional dict mapping compartment names to
            dependency variable names for temperature (e.g., {"RIVER_NETWORK_REACHES": "air_temperature"})
        phreeqc_pressure_mapping: Optional dict mapping compartment names to
            dependency variable names for pressure
    """

    # Create the destination directory path if it doesn't exist
    dest_path = Path(bgc_config_filepath)
    dest_dir = dest_path.parent
    dest_dir.mkdir(parents=True, exist_ok=True)

    # Validate source files exist
    pqi_source = Path(phreeqc_input_filepath)
    dat_source = Path(phreeqc_database_filepath)

    if not pqi_source.exists():
        raise FileNotFoundError(
            f"PHREEQC input file not found: {phreeqc_input_filepath}\n"
            f"Please provide a valid .pqi file with SOLUTION and reaction definitions."
        )

    if not dat_source.exists():
        raise FileNotFoundError(
            f"PHREEQC database file not found: {phreeqc_database_filepath}\n"
            f"Please download phreeqc.dat from https://www.usgs.gov/software/phreeqc-version-3"
        )

    # Copy PHREEQC files to the output directory
    pqi_dest = dest_dir / pqi_source.name
    dat_dest = dest_dir / dat_source.name

    if pqi_source.resolve() != pqi_dest.resolve():
        shutil.copy2(phreeqc_input_filepath, pqi_dest)
        print(f"  Copied PHREEQC input file to: {pqi_dest}")

    if dat_source.resolve() != dat_dest.resolve():
        shutil.copy2(phreeqc_database_filepath, dat_dest)
        print(f"  Copied PHREEQC database file to: {dat_dest}")

    # Build relative paths for JSON (relative to where master JSON lives,
    # which is one level above openwq_in/)
    pqi_relpath = f"openwq_in/{pqi_source.name}"
    dat_relpath = f"openwq_in/{dat_source.name}"

    # Build the JSON structure
    config = {
        "MODULE_NAME": "PHREEQC",
        "FILEPATH": pqi_relpath,
        "DATABASE": dat_relpath,
        "COMPONENT_H2O": phreeqc_component_h2o,
        "BGC_GENERAL_MOBILE_SPECIES": phreeqc_mobile_species
    }

    # Add optional temperature mapping
    if phreeqc_temperature_mapping:
        config["TEMPERATURE"] = phreeqc_temperature_mapping

    # Add optional pressure mapping
    if phreeqc_pressure_mapping:
        config["PRESSURE"] = phreeqc_pressure_mapping

    # Convert to formatted JSON string
    json_string = json.dumps(config, indent=4)

    # Compact arrays onto single lines
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
    with open(bgc_config_filepath, 'w') as f:
        for comment in json_header_comment:
            f.write(comment + "\n")
        f.write(json_string)

    print(f"  PHREEQC module config saved to: {bgc_config_filepath}")
