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
OpenWQ Transport Configuration Generator
Simple function to create transport JSON configuration from individual arguments
"""

import json
from pathlib import Path
from typing import List, Union


def create_td_module_json(

        td_config_filepath: str,

        json_header_comment: List[str],

        # Module name
        td_module_name: str,

        # Transport configuration parameters
        td_module_dispersion_xyz: List[Union[int, float]],

        # Characteristic length between cell centers [m]
        # Used to compute effective dispersion rate: D_eff = D_avg / L^2
        td_module_characteristic_length_m: Union[int, float] = 1.0

) -> None:
    """
    Create OpenWQ Transport configuration JSON file from individual arguments.

    Parameters:
        td_config_filepath: Full path where the JSON file will be saved (directories will be created if needed)
        json_header_comment: List of comment lines to add at the top of the file
        td_module_name: Name of the transport module (e.g., "OPENWQ_NATIVE_TD_ADVDISP")
        td_module_dispersion_xyz: Dispersion coefficients [Dx, Dy, Dz] in m2/s
        td_module_characteristic_length_m: Characteristic distance between cell centers [m].
            Used with OPENWQ_NATIVE_TD_ADVDISP to compute the effective dispersion rate
            D_eff = D_avg / L^2, where D_avg = (Dx + Dy + Dz) / 3.
    """

    dispersion_x = td_module_dispersion_xyz[0]
    dispersion_y = td_module_dispersion_xyz[1]
    dispersion_z = td_module_dispersion_xyz[2]

    # Create the directory path if it doesn't exist
    output_path = Path(td_config_filepath)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Build the complete JSON structure
    config = {
        "MODULE_NAME": td_module_name,
        "TRANSPORT_CONFIGURATION": {
            "dispersion_x_m2/s": dispersion_x,
            "dispersion_y_m2/s": dispersion_y,
            "dispersion_z_m2/s": dispersion_z,
            "characteristic_length_m": td_module_characteristic_length_m
        }
    }

    # Convert to JSON string with standard formatting
    json_string = json.dumps(config, indent=4)

    # Write to file
    with open(td_config_filepath, 'w') as f:
        # Write comment lines first
        for comment in json_header_comment:
            f.write(comment + "\n")
        # Write the JSON content
        f.write(json_string)

    print(f"  Transport config file saved to: {td_config_filepath}")