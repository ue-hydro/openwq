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
OpenWQ Sorption Isotherm Configuration Generator
Creates JSON configuration files for Freundlich or Langmuir sorption isotherms.
"""

import json
from pathlib import Path
from typing import Dict, List, Optional, Union


def create_si_module_json(

        si_config_filepath: str,

        json_header_comment: List[str],

        # Module name: "FREUNDLICH", "LANGMUIR", or "NONE"
        si_module_name: str,

        # Soil properties (shared by both isotherms)
        si_bulk_density_kg_m3: float = 1500.0,
        si_layer_thickness_m: float = 1.0,

        # Per-species isotherm parameters
        # For Freundlich: dict of species_name -> {"Kfr": float, "Nfr": float, "Kadsdes_1_per_s": float}
        # For Langmuir:   dict of species_name -> {"qmax_mg_per_kg": float, "KL_L_per_mg": float, "Kadsdes_1_per_s": float}
        si_species_params: Optional[Dict[str, Dict[str, float]]] = None

) -> None:
    """
    Create OpenWQ Sorption Isotherm configuration JSON file.

    Parameters:
        si_config_filepath: Full path where the JSON file will be saved
        json_header_comment: List of comment lines to add at the top
        si_module_name: "FREUNDLICH", "LANGMUIR", or "NONE"
        si_bulk_density_kg_m3: Bulk density of the soil/medium [kg/m3]
        si_layer_thickness_m: Representative layer thickness [m]
        si_species_params: Per-species isotherm parameters (see format above)
    """

    if si_module_name == "NONE":
        # Create minimal config file
        config = {
            "MODULE_NAME": "NONE"
        }
    elif si_module_name == "FREUNDLICH":
        species_block = {}
        if si_species_params:
            for species_name, params in si_species_params.items():
                species_block[species_name] = {
                    "Kfr": params.get("Kfr", 0.5),
                    "Nfr": params.get("Nfr", 0.7),
                    "Kadsdes_1/s": params.get("Kadsdes_1_per_s", 0.001)
                }

        config = {
            "MODULE_NAME": "FREUNDLICH",
            "SOIL_PROPERTIES": {
                "bulk_density_kg/m3": si_bulk_density_kg_m3,
                "layer_thickness_m": si_layer_thickness_m
            },
            "SPECIES": species_block
        }

    elif si_module_name == "LANGMUIR":
        species_block = {}
        if si_species_params:
            for species_name, params in si_species_params.items():
                species_block[species_name] = {
                    "qmax_mg/kg": params.get("qmax_mg_per_kg", 100.0),
                    "KL_L/mg": params.get("KL_L_per_mg", 0.01),
                    "Kadsdes_1/s": params.get("Kadsdes_1_per_s", 0.001)
                }

        config = {
            "MODULE_NAME": "LANGMUIR",
            "SOIL_PROPERTIES": {
                "bulk_density_kg/m3": si_bulk_density_kg_m3,
                "layer_thickness_m": si_layer_thickness_m
            },
            "SPECIES": species_block
        }
    else:
        print(f"WARNING: Unknown SI module name '{si_module_name}'. "
              f"Valid options: FREUNDLICH, LANGMUIR, NONE")
        return

    # Create the directory path if it doesn't exist
    output_path = Path(si_config_filepath)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Convert to JSON string with standard formatting
    json_string = json.dumps(config, indent=4)

    # Write to file
    with open(si_config_filepath, 'w') as f:
        # Write comment lines first
        for comment in json_header_comment:
            f.write(comment + "\n")
        # Write the JSON content
        f.write(json_string)

    print(f"  Sorption Isotherm config file saved to: {si_config_filepath}")
