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
OpenWQ Sorption Isotherm (model_SI) Configuration Generator.

Writes the SI module JSON (Freundlich or Langmuir). The module file holds only
the soil properties and a pointer to the sorption-parameter database
(SI_param_database.json). The per-species isotherm parameters live in that
database, INDEPENDENT of the BGC/model_CH templates, and the database is copied
next to the module file (into openwq_in/) so each run keeps its own editable
copy.
"""

import json
import shutil
from pathlib import Path
from typing import List


def create_si_module_json(

        si_config_filepath: str,

        json_header_comment: List[str],

        # Module name: "FREUNDLICH", "LANGMUIR", or "NONE"
        si_module_name: str,

        # Soil properties (shared by both isotherms)
        si_bulk_density_kg_m3: float = 1500.0,
        si_layer_thickness_m: float = 1.0,

        # Source path to the sorption-parameter database to ship with the run.
        # The caller (Gen_Input_Driver) defaults this to the shipped
        # config_support_lib/SI_model_parameters/SI_param_database.json; a user
        # may point it to a custom database with the same structure.
        si_param_database_filepath: str = ""

) -> None:
    """
    Create the OpenWQ Sorption Isotherm module JSON file and copy the
    sorption-parameter database into the run input folder (openwq_in/).

    Parameters:
        si_config_filepath: Full path where the SI module JSON will be saved
        json_header_comment: List of comment lines to add at the top
        si_module_name: "FREUNDLICH", "LANGMUIR", or "NONE"
        si_bulk_density_kg_m3: Bulk density of the soil/medium [kg/m3]
        si_layer_thickness_m: Representative layer thickness [m]
        si_param_database_filepath: Path to the SI parameter database
            (SI_param_database.json) to copy into the run and point to.
    """

    if si_module_name == "NONE":
        return

    if si_module_name not in ("FREUNDLICH", "LANGMUIR"):
        print(f"WARNING: Unknown SI module name '{si_module_name}'. "
              f"Valid options: FREUNDLICH, LANGMUIR, NONE")
        return

    # ------------------------------------------------------------------
    # Copy the parameter database next to the SI module JSON (openwq_in/)
    # and reference it with a path relative to the run CWD — mirroring how
    # the PHREEQC module ships its thermodynamic database.
    # ------------------------------------------------------------------
    db_source = Path(si_param_database_filepath)
    if not db_source.is_file():
        raise FileNotFoundError(
            f"SI parameter database not found: '{si_param_database_filepath}'.\n"
            f"  Expected the shipped config_support_lib/SI_model_parameters/"
            f"SI_param_database.json or a user-provided path with the same "
            f"structure."
        )

    out_dir = Path(si_config_filepath).parent
    out_dir.mkdir(parents=True, exist_ok=True)
    db_dest = out_dir / db_source.name
    if db_source.resolve() != db_dest.resolve():
        shutil.copy2(str(db_source), str(db_dest))
        print(f"  Copied SI parameter database to: {db_dest}")

    # Path embedded in the module JSON is relative to the run CWD (the master
    # JSON lives one level above openwq_in/). The key name contains
    # "DATABASE"/"FILEPATH", so openWQ preserves its case on load.
    db_relpath = f"openwq_in/{db_source.name}"

    config = {
        "MODULE_NAME": si_module_name,
        "SOIL_PROPERTIES": {
            "bulk_density_kg/m3": si_bulk_density_kg_m3,
            "layer_thickness_m": si_layer_thickness_m
        },
        "SI_PARAMETER_DATABASE_FILEPATH": db_relpath
    }

    # Convert to JSON string with standard formatting
    json_string = json.dumps(config, indent=4)

    # Write to file (comment lines first, then JSON content)
    with open(si_config_filepath, 'w') as f:
        for comment in json_header_comment:
            f.write(comment + "\n")
        f.write(json_string)

    print(f"  Sorption Isotherm config file saved to: {si_config_filepath}")
