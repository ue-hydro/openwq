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

=============================================================================
PHREEQC SUPPORTED CHEMICAL SPECIES (phreeqc.dat database)
=============================================================================

MASTER SPECIES (Elements/Components):
-------------------------------------
These are the primary elements that PHREEQC tracks. When you define mobile
species in OpenWQ, use these element names (case-insensitive).

| Element    | Description                  | Default Species | Gram Formula Weight |
|------------|------------------------------|-----------------|---------------------|
| H          | Hydrogen                     | H+              | 1.008               |
| H(0)       | Hydrogen gas (H2)            | H2              | -                   |
| H(1)       | Hydrogen ion                 | H+              | -                   |
| O          | Oxygen                       | H2O             | 16.0                |
| O(0)       | Dissolved oxygen (O2)        | O2              | -                   |
| O(-2)      | Oxygen in water              | H2O             | -                   |
| Ca         | Calcium                      | Ca+2            | 40.08               |
| Mg         | Magnesium                    | Mg+2            | 24.312              |
| Na         | Sodium                       | Na+             | 22.9898             |
| K          | Potassium                    | K+              | 39.102              |
| Fe         | Iron                         | Fe+2            | 55.847              |
| Fe(+2)     | Ferrous iron                 | Fe+2            | -                   |
| Fe(+3)     | Ferric iron                  | Fe+3            | -                   |
| Mn         | Manganese                    | Mn+2            | 54.938              |
| Mn(+2)     | Manganous                    | Mn+2            | -                   |
| Mn(+3)     | Manganic                     | Mn+3            | -                   |
| Al         | Aluminum                     | Al+3            | 26.9815             |
| Ba         | Barium                       | Ba+2            | 137.34              |
| Sr         | Strontium                    | Sr+2            | 87.62               |
| Si         | Silicon (silica)             | H4SiO4          | 28.0843             |
| Cl         | Chloride                     | Cl-             | 35.453              |
| C          | Carbon (inorganic)           | CO3-2           | 12.0111             |
| C(+4)      | Carbonate carbon             | CO3-2           | -                   |
| C(-4)      | Methane carbon               | CH4             | -                   |
| S          | Sulfur                       | SO4-2           | 32.064              |
| S(6)       | Sulfate                      | SO4-2           | -                   |
| S(-2)      | Sulfide                      | HS-             | -                   |
| N          | Nitrogen                     | NO3-            | 14.0067             |
| N(+5)      | Nitrate nitrogen             | NO3-            | -                   |
| N(+3)      | Nitrite nitrogen             | NO2-            | -                   |
| N(0)       | Dissolved N2                 | N2              | -                   |
| N(-3)      | Ammonium nitrogen            | NH4+            | 14.0067             |
| B          | Boron                        | H3BO3           | 10.81               |
| P          | Phosphorus                   | PO4-3           | 30.9738             |
| F          | Fluoride                     | F-              | 18.9984             |
| Li         | Lithium                      | Li+             | 6.939               |
| Br         | Bromide                      | Br-             | 79.904              |
| Zn         | Zinc                         | Zn+2            | 65.37               |
| Cd         | Cadmium                      | Cd+2            | 112.4               |
| Pb         | Lead                         | Pb+2            | 207.19              |
| Cu         | Copper                       | Cu+2            | 63.546              |
| Cu(+2)     | Cupric copper                | Cu+2            | -                   |
| Cu(+1)     | Cuprous copper               | Cu+1            | -                   |
| Alkalinity | Alkalinity (as CaCO3)        | CO3-2           | 50.05               |
| Charge     | Electrical charge balance    | -               | -                   |

REDOX-UNCOUPLED GAS SPECIES:
----------------------------
These can be used for gas-water equilibrium without affecting other redox couples.

| Species | Description     | Formula Weight |
|---------|-----------------|----------------|
| Hdg     | Hydrogen gas    | 2.016          |
| Oxg     | Oxygen gas      | 32.0           |
| Mtg     | Methane gas     | 16.032         |
| Sg      | Hydrogen sulfide| 32.064         |
| Ntg     | Nitrogen gas    | 28.0134        |

MINERAL PHASES (for EQUILIBRIUM_PHASES):
----------------------------------------
Carbonate Minerals:
  - Calcite (CaCO3)
  - Aragonite (CaCO3)
  - Dolomite (CaMg(CO3)2)
  - Siderite (FeCO3)
  - Rhodochrosite (MnCO3)
  - Strontianite (SrCO3)
  - Witherite (BaCO3)
  - Cerussite (PbCO3)
  - Smithsonite (ZnCO3)
  - Otavite (CdCO3)

Sulfate Minerals:
  - Gypsum (CaSO4:2H2O)
  - Anhydrite (CaSO4)
  - Celestite (SrSO4)
  - Barite (BaSO4)
  - Anglesite (PbSO4)
  - Arcanite (K2SO4)
  - Mirabilite (Na2SO4:10H2O)
  - Thenardite (Na2SO4)
  - Epsomite (MgSO4:7H2O)
  - Melanterite (FeSO4:7H2O)

Silicate Minerals:
  - SiO2(a) (amorphous silica)
  - Chalcedony (SiO2)
  - Quartz (SiO2)
  - Gibbsite (Al(OH)3)
  - Kaolinite (Al2Si2O5(OH)4)
  - Albite (NaAlSi3O8)
  - Anorthite (CaAl2Si2O8)
  - K-feldspar (KAlSi3O8)
  - K-mica (muscovite)
  - Illite
  - Ca-Montmorillonite
  - Chlorite(14A)
  - Talc
  - Chrysotile
  - Sepiolite

Iron/Manganese Minerals:
  - Hematite (Fe2O3)
  - Goethite (FeOOH)
  - Fe(OH)3(a) (ferrihydrite)
  - Pyrite (FeS2)
  - FeS(ppt) (amorphous FeS)
  - Mackinawite (FeS)
  - Hausmannite (Mn3O4)
  - Manganite (MnOOH)
  - Pyrochroite (Mn(OH)2)

Other Minerals:
  - Hydroxyapatite (Ca5(PO4)3OH)
  - Vivianite (Fe3(PO4)2:8H2O)
  - Fluorite (CaF2)
  - Halite (NaCl)
  - Sylvite (KCl)
  - Sulfur (S)
  - Sphalerite (ZnS)
  - Willemite (Zn2SiO4)
  - Alunite (KAl3(SO4)2(OH)6)
  - Jarosite-K (KFe3(SO4)2(OH)6)

GAS PHASES (for EQUILIBRIUM_PHASES or GAS_PHASE):
-------------------------------------------------
  - CO2(g) - Carbon dioxide
  - O2(g) - Oxygen
  - H2(g) - Hydrogen
  - N2(g) - Nitrogen
  - H2S(g) - Hydrogen sulfide
  - CH4(g) - Methane
  - NH3(g) - Ammonia
  - H2O(g) - Water vapor

=============================================================================
USAGE NOTES:
=============================================================================

1. In BGC_GENERAL_MOBILE_SPECIES, use the ELEMENT names (e.g., "Ca", "N", "S"),
   not the species names (e.g., NOT "Ca+2", "NO3-", "SO4-2").

2. Species names are CASE-INSENSITIVE in OpenWQ. Both "Ca" and "CA" work.

3. For nitrogen species, use "N" as the mobile species. PHREEQC will track
   all nitrogen forms (NO3-, NO2-, NH4+, N2) internally.

4. For redox-sensitive elements (Fe, Mn, S, N, C), PHREEQC automatically
   handles speciation based on pe (redox potential).

5. Temperature is auto-converted: OpenWQ detects Kelvin (>100) and converts
   to Celsius for PHREEQC.

Example mobile species for water quality modeling:
  - Major ions: ["Ca", "Mg", "Na", "K", "Cl", "S", "C"]
  - Nutrients: ["N", "P"]
  - Metals: ["Fe", "Mn", "Zn", "Cu", "Pb"]
  - Full suite: ["Ca", "Mg", "Na", "K", "Cl", "S", "C", "N", "P", "Si", "Fe", "Mn"]

=============================================================================
"""

import json
import re
import shutil
from pathlib import Path
from typing import Dict, List, Optional, Union


# Complete list of PHREEQC master species from phreeqc.dat database
PHREEQC_MASTER_SPECIES = {
    # Primary elements
    "H": {"description": "Hydrogen", "default_species": "H+", "gfw": 1.008},
    "O": {"description": "Oxygen", "default_species": "H2O", "gfw": 16.0},
    "Ca": {"description": "Calcium", "default_species": "Ca+2", "gfw": 40.08},
    "Mg": {"description": "Magnesium", "default_species": "Mg+2", "gfw": 24.312},
    "Na": {"description": "Sodium", "default_species": "Na+", "gfw": 22.9898},
    "K": {"description": "Potassium", "default_species": "K+", "gfw": 39.102},
    "Fe": {"description": "Iron", "default_species": "Fe+2", "gfw": 55.847},
    "Mn": {"description": "Manganese", "default_species": "Mn+2", "gfw": 54.938},
    "Al": {"description": "Aluminum", "default_species": "Al+3", "gfw": 26.9815},
    "Ba": {"description": "Barium", "default_species": "Ba+2", "gfw": 137.34},
    "Sr": {"description": "Strontium", "default_species": "Sr+2", "gfw": 87.62},
    "Si": {"description": "Silicon (silica)", "default_species": "H4SiO4", "gfw": 28.0843},
    "Cl": {"description": "Chloride", "default_species": "Cl-", "gfw": 35.453},
    "C": {"description": "Carbon (inorganic)", "default_species": "CO3-2", "gfw": 12.0111},
    "S": {"description": "Sulfur", "default_species": "SO4-2", "gfw": 32.064},
    "N": {"description": "Nitrogen", "default_species": "NO3-", "gfw": 14.0067},
    "B": {"description": "Boron", "default_species": "H3BO3", "gfw": 10.81},
    "P": {"description": "Phosphorus", "default_species": "PO4-3", "gfw": 30.9738},
    "F": {"description": "Fluoride", "default_species": "F-", "gfw": 18.9984},
    "Li": {"description": "Lithium", "default_species": "Li+", "gfw": 6.939},
    "Br": {"description": "Bromide", "default_species": "Br-", "gfw": 79.904},
    "Zn": {"description": "Zinc", "default_species": "Zn+2", "gfw": 65.37},
    "Cd": {"description": "Cadmium", "default_species": "Cd+2", "gfw": 112.4},
    "Pb": {"description": "Lead", "default_species": "Pb+2", "gfw": 207.19},
    "Cu": {"description": "Copper", "default_species": "Cu+2", "gfw": 63.546},
}

# Redox states for elements with multiple oxidation states
PHREEQC_REDOX_STATES = {
    "H": ["H(0)", "H(1)"],
    "O": ["O(0)", "O(-2)"],
    "Fe": ["Fe(+2)", "Fe(+3)"],
    "Mn": ["Mn(+2)", "Mn(+3)"],
    "Cu": ["Cu(+1)", "Cu(+2)"],
    "C": ["C(+4)", "C(-4)"],
    "S": ["S(6)", "S(-2)"],
    "N": ["N(+5)", "N(+3)", "N(0)", "N(-3)"],
}

# Gas species (redox-uncoupled)
PHREEQC_GAS_SPECIES = {
    "Hdg": {"description": "Hydrogen gas (H2)", "gfw": 2.016},
    "Oxg": {"description": "Oxygen gas (O2)", "gfw": 32.0},
    "Mtg": {"description": "Methane gas (CH4)", "gfw": 16.032},
    "Sg": {"description": "Hydrogen sulfide gas (H2S)", "gfw": 32.064},
    "Ntg": {"description": "Nitrogen gas (N2)", "gfw": 28.0134},
}

# Common mineral phases
PHREEQC_MINERAL_PHASES = [
    # Carbonates
    "Calcite", "Aragonite", "Dolomite", "Siderite", "Rhodochrosite",
    "Strontianite", "Witherite", "Cerussite", "Smithsonite", "Otavite",
    # Sulfates
    "Gypsum", "Anhydrite", "Celestite", "Barite", "Anglesite",
    "Arcanite", "Mirabilite", "Thenardite", "Epsomite", "Melanterite",
    # Silicates
    "SiO2(a)", "Chalcedony", "Quartz", "Gibbsite", "Kaolinite",
    "Albite", "Anorthite", "K-feldspar", "K-mica", "Illite",
    "Ca-Montmorillonite", "Chlorite(14A)", "Talc", "Chrysotile", "Sepiolite",
    # Iron/Manganese
    "Hematite", "Goethite", "Fe(OH)3(a)", "Pyrite", "FeS(ppt)",
    "Mackinawite", "Hausmannite", "Manganite", "Pyrochroite",
    # Other
    "Hydroxyapatite", "Vivianite", "Fluorite", "Halite", "Sylvite",
    "Sulfur", "Sphalerite", "Willemite", "Alunite", "Jarosite-K",
]

# Gas phases for equilibrium
PHREEQC_GAS_PHASES = [
    "CO2(g)", "O2(g)", "H2(g)", "N2(g)", "H2S(g)", "CH4(g)", "NH3(g)", "H2O(g)",
]


def get_phreeqc_species_list() -> List[str]:
    """
    Return a list of all valid PHREEQC master species names.
    These are the element names that can be used in BGC_GENERAL_MOBILE_SPECIES.
    """
    return list(PHREEQC_MASTER_SPECIES.keys())


def get_phreeqc_species_info(species: str) -> Optional[Dict]:
    """
    Get information about a PHREEQC master species.

    Parameters:
        species: Element name (case-insensitive)

    Returns:
        Dictionary with description, default_species, and gfw, or None if not found
    """
    return PHREEQC_MASTER_SPECIES.get(species) or PHREEQC_MASTER_SPECIES.get(species.capitalize())


def validate_phreeqc_species(species_list: List[str]) -> List[str]:
    """
    Validate a list of species names against known PHREEQC master species.

    Parameters:
        species_list: List of species names to validate

    Returns:
        List of invalid species names (empty if all valid)
    """
    valid_species = set(s.lower() for s in PHREEQC_MASTER_SPECIES.keys())
    valid_species.add("charge")  # Special case
    invalid = []
    for species in species_list:
        if species.lower() not in valid_species:
            invalid.append(species)
    return invalid


def print_phreeqc_species_table():
    """Print a formatted table of all PHREEQC master species."""
    print("\nPHREEQC Master Species (phreeqc.dat database)")
    print("=" * 70)
    print(f"{'Element':<12} {'Description':<25} {'Default Species':<15} {'GFW':>10}")
    print("-" * 70)
    for element, info in PHREEQC_MASTER_SPECIES.items():
        print(f"{element:<12} {info['description']:<25} {info['default_species']:<15} {info['gfw']:>10.4f}")
    print("=" * 70)
    print("\nNote: Use element names (e.g., 'Ca', 'N', 'S') in BGC_GENERAL_MOBILE_SPECIES")
    print("Species names are case-insensitive in OpenWQ.\n")


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
            dependency variable names for temperature (e.g., {"RIVER_NETWORK_REACHES": "Treach_K"})
            NOTE: OpenWQ auto-detects Kelvin (>100) and converts to Celsius for PHREEQC.
            For MIZUROUTE, use "Treach_K" (Kelvin). For SUMMA, use "Tair_K" or "Tsoil_K".
        phreeqc_pressure_mapping: Optional dict mapping compartment names to
            dependency variable names for pressure (in atm)
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
