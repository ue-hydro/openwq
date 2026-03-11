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

"""
Parameter Extraction Module
============================

Auto-extracts calibration parameters from BGC template _PARAMETERS_INFO blocks.
All parameters are included; those with an explicit RANGE field use those bounds,
while others get auto-generated bounds based on their VALUE.

Also provides helper functions for building non-BGC parameter entries
(transport, sorption, sediment, lateral exchange, source/sink).
"""

import json
import os
import logging
from typing import List, Dict, Any, Optional, Tuple

logger = logging.getLogger(__name__)


def _round_sig(x: float, sig: int = 4) -> float:
    """Round *x* to *sig* significant figures to avoid float artefacts."""
    if x == 0:
        return 0.0
    from math import log10, floor
    return round(x, sig - 1 - int(floor(log10(abs(x)))))


def extract_calibration_parameters(bgc_template_path: str) -> List[Dict]:
    """
    Auto-extract calibration parameters from a NATIVE_BGC_FLEX template.

    Walks the BGC JSON structure:
        CYCLING_FRAMEWORKS -> framework -> reaction_number -> _PARAMETERS_INFO

    ALL parameters are extracted. Those with an explicit RANGE field use those
    bounds; others get auto-generated bounds based on their VALUE (+/- factor).
    Parameters without RANGE are marked ``has_explicit_range=False``.

    Parameters
    ----------
    bgc_template_path : str
        Path to the BGC template JSON file (e.g., SWAT_full_nutrients.json)

    Returns
    -------
    List[Dict]
        List of calibration parameter dicts, each with:
        - name: unique identifier (framework_reactionName_paramName)
        - file_type: "bgc_json"
        - path: list of JSON keys to reach PARAMETER_VALUES
        - initial: default value from the template
        - bounds: (min, max) tuple from RANGE or auto-generated
        - transform: "log" or "linear"
        - units: parameter units (if available)
        - description: parameter description (if available)
        - source: "auto-extracted"
        - has_explicit_range: True if RANGE was in the template
    """
    if not os.path.isfile(bgc_template_path):
        logger.warning(f"BGC template not found: {bgc_template_path}")
        return []

    # Read JSON (may have comment header lines)
    with open(bgc_template_path, 'r') as f:
        content = f.read()

    json_start = content.find('{')
    if json_start < 0:
        logger.warning(f"No JSON content found in {bgc_template_path}")
        return []

    try:
        data = json.loads(content[json_start:])
    except json.JSONDecodeError as e:
        logger.error(f"Failed to parse BGC template JSON: {e}")
        return []

    parameters = []

    cycling_frameworks = data.get("CYCLING_FRAMEWORKS", {})

    for framework_name, framework in cycling_frameworks.items():
        if not isinstance(framework, dict):
            continue

        for rxn_num, reaction in framework.items():
            if not isinstance(reaction, dict):
                continue

            # Get the reaction name for readable parameter naming
            rxn_name = reaction.get("_NAME", f"rxn{rxn_num}")
            # Clean up reaction name for use in parameter identifier
            rxn_name_clean = rxn_name.replace(" ", "_").replace("-", "_").upper()

            params_info = reaction.get("_PARAMETERS_INFO", {})
            if not isinstance(params_info, dict):
                continue

            for param_name, info in params_info.items():
                if not isinstance(info, dict):
                    continue

                value = float(info.get("VALUE", 0.0))

                # Use explicit RANGE if available, otherwise auto-generate
                param_range = info.get("RANGE")
                has_explicit_range = (
                    param_range is not None
                    and isinstance(param_range, (list, tuple))
                    and len(param_range) == 2
                )

                if has_explicit_range:
                    range_min = float(param_range[0])
                    range_max = float(param_range[1])
                else:
                    # Auto-generate bounds: ±50% of value, with floor
                    abs_val = abs(value) if value != 0 else 0.01
                    range_min = _round_sig(abs_val * 0.5)
                    range_max = _round_sig(abs_val * 2.0)
                    # Ensure min < max
                    if range_min >= range_max:
                        range_min, range_max = range_max, range_min
                    if range_min == range_max:
                        range_max = range_min * 2.0

                # Determine transform: use log if range spans >2 orders of
                # magnitude and value is positive
                if value > 0 and range_min > 0 and range_max / range_min > 100:
                    transform = "log"
                else:
                    transform = "linear"

                # Build unique parameter name
                param_id = f"{framework_name}_{rxn_name_clean}_{param_name}"

                param_entry = {
                    "name": param_id,
                    "file_type": "bgc_json",
                    "path": [
                        "CYCLING_FRAMEWORKS",
                        framework_name,
                        rxn_num,
                        "PARAMETER_VALUES",
                        param_name
                    ],
                    "initial": float(value),
                    "bounds": (range_min, range_max),
                    "transform": transform,
                    "units": info.get("UNITS", ""),
                    "description": info.get("DESCRIPTION", ""),
                    "source": "auto-extracted",
                    "has_explicit_range": has_explicit_range,
                    # Extra metadata for the report
                    "_framework": framework_name,
                    "_reaction": rxn_name,
                    "_reaction_num": rxn_num,
                }

                parameters.append(param_entry)
                logger.debug(
                    f"Extracted: {param_id} = {value} "
                    f"[{range_min}, {range_max}] ({transform})"
                )

    logger.info(
        f"Auto-extracted {len(parameters)} calibration parameters "
        f"from {os.path.basename(bgc_template_path)}"
    )
    return parameters


def apply_overrides(
    parameters: List[Dict],
    overrides: Dict[str, Optional[Dict]]
) -> List[Dict]:
    """
    Apply user overrides to auto-extracted parameters.

    Parameters
    ----------
    parameters : List[Dict]
        Auto-extracted parameter list
    overrides : Dict[str, Optional[Dict]]
        Mapping of parameter name -> override dict or None.
        - If value is a dict: merge into the parameter entry (e.g., new bounds)
        - If value is None: remove the parameter from calibration

    Returns
    -------
    List[Dict]
        Updated parameter list with overrides applied
    """
    result = []
    override_names = set(overrides.keys())
    applied = set()

    for param in parameters:
        name = param["name"]
        if name in overrides:
            applied.add(name)
            override = overrides[name]
            if override is None:
                # User wants to exclude this parameter
                logger.info(f"Excluded parameter: {name}")
                continue
            else:
                # Merge override into parameter
                merged = dict(param)
                merged.update(override)
                merged["source"] = "auto-extracted (user override)"
                result.append(merged)
                logger.info(f"Overridden parameter: {name} -> {override}")
        else:
            result.append(param)

    # Warn about overrides that didn't match any parameter
    unmatched = override_names - applied
    for name in unmatched:
        logger.warning(
            f"Override for '{name}' did not match any auto-extracted parameter"
        )

    return result


def merge_additional(
    parameters: List[Dict],
    additional: List[Dict]
) -> List[Dict]:
    """
    Merge additional user-defined parameters into the parameter list.

    Parameters
    ----------
    parameters : List[Dict]
        Existing parameter list (auto-extracted + overrides)
    additional : List[Dict]
        Additional parameter dicts provided by the user

    Returns
    -------
    List[Dict]
        Combined parameter list
    """
    existing_names = {p["name"] for p in parameters}

    for param in additional:
        name = param.get("name", "")
        if not name:
            logger.warning("Skipping additional parameter with no name")
            continue
        if name in existing_names:
            logger.warning(
                f"Additional parameter '{name}' conflicts with existing; skipping"
            )
            continue

        # Mark source
        param_copy = dict(param)
        if "source" not in param_copy:
            param_copy["source"] = "user-defined"
        parameters.append(param_copy)
        existing_names.add(name)

    return parameters


def print_parameter_table(parameters: List[Dict]) -> str:
    """
    Print a formatted table of calibration parameters.

    Parameters
    ----------
    parameters : List[Dict]
        Parameter list

    Returns
    -------
    str
        Formatted table string
    """
    if not parameters:
        return "No calibration parameters found."

    # Column widths
    name_w = max(len(p["name"]) for p in parameters)
    name_w = max(name_w, 4)  # min width for header

    lines = []
    lines.append("")
    lines.append("=" * (name_w + 80))
    lines.append("CALIBRATION PARAMETERS")
    lines.append("=" * (name_w + 80))
    lines.append(
        f"{'Name':<{name_w}}  {'Initial':>10}  {'Min':>10}  {'Max':>10}  "
        f"{'Transform':>9}  {'Units':>10}  {'Source':>16}"
    )
    lines.append("-" * (name_w + 80))

    for p in parameters:
        bounds = p.get("bounds", (0, 0))
        lines.append(
            f"{p['name']:<{name_w}}  {p.get('initial', 0):>10.4g}  "
            f"{bounds[0]:>10.4g}  {bounds[1]:>10.4g}  "
            f"{p.get('transform', 'linear'):>9}  "
            f"{p.get('units', ''):>10}  "
            f"{p.get('source', ''):>16}"
        )

    lines.append("-" * (name_w + 80))
    lines.append(f"Total: {len(parameters)} parameters")
    lines.append("")

    table_str = "\n".join(lines)
    print(table_str)
    return table_str


# =========================================================================
# Helper functions for building non-BGC parameter entries
# =========================================================================

def make_transport_param(
    name: str,
    param_key: str,
    initial: float,
    bounds: Tuple[float, float],
    transform: str = "log",
    units: str = "m2/s",
    description: str = ""
) -> Dict:
    """Build a transport dissolved parameter entry."""
    return {
        "name": name,
        "file_type": "transport_json",
        "path": {"param": param_key},
        "initial": initial,
        "bounds": bounds,
        "transform": transform,
        "units": units,
        "description": description,
        "source": "user-defined",
    }


def make_lateral_exchange_param(
    name: str,
    exchange_id: int,
    initial: float,
    bounds: Tuple[float, float],
    transform: str = "log",
    description: str = ""
) -> Dict:
    """Build a lateral exchange K_val parameter entry."""
    return {
        "name": name,
        "file_type": "lateral_exchange_json",
        "path": {"exchange_id": exchange_id, "param": "K_val"},
        "initial": initial,
        "bounds": bounds,
        "transform": transform,
        "units": "1/s",
        "description": description,
        "source": "user-defined",
    }


def make_sediment_param(
    name: str,
    module: str,
    param_key: str,
    initial: float,
    bounds: Tuple[float, float],
    transform: str = "linear",
    units: str = "",
    description: str = "",
    index: Optional[int] = None
) -> Dict:
    """Build a sediment transport parameter entry."""
    path = {"module": module, "param": param_key}
    if index is not None:
        path["index"] = index
    return {
        "name": name,
        "file_type": "sediment_json",
        "path": path,
        "initial": initial,
        "bounds": bounds,
        "transform": transform,
        "units": units,
        "description": description,
        "source": "user-defined",
    }


def make_sorption_param(
    name: str,
    module: str,
    species: Optional[str],
    param_key: str,
    initial: float,
    bounds: Tuple[float, float],
    transform: str = "log",
    units: str = "",
    description: str = ""
) -> Dict:
    """Build a sorption isotherm parameter entry."""
    path = {"module": module, "param": param_key}
    if species:
        path["species"] = species
    return {
        "name": name,
        "file_type": "sorption_json",
        "path": path,
        "initial": initial,
        "bounds": bounds,
        "transform": transform,
        "units": units,
        "description": description,
        "source": "user-defined",
    }


def make_ss_csv_scale_param(
    name: str,
    species: str = "all",
    initial: float = 1.0,
    bounds: Tuple[float, float] = (0.1, 5.0),
    description: str = "Source/sink load scaling factor"
) -> Dict:
    """Build a source/sink CSV scaling parameter entry."""
    return {
        "name": name,
        "file_type": "ss_csv_scale",
        "path": {"species": species},
        "initial": initial,
        "bounds": bounds,
        "transform": "linear",
        "units": "multiplier",
        "description": description,
        "source": "user-defined",
    }


def make_ss_copernicus_param(
    name: str,
    lulc_class: int,
    species: str,
    initial: float,
    bounds: Tuple[float, float],
    dynamic: bool = False,
    description: str = ""
) -> Dict:
    """Build a Copernicus LULC export coefficient parameter entry."""
    return {
        "name": name,
        "file_type": "ss_copernicus_dynamic" if dynamic else "ss_copernicus_static",
        "path": {"lulc_class": lulc_class, "species": species},
        "initial": initial,
        "bounds": bounds,
        "transform": "linear",
        "units": "kg/ha/yr",
        "description": description,
        "source": "user-defined",
    }


# =========================================================================
# Bounds lookup tables for auto-extraction
# =========================================================================

_HYPE_MMF_BOUNDS = {
    "COHESION": (1.0, 30.0, "linear", "kPa", "Soil cohesion"),
    "ERODIBILITY": (0.5, 10.0, "linear", "g/J", "Soil erodibility coefficient"),
    "SREROEXP": (0.5, 2.5, "linear", "", "Splash erosion exponent"),
    "CROPCOVER": (0.0, 1.0, "linear", "fraction", "Crop canopy cover fraction"),
    "GROUNDCOVER": (0.0, 1.0, "linear", "fraction", "Ground cover fraction"),
    "SLOPE": (0.001, 0.5, "linear", "m/m", "Terrain slope"),
    "TRANSPORT_FACTOR_1": (0.05, 0.5, "linear", "", "Sediment transport capacity factor 1"),
    "TRANSPORT_FACTOR_2": (0.1, 1.0, "linear", "", "Sediment transport capacity factor 2"),
}

_HYPE_HBVSED_BOUNDS = {
    "SLOPE": (0.001, 0.5, "linear", "m/m", "Terrain slope"),
    "EROSION_INDEX": (0.1, 2.0, "linear", "", "Soil erosion susceptibility"),
    "SOIL_EROSION_FACTOR_LAND_DEPENDENCE": (0.1, 2.0, "linear", "", "Land-use erosion factor"),
    "SOIL_EROSION_FACTOR_SOIL_DEPENDENCE": (0.1, 2.0, "linear", "", "Soil-type erosion factor"),
    "SLOPE_EROSION_FACTOR_EXPONENT": (0.5, 3.0, "linear", "", "Nonlinear slope effect"),
    "PRECIP_EROSION_FACTOR_EXPONENT": (0.5, 3.0, "linear", "", "Nonlinear rainfall effect"),
    "PARAM_SCALING_EROSION_INDEX": (0.1, 2.0, "linear", "", "Scaling factor for erosion index"),
}

_FREUNDLICH_BOUNDS = {
    "Kfr": (0.01, 100.0, "log", "L/kg", "Freundlich partition coefficient"),
    "Nfr": (0.3, 1.0, "linear", "", "Freundlich exponent"),
    "Kadsdes_1_per_s": (1e-6, 0.1, "log", "1/s", "Adsorption/desorption rate"),
}

_LANGMUIR_BOUNDS = {
    "qmax_mg_per_kg": (10.0, 5000.0, "log", "mg/kg", "Maximum sorption capacity"),
    "KL_L_per_mg": (0.001, 5.0, "log", "L/mg", "Langmuir affinity constant"),
    "Kadsdes_1_per_s": (1e-7, 0.01, "log", "1/s", "Adsorption/desorption rate"),
}

# LULC class names for Copernicus
_LULC_CLASS_NAMES = {
    10: "Cropland", 20: "Cropland_irrigated", 30: "Mosaic_crop",
    40: "Mosaic_natural", 50: "Broadleaf_evergreen", 60: "Broadleaf_deciduous",
    70: "Needleleaf_evergreen", 80: "Needleleaf_deciduous", 90: "Mixed_forest",
    100: "Mosaic_tree_shrub", 110: "Mosaic_herb", 120: "Shrubland",
    130: "Grassland", 140: "Lichens_mosses", 150: "Sparse_vegetation",
    160: "Freshwater_flood", 170: "Saltwater_flood", 180: "Shrub_herb_flood",
    190: "Urban", 200: "Bare", 210: "Water", 220: "Ice_snow",
}


# =========================================================================
# Auto-extract parameters for ALL active modules
# =========================================================================

def extract_all_module_parameters(
    model_config: Dict[str, Any],
    bgc_params: Optional[List[Dict]] = None,
    ss_load_species: Optional[set] = None,
) -> Dict[str, List[Dict]]:
    """
    Auto-extract calibration parameters for ALL active modules.

    Reads module selection variables from model_config and generates
    parameter entries for every active module. Only groups for active
    modules are included in the result.

    Parameters
    ----------
    model_config : Dict[str, Any]
        Loaded model configuration (from config_integration.load_model_config).
    bgc_params : List[Dict], optional
        Pre-extracted BGC parameters. If None, BGC group is skipped.
    ss_load_species : set, optional
        Set of model species names that have SS load entries
        (e.g. {"NH4-N", "NO3-N"}), as resolved by the config template.
        If None, falls back to custom coefficient species.

    Returns
    -------
    Dict[str, List[Dict]]
        Mapping of group_key -> list of parameter dicts.
        Keys: "bgc", "transport_dissolved", "sediment_transport",
              "lateral_exchange", "sorption_isotherm", "source_sink"
    """
    groups = {}

    # ── BGC parameters ──
    if bgc_params:
        groups["bgc"] = bgc_params

    # ── Transport Dissolved ──
    td_module = model_config.get("td_module_name", "NONE")
    if td_module == "OPENWQ_NATIVE_TD_ADVDISP":
        td_disp = model_config.get("td_module_dispersion_xyz", [1.0, 0.1, 0.001])
        char_len = model_config.get("td_module_characteristic_length_m", 100.0)
        td_params = [
            make_transport_param(
                "TD_dispersion_x", "dispersion_x",
                initial=float(td_disp[0]) if len(td_disp) > 0 else 1.0,
                bounds=(0.01, 100.0), transform="log", units="m2/s",
                description="Longitudinal dispersion coefficient"),
            make_transport_param(
                "TD_dispersion_y", "dispersion_y",
                initial=float(td_disp[1]) if len(td_disp) > 1 else 0.1,
                bounds=(0.001, 10.0), transform="log", units="m2/s",
                description="Transverse dispersion coefficient"),
            make_transport_param(
                "TD_dispersion_z", "dispersion_z",
                initial=float(td_disp[2]) if len(td_disp) > 2 else 0.001,
                bounds=(1e-5, 0.1), transform="log", units="m2/s",
                description="Vertical dispersion coefficient"),
            make_transport_param(
                "TD_characteristic_length", "characteristic_length",
                initial=float(char_len), bounds=(10.0, 10000.0),
                transform="log", units="m",
                description="Characteristic mixing length"),
        ]
        for p in td_params:
            p["source"] = "auto-extracted"
        groups["transport_dissolved"] = td_params

    # ── Sediment Transport ──
    ts_module = model_config.get("ts_module_name", "NONE")
    if ts_module == "HYPE_MMF":
        defaults = model_config.get("ts_mmf_defaults", {})
        bounds_tbl = _HYPE_MMF_BOUNDS
        ts_params = []
        for key, (bmin, bmax, transform, units, desc) in bounds_tbl.items():
            initial = float(defaults.get(key, (bmin + bmax) / 2))
            ts_params.append(make_sediment_param(
                f"TS_MMF_{key}", "HYPE_MMF", key,
                initial=initial, bounds=(bmin, bmax),
                transform=transform, units=units, description=desc))
        for p in ts_params:
            p["source"] = "auto-extracted"
        groups["sediment_transport"] = ts_params

    elif ts_module == "HYPE_HBVSED":
        defaults = model_config.get("ts_hbvsed_defaults", {})
        bounds_tbl = _HYPE_HBVSED_BOUNDS
        ts_params = []
        for key, (bmin, bmax, transform, units, desc) in bounds_tbl.items():
            initial = float(defaults.get(key, (bmin + bmax) / 2))
            ts_params.append(make_sediment_param(
                f"TS_HBVSED_{key}", "HYPE_HBVSED", key,
                initial=initial, bounds=(bmin, bmax),
                transform=transform, units=units, description=desc))
        for p in ts_params:
            p["source"] = "auto-extracted"
        groups["sediment_transport"] = ts_params

    # ── Lateral Exchange ──
    le_module = model_config.get("le_module_name", "NONE")
    if le_module == "NATIVE_LE_BOUNDMIX":
        le_config = model_config.get("le_module_config", [])
        if isinstance(le_config, list) and le_config:
            le_params = []
            for i, entry in enumerate(le_config):
                if not isinstance(entry, dict):
                    continue
                upper = entry.get("upper_compartment", f"COMP{i}")
                lower = entry.get("lower_compartment", f"COMP{i+1}")
                k_val = float(entry.get("K_val", 1e-9))
                name = f"LE_Kval_{upper}_{lower}"
                le_params.append(make_lateral_exchange_param(
                    name=name, exchange_id=i, initial=k_val,
                    bounds=(1e-14, 1e-6), transform="log",
                    description=f"Exchange: {upper} <-> {lower}"))
            for p in le_params:
                p["source"] = "auto-extracted"
            if le_params:
                groups["lateral_exchange"] = le_params

    # ── Sorption Isotherm ──
    si_module = model_config.get("si_module_name", "NONE")
    si_species_params = model_config.get("si_species_params", None)

    if si_module == "FREUNDLICH" and isinstance(si_species_params, dict):
        si_params = []
        for species, sp_vals in si_species_params.items():
            if not isinstance(sp_vals, dict):
                continue
            sp_clean = species.replace("-", "_").replace(" ", "_")
            for key, (bmin, bmax, transform, units, desc) in _FREUNDLICH_BOUNDS.items():
                initial = float(sp_vals.get(key, (bmin + bmax) / 2))
                si_params.append(make_sorption_param(
                    f"SI_FR_{sp_clean}_{key}", "FREUNDLICH", species, key,
                    initial=initial, bounds=(bmin, bmax),
                    transform=transform, units=units,
                    description=f"{desc} ({species})"))
        for p in si_params:
            p["source"] = "auto-extracted"
        if si_params:
            groups["sorption_isotherm"] = si_params

    elif si_module == "LANGMUIR" and isinstance(si_species_params, dict):
        si_params = []
        for species, sp_vals in si_species_params.items():
            if not isinstance(sp_vals, dict):
                continue
            sp_clean = species.replace("-", "_").replace(" ", "_")
            for key, (bmin, bmax, transform, units, desc) in _LANGMUIR_BOUNDS.items():
                initial = float(sp_vals.get(key, (bmin + bmax) / 2))
                si_params.append(make_sorption_param(
                    f"SI_LM_{sp_clean}_{key}", "LANGMUIR", species, key,
                    initial=initial, bounds=(bmin, bmax),
                    transform=transform, units=units,
                    description=f"{desc} ({species})"))
        for p in si_params:
            p["source"] = "auto-extracted"
        if si_params:
            groups["sorption_isotherm"] = si_params

    # ── Source/Sink Loads ──
    ss_method = model_config.get("ss_method", "none")
    species_list = model_config.get("chemical_species", [])
    if isinstance(species_list, str):
        species_list = []

    if ss_method == "load_from_csv":
        ss_params = [
            make_ss_csv_scale_param(
                "SS_CSV_scale_global", species="all",
                initial=1.0, bounds=(0.1, 5.0),
                description="Global source/sink load scaling factor"),
        ]
        for sp in species_list:
            sp_clean = sp.replace("-", "_").replace(" ", "_")
            ss_params.append(make_ss_csv_scale_param(
                f"SS_CSV_scale_{sp_clean}", species=sp,
                initial=1.0, bounds=(0.1, 5.0),
                description=f"Load scaling factor for {sp}"))
        for p in ss_params:
            p["source"] = "auto-extracted"
        groups["source_sink"] = ss_params

    elif ss_method in ("using_copernicus_lulc_with_static_coeff",
                       "using_copernicus_lulc_with_dynamic_coeff"):
        # ss_load_species contains the actual model species that have
        # loads (already resolved by the config template, including
        # stoichiometric conversions).  Generate a scaling factor param
        # per SS load species.
        ss_params = []
        if ss_load_species:
            for sp in sorted(ss_load_species):
                sp_clean = sp.replace("-", "_").replace(" ", "_")
                ss_params.append(make_ss_csv_scale_param(
                    f"SS_COP_scale_{sp_clean}",
                    species=sp,
                    initial=1.0,
                    bounds=(0.1, 5.0),
                    description=f"Load scaling factor for {sp}"))
        else:
            # Fallback: use custom coefficient species if ss_load_species
            # was not provided (SS JSON not generated yet)
            custom_coeffs = model_config.get(
                "ss_method_copernicus_optional_custom_annual_load_coeffs_per_lulc_class", {}
            )
            if isinstance(custom_coeffs, dict) and custom_coeffs:
                seen_species = set()
                for lulc_class, sp_loads in custom_coeffs.items():
                    if isinstance(sp_loads, dict):
                        for sp in sp_loads:
                            if sp not in seen_species:
                                seen_species.add(sp)
                                sp_clean = sp.replace("-", "_").replace(" ", "_")
                                ss_params.append(make_ss_csv_scale_param(
                                    f"SS_COP_scale_{sp_clean}",
                                    species=sp,
                                    initial=1.0,
                                    bounds=(0.1, 5.0),
                                    description=f"Load scaling factor for {sp}"))

        # Climate response params (dynamic SS only)
        is_dynamic = (ss_method == "using_copernicus_lulc_with_dynamic_coeff")
        if is_dynamic:
            psp = float(model_config.get("ss_climate_precip_scaling_power", 1.0))
            q10 = float(model_config.get("ss_climate_temp_q10", 2.0))
            tref = float(model_config.get("ss_climate_temp_reference_c", 15.0))
            ss_params.extend([
                {"name": "SS_climate_precip_power", "file_type": "ss_climate",
                 "path": {"param": "precip_scaling_power"},
                 "initial": psp, "bounds": (0.5, 2.5), "transform": "linear",
                 "units": "", "description": "Precipitation load exponent",
                 "source": "auto-extracted"},
                {"name": "SS_climate_Q10", "file_type": "ss_climate",
                 "path": {"param": "Q10_biological"},
                 "initial": q10, "bounds": (1.5, 4.0), "transform": "linear",
                 "units": "", "description": "Temperature Q10 factor",
                 "source": "auto-extracted"},
                {"name": "SS_climate_T_reference", "file_type": "ss_climate",
                 "path": {"param": "T_reference"},
                 "initial": tref, "bounds": (10.0, 25.0), "transform": "linear",
                 "units": "C", "description": "Reference temperature",
                 "source": "auto-extracted"},
            ])

        for p in ss_params:
            if "source" not in p:
                p["source"] = "auto-extracted"
        if ss_params:
            groups["source_sink"] = ss_params

    n_total = sum(len(v) for v in groups.values())
    logger.info(f"Extracted {n_total} parameters across {len(groups)} module groups")
    return groups
