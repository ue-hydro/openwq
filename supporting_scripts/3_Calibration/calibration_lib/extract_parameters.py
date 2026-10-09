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

try:                                   # works as a package module or standalone
    from . import ml_regionalization
except ImportError:                    # pragma: no cover
    import ml_regionalization

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

                # Regionalization (Layer 1): if this parameter carries an
                # ML_REGIONALIZE block, expand it into low-dimensional DDS
                # sub-parameters (per-class values or regression coefficients)
                # instead of a single global scalar. The parameter handler
                # regroups them per parameter and writes the per-cell
                # {"DEFAULT","CELLS"} map that openWQ's OpenWQ_load_param reads.
                ml_spec = info.get("ML_REGIONALIZE")
                if isinstance(ml_spec, dict):
                    json_path = ["CYCLING_FRAMEWORKS", framework_name, rxn_num,
                                 "PARAMETER_VALUES", param_name]
                    subs = ml_regionalization.expand_regionalized_param(
                        param_id, json_path, ml_spec)
                    _rung = ml_spec.get("rung", "per_class")
                    _attr = (ml_spec.get("attribute")
                             or ", ".join(ml_spec.get("attributes", []) or []))
                    for s in subs:
                        s["units"] = info.get("UNITS", "")
                        s["description"] = (
                            f"regionalized ({_rung}"
                            + (f" by {_attr}" if _attr else "")
                            + f") · '{s['subparam_key']}'"
                            + (f" — {info.get('DESCRIPTION','')}"
                               if info.get("DESCRIPTION") else ""))
                        s["regionalize_of"] = param_id     # parent physical param
                        s["_framework"] = framework_name
                        s["_reaction"] = rxn_name
                        s["_reaction_num"] = rxn_num
                    parameters.extend(subs)
                    logger.info(
                        f"Regionalized {param_id}: {len(subs)} sub-params "
                        f"(rung={ml_spec.get('rung', 'per_class')})")
                    continue

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

    # Module-level shared parameters: one value used by several reactions
    # (GLOBAL_PARAMETERS, read by the engine when a reaction's PARAMETER_VALUES
    # does not list the name). Calibration metadata sits in
    # _GLOBAL_PARAMETERS_INFO, mirroring _PARAMETERS_INFO of a reaction; only
    # the names listed there are extracted.
    global_values = data.get("GLOBAL_PARAMETERS", {}) or {}
    for param_name, info in (data.get("_GLOBAL_PARAMETERS_INFO", {}) or {}).items():
        if not isinstance(info, dict) or param_name not in global_values:
            continue
        value = float(global_values[param_name])
        param_range = info.get("RANGE")
        has_explicit_range = (isinstance(param_range, (list, tuple))
                              and len(param_range) == 2)
        if has_explicit_range:
            range_min, range_max = float(param_range[0]), float(param_range[1])
        else:
            abs_val = abs(value) if value != 0 else 0.01
            range_min, range_max = _round_sig(abs_val * 0.5), _round_sig(abs_val * 2.0)
        transform = ("log" if value > 0 and range_min > 0
                     and range_max / range_min > 100 else "linear")
        parameters.append({
            "name": f"GLOBAL_{param_name}",
            "file_type": "bgc_json",
            "path": ["GLOBAL_PARAMETERS", param_name],
            "initial": value,
            "bounds": (range_min, range_max),
            "transform": transform,
            "units": info.get("UNITS", ""),
            "description": info.get("DESCRIPTION", ""),
            "source": "auto-extracted",
            "has_explicit_range": has_explicit_range,
            "_framework": "",
            "_reaction": "",
            "_reaction_num": "",
        })

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


def make_sorption_regionalize_param(
    base_name: str,
    module: str,
    species: str,
    param_key: str,
    ml_spec: Dict,
    database_file: str = "openwq_in/SI_param_database.json"
) -> List[Dict]:
    """Declare a REGIONALIZED sorption parameter (ML Layer 1).

    Returns a LIST of low-dimensional DDS sub-parameters (per-class values or
    regression coefficients) — use ``parameters.extend(...)`` in the template.
    Each carries ``file_type="sorption_regionalize"`` and a ``path`` dict
    ``{module, species, param, database_file}``; the parameter handler groups
    them and writes a per-cell ``{"DEFAULT","CELLS"}`` map to the SI parameter
    database at ``[species][module][param_key]`` (which openWQ's SI loader,
    already routed through OpenWQ_load_param, reads as a spatial field).

    ``ml_spec`` is the same block as BGC regionalization, e.g.::

        {"rung":"per_class","attribute":"soil_class","default":50.0,
         "attribute_table":"openwq_in/attributes.csv",
         "mapping_source":"openwq_out/HDF5/...main.h5",
         "classes":{"sand":[10,80],"clay":[10,80],"peat":[30,120]}}
    """
    path = {"module": module, "species": species,
            "param": param_key, "database_file": database_file}
    return ml_regionalization.expand_regionalized_param(
        base_name, path, ml_spec, file_type="sorption_regionalize")


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


def make_ss_seasonal_param(
    name: str,
    species: str,
    file_type: str,
    initial: float,
    bounds: Tuple[float, float],
    month: int = None,
    description: str = ""
) -> Dict:
    """Build a seasonal source/sink scaling parameter.

    ``file_type`` is one of ``ss_seasonal_amp`` / ``ss_seasonal_phase``
    (harmonic mode) or ``ss_seasonal_month`` (per-month mode). The month-mode
    parameters carry their month (1-12) in ``path`` so the handler can scale
    only that month's load rows."""
    path = {"species": species}
    if month is not None:
        path["month"] = int(month)
    return {
        "name": name,
        "file_type": file_type,
        "path": path,
        "initial": initial,
        "bounds": bounds,
        "transform": "linear",
        "units": "multiplier" if file_type != "ss_seasonal_phase" else "months",
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


def apply_regionalize_to_params(params: List[Dict],
                                ml_regionalize: Optional[Dict]) -> List[Dict]:
    """Expand any calibration parameter named in ``ml_regionalize`` into its
    low-dimensional regionalization sub-parameters (hybrid-ML Layer 1A), in
    place of the single scalar. Parameters not listed pass through unchanged.

    The generated run script calls this (only when the report's ML tab activated
    regionalization) right after defining ``calibration_parameters``. Keys in
    ``ml_regionalize`` are the parameter ``name`` as shown in the report — which
    is model-tagged (``m{i}:``) in a chain; the un-tagged name is also accepted.
    Each expanded sub-parameter (per-class value / regression coefficient) is a
    normal low-dim DDS knob the parameter handler already knows how to apply.
    """
    if not ml_regionalize:
        return params

    _MODULE_FT = {"td": "module_regionalize", "le": "module_regionalize",
                  "ts": "ts_regionalize", "ss": "ss_regionalize"}  # else -> bgc_regionalize

    def _untag(nm):
        head, sep, rest = nm.partition(":")
        if sep and head.startswith("m") and head[1:].isdigit():
            return head + ":", rest, int(head[1:])
        return "", nm, None

    def _expand(name, spec, path, model_index=None, units=None):
        """Expand one regionalize entry into its DDS sub-params, tagged with the
        right file_type / module so the parameter handler writes to the correct
        module file. Works whether the target param is a BGC param (path from
        calibration_parameters) or a non-BGC module param (path baked in spec)."""
        tag, base_id, mi = _untag(name)
        if model_index is not None:
            mi = model_index
        subs = ml_regionalization.expand_regionalized_param(base_id, path, spec)
        module = spec.get("module", "bgc")
        ft = _MODULE_FT.get(module, "bgc_regionalize")
        _rung = spec.get("rung", "per_class")
        _attr = (spec.get("attribute")
                 or ", ".join(spec.get("attributes", []) or []))
        for s in subs:
            if tag:
                s["name"] = tag + s["name"]
            if mi is not None:
                s["model_index"] = mi
            if units:
                s.setdefault("units", units)
            s["file_type"] = ft
            if module in ("td", "le"):
                s["module_key"] = spec.get("module_key")
            if module == "ts":
                s["ts_param"] = spec.get("ts_param") or (path[-1] if path else None)
            s.setdefault("description",
                         f"regionalized ({_rung}"
                         + (f" by {_attr}" if _attr else "")
                         + f") · '{s.get('subparam_key')}'")
            s["regionalize_of"] = name
        return subs

    out: List[Dict] = []
    consumed = set()
    # Pass 1: BGC (and any) params ALREADY in calibration_parameters — replace
    # the lumped scalar with its regionalized sub-params (a param NOT listed here
    # stays a normal lumped scalar knob, so lumped calibration still works).
    for p in params:
        name = p.get("name", "")
        _, base, _ = _untag(name)
        key = name if name in ml_regionalize else (base if base in ml_regionalize else None)
        spec = ml_regionalize.get(key) if key else None
        if isinstance(spec, dict):
            consumed.add(key)
            path = spec.get("path") or p.get("path")
            subs = _expand(name, spec, path,
                           model_index=p.get("model_index"), units=p.get("units"))
            out.extend(subs)
            logger.info(f"Regionalized {name}: {len(subs)} sub-params")
        else:
            out.append(p)
    # Pass 2: non-BGC module/TS params (TD/LE/TS) that are not calibration
    # parameters — expand from the entry's own baked path/module.
    for key, spec in ml_regionalize.items():
        if key in consumed or not isinstance(spec, dict):
            continue
        subs = _expand(key, spec, spec.get("path"),
                       model_index=spec.get("model_index"))
        out.extend(subs)
        logger.info(f"Regionalized {key} (module): {len(subs)} sub-params")
    return out


# =========================================================================
# Layer-2 DDS-CALIBRATED per-species DERIVATIVE closures (no training / file)
# =========================================================================


def apply_calibrated_closures_to_params(params: List[Dict],
                                        ml_closures: Optional[Dict]) -> List[Dict]:
    """Expand every Layer-2 closure declared with ``mode == "calibrated"`` into
    its DDS sub-parameters — the weight ``cw`` and bias ``cb`` of a small bounded
    correction ``g = tanh(cw*x + cb)`` where ``x`` is the target species' own
    state (the solver passes a 1-vector, so ``n_in`` is always 1 — uniform for
    CHEM / SORPT / SS, no per-reaction width machinery). No offline training and
    no weights file: DDS fits the closure against the same observations as every
    other parameter, and the parameter handler serializes the calibrated
    ``(cw, cb)`` to the weights JSON openWQ reads each evaluation (see
    ``ParameterHandler._apply_calibrated_closures``).

    The generated run script calls this right after ``calibration_parameters`` is
    defined (only when the ML tab activated a calibrated closure). Closures with
    ``mode != "calibrated"`` (pretrained weights-file) pass through untouched.
    """
    if not ml_closures:
        return params
    out = list(params)
    for key, spec in ml_closures.items():
        if not isinstance(spec, dict) or spec.get("mode") != "calibrated":
            continue
        mi = spec.get("model_index", 0)
        _term = spec.get("term", "?")
        _sp = spec.get("species", "?")
        # cw (state sensitivity) and cb (offset) both start at 0 -> g=tanh(0)=0
        # -> factor=1 (exact physics) at the DDS start; the optimizer moves them.
        for sk, (lo, hi) in (("cw", (-5.0, 5.0)), ("cb", (-3.0, 3.0))):
            out.append({
                "name": f"{key}@{sk}",
                "file_type": "closure_calib",
                "path": [],
                "bounds": (lo, hi),
                "initial": 0.0,
                "transform": "linear",
                "closure_group": key,
                "subparam_key": sk,
                "model_index": mi,
                "source": "ml-closure-calib",
                "description": f"calibrated closure ({_term} on {_sp}) · '{sk}'",
            })
        logger.info(f"Calibrated closure {key}: 2 DDS sub-params (cw, cb)")
    return out


# =========================================================================
# Auto-extract parameters for ALL active modules
# =========================================================================

def extract_phreeqc_parameters(pqi_path: str) -> List[Dict]:
    """
    Extract calibration parameters from a PHREEQC input file (.pqi).

    A value is calibratable when its line carries a tag
    ``# CALIBRATE [min, max]``. Supported places:

    * CALCULATE_VALUES: the tag sits on the function-name line and the value is
      the number in its ``SAVE`` line, e.g.::

          CALCULATE_VALUES
          k_nit   # maximum nitrification rate [mg N/L/day]  # CALIBRATE [0.01, 20]
          -start
          10 SAVE 1.0
          -end

      RATES read it with CALC_VALUE("k_nit"), so one value is shared by every
      KINETICS block that uses the rate.
    * EXCHANGE: ``X  0.05   # CALIBRATE [0.001, 0.5]`` (moles of exchanger).
    * EQUILIBRIUM_PHASES: ``CO2(g)  -2.0  10   # CALIBRATE [-3.5, -1.0]`` (the
      target saturation index).
    """
    import re
    if not pqi_path or not os.path.isfile(pqi_path):
        logger.warning(f"PHREEQC input file not found: {pqi_path}")
        return []
    tag = re.compile(r"#\s*CALIBRATE\s*\[\s*([-+0-9.eE]+)\s*,\s*([-+0-9.eE]+)\s*\]")
    blocks = ("SOLUTION_MASTER_SPECIES", "SOLUTION_SPECIES", "EXCHANGE_SPECIES",
              "SURFACE_SPECIES", "PHASES", "CALCULATE_VALUES", "RATES", "SOLUTION",
              "EQUILIBRIUM_PHASES", "EXCHANGE", "SURFACE", "KINETICS", "GAS_PHASE",
              "SOLID_SOLUTIONS", "SELECTED_OUTPUT", "USER_PUNCH", "TITLE", "END")
    lines = open(pqi_path).read().splitlines()
    params, block = [], ""
    for i, raw in enumerate(lines):
        code = raw.split("#", 1)[0].strip()
        first = code.split()[0].upper() if code else ""
        if first in blocks:
            block = first
            continue
        m = tag.search(raw)
        if not m or not code:
            continue
        lo, hi = float(m.group(1)), float(m.group(2))
        desc = raw.split("#", 1)[1].split("CALIBRATE")[0].strip(" #") if "#" in raw else ""
        entry = None
        if block == "CALCULATE_VALUES":
            name = code.split()[0]
            value = None
            for nxt in lines[i + 1:]:
                s = nxt.split("#", 1)[0].strip()
                if s.lower() == "-end":
                    break
                mm = re.match(r"^\d+\s+SAVE\s+([-+0-9.eE]+)\s*$", s, re.I)
                if mm:
                    value = float(mm.group(1))
            if value is not None:
                entry = (f"PHREEQC_{name}", {"block": "CALCULATE_VALUES", "name": name}, value)
        elif block == "EXCHANGE":
            parts = code.split()
            if len(parts) >= 2:
                entry = (f"PHREEQC_EXCHANGE_{parts[0]}",
                         {"block": "EXCHANGE", "species": parts[0]}, float(parts[1]))
        elif block == "EQUILIBRIUM_PHASES":
            parts = code.split()
            if len(parts) >= 2:
                entry = (f"PHREEQC_SI_{parts[0]}",
                         {"block": "EQUILIBRIUM_PHASES", "phase": parts[0], "field": "si"},
                         float(parts[1]))
        if entry is None:
            logger.warning(f"CALIBRATE tag not understood in {os.path.basename(pqi_path)} "
                           f"line {i + 1}: {raw.strip()}")
            continue
        name, path, value = entry
        params.append({
            "name": name, "file_type": "phreeqc_pqi", "path": path,
            "initial": value, "bounds": (lo, hi),
            "transform": "log" if lo > 0 and hi / lo > 100 else "linear",
            "units": "", "description": desc, "source": "auto-extracted",
            "has_explicit_range": True,
            "_framework": "", "_reaction": "", "_reaction_num": "",
        })
    logger.info(f"Extracted {len(params)} calibration parameters from "
                f"{os.path.basename(pqi_path)}")
    return params


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
    elif model_config.get("bgc_module_name", "") == "PHREEQC":
        _pqi = extract_phreeqc_parameters(model_config.get("phreeqc_input_filepath", ""))
        if _pqi:
            groups["bgc"] = _pqi

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

    elif ss_method == "based_on_lulc":
        # ss_load_species contains the actual model species that have
        # loads (already resolved by the config template, including
        # stoichiometric conversions).
        #
        # Calibration knobs depend on the based_on_lulc sub-schema:
        #   lulc_loads='static'  (or 'dynamic'+shape='uniform') -> one constant
        #       SS_COP_scale (no within-year shape)
        #   lulc_loads='dynamic', lulc_loads_dynamic_shape='harmonic' ->
        #       SS_COP_scale (S0) + SS_COP_seasamp + SS_COP_seasphase
        #       factor(m)=S0*(1+A*cos(2pi(m-phi)/12))  [flat base]
        #   lulc_loads='dynamic', lulc_loads_dynamic_shape='monthly' ->
        #       12 free multipliers SS_COP_scale_mNN  [flat base]
        #   climate_dependency=True -> also calibrate the climate response
        #       coefficients (precip power, Q10, T_ref)
        #   species_dependency=False -> ONE shared knob set for all species
        #       (species='all' broadcasts); True -> a set per species.
        _lulc_loads = str(model_config.get("lulc_loads", "static")).lower()
        _shape = str(model_config.get("lulc_loads_dynamic_shape", "uniform")).lower()
        _climate_dep = bool(model_config.get("climate_dependency", False))
        _per_species = bool(model_config.get("species_dependency", True))
        _seas = _shape if _lulc_loads == "dynamic" else "uniform"

        def _ss_params_for(sp):
            sp_clean = sp.replace("-", "_").replace(" ", "_")
            if _seas == "monthly":
                _mp = [make_ss_seasonal_param(
                    f"SS_COP_scale_m{m:02d}_{sp_clean}", sp, "ss_seasonal_month",
                    initial=1.0, bounds=(0.001, 10.0), month=m,
                    description=f"Month-{m:02d} load multiplier for {sp}")
                    for m in range(1, 13)]
                # Log transform: multiplicative steps resolve near-zero loads.
                for _p in _mp:
                    _p["transform"] = "log"
                return _mp
            # Multiplier range kept physically plausible: a very large multiplier
            # both runs the model far slower (huge in-stream loads stiffen the
            # solver) and scores poorly, so it only wastes long evaluations.
            out = [make_ss_csv_scale_param(
                f"SS_COP_scale_{sp_clean}", species=sp, initial=1.0,
                bounds=(0.001, 10.0), description=f"Load scaling factor for {sp}")]
            # Log transform: multiplicative DDS steps so it can resolve loads
            # right down toward zero (linear steps over [0,10] are too coarse near 0).
            out[0]["transform"] = "log"
            if _seas == "harmonic":
                out.append(make_ss_seasonal_param(
                    f"SS_COP_seasamp_{sp_clean}", sp, "ss_seasonal_amp",
                    initial=0.0, bounds=(0.0, 0.95),
                    description=f"Seasonal amplitude (0=flat) for {sp}"))
                out.append(make_ss_seasonal_param(
                    f"SS_COP_seasphase_{sp_clean}", sp, "ss_seasonal_phase",
                    initial=0.0, bounds=(0.0, 12.0),
                    description=f"Seasonal phase offset (months) for {sp}"))
            return out

        ss_params = []
        if not _per_species:
            # species_dependency=False -> one shared set of knobs for ALL
            # species (species='all' broadcasts across every SS load row).
            ss_params.extend(_ss_params_for("all"))
        elif ss_load_species:
            for sp in sorted(ss_load_species):
                ss_params.extend(_ss_params_for(sp))
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
                                ss_params.extend(_ss_params_for(sp))

        # Climate response params (only when climate_dependency=True)
        is_dynamic = _climate_dep
        if is_dynamic:
            psp = float(model_config.get("ss_climate_precip_scaling_power", 1.0))
            q10 = float(model_config.get("ss_climate_temp_q10", 2.0))
            tref = float(model_config.get("ss_climate_temp_reference_c", 15.0))
            ss_params.extend([
                {"name": "SS_climate_precip_power", "file_type": "ss_climate",
                 "path": {"param": "precip_scaling_power"},
                 "initial": psp, "bounds": (0.1667, 5.0), "transform": "linear",
                 "units": "", "description": "Precipitation load exponent",
                 "source": "auto-extracted"},
                {"name": "SS_climate_Q10", "file_type": "ss_climate",
                 "path": {"param": "Q10_biological"},
                 "initial": q10, "bounds": (0.5, 8.0), "transform": "linear",
                 "units": "", "description": "Temperature Q10 factor",
                 "source": "auto-extracted"},
                {"name": "SS_climate_T_reference", "file_type": "ss_climate",
                 "path": {"param": "T_reference"},
                 "initial": tref, "bounds": (3.3333, 50.0), "transform": "linear",
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
