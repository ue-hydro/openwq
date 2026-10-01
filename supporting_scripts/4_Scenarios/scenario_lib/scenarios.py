# Copyright 2026, Diogo Costa, diogo.costa@uevora.pt
# This file is part of OpenWQ model.
#
# This program, openWQ, is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
Scenario analysis for openWQ — the lever catalogue and how each lever is
APPLIED to a run folder.

A *scenario* is a named list of *levers* (management practices, land-use
change, point sources, governance targets, climate change, in-stream measures,
initial stores).  The scenario runner (``run_scenarios``) builds one run folder
per scenario on top of the CALIBRATED parameter set (or the baseline when no
calibration has run), applies the levers to the openWQ input files (and, for
climate, to the host-model forcing), runs the coupled model and writes a
comparison report.

Every lever is expressed through the INPUT FILES the model reads — nothing in
the engine changes — so the same catalogue works for every host model, load
source and module combination.  Mechanisms:

* ``ss_factor``   per-unit × per-species × per-month multiplicative factors on
                  the source/sink loads (computed from the per-LULC-class load
                  breakdown the load generator writes to
                  ``ss_copernicus_files/nutrient_loads_detailed.csv``);
* ``ss_shift``    move loads between months (application windows / bans);
* ``ss_add``      append a new source/sink file (point sources);
* ``ss_entry``    scale an existing SS entry (point-source treatment / growth);
* ``ic_factor``   scale initial conditions (legacy stores);
* ``bgc_map``     a per-cell BGC parameter map on selected units (in-stream
                  measures) — the Layer-1 spatial-parameter machinery;
* ``forcing``     perturb the host forcing (precipitation / temperature) and
                  scale the loads with the generator's own climate response.

The catalogue (``LEVERS``) is the single source of truth for the setup
report's Scenarios tab (labels, parameters, defaults) and for the runner.
This package (``supporting_scripts/4_Scenarios/scenario_lib``) builds on the
calibration framework in ``../../3_Calibration/calibration_lib``.
"""
from __future__ import annotations

import os
import re
import json
import math
import copy
import glob
import shutil
import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Species groups (used for defaults; matched case-insensitively on the model's
# own species names)
# ---------------------------------------------------------------------------
_N_SPECIES = ("NO3-N", "NH4-N", "NO2-N", "TN", "ORG_N_ACTIVE", "ORG_N_FRESH", "DON", "TKN")
_P_SPECIES = ("PO4-P", "PO4-P_SOL", "TP", "ORG_P_HUMUS", "ORG_P_FRESH", "PO4-P_ACTIVE",
              "PO4-P_STABLE", "PP", "SRP", "DRP")
_NUTRIENT_SPECIES = _N_SPECIES + _P_SPECIES

# ESA-CCI / Copernicus land-cover code families (the default LULC source);
# other sources are namespaced (source_id*1e6 + native) and named by the
# generator — the report passes real class ids/names, these are only fallbacks.
_CROP_CODES = {10, 11, 12, 20, 30, 40}
_GRASS_CODES = {110, 120, 121, 122, 130, 140, 150, 151, 152, 153, 180}
_FOREST_CODES = {50, 60, 61, 62, 70, 71, 72, 80, 81, 82, 90, 100, 160, 170}
_URBAN_CODES = {190}


def _canon(s) -> str:
    return re.sub(r"[^A-Z0-9]", "", str(s).upper())


def resolve_model_species(model_config: Dict[str, Any]) -> List[str]:
    """The model's species as a LIST (``chemical_species`` may be ``"all"``);
    see ``calibration_lib.config_integration.resolve_model_species``."""
    from calibration_lib import config_integration as _ci
    return _ci.resolve_model_species(model_config)


_N_TOKENS = {"NO3", "NO2", "NH4", "NH3", "TN", "TKN", "DON", "DIN", "PON", "TON", "N"}
_P_TOKENS = {"PO4", "TP", "SRP", "DRP", "PP", "DOP", "DIP", "TDP", "P"}


def _species_group(sp) -> str:
    """'N' | 'P' | 'other' for a model species name (NO3-N, ORG_N_active,
    PO4-P_sol, PP, ... — matched on the name's tokens, then on the known
    lists). Gas-loss trackers (``*_loss``) are 'other'."""
    toks = [x for x in re.split(r"[^A-Za-z0-9]+", str(sp).upper()) if x]
    if "LOSS" in toks:
        return "other"
    if any(x in _N_TOKENS for x in toks):
        return "N"
    if any(x in _P_TOKENS for x in toks):
        return "P"
    c = _canon(sp)
    if c in {_canon(x) for x in _N_SPECIES} or c.startswith(("NITRAT", "NITRIT", "AMMONI", "NITROGEN")):
        return "N"
    if c in {_canon(x) for x in _P_SPECIES} or c.startswith(("PHOSPH", "ORTHOP")):
        return "P"
    return "other"


def _species_default(model_species, group="all"):
    """Model species of a default group: 'N', 'P' (exactly those — empty when
    the model has none) or 'all' (the nutrient species; every species when the
    model has no N / P species at all)."""
    ms = list(model_species or [])
    if group in ("N", "P"):
        return [s for s in ms if _species_group(s) == group]
    out = [s for s in ms if _species_group(s) in ("N", "P")]
    return out or ms


# ---------------------------------------------------------------------------
# The lever catalogue
# ---------------------------------------------------------------------------
# Two top-level GROUPS — "Climate" and "Management & policy" — each split into
# categories. Parameter kinds: number | percent | int | select | species |
# classes (LULC classes) | class | units (spatial units = HRU / reach ids:
# 'all', a list, or blank = the scenario's own units; selectable on the report
# map) | unit (one unit) | months | text | bool | param (a BGC parameter).
# Report hints: a parameter with "adv": True is shown under the option's
# "details"; a lever with "multi": True can be added several times
# ("add_label" names the button), the others are single tick-boxes.
GROUPS = [
    ("climate", "Climate"),
    ("management", "Management & policy"),
]

CATEGORIES = [   # (id, label, group, one-line hint shown in the report)
    ("climate", "Climate change", "climate", "delta-change on the forcing and the loads"),
    ("landuse", "Land-use change & disturbance", "management",
     "convert a share of one class into another; forest harvesting, fire"),
    ("agri", "Agricultural practices", "management",
     "fertilizer, tillage, cover crops, buffers, drainage, grazing, biochar"),
    ("urban_point", "Urban areas & point sources", "management",
     "stormwater BMPs, septic and sewer systems, point sources"),
    ("policy", "Policy targets & regulation", "management",
     "load targets, application limits, wastewater standards, deposition, detergents"),
    ("instream", "In-stream measures & legacy stores", "management",
     "stream restoration, in-stream processing rates, legacy nutrient stores"),
    ("generic", "Generic adjustment", "management", "a free load multiplier"),
]

# "Where": every management / policy option can be limited to a set of spatial
# units (HRUs / reaches). Blank = the units chosen for the whole scenario.
_UNITS_PARAM = {"key": "units", "label": "Where (unit ids; blank = the scenario's units)",
                "kind": "units", "default": ""}
# "When": an option can start in a given year (phased adoption, as SWAT's
# dated scheduled operations); blank = the whole run.
_START_YEAR_PARAM = {"key": "start_year", "label": "Start year (blank = whole run)", "kind": "int",
                     "default": "", "adv": True}


def _cr(classes: str, n: float, p: float, other: float = 0.0, extra=None) -> List[Dict[str, Any]]:
    """Parameters of a CLASS-REDUCTION practice: the export of the selected
    land-use classes is cut by a percentage per nutrient group (N species, P
    species, every other species), on a share of the class area (adoption), in
    the selected units. Negative values increase the export."""
    return [
        {"key": "classes", "label": "Land-use classes", "kind": "classes", "default": classes},
        {"key": "reduction_N_pct", "label": "N export reduction (%)", "kind": "percent", "default": n, "min": -100, "max": 100},
        {"key": "reduction_P_pct", "label": "P export reduction (%)", "kind": "percent", "default": p, "min": -100, "max": 100},
        {"key": "reduction_other_pct", "label": "Other species reduction (%)", "kind": "percent", "default": other,
         "min": -100, "max": 100, "adv": True},
        {"key": "adoption_pct", "label": "Adoption (% of the class area)", "kind": "percent", "default": 100, "min": 0, "max": 100},
    ] + [dict(e, adv=True) for e in (extra or [])] + [dict(_START_YEAR_PARAM), dict(_UNITS_PARAM)]


def _cr_inc(classes: str, n: float, p: float, other: float = 0.0) -> List[Dict[str, Any]]:
    """Parameters of a DISTURBANCE (``"sign": 1``): the export of the selected
    classes INCREASES by a percentage per nutrient group, on a share of the
    class area, in the selected units."""
    return [
        {"key": "classes", "label": "Land-use classes", "kind": "classes", "default": classes},
        {"key": "increase_N_pct", "label": "N export increase (%)", "kind": "percent", "default": n, "min": 0, "max": 2000},
        {"key": "increase_P_pct", "label": "P export increase (%)", "kind": "percent", "default": p, "min": 0, "max": 2000},
        {"key": "increase_other_pct", "label": "Other species increase (%)", "kind": "percent", "default": other,
         "min": 0, "max": 2000, "adv": True},
        {"key": "adoption_pct", "label": "Area affected (% of the class area)", "kind": "percent", "default": 100, "min": 0, "max": 100},
        dict(_START_YEAR_PARAM),
        dict(_UNITS_PARAM),
    ]


# Default efficiencies are MID-RANGE literature values (Iowa Nutrient Reduction
# Strategy science assessment; Chesapeake Bay Program BMP efficiencies; the
# SWAT conservation-practice guides) — illustrative, to be replaced by regional
# numbers. They act on the nutrient EXPORT of the selected classes.
LEVERS: List[Dict[str, Any]] = [
    # ═════════════════════════ CLIMATE ═════════════════════════════════════
    {"id": "climate_delta", "cat": "climate",
     "label": "Climate change (delta-change on forcing + load response)",
     "desc": "Perturb the host-model forcing (precipitation × factor, temperature + "
             "ΔT; SUMMA forcing NetCDF; for mizuRoute-only runs the runoff input is "
             "scaled by the precipitation factor as a proxy) and scale the loads with "
             "the load generator's own climate response "
             "(precip^power × Q10^(ΔT/10)). Annual or 12 monthly values. Presets are "
             "ILLUSTRATIVE global-mean CMIP6 ranges — replace with regional deltas. "
             "Under 'details': shortwave radiation, humidity and wind changes on the "
             "SUMMA forcing (SWAT's RADINC / HUMINC equivalents; no effect on a "
             "mizuRoute-only run). A CO2 change is not represented. Applies to the "
             "whole domain.",
     "ref": "Delta-change method (Hay et al. 2000); IPCC AR6 SSP ranges",
     "mech": "forcing",
     "params": [
         {"key": "dP_pct", "label": "Precipitation change (%, one value or 12 monthly)", "kind": "text", "default": "5"},
         {"key": "dT_c", "label": "Temperature change (°C, one value or 12 monthly)", "kind": "text", "default": "2.0"},
         {"adv": True, "key": "apply_forcing", "label": "Apply to host forcing (re-runs hydrology)", "kind": "bool", "default": True},
         {"adv": True, "key": "apply_loads", "label": "Apply to diffuse loads (climate response)", "kind": "bool", "default": True},
         {"adv": True, "key": "precip_power", "label": "Load precipitation exponent", "kind": "number", "default": 1.0, "min": 0, "max": 5},
         {"adv": True, "key": "q10", "label": "Load Q10 (temperature response)", "kind": "number", "default": 2.0, "min": 0.5, "max": 8},
         {"adv": True, "key": "dSW_pct", "label": "Shortwave radiation change (%)", "kind": "number", "default": 0.0, "min": -50, "max": 50},
         {"adv": True, "key": "dHum_pct", "label": "Specific humidity change (%)", "kind": "number", "default": 0.0, "min": -50, "max": 50},
         {"adv": True, "key": "dWind_pct", "label": "Wind speed change (%)", "kind": "number", "default": 0.0, "min": -50, "max": 50},
     ],
     "presets": {
         "SSP1-2.6 · 2050": {"dP_pct": "3", "dT_c": "1.0"},
         "SSP2-4.5 · 2050": {"dP_pct": "4", "dT_c": "1.5"},
         "SSP2-4.5 · 2080": {"dP_pct": "6", "dT_c": "2.2"},
         "SSP5-8.5 · 2050": {"dP_pct": "5", "dT_c": "2.0"},
         "SSP5-8.5 · 2080": {"dP_pct": "8", "dT_c": "3.7"},
     }},
    {"id": "precip_intensification", "cat": "climate",
     "label": "Precipitation intensification (more extreme events)",
     "desc": "Scale the wettest days (above a percentile of the forcing) by a factor "
             "and rescale the rest so the annual total is unchanged: same rain, more "
             "of it in storms. Host forcing only (SUMMA); loads unchanged. Applies to "
             "the whole domain.",
     "ref": "Clausius–Clapeyron scaling of extremes (~7 %/°C); Fowler et al. 2021",
     "mech": "forcing",
     "params": [
         {"key": "percentile", "label": "Wet-day percentile threshold", "kind": "number", "default": 95, "min": 50, "max": 99.9},
         {"key": "factor", "label": "Multiplier on days above it", "kind": "number", "default": 1.2, "min": 1, "max": 3},
     ]},

    # ═══════════════ MANAGEMENT & POLICY · Land-use change ═════════════════
    {"id": "landuse_change", "cat": "landuse", "multi": True, "add_label": "Add a land-use conversion",
     "label": "Land-use conversion (class A → class B)",
     "desc": "Convert a fraction of one land-use class into another in the selected "
             "units (afforestation, land retirement, agricultural expansion, "
             "urbanization). The unit's load is recomputed from the per-class export "
             "coefficients: the converted area exports at the destination class's "
             "coefficient.",
     "ref": "SWAT land-use update (lup.dat) / export-coefficient models",
     "mech": "ss_factor",
     "params": [
         {"key": "from_class", "label": "From class", "kind": "class", "default": "crop"},
         {"key": "to_class", "label": "To class", "kind": "class", "default": "forest"},
         {"key": "fraction_pct", "label": "Converted fraction (% of the 'from' area)", "kind": "percent", "default": 25, "min": 0, "max": 100},
         {"adv": True, "key": "species", "label": "Species", "kind": "species", "default": "all"},
         dict(_START_YEAR_PARAM),
         dict(_UNITS_PARAM),
     ],
     "presets": {
         "Afforestation (crop → forest)": {"from_class": "crop", "to_class": "forest", "fraction_pct": 25},
         "Reforest pasture (grass → forest)": {"from_class": "grass", "to_class": "forest", "fraction_pct": 25},
         "Land retirement (crop → grass)": {"from_class": "crop", "to_class": "grass", "fraction_pct": 20},
         "Agricultural expansion (forest → crop)": {"from_class": "forest", "to_class": "crop", "fraction_pct": 10},
         "Urbanization (crop → urban)": {"from_class": "crop", "to_class": "urban", "fraction_pct": 10},
     }},

    # ═══════════ MANAGEMENT & POLICY · Agricultural practices ══════════════
    {"id": "fert_reduction", "cat": "agri", "model": "class_reduction",
     "label": "Fertilizer / manure rate reduction (nutrient management plan)",
     "desc": "Lower fertilizer or manure application rate on the selected classes, "
             "expressed as a cut of their nutrient export (SWAT: fertilizer operation "
             "FRT_KG, auto-fertilization AUTO_NSTRS / AUTO_NAPP).",
     "ref": "SWAT management operations: fertilizer / auto-fertilization (Arnold et al. "
            "2012); Arabi et al. 2008",
     "mech": "ss_factor", "params": _cr("crop", 30, 30)},
    {"id": "fert_timing", "cat": "agri",
     "label": "Application timing / spreading window",
     "desc": "No application in the selected months (e.g. a winter spreading ban). "
             "The loads of those months are either moved to the first allowed month "
             "(same annual total, different timing) or removed.",
     "ref": "SWAT scheduled management (date-based operations); EU Nitrates Directive "
            "closed periods",
     "mech": "ss_shift",
     "params": [
         {"key": "classes", "label": "Land-use classes", "kind": "classes", "default": "crop"},
         {"adv": True, "key": "species", "label": "Species", "kind": "species", "default": "all"},
         {"key": "months", "label": "Closed months", "kind": "months", "default": [11, 12, 1, 2]},
         {"key": "mode", "label": "Closed-month loads", "kind": "select",
          "options": ["move_to_next_open_month", "remove"], "default": "move_to_next_open_month"},
         dict(_UNITS_PARAM),
     ]},
    {"id": "fert_placement", "cat": "agri", "model": "class_reduction",
     "label": "Fertilizer placement & enhanced-efficiency products (4R)",
     "desc": "Incorporation / injection instead of surface broadcast, nitrification or "
             "urease inhibitors, precision (variable-rate) application (SWAT: "
             "FRT_SURFACE, the fraction applied to the top 10 mm).",
     "ref": "SWAT fertilizer operation (FRT_SURFACE); 4R nutrient stewardship; Iowa "
            "Nutrient Reduction Strategy",
     "mech": "ss_factor", "params": _cr("crop", 10, 25),
     "presets": {
         "Incorporation / injection": {"reduction_N_pct": 10, "reduction_P_pct": 25},
         "Nitrification / urease inhibitors": {"reduction_N_pct": 15, "reduction_P_pct": 0},
         "Precision (variable-rate) application": {"reduction_N_pct": 10, "reduction_P_pct": 10},
     }},
    {"id": "manure_management", "cat": "agri", "model": "class_reduction",
     "label": "Manure management (storage, treatment, incorporation)",
     "desc": "Covered storage, no spreading on frozen or saturated soil, incorporation "
             "after spreading (SWAT: continuous fertilization / manure operations, "
             "CFRT_KG, MANURE_KG).",
     "ref": "SWAT continuous-fertilization and grazing-manure operations; Chesapeake "
            "Bay Program BMP efficiencies",
     "mech": "ss_factor", "params": _cr("crop", 20, 20)},
    {"id": "cover_crops", "cat": "agri", "model": "class_reduction",
     "label": "Cover crops / catch crops",
     "desc": "Winter cover takes up residual nitrogen and protects the soil: −30 to "
             "−50 % nitrate leaching, a smaller effect on P (SWAT: plant / kill "
             "operations of a cover crop in the rotation).",
     "ref": "SWAT plant-growth & management operations; Abdalla et al. 2019 "
            "meta-analysis; Iowa Nutrient Reduction Strategy",
     "mech": "ss_factor", "params": _cr("crop", 35, 10)},
    {"id": "crop_rotation", "cat": "agri", "model": "class_reduction",
     "label": "Crop rotation / diversification",
     "desc": "Longer rotations with legumes or perennials (less fertilizer N, more "
             "ground cover) instead of continuous row crops (SWAT: the plant / harvest "
             "schedule of the rotation).",
     "ref": "SWAT management schedule (.mgt rotations); Iowa Nutrient Reduction Strategy "
            "(extended rotations)",
     "mech": "ss_factor", "params": _cr("crop", 20, 10)},
    {"id": "conservation_tillage", "cat": "agri", "model": "class_reduction",
     "label": "Conservation tillage / no-till",
     "desc": "Less soil disturbance: lower sediment and particulate-P export (SWAT "
             "tillage operations: mixing efficiency, CN, USLE). Dissolved P can "
             "increase under no-till — use a negative value if relevant. Optionally "
             "scales the sediment erosion index of the transport module (whole domain).",
     "ref": "SWAT tillage operations (till.dat); Arabi et al. 2008",
     "mech": "ss_factor",
     "params": _cr("crop", 10, 30, extra=[
         {"key": "erosion_factor", "label": "Sediment erosion-index multiplier (TS module, whole domain)",
          "kind": "number", "default": 0.7, "min": 0, "max": 2}])},
    {"id": "residue_management", "cat": "agri", "model": "class_reduction",
     "label": "Residue management",
     "desc": "Crop residue left on the field after harvest protects the surface and "
             "slows runoff (SWAT scheduled operation: residue management).",
     "ref": "SWAT scheduled operations (.ops): residue management; Waidler et al. 2011",
     "mech": "ss_factor", "params": _cr("crop", 5, 20)},
    {"id": "contour_farming", "cat": "agri", "model": "class_reduction",
     "label": "Contour farming",
     "desc": "Tillage and planting along the contour: less runoff and erosion (SWAT "
             "scheduled operation: contouring, CONT_CN / CONT_P).",
     "ref": "SWAT scheduled operations (.ops): contouring; Arabi et al. 2008",
     "mech": "ss_factor", "params": _cr("crop", 10, 30)},
    {"id": "strip_cropping", "cat": "agri", "model": "class_reduction",
     "label": "Strip cropping",
     "desc": "Alternating strips of row crops and close-growing crops across the "
             "slope (SWAT scheduled operation: strip cropping, STRIP_N / STRIP_CN / "
             "STRIP_C / STRIP_P).",
     "ref": "SWAT scheduled operations (.ops): strip cropping; Waidler et al. 2011",
     "mech": "ss_factor", "params": _cr("crop", 20, 35)},
    {"id": "terracing", "cat": "agri", "model": "class_reduction",
     "label": "Terracing",
     "desc": "Terraces shorten the slope length and pond runoff (SWAT scheduled "
             "operation: terracing, TERR_P / TERR_CN / TERR_SL).",
     "ref": "SWAT scheduled operations (.ops): terracing; Arabi et al. 2008",
     "mech": "ss_factor", "params": _cr("crop", 20, 50)},
    {"id": "grassed_waterways", "cat": "agri", "model": "class_reduction",
     "label": "Grassed waterways",
     "desc": "Vegetated channels that carry concentrated runoff without gully erosion "
             "and trap sediment-bound nutrients (SWAT scheduled operation: grassed "
             "waterways, GWATN / GWATW / GWATL).",
     "ref": "SWAT scheduled operations (.ops): grassed waterways; Waidler et al. 2011",
     "mech": "ss_factor", "params": _cr("crop", 15, 35)},
    {"id": "filter_strip", "cat": "agri", "model": "class_reduction",
     "label": "Vegetative filter strips (edge of field)",
     "desc": "Grass strips at the field edge trap the exported load (SWAT VFS: "
             "VFSRATIO / VFSCON / VFSCH). Typical efficiencies: sediment 65–75 %, total "
             "N 30–50 %, total P 40–60 %, dissolved N 20–40 %.",
     "ref": "SWAT VFS sub-model (White & Arnold 2009); Liu et al. 2008 review",
     "mech": "ss_factor", "params": _cr("crop", 35, 50)},
    {"id": "riparian_buffer", "cat": "agri", "model": "class_reduction",
     "label": "Riparian buffer zones",
     "desc": "Forested or grassed buffers along the streams intercept surface and "
             "shallow subsurface flow from the upslope land: nitrate removal by uptake "
             "and denitrification, P by sedimentation. 'Adoption' is the share of the "
             "stream length that is buffered.",
     "ref": "Mayer et al. 2007 (N removal meta-analysis); Hoffmann et al. 2009 (P); SWAT "
            "VFS / REMM riparian model",
     "mech": "ss_factor", "params": _cr("all", 50, 45)},
    {"id": "drainage_water_management", "cat": "agri", "model": "class_reduction",
     "label": "Drainage water management (controlled drainage, bioreactors, saturated buffers)",
     "desc": "Measures on tile-drained land: raising the drain outlet, woodchip "
             "denitrifying bioreactors, saturated buffers. They act mainly on nitrate "
             "(SWAT tile drainage: DDRAIN / TDRAIN / GDRAIN, DRAINMOD routines).",
     "ref": "SWAT tile-drainage operations; Iowa Nutrient Reduction Strategy (controlled "
            "drainage −33 %, bioreactors −43 %, saturated buffers −50 % nitrate)",
     "mech": "ss_factor", "params": _cr("crop", 35, 5),
     "presets": {
         "Controlled drainage": {"reduction_N_pct": 33, "reduction_P_pct": 5},
         "Denitrifying bioreactor": {"reduction_N_pct": 43, "reduction_P_pct": 0},
         "Saturated buffer": {"reduction_N_pct": 50, "reduction_P_pct": 0},
     }},
    {"id": "irrigation_management", "cat": "agri", "model": "class_reduction",
     "label": "Irrigation management",
     "desc": "Efficient scheduling and application (deficit irrigation, drip instead "
             "of furrow) cut return flow and leaching (SWAT irrigation / auto-irrigation "
             "operations: IRR_AMT, IRR_EFM, AUTO_WSTRS).",
     "ref": "SWAT irrigation and auto-irrigation operations",
     "mech": "ss_factor", "params": _cr("crop", 15, 10)},
    {"id": "grazing_management", "cat": "agri", "model": "class_reduction",
     "label": "Grazing management / livestock exclusion",
     "desc": "Reduced stocking, rotational grazing or fencing livestock out of the "
             "streams on grassland / pasture classes (SWAT grazing operation: GRZ_DAYS, "
             "BIO_EAT, MANURE_KG).",
     "ref": "SWAT grazing operation; Chesapeake Bay Program BMP efficiencies",
     "mech": "ss_factor", "params": _cr("grass", 25, 25)},
    {"id": "biochar", "cat": "agri", "model": "class_reduction",
     "label": "Biochar soil amendment",
     "desc": "Biochar worked into cropland soil retains nitrate and ammonium: "
             "meta-analyses report about −13 % nitrate leaching on average and −26 % "
             "in experiments longer than a month. The P response is variable (available "
             "P can increase), hence 0 by default. Not a native SWAT operation: "
             "represented here as an export reduction.",
     "ref": "Borchard et al. 2019 (biochar meta-analysis: nitrate leaching and N2O)",
     "mech": "ss_factor", "params": _cr("crop", 15, 0)},
    {"id": "edge_of_field_retention", "cat": "agri", "model": "class_reduction",
     "label": "Constructed wetlands / retention ponds",
     "desc": "Intercept the load leaving the selected units before it reaches the "
             "water body (SWAT ponds and wetlands: .pnd). 'Adoption' is the share of "
             "the unit's drainage that passes through the wetland or pond.",
     "ref": "SWAT ponds & wetlands; Land et al. 2016 (wetland N / P removal review); "
            "Iowa Nutrient Reduction Strategy (wetlands −52 % nitrate)",
     "mech": "ss_factor", "params": _cr("all", 40, 40)},
    {"id": "conservation_practice_user", "cat": "agri", "model": "class_reduction",
     "label": "User-defined conservation practice (generic BMP)",
     "desc": "Any other practice, given directly as removal efficiencies per nutrient "
             "group (SWAT generic conservation practice: BMP_FLAG with BMP_SED / BMP_PP "
             "/ BMP_SP / BMP_PN / BMP_SN).",
     "ref": "SWAT scheduled operations (.ops): generic conservation practice",
     "mech": "ss_factor", "params": _cr("crop", 0, 0)},

    # ═════════ MANAGEMENT & POLICY · Urban areas & point sources ═══════════
    {"id": "urban_bmp", "cat": "urban_point", "model": "class_reduction",
     "label": "Urban stormwater BMPs (detention, LID, street sweeping)",
     "desc": "Reduced export from built-up classes through detention and wet ponds, "
             "bioretention / green infrastructure (rain gardens, green roofs, porous "
             "pavement) or street sweeping (SWAT urban BMPs and the sweep operation).",
     "ref": "SWAT urban BMP / LID modules; Jia et al. 2012",
     "mech": "ss_factor", "params": _cr("urban", 40, 40)},
    {"id": "point_source_add", "cat": "urban_point", "multi": True, "add_label": "Add a point source",
     "label": "New point source (WWTP, industry, aquaculture)",
     "desc": "A constant daily load at one spatial unit (reach / HRU id), optionally "
             "growing every year (population / production growth).",
     "ref": "SWAT point sources (.pnt / recday) ",
     "mech": "ss_add",
     "params": [
         {"key": "unit", "label": "Unit id (reach / HRU)", "kind": "unit", "default": ""},
         {"key": "species", "label": "Species", "kind": "species", "default": "all"},
         {"key": "load_kg_day", "label": "Load (kg/day per species)", "kind": "number", "default": 5.0, "min": 0},
         {"adv": True, "key": "growth_pct_yr", "label": "Growth (% per year)", "kind": "number", "default": 0.0, "min": -50, "max": 50},
         {"adv": True, "key": "start_year", "label": "Start year (blank = whole run)", "kind": "int", "default": ""},
     ]},
    {"id": "point_source_scale", "cat": "urban_point",
     "label": "Scale existing point-source loads (treatment upgrade / growth)",
     "desc": "Multiply the loads of the existing CSV / point-source entries of the "
             "master file (all entries except the land-use-based one) in the selected "
             "units — e.g. a tertiary-treatment upgrade (N −70 %, P −85 %) or "
             "population growth (+1.5 %/yr). Does nothing when the run has no such "
             "entries.",
     "ref": "EU Urban Waste Water Treatment Directive limits; SWAT point-source editing",
     "mech": "ss_entry",
     "params": [
         {"key": "species", "label": "Species", "kind": "species", "default": "all"},
         {"key": "factor", "label": "Multiplier (1 = unchanged)", "kind": "number", "default": 0.3, "min": 0},
         {"adv": True, "key": "growth_pct_yr", "label": "Additional growth (% per year)", "kind": "number", "default": 0.0, "min": -50, "max": 50},
         dict(_UNITS_PARAM),
     ]},

    # ═════════ MANAGEMENT & POLICY · Policy targets & regulation ═══════════
    {"id": "load_cap", "cat": "policy",
     "label": "Load-reduction target",
     "desc": "A policy target expressed as a uniform reduction of every diffuse load "
             "of the selected species in the selected units (all classes) — the 'what "
             "if the target is met' bookend, independent of how it is met.",
     "ref": "TMDL / Water Framework Directive programme-of-measures targets",
     "mech": "ss_factor",
     "params": [
         {"key": "species", "label": "Species", "kind": "species", "default": "all"},
         {"key": "reduction_pct", "label": "Reduction target (%)", "kind": "percent", "default": 20, "min": 0, "max": 100},
         dict(_START_YEAR_PARAM),
         dict(_UNITS_PARAM),
     ]},
    {"id": "application_rate_limit", "cat": "policy",
     "label": "Maximum application rate per hectare",
     "desc": "Cap the export coefficient of the selected classes at a maximum "
             "kg/ha/yr for each species (e.g. a 170 kg N/ha manure limit expressed "
             "as an export cap). Classes already below the cap are unchanged.",
     "ref": "EU Nitrates Directive 170 kg N/ha; provincial nutrient-management regulations",
     "mech": "ss_factor",
     "params": [
         {"key": "classes", "label": "Land-use classes", "kind": "classes", "default": "crop"},
         {"key": "species", "label": "Species", "kind": "species", "default": "N"},
         {"key": "max_kg_ha_yr", "label": "Maximum export (kg/ha/yr)", "kind": "number", "default": 10.0, "min": 0},
         dict(_UNITS_PARAM),
     ]},
    {"id": "wastewater_standard", "cat": "policy",
     "label": "Wastewater treatment standard",
     "desc": "An effluent standard or a treatment-class upgrade, as the required "
             "removal of N and of P from the existing point loads (CSV / point-source "
             "entries) in the selected units. MapShed / PRedICT, MONERIS and the MIKE "
             "Load Calculator step the treatment class (primary, secondary, tertiary); "
             "INCA and HYPE change the effluent concentrations; Balt-HYPE used 70–80 % N "
             "and 90 % P removal.",
     "ref": "UWWTD sensitive-area limits (N 70–80 %, P 80 % removal); Arheimer et al. 2012 (Balt-HYPE)",
     "mech": "ss_entry",
     "params": [
         {"key": "reduction_N_pct", "label": "N removal (%)", "kind": "percent", "default": 75, "min": 0, "max": 100},
         {"key": "reduction_P_pct", "label": "P removal (%)", "kind": "percent", "default": 80, "min": 0, "max": 100},
         dict(_UNITS_PARAM),
     ],
     "presets": {
         "Tertiary treatment (UWWTD sensitive areas)": {"reduction_N_pct": 75, "reduction_P_pct": 80},
         "Balt-HYPE wastewater scenario": {"reduction_N_pct": 75, "reduction_P_pct": 90},
         "Phosphorus stripping only": {"reduction_N_pct": 0, "reduction_P_pct": 80},
     }},

    # ═══════ MANAGEMENT & POLICY · In-stream measures & legacy stores ══════
    {"id": "instream_processing", "cat": "instream", "multi": True, "add_label": "Add an in-stream rate change",
     "label": "Enhanced in-stream processing (restoration, wetlands, beaver dams)",
     "desc": "Multiply one biogeochemical rate constant (e.g. denitrification k_den) "
             "in the selected units only — a per-cell parameter map through the "
             "Layer-1 spatial-parameter machinery. Other units keep the calibrated "
             "value.",
     "ref": "SWAT in-stream (QUAL2E) rate coefficients; restoration studies",
     "mech": "bgc_map",
     "params": [
         {"key": "parameter", "label": "Rate constant", "kind": "param", "default": ""},
         {"key": "factor", "label": "Multiplier", "kind": "number", "default": 2.0, "min": 0},
         dict(_UNITS_PARAM),
     ]},
    {"id": "legacy_stores", "cat": "instream",
     "label": "Legacy nutrient stores (initial conditions)",
     "desc": "Scale the initial concentration of the selected species in the "
             "selected compartments — legacy soil / groundwater nutrient stores "
             "that keep exporting after loads are cut. Applies to the whole domain.",
     "ref": "Van Meter et al. 2016 (legacy nitrogen)",
     "mech": "ic_factor",
     "params": [
         {"key": "species", "label": "Species", "kind": "species", "default": "all"},
         {"key": "compartments", "label": "Compartments (names, or 'all')", "kind": "text", "default": "all"},
         {"key": "factor", "label": "Multiplier", "kind": "number", "default": 2.0, "min": 0},
     ]},

    # ═══════════════ MANAGEMENT & POLICY · Generic adjustment ══════════════
    {"id": "load_scale", "cat": "generic", "multi": True, "add_label": "Add a load multiplier",
     "label": "Generic load multiplier",
     "desc": "Multiply the diffuse loads of the selected species by a factor — in "
             "the selected classes, units and months (atmospheric deposition trends, "
             "sensitivity bookends, anything not covered above).",
     "ref": "—",
     "mech": "ss_factor",
     "params": [
         {"key": "species", "label": "Species", "kind": "species", "default": "all"},
         {"key": "factor", "label": "Multiplier", "kind": "number", "default": 1.5, "min": 0},
         {"key": "classes", "label": "Land-use classes (or all)", "kind": "classes", "default": "all"},
         {"key": "months", "label": "Months (empty = all)", "kind": "months", "default": []},
         dict(_UNITS_PARAM),
     ]},
]

# ---------------------------------------------------------------------------
# Options added after a review of what other catchment models offer for
# scenario analysis (SWAT / SWAT+, HYPE, INCA, MIKE HYDRO Basin and MIKE SHE-
# DAISY, MONERIS, GWLF-E / MapShed, the Chesapeake Bay Program BMP list, HSPF,
# APEX, AnnAGNPS, SPARROW, the Iowa Nutrient Reduction Strategy).
# "INRS" = Iowa Nutrient Reduction Strategy science assessment (ISU Extension
# SP 435, 2013); "CBP" = Chesapeake Bay Program BMP efficiencies (Scenario
# Builder table, February 2011). A default with no such source is marked
# ILLUSTRATIVE in the description.
# ---------------------------------------------------------------------------
LEVERS += [
    # ── land use: disturbances (the export increases) ─────────────────────
    {"id": "forest_harvest", "cat": "landuse", "model": "class_reduction", "sign": 1,
     "label": "Forest harvesting (clear-felling)",
     "desc": "Harvesting raises the nutrient export of the felled area for some years "
             "(strongly on nitrogen-rich sites, hardly on boreal ones: in an INCA-N study "
             "a 20 % larger harvested area gave no change in the N load). 'Area' is the "
             "share of the forest classes harvested. Defaults are ILLUSTRATIVE. With "
             "harvesting BMPs the CBP credits −50 % N and −60 % P of the harvest load.",
     "ref": "INCA forestry scenarios (Rankinen et al. 2006); CBP forest harvesting practices",
     "mech": "ss_factor", "params": _cr_inc("forest", 50, 20),
     "presets": {"Clear-felling": {"increase_N_pct": 50, "increase_P_pct": 20},
                 "With harvesting BMPs (CBP)": {"increase_N_pct": 25, "increase_P_pct": 8}}},
    {"id": "wildfire", "cat": "landuse", "model": "class_reduction", "sign": 1,
     "label": "Wildfire / prescribed burning",
     "desc": "Burning removes cover and releases nutrients in ash: higher nutrient and "
             "sediment export from the burnt area for some years (SWAT: burn operation "
             "and the fire scheduled operation). 'Area' is the share of the classes "
             "burnt. Defaults are ILLUSTRATIVE.",
     "ref": "SWAT burn operation (BURN_FRLB) and scheduled fire operation (FIRE_CN)",
     "mech": "ss_factor", "params": _cr_inc("forest", 100, 150)},

    # ── agricultural practices ────────────────────────────────────────────
    {"id": "soil_test_p", "cat": "agri", "model": "class_reduction",
     "label": "Soil-test-based phosphorus application",
     "desc": "No P applied until the soil-test P has dropped to the agronomic optimum "
             "(INRS: −17 % P loss).",
     "ref": "INRS phosphorus practices",
     "mech": "ss_factor", "params": _cr("crop", 0, 17)},
    {"id": "pesticide_reduction", "cat": "agri", "model": "class_reduction",
     "label": "Pesticide / other contaminant application reduction",
     "desc": "Lower application of pesticides or of any other modelled substance that "
             "is neither N nor P (SWAT pesticide and continuous-pesticide operations). "
             "Acts on the species that are not nitrogen or phosphorus.",
     "ref": "SWAT pesticide operations (PEST_KG, CPST_KG)",
     "mech": "ss_factor", "table": False,     # shown as a line with its own fields, not in the N / P table
     "params": [dict(q, adv=(q["key"] in ("reduction_N_pct", "reduction_P_pct", "start_year")))
                for q in _cr("crop", 0, 0, other=30)]},
    {"id": "feed_management", "cat": "agri", "model": "class_reduction",
     "label": "Animal feed management (phytase, precision feeding)",
     "desc": "Phytase in poultry / pig feed and precision feeding of dairy cattle lower "
             "the nutrient content of manure. The CBP and MapShed / PRedICT credit it as "
             "an application reduction without a fixed efficiency: defaults are "
             "ILLUSTRATIVE.",
     "ref": "MapShed / PRedICT rural BMPs; CBP dairy precision feeding, poultry phytase",
     "mech": "ss_factor", "params": _cr("grass", 5, 15)},
    {"id": "tillage_timing", "cat": "agri", "model": "class_reduction",
     "label": "Tillage timing (spring instead of autumn ploughing)",
     "desc": "Ploughing in spring keeps the soil covered over winter and delays "
             "mineralisation (HYPE crop calendar: spring and autumn ploughing dates). No "
             "documented efficiency range was found: defaults are ILLUSTRATIVE.",
     "ref": "HYPE CropData (bd1 spring ploughing, bd4 autumn ploughing)",
     "mech": "ss_factor", "params": _cr("crop", 10, 0)},
    {"id": "low_input_farming", "cat": "agri", "model": "class_reduction",
     "label": "Organic / low-input farming",
     "desc": "Conversion to organic or low-input systems with lower nutrient surpluses "
             "(MONERIS evaluates it through the nitrogen surplus; MIKE SHE-DAISY through "
             "the crop rotation and fertilisation plan). Defaults are ILLUSTRATIVE.",
     "ref": "MONERIS nitrogen-surplus measures; MIKE SHE-DAISY agricultural management",
     "mech": "ss_factor", "params": _cr("crop", 20, 10)},
    {"id": "sediment_basins", "cat": "agri", "model": "class_reduction",
     "label": "Water and sediment control basins / check dams",
     "desc": "Small basins, check dams and grade-stabilisation structures in the field "
             "drainage pond the runoff and trap sediment-bound P (INRS: −85 % P for the "
             "area draining to the basin; APEX / SWAT: sediment basins, check dams).",
     "ref": "INRS sedimentation basins or ponds; Waidler et al. 2011 (APEX / SWAT)",
     "mech": "ss_factor", "params": _cr("crop", 0, 85)},
    {"id": "livestock_density", "cat": "agri", "model": "class_reduction",
     "label": "Livestock numbers (stocking density)",
     "desc": "Fewer animals per hectare: less manure and less grazing pressure (MIKE "
             "Load Calculator and MapShed: animal numbers; INCA: livestock input per "
             "land class). Defaults are ILLUSTRATIVE.",
     "ref": "MIKE HYDRO Basin Load Calculator; GWLF-E / MapShed farm animals; INCA-N",
     "mech": "ss_factor", "params": _cr("grass", 20, 20)},
    {"id": "barnyard_runoff", "cat": "agri", "model": "class_reduction",
     "label": "Barnyard / feedlot runoff control",
     "desc": "Roof gutters, diversions and treatment of the runoff from barnyards, "
             "loafing lots and confined feeding areas (CBP: −20 % N, −20 % P, −40 % "
             "sediment).",
     "ref": "CBP barnyard runoff control and loafing-lot management; MapShed; AnnAGNPS feedlots",
     "mech": "ss_factor", "params": _cr("grass", 20, 20, other=40)},
    {"id": "structure_liming", "cat": "agri", "model": "class_reduction",
     "label": "Structure liming of clay soils",
     "desc": "Quicklime or slaked lime worked into clay soils improves the aggregate "
             "stability and cuts P losses (Swedish programme of measures, which uses "
             "S-HYPE loads: −30 % P on clay soils).",
     "ref": "Vattenmyndigheterna 2016:19 (Swedish measures programme)",
     "mech": "ss_factor", "params": _cr("crop", 0, 30)},

    # ── urban areas & point sources ───────────────────────────────────────
    {"id": "urban_nutrient_mgmt", "cat": "urban_point", "model": "class_reduction",
     "label": "Urban nutrient management (lawn fertilizer)",
     "desc": "Restrictions and advice on lawn and turf fertilizer (CBP: −17 % N, "
             "−22 % P).",
     "ref": "CBP urban nutrient management",
     "mech": "ss_factor", "params": _cr("urban", 17, 22)},
    {"id": "construction_erosion", "cat": "urban_point", "model": "class_reduction",
     "label": "Construction-site erosion and sediment control",
     "desc": "Silt fences, sediment traps and stabilised exits on construction sites "
             "(CBP: −25 % N, −40 % P, −40 % sediment).",
     "ref": "CBP erosion and sediment control; Waidler et al. 2011 (APEX / SWAT)",
     "mech": "ss_factor", "params": _cr("urban", 25, 40, other=40)},
    {"id": "septic_upgrade", "cat": "urban_point", "model": "class_reduction",
     "label": "Septic systems: upgrade or connection to sewers",
     "desc": "On-site wastewater of unsewered households: denitrifying units (CBP: "
             "−50 % N), regular pumping (CBP: −5 % N), or connection to a sewer with "
             "treatment (SWAT septic systems; HYPE rural households; MONERIS and the "
             "MIKE Load Calculator: connection rates). 'Adoption' is the share of the "
             "unsewered households concerned.",
     "ref": "CBP septic BMPs; SWAT .sep; HYPE rural household sources; MONERIS",
     "mech": "ss_factor", "params": _cr("urban", 50, 0),
     "presets": {"Denitrifying units (CBP)": {"reduction_N_pct": 50, "reduction_P_pct": 0},
                 "Regular pumping (CBP)": {"reduction_N_pct": 5, "reduction_P_pct": 0}}},
    {"id": "sewer_measures", "cat": "urban_point", "model": "class_reduction",
     "label": "Sewer system: overflow storage and rehabilitation",
     "desc": "More storage for combined-sewer overflows and repair of leaking sewers "
             "(MONERIS urban-systems pathway). Sewer exfiltration was 9.8 % of the "
             "nitrate and 17.2 % of the phosphate urban-system loads in Germany, which "
             "bounds what rehabilitation can gain; no efficiency range is documented for "
             "overflow storage: defaults are ILLUSTRATIVE.",
     "ref": "MONERIS urban systems; Nguyen & Venohr 2021 (sewer exfiltration)",
     "mech": "ss_factor", "params": _cr("urban", 10, 17)},

    # ── policy targets & regulation ───────────────────────────────────────
    {"id": "atmospheric_deposition", "cat": "policy", "model": "class_reduction",
     "label": "Atmospheric nitrogen deposition change",
     "desc": "Air-quality policy changes the N deposited on every land-use class "
             "(HYPE AtmdepData; INCA wet and dry deposition; MONERIS and SPARROW source). "
             "Only part of the export responds: a 30 % cut in deposition gave about −7 % "
             "N load to the Baltic Sea, a 20 % cut about −6 % in an INCA-N study.",
     "ref": "Bartosova et al. 2019 (E-HYPE); Flynn et al. 2002 (INCA-N); MONERIS; SPARROW",
     "mech": "ss_factor", "params": _cr("all", 7, 0)},
    {"id": "p_free_detergents", "cat": "policy", "model": "class_reduction",
     "label": "Phosphate-free detergents",
     "desc": "A ban on phosphates in laundry and dishwasher detergents lowers the P in "
             "wastewater (MONERIS input; detergents were about 15 % of the P emissions "
             "of the Danube basin).",
     "ref": "MONERIS; ICPDR Danube nutrient assessments",
     "mech": "ss_factor", "params": _cr("urban", 0, 15)},

    # ── in-stream measures ────────────────────────────────────────────────
    {"id": "stream_restoration", "cat": "instream", "model": "class_reduction",
     "label": "Stream restoration / floodplain reconnection",
     "desc": "Re-meandering, floodplain reconnection and streambank stabilisation: more "
             "retention and less bank erosion. An E-HYPE study found −30 to −52 % P and "
             "−29 to −74 % inorganic N for floodplain 'stream mitigation' in two "
             "catchments; the CBP credits stream restoration per length of stream. "
             "'Adoption' is the share of the drainage that passes through the restored "
             "reaches.",
     "ref": "Wynants et al. 2024 (HYPE); CBP stream restoration; MapShed streambank stabilization",
     "mech": "ss_factor", "params": _cr("all", 30, 30)},
]


def _retune(lever_id: str, defaults: Dict[str, Any] = None, presets: Dict[str, Any] = None,
            label: str = None, desc_add: str = None, ref_add: str = None) -> None:
    """Defaults / presets of an option above, set from the reviewed sources."""
    lv = next(x for x in LEVERS if x["id"] == lever_id)
    for p in lv["params"]:
        if defaults and p["key"] in defaults:
            p["default"] = defaults[p["key"]]
    if presets:
        lv.setdefault("presets", {}).update(presets)
    if label:
        lv["label"] = label
    if desc_add:
        lv["desc"] = lv["desc"].rstrip() + " " + desc_add
    if ref_add:
        lv["ref"] = (lv["ref"].rstrip() + "; " + ref_add) if lv.get("ref") and lv["ref"] != "—" else ref_add


def _np(n, p):
    return {"reduction_N_pct": n, "reduction_P_pct": p}


_retune("cover_crops", _np(31, 10),
        desc_add="INRS: −31 % nitrate for a rye cover crop; CBP: −9 to −45 % N and −7 to −20 % P "
                 "depending on species, planting date and region.", ref_add="INRS; CBP")
_retune("crop_rotation", presets={"Living mulch (INRS)": _np(41, 0),
                                  "Perennial energy crops (INRS)": _np(72, 34)}, ref_add="INRS")
_retune("conservation_tillage",
        desc_add="CBP continuous no-till (pending approval): −10 to −15 % N, −20 to −40 % P, −70 % sediment.",
        ref_add="CBP")
_retune("filter_strip", _np(30, 40),
        desc_add="CBP grass buffers: −13 to −46 % N and −30 to −45 % P depending on the region.", ref_add="CBP")
_retune("riparian_buffer", _np(40, 40),
        desc_add="CBP forest buffers: −19 to −65 % N and −30 to −45 % P by region; MapShed default "
                 "−40 % N, −40 % P; Swedish programme of measures: −13 to −72 % P.",
        ref_add="CBP; MapShed; Vattenmyndigheterna 2016:19")
_retune("drainage_water_management", _np(33, 0),
        presets={"Controlled drainage (INRS, CBP)": _np(33, 0), "Denitrifying bioreactor (INRS)": _np(43, 0),
                 "Saturated buffer": _np(50, 0), "Shallow drainage (INRS)": _np(32, 0),
                 "Lime-filter drains (P)": _np(0, 25)},
        ref_add="INRS; CBP water control structures; Vattenmyndigheterna 2016:19 (lime-filter drains)")
_retune("edge_of_field_retention",
        presets={"Targeted wetland (INRS)": _np(52, 40), "Wetland restoration (CBP, mid-range)": _np(16, 30)},
        desc_add="INRS: −52 % nitrate for targeted wetlands; CBP wetland restoration: −7 to −25 % N and "
                 "−12 to −50 % P by region. HYPE represents them as N and P retention in wetland area.",
        ref_add="INRS; CBP; HYPE wetlands")
_retune("grazing_management", _np(10, 25),
        presets={"Prescribed grazing (CBP)": _np(9, 24), "Off-stream watering (CBP)": _np(5, 8)},
        desc_add="CBP: prescribed grazing −9 % N, −24 % P; off-stream watering −5 % N, −8 % P.", ref_add="CBP")
_retune("fert_placement", presets={"Manure / litter injection (CBP, interim)": _np(25, 0)}, ref_add="CBP")
_retune("manure_management",
        presets={"Animal waste management system (CBP: of the facility load)": _np(75, 75)},
        desc_add="CBP animal waste management systems: −75 % N and P of the load of the facility itself. "
                 "Manure or litter export out of the basin is credited as an application reduction.",
        ref_add="CBP; MapShed animal waste management systems")
_retune("biochar", _np(13, 0))
_retune("fert_reduction",
        desc_add="HYPE sets the fertilizer and manure amounts per crop (CropData); SPARROW and MONERIS "
                 "change the fertilizer source or the nitrogen surplus; Balt-HYPE used a flat −20 % arable "
                 "leaching scenario.", ref_add="HYPE CropData; MONERIS; SPARROW")
_retune("urban_bmp", _np(20, 45),
        presets={"Wet ponds and wetlands (CBP)": _np(20, 45), "Dry detention ponds (CBP)": _np(5, 10),
                 "Dry extended detention (CBP)": _np(20, 20), "Infiltration practices (CBP)": _np(80, 85),
                 "Filtering practices (CBP)": _np(40, 60), "Bioretention, C/D soils (CBP)": _np(25, 45),
                 "Bioretention, A/B soils (CBP)": _np(70, 75), "Vegetated open channels, A/B soils (CBP)": _np(45, 45),
                 "Bioswale (CBP)": _np(70, 75), "Permeable pavement, A/B soils (CBP)": _np(45, 50),
                 "Urban forest buffers (CBP)": _np(25, 50), "Street sweeping (CBP)": _np(3, 3)},
        ref_add="CBP urban BMP efficiencies")
_retune("landuse_change",
        presets={"Perennial energy crops (crop → grass)": {"from_class": "crop", "to_class": "grass", "fraction_pct": 20}})

_retune("conservation_practice_user",
        presets={"Diversion (SWAT table)": _np(10, 30), "Alum treatment (SWAT table)": _np(60, 80),
                 "Waste management system (SWAT table)": _np(80, 90),
                 "Solids separation basin (SWAT table)": _np(35, 31),
                 "Feedlot waste storage (SWAT table)": _np(65, 60),
                 "Pet waste management (STEPL)": _np(80, 90)},
        desc_add="The SWAT documentation tabulates typical removals for such practices: diversion 10 % N / "
                 "30 % P, alum treatment 60 / 80, waste management system 80 / 90, solids separation basin "
                 "35 / 31, feedlot waste storage 65 / 60.",
        ref_add="SWAT I/O documentation ch. 33 (conservation practice table)")
_retune("barnyard_runoff", presets={"Feedlot waste storage (SWAT table)": _np(65, 60),
                                    "Solids separation basin (SWAT table)": _np(35, 31)},
        ref_add="SWAT I/O documentation ch. 33")
_retune("septic_upgrade",
        presets={"Treatment technologies (CBP panel: 20–50 %)": _np(35, 0), "Connection to a sewer": _np(100, 100)},
        desc_add="CBP on-site wastewater panel: treatment technologies −20 to −50 % N, sewer connection "
                 "removes the on-site load.", ref_add="CBP OWTS Expert Panel 2014")
_retune("stream_restoration",
        presets={"Streambank stabilization and fencing (SWAT table)": _np(75, 75),
                 "Floodplain reconnection (Wynants et al. 2024, low end)": _np(30, 30)},
        desc_add="The SWAT documentation lists 75 % N, P and sediment for streambank stabilization and "
                 "fencing (of the bank-derived load).", ref_add="SWAT I/O documentation ch. 33")
_retune("drainage_water_management",
        desc_add="Installing NEW tile drains works the other way (more nitrate reaches the stream): enter "
                 "a negative N reduction for that case. SWAT sets the drain depth, time and lag "
                 "(DDRAIN / TDRAIN / GDRAIN) and can switch drains on at a date.")
_retune("residue_management",
        desc_add="Removing residue or stover (SWAT harvest operations: harvest efficiency, stover fraction) "
                 "works the other way: enter negative values.")
_retune("irrigation_management",
        desc_add="SWAT+ can also carry nitrate and phosphate in the irrigation water; add that load with a "
                 "point source or the load multiplier.")
_retune("pesticide_reduction",
        desc_add="Also for bacteria / pathogens carried by manure (SWAT bacteria options).")

# Sub-groups of the agricultural practices (headings of the practices table)
# and the display order of every option inside its category.
_AGRI_SUBS = [
    ("Nutrient & manure management", ["fert_reduction", "soil_test_p", "fert_timing", "fert_placement",
                                      "manure_management", "feed_management"]),
    ("Crop & soil management", ["cover_crops", "crop_rotation", "conservation_tillage", "tillage_timing",
                                "residue_management", "irrigation_management", "low_input_farming"]),
    ("Erosion & runoff control", ["contour_farming", "strip_cropping", "terracing", "grassed_waterways",
                                  "sediment_basins"]),
    ("Edge of field & drainage", ["filter_strip", "riparian_buffer", "drainage_water_management",
                                  "edge_of_field_retention"]),
    ("Livestock", ["grazing_management", "livestock_density", "barnyard_runoff"]),
    ("Soil amendments", ["biochar", "structure_liming"]),
    ("Other", ["conservation_practice_user", "pesticide_reduction"]),
]
_ORDER = ["climate_delta", "precip_intensification",
          "landuse_change", "forest_harvest", "wildfire"]
for _sub, _ids in _AGRI_SUBS:
    for _i in _ids:
        next(x for x in LEVERS if x["id"] == _i)["sub"] = _sub
    _ORDER += _ids
_ORDER += ["urban_bmp", "urban_nutrient_mgmt", "construction_erosion", "septic_upgrade", "sewer_measures",
           "point_source_scale", "point_source_add",
           "load_cap", "application_rate_limit", "wastewater_standard", "atmospheric_deposition",
           "p_free_detergents",
           "stream_restoration", "instream_processing", "legacy_stores", "load_scale"]
_CAT_IX = {c[0]: i for i, c in enumerate(CATEGORIES)}
LEVERS.sort(key=lambda lv: (_CAT_IX.get(lv["cat"], 99),
                            _ORDER.index(lv["id"]) if lv["id"] in _ORDER else len(_ORDER)))

LEVER_BY_ID = {lv["id"]: lv for lv in LEVERS}
_CAT_GROUP = {c[0]: c[2] for c in CATEGORIES}


def lever_group(lever_id: str) -> str:
    """'climate' | 'management' for a lever id."""
    lv = LEVER_BY_ID.get(lever_id) or {}
    return _CAT_GROUP.get(lv.get("cat"), "management")


def lever_spatial(lever: Dict[str, Any]) -> str:
    """'units' (a set of units), 'unit' (one unit) or '' (whole domain)."""
    kinds = {p.get("kind") for p in lever.get("params") or []}
    return "units" if "units" in kinds else ("unit" if "unit" in kinds else "")


# ---------------------------------------------------------------------------
# LULC / load context of a run (from the generator's intermediate files)
# ---------------------------------------------------------------------------
def _read_json_hdr(path):
    txt = open(path).read()
    hdr = [l for l in txt.splitlines() if l.strip().startswith("//")]
    body = "\n".join(l for l in txt.splitlines() if not l.strip().startswith("//"))
    return json.loads(body), hdr


def _write_json_hdr(path, data, hdr):
    with open(path, "w") as f:
        if hdr:
            f.write("\n".join(hdr) + "\n")
        f.write(json.dumps(data, indent=2) + "\n")


def load_lulc_context(run_dir: str) -> Dict[str, Any]:
    """Per-unit × per-class × per-species annual loads and class areas of a run
    (``<run_dir>/openwq_in/ss_copernicus_files/``). Returns::

        {"loads": {unit_id: {species: {class_id: kg_yr}}},   # one representative year
         "areas": {unit_id: {class_id: ha}},
         "coef":  {species: {class_id: kg_ha_yr}},
         "classes": {class_id: {"name": str, "area_ha": total}},
         "years": [..]}
    Empty dict when the run has no land-use-based loads."""
    import pandas as pd
    d = os.path.join(run_dir, "openwq_in", "ss_copernicus_files")
    p = os.path.join(d, "nutrient_loads_detailed.csv")
    if not os.path.isfile(p):
        return {}
    df = pd.read_csv(p)
    need = {"Year", "LC_Class", "Area_ha", "Nutrient", "Coefficient_kg_ha_yr", "Load_kg_yr"}
    if not need.issubset(df.columns):
        return {}
    # the unit id column is the shapefile mapping key (GRU_ID, SubId, hruId, ...):
    # the generator writes it FIRST
    id_col = next((c for c in df.columns if c not in need), df.columns[0])
    df["unit"] = df[id_col].map(_norm_id)
    df["cls"] = df["LC_Class"].map(_norm_id)
    df["sp"] = df["Nutrient"].astype(str)
    years = sorted(int(y) for y in df["Year"].dropna().unique())
    y0 = years[-1] if years else None
    dfy = df[df["Year"] == y0] if y0 is not None else df
    loads: Dict[str, Dict[str, Dict[str, float]]] = {}
    areas: Dict[str, Dict[str, float]] = {}
    coef: Dict[str, Dict[str, float]] = {}
    for r in dfy.itertuples(index=False):
        u, c, sp = str(r.unit), str(r.cls), str(r.sp)
        loads.setdefault(u, {}).setdefault(sp, {})[c] = float(r.Load_kg_yr or 0.0)
        areas.setdefault(u, {})[c] = float(r.Area_ha or 0.0)
        if r.Coefficient_kg_ha_yr == r.Coefficient_kg_ha_yr:
            coef.setdefault(sp, {})[c] = float(r.Coefficient_kg_ha_yr)
    classes: Dict[str, Dict[str, Any]] = {}
    for u, cs in areas.items():
        for c, a in cs.items():
            classes.setdefault(c, {"name": "", "area_ha": 0.0})["area_ha"] += a
    # class names: the generator's reference / a names CSV if present
    names = _lulc_class_names(d)
    for c in classes:
        classes[c]["name"] = names.get(c, "")
    return {"loads": loads, "areas": areas, "coef": coef, "classes": classes,
            "years": years, "dir": d}


def _lulc_class_names(ss_dir: str) -> Dict[str, str]:
    """Class id -> name. Looks for a names table next to the loads; falls back
    to the ESA-CCI legend for the plain numeric codes."""
    names: Dict[str, str] = {}
    leg = os.path.join(ss_dir, "lulc_source_class_legend.json")
    if os.path.isfile(leg):
        try:
            for k, v in json.load(open(leg)).items():
                names[_norm_id(k)] = str(v).split(":", 1)[-1].strip()
        except Exception:
            pass
    for cand in ("lulc_class_names.csv", "lulc_classes.csv", "class_names.csv"):
        p = os.path.join(ss_dir, cand)
        if os.path.isfile(p):
            try:
                import pandas as pd
                t = pd.read_csv(p)
                cid = next((c for c in t.columns if c.lower() in ("class", "lc_class", "code", "id", "class_id")), None)
                cnm = next((c for c in t.columns if "name" in c.lower() or "label" in c.lower()), None)
                if cid and cnm:
                    for _, r in t.iterrows():
                        names[_norm_id(r[cid])] = str(r[cnm])
            except Exception:
                pass
    esa = {"0": "no data", "10": "cropland rainfed", "11": "cropland herbaceous", "12": "cropland tree/shrub",
           "20": "cropland irrigated", "30": "mosaic cropland (>50%)", "40": "mosaic natural veg (>50%)",
           "50": "tree broadleaf evergreen", "60": "tree broadleaf deciduous", "61": "tree broadleaf deciduous closed",
           "62": "tree broadleaf deciduous open", "70": "tree needleleaf evergreen", "71": "tree needleleaf evergreen closed",
           "72": "tree needleleaf evergreen open", "80": "tree needleleaf deciduous", "81": "tree needleleaf deciduous closed",
           "82": "tree needleleaf deciduous open", "90": "tree mixed", "100": "mosaic tree/shrub", "110": "mosaic herbaceous",
           "120": "shrubland", "121": "shrubland evergreen", "122": "shrubland deciduous", "130": "grassland",
           "140": "lichens and mosses", "150": "sparse vegetation", "151": "sparse tree", "152": "sparse shrub",
           "153": "sparse herbaceous", "160": "tree flooded fresh", "170": "tree flooded saline", "180": "shrub/herb flooded",
           "190": "urban", "200": "bare areas", "201": "bare consolidated", "202": "bare unconsolidated",
           "210": "water", "220": "snow and ice"}
    for k, v in esa.items():
        names.setdefault(k, v)
    return names


def class_group_ids(ctx: Dict[str, Any], group: str) -> List[str]:
    """Class ids of a default group ('crop'|'grass'|'forest'|'urban'|'all')
    among the classes present in the run."""
    present = list((ctx.get("classes") or {}).keys())
    if group == "all":
        return present
    fam = {"crop": _CROP_CODES, "grass": _GRASS_CODES, "forest": _FOREST_CODES, "urban": _URBAN_CODES}.get(group, set())

    def _native(c):
        try:
            v = int(float(c))
            return v % 1000000 if v >= 1000000 else v
        except (TypeError, ValueError):
            return None
    out = [c for c in present if _native(c) in fam]
    if not out:   # fall back to the names
        kw = {"crop": ("crop", "agric", "arable", "cultiv"), "grass": ("grass", "pasture", "herb"),
              "forest": ("tree", "forest", "wood"), "urban": ("urban", "built", "artificial")}.get(group, ())
        out = [c for c in present if any(k in (ctx["classes"][c].get("name") or "").lower() for k in kw)]
    return out


def _norm_id(v) -> str:
    try:
        f = float(v)
        return str(int(f)) if f.is_integer() else str(f)
    except (TypeError, ValueError):
        return str(v).strip()


def _as_list(v) -> List[str]:
    if v is None:
        return []
    if isinstance(v, (list, tuple, set)):
        return [str(x).strip() for x in v if str(x).strip()]
    s = str(v).strip()
    if not s:
        return []
    return [x.strip() for x in re.split(r"[,\s;]+", s) if x.strip()]


def _resolve_classes(val, ctx) -> Optional[List[str]]:
    """'crop'|'grass'|'forest'|'urban'|'all'|list of ids -> list of class ids
    (None = every class)."""
    if val is None or (isinstance(val, str) and val.strip().lower() in ("", "all")):
        return None
    if isinstance(val, str) and val.strip().lower() in ("crop", "grass", "forest", "urban"):
        return class_group_ids(ctx, val.strip().lower())
    return [_norm_id(x) for x in _as_list(val)]


def _resolve_species(val, model_species) -> List[str]:
    if val is None or (isinstance(val, str) and val.strip().lower() in ("", "all", "n", "p")):
        g = (val or "all").strip().lower() if isinstance(val, str) and val.strip() else "all"
        return _species_default(model_species, {"all": "all", "n": "N", "p": "P"}[g])
    want = {_canon(x) for x in _as_list(val)}
    return [s for s in (model_species or []) if _canon(s) in want] or list(_as_list(val))


def _resolve_units(val) -> Optional[set]:
    """'all' / blank / None -> None (every unit); else the set of unit ids."""
    if val is None or (isinstance(val, str) and val.strip().lower() in ("", "all")):
        return None
    ids = {_norm_id(x) for x in _as_list(val)}
    if not ids or any(str(x).lower() == "all" for x in ids):
        return None
    return ids


def _monthly(val, n=12) -> List[float]:
    """'5' -> [5]*12 ; '1,2,...,12' -> 12 values."""
    vals = [float(x) for x in _as_list(val)] if not isinstance(val, (int, float)) else [float(val)]
    if not vals:
        return [0.0] * n
    if len(vals) == 1:
        return vals * n
    if len(vals) != n:
        raise ValueError(f"expected 1 or {n} values, got {len(vals)}: {val}")
    return vals


# ---------------------------------------------------------------------------
# Factor field: {unit: {species: {month(1..12) or 0 (=all months): factor}}}
# Every ss_factor lever contributes multiplicatively into this field; it is
# applied once to the SS files at the end (so levers compose).
# ---------------------------------------------------------------------------
def _mul(field, unit, species, month, f):
    field.setdefault(unit, {}).setdefault(species, {})
    field[unit][species][month] = field[unit][species].get(month, 1.0) * float(f)


def _lever_factor_field(lever: Dict, params: Dict, ctx: Dict, model_species: List[str],
                        field: Dict, notes: List[str]) -> None:
    lid = lever["id"]
    species = _resolve_species(params.get("species"), model_species)
    loads, areas, coef = ctx.get("loads", {}), ctx.get("areas", {}), ctx.get("coef", {})
    units_all = list(loads.keys())

    def _per(u, sp):
        """Per-class loads of species sp in unit u (the breakdown's nutrient
        that corresponds to the model species)."""
        d = loads.get(u) or {}
        k = _load_key(sp, d)
        return (d.get(k) if k is not None else None) or {}

    def _share(u, sp, classes):
        """Share of unit u's species-sp load that comes from `classes`."""
        per = _per(u, sp)
        tot = sum(per.values())
        if tot <= 0:
            return 0.0
        return sum(v for c, v in per.items() if classes is None or c in classes) / tot

    units = _resolve_units(params.get("units"))
    if units is not None and loads:
        _missing = sorted(units - set(units_all))
        if _missing and lever.get("mech") == "ss_factor":
            notes.append(f"{lid}: {len(_missing)} of the selected unit(s) have no land-use-based "
                         f"loads (e.g. {', '.join(_missing[:3])}) — no effect there.")
    targets = [u for u in units_all if units is None or u in units]

    if lever.get("model") == "class_reduction":
        classes = _resolve_classes(params.get("classes"), ctx)
        adopt = float(params.get("adoption_pct", 100) or 0) / 100.0
        # reduction per species: N / P / other groups (legacy form: one
        # `reduction_pct` on the `species` list)
        red_of: Dict[str, float] = {}
        if params.get("reduction_pct") not in (None, ""):
            for sp in species:
                red_of[sp] = float(params["reduction_pct"]) / 100.0
        else:
            if lever.get("sign") == 1:      # a disturbance: the export INCREASES
                rn = -float(params.get("increase_N_pct") or 0.0) / 100.0
                rp = -float(params.get("increase_P_pct") or 0.0) / 100.0
                ro = -float(params.get("increase_other_pct") or 0.0) / 100.0
            else:
                rn = float(params.get("reduction_N_pct") or 0.0) / 100.0
                rp = float(params.get("reduction_P_pct") or 0.0) / 100.0
                ro = float(params.get("reduction_other_pct") or 0.0) / 100.0
            for sp in model_species:
                r = {"N": rn, "P": rp}.get(_species_group(sp), ro)
                if r:
                    red_of[sp] = r
        if not red_of:
            notes.append(f"{lid}: no model species is affected (all reductions are 0, or the "
                         "model has no species of the reduced group) — no effect.")
            return
        if not loads:
            notes.append(f"{lid}: no per-class load breakdown for this run "
                         "(ss_copernicus_files/nutrient_loads_detailed.csv) — applied as a "
                         "uniform reduction of every diffuse load"
                         + ("" if units is None else f" in the {len(units)} selected unit(s)") + ".")
            for sp, r in red_of.items():
                for u in (units if units is not None else ["*"]):
                    _mul(field, u, sp, 0, max(0.0, 1.0 - r * adopt))
            return
        if classes is not None and not classes:
            notes.append(f"{lid}: none of the requested classes exist in this basin — no effect.")
            return
        n_hit = 0
        for u in targets:
            for sp, r in red_of.items():
                s = _share(u, sp, classes)
                if s > 0:
                    _mul(field, u, sp, 0, max(0.0, 1.0 - s * r * adopt))
                    n_hit += 1
        if n_hit == 0:
            notes.append(f"{lid}: the selected classes export nothing in the selected units — no effect.")
        return

    if lid == "landuse_change":
        if not loads:
            notes.append("landuse_change: needs the per-class load breakdown — skipped.")
            return
        fr = _resolve_classes(params.get("from_class"), ctx) or []
        to = _resolve_classes(params.get("to_class"), ctx) or []
        if not fr:
            notes.append("landuse_change: 'from' class not present in this basin — no effect.")
            return
        to_c = to[0] if to else None
        frac = float(params.get("fraction_pct", 0)) / 100.0
        for u in targets:
            for sp in species:
                per = _per(u, sp)
                tot = sum(per.values())
                if tot <= 0:
                    continue
                spc = _load_key(sp, coef)
                c_to = (coef.get(spc, {}) or {}).get(to_c, None) if to_c else 0.0
                if c_to is None:
                    c_to = 0.0
                    notes.append(f"landuse_change: no export coefficient for class {to_c} "
                                 f"({sp}) in this run — converted area assumed to export 0.")
                removed = sum(per.get(c, 0.0) for c in fr) * frac
                added = sum((areas.get(u, {}).get(c, 0.0)) for c in fr) * frac * float(c_to)
                new = max(tot - removed + added, 0.0)
                _mul(field, u, sp, 0, new / tot)
        return

    if lid == "load_cap":
        red = float(params.get("reduction_pct", 0)) / 100.0
        for sp in species:
            if units is None:
                _mul(field, "*", sp, 0, 1.0 - red)
            else:
                for u in units:
                    _mul(field, u, sp, 0, 1.0 - red)
        return

    if lid == "application_rate_limit":
        classes = _resolve_classes(params.get("classes"), ctx)
        cap = float(params.get("max_kg_ha_yr", 0))
        if not loads:
            notes.append("application_rate_limit: needs the per-class load breakdown — skipped.")
            return
        for u in targets:
            for sp in species:
                per = _per(u, sp)
                tot = sum(per.values())
                if tot <= 0:
                    continue
                spc = _load_key(sp, coef)
                new = 0.0
                for c, v in per.items():
                    cf = (coef.get(spc, {}) or {}).get(c, None)
                    if (classes is None or c in classes) and cf and cf > cap:
                        new += v * cap / cf
                    else:
                        new += v
                _mul(field, u, sp, 0, new / tot)
        return

    if lid == "load_scale":
        f = float(params.get("factor", 1.0))
        classes = _resolve_classes(params.get("classes"), ctx)
        months = [int(m) for m in _as_list(params.get("months"))] or [0]
        for sp in species:
            if classes is None or not loads:
                for u in (units if units is not None else ["*"]):
                    for m in months:
                        _mul(field, u, sp, m, f)
            else:
                for u in (units if units is not None else units_all):
                    s = _share(u, sp, classes)
                    for m in months:
                        _mul(field, u, sp, m, 1.0 + s * (f - 1.0))
        return

    if lid == "climate_delta" and params.get("apply_loads", True):
        dP = _monthly(params.get("dP_pct", 0))
        dT = _monthly(params.get("dT_c", 0))
        pw = float(params.get("precip_power", 1.0))
        q10 = float(params.get("q10", 2.0))
        for m in range(1, 13):
            f = max(0.0, 1.0 + dP[m - 1] / 100.0) ** pw * (q10 ** (dT[m - 1] / 10.0))
            for sp in model_species:
                _mul(field, "*", sp, m, f)
        return


def _canon_match(name, d):
    """Key of dict d whose canonical form equals name's (or None)."""
    if not d:
        return None
    c = _canon(name)
    for k in d:
        if _canon(k) == c:
            return k
    return None


def _load_key(species, d):
    """Key of the load-breakdown dict ``d`` (keyed by the generator's nutrient
    names: NO3-N, NH4-N, PO4-P, TN, TP, ...) that corresponds to a MODEL
    species: the same name, else a name one is a prefix of (PO4-P_sol ->
    PO4-P), else the total of the species' group (TN / TP). None if nothing
    fits."""
    if not d:
        return None
    k = species if species in d else _canon_match(species, d)
    if k is not None:
        return k
    c = _canon(species)
    pref = [x for x in d if _canon(x) and (c.startswith(_canon(x)) or _canon(x).startswith(c))]
    if pref:
        return max(pref, key=lambda x: len(_canon(x)))
    tot = {"N": "TN", "P": "TP"}.get(_species_group(species))
    return _canon_match(tot, d) if tot else None


# ---------------------------------------------------------------------------
# Applying the field / shifts / additions to the run folder's SS files
# ---------------------------------------------------------------------------
def _ss_files(eval_dir: Path) -> List[Tuple[str, Path]]:
    """(label, path) of every SS file listed in the master (in order)."""
    master = eval_dir / "openWQ_master.json"
    if not master.exists():
        return []
    data, _ = _read_json_hdr(master)
    out = []
    ss = (data.get("OPENWQ_INPUT") or {}).get("SINK_SOURCE") or {}
    for k in sorted(ss, key=lambda x: int(x) if str(x).isdigit() else 0):
        e = ss[k]
        if isinstance(e, dict) and e.get("FILEPATH"):
            p = Path(e["FILEPATH"])
            out.append((str(e.get("LABEL", "")), p if p.is_absolute() else eval_dir / p))
    return out


def _is_diffuse_file(path: Path) -> bool:
    return "based_on_lulc" in path.name.lower() or "lulc" in path.name.lower()


def _row_fields(row):
    """(year, month, day, ix, value_index) of a JSON SS row (ints or 'all')."""
    def _i(v):
        try:
            return int(v)
        except (TypeError, ValueError):
            return None
    return _i(row[0]), _i(row[1]), _i(row[2]), _norm_id(row[6]), 9


def _unit_to_ix(unit_ids, mapper) -> Dict[str, str]:
    """unit id -> row ix (string). Uses the reach/HRU mapping when available;
    otherwise identity."""
    out = {}
    for u in unit_ids:
        ix = None
        if mapper is not None:
            try:
                from calibration_lib.ml_regionalization import resolve_cell_columns
                cols = resolve_cell_columns(mapper, u)
                if cols:
                    ix = str(cols[0][0])
            except Exception:
                ix = None
        out[str(u)] = ix if ix is not None else str(u)
    return out


def apply_factor_field(eval_dir: Path, field: Dict, mapper=None, diffuse_only=True,
                       year_factor: Optional[Dict[int, float]] = None,
                       year_min: Optional[int] = None) -> Dict[str, Any]:
    """Multiply SS rows in place. ``field[unit][species][month]`` with unit '*'
    = every unit and month 0 = every month. ``year_min``: only the rows of
    that year and later (an option with a start year). Returns a small summary."""
    if not field:
        return {"files": 0, "rows": 0}
    units = {u for u in field if u != "*"}
    u2ix = _unit_to_ix(units, mapper)
    ix2u = {}
    for u, ix in u2ix.items():
        ix2u.setdefault(ix, u)
    n_rows = n_files = 0
    for label, p in _ss_files(eval_dir):
        if not p.exists():
            continue
        if diffuse_only and not _is_diffuse_file(p):
            continue
        data, hdr = _read_json_hdr(p)
        changed = False
        for blk in data.values():
            if not isinstance(blk, dict):
                continue
            sp = str(blk.get("CHEMICAL_NAME", ""))
            rows = blk.get("DATA")
            if not rows:
                continue
            it = rows.items() if isinstance(rows, dict) else enumerate(rows)
            for k, row in it:
                if not isinstance(row, list) or len(row) < 10:
                    continue
                y, m, d, ix, vi = _row_fields(row)
                if year_min is not None and (y is None or y < year_min):
                    continue
                f = 1.0
                # rows carry the unit id itself (generator rows); the internal
                # index is only a fallback for hand-written rows
                for u in ("*", ix if ix in field else ix2u.get(ix, ix)):
                    per = field.get(u, {})
                    spk = _canon_match(sp, per)
                    if spk is None:
                        continue
                    fm = per[spk]
                    f *= fm.get(0, 1.0)
                    if m is not None:
                        f *= fm.get(m, 1.0)
                    else:
                        # a row spanning every month: use the mean monthly factor
                        mm = [fm.get(i, 1.0) for i in range(1, 13)]
                        f *= sum(mm) / 12.0
                if year_factor and y in year_factor:
                    f *= year_factor[y]
                if f != 1.0:
                    try:
                        row[vi] = float(row[vi]) * f
                        changed = True
                        n_rows += 1
                    except (TypeError, ValueError):
                        pass
        if changed:
            _write_json_hdr(p, data, hdr)
            n_files += 1
    return {"files": n_files, "rows": n_rows}


def apply_month_shift(eval_dir: Path, closed_months: List[int], mode: str,
                      species: List[str], class_share: Dict[str, Dict[str, float]],
                      mapper=None, units: Optional[set] = None) -> Dict[str, Any]:
    """Application-window lever on the diffuse SS files: for each unit ×
    species × year, the closed months' loads (× the share of the load that
    comes from the affected classes) are removed or moved to the first open
    month after the closed block. Daily 'continuous' rows are assumed (the
    generator's format); other row types are left untouched."""
    closed = {int(m) for m in closed_months}
    if not closed or len(closed) >= 12:
        return {"rows": 0}
    order = list(range(1, 13))

    def _next_open(m):
        i = order.index(m)
        for k in range(1, 13):
            mm = order[(i + k) % 12]
            if mm not in closed:
                return mm
        return m
    u2ix = _unit_to_ix(class_share.keys(), mapper)
    ix2u = {ix: u for u, ix in u2ix.items()}
    n = 0
    for label, p in _ss_files(eval_dir):
        if not p.exists() or not _is_diffuse_file(p):
            continue
        data, hdr = _read_json_hdr(p)
        changed = False
        for blk in data.values():
            if not isinstance(blk, dict):
                continue
            sp = str(blk.get("CHEMICAL_NAME", ""))
            if _canon(sp) not in {_canon(s) for s in species}:
                continue
            rows = blk.get("DATA")
            if not rows:
                continue
            items = list(rows.items()) if isinstance(rows, dict) else list(enumerate(rows))
            # group by (year, ix): move the closed months' mass to the next open month
            by_key: Dict[tuple, List[tuple]] = {}
            for k, row in items:
                if isinstance(row, list) and len(row) >= 10:
                    y, m, d, ix, vi = _row_fields(row)
                    if y is not None and m is not None:
                        by_key.setdefault((y, ix), []).append((k, row, m, d))
            for (y, ix), lst in by_key.items():
                u = ix if (ix in class_share or (units is not None and ix in units)) else ix2u.get(ix, ix)
                if units is not None and u not in units:
                    continue
                share = 1.0
                if class_share:
                    sh = class_share.get(u, {})
                    _k = _load_key(sp, sh) if sh else None
                    share = float(sh[_k]) if _k is not None else float(sh.get("*", 1.0) if sh else 1.0)
                if share <= 0:
                    continue
                moved: Dict[int, float] = {}
                for k, row, m, d in lst:
                    if m in closed:
                        v = float(row[9]) * share
                        row[9] = float(row[9]) - v
                        if mode != "remove":
                            moved[_next_open(m)] = moved.get(_next_open(m), 0.0) + v
                        changed = True
                        n += 1
                for mo, mass in moved.items():
                    tgt = [r for (k, r, m, d) in lst if m == mo]
                    if tgt:
                        add = mass / len(tgt)
                        for r in tgt:
                            r[9] = float(r[9]) + add
        if changed:
            _write_json_hdr(p, data, hdr)
    return {"rows": n}


def add_point_source(eval_dir: Path, unit: str, species: List[str], load_kg_day: float,
                     growth_pct_yr: float, start_year, compartment: str, years: List[int],
                     mapper=None, label: str = "scenario point source") -> Dict[str, Any]:
    """Append a constant (optionally growing) daily load at one unit as a new
    SS file + master entry."""
    # rows carry the unit id itself (the cell id openWQ resolves through its
    # hruId/reachID mapping, like the generator's own rows) — never the
    # internal index
    ix = str(unit).strip()
    y_first = int(start_year) if str(start_year).strip() else (years[0] if years else None)
    if not years:
        years = [y_first] if y_first else []
    blocks = {}
    n = 1
    for sp in species:
        rows = {}
        r = 1
        for y in years:
            if y_first and y < y_first:
                continue
            g = (1.0 + growth_pct_yr / 100.0) ** (y - (y_first or y))
            rows[str(r)] = [int(y), "all", "all", "all", "all", "all", str(ix), "all", "all",
                            float(load_kg_day) * g, "continuous", "day"]
            r += 1
        blocks[str(n)] = {"CHEMICAL_NAME": sp, "COMPARTMENT_NAME": compartment,
                          "COMMENT": f"{label}: {load_kg_day} kg/day at unit {unit}"
                                     + (f", +{growth_pct_yr}%/yr" if growth_pct_yr else ""),
                          "TYPE": "source", "UNITS": "kg", "DATA_FORMAT": "JSON", "DATA": rows}
        n += 1
    fn = eval_dir / "openwq_in" / "openWQ_SS_scenario_point_sources.json"
    if fn.exists():
        data, hdr = _read_json_hdr(fn)
        k0 = max([int(k) for k in data if str(k).isdigit()] or [0])
        for i, b in enumerate(blocks.values(), start=k0 + 1):
            data[str(i)] = b
    else:
        data, hdr = blocks, ["// Scenario point sources (added by the scenario runner)"]
    _write_json_hdr(fn, data, hdr)
    master = eval_dir / "openWQ_master.json"
    md, mh = _read_json_hdr(master)
    ss = md.setdefault("OPENWQ_INPUT", {}).setdefault("SINK_SOURCE", {})
    if not any(isinstance(e, dict) and str(e.get("FILEPATH", "")).endswith(fn.name) for e in ss.values()):
        k = max([int(x) for x in ss if str(x).isdigit()] or [0]) + 1
        ss[str(k)] = {"LABEL": "Scenario point sources", "FILEPATH": f"openwq_in/{fn.name}"}
        _write_json_hdr(master, md, mh)
    return {"file": str(fn), "species": species, "ix": ix}


def scale_point_entries(eval_dir: Path, species: List[str], factor: float,
                        growth_pct_yr: float, years: List[int],
                        units: Optional[set] = None) -> Dict[str, Any]:
    """Multiply the non-diffuse SS entries (CSV / point sources) by ``factor``
    (× growth compounding per year), in every unit or only in ``units``."""
    yf = {y: (1.0 + growth_pct_yr / 100.0) ** (i) for i, y in enumerate(years)} if growth_pct_yr else None
    n_rows = n_files = 0
    spc = {_canon(s) for s in species}
    for label, p in _ss_files(eval_dir):
        if not p.exists() or _is_diffuse_file(p) or "scenario_point_sources" in p.name:
            continue
        try:
            data, hdr = _read_json_hdr(p)
        except Exception:
            continue      # CSV-format entries are not editable here
        changed = False
        for blk in data.values():
            if not isinstance(blk, dict) or _canon(blk.get("CHEMICAL_NAME", "")) not in spc:
                continue
            rows = blk.get("DATA") or {}
            for row in (rows.values() if isinstance(rows, dict) else rows):
                if isinstance(row, list) and len(row) >= 10:
                    y = _row_fields(row)[0]
                    if units is not None and _norm_id(row[6]) not in units:
                        continue
                    f = factor * ((yf or {}).get(y, 1.0) if y is not None else 1.0)
                    try:
                        row[9] = float(row[9]) * f
                        changed = True
                        n_rows += 1
                    except (TypeError, ValueError):
                        pass
        if changed:
            _write_json_hdr(p, data, hdr)
            n_files += 1
    return {"files": n_files, "rows": n_rows}


def scale_initial_conditions(eval_dir: Path, species: List[str], compartments, factor: float) -> Dict[str, Any]:
    cfg = eval_dir / "openwq_in" / "openWQ_config.json"
    if not cfg.exists():
        return {"n": 0}
    data, hdr = _read_json_hdr(cfg)
    bgc = data.get("BIOGEOCHEMISTRY_CONFIGURATION") or {}
    comps = None if (compartments is None or str(compartments).strip().lower() in ("", "all")) \
        else {_canon(c) for c in _as_list(compartments)}
    spc = {_canon(s) for s in species}
    n = 0
    for cname, cblk in bgc.items():
        if comps is not None and _canon(cname) not in comps:
            continue
        ic = ((cblk or {}).get("INITIAL_CONDITIONS") or {}).get("DATA") or {}
        for sp, rows in ic.items():
            if _canon(sp) not in spc or not isinstance(rows, dict):
                continue
            for r in rows.values():
                if isinstance(r, list) and len(r) >= 4:
                    try:
                        r[3] = float(r[3]) * factor
                        n += 1
                    except (TypeError, ValueError):
                        pass
    if n:
        _write_json_hdr(cfg, data, hdr)
    return {"n": n}


def apply_bgc_map(eval_dir: Path, path: List[str], units, factor: float, mapper) -> Dict[str, Any]:
    """Multiply a BGC parameter in the selected units only (per-cell map); the
    other cells keep the (calibrated) scalar as DEFAULT."""
    from calibration_lib.ml_regionalization import build_spatial_param
    bgc = eval_dir / "openwq_in" / "openWQ_MODULE_NATIVE_BGC_FLEX.json"
    if not bgc.exists() or not path:
        return {"n": 0}
    data, hdr = _read_json_hdr(bgc)
    obj = data
    for k in path[:-1]:
        if not isinstance(obj, dict) or k not in obj:
            return {"n": 0, "error": f"path {path} not found"}
        obj = obj[k]
    cur = obj.get(path[-1])
    if isinstance(cur, dict):
        base = float(cur.get("DEFAULT", 0.0))
    else:
        base = float(cur)
    us = _resolve_units(units)
    if us is None:
        obj[path[-1]] = base * factor
        _write_json_hdr(bgc, data, hdr)
        return {"n": "all", "value": base * factor}
    if mapper is None:
        return {"n": 0, "error": "no unit mapping available"}
    m = build_spatial_param({u: base * factor for u in us}, mapper, default=base)
    obj[path[-1]] = m
    _write_json_hdr(bgc, data, hdr)
    return {"n": len(m.get("CELLS", [])), "default": base, "value": base * factor}


def scale_ts_param(eval_dir: Path, param_name_fragment: str, factor: float) -> Dict[str, Any]:
    """Multiply a sediment-transport parameter (any module file whose PARAMETERS
    hold a scalar with that name)."""
    n = 0
    for p in (eval_dir / "openwq_in").glob("openWQ_MODULE_*.json"):
        if "BGC" in p.name.upper():
            continue
        try:
            data, hdr = _read_json_hdr(p)
        except Exception:
            continue
        changed = False
        for sect in ("PARAMETERS", "PARAMETER_DEFAULTS"):
            d = data.get(sect)
            if isinstance(d, dict):
                for k, v in list(d.items()):
                    if param_name_fragment.upper() in str(k).upper() and isinstance(v, (int, float)):
                        d[k] = float(v) * factor
                        changed = True
                        n += 1
        if changed:
            _write_json_hdr(p, data, hdr)
    return {"n": n}


# ---------------------------------------------------------------------------
# Host forcing perturbation (climate levers)
# ---------------------------------------------------------------------------
def perturb_forcing(eval_dir: Path, model_config: Dict[str, Any], hostmodel: str,
                    dP_pct_monthly: List[float], dT_monthly: List[float],
                    intensify: Optional[Tuple[float, float]] = None,
                    scale_vars: Optional[Dict[str, float]] = None) -> Dict[str, Any]:
    """Write a perturbed copy of the host forcing under ``<eval>/forcing_scenario/``
    and return the override the runner needs (``forcing_path`` for SUMMA,
    ``input_dir`` for mizuRoute). ``intensify`` = (percentile, factor) scales
    the wettest days and rescales the rest to keep the total."""
    import numpy as np
    try:
        import netCDF4
    except ImportError as e:
        raise RuntimeError("climate levers need the netCDF4 package: pip install netCDF4") from e
    out_dir = eval_dir / "forcing_scenario"
    out_dir.mkdir(parents=True, exist_ok=True)
    pf = np.array([1.0 + x / 100.0 for x in dP_pct_monthly], dtype=float)
    dt = np.array(dT_monthly, dtype=float)

    def _months(tvar):
        units = getattr(tvar, "units", "")
        cal = getattr(tvar, "calendar", "standard")
        try:
            dates = netCDF4.num2date(tvar[:], units, cal)
            return np.array([d.month for d in dates], dtype=int)
        except Exception:
            return None

    def _perturb_var(ds, name, kind, months):
        v = ds.variables[name]
        tdim = v.dimensions.index("time") if "time" in v.dimensions else 0
        n = v.shape[tdim]
        step = max(1, min(n, 200000))
        for i0 in range(0, n, step):
            sl = [slice(None)] * v.ndim
            sl[tdim] = slice(i0, min(n, i0 + step))
            arr = np.array(v[tuple(sl)], dtype=float)
            mm = months[i0:i0 + arr.shape[tdim]] if months is not None else None
            shp = [1] * arr.ndim
            shp[tdim] = arr.shape[tdim]
            if kind == "P":
                fac = (pf[mm - 1] if mm is not None else np.full(arr.shape[tdim], pf.mean())).reshape(shp)
                arr = arr * fac
            else:
                add = (dt[mm - 1] if mm is not None else np.full(arr.shape[tdim], dt.mean())).reshape(shp)
                arr = arr + add
            v[tuple(sl)] = arr

    def _intensify(ds, name, pct, fac):
        v = ds.variables[name]
        arr = np.array(v[:], dtype=float)
        pos = arr > 0
        if not pos.any():
            return
        thr = np.percentile(arr[pos], pct)
        big = arr >= thr
        tot = arr.sum()
        arr[big] *= fac
        rest = (~big) & pos
        excess = arr.sum() - tot
        if rest.any() and arr[rest].sum() > excess:
            arr[rest] *= (arr[rest].sum() - excess) / arr[rest].sum()
        v[:] = arr

    result: Dict[str, Any] = {"dir": str(out_dir)}
    if hostmodel == "summa":
        fm = model_config.get("file_manager_path")
        if not fm or not os.path.isfile(fm):
            raise RuntimeError("SUMMA fileManager not found for the forcing perturbation")
        keys = {}
        for line in open(fm):
            parts = line.split("!")[0].strip().split(None, 1)
            if len(parts) == 2:
                keys[parts[0]] = parts[1].strip().strip("'\"")
        fpath = keys.get("forcingPath", "")
        flist = keys.get("forcingListFile") or keys.get("forcingList") or "forcingFileList.txt"
        spath = keys.get("settingsPath", "")
        # host-side paths (the fileManager holds container paths)
        host_root, cont_root = _docker_roots(model_config)
        fpath_h = fpath.replace(cont_root, host_root, 1) if cont_root and fpath.startswith(cont_root) else fpath
        spath_h = spath.replace(cont_root, host_root, 1) if cont_root and spath.startswith(cont_root) else spath
        lst = os.path.join(spath_h, flist)
        files = [l.strip().strip("'\"") for l in open(lst) if l.strip() and not l.strip().startswith("!")] \
            if os.path.isfile(lst) else sorted(os.path.basename(x) for x in glob.glob(os.path.join(fpath_h, "*.nc")))
        pvar = str(model_config.get("forcing_precip_var") or "pptrate")
        tvar_name = str(model_config.get("forcing_temp_var") or "airtemp")
        for f in files:
            src = os.path.join(fpath_h, f)
            dst = out_dir / f
            shutil.copyfile(src, dst)
            with netCDF4.Dataset(dst, "r+") as ds:
                months = _months(ds.variables["time"]) if "time" in ds.variables else None
                if pvar in ds.variables and np.any(pf != 1.0):
                    _perturb_var(ds, pvar, "P", months)
                if tvar_name in ds.variables and np.any(dt != 0.0):
                    _perturb_var(ds, tvar_name, "T", months)
                if intensify and pvar in ds.variables:
                    _intensify(ds, pvar, float(intensify[0]), float(intensify[1]))
                # other forcing variables scaled by a constant factor
                # (shortwave radiation, specific humidity, wind speed)
                for _vn, _fac in (scale_vars or {}).items():
                    if _vn in ds.variables and _fac != 1.0:
                        _v = ds.variables[_vn]
                        _v[:] = np.array(_v[:], dtype=float) * float(_fac)
                        result.setdefault("scaled", {})[_vn] = float(_fac)
        result["forcing_path"] = str(out_dir)
        result["files"] = files
    else:
        # mizuRoute: scale the runoff input by the precipitation factor (proxy)
        ctl = model_config.get("file_manager_path")
        keys = {}
        for line in open(ctl):
            m = re.match(r"\s*<(\w+)>\s+(\S+)", line)
            if m:
                keys[m.group(1)] = m.group(2)
        host_root, cont_root = _docker_roots(model_config)
        idir = keys.get("input_dir", "")
        idir_h = re.sub("^/+", "/", idir)
        if cont_root and idir_h.startswith(cont_root):
            idir_h = idir_h.replace(cont_root, host_root, 1)
        qfile = keys.get("fname_qsim", "")
        qvar = keys.get("vname_qsim", "ROF")
        src = os.path.join(idir_h, qfile)
        dst = out_dir / qfile
        shutil.copyfile(src, dst)
        with netCDF4.Dataset(dst, "r+") as ds:
            tname = keys.get("vname_time", "time")
            months = _months(ds.variables[tname]) if tname in ds.variables else None
            if qvar in ds.variables and np.any(pf != 1.0):
                _perturb_var(ds, qvar, "P", months)
        # everything else the control expects from input_dir is linked/copied
        for other in glob.glob(os.path.join(idir_h, "*")):
            dst_o = out_dir / os.path.basename(other)
            if os.path.basename(other) != qfile and not dst_o.exists():
                if os.path.isdir(other):
                    shutil.copytree(other, dst_o)
                else:
                    shutil.copyfile(other, dst_o)   # a symlink would not resolve inside the container
        result["input_dir"] = str(out_dir)
        result["files"] = [qfile]
        if np.any(dt != 0.0):
            result["note"] = ("temperature change has no effect on a mizuRoute-only run "
                              "(runoff is precomputed); only the load response uses it")
    return result


def _docker_roots(model_config) -> Tuple[str, str]:
    """(host_root, container_root) of the docker bind mount, best effort."""
    try:
        from calibration_lib import config_integration as ci
        cc = ci.get_container_config(model_config) or {}
        h, c = cc.get("docker_host_path"), cc.get("docker_container_path")
        if h and c:
            return str(h), str(c)
    except Exception:
        pass
    try:
        import subprocess
        r = subprocess.run(["docker", "inspect", "docker_openwq", "--format",
                            "{{range .Mounts}}{{.Source}}||{{.Destination}}{{end}}"],
                           capture_output=True, text=True, timeout=5)
        if r.returncode == 0 and "||" in r.stdout:
            h, c = r.stdout.strip().split("||")[:2]
            return h, c
    except Exception:
        pass
    return str(Path.home()), "/code"


# ---------------------------------------------------------------------------
# Apply a whole scenario to a run folder
# ---------------------------------------------------------------------------
def scenario_levers(scenario: Dict[str, Any]) -> List[Dict[str, Any]]:
    """The flat list of options of a scenario, each ``{"id", "params", "where"}``.

    A scenario is written as::

        {"name": ..., "description": ...,
         "climate":    [ {"id": ..., "params": {...}}, ... ],          # whole domain
         "management": [ {"where": "all", "options": [ {...}, ... ]},   # every HRU
                         {"name": "headwaters", "where": ["12", "15"],   # a group of HRUs
                          "options": [ {...}, ... ]} ]}

    ``management`` is a list of ZONES. The zone ``where: "all"`` holds the
    options applied everywhere; a group zone holds the options tailored to its
    HRUs. In a group, an option REPLACES the same option of the ``all`` zone
    (the ``all`` instance skips the group's HRUs — ``exclude_units``);
    different options add up. Repeatable options (``multi``) always add up.

    The earlier format (``"levers": [...]`` with an optional scenario-level
    ``"units"`` and per-option ``"units"``) is still accepted."""
    out: List[Dict[str, Any]] = []
    if "climate" not in scenario and "management" not in scenario:
        scn_units = scenario.get("units")
        for lv in scenario.get("levers") or []:
            lever = LEVER_BY_ID.get(lv.get("id")) or {}
            p = dict(lv.get("params") or {})
            if lever and lever_spatial(lever) == "units" and p.get("units") in (None, "", []):
                p["units"] = scn_units if scn_units not in (None, "", []) else "all"
            u = _resolve_units(p.get("units")) if lever and lever_spatial(lever) == "units" else None
            out.append({"id": lv.get("id"), "params": p,
                        "where": ("whole domain" if not lever or not lever_spatial(lever)
                                  else ("all units" if u is None else f"{len(u)} unit(s)"))})
        return out
    for lv in scenario.get("climate") or []:
        out.append({"id": lv.get("id"), "params": dict(lv.get("params") or {}), "where": "whole domain"})
    zones = [z for z in (scenario.get("management") or []) if isinstance(z, dict)]
    groups = [(z, _resolve_units(z.get("where"))) for z in zones]
    for z, zu in groups:
        zname = str(z.get("name") or "").strip()
        for lv in z.get("options") or []:
            lid = lv.get("id")
            lever = LEVER_BY_ID.get(lid) or {}
            p = dict(lv.get("params") or {})
            sp = lever_spatial(lever) if lever else ""
            if zu is None:
                where = "all units"
                if sp == "units":
                    p["units"] = "all"
                    if not lever.get("multi"):      # a group's own setting of this option replaces it there
                        excl = set()
                        for z2, zu2 in groups:
                            if zu2 is not None and any(o.get("id") == lid for o in z2.get("options") or []):
                                excl |= zu2
                        if excl:
                            p["exclude_units"] = sorted(excl)
                            where = f"all units except {len(excl)} with their own setting"
                elif not sp:
                    where = "whole domain"
            else:
                where = (f"group '{zname}'" if zname else "group") + f" ({len(zu)} unit(s))"
                if sp == "units":
                    p["units"] = sorted(zu)
                elif not sp:
                    where = "whole domain (set in " + (f"group '{zname}'" if zname else "a group") + ")"
            if sp == "unit":
                where = f"unit {p.get('unit') or '?'}"
            out.append({"id": lid, "params": p, "where": where})
    return out


def scenario_summary_counts(scenario: Dict[str, Any]) -> Tuple[int, int]:
    """(number of climate options, number of management & policy options)."""
    lv = scenario_levers(scenario)
    nc = sum(1 for x in lv if lever_group(x.get("id")) == "climate")
    return nc, len(lv) - nc


def apply_scenario(eval_dir, scenario: Dict[str, Any], model_config: Dict[str, Any],
                   ctx: Dict[str, Any], model_species: List[str], years: List[int],
                   mapper=None, param_info: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    """Apply every lever of ``scenario`` to the run folder. Returns a log dict
    (written to ``<eval>/scenario_applied.json`` by the runner) including the
    forcing override for the model runner."""
    eval_dir = Path(eval_dir)
    hostmodel = str(model_config.get("hostmodel") or "mizuroute").lower()
    field: Dict = {}
    notes: List[str] = []
    log: Dict[str, Any] = {"name": scenario.get("name"), "levers": [], "notes": notes}
    dated: Dict[int, Dict] = {}     # start year -> factor field of the options that start then
    # every unit id of the run (to turn "all units except ..." into a list)
    all_units = set((ctx.get("loads") or {}).keys())
    if not all_units and mapper is not None:
        try:
            all_units = {_norm_id(u) for u in mapper.get_all_reach_ids()}
        except Exception:
            all_units = set()
    forcing_override: Dict[str, Any] = {}
    compartment = str(model_config.get("ss_method_copernicus_compartment_name_for_load")
                      or ("ILAYERVOLFRACWAT_SOIL" if hostmodel == "summa" else "RIVER_NETWORK_REACHES"))
    for lv in scenario_levers(scenario):
        lid = lv.get("id")
        lever = LEVER_BY_ID.get(lid)
        if lever is None:
            notes.append(f"unknown option '{lid}' skipped")
            continue
        params = {p["key"]: p.get("default") for p in lever["params"]}
        params.update({k: v for k, v in (lv.get("params") or {}).items() if v is not None and v != ""})
        if lever_spatial(lever) == "units" and params.get("units") in (None, "", []):
            params["units"] = "all"
        # "all units except the groups that set this option themselves"
        excl = _resolve_units(params.pop("exclude_units", None))
        if excl and lever_spatial(lever) == "units":
            base = _resolve_units(params.get("units"))
            if base is not None:
                params["units"] = sorted(base - excl)
            elif all_units:
                params["units"] = sorted(all_units - excl)
            else:
                notes.append(f"{lid}: the unit ids of this run are not known, so the option could not be "
                             "limited to the units outside the groups — applied to every unit.")
            if isinstance(params.get("units"), list) and not params["units"]:
                notes.append(f"{lid}: every unit has its own group setting — the 'all units' setting has no effect.")
                log["levers"].append({"id": lid, "label": lever["label"], "params": params, "where": lv.get("where")})
                continue
        entry = {"id": lid, "label": lever["label"], "params": params, "where": lv.get("where")}
        try:
            # an option with a start year writes into its own dated factor field
            _sy = params.get("start_year")
            _fld = field
            if lever["mech"] in ("ss_factor", "forcing") and str(_sy if _sy is not None else "").strip() not in ("", "0"):
                try:
                    _fld = dated.setdefault(int(float(_sy)), {})
                except (TypeError, ValueError):
                    notes.append(f"{lid}: start year '{_sy}' not understood — applied to the whole run.")
            if lever["mech"] == "ss_factor":
                _lever_factor_field(lever, params, ctx, model_species, _fld, notes)
                if lid == "conservation_tillage" and params.get("erosion_factor") not in (None, "", 1, 1.0):
                    entry["ts"] = scale_ts_param(eval_dir, "EROSION_INDEX", float(params["erosion_factor"]))
            elif lever["mech"] == "ss_shift":
                species = _resolve_species(params.get("species"), model_species)
                classes = _resolve_classes(params.get("classes"), ctx)
                share = {}
                for u, per_sp in (ctx.get("loads") or {}).items():
                    share[u] = {}
                    for sp, per in per_sp.items():
                        tot = sum(per.values())
                        share[u][sp] = (sum(v for c, v in per.items() if classes is None or c in classes) / tot) if tot > 0 else 0.0
                entry["result"] = apply_month_shift(eval_dir, [int(m) for m in _as_list(params.get("months"))],
                                                    str(params.get("mode") or "move_to_next_open_month"),
                                                    species, share, mapper,
                                                    units=_resolve_units(params.get("units")))
            elif lever["mech"] == "ss_add":
                species = _resolve_species(params.get("species"), model_species)
                entry["result"] = add_point_source(
                    eval_dir, str(params.get("unit") or "1"), species, float(params.get("load_kg_day") or 0.0),
                    float(params.get("growth_pct_yr") or 0.0), params.get("start_year"),
                    compartment, years, mapper)
            elif lever["mech"] == "ss_entry":
                species = _resolve_species(params.get("species"), model_species)
                _pu = _resolve_units(params.get("units"))
                if lid == "wastewater_standard" and params.get("reduction_pct") in (None, ""):
                    res = {"files": 0, "rows": 0}
                    for _g, _k in (("N", "reduction_N_pct"), ("P", "reduction_P_pct")):
                        _r = float(params.get(_k) or 0.0) / 100.0
                        _sp = _species_default(model_species, _g)
                        if _r and _sp:
                            _x = scale_point_entries(eval_dir, _sp, 1.0 - _r, 0.0, years, units=_pu)
                            res = {"files": max(res["files"], _x.get("files", 0)), "rows": res["rows"] + _x.get("rows", 0)}
                    entry["result"] = res
                else:
                    if lid == "wastewater_standard":      # earlier format: one reduction on a species list
                        f, g = 1.0 - float(params.get("reduction_pct", 0)) / 100.0, 0.0
                    else:
                        f, g = float(params.get("factor", 1.0)), float(params.get("growth_pct_yr") or 0.0)
                    entry["result"] = scale_point_entries(eval_dir, species, f, g, years, units=_pu)
                if entry["result"].get("files", 0) == 0:
                    notes.append(f"{lid}: this run has no point-source / CSV load entries — no effect.")
            elif lever["mech"] == "ic_factor":
                species = _resolve_species(params.get("species"), model_species)
                entry["result"] = scale_initial_conditions(eval_dir, species, params.get("compartments"),
                                                           float(params.get("factor", 1.0)))
            elif lever["mech"] == "bgc_map":
                pname = str(params.get("parameter") or "")
                pinfo = (param_info or {}).get(pname) or {}
                path = pinfo.get("path")
                if not path:
                    notes.append(f"instream_processing: unknown parameter '{pname}' — skipped.")
                else:
                    entry["result"] = apply_bgc_map(eval_dir, list(path), params.get("units"),
                                                    float(params.get("factor", 1.0)), mapper)
            elif lever["mech"] == "forcing":
                if lid == "climate_delta":
                    dP = _monthly(params.get("dP_pct", 0))
                    dT = _monthly(params.get("dT_c", 0))
                    if params.get("apply_forcing", True) in (True, "true", "True", 1, "1"):
                        _sv = {}
                        for _k, _vn in (("dSW_pct", str(model_config.get("forcing_swrad_var") or "SWRadAtm")),
                                        ("dHum_pct", str(model_config.get("forcing_humidity_var") or "spechum")),
                                        ("dWind_pct", str(model_config.get("forcing_wind_var") or "windspd"))):
                            try:
                                _pc = float(params.get(_k) or 0.0)
                            except (TypeError, ValueError):
                                _pc = 0.0
                            if _pc:
                                _sv[_vn] = 1.0 + _pc / 100.0
                        forcing_override.update(perturb_forcing(eval_dir, model_config, hostmodel, dP, dT,
                                                                scale_vars=_sv))
                    if params.get("apply_loads", True) in (True, "true", "True", 1, "1"):
                        _lever_factor_field(lever, params, ctx, model_species, field, notes)
                else:   # precip_intensification
                    forcing_override.update(perturb_forcing(
                        eval_dir, model_config, hostmodel, [0.0] * 12, [0.0] * 12,
                        intensify=(float(params.get("percentile", 95)), float(params.get("factor", 1.2)))))
                entry["forcing"] = dict(forcing_override)   # incl. the perturbed file list
        except Exception as e:
            notes.append(f"{lid}: FAILED — {e}")
            entry["error"] = str(e)
            logger.warning(f"scenario '{scenario.get('name')}' lever {lid} failed: {e}")
        log["levers"].append(entry)
    log["ss_factor"] = apply_factor_field(eval_dir, field, mapper)
    for _y in sorted(dated):
        _r = apply_factor_field(eval_dir, dated[_y], mapper, year_min=_y)
        log["ss_factor"] = {"files": max(log["ss_factor"].get("files", 0), _r.get("files", 0)),
                            "rows": log["ss_factor"].get("rows", 0) + _r.get("rows", 0)}
    log["forcing_override"] = forcing_override
    return log


# ===========================================================================
# Running scenarios (built on the calibration machinery)
# ===========================================================================
def _slug(name: str) -> str:
    s = re.sub(r"[^A-Za-z0-9._-]+", "_", str(name).strip()).strip("_")
    return s[:60] or "scenario"


def _h5_reader_path() -> str:
    _this = Path(__file__).resolve().parent
    return os.environ.get("OPENWQ_H5_SUPPORT_LIB") or str(
        _this.parent.parent / "2_Read_Outputs" / "hdf5_support_lib")


def extract_series(eval_dir, species: List[str], compartments: List[str],
                   mapping_key: str, units: str = "MG/L") -> "pd.DataFrame":
    """Simulated series of a run as a long frame
    ``datetime, unit, species, compartment, value`` (SUMMA soil layers
    ``<hru>_z<k>`` averaged per HRU). Uses the hdf5_support_lib reader like the
    calibration objective does."""
    import pandas as pd
    import sys as _sys
    hp = _h5_reader_path()
    for p in (hp, os.path.dirname(hp)):
        if p not in _sys.path:
            _sys.path.insert(0, p)
    import Read_h5_driver as h5_lib
    out_dir = Path(eval_dir) / "openwq_out"
    res = h5_lib.Read_h5_driver(
        openwq_info={"path_to_results": str(out_dir), "mapping_key": mapping_key},
        output_format="HDF5", debugmode=False, cmp=list(compartments),
        space_elem="all", chemSpec=list(species), chemUnits=units, noDataFlag=-9999)
    frames = []
    for key, extensions in (res or {}).items():
        comp = key.split("@")[0]
        sp = key.split("@")[1].split("#")[0] if "@" in key else key
        for ext_name, data_list in extensions:
            if ext_name != "main":
                continue
            for filename, df, coords in data_list:
                if df is None or df.empty:
                    continue
                long = df.copy()
                long.index.name = "datetime"
                long = long.reset_index().melt(id_vars="datetime", var_name="col", value_name="value")
                long["value"] = pd.to_numeric(long["value"], errors="coerce")
                long = long[long["value"].notna() & (long["value"] > -9998)]
                long["unit"] = long["col"].astype(str).str.replace(r"^(reachID|hruId)_", "", regex=True)\
                    .str.replace(r"_z\d+$", "", regex=True)
                g = long.groupby(["datetime", "unit"], as_index=False)["value"].mean()
                g["species"] = sp
                g["compartment"] = comp
                frames.append(g)
    if not frames:
        return pd.DataFrame(columns=["datetime", "unit", "species", "compartment", "value"])
    df = pd.concat(frames, ignore_index=True)
    df["datetime"] = pd.to_datetime(df["datetime"], errors="coerce")
    return df.dropna(subset=["datetime"])


class _ScenarioModelRunner:
    """Factory for a ModelRunner whose per-run SUMMA fileManager / mizuRoute
    control can point at a scenario's perturbed forcing (``forcing_scenario/``
    inside the run folder). Everything else is the calibration ModelRunner."""

    @staticmethod
    def build(**kwargs):
        from calibration_lib.model_runner import ModelRunner

        class ScenarioModelRunner(ModelRunner):
            forcing_overrides: Dict[str, Dict[str, Any]] = {}

            def _summa_eval_filemanager(self, eval_dir, container_eval_dir):
                out = super()._summa_eval_filemanager(eval_dir, container_eval_dir)
                ov = self.forcing_overrides.get(str(Path(eval_dir).resolve()))
                if out and ov and ov.get("forcing_path"):
                    try:
                        p = Path(eval_dir) / "fileManager_eval.txt"
                        txt = p.read_text()
                        cpath = f"{container_eval_dir.rstrip('/')}/forcing_scenario/"
                        txt = re.sub(r"(forcingPath\s+')[^']*(')", lambda m: m.group(1) + cpath + m.group(2), txt)
                        p.write_text(txt)
                    except OSError as e:
                        logger.warning(f"scenario forcing path not applied: {e}")
                return out

            def _mizuroute_eval_control(self, eval_dir, container_eval_dir):
                ov = self.forcing_overrides.get(str(Path(eval_dir).resolve()))
                out = super()._mizuroute_eval_control(eval_dir, container_eval_dir)
                if not (ov and ov.get("input_dir")):
                    return out
                try:
                    src = Path(eval_dir) / "control_eval.txt"
                    txt = src.read_text() if out and src.exists() else Path(self.file_manager_path).read_text()
                    cpath = f"{container_eval_dir.rstrip('/')}/forcing_scenario/"
                    txt = re.sub(r"(<input_dir>\s+)[^!\n]*", lambda m: m.group(1) + cpath + "    ", txt)
                    src.write_text(txt)
                    return f"{container_eval_dir.rstrip('/')}/control_eval.txt"
                except OSError as e:
                    logger.warning(f"scenario input_dir not applied: {e}")
                    return out

        return ScenarioModelRunner(**kwargs)


def _base_values(calibration_parameters: List[Dict], work_dir: Path, param_source: str):
    """Parameter vector the scenarios start from: the calibrated best (from
    results/best_parameters.json, falling back per parameter to 'initial') or
    the initial values."""
    import numpy as np
    best = {}
    src = "initial"
    if param_source != "initial":
        bp = work_dir / "results" / "best_parameters.json"
        if bp.is_file():
            try:
                best = json.load(open(bp))
                src = "calibrated best"
            except Exception:
                best = {}
    vals = []
    n_hit = 0
    for p in calibration_parameters:
        v = best.get(p["name"])
        if v is None:
            v = p.get("initial", 0.0)
        else:
            n_hit += 1
        vals.append(float(v))
    if param_source != "initial" and calibration_parameters and n_hit == 0:
        src = "initial (no calibrated values found)"
    return np.array(vals, dtype=float), src, n_hit


def run_scenarios(*, model_config: Dict[str, Any], scenarios: List[Dict[str, Any]],
                  work_dir: Optional[str] = None, calibration_work_dir: Optional[str] = None,
                  calibration_parameters: Optional[List[Dict]] = None,
                  period: Optional[Tuple[str, str]] = None, spinup_start: Optional[str] = None,
                  container_runtime: str = "docker", docker_container_name: str = "docker_openwq",
                  docker_compose_path: str = "", executable_full_path: str = "",
                  file_manager_path: str = "", apptainer_sif_path=None, apptainer_bind_path=None,
                  n_parallel: int = 1, species: Optional[List[str]] = None,
                  compartments: Optional[List[str]] = None,
                  thresholds: Optional[Dict[str, float]] = None,
                  param_source: str = "best", include_baseline: bool = True,
                  run_mode_debug: bool = False, ml_closures=None, ml_runtime=None,
                  report_stem: str = "scenarios", clean: bool = True,
                  param_info: Optional[Dict[str, Any]] = None,
                  timeout_seconds: int = 7200) -> Dict[str, Any]:
    """Run every scenario (plus the baseline) on top of the base model config
    — normally the CALIBRATED model config (``<template>_config_run.py``,
    whose inputs already hold the calibrated values) — extract the simulated
    series and write the comparison report.
    ``calibration_parameters`` + ``results/best_parameters.json`` under
    ``work_dir`` are only used when the base is a calibration folder (legacy).
    Returns ``{"scenarios": {name: {...}}, "report": path, ...}``."""
    import numpy as np
    import pandas as pd
    from calibration_lib.parameter_handler import ParameterHandler
    from calibration_lib import config_integration as ci
    from calibration_lib.calibrated_config import read_calibrated_setup

    calibration_parameters = list(calibration_parameters or [])
    work_dir = Path(work_dir or calibration_work_dir)
    scen_root = work_dir / "scenarios"
    scen_root.mkdir(parents=True, exist_ok=True)
    hostmodel = str(model_config.get("hostmodel") or "mizuroute").lower()
    spatial = ci.get_spatial_mapping(model_config) if hasattr(ci, "get_spatial_mapping") else {}
    mapping_key = spatial.get("h5_mapping_key") or ("hruId" if hostmodel == "summa" else "reachID")
    model_species = resolve_model_species(model_config)
    species = list(species or model_species)
    compartments = list(compartments or [str(model_config.get("ss_method_copernicus_compartment_name_for_load")
                                              or ("ILAYERVOLFRACWAT_SOIL" if hostmodel == "summa"
                                                  else "RIVER_NETWORK_REACHES"))])
    sim_window = None
    if period and period[1]:
        sim_window = (spinup_start or period[0], period[1])

    base_dir = str(model_config.get("dir2save_input_files") or "")
    ctx = load_lulc_context(base_dir) if base_dir else {}
    years = ctx.get("years") or []
    if period:
        try:
            y0, y1 = int(str(period[0])[:4]), int(str(period[1])[:4])
            years = list(range(y0, y1 + 1))
        except ValueError:
            pass

    # unit mapping (reach/HRU id -> cell): the base run's output, the
    # calibration's mapping.json (via the calibrated setup provenance), or
    # this work dir's own ml_attributes
    prov = read_calibrated_setup(base_dir) if base_dir else None
    mapper = None
    try:
        from calibration_lib.reach_mapping import ReachMapper
        mapper = ReachMapper(hostmodel=hostmodel)
        cands = [Path(base_dir) / "openwq_out" / "HDF5" if base_dir else None,
                 Path(prov["mapping_json"]) if prov and prov.get("mapping_json") else None,
                 Path(prov["calibration_dir"]) / "ml_attributes" / "mapping.json" if prov and prov.get("calibration_dir") else None,
                 work_dir / "ml_attributes" / "mapping.json"]
        ok = False
        for cand in cands:
            if cand is not None and cand.exists() and mapper.load_mapping(str(cand)):
                ok = True
                break
        if not ok:
            mapper = None
    except Exception:
        mapper = None

    if calibration_parameters:
        values, src_label, n_hit = _base_values(calibration_parameters, work_dir, param_source)
        logger.info(f"Scenario parameter set: {src_label} ({n_hit}/{len(calibration_parameters)} parameters "
                    f"from the calibration results)")
    else:
        values, n_hit = np.array([], dtype=float), 0
        if prov:
            src_label = (f"calibrated model config ({prov.get('n_from_best', 0)} calibrated parameter(s)"
                         + (f", best eval {prov['best_eval']}" if prov.get("best_eval") else "") + ")")
        else:
            src_label = "model config as is (no calibrated setup found)"
        logger.info(f"Scenario base: {src_label}")

    runner = _ScenarioModelRunner.build(
        runtime=container_runtime, docker_container_name=docker_container_name,
        docker_compose_path=docker_compose_path, apptainer_sif_path=apptainer_sif_path,
        apptainer_bind_path=apptainer_bind_path, file_manager_path=file_manager_path,
        executable_full_path=executable_full_path, hostmodel=hostmodel,
        calibration_work_dir=str(scen_root), calibration_period=sim_window,
        total_evaluations=len(scenarios) + (1 if include_baseline else 0),
        timeout_seconds=timeout_seconds)
    ok_pf, msg_pf = runner.preflight()
    if not ok_pf:
        raise RuntimeError("Pre-flight check failed - " + str(msg_pf).splitlines()[0])

    todo: List[Dict[str, Any]] = []
    if include_baseline:
        todo.append({"name": "baseline", "description": "the base model, no option", "levers": []})
    todo += [s for s in scenarios if _slug(s.get("name", "")) != "baseline"]

    results: Dict[str, Any] = {"scenarios": {}, "order": [], "param_source": src_label,
                               "base_dir": base_dir, "calibrated_setup": prov,
                               "period": list(period) if period else None, "species": species,
                               "compartments": compartments, "thresholds": thresholds or {}}
    all_series = []
    for i, sc in enumerate(todo):
        name = sc.get("name") or f"scenario_{i}"
        slug = _slug(name)
        sdir = scen_root / slug
        if clean and sdir.exists():
            shutil.rmtree(sdir)
        sdir.mkdir(parents=True, exist_ok=True)
        cfg = copy.deepcopy(model_config)
        cfg["run_mode_debug"] = bool(run_mode_debug)
        handler = ParameterHandler(calibration_work_dir=str(sdir), model_config=cfg,
                                   running_on_docker=(container_runtime == "docker"),
                                   calibration_period=sim_window, ml_closures=ml_closures,
                                   ml_runtime=ml_runtime)
        entry: Dict[str, Any] = {"name": name, "slug": slug, "dir": str(sdir), "levers": scenario_levers(sc)}
        results["order"].append(name)
        results["scenarios"][name] = entry
        try:
            eval_dir = handler.setup_working_directory(0, calibration_parameters, values)
            if calibration_parameters:
                handler.apply_parameters(eval_dir, calibration_parameters, values)
            json.dump({p["name"]: float(v) for p, v in zip(calibration_parameters, values)},
                      open(Path(eval_dir) / "parameters.json", "w"), indent=2)
            log = apply_scenario(eval_dir, sc, cfg, ctx, model_species, years, mapper, param_info)
            json.dump(log, open(Path(eval_dir) / "scenario_applied.json", "w"), indent=2, default=str)
            entry["applied"] = log
            entry["eval_dir"] = str(eval_dir)
            if log.get("forcing_override"):
                runner.forcing_overrides[str(Path(eval_dir).resolve())] = log["forcing_override"]
            logger.info(f"--- Scenario '{name}' ({i + 1}/{len(todo)}) — running the model ---")
            master = str(Path(eval_dir) / "openWQ_master.json")
            ok, rt, err = runner.run_single_evaluation(Path(eval_dir), master, 900000 + i)
            entry.update({"success": bool(ok), "runtime_s": float(rt or 0), "error": err or ""})
            if not ok:
                logger.error(f"scenario '{name}' failed: {err}")
                continue
            df = extract_series(eval_dir, species, compartments, mapping_key)
            if period:
                df = df[(df["datetime"] >= pd.to_datetime(period[0])) & (df["datetime"] <= pd.to_datetime(period[1]))]
            df.insert(0, "scenario", name)
            df.to_csv(sdir / "simulated.csv", index=False)
            entry["n_rows"] = int(len(df))
            all_series.append(df)
        except Exception as e:
            entry.update({"success": False, "error": str(e)})
            logger.exception(f"scenario '{name}' failed")

    if all_series:
        big = pd.concat(all_series, ignore_index=True)
        big.to_csv(scen_root / "scenarios_simulated.csv", index=False)
        results["summary"] = summarize(big, thresholds or {}, baseline="baseline" if include_baseline else None)
        results["summary"].to_csv(scen_root / "scenarios_summary.csv", index=False)
    json.dump({k: v for k, v in results.items() if k != "summary"},
              open(scen_root / "scenarios_index.json", "w"), indent=2, default=str)
    try:
        from . import Gen_Scenario_Results_Report as G
        results["report"] = G.generate_scenario_report(
            work_dir=str(work_dir), model_config=model_config, results=results,
            scenarios=todo, report_stem=report_stem)
    except Exception as e:
        logger.exception("scenario report failed")
        results["report_error"] = str(e)
    return results


def summarize(df, thresholds: Dict[str, float], baseline: Optional[str]) -> "pd.DataFrame":
    """Per (scenario, species, unit) statistics + basin aggregate + % change
    versus the baseline scenario."""
    import pandas as pd
    import numpy as np
    rows = []
    thr = {_canon(k): float(v) for k, v in (thresholds or {}).items()}
    for (sc, sp, u), g in df.groupby(["scenario", "species", "unit"]):
        v = g["value"].to_numpy(dtype=float)
        t = thr.get(_canon(sp))
        rows.append({"scenario": sc, "species": sp, "unit": str(u), "n": len(v),
                     "mean": float(np.nanmean(v)), "median": float(np.nanmedian(v)),
                     "p90": float(np.nanpercentile(v, 90)), "max": float(np.nanmax(v)),
                     "min": float(np.nanmin(v)),
                     "exceed_frac": float(np.mean(v > t)) if t is not None else np.nan,
                     "threshold": t if t is not None else np.nan})
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    # basin aggregate = mean over units of the daily means (unweighted)
    agg = (df.groupby(["scenario", "species", "datetime"], as_index=False)["value"].mean())
    for (sc, sp), g in agg.groupby(["scenario", "species"]):
        v = g["value"].to_numpy(dtype=float)
        t = thr.get(_canon(sp))
        out.loc[len(out)] = {"scenario": sc, "species": sp, "unit": "ALL", "n": len(v),
                             "mean": float(np.nanmean(v)), "median": float(np.nanmedian(v)),
                             "p90": float(np.nanpercentile(v, 90)), "max": float(np.nanmax(v)),
                             "min": float(np.nanmin(v)),
                             "exceed_frac": float(np.mean(v > t)) if t is not None else np.nan,
                             "threshold": t if t is not None else np.nan}
    if baseline and baseline in set(out["scenario"]):
        b = out[out["scenario"] == baseline].set_index(["species", "unit"])
        for col in ("mean", "median", "p90", "max"):
            out[f"{col}_pct_change"] = [
                (100.0 * (r[col] - b.loc[(r["species"], r["unit"]), col]) / b.loc[(r["species"], r["unit"]), col])
                if (r["species"], r["unit"]) in b.index and b.loc[(r["species"], r["unit"]), col] not in (0, 0.0)
                else np.nan for _, r in out.iterrows()]
    return out


def regenerate_scenario_report(*, model_config: Dict[str, Any], scenarios: List[Dict[str, Any]],
                               work_dir: Optional[str] = None, calibration_work_dir: Optional[str] = None,
                               report_stem: str = "scenarios") -> Optional[str]:
    """Rebuild the scenario report from what is on disk (``scenarios/``):
    the index json + summary CSV written by :func:`run_scenarios`."""
    import pandas as pd
    work_dir = Path(work_dir or calibration_work_dir)
    idx = work_dir / "scenarios" / "scenarios_index.json"
    if not idx.is_file():
        logger.warning("no scenario results on disk yet (scenarios/scenarios_index.json)")
        return None
    results = json.load(open(idx))
    p = work_dir / "scenarios" / "scenarios_summary.csv"
    results["summary"] = pd.read_csv(p) if p.is_file() else pd.DataFrame()
    todo = [{"name": "baseline", "description": "the base model, no option", "levers": []}] + \
        [s for s in scenarios if _slug(s.get("name", "")) != "baseline"]
    from . import Gen_Scenario_Results_Report as G
    return G.generate_scenario_report(work_dir=str(work_dir), model_config=model_config,
                                      results=results, scenarios=todo, report_stem=report_stem)


def scenario_info_for_report(model_config: Dict[str, Any], work_dir: str,
                             param_info: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    """Everything the setup report's Scenarios tab needs to offer real choices:
    the lever catalogue, the LULC classes of this run (id, name, area), the
    spatial units, the model species, the BGC parameters usable by the
    in-stream lever, and whether point-source entries exist."""
    base = str(model_config.get("dir2save_input_files") or "")
    ctx = load_lulc_context(base) if base else {}
    classes = [{"id": c, "name": v.get("name", ""), "area_ha": round(float(v.get("area_ha", 0.0)), 1)}
               for c, v in sorted((ctx.get("classes") or {}).items(),
                                  key=lambda kv: -float(kv[1].get("area_ha", 0.0)))]
    groups = {g: class_group_ids(ctx, g) for g in ("crop", "grass", "forest", "urban")} if ctx else {}
    units: List[str] = list((ctx.get("loads") or {}).keys())
    if not units:
        from calibration_lib.calibrated_config import read_calibrated_setup
        prov = read_calibrated_setup(base) if base else None
        for mp in (Path(work_dir) / "ml_attributes" / "mapping.json",
                   Path(prov["mapping_json"]) if prov and prov.get("mapping_json") else None):
            try:
                if mp is not None and mp.is_file():
                    units = [str(k) for k in json.load(open(mp)).keys()]
                    break
            except Exception:
                continue
    bgc = {k: {"path": v.get("path")} for k, v in (param_info or {}).items()
           if isinstance(v, dict) and v.get("path") and (v.get("module") in (None, "bgc"))
           and str(v.get("path", [""])[0]).upper() == "CYCLING_FRAMEWORKS"}
    n_point = 0
    try:
        for label, p in _ss_files(Path(base)):
            if not _is_diffuse_file(p):
                n_point += 1
    except Exception:
        pass
    levers = [dict(lv, group=lever_group(lv["id"]), spatial=lever_spatial(lv)) for lv in LEVERS]
    # export coefficients (kg/ha/yr) and basin loads (kg/yr) per class and per
    # class group, for the N and P groups: the report's "estimate from
    # measurements" calculators express a physical change relative to them
    coef_cls: Dict[str, Dict[str, float]] = {}
    loads_grp: Dict[str, Dict[str, float]] = {}
    coef_grp: Dict[str, Dict[str, float]] = {}
    try:
        cf = ctx.get("coef") or {}
        def _grp_coef(c, g):
            tot = "TN" if g == "N" else "TP"
            k = _canon_match(tot, cf)
            if k is not None and c in cf[k]:
                return float(cf[k][c])
            return float(sum(v.get(c, 0.0) for sp, v in cf.items()
                             if _species_group(sp) == g and _canon(sp) not in ("TN", "TP", "TKN")))
        for c in (ctx.get("classes") or {}):
            coef_cls[c] = {"N": round(_grp_coef(c, "N"), 4), "P": round(_grp_coef(c, "P"), 4)}
        for g in ("crop", "grass", "forest", "urban", "all"):
            ids = class_group_ids(ctx, g)
            area = sum(float((ctx["classes"].get(c) or {}).get("area_ha", 0.0)) for c in ids)
            if area > 0:
                coef_grp[g] = {k: round(sum(coef_cls[c][k] * float(ctx["classes"][c]["area_ha"]) for c in ids) / area, 4)
                               for k in ("N", "P")}
                loads_grp[g] = {k: round(sum(coef_cls[c][k] * float(ctx["classes"][c]["area_ha"]) for c in ids), 1)
                                for k in ("N", "P")}
    except Exception as _e:
        logger.info(f"scenario report: class coefficients not summarised ({_e})")
    return {"levers": levers, "categories": CATEGORIES, "groups_top": GROUPS, "classes": classes[:60],
            "groups": groups, "coef": coef_cls, "coef_groups": coef_grp, "loads_groups": loads_grp,
            "units": units[:20000], "n_units": len(units), "species": resolve_model_species(model_config),
            "bgc_params": bgc, "n_point_entries": n_point,
            "hostmodel": str(model_config.get("hostmodel") or "mizuroute").lower(),
            "compartment": str(model_config.get("ss_method_copernicus_compartment_name_for_load") or ""),
            "has_lulc_breakdown": bool(ctx.get("loads"))}
