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
Gen_LULC_Sources.py — multi-source land-use/land-cover (LULC) orchestrator for
the `based_on_lulc` source/sink (SS) method.

The user selects ONE land-cover dataset in the model config (`lulc_source`).
For SS loading, LULC maps are ALTERNATIVES (not additive like observations):
you use a single map per basin, so exactly one source is used. (A list is still
accepted for back-compatibility, but only its first valid entry is used.)

Acquisition is done with Google Earth Engine (GEE): GEE hosts all of these
datasets and computes the AREA per HRU per land-cover class *server-side*
(`pixelArea().addBands(class).reduceRegions(hrus, sum().group())`), returning a
tiny table — no multi-GB local raster downloads. The output conforms to the
exact schema the existing SS pipeline consumes (Gen_SS_Driver.calc_copernicus_
lulc):
    [<mapping_key>, Year, LC_Class, Pixel_Count, Area_m2, Area_ha, Area_km2]
plus a `Source` column.

Each source keeps its NATIVE class codes, mapped to a canonical land-use
category whose export coefficients (kg/ha/yr per species) drive the SS loads —
so e.g. USDA-CDL 'corn' vs 'soybeans' get different N/P loads (crop-resolved
loading, the reason these sources are worth adding). To let one merged
coefficient table flow through the unchanged downstream, class codes are
SOURCE-NAMESPACED: LC_Class = source_id * 1_000_000 + native_code.

Requires the `earthengine-api` package and an authenticated Earth Engine
account (`earthengine authenticate`). The legacy ESA CCI path (lulc_sources =
["copernicus"] with a local/CDS raster dir) is unaffected and does NOT need GEE.
"""

import os
import sys
import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd

# Class-code namespacing: LC_Class = source_id * _NS + native_code. 1e6 is far
# larger than any native legend (CDL tops out ~256), so codes never collide.
_NS = 1_000_000

# Species the coefficient tables provide (superset; downstream filters to the
# BGC template's species). kg/ha/yr.
SPECIES = ["TN", "TP", "NO3-N", "NO2-N", "NH4-N", "PO4-P"]

# ─────────────────────────────────────────────────────────────────────────────
#  Canonical land-use categories → export coefficients (kg/ha/yr).
#  Magnitudes are consistent with the ESA CCI defaults in Gen_SS_Driver
#  (get_default_copernicus_load_coefficients); see SS_lulc_loads_reference.md.
#  Every source maps its native classes onto ONE of these categories, so the
#  coefficients stay internally consistent across datasets while crop-specific
#  categories preserve the fertiliser/manure contrast that matters for loading.
# ─────────────────────────────────────────────────────────────────────────────
BASE_COEFFICIENTS: Dict[str, Dict[str, float]] = {
    # ── Agriculture (differentiated — this is the crop-resolved value) ──
    "cropland_rowcrop":   {"TN": 18.0, "TP": 2.5, "NO3-N": 13.0, "NO2-N": 0.06, "NH4-N": 0.8, "PO4-P": 0.40},
    "cropland_generic":   {"TN": 14.0, "TP": 2.0, "NO3-N": 10.0, "NO2-N": 0.05, "NH4-N": 0.6, "PO4-P": 0.30},
    "cropland_smallgrain":{"TN": 11.0, "TP": 1.4, "NO3-N":  8.0, "NO2-N": 0.04, "NH4-N": 0.5, "PO4-P": 0.22},
    "cropland_legume":    {"TN":  9.0, "TP": 1.2, "NO3-N":  6.0, "NO2-N": 0.03, "NH4-N": 0.4, "PO4-P": 0.20},
    "rice":               {"TN": 16.0, "TP": 1.8, "NO3-N":  9.0, "NO2-N": 0.05, "NH4-N": 1.2, "PO4-P": 0.35},
    "orchard_vineyard":   {"TN": 12.0, "TP": 1.5, "NO3-N":  8.0, "NO2-N": 0.04, "NH4-N": 0.5, "PO4-P": 0.25},
    "pasture_hay":        {"TN":  6.0, "TP": 0.5, "NO3-N":  3.0, "NO2-N": 0.02, "NH4-N": 0.4, "PO4-P": 0.08},
    # ── Semi-natural / natural ──
    "grassland":          {"TN":  6.0, "TP": 0.2, "NO3-N":  3.0, "NO2-N": 0.01, "NH4-N": 0.25, "PO4-P": 0.03},
    "shrubland":          {"TN":  4.0, "TP": 0.15, "NO3-N": 2.0, "NO2-N": 0.01, "NH4-N": 0.15, "PO4-P": 0.02},
    "forest":             {"TN":  3.0, "TP": 0.15, "NO3-N": 1.5, "NO2-N": 0.005, "NH4-N": 0.15, "PO4-P": 0.01},
    "wetland":            {"TN":  5.0, "TP": 0.3, "NO3-N":  2.5, "NO2-N": 0.02, "NH4-N": 0.3, "PO4-P": 0.05},
    "sparse_veg":         {"TN":  1.0, "TP": 0.05, "NO3-N": 0.5, "NO2-N": 0.005, "NH4-N": 0.05, "PO4-P": 0.005},
    "barren":             {"TN":  0.5, "TP": 0.05, "NO3-N": 0.25, "NO2-N": 0.002, "NH4-N": 0.025, "PO4-P": 0.005},
    # ── Developed ──
    "urban":              {"TN":  8.0, "TP": 1.0, "NO3-N":  5.0, "NO2-N": 0.03, "NH4-N": 0.4, "PO4-P": 0.15},
    "urban_open":         {"TN":  6.0, "TP": 0.7, "NO3-N":  3.5, "NO2-N": 0.02, "NH4-N": 0.3, "PO4-P": 0.10},
    # ── No load ──
    "water":              {"TN": 0.0, "TP": 0.0, "NO3-N": 0.0, "NO2-N": 0.0, "NH4-N": 0.0, "PO4-P": 0.0},
    "snow_ice":           {"TN": 0.0, "TP": 0.0, "NO3-N": 0.0, "NO2-N": 0.0, "NH4-N": 0.0, "PO4-P": 0.0},
}


# ─────────────────────────────────────────────────────────────────────────────
#  Per-source native class → (name, canonical category) maps.
#  Only the categories are used for coefficients; names are for the legend.
# ─────────────────────────────────────────────────────────────────────────────
# MODIS MCD12Q1, band LC_Type1 (IGBP 17-class).
_MODIS_IGBP = {
    1: ("Evergreen needleleaf forest", "forest"),
    2: ("Evergreen broadleaf forest", "forest"),
    3: ("Deciduous needleleaf forest", "forest"),
    4: ("Deciduous broadleaf forest", "forest"),
    5: ("Mixed forest", "forest"),
    6: ("Closed shrublands", "shrubland"),
    7: ("Open shrublands", "shrubland"),
    8: ("Woody savannas", "grassland"),
    9: ("Savannas", "grassland"),
    10: ("Grasslands", "grassland"),
    11: ("Permanent wetlands", "wetland"),
    12: ("Croplands", "cropland_generic"),
    13: ("Urban and built-up", "urban"),
    14: ("Cropland/natural mosaic", "cropland_generic"),
    15: ("Snow and ice", "snow_ice"),
    16: ("Barren", "barren"),
    17: ("Water bodies", "water"),
}

# ESA WorldCover v100(2020)/v200(2021), band Map.
_WORLDCOVER = {
    10: ("Tree cover", "forest"),
    20: ("Shrubland", "shrubland"),
    30: ("Grassland", "grassland"),
    40: ("Cropland", "cropland_generic"),
    50: ("Built-up", "urban"),
    60: ("Bare / sparse vegetation", "barren"),
    70: ("Snow and ice", "snow_ice"),
    80: ("Permanent water bodies", "water"),
    90: ("Herbaceous wetland", "wetland"),
    95: ("Mangroves", "wetland"),
    100: ("Moss and lichen", "sparse_veg"),
}

# Google Dynamic World V1, band label (annual mode composite).
_DYNAMIC_WORLD = {
    0: ("Water", "water"),
    1: ("Trees", "forest"),
    2: ("Grass", "grassland"),
    3: ("Flooded vegetation", "wetland"),
    4: ("Crops", "cropland_generic"),
    5: ("Shrub and scrub", "shrubland"),
    6: ("Built area", "urban"),
    7: ("Bare ground", "barren"),
    8: ("Snow and ice", "snow_ice"),
}

# ESRI/Impact Observatory 10 m Annual LC time series, band b1 (remapped codes).
_ESRI = {
    1: ("Water", "water"),
    2: ("Trees", "forest"),
    4: ("Flooded vegetation", "wetland"),
    5: ("Crops", "cropland_generic"),
    7: ("Built area", "urban"),
    8: ("Bare ground", "barren"),
    9: ("Snow and ice", "snow_ice"),
    10: ("Clouds", "barren"),
    11: ("Rangeland", "grassland"),
}

# Copernicus CGLS-LC100, band discrete_classification (23-class).
_CGLS = {
    0: ("Unknown", "barren"),
    20: ("Shrubs", "shrubland"),
    30: ("Herbaceous vegetation", "grassland"),
    40: ("Cropland", "cropland_generic"),
    50: ("Urban / built-up", "urban"),
    60: ("Bare / sparse vegetation", "barren"),
    70: ("Snow and ice", "snow_ice"),
    80: ("Permanent water bodies", "water"),
    90: ("Herbaceous wetland", "wetland"),
    100: ("Moss and lichen", "sparse_veg"),
    111: ("Closed forest, evergreen needleleaf", "forest"),
    112: ("Closed forest, evergreen broadleaf", "forest"),
    113: ("Closed forest, deciduous needleleaf", "forest"),
    114: ("Closed forest, deciduous broadleaf", "forest"),
    115: ("Closed forest, mixed", "forest"),
    116: ("Closed forest, other", "forest"),
    121: ("Open forest, evergreen needleleaf", "forest"),
    122: ("Open forest, evergreen broadleaf", "forest"),
    123: ("Open forest, deciduous needleleaf", "forest"),
    124: ("Open forest, deciduous broadleaf", "forest"),
    125: ("Open forest, mixed", "forest"),
    126: ("Open forest, other", "forest"),
    200: ("Open sea", "water"),
}

# CORINE Land Cover (EEA), band landcover (44-class; grouped).
_CORINE = {
    111: ("Continuous urban fabric", "urban"),
    112: ("Discontinuous urban fabric", "urban"),
    121: ("Industrial/commercial units", "urban"),
    122: ("Road and rail networks", "urban"),
    123: ("Port areas", "urban"),
    124: ("Airports", "urban"),
    131: ("Mineral extraction sites", "barren"),
    132: ("Dump sites", "barren"),
    133: ("Construction sites", "barren"),
    141: ("Green urban areas", "urban_open"),
    142: ("Sport and leisure facilities", "urban_open"),
    211: ("Non-irrigated arable land", "cropland_generic"),
    212: ("Permanently irrigated land", "cropland_rowcrop"),
    213: ("Rice fields", "rice"),
    221: ("Vineyards", "orchard_vineyard"),
    222: ("Fruit trees and berry plantations", "orchard_vineyard"),
    223: ("Olive groves", "orchard_vineyard"),
    231: ("Pastures", "pasture_hay"),
    241: ("Annual + permanent crops", "cropland_generic"),
    242: ("Complex cultivation patterns", "cropland_generic"),
    243: ("Agriculture with natural vegetation", "cropland_generic"),
    244: ("Agro-forestry areas", "orchard_vineyard"),
    311: ("Broad-leaved forest", "forest"),
    312: ("Coniferous forest", "forest"),
    313: ("Mixed forest", "forest"),
    321: ("Natural grasslands", "grassland"),
    322: ("Moors and heathland", "shrubland"),
    323: ("Sclerophyllous vegetation", "shrubland"),
    324: ("Transitional woodland-shrub", "shrubland"),
    331: ("Beaches, dunes, sands", "barren"),
    332: ("Bare rocks", "barren"),
    333: ("Sparsely vegetated areas", "sparse_veg"),
    334: ("Burnt areas", "barren"),
    335: ("Glaciers and perpetual snow", "snow_ice"),
    411: ("Inland marshes", "wetland"),
    412: ("Peat bogs", "wetland"),
    421: ("Salt marshes", "wetland"),
    422: ("Salines", "wetland"),
    423: ("Intertidal flats", "wetland"),
    511: ("Water courses", "water"),
    512: ("Water bodies", "water"),
    521: ("Coastal lagoons", "water"),
    522: ("Estuaries", "water"),
    523: ("Sea and ocean", "water"),
}

# USGS NLCD, band landcover (16-class).
_NLCD = {
    11: ("Open water", "water"),
    12: ("Perennial ice/snow", "snow_ice"),
    21: ("Developed, open space", "urban_open"),
    22: ("Developed, low intensity", "urban"),
    23: ("Developed, medium intensity", "urban"),
    24: ("Developed, high intensity", "urban"),
    31: ("Barren land", "barren"),
    41: ("Deciduous forest", "forest"),
    42: ("Evergreen forest", "forest"),
    43: ("Mixed forest", "forest"),
    51: ("Dwarf scrub", "shrubland"),
    52: ("Shrub/scrub", "shrubland"),
    71: ("Grassland/herbaceous", "grassland"),
    72: ("Sedge/herbaceous", "grassland"),
    73: ("Lichens", "sparse_veg"),
    74: ("Moss", "sparse_veg"),
    81: ("Pasture/hay", "pasture_hay"),
    82: ("Cultivated crops", "cropland_generic"),
    90: ("Woody wetlands", "wetland"),
    95: ("Emergent herbaceous wetlands", "wetland"),
}

# CEC NALCMS North American Land Cover, 19-class.
_NALCMS = {
    1: ("Temperate/sub-polar needleleaf forest", "forest"),
    2: ("Sub-polar taiga needleleaf forest", "forest"),
    3: ("Tropical/sub-tropical broadleaf evergreen forest", "forest"),
    4: ("Tropical/sub-tropical broadleaf deciduous forest", "forest"),
    5: ("Temperate/sub-polar broadleaf deciduous forest", "forest"),
    6: ("Mixed forest", "forest"),
    7: ("Tropical/sub-tropical shrubland", "shrubland"),
    8: ("Temperate/sub-polar shrubland", "shrubland"),
    9: ("Tropical/sub-tropical grassland", "grassland"),
    10: ("Temperate/sub-polar grassland", "grassland"),
    11: ("Sub-polar/polar shrubland-lichen-moss", "shrubland"),
    12: ("Sub-polar/polar grassland-lichen-moss", "grassland"),
    13: ("Sub-polar/polar barren-lichen-moss", "sparse_veg"),
    14: ("Wetland", "wetland"),
    15: ("Cropland", "cropland_generic"),
    16: ("Barren land", "barren"),
    17: ("Urban and built-up", "urban"),
    18: ("Water", "water"),
    19: ("Snow and ice", "snow_ice"),
}

# GLC_FCS30D fine land cover (35-class; grouped to canonical categories).
_GLC_FCS30D = {
    10: ("Rainfed cropland", "cropland_generic"),
    11: ("Herbaceous cover cropland", "cropland_generic"),
    12: ("Tree/shrub cover cropland", "orchard_vineyard"),
    20: ("Irrigated cropland", "cropland_rowcrop"),
    51: ("Open evergreen broadleaved forest", "forest"),
    52: ("Closed evergreen broadleaved forest", "forest"),
    61: ("Open deciduous broadleaved forest", "forest"),
    62: ("Closed deciduous broadleaved forest", "forest"),
    71: ("Open evergreen needleleaved forest", "forest"),
    72: ("Closed evergreen needleleaved forest", "forest"),
    81: ("Open deciduous needleleaved forest", "forest"),
    82: ("Closed deciduous needleleaved forest", "forest"),
    91: ("Open mixed-leaf forest", "forest"),
    92: ("Closed mixed-leaf forest", "forest"),
    120: ("Shrubland", "shrubland"),
    121: ("Evergreen shrubland", "shrubland"),
    122: ("Deciduous shrubland", "shrubland"),
    130: ("Grassland", "grassland"),
    140: ("Lichens and mosses", "sparse_veg"),
    150: ("Sparse vegetation", "sparse_veg"),
    152: ("Sparse shrubland", "sparse_veg"),
    153: ("Sparse herbaceous", "sparse_veg"),
    181: ("Swamp", "wetland"),
    182: ("Marsh", "wetland"),
    183: ("Flooded flat", "wetland"),
    184: ("Saline", "wetland"),
    185: ("Mangrove", "wetland"),
    186: ("Salt marsh", "wetland"),
    187: ("Tidal flat", "wetland"),
    190: ("Impervious surfaces", "urban"),
    200: ("Bare areas", "barren"),
    201: ("Consolidated bare areas", "barren"),
    202: ("Unconsolidated bare areas", "barren"),
    210: ("Water body", "water"),
    220: ("Permanent ice and snow", "snow_ice"),
}

# USDA Cropland Data Layer, band cropland (crop-resolved; major codes).
# Non-crop CDL codes (>=111) mirror NLCD.  Unlisted crop codes default to
# cropland_generic via the source's `default_category`.
_CDL = {
    1: ("Corn", "cropland_rowcrop"),
    2: ("Cotton", "cropland_rowcrop"),
    3: ("Rice", "rice"),
    4: ("Sorghum", "cropland_rowcrop"),
    5: ("Soybeans", "cropland_legume"),
    6: ("Sunflower", "cropland_rowcrop"),
    10: ("Peanuts", "cropland_legume"),
    11: ("Tobacco", "cropland_rowcrop"),
    12: ("Sweet corn", "cropland_rowcrop"),
    13: ("Pop/orn corn", "cropland_rowcrop"),
    21: ("Barley", "cropland_smallgrain"),
    22: ("Durum wheat", "cropland_smallgrain"),
    23: ("Spring wheat", "cropland_smallgrain"),
    24: ("Winter wheat", "cropland_smallgrain"),
    25: ("Other small grains", "cropland_smallgrain"),
    26: ("Winter wheat/soybeans", "cropland_generic"),
    27: ("Rye", "cropland_smallgrain"),
    28: ("Oats", "cropland_smallgrain"),
    29: ("Millet", "cropland_smallgrain"),
    36: ("Alfalfa", "pasture_hay"),
    37: ("Other hay/non-alfalfa", "pasture_hay"),
    41: ("Sugarbeets", "cropland_rowcrop"),
    42: ("Dry beans", "cropland_legume"),
    43: ("Potatoes", "cropland_rowcrop"),
    49: ("Onions", "cropland_rowcrop"),
    54: ("Tomatoes", "cropland_rowcrop"),
    61: ("Fallow/idle cropland", "grassland"),
    111: ("Open water", "water"),
    112: ("Perennial ice/snow", "snow_ice"),
    121: ("Developed, open space", "urban_open"),
    122: ("Developed, low intensity", "urban"),
    123: ("Developed, medium intensity", "urban"),
    124: ("Developed, high intensity", "urban"),
    131: ("Barren", "barren"),
    141: ("Deciduous forest", "forest"),
    142: ("Evergreen forest", "forest"),
    143: ("Mixed forest", "forest"),
    152: ("Shrubland", "shrubland"),
    176: ("Grassland/pasture", "pasture_hay"),
    190: ("Woody wetlands", "wetland"),
    195: ("Herbaceous wetlands", "wetland"),
    204: ("Pistachios", "orchard_vineyard"),
    211: ("Olives", "orchard_vineyard"),
    212: ("Oranges", "orchard_vineyard"),
    215: ("Avocados", "orchard_vineyard"),
    217: ("Pomegranates", "orchard_vineyard"),
    218: ("Nectarines", "orchard_vineyard"),
    220: ("Plums", "orchard_vineyard"),
    221: ("Strawberries", "cropland_rowcrop"),
    225: ("Winter wheat/corn", "cropland_generic"),
    229: ("Pumpkins", "cropland_rowcrop"),
    236: ("Winter wheat/sorghum", "cropland_generic"),
    242: ("Blueberries", "orchard_vineyard"),
    243: ("Cabbage", "cropland_rowcrop"),
    244: ("Cauliflower", "cropland_rowcrop"),
    246: ("Radishes", "cropland_rowcrop"),
    247: ("Turnips", "cropland_rowcrop"),
    250: ("Cranberries", "orchard_vineyard"),
}

# AAFC Annual Crop Inventory (Canada), band landcover (crop-resolved; major).
_AAFC = {
    10: ("Cloud", "barren"),
    20: ("Water", "water"),
    30: ("Exposed land / barren", "barren"),
    34: ("Urban / developed", "urban"),
    35: ("Greenhouses", "urban"),
    50: ("Shrubland", "shrubland"),
    80: ("Wetland", "wetland"),
    110: ("Grassland", "grassland"),
    120: ("Agriculture (undifferentiated)", "cropland_generic"),
    122: ("Pasture and forages", "pasture_hay"),
    130: ("Too wet to be seeded", "cropland_generic"),
    131: ("Fallow", "grassland"),
    133: ("Barley", "cropland_smallgrain"),
    134: ("Other grains", "cropland_smallgrain"),
    135: ("Millet", "cropland_smallgrain"),
    136: ("Oats", "cropland_smallgrain"),
    137: ("Rye", "cropland_smallgrain"),
    138: ("Spelt", "cropland_smallgrain"),
    139: ("Triticale", "cropland_smallgrain"),
    140: ("Wheat", "cropland_smallgrain"),
    141: ("Switchgrass", "grassland"),
    142: ("Sorghum", "cropland_rowcrop"),
    145: ("Winter wheat", "cropland_smallgrain"),
    146: ("Spring wheat", "cropland_smallgrain"),
    147: ("Corn", "cropland_rowcrop"),
    148: ("Tobacco", "cropland_rowcrop"),
    149: ("Ginseng", "cropland_rowcrop"),
    150: ("Oilseeds", "cropland_rowcrop"),
    151: ("Borage", "cropland_rowcrop"),
    153: ("Canola / rapeseed", "cropland_rowcrop"),
    154: ("Flaxseed", "cropland_rowcrop"),
    155: ("Mustard", "cropland_rowcrop"),
    156: ("Sunflower", "cropland_rowcrop"),
    157: ("Soybeans", "cropland_legume"),
    158: ("Peas", "cropland_legume"),
    160: ("Beans", "cropland_legume"),
    162: ("Lentils", "cropland_legume"),
    167: ("Sugarbeets", "cropland_rowcrop"),
    174: ("Potatoes", "cropland_rowcrop"),
    177: ("Vegetables", "cropland_rowcrop"),
    188: ("Orchards", "orchard_vineyard"),
    189: ("Berries", "orchard_vineyard"),
    190: ("Nursery", "orchard_vineyard"),
    191: ("Vineyards", "orchard_vineyard"),
    200: ("Forest (undifferentiated)", "forest"),
    210: ("Coniferous forest", "forest"),
    220: ("Broadleaf forest", "forest"),
    230: ("Mixed forest", "forest"),
}


# ─────────────────────────────────────────────────────────────────────────────
#  Source registry.
#   source_id : small unique int used for class-code namespacing
#   access    : "gee" (official EE dataset) | "gee_community" (sat-io etc.)
#               | "local" (legacy ESA-CCI via CDS; NOT GEE)
#   gee.strategy : how to get the year's single-band class image (see engine)
# ─────────────────────────────────────────────────────────────────────────────
LULC_SOURCES: Dict[str, dict] = {
    # ── Legacy local ESA CCI (kept for back-compat; handled by Gen_SS_Driver,
    #    NOT this GEE engine). Global 300 m annual 1992-2022. ──
    "copernicus": {
        "name": "ESA CCI Land Cover (Copernicus CDS) — legacy local path",
        "description": ("Global 300 m annual land cover 1992-2022 (UN-LCCS 22 "
                        "classes). Downloaded via the Copernicus Climate Data "
                        "Store and clipped locally — the original path; does NOT "
                        "use Earth Engine."),
        "coverage": "Global",
        "years": (1992, 2022),
        "res_m": 300,
        "access": "local",
        "source_id": 0,   # native ESA CCI codes are used as-is on this path
    },
    # ── Global, annual multi-year (drop-in ESA-CCI alternatives) ──
    "modis": {
        "name": "MODIS MCD12Q1 (IGBP)",
        "description": ("Global 500 m annual land cover 2001-present, IGBP "
                        "17-class."),
        "coverage": "Global", "years": (2001, 2023), "res_m": 500,
        "access": "gee", "source_id": 1,
        "gee": {"asset": "MODIS/061/MCD12Q1", "type": "collection",
                "band": "LC_Type1", "strategy": "collection_year", "scale": 500},
        "classes": _MODIS_IGBP, "default_category": "grassland",
    },
    "glc_fcs30d": {
        "name": "GLC_FCS30D — fine global (30 m, annual 1985-2022)",
        "description": ("Global 30 m fine land cover, 35 classes, annual "
                        "1985-2022 — high-res long historical record."),
        "coverage": "Global", "years": (1985, 2022), "res_m": 30,
        "access": "gee_community", "source_id": 2,
        "gee": {"asset": "projects/sat-io/open-datasets/GLC-FCS30D/annual",
                "type": "collection", "band": None, "strategy": "glcfcs_band",
                "scale": 30},
        "classes": _GLC_FCS30D, "default_category": "grassland",
    },
    "cgls_lc100": {
        "name": "Copernicus Global Land Service (CGLS-LC100)",
        "description": ("Global 100 m annual 2015-2019, 23 discrete classes "
                        "(+ fractional cover)."),
        "coverage": "Global", "years": (2015, 2019), "res_m": 100,
        "access": "gee", "source_id": 3,
        "gee": {"asset": "COPERNICUS/Landcover/100m/Proba-V-C3/Global",
                "type": "collection", "band": "discrete_classification",
                "strategy": "collection_year", "scale": 100},
        "classes": _CGLS, "default_category": "grassland",
    },
    # ── Global, high-res recent (short period) ──
    "worldcover": {
        "name": "ESA WorldCover (10 m, 2020/2021)",
        "description": ("Global 10 m land cover, 11 classes; two epochs (2020, "
                        "2021)."),
        "coverage": "Global", "years": (2020, 2021), "res_m": 10,
        "access": "gee", "source_id": 4,
        "gee": {"type": "image_per_year", "band": "Map",
                "assets": {2020: "ESA/WorldCover/v100/2020",
                           2021: "ESA/WorldCover/v200/2021"},
                "strategy": "image_per_year", "scale": 10},
        "classes": _WORLDCOVER, "default_category": "grassland",
    },
    "dynamic_world": {
        "name": "Google Dynamic World (10 m, 2015-present)",
        "description": ("Global 10 m near-real-time land cover, 9 classes, "
                        "2015-present (annual mode composite)."),
        "coverage": "Global", "years": (2015, 2024), "res_m": 10,
        "access": "gee", "source_id": 5,
        "gee": {"asset": "GOOGLE/DYNAMICWORLD/V1", "type": "collection",
                "band": "label", "strategy": "dynamic_world", "scale": 10},
        "classes": _DYNAMIC_WORLD, "default_category": "grassland",
    },
    "esri_lulc": {
        "name": "ESRI/Impact Observatory 10 m Annual LC",
        "description": ("Global 10 m annual land cover, 9 classes, 2017-present."),
        "coverage": "Global", "years": (2017, 2023), "res_m": 10,
        "access": "gee_community", "source_id": 6,
        "gee": {"asset": "projects/sat-io/open-datasets/landcover/"
                         "ESRI_Global-LULC_10m_TS",
                "type": "collection", "band": "b1", "strategy": "collection_year",
                "scale": 10},
        "classes": _ESRI, "default_category": "grassland",
    },
    # ── Regional, high-detail ──
    "nlcd": {
        "name": "NLCD — US National Land Cover Database",
        "description": ("CONUS 30 m, 16 classes, 2001-2021 (~2-3 yr cadence)."),
        "coverage": "United States (CONUS)", "years": (2001, 2021), "res_m": 30,
        "access": "gee", "source_id": 7,
        "gee": {"asset": "USGS/NLCD_RELEASES/2021_REL/NLCD", "type": "collection",
                "band": "landcover", "strategy": "nlcd_year", "scale": 30},
        "classes": _NLCD, "default_category": "grassland",
    },
    "nalcms": {
        "name": "NALCMS — North American Land Cover (CEC)",
        "description": ("Canada+US+Mexico 30 m, 19 harmonized classes; epochs "
                        "2005/2010/2015/2020."),
        "coverage": "North America", "years": (2005, 2020), "res_m": 30,
        "access": "gee_community", "source_id": 8,
        "gee": {"asset": "projects/sat-io/open-datasets/CEC_NALCMS/"
                         "landcover_30m",
                "type": "collection", "band": None, "strategy": "nearest_image",
                "scale": 30},
        "classes": _NALCMS, "default_category": "grassland",
    },
    "corine": {
        "name": "CORINE Land Cover (EEA, Europe)",
        "description": ("Europe 100 m, 44 classes; epochs 1990/2000/2006/2012/"
                        "2018."),
        "coverage": "Europe", "years": (1990, 2018), "res_m": 100,
        "access": "gee", "source_id": 9,
        "gee": {"asset": "COPERNICUS/CORINE/V20/100m", "type": "collection",
                "band": "landcover", "strategy": "corine_year", "scale": 100},
        "classes": _CORINE, "default_category": "cropland_generic",
    },
    # ── Crop-specific (best for nutrient loading) ──
    "usda_cdl": {
        "name": "USDA Cropland Data Layer (CDL)",
        "description": ("CONUS 30 m, annual 2008-present, 100+ crop-specific "
                        "classes — crop-resolved N/P loading."),
        "coverage": "United States (CONUS)", "years": (2008, 2023), "res_m": 30,
        "access": "gee", "source_id": 10,
        "gee": {"asset": "USDA/NASS/CDL", "type": "collection",
                "band": "cropland", "strategy": "collection_year", "scale": 30},
        "classes": _CDL, "default_category": "cropland_generic",
    },
    "aafc_aci": {
        "name": "AAFC Annual Crop Inventory (Canada)",
        "description": ("Canada 30 m, annual 2009-present, crop-specific "
                        "classes — crop-resolved N/P loading."),
        "coverage": "Canada", "years": (2009, 2023), "res_m": 30,
        "access": "gee", "source_id": 11,
        "gee": {"asset": "AAFC/ACI", "type": "collection", "band": "landcover",
                "strategy": "collection_year", "scale": 30},
        "classes": _AAFC, "default_category": "cropland_generic",
    },
    # ── VECTOR source (polygons, exact areas; no Earth Engine, no login) ──
    "corine_vector": {
        "name": "CORINE Land Cover (EEA, Europe) — VECTOR (exact areas)",
        "description": ("Europe, CORINE 2018 VECTOR polygons (44 classes, 25 ha "
                        "MMU) via the EEA public ArcGIS REST service. Areas per "
                        "HRU are computed by EXACT polygon overlay (no pixels) — "
                        "needs neither Earth Engine nor a login; the 2018 epoch "
                        "is applied across the simulation period."),
        "coverage": "Europe", "years": (2018, 2018), "res_m": None,
        "access": "vector", "source_id": 12,
        "vector": {
            "service": ("https://image.discomap.eea.europa.eu/arcgis/rest/"
                        "services/Corine/CLC2018_LAEA/MapServer/0"),
            "class_field": "Code_18",
            "area_crs": 3035,   # ETRS89-LAEA Europe — equal-area, for exact m²
        },
        "classes": _CORINE, "default_category": "cropland_generic",
    },
}

SOURCE_ALIASES = {
    "esa_cci": "copernicus", "esacci": "copernicus", "cci": "copernicus",
    "cds": "copernicus", "mcd12q1": "modis", "modis_lc": "modis",
    "worldcover": "worldcover", "esa_worldcover": "worldcover",
    "dynamicworld": "dynamic_world", "dw": "dynamic_world",
    "esri": "esri_lulc", "io_lulc": "esri_lulc",
    "cgls": "cgls_lc100", "lc100": "cgls_lc100",
    "glc_fcs30": "glc_fcs30d", "glcfcs": "glc_fcs30d",
    "cdl": "usda_cdl", "cropland_data_layer": "usda_cdl", "cropscape": "usda_cdl",
    "aci": "aafc_aci", "aafc": "aafc_aci",
    "clc": "corine", "nlcd_us": "nlcd",
    "corinevector": "corine_vector", "clc_vector": "corine_vector",
    "corine_shp": "corine_vector", "corine_shapefile": "corine_vector",
}

# Access types handled by the Gen_LULC_Sources engine (compute_areas_per_hru),
# as opposed to the legacy local ESA-CCI path ("local"). GEE = raster pixel
# summation; "vector" = polygon overlay (exact areas, no Earth Engine).
_ENGINE_ACCESS = ("gee", "gee_community", "vector")

_BROWSER_UA = ("Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
               "AppleWebKit/537.36")


# ─────────────────────────────────────────────────────────────────────────────
#  Selection + priority
# ─────────────────────────────────────────────────────────────────────────────
def normalize_selection(selected) -> List[str]:
    """Resolve the config value to the SINGLE LULC source to use (returned as a
    one-item list for the internal helpers).

    Land-cover maps are alternatives — one map per basin — so exactly one source
    is used. Accepts a string ('nlcd') or, for back-compatibility, a list; if a
    list with more than one entry is given, the FIRST valid one is used and the
    rest are ignored with a note. Unknown names are dropped with a warning.
    """
    if selected is None:
        return ["copernicus"]
    if isinstance(selected, str):
        selected = [selected]
    out: List[str] = []
    for s in selected:
        key = str(s).strip().lower()
        if not key:
            continue
        key = SOURCE_ALIASES.get(key, key)
        if key in LULC_SOURCES:
            if key not in out:
                out.append(key)
        else:
            print(f"  WARNING: unknown LULC source '{s}' — ignored. Valid: "
                  f"{', '.join(LULC_SOURCES)}")
    if not out:
        return ["copernicus"]
    if len(out) > 1:
        print(f"  NOTE: LULC uses a SINGLE land-cover source (maps are "
              f"alternatives, not combined) — using '{out[0]}', ignoring "
              f"{out[1:]}.")
    return out[:1]


def engine_sources(sources) -> List[str]:
    """The selection's source(s) handled by this module's engine (GEE raster OR
    vector overlay) — i.e. everything except the legacy local ESA-CCI path."""
    return [s for s in normalize_selection(sources)
            if LULC_SOURCES[s].get("access") in _ENGINE_ACCESS]


def is_engine_selection(sources) -> bool:
    """True when the selection is served by this module (GEE or vector), not the
    legacy local ESA-CCI 'copernicus' path."""
    return bool(engine_sources(sources))


def is_gee_selection(sources) -> bool:
    """True when the selection specifically needs Google Earth Engine (raster
    sources). Vector sources are engine sources but NOT GEE."""
    norm = normalize_selection(sources)
    return any(LULC_SOURCES[s].get("access") in ("gee", "gee_community")
               for s in norm)


def gee_sources(sources) -> List[str]:
    """The GEE-capable (raster) sources of the selection — excludes the legacy
    'copernicus' local path AND any vector source."""
    return [s for s in normalize_selection(sources)
            if LULC_SOURCES[s].get("access") in ("gee", "gee_community")]


# ─────────────────────────────────────────────────────────────────────────────
#  Coefficients (native codes, source-namespaced so one merged table serves the
#  unchanged downstream loads step).
# ─────────────────────────────────────────────────────────────────────────────
def namespaced_code(source_key: str, native_code: int) -> int:
    return LULC_SOURCES[source_key]["source_id"] * _NS + int(native_code)


def source_and_native(ns_code: int) -> Tuple[Optional[str], int]:
    sid, native = divmod(int(ns_code), _NS)
    for k, cfg in LULC_SOURCES.items():
        if cfg.get("source_id") == sid:
            return k, native
    return None, native


def build_load_coefficients(sources, user_overrides=None,
                            observed_classes=None
                            ) -> Dict[int, Dict[str, float]]:
    """Merged export-coefficient table keyed by SOURCE-NAMESPACED class code.

    For the selected engine source (GEE raster or vector) and each of its native
    classes, the coefficient is BASE_COEFFICIENTS[category]. `user_overrides`
    (from the config's custom per-class dict) win: keys may be given either as
    namespaced codes or, for a single-source run, as that source's native codes.

    `observed_classes` (namespaced codes actually present in the LULC areas) is
    the SAFETY NET: any observed class NOT in the source's explicit map is given
    the source's `default_category` coefficient — so no land-use area is ever
    silently dropped to zero load (important for long-tail legends like USDA-CDL
    where only the major crops are enumerated). The defaulted classes are logged
    so they can be mapped/overridden for exact loads.
    """
    coeffs: Dict[int, Dict[str, float]] = {}
    srcs = engine_sources(sources)
    for skey in srcs:
        cfg = LULC_SOURCES[skey]
        classes = cfg.get("classes", {})
        default_cat = cfg.get("default_category", "grassland")
        for native, (_name, category) in classes.items():
            base = BASE_COEFFICIENTS.get(category, BASE_COEFFICIENTS[default_cat])
            coeffs[namespaced_code(skey, native)] = dict(base)
    # ── Safety net: cover every observed-but-unmapped class ──
    if observed_classes:
        _defaulted: Dict[str, list] = {}
        for ns in observed_classes:
            try:
                ns = int(ns)
            except (TypeError, ValueError):
                continue
            if ns in coeffs:
                continue
            skey, native = source_and_native(ns)
            if skey is None or skey not in LULC_SOURCES:
                continue
            dcat = LULC_SOURCES[skey].get("default_category", "grassland")
            coeffs[ns] = dict(BASE_COEFFICIENTS[dcat])
            _defaulted.setdefault(skey, []).append(native)
        for skey, natives in _defaulted.items():
            dcat = LULC_SOURCES[skey].get("default_category", "grassland")
            print(f"  NOTE: {len(natives)} land-cover class(es) present in the "
                  f"data but not explicitly mapped for '{skey}' "
                  f"{sorted(set(natives))[:25]} → assigned the default category "
                  f"'{dcat}'. Add them to the class map or the custom-coefficient "
                  "dict for exact loads (they are NOT zeroed).")
    # Apply user overrides
    if user_overrides:
        single = srcs[0] if len(srcs) == 1 else None
        for k, v in user_overrides.items():
            try:
                code = int(k)
            except (TypeError, ValueError):
                continue
            # Heuristic: a small code on a single-source run is a native code.
            if single is not None and code < _NS:
                code = namespaced_code(single, code)
            coeffs[code] = dict(v)
    return coeffs


def class_legend(sources) -> Dict[int, str]:
    """{namespaced_code: 'source: ClassName'} for reports/reference."""
    legend: Dict[int, str] = {}
    for skey in engine_sources(sources):
        cfg = LULC_SOURCES[skey]
        short = cfg["name"].split("—")[0].strip()
        for native, (name, _cat) in cfg.get("classes", {}).items():
            legend[namespaced_code(skey, native)] = f"{short}: {name}"
    return legend


# ─────────────────────────────────────────────────────────────────────────────
#  Google Earth Engine acquisition engine
# ─────────────────────────────────────────────────────────────────────────────
def _ee_interactive() -> bool:
    """True when we may prompt the user (a real terminal, prompts not suppressed
    by the calibration/HPC non-interactive path)."""
    return bool(getattr(sys, "stdin", None) and sys.stdin.isatty()
                and os.environ.get("OPENWQ_SUPPRESS_PROMPTS") != "1")


def _ee_pip_spec() -> str:
    """The earthengine-api version spec appropriate for this Python (>= 1.0
    needs Python >= 3.10; on 3.9 the last 0.1.x line)."""
    return ("earthengine-api>=1.0.0" if sys.version_info >= (3, 10)
            else "earthengine-api>=0.1.380,<1.0.0")


def _pip_install(spec: str) -> bool:
    import subprocess
    print(f"    → pip install '{spec}' ...")
    r = subprocess.run([sys.executable, "-m", "pip", "install", spec],
                       capture_output=True, text=True)
    if r.returncode != 0:
        print((r.stdout or "")[-1200:])
        print((r.stderr or "")[-1200:])
        print("    ✗ install failed.")
    return r.returncode == 0


def _ensure_ee():
    """Import + initialize Earth Engine, DRIVING the one-time setup in the
    terminal so the user only has to select a GEE source in the config:

      1. earthengine-api missing → offer to pip-install it (right version for
         this Python);
      2. not authenticated → run the official Google OAuth flow (ee.Authenticate,
         opens a browser / prints a URL);
      3. Cloud project needed → use $EARTHENGINE_PROJECT or prompt for it.

    Returns the initialized `ee` module, or None (with clear guidance) when the
    run is non-interactive (calibration/HPC) or the user declines. Credentials
    are cached (~/.config/earthengine) so this only happens once.
    """
    # 1) import — offer to install when missing
    try:
        import ee
    except Exception:
        spec = _ee_pip_spec()
        print("\n  ⚙  Earth Engine API (earthengine-api) is needed for the GEE "
              "LULC source(s) you selected.")
        if _ee_interactive():
            try:
                ans = input(f"     Install it now?  ({spec})  [Y/n]: ").strip().lower()
            except (EOFError, KeyboardInterrupt):
                ans = "n"
            if ans in ("", "y", "yes") and _pip_install(spec):
                try:
                    import ee
                except Exception as exc:
                    print(f"     ✗ still cannot import ee ({exc}).")
                    return None
            else:
                print(f"     Skipped. Install manually:  pip install '{spec}'\n"
                      "     (Or use lulc_sources=['copernicus'] — no GEE needed.)")
                return None
        else:
            print(f"     Non-interactive run — install first:  pip install '{spec}'")
            return None

    # 2) initialize (already-authenticated fast path)
    proj = (os.environ.get("EARTHENGINE_PROJECT")
            or os.environ.get("GOOGLE_CLOUD_PROJECT") or None)

    def _try_init(p):
        try:
            ee.Initialize(project=p) if p else ee.Initialize()
            return True, None
        except Exception as e:
            return False, e

    ok, err = _try_init(proj)
    if ok:
        return ee

    # Not initialized → needs a one-time OAuth (and maybe a project).
    if not _ee_interactive():
        print(f"\n  ✗ Earth Engine is not authenticated ({err}).\n"
              "    This run is non-interactive (calibration/HPC). Authenticate "
              "ONCE beforehand on a machine with a browser:\n"
              "        earthengine authenticate\n"
              "    and set EARTHENGINE_PROJECT to your Google Cloud project id.")
        return None

    print("\n  ⚙  Earth Engine needs a one-time sign-in (opens your browser; "
          "paste the token back here if prompted)...")
    try:
        ee.Authenticate()          # official Google OAuth flow — user-driven
    except Exception as exc:
        print(f"     ✗ authentication failed ({exc}). You can also run "
              "`earthengine authenticate` in a terminal, then re-run.")
        return None

    ok, err = _try_init(proj)
    if ok:
        return ee

    # Still failing → almost always a missing/invalid Cloud project id.
    print("     Earth Engine needs a Google Cloud project (free to create at "
          "https://console.cloud.google.com/ and register for Earth Engine).")
    for _ in range(3):
        try:
            entered = input("     Enter your Earth Engine Cloud project id "
                            "(blank to abort): ").strip()
        except (EOFError, KeyboardInterrupt):
            entered = ""
        if not entered:
            break
        ok, err = _try_init(entered)
        if ok:
            print(f"     ✓ initialized (project '{entered}'). Tip: "
                  f"export EARTHENGINE_PROJECT={entered}  to skip this next time.")
            return ee
        print(f"       ✗ '{entered}' didn't work ({err}).")
    print("  ✗ Could not initialize Earth Engine — using no GEE source.")
    return None


def _hru_featurecollection(ee, basins_hrus: Dict[str, str], id_col: str):
    """Read the basin/HRU shapefile → EE FeatureCollection in EPSG:4326, keeping
    only the id property (as a string). Uses geopandas locally (no upload)."""
    import geopandas as gpd
    shp = basins_hrus["path_to_shp"]
    gdf = gpd.read_file(shp)
    if id_col not in gdf.columns:
        raise KeyError(f"mapping_key '{id_col}' not in {shp} "
                       f"(columns: {list(gdf.columns)})")
    try:
        if gdf.crs is not None and str(gdf.crs).upper() not in ("EPSG:4326",):
            gdf = gdf.to_crs("EPSG:4326")
    except Exception:
        pass
    gdf = gdf[[id_col, "geometry"]].copy()
    gdf[id_col] = gdf[id_col].astype(str)
    feats = []
    for _, row in gdf.iterrows():
        if row.geometry is None or row.geometry.is_empty:
            continue
        geom = ee.Geometry(row.geometry.__geo_interface__, proj="EPSG:4326",
                           geodesic=False)
        feats.append(ee.Feature(geom, {id_col: row[id_col]}))
    return ee.FeatureCollection(feats), gdf


def _class_image(ee, skey: str, year: int):
    """Return the single-band integer class ee.Image for `skey`+`year`, or None
    if the source has no image for that year. Handles each source's layout."""
    cfg = LULC_SOURCES[skey]
    g = cfg["gee"]
    strat = g["strategy"]
    band = g.get("band")
    try:
        if strat == "image_per_year":
            asset = g["assets"].get(int(year))
            if asset is None:
                return None
            return ee.Image(asset).select(band)
        if strat == "collection_year":
            col = ee.ImageCollection(g["asset"]).filter(
                ee.Filter.calendarRange(int(year), int(year), "year"))
            img = col.first()
            return ee.Image(img).select(band)
        if strat == "dynamic_world":
            col = (ee.ImageCollection(g["asset"])
                   .filterDate(f"{int(year)}-01-01", f"{int(year) + 1}-01-01")
                   .select(band))
            return col.reduce(ee.Reducer.mode()).rename(band).toInt()
        if strat == "nlcd_year":
            col = ee.ImageCollection(g["asset"]).filter(
                ee.Filter.eq("system:index", str(int(year))))
            img = col.first()
            return ee.Image(img).select(band)
        if strat == "corine_year":
            # CORINE epochs: 1990,2000,2006,2012,2018 — pick nearest ≤/≈ year.
            epochs = [1990, 2000, 2006, 2012, 2018]
            ey = min(epochs, key=lambda e: abs(e - int(year)))
            col = ee.ImageCollection(g["asset"]).filter(
                ee.Filter.eq("system:index", f"Y{ey}"))
            img = col.first()
            return ee.Image(img).select(band)
        if strat == "nearest_image":
            # A collection of epoch images without a year index: mosaic all.
            return ee.ImageCollection(g["asset"]).mosaic().toInt()
        if strat == "glcfcs_band":
            # GLC_FCS30D annual: an ImageCollection of continental tiles, each a
            # multiband image with one band per year (b1=1985, then annual from
            # 2000). Mosaic tiles, then select the band for this year.
            base_years = [1985, 1990, 1995] + list(range(2000, 2023))
            if int(year) < base_years[0]:
                idx = 0
            elif int(year) in base_years:
                idx = base_years.index(int(year))
            else:
                idx = min(range(len(base_years)),
                          key=lambda i: abs(base_years[i] - int(year)))
            mosaic = ee.ImageCollection(g["asset"]).mosaic()
            return mosaic.select([idx]).toInt()
    except Exception as exc:
        print(f"    WARNING: could not build {skey} image for {year} ({exc}).")
        return None
    return None


def _areas_for_hru_chunk(ee, class_img, hru_fc, scale, id_col):
    """Server-side area (m²) per HRU per class for one FeatureCollection chunk.
    Returns a list of {id_col, 'LC_Class', 'Area_m2'} dicts."""
    area_img = ee.Image.pixelArea().addBands(class_img.rename("lc_class"))
    reduced = area_img.reduceRegions(
        collection=hru_fc,
        reducer=ee.Reducer.sum().group(groupField=1, groupName="lc_class"),
        scale=scale,
    )
    rows = []
    info = reduced.getInfo()
    for feat in info.get("features", []):
        props = feat.get("properties", {})
        hid = props.get(id_col)
        for grp in props.get("groups", []):
            try:
                cls = int(grp.get("lc_class"))
                area_m2 = float(grp.get("sum", 0.0))
            except (TypeError, ValueError):
                continue
            if area_m2 <= 0:
                continue
            rows.append({id_col: hid, "LC_Class": cls, "Area_m2": area_m2})
    return rows


def _resolve_years(skey: str, requested_years: List[int]) -> Dict[int, int]:
    """Map each requested year → the nearest year the source actually provides
    (its (min,max) range; epoch sources still return a valid image via nearest
    logic in _class_image). Returns {requested_year: source_year}."""
    lo, hi = LULC_SOURCES[skey].get("years", (None, None))
    out = {}
    for y in requested_years:
        sy = y
        if lo is not None and y < lo:
            sy = lo
        elif hi is not None and y > hi:
            sy = hi
        out[y] = sy
    return out


# ─────────────────────────────────────────────────────────────────────────────
#  Vector acquisition engine (polygon overlay — exact areas, no Earth Engine)
# ─────────────────────────────────────────────────────────────────────────────
def _query_arcgis_vector(service, class_field, bbox, page=1000):
    """Fetch polygons intersecting bbox=(minx,miny,maxx,maxy) [WGS84] from an
    ArcGIS REST MapServer/FeatureServer layer, paginated. Returns a list of
    (class_code:int, geometry:geojson-dict)."""
    import urllib.parse
    minx, miny, maxx, maxy = bbox
    out, offset = [], 0
    while True:
        params = {
            "where": "1=1",
            "geometry": f"{minx},{miny},{maxx},{maxy}",
            "geometryType": "esriGeometryEnvelope",
            "inSR": "4326", "outSR": "4326",
            "spatialRel": "esriSpatialRelIntersects",
            "outFields": class_field, "returnGeometry": "true", "f": "geojson",
            "resultRecordCount": page, "resultOffset": offset,
        }
        url = service + "/query?" + urllib.parse.urlencode(params)
        try:
            req = urllib.request.Request(url, headers={"User-Agent": _BROWSER_UA})
            with urllib.request.urlopen(req, timeout=180) as r:   # noqa: S310
                j = json.loads(r.read().decode("utf-8", "replace"))
        except Exception as exc:
            print(f"    WARNING: vector query failed ({exc}).")
            break
        feats = j.get("features", []) if isinstance(j, dict) else []
        for ft in feats:
            code = (ft.get("properties") or {}).get(class_field)
            geom = ft.get("geometry")
            try:
                code = int(str(code).strip())
            except (TypeError, ValueError):
                continue
            if geom:
                out.append((code, geom))
        if len(feats) < page:
            break
        offset += page
        if offset > 200000:
            print(f"    NOTE: vector query hit the {offset}-feature cap; stopping.")
            break
    return out


def _compute_areas_vector(skey, basins_hrus, years, output_dir):
    """Exact area per HRU per class by polygon overlay (HRUs × land-cover
    polygons fetched from an ArcGIS REST service). No Earth Engine required."""
    cfg = LULC_SOURCES[skey]
    vec = cfg["vector"]
    id_col = basins_hrus.get("mapping_key", "HRU_ID")
    try:
        import geopandas as gpd
        from shapely.geometry import shape as _shape
    except Exception as exc:
        print(f"  ✗ geopandas/shapely are required for vector LULC ({exc}).")
        return None
    try:
        gdf = gpd.read_file(basins_hrus["path_to_shp"])
    except Exception as exc:
        print(f"  ✗ could not read HRUs ({exc}).")
        return None
    if id_col not in gdf.columns:
        print(f"  ✗ mapping_key '{id_col}' not in HRU shapefile "
              f"({list(gdf.columns)}).")
        return None
    gdf = gdf[[id_col, "geometry"]].copy()
    gdf[id_col] = gdf[id_col].astype(str)
    if gdf.crs is None:
        gdf = gdf.set_crs("EPSG:4326")

    minx, miny, maxx, maxy = gdf.to_crs("EPSG:4326").total_bounds
    print(f"\n  → LULC source '{skey}' ({cfg['name']}) via EEA vector service "
          f"[bbox {minx:.3f},{miny:.3f},{maxx:.3f},{maxy:.3f}]")
    feats = _query_arcgis_vector(vec["service"], vec["class_field"],
                                 (minx, miny, maxx, maxy))
    if not feats:
        print("  vector source: no polygons returned here (outside coverage, "
              "or the service is unreachable).")
        return None
    clc = gpd.GeoDataFrame({"_class": [c for c, _ in feats]},
                           geometry=[_shape(g) for _, g in feats],
                           crs="EPSG:4326")
    ea = vec.get("area_crs", 3035)
    try:
        inter = gpd.overlay(gdf.to_crs(ea), clc.to_crs(ea),
                            how="intersection", keep_geom_type=True)
    except Exception as exc:
        print(f"  ✗ HRU × land-cover overlay failed ({exc}).")
        return None
    if len(inter) == 0:
        print("  vector source: no overlap between HRUs and land-cover polygons.")
        return None
    inter["Area_m2"] = inter.geometry.area
    grp = inter.groupby([id_col, "_class"])["Area_m2"].sum().reset_index()

    req_years = sorted({int(y) for y in years}) or [max(int(y) for y in years)]
    rows = []
    for _, r in grp.iterrows():
        for yy in req_years:      # single epoch → applied to every sim year
            rows.append({id_col: str(r[id_col]), "Year": yy,
                         "LC_Class": namespaced_code(skey, int(r["_class"])),
                         "Area_m2": float(r["Area_m2"])})
    df = pd.DataFrame(rows)
    df["Area_ha"] = df["Area_m2"] / 10_000.0
    df["Area_km2"] = df["Area_m2"] / 1_000_000.0
    df["Pixel_Count"] = 0        # vector — no pixels (areas are exact)
    df["Source"] = skey
    df = df[[id_col, "Year", "LC_Class", "Pixel_Count", "Area_m2",
             "Area_ha", "Area_km2", "Source"]]
    try:
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        df.to_csv(out / "lulc_areas_all.csv", index=False)
        leg = class_legend([skey])
        with open(out / "lulc_source_class_legend.json", "w") as fh:
            json.dump({str(k): v for k, v in leg.items()}, fh, indent=2)
    except Exception:
        pass
    print(f"    ✓ '{skey}': {len(grp)} HRU×class records (EXACT polygon areas) "
          f"across {df[id_col].nunique()} HRU(s); 2018 epoch applied to "
          f"{len(req_years)} year(s).")
    return df


def compute_areas_per_hru(sources, basins_hrus: Dict[str, str],
                          years: List[int], output_dir,
                          interactive: bool = True,
                          hru_chunk: int = 40) -> Optional[pd.DataFrame]:
    """Compute area per HRU per land-cover class per year for the selected
    source, returning a DataFrame in the Gen_SS_Driver schema (+ 'Source'), or
    None if unavailable.

    Dispatches by the source's access type: a "vector" source uses an exact
    polygon overlay (no Earth Engine); a GEE source uses server-side pixel-area
    summation. `sources` normalizes to a single source (LULC maps are
    alternatives). The source is used for ALL years (nearest-year proxying
    inside its temporal range).
    """
    srcs = engine_sources(sources)
    if not srcs:
        print("  No engine LULC source selected (only the legacy 'copernicus' "
              "local path).")
        return None
    # ── Vector source: polygon overlay, no Earth Engine ──
    if LULC_SOURCES[srcs[0]].get("access") == "vector":
        return _compute_areas_vector(srcs[0], basins_hrus, years, output_dir)
    # ── GEE (raster) sources ──
    ee = _ensure_ee()
    if ee is None:
        return None

    id_col = basins_hrus.get("mapping_key", "HRU_ID")
    try:
        hru_fc, gdf = _hru_featurecollection(ee, basins_hrus, id_col)
    except Exception as exc:
        print(f"  ✗ Could not load HRUs for GEE ({exc}).")
        return None
    n_hru = len(gdf)
    hru_list = ee.FeatureCollection(hru_fc).toList(n_hru)

    req_years = sorted({int(y) for y in years}) or [max(y for y in years)]

    for skey in srcs:
        cfg = LULC_SOURCES[skey]
        scale = cfg["gee"].get("scale", 100)
        ymap = _resolve_years(skey, req_years)
        print(f"\n  → LULC source '{skey}' ({cfg['name']}) via GEE "
              f"[{scale} m]; years {req_years} → source years "
              f"{sorted(set(ymap.values()))}")
        all_rows = []
        source_years_done = {}
        try:
            for ry in req_years:
                sy = ymap[ry]
                if sy in source_years_done:
                    # reuse the already-computed source-year areas (proxy)
                    for r in source_years_done[sy]:
                        all_rows.append({**r, "Year": ry})
                    continue
                img = _class_image(ee, skey, sy)
                if img is None:
                    print(f"    · no image for {skey} {sy}; skipping year {ry}.")
                    continue
                yr_rows = []
                for i in range(0, n_hru, hru_chunk):
                    chunk = ee.FeatureCollection(
                        hru_list.slice(i, min(i + hru_chunk, n_hru)))
                    yr_rows.extend(
                        _areas_for_hru_chunk(ee, img, chunk, scale, id_col))
                source_years_done[sy] = yr_rows
                for r in yr_rows:
                    all_rows.append({**r, "Year": ry})
                print(f"    · year {ry} (map {sy}): "
                      f"{len(yr_rows)} HRU×class area records")
        except Exception as exc:
            print(f"    WARNING: GEE query failed for '{skey}' ({exc}); "
                  "trying next source.")
            continue

        if not all_rows:
            print(f"    '{skey}' returned no area in the basin; trying next "
                  "source.")
            continue

        df = pd.DataFrame(all_rows)
        # Namespace class codes + fill the schema the SS pipeline expects.
        df["LC_Class"] = df["LC_Class"].apply(
            lambda c: namespaced_code(skey, int(c)))
        df["Area_ha"] = df["Area_m2"] / 10_000.0
        df["Area_km2"] = df["Area_m2"] / 1_000_000.0
        px = float(scale) * float(scale)
        df["Pixel_Count"] = (df["Area_m2"] / px).round().astype(int)
        df["Source"] = skey
        df = df[[id_col, "Year", "LC_Class", "Pixel_Count", "Area_m2",
                 "Area_ha", "Area_km2", "Source"]]

        # Persist (same filename/dir the local pipeline uses, so caching + the
        # rest of Gen_SS_Driver find it).
        try:
            out = Path(output_dir)
            out.mkdir(parents=True, exist_ok=True)
            df.to_csv(out / "lulc_areas_all.csv", index=False)
            # A small legend for reference/report.
            leg = class_legend([skey])
            with open(out / "lulc_source_class_legend.json", "w") as fh:
                json.dump({str(k): v for k, v in leg.items()}, fh, indent=2)
        except Exception:
            pass
        print(f"    ✓ '{skey}': {len(df)} area records across {n_hru} HRU(s), "
              f"{df['Year'].nunique()} year(s). Using this source.")
        return df

    print("  ✗ None of the selected GEE LULC sources returned data for this "
          "basin.")
    return None
