"""
Hybrid physics-ML LAYER 1 helper — build the per-unit ATTRIBUTE TABLE that a
regionalized parameter needs (``theta = g(attributes)``), automatically, from
what a domain already has:

* the model's spatial units  — the basin / catchment shapefile of the model
  config (one polygon per reach or HRU, keyed by the mapping id), with a
  fallback to a per-reach *delineated* catchment layer when the configured
  shapefile is a single lumped polygon;
* the domain's attribute rasters — ``<domain>/attributes/**/*.tif`` (soil
  class, land class, elevation, slope), summarised per polygon (majority class
  + its areal fraction for categorical rasters; mean / std / min / max for
  continuous ones);
* the river-network shapefile — per-reach slope, upstream area, length,
  sinuosity when the ids match.

Output: ``attributes.csv`` (``id`` + one column per attribute) in the schema
``ml_regionalization.load_attribute_table`` reads, plus a JSON summary the
interactive setup report bakes into the Layer-1 row (attribute dropdown,
auto-filled per-class bounds, "N spatial units" check).

Only geopandas + rasterio are required (both already used by the config
support library); nothing here touches the model configuration.
"""
from __future__ import annotations

import csv
import glob
import json
import logging
import math
import os
from typing import Any, Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)

# MODIS MCD12Q1 IGBP land-cover classes (the century-basins land_classes.tif)
MODIS_IGBP = {
    1: "Evergreen needleleaf forest", 2: "Evergreen broadleaf forest",
    3: "Deciduous needleleaf forest", 4: "Deciduous broadleaf forest",
    5: "Mixed forest", 6: "Closed shrubland", 7: "Open shrubland",
    8: "Woody savanna", 9: "Savanna", 10: "Grassland", 11: "Permanent wetland",
    12: "Cropland", 13: "Urban", 14: "Cropland/natural mosaic",
    15: "Snow and ice", 16: "Barren", 17: "Water",
}
# USDA soil texture classes as coded by SoilGrids (the soil_classes.tif)
USDA_TEXTURE = {
    1: "Clay", 2: "Silty clay", 3: "Sandy clay", 4: "Clay loam",
    5: "Silty clay loam", 6: "Sandy clay loam", 7: "Loam", 8: "Silt loam",
    9: "Sandy loam", 10: "Silt", 11: "Loamy sand", 12: "Sand",
}

_RASTER_KINDS = (
    # (column, keywords in the path, categorical?, class-name table)
    ("soil_class", ("soil",), True, USDA_TEXTURE),
    ("land_class", ("land", "lulc", "landcover", "land_cover", "lc_"), True, MODIS_IGBP),
    ("elevation", ("elev", "dem"), False, None),
    ("slope", ("slope",), False, None),
)
_RIVER_ID_CANDIDATES = ("LINKNO", "SegId", "segId", "seg_id", "reachID", "reachId",
                        "SubId", "COMID", "ID", "id")
_RIVER_NUMERIC = {"Slope": "slope", "slope": "slope", "uparea": "uparea",
                  "lengthdir": "length", "Length": "length", "length": "length",
                  "sinuosity": "sinuosity", "strmOrder": "stream_order",
                  "order": "stream_order",
                  # BasinMaker / mizuRoute-ready networks (SubId-keyed)
                  "RivSlope": "slope", "RivLength": "length", "BasSlope": "basin_slope",
                  "BasAspect": "basin_aspect", "BasArea": "basin_area",
                  "BkfWidth": "bankfull_width", "BkfDepth": "bankfull_depth",
                  "Lake_Cat": "lake_cat", "LakeArea": "lake_area", "LakeDepth": "lake_depth",
                  "MeanElev": "elevation", "DrainArea": "drain_area", "Strahler": "stream_order",
                  "Q_Mean": "q_mean", "Ch_n": "channel_n", "FloodP_n": "floodplain_n",
                  "DA_Slope": "drain_slope", "Max_DEM": "elev_max", "Min_DEM": "elev_min"}


# ---------------------------------------------------------------------------
# Config helpers
# ---------------------------------------------------------------------------
def _basin_shapefile(model_config: Dict[str, Any]) -> Tuple[Optional[str], Optional[str]]:
    """(path, mapping_key) of the basin/catchment shapefile in the model config."""
    for key in ("ss_method_copernicus_basins_hrus", "basins_hrus", "basin_shapefile_cfg"):
        blk = model_config.get(key)
        if isinstance(blk, dict) and blk.get("path_to_shp"):
            return str(blk["path_to_shp"]), (blk.get("mapping_key") or None)
    p = model_config.get("basin_shapefile")
    if p:
        return str(p), (model_config.get("basin_mapping_key") or None)
    return None, None


def _domain_dirs(model_config: Dict[str, Any]) -> List[str]:
    """Candidate domain roots in which to look for ``attributes/``."""
    cands: List[str] = []
    d2s = model_config.get("dir2save_input_files")
    if d2s:
        cands.append(os.path.dirname(os.path.abspath(str(d2s))))
    shp, _ = _basin_shapefile(model_config)
    if shp:
        p = os.path.abspath(shp)
        for _ in range(3):                     # shapefiles/catchment/x.shp -> domain
            p = os.path.dirname(p)
            cands.append(p)
    rn = model_config.get("river_network_shapefile")
    if rn:
        p = os.path.abspath(str(rn))
        for _ in range(3):
            p = os.path.dirname(p)
            cands.append(p)
    out: List[str] = []
    for c in cands:
        if c and c not in out and os.path.isdir(c):
            out.append(c)
    return out


def discover_attribute_rasters(model_config: Dict[str, Any]) -> Dict[str, str]:
    """Find ``<domain>/attributes/**/*.tif`` and label them by path keywords
    -> ``{"soil_class": path, "land_class": path, "elevation": path, ...}``."""
    found: Dict[str, str] = {}
    for root in _domain_dirs(model_config):
        adir = os.path.join(root, "attributes")
        if not os.path.isdir(adir):
            continue
        tifs = sorted(glob.glob(os.path.join(adir, "**", "*.tif"), recursive=True)
                      + glob.glob(os.path.join(adir, "**", "*.tiff"), recursive=True))
        for col, kws, _cat, _names in _RASTER_KINDS:
            if col in found:
                continue
            for t in tifs:
                rel = os.path.relpath(t, adir).lower()
                if any(k in rel for k in kws):
                    found[col] = t
                    break
        if found:
            break
    return found


# ---------------------------------------------------------------------------
# Spatial units
# ---------------------------------------------------------------------------
def _model_output_ids(model_config: Dict[str, Any], hostmodel: Optional[str] = None):
    """Ids present in the baseline model output (``<dir2save>/openwq_out/HDF5``)
    via ReachMapper, or None when no baseline output exists."""
    try:
        d2s = model_config.get("dir2save_input_files")
        h5dir = os.path.join(str(d2s), "openwq_out", "HDF5") if d2s else ""
        if not (h5dir and os.path.isdir(h5dir)):
            return None
        try:
            from .reach_mapping import ReachMapper
        except ImportError:                          # pragma: no cover
            from reach_mapping import ReachMapper
        mp = ReachMapper(hostmodel=(hostmodel or model_config.get("hostmodel") or "mizuroute"))
        if not mp.load_mapping(h5dir):
            return None
        return mp
    except Exception:
        return None


def resolve_spatial_units(model_config: Dict[str, Any],
                          model_ids: Optional[Any] = None,
                          n_model_units: Optional[int] = None):
    """Return ``(gdf, id_col, notes)`` — one polygon per model spatial unit.

    Uses the model config's basin shapefile + mapping key. If that layer is a
    SINGLE lumped polygon while the MODEL has more than one unit (ids in its
    baseline output) — or the model is mizuRoute and no output exists yet — a
    sibling ``*delineated*.shp`` (per-reach catchments) carrying an id column
    is used instead: the lumped layer cannot regionalize anything. A lumped
    model (1 unit in its output) keeps its single polygon and gets a note."""
    import geopandas as gpd
    notes: List[str] = []
    n_model = n_model_units
    if n_model is None and model_ids is not None:
        try:
            from .ml_regionalization import model_unit_ids
        except ImportError:                          # pragma: no cover
            from ml_regionalization import model_unit_ids
        n_model = len(model_unit_ids(model_ids))
    hostmodel = str(model_config.get("hostmodel") or "mizuroute").lower()
    shp, key = _basin_shapefile(model_config)
    if not shp or not os.path.isfile(shp):
        raise FileNotFoundError("model config has no basin/catchment shapefile "
                                "(ss_method_copernicus_basins_hrus.path_to_shp)")
    gdf = gpd.read_file(shp)
    if not key or key not in gdf.columns:
        for c in ("GRU_ID", "HRU_ID", "hruId", "LINKNO", "SegId", "segId", "COMID", "ID"):
            if c in gdf.columns:
                notes.append(f"mapping key '{key}' not in shapefile -> using '{c}'")
                key = c
                break
    if not key or key not in gdf.columns:
        raise KeyError(f"no id column found in {os.path.basename(shp)}")
    if len(gdf) <= 1 and n_model == 1:
        notes.append("lumped model (1 spatial unit in its output): regionalization "
                     "has no effect on this domain")
    elif len(gdf) <= 1 and (n_model is None and hostmodel != "mizuroute"):
        notes.append("configured basin shapefile has 1 polygon and no baseline model "
                     "output to check the unit count against - kept as is")
    elif len(gdf) <= 1:
        sib = sorted(glob.glob(os.path.join(os.path.dirname(shp), "*delineated*.shp")))
        for s in sib:
            g2 = gpd.read_file(s)
            k2 = key if key in g2.columns else next(
                (c for c in ("GRU_ID", "LINKNO", "SegId", "segId", "COMID", "ID")
                 if c in g2.columns), None)
            if k2 and len(g2) > len(gdf):
                notes.append(f"configured basin shapefile has {len(gdf)} polygon "
                             f"-> using per-reach layer {os.path.basename(s)} "
                             f"({len(g2)} units, id '{k2}')")
                gdf, key, shp = g2, k2, s
                break
    # keep the layer's own attribute columns (BasinMaker / hydrofabric polygons
    # carry slope, elevation, drainage area, ...) — build_attribute_table reads
    # the known ones; the id column is moved first for readability
    gdf = gdf[[key] + [c for c in gdf.columns if c != key]].copy()
    gdf.attrs["source_path"] = shp
    return gdf, key, notes


def _norm_id(v):
    try:
        f = float(v)
        return int(f) if f.is_integer() else f
    except (TypeError, ValueError):
        return str(v)


# ---------------------------------------------------------------------------
# Zonal statistics (rasterio only — no rasterstats dependency)
# ---------------------------------------------------------------------------
def _zonal_arrays(src, geoms):
    """Yield the masked pixel array of ``src`` inside each geometry (falls back
    to all-touched for polygons smaller than a pixel)."""
    import numpy as np
    from rasterio.mask import mask as _mask
    nodata = src.nodata
    for geom in geoms:
        arr = None
        for all_touched in (False, True):
            try:
                out, _ = _mask(src, [geom.__geo_interface__], crop=True,
                               all_touched=all_touched, filled=True,
                               nodata=(nodata if nodata is not None else -9999))
            except Exception:
                out = None
            if out is None:
                continue
            a = out[0].astype(float)
            nd = nodata if nodata is not None else -9999
            a = a[np.isfinite(a) & (a != nd)]
            if a.size:
                arr = a
                break
        yield (arr if arr is not None else np.array([], dtype=float))


def _categorical_summary(arr, treat_zero_as_nodata=True):
    import numpy as np
    a = arr
    if treat_zero_as_nodata:
        a = a[a != 0]
    if a.size == 0:
        return None, 0.0, 0
    vals, counts = np.unique(np.round(a).astype(int), return_counts=True)
    i = int(np.argmax(counts))
    return int(vals[i]), float(counts[i]) / float(a.size), int(a.size)


def _numeric_summary(arr):
    import numpy as np
    if arr.size == 0:
        return {}
    return {"mean": float(np.mean(arr)), "std": float(np.std(arr)),
            "min": float(np.min(arr)), "max": float(np.max(arr)), "n": int(arr.size)}


# ---------------------------------------------------------------------------
# Main entry
# ---------------------------------------------------------------------------
def build_attribute_table(model_config: Dict[str, Any],
                          out_dir: str,
                          filename: str = "attributes.csv",
                          hostmodel: Optional[str] = None) -> Dict[str, Any]:
    """Build ``<out_dir>/attributes.csv`` for the model's spatial units and
    return a summary::

        {"path", "id_col", "n_units", "columns": [...],
         "categorical": {col: {"classes": {code: {"n": int, "name": str}}}},
         "numeric": {col: {"min", "max"}},
         "n_units_in_model_output": int|None, "notes": [...], "sources": {...}}

    Raises only on a missing/unreadable basin shapefile; every optional source
    (rasters, river network, model output) degrades to a note."""
    import geopandas as gpd
    import rasterio

    hostmodel = (hostmodel or model_config.get("hostmodel") or "mizuroute")
    _mp = _model_output_ids(model_config, hostmodel)
    gdf, id_col, notes = resolve_spatial_units(model_config, model_ids=_mp)
    ids = [_norm_id(v) for v in gdf[id_col].tolist()]
    rows: Dict[Any, Dict[str, Any]] = {i: {} for i in ids}
    sources: Dict[str, str] = {"units": gdf.attrs.get("source_path", "")}
    categorical: Dict[str, Dict[str, Any]] = {}
    numeric: Dict[str, Dict[str, float]] = {}

    # area (km2) in an equal-area projection
    try:
        cea = gdf.to_crs("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
        for i, a in zip(ids, cea.geometry.area.tolist()):
            rows[i]["area_km2"] = float(a) / 1e6
        numeric["area_km2"] = {}
    except Exception as e:
        notes.append(f"area not computed ({e})")

    # numeric attributes carried by the spatial-unit layer itself (BasinMaker /
    # TauDEM / hydrofabric polygons ship slope, elevation, drainage area, ...)
    _poly_cols = []
    for src_col, dst in _RIVER_NUMERIC.items():
        if src_col in gdf.columns and dst not in numeric:
            try:
                vals = [float(v) for v in gdf[src_col].tolist()]
            except (TypeError, ValueError):
                continue
            for i, v in zip(ids, vals):
                if v == v and v not in (-9999.0, -999.0):   # skip NaN / nodata
                    rows[i][dst] = v
            if any(dst in rows[i] for i in ids):
                numeric[dst] = {}
                _poly_cols.append(dst)
    if _poly_cols:
        notes.append(f"unit-layer attributes read from polygon columns: {', '.join(_poly_cols)}")

    # rasters
    rasters = discover_attribute_rasters(model_config)
    for col, kws, is_cat, names in _RASTER_KINDS:
        path = rasters.get(col)
        if not path:
            continue
        try:
            with rasterio.open(path) as src:
                g = gdf.to_crs(src.crs) if (src.crs and gdf.crs and src.crs != gdf.crs) else gdf
                for i, arr in zip(ids, _zonal_arrays(src, g.geometry.tolist())):
                    if is_cat:
                        code, frac, n = _categorical_summary(arr)
                        rows[i][col] = code
                        rows[i][f"{col}_frac"] = round(frac, 4)
                        if names is not None:
                            rows[i][f"{col}_name"] = names.get(code, "") if code is not None else ""
                    else:
                        st = _numeric_summary(arr)
                        rows[i][f"{col}_mean"] = st.get("mean")
                        rows[i][f"{col}_std"] = st.get("std")
            sources[col] = path
            if is_cat:
                cls: Dict[str, Dict[str, Any]] = {}
                for i in ids:
                    c = rows[i].get(col)
                    if c is None:
                        continue
                    k = str(c)
                    cls.setdefault(k, {"n": 0, "name": (names or {}).get(c, "")})
                    cls[k]["n"] += 1
                categorical[col] = {"classes": dict(sorted(cls.items(), key=lambda kv: int(kv[0]) if kv[0].lstrip('-').isdigit() else kv[0]))}
                numeric[f"{col}_frac"] = {}
            else:
                numeric[f"{col}_mean"] = {}
                numeric[f"{col}_std"] = {}
        except Exception as e:
            notes.append(f"{col}: raster {os.path.basename(path)} not summarised ({e})")

    # river network attributes (per reach)
    rn = model_config.get("river_network_shapefile")
    if rn and os.path.isfile(str(rn)):
        try:
            r = gpd.read_file(str(rn))
            # candidate id columns: the configured key first, then the known
            # names — keep the one whose values overlap the unit ids best (a
            # BasinMaker layer may carry both SubId and a re-coded SegId).
            _cands = [model_config.get("river_network_mapping_key")] + list(_RIVER_ID_CANDIDATES)
            _cands = [c for c in _cands if c and c in r.columns]
            rkey, rid, common = None, [], set()
            for c in dict.fromkeys(_cands):
                _rid = [_norm_id(v) for v in r[c].tolist()]
                _common = set(_rid) & set(ids)
                if len(_common) > len(common):
                    rkey, rid, common = c, _rid, _common
            if rkey is None and _cands:
                rkey = _cands[0]
            if rkey:
                if common:
                    lut = {k: idx for idx, k in enumerate(rid)}
                    for src_col, dst in _RIVER_NUMERIC.items():
                        if src_col in r.columns and dst not in numeric:
                            for i in ids:
                                j = lut.get(i)
                                if j is not None:
                                    try:
                                        rows[i][dst] = float(r[src_col].iloc[j])
                                    except (TypeError, ValueError):
                                        pass
                            numeric[dst] = {}
                    sources["river_network"] = str(rn)
                    notes.append(f"river-network attributes joined on '{rkey}' "
                                 f"({len(common)}/{len(ids)} units)")
                else:
                    notes.append(f"river-network ids ('{rkey}') do not match the unit ids")
        except Exception as e:
            notes.append(f"river network not joined ({e})")

    # do the unit ids exist in the model's baseline output? (mapping check)
    n_in_model = None
    n_model_units = None
    try:
        if _mp is not None:
            try:
                from .ml_regionalization import model_unit_ids, resolve_cell_xyz
            except ImportError:                      # pragma: no cover
                from ml_regionalization import model_unit_ids, resolve_cell_xyz
            n_model_units = len(model_unit_ids(_mp))
            n_in_model = sum(1 for i in ids if resolve_cell_xyz(_mp, i))
            if n_in_model < len(ids):
                notes.append(f"{len(ids) - n_in_model} unit id(s) not found in the "
                             f"baseline model output (model has {n_model_units} unit(s))")
    except Exception as e:
        notes.append(f"model-output mapping check skipped ({e})")

    # numeric ranges for the summary
    for col in list(numeric.keys()):
        vals = [rows[i][col] for i in ids if isinstance(rows[i].get(col), (int, float))
                and math.isfinite(rows[i][col])]
        if vals:
            numeric[col] = {"min": float(min(vals)), "max": float(max(vals))}
        else:
            numeric.pop(col, None)

    # write CSV (id + stable column order)
    os.makedirs(out_dir, exist_ok=True)
    columns: List[str] = []
    for i in ids:
        for k in rows[i]:
            if k not in columns:
                columns.append(k)
    out_path = os.path.join(out_dir, filename)
    with open(out_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["id"] + columns)
        for i in ids:
            w.writerow([i] + [("" if rows[i].get(c) is None else rows[i].get(c)) for c in columns])

    # Save the id -> (ix,iy,iz) mapping next to the table so every evaluation
    # can regionalize BEFORE it has produced any output of its own.
    mapping_file = None
    try:
        if _mp is not None:
            mapping_file = os.path.join(out_dir, "mapping.json")
            _mp.save_mapping(mapping_file)
            mapping_file = os.path.abspath(mapping_file)
    except Exception as e:
        notes.append(f"mapping not saved ({e})")
        mapping_file = None

    # Prepared Layer-1B nets in this folder (ml_regionalization.prepare_runtime_param)
    runtime_files = []
    for wf in sorted(glob.glob(os.path.join(out_dir, "_ml_runtime_*_weights.json"))):
        af = wf.replace("_weights.json", "_attributes.json")
        try:
            w = json.load(open(wf))
        except Exception:
            continue
        nm = os.path.basename(wf)[len("_ml_runtime_"):-len("_weights.json")]
        tr = w.get("_training") or {}
        runtime_files.append({"name": nm, "weights": os.path.abspath(wf),
                              "attributes": os.path.abspath(af) if os.path.isfile(af) else "",
                              "features": w.get("_features") or [],
                              "default": w.get("_default"),
                              "r2": tr.get("r2")})

    summary = {
        "path": os.path.abspath(out_path), "id_col": id_col, "n_units": len(ids),
        "mapping_file": mapping_file, "runtime_files": runtime_files,
        "columns": columns, "categorical": categorical, "numeric": numeric,
        "n_units_in_model_output": n_in_model, "n_model_units": n_model_units,
        "notes": notes, "sources": sources, "hostmodel": hostmodel,
    }
    with open(os.path.join(out_dir, "attributes_summary.json"), "w") as f:
        json.dump(summary, f, indent=1)
    logger.info(f"Layer-1 attribute table: {len(ids)} units x {len(columns)} attributes -> {out_path}")
    return summary
