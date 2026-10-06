# Copyright 2026, Diogo Costa
# This file is part of OpenWQ model.
#
# openWQ hybrid physics-ML — LAYER 1 (parameter learning / regionalization),
# Mode A (offline bake).
"""
Turn a LOW-DIMENSIONAL attribute->parameter mapping into the per-cell
``{"DEFAULT": d, "CELLS": [[icmp, ix, iy, iz, value], ...]}`` spatial-parameter
maps that openWQ's ``OpenWQ_load_param`` (Phase 0) already reads. No extra C++
is needed for Mode A — this module IS the ingestion path.

Why a mapping and not free per-cell values
-------------------------------------------
Calibrating one free value per (HRU x parameter) is an ``N_HRU x N_params``
explosion — ill-posed and untrainable. Instead the TRAINABLE object is a small
shared mapping ``theta = g(attributes)``; its dimension is FIXED and independent
of the number of HRUs/reaches. Adding cells adds constraint (data), not unknowns.

Rungs (increasing flexibility; all calibratable with the existing DDS optimizer)
-------------------------------------------------------------------------------
* ``per_class``  : ``theta_cell = class_values[class_of(cell)]``
                   trainable = one value per (soil x land-use ...) class.
* ``regression`` : ``theta_cell = intercept + sum_k coeff_k * attribute_k(cell)``
                   trainable = a handful of global coefficients.
(A ``dpl`` rung — a small NN — trains in Python by gradient descent and writes
the same maps; it is a later add. The map format below is identical for it.)

The generator is attribute-source-agnostic: it consumes plain dicts of
per-cell attributes, so the same code serves SUMMA HRUs, mizuRoute reaches, or
a gridded host. Internal (ix,iy,iz) indices come from the calibration
``ReachMapper`` (reach/HRU id -> internal index), so callers work in the
familiar host ids, never internal array indices.
"""

from typing import Dict, List, Optional, Any, Tuple
import csv
import json
import os


# ---------------------------------------------------------------------------
# Rung 2 — per-class lookup
# ---------------------------------------------------------------------------
def canon_class_key(v: Any) -> str:
    """Canonical string for a class code so ``3``, ``3.0`` and ``"3"`` (CSV
    round-trips turn integer codes into floats) all address the same class."""
    try:
        f = float(v)
        if f.is_integer():
            return str(int(f))
        return str(v)
    except (TypeError, ValueError):
        return str(v).strip()


def values_per_class(cell_classes: Dict[Any, Any],
                     class_values: Dict[Any, float]) -> Dict[Any, float]:
    """rung 2: every cell takes its class's value.

    Parameters
    ----------
    cell_classes : dict ``cell_id -> class_key``
    class_values : dict ``class_key -> value`` (the trainable object)

    Returns ``cell_id -> value``. A cell whose class is not in ``class_values``
    is omitted (it will inherit the map's DEFAULT).
    """
    cv = {canon_class_key(k): float(v) for k, v in class_values.items()}
    out: Dict[Any, float] = {}
    for cid, ck in cell_classes.items():
        key = canon_class_key(ck)
        if key in cv:
            out[cid] = cv[key]
    return out


# ---------------------------------------------------------------------------
# Rung 3 — attribute regression (pedotransfer)
# ---------------------------------------------------------------------------
def values_regression(cell_attributes: Dict[Any, Dict[str, float]],
                      coeffs: Dict[str, float],
                      intercept: float = 0.0,
                      lower: Optional[float] = None,
                      upper: Optional[float] = None) -> Dict[Any, float]:
    """rung 3: ``theta = intercept + sum_k coeff_k * attr_k``, optionally clamped.

    Parameters
    ----------
    cell_attributes : dict ``cell_id -> {attr_name: value}``
    coeffs : dict ``attr_name -> coefficient`` (trainable)
    intercept : float (trainable)
    lower, upper : optional physical bounds applied per cell after the sum.

    Returns ``cell_id -> value``. Missing attributes are treated as 0.
    """
    out: Dict[Any, float] = {}
    for cid, attrs in cell_attributes.items():
        v = float(intercept)
        for name, c in coeffs.items():
            v += float(c) * float(attrs.get(name, 0.0))
        if lower is not None:
            v = max(float(lower), v)
        if upper is not None:
            v = min(float(upper), v)
        out[cid] = v
    return out


# ---------------------------------------------------------------------------
# Rung 3 helper — standardize the regression attributes across cells
# ---------------------------------------------------------------------------
def standardize_attributes(cell_attributes: Dict[Any, Dict[str, Any]],
                           names: List[str]):
    """z-score the named attributes across the cells that carry them:
    ``z = (x - mean) / std`` (std == 0 -> 1, so a constant attribute maps to 0).

    Without this the regression ``theta = a + sum_k c_k * attr_k`` mixes raw
    units (an elevation of ~1000 m next to a slope of ~0.01), the coefficient
    bounds mean something different for every attribute, and the per-cell
    value pins to the clamp bounds — the same saturation the Layer-2 closure
    showed with a raw-mass input. With z-scores a coefficient is "change in
    theta per one standard deviation of the attribute", the intercept is theta
    at the mean attributes, and c = 0 recovers the lumped value exactly.

    Returns ``(z_attributes, stats)`` with ``stats = {name: {"mean", "std",
    "n"}}``; non-numeric values are ignored (treated as missing -> z = 0).
    """
    import math as _math
    stats: Dict[str, Dict[str, float]] = {}
    for nm in names:
        xs = []
        for attrs in cell_attributes.values():
            v = attrs.get(nm)
            try:
                fv = float(v)
            except (TypeError, ValueError):
                continue
            if _math.isfinite(fv):
                xs.append(fv)
        if not xs:
            stats[nm] = {"mean": 0.0, "std": 1.0, "n": 0}
            continue
        mean = sum(xs) / len(xs)
        var = sum((x - mean) ** 2 for x in xs) / len(xs)
        std = _math.sqrt(var)
        stats[nm] = {"mean": mean, "std": (std if std > 0 else 1.0), "n": len(xs)}
    z: Dict[Any, Dict[str, float]] = {}
    for cid, attrs in cell_attributes.items():
        row: Dict[str, float] = {}
        for nm in names:
            try:
                fv = float(attrs.get(nm))
            except (TypeError, ValueError):
                fv = None
            if fv is None or not _math.isfinite(fv):
                row[nm] = 0.0                      # missing -> at the mean
            else:
                st = stats[nm]
                row[nm] = (fv - st["mean"]) / st["std"]
        z[cid] = row
    return z, stats


# ---------------------------------------------------------------------------
# Build the openWQ spatial-parameter object
# ---------------------------------------------------------------------------
def build_spatial_param(cell_values: Dict[Any, float],
                        reach_mapper,
                        default: float,
                        icmp: int = -1,
                        skip_default_equal: bool = True,
                        tol: float = 1e-15) -> Dict[str, Any]:
    """Build the ``{"DEFAULT": d, "CELLS": [[icmp, ix, iy, iz, value], ...]}``
    object that ``OpenWQ_load_param`` reads.

    Index convention (matches every other openWQ JSON input and the
    ``xyz_elements`` written to the HDF5 output, which is where ``ReachMapper``
    takes them from): ``ix, iy, iz`` are ONE-based cell indices; ``icmp`` is the
    0-based index into the host model's compartment list, or -1 = EVERY
    compartment (the default: a BGC/SS parameter acts wherever its framework
    runs). Rows are emitted per unit COLUMN with ``iz = -1`` (every vertical
    element), so a SUMMA HRU gets the value in all its soil layers and in the
    runoff/aquifer/snow compartments too, whatever their layer counts.

    Parameters
    ----------
    cell_values : dict ``reach/HRU id -> value``
    reach_mapper : a ``ReachMapper`` (``get_xyz(id) -> (ix, iy, iz)``)
    default : baseline value for every unlisted cell
    icmp : target compartment index; -1 (default) = every compartment
    skip_default_equal : omit cells whose value equals ``default`` (compact map)

    Cells whose id cannot be mapped are skipped (they inherit ``default``).
    """
    # ``icmp`` may also be a list of compartment indices: one row per index
    # (coupled SUMMA + mizuRoute: a land value addresses every land compartment
    # explicitly, never the river compartment whose cell indices overlap)
    icmps = [int(i) for i in icmp] if isinstance(icmp, (list, tuple, set)) else [int(icmp)]
    if not icmps:
        icmps = [-1]
    cells = []
    for rid, val in cell_values.items():
        cols = resolve_cell_columns(reach_mapper, rid)
        if not cols:
            continue                                 # unmapped -> inherits DEFAULT
        if skip_default_equal and abs(float(val) - float(default)) <= tol:
            continue
        for ix, iy in cols:
            for ic in icmps:
                cells.append([int(ic), int(ix), int(iy), -1, float(val)])
    return {"DEFAULT": float(default), "CELLS": cells}


def resolve_cell_xyz(reach_mapper, rid) -> List[Tuple[int, int, int]]:
    """All internal (ix, iy, iz) cells a unit id maps to. mizuRoute ids map
    1:1; a SUMMA HRU is several vertical elements the mapper keys as
    ``<hruId>_z<k>`` (one per soil layer), so a unit id expands to ALL its
    layers (the whole HRU column gets the regionalized value)."""
    xyz = reach_mapper.get_xyz(rid)
    if xyz is None:
        xyz = reach_mapper.get_xyz(str(rid))
    if xyz is not None:
        return [xyz]
    cache = getattr(reach_mapper, "_owq_layer_index", None)
    if cache is None:
        cache = {}
        try:
            import re as _re
            for full in reach_mapper.get_all_reach_ids():
                base = _re.sub(r"_z\d+$", "", str(full))
                cache.setdefault(base, []).append(full)
        except Exception:
            pass
        try:
            setattr(reach_mapper, "_owq_layer_index", cache)
        except Exception:
            pass
    out = []
    for full in cache.get(str(rid), []):
        x = reach_mapper.get_xyz(full)
        if x is not None:
            out.append(x)
    return out


def resolve_cell_columns(reach_mapper, rid) -> List[Tuple[int, int]]:
    """The distinct (ix, iy) columns a unit id occupies (a SUMMA HRU's soil
    layers share one column; a mizuRoute reach is one column). Used to emit
    wildcard rows ``[-1, ix, iy, -1, value]`` = every compartment and every
    vertical element of that column — so the regionalized value follows the
    unit through soil/runoff/aquifer alike, whatever their layer counts."""
    cols: List[Tuple[int, int]] = []
    for ix, iy, _iz in resolve_cell_xyz(reach_mapper, rid):
        if (int(ix), int(iy)) not in cols:
            cols.append((int(ix), int(iy)))
    return cols


def model_unit_ids(reach_mapper) -> List[str]:
    """Distinct spatial UNITS in a loaded mapper (SUMMA's ``<hru>_z<k>`` layer
    elements collapse to their HRU id)."""
    import re as _re
    seen = []
    for full in reach_mapper.get_all_reach_ids():
        base = _re.sub(r"_z\d+$", "", str(full))
        if base not in seen:
            seen.append(base)
    return seen


# ---------------------------------------------------------------------------
# One-call bridge used by the calibration parameter handler
# ---------------------------------------------------------------------------
def make_spatial_param(spec: Dict[str, Any],
                       subparam_values: Dict[str, float],
                       cell_attributes: Dict[Any, Dict[str, Any]],
                       reach_mapper,
                       icmp=-1,
                       diagnostics: Optional[Dict[str, Any]] = None,
                       river_attributes: Optional[Dict[Any, Dict[str, Any]]] = None,
                       river_mapper=None,
                       river_icmp: Optional[int] = None) -> Dict[str, Any]:
    """Produce the spatial-parameter object for one regionalized parameter.

    ``icmp`` is an int or a list of compartment indices.  For a coupled SUMMA +
    mizuRoute run, ``river_attributes`` (reach id -> attributes), ``river_mapper``
    and ``river_icmp`` add the rows of the river compartment: the same rung and
    calibrated values are applied to the reaches through their own id space.

    Parameters
    ----------
    spec : the regionalization declaration, e.g.::

        {"rung": "per_class", "default": 0.05, "attribute": "soil_class"}
        {"rung": "regression", "default": 0.05,
         "attributes": ["clay_frac", "soc"], "lower": 0.001, "upper": 0.2}

    subparam_values : the calibrated low-dimensional values —
        rung 2: ``{class_key: value, ...}``;
        rung 3: ``{"intercept": v, attr_name: coeff, ...}``.
    cell_attributes : ``cell_id -> {attr_name: value}`` (rung 2 must include the
        class attribute; rung 3 the regression attributes).
    reach_mapper : ``ReachMapper`` for id -> internal index.
    icmp : target compartment index.
    diagnostics : optional dict, FILLED with what the mapping did (per-cell
        values / classes / z-scores, standardization stats, clamp counts,
        mapped vs unmapped ids) — the parameter handler writes it next to the
        eval config so the results report can show the regionalization.
    """
    rung = spec.get("rung", "per_class")
    default = float(spec.get("default", 0.0))
    diag = diagnostics if isinstance(diagnostics, dict) else {}
    diag.update({"rung": rung, "default": default, "n_cells": len(cell_attributes)})

    if rung == "per_class":
        attr = spec["attribute"]
        cell_classes = {cid: attrs[attr]
                        for cid, attrs in cell_attributes.items()
                        if attr in attrs}
        vals = values_per_class(cell_classes, subparam_values)
        diag["attribute"] = attr
        diag["cell_classes"] = {str(k): canon_class_key(v) for k, v in cell_classes.items()}
        _calib_keys = {canon_class_key(k) for k in subparam_values}
        diag["n_cells_class_not_calibrated"] = sum(
            1 for ck in cell_classes.values() if canon_class_key(ck) not in _calib_keys)

    elif rung == "regression":
        intercept = float(subparam_values.get("intercept", spec.get("intercept", 0.0)))
        coeffs = {k: v for k, v in subparam_values.items() if k != "intercept"}
        names = list(coeffs.keys())
        standardize = bool(spec.get("standardize", True))
        attr_in, stats = cell_attributes, {}
        if standardize and names:
            attr_in, stats = standardize_attributes(cell_attributes, names)
        lower, upper = spec.get("lower"), spec.get("upper")
        raw = values_regression(attr_in, coeffs, intercept, None, None)
        vals = values_regression(attr_in, coeffs, intercept, lower, upper)
        diag.update({
            "attributes": names, "standardize": standardize,
            "standardization": stats,
            "intercept": intercept,
            "coeffs": {k: float(v) for k, v in coeffs.items()},
            "lower": lower, "upper": upper,
            "n_cells_clamped": sum(1 for cid in vals
                                   if abs(vals[cid] - raw.get(cid, vals[cid])) > 0),
            "cell_attributes": {str(cid): {nm: cell_attributes[cid].get(nm) for nm in names}
                                for cid in cell_attributes},
            "cell_attributes_z": ({str(cid): {nm: float(attr_in[cid].get(nm, 0.0)) for nm in names}
                                   for cid in attr_in} if standardize else {}),
        })

    else:
        raise ValueError(f"unknown regionalization rung: {rung!r} "
                         "(expected 'per_class' or 'regression')")

    diag["cell_values"] = {str(k): float(v) for k, v in vals.items()}
    # ids that resolve to a model cell (the rest silently inherit DEFAULT)
    _mapped = {cid for cid in vals if resolve_cell_xyz(reach_mapper, cid)}
    diag["n_mapped"] = len(_mapped)
    diag["ids_unmapped"] = [str(cid) for cid in vals if cid not in _mapped][:50]
    out = build_spatial_param(vals, reach_mapper, default, icmp)
    # river compartment of a coupled run: same rung on the reach attribute table
    if river_attributes and river_mapper is not None and river_icmp is not None:
        if rung == "per_class":
            r_classes = {cid: attrs[attr] for cid, attrs in river_attributes.items() if attr in attrs}
            r_vals = values_per_class(r_classes, subparam_values)
        else:
            r_in = river_attributes
            if standardize and names:
                r_in, _ = standardize_attributes(river_attributes, names)
            r_vals = values_regression(r_in, coeffs, intercept, lower, upper)
        r_mapped = {cid for cid in r_vals if resolve_cell_xyz(river_mapper, cid)}
        diag["river_cell_values"] = {str(k): float(v) for k, v in r_vals.items()}
        diag["n_mapped_river"] = len(r_mapped)
        diag["river_icmp"] = int(river_icmp)
        out["CELLS"].extend(build_spatial_param(r_vals, river_mapper, default, int(river_icmp))["CELLS"])
    return out


# ---------------------------------------------------------------------------
# Calibration integration: expand an ML_REGIONALIZE block into DDS sub-params
# ---------------------------------------------------------------------------
def _canon_spec(ml_spec: Dict[str, Any]) -> Dict[str, Any]:
    """Normalize an ML_REGIONALIZE block into the compact ``spec`` that
    ``make_spatial_param`` consumes (carried on every sub-param)."""
    return {
        "rung":            ml_spec.get("rung", "per_class"),
        "default":         float(ml_spec.get("default", 0.0)),
        "attribute":       ml_spec.get("attribute"),
        "attributes":      ml_spec.get("attributes"),
        "lower":           ml_spec.get("lower"),
        "upper":           ml_spec.get("upper"),
        "attribute_table": ml_spec.get("attribute_table"),
        "mapping_source":  ml_spec.get("mapping_source"),
        "icmp":            int(ml_spec.get("icmp", -1)),   # -1 = every compartment
        # coupled SUMMA + mizuRoute: attribute table of the reaches (river compartment)
        "attribute_table_river": ml_spec.get("attribute_table_river"),
        # rung 3: z-score the attributes across cells (default ON; see
        # standardize_attributes) so coefficient bounds are comparable.
        "standardize":     bool(ml_spec.get("standardize", True)),
        # SS load scale of ONE species (None = every species); the handler
        # merges per-species groups into the species-keyed ML_SCALE.
        "species":         ml_spec.get("species"),
    }


def expand_regionalized_param(base_name: str,
                              json_path: List[str],
                              ml_spec: Dict[str, Any],
                              file_type: str = "bgc_regionalize") -> List[Dict[str, Any]]:
    """Expand one ``ML_REGIONALIZE`` block into a list of DDS sub-parameters.

    Each returned dict is a NORMAL calibration parameter (it has ``bounds`` so
    the optimizer varies it), tagged so the parameter handler can regroup and
    apply them together:

    * ``regionalize_group`` = ``base_name`` (the physical parameter they build)
    * ``subparam_key``      = the class name (rung 2) or coefficient name (rung 3)
    * ``regionalize_spec``  = the compact spec for ``make_spatial_param``
    * ``path``              = the shared JSON path of the physical parameter

    rung 2 (``classes``): ``{class_key: [lo, hi], ...}`` -> one sub-param per class.
    rung 3 (``coeffs``):  ``{coeff_name: [lo, hi], ...}`` (may include
    ``"intercept"``) -> one sub-param per coefficient.
    """
    spec = _canon_spec(ml_spec)
    rung = spec["rung"]

    if rung == "per_class":
        items = ml_spec.get("classes", {})
    elif rung == "regression":
        items = ml_spec.get("coeffs", {})
    else:
        raise ValueError(f"unknown regionalization rung: {rung!r}")

    subs: List[Dict[str, Any]] = []
    for key, rng in items.items():
        lo, hi = float(rng[0]), float(rng[1])
        if lo > hi:
            lo, hi = hi, lo
        # Start = the LUMPED model exactly (like the Layer-2 closure's cw=cb=0):
        # rung-2 class values start at the block default; rung-3 intercept
        # starts at the default (= theta at the mean attributes once the
        # attributes are standardized) and every coefficient at 0 (no
        # attribute dependence). Midpoint fallback if 0 / the default is
        # outside the user's bounds.
        if rung == "per_class" or key == "intercept":
            initial = spec["default"]
        else:
            initial = 0.0 if lo <= 0.0 <= hi else 0.5 * (lo + hi)
        # keep the initial inside the bounds
        initial = min(max(initial, lo), hi)
        subs.append({
            "name":              f"{base_name}@{key}",
            "file_type":         file_type,
            "path":              json_path,
            "bounds":            (lo, hi),
            "initial":           float(initial),
            "transform":         "linear",
            "regionalize_group": base_name,
            "subparam_key":      key,
            "regionalize_spec":  spec,
            "source":            "ml-regionalize",
        })
    return subs


# ---------------------------------------------------------------------------
# Layer 1, Mode B: export a trained MLP for RUNTIME inference inside openWQ
# ---------------------------------------------------------------------------
def export_mlp_weights(layers,
                       in_mean=None, in_std=None,
                       out_scale=None, out_offset=None,
                       in_transform=None) -> Dict[str, Any]:
    """Serialize a trained MLP to the JSON that openWQ's ``OpenWQ_ML_from_json``
    reads.

    Parameters
    ----------
    layers : list of ``(W, b, activation)`` — ``W`` is ``[n_out][n_in]`` (nested
        list or numpy array), ``b`` is ``[n_out]``, ``activation`` in
        ``{"tanh","relu","sigmoid","linear"}``.
    in_mean, in_std : optional per-feature input standardization (must match
        training): openWQ applies ``(x - mean)/std`` before the first layer.
    out_scale, out_offset : optional output affine de-normalization
        (``y*scale + offset`` after the last layer).
    in_transform : optional element-wise input transform applied BEFORE the
        standardization (``"log1p"`` = log(1+max(x,0)); None = identity).
    """
    def _l(a):
        return a.tolist() if hasattr(a, "tolist") else a
    w: Dict[str, Any] = {"layers": [
        {"W": _l(W), "b": _l(b), "activation": act} for (W, b, act) in layers]}
    if in_transform:
        w["in_transform"] = str(in_transform)
    if in_mean is not None and in_std is not None:
        w["in_mean"] = _l(in_mean)
        w["in_std"] = _l(in_std)
    if out_scale is not None:
        w["out_scale"] = _l(out_scale)
    if out_offset is not None:
        w["out_offset"] = _l(out_offset)
    return w


def build_ml_runtime_param(weights: Dict[str, Any],
                           cell_attributes: Dict[Any, Dict[str, Any]],
                           reach_mapper,
                           default: float,
                           feature_order: Optional[List[str]] = None,
                           icmp: int = -1,
                           weights_file: Optional[str] = None,
                           attributes_file: Optional[str] = None) -> Dict[str, Any]:
    """Build the ``{"ML_RUNTIME": {...}}`` config for a parameter (Mode B).

    openWQ evaluates the net per cell AT LOAD TIME, so the config carries the
    net + each unit's attribute vector, keyed by unit column (``ix, iy``
    ONE-based = the ``xyz_elements`` convention; ``iz = -1`` every vertical
    element; ``icmp`` 0-based compartment index or -1 = every compartment)::

        {"ML_RUNTIME": {"weights": ...,
                        "default": d,
                        "attributes": [[icmp, ix, iy, iz, a1, a2, ...], ...]}}

    ``feature_order`` fixes the attribute column order fed to the net (it MUST
    match the order the net was trained on). Unmapped ids are skipped.
    """
    rows = []
    for rid, feats in cell_attributes.items():
        if feature_order:
            vec = [float(feats.get(k, 0.0) or 0.0) for k in feature_order]
        else:
            vec = [float(v) for v in feats.values()]
        # one row per unit column; iz=-1 -> every vertical element (a SUMMA
        # HRU = all its soil layers, in every compartment when icmp=-1)
        for ix, iy in resolve_cell_columns(reach_mapper, rid):
            rows.append([int(icmp), int(ix), int(iy), -1] + vec)
    # For a REAL openWQ config, write weights/attributes to separate files and
    # reference them (openWQ's config normalizer upper-cases keys/values and
    # cannot hold the nested inline forms; *FILEPATH values are left untouched).
    if weights_file or attributes_file:
        block: Dict[str, Any] = {"DEFAULT": float(default)}
        if weights_file:
            json.dump(weights, open(weights_file, "w"))
            block["WEIGHTS_FILEPATH"] = weights_file
        else:
            block["weights"] = weights
        if attributes_file:
            json.dump(rows, open(attributes_file, "w"))
            block["ATTRIBUTES_FILEPATH"] = attributes_file
        else:
            block["attributes"] = rows
        return {"ML_RUNTIME": block}
    return {"ML_RUNTIME": {"weights": weights,
                           "default": float(default),
                           "attributes": rows}}




# ---------------------------------------------------------------------------
# Layer 1, Mode B — a small numpy MLP trainer + the matching forward pass
# ---------------------------------------------------------------------------
def forward_mlp(weights: Dict[str, Any], x) -> float:
    """Evaluate an exported net (``export_mlp_weights`` JSON) on one feature
    vector — the exact Python twin of ``OpenWQ_ML::forward`` (input transform,
    standardization, layers, output de-normalization)."""
    import numpy as _np
    h = _np.asarray(x, dtype=float).reshape(-1)
    if str(weights.get("in_transform", "")).lower() == "log1p":
        h = _np.log1p(_np.maximum(h, 0.0))
    if weights.get("in_mean") is not None and weights.get("in_std") is not None:
        h = (h - _np.asarray(weights["in_mean"], float)) / _np.asarray(weights["in_std"], float)
    for L in weights.get("layers", []):
        h = _np.asarray(L["W"], float) @ h + _np.asarray(L["b"], float)
        act = L.get("activation", "linear")
        if act == "tanh":
            h = _np.tanh(h)
        elif act == "relu":
            h = _np.maximum(h, 0.0)
        elif act == "sigmoid":
            h = 1.0 / (1.0 + _np.exp(-h))
    if weights.get("out_scale") is not None:
        h = h * _np.asarray(weights["out_scale"], float) + _np.asarray(
            weights.get("out_offset", [0.0] * len(h)), float)
    return float(h[0])


def train_mlp(X, y, hidden: int = 6, iters: int = 4000, lr: float = 0.02,
              seed: int = 0, lower: Optional[float] = None,
              upper: Optional[float] = None) -> Dict[str, Any]:
    """Fit a 1-hidden-layer tanh MLP ``y = g(x)`` by full-batch Adam on the
    mean-squared error and export it in the openWQ weights JSON.

    Inputs are z-scored (baked into ``in_mean/in_std``) and the output is
    de-normalized (``out_scale/out_offset`` = the target std/mean), so the net
    is trained on unit-scale numbers but openWQ feeds it RAW attributes and
    gets the parameter in its physical units. Dependency-free (numpy only)."""
    import numpy as _np
    X = _np.asarray(X, float)
    y = _np.asarray(y, float).reshape(-1)
    if X.ndim == 1:
        X = X.reshape(-1, 1)
    n, d = X.shape
    xm, xs = X.mean(axis=0), X.std(axis=0)
    xs = _np.where(xs > 0, xs, 1.0)
    ym, ys = float(y.mean()), float(y.std()) if float(y.std()) > 0 else 1.0
    Xz = (X - xm) / xs
    yz = (y - ym) / ys
    rng = _np.random.default_rng(seed)
    W1 = rng.normal(0, 0.5, (hidden, d)); b1 = _np.zeros(hidden)
    W2 = rng.normal(0, 0.5, (1, hidden)); b2 = _np.zeros(1)
    params = [W1, b1, W2, b2]
    m = [_np.zeros_like(p) for p in params]; v = [_np.zeros_like(p) for p in params]
    b1a, b2a, eps = 0.9, 0.999, 1e-8
    for it in range(1, iters + 1):
        H = _np.tanh(Xz @ W1.T + b1)                  # n x hidden
        out = H @ W2.T + b2                            # n x 1
        err = out.reshape(-1) - yz                     # n
        loss = float(_np.mean(err ** 2)) + 1e-4 * float(_np.sum(W1 ** 2) + _np.sum(W2 ** 2))
        g_out = (2.0 / n) * err.reshape(-1, 1)         # n x 1
        gW2 = g_out.T @ H + 2e-4 * W2; gb2 = g_out.sum(axis=0)
        gH = g_out @ W2 * (1.0 - H ** 2)               # n x hidden
        gW1 = gH.T @ Xz + 2e-4 * W1; gb1 = gH.sum(axis=0)
        for i, (p, g) in enumerate(zip(params, (gW1, gb1, gW2, gb2))):
            m[i] = b1a * m[i] + (1 - b1a) * g
            v[i] = b2a * v[i] + (1 - b2a) * g ** 2
            mh = m[i] / (1 - b1a ** it); vh = v[i] / (1 - b2a ** it)
            p -= lr * mh / (_np.sqrt(vh) + eps)
    W1, b1, W2, b2 = params
    w = export_mlp_weights([(W1, b1, "tanh"), (W2, b2, "linear")],
                           in_mean=xm, in_std=xs, out_scale=[ys], out_offset=[ym])
    # report fit quality + range on the training cells
    pred = _np.array([forward_mlp(w, xi) for xi in X])
    if lower is not None or upper is not None:
        pred = _np.clip(pred, lower if lower is not None else -_np.inf,
                        upper if upper is not None else _np.inf)
    ss_res = float(_np.sum((y - pred) ** 2)); ss_tot = float(_np.sum((y - y.mean()) ** 2)) or 1.0
    w["_training"] = {"n": int(n), "features": int(d), "hidden": int(hidden),
                      "iters": int(iters), "loss": loss, "r2": 1.0 - ss_res / ss_tot,
                      "pred_min": float(pred.min()), "pred_max": float(pred.max())}
    return w


def prepare_runtime_param(attribute_table: str,
                          targets: Dict[Any, float],
                          features: List[str],
                          reach_mapper,
                          default: float,
                          out_dir: str,
                          name: str = "param",
                          hidden: int = 6,
                          iters: int = 4000,
                          icmp: int = -1,
                          lower: Optional[float] = None,
                          upper: Optional[float] = None) -> Dict[str, Any]:
    """One-call Layer-1B preparation: train ``theta = g(features)`` on the
    cells that have a target value, then write the two files a runtime-NN
    row in the calibration report needs::

        <out_dir>/_ml_runtime_<name>_weights.json     (export_mlp_weights)
        <out_dir>/_ml_runtime_<name>_attributes.json  ([[icmp,ix,iy,iz,f1,..],...])

    Returns ``{"weights_file", "attributes_file", "feature_order", "training",
    "field": {id: value}, "n_cells"}``. Typical source of ``targets``: the
    per-cell field of a Layer-1A calibration (``_ml_regionalize_<param>.json``
    -> ``cell_values``), so the network generalizes what DDS found per class
    into a smooth attribute relation."""
    attrs = load_attribute_table(attribute_table)
    ids = [cid for cid in attrs if canon_class_key(cid) in
           {canon_class_key(t) for t in targets}]
    tmap = {canon_class_key(t): float(v) for t, v in targets.items()}
    X = [[float(attrs[cid].get(f, 0.0) or 0.0) for f in features] for cid in ids]
    y = [tmap[canon_class_key(cid)] for cid in ids]
    if len(ids) < 3:
        raise ValueError("need at least 3 cells with a target value to train")
    w = train_mlp(X, y, hidden=hidden, iters=iters, lower=lower, upper=upper)
    w["_features"] = list(features)
    w["_default"] = float(default)
    os.makedirs(out_dir, exist_ok=True)
    safe = "".join(ch if ch.isalnum() or ch in "_-." else "_" for ch in str(name))
    wfile = os.path.join(out_dir, f"_ml_runtime_{safe}_weights.json")
    afile = os.path.join(out_dir, f"_ml_runtime_{safe}_attributes.json")
    build_ml_runtime_param(w, attrs, reach_mapper, default, feature_order=list(features),
                           icmp=icmp, weights_file=wfile, attributes_file=afile)
    field = {str(cid): forward_mlp(w, [float(attrs[cid].get(f, 0.0) or 0.0) for f in features])
             for cid in attrs}
    return {"weights_file": os.path.abspath(wfile), "attributes_file": os.path.abspath(afile),
            "feature_order": list(features), "training": w.get("_training", {}),
            "field": field, "n_cells": len(field)}


# ---------------------------------------------------------------------------
# Attribute-table loading (per-cell attributes: id -> {attr: value})
# ---------------------------------------------------------------------------
def load_attribute_table(path: str,
                         id_column: str = "id") -> Dict[Any, Dict[str, Any]]:
    """Load a per-cell attribute table into ``{id -> {attr_name: value}}``.

    Supports CSV (header row, one column is the id) and JSON (either a mapping
    ``{id: {attr: value}}`` or a list of row dicts). Numeric-looking CSV values
    are converted to float; everything else stays a string (so categorical
    class labels for rung 2 are preserved).
    """
    if not path or not os.path.isfile(path):
        return {}

    if str(path).lower().endswith(".json"):
        with open(path) as f:
            data = json.load(f)
        if isinstance(data, dict):
            return {k: dict(v) for k, v in data.items()}
        out: Dict[Any, Dict[str, Any]] = {}
        for row in data:                       # list of row dicts
            rid = row.get(id_column)
            out[rid] = {k: v for k, v in row.items() if k != id_column}
        return out

    # CSV
    out = {}
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            rid = row.get(id_column)
            if rid is None:
                continue
            attrs: Dict[str, Any] = {}
            for k, v in row.items():
                if k == id_column:
                    continue
                try:
                    attrs[k] = float(v)
                except (TypeError, ValueError):
                    attrs[k] = v               # categorical (e.g. soil_class)
            # ids are often integers; keep both the raw and int-normalized form
            try:
                out[int(rid)] = attrs
            except (TypeError, ValueError):
                out[rid] = attrs
    return out
