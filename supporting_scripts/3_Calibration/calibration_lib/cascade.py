# Copyright 2026, Diogo Costa
# This file is part of OpenWQ model.
#
# Sub-basin ("cascade") calibration — stage driver.
"""Upstream-to-downstream calibration of water-quality zones.

``run_cascade`` receives the arguments of :func:`calibration_driver.run_calibration`
plus ``cascade={...}`` and runs a sequence of ordinary calibrations in
sub-folders of the work directory:

* ``stage_00_global``: every selected parameter lumped (one value for the whole
  domain), scored at all target stations. Baseline and starting point.
* ``stage_NN_levelL``: one stage per topological level of the gauged zones. The
  zone-specific parameters of the zones in that level are free (one value per
  zone, starting at the global value); zones already calibrated keep their
  values, zones not yet calibrated and the ungauged zone use the inherited
  value; the non-zone parameters keep their global values. Scored at the
  stations of the level.
* ``stage_NN_polish`` (optional): every zone value free, scored at all stations.

Zone-specific values travel as Layer-1 ``per_class`` sub-parameters on the
``wq_zone`` attribute (see :mod:`wq_zones`), so the parameter handler writes
them as the usual per-cell maps and the ML tab of the report can combine them
with any other regionalization. Results of the last stage are copied to
``<work_dir>/results`` so the standard results report works, and
``<work_dir>/wq_zones/cascade_summary.json`` holds the per-zone values and fits.

``cascade`` keys (all optional except ``enabled``):
    min_obs (20), target_reaches [..], excluded_reaches [..],
    zone_parameters [names], global_evaluations, stage_evaluations,
    polish_evaluations (0 = none), inherit ("global" | "upstream" | "downstream_most")
"""

import copy
import json
import logging
import os
import shutil
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np

try:
    from . import wq_zones
    from . import extract_parameters as _ep
except ImportError:                                  # pragma: no cover
    import wq_zones
    import extract_parameters as _ep

logger = logging.getLogger(__name__)

STAGE_PARAMS_FILE = "cascade_stage_params.json"


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
# Parameter types that can carry one value per zone: reaction-network rate
# parameters and the sink/source load scale (written as per-cell maps). The
# other types (climate response and shape of the loads, sediment, lateral
# exchange) stay global.
_ZONE_ELIGIBLE_PREFIXES = ("bgc_json", "ss_csv_scale")


def zone_eligible(param: Dict[str, Any]) -> bool:
    ft = str(param.get("file_type", ""))
    return any(ft.startswith(x) for x in _ZONE_ELIGIBLE_PREFIXES)


def _untag(name: str):
    head, sep, rest = name.partition(":")
    if sep and head.startswith("m") and head[1:].isdigit():
        return int(head[1:]), rest
    return None, name


def _zone_spec(param: Dict[str, Any], classes: Dict[str, List[float]], default: float,
               table_path: str, river_table_path: Optional[str] = None) -> Dict[str, Any]:
    """The ML_REGIONALIZE block that makes ``param`` zone-specific.

    ``river_table_path`` (coupled SUMMA + mizuRoute) is the reach table: the
    land units follow ``table_path`` and the reaches this one, each in its own
    compartment family."""
    ft = str(param.get("file_type", ""))
    spec: Dict[str, Any] = {
        "rung": "per_class",
        "attribute": wq_zones.ZONE_ATTRIBUTE,
        "attribute_table": table_path,
        "id_column": "id",
        "default": float(default),
        "classes": classes,
        "icmp": -1,
    }
    if river_table_path:
        spec["attribute_table_river"] = river_table_path
    if ft.startswith("ss_"):
        spec["module"] = "ss"
        if param.get("species"):
            spec["species"] = param["species"]
    elif ft.startswith("ts_"):
        spec["module"] = "ts"
        spec["ts_param"] = param.get("ts_param")
    elif ft.startswith("module_") or param.get("module_key"):
        spec["module"] = param.get("module", "td")
        spec["module_key"] = param.get("module_key")
    else:
        spec["module"] = "bgc"
    if param.get("path"):
        spec["path"] = param["path"]
    return spec


def _station_fit(stage_dir: Path, metric_name: str) -> Dict[str, Dict[str, float]]:
    """Per-station fit of the best evaluation of a stage from its matched data."""
    out: Dict[str, Dict[str, float]] = {}
    f = stage_dir / "results" / "matched_data.csv"
    if not f.is_file():
        return out
    try:
        import pandas as pd
        try:
            from .objective_functions import ObjectiveFunction as _OF
        except ImportError:                          # pragma: no cover
            from objective_functions import ObjectiveFunction as _OF
        df = pd.read_csv(f)
        for rid, g in df.groupby("reach_id"):
            o, s = g["observed"].values.astype(float), g["simulated"].values.astype(float)
            m = ~(np.isnan(o) | np.isnan(s))
            if m.sum() < 2:
                continue
            key = str(int(float(rid)))
            out[key] = {"n": int(m.sum()), "kge": float(_OF.kge(o[m], s[m])),
                        "pbias": float(_OF.pbias(o[m], s[m])) if hasattr(_OF, "pbias") else float("nan")}
    except Exception as e:                           # the fit table is a convenience only
        logger.debug(f"station fit not computed for {stage_dir}: {e}")
    return out


def _best_params(stage_dir: Path) -> Dict[str, float]:
    f = stage_dir / "results" / "best_parameters.json"
    if f.is_file():
        with open(f) as fh:
            return {k: float(v) for k, v in json.load(fh).items()}
    return {}


def _stage_done(stage_dir: Path) -> bool:
    return (stage_dir / "results" / "best_parameters.json").is_file()


# ---------------------------------------------------------------------------
# main entry
# ---------------------------------------------------------------------------
def run_cascade(**call) -> Dict[str, Any]:
    try:
        from .calibration_driver import run_calibration
    except ImportError:                              # pragma: no cover
        from calibration_driver import run_calibration

    cfg = dict(call.pop("cascade"))
    work_dir = Path(call["calibration_work_dir"])
    work_dir.mkdir(parents=True, exist_ok=True)
    resume = bool(call.get("resume"))
    params: List[Dict[str, Any]] = list(call.get("calibration_parameters") or [])
    targets = dict(call.get("calibration_targets") or {})
    species = list(targets.get("species") or [])
    period = call.get("calibration_period")

    log_file = work_dir / "cascade.log"
    fh = logging.FileHandler(log_file)
    fh.setFormatter(logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s"))
    logging.getLogger().addHandler(fh)

    # ── model with the river network: the last model of a chain ──
    chain = call.get("model_chain")
    model_cfg = call.get("model_config")
    network_cfg = chain[-1] if chain else model_cfg
    network_index = (len(chain) - 1) if chain else 0

    # ── zones ──
    obs_csv = str(work_dir / "calibration_observations.csv")
    if not os.path.isfile(obs_csv):
        try:
            try:
                from . import config_integration as _ci
            except ImportError:                      # pragma: no cover
                import config_integration as _ci
            _ci.prepare_calibration_observations_csv(network_cfg, str(work_dir),
                                                     log=lambda *a, **k: None)
        except Exception as e:
            logger.warning(f"cascade: could not prepare the observation CSV: {e}")
    zones = wq_zones.load_zones(str(work_dir)) if resume else None
    if zones is None:
        zones = wq_zones.delineate(
            network_cfg, obs_csv, species,
            min_obs=int(cfg.get("min_obs", 20)),
            period=tuple(period) if period else None,
            target_reaches=cfg.get("target_reaches") or None,
            excluded_reaches=cfg.get("excluded_reaches") or None)
    tables = wq_zones.write_zone_tables(zones, str(work_dir / wq_zones.ZONES_DIRNAME))
    if not zones["zones"]:
        raise ValueError("cascade: no station qualifies as a calibration target "
                         "(check min_obs and the target reaches)")
    station_of = {z["id"]: int(z["station_reach"]) for z in zones["zones"]}
    all_stations = sorted(station_of.values())
    logger.info("Cascade: %d zones in %d levels: %s", len(zones["zones"]), len(zones["levels"]),
                "; ".join(f"L{i+1}: zones {lv}" for i, lv in enumerate(zones["levels"])))

    # ── which parameters become zone-specific ──
    wanted = [str(n) for n in (cfg.get("zone_parameters") or [])]
    by_name = {p["name"]: p for p in params}
    zone_params: List[Dict[str, Any]] = []
    for n in wanted:
        p = by_name.get(n)
        if p is None:
            # accept an un-tagged name for a chain-tagged parameter
            p = next((q for q in params if _untag(q["name"])[1] == n), None)
        if p is None:
            logger.warning(f"cascade: zone parameter '{n}' is not a calibration parameter, ignored")
        elif p.get("regionalize_group"):
            logger.warning(f"cascade: '{n}' is already regionalized in the ML tab, kept as is")
        elif not zone_eligible(p):
            logger.warning(f"cascade: '{n}' ({p.get('file_type')}) has one value per domain and "
                           "cannot vary by zone, kept global")
        else:
            zone_params.append(p)
    if not zone_params:
        logger.warning("cascade: no zone-specific parameter selected; every stage after the "
                       "global one would be empty. Running the global stage only.")

    # SUMMA with internally coupled mizuRoute: one model holds land units (HRU
    # ids) and reaches (segIds); zone values go to both through their own tables
    _coupled = (chain is None and str((network_cfg or {}).get("hostmodel") or "").lower() == "summa"
                and bool((network_cfg or {}).get("mizuroute_config_path")))
    if _coupled:
        logger.info("Cascade: coupled land + river model, zone values applied to the HRUs "
                    "(land compartments) and to the reaches (river compartment)")

    def _table_for(p):
        mi = p.get("model_index")
        if mi is None:
            mi = _untag(p["name"])[0] or 0
        if _coupled:
            return tables["hrus"]
        return tables["reaches"] if (chain is None or mi == network_index) else tables["hrus"]

    def _river_table_for(p):
        return tables["reaches"] if _coupled else None

    # bookkeeping across stages
    zone_values: Dict[str, Dict[int, float]] = {p["name"]: {} for p in zone_params}
    summary: Dict[str, Any] = {"schema": "openwq_cascade/1", "zones": zones["zones"],
                               "ungauged": {k: v for k, v in zones["ungauged"].items()
                                            if k not in ("reaches", "hrus")},
                               "levels": zones["levels"], "stages": [],
                               "zone_parameters": [p["name"] for p in zone_params],
                               "inherit": cfg.get("inherit", "global")}
    summary_path = work_dir / wq_zones.ZONES_DIRNAME / "cascade_summary.json"

    def _save_summary():
        with open(summary_path, "w") as f:
            json.dump(summary, f, indent=2, default=float)

    # number of stages of this run, for the "cascade k/N" progress text
    n_stages_total = 1 + (len(zones["levels"]) if zone_params else 0) \
        + (1 if (zone_params and int(cfg.get("polish_evaluations") or 0) > 0) else 0)

    def _progress(stage_dir: Path, label: str, reach_ids) -> str:
        k = int(stage_dir.name.split("_")[1]) + 1
        if label == "global":
            what = "global"
        elif label.startswith("level"):
            zs = [z for z in zones["zones"] if z["station_reach"] in set(int(r) for r in reach_ids)]
            what = "zone " + ", ".join(str(z["id"]) for z in zs)
        else:
            what = "polish (all zones)"
        return (f"cascade stage {k}/{n_stages_total} · {what} · stations "
                + ", ".join(str(r) for r in reach_ids))

    def _stage_call(stage_dir: Path, stage_params, fixed, reach_ids, n_evals, label, start_values=None):
        # the reports rebuild a stage from this file (its own parameter set)
        stage_dir.mkdir(parents=True, exist_ok=True)
        try:
            with open(stage_dir / STAGE_PARAMS_FILE, "w") as f:
                json.dump({"label": label, "stations": list(reach_ids),
                           "calibration_parameters": stage_params}, f, indent=1, default=str)
        except Exception as e:
            logger.debug(f"cascade: stage parameter file not written: {e}")
        summary["active_stage"] = stage_dir.name
        _save_summary()
        c = copy.copy(call)
        c["calibration_work_dir"] = str(stage_dir)
        c["calibration_parameters"] = stage_params
        c["fixed_parameters"] = fixed
        c["calibration_targets"] = {**targets, "reach_ids": list(reach_ids)}
        c["max_evaluations"] = int(n_evals)
        c["_cascade_stage"] = label
        c["cascade_progress"] = _progress(stage_dir, label, reach_ids)
        logger.info("Cascade progress: %s", c["cascade_progress"])
        c["resume"] = resume and stage_dir.exists()
        c.pop("cascade", None)
        return run_calibration(**c)

    # ── stage 0: global ──
    n_global = int(cfg.get("global_evaluations") or call.get("max_evaluations") or 100)
    g_dir = work_dir / "stage_00_global"
    if resume and _stage_done(g_dir):
        logger.info("Cascade stage 00 (global): already complete, reusing")
    else:
        logger.info("=" * 60)
        logger.info("CASCADE STAGE 00: global calibration (%d parameters, stations %s)",
                    len(params), all_stations)
        logger.info("=" * 60)
        _stage_call(g_dir, params, [], all_stations, n_global, "global")
    global_best = _best_params(g_dir)
    summary["stages"].append({"stage": "stage_00_global", "level": 0, "zones": [],
                              "stations": all_stations, "best": global_best,
                              "fit": _station_fit(g_dir, call.get("objective_function", "KGE"))})
    _save_summary()
    global_fixed = [{**p, "value": global_best.get(p["name"], p["initial"])} for p in params]

    # ── one stage per level ──
    n_stage = int(cfg.get("stage_evaluations") or call.get("max_evaluations") or 100)
    inherit = cfg.get("inherit", "global")
    stage_no = 1
    for level, zone_ids in enumerate(zones["levels"], start=1):
        if not zone_params:
            break
        stage_dir = work_dir / f"stage_{stage_no:02d}_level{level}"
        stations = sorted(station_of[z] for z in zone_ids)
        free: List[Dict[str, Any]] = []
        fixed: List[Dict[str, Any]] = [dict(f) for f in global_fixed if f["name"] not in zone_values]
        for p in zone_params:
            gval = float(global_best.get(p["name"], p["initial"]))
            lo, hi = float(p["bounds"][0]), float(p["bounds"][1])
            classes = {str(z): [lo, hi] for z in zone_ids}
            spec = _zone_spec(p, classes, gval, _table_for(p), _river_table_for(p))
            subs = _ep.apply_regionalize_to_params([dict(p)], {p["name"]: spec})
            subs = [s for s in subs if s.get("regionalize_group")]
            for s in subs:
                s["initial"] = gval
            free.extend(subs)
            # frozen classes: calibrated zones at their values, zone 0 at the inherited value
            frozen = dict(zone_values[p["name"]])
            frozen[0] = wq_zones.inherit_value(zones, zone_values[p["name"]], gval, inherit)
            fclasses = {str(z): [v, v] for z, v in frozen.items()}
            if fclasses:
                fspec = _zone_spec(p, fclasses, gval, _table_for(p), _river_table_for(p))
                for s in _ep.apply_regionalize_to_params([dict(p)], {p["name"]: fspec}):
                    if s.get("regionalize_group"):
                        s["value"] = float(frozen[int(s["subparam_key"])])
                        fixed.append(s)
        if resume and _stage_done(stage_dir):
            logger.info("Cascade stage %02d (level %d): already complete, reusing", stage_no, level)
        else:
            logger.info("=" * 60)
            logger.info("CASCADE STAGE %02d: level %d, zones %s, stations %s, %d free values",
                        stage_no, level, zone_ids, stations, len(free))
            logger.info("=" * 60)
            _stage_call(stage_dir, free, fixed, stations, n_stage, f"level{level}")
        best = _best_params(stage_dir)
        for s in free:
            if s["name"] in best:
                zone_values[s["regionalize_of"]][int(s["subparam_key"])] = best[s["name"]]
        summary["stages"].append({"stage": stage_dir.name, "level": level, "zones": zone_ids,
                                  "stations": stations, "best": best,
                                  "fit": _station_fit(stage_dir, call.get("objective_function", "KGE"))})
        summary["zone_values"] = {k: {str(z): v for z, v in d.items()} for k, d in zone_values.items()}
        _save_summary()
        stage_no += 1

    # ── optional polish: every zone value free, all stations ──
    n_polish = int(cfg.get("polish_evaluations") or 0)
    last_dir = work_dir / (f"stage_{stage_no-1:02d}_level{len(zones['levels'])}" if zone_params else "stage_00_global")
    if zone_params and n_polish > 0:
        stage_dir = work_dir / f"stage_{stage_no:02d}_polish"
        free, fixed = [], [dict(f) for f in global_fixed if f["name"] not in zone_values]
        for p in zone_params:
            gval = float(global_best.get(p["name"], p["initial"]))
            lo, hi = float(p["bounds"][0]), float(p["bounds"][1])
            classes = {str(z): [lo, hi] for z in zone_values[p["name"]]}
            classes["0"] = [lo, hi]
            spec = _zone_spec(p, classes, gval, _table_for(p), _river_table_for(p))
            for s in _ep.apply_regionalize_to_params([dict(p)], {p["name"]: spec}):
                if not s.get("regionalize_group"):
                    continue
                z = int(s["subparam_key"])
                s["initial"] = float(zone_values[p["name"]].get(z,
                                     wq_zones.inherit_value(zones, zone_values[p["name"]], gval, inherit)))
                free.append(s)
        if resume and _stage_done(stage_dir):
            logger.info("Cascade polish stage: already complete, reusing")
        else:
            logger.info("CASCADE STAGE %02d: polish, %d free values, stations %s", stage_no, len(free), all_stations)
            _stage_call(stage_dir, free, fixed, all_stations, n_polish, "polish")
        best = _best_params(stage_dir)
        for s in free:
            if s["name"] in best:
                zone_values[s["regionalize_of"]][int(s["subparam_key"])] = best[s["name"]]
        summary["stages"].append({"stage": stage_dir.name, "level": "polish", "zones": [z["id"] for z in zones["zones"]],
                                  "stations": all_stations, "best": best,
                                  "fit": _station_fit(stage_dir, call.get("objective_function", "KGE"))})
        summary["zone_values"] = {k: {str(z): v for z, v in d.items()} for k, d in zone_values.items()}
        last_dir = stage_dir
    _save_summary()

    # ── expose the final stage as the run's results ──
    res_dir = work_dir / "results"
    if (last_dir / "results").is_dir():
        res_dir.mkdir(exist_ok=True)
        for f in (last_dir / "results").iterdir():
            if f.is_file():
                shutil.copy2(f, res_dir / f.name)
    combined = dict(global_best)
    for pname, d in zone_values.items():
        for z, v in d.items():
            combined[f"{pname}@zone{z}"] = v
    with open(res_dir / "best_parameters_cascade.json", "w") as f:
        json.dump(combined, f, indent=2)
    (work_dir / "evaluations").mkdir(exist_ok=True)
    logger.info("Cascade complete: %d stages, summary in %s", len(summary["stages"]), summary_path)
    logging.getLogger().removeHandler(fh)

    last_best = _best_params(last_dir)
    return {"best_params": combined, "best_objective": None, "results_dir": str(res_dir),
            "cascade": summary, "last_stage": str(last_dir),
            "n_evaluations": sum(int(cfg.get("stage_evaluations") or 0) for _ in zones["levels"]) + n_global,
            "last_stage_best": last_best}
