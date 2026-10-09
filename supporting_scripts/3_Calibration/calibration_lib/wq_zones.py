# Copyright 2026, Diogo Costa
# This file is part of OpenWQ model.
#
# Sub-basin ("cascade") calibration — zone delineation.
"""Delineate water-quality calibration zones from the river network and the
observation stations, for the upstream-to-downstream cascade calibration.

A *zone* is the set of reaches (and the land units draining to them) that lie
upstream of a station but not upstream of any other station. Zones nest along
the network, so they can be calibrated in topological order: the zones with no
gauged zone upstream first (level 1), then the zones whose upstream zones are
all calibrated (level 2), and so on. Reaches below the last station, or on
tributaries without a station, form the *ungauged* zone 0, whose parameters are
inherited (see :func:`inherit_value`).

The hydrology sub-basins need not match the zones: a zone is simply a union of
model units (reaches and the HRUs that drain to them, from ``hruToSegId``).
Zones are written as per-unit tables (``zone_reaches.csv`` for the routing
model, ``zone_hrus.csv`` for the land model) with an ``id`` and a ``wq_zone``
column, i.e. an attribute table that the Layer-1 ``per_class`` regionalization
can consume directly. A zone-specific parameter is therefore nothing but a
per-class parameter on the ``wq_zone`` attribute.
"""

import json
import logging
import os
import re
from typing import Any, Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)

ZONES_DIRNAME = "wq_zones"
ZONE_ATTRIBUTE = "wq_zone"


# ---------------------------------------------------------------------------
# River-network topology
# ---------------------------------------------------------------------------
def _control_file(model_config: Dict[str, Any]) -> str:
    return (model_config.get("file_manager_path")
            or model_config.get("control_file_path") or "")


def find_topology(model_config: Dict[str, Any]) -> Optional[Dict[str, str]]:
    """Locate the mizuRoute network topology of a model configuration.

    Reads ``<fname_ntopOld>`` and the ``<varname_*>`` entries of the control
    file; the topology lives in the control file's own directory (the
    ``<ancil_dir>`` may be a container path). Returns ``None`` when the model
    has no river network (e.g. a SUMMA-only configuration).
    """
    # SUMMA with internally coupled mizuRoute: the network is named by the TOML
    toml = model_config.get("mizuroute_config_path") or ""
    if toml and str(model_config.get("hostmodel") or "").lower() == "summa" and os.path.isfile(toml):
        return _find_topology_toml(model_config, toml)

    fm = _control_file(model_config)
    if not fm or not os.path.isfile(fm):
        return None
    text = open(fm, encoding="utf-8", errors="ignore").read()
    if "<fname_ntopOld>" not in text and "<ancil_dir>" not in text:
        return None

    def _tag(name, default):
        m = re.search(r"<%s>\s*([^\n!<]+)" % name, text)
        return m.group(1).strip() if m else default

    topo_name = _tag("fname_ntopOld", "topology.nc")
    cands = [os.path.join(os.path.dirname(os.path.abspath(fm)), topo_name)]
    ancil = _tag("ancil_dir", "")
    if ancil:
        cands.append(os.path.join(ancil, topo_name))
        # the control file may hold the container path of the ancillary folder;
        # try the same folder name under the domain of the model configuration
        _dom = model_config.get("dir2save_input_files") or ""
        if _dom:
            _domain_dir = os.path.dirname(os.path.abspath(_dom))
            _tail = os.path.basename(ancil.rstrip("/"))
            cands.append(os.path.join(_domain_dir, "settings", _tail, topo_name))
            cands.append(os.path.join(_domain_dir, _tail, topo_name))
    topo = next((c for c in cands if os.path.isfile(c)), None)
    if topo is None:
        _dom = model_config.get("dir2save_input_files") or ""
        if _dom:
            _domain_dir = os.path.dirname(os.path.abspath(_dom))
            for root, _dirs, files in os.walk(_domain_dir):
                if root.count(os.sep) - _domain_dir.count(os.sep) > 3:
                    continue
                if topo_name in files:
                    topo = os.path.join(root, topo_name)
                    break
    if topo is None:
        return None
    return {
        "path": topo,
        "seg_id": _tag("varname_segId", "segId"),
        "down_id": _tag("varname_downSegId", "downSegId"),
        "hru_id": _tag("varname_hruid", "hruId"),
        "hru_seg": _tag("varname_hruSegId", "hruToSegId"),
        "hru_area": _tag("varname_area", "area"),
    }


def _find_topology_toml(model_config: Dict[str, Any], toml: str) -> Optional[Dict[str, str]]:
    """Topology of a coupled SUMMA + mizuRoute run from its TOML
    (``[hydrofabric] hfabric_path / hfabric_file / varname_*``); the path in the
    TOML is usually a container path, so the same folder is also looked for
    under the domain of the model configuration."""
    text = open(toml, encoding="utf-8", errors="ignore").read()

    def _key(name, default):
        m = re.search(r'^\s*%s\s*=\s*"?([^"\n#]+)"?' % re.escape(name), text, re.M)
        return m.group(1).strip() if m else default

    fab_dir = _key("hfabric_path", "")
    fab_file = _key("hfabric_file", "topology.nc")
    cands = [os.path.join(fab_dir, fab_file) if fab_dir else fab_file,
             os.path.join(os.path.dirname(os.path.abspath(toml)), fab_file)]
    _dom = model_config.get("dir2save_input_files") or ""
    if _dom:
        _domain_dir = os.path.dirname(os.path.abspath(_dom))
        _tail = os.path.basename(fab_dir.rstrip("/")) if fab_dir else ""
        if _tail:
            cands.append(os.path.join(_domain_dir, "settings", _tail, fab_file))
            cands.append(os.path.join(_domain_dir, _tail, fab_file))
    topo = next((c for c in cands if os.path.isfile(c)), None)
    if topo is None and _dom:
        _domain_dir = os.path.dirname(os.path.abspath(_dom))
        for root, _dirs, files in os.walk(_domain_dir):
            if root.count(os.sep) - _domain_dir.count(os.sep) > 3:
                continue
            if fab_file in files:
                topo = os.path.join(root, fab_file)
                break
    if topo is None:
        return None
    return {
        "path": topo,
        "seg_id": _key("varname_segId", "segId"),
        "down_id": _key("varname_downSegId", "downSegId"),
        "hru_id": _key("varname_HRUid", "hruId"),
        "hru_seg": _key("varname_hruSegId", "hruToSegId"),
        "hru_area": _key("varname_area", "area"),
    }


def read_topology(info: Dict[str, str]) -> Dict[str, Any]:
    """Read reach ids, downstream ids and the HRU -> reach mapping."""
    import netCDF4 as nc
    import numpy as np
    ds = nc.Dataset(info["path"])
    try:
        def _arr(name):
            return np.asarray(ds[name][:]).ravel() if name in ds.variables else None
        seg = _arr(info["seg_id"]).astype("int64")
        down = _arr(info["down_id"]).astype("int64")
        hru = _arr(info["hru_id"])
        hru_seg = _arr(info["hru_seg"])
        area = _arr(info["hru_area"])
        out = {
            "seg_ids": [int(s) for s in seg],
            "down_ids": {int(s): int(d) for s, d in zip(seg, down)},
            "hru_ids": [int(h) for h in hru] if hru is not None else [],
            "hru_seg": ({int(h): int(s) for h, s in zip(hru, hru_seg)}
                        if hru is not None and hru_seg is not None else {}),
            "hru_area": ({int(h): float(a) for h, a in zip(hru, area)}
                         if hru is not None and area is not None else {}),
        }
    finally:
        ds.close()
    return out


def upstream_reaches(down_ids: Dict[int, int], outlet: int) -> set:
    """All reaches draining to ``outlet`` (inclusive)."""
    ups = {}
    for s, d in down_ids.items():
        ups.setdefault(d, []).append(s)
    seen, stack = {outlet}, [outlet]
    while stack:
        r = stack.pop()
        for u in ups.get(r, []):
            if u not in seen:
                seen.add(u)
                stack.append(u)
    return seen


# ---------------------------------------------------------------------------
# Observations
# ---------------------------------------------------------------------------
def station_obs_counts(obs_csv: str, species: Optional[List[str]] = None,
                       period: Optional[Tuple[str, str]] = None,
                       primary_only: bool = True) -> Dict[int, Dict[str, Any]]:
    """Observations per station reach: ``{reach_id: {"n_obs", "first", "last", "n_stations"}}``."""
    import pandas as pd
    if not obs_csv or not os.path.isfile(obs_csv):
        return {}
    df = pd.read_csv(obs_csv)
    if "reach_id" not in df.columns or "datetime" not in df.columns:
        return {}
    df["datetime"] = pd.to_datetime(df["datetime"], errors="coerce")
    df = df.dropna(subset=["datetime", "reach_id"])
    if species and "species" in df.columns:
        df = df[df["species"].astype(str).isin([str(s) for s in species])]
    if primary_only and "is_primary" in df.columns:
        prim = df["is_primary"].fillna(True).astype(bool)
        if prim.any():
            df = df[prim]
    if period and period[0] and period[1]:
        df = df[(df["datetime"] >= pd.to_datetime(period[0]))
                & (df["datetime"] <= pd.to_datetime(period[1]))]
    out = {}
    for rid, g in df.groupby("reach_id"):
        try:
            key = int(float(rid))
        except (TypeError, ValueError):
            continue
        out[key] = {
            "n_obs": int(len(g)),
            "first": str(g["datetime"].min().date()),
            "last": str(g["datetime"].max().date()),
            "n_stations": int(g["source"].nunique()) if "source" in g.columns else 1,
        }
    return out


# ---------------------------------------------------------------------------
# Delineation
# ---------------------------------------------------------------------------
def delineate(model_config: Dict[str, Any],
              obs_csv: str,
              species: Optional[List[str]] = None,
              min_obs: int = 20,
              period: Optional[Tuple[str, str]] = None,
              target_reaches: Optional[List[int]] = None,
              excluded_reaches: Optional[List[int]] = None) -> Dict[str, Any]:
    """Build the zones of the cascade calibration.

    ``model_config`` is the configuration of the model that holds the river
    network (the last model of a chain). Stations are the reaches with
    observations; a station is a calibration target when it has at least
    ``min_obs`` observations (in ``period``, when given), unless the caller
    fixes the targets with ``target_reaches`` / ``excluded_reaches``.
    """
    info = find_topology(model_config)
    if info is None:
        raise ValueError("the model has no river-network topology (mizuRoute control file "
                         "with <fname_ntopOld>); the cascade needs a routed network")
    topo = read_topology(info)
    down = topo["down_ids"]
    counts = station_obs_counts(obs_csv, species, period)

    # candidate stations
    stations = [r for r in counts if r in down]
    rejected = [{"station_reach": int(r), "n_obs": counts[r]["n_obs"],
                 "reason": "reach not in the river network"}
                for r in counts if r not in down]
    if target_reaches:
        targets = [int(r) for r in target_reaches if int(r) in down]
        for r in stations:
            if r not in targets:
                rejected.append({"station_reach": int(r), "n_obs": counts[r]["n_obs"],
                                 "reason": "not selected as a target"})
    else:
        targets = []
        for r in stations:
            if excluded_reaches and int(r) in [int(x) for x in excluded_reaches]:
                rejected.append({"station_reach": int(r), "n_obs": counts[r]["n_obs"],
                                 "reason": "excluded by the user"})
            elif counts[r]["n_obs"] < int(min_obs):
                rejected.append({"station_reach": int(r), "n_obs": counts[r]["n_obs"],
                                 "reason": f"fewer than {int(min_obs)} observations"})
            else:
                targets.append(int(r))

    ups = {r: upstream_reaches(down, r) for r in targets}
    # reaches of each zone: upstream of the station minus upstream of the stations above it
    zones = []
    for r in targets:
        above = [c for c in targets if c != r and c in ups[r]]
        own = set(ups[r])
        for c in above:
            own -= ups[c]
        zones.append({"station_reach": r, "reaches": sorted(own),
                      "upstream_stations": sorted(above)})
    # nesting: immediate upstream zones and level
    for z in zones:
        # a station c is immediately upstream when no other upstream station lies between
        z["upstream_zones_stations"] = [
            c for c in z["upstream_stations"]
            if not any((c in ups[d]) and (d != c) for d in z["upstream_stations"])]
    level_of: Dict[int, int] = {}

    def _level(z):
        r = z["station_reach"]
        if r in level_of:
            return level_of[r]
        if not z["upstream_stations"]:
            level_of[r] = 1
        else:
            level_of[r] = 1 + max(_level(next(y for y in zones if y["station_reach"] == c))
                                  for c in z["upstream_zones_stations"])
        return level_of[r]

    for z in zones:
        _level(z)
    # order zones by level, then by drainage area
    zones.sort(key=lambda z: (level_of[z["station_reach"]], len(ups[z["station_reach"]])))
    zone_id = {z["station_reach"]: i + 1 for i, z in enumerate(zones)}
    hru_area = topo["hru_area"]
    hru_seg = topo["hru_seg"]
    seg_hrus: Dict[int, List[int]] = {}
    for h, s in hru_seg.items():
        seg_hrus.setdefault(s, []).append(h)

    reach_zone = {int(s): 0 for s in topo["seg_ids"]}
    for z in zones:
        zid = zone_id[z["station_reach"]]
        z["id"] = zid
        z["level"] = level_of[z["station_reach"]]
        z["upstream_zones"] = [zone_id[c] for c in z["upstream_zones_stations"]]
        z["n_reaches"] = len(z["reaches"])
        z["hrus"] = sorted(h for r in z["reaches"] for h in seg_hrus.get(r, []))
        z["area_m2"] = float(sum(hru_area.get(h, 0.0) for h in z["hrus"]))
        z["drainage_area_m2"] = float(sum(hru_area.get(h, 0.0)
                                          for r in ups[z["station_reach"]]
                                          for h in seg_hrus.get(r, [])))
        z["n_obs"] = counts[z["station_reach"]]["n_obs"]
        z["obs_first"] = counts[z["station_reach"]]["first"]
        z["obs_last"] = counts[z["station_reach"]]["last"]
        z["role"] = "target"
        for r in z["reaches"]:
            reach_zone[int(r)] = zid
        z.pop("upstream_zones_stations", None)
    for z in zones:
        d = down.get(z["station_reach"])
        z["downstream_zone"] = reach_zone.get(d, 0) if d in reach_zone else 0
    hru_zone = {int(h): reach_zone.get(int(s), 0) for h, s in hru_seg.items()}

    ung_reaches = sorted(r for r, zid in reach_zone.items() if zid == 0)
    ungauged = {
        "id": 0,
        "reaches": ung_reaches,
        "hrus": sorted(h for r in ung_reaches for h in seg_hrus.get(r, [])),
        "n_reaches": len(ung_reaches),
    }
    ungauged["area_m2"] = float(sum(hru_area.get(h, 0.0) for h in ungauged["hrus"]))
    levels = sorted({z["level"] for z in zones})
    return {
        "schema": "openwq_wq_zones/1",
        "topology": info["path"],
        "species": list(species or []),
        "min_obs": int(min_obs),
        "period": list(period) if period else None,
        "n_reaches": len(topo["seg_ids"]),
        "n_hrus": len(topo["hru_ids"]),
        "zones": zones,
        "ungauged": ungauged,
        "levels": [[z["id"] for z in zones if z["level"] == lv] for lv in levels],
        "reach_zone": {str(k): v for k, v in reach_zone.items()},
        "hru_zone": {str(k): v for k, v in hru_zone.items()},
        "rejected": rejected,
    }


# ---------------------------------------------------------------------------
# Tables and inheritance
# ---------------------------------------------------------------------------
def write_zone_tables(zones: Dict[str, Any], out_dir: str) -> Dict[str, str]:
    """Write ``zones.json`` and the per-unit attribute tables; returns their paths."""
    os.makedirs(out_dir, exist_ok=True)
    paths = {"zones": os.path.join(out_dir, "zones.json"),
             "reaches": os.path.join(out_dir, "zone_reaches.csv"),
             "hrus": os.path.join(out_dir, "zone_hrus.csv"),
             "summary": os.path.join(out_dir, "zones_summary.csv")}
    with open(paths["zones"], "w") as f:
        json.dump(zones, f, indent=2)
    with open(paths["reaches"], "w") as f:
        f.write("id,%s\n" % ZONE_ATTRIBUTE)
        for k, v in zones["reach_zone"].items():
            f.write(f"{k},{v}\n")
    with open(paths["hrus"], "w") as f:
        f.write("id,%s\n" % ZONE_ATTRIBUTE)
        for k, v in zones["hru_zone"].items():
            f.write(f"{k},{v}\n")
    with open(paths["summary"], "w") as f:
        f.write("zone,station_reach,level,n_obs,obs_first,obs_last,n_reaches,area_km2,drainage_area_km2,upstream_zones,downstream_zone\n")
        for z in zones["zones"]:
            f.write("%d,%d,%d,%d,%s,%s,%d,%.1f,%.1f,%s,%s\n" % (
                z["id"], z["station_reach"], z["level"], z["n_obs"], z["obs_first"],
                z["obs_last"], z["n_reaches"], z["area_m2"] / 1e6,
                z["drainage_area_m2"] / 1e6, ";".join(map(str, z["upstream_zones"])),
                z["downstream_zone"]))
        u = zones["ungauged"]
        f.write("0,,,,,,%d,%.1f,,,\n" % (u["n_reaches"], u["area_m2"] / 1e6))
    return paths


def load_zones(work_dir: str) -> Optional[Dict[str, Any]]:
    p = os.path.join(work_dir, ZONES_DIRNAME, "zones.json")
    if os.path.isfile(p):
        with open(p) as f:
            return json.load(f)
    return None


def inherit_value(zones: Dict[str, Any], zone_values: Dict[int, float],
                  global_value: float, rule: str = "global") -> float:
    """Value for the ungauged zone 0: ``global``, ``upstream`` (area-weighted
    mean of the gauged zones that drain into it) or ``downstream_most`` (the
    zone of the last station, the one closest to the outlet)."""
    if rule == "upstream":
        num = den = 0.0
        for z in zones["zones"]:
            if z.get("downstream_zone") == 0 and z["id"] in zone_values:
                num += zone_values[z["id"]] * z["drainage_area_m2"]
                den += z["drainage_area_m2"]
        if den > 0:
            return num / den
    if rule == "downstream_most":
        last = max(zones["zones"], key=lambda z: z["drainage_area_m2"], default=None)
        if last and last["id"] in zone_values:
            return zone_values[last["id"]]
    return global_value
