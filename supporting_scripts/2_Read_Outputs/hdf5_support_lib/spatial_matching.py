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
Shared spatial-matching helpers for the OpenWQ supporting scripts.

These functions are the single source of truth for matching observation
stations to the model's spatial elements (river reaches for mizuRoute,
basin polygons / HRUs for SUMMA).  They are used by:

* ``2_Read_Outputs/hdf5_support_lib/Plot_h5_driver.py`` — config-side
  interactive viewer (post-run report).
* ``3_Calibration/calibration_lib/observation_data/prepare_calibration_observations.py``
  — bakes ``reach_id`` and ``is_primary`` into the calibration obs CSV.
* ``3_Calibration/calibration_lib/Gen_Calibration_Results_Report.py``
  — renders the same map visualisation on the calibration side.

Keeping a single implementation here means primary-vs-secondary semantics
(pouring-point obs vs. extra-inside-basin obs) and the basin-outlet logic
(river-network leaf nearest the polygon boundary) stay consistent across
the manual, automatic, and report-rendering paths.
"""

from __future__ import annotations

import math
from typing import Dict, Iterable, List, Optional, Set, Tuple


# ---------------------------------------------------------------------------
# Distance + geometry primitives
# ---------------------------------------------------------------------------

def haversine(lat1: float, lon1: float,
              lat2: float, lon2: float) -> float:
    """Great-circle distance in km between two ``(lat, lon)`` points."""
    R = 6371.0
    dlat = math.radians(lat2 - lat1)
    dlon = math.radians(lon2 - lon1)
    a = (math.sin(dlat / 2) ** 2 +
         math.cos(math.radians(lat1)) * math.cos(math.radians(lat2)) *
         math.sin(dlon / 2) ** 2)
    return R * 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))


def extract_coords(geom: dict) -> List[List[float]]:
    """Flatten every coordinate pair from a GeoJSON geometry."""
    gtype = geom.get('type', '')
    coords = geom.get('coordinates', [])
    if gtype == 'Point':
        return [coords]
    if gtype == 'LineString':
        return coords
    if gtype in ('MultiLineString', 'Polygon'):
        return [c for ring in coords for c in ring]
    if gtype == 'MultiPolygon':
        return [c for poly in coords for ring in poly for c in ring]
    return []


def polygon_rings(geom: dict) -> List[List[List[float]]]:
    """Return outer rings of a Polygon / MultiPolygon as ``[[lon, lat], ...]``.

    GeoJSON Polygon stores its outer ring at index 0 of ``coordinates`` and
    any interior holes at subsequent indices; we ignore holes here.
    """
    gtype = geom.get('type', '')
    coords = geom.get('coordinates', [])
    if gtype == 'Polygon':
        return [coords[0]] if coords else []
    if gtype == 'MultiPolygon':
        return [poly[0] for poly in coords if poly]
    return []


def point_in_polygon(lon: float, lat: float,
                     ring: List[List[float]]) -> bool:
    """Ray-casting point-in-polygon test.

    ``ring`` is a list of ``[lon, lat]`` vertices.  Returns ``False`` for
    degenerate rings of fewer than three vertices.
    """
    n = len(ring)
    if n < 3:
        return False
    inside = False
    j = n - 1
    for i in range(n):
        xi, yi = ring[i][0], ring[i][1]
        xj, yj = ring[j][0], ring[j][1]
        if ((yi > lat) != (yj > lat)) and \
           (lon < (xj - xi) * (lat - yi) / ((yj - yi) or 1e-30) + xi):
            inside = not inside
        j = i
    return inside


def polygon_centroid(rings: List[List[List[float]]]
                     ) -> Optional[Tuple[float, float]]:
    """Area-weighted centroid of a ``(Multi)Polygon`` outer rings.

    Falls back to the vertex mean if every ring is degenerate.  Returns
    ``(lon, lat)`` or ``None`` if no vertices exist.
    """
    sx = sy = sa = 0.0
    for ring in rings:
        n = len(ring)
        if n < 3:
            continue
        ring_closed = ring if ring[0] == ring[-1] else ring + [ring[0]]
        for i in range(len(ring_closed) - 1):
            x0, y0 = ring_closed[i][0], ring_closed[i][1]
            x1, y1 = ring_closed[i + 1][0], ring_closed[i + 1][1]
            a = x0 * y1 - x1 * y0
            sa += a
            sx += (x0 + x1) * a
            sy += (y0 + y1) * a
    if sa != 0:
        sa *= 0.5
        return (sx / (6 * sa), sy / (6 * sa))
    pts = [p for ring in rings for p in ring]
    if not pts:
        return None
    return (sum(p[0] for p in pts) / len(pts),
            sum(p[1] for p in pts) / len(pts))


# ---------------------------------------------------------------------------
# River-network endpoint graph (used to locate basin outlets)
# ---------------------------------------------------------------------------

def river_endpoint_degrees(river_geojson: Optional[dict]
                           ) -> Dict[Tuple[float, float], int]:
    """Return ``{(lon, lat)-rounded: degree}`` for the river-network graph.

    Degree = number of LineString segments that share this endpoint.  A
    *leaf* (degree==1) is topologically either a headwater or the basin
    outlet — never an inflow junction.
    """
    if not river_geojson or not river_geojson.get('features'):
        return {}
    deg: Dict[Tuple[float, float], int] = {}

    def _key(pt):
        return (round(float(pt[0]), 7), round(float(pt[1]), 7))

    for feat in river_geojson['features']:
        geom = feat.get('geometry') or {}
        gtype = geom.get('type')
        coords = geom.get('coordinates') or []
        lines = []
        if gtype == 'LineString' and len(coords) >= 2:
            lines = [coords]
        elif gtype == 'MultiLineString':
            lines = [ln for ln in coords if len(ln) >= 2]
        for line in lines:
            for pt in (line[0], line[-1]):
                k = _key(pt)
                deg[k] = deg.get(k, 0) + 1
    return deg


def find_basin_outlet(rings: List[List[List[float]]],
                      river_geojson: Optional[dict],
                      endpoint_deg:
                      Optional[Dict[Tuple[float, float], int]] = None
                      ) -> Optional[Tuple[float, float]]:
    """Return ``(lon, lat)`` of the basin's outlet via the river network.

    Strategy — every river *segment* with one endpoint inside the basin and
    one outside crosses the boundary.  The outside endpoint of that segment
    is a candidate outlet.  We then prefer endpoints that are LEAVES of the
    river graph (degree==1), i.e. terminal nodes that cannot be inflow
    junctions — they're almost always the basin's outlet rather than a
    transient crossing.  Ties broken by southernmost (rough downstream
    proxy for northern-hemisphere basins).

    If no crossing segments are found (lumped basin where the river is
    fully interior), fall back to the river endpoint inside the basin that
    is closest to the polygon boundary — still a better outlet estimate
    than the polygon centroid.
    """
    if not rings or not river_geojson:
        return None
    feats = river_geojson.get('features') or []
    if not feats:
        return None
    if endpoint_deg is None:
        endpoint_deg = river_endpoint_degrees(river_geojson)

    def _in_polygon(lon, lat):
        return any(point_in_polygon(lon, lat, r) for r in rings)

    def _ekey(pt):
        return (round(float(pt[0]), 7), round(float(pt[1]), 7))

    # 1) Crossing segments -> outside endpoint is the candidate outlet.
    crossings = []  # list of (degree, lat, lon)
    for feat in feats:
        geom = feat.get('geometry') or {}
        gtype = geom.get('type')
        coords = geom.get('coordinates') or []
        lines = []
        if gtype == 'LineString' and len(coords) >= 2:
            lines = [coords]
        elif gtype == 'MultiLineString':
            lines = [ln for ln in coords if len(ln) >= 2]
        for line in lines:
            sx, sy = float(line[0][0]), float(line[0][1])
            ex, ey = float(line[-1][0]), float(line[-1][1])
            si = _in_polygon(sx, sy)
            ei = _in_polygon(ex, ey)
            if si == ei:
                continue
            out_pt = (ex, ey) if si else (sx, sy)
            d = endpoint_deg.get(_ekey(out_pt), 0)
            crossings.append((d, out_pt[1], out_pt[0]))

    if crossings:
        leaves = [c for c in crossings if c[0] == 1]
        pool = leaves if leaves else crossings
        pool.sort(key=lambda c: (c[1], c[0]))
        return (pool[0][2], pool[0][1])

    # 2) No crossings — closest-interior-endpoint fallback.
    boundary_pts = [(v[0], v[1]) for r in rings for v in r]
    if not boundary_pts:
        return None
    candidates = []  # (dist_to_boundary_km, degree, lat, lon)
    for (lon, lat), degree in endpoint_deg.items():
        if not _in_polygon(lon, lat):
            continue
        min_d = min(haversine(lat, lon, v[1], v[0]) for v in boundary_pts)
        candidates.append((min_d, degree, lat, lon))
    if not candidates:
        return None
    leaves = [c for c in candidates if c[1] == 1]
    pool = leaves if leaves else candidates
    pool.sort(key=lambda c: (c[0], c[2]))
    return (pool[0][3], pool[0][2])


# ---------------------------------------------------------------------------
# Station ↔ feature matching (the public API)
# ---------------------------------------------------------------------------

def match_stations_to_features(
        station_locations: Dict[str, Tuple[float, float]],
        river_geojson: dict,
        mapping_key: str
) -> Dict[str, str]:
    """Match each station to its NEAREST river reach (mizuRoute case).

    Parameters
    ----------
    station_locations
        ``{station_id: (lat, lon)}``.
    river_geojson
        GeoJSON FeatureCollection for the river network.
    mapping_key
        Property key holding the reach identifier in each feature.

    Returns
    -------
    Dict[str, str]
        ``{station_id: feature_id}``.  Every station is matched (nearest
        wins) as long as the river network has at least one feature with
        coordinates.
    """
    matches: Dict[str, str] = {}
    if not station_locations or not river_geojson:
        return matches
    for stn_id, (slat, slon) in station_locations.items():
        best_fid: Optional[str] = None
        best_dist = float('inf')
        for feat in river_geojson.get('features', []):
            fid = str(feat['properties'].get(mapping_key, ''))
            for c in extract_coords(feat.get('geometry') or {}):
                d = haversine(slat, slon, c[1], c[0])
                if d < best_dist:
                    best_dist = d
                    best_fid = fid
        if best_fid is not None:
            matches[stn_id] = best_fid
    return matches


def match_stations_to_basins(
        station_locations: Dict[str, Tuple[float, float]],
        basin_geojson: dict,
        mapping_key: str,
        river_geojson: Optional[dict] = None,
        log: Optional[callable] = None,
        max_distance_km: Optional[float] = 5.0,
) -> Tuple[Dict[str, str], Set[str]]:
    """Match each station to a basin and designate the pouring-point obs.

    Matching is robust:

    1. Point-in-polygon — for stations actually inside a basin polygon.
    2. Nearest-basin-centroid fallback — for stations digitised just
       outside the polygon boundary (the common reason every observation
       ends up unmatched in practice when the user's basin shapefile is
       coarser than the gauge location).

    The pouring-point station of a basin is the station closest to the
    basin's RIVER-NETWORK OUTLET (see :func:`find_basin_outlet`).  When no
    river network is supplied we fall back to the polygon centroid — less
    accurate but the only sensible default in that case.

    Parameters
    ----------
    station_locations
        ``{station_id: (lat, lon)}``.
    basin_geojson
        GeoJSON FeatureCollection of basin / HRU polygons.
    mapping_key
        Property key holding the basin identifier in each feature.
    river_geojson
        Optional river-network FeatureCollection used to locate the
        outlet of each basin (highly recommended for SUMMA).
    log
        Optional ``print``-like callable used to surface step counts; if
        omitted nothing is printed.

    Returns
    -------
    Tuple[Dict[str, str], Set[str]]
        ``(station_to_basin, pouring_point_stations)``.  The set is
        always a subset of ``station_to_basin.keys()``.
    """
    _log = log or (lambda *a, **kw: None)

    if not station_locations or not basin_geojson:
        return {}, set()

    # Pre-cache per-basin geometry helpers.
    rings_by_fid: Dict[str, List[List[List[float]]]] = {}
    centroids: Dict[str, Tuple[float, float]] = {}
    for feat in basin_geojson.get('features', []):
        fid = str(feat['properties'].get(mapping_key, ''))
        if not fid:
            continue
        rings = polygon_rings(feat.get('geometry') or {})
        if not rings:
            continue
        rings_by_fid[fid] = rings
        c = polygon_centroid(rings)
        if c is not None:
            centroids[fid] = c

    # 1) Point-in-polygon pass.
    station_to_basin: Dict[str, str] = {}
    for fid, rings in rings_by_fid.items():
        for sid, (slat, slon) in station_locations.items():
            if sid in station_to_basin:
                continue
            for ring in rings:
                if point_in_polygon(slon, slat, ring):
                    station_to_basin[sid] = fid
                    break
    _n_pip = len(station_to_basin)
    _log(f"  [match] {_n_pip}/{len(station_locations)} stations matched "
         f"via point-in-polygon")

    # 2) Nearest-basin-centroid fallback — BOUNDED by max_distance_km.
    # Stations beyond the threshold stay unmatched (gray on the map),
    # which is the right behaviour for gauges that happen to be in the
    # bounding box but actually sit outside the modelled domain.  When
    # the caller passes ``max_distance_km=None`` we revert to the
    # legacy unbounded behaviour for backward compatibility.
    if _n_pip < len(station_locations) and centroids:
        _n_within = 0
        _n_out = 0
        for sid, (slat, slon) in station_locations.items():
            if sid in station_to_basin:
                continue
            best_fid, best_dist = None, float('inf')
            for fid, (clon, clat) in centroids.items():
                d = haversine(slat, slon, clat, clon)
                if d < best_dist:
                    best_dist = d
                    best_fid = fid
            if best_fid is not None and (
                    max_distance_km is None
                    or best_dist <= max_distance_km):
                station_to_basin[sid] = best_fid
                _n_within += 1
            else:
                _n_out += 1
        _msg = (f"  [match] {_n_within} stations matched via nearest-basin "
                f"fallback")
        if max_distance_km is not None:
            _msg += f" (within {max_distance_km:.1f} km)"
        if _n_out:
            _msg += f"; {_n_out} stations beyond threshold left unmatched"
        _log(_msg)

    # 3) Determine each basin's outlet via the river network.
    endpoint_deg = (river_endpoint_degrees(river_geojson)
                    if river_geojson else {})
    _log(f"  [match] River-network endpoints in graph: {len(endpoint_deg)}")
    basin_outlets: Dict[str, Tuple[float, float]] = {}
    for fid, rings in rings_by_fid.items():
        outlet = find_basin_outlet(rings, river_geojson, endpoint_deg)
        if outlet is None and fid in centroids:
            outlet = centroids[fid]
        if outlet is not None:
            basin_outlets[fid] = outlet
    _n_river_outlets = sum(1 for fid in basin_outlets
                           if fid in rings_by_fid
                           and basin_outlets[fid] != centroids.get(fid))
    _log(f"  [match] Basin outlets via river network: "
         f"{_n_river_outlets}/{len(rings_by_fid)} basins "
         f"(rest fell back to centroid)")

    # 4) Pouring-point designation: closest station to each basin's outlet.
    by_basin: Dict[str, List[str]] = {}
    for sid, fid in station_to_basin.items():
        by_basin.setdefault(fid, []).append(sid)
    pouring: Set[str] = set()
    for fid, sids in by_basin.items():
        if len(sids) == 1:
            pouring.add(sids[0])
            continue
        outlet = basin_outlets.get(fid)
        if outlet is None:
            pouring.add(sids[0])
            continue
        olon, olat = outlet
        best = min(sids, key=lambda s: haversine(
            station_locations[s][0], station_locations[s][1], olat, olon))
        pouring.add(best)

    return station_to_basin, pouring


# ---------------------------------------------------------------------------
# Host-aware top-level dispatcher
# ---------------------------------------------------------------------------

def match_stations(
        station_locations: Dict[str, Tuple[float, float]],
        hostmodel: str,
        river_geojson: Optional[dict] = None,
        basin_geojson: Optional[dict] = None,
        river_mapping_key: Optional[str] = None,
        basin_mapping_key: Optional[str] = None,
        log: Optional[callable] = None,
        max_distance_km: Optional[float] = 5.0,
) -> Tuple[Dict[str, str], Set[str]]:
    """One-stop entry point for "match my stations to the model elements".

    Dispatches to the right matcher based on ``hostmodel``:

    * ``"summa"`` -> :func:`match_stations_to_basins` using the basin
      GeoJSON + (optionally) the river network for outlet detection.
    * anything else (``"mizuroute"``, ``"mizuroute_cslm"``, etc.) ->
      :func:`match_stations_to_features` using the river network.

    ``max_distance_km`` bounds the matching distance: stations farther
    than this from any basin (SUMMA) or river segment (mizuRoute) stay
    unmatched and render gray on the map.  ``None`` disables the bound
    (legacy behaviour: nearest-feature wins regardless of distance).

    Returns ``(station_to_feature, primary_stations)``.  For the line-
    geometry hosts every matched station is treated as primary (set
    equals ``station_to_feature.keys()``).
    """
    h = (hostmodel or '').lower()
    _log = log or (lambda *a, **kw: None)
    if h == 'summa':
        if not basin_geojson:
            (log or print)(
                "  [match] hostmodel='summa' but no basin_geojson supplied; "
                "returning empty mapping."
            )
            return {}, set()
        return match_stations_to_basins(
            station_locations, basin_geojson,
            mapping_key=basin_mapping_key or 'HRU_ID',
            river_geojson=river_geojson, log=log,
            max_distance_km=max_distance_km)
    # mizuRoute family — bound the nearest-reach match by distance too.
    if not river_geojson:
        (log or print)(
            "  [match] hostmodel='%s' but no river_geojson supplied; "
            "returning empty mapping." % h
        )
        return {}, set()
    raw = match_stations_to_features(
        station_locations, river_geojson,
        mapping_key=river_mapping_key or 'SegId')
    if max_distance_km is None:
        return raw, set(raw.keys())
    # Re-compute the nearest distance for each matched station and drop
    # any beyond max_distance_km.  This walks the same coordinate set
    # the original matcher used so the threshold has the same units.
    feats = river_geojson.get('features') or []
    feat_coords: Dict[str, List[List[float]]] = {}
    for f in feats:
        fid = str(f['properties'].get(river_mapping_key or 'SegId', ''))
        if fid:
            feat_coords[fid] = list(extract_coords(f.get('geometry') or {}))
    bounded: Dict[str, str] = {}
    n_out = 0
    for sid, fid in raw.items():
        slat, slon = station_locations[sid]
        best_d = float('inf')
        for c in feat_coords.get(fid, []):
            d = haversine(slat, slon, c[1], c[0])
            if d < best_d:
                best_d = d
                if best_d <= max_distance_km:
                    break
        if best_d <= max_distance_km:
            bounded[sid] = fid
        else:
            n_out += 1
    if n_out:
        _log(f"  [match] {n_out} stations beyond {max_distance_km:.1f} km "
             f"from any reach — left unmatched")
    return bounded, set(bounded.keys())


__all__ = [
    'haversine', 'extract_coords', 'polygon_rings',
    'point_in_polygon', 'polygon_centroid',
    'river_endpoint_degrees', 'find_basin_outlet',
    'match_stations_to_features', 'match_stations_to_basins',
    'match_stations',
]
