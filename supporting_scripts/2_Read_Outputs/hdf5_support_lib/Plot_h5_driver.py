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
OpenWQ Interactive Time-Series Plotting Tool
Generates a self-contained HTML report with interactive Plotly.js charts.
Supports multiple chemical species (all plots in a single HTML file).
"""

import numpy as np
import pandas as pd
import netCDF4
import os
import json
import re
import datetime
import math
from collections import defaultdict


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def _sanitize_id(name):
    """Convert a species/extension name to a valid HTML element id."""
    return re.sub(r'[^a-zA-Z0-9_]', '_', str(name))


def _haversine(lat1, lon1, lat2, lon2):
    """Haversine distance in km between two lat/lon points."""
    R = 6371
    dlat = math.radians(lat2 - lat1)
    dlon = math.radians(lon2 - lon1)
    a = (math.sin(dlat / 2) ** 2 +
         math.cos(math.radians(lat1)) * math.cos(math.radians(lat2)) *
         math.sin(dlon / 2) ** 2)
    return R * 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))


def _extract_coords(geom):
    """Flatten all coordinate pairs from a GeoJSON geometry dict."""
    gtype = geom.get('type', '')
    coords = geom.get('coordinates', [])
    if gtype == 'Point':
        return [coords]
    elif gtype == 'LineString':
        return coords
    elif gtype in ('MultiLineString', 'Polygon'):
        return [c for ring in coords for c in ring]
    elif gtype == 'MultiPolygon':
        return [c for poly in coords for ring in poly for c in ring]
    return []


def _match_stations_to_features(station_locations, river_geojson, mapping_key):
    """Match observation stations to their nearest river feature.

    Parameters
    ----------
    station_locations : dict
        {station_id: (lat, lon)}
    river_geojson : dict
        GeoJSON FeatureCollection for the river network.
    mapping_key : str
        Property key for the feature ID in GeoJSON features.

    Returns
    -------
    dict
        {station_id: feature_id}
    """
    matches = {}
    for stn_id, (slat, slon) in station_locations.items():
        best_fid = None
        best_dist = float('inf')
        for feat in river_geojson['features']:
            fid = str(feat['properties'].get(mapping_key, ''))
            for c in _extract_coords(feat['geometry']):
                d = _haversine(slat, slon, c[1], c[0])
                if d < best_dist:
                    best_dist = d
                    best_fid = fid
        if best_fid is not None:
            matches[stn_id] = best_fid
    return matches


def _load_observation_data(obs_dir=None, obs_csv=None):
    """Load observation data from a GRQA clipped directory or user CSV.

    Returns
    -------
    obs_by_species : dict or None
        {species_name: list of {station_id, datetime, value}}
    station_locations : dict or None
        {station_id: (lat, lon)}
    """
    obs_by_species = {}
    station_locations = {}

    if obs_dir and os.path.isdir(obs_dir):
        obs_csv_path = os.path.join(obs_dir, 'grqa_clipped_observations.csv')
        stn_csv_path = os.path.join(obs_dir, 'grqa_clipped_stations.csv')
        if not os.path.isfile(obs_csv_path):
            print("  No GRQA observations file found in obs dir.")
            return None, None
        try:
            df = pd.read_csv(obs_csv_path)
            # Station locations
            if os.path.isfile(stn_csv_path):
                stn_df = pd.read_csv(stn_csv_path)
                for _, row in stn_df.iterrows():
                    sid = str(row.get('site_id', ''))
                    lat = float(row.get('lat_wgs84', 0))
                    lon = float(row.get('lon_wgs84', 0))
                    if sid:
                        station_locations[sid] = (lat, lon)
            else:
                for _, row in df.drop_duplicates('site_id').iterrows():
                    sid = str(row.get('site_id', ''))
                    lat = float(row.get('lat_wgs84', 0))
                    lon = float(row.get('lon_wgs84', 0))
                    if sid:
                        station_locations[sid] = (lat, lon)
            # Group by model_species
            if 'model_species' in df.columns:
                for species, grp in df.groupby('model_species'):
                    records = []
                    for _, row in grp.iterrows():
                        try:
                            records.append({
                                'station_id': str(row['site_id']),
                                'datetime': str(row.get('obs_date', '')),
                                'value': float(row.get('obs_value', 0)),
                            })
                        except (ValueError, TypeError):
                            continue
                    if records:
                        obs_by_species[species] = records
            print(f"  ✓ Loaded GRQA observations: {len(df)} records, "
                  f"{len(obs_by_species)} species, "
                  f"{len(station_locations)} stations")
        except Exception as e:
            print(f"  WARNING: Failed to load GRQA observations: {e}")
            return None, None

    elif obs_csv and os.path.isfile(obs_csv):
        try:
            df = pd.read_csv(obs_csv)
            required = ['station_id', 'lat', 'lon', 'parameter', 'value']
            missing = [c for c in required if c not in df.columns]
            if missing:
                print(f"  WARNING: User CSV missing columns: {missing}")
                return None, None
            for _, row in df.drop_duplicates('station_id').iterrows():
                sid = str(row['station_id'])
                station_locations[sid] = (float(row['lat']), float(row['lon']))
            date_cols = ['year', 'month', 'day']
            if all(c in df.columns for c in date_cols):
                df['_datetime'] = pd.to_datetime(
                    df[date_cols].assign(hour=0,
                                         minute=df.get('minute', 0)),
                    errors='coerce')
            else:
                df['_datetime'] = pd.NaT
            for species, grp in df.groupby('parameter'):
                records = []
                for _, row in grp.iterrows():
                    try:
                        dt = row['_datetime']
                        dt_str = (dt.strftime('%Y-%m-%d')
                                  if pd.notna(dt) else '')
                        records.append({
                            'station_id': str(row['station_id']),
                            'datetime': dt_str,
                            'value': float(row['value']),
                        })
                    except (ValueError, TypeError):
                        continue
                if records:
                    obs_by_species[species] = records
            print(f"  ✓ Loaded user observations: {len(df)} records, "
                  f"{len(obs_by_species)} species, "
                  f"{len(station_locations)} stations")
        except Exception as e:
            print(f"  WARNING: Failed to load user observations: {e}")
            return None, None
    else:
        return None, None

    return obs_by_species, station_locations


def _build_traces(feature_data, n_visible=10, feature_label='Feature'):
    """Build Plotly trace dicts from a {fid: pd.Series} dictionary.

    Parameters
    ----------
    feature_data : dict
        {feature_id: pd.Series} — each Series has a time index.
    n_visible : int
        When more than 20 features, only the first *n_visible* are shown
        by default; the rest are set to ``'legendonly'``.
    feature_label : str
        Label prefix for legend entries (e.g. 'SubId', 'SegId').

    Returns
    -------
    list[dict]
        Plotly trace dicts ready for ``json.dumps()``.
    """
    traces = []
    n_features = len(feature_data)

    for idx, (fid, series) in enumerate(feature_data.items()):
        # Convert time index to ISO strings for Plotly
        if isinstance(series.index, pd.DatetimeIndex):
            x_vals = series.index.strftime('%Y-%m-%dT%H:%M:%S').tolist()
        else:
            x_vals = series.index.tolist()

        # Convert NaN to None (JSON null) so Plotly renders gaps
        y_vals = [None if (isinstance(v, float) and np.isnan(v)) else float(v)
                  for v in series.values]

        trace = {
            'x': x_vals,
            'y': y_vals,
            'mode': 'lines',
            'name': f'{feature_label} {fid}',
            'hovertemplate': f'{feature_label} {fid}<br>%{{x}}<br>%{{y:.4g}}<extra></extra>',
        }

        # For >20 features, hide extras behind legend toggle
        if n_features > 20 and idx >= n_visible:
            trace['visible'] = 'legendonly'

        traces.append(trace)

    return traces


# ── Debug-mode helpers ──────────────────────────────────────────────────
# Colorway matching Plotly layout — explicit assignment lets debug traces
# reuse the same colour as their corresponding main trace.
_COLORWAY = ['#0066cc', '#00a86b', '#ff6b35', '#004499', '#34d399', '#fb923c',
             '#667eea', '#764ba2', '#e63946', '#2ec4b6', '#e9c46a', '#264653']

# Metadata for each debug file extension
_DEBUG_EXT_META = {
    'd_output_dt_chemistry': {'label': 'dC/dt chem',   'dash': 'dash'},
    'd_output_dt_transport': {'label': 'dC/dt transp', 'dash': 'dot'},
    'd_output_ss':           {'label': 'SS',            'dash': 'dashdot'},
    'd_output_ewf':          {'label': 'EWF',           'dash': 'longdash'},
    'd_output_ic':           {'label': 'IC',            'dash': 'longdashdot'},
}


def _hru_from_fid(fid):
    """Extract the HRU/reach ID from a display fid like ``"1 (z1)"`` → ``"1"``."""
    import re
    m = re.match(r'^(.+?) \(z\d+\)$', str(fid))
    return m.group(1) if m else str(fid)


def _build_global_color_map(all_fids):
    """Build a consistent HRU→colour mapping from all fids across all plots.

    Colours are assigned by sorted HRU ID so they stay stable regardless
    of compartment or data type.
    """
    hru_set = set()
    for fid in all_fids:
        hru_set.add(_hru_from_fid(fid))
    hru_sorted = sorted(hru_set)
    hru_color = {hru: _COLORWAY[i % len(_COLORWAY)]
                 for i, hru in enumerate(hru_sorted)}
    # Expand to full fids (with layer suffixes)
    color_map = {}
    for fid in all_fids:
        hru = _hru_from_fid(fid)
        color_map[fid] = hru_color[hru]
    return color_map


def _build_traces_with_colors(feature_data, n_visible=10, feature_label='Feature',
                              global_color_map=None):
    """Like ``_build_traces`` but assigns explicit colours from ``_COLORWAY``.

    Parameters
    ----------
    global_color_map : dict or None
        Pre-built fid→colour mapping.  When provided, colours are taken
        from this map so they stay consistent across all compartments.

    Returns
    -------
    (list[dict], dict)
        (traces, color_map) where *color_map* maps ``fid → hex colour``.
    """
    traces = []
    color_map = {}
    n_features = len(feature_data)

    for idx, (fid, series) in enumerate(feature_data.items()):
        if global_color_map and fid in global_color_map:
            color = global_color_map[fid]
        else:
            color = _COLORWAY[idx % len(_COLORWAY)]
        color_map[fid] = color

        if isinstance(series.index, pd.DatetimeIndex):
            x_vals = series.index.strftime('%Y-%m-%dT%H:%M:%S').tolist()
        else:
            x_vals = series.index.tolist()

        y_vals = [None if (isinstance(v, float) and np.isnan(v)) else float(v)
                  for v in series.values]

        trace = {
            'x': x_vals,
            'y': y_vals,
            'mode': 'lines',
            'name': f'{feature_label} {fid}',
            'hovertemplate': (f'{feature_label} {fid}<br>'
                              f'%{{x}}<br>%{{y:.4g}}<extra></extra>'),
            'line': {'color': color},
        }

        if n_features > 20 and idx >= n_visible:
            trace['visible'] = 'legendonly'

        traces.append(trace)

    return traces, color_map


def _build_debug_traces(feature_data, color_map, ext_label, dash_style,
                        feature_label='Feature'):
    """Build traces for a single debug extension.

    Each trace gets its own legend entry showing the feature ID and the
    derivative label, e.g. ``HRU 1 (z1) dC/dt chem``.  Traces share
    the colour of the matching main-trace but use a different dash pattern.
    They start visible and can be toggled with the Debug button.
    """
    traces = []

    for fid, series in feature_data.items():
        color = color_map.get(fid, '#888')

        if isinstance(series.index, pd.DatetimeIndex):
            x_vals = series.index.strftime('%Y-%m-%dT%H:%M:%S').tolist()
        else:
            x_vals = series.index.tolist()

        y_vals = [None if (isinstance(v, float) and np.isnan(v)) else float(v)
                  for v in series.values]

        trace = {
            'x': x_vals,
            'y': y_vals,
            'mode': 'lines',
            'name': f'{feature_label} {fid} {ext_label}',
            'hovertemplate': (f'{feature_label} {fid} {ext_label}<br>'
                              f'%{{x}}<br>%{{y:.4g}}<extra></extra>'),
            'line': {'color': color, 'dash': dash_style, 'width': 1.5},
            'legendgroup': f'debug_{ext_label}_{fid}',
            'showlegend': True,
            'visible': True,
        }

        traces.append(trace)

    return traces


def _build_html(plots, what2map, hostmodel, river_geojson=None,
                map_center=None, map_bounds=None, mapping_key='SegId',
                observation_data=None, feature_label=None,
                map_geom_type='line', station_locations=None,
                station_to_feature=None, separator=' | ',
                basin_geojson=None, river_line_geojson=None):
    """Build a self-contained HTML string with interactive Plotly.js charts.

    Parameters
    ----------
    plots : list[dict]
        Each dict has keys: id, title, traces, ylabel, xtype, caption.
    what2map : str
        'hostmodel' or 'openwq'.
    hostmodel : str
        Host model name (e.g. 'mizuroute').
    river_geojson : dict or None
        GeoJSON FeatureCollection for the PRIMARY interactive layer
        (the one whose features match plot data and are clickable).
    map_center : list or None
        [lat, lng] for the initial map center.
    map_bounds : tuple or None
        (minx, miny, maxx, maxy) for auto-fitting the map.
    mapping_key : str
        Property name in the GeoJSON features that matches feature IDs in plots.
    feature_label : str or None
        Display label for features in legends, tooltips, etc.
        If None, defaults to mapping_key.
    map_geom_type : str
        'line' for river segments, 'polygon' for basin/HRU areas.
        Describes the geometry type of the PRIMARY interactive layer.
    basin_geojson : dict or None
        GeoJSON FeatureCollection for basin/catchment polygons.
        Always rendered at the back. Interactive if map_geom_type=='polygon'.
    river_line_geojson : dict or None
        GeoJSON FeatureCollection for river network lines.
        Rendered on top of basins. Interactive if map_geom_type=='line'.

    Returns
    -------
    str
        Complete HTML document as a string.
    """
    if feature_label is None:
        feature_label = mapping_key
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    n_plots = len(plots)

    # Build a stable fid → color mapping so map and plots use identical colors.
    _colorway = ['#0066cc','#00a86b','#ff6b35','#004499','#34d399','#fb923c',
                 '#667eea','#764ba2','#e63946','#2ec4b6','#e9c46a','#264653']
    _label_prefix = f'{feature_label} '
    _all_fids_set = set()
    for p in plots:
        for t in p['traces']:
            # Only collect feature IDs from main traces (not debug)
            if t.get('legendgroup', '').startswith('debug_'):
                continue
            fid = t['name'].replace(_label_prefix, '')
            _all_fids_set.add(fid)
    _all_fids_sorted = sorted(_all_fids_set)

    # Detect layers: fids like "1 (z1)", "1 (z2)" have a " (z<n>)" suffix
    _all_layers = set()
    _all_hru_ids = set()
    _LAYER_RE = re.compile(r'^(.+?) \((z\d+)\)$')
    for fid in _all_fids_sorted:
        m = _LAYER_RE.match(fid)
        if m:
            _all_layers.add(m.group(2))
            _all_hru_ids.add(m.group(1))
        else:
            _all_hru_ids.add(fid)
    _all_layers_sorted = sorted(_all_layers, key=lambda x: int(x[1:]) if x[1:].isdigit() else x)
    _has_layers = len(_all_layers) > 0

    # Assign colours based on HRU ID (not layer) so same HRU gets same
    # colour across all layers
    _hru_ids_sorted = sorted(_all_hru_ids)
    _hru_color = {hid: _colorway[i % len(_colorway)]
                  for i, hid in enumerate(_hru_ids_sorted)}

    # Build the full fid→color mapping.  For layered fids like "1 (z1)",
    # use the HRU ID colour so all layers of the same HRU share the colour.
    _fid_color = {}
    for fid in _all_fids_sorted:
        m = _LAYER_RE.match(fid)
        if m:
            _fid_color[fid] = _hru_color.get(m.group(1), '#888')
        else:
            _fid_color[fid] = _hru_color.get(fid, '#888')

    # Build a separate color mapping for map GeoJSON features.
    # When plot IDs and map IDs are in the same space (mizuRoute),
    # _fid_color already covers them.  When they differ (SUMMA: HDF5
    # uses sequential 1,2,3 but shapefile uses GRU_ID 740457190 etc.),
    # we build an independent map color dict from the GeoJSON features.
    _map_fid_color = {}
    if river_geojson and river_geojson.get('features'):
        _map_fid_set = set()
        for feat in river_geojson['features']:
            mfid = str(feat['properties'].get(mapping_key, ''))
            if mfid:
                _map_fid_set.add(mfid)
        _map_fids_sorted = sorted(_map_fid_set)
        # Check if map IDs overlap with plot IDs — if so, reuse plot colors
        if _map_fid_set & (_all_fids_set | _all_hru_ids):
            # IDs match (mizuRoute case): reuse _fid_color + _hru_color
            for mfid in _map_fids_sorted:
                _map_fid_color[mfid] = (_fid_color.get(mfid) or
                                        _hru_color.get(mfid) or '#888')
        else:
            # IDs don't match (SUMMA case): assign independent colors
            for i, mfid in enumerate(_map_fids_sorted):
                _map_fid_color[mfid] = _colorway[i % len(_colorway)]

    # Stamp explicit line colors onto every trace so Plotly matches the map.
    for p in plots:
        for t in p['traces']:
            if t.get('legendgroup', '').startswith('debug_'):
                # Debug trace: extract fid from hovertemplate (format:
                # "{ext}: {feature_label} {fid}<br>...")
                _ht = t.get('hovertemplate', '')
                _after = _ht.split(_label_prefix)[-1] if _label_prefix in _ht else ''
                _dfid = _after.split('<br>')[0] if '<br>' in _after else ''
                if _dfid and _dfid in _fid_color:
                    _existing = t.get('line', {})
                    _existing['color'] = _fid_color[_dfid]
                    t['line'] = _existing
            else:
                fid = t['name'].replace(_label_prefix, '')
                t['line'] = {'color': _fid_color.get(fid, '#888')}

    # --- Merge observation traces (if any) into each plot ---
    _obs_meta = {}   # {div_id: {trace_index_str: feature_id}}
    _has_obs = set()  # div_ids that have observation traces
    if observation_data:
        for p in plots:
            plot_obs = observation_data.get(p['id'], [])
            if not plot_obs:
                continue
            div_id = f'{p["id"]}_div'
            p_meta = {}
            base_idx = len(p['traces'])
            _species_name = p.get('species', '')
            for i, od in enumerate(plot_obs):
                fid = str(od['feature_id'])
                color = _fid_color.get(fid, '#888')
                trace = {
                    'x': od['x'],
                    'y': od['y'],
                    'mode': 'markers',
                    'name': f'Obs: {od["station_id"]}',
                    'marker': {
                        'color': color,
                        'size': 5,
                        'symbol': 'circle',
                        'line': {'width': 0.5, 'color': '#333'},
                    },
                    'hovertemplate': (
                        f'Station {od["station_id"]}<br>'
                        f'Species: {_species_name}<br>'
                        f'%{{x}}<br>%{{y:.4g}}'
                        f'<extra>\u2192 {feature_label} {fid}</extra>'
                    ),
                    'showlegend': False,
                }
                p['traces'].append(trace)
                p_meta[str(base_idx + i)] = fid
            if p_meta:
                _obs_meta[div_id] = p_meta
                _has_obs.add(div_id)

    # --- Compute summary stats for header KPI boxes ---
    _species_list = []
    _seen_species = set()
    for p in plots:
        s = p.get('species', p.get('title', ''))
        if s and s not in _seen_species:
            _seen_species.add(s)
            _species_list.append(s)
    _n_features = len(_all_fids_sorted)

    # Time range from traces
    _date_strs = []
    for p in plots:
        for t in p['traces']:
            if t.get('x') and t.get('mode') == 'lines':
                _date_strs.extend([str(x) for x in t['x'] if x])
    if _date_strs:
        _date_strs_sorted = sorted(_date_strs)
        _time_range_str = (f"{_date_strs_sorted[0][:10]} &rarr; "
                           f"{_date_strs_sorted[-1][:10]}")
    else:
        _time_range_str = "N/A"

    # Observation count
    _n_obs_stations = 0
    _n_obs_records = 0
    if observation_data:
        _obs_stn_set = set()
        for plot_obs_list in observation_data.values():
            for od in plot_obs_list:
                _obs_stn_set.add(od.get('station_id'))
                _n_obs_records += len(od.get('y', []))
        _n_obs_stations = len(_obs_stn_set)

    H = []

    # --- HEAD ---
    H.append("""<!DOCTYPE html>
<html lang="en" data-theme="dark">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>OpenWQ Results</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link href="https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&family=JetBrains+Mono:wght@400;500&display=swap" rel="stylesheet">
<script src="https://cdn.plot.ly/plotly-basic-2.35.2.min.js" defer></script>
<link rel="stylesheet" href="https://unpkg.com/leaflet@1.9.4/dist/leaflet.css"/>
<script src="https://unpkg.com/leaflet@1.9.4/dist/leaflet.js" defer></script>
<style>
:root,[data-theme="light"]{
  --primary:#0066cc;--primary-dark:#004499;
  --secondary:#00a86b;--accent:#ff6b35;
  --bg:#f0f2f5;--surface:#fff;--text:#2d3436;--text2:#636e72;
  --border:#e2e8f0;--border2:#d0d4da;
  --shadow:0 4px 20px rgba(0,0,0,.06);--shadow-lg:0 10px 40px rgba(0,0,0,.08);
  --glass:rgba(255,255,255,.85);
  --plot-bg:#fff;--plot-grid:#f0f0f0;--plot-font:#333;
}
[data-theme="dark"]{
  --primary:#4d9ee8;--primary-dark:#3a8bd4;
  --secondary:#34d399;--accent:#fb923c;
  --bg:#0f1117;--surface:#1a1b2e;--text:#e2e8f0;--text2:#a0aec0;
  --border:#2a2b3d;--border2:#3a3b4d;
  --shadow:0 4px 20px rgba(0,0,0,.3);--shadow-lg:0 10px 40px rgba(0,0,0,.4);
  --glass:rgba(15,17,23,.92);
  --plot-bg:#1e1f2e;--plot-grid:#2d2e42;--plot-font:#ccc;
}
*{margin:0;padding:0;box-sizing:border-box}
body{font-family:'Inter',-apple-system,BlinkMacSystemFont,sans-serif;background:var(--bg);
  color:var(--text);line-height:1.65}
html,.main{scroll-behavior:smooth}
a{color:var(--primary);text-decoration:none}
.layout{display:flex;height:100vh;overflow:hidden}
.sidebar{width:220px;min-width:120px;max-width:50vw;height:100vh;background:var(--glass);
  backdrop-filter:blur(12px);border-right:none;padding:1.2rem .8rem;
  display:flex;flex-direction:column;z-index:10;overflow-y:auto;flex-shrink:0;
  transition:background .3s;position:relative}
.sidebar-handle{position:absolute;top:0;right:0;width:5px;height:100%;
  cursor:col-resize;background:var(--border);transition:background .15s;z-index:11}
.sidebar-handle:hover,.sidebar-handle.dragging{background:var(--primary);width:5px}
.sidebar .logo{font-weight:700;font-size:1.1rem;margin-bottom:1.2rem;
  padding-left:.7rem;letter-spacing:-.3px;color:var(--text)}
.sidebar .logo span{color:var(--primary)}
.sidebar nav{display:flex;flex-direction:column;gap:2px}
.sidebar nav a{display:block;padding:.4rem .7rem;border-radius:6px;font-size:.82rem;
  color:var(--text2);transition:all .15s;white-space:nowrap;overflow:hidden;text-overflow:ellipsis}
.sidebar nav a:hover{background:rgba(0,102,204,.06);color:var(--primary)}
.sidebar nav a.active{background:rgba(0,102,204,.1);color:var(--primary);font-weight:600}
.theme-toggle{margin-top:auto;padding-top:1rem;border-top:1px solid var(--border)}
.theme-btn{display:flex;align-items:center;gap:.5rem;width:100%;padding:.45rem .7rem;
  border:1px solid var(--border);border-radius:8px;background:var(--surface);
  color:var(--text2);font-size:.78rem;cursor:pointer;transition:all .15s;font-family:inherit}
.theme-btn:hover{border-color:var(--primary);color:var(--primary)}
.theme-btn svg{width:16px;height:16px;flex-shrink:0}
.main{flex:1;min-width:0;overflow-y:auto;overflow-x:hidden;height:100vh}
.top-fixed{position:sticky;top:0;z-index:20;flex-shrink:0}
.header{background:linear-gradient(135deg,#1a1a2e 0%,#16213e 50%,#0f3460 100%);
  color:#fff;padding:1.2rem 2.5rem 1rem;position:relative;overflow:hidden}
.header::before{content:'';position:absolute;inset:0;
  background:radial-gradient(circle at 80% 20%,rgba(255,255,255,.08) 0%,transparent 60%)}
.header h1{font-size:2rem;font-weight:700;letter-spacing:-.5px;position:relative}
.header h1 span{color:var(--primary)}
.header .subtitle{opacity:.8;font-size:.92rem;margin-top:.4rem;position:relative}
.header .meta{display:flex;gap:1.2rem;margin-top:1rem;font-size:.78rem;opacity:.7;
  flex-wrap:wrap;position:relative}
.header .meta span{background:rgba(255,255,255,.1);padding:.2rem .6rem;border-radius:20px}
.container{max-width:1100px;margin:0 auto;padding:2rem 2rem 3rem}
.section{margin-bottom:2.5rem}
.section h2{font-size:1.2rem;font-weight:700;color:var(--text);margin-bottom:1rem;
  display:flex;align-items:center;gap:.5rem;letter-spacing:-.3px}
.section h2::before{content:'';width:5px;height:1.2em;border-radius:3px;
  background:linear-gradient(180deg,var(--primary),var(--primary-dark));flex-shrink:0}
.card{background:var(--surface);border:1px solid var(--border);border-radius:16px;
  padding:1.5rem;box-shadow:var(--shadow);transition:transform .3s,box-shadow .3s}
.card:hover{box-shadow:var(--shadow-lg);border-color:var(--border2)}
.plot-caption{font-size:.82rem;color:var(--text2);margin-top:.5rem;text-align:right}
.plot-toolbar{display:flex;align-items:center;gap:6px;padding:6px 10px;margin-bottom:4px;
  font-size:.8rem;flex-wrap:wrap}
.plot-toolbar label{color:var(--text2);font-weight:600}
.tb-btn{padding:3px 10px;border:1px solid var(--border);border-radius:4px;
  background:var(--surface);color:var(--text);cursor:pointer;font-size:.78rem;
  transition:background .15s,color .15s}
.tb-btn:hover{background:var(--primary);color:#fff}
.tb-btn.active{background:var(--primary);color:#fff}
.tb-input{width:80px;padding:2px 6px;border:1px solid var(--border);border-radius:4px;
  font-size:.78rem;background:var(--surface);color:var(--text)}
.tb-date{width:120px}
.tb-sep{width:1px;height:18px;background:var(--border);margin:0 4px}
.layer-bar{background:var(--surface);
  border-bottom:1px solid var(--border);padding:6px 12px;display:flex;
  align-items:center;gap:6px;font-size:.82rem;box-shadow:0 2px 8px rgba(0,0,0,.08)}
.layer-bar label{font-weight:700;margin-right:4px;color:var(--text)}
#mapContainer{width:100%;height:420px;border-radius:12px;z-index:1}
.map-legend{background:var(--glass);
  backdrop-filter:blur(8px);border:1px solid var(--border);border-radius:8px;
  padding:.5rem .7rem;font-size:.72rem;overflow-y:auto;
  color:var(--text);min-width:90px}
.map-legend .ml-title{font-weight:600;margin-bottom:.3rem;font-size:.8rem}
.map-legend .ml-item{display:flex;align-items:center;gap:.4rem;padding:2px 4px;
  border-radius:4px;cursor:pointer;transition:background .15s}
.map-legend .ml-item:hover{background:rgba(0,102,204,.08)}
.map-legend .ml-item.active{background:rgba(0,102,204,.18);font-weight:600;border-left:3px solid var(--primary,#0066cc);padding-left:3px}
.map-legend.has-selection .ml-item:not(.active){opacity:.35}
.map-legend .ml-swatch{width:14px;height:4px;border-radius:2px;flex-shrink:0}
.map-legend .ml-toggles{display:flex;flex-direction:column;gap:2px;margin-bottom:2px}
.map-legend .ml-toggle{display:flex;align-items:center;gap:4px;font-size:.72rem;cursor:pointer}
.map-legend .ml-toggle input{margin:0;cursor:pointer}
.map-legend .ml-toggle label{cursor:pointer;user-select:none}
</style>
</head>
<body>""")

    # --- SIDEBAR ---
    H.append('<div class="layout">')
    H.append('<aside class="sidebar">')
    H.append('<div class="logo">Open<span>WQ</span> Results</div>')
    H.append('<nav>')
    _map_title = 'HRU / Basin Map' if map_geom_type == 'polygon' else 'River Network Map'
    if river_geojson:
        H.append(f'<a href="#rivermap">{_map_title}</a>')
    _nav_prev_comp = None
    for p in plots:
        _nc = p.get('compartment', '')
        if _nc and _nc != _nav_prev_comp:
            _cid = f'comp_{re.sub(r"[^a-zA-Z0-9]", "_", _nc)}'
            H.append(f'<a href="#{_cid}" style="font-weight:700;'
                     f'margin-top:.5rem;color:var(--text)">{_nc}</a>')
            _nav_prev_comp = _nc
        _nav_label = p.get('nav_title', p['title'])
        H.append(f'<a href="#{p["id"]}" style="padding-left:1.2rem">'
                 f'{_nav_label}</a>')
    H.append('</nav>')
    H.append("""<div class="theme-toggle">
<button class="theme-btn" onclick="toggleTheme()">
<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"
  stroke-linecap="round" stroke-linejoin="round">
<circle cx="12" cy="12" r="5"/><line x1="12" y1="1" x2="12" y2="3"/>
<line x1="12" y1="21" x2="12" y2="23"/><line x1="4.22" y1="4.22" x2="5.64" y2="5.64"/>
<line x1="18.36" y1="18.36" x2="19.78" y2="19.78"/><line x1="1" y1="12" x2="3" y2="12"/>
<line x1="21" y1="12" x2="23" y2="12"/><line x1="4.22" y1="19.78" x2="5.64" y2="18.36"/>
<line x1="18.36" y1="5.64" x2="19.78" y2="4.22"/>
</svg>
<span>Toggle Theme</span>
</button>
</div>""")
    H.append('<div class="sidebar-handle" id="sidebarHandle"></div>')
    H.append('</aside>')

    # --- MAIN ---
    H.append('<div class="main">')

    # --- FIXED TOP BAR (header + layer selector, always pinned) ---
    H.append('<div class="top-fixed">')

    # Header
    mode_label = what2map.upper() if what2map else 'N/A'
    host_label = (hostmodel or 'N/A').upper()
    _obs_label = (f'{_n_obs_stations} station(s), {_n_obs_records} records'
                  if _n_obs_stations else 'None')
    _species_str = ', '.join(_species_list) if _species_list else 'N/A'

    H.append(f"""<div class="header">
<h1>Open<span>WQ</span> &mdash; Results</h1>
<p class="subtitle">{mode_label} mode &mdash; {n_plots} time-series plot(s)</p>
<div class="meta">
<span>Host Model: {host_label}</span>
<span>Generated: {now}</span>
<span>Species: {_species_str}</span>
<span>{feature_label}s: {_n_features}</span>
<span>Period: {_time_range_str}</span>
<span>Observations: {_obs_label}</span>
</div>
</div>""")

    # Layer selector removed from top bar — now per-plot (see plot toolbar below)

    H.append('</div>')  # close top-fixed

    H.append('<div class="container">')

    # --- RIVER NETWORK MAP ---
    if river_geojson:
        _feat_name = feature_label or mapping_key or 'feature'
        _map_click_hint = f'Click a {_feat_name} to select/deselect it. Selected features are isolated in all plots below.'
        H.append(f"""<div class="section" id="rivermap">
<h2>{_map_title}</h2>
<div class="card" style="padding:0;overflow:hidden;position:relative">
<div id="mapContainer"></div>
</div>
<p class="plot-caption" style="margin-top:.4rem;font-size:.78rem;color:var(--text2)">
{_map_click_hint}</p>
</div>""")

    # --- PLOT SECTIONS (grouped by compartment) ---
    _prev_comp = None
    for p in plots:
        _cur_comp = p.get('compartment', '')
        if _cur_comp and _cur_comp != _prev_comp:
            if _prev_comp is not None:
                H.append('</div>')  # close previous compartment group
            _comp_id = f'comp_{re.sub(r"[^a-zA-Z0-9]", "_", _cur_comp)}'
            H.append(f'<div class="section" id="{_comp_id}">'
                     f'<h2 style="font-size:1.3rem;margin-bottom:.5rem">'
                     f'{_cur_comp}</h2></div>')
            H.append(f'<div style="margin-left:.5rem">')
            _prev_comp = _cur_comp
        div_id = f'{p["id"]}_div'
        # Build second row: "Active:" + Concentrations toggle + per-derivative toggles
        _trace_btns = (
            f'  <button class="tb-btn active" '
            f'data-deriv="conc" data-plot="{div_id}" '
            f'onclick="_owqToggleDeriv(this)">'
            f'Concentrations \u2713</button>'
        )
        for _dlbl in p.get('debug_labels', []):
            _btn_id = f'{div_id}_dbg_{_dlbl.replace("/","").replace(" ","_")}'
            _trace_btns += (
                f'  <button class="tb-btn active" id="{_btn_id}" '
                f'data-deriv="{_dlbl}" data-plot="{div_id}" '
                f'onclick="_owqToggleDeriv(this)">'
                f'{_dlbl} \u2713</button>'
            )
        # Observation toggle button (only shown when obs data exists for this plot)
        if div_id in _has_obs:
            _trace_btns += (
                f'  <span class="tb-sep"></span>'
                f'  <button class="tb-btn active" '
                f'data-obs="1" data-plot="{div_id}" '
                f'onclick="_owqToggleObs(this)">'
                f'Observations \u2713</button>'
            )
        # Layer range selector (per-plot, only when THIS plot has layers)
        _layer_row = ''
        _plot_layers = set()
        for _t in p['traces']:
            _tm = _LAYER_RE.match(
                (_t.get('name', '').replace(_label_prefix, '')))
            if _tm:
                _plot_layers.add(_tm.group(2))
        _plot_layers_sorted = sorted(
            _plot_layers,
            key=lambda x: int(x[1:]) if x[1:].isdigit() else x)
        if len(_plot_layers_sorted) > 1:
            _opts = ''.join(
                f'<option value="{lyr}">{lyr}</option>'
                for lyr in _plot_layers_sorted
            )
            _layer_row = (
                f'<div class="plot-toolbar" data-plot="{div_id}" style="padding-top:0">'
                f'  <label style="font-weight:700">Layers:</label>'
                f'  <select class="tb-input" id="{div_id}_lyrFrom" '
                f'style="width:auto;min-width:60px;padding:2px 4px" '
                f'onchange="_owqApplyLayerRange(\'{div_id}\')">{_opts}</select>'
                f'  <label style="margin:0 2px">to</label>'
                f'  <select class="tb-input" id="{div_id}_lyrTo" '
                f'style="width:auto;min-width:60px;padding:2px 4px" '
                f'onchange="_owqApplyLayerRange(\'{div_id}\')">{_opts}</select>'
                f'  <button class="tb-btn" '
                f'onclick="_owqLayerAll(\'{div_id}\')">All</button>'
                f'</div>'
            )
        # Display title: use short species name when inside a compartment group
        _display_title = p['title']
        _nav_t = p.get('nav_title', '')
        if _nav_t and _cur_comp:
            _display_title = (f'<span style="color:var(--primary);font-weight:400">'
                              f'{separator.strip()}</span> {_nav_t}')
        H.append(f"""<div class="section" id="{p['id']}">
<h2>{_display_title}</h2>
<div class="plot-toolbar" data-plot="{div_id}">
  <button class="tb-btn" onclick="_owqToggleLog('{div_id}',this)">Y Log</button>
  <span class="tb-sep"></span>
  <label>Y:</label>
  <input type="number" class="tb-input" id="{div_id}_ymin" placeholder="min" step="any">
  <input type="number" class="tb-input" id="{div_id}_ymax" placeholder="max" step="any">
  <span class="tb-sep"></span>
  <label>X:</label>
  <input type="date" class="tb-input tb-date" id="{div_id}_xmin">
  <input type="date" class="tb-input tb-date" id="{div_id}_xmax">
  <span class="tb-sep"></span>
  <button class="tb-btn" onclick="_owqApplyRange('{div_id}')">Apply</button>
  <button class="tb-btn" onclick="_owqAxisTight('{div_id}')">Tight</button>
  <button class="tb-btn" onclick="_owqResetRange('{div_id}',this)">Reset</button>
</div>
<div class="plot-toolbar" data-plot="{div_id}" style="padding-top:0">
  <label style="font-weight:700">Active:</label>
{_trace_btns}
</div>
{_layer_row}
<div class="card">
<div id="{div_id}" style="width:100%;min-height:500px"></div>
<p class="plot-caption">{p.get('caption', '')}</p>
</div>
</div>""")

    if _prev_comp is not None:
        H.append('</div>')  # close last compartment group

    H.append('</div>')  # container
    H.append('</div>')  # main
    H.append('</div>')  # layout

    # --- JAVASCRIPT ---
    # Theme restore + sidebar resize run immediately (before deferred scripts)
    H.append("""<script>
(function(){
  var saved = localStorage.getItem('owq_theme');
  if(saved) document.documentElement.setAttribute('data-theme', saved);
})();
// Sidebar drag-to-resize
(function(){
  var handle=document.getElementById('sidebarHandle');
  if(!handle) return;
  var sidebar=handle.parentElement;
  var dragging=false;
  handle.addEventListener('mousedown',function(e){
    e.preventDefault();
    dragging=true;
    handle.classList.add('dragging');
    document.body.style.cursor='col-resize';
    document.body.style.userSelect='none';
  });
  document.addEventListener('mousemove',function(e){
    if(!dragging) return;
    var newW=e.clientX;
    if(newW<120) newW=120;
    if(newW>window.innerWidth*0.5) newW=window.innerWidth*0.5;
    sidebar.style.width=newW+'px';
  });
  document.addEventListener('mouseup',function(){
    if(!dragging) return;
    dragging=false;
    handle.classList.remove('dragging');
    document.body.style.cursor='';
    document.body.style.userSelect='';
  });
})();
</script>""")

    # Everything else waits for deferred scripts via DOMContentLoaded
    H.append('<script>')
    H.append('document.addEventListener("DOMContentLoaded",function(){')

    # Theme toggle
    H.append("""
window.toggleTheme = function(){
  var el = document.documentElement;
  var t = el.getAttribute('data-theme')==='dark' ? 'light' : 'dark';
  el.setAttribute('data-theme', t);
  localStorage.setItem('owq_theme', t);
  _owqRelayoutAll();
};
""")

    # Feature label prefix (from mapping_key) used by JS for name parsing
    H.append(f"\nvar _owqFeatureLabel = '{feature_label} ';\n")

    # Plotly layout helper
    H.append("""
function _owqLayout(extra){
  var t = document.documentElement.getAttribute('data-theme')==='dark';
  var base = {
    margin:{l:60,r:25,t:30,b:45},
    plot_bgcolor: t?'#1e1f2e':'#fff',
    paper_bgcolor: t?'#1e1f2e':'#fff',
    font:{color:t?'#ccc':'#333',family:'Inter,sans-serif',size:12},
    colorway:['#0066cc','#00a86b','#ff6b35','#004499','#34d399','#fb923c',
              '#667eea','#764ba2','#e63946','#2ec4b6','#e9c46a','#264653'],
    xaxis:{gridcolor:t?'#2a2b3d':'#eee'},
    yaxis:{gridcolor:t?'#2a2b3d':'#eee'},
    legend:{font:{size:11}}
  };
  if(extra){for(var k in extra){base[k]=extra[k];}}
  return base;
}

var PCFG = {responsive:true, displayModeBar:true,
  modeBarButtonsToRemove:['lasso2d','select2d'],
  toImageButtonOptions:{format:'svg',filename:'openwq_plot'}};

window._owqPlotIds = [];
function _owqRegister(id){ window._owqPlotIds.push(id); }
var _owqObsMeta = """ + json.dumps(_obs_meta) + """;

// Toggle observation traces on/off for a specific plot
window._owqToggleObs = function(btn){
  var plotId = btn.getAttribute('data-plot');
  var isActive = btn.classList.toggle('active');
  btn.textContent = isActive ? 'Observations \\u2713' : 'Observations';
  var gd = document.getElementById(plotId);
  if(!gd || !gd.data) return;
  var om = _owqObsMeta[plotId] || {};
  var selFids = window._owqSelectedFids || {};
  var hasSel = Object.keys(selFids).length > 0;
  var newVis = gd.data.map(function(t, idx){
    if(om[idx] !== undefined){
      if(!isActive) return false;
      if(hasSel) return selFids[om[idx]] ? true : false;
      return true;
    }
    return t.visible === undefined ? true : t.visible;
  });
  Plotly.restyle(gd, {visible: newVis});
};

function _owqRelayoutAll(){
  var t=document.documentElement.getAttribute('data-theme')==='dark';
  window._owqPlotIds.forEach(function(id){
    try{Plotly.relayout(id,{
      'plot_bgcolor':t?'#1e1f2e':'#fff',
      'paper_bgcolor':t?'#1e1f2e':'#fff',
      'font.color':t?'#ccc':'#333',
      'xaxis.gridcolor':t?'#2a2b3d':'#eee',
      'yaxis.gridcolor':t?'#2a2b3d':'#eee'
    });}catch(e){}
  });
}
""")

    # Helper: get active layers for a specific plot (reads its range dropdowns)
    H.append(f"""
var _owqAllLayers={json.dumps(_all_layers_sorted)};
function _owqGetPlotLayers(plotId){{
  var active={{}};
  var fromEl=document.getElementById(plotId+'_lyrFrom');
  var toEl=document.getElementById(plotId+'_lyrTo');
  if(!fromEl||!toEl) return active;
  function _idx(n){{ var m=n.match(/z(\\d+)/); return m?parseInt(m[1],10):0; }}
  var fi=_idx(fromEl.value), ti=_idx(toEl.value);
  var lo=Math.min(fi,ti), hi=Math.max(fi,ti);
  _owqAllLayers.forEach(function(lyr){{ var li=_idx(lyr); if(li>=lo&&li<=hi) active[lyr]=true; }});
  return active;
}}
""")

    # Plot toolbar controls: log scale toggle, axis range, reset
    H.append("""
window._owqToggleLog = function(plotId, btn){
  var gd=document.getElementById(plotId);
  var cur=(gd.layout&&gd.layout.yaxis&&gd.layout.yaxis.type)||'linear';
  var nxt=cur==='log'?'linear':'log';
  Plotly.relayout(plotId,{'yaxis.type':nxt});
  btn.classList.toggle('active',nxt==='log');
  btn.textContent=nxt==='log'?'Y Log \\u2713':'Y Log';
};
window._owqApplyRange = function(plotId){
  var u={};
  var ymin=document.getElementById(plotId+'_ymin').value;
  var ymax=document.getElementById(plotId+'_ymax').value;
  var xmin=document.getElementById(plotId+'_xmin').value;
  var xmax=document.getElementById(plotId+'_xmax').value;
  if(ymin!==''&&ymax!==''){u['yaxis.range']=[parseFloat(ymin),parseFloat(ymax)];u['yaxis.autorange']=false;}
  if(xmin&&xmax){u['xaxis.range']=[xmin,xmax];u['xaxis.autorange']=false;}
  if(Object.keys(u).length)Plotly.relayout(plotId,u);
};
window._owqResetRange = function(plotId,btn){
  Plotly.relayout(plotId,{'yaxis.autorange':true,'xaxis.autorange':true,'yaxis.type':'linear'});
  var tb=btn.parentElement;
  var lb=tb.querySelector('.tb-btn');
  if(lb){lb.classList.remove('active');lb.textContent='Y Log';}
  tb.querySelectorAll('.tb-input').forEach(function(inp){inp.value='';});
};
window._owqAxisTight = function(plotId){
  var gd=document.getElementById(plotId);
  if(!gd||!gd.data) return;
  var xmin=null,xmax=null,ymin=Infinity,ymax=-Infinity;
  gd.data.forEach(function(t){
    if(t.visible==='legendonly'||t.visible===false) return;
    if(!t.x||!t.y) return;
    t.y.forEach(function(v,i){
      if(v===null||v===undefined) return;
      v=parseFloat(v); if(isNaN(v)) return;
      if(v<ymin) ymin=v; if(v>ymax) ymax=v;
      var xv=t.x[i];
      if(xv!==null&&xv!==undefined){
        if(typeof xv==='string'&&xv.match(/^\\d{4}/)){
          if(!xmin||xv<xmin) xmin=xv;
          if(!xmax||xv>xmax) xmax=xv;
        }else{
          var xn=parseFloat(xv);
          if(!isNaN(xn)){if(xmin===null||xn<xmin)xmin=xn;if(xmax===null||xn>xmax)xmax=xn;}
        }
      }
    });
  });
  if(ymin===Infinity) return;
  var ypad=(ymax-ymin)*0.05; if(ypad===0) ypad=Math.abs(ymax)*0.05||1;
  var upd={'yaxis.range':[ymin-ypad,ymax+ypad],'yaxis.autorange':false};
  if(xmin!==null&&xmax!==null){upd['xaxis.range']=[xmin,xmax];upd['xaxis.autorange']=false;}
  Plotly.relayout(plotId,upd);
};
// Per-trace-type toggle: "conc" for concentrations, or a derivative label
// Uses a full-trace restyle (like _owqRefreshSelection) to keep all states consistent.
window._owqToggleDeriv = function(btn){
  var plotId=btn.getAttribute('data-plot');
  var derivLabel=btn.getAttribute('data-deriv');
  var gd=document.getElementById(plotId);
  if(!gd||!gd.data) return;
  var show=!btn.classList.contains('active');
  var displayText=derivLabel==='conc'?'Concentrations':derivLabel;
  btn.classList.toggle('active',show);
  btn.textContent=show?displayText+' \u2713':displayText;
  // Full restyle of ALL traces — respects all button states, feature selection, layers
  var active=_owqGetPlotLayers(plotId);
  var hasActiveLyr=Object.keys(active).length>0;
  var selFids=window._owqSelectedFids||{};
  var hasSel=Object.keys(selFids).length>0;
  var om=_owqObsMeta[plotId]||{};
  var obsBtn=document.querySelector('[data-plot="'+plotId+'"][data-obs]');
  var obsOn=obsBtn?obsBtn.classList.contains('active'):true;
  var derivBtns=document.querySelectorAll('[data-plot="'+plotId+'"][data-deriv]');
  var newVis=gd.data.map(function(t,tidx){
    // Observation traces
    if(om[tidx]!==undefined){
      if(!obsOn) return false;
      if(hasSel) return selFids[om[tidx]]?true:false;
      return true;
    }
    var lg=t.legendgroup||'';
    var isDebug=lg.indexOf('debug_')===0;
    var tname=(t.name||'').replace(_owqFeatureLabel,'');
    // Check trace-type toggle (conc or derivative)
    var typeOn=true;
    if(isDebug){
      typeOn=false;
      derivBtns.forEach(function(b){
        var dd=b.getAttribute('data-deriv');
        if(dd!=='conc' && lg.indexOf('debug_'+dd+'_')===0
           && b.classList.contains('active')) typeOn=true;
      });
    }else{
      var concBtn=document.querySelector('[data-plot="'+plotId+'"][data-deriv="conc"]');
      if(concBtn && !concBtn.classList.contains('active')) typeOn=false;
    }
    if(!typeOn) return false;
    // Feature selection filter
    if(hasSel){
      var info=_owqExtractHru(tname,lg);
      if(!selFids[info.hru]) return false;
    }
    // Layer filter
    if(hasActiveLyr){
      var m=tname.match(/\(z(\d+)\)/);
      if(m && !active['z'+m[1]]) return false;
    }
    return true;
  });
  Plotly.restyle(gd,{visible:newVis});
};
""")

    # Legend click: disabled — the map is the master controller for feature selection.
    # Clicking a legend item in a plot does nothing; use the map or map legend
    # to select/deselect features across all plots.
    H.append("""
function _owqExtractHru(tname, legendgroup){
  // Extract feature ID (hru) and optional layer from trace name.
  // Trace names: "ID", "ID (z1)", "ID dC/dt chem", "ID (z1) dC/dt chem"
  // For debug traces, the fid is embedded in legendgroup: "debug_{label}_{fid}"
  // 1) Try layer pattern first
  var m=tname.match(/^(.+?) \((z\d+)\)/);
  if(m) return {hru:m[1],layer:m[2]};
  // 2) For debug traces, extract fid from legendgroup
  if(legendgroup && legendgroup.indexOf('debug_')===0){
    // legendgroup format: "debug_{derivLabel}_{fid}"
    // derivLabel can contain '/' so we can't just split on '_'
    // Instead, find the fid at the end: it's the last '_'-separated segment
    var lastUs=legendgroup.lastIndexOf('_');
    if(lastUs>6) return {hru:legendgroup.substring(lastUs+1),layer:null};
  }
  return {hru:tname,layer:null};
}
function _owqLegendClick(gd){
  // Disable legend click — map is the master controller
  gd.on('plotly_legendclick',function(e){ return false; });
}
""")

    # Per-plot data — staggered via requestAnimationFrame for responsiveness
    H.append("""
var _owqPlotQueue = [];
function _owqFlushQueue(){
  if(!_owqPlotQueue.length) return;
  var job = _owqPlotQueue.shift();
  var gd = document.getElementById(job.id);
  Plotly.react(gd, job.traces, job.layout, PCFG);
  _owqRegister(job.id);
  _owqLegendClick(gd);
  if(_owqPlotQueue.length) requestAnimationFrame(_owqFlushQueue);
}
""")
    for p in plots:
        div_id = f'{p["id"]}_div'
        traces_json = json.dumps(p['traces'])
        xaxis_cfg = "title:'Time'"
        if p.get('xtype') == 'date':
            xaxis_cfg += ",type:'date'"

        H.append(f"""
_owqPlotQueue.push({{id:'{div_id}',traces:{traces_json},
  layout:_owqLayout({{xaxis:{{{xaxis_cfg}}},yaxis:{{title:'{p["ylabel"]}'}},
    hovermode:'closest'}})}});
""")
    H.append("requestAnimationFrame(_owqFlushQueue);")

    # --- LEAFLET MAP INITIALIZATION ---
    if river_geojson:
        # river_geojson is the primary interactive layer data (used for FID extraction)
        pass
        _fid_color_json = json.dumps(_fid_color)
        _fid_list_json = json.dumps(_all_fids_sorted)
        _hru_color_json = json.dumps(_hru_color)
        _hru_ids_json = json.dumps(_hru_ids_sorted)

        H.append(f"""
// --- River Network Map ---
(function(){{
  var mapEl=document.getElementById('mapContainer');
  if(!mapEl) return;
  // Basemap tile layers (Satellite, Light, Topo)
  var satTile=L.tileLayer('https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{{z}}/{{y}}/{{x}}',{{
    attribution:'Esri, Maxar, Earthstar Geographics',maxZoom:19}});
  var lightTile=L.tileLayer('https://{{s}}.basemaps.cartocdn.com/light_all/{{z}}/{{x}}/{{y}}{{r}}.png',{{
    attribution:'&copy; <a href="https://carto.com">CARTO</a>',maxZoom:18}});
  var topoTile=L.tileLayer('https://{{s}}.tile.opentopomap.org/{{z}}/{{x}}/{{y}}.png',{{
    attribution:'OpenTopoMap',maxZoom:17}});
  // Start with interactions disabled (locked view)
  var map=L.map('mapContainer',{{
    zoomControl:false,dragging:false,scrollWheelZoom:false,
    doubleClickZoom:false,touchZoom:false,boxZoom:false,keyboard:false
  }});
  satTile.addTo(map);  // default basemap
  L.control.layers({{'Satellite':satTile,'Light':lightTile,'Topo':topoTile}},{{}},{{collapsed:true,position:'bottomleft'}}).addTo(map);

  // Lock / unlock control (top-left)
  var mapLocked=true;
  var LockCtrl=L.Control.extend({{
    options:{{position:'topleft'}},
    onAdd:function(){{
      var btn=L.DomUtil.create('button','leaflet-bar');
      btn.style.cssText='width:34px;height:34px;background:var(--surface,#fff);border:none;'
        +'border-radius:4px;cursor:pointer;font-size:16px;line-height:34px;text-align:center;'
        +'box-shadow:0 1px 5px rgba(0,0,0,.3);';
      btn.title='Toggle map lock';
      btn.innerHTML='&#x1F512;';  // locked icon
      L.DomEvent.disableClickPropagation(btn);
      btn.addEventListener('click',function(){{
        mapLocked=!mapLocked;
        if(mapLocked){{
          map.dragging.disable();map.scrollWheelZoom.disable();
          map.doubleClickZoom.disable();map.touchZoom.disable();
          map.boxZoom.disable();map.keyboard.disable();
          btn.innerHTML='&#x1F512;';btn.title='Unlock map';
        }}else{{
          map.dragging.enable();map.scrollWheelZoom.enable();
          map.doubleClickZoom.enable();map.touchZoom.enable();
          map.boxZoom.enable();map.keyboard.enable();
          btn.innerHTML='&#x1F513;';btn.title='Lock map';
        }}
      }});
      return btn;
    }}
  }});
  new LockCtrl().addTo(map);
  // Show-all control (top-right) — restores all features in map + plots
  var ShowAllCtrl=L.Control.extend({{
    options:{{position:'topright'}},
    onAdd:function(){{
      var btn=L.DomUtil.create('button','leaflet-bar');
      btn.style.cssText='padding:4px 10px;background:var(--surface,#fff);border:none;'
        +'border-radius:4px;cursor:pointer;font-size:12px;line-height:20px;text-align:center;'
        +'box-shadow:0 1px 5px rgba(0,0,0,.3);color:var(--text,#333);white-space:nowrap;font-weight:600;';
      btn.title='Restore all {feature_label}s on map and plots';
      btn.innerHTML='Restore All';
      L.DomEvent.disableClickPropagation(btn);
      btn.addEventListener('click',function(){{
        if(typeof _owqSelectFeature==='function') _owqSelectFeature(null);
      }});
      return btn;
    }}
  }});
  new ShowAllCtrl().addTo(map);
  L.control.zoom({{position:'topleft'}}).addTo(map);

  var mapKey='{mapping_key}';
  var fidColor={_fid_color_json};
  var plotFids={_fid_list_json};
  var mapFidColor={json.dumps(_map_fid_color)};

  var featureLayers={{}};  // fid → Leaflet layer
  var selectedFids={{}};   // fid → true (multi-select set)
  window._owqSelectedFids=selectedFids;  // expose for layer selector

  var isPolygon={json.dumps(map_geom_type == 'polygon')};
  var hasLayers={json.dumps(_has_layers)};
  var allLayers={json.dumps(_all_layers_sorted)};
  var hruColor={json.dumps(_hru_color)};
  var currentLayer=hasLayers&&allLayers.length?allLayers[0]:null;

  // Map color lookup: try map-specific dict first, then plot dicts, then gray
  function _mapColor(fid){{
    return mapFidColor[fid]||fidColor[fid]||hruColor[fid]||'#888';
  }}

  function _defaultStyle(fid){{
    var c=_mapColor(fid);
    if(isPolygon){{
      return {{color:'#444',weight:1.5,opacity:0.7,fillColor:c,fillOpacity:0.5}};
    }}
    return {{color:c,weight:3,opacity:0.85}};
  }}
  function _highlightStyle(fid){{
    var c=_mapColor(fid);
    if(isPolygon){{
      return {{color:'#222',weight:2.5,opacity:1,fillColor:c,fillOpacity:0.75}};
    }}
    return {{color:c,weight:5,opacity:1}};
  }}
  function _dimStyle(fid){{
    var c=_mapColor(fid);
    if(isPolygon){{
      return {{color:'#999',weight:0.5,opacity:0.3,fillColor:c,fillOpacity:0.15}};
    }}
    return {{color:c,weight:1.5,opacity:0.3}};
  }}

  // === LAYER ORDER: Basin (back) → River Network (middle) → Stations (front) ===

  // Helper: make a layer interactive (clickable, with selection logic)
  function _makeInteractive(gjData){{
    return L.geoJSON(gjData,{{
      style:function(feature){{
        var fid=String(feature.properties[mapKey]||'');
        return _defaultStyle(fid);
      }},
      onEachFeature:function(feature,layer){{
        var fid=String(feature.properties[mapKey]||'');
        featureLayers[fid]=layer;
        layer.bindTooltip(_owqFeatureLabel+fid,{{sticky:true,direction:'top',opacity:0.9}});
        layer.on('click',function(){{ _owqToggleFeature(fid); }});
        layer.on('mouseover',function(){{ layer.setStyle(_highlightStyle(fid)); }});
        layer.on('mouseout',function(){{
          var hasSel=Object.keys(selectedFids).length>0;
          if(hasSel&&!selectedFids[fid]){{ layer.setStyle(_dimStyle(fid)); }}
          else{{ layer.setStyle(_defaultStyle(fid)); }}
        }});
      }}
    }});
  }}

  // 1) BASIN layer (always at the back)
  var _basinData={json.dumps(basin_geojson) if basin_geojson else 'null'};
  var _basinLayer=null;
  var _basinVisible=true;
  var _basinIsInteractive=isPolygon;  // interactive when primary is polygon (SUMMA)
  if(_basinData){{
    if(_basinIsInteractive){{
      _basinLayer=_makeInteractive(_basinData);
    }}else{{
      _basinLayer=L.geoJSON(_basinData,{{
        interactive:false,
        style:function(){{
          return {{color:'#555',weight:1.5,opacity:0.7,fillColor:'#b0b0b0',fillOpacity:0.25}};
        }}
      }});
    }}
    _basinLayer.addTo(map);
  }}

  // 2) RIVER NETWORK layer (on top of basins)
  var _riverData={json.dumps(river_line_geojson) if river_line_geojson else 'null'};
  var _riverLayer=null;
  var _riverVisible=true;
  var _riverIsInteractive=!isPolygon;  // interactive when primary is line (mizuRoute)
  if(_riverData){{
    if(_riverIsInteractive){{
      _riverLayer=_makeInteractive(_riverData);
    }}else{{
      _riverLayer=L.geoJSON(_riverData,{{
        interactive:false,
        style:function(){{
          return {{color:'#fff',weight:2.5,opacity:0.9}};
        }}
      }});
    }}
    _riverLayer.addTo(map);
  }}

  // Reference to the primary interactive layer (for fitBounds, bringToFront, etc.)
  var geoLayer=_basinIsInteractive?_basinLayer:_riverLayer;

  // --- Observation station markers ---
  var stationData={json.dumps({sid: list(loc) for sid, loc in (station_locations or {}).items()})};
  var stationToFeature={json.dumps(station_to_feature or {})};
  var stationMarkers={{}};  // sid → L.circleMarker
  var stationLayer=L.layerGroup();
  var stationsVisible=true;
  Object.keys(stationData).forEach(function(sid){{
    var ll=stationData[sid];
    var matchedFid=stationToFeature[sid]||null;
    var marker=L.circleMarker([ll[0],ll[1]],{{
      radius:4,fillColor:'#fff',color:'#e00',weight:1.5,opacity:1,fillOpacity:0.9
    }});
    var popupHtml='<b>Observation Station</b><br>ID: '+sid
      +'<br>Lat: '+ll[0].toFixed(4)+'<br>Lon: '+ll[1].toFixed(4);
    if(matchedFid) popupHtml+='<br>Matched {feature_label}: '+matchedFid;
    marker.bindTooltip('Station: '+sid+(matchedFid?' \\u2192 {feature_label} '+matchedFid:''),{{direction:'top',offset:[0,-6]}});
    marker.bindPopup(popupHtml);
    stationMarkers[sid]=marker;
    stationLayer.addLayer(marker);
  }});
  // Filter station markers when a feature is selected
  window._owqFilterStations=function(selSet){{
    if(!stationsVisible) return;
    var hasSel=Object.keys(selSet).length>0;
    Object.keys(stationMarkers).forEach(function(sid){{
      var mk=stationMarkers[sid];
      var matchedFid=stationToFeature[sid]||null;
      if(!hasSel){{
        mk.setStyle({{fillOpacity:0.85,opacity:1}});
      }}else if(matchedFid && selSet[String(matchedFid)]){{
        mk.setStyle({{fillOpacity:0.85,opacity:1}});
      }}else{{
        mk.setStyle({{fillOpacity:0.15,opacity:0.2}});
      }}
    }});
  }};
  if(Object.keys(stationData).length>0){{
    stationLayer.addTo(map);
  }}

  // Fit map tightly to data bounds (zero padding, delayed to ensure container is sized)
  var _fitBounds=null;
  try{{
    _fitBounds=geoLayer.getBounds();
    if(_basinLayer && _basinLayer!==geoLayer){{ _fitBounds.extend(_basinLayer.getBounds()); }}
    if(_riverLayer && _riverLayer!==geoLayer){{ _fitBounds.extend(_riverLayer.getBounds()); }}
    // Delay fitBounds so the map container has its final dimensions
    setTimeout(function(){{
      map.invalidateSize();
      map.fitBounds(_fitBounds,{{padding:[0,0]}});
    }}, 100);
  }}catch(e){{}}

  // Re-center control (top-left)
  var RecenterCtrl=L.Control.extend({{
    options:{{position:'topleft'}},
    onAdd:function(){{
      var btn=L.DomUtil.create('button','leaflet-bar');
      btn.style.cssText='width:34px;height:34px;background:var(--surface,#fff);border:none;'
        +'border-radius:4px;cursor:pointer;font-size:18px;line-height:34px;text-align:center;'
        +'box-shadow:0 1px 5px rgba(0,0,0,.3);color:var(--text,#333);';
      btn.title='Re-center map';
      btn.innerHTML='&#x2316;';  // ⌖ crosshair
      L.DomEvent.disableClickPropagation(btn);
      btn.addEventListener('click',function(){{
        if(_fitBounds) map.fitBounds(_fitBounds,{{padding:[0,0]}});
      }});
      return btn;
    }}
  }});
  new RecenterCtrl().addTo(map);

  // Helper: re-stack layers in correct z-order after toggling
  function _restackLayers(){{
    if(_basinLayer && _basinVisible) _basinLayer.bringToBack();
    if(_riverLayer && _riverVisible) _riverLayer.bringToFront();
    if(stationsVisible && Object.keys(stationData).length>0) stationLayer.bringToFront();
  }}

  // Build map legend as a Leaflet control (top-right, below Restore All)
  var _mapFidsSorted={json.dumps(sorted(_map_fid_color.keys()) if _map_fid_color else [])};
  var legendFids=_mapFidsSorted.length?_mapFidsSorted:(hasLayers?{_hru_ids_json}:plotFids);
  var legendEl=null;
  var LegendCtrl=L.Control.extend({{
    options:{{position:'topright'}},
    onAdd:function(){{
      var div=L.DomUtil.create('div','map-legend');
      L.DomEvent.disableClickPropagation(div);
      L.DomEvent.disableScrollPropagation(div);
      // Compute max-height: from control top to just above attribution (~30px from bottom)
      var mapH=mapEl.offsetHeight||420;
      // Leaflet topright container starts ~10px from top; Restore All is ~30px tall
      // Leave ~40px at bottom for attribution
      div.style.maxHeight=(mapH - 80)+'px';
      var html='';
      // Layer toggles section
      var _toggles=[];
      if(_basinLayer && !_basinIsInteractive){{
        _toggles.push('<div class="ml-toggle" title="Toggle basin polygons">'
          +'<input type="checkbox" checked id="mlChkBasin"><label for="mlChkBasin">Basins</label></div>');
      }}
      if(_riverLayer && !_riverIsInteractive){{
        _toggles.push('<div class="ml-toggle" title="Toggle river network">'
          +'<input type="checkbox" checked id="mlChkRiver"><label for="mlChkRiver">River Network</label></div>');
      }}
      if(Object.keys(stationData).length>0){{
        _toggles.push('<div class="ml-toggle" title="Toggle observation stations">'
          +'<input type="checkbox" checked id="mlChkSta"><label for="mlChkSta">Stations</label></div>');
      }}
      if(_toggles.length){{
        html+='<div class="ml-toggles">'+_toggles.join('')+'</div>'
             +'<div style="border-top:1px solid var(--border,#ddd);margin:4px 0"></div>';
      }}
      // Feature ID list
      html+='<div class="ml-title">{feature_label}s</div>';
      legendFids.forEach(function(fid){{
        var c=_mapColor(fid);
        html+='<div class="ml-item" data-fid="'+fid+'">'
             +'<span class="ml-swatch" style="background:'+c+';'
             +(isPolygon?'height:10px;width:16px;border-radius:3px;border:1px solid #666':'')
             +'"></span>'
             +'<span>'+fid+'</span></div>';
      }});
      div.innerHTML=html;
      return div;
    }}
  }});
  var _legendCtrl=new LegendCtrl();
  _legendCtrl.addTo(map);
  legendEl=_legendCtrl.getContainer();

  // Bind feature click handlers
  if(legendEl){{
    legendEl.querySelectorAll('.ml-item').forEach(function(el){{
      el.addEventListener('click',function(){{
        var fid=el.getAttribute('data-fid');
        _owqToggleFeature(fid);
      }});
    }});
    // Bind layer toggle handlers
    var chkBasin=document.getElementById('mlChkBasin');
    if(chkBasin){{
      chkBasin.addEventListener('change',function(){{
        _basinVisible=chkBasin.checked;
        if(_basinVisible){{ _basinLayer.addTo(map); _restackLayers(); }}
        else{{ map.removeLayer(_basinLayer); }}
      }});
    }}
    var chkRiver=document.getElementById('mlChkRiver');
    if(chkRiver){{
      chkRiver.addEventListener('change',function(){{
        _riverVisible=chkRiver.checked;
        if(_riverVisible){{ _riverLayer.addTo(map); _restackLayers(); }}
        else{{ map.removeLayer(_riverLayer); }}
      }});
    }}
    var chkSta=document.getElementById('mlChkSta');
    if(chkSta){{
      chkSta.addEventListener('change',function(){{
        stationsVisible=chkSta.checked;
        if(stationsVisible){{ stationLayer.addTo(map); _restackLayers(); }}
        else{{ map.removeLayer(stationLayer); }}
      }});
    }}
  }}

  // Toggle a feature in/out of the multi-select set
  window._owqToggleFeature=function(fid){{
    if(selectedFids[fid]){{ delete selectedFids[fid]; }}
    else{{ selectedFids[fid]=true; }}
    _owqRefreshSelection();
  }};

  // Clear all selections (called by Restore All)
  window._owqSelectFeature=function(arg){{
    // arg===null means clear all
    if(arg===null||arg===undefined){{
      Object.keys(selectedFids).forEach(function(k){{ delete selectedFids[k]; }});
    }}
    _owqRefreshSelection();
  }};

  // Master refresh — updates map styles, legend, stations, and all plot traces
  var _selTimer=null;
  function _owqRefreshSelection(){{
    var hasSel=Object.keys(selectedFids).length>0;
    // 1) Update map feature styles
    for(var f in featureLayers){{
      var ly=featureLayers[f];
      if(!hasSel){{ ly.setStyle(_defaultStyle(f)); }}
      else if(selectedFids[f]){{ ly.setStyle(_highlightStyle(f)); }}
      else{{ ly.setStyle(_dimStyle(f)); }}
    }}
    // 2) Filter station markers
    if(typeof _owqFilterStations==='function') _owqFilterStations(selectedFids);
    // 3) Update map legend highlights
    if(legendEl){{
      legendEl.classList.toggle('has-selection',hasSel);
      legendEl.querySelectorAll('.ml-item').forEach(function(el){{
        el.classList.toggle('active',!!selectedFids[el.getAttribute('data-fid')]);
      }});
    }}
    // 4) Debounce Plotly restyles (expensive)
    if(_selTimer) cancelAnimationFrame(_selTimer);
    _selTimer=requestAnimationFrame(function(){{
      var ids=window._owqPlotIds, i=0;
      function _nextRestyle(){{
        if(i>=ids.length) return;
        var gd=document.getElementById(ids[i++]);
        if(!gd||!gd.data){{ _nextRestyle(); return; }}
        var om=_owqObsMeta[gd.id]||{{}};
        var obsBtn=document.querySelector('[data-plot="'+gd.id+'"][data-obs]');
        var obsOn=obsBtn?obsBtn.classList.contains('active'):true;
        var active=_owqGetPlotLayers(gd.id);
        var hasActiveLyr=Object.keys(active).length>0;
        var derivBtns=document.querySelectorAll('[data-plot="'+gd.id+'"][data-deriv]');
        var newVis=gd.data.map(function(t,tidx){{
          var isDebug=t.legendgroup&&t.legendgroup.indexOf('debug_')===0;
          // Observation traces
          if(om[tidx]!==undefined){{
            if(!obsOn) return false;
            if(!hasSel) return true;
            return selectedFids[om[tidx]]?true:false;
          }}
          var lg=t.legendgroup||'';
          var tname=(t.name||'').replace(_owqFeatureLabel,'');
          var info=_owqExtractHru(tname,lg);
          // Hide traces from inactive layers
          if(hasActiveLyr && info.layer && !active[info.layer]) return false;
          // Check trace-type toggle (conc or derivative)
          var typeOn=true;
          if(isDebug){{
            typeOn=false;
            derivBtns.forEach(function(b){{
              var dd=b.getAttribute('data-deriv');
              if(dd!=='conc' && lg.indexOf('debug_'+dd+'_')===0
                 && b.classList.contains('active')) typeOn=true;
            }});
          }}else{{
            var concBtn=document.querySelector('[data-plot="'+gd.id+'"][data-deriv="conc"]');
            if(concBtn && !concBtn.classList.contains('active')) typeOn=false;
          }}
          if(!typeOn) return false;
          if(!hasSel) return true;
          // Multi-select: show only selected features, completely hide others
          return selectedFids[info.hru]?true:false;
        }});
        Plotly.restyle(gd,{{visible:newVis}}).then(_nextRestyle);
      }}
      _nextRestyle();
    }});
  }}

}})();
""")
    else:
        # No map — provide a no-op so legend click broadcasts don't error
        H.append("""
window._owqSelectFeature = window._owqSelectFeature || function(){};
""")

    # --- LAYER SELECTOR LOGIC (per-plot range dropdowns) ---
    H.append(f"""
// Layer selector (SUMMA multi-layer support — per-plot range)
(function(){{

  // Apply layer range for a single plot
  window._owqApplyLayerRange=function(plotId){{
    var active=_owqGetPlotLayers(plotId);
    var gd=document.getElementById(plotId);
    if(!gd||!gd.data) return;
    var btns=document.querySelectorAll('[data-plot="'+plotId+'"][data-deriv]');
    var om=_owqObsMeta[plotId]||{{}};
    var obsBtn=document.querySelector('[data-plot="'+plotId+'"][data-obs]');
    var obsOn=obsBtn?obsBtn.classList.contains('active'):true;
    var selFids=window._owqSelectedFids||{{}};
    var hasSel=Object.keys(selFids).length>0;
    var newVis=gd.data.map(function(t,tidx){{
      // Observation traces
      if(om[tidx]!==undefined){{
        if(!obsOn) return false;
        if(hasSel) return selFids[om[tidx]]?true:false;
        return true;
      }}
      var lg=t.legendgroup||'';
      var isDebug=lg.indexOf('debug_')===0;
      var tname=(t.name||'').replace(_owqFeatureLabel,'');
      var typeOn=true;
      if(isDebug){{
        typeOn=false;
        btns.forEach(function(b){{
          var d=b.getAttribute('data-deriv');
          if(d!=='conc' && lg.indexOf('debug_'+d+'_')===0
             && b.classList.contains('active')) typeOn=true;
        }});
      }}else{{
        var concBtn=document.querySelector('[data-plot="'+plotId+'"][data-deriv="conc"]');
        if(concBtn && !concBtn.classList.contains('active')) typeOn=false;
      }}
      if(!typeOn) return false;
      // Feature selection filter
      if(hasSel){{
        var info=_owqExtractHru(tname,lg);
        if(!selFids[info.hru]) return false;
      }}
      var m=tname.match(/\\(z(\\d+)\\)/);
      if(m){{
        var traceLayer='z'+m[1];
        return active[traceLayer]===true?true:false;
      }}
      return true;
    }});
    Plotly.restyle(gd,{{visible:newVis}});
  }};

  // "All" button — set range to first..last option in this plot's dropdowns
  window._owqLayerAll=function(plotId){{
    var fromEl=document.getElementById(plotId+'_lyrFrom');
    var toEl=document.getElementById(plotId+'_lyrTo');
    if(!fromEl||!toEl) return;
    fromEl.selectedIndex=0;
    toEl.selectedIndex=toEl.options.length-1;
    _owqApplyLayerRange(plotId);
  }};

  // Backwards compat stubs
  window._owqToggleLayer=function(){{}};
  window._owqSelectLayer=function(){{}};
  window._owqActiveLayers={{}};
}})();
""")

    # Sidebar active state on scroll (throttled) + click-to-scroll
    # NOTE: .main is the scroll container (overflow-y:auto), not window.
    H.append("""
(function(){
  var mainEl = document.querySelector('.main');
  if(!mainEl) return;
  var links = document.querySelectorAll('.sidebar nav a');
  // Cache section elements once
  var sectionMap = [];
  links.forEach(function(link){
    var id = link.getAttribute('href');
    if(!id||!id.startsWith('#')) return;
    var sec = document.querySelector(id);
    if(sec) sectionMap.push({link:link, sec:sec});
  });

  // Click handler: smooth-scroll to section inside .main
  links.forEach(function(link){
    link.addEventListener('click', function(e){
      var id = link.getAttribute('href');
      if(!id||!id.startsWith('#')) return;
      e.preventDefault();
      var target = document.querySelector(id);
      if(target){
        var topFixed = document.querySelector('.top-fixed');
        var offset = topFixed ? topFixed.offsetHeight + 10 : 0;
        var y = target.offsetTop - offset;
        mainEl.scrollTo({top: y, behavior: 'smooth'});
      }
    });
  });

  // Throttled scroll listener for active state (on .main, not window)
  var _scrollTick = false;
  mainEl.addEventListener('scroll', function(){
    if(_scrollTick) return;
    _scrollTick = true;
    requestAnimationFrame(function(){
      var topFixed = document.querySelector('.top-fixed');
      var fixedH = topFixed ? topFixed.offsetHeight : 0;
      var fromTop = mainEl.scrollTop + fixedH + 20;
      sectionMap.forEach(function(item){
        // Skip hidden sections (filtered out by layer selector)
        if(item.sec.style.display==='none'||item.link.style.display==='none'){
          item.link.classList.remove('active');
          return;
        }
        var t = item.sec.offsetTop, h = item.sec.offsetHeight;
        if(t <= fromTop && t + h > fromTop){
          item.link.classList.add('active');
        } else {
          item.link.classList.remove('active');
        }
      });
      _scrollTick = false;
    });
  });
})();
""")

    # Auto-apply first layer on page load (after plots are rendered)
    if _has_layers:
        H.append("""
// Auto-apply first layer to plots that have layer dropdowns
setTimeout(function(){
  window._owqPlotIds.forEach(function(id){
    var fromEl=document.getElementById(id+'_lyrFrom');
    var toEl=document.getElementById(id+'_lyrTo');
    if(!fromEl||!toEl) return;  // no layers for this plot
    // Set both dropdowns to the first available option
    fromEl.selectedIndex=0;
    toEl.selectedIndex=0;
    if(typeof _owqApplyLayerRange==='function') _owqApplyLayerRange(id);
  });
}, 500);
""")

    # Close DOMContentLoaded wrapper
    H.append('});')  # end DOMContentLoaded
    H.append('</script>')
    H.append('</body></html>')

    return '\n'.join(H)


# ---------------------------------------------------------------------------
# Main driver
# ---------------------------------------------------------------------------

def Plot_h5_driver(what2map=None,
                   hostmodel=None,
                   mapping_key_values=None,
                   openwq_results=None,
                   chemSpec=None,
                   hydromodel_info=None,
                   hydromodel_var2print=None,
                   output_path=None,
                   debugmode=False,
                   sediment_as_well=False,
                   combine_features=True,
                   figsize=(12, 6),
                   river_network_shp=None,
                   basin_shapefile=None,
                   basin_mapping_key=None,
                   mapping_key='SegId',
                   feature_label=None,
                   observation_dir=None,
                   observation_csv=None,
                   observation_compartments=None,
                   separator=' | '):
    """
    Generate interactive HTML time-series plots (Plotly.js).

    Produces a single self-contained HTML file with one interactive chart
    per chemical species (or host-model variable).  The file uses Plotly.js
    loaded from CDN — no extra Python packages are required.

    Parameters
    ----------
    what2map : str
        'hostmodel' or 'openwq'
    hostmodel : str
        'mizuroute' or 'summa'
    mapping_key_values : list
        List of feature IDs to plot (e.g., [200005899, 200002273])
        Use ["all"] to plot every feature in the data.
    openwq_results : dict
        OpenWQ results dictionary from Read_h5_driver.
        Used if what2map='openwq'.
    chemSpec : list or str
        Chemical species to plot (e.g., ['NO3-N', 'NH4-N'] or 'NO3-N').
        Used if what2map='openwq'.
    hydromodel_info : dict
        Dictionary with 'path_to_results' and 'mapping_key'.
        Used if what2map='hostmodel'.
    hydromodel_var2print : str
        Variable to plot (e.g., 'basRunoff', 'DWroutedRunoff').
        Used if what2map='hostmodel'.
    output_path : str
        Path where the HTML file will be saved.
        Legacy .png extensions are automatically converted to .html.
    debugmode : bool
        If True, process all debug extensions in addition to 'main'.
    sediment_as_well : bool
        If True, also plot sediment transport results.
    combine_features : bool
        Kept for backward compatibility (ignored - features are always combined).
    figsize : tuple
        Kept for backward compatibility (ignored - Plotly charts are responsive).
    separator : str
        Separator between compartment and species names in plot titles.
        Default is ' | '. Examples: ' | ', ' - ', ' :: '.

    Returns
    -------
    str or None
        Path to the generated HTML file, or None if no plots were created.
    """

    # Default feature_label to mapping_key if not provided
    if feature_label is None:
        feature_label = mapping_key

    print("=" * 70)
    print("TIME-SERIES PLOTTING (Interactive HTML)")
    print(f"Mode: {what2map.upper()}")
    if what2map == 'openwq':
        print(f"Debug mode: {debugmode}")
    print("=" * 70)

    # Validate inputs
    if what2map not in ['hostmodel', 'openwq']:
        print(f"✗ Error: what2map must be 'hostmodel' or 'openwq', got '{what2map}'")
        return None

    if mapping_key_values is None or len(mapping_key_values) == 0:
        print("✗ Error: mapping_key_values must be provided")
        return None

    # Convert to list if single value
    if not isinstance(mapping_key_values, list):
        mapping_key_values = [mapping_key_values]

    print(f"Features to plot: {mapping_key_values}")

    # Resolve output path (.png → .html migration)
    if output_path is None:
        output_path = 'timeseries_results.html'
    else:
        base, ext = os.path.splitext(output_path)
        if ext.lower() == '.png':
            output_path = base + '.html'
            print(f"  NOTE: Output format changed from .png to .html")
            print(f"        Interactive plots are now saved as: {output_path}")
        elif ext.lower() != '.html':
            output_path = base + '.html'

    # Create output directory if needed
    output_dir = os.path.dirname(output_path)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)

    # Collect all plot data
    plots = []

    # =====================================================================
    # HOSTMODEL MODE
    # =====================================================================
    if what2map == 'hostmodel':
        print("\n" + "=" * 70)
        print("Loading hostmodel data...")
        print("=" * 70)

        if hydromodel_info is None:
            print("✗ Error: hydromodel_info must be provided for hostmodel mode")
            return None

        if hydromodel_var2print is None:
            print("✗ Error: hydromodel_var2print must be specified for hostmodel mode")
            return None

        nc_path = hydromodel_info['path_to_results']
        mapping_key = hydromodel_info['mapping_key']

        print(f"  NetCDF: {nc_path}")
        print(f"  Variable: {hydromodel_var2print}")
        print(f"  Mapping key: {mapping_key}")

        nc = netCDF4.Dataset(nc_path, 'r')

        try:
            all_ids = np.array(nc.variables[mapping_key][:])
            print(f"  Total features in file: {len(all_ids)}")

            data_var = nc.variables[hydromodel_var2print]

            # Get time index
            if 'time' in nc.variables:
                time_var = nc.variables['time']
                time_data = time_var[:]
                if hasattr(time_var, 'units'):
                    from netCDF4 import num2date
                    dates = num2date(time_data, time_var.units,
                                     only_use_cftime_datetimes=False)
                    time_index = pd.DatetimeIndex(dates)
                else:
                    time_index = pd.RangeIndex(len(time_data))
            else:
                time_index = pd.RangeIndex(data_var.shape[0])

            print(f"  Time range: {time_index[0]} to {time_index[-1]}")
            print(f"  Time steps: {len(time_index)}")

            units = (nc.variables[hydromodel_var2print].units
                     if hasattr(nc.variables[hydromodel_var2print], 'units')
                     else 'units')

            # Extract feature data
            feature_data = {}
            missing_features = []

            for fid in mapping_key_values:
                idx = np.where(all_ids == fid)[0]
                if len(idx) > 0:
                    data = data_var[:, idx[0]]
                    feature_data[fid] = pd.Series(data, index=time_index)
                    print(f"  ✓ Found data for feature {fid}")
                else:
                    missing_features.append(fid)
                    print(f"  ✗ Feature {fid} not found")

            nc.close()

            if len(feature_data) == 0:
                print("✗ Error: No data found for any requested features")
                return None

            if missing_features:
                print(f"\n⚠ Warning: Missing features: {missing_features}")

            # Build plot entry
            traces = _build_traces(feature_data, feature_label=feature_label)
            plot_id = f'plot_{_sanitize_id(hydromodel_var2print)}'
            time_range = (f'{time_index[0]} to {time_index[-1]}'
                          if len(time_index) > 0 else '')
            plots.append({
                'id': plot_id,
                'title': f'{hostmodel.upper()}{separator}{hydromodel_var2print}',
                'traces': traces,
                'ylabel': f'{hydromodel_var2print} ({units})',
                'xtype': 'date' if isinstance(time_index, pd.DatetimeIndex) else 'linear',
                'caption': f'{len(feature_data)} {feature_label}s | {time_range}',
            })

        except Exception as e:
            print(f"✗ Error loading hostmodel data: {e}")
            nc.close()
            return None

    # =====================================================================
    # OPENWQ MODE
    # =====================================================================
    else:
        print("\n" + "=" * 70)
        print("Loading OpenWQ data...")
        print("=" * 70)

        if openwq_results is None:
            print("✗ Error: openwq_results must be provided for openwq mode")
            return None

        if chemSpec is None:
            print("✗ Error: chemSpec must be specified for openwq mode")
            return None

        if isinstance(chemSpec, str):
            chemSpec = [chemSpec]

        if sediment_as_well and 'Sediment' not in chemSpec:
            chemSpec = list(chemSpec) + ['Sediment']

        # Extensions to process
        if debugmode:
            file_extensions = [
                'main',
                'd_output_dt_chemistry',
                'd_output_dt_transport',
                'd_output_ss',
                'd_output_ewf',
                'd_output_ic'
            ]
        else:
            file_extensions = ['main']

        # Helper: parse column name into (feature_id, layer) tuple.
        # SUMMA columns: "hruId_1_z1" → fid="1", layer="z1"
        # mizuRoute columns: "reachID_740493340" → fid="740493340", layer=None
        def _parse_col(col):
            """Parse a DataFrame column name into (fid, layer).

            For SUMMA-style columns like ``hruId_1_z1`` the value part
            is ``1_z1``.  We detect the ``_z\d+`` suffix to split out
            the layer; everything before it is the feature ID.

            For mizuRoute-style columns like ``reachID_740493340`` there
            is no ``_z\d+`` suffix, so we return ``(740493340, None)``.
            """
            # Strip the mapping-key prefix (e.g. "hruId_" or "reachID_")
            parts = col.split('_', 1)  # ['hruId', '1_z1'] or ['reachID', '740493340']
            if len(parts) < 2:
                return (col, None)
            value = parts[1]  # '1_z1' or '740493340'
            # Check for layer suffix _z<digits>
            m = re.match(r'^(.+)_(z\d+)$', value)
            if m:
                return (m.group(1), m.group(2))
            return (value, None)

        # Helper: extract feature data from a DataFrame
        def _extract_features(ttdata, mkv):
            """Extract feature data, returning {display_id: Series}.

            For SUMMA (with layers), display_id is ``"{hru} (z{n})"``
            so that each HRU+layer pair gets its own trace with a
            readable legend name.
            For mizuRoute (no layers), display_id is the raw feature ID.

            Also returns a set of unique layers found.
            """
            fd, miss = {}, []
            _layers_found = set()

            if mkv == ["all"] or mkv == "all":
                for col in ttdata.columns:
                    fid, layer = _parse_col(col)
                    if layer:
                        _layers_found.add(layer)
                        display_id = f'{fid} ({layer})'
                    else:
                        display_id = fid
                    fd[display_id] = ttdata[col]
            else:
                for target_fid in mkv:
                    found_any = False
                    for col in ttdata.columns:
                        fid, layer = _parse_col(col)
                        if str(fid) == str(target_fid):
                            if layer:
                                _layers_found.add(layer)
                                display_id = f'{fid} ({layer})'
                            else:
                                display_id = str(fid)
                            fd[display_id] = ttdata[col]
                            found_any = True
                    if not found_any:
                        miss.append(target_fid)
            return fd, miss, _layers_found

        # ── Discover compartments from openwq_results keys ──────────
        # Keys are "COMPARTMENT@CHEMICAL#UNITS" or "COMPARTMENT@Sediment"
        _compartments_seen = []
        _comp_set = set()
        for key in openwq_results.keys():
            comp = key.split('@')[0]
            if comp not in _comp_set:
                _comp_set.add(comp)
                _compartments_seen.append(comp)

        print(f"  Compartments: {_compartments_seen}")
        print(f"  Chemical species: {chemSpec}")
        print(f"  Extensions: {file_extensions}")
        if debugmode:
            print(f"  Debug mode ON — debug traces overlaid on each plot")

        # ── Pre-pass: collect all fids to build a global colour map ──
        _all_fids_for_colors = set()
        for _pk in openwq_results.keys():
            for _pext, _pdata_list in openwq_results[_pk]:
                if _pext == 'main' and _pdata_list and len(_pdata_list) > 0:
                    _, _ptt, _ = _pdata_list[0]
                    _pfd, _, _ = _extract_features(_ptt, mapping_key_values)
                    _all_fids_for_colors.update(_pfd.keys())
        _global_color_map = _build_global_color_map(_all_fids_for_colors)
        print(f"  Global colour map: {len(_global_color_map)} fids")

        # ── Loop: compartment → species ─────────────────────────────
        for comp in _compartments_seen:

            for spec in chemSpec:
                # Build the result key
                # Normal species: "COMP@CHEM#UNITS"
                # Sediment:       "COMP@Sediment"
                result_key = None
                for key in openwq_results.keys():
                    kcomp = key.split('@')[0]
                    if kcomp != comp:
                        continue
                    if spec == 'Sediment':
                        if key == f'{comp}@Sediment':
                            result_key = key
                            break
                    elif spec in key:
                        result_key = key
                        break

                if result_key is None:
                    continue

                print(f"\n{'=' * 70}")
                print(f"  {comp} / {spec}  (key={result_key})")
                print("=" * 70)

                try:
                    # Parse units & chemical name from key
                    parts = result_key.split('#')
                    units = parts[1] if len(parts) > 1 else (
                        'kg' if 'Sediment' in result_key else 'concentration')
                    chem_parts = result_key.split('@')
                    chem_name = (chem_parts[1].split('#')[0]
                                 if len(chem_parts) > 1 else spec)

                    # ── Extract MAIN data ─────────────────────────
                    main_ttdata = None
                    for ext, data_list in openwq_results[result_key]:
                        if ext == 'main' and data_list and len(data_list) > 0:
                            _, main_ttdata, _ = data_list[0]
                            print(f"  [main] shape={main_ttdata.shape}  "
                                  f"{main_ttdata.index[0]} – "
                                  f"{main_ttdata.index[-1]}")
                            break

                    if main_ttdata is None:
                        print(f"  ✗ No 'main' data — skipped")
                        continue

                    main_fd, main_miss, layers_found = _extract_features(
                        main_ttdata, mapping_key_values)
                    if not main_fd:
                        print(f"  ✗ No feature data — skipped")
                        continue
                    if main_miss:
                        print(f"  ⚠ Missing features: {main_miss}")
                    if layers_found:
                        print(f"  Layers detected: {sorted(layers_found)}")

                    # Build main traces with global colours
                    all_traces, color_map = _build_traces_with_colors(
                        main_fd, feature_label=feature_label,
                        global_color_map=_global_color_map)

                    # ── Overlay debug extensions ──────────────────
                    _debug_labels_found = []  # list of labels present
                    if debugmode:
                        for dbg_ext in file_extensions[1:]:
                            meta = _DEBUG_EXT_META.get(dbg_ext)
                            if meta is None:
                                continue

                            dbg_ttdata = None
                            for ext, data_list in openwq_results[result_key]:
                                if (ext == dbg_ext and data_list
                                        and len(data_list) > 0):
                                    _, dbg_ttdata, _ = data_list[0]
                                    break

                            if dbg_ttdata is None:
                                continue

                            dbg_fd, _, _ = _extract_features(
                                dbg_ttdata, mapping_key_values)
                            dbg_fd = {k: v for k, v in dbg_fd.items()
                                      if k in color_map}
                            if not dbg_fd:
                                continue

                            dbg_traces = _build_debug_traces(
                                dbg_fd, color_map,
                                meta['label'], meta['dash'],
                                feature_label=feature_label)
                            all_traces.extend(dbg_traces)
                            _debug_labels_found.append(meta['label'])
                            print(f"  [debug] {meta['label']}: "
                                  f"{len(dbg_fd)} feature(s)")

                    # ── Append plot entry ─────────────────────────
                    plot_id = (f'plot_{_sanitize_id(comp)}_'
                               f'{_sanitize_id(spec)}_main')
                    time_range = (
                        f'{main_ttdata.index[0]} to '
                        f'{main_ttdata.index[-1]}'
                        if len(main_ttdata.index) > 0 else '')
                    plots.append({
                        'id': plot_id,
                        'title': f'{comp}{separator}{chem_name}',
                        'nav_title': chem_name,
                        'traces': all_traces,
                        'ylabel': f'{chem_name} ({units})',
                        'xtype': ('date'
                                  if isinstance(main_ttdata.index,
                                                pd.DatetimeIndex)
                                  else 'linear'),
                        'caption': (f'{len(main_fd)} {feature_label}s'
                                    f' | {time_range}'),
                        'species': spec,
                        'compartment': comp,
                        'debug_labels': _debug_labels_found,
                    })

                    print(f"  ✓ Plot: {comp}{separator}{chem_name}"
                          + (f" (+{', '.join(_debug_labels_found)})"
                             if _debug_labels_found else ""))

                except Exception as e:
                    print(f"  ✗ Error: {comp}/{spec}: {e}")
                    import traceback
                    traceback.print_exc()
                    continue

    # =====================================================================
    # BUILD HTML AND WRITE
    # =====================================================================
    if len(plots) == 0:
        print("\n✗ No plots were created.")
        print("  Check that the model has been run and HDF5 output files exist")
        print(f"  for the requested species: {chemSpec}")
        # Remove stale report from a previous run so it cannot be
        # accidentally opened and mistaken for current results.
        if output_path and os.path.isfile(output_path):
            try:
                os.remove(output_path)
                print(f"  Removed stale report: {output_path}")
            except OSError:
                pass
        return None

    print(f"\n{'=' * 70}")
    print(f"Building interactive HTML with {len(plots)} plot(s)...")
    print("=" * 70)

    # Load shapefile for interactive map (optional)
    # For SUMMA: prefer basin_shapefile (HRU polygons) over river_network_shp
    # For mizuRoute: use river_network_shp (line segments)
    _river_geojson = None
    _river_bounds = None
    _basin_geojson = None
    _basin_bounds = None
    _map_center = None
    _map_bounds = None
    _map_geom_type = 'line'  # 'line' for rivers, 'polygon' for basins

    def _load_shapefile(shp_path, label='shapefile'):
        """Load a shapefile and return (geojson_dict, bounds) or (None, None)."""
        try:
            import fiona
            from shapely.geometry import shape, mapping as shp_mapping

            features = []
            src_proj4 = None
            with fiona.open(shp_path) as src:
                if src.crs:
                    if isinstance(src.crs, dict):
                        parts = []
                        for k, v in src.crs.items():
                            parts.append(f'+{k}' if v is True else f'+{k}={v}')
                        src_proj4 = ' '.join(parts)
                    else:
                        try:
                            src_proj4 = src.crs.to_proj4()
                        except Exception:
                            src_proj4 = None
                for feat in src:
                    geom = shape(feat['geometry'])
                    features.append({
                        'type': 'Feature',
                        'geometry': geom,
                        'properties': dict(feat['properties']),
                    })

            if features:
                sb = features[0]['geometry'].bounds
                if (abs(sb[0]) > 360 or abs(sb[1]) > 360) and src_proj4:
                    print(f"  Reprojecting {label} to WGS84...")
                    try:
                        from pyproj import Proj
                        from shapely.ops import transform as shp_transform
                        src_p = Proj(src_proj4)
                        def _to_wgs(x, y):
                            return src_p(x, y, inverse=True)
                        for f in features:
                            f['geometry'] = shp_transform(_to_wgs, f['geometry'])
                    except Exception as _re:
                        print(f"  WARNING: Reprojection failed: {_re}")

                gj_features = []
                for f in features:
                    gj_features.append({
                        'type': 'Feature',
                        'geometry': shp_mapping(f['geometry']),
                        'properties': {k: (str(v) if v is not None else '')
                                       for k, v in f['properties'].items()},
                    })
                geojson = {'type': 'FeatureCollection',
                           'features': gj_features}
                all_b = [f['geometry'].bounds for f in features]
                bounds = (min(b[0] for b in all_b), min(b[1] for b in all_b),
                          max(b[2] for b in all_b), max(b[3] for b in all_b))
                print(f"  ✓ Loaded {label}: {len(gj_features)} features")
                return geojson, bounds
        except ImportError:
            print(f"  WARNING: fiona/shapely not installed — skipping {label}")
        except Exception as _e:
            print(f"  WARNING: Failed to load {label}: {_e}")
        return None, None

    # Load basin shapefile (polygons)
    if basin_shapefile and os.path.isfile(basin_shapefile):
        _basin_geojson, _basin_bounds = _load_shapefile(
            basin_shapefile, 'basin/HRU shapefile')

    # Load river network shapefile (lines)
    if river_network_shp and os.path.isfile(river_network_shp):
        _river_geojson, _river_bounds = _load_shapefile(
            river_network_shp, 'river network')

    # Determine primary interactive layer based on host model
    # mizuRoute: river network is clickable (line)
    # SUMMA: basin is clickable (polygon)
    # Both layers are always shown (when available) in order: basin → river → stations
    _hm = (hostmodel or '').lower()
    if _hm in ('mizuroute', 'mizuroute_cslm'):
        _primary_geojson = _river_geojson or _basin_geojson
        _map_geom_type = 'line' if _river_geojson else 'polygon'
    else:
        _primary_geojson = _basin_geojson or _river_geojson
        _map_geom_type = 'polygon' if _basin_geojson else 'line'

    # Compute map bounds from whichever layers are available
    _all_bounds = []
    if _basin_geojson and _basin_bounds:
        _all_bounds.append(_basin_bounds)
    if _river_geojson and _river_bounds:
        _all_bounds.append(_river_bounds)
    if _all_bounds:
        _map_bounds = (min(b[0] for b in _all_bounds),
                       min(b[1] for b in _all_bounds),
                       max(b[2] for b in _all_bounds),
                       max(b[3] for b in _all_bounds))
        _map_center = [(_map_bounds[1]+_map_bounds[3])/2,
                       (_map_bounds[0]+_map_bounds[2])/2]

    # Auto-detect basin mapping key from actual GeoJSON properties.
    # Strategy: find a property whose values overlap with the HRU IDs
    # extracted from the HDF5 data columns.
    _basin_mapping_key = basin_mapping_key or 'HRU_ID'
    if _basin_geojson and _basin_geojson.get('features'):
        _props = _basin_geojson['features'][0].get('properties', {})

        # Collect HRU IDs from plot data for matching
        _plot_hru_ids = set()
        if plots:
            for p in plots:
                for t in p.get('traces', []):
                    if t.get('legendgroup', '').startswith('debug_'):
                        continue
                    tname = t.get('name', '')
                    if feature_label and tname.startswith(feature_label + ' '):
                        tfid = tname[len(feature_label) + 1:]
                    else:
                        tfid = tname
                    m = re.match(r'^(.+?) \(z\d+\)$', tfid)
                    _plot_hru_ids.add(m.group(1) if m else tfid)

        # First try the requested key
        _found_match = False
        if _basin_mapping_key in _props:
            _map_vals = set()
            for feat in _basin_geojson['features']:
                _map_vals.add(str(feat['properties'].get(_basin_mapping_key, '')))
            if _plot_hru_ids & _map_vals:
                _found_match = True
                print(f"  Basin mapping key '{_basin_mapping_key}' matches plot data ✓")

        # If no match, try all properties to find one that overlaps
        if not _found_match:
            for _try_key in list(_props.keys()):
                _map_vals = set()
                for feat in _basin_geojson['features']:
                    _map_vals.add(str(feat['properties'].get(_try_key, '')))
                if _plot_hru_ids and (_plot_hru_ids & _map_vals):
                    _basin_mapping_key = _try_key
                    _found_match = True
                    print(f"  Auto-detected basin mapping key: '{_try_key}' "
                          f"(matches plot HRU IDs)")
                    break

            if not _found_match and _basin_mapping_key not in _props:
                # Fallback: use first available property
                for _try_key in ['HRU_ID', 'GRU_ID', 'hru_id', 'gru_id']:
                    if _try_key in _props:
                        _basin_mapping_key = _try_key
                        break
                print(f"  WARNING: No basin property matches plot HRU IDs. "
                      f"Using '{_basin_mapping_key}'. Map colours may be wrong.")

    # --- Load observation data and match to features ---
    _observation_data = None    # {plot_id: [obs_dict, ...]}
    station_locations = None     # {station_id: (lat, lon)}
    station_to_feature = None    # {station_id: feature_id}
    if (observation_dir or observation_csv) and _primary_geojson:
        print("\n" + "-" * 40)
        print("Loading observation data...")
        obs_by_species, station_locations = _load_observation_data(
            obs_dir=observation_dir, obs_csv=observation_csv)
        if obs_by_species and station_locations:
            _obs_geojson = _primary_geojson
            _obs_map_key = _basin_mapping_key if _basin_geojson else mapping_key
            station_to_feature = _match_stations_to_features(
                station_locations, _obs_geojson, _obs_map_key)
            if station_to_feature:
                print(f"  Matched {len(station_to_feature)} stations to features:")
                for sid, fid in station_to_feature.items():
                    print(f"    Station {sid} → Feature {fid}")
                _observation_data = {}
                # Normalise compartment filter to a set of lowercase names
                _obs_comps = None
                if observation_compartments:
                    if isinstance(observation_compartments, str):
                        _obs_comps = {observation_compartments.strip().lower()}
                    else:
                        _obs_comps = {c.strip().lower()
                                      for c in observation_compartments}
                for p in plots:
                    # Filter by compartment if specified
                    if _obs_comps:
                        p_comp = (p.get('compartment') or '').strip().lower()
                        if p_comp not in _obs_comps:
                            continue
                    p_species = p.get('species')
                    if not p_species:
                        continue
                    species_obs = obs_by_species.get(p_species)
                    if not species_obs:
                        continue
                    # Group by station
                    stn_groups = defaultdict(lambda: {'x': [], 'y': []})
                    for rec in species_obs:
                        stn_groups[rec['station_id']]['x'].append(rec['datetime'])
                        stn_groups[rec['station_id']]['y'].append(rec['value'])
                    plot_obs = []
                    for stn_id, data in stn_groups.items():
                        fid = station_to_feature.get(stn_id)
                        if fid is None:
                            continue
                        plot_obs.append({
                            'x': data['x'],
                            'y': data['y'],
                            'station_id': stn_id,
                            'feature_id': str(fid),
                        })
                    if plot_obs:
                        _observation_data[p['id']] = plot_obs
                        print(f"  ✓ {p_species}: {len(plot_obs)} station(s) "
                              f"with observations")
                if not _observation_data:
                    _observation_data = None
                    print("  No observations matched any plotted species.")
            else:
                print("  No stations could be matched to features.")
        print("-" * 40)
    elif (observation_dir or observation_csv) and not _primary_geojson:
        print("  WARNING: Observation data provided but no shapefile "
              "loaded — cannot match stations to features.")

    # Determine the map key for the geojson features
    _map_key = _basin_mapping_key if _basin_geojson else mapping_key

    # ── Remap sequential HRU IDs → actual shapefile IDs in traces ──
    # SUMMA HDF5 stores sequential IDs (1, 2, 3...) but the shapefile
    # has the real GRU_IDs (740457190, ...).  Build a lookup from the
    # basin shapefile feature order and rewrite trace names/hovertemplates.
    if _basin_geojson and _basin_geojson.get('features'):
        # Build mapping: "1" → first feature's GRU_ID, "2" → second, etc.
        _id_remap = {}
        for i, feat in enumerate(_basin_geojson['features']):
            seq_id = str(i + 1)
            real_id = str(feat['properties'].get(_basin_mapping_key, seq_id))
            # Strip trailing ".0" from float-converted ints
            if real_id.endswith('.0'):
                real_id = real_id[:-2]
            if seq_id != real_id:
                _id_remap[seq_id] = real_id

        if _id_remap:
            print(f"  ID remapping ({len(_id_remap)} entries): "
                  f"e.g. {list(_id_remap.items())[:3]}")
            _fl = feature_label or ''

            def _remap_fid(display_id):
                """Remap '1 (z1)' → '740457190 (z1)' or '1' → '740457190'."""
                m = re.match(r'^(.+?) \((z\d+)\)$', display_id)
                if m:
                    hru = m.group(1)
                    lyr = m.group(2)
                    return f'{_id_remap.get(hru, hru)} ({lyr})'
                return _id_remap.get(display_id, display_id)

            for p in plots:
                for t in p['traces']:
                    # Rewrite trace name
                    tname = t.get('name', '')
                    if _fl and tname.startswith(_fl + ' '):
                        old_fid = tname[len(_fl) + 1:]
                        new_fid = _remap_fid(old_fid)
                        if new_fid != old_fid:
                            t['name'] = f'{_fl} {new_fid}'
                    elif not _fl:
                        t['name'] = _remap_fid(tname)

                    # Rewrite hovertemplate
                    ht = t.get('hovertemplate', '')
                    if _fl and (_fl + ' ') in ht:
                        # Extract fid from hovertemplate
                        prefix = _fl + ' '
                        idx = ht.index(prefix) + len(prefix)
                        rest = ht[idx:]
                        br = rest.find('<br>')
                        if br >= 0:
                            old_fid = rest[:br]
                            new_fid = _remap_fid(old_fid)
                            if new_fid != old_fid:
                                t['hovertemplate'] = (ht[:idx] + new_fid +
                                                      rest[br:])

    html_content = _build_html(plots, what2map, hostmodel,
                               river_geojson=_primary_geojson,
                               map_center=_map_center,
                               map_bounds=_map_bounds,
                               mapping_key=_map_key,
                               observation_data=_observation_data,
                               feature_label=feature_label,
                               map_geom_type=_map_geom_type,
                               station_locations=station_locations,
                               station_to_feature=station_to_feature,
                               separator=separator,
                               basin_geojson=_basin_geojson,
                               river_line_geojson=_river_geojson)

    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(html_content)

    # Summary
    print(f"\n{'=' * 70}")
    print("✓ COMPLETE!")
    print("=" * 70)
    print(f"  Plots created: {len(plots)}")
    print(f"  Output: {output_path}")
    for p in plots:
        print(f"    - {p['title']}")
    print("=" * 70)
    print("  Open the HTML file in a browser for interactive charts.")

    return output_path
