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


def _build_html(plots, what2map, hostmodel, river_geojson=None,
                map_center=None, map_bounds=None, mapping_key='SegId',
                observation_data=None, feature_label=None):
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
        GeoJSON FeatureCollection for the river network.
    map_center : list or None
        [lat, lng] for the initial map center.
    map_bounds : tuple or None
        (minx, miny, maxx, maxy) for auto-fitting the map.
    mapping_key : str
        Property name in the GeoJSON features that matches feature IDs in plots.
    feature_label : str or None
        Display label for features in legends, tooltips, etc.
        If None, defaults to mapping_key.

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
            _all_fids_set.add(t['name'].replace(_label_prefix, ''))
    _all_fids_sorted = sorted(_all_fids_set)
    _fid_color = {fid: _colorway[i % len(_colorway)]
                  for i, fid in enumerate(_all_fids_sorted)}

    # Stamp explicit line colors onto every trace so Plotly matches the map
    for p in plots:
        for t in p['traces']:
            fid = t['name'].replace(_label_prefix, '')
            t['line'] = {'color': _fid_color.get(fid, '#888')}

    # --- Merge observation traces (if any) into each plot ---
    _obs_meta = {}   # {div_id: {trace_index_str: feature_id}}
    if observation_data:
        for p in plots:
            plot_obs = observation_data.get(p['id'], [])
            if not plot_obs:
                continue
            div_id = f'{p["id"]}_div'
            p_meta = {}
            base_idx = len(p['traces'])
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
                        'size': 7,
                        'symbol': 'circle',
                        'line': {'width': 1, 'color': '#333'},
                    },
                    'hovertemplate': (
                        f'Station {od["station_id"]}<br>'
                        f'%{{x}}<br>%{{y:.4g}}'
                        f'<extra>\u2192 {feature_label} {fid}</extra>'
                    ),
                    'showlegend': False,
                }
                p['traces'].append(trace)
                p_meta[str(base_idx + i)] = fid
            if p_meta:
                _obs_meta[div_id] = p_meta

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
html{scroll-behavior:smooth}
a{color:var(--primary);text-decoration:none}
.layout{display:flex;min-height:100vh}
.sidebar{width:220px;position:sticky;top:0;height:100vh;background:var(--glass);
  backdrop-filter:blur(12px);border-right:1px solid var(--border);padding:1.2rem .8rem;
  display:flex;flex-direction:column;z-index:10;overflow-y:auto;transition:background .3s,border .3s}
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
.main{flex:1;min-width:0;overflow-x:hidden}
.header{background:linear-gradient(135deg,#1a1a2e 0%,#16213e 50%,#0f3460 100%);
  color:#fff;padding:3rem 2.5rem 2rem;position:relative;overflow:hidden}
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
#mapContainer{width:100%;height:420px;border-radius:12px;z-index:1}
.map-legend{position:absolute;bottom:12px;right:12px;background:var(--glass);
  backdrop-filter:blur(8px);border:1px solid var(--border);border-radius:10px;
  padding:.6rem .8rem;font-size:.75rem;max-height:300px;overflow-y:auto;z-index:800;
  color:var(--text)}
.map-legend .ml-title{font-weight:600;margin-bottom:.3rem;font-size:.8rem}
.map-legend .ml-item{display:flex;align-items:center;gap:.4rem;padding:2px 4px;
  border-radius:4px;cursor:pointer;transition:background .15s}
.map-legend .ml-item:hover{background:rgba(0,102,204,.08)}
.map-legend .ml-item.active{background:rgba(0,102,204,.15);font-weight:600}
.map-legend .ml-swatch{width:14px;height:4px;border-radius:2px;flex-shrink:0}
</style>
</head>
<body>""")

    # --- SIDEBAR ---
    H.append('<div class="layout">')
    H.append('<aside class="sidebar">')
    H.append('<div class="logo">Open<span>WQ</span> Results</div>')
    H.append('<nav>')
    if river_geojson:
        H.append('<a href="#rivermap">River Network Map</a>')
    for p in plots:
        H.append(f'<a href="#{p["id"]}">{p["title"]}</a>')
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
    H.append('</aside>')

    # --- MAIN ---
    H.append('<div class="main">')

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

    H.append('<div class="container">')

    # --- RIVER NETWORK MAP ---
    if river_geojson:
        H.append("""<div class="section" id="rivermap">
<h2>River Network Map</h2>
<div class="card" style="padding:0;overflow:hidden;position:relative">
<div id="mapContainer"></div>
<div class="map-legend" id="mapLegend"></div>
</div>
<p class="plot-caption" style="margin-top:.4rem;font-size:.78rem;color:var(--text2)">
Click a river segment to isolate it in all plots below.</p>
</div>""")

    # --- PLOT SECTIONS ---
    for p in plots:
        div_id = f'{p["id"]}_div'
        H.append(f"""<div class="section" id="{p['id']}">
<h2>{p['title']}</h2>
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
<div class="card">
<div id="{div_id}" style="width:100%;min-height:500px"></div>
<p class="plot-caption">{p.get('caption', '')}</p>
</div>
</div>""")

    H.append('</div>')  # container
    H.append('</div>')  # main
    H.append('</div>')  # layout

    # --- JAVASCRIPT ---
    # Theme restore runs immediately (before deferred scripts)
    H.append("""<script>
(function(){
  var saved = localStorage.getItem('owq_theme');
  if(saved) document.documentElement.setAttribute('data-theme', saved);
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
""")

    # Legend click: single-click isolates a trace (others become 'legendonly',
    # i.e. grayed out in legend but still clickable); click same trace again
    # to restore all.  Also broadcasts to map via _owqSelectFeature.
    H.append("""
function _owqLegendClick(gd){
  gd.on('plotly_legendclick',function(e){
    var d=gd.data, n=d.length, ci=e.curveNumber;
    var om=_owqObsMeta[gd.id]||{};
    // Count visible non-obs feature traces
    var vis=d.map(function(t,idx){
      if(om[idx]!==undefined) return false;
      return t.visible===undefined||t.visible===true;
    });
    var onlyMe=vis.filter(Boolean).length===1 && vis[ci];
    var clickedFid=(d[ci].name||'').replace(_owqFeatureLabel,'');
    var newVis=d.map(function(t,idx){
      if(om[idx]!==undefined){
        if(onlyMe) return true;
        return om[idx]===clickedFid?true:false;
      }
      return onlyMe?true:(idx===ci?true:'legendonly');
    });
    Plotly.restyle(gd,{visible:newVis});
    // Broadcast to map
    if(typeof _owqSelectFeature==='function'){
      _owqSelectFeature(onlyMe?null:clickedFid);
    }
    return false;   // prevent Plotly default toggle
  });
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
        geojson_str = json.dumps(river_geojson)
        _fid_color_json = json.dumps(_fid_color)
        _fid_list_json = json.dumps(_all_fids_sorted)

        H.append(f"""
// --- River Network Map ---
(function(){{
  var mapEl=document.getElementById('mapContainer');
  if(!mapEl) return;
  var isDark=document.documentElement.getAttribute('data-theme')==='dark';
  var tileUrl=isDark
    ?'https://{{s}}.basemaps.cartocdn.com/dark_all/{{z}}/{{x}}/{{y}}{{r}}.png'
    :'https://{{s}}.basemaps.cartocdn.com/light_all/{{z}}/{{x}}/{{y}}{{r}}.png';
  // Start with interactions disabled (locked view)
  var map=L.map('mapContainer',{{
    zoomControl:false,dragging:false,scrollWheelZoom:false,
    doubleClickZoom:false,touchZoom:false,boxZoom:false,keyboard:false
  }});
  var tileLayer=L.tileLayer(tileUrl,{{
    attribution:'&copy; <a href="https://carto.com">CARTO</a>',maxZoom:18}}).addTo(map);

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

  var geojsonData={geojson_str};
  var mapKey='{mapping_key}';
  var fidColor={_fid_color_json};
  var plotFids={_fid_list_json};

  var featureLayers={{}};  // fid → Leaflet layer
  var selectedFid=null;

  var geoLayer=L.geoJSON(geojsonData,{{
    style:function(feature){{
      var fid=String(feature.properties[mapKey]||'');
      var c=fidColor[fid]||'#888';
      return {{color:c,weight:3,opacity:0.85}};
    }},
    onEachFeature:function(feature,layer){{
      var fid=String(feature.properties[mapKey]||'');
      featureLayers[fid]=layer;
      // Tooltip with feature name on hover
      layer.bindTooltip(_owqFeatureLabel+fid,{{sticky:true,direction:'top',opacity:0.9}});
      layer.on('click',function(){{
        if(selectedFid===fid){{
          _owqSelectFeature(null);  // deselect → restore all
        }}else{{
          _owqSelectFeature(fid);   // select this feature
        }}
      }});
      layer.on('mouseover',function(){{ layer.setStyle({{weight:5,opacity:1}}); }});
      layer.on('mouseout',function(){{
        var w=(selectedFid&&selectedFid!==fid)?1.5:3;
        var o=(selectedFid&&selectedFid!==fid)?0.3:0.85;
        layer.setStyle({{weight:w,opacity:o}});
      }});
    }}
  }}).addTo(map);

  // Fit map to data bounds
  var _fitBounds=null;
  try{{ _fitBounds=geoLayer.getBounds(); map.fitBounds(_fitBounds,{{padding:[20,20]}}); }}catch(e){{}}

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
        if(_fitBounds) map.fitBounds(_fitBounds,{{padding:[20,20]}});
      }});
      return btn;
    }}
  }});
  new RecenterCtrl().addTo(map);

  // Build map legend
  var legendEl=document.getElementById('mapLegend');
  if(legendEl){{
    var html='<div class="ml-title">{feature_label}s</div>';
    plotFids.forEach(function(fid){{
      var c=fidColor[fid]||'#888';
      html+='<div class="ml-item" data-fid="'+fid+'">'
           +'<span class="ml-swatch" style="background:'+c+'"></span>'
           +'<span>'+fid+'</span></div>';
    }});
    legendEl.innerHTML=html;
    legendEl.querySelectorAll('.ml-item').forEach(function(el){{
      el.addEventListener('click',function(){{
        var fid=el.getAttribute('data-fid');
        if(selectedFid===fid){{ _owqSelectFeature(null); }}
        else{{ _owqSelectFeature(fid); }}
      }});
    }});
  }}

  // Global selection function — debounced to avoid rapid Plotly restyles
  var _selTimer=null;
  window._owqSelectFeature=function(fid){{
    selectedFid=fid;
    // 1) Update map styles (cheap — do immediately)
    for(var f in featureLayers){{
      var ly=featureLayers[f];
      if(!fid){{ ly.setStyle({{weight:3,opacity:0.85}}); }}
      else if(f===fid){{ ly.setStyle({{weight:5,opacity:1}}); }}
      else{{ ly.setStyle({{weight:1.5,opacity:0.3}}); }}
    }}
    // 2) Update map legend highlights (cheap — do immediately)
    if(legendEl){{
      legendEl.querySelectorAll('.ml-item').forEach(function(el){{
        el.classList.toggle('active',el.getAttribute('data-fid')===fid);
      }});
    }}
    // 3) Debounce Plotly restyles (expensive)
    if(_selTimer) cancelAnimationFrame(_selTimer);
    _selTimer=requestAnimationFrame(function(){{
      var ids=window._owqPlotIds, i=0;
      function _nextRestyle(){{
        if(i>=ids.length) return;
        var gd=document.getElementById(ids[i++]);
        if(!gd||!gd.data){{ _nextRestyle(); return; }}
        var om=_owqObsMeta[gd.id]||{{}};
        var newVis=gd.data.map(function(t,tidx){{
          if(om[tidx]!==undefined){{
            if(!fid) return true;
            return om[tidx]===fid?true:false;
          }}
          var tfid=(t.name||'').replace(_owqFeatureLabel,'');
          if(!fid) return true;
          return tfid===fid?true:'legendonly';
        }});
        Plotly.restyle(gd,{{visible:newVis}}).then(_nextRestyle);
      }}
      _nextRestyle();
    }});
  }};

  // Theme toggle: swap tile layer
  var origToggle=window.toggleTheme;
  window.toggleTheme=function(){{
    origToggle();
    var d=document.documentElement.getAttribute('data-theme')==='dark';
    var url=d
      ?'https://{{s}}.basemaps.cartocdn.com/dark_all/{{z}}/{{x}}/{{y}}{{r}}.png'
      :'https://{{s}}.basemaps.cartocdn.com/light_all/{{z}}/{{x}}/{{y}}{{r}}.png';
    tileLayer.setUrl(url);
  }};
}})();
""")
    else:
        # No map — provide a no-op so legend click broadcasts don't error
        H.append("""
window._owqSelectFeature = window._owqSelectFeature || function(){};
""")

    # Sidebar active state on scroll (throttled) + click-to-scroll
    H.append("""
(function(){
  var links = document.querySelectorAll('.sidebar nav a');
  // Cache section elements once
  var sectionMap = [];
  links.forEach(function(link){
    var id = link.getAttribute('href');
    if(!id||!id.startsWith('#')) return;
    var sec = document.querySelector(id);
    if(sec) sectionMap.push({link:link, sec:sec});
  });

  // Click handler: smooth-scroll to section
  links.forEach(function(link){
    link.addEventListener('click', function(e){
      var id = link.getAttribute('href');
      if(!id||!id.startsWith('#')) return;
      e.preventDefault();
      var target = document.querySelector(id);
      if(target) target.scrollIntoView({behavior:'smooth', block:'start'});
    });
  });

  // Throttled scroll listener for active state
  var _scrollTick = false;
  window.addEventListener('scroll', function(){
    if(_scrollTick) return;
    _scrollTick = true;
    requestAnimationFrame(function(){
      var fromTop = window.scrollY + 80;
      sectionMap.forEach(function(item){
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
                   mapping_key='SegId',
                   feature_label=None,
                   observation_dir=None,
                   observation_csv=None):
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
        Kept for backward compatibility (ignored — features are always combined).
    figsize : tuple
        Kept for backward compatibility (ignored — Plotly charts are responsive).

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
                'title': f'{hostmodel.upper()} — {hydromodel_var2print}',
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

        print(f"  Chemical species: {chemSpec}")
        print(f"  Extensions to process: {file_extensions}")
        print(f"  Will create {len(chemSpec) * len(file_extensions)} plot(s)")

        # Loop through species and extensions
        for spec_idx, spec in enumerate(chemSpec):
            print(f"\n{'=' * 70}")
            print(f"Processing species {spec_idx + 1}/{len(chemSpec)}: {spec}")
            print("=" * 70)

            try:
                # Find matching species in openwq_results
                species_key = None
                for key in openwq_results.keys():
                    if spec in key:
                        species_key = key
                        print(f"  Found species key: {species_key}")
                        break

                if species_key is None:
                    print(f"  ✗ No data found for species: {spec}")
                    print(f"  Available species:")
                    for key in list(openwq_results.keys())[:5]:
                        print(f"    - {key}")
                    continue

                # Parse units
                parts = species_key.split('#')
                if len(parts) > 1:
                    units = parts[1]
                else:
                    units = 'kg' if 'Sediment' in species_key else 'concentration'

                # Get chemical name
                chem_parts = species_key.split('@')
                chem_name = chem_parts[1].split('#')[0] if len(chem_parts) > 1 else spec

                # Loop through extensions
                for ext_idx, file_extension in enumerate(file_extensions):
                    print(f"\n  {'- ' * 33}")
                    print(f"  Processing extension {ext_idx + 1}/{len(file_extensions)}: {file_extension}")
                    print(f"  {'- ' * 33}")

                    # Extract data
                    data_found = False
                    ttdata = None
                    for ext, data_list in openwq_results[species_key]:
                        if ext == file_extension:
                            if data_list and len(data_list) > 0:
                                filename, ttdata, xyz_coords = data_list[0]
                                data_found = True
                                print(f"    Data shape: {ttdata.shape}")
                                print(f"    Time range: {ttdata.index[0]} to {ttdata.index[-1]}")
                                break

                    if not data_found:
                        print(f"    ✗ No data found for extension '{file_extension}'")
                        continue

                    # Extract feature data
                    feature_data = {}
                    missing_features = []

                    if mapping_key_values == ["all"] or mapping_key_values == "all":
                        for col in ttdata.columns:
                            fid = col.split('_')[-1]
                            feature_data[fid] = ttdata[col]
                            print(f"    ✓ Found data for feature {fid}")
                    else:
                        for fid in mapping_key_values:
                            matching_cols = [col for col in ttdata.columns
                                             if str(fid) in col.split('_')[-1]]
                            if len(matching_cols) > 0:
                                feature_data[fid] = ttdata[matching_cols[0]]
                                print(f"    ✓ Found data for feature {fid}")
                            else:
                                missing_features.append(fid)
                                print(f"    ✗ Feature {fid} not found")

                    if len(feature_data) == 0:
                        print(f"    ✗ No data found for any requested features")
                        continue

                    if missing_features:
                        print(f"    ⚠ Missing features: {missing_features}")

                    # Build plot entry
                    traces = _build_traces(feature_data, feature_label=feature_label)
                    plot_id = f'plot_{_sanitize_id(spec)}_{_sanitize_id(file_extension)}'
                    time_range = (f'{ttdata.index[0]} to {ttdata.index[-1]}'
                                  if len(ttdata.index) > 0 else '')
                    plots.append({
                        'id': plot_id,
                        'title': f'{chem_name} ({file_extension})',
                        'traces': traces,
                        'ylabel': f'{chem_name} ({units})',
                        'xtype': ('date' if isinstance(ttdata.index, pd.DatetimeIndex)
                                  else 'linear'),
                        'caption': f'{len(feature_data)} {feature_label}s | {time_range}',
                        'species': spec,
                    })

                    print(f"    ✓ Plot data collected for {chem_name} ({file_extension})")

            except Exception as e:
                print(f"  ✗ Error processing species {spec}: {e}")
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

    # Load river network shapefile for interactive map (optional)
    _river_geojson = None
    _map_center = None
    _map_bounds = None
    if river_network_shp and os.path.isfile(river_network_shp):
        try:
            import fiona
            from shapely.geometry import shape, mapping as shp_mapping

            features = []
            src_proj4 = None
            with fiona.open(river_network_shp) as src:
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
                # Check if reprojection needed
                sb = features[0]['geometry'].bounds
                if (abs(sb[0]) > 360 or abs(sb[1]) > 360) and src_proj4:
                    print("  Reprojecting shapefile to WGS84...")
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
                _river_geojson = {'type': 'FeatureCollection',
                                  'features': gj_features}
                all_b = [f['geometry'].bounds for f in features]
                _map_bounds = (min(b[0] for b in all_b), min(b[1] for b in all_b),
                               max(b[2] for b in all_b), max(b[3] for b in all_b))
                _map_center = [(_map_bounds[1]+_map_bounds[3])/2,
                               (_map_bounds[0]+_map_bounds[2])/2]
                print(f"  ✓ Loaded river network: {len(gj_features)} features")
        except ImportError:
            print("  WARNING: fiona/shapely not installed — skipping map")
        except Exception as _e:
            print(f"  WARNING: Failed to load shapefile: {_e}")

    # --- Load observation data and match to features ---
    _observation_data = None    # {plot_id: [obs_dict, ...]}
    if (observation_dir or observation_csv) and _river_geojson:
        print("\n" + "-" * 40)
        print("Loading observation data...")
        obs_by_species, station_locations = _load_observation_data(
            obs_dir=observation_dir, obs_csv=observation_csv)
        if obs_by_species and station_locations:
            station_to_feature = _match_stations_to_features(
                station_locations, _river_geojson, mapping_key)
            if station_to_feature:
                print(f"  Matched {len(station_to_feature)} stations to features:")
                for sid, fid in station_to_feature.items():
                    print(f"    Station {sid} → Feature {fid}")
                _observation_data = {}
                for p in plots:
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
    elif (observation_dir or observation_csv) and not _river_geojson:
        print("  WARNING: Observation data provided but no river network "
              "shapefile loaded — cannot match stations to features.")

    html_content = _build_html(plots, what2map, hostmodel,
                               river_geojson=_river_geojson,
                               map_center=_map_center,
                               map_bounds=_map_bounds,
                               mapping_key=mapping_key,
                               observation_data=_observation_data,
                               feature_label=feature_label)

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
