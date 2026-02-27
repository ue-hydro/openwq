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
Gen_Report.py — Self-contained HTML report for OpenWQ simulation results.

Generates an interactive HTML report after model execution, following the
design patterns from the calibration postprocessing module (results_analysis.py):
  - Dark/light theme toggle with localStorage persistence
  - Plotly.js interactive charts (CDN)
  - Leaflet.js interactive basin map (CDN)
  - Responsive sidebar navigation
"""

import os
import sys
import json
import datetime
import glob as _glob
import numpy as np

# HDF5 reader — add 2_Read_Outputs to path
_this_dir = os.path.dirname(os.path.abspath(__file__))
_read_outputs_dir = os.path.normpath(
    os.path.join(_this_dir, '..', '..', '2_Read_Outputs', 'hdf5_support_lib'))
if _read_outputs_dir not in sys.path:
    sys.path.insert(0, _read_outputs_dir)

# GRQA extraction tools — add 3_Calibration path
_calibration_dir = os.path.normpath(
    os.path.join(_this_dir, '..', '..', '3_Calibration'))
if _calibration_dir not in sys.path:
    sys.path.insert(0, _calibration_dir)

try:
    from calibration_lib.observation_data.grqa_extract_stations import (
        GRQACalibrationExtractor, SpeciesMapper, GRQA_PARAMETERS)
    _GRQA_AVAILABLE = True
except ImportError:
    _GRQA_AVAILABLE = False

# Auto-mapping: model species name → GRQA parameter code
_MODEL_SPECIES_TO_GRQA = {
    'NO3-N': 'NO3', 'NO3': 'NO3', 'Nitrate': 'NO3',
    'NH4-N': 'NH4', 'NH4': 'NH4', 'Ammonium': 'NH4',
    'NO2-N': 'NO2', 'NO2': 'NO2',
    'TN': 'TN', 'Total Nitrogen': 'TN',
    'DON': 'DON', 'DIN': 'DIN', 'TKN': 'TKN',
    'PO4-P': 'PO4', 'PO4': 'PO4', 'Phosphate': 'PO4',
    'TP': 'TP', 'Total Phosphorus': 'TP', 'DP': 'DP',
    'DO': 'DO', 'Dissolved Oxygen': 'DO',
    'BOD': 'BOD', 'BOD5': 'BOD5', 'COD': 'COD',
    'DOC': 'DOC', 'TOC': 'TOC', 'DIC': 'DIC',
    'TSS': 'TSS', 'TDS': 'TDS',
    'N_ORG_active': 'DON', 'N_ORG_fresh': 'DON', 'N_ORG_stable': 'TKN',
}


def _js(obj):
    """Convert Python object to safe JSON string for embedding in HTML/JS."""
    if isinstance(obj, (np.integer,)):
        return json.dumps(int(obj))
    if isinstance(obj, (np.floating,)):
        if np.isnan(obj) or np.isinf(obj):
            return 'null'
        return json.dumps(float(obj))
    if isinstance(obj, np.ndarray):
        return json.dumps(obj.tolist())
    return json.dumps(obj)


def _read_module_json(output_dir, module_name):
    """Read a module JSON file from openwq_in/, stripping comment lines."""
    if module_name.upper() == 'NONE':
        return None
    fpath = os.path.join(output_dir, 'openwq_in',
                         f'openWQ_MODULE_{module_name}.json')
    if not os.path.isfile(fpath):
        return None
    try:
        with open(fpath, 'r') as f:
            lines = [l for l in f.readlines()
                     if not l.strip().startswith('//')]
        return json.loads(''.join(lines))
    except Exception:
        return None


def _render_module_params(module_data, module_name):
    """Generate HTML list for module parameter display inside a <details> block."""
    if module_data is None or module_name.upper() == 'NONE':
        return []

    H = []
    mn = module_name.upper()

    # --- Transport Dissolved ---
    if 'TRANSPORT_CONFIGURATION' in module_data:
        tc = module_data['TRANSPORT_CONFIGURATION']
        H.append('<table class="param-table">')
        H.append('<tr><th>Parameter</th><th>Value</th></tr>')
        for k, v in tc.items():
            H.append(f'<tr><td>{k}</td><td class="num">{v}</td></tr>')
        H.append('</table>')

    # --- Lateral Exchange ---
    elif 'CONFIGURATION' in module_data and 'LE' in mn:
        cfg = module_data['CONFIGURATION']
        H.append('<table class="param-table">')
        H.append('<tr><th>#</th><th>Direction</th><th>Upper Compartment</th>'
                 '<th>Lower Compartment</th><th>K<sub>val</sub></th></tr>')
        for idx, entry in cfg.items():
            if isinstance(entry, dict):
                H.append(
                    f'<tr><td>{idx}</td>'
                    f'<td>{entry.get("direction", "")}</td>'
                    f'<td>{entry.get("upper_compartment", "")}</td>'
                    f'<td>{entry.get("lower_compartment", "")}</td>'
                    f'<td class="num">{entry.get("K_val", "")}</td></tr>')
        H.append('</table>')

    # --- Sediment Transport (HYPE_MMF / HYPE_HBVSED) ---
    elif mn.startswith('HYPE'):
        if 'CONFIGURATION' in module_data:
            cfg = module_data['CONFIGURATION']
            H.append('<h4>Configuration</h4>')
            H.append('<table class="param-table">')
            H.append('<tr><th>Setting</th><th>Value</th></tr>')
            for k, v in cfg.items():
                H.append(f'<tr><td>{k}</td><td>{v}</td></tr>')
            H.append('</table>')

        if 'PARAMETER_DEFAULTS' in module_data:
            pd = module_data['PARAMETER_DEFAULTS']
            H.append('<h4>Parameter Defaults</h4>')
            H.append('<table class="param-table">')
            H.append('<tr><th>Parameter</th><th>Value</th></tr>')
            for k, v in pd.items():
                H.append(f'<tr><td>{k}</td><td class="num">{v}</td></tr>')
            H.append('</table>')

        if 'MONTHLY_EROSION_FACTOR' in module_data:
            mef = module_data['MONTHLY_EROSION_FACTOR']
            months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
                      'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
            H.append('<h4>Monthly Erosion Factor</h4>')
            H.append('<table class="param-table">')
            H.append('<tr>' + ''.join(f'<th>{m}</th>' for m in months) + '</tr>')
            H.append('<tr>' + ''.join(f'<td class="num">{v}</td>' for v in mef) + '</tr>')
            H.append('</table>')

    # --- BGC (NATIVE_BGC_FLEX) ---
    elif 'BIOGEOCHEMISTRY_CONFIGURATION' in module_data:
        bgc = module_data['BIOGEOCHEMISTRY_CONFIGURATION']

        if 'CHEMICAL_SPECIES' in bgc:
            species = bgc['CHEMICAL_SPECIES']
            H.append('<h4>Chemical Species</h4>')
            H.append('<table class="param-table">')
            H.append('<tr><th>Species</th><th>Description</th>'
                     '<th>Units</th><th>Initial Condition</th></tr>')
            for k, v in species.items():
                if k.startswith('_') or not isinstance(v, dict):
                    continue
                H.append(
                    f'<tr><td><code>{k}</code></td>'
                    f'<td>{v.get("DESCRIPTION", "")}</td>'
                    f'<td>{v.get("UNITS", "")}</td>'
                    f'<td class="num">{v.get("INITIAL_CONDITION", "")}</td></tr>')
            H.append('</table>')

        if 'CYCLING_FRAMEWORKS' in bgc:
            frameworks = bgc['CYCLING_FRAMEWORKS']
            H.append('<h4>Cycling Frameworks</h4>')
            for fw_name, fw_data in frameworks.items():
                if not isinstance(fw_data, dict):
                    continue
                desc = fw_data.get('_DESCRIPTION', fw_name)
                H.append(f'<details class="nested-details">'
                         f'<summary>{fw_name} &mdash; {desc}</summary>')
                H.append('<table class="param-table">')
                H.append('<tr><th>Transformation</th><th>Consumed</th>'
                         '<th>Produced</th><th>Kinetics</th></tr>')
                for t_name, t_data in fw_data.items():
                    if t_name.startswith('_') or not isinstance(t_data, dict):
                        continue
                    if 'KINETICS' not in t_data and 'FORMULA' not in t_data:
                        continue
                    consumed = t_data.get('CONSUMED', '')
                    produced = t_data.get('PRODUCED', '')
                    kinetics = t_data.get('KINETICS', t_data.get('FORMULA', ''))
                    H.append(
                        f'<tr><td><code>{t_name}</code></td>'
                        f'<td>{consumed}</td>'
                        f'<td>{produced}</td>'
                        f'<td><code style="font-size:.75rem">{kinetics}</code></td></tr>')
                H.append('</table></details>')

        if 'COMPARTMENT_ASSIGNMENT' in bgc:
            ca = bgc['COMPARTMENT_ASSIGNMENT']
            H.append('<h4>Compartment Assignment</h4>')
            H.append('<table class="param-table">')
            H.append('<tr><th>Compartment</th><th>Frameworks</th></tr>')
            for comp, fws in ca.items():
                fws_str = ", ".join(fws) if isinstance(fws, list) else str(fws)
                H.append(f'<tr><td>{comp}</td><td>{fws_str}</td></tr>')
            H.append('</table>')

    # --- PHREEQC ---
    elif mn == 'PHREEQC':
        if 'PHREEQC_CONFIGURATION' in module_data:
            pc = module_data['PHREEQC_CONFIGURATION']
            H.append('<table class="param-table">')
            H.append('<tr><th>Setting</th><th>Value</th></tr>')
            for k, v in pc.items():
                val = ", ".join(v) if isinstance(v, list) else str(v)
                H.append(f'<tr><td>{k}</td><td>{val}</td></tr>')
            H.append('</table>')

    # --- Sorption Isotherm (Freundlich / Langmuir) ---
    elif mn in ('FREUNDLICH', 'LANGMUIR'):
        for key, val in module_data.items():
            if key == 'MODULE_NAME':
                continue
            if isinstance(val, dict):
                H.append(f'<h4>{key}</h4>')
                H.append('<table class="param-table">')
                H.append('<tr><th>Parameter</th><th>Value</th></tr>')
                for pk, pv in val.items():
                    if isinstance(pv, dict):
                        for spk, spv in pv.items():
                            H.append(f'<tr><td>{pk}.{spk}</td>'
                                     f'<td class="num">{spv}</td></tr>')
                    else:
                        H.append(f'<tr><td>{pk}</td>'
                                 f'<td class="num">{pv}</td></tr>')
                H.append('</table>')

    # --- Generic fallback ---
    elif module_data:
        for key, val in module_data.items():
            if key == 'MODULE_NAME':
                continue
            if isinstance(val, dict):
                H.append(f'<h4>{key}</h4>')
                H.append('<table class="param-table">')
                H.append('<tr><th>Parameter</th><th>Value</th></tr>')
                for pk, pv in val.items():
                    H.append(f'<tr><td>{pk}</td><td class="num">{pv}</td></tr>')
                H.append('</table>')
            elif isinstance(val, list):
                H.append(f'<p><strong>{key}:</strong> '
                         f'[{", ".join(str(x) for x in val)}]</p>')
            else:
                H.append(f'<p><strong>{key}:</strong> {val}</p>')

    return H


def _read_hdf5_outputs(output_dir, compartments_and_cells, chemical_species, units):
    """Read HDF5 simulation outputs.

    OpenWQ HDF5 format: each file has timestamp-named datasets (e.g.
    '1993JAN01-01:00:00') with shape (1, N_cells), plus metadata datasets
    'reachID', 'xyz_elements', 'xyz_elements_size'.

    Returns dict: {species_name: pandas.DataFrame} where DataFrame has
    DatetimeIndex and columns per reach/cell.
    """
    h5_dir = os.path.join(output_dir, 'openwq_out', 'HDF5')
    if not os.path.isdir(h5_dir):
        print(f"  WARNING: HDF5 directory not found: {h5_dir}")
        return {}

    h5_files = [f for f in os.listdir(h5_dir) if f.endswith('.h5')]
    if not h5_files:
        print(f"  WARNING: No .h5 files in {h5_dir}")
        return {}

    # Metadata keys to skip when collecting timestamp datasets
    _META_KEYS = {'reachID', 'hruId', 'xyz_elements', 'xyz_elements_size',
                  'mapping_key', 'coordinates', 'timestamps', 'data', 'concentration'}

    compartments = list(compartments_and_cells.keys())
    units_clean = units.replace('/', '|')

    results = {}
    try:
        import h5py
        import pandas as pd

        for species in chemical_species:
            species_upper = species.upper().replace('-', '_')
            for cmp in compartments:
                # Filename pattern: COMPARTMENT@SPECIES#UNITS-main.h5
                fname = f"{cmp}@{species_upper}#{units_clean}-main.h5"
                fpath = os.path.join(h5_dir, fname)

                if not os.path.isfile(fpath):
                    # Try alternate naming
                    candidates = _glob.glob(os.path.join(h5_dir, f"*{species_upper}*main.h5"))
                    if candidates:
                        fpath = candidates[0]
                        fname = os.path.basename(fpath)
                    else:
                        continue

                try:
                    with h5py.File(fpath, 'r') as f:
                        all_keys = list(f.keys())

                        # Get reach/cell IDs for column names
                        if 'reachID' in f:
                            ids_raw = f['reachID'][:]
                            ids = [x.decode() if isinstance(x, bytes) else str(x) for x in ids_raw]
                            col_names = [f"reach_{rid}" for rid in ids]
                        elif 'hruId' in f:
                            ids_raw = f['hruId'][:]
                            ids = [x.decode() if isinstance(x, bytes) else str(x) for x in ids_raw]
                            col_names = [f"hru_{hid}" for hid in ids]
                        else:
                            col_names = None

                        # Collect timestamp datasets (keys not in metadata)
                        ts_keys = [k for k in all_keys if k not in _META_KEYS]

                        # Parse timestamps from key names (format: YYYYMmmDD-HH:MM:SS)
                        rows = []
                        timestamps = []
                        for ts_key in ts_keys:
                            try:
                                # Parse format like '1993JAN01-01:00:00'
                                dt = pd.to_datetime(ts_key, format='%Y%b%d-%H:%M:%S')
                            except (ValueError, TypeError):
                                try:
                                    dt = pd.to_datetime(ts_key, errors='coerce')
                                except Exception:
                                    continue
                            if pd.isna(dt):
                                continue

                            ds = f[ts_key][:]
                            # Squeeze from (1, N) to (N,)
                            if ds.ndim > 1:
                                ds = ds.squeeze(axis=0)
                            rows.append(ds)
                            timestamps.append(dt)

                        if not rows:
                            continue

                        data = np.array(rows)
                        if col_names is None:
                            col_names = [f"cell_{j}" for j in range(data.shape[1])]
                        elif len(col_names) != data.shape[1]:
                            col_names = [f"cell_{j}" for j in range(data.shape[1])]

                        df = pd.DataFrame(data, index=pd.DatetimeIndex(timestamps),
                                          columns=col_names)
                        df.sort_index(inplace=True)

                        # Replace no-data flags
                        df.replace(-9999, np.nan, inplace=True)
                        df.replace(-9999.0, np.nan, inplace=True)

                        results[species] = df
                        print(f"  Loaded: {fname} ({len(df)} timesteps, {len(df.columns)} cells)")

                except Exception as e:
                    print(f"  WARNING: Failed to read {fpath}: {e}")

    except ImportError as e:
        print(f"  WARNING: Cannot read HDF5 — missing dependency: {e}")

    return results


def _reproject_geom_to_wgs84(geom, src_proj4):
    """Reproject a shapely geometry from source CRS (proj4 string) to WGS84.

    Uses pyproj.Proj with inverse=True — robust against PROJ database issues.
    Returns reprojected geometry, or original if reprojection fails.
    """
    try:
        from pyproj import Proj
        from shapely.ops import transform as shp_transform
        src_proj = Proj(src_proj4)

        def _to_wgs84(x, y):
            return src_proj(x, y, inverse=True)

        return shp_transform(_to_wgs84, geom)
    except Exception:
        return geom


def _load_single_geojson(path):
    """Load a single shapefile/GeoPackage as a GeoJSON dict in WGS84.

    Attempts reprojection to WGS84 (EPSG:4326) so Leaflet can render it.
    Uses fiona + shapely with pyproj for CRS reprojection.
    Returns (geojson_dict, bounds_wgs84) or (None, None) on failure.
    """
    if path is None or not os.path.isfile(path):
        return None, None

    try:
        import fiona
        from shapely.geometry import shape, mapping

        features = []
        src_proj4 = None

        with fiona.open(path) as src:
            # Get source CRS as proj4 string for reprojection
            if src.crs:
                if isinstance(src.crs, dict):
                    parts = []
                    for k, v in src.crs.items():
                        if v is True:
                            parts.append(f'+{k}')
                        else:
                            parts.append(f'+{k}={v}')
                    src_proj4 = ' '.join(parts)
                else:
                    try:
                        src_proj4 = src.crs.to_proj4()
                    except Exception:
                        src_proj4 = None

            for feat in src:
                geom = shape(feat['geometry'])
                properties = dict(feat['properties'])
                features.append({
                    'type': 'Feature',
                    'geometry': geom,
                    'properties': properties,
                })

        if not features:
            print(f"  WARNING: Shapefile has no features: {path}")
            return None, None

        # Check if reprojection is needed (coordinates > 360 ⇒ projected)
        sample_bounds = features[0]['geometry'].bounds
        needs_reproj = (abs(sample_bounds[0]) > 360 or abs(sample_bounds[1]) > 360)

        if needs_reproj and src_proj4:
            print(f"  Reprojecting to WGS84...")
            for f in features:
                f['geometry'] = _reproject_geom_to_wgs84(f['geometry'], src_proj4)

        # Convert to GeoJSON dicts
        geojson_features = []
        for f in features:
            geojson_features.append({
                'type': 'Feature',
                'geometry': mapping(f['geometry']),
                'properties': {k: (str(v) if v is not None else '')
                               for k, v in f['properties'].items()},
            })

        geojson = {'type': 'FeatureCollection', 'features': geojson_features}
        all_bounds = [f['geometry'].bounds for f in features]
        bounds = (
            min(b[0] for b in all_bounds),
            min(b[1] for b in all_bounds),
            max(b[2] for b in all_bounds),
            max(b[3] for b in all_bounds),
        )
        return geojson, bounds

    except ImportError:
        print("  WARNING: fiona/shapely not installed — skipping shapefile")
        return None, None
    except Exception as e:
        print(f"  WARNING: Failed to load shapefile {path}: {e}")
        return None, None


def _build_map_layers(river_network_shapefile, basin_shapefile):
    """Build map layers from the river network and basin shapefiles.

    Returns:
        (layers, map_center) where layers is a list of dicts with
        'geojson_str', 'label', 'style', 'type', and map_center is [lat, lng].
        Returns ([], None) if no layers loaded.
    """
    layers = []
    all_bounds = []

    # Layer 1: Basin polygons (from Copernicus LULC shapefile)
    if basin_shapefile:
        geojson, bounds = _load_single_geojson(basin_shapefile)
        if geojson is not None:
            layers.append({
                'geojson_str': json.dumps(geojson),
                'label': 'Basin',
                'style': {'color': '#ffffff', 'weight': 1.5, 'opacity': 0.9,
                          'fillOpacity': 0.08},
                'type': 'polygon',
            })
            all_bounds.append(bounds)

    # Layer 2: River network
    if river_network_shapefile:
        geojson, bounds = _load_single_geojson(river_network_shapefile)
        if geojson is not None:
            layers.append({
                'geojson_str': json.dumps(geojson),
                'label': 'River Network',
                'style': {'color': '#0066cc', 'weight': 2.5, 'opacity': 0.85},
                'type': 'line',
            })
            all_bounds.append(bounds)

    if not layers:
        return [], None

    minx = min(b[0] for b in all_bounds)
    miny = min(b[1] for b in all_bounds)
    maxx = max(b[2] for b in all_bounds)
    maxy = max(b[3] for b in all_bounds)
    map_center = [(miny + maxy) / 2, (minx + maxx) / 2]

    return layers, map_center


def _extract_grqa_for_report(river_network_shapefile, chemical_species,
                              output_dir, grqa_local_data_path=None):
    """Extract GRQA observation stations for the simulated species.

    Uses the GRQACalibrationExtractor to find monitoring stations near the
    river network that have observations for the species being simulated.

    Parameters:
        river_network_shapefile: path to river network shapefile
        chemical_species: list of model species names (e.g. ['NO3-N', 'NH4-N'])
        output_dir: directory to cache GRQA extraction results
        grqa_local_data_path: path to pre-downloaded GRQA data (or None)

    Returns:
        (stations_geojson_str, grqa_stats) or (None, None) on failure.
        grqa_stats is a dict with keys: n_stations, n_observations,
        species_stats (list), year_range, etc.
    """
    if not _GRQA_AVAILABLE:
        print("  GRQA tools not available (missing geopandas/shapely). Skipping.")
        return None, None

    if not river_network_shapefile or not os.path.isfile(river_network_shapefile):
        print("  No river network shapefile — skipping GRQA extraction.")
        return None, None

    # Build species mapping: GRQA code → model species name
    species_mapping = {}
    for model_name in chemical_species:
        grqa_code = _MODEL_SPECIES_TO_GRQA.get(model_name)
        if grqa_code and grqa_code not in species_mapping:
            species_mapping[grqa_code] = model_name

    if not species_mapping:
        print("  No GRQA-compatible species found. Skipping GRQA extraction.")
        return None, None

    grqa_params = list(species_mapping.keys())
    print(f"  GRQA species mapping: {species_mapping}")
    print(f"  Will extract {len(grqa_params)} GRQA parameters: {', '.join(grqa_params)}")
    if grqa_local_data_path is None:
        print("  No local GRQA data provided — will download from Zenodo (one-time, ~1.2 GB).")
        print("  Tip: set grqa_local_data_path to skip future downloads.")

    grqa_output_dir = os.path.join(output_dir, 'openwq_in', 'grqa_report_data')

    try:
        mapper = SpeciesMapper(mapping=species_mapping)
        extractor = GRQACalibrationExtractor(
            output_dir=grqa_output_dir,
            species_mapper=mapper,
            local_data_path=grqa_local_data_path,
            buffer_distance_m=1000,
        )

        # Load river network — CRS reprojection to WGS84 is handled
        # automatically by the extractor's load_river_network() method.
        river_gdf = extractor.load_river_network(river_network_shapefile)

        buffer_gdf = extractor.create_buffer(river_gdf)
        stations, observations = extractor.extract_stations_and_observations(buffer_gdf)

        if stations is None or stations.empty:
            print("  No GRQA stations found near the river network.")
            return None, None

        # Generate files (shapefile + CSV)
        extractor.generate_calibration_files(prefix="report_grqa")

        # Build GeoJSON for map
        from shapely.geometry import mapping
        features = []
        for _, row in stations.iterrows():
            props = {k: (str(v) if v is not None else '')
                     for k, v in row.items() if k != 'geometry'}
            features.append({
                'type': 'Feature',
                'geometry': mapping(row.geometry),
                'properties': props,
            })
        stations_geojson = {'type': 'FeatureCollection', 'features': features}

        # Compute statistics
        import pandas as pd
        site_col = extractor.column_map.get('site_id', 'site_id')
        date_col = extractor.column_map.get('obs_date', 'obs_date')

        dates = pd.to_datetime(observations[date_col], errors='coerce')
        species_stats = []
        for model_sp in observations['model_species'].unique():
            sp_obs = observations[observations['model_species'] == model_sp]
            sp_dates = pd.to_datetime(sp_obs[date_col], errors='coerce')
            sp_stations = sp_obs[site_col].nunique()
            species_stats.append({
                'species': model_sp,
                'n_stations': sp_stations,
                'n_observations': len(sp_obs),
                'year_start': int(sp_dates.min().year) if not sp_dates.isna().all() else None,
                'year_end': int(sp_dates.max().year) if not sp_dates.isna().all() else None,
            })

        grqa_stats = {
            'n_stations': len(stations),
            'n_observations': len(observations),
            'year_start': int(dates.min().year) if not dates.isna().all() else None,
            'year_end': int(dates.max().year) if not dates.isna().all() else None,
            'species_stats': species_stats,
            'output_dir': grqa_output_dir,
        }

        return json.dumps(stations_geojson), grqa_stats

    except Exception as e:
        print(f"  WARNING: GRQA extraction failed: {e}")
        import traceback
        traceback.print_exc()
        return None, None


def _compute_spatial_stats(results):
    """Compute spatial statistics (min/max/mean/std) per species across all cells."""
    stats = []
    for species, df in results.items():
        numeric = df.select_dtypes(include=[np.number])
        if numeric.empty:
            continue
        # Statistics across all cells and timesteps
        all_vals = numeric.values.flatten()
        all_vals = all_vals[~np.isnan(all_vals)]
        if len(all_vals) == 0:
            continue
        stats.append({
            'species': species,
            'min': float(np.min(all_vals)),
            'max': float(np.max(all_vals)),
            'mean': float(np.mean(all_vals)),
            'std': float(np.std(all_vals)),
            'median': float(np.median(all_vals)),
            'n_cells': numeric.shape[1],
            'n_timesteps': numeric.shape[0],
        })
    return stats


def generate_simulation_report(
        output_dir,
        project_name,
        authors,
        date,
        comment,
        hostmodel,
        bgc_module_name,
        td_module_name,
        le_module_name,
        ts_module_name,
        si_module_name,
        ss_method,
        solver,
        chemical_species,
        units,
        compartments_and_cells,
        timestep,
        river_network_shapefile=None,
        basin_shapefile=None,
        grqa_local_data_path=None,
        container_runtime="docker",
        mpi_np=2,
        run_time_seconds=0.0,
        model_was_run=True,
):
    """Generate a self-contained HTML simulation report.

    Parameters:
        river_network_shapefile: path to river network shapefile (.shp/.gpkg)
        basin_shapefile: path to basin/catchment shapefile (auto-loaded from
            ss_method_copernicus_basin_info if available)
        grqa_local_data_path: "auto" or None to download from Zenodo,
            a path to pre-downloaded GRQA folder, or "skip" to disable.
            When enabled and river_network_shapefile is set, GRQA observation
            stations are shown on the map and a data availability report is generated.

    Returns the path to the generated HTML file.
    """
    report_path = os.path.join(output_dir, "openwq_report.html")
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")

    results = {}
    if model_was_run:
        print(f"  Reading HDF5 outputs...")
        results = _read_hdf5_outputs(output_dir, compartments_and_cells, chemical_species, units)

    # Spatial statistics
    spatial_stats = _compute_spatial_stats(results)

    # Basin map layers (river network + basin polygons)
    map_layers, map_center = _build_map_layers(river_network_shapefile, basin_shapefile)

    # GRQA observation stations
    grqa_geojson_str = None
    grqa_stats = None
    _grqa_skip = (isinstance(grqa_local_data_path, str)
                  and grqa_local_data_path.strip().lower() == "skip")
    if river_network_shapefile and not _grqa_skip:
        print("  Extracting GRQA observation stations...")
        # "auto" or None → download from Zenodo; anything else → treat as path
        _grqa_path = None
        if isinstance(grqa_local_data_path, str) and \
                grqa_local_data_path.strip().lower() not in ("auto", ""):
            _grqa_path = grqa_local_data_path
        grqa_geojson_str, grqa_stats = _extract_grqa_for_report(
            river_network_shapefile, chemical_species,
            output_dir, _grqa_path)
    elif _grqa_skip:
        print("  GRQA: skipped by user (grqa_local_data_path = 'skip').")

    # Read source/sink config for summary
    ss_config_path = os.path.join(output_dir, 'openwq_in',
                                  f'openWQ_SS_{ss_method}.json')
    ss_summary = []
    ss_metadata = {}
    if os.path.isfile(ss_config_path):
        try:
            with open(ss_config_path, 'r') as f:
                lines = [l for l in f.readlines()
                         if not l.strip().startswith('//')]
                ss_data = json.loads(''.join(lines))

            # Extract metadata block
            if 'METADATA' in ss_data:
                ss_metadata = ss_data['METADATA']

            # Entries are stored under numeric keys ("1", "2", ...)
            # or under a "SINKSOURCE" wrapper (legacy format)
            entries_dict = ss_data
            if 'SINKSOURCE' in ss_data:
                entries_dict = ss_data['SINKSOURCE']

            for key, entry in entries_dict.items():
                if key == 'METADATA' or not isinstance(entry, dict):
                    continue
                if 'CHEMICAL_NAME' not in entry:
                    continue

                # Parse DATA block for statistics
                data_block = entry.get('DATA', {})
                n_entries = len(data_block)
                cells = set()
                years = set()
                total_load = 0.0
                for _dk, row in data_block.items():
                    if isinstance(row, list) and len(row) >= 10:
                        years.add(row[0])
                        cells.add(str(row[6]))  # cell_id or ix
                        try:
                            total_load += float(row[9])
                        except (ValueError, TypeError):
                            pass

                ss_summary.append({
                    'chemical': entry.get('CHEMICAL_NAME', 'N/A'),
                    'compartment': entry.get('COMPARTMENT_NAME', 'N/A'),
                    'type': entry.get('TYPE', 'N/A'),
                    'units': entry.get('UNITS', 'N/A'),
                    'data_format': entry.get('DATA_FORMAT', ''),
                    'comment': entry.get('COMMENT', ''),
                    'n_entries': n_entries,
                    'n_cells': len(cells),
                    'year_range': (min(years), max(years)) if years else None,
                    'total_load': total_load,
                })
        except Exception as e:
            print(f"  WARNING: Could not parse SS config: {e}")

    # Format runtime
    if not model_was_run:
        runtime_str = "Not executed"
    elif run_time_seconds >= 3600:
        runtime_str = f"{run_time_seconds/3600:.1f} hours"
    elif run_time_seconds >= 60:
        runtime_str = f"{run_time_seconds/60:.1f} min"
    else:
        runtime_str = f"{run_time_seconds:.1f}s"

    # =========================================================================
    # Build HTML
    # =========================================================================
    H = []

    # --- HEAD ---
    H.append("""<!DOCTYPE html>
<html lang="en" data-theme="dark">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>OpenWQ Simulation Report</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link href="https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&family=JetBrains+Mono:wght@400;500&display=swap" rel="stylesheet">
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
<link rel="stylesheet" href="https://unpkg.com/leaflet@1.9.4/dist/leaflet.css">
<script src="https://unpkg.com/leaflet@1.9.4/dist/leaflet.js"></script>
<style>
:root,[data-theme="light"]{
  --primary:#0066cc;--primary-dark:#004499;
  --secondary:#00a86b;--accent:#ff6b35;
  --dark:#1a1a2e;--light:#f8f9fa;
  --bg:#f0f2f5;--surface:#fff;--text:#2d3436;--text2:#636e72;--text3:#a0aec0;
  --border:#e2e8f0;--border2:#d0d4da;
  --shadow:0 4px 20px rgba(0,0,0,.06);--shadow-lg:0 10px 40px rgba(0,0,0,.08);
  --plot-bg:#fff;--plot-grid:#f0f0f0;--plot-font:#333;
  --glass:rgba(255,255,255,.85);
}
[data-theme="dark"]{
  --primary:#4d9ee8;--primary-dark:#3a8bd4;
  --secondary:#34d399;--accent:#fb923c;
  --dark:#1a1a2e;--light:#1e1f2e;
  --bg:#0f1117;--surface:#1a1b2e;--text:#e2e8f0;--text2:#a0aec0;--text3:#636e72;
  --border:#2a2b3d;--border2:#3a3b4d;
  --shadow:0 4px 20px rgba(0,0,0,.3);--shadow-lg:0 10px 40px rgba(0,0,0,.4);
  --plot-bg:#1e1f2e;--plot-grid:#2d2e42;--plot-font:#ccc;
  --glass:rgba(15,17,23,.92);
}
*{margin:0;padding:0;box-sizing:border-box}
body{font-family:'Inter',-apple-system,BlinkMacSystemFont,sans-serif;background:var(--bg);
  color:var(--text);line-height:1.65;transition:background .3s,color .3s}
a{color:var(--primary);text-decoration:none}

/* Layout */
.layout{display:flex;min-height:100vh}
.sidebar{width:220px;position:sticky;top:0;height:100vh;background:var(--glass);
  backdrop-filter:blur(12px);border-right:1px solid var(--border);padding:1.2rem .8rem;
  display:flex;flex-direction:column;z-index:10;overflow-y:auto;transition:background .3s,border .3s}
.sidebar .logo{font-weight:700;font-size:1.1rem;color:var(--dark);margin-bottom:1.2rem;
  padding-left:.7rem;letter-spacing:-.3px}
[data-theme="dark"] .sidebar .logo{color:var(--text)}
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

/* Main + Header */
.main{flex:1;min-width:0;overflow-x:hidden}
.header{background:linear-gradient(135deg,var(--dark) 0%,#16213e 50%,#0f3460 100%);
  color:#fff;padding:3rem 2.5rem 2rem;position:relative;overflow:hidden}
.header::before{content:'';position:absolute;inset:0;
  background:radial-gradient(circle at 80% 20%,rgba(255,255,255,.08) 0%,transparent 60%)}
.header h1{font-size:2rem;font-weight:700;letter-spacing:-.5px;position:relative}
.header h1 span{color:var(--primary)}
.header .subtitle{opacity:.8;font-size:.92rem;margin-top:.4rem;position:relative}
.header .meta{display:flex;gap:1.2rem;margin-top:1rem;font-size:.78rem;opacity:.7;
  flex-wrap:wrap;position:relative}
.header .meta span{background:rgba(255,255,255,.1);padding:.2rem .6rem;border-radius:20px}

/* Container + Sections */
.container{max-width:1100px;margin:0 auto;padding:2rem 2rem 3rem}
.section{margin-bottom:2.5rem;opacity:0;transform:translateY(16px);
  transition:opacity .5s ease,transform .5s ease}
.section.visible{opacity:1;transform:none}
.section h2{font-size:1.2rem;font-weight:700;color:var(--text);margin-bottom:1rem;
  display:flex;align-items:center;gap:.5rem;letter-spacing:-.3px}
.section h2::before{content:'';width:5px;height:1.2em;border-radius:3px;
  background:linear-gradient(180deg,var(--primary),var(--primary-dark));flex-shrink:0}

/* Cards */
.card{background:var(--surface);border:1px solid var(--border);border-radius:16px;
  padding:1.5rem;box-shadow:var(--shadow);transition:transform .3s,box-shadow .3s}
.card:hover{box-shadow:var(--shadow-lg);border-color:var(--border2)}
.card.primary{border-left:5px solid var(--primary)}
.card.secondary{border-left:5px solid var(--secondary)}
.card.accent{border-left:5px solid var(--accent)}
.card h3{font-size:1.1rem;color:var(--text);margin-bottom:.8rem;font-weight:600;
  display:flex;align-items:center;gap:10px}

/* KPI Grid */
.kpi-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(160px,1fr));gap:1rem;margin-bottom:1.5rem}
.kpi{background:var(--surface);border:1px solid var(--border);border-radius:16px;
  padding:1.2rem 1rem;text-align:center;box-shadow:var(--shadow);
  transition:transform .3s,box-shadow .3s;position:relative;overflow:hidden}
.kpi:hover{transform:translateY(-3px);box-shadow:var(--shadow-lg)}
.kpi::before{content:'';position:absolute;top:0;left:0;right:0;height:3px;
  background:linear-gradient(90deg,var(--primary),var(--secondary))}
.kpi .icon{font-size:1.6rem;margin-bottom:.4rem}
.kpi .value{font-size:1.5rem;font-weight:700;color:var(--primary);
  font-family:'JetBrains Mono',monospace;overflow-wrap:break-word;word-break:break-word;
  font-size:clamp(.85rem,3.5vw,1.5rem)}
.kpi .label{font-size:.7rem;color:var(--text3);text-transform:uppercase;
  letter-spacing:.8px;margin-top:.3rem}

/* Tables */
.table-wrap{overflow-x:auto;border-radius:12px;border:1px solid var(--border);
  box-shadow:0 4px 15px rgba(0,0,0,.05)}
table{width:100%;border-collapse:collapse;font-size:.9rem;background:var(--surface)}
th{background:var(--dark);color:#fff;padding:.7rem .9rem;text-align:left;font-weight:600;
  font-size:.72rem;text-transform:uppercase;letter-spacing:.5px;position:sticky;top:0}
[data-theme="dark"] th{background:#252640}
td{padding:.6rem .9rem;border-bottom:1px solid var(--border)}
tr:nth-child(even){background:rgba(0,102,204,.03)}
[data-theme="dark"] tr:nth-child(even){background:rgba(77,158,232,.05)}
tr:hover{background:rgba(0,102,204,.05)}
td.num{text-align:right;font-family:'JetBrains Mono',monospace;font-size:.82rem}

/* Badges */
.badge{display:inline-flex;align-items:center;gap:.25rem;padding:.2rem .6rem;
  border-radius:20px;font-size:.72rem;font-weight:600}
.badge-primary{background:rgba(0,102,204,.1);color:var(--primary)}
.badge-secondary{background:rgba(0,168,107,.1);color:var(--secondary)}
.badge-accent{background:rgba(255,107,53,.1);color:var(--accent)}
.badge-none{background:rgba(0,0,0,.04);color:var(--text3)}
[data-theme="dark"] .badge-none{background:rgba(255,255,255,.06)}

/* Plot + Map */
.plot-container{width:100%;min-height:400px;margin:1rem 0}
.map-container{border-radius:12px;overflow:hidden;border:1px solid var(--border)}
#basinMap{height:480px;width:100%;background:var(--surface);z-index:1}
[data-theme="dark"] .leaflet-tile{filter:brightness(.8) contrast(1.1) saturate(.85)}
[data-theme="dark"] .leaflet-popup-content-wrapper{background:#1e1f2e !important;color:#e2e8f0 !important}
[data-theme="dark"] .leaflet-popup-tip{border-top-color:#1e1f2e !important}
.leaflet-popup-content-wrapper{border-radius:8px !important;font-family:'Inter',sans-serif;font-size:.82rem}
.leaflet-popup-content{margin:10px 14px !important}
.leaflet-popup-content b{color:var(--primary)}

/* Highlight boxes */
.highlight-box{padding:1.2rem 1.5rem;border-radius:0 12px 12px 0;margin:1rem 0}
.highlight-box.info{background:linear-gradient(135deg,#eff6ff,#dbeafe);border-left:5px solid var(--primary)}
.highlight-box.success{background:linear-gradient(135deg,#ecfdf5,#d1fae5);border-left:5px solid var(--secondary)}
.highlight-box.warning{background:linear-gradient(135deg,#fff7ed,#ffedd5);border-left:5px solid var(--accent)}
[data-theme="dark"] .highlight-box.info{background:rgba(0,102,204,.08)}
[data-theme="dark"] .highlight-box.success{background:rgba(0,168,107,.08)}
[data-theme="dark"] .highlight-box.warning{background:rgba(255,107,53,.08)}

/* Module details (collapsible) */
details.module-details{margin-top:.8rem}
details.module-details>summary{cursor:pointer;padding:.6rem 1rem;border-radius:8px;
  background:var(--surface);border:1px solid var(--border);font-weight:600;font-size:.88rem;
  transition:all .15s;list-style:none;display:flex;align-items:center;gap:.5rem}
details.module-details>summary::before{content:'\25B8';font-size:.9rem;transition:transform .2s}
details[open].module-details>summary::before{transform:rotate(90deg)}
details.module-details>summary:hover{border-color:var(--primary);color:var(--primary)}
.module-content{padding:1rem;border:1px solid var(--border);border-top:none;
  border-radius:0 0 12px 12px;background:var(--surface)}
.module-content h4{font-size:.85rem;font-weight:600;color:var(--text2);margin:1rem 0 .5rem;
  text-transform:uppercase;letter-spacing:.5px}
.module-content h4:first-child{margin-top:0}
table.param-table{width:100%;border-collapse:collapse;font-size:.82rem;margin-bottom:.5rem}
table.param-table th{background:var(--dark);color:#fff;padding:.45rem .7rem;text-align:left;
  font-size:.7rem;text-transform:uppercase;letter-spacing:.5px}
[data-theme="dark"] table.param-table th{background:#252640}
table.param-table td{padding:.4rem .7rem;border-bottom:1px solid var(--border)}
table.param-table code{font-family:'JetBrains Mono',monospace;font-size:.78rem;
  background:rgba(0,102,204,.05);padding:.1rem .3rem;border-radius:4px}
details.nested-details{margin:.5rem 0}
details.nested-details>summary{cursor:pointer;padding:.4rem .7rem;border-radius:6px;
  font-size:.82rem;font-weight:500;color:var(--primary);background:rgba(0,102,204,.04);
  border:1px solid transparent;transition:all .15s;list-style:none;
  display:flex;align-items:center;gap:.4rem}
details.nested-details>summary::before{content:'\25B8';font-size:.8rem;transition:transform .2s}
details[open].nested-details>summary::before{transform:rotate(90deg)}
details.nested-details>summary:hover{border-color:var(--primary);background:rgba(0,102,204,.08)}

/* Footer */
.footer{text-align:center;padding:2rem;color:var(--text3);font-size:.78rem;
  border-top:1px solid var(--border)}
.footer .logo{font-weight:700;font-size:.9rem;color:var(--text);margin-bottom:.3rem}
.footer .logo span{color:var(--primary)}
.footer a{color:var(--primary)}
.footer a:hover{text-decoration:underline}

/* Animations */
@keyframes fadeInUp{from{opacity:0;transform:translateY(30px)}to{opacity:1;transform:none}}

/* Responsive */
@media(max-width:900px){
  .sidebar{display:none}
  .container{padding:1.2rem}
  .kpi-grid{grid-template-columns:repeat(2,1fr)}
  .header{padding:2rem 1.5rem 1.5rem}
}
</style>
</head>
<body>
""")

    # --- SIDEBAR ---
    nav_items = [
        ('project', 'Project'),
        ('summary', 'Summary'),
    ]
    if map_layers or grqa_geojson_str:
        nav_items.append(('map', 'Basin Map'))
    if grqa_stats:
        nav_items.append(('grqa', 'GRQA Observations'))
    nav_items.append(('config', 'Configuration'))
    if results:
        nav_items.append(('timeseries', 'Time Series'))
    if spatial_stats:
        nav_items.append(('spatial', 'Spatial Statistics'))
    if ss_summary:
        nav_items.append(('sources', 'Source/Sink Setup'))
    nav_items.append(('metadata', 'Run Metadata'))

    H.append('<div class="layout">')
    H.append('<aside class="sidebar">')
    H.append('<div class="logo">Open<span>WQ</span> Report</div>')
    H.append('<nav>')
    for sid, label in nav_items:
        H.append(f'<a href="#{sid}">{label}</a>')
    H.append('</nav>')
    H.append('<div class="theme-toggle">')
    H.append("""<button class="theme-btn" onclick="toggleTheme()">
<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
<circle cx="12" cy="12" r="5"/><line x1="12" y1="1" x2="12" y2="3"/>
<line x1="12" y1="21" x2="12" y2="23"/><line x1="4.22" y1="4.22" x2="5.64" y2="5.64"/>
<line x1="18.36" y1="18.36" x2="19.78" y2="19.78"/><line x1="1" y1="12" x2="3" y2="12"/>
<line x1="21" y1="12" x2="23" y2="12"/><line x1="4.22" y1="19.78" x2="5.64" y2="18.36"/>
<line x1="18.36" y1="5.64" x2="19.78" y2="4.22"/>
</svg>
<span>Toggle Theme</span>
</button>""")
    H.append('</div>')
    H.append('</aside>')

    # --- MAIN ---
    H.append('<div class="main">')

    # Header
    H.append(f"""<div class="header">
<h1>Open<span>WQ</span> &mdash; {project_name}</h1>
<p class="subtitle">{comment}</p>
<div class="meta">
<span>Authors: {authors}</span>
<span>Date: {date}</span>
<span>Host Model: {hostmodel}</span>
<span>Report: {now}</span>
</div>
</div>""")

    H.append('<div class="container">')

    # --- SECTION: Project Information ---
    H.append(f"""<div class="section" id="project">
<h2>Project Information</h2>
<div class="card primary"><div class="table-wrap"><table>
<tr><th>Property</th><th>Value</th></tr>
<tr><td>Project Name</td><td><strong>{project_name}</strong></td></tr>
<tr><td>Authors</td><td>{authors}</td></tr>
<tr><td>Description</td><td>{comment}</td></tr>
<tr><td>Host Model</td><td><span class="badge badge-primary">{hostmodel}</span></td></tr>
<tr><td>Date</td><td>{date}</td></tr>
</table></div></div>
</div>""")

    # --- SECTION: Summary KPIs ---
    n_species = len(chemical_species)
    n_compartments = len(compartments_and_cells)
    n_cells = sum(1 for cmp in compartments_and_cells.values()
                  for _ in cmp.values())
    n_outputs = len(results) if results else n_species

    H.append(f"""<div class="section" id="summary">
<h2>Summary</h2>
<div class="kpi-grid">
<div class="kpi"><div class="icon">&#x1f9ea;</div><div class="value">{n_species}</div><div class="label">Chemical Species</div></div>
<div class="kpi"><div class="icon">&#x1f4e6;</div><div class="value">{n_compartments}</div><div class="label">Compartments</div></div>
<div class="kpi"><div class="icon">&#x2699;</div><div class="value">{solver}</div><div class="label">Solver</div></div>
<div class="kpi"><div class="icon">&#x23f1;</div><div class="value">{runtime_str}</div><div class="label">Runtime</div></div>
<div class="kpi"><div class="icon">&#x1f4ca;</div><div class="value">{n_outputs}</div><div class="label">Output Species</div></div>
<div class="kpi"><div class="icon">&#x1f552;</div><div class="value">{timestep[0]} {timestep[1]}</div><div class="label">Output Timestep</div></div>
</div>
</div>""")

    # Notice when model was not run
    if not model_was_run:
        H.append("""<div class="section visible" style="margin-bottom:1.5rem">
<div class="highlight-box warning">
<strong>Configuration-only report</strong> &mdash; The model was not executed.
This report shows the selected modules and their parameters.
Time series, spatial statistics, and runtime metadata are not available.
</div></div>""")

    # --- SECTION: Basin Map ---
    if map_layers or grqa_geojson_str:
        H.append(f"""<div class="section" id="map">
<h2>Basin Map</h2>
<div class="card secondary">
<div class="map-container"><div id="basinMap"></div></div>
</div>
</div>""")

    # --- SECTION: GRQA Observations ---
    if grqa_stats:
        H.append('<div class="section" id="grqa"><h2>GRQA Observation Data</h2>')
        H.append('<div class="highlight-box info">'
                 '<strong>Global River Water Quality Archive (GRQA)</strong> — '
                 'Monitoring stations near the river network with observations '
                 'for the simulated species. '
                 '<a href="https://zenodo.org/records/15335450" target="_blank">'
                 'GRQA Dataset</a></div>')

        # Summary KPIs
        yr_range = ''
        if grqa_stats.get('year_start') and grqa_stats.get('year_end'):
            yr_range = f"{grqa_stats['year_start']}&ndash;{grqa_stats['year_end']}"
        H.append(f"""<div class="kpi-grid" style="margin:1rem 0">
<div class="kpi"><div class="icon">&#x1f4cd;</div><div class="value">{grqa_stats['n_stations']}</div><div class="label">Stations</div></div>
<div class="kpi"><div class="icon">&#x1f4c8;</div><div class="value">{grqa_stats['n_observations']:,}</div><div class="label">Observations</div></div>
<div class="kpi"><div class="icon">&#x1f4c5;</div><div class="value">{yr_range}</div><div class="label">Period</div></div>
<div class="kpi"><div class="icon">&#x1f9ea;</div><div class="value">{len(grqa_stats.get('species_stats', []))}</div><div class="label">Species Matched</div></div>
</div>""")

        # Per-species table
        if grqa_stats.get('species_stats'):
            H.append('<div class="card primary"><div class="table-wrap"><table>')
            H.append('<tr><th>Species</th><th>Stations</th><th>Observations</th>'
                     '<th>Period</th></tr>')
            for sp in grqa_stats['species_stats']:
                yr = ''
                if sp.get('year_start') and sp.get('year_end'):
                    yr = f"{sp['year_start']}&ndash;{sp['year_end']}"
                H.append(
                    f'<tr><td><span class="badge badge-primary">{sp["species"]}</span></td>'
                    f'<td class="num">{sp["n_stations"]}</td>'
                    f'<td class="num">{sp["n_observations"]:,}</td>'
                    f'<td>{yr}</td></tr>')
            H.append('</table></div></div>')

        H.append('</div>')

    # --- SECTION: Configuration ---
    modules = [
        ('Biogeochemistry', bgc_module_name),
        ('Transport Dissolved', td_module_name),
        ('Lateral Exchange', le_module_name),
        ('Sediment Transport', ts_module_name),
        ('Sorption Isotherm', si_module_name),
    ]

    H.append('<div class="section" id="config"><h2>Model Configuration</h2>')
    H.append('<div class="card primary"><div class="table-wrap"><table>')
    H.append('<tr><th>Module</th><th>Selection</th></tr>')
    for mod_label, mod_name in modules:
        badge = 'badge-none' if mod_name.upper() == 'NONE' else 'badge-primary'
        H.append(f'<tr><td>{mod_label}</td><td><span class="badge {badge}">{mod_name}</span></td></tr>')
    H.append(f'<tr><td>Source/Sink Method</td><td><span class="badge badge-secondary">{ss_method}</span></td></tr>')
    H.append(f'<tr><td>Output Format</td><td>{units}</td></tr>')
    H.append(f'<tr><td>Species</td><td>{", ".join(chemical_species)}</td></tr>')
    cmp_names = ", ".join(compartments_and_cells.keys())
    H.append(f'<tr><td>Compartments</td><td>{cmp_names}</td></tr>')
    H.append('</table></div></div>')

    # Module parameter details (collapsible)
    for mod_label, mod_name in modules:
        mod_data = _read_module_json(output_dir, mod_name)
        if mod_data is None:
            continue
        param_html = _render_module_params(mod_data, mod_name)
        if not param_html:
            continue
        H.append(f'<details class="module-details">')
        H.append(f'<summary>{mod_label}: '
                 f'<span class="badge badge-primary">{mod_name}</span></summary>')
        H.append('<div class="module-content">')
        H.extend(param_html)
        H.append('</div></details>')

    H.append('</div>')  # close section

    # --- SECTION: Time Series ---
    if results:
        H.append('<div class="section" id="timeseries"><h2>Time Series</h2>')
        for i, (species, df) in enumerate(results.items()):
            plot_id = f"plotTS_{i}"
            H.append(f'<div class="card accent" style="margin-bottom:1rem">')
            H.append(f'<h3>{species}</h3>')
            H.append(f'<div id="{plot_id}" class="plot-container"></div>')
            H.append('</div>')
        H.append('</div>')

    # --- SECTION: Spatial Statistics ---
    if spatial_stats:
        H.append('<div class="section" id="spatial"><h2>Spatial Statistics</h2>')
        H.append('<div class="card secondary"><div class="table-wrap"><table>')
        H.append('<tr><th>Species</th><th>Min</th><th>Max</th><th>Mean</th>'
                 '<th>Std</th><th>Median</th><th>Cells</th><th>Timesteps</th></tr>')
        for s in spatial_stats:
            H.append(f'<tr><td>{s["species"]}</td>'
                     f'<td class="num">{s["min"]:.4g}</td>'
                     f'<td class="num">{s["max"]:.4g}</td>'
                     f'<td class="num">{s["mean"]:.4g}</td>'
                     f'<td class="num">{s["std"]:.4g}</td>'
                     f'<td class="num">{s["median"]:.4g}</td>'
                     f'<td class="num">{s["n_cells"]}</td>'
                     f'<td class="num">{s["n_timesteps"]}</td></tr>')
        H.append('</table></div></div></div>')

    # --- SECTION: Source/Sink Summary ---
    if ss_summary:
        H.append('<div class="section" id="sources"><h2>Source/Sink Setup</h2>')

        # Metadata card
        if ss_metadata:
            H.append('<div class="card primary" style="margin-bottom:1rem">'
                     '<div class="table-wrap"><table>')
            H.append('<tr><th>Property</th><th>Value</th></tr>')
            H.append(f'<tr><td>Method</td><td><span class="badge badge-secondary">'
                     f'{ss_method}</span></td></tr>')
            for mk, mv in ss_metadata.items():
                H.append(f'<tr><td>{mk}</td><td>{mv}</td></tr>')
            H.append('</table></div></div>')

        # Per-species detail table
        H.append('<div class="card accent"><div class="table-wrap"><table>')
        H.append('<tr><th>Chemical</th><th>Compartment</th><th>Type</th>'
                 '<th>Units</th><th>Cells</th><th>Period</th>'
                 '<th>Entries</th><th>Total Load</th></tr>')
        for entry in ss_summary:
            yr = entry.get('year_range')
            period = f"{yr[0]}&ndash;{yr[1]}" if yr and yr[0] != yr[1] else (
                str(yr[0]) if yr else 'N/A')
            total = entry.get('total_load', 0)
            total_str = f"{total:.4g}" if total else "0"
            H.append(
                f'<tr><td>{entry["chemical"]}</td>'
                f'<td>{entry["compartment"]}</td>'
                f'<td>{entry["type"]}</td>'
                f'<td>{entry["units"]}</td>'
                f'<td class="num">{entry.get("n_cells", "N/A")}</td>'
                f'<td>{period}</td>'
                f'<td class="num">{entry.get("n_entries", "N/A")}</td>'
                f'<td class="num">{total_str} {entry["units"]}</td></tr>')
        H.append('</table></div></div>')

        H.append('</div>')

    # --- SECTION: Run Metadata ---
    H.append('<div class="section" id="metadata"><h2>Run Metadata</h2>')
    H.append('<div class="card primary"><div class="table-wrap"><table>')
    H.append('<tr><th>Property</th><th>Value</th></tr>')
    if model_was_run:
        H.append(f'<tr><td>Container Runtime</td><td>{container_runtime}</td></tr>')
        H.append(f'<tr><td>MPI Processes</td><td>{mpi_np}</td></tr>')
        H.append(f'<tr><td>Runtime</td><td>{runtime_str}</td></tr>')
    else:
        H.append('<tr><td>Model Execution</td>'
                 '<td><span class="badge badge-none">Not executed</span></td></tr>')
    H.append(f'<tr><td>Output Directory</td>'
             f'<td style="font-family:monospace;font-size:.8rem">{output_dir}</td></tr>')
    H.append(f'<tr><td>Report Generated</td><td>{now}</td></tr>')
    H.append(f'<tr><td>Solver</td><td>{solver}</td></tr>')
    H.append(f'<tr><td>Output Timestep</td><td>{timestep[0]} {timestep[1]}</td></tr>')
    H.append('</table></div></div></div>')

    H.append('</div>')  # container

    # Footer
    H.append(f"""<div class="footer">
<div class="logo">Open<span>WQ</span></div>
<p>Simulation report generated on {now}</p>
<p><a href="https://github.com/DigitalWaterAnalytics/openwq" target="_blank">github.com/DigitalWaterAnalytics/openwq</a></p>
</div>""")

    H.append('</div>')  # main
    H.append('</div>')  # layout

    # =========================================================================
    # JAVASCRIPT
    # =========================================================================
    H.append('<script>')

    # Theme toggle
    H.append("""
(function(){
  var saved = localStorage.getItem('owq_theme');
  if(saved) document.documentElement.setAttribute('data-theme', saved);
})();
function toggleTheme(){
  var el = document.documentElement;
  var t = el.getAttribute('data-theme')==='dark' ? 'light' : 'dark';
  el.setAttribute('data-theme', t);
  localStorage.setItem('owq_theme', t);
  _owqRelayoutAll();
}
""")

    # Plotly layout helper
    H.append("""
function _owqLayout(extra){
  var t = document.documentElement.getAttribute('data-theme')==='dark';
  var base = {
    margin:{l:60,r:25,t:30,b:45},
    plot_bgcolor: t?'#1e1f2e':'#fff',
    paper_bgcolor: t?'#1e1f2e':'#fff',
    font:{color:t?'#ccc':'#333',family:'Inter,sans-serif',size:12},
    colorway:['#0066cc','#00a86b','#ff6b35','#004499','#34d399','#fb923c','#667eea','#764ba2'],
    xaxis:{gridcolor:t?'#2a2b3d':'#eee'},
    yaxis:{gridcolor:t?'#2a2b3d':'#eee'}
  };
  return Object.assign({},base,extra||{});
}
var PCFG={responsive:true,displayModeBar:true,
  modeBarButtonsToRemove:['lasso2d','select2d'],
  toImageButtonOptions:{format:'svg',filename:'openwq_plot'}};
window._owqPlotIds=[];
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

    # Time series plots
    if results:
        for i, (species, df) in enumerate(results.items()):
            plot_id = f"plotTS_{i}"
            H.append(f'window._owqPlotIds.push("{plot_id}");')

            # Select up to 10 representative columns
            cols = list(df.columns)
            if len(cols) > 10:
                step = max(1, len(cols) // 10)
                cols = cols[::step][:10]

            timestamps = df.index.tolist()
            ts_strings = [str(t) for t in timestamps]

            traces = []
            for col in cols:
                vals = df[col].tolist()
                # Replace NaN with None for JSON
                vals = [None if (isinstance(v, float) and np.isnan(v)) else v for v in vals]
                traces.append({
                    'x': ts_strings,
                    'y': vals,
                    'name': col,
                    'type': 'scatter',
                    'mode': 'lines',
                })

            H.append(f"""
Plotly.newPlot('{plot_id}',{json.dumps(traces)},
  _owqLayout({{
    title:'{species} ({units})',
    xaxis:{{title:'Time'}},
    yaxis:{{title:'{units}'}}
  }}),PCFG);
""")

    # Basin map — river network + basin + GRQA stations
    if (map_layers and map_center) or grqa_geojson_str:
        layer_js_blocks = []
        overlay_entries = []
        fit_var = None

        # Shapefile layers (basin polygons + river network)
        for idx, layer in enumerate(map_layers):
            var_name = f"lyr_{idx}"
            geojson_var = f"gj_{idx}"
            style = layer['style']
            label = layer['label'].replace("'", "\\'")
            layer_type = layer['type']

            layer_js_blocks.append(f"  var {geojson_var} = {layer['geojson_str']};")

            color = style.get('color', '#0066cc')
            weight = style.get('weight', 2.5)
            opacity = style.get('opacity', 0.85)
            fill_opacity = style.get('fillOpacity', 0.2)
            layer_js_blocks.append(f"""  var {var_name} = L.geoJSON({geojson_var}, {{
    style: function(f){{ return {{color:'{color}',weight:{weight},opacity:{opacity},fillOpacity:{fill_opacity}}}; }},
    onEachFeature: function(f,l){{
      var props = f.properties;
      var html = '<b>{label}</b><br>';
      for(var k in props){{ if(props[k]!==null && props[k]!=='') html += '<b>'+k+'</b>: '+props[k]+'<br>'; }}
      l.bindPopup(html);
    }}
  }}).addTo(map);""")
            overlay_entries.append(f"'{label}': {var_name}")
            if fit_var is None:
                fit_var = var_name

        # GRQA stations layer (circle markers)
        if grqa_geojson_str:
            layer_js_blocks.append(f"  var gj_grqa = {grqa_geojson_str};")
            layer_js_blocks.append("""  var lyr_grqa = L.geoJSON(gj_grqa, {
    pointToLayer: function(f, latlng){
      return L.circleMarker(latlng, {
        radius: 7, fillColor: '#e63946', color: '#333',
        weight: 1.5, opacity: 0.9, fillOpacity: 0.85
      });
    },
    onEachFeature: function(f,l){
      var props = f.properties;
      var html = '<b>GRQA Station</b><br>';
      for(var k in props){ if(props[k]!==null && props[k]!=='') html += '<b>'+k+'</b>: '+props[k]+'<br>'; }
      l.bindPopup(html);
    }
  }).addTo(map);""")
            overlay_entries.append("'GRQA Stations': lyr_grqa")
            if fit_var is None:
                fit_var = "lyr_grqa"

        overlay_obj = "{" + ", ".join(overlay_entries) + "}"
        center_js = _js(map_center) if map_center else "[0, 0]"
        zoom = 9 if map_center else 2

        H.append(f"""
(function(){{
  var map = L.map('basinMap').setView({center_js}, {zoom});
  var satTile = L.tileLayer('https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{{z}}/{{y}}/{{x}}',{{
    attribution:'Esri, Maxar, Earthstar Geographics', maxZoom:19}}).addTo(map);
  var lightTile = L.tileLayer('https://{{s}}.basemaps.cartocdn.com/light_all/{{z}}/{{x}}/{{y}}{{r}}.png',{{
    attribution:'&copy; OpenStreetMap &copy; CARTO', maxZoom:19
  }});
  var topoTile = L.tileLayer('https://{{s}}.tile.opentopomap.org/{{z}}/{{x}}/{{y}}.png',{{
    attribution:'OpenTopoMap', maxZoom:17}});
{chr(10).join(layer_js_blocks)}
  L.control.layers({{'Satellite':satTile,'Light':lightTile,'Topo':topoTile}}, {overlay_obj}, {{collapsed:false}}).addTo(map);
  var _fitBounds = null;
  try{{ _fitBounds = {fit_var}.getBounds(); map.fitBounds(_fitBounds); }}catch(e){{}}
  L.control.scale().addTo(map);
  var recenterBtn = L.control({{position:'topleft'}});
  recenterBtn.onAdd = function(){{
    var div = L.DomUtil.create('div','leaflet-bar leaflet-control');
    div.innerHTML = '<a href="#" title="Re-center map" style="font-size:18px;line-height:30px;width:30px;height:30px;display:block;text-align:center;text-decoration:none;color:#333">⌖</a>';
    div.firstChild.onclick = function(e){{
      e.preventDefault(); e.stopPropagation();
      if(_fitBounds) map.fitBounds(_fitBounds);
      else map.setView({center_js}, {zoom});
    }};
    return div;
  }};
  recenterBtn.addTo(map);
}})();
""")

    # Scroll animation + sidebar active state
    H.append("""
(function(){
  var sections = document.querySelectorAll('.section');
  var links = document.querySelectorAll('.sidebar a');
  var obs = new IntersectionObserver(function(entries){
    entries.forEach(function(e){
      if(e.isIntersecting) e.target.classList.add('visible');
    });
  },{threshold:0.1});
  sections.forEach(function(s){obs.observe(s);});

  window.addEventListener('scroll',function(){
    var fromTop=window.scrollY+80;
    links.forEach(function(link){
      var id=link.getAttribute('href');
      if(!id||!id.startsWith('#'))return;
      var sec=document.querySelector(id);
      if(!sec)return;
      if(sec.offsetTop<=fromTop && sec.offsetTop+sec.offsetHeight>fromTop){
        link.classList.add('active');
      }else{
        link.classList.remove('active');
      }
    });
  });
})();
""")

    H.append('</script>')
    H.append('</body></html>')

    # Write report
    html_content = '\n'.join(H)
    with open(report_path, 'w', encoding='utf-8') as f:
        f.write(html_content)

    return report_path


def generate_report(
        dir2save_input_files,
        project_name,
        authors,
        date,
        comment,
        hostmodel,
        bgc_module_name,
        td_module_name,
        le_module_name,
        ts_module_name,
        si_module_name,
        ss_method,
        solver,
        chemical_species,
        units,
        compartments_and_cells,
        timestep,
        river_network_shapefile=None,
        basin_shapefile=None,
        grqa_local_data_path=None,
        container_runtime="docker",
        mpi_np=2,
        run_time_seconds=0.0,
        model_was_run=True,
):
    """Generate an HTML simulation report (entry point for template).

    Wraps generate_simulation_report() with logging and error handling.

    Parameters:
        river_network_shapefile: path to river network shapefile
        basin_shapefile: path to basin/catchment shapefile
        grqa_local_data_path: path to pre-downloaded GRQA data (or None)

    Returns:
        str: Path to the generated HTML report, or None on failure
    """
    print("\n" + "=" * 60)
    print("GENERATING REPORT")
    print("=" * 60)

    try:
        report_path = generate_simulation_report(
            output_dir=dir2save_input_files,
            project_name=project_name,
            authors=authors,
            date=date,
            comment=comment,
            hostmodel=hostmodel,
            bgc_module_name=bgc_module_name,
            td_module_name=td_module_name,
            le_module_name=le_module_name,
            ts_module_name=ts_module_name,
            si_module_name=si_module_name,
            ss_method=ss_method,
            solver=solver,
            chemical_species=chemical_species,
            units=units,
            compartments_and_cells=compartments_and_cells,
            timestep=timestep,
            river_network_shapefile=river_network_shapefile,
            basin_shapefile=basin_shapefile,
            grqa_local_data_path=grqa_local_data_path,
            container_runtime=container_runtime,
            mpi_np=mpi_np,
            run_time_seconds=run_time_seconds,
            model_was_run=model_was_run,
        )

        print(f"  Report saved: {report_path}")
        print("  Open in browser to view interactive charts and maps.")
        return report_path

    except Exception as e:
        print(f"  WARNING: Report generation failed: {e}")
        import traceback
        traceback.print_exc()
        print("  Model outputs are still available in openwq_out/")
        return None
