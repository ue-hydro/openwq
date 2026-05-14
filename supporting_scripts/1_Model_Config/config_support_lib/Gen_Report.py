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
import platform
import numpy as np

# Ensure config_support_lib is on path (for local imports like Gen_BGC_Diagram)
_this_dir = os.path.dirname(os.path.abspath(__file__))
if _this_dir not in sys.path:
    sys.path.insert(0, _this_dir)

# HDF5 reader — add 2_Read_Outputs to path
_read_outputs_dir = os.path.normpath(
    os.path.join(_this_dir, '..', '..', '2_Read_Outputs', 'hdf5_support_lib'))
if _read_outputs_dir not in sys.path:
    sys.path.insert(0, _read_outputs_dir)

# GRQA extraction tools (standalone — no calibration_lib dependency)
try:
    from Gen_GRQA_Extract import (
        GRQACalibrationExtractor, SpeciesMapper, GRQA_PARAMETERS)
    _GRQA_AVAILABLE = True
except ImportError:
    try:
        from .Gen_GRQA_Extract import (
            GRQACalibrationExtractor, SpeciesMapper, GRQA_PARAMETERS)
        _GRQA_AVAILABLE = True
    except (ImportError, SystemError):
        _GRQA_AVAILABLE = False

try:
    from Gen_BGC_Diagram import generate_bgc_reaction_diagram
    _DIAGRAM_AVAILABLE = True
except ImportError:
    try:
        from .Gen_BGC_Diagram import generate_bgc_reaction_diagram
        _DIAGRAM_AVAILABLE = True
    except (ImportError, SystemError):
        _DIAGRAM_AVAILABLE = False

# Auto-mapping: model species name → GRQA parameter code.
#
# BGC templates use the "-N" / "-P" / "-S" convention to indicate that the
# species is reported "as element" using stoichiometry:
#   NO3-N  = nitrate as nitrogen   (mg-N/L)
#   NH4-N  = ammonium as nitrogen  (mg-N/L)
#   PO4-P  = phosphate as phosphorus (mg-P/L)
#   SO4-S  = sulfate as sulfur     (mg-S/L)
#
# GRQA reports the full ion (e.g. NO3 in mg-NO3/L).
# Conversion:  model_value = grqa_value × MW(element) / MW(ion)
_GRQA_TO_MODEL_FACTOR = {
    # GRQA code    factor          explanation
    'NO3':       14.007 / 62.004,  # N / NO3  ≈ 0.226
    'NH4':       14.007 / 18.039,  # N / NH4  ≈ 0.776
    'NO2':       14.007 / 46.006,  # N / NO2  ≈ 0.304
    'PO4':       30.974 / 94.971,  # P / PO4  ≈ 0.326
    'SO4':       32.065 / 96.06,   # S / SO4  ≈ 0.334
    # All other GRQA parameters: factor = 1.0 (no conversion needed)
}

_MODEL_SPECIES_TO_GRQA = {
    # ── "as N" species (templates use -N suffix) ──
    'NO3-N': 'NO3', 'NH4-N': 'NH4', 'NO2-N': 'NO2',
    'NO3-N_soil': 'NO3', 'NO3-N_tile': 'NO3',
    'NH4-N_soil': 'NH4', 'SED_NH4': 'NH4',
    # ── "as P" species (templates use -P suffix) ──
    'PO4-P': 'PO4',
    'PO4-P_sol': 'PO4', 'PO4-P_active': 'PO4', 'PO4-P_stable': 'PO4',
    'PO4-P_soil': 'PO4', 'PO4-P_tile': 'PO4',
    # ── "as S" species ──
    'SO4-S': 'SO4',
    # ── Direct GRQA codes (WASP8 and similar use full-ion names) ──
    'NO3': 'NO3', 'NH4': 'NH4', 'NO2': 'NO2', 'PO4': 'PO4', 'SO4': 'SO4',
    'DO': 'DO', 'DON': 'DON', 'DIN': 'DIN', 'TKN': 'TKN',
    'TN': 'TN', 'TP': 'TP', 'DP': 'DP', 'DOP': 'DOP', 'PP': 'PP',
    'BOD': 'BOD', 'BOD5': 'BOD5', 'COD': 'COD',
    'DOC': 'DOC', 'TOC': 'TOC', 'DIC': 'DIC', 'POC': 'POC',
    'TSS': 'TSS', 'TDS': 'TDS', 'SS': 'SS',
    # ── Sub-pool aggregations ──
    'BOD_fast': 'BOD', 'BOD_slow': 'BOD',
    'CBOD': 'BOD', 'CBOD_fast': 'BOD', 'CBOD_slow': 'BOD',
    'DOC_labile': 'DOC', 'DOC_refrac': 'DOC',
    'DOC_soil': 'DOC', 'DOC_tile': 'DOC',
    'POC_labile': 'POC', 'POC_refrac': 'POC',
    'TSS_tile': 'TSS', 'SED': 'TSS',
    'PON': 'PN', 'SED_PON1': 'PN', 'SED_PON2': 'PN',
    'DSi': 'Si', 'BSi': 'Si',
    'ALK': 'Alk', 'Chla': 'Chl_a', 'EC_generic': 'EC',
    # ── Organic-N / organic-P pool names ──
    'ORG_N': 'DON', 'ON': 'DON', 'ON_labile': 'DON', 'ON_refrac': 'TKN',
    'N_ORG': 'DON', 'N_ORG_active': 'DON', 'N_ORG_fresh': 'DON',
    'N_ORG_stable': 'TKN', 'N_ORG_soil': 'TKN',
    'ORG_N_active': 'DON', 'ORG_N_fresh': 'DON', 'ORG_N_stable': 'TKN',
    'ORG_P': 'DP', 'ORG_P_humus': 'DP', 'ORG_P_fresh': 'DP',
    'P_ORG': 'DP', 'P_ORG_active': 'DP',
    'OP_labile': 'DOP', 'OP_refrac': 'DP',
    'ORG_C_active': 'DOC', 'ORG_C_fresh': 'DOC', 'ORG_C_stable': 'TOC',
    'PART_P': 'PP', 'P_PART': 'PP', 'PartP': 'PP',
    # ── HYPE model names ──
    'IN': 'DIN', 'SP': 'PO4',
    'fastN': 'DIN', 'humusN': 'TKN', 'fastP': 'PO4', 'humusP': 'DP',
    # ── QUAL2E model names ──
    'DISS_P': 'DP',
    # ── Legacy / alternative names (user convenience) ──
    'NO3_N': 'NO3', 'NH4_N': 'NH4', 'NO2_N': 'NO2', 'PO4_P': 'PO4',
    'SO4_S': 'SO4',
    'Nitrate': 'NO3', 'Ammonium': 'NH4', 'Phosphate': 'PO4',
    'Dissolved Oxygen': 'DO', 'Total Nitrogen': 'TN',
    'Total Phosphorus': 'TP',
}

# Full set of GRQA parameters available in the archive (for the report).
_GRQA_AVAILABLE_PARAMS = {
    'TN': 'Total Nitrogen', 'DN': 'Dissolved Nitrogen',
    'NO3': 'Nitrate', 'NO2': 'Nitrite', 'NH4': 'Ammonium',
    'NO3_NO2': 'Nitrate + Nitrite', 'DIN': 'Dissolved Inorganic N',
    'DON': 'Dissolved Organic N', 'TKN': 'Total Kjeldahl N',
    'PN': 'Particulate N',
    'TP': 'Total Phosphorus', 'DP': 'Dissolved Phosphorus',
    'DIP': 'Dissolved Inorganic P', 'PO4': 'Phosphate',
    'DOP': 'Dissolved Organic P', 'PP': 'Particulate P',
    'DO': 'Dissolved Oxygen', 'DO_sat': 'DO Saturation',
    'BOD': 'Biochemical O Demand', 'BOD5': 'BOD (5-day)',
    'COD': 'Chemical O Demand', 'CODMn': 'COD (Permanganate)',
    'DOC': 'Dissolved Organic C', 'TOC': 'Total Organic C',
    'TC': 'Total Carbon', 'TIC': 'Total Inorganic C',
    'DIC': 'Dissolved Inorganic C', 'POC': 'Particulate Organic C',
    'PC': 'Particulate Carbon',
    'TSS': 'Total Suspended Solids', 'TDS': 'Total Dissolved Solids',
    'SS': 'Suspended Solids', 'Turbidity': 'Turbidity',
    'Temp': 'Water Temperature', 'pH': 'pH',
    'EC': 'Electrical Conductivity', 'SC': 'Specific Conductance',
    'Chl_a': 'Chlorophyll-a',
    'Si': 'Silica', 'SO4': 'Sulfate', 'Cl': 'Chloride', 'Alk': 'Alkalinity',
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


def _render_module_params(module_data, module_name, output_dir=None,
                          species_obs_availability=None):
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
    elif ('CHEMICAL_SPECIES' in module_data
          or 'BIOGEOCHEMISTRY_CONFIGURATION' in module_data):
        # Support both runtime format (top-level CHEMICAL_SPECIES)
        # and legacy format (nested under BIOGEOCHEMISTRY_CONFIGURATION)
        if 'BIOGEOCHEMISTRY_CONFIGURATION' in module_data:
            bgc = module_data['BIOGEOCHEMISTRY_CONFIGURATION']
        else:
            bgc = module_data

        if 'CHEMICAL_SPECIES' in bgc:
            chem_section = bgc['CHEMICAL_SPECIES']
            # Collect rows first, then decide which columns to render
            _list_descs = chem_section.get('_LIST_DESCRIPTIONS', {})
            _rows = []  # list of (name, description, units, ic)

            if 'LIST' in chem_section:
                for idx in sorted(chem_section['LIST'].keys(),
                                  key=lambda x: int(x)):
                    name = chem_section['LIST'][idx]
                    info = chem_section.get(f'_{name}_INFO', {})
                    desc = (_list_descs.get(name, '')
                            or info.get('DESCRIPTION', ''))
                    units = info.get('UNITS', '')
                    ic = info.get('INITIAL_CONDITION', '')
                    _rows.append((name, desc, units, ic))
            else:
                for k, v in chem_section.items():
                    if k.startswith('_') or not isinstance(v, dict):
                        continue
                    desc = (_list_descs.get(k, '')
                            or v.get('DESCRIPTION', ''))
                    _rows.append((k, desc,
                                  v.get('UNITS', ''),
                                  v.get('INITIAL_CONDITION', '')))

            # Determine which optional columns have data
            _has_desc = any(r[1] for r in _rows)
            _has_units = any(r[2] for r in _rows)
            _has_ic = any(r[3] for r in _rows)

            H.append('<h4>Chemical Species</h4>')
            H.append('<table class="param-table">')
            _hdr = '<tr><th>Species</th>'
            if _has_desc:
                _hdr += '<th>Description</th>'
            if _has_units:
                _hdr += '<th>Units</th>'
            if _has_ic:
                _hdr += '<th>Initial Condition</th>'
            _hdr += '</tr>'
            H.append(_hdr)

            for name, desc, units, ic in _rows:
                _row = f'<tr><td><code>{name}</code></td>'
                if _has_desc:
                    _row += f'<td>{desc}</td>'
                if _has_units:
                    _row += f'<td>{units}</td>'
                if _has_ic:
                    _row += f'<td class="num">{ic}</td>'
                _row += '</tr>'
                H.append(_row)
            H.append('</table>')

        # Reaction network diagram (right after Chemical Species)
        if _DIAGRAM_AVAILABLE:
            try:
                _diag = generate_bgc_reaction_diagram(
                    bgc_data=bgc,
                    species_obs_availability=species_obs_availability,
                )
                if _diag:
                    H.append('<h4>Reaction Network</h4>')
                    H.append(_diag)
            except Exception:
                pass

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

                # Runtime format: numbered entries with LIST_TRANSFORMATIONS
                if 'LIST_TRANSFORMATIONS' in fw_data:
                    lt = fw_data['LIST_TRANSFORMATIONS']
                    for idx in sorted(lt.keys(), key=lambda x: int(x)):
                        t_data = fw_data.get(idx, {})
                        t_name = t_data.get('_NAME', lt[idx])
                        consumed = t_data.get('CONSUMED', '')
                        produced = t_data.get('PRODUCED', '')
                        kin = t_data.get('KINETICS', '')
                        if isinstance(kin, list):
                            kin = f'{kin[0]} [{kin[1]}]'
                        H.append(
                            f'<tr><td><code>{t_name}</code></td>'
                            f'<td>{consumed}</td>'
                            f'<td>{produced}</td>'
                            f'<td><code style="font-size:.75rem">'
                            f'{kin}</code></td></tr>')
                else:
                    # Legacy format
                    for t_name, t_data in fw_data.items():
                        if t_name.startswith('_'):
                            continue
                        if not isinstance(t_data, dict):
                            continue
                        if ('KINETICS' not in t_data
                                and 'FORMULA' not in t_data):
                            continue
                        consumed = t_data.get('CONSUMED', '')
                        produced = t_data.get('PRODUCED', '')
                        kinetics = t_data.get(
                            'KINETICS', t_data.get('FORMULA', ''))
                        H.append(
                            f'<tr><td><code>{t_name}</code></td>'
                            f'<td>{consumed}</td>'
                            f'<td>{produced}</td>'
                            f'<td><code style="font-size:.75rem">'
                            f'{kinetics}</code></td></tr>')
                H.append('</table></details>')

        # Compartment assignment (runtime: _COMPARTMENT_ASSIGNMENT, legacy: COMPARTMENT_ASSIGNMENT)
        ca = bgc.get('COMPARTMENT_ASSIGNMENT',
                      bgc.get('_COMPARTMENT_ASSIGNMENT'))
        if ca:
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

        # PHREEQC reaction network diagram
        if _DIAGRAM_AVAILABLE:
            pqi_path = os.path.join(output_dir, 'openwq_in',
                                    'openWQ_PHREEQC.pqi')
            if os.path.isfile(pqi_path):
                try:
                    _diag = generate_bgc_reaction_diagram(
                        pqi_filepath=pqi_path, module_name='PHREEQC')
                    if _diag:
                        H.append('<h4>Reaction Network</h4>')
                        H.append(_diag)
                except Exception:
                    pass

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


def _build_map_layers(river_network_shapefile, basin_shapefile, hostmodel="mizuroute"):
    """Build map layers from the river network and basin shapefiles.

    Both layers are always included when available.  The *hostmodel* parameter
    controls which layer is **interactive** (clickable, linked to plots) and
    which is a background reference layer (lighter, transparent, non-clickable).

    - mizuRoute: river network is interactive, basin is background.
    - SUMMA:     basin is interactive, river network is background.

    Returns:
        (layers, map_center) where layers is a list of dicts with
        'geojson_str', 'label', 'style', 'type', 'interactive', and
        map_center is [lat, lng].
        Returns ([], None) if no layers loaded.
    """
    is_summa = hostmodel.lower() == "summa"
    layers = []
    all_bounds = []

    # Layer 1: Basin polygons
    if basin_shapefile:
        geojson, bounds = _load_single_geojson(basin_shapefile)
        if geojson is not None:
            if is_summa:
                # SUMMA: basin is the interactive layer
                style = {'color': '#ffffff', 'weight': 1.5, 'opacity': 0.9,
                         'fillColor': '#4dabf7', 'fillOpacity': 0.25}
            else:
                # mizuRoute: basin is a background reference layer
                style = {'color': '#ffffff', 'weight': 0.8, 'opacity': 0.4,
                         'fillColor': '#aaaaaa', 'fillOpacity': 0.06}
            layers.append({
                'geojson_str': json.dumps(geojson),
                'label': 'Basin',
                'style': style,
                'type': 'polygon',
                'interactive': is_summa,
            })
            all_bounds.append(bounds)

    # Layer 2: River network
    if river_network_shapefile:
        geojson, bounds = _load_single_geojson(river_network_shapefile)
        if geojson is not None:
            if not is_summa:
                # mizuRoute: river network is the interactive layer
                style = {'color': '#0066cc', 'weight': 2.5, 'opacity': 0.85}
            else:
                # SUMMA: river network is a background reference layer
                style = {'color': '#6699cc', 'weight': 1.2, 'opacity': 0.35}
            layers.append({
                'geojson_str': json.dumps(geojson),
                'label': 'River Network',
                'style': style,
                'type': 'line',
                'interactive': not is_summa,
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


def _extract_grqa_for_report(river_network_shapefile, basin_shapefile,
                              chemical_species, output_dir,
                              grqa_local_data_path=None,
                              grqa_buffer_km=10):
    """Extract GRQA observation stations for the simulated species.

    Uses the standalone Gen_GRQA_Extract module (no calibration_lib dependency).

    Search strategy:
      1. If a basin shapefile is provided, use it as the search area
         (with a user-defined buffer around the basin boundary).
      2. Otherwise use the buffer around the river network.

    After extraction the clipped observations are saved to
    ``openwq_in/grqa_clipped_data/`` and the raw cache (full GRQA CSVs)
    is deleted to free disk space.

    Returns:
        (stations_geojson_str, grqa_stats) or (None, None) on failure.
    """
    if not _GRQA_AVAILABLE:
        print("  GRQA tools not available (missing geopandas/shapely). Skipping.")
        return None, None

    import geopandas as _gpd
    from shapely.geometry import mapping as _shp_mapping
    from shapely.ops import unary_union as _unary_union
    import pandas as _pd

    # Need at least one shapefile to define the search area
    shp_for_search = basin_shapefile or river_network_shapefile
    if not shp_for_search or not os.path.isfile(shp_for_search):
        print("  No shapefile available for GRQA spatial search. Skipping.")
        return None, None

    # Build species mapping: GRQA code -> model species name
    species_mapping = {}
    unmapped_species = []
    for model_name in chemical_species:
        grqa_code = _MODEL_SPECIES_TO_GRQA.get(model_name)
        if grqa_code and grqa_code not in species_mapping:
            species_mapping[grqa_code] = model_name
        elif not grqa_code:
            unmapped_species.append(model_name)

    if not species_mapping:
        print("  No GRQA-compatible species found. Skipping GRQA extraction.")
        return None, None

    grqa_params = list(species_mapping.keys())
    buffer_m = grqa_buffer_km * 1000
    print(f"  GRQA species mapping: {species_mapping}")
    print(f"  Will extract {len(grqa_params)} GRQA parameters: {', '.join(grqa_params)}")
    print(f"  Search buffer: {grqa_buffer_km} km")
    if grqa_local_data_path is None:
        print("  No local GRQA data provided — will download from Zenodo (one-time, ~1.2 GB).")
        print("  Tip: set grqa_local_data_path to skip future downloads.")

    grqa_output_dir = os.path.join(output_dir, 'openwq_in', 'grqa_clipped_data')

    try:
        import shutil

        # Determine search buffer
        if basin_shapefile and os.path.isfile(basin_shapefile):
            print(f"  Using basin shapefile for GRQA search area (+ {grqa_buffer_km} km buffer).")
            basin_gdf = _gpd.read_file(basin_shapefile)
            if basin_gdf.crs and not basin_gdf.crs.is_geographic:
                basin_gdf = basin_gdf.to_crs(epsg=4326)
            centroid = basin_gdf.geometry.unary_union.centroid
            utm_zone = int((centroid.x + 180) / 6) + 1
            hemi = 'north' if centroid.y >= 0 else 'south'
            utm_epsg = 32600 + utm_zone if hemi == 'north' else 32700 + utm_zone
            basin_proj = basin_gdf.to_crs(epsg=utm_epsg)
            search_geom = _unary_union(basin_proj.geometry).buffer(buffer_m)
            search_gdf = _gpd.GeoDataFrame(geometry=[search_geom],
                                           crs=f'EPSG:{utm_epsg}')
            search_gdf = search_gdf.to_crs(epsg=4326)
        else:
            print(f"  Using river network with {grqa_buffer_km} km buffer for GRQA search area.")

        mapper = SpeciesMapper(mapping=species_mapping)
        extractor = GRQACalibrationExtractor(
            output_dir=grqa_output_dir,
            species_mapper=mapper,
            local_data_path=grqa_local_data_path,
            buffer_distance_m=buffer_m,
        )

        if basin_shapefile and os.path.isfile(basin_shapefile):
            stations, observations = extractor.extract_stations_and_observations(
                search_gdf)
        else:
            river_gdf = extractor.load_river_network(river_network_shapefile)
            buffer_gdf = extractor.create_buffer(river_gdf)
            stations, observations = extractor.extract_stations_and_observations(
                buffer_gdf)

        # Clean up: delete the large raw cache files
        cache_dir = os.path.join(grqa_output_dir, 'grqa_cache')
        if os.path.isdir(cache_dir):
            print(f"  Cleaning up GRQA cache ({cache_dir})...")
            shutil.rmtree(cache_dir, ignore_errors=True)

        if stations is None or stations.empty:
            print("  No GRQA stations found in the search area.")
            grqa_stats = {
                'n_stations': 0, 'n_observations': 0,
                'year_start': None, 'year_end': None,
                'species_stats': [],
                'output_dir': grqa_output_dir,
                'searched_species': list(species_mapping.values()),
                'unmapped_species': unmapped_species,
                'buffer_km': grqa_buffer_km,
            }
            return None, grqa_stats

        # Generate clipped observation files (CSV)
        extractor.generate_calibration_files(prefix="grqa_clipped")

        # Build GeoJSON for map
        features = []
        for _, row in stations.iterrows():
            props = {k: (str(v) if v is not None else '')
                     for k, v in row.items() if k != 'geometry'}
            features.append({
                'type': 'Feature',
                'geometry': _shp_mapping(row.geometry),
                'properties': props,
            })
        stations_geojson = {'type': 'FeatureCollection', 'features': features}

        # Compute statistics
        site_col = extractor.column_map.get('site_id', 'site_id')
        date_col = extractor.column_map.get('obs_date', 'obs_date')
        dates = _pd.to_datetime(observations[date_col], errors='coerce')
        species_stats = []
        found_model_species = set()
        for model_sp in observations['model_species'].unique():
            sp_obs = observations[observations['model_species'] == model_sp]
            sp_dates = _pd.to_datetime(sp_obs[date_col], errors='coerce')
            sp_stations = sp_obs[site_col].nunique()
            species_stats.append({
                'species': model_sp,
                'n_stations': sp_stations,
                'n_observations': len(sp_obs),
                'year_start': int(sp_dates.min().year) if not sp_dates.isna().all() else None,
                'year_end': int(sp_dates.max().year) if not sp_dates.isna().all() else None,
            })
            found_model_species.add(model_sp)

        no_data_species = [s for s in species_mapping.values()
                           if s not in found_model_species]

        grqa_stats = {
            'n_stations': len(stations),
            'n_observations': len(observations),
            'year_start': int(dates.min().year) if not dates.isna().all() else None,
            'year_end': int(dates.max().year) if not dates.isna().all() else None,
            'species_stats': species_stats,
            'output_dir': grqa_output_dir,
            'buffer_km': grqa_buffer_km,
            'no_data_species': no_data_species,
            'unmapped_species': unmapped_species,
        }

        return json.dumps(stations_geojson), grqa_stats

    except Exception:
        raise


def _load_clipped_observations_for_report(output_dir, chemical_species):
    """Load pre-extracted observation data from ``openwq_in/grqa_clipped_data/``.

    This is a lightweight loader that reads the CSV files produced by the
    calibration tools (GRQA extraction).  It does NOT depend on
    ``calibration_lib`` — only standard library + pandas.

    Expected files:
        openwq_in/grqa_clipped_data/grqa_clipped_observations.csv
        openwq_in/grqa_clipped_data/grqa_clipped_stations.csv

    Returns:
        (stations_geojson_str, obs_stats) or (None, None) when data is absent.
    """
    clipped_dir = os.path.join(output_dir, 'openwq_in', 'grqa_clipped_data')
    obs_csv = os.path.join(clipped_dir, 'grqa_clipped_observations.csv')
    stn_csv = os.path.join(clipped_dir, 'grqa_clipped_stations.csv')

    if not os.path.isfile(obs_csv):
        print(f"  No pre-extracted observation data found at {clipped_dir}")
        print("  To generate it, run the GRQA extraction from the calibration tools first,")
        print("  or set observation_data_source = 'user_csv' with your own CSV.")
        return None, None

    try:
        import pandas as pd

        observations = pd.read_csv(obs_csv)
        if observations.empty:
            print("  Clipped observations CSV is empty.")
            return None, None

        # Detect column names (GRQA clipped format)
        site_col = 'site_id'
        date_col = 'obs_date'
        lat_col = 'lat_wgs84' if 'lat_wgs84' in observations.columns else 'lat'
        lon_col = 'lon_wgs84' if 'lon_wgs84' in observations.columns else 'lon'

        # Build station GeoJSON from the stations CSV or from observations
        features = []
        if os.path.isfile(stn_csv):
            stations = pd.read_csv(stn_csv)
            stn_lat = 'lat_wgs84' if 'lat_wgs84' in stations.columns else 'lat'
            stn_lon = 'lon_wgs84' if 'lon_wgs84' in stations.columns else 'lon'
            for _, row in stations.iterrows():
                props = {k: str(v) for k, v in row.items()
                         if k not in (stn_lat, stn_lon)}
                features.append({
                    'type': 'Feature',
                    'geometry': {
                        'type': 'Point',
                        'coordinates': [float(row[stn_lon]), float(row[stn_lat])],
                    },
                    'properties': props,
                })
        else:
            # Fallback: derive stations from the observations CSV
            stn_grp = observations.groupby(site_col).first()
            for sid, row in stn_grp.iterrows():
                features.append({
                    'type': 'Feature',
                    'geometry': {
                        'type': 'Point',
                        'coordinates': [float(row[lon_col]), float(row[lat_col])],
                    },
                    'properties': {'site_id': str(sid)},
                })

        stations_geojson = {'type': 'FeatureCollection', 'features': features}

        # Per-species statistics
        dates = pd.to_datetime(observations[date_col], errors='coerce')
        species_col = 'model_species' if 'model_species' in observations.columns else 'parameter'
        species_stats = []
        found_species = set()
        for model_sp in observations[species_col].unique():
            sp_obs = observations[observations[species_col] == model_sp]
            sp_dates = pd.to_datetime(sp_obs[date_col], errors='coerce')
            sp_stations = sp_obs[site_col].nunique()
            species_stats.append({
                'species': model_sp,
                'n_stations': sp_stations,
                'n_observations': len(sp_obs),
                'year_start': int(sp_dates.min().year) if not sp_dates.isna().all() else None,
                'year_end': int(sp_dates.max().year) if not sp_dates.isna().all() else None,
            })
            found_species.add(model_sp)

        # Species in the model that have no observations
        chem_set = set(chemical_species) if chemical_species else set()
        no_data_species = sorted(chem_set - found_species)

        n_stations = len(features)
        obs_stats = {
            'n_stations': n_stations,
            'n_observations': len(observations),
            'year_start': int(dates.min().year) if not dates.isna().all() else None,
            'year_end': int(dates.max().year) if not dates.isna().all() else None,
            'species_stats': species_stats,
            'output_dir': clipped_dir,
            'no_data_species': no_data_species,
            'unmapped_species': [],
        }

        print(f"  Loaded {len(observations)} observations from {n_stations} stations.")
        print(f"  Species with data: {sorted(found_species)}")
        if no_data_species:
            print(f"  Species without observations: {no_data_species}")

        return json.dumps(stations_geojson), obs_stats

    except Exception as e:
        print(f"  WARNING: Failed to load clipped observations: {e}")
        return None, None


def _load_user_observations_for_report(user_observation_csv, chemical_species):
    """Load user-provided observation data from a CSV file.

    Expected CSV columns:
        station_id, lat, lon, parameter, year, month, day, minute, value, units

    The ``parameter`` column must contain species names that match the BGC
    template (e.g. "NO3-N", "PO4-P", "DO", "TSS").

    Parameters:
        user_observation_csv: path to the CSV file
        chemical_species: list of model species names (used for filtering/matching)

    Returns:
        (stations_geojson_str, obs_stats) or (None, None) on failure.
        obs_stats has the same structure as grqa_stats for interoperability.
    """
    if not user_observation_csv or not os.path.isfile(user_observation_csv):
        print(f"  User observation CSV not found: {user_observation_csv}")
        return None, None

    try:
        import pandas as pd

        required_cols = {'station_id', 'lat', 'lon', 'parameter',
                         'year', 'month', 'day', 'minute', 'value', 'units'}

        df = pd.read_csv(user_observation_csv)
        df.columns = df.columns.str.strip().str.lower()

        missing = required_cols - set(df.columns)
        if missing:
            print(f"  ERROR: User CSV missing required columns: {missing}")
            return None, None

        # Strip whitespace from string columns
        for col in ('station_id', 'parameter', 'units'):
            df[col] = df[col].astype(str).str.strip()

        # Filter to species present in the model configuration
        chem_set = set(chemical_species)
        df_matched = df[df['parameter'].isin(chem_set)]
        unmapped_csv = sorted(set(df['parameter'].unique()) - chem_set)
        unmatched_model = sorted(chem_set - set(df['parameter'].unique()))

        if df_matched.empty:
            print("  No matching species found between CSV and model configuration.")
            obs_stats = {
                'n_stations': 0,
                'n_observations': 0,
                'year_start': None,
                'year_end': None,
                'species_stats': [],
                'buffer_km': None,
                'no_data_species': list(chem_set),
                'unmapped_species': unmapped_csv,
                'searched_species': list(chem_set),
            }
            return None, obs_stats

        n_observations = len(df_matched)
        stations = df_matched.groupby('station_id').first()[['lat', 'lon']]
        n_stations = len(stations)

        year_min = int(df_matched['year'].min())
        year_max = int(df_matched['year'].max())

        # Per-species statistics
        species_stats = []
        found_species = set()
        for param, grp in df_matched.groupby('parameter'):
            sp_stations = grp['station_id'].nunique()
            sp_obs = len(grp)
            sp_yr_min = int(grp['year'].min())
            sp_yr_max = int(grp['year'].max())
            species_stats.append({
                'species': param,
                'n_stations': sp_stations,
                'n_observations': sp_obs,
                'year_start': sp_yr_min,
                'year_end': sp_yr_max,
            })
            found_species.add(param)

        no_data_species = sorted(chem_set - found_species)

        # Build GeoJSON for station markers
        features = []
        for stn_id, row in stations.iterrows():
            # Collect parameters at this station
            stn_params = df_matched[df_matched['station_id'] == stn_id]['parameter'].unique()
            features.append({
                'type': 'Feature',
                'geometry': {
                    'type': 'Point',
                    'coordinates': [float(row['lon']), float(row['lat'])],
                },
                'properties': {
                    'station_id': str(stn_id),
                    'parameters': ', '.join(sorted(stn_params)),
                    'n_observations': int(
                        df_matched[df_matched['station_id'] == stn_id].shape[0]),
                },
            })

        stations_geojson = {
            'type': 'FeatureCollection',
            'features': features,
        }

        obs_stats = {
            'n_stations': n_stations,
            'n_observations': n_observations,
            'year_start': year_min,
            'year_end': year_max,
            'species_stats': species_stats,
            'buffer_km': None,
            'no_data_species': no_data_species,
            'unmapped_species': unmapped_csv,
            'searched_species': list(chem_set),
        }

        print(f"  Loaded {n_observations} observations from {n_stations} stations "
              f"({len(species_stats)} species matched).")

        return json.dumps(stations_geojson), obs_stats

    except Exception:
        raise


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


def _parse_hostmodel_output_info(file_manager_path, hostmodel, host_exe_dir):
    """Parse the host-model control file to find the output directory and file prefix.

    For mizuRoute: reads ``<output_dir>`` and ``<case_name>`` tags.
    For SUMMA:     reads ``outputPath`` and ``outFilePrefix`` keys.

    Container paths (starting with ``/code/...``) are converted to host paths
    by matching path components with the known *host_exe_dir*.

    Returns (host_output_dir, file_prefix) or (None, None) on failure.
    """
    if not file_manager_path or not os.path.isfile(file_manager_path):
        return None, None

    try:
        with open(file_manager_path) as f:
            lines = f.readlines()
    except Exception:
        return None, None

    output_dir_raw = None
    file_prefix = None

    if hostmodel.lower() == "mizuroute":
        for line in lines:
            s = line.strip()
            if s.startswith('<output_dir>'):
                output_dir_raw = s.split('<output_dir>')[1].split('!')[0].strip()
            elif s.startswith('<case_name>'):
                file_prefix = s.split('<case_name>')[1].split('!')[0].strip()
    elif hostmodel.lower() == "summa":
        for line in lines:
            parts = line.strip().split()
            if len(parts) >= 2:
                if parts[0] == 'outputPath':
                    output_dir_raw = parts[1].strip("'\"")
                elif parts[0] == 'outFilePrefix':
                    file_prefix = parts[1].strip("'\"")

    if not output_dir_raw:
        return None, None

    # --- Resolve container path to host path ---
    norm = output_dir_raw.replace('//', '/').rstrip('/')

    # 1. Already a valid host path?
    if os.path.isdir(norm):
        return norm, file_prefix or ""

    # 2. Match path components against host_exe_dir to find the mount mapping.
    #    e.g. container "/code/proj/bin/mizuroute_out" vs host "/Users/.../proj/bin"
    #    → shared tail "proj/bin", extra = "mizuroute_out"
    #    → host result = "/Users/.../proj/bin/mizuroute_out"
    c_parts = [p for p in norm.split('/') if p]
    h_parts = [p for p in os.path.normpath(host_exe_dir).split(os.sep) if p]

    best_host_path = None
    best_match = 0
    for ci in range(len(c_parts)):
        for hi in range(len(h_parts)):
            ml = 0
            while (ci + ml < len(c_parts) and hi + ml < len(h_parts)
                   and c_parts[ci + ml] == h_parts[hi + ml]):
                ml += 1
            if ml > best_match:
                best_match = ml
                remaining = c_parts[ci + ml:]
                best_host_path = os.path.join(host_exe_dir, *remaining) if remaining else host_exe_dir

    if best_host_path and best_match >= 2 and os.path.isdir(best_host_path):
        return best_host_path, file_prefix or ""

    # 3. Fallback: try appending the last 1-3 components to host_exe_dir.
    for n in range(1, min(4, len(c_parts) + 1)):
        candidate = os.path.join(host_exe_dir, *c_parts[-n:])
        if os.path.isdir(candidate):
            return candidate, file_prefix or ""

    # 4. Return raw path as last resort.
    return norm, file_prefix or ""


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
        observation_data_source="grqa",
        grqa_local_data_path=None,
        grqa_buffer_km=10,
        user_observation_csv=None,
        observation_compartments=None,
        container_runtime="docker",
        mpi_np=2,
        docker_container_name="docker_openwq",
        executable_path=None,
        file_manager_path=None,
        config_errors=None,
        bgc_template_path=None,
        openwq_h5_mapping_key=None,
        seasonal_loads_method=None,
        plot_separator=' | ',
):
    """Generate a self-contained HTML simulation report.

    Parameters:
        river_network_shapefile: path to river network shapefile (.shp/.gpkg)
        basin_shapefile: path to basin/catchment shapefile (auto-loaded from
            ss_method_copernicus_basin_info if available)
        observation_data_source: "grqa", "user_csv", or "skip"
        grqa_local_data_path: "auto" or None to download from Zenodo,
            a path to pre-downloaded GRQA folder.
        grqa_buffer_km: buffer distance in km around basin/river for GRQA search
        user_observation_csv: path to user-provided CSV observation file
        seasonal_loads_method: temporal distribution method for annual loads
            (e.g. 'uniform', 'climate_adjusted')

    Returns the path to the generated HTML file.
    """
    report_path = os.path.join(output_dir, "openwq_config_report.html")
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")

    # Track errors for each section so the report can still be generated
    _section_errors = {}   # section_id → error message

    # Merge any config-generation errors passed from the caller
    if config_errors:
        _section_errors.update(config_errors)

    # (Time Series and Spatial Statistics sections removed — they require
    #  a completed model run which is not available at config time.)

    # Basin map layers (river network + basin polygons)
    map_layers, map_center = [], None
    try:
        map_layers, map_center = _build_map_layers(river_network_shapefile, basin_shapefile, hostmodel)
    except Exception as _e:
        print(f"  WARNING: Map layers failed: {_e}")
        _section_errors['map'] = str(_e)

    # Observation stations (GRQA or user-provided CSV)
    obs_geojson_str = None
    obs_stats = None
    _obs_source = (observation_data_source or "grqa").strip().lower()
    _valid_obs_sources = {"grqa", "user_csv", "skip"}
    if _obs_source not in _valid_obs_sources:
        print(f"  WARNING: Unrecognised observation_data_source = '{observation_data_source}'. "
              f"Valid options: {', '.join(sorted(_valid_obs_sources))}. Defaulting to 'skip'.")
        _obs_source = "skip"
    _obs_label = "GRQA"   # human-readable label for the observation source

    if _obs_source == "skip":
        print("  Observations: skipped by user (observation_data_source = 'skip').")
    elif _obs_source == "user_csv":
        _obs_label = "User CSV"
        try:
            print("  Loading user-provided observation CSV...")
            obs_geojson_str, obs_stats = _load_user_observations_for_report(
                user_observation_csv, chemical_species)
        except Exception as _e:
            print(f"  WARNING: User CSV observation loading failed: {_e}")
            _section_errors['observations'] = str(_e)
    else:
        # Default: GRQA extraction (or fallback to pre-extracted CSVs)
        _obs_label = "GRQA"
        if river_network_shapefile or basin_shapefile:
            try:
                print("  Extracting GRQA observation stations...")
                _grqa_path = None
                if isinstance(grqa_local_data_path, str) and \
                        grqa_local_data_path.strip().lower() not in ("auto", ""):
                    _grqa_path = grqa_local_data_path
                obs_geojson_str, obs_stats = _extract_grqa_for_report(
                    river_network_shapefile, basin_shapefile, chemical_species,
                    output_dir, _grqa_path, grqa_buffer_km=grqa_buffer_km)
            except Exception as _e:
                print(f"  WARNING: GRQA extraction failed: {_e}")
                print("  Falling back to pre-extracted observation data...")
                try:
                    obs_geojson_str, obs_stats = _load_clipped_observations_for_report(
                        output_dir, chemical_species)
                except Exception as _e2:
                    print(f"  WARNING: Fallback also failed: {_e2}")
                    _section_errors['observations'] = str(_e)
        else:
            # No shapefiles — try pre-extracted data
            try:
                print("  No shapefile for GRQA search. Looking for pre-extracted data...")
                obs_geojson_str, obs_stats = _load_clipped_observations_for_report(
                    output_dir, chemical_species)
            except Exception as _e:
                print(f"  WARNING: Observation data loading failed: {_e}")
                _section_errors['observations'] = str(_e)

    # Build species observation availability dict from obs_stats
    _species_obs_avail = None
    if obs_stats and chemical_species:
        _species_with_data = set(
            s.get('species', '') for s in obs_stats.get('species_stats', [])
        )
        _species_obs_avail = {
            sp: {"has_obs": sp in _species_with_data}
            for sp in chemical_species
        }

    # Read source/sink config for summary
    ss_config_path = os.path.join(output_dir, 'openwq_in',
                                  f'openWQ_SS_{ss_method}.json')
    ss_summary = []
    ss_metadata = {}
    if os.path.isfile(ss_config_path):
        try:
            with open(ss_config_path, 'r') as f:
                all_lines = f.readlines()

            # Extract metadata from // comment lines
            for l in all_lines:
                stripped = l.strip()
                if stripped.startswith('// METADATA - '):
                    # Parse "// METADATA - Key: Value"
                    meta_part = stripped[len('// METADATA - '):]
                    if ': ' in meta_part:
                        mk, mv = meta_part.split(': ', 1)
                        ss_metadata[mk.strip()] = mv.strip()

            # Also support legacy METADATA key inside JSON body
            json_lines = [l for l in all_lines
                          if not l.strip().startswith('//')]
            ss_data = json.loads(''.join(json_lines))

            if 'METADATA' in ss_data:
                ss_metadata.update(ss_data['METADATA'])

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
                data_format = entry.get('DATA_FORMAT', '')
                n_entries = 0
                cells = set()
                years = set()
                total_load = 0.0

                if data_format == 'ASCII' and 'FILEPATH' in data_block:
                    # ASCII/CSV: read the actual CSV file to extract stats
                    csv_path = data_block['FILEPATH']
                    csv_delim = data_block.get('DELIMITER', ',')

                    # The JSON stores Docker paths (/code/...) but the report
                    # runs on the host.  Reverse the mapping if needed.
                    if not os.path.isfile(csv_path):
                        try:
                            _rpt_dir = os.path.dirname(os.path.abspath(__file__))
                            _dc_yml = os.path.normpath(os.path.join(
                                _rpt_dir, '..', '..', '..',
                                'containers', 'docker-compose.yml'))
                            if os.path.isfile(_dc_yml):
                                import re as _re_mod
                                with open(_dc_yml, 'r') as _dcf:
                                    _dc_content = _dcf.read()
                                _vol_pat = _re_mod.compile(
                                    r'^\s*-\s+([^:]+):([^:]+)(?::.*)?$',
                                    _re_mod.MULTILINE)
                                _vol_matches = _vol_pat.findall(_dc_content)
                                if _vol_matches:
                                    _hp, _cp = _vol_matches[0]
                                    _hp = os.path.normpath(os.path.join(
                                        os.path.dirname(_dc_yml),
                                        _hp.strip()))
                                    if not _hp.endswith('/'):
                                        _hp += '/'
                                    _cp = _cp.strip()
                                    if not _cp.endswith('/'):
                                        _cp += '/'
                                    if csv_path.startswith(_cp):
                                        csv_path = _hp + csv_path[len(_cp):]
                        except Exception:
                            pass

                    try:
                        with open(csv_path, 'r') as _cf:
                            csv_lines = _cf.readlines()
                        hdr_cols = None
                        data_rows = []
                        for _cl in csv_lines:
                            if hdr_cols is None and 'YYYY' in _cl.upper():
                                hdr_cols = [c.strip().upper()
                                            for c in _cl.strip().split(csv_delim)]
                            elif hdr_cols is not None and _cl.strip():
                                data_rows.append(_cl.strip().split(csv_delim))
                        if hdr_cols and data_rows:
                            col_idx = {name: i for i, name in enumerate(hdr_cols)}
                            i_yyyy = col_idx.get('YYYY')
                            i_ix   = col_idx.get('IX')
                            i_load = col_idx.get('LOAD')
                            for dr in data_rows:
                                if i_yyyy is not None and i_yyyy < len(dr):
                                    try:
                                        years.add(int(dr[i_yyyy].strip()))
                                    except ValueError:
                                        pass
                                if i_ix is not None and i_ix < len(dr):
                                    cells.add(dr[i_ix].strip())
                                if i_load is not None and i_load < len(dr):
                                    try:
                                        total_load += float(dr[i_load].strip())
                                    except ValueError:
                                        pass
                            n_entries = len(data_rows)
                    except Exception:
                        pass  # file not found — leave stats at zero
                else:
                    # JSON format: DATA contains numbered entries
                    n_entries = len(data_block)
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
            _section_errors['sources'] = str(e)

    # =========================================================================
    # Build HTML
    # =========================================================================
    import html as _html_mod
    import traceback as _tb_mod

    def _error_card(section_title, error_msg, tb_str=None):
        """Return HTML for an inline error card shown when a section fails."""
        safe_msg = _html_mod.escape(str(error_msg))
        card = (f'<div class="card" style="border-left:4px solid #e74c3c;'
                f'background:rgba(231,76,60,.06);margin-bottom:1rem">'
                f'<p style="color:#e74c3c;font-weight:600;margin:0 0 .3rem 0">'
                f'&#x26a0; Error generating {_html_mod.escape(section_title)}</p>'
                f'<p style="font-size:.85rem;margin:0">{safe_msg}</p>')
        if tb_str:
            safe_tb = _html_mod.escape(tb_str)
            card += (f'<details style="margin-top:.5rem"><summary style="font-size:.8rem;'
                     f'cursor:pointer;color:var(--muted)">Show traceback</summary>'
                     f'<pre style="font-size:.75rem;overflow-x:auto;margin-top:.3rem;'
                     f'background:var(--glass);padding:.5rem;border-radius:4px">'
                     f'{safe_tb}</pre></details>')
        card += '</div>'
        return card

    def _safe_section(H, section_id, section_title, builder_fn):
        """Run builder_fn(H) inside a try/except. On failure, render an error
        card inside the section div so the rest of the report still works."""
        try:
            builder_fn(H)
        except Exception as _e:
            print(f"  WARNING: Section '{section_title}' failed: {_e}")
            H.append(f'<div class="section" id="{section_id}">'
                     f'<h2>{_html_mod.escape(section_title)}</h2>'
                     f'{_error_card(section_title, _e, _tb_mod.format_exc())}'
                     f'</div>')

    class _SectionGuard:
        """Context manager that rolls back H on exception and emits an error card."""
        def __init__(self, h_list, section_id, section_title):
            self.h_list = h_list
            self.section_id = section_id
            self.section_title = section_title
            self.snapshot_len = 0
            self.failed = False
        def __enter__(self):
            self.snapshot_len = len(self.h_list)
            return self
        def __exit__(self, exc_type, exc_val, exc_tb):
            if exc_type is not None:
                # Roll back any partial HTML appended by the failed section
                del self.h_list[self.snapshot_len:]
                print(f"  WARNING: Section '{self.section_title}' failed: {exc_val}")
                self.h_list.append(
                    f'<div class="section" id="{self.section_id}">'
                    f'<h2>{_html_mod.escape(self.section_title)}</h2>'
                    f'{_error_card(self.section_title, exc_val, _tb_mod.format_exc())}'
                    f'</div>')
                self.failed = True
                return True  # suppress exception

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
.section{margin-bottom:2.5rem}
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

/* BGC Reaction Network Diagram */
.bgc-diagram-container{overflow-x:auto;margin:1rem 0;border:1px solid var(--border);
  border-radius:12px;padding:1rem 1rem .6rem;background:var(--surface)}
.bgc-diagram-container svg{display:block;margin:0 auto;max-width:100%;
  font-family:'Inter',-apple-system,BlinkMacSystemFont,sans-serif}
.bgc-diagram-header{display:flex;justify-content:space-between;align-items:center;margin-bottom:.4rem}
.bgc-diagram-title{font-size:.9rem;font-weight:600;color:var(--text)}
.bgc-zoom-bar{display:flex;align-items:center;gap:4px}
.zoom-btn{width:26px;height:26px;border-radius:6px;border:1px solid var(--border);
  background:var(--surface);color:var(--text);font-size:.85rem;cursor:pointer;
  display:inline-flex;align-items:center;justify-content:center;transition:background .15s;
  font-family:inherit;line-height:1;padding:0}
.zoom-btn:hover{background:var(--border)}
.zoom-hint{font-size:.62rem;color:var(--text3);margin-left:4px}
/* Framework toggle bar */
.bgc-fw-toggle-bar{display:flex;flex-wrap:wrap;gap:6px;margin-bottom:.6rem;align-items:center}
.fw-toggle-chip{display:inline-flex;align-items:center;gap:4px;padding:3px 10px;
  font-size:.72rem;font-weight:500;border:1.5px solid var(--chip-color,var(--primary));
  border-radius:14px;background:transparent;color:var(--chip-color,var(--primary));
  cursor:pointer;transition:background .2s,opacity .2s;font-family:inherit;line-height:1.4}
.fw-toggle-chip.active{background:var(--chip-color,var(--primary));color:var(--surface,#fff)}
.fw-toggle-chip:not(.active){opacity:.45}
.fw-toggle-chip[data-fw="__all__"]{--chip-color:var(--text2)}
.chip-swatch{width:10px;height:10px;border-radius:50%;display:inline-block}
.fw-toggle-chip.active .chip-swatch{border:1.5px solid var(--surface,#fff)}
.chip-tick{font-size:.7rem;line-height:1;margin-right:1px;display:inline-block;transition:opacity .15s}
.fw-toggle-chip:not(.active) .chip-tick{opacity:0;width:0;margin:0;overflow:hidden}
/* Hover tooltip */
.bgc-tooltip{position:fixed;z-index:9999;pointer-events:none;max-width:420px;padding:8px 12px;
  background:var(--surface,#1a1b2e);color:var(--text,#e2e8f0);border:1px solid var(--border,#2a2b3d);
  border-radius:8px;font-size:.78rem;line-height:1.55;box-shadow:0 6px 24px rgba(0,0,0,.25);
  font-family:'Inter',-apple-system,BlinkMacSystemFont,sans-serif}
.bgc-tooltip strong{display:block;margin-bottom:3px;font-weight:600}
.bgc-tip-reaction{display:block;color:var(--text2,#a0aec0);font-size:.72rem;margin-bottom:4px}
.bgc-tip-expr{display:block;font-family:'JetBrains Mono',monospace;font-size:.7rem;
  color:var(--secondary,#34d399);background:rgba(0,0,0,.15);padding:3px 6px;
  border-radius:4px;margin:3px 0;word-break:break-all}
.bgc-tip-params{display:block;font-size:.68rem;color:var(--text3,#636e72);margin-top:2px}
/* Diagram interactivity */
.bgc-diagram-container .reaction-arrow{transition:opacity .25s ease,filter .25s ease;cursor:pointer}
.bgc-diagram-container .species-node{transition:opacity .25s ease,filter .25s ease;cursor:pointer}
.bgc-diagram-container svg{cursor:crosshair}
.bgc-diagram-container .reaction-arrow path.hit-area{pointer-events:stroke}
.bgc-diagram-container .reaction-arrow path.vis-path{pointer-events:none}
.bgc-diagram-container .reaction-arrow text{pointer-events:none}
.bgc-diagram-container .label-bg{pointer-events:none}

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
    if map_layers or obs_geojson_str:
        nav_items.append(('map', 'Basin Map'))
    if obs_stats:
        _obs_nav = f'{_obs_label} Observations' if _obs_label != "GRQA" else 'GRQA Observations'
        nav_items.append(('observations', _obs_nav))
    nav_items.append(('config', 'Configuration'))
    if ss_summary:
        nav_items.append(('sources', 'Source/Sink Setup'))
    nav_items.append(('metadata', 'Run Metadata'))
    _any_errors = bool(_section_errors)
    if _any_errors:
        nav_items.append(('consoleerrors', 'Errors'))
    if not _any_errors:
        nav_items.append(('nextsteps', 'Next Steps'))

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
    try:
        n_species = len(chemical_species)
        n_compartments = len(compartments_and_cells)
        # Count spatial elements from river network shapefile if available,
        # otherwise fall back to counting config output entries.
        n_cells = 0
        for layer in map_layers:
            if layer.get('label') == 'River Network':
                try:
                    _rn_geojson = json.loads(layer['geojson_str'])
                    n_cells = len(_rn_geojson.get('features', []))
                except Exception:
                    pass
                break
        if n_cells == 0:
            n_cells = sum(1 for cmp in compartments_and_cells.values()
                          for _ in cmp.values())
        n_outputs = n_species
        comp_names = ', '.join(compartments_and_cells.keys())

        # Extract simulation period from host-model control file
        sim_period_str = 'N/A'
        try:
            from Gen_Input_Driver import parse_sim_period_from_control_file
            sp = parse_sim_period_from_control_file(file_manager_path, hostmodel)
            sim_period_str = f"{sp[0]}&ndash;{sp[1]}" if sp[0] != sp[1] else str(sp[0])
        except Exception:
            pass

        # Count active (non-NONE) modules
        _all_modules = [bgc_module_name, td_module_name, le_module_name,
                        ts_module_name, si_module_name]
        _n_active = sum(1 for m in _all_modules
                        if m and m.upper() != 'NONE')
        _active_str = f'{_n_active}/{len(_all_modules)}'

        # Source/sink load count
        _n_ss = len(ss_summary) if ss_summary else 0

        # Runtime info
        _runtime_str = (f'{container_runtime}'
                        if container_runtime else 'N/A')

        # Observation stations
        _n_obs_stn = obs_stats['n_stations'] if obs_stats else 0

        H.append(f"""<div class="section" id="summary">
<h2>Summary</h2>
<div class="kpi-grid">
<div class="kpi"><div class="icon">&#x1f9ea;</div><div class="value">{n_species}</div><div class="label">Chemical Species</div><div style="font-size:0.7rem;color:#888;margin-top:2px;">({bgc_module_name})</div>{f'<div style="font-size:0.65rem;color:#888;margin-top:1px;">{os.path.basename(bgc_template_path)}</div>' if bgc_template_path else ''}</div>
<div class="kpi"><div class="icon">&#x1f4e6;</div><div class="value">{n_compartments}</div><div class="label">Compartments</div><div style="font-size:0.7rem;color:#888;margin-top:2px;">{comp_names}</div></div>
<div class="kpi"><div class="icon">&#x2699;</div><div class="value">{solver}</div><div class="label">Solver</div></div>
<div class="kpi"><div class="icon">&#x1f4ca;</div><div class="value">{n_outputs}</div><div class="label">Output Species</div></div>
<div class="kpi"><div class="icon">&#x1f552;</div><div class="value">{timestep[0]} {timestep[1]}</div><div class="label">Output Timestep</div></div>
<div class="kpi"><div class="icon">&#x1f4c5;</div><div class="value">{sim_period_str}</div><div class="label">Simulation Period</div></div>
<div class="kpi"><div class="icon">&#x1f310;</div><div class="value">{n_cells}</div><div class="label">Spatial Elements</div></div>
<div class="kpi"><div class="icon">&#x1f9e9;</div><div class="value">{_active_str}</div><div class="label">Active Modules</div></div>
<div class="kpi"><div class="icon">&#x1f4a7;</div><div class="value">{_n_ss}</div><div class="label">Source/Sink Loads</div></div>
<div class="kpi"><div class="icon">&#x1f4cd;</div><div class="value">{_n_obs_stn}</div><div class="label">Obs Stations</div></div>
</div>
</div>""")
    except Exception as _e:
        H.append(f'<div class="section" id="summary"><h2>Summary</h2>'
                 f'{_error_card("Summary", _e, _tb_mod.format_exc())}</div>')

    # --- SECTION: Basin Map ---
    try:
        if 'map' in _section_errors:
            H.append(f'<div class="section" id="map"><h2>Basin Map</h2>'
                     f'{_error_card("Basin Map", _section_errors["map"])}</div>')
        elif map_layers or obs_geojson_str:
            H.append(f"""<div class="section" id="map">
<h2>Basin Map</h2>
<div class="card secondary">
<div class="map-container"><div id="basinMap"></div></div>
</div>
</div>""")
    except Exception as _e:
        H.append(f'<div class="section" id="map"><h2>Basin Map</h2>'
                 f'{_error_card("Basin Map", _e, _tb_mod.format_exc())}</div>')

    # --- SECTION: Observation Data (GRQA or User CSV) ---
    _obs_section_title = f'{_obs_label} Observation Data'
    with _SectionGuard(H, 'observations', _obs_section_title):
        if obs_stats:
            H.append(f'<div class="section" id="observations">'
                     f'<h2>{_obs_section_title}</h2>')

            # Source-specific intro box
            if _obs_source == "user_csv":
                _csv_name = os.path.basename(user_observation_csv) \
                    if user_observation_csv else 'N/A'
                H.append(f'<div class="highlight-box info">'
                         f'<strong>User-Provided Observations</strong> — '
                         f'Loaded from <code>{_csv_name}</code>.</div>')
            else:
                _buf_km = obs_stats.get('buffer_km', '?')
                H.append(f'<div class="highlight-box info">'
                         f'<strong>Global River Water Quality Archive '
                         f'(GRQA)</strong> — '
                         f'Monitoring stations within {_buf_km} km of the '
                         f'basin/river network. '
                         f'<a href="https://zenodo.org/records/15335450" '
                         f'target="_blank">GRQA Dataset</a></div>')

            _no_data = obs_stats.get('no_data_species', [])
            _unmapped = obs_stats.get('unmapped_species', [])

            if obs_stats['n_stations'] == 0:
                # No stations/observations found
                searched = obs_stats.get('searched_species', [])
                sp_list = ', '.join(f'<code>{s}</code>' for s in searched) \
                    if searched else 'N/A'
                if _obs_source == "user_csv":
                    H.append(
                        f'<div class="highlight-box" '
                        f'style="border-left-color:var(--accent);">'
                        f'<strong>No matching observations found</strong> in '
                        f'the CSV for the following model species: {sp_list}. '
                        f'Make sure the <code>parameter</code> column in '
                        f'the CSV matches your BGC template species names.'
                        f'</div>')
                else:
                    _buf_km = obs_stats.get('buffer_km', '?')
                    H.append(
                        f'<div class="highlight-box" '
                        f'style="border-left-color:var(--accent);">'
                        f'<strong>No GRQA monitoring stations found</strong> '
                        f'within {_buf_km} km of the basin/river network for '
                        f'the following GRQA-supported species: {sp_list}. '
                        f'Try increasing <code>grqa_buffer_km</code> to widen '
                        f'the search area.</div>')
                if _unmapped:
                    _um_list = ', '.join(f'<code>{s}</code>' for s in _unmapped)
                    if _obs_source == "user_csv":
                        H.append(
                            f'<div class="highlight-box" '
                            f'style="border-left-color:var(--text-light);">'
                            f'The following species in the CSV are <strong>not '
                            f'in the model configuration</strong> and were '
                            f'excluded: {_um_list}.</div>')
                    else:
                        H.append(
                            f'<div class="highlight-box" '
                            f'style="border-left-color:var(--text-light);">'
                            f'The following model species are <strong>not '
                            f'available</strong> in the GRQA database and were '
                            f'excluded from the search: {_um_list}.</div>')
            else:
                # Summary KPIs
                yr_range = ''
                if obs_stats.get('year_start') and obs_stats.get('year_end'):
                    yr_range = (f"{obs_stats['year_start']}&ndash;"
                                f"{obs_stats['year_end']}")
                H.append(f"""<div class="kpi-grid" style="margin:1rem 0">
<div class="kpi"><div class="icon">&#x1f4cd;</div><div class="value">{obs_stats['n_stations']}</div><div class="label">Stations</div></div>
<div class="kpi"><div class="icon">&#x1f4c8;</div><div class="value">{obs_stats['n_observations']:,}</div><div class="label">Observations</div></div>
<div class="kpi"><div class="icon">&#x1f4c5;</div><div class="value">{yr_range}</div><div class="label">Period</div></div>
<div class="kpi"><div class="icon">&#x1f9ea;</div><div class="value">{len(obs_stats.get('species_stats', []))}</div><div class="label">Species Matched</div></div>
</div>""")

                # Per-species table
                if obs_stats.get('species_stats'):
                    H.append('<div class="card primary"><div class="table-wrap">'
                             '<table>')
                    H.append('<tr><th>Species</th><th>Stations</th>'
                             '<th>Observations</th><th>Period</th></tr>')
                    for sp in obs_stats['species_stats']:
                        yr = ''
                        if sp.get('year_start') and sp.get('year_end'):
                            yr = f"{sp['year_start']}&ndash;{sp['year_end']}"
                        H.append(
                            f'<tr><td><span class="badge badge-primary">'
                            f'{sp["species"]}</span></td>'
                            f'<td class="num">{sp["n_stations"]}</td>'
                            f'<td class="num">{sp["n_observations"]:,}</td>'
                            f'<td>{yr}</td></tr>')
                    H.append('</table></div></div>')

                # Species with no data
                if _no_data:
                    _nd_list = ', '.join(f'<code>{s}</code>' for s in _no_data)
                    if _obs_source == "user_csv":
                        H.append(
                            f'<div class="highlight-box" '
                            f'style="border-left-color:var(--accent);">'
                            f'The following model species have <strong>no '
                            f'observations</strong> in the CSV: {_nd_list}.'
                            f'</div>')
                    else:
                        _buf_km = obs_stats.get('buffer_km', '?')
                        H.append(
                            f'<div class="highlight-box" '
                            f'style="border-left-color:var(--accent);">'
                            f'The following species exist in the GRQA database '
                            f'but <strong>no stations were found</strong> '
                            f'within {_buf_km} km: {_nd_list}. '
                            f'Try increasing <code>grqa_buffer_km</code>.'
                            f'</div>')

                # Unmapped species
                if _unmapped:
                    _um_list = ', '.join(f'<code>{s}</code>' for s in _unmapped)
                    if _obs_source == "user_csv":
                        H.append(
                            f'<div class="highlight-box" '
                            f'style="border-left-color:var(--text-light);">'
                            f'The following species in the CSV are <strong>not '
                            f'in the model configuration</strong>: '
                            f'{_um_list}.</div>')
                    else:
                        H.append(
                            f'<div class="highlight-box" '
                            f'style="border-left-color:var(--text-light);">'
                            f'The following model species are <strong>not '
                            f'available</strong> in the GRQA database: '
                            f'{_um_list}.</div>')

            # GRQA-specific extras (stoichiometric conversion + parameter list)
            if _obs_source != "user_csv":
                H.append(
                    '<details style="margin-top:1rem;">'
                    '<summary style="cursor:pointer;font-weight:600;'
                    'color:var(--primary);font-size:0.95rem;">'
                    'Stoichiometric Conversion: GRQA &rarr; Model Units'
                    '</summary>'
                    '<div class="card primary" style="margin-top:0.5rem;">'
                    '<p style="margin-bottom:0.5rem;">GRQA reports '
                    'concentrations of the full ion '
                    '(e.g.&nbsp;mg&#8209;NO<sub>3</sub>/L), while '
                    'the model uses the <strong>&ldquo;as element&rdquo;'
                    '</strong> convention '
                    '(e.g.&nbsp;NO3&#8209;N in mg&#8209;N/L). The '
                    'conversion uses stoichiometry:</p>'
                    '<p style="text-align:center;font-size:0.95rem;'
                    'margin:0.5rem 0;">'
                    '<code>model_value = GRQA_value &times; MW(element) / '
                    'MW(ion)</code></p>'
                    '<div class="table-wrap"><table>'
                    '<tr><th>GRQA</th><th>Model</th>'
                    '<th>Factor</th><th>Formula</th></tr>')
                for grqa_code, factor in _GRQA_TO_MODEL_FACTOR.items():
                    elem = {'NO3': 'N', 'NH4': 'N', 'NO2': 'N',
                            'PO4': 'P', 'SO4': 'S'}[grqa_code]
                    model_name = f'{grqa_code}-{elem}'
                    H.append(
                        f'<tr><td><code>{grqa_code}</code> '
                        f'(mg&#8209;{grqa_code}/L)</td>'
                        f'<td><code>{model_name}</code> '
                        f'(mg&#8209;{elem}/L)</td>'
                        f'<td class="num">{factor:.4f}</td>'
                        f'<td>MW({elem}) / MW({grqa_code})</td></tr>')
                H.append('</table></div></div></details>')

                # Collapsible table of all GRQA parameters available
                _grqa_groups = {}
                for code, name in _GRQA_AVAILABLE_PARAMS.items():
                    if code in ('TN', 'DN', 'NO3', 'NO2', 'NH4', 'NO3_NO2',
                                'DIN', 'DON', 'TKN', 'PN'):
                        grp = 'Nitrogen'
                    elif code in ('TP', 'DP', 'DIP', 'PO4', 'DOP', 'PP'):
                        grp = 'Phosphorus'
                    elif code in ('DO', 'DO_sat', 'BOD', 'BOD5', 'COD',
                                  'CODMn'):
                        grp = 'Oxygen'
                    elif code in ('DOC', 'TOC', 'TC', 'TIC', 'DIC', 'POC',
                                  'PC'):
                        grp = 'Carbon'
                    elif code in ('TSS', 'TDS', 'SS', 'Turbidity'):
                        grp = 'Sediment'
                    elif code in ('Temp', 'pH', 'EC', 'SC'):
                        grp = 'Physical'
                    else:
                        grp = 'Other'
                    _grqa_groups.setdefault(grp, []).append((code, name))

                H.append(
                    '<details style="margin-top:1rem;">'
                    '<summary style="cursor:pointer;font-weight:600;'
                    'color:var(--primary);font-size:0.95rem;">'
                    'All GRQA Parameters Available in the Archive '
                    f'({len(_GRQA_AVAILABLE_PARAMS)} parameters)</summary>'
                    '<div class="card teal" style="margin-top:0.5rem;">'
                    '<div class="table-wrap"><table>'
                    '<tr><th>Group</th><th>Code</th><th>Parameter</th></tr>')
                for grp in ('Nitrogen', 'Phosphorus', 'Oxygen', 'Carbon',
                            'Sediment', 'Physical', 'Other'):
                    params = _grqa_groups.get(grp, [])
                    for i, (code, name) in enumerate(params):
                        grp_cell = (f'<td rowspan="{len(params)}" '
                                    f'style="font-weight:600;">{grp}</td>'
                                    if i == 0 else '')
                        H.append(f'<tr>{grp_cell}<td><code>{code}</code></td>'
                                 f'<td>{name}</td></tr>')
                H.append('</table></div></div></details>')

            H.append('</div>')

    # --- SECTION: Configuration ---
    with _SectionGuard(H, 'config', 'Model Configuration'):
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
            # Pass obs availability only for BGC modules
            _obs_for_mod = (_species_obs_avail
                            if 'BGC' in mod_name.upper() else None)
            param_html = _render_module_params(
                mod_data, mod_name, output_dir=output_dir,
                species_obs_availability=_obs_for_mod)
            if not param_html:
                continue
            H.append(f'<details class="module-details">')
            H.append(f'<summary>{mod_label}: '
                     f'<span class="badge badge-primary">{mod_name}</span></summary>')
            H.append('<div class="module-content">')
            H.extend(param_html)
            H.append('</div></details>')

        H.append('</div>')  # close section

    # --- SECTION: Source/Sink Summary ---
    with _SectionGuard(H, 'sources', 'Source/Sink Setup'):
        if ss_summary:
            H.append('<div class="section" id="sources"><h2>Source/Sink Setup</h2>')

            # Metadata card
            if ss_metadata:
                H.append('<div class="card primary" style="margin-bottom:1rem">'
                         '<div class="table-wrap"><table>')
                H.append('<tr><th>Property</th><th>Value</th></tr>')
                H.append(f'<tr><td>Method</td><td><span class="badge badge-secondary">'
                         f'{ss_method}</span></td></tr>')
                if seasonal_loads_method:
                    H.append(f'<tr><td>Seasonal Distribution</td><td>'
                             f'<span class="badge badge-secondary">'
                             f'{seasonal_loads_method}</span></td></tr>')
                # Observation compartments
                if observation_compartments:
                    if isinstance(observation_compartments, str):
                        _obs_comp_str = observation_compartments
                    else:
                        _obs_comp_str = ', '.join(str(c) for c in observation_compartments)
                    H.append(f'<tr><td>Observation Compartments</td>'
                             f'<td><span class="badge badge-secondary">'
                             f'{_obs_comp_str}</span></td></tr>')
                else:
                    H.append(f'<tr><td>Observation Compartments</td>'
                             f'<td><span class="badge badge-secondary">'
                             f'All (no filter)</span></td></tr>')
                for mk, mv in ss_metadata.items():
                    H.append(f'<tr><td>{mk}</td><td>{mv}</td></tr>')
                H.append('</table></div></div>')

            # Proxy year warning (when Copernicus proxy data is used)
            if ss_metadata.get('Proxy_Years'):
                H.append(
                    '<div class="card" style="border-left:4px solid #f39c12;'
                    'background:rgba(243,156,18,.08);margin-bottom:1rem;'
                    'padding:.8rem 1rem">'
                    '<p style="color:#e67e22;font-weight:600;margin:0 0 .4rem 0">'
                    '&#x26a0; Proxy Year Warning</p>'
                    f'<p style="margin:0;font-size:.85rem">'
                    f'{_html_mod.escape(str(ss_metadata["Proxy_Years"]))}</p>'
                    '</div>')

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
    with _SectionGuard(H, 'metadata', 'Run Metadata'):
        H.append('<div class="section" id="metadata"><h2>Run Metadata</h2>')
        H.append('<div class="card primary"><div class="table-wrap"><table>')
        H.append('<tr><th>Property</th><th>Value</th></tr>')
        H.append(f'<tr><td>Output Directory</td>'
                 f'<td style="font-family:monospace;font-size:.8rem">{output_dir}</td></tr>')
        H.append(f'<tr><td>Report Generated</td><td>{now}</td></tr>')
        H.append(f'<tr><td>Solver</td><td>{solver}</td></tr>')
        H.append(f'<tr><td>Output Timestep</td><td>{timestep[0]} {timestep[1]}</td></tr>')
        H.append('</table></div></div></div>')

    # --- SECTION: Errors (when any errors occurred) ---
    if _any_errors:
        with _SectionGuard(H, 'consoleerrors', 'Errors'):
            H.append('<div class="section" id="consoleerrors">'
                     '<h2>Errors</h2>')
            H.append('<div class="card" style="border-left:4px solid #e74c3c;'
                     'background:rgba(231,76,60,.06)">')
            H.append('<p style="color:#e74c3c;font-weight:600;margin:0 0 .8rem 0">'
                     '&#x26a0; Errors were encountered. '
                     'Please fix the issues and re-run your copy of <code>model_config_template.py</code> '
                     'to generate all the necessary input files.</p>')
            # Show per-section errors
            if _section_errors:
                for _sec_id, _sec_msg in _section_errors.items():
                    _safe_msg = _html_mod.escape(str(_sec_msg))
                    H.append(f'<div style="margin:.6rem 0;padding:.5rem .8rem;'
                             f'background:var(--glass);border-radius:6px;'
                             f'font-size:.85rem">'
                             f'<strong>{_html_mod.escape(_sec_id)}:</strong> '
                             f'{_safe_msg}</div>')
            H.append('</div></div>')

    # --- SECTION: Next Steps (skipped when errors occurred) ---
    if _any_errors:
        pass  # Next Steps omitted — errors must be fixed first
    else:
      with _SectionGuard(H, 'nextsteps', 'Next Steps'):
        # Helper: wrap a code snippet in a <pre> with a copy button.
        # The raw text is stored in a hidden <template> tag to avoid any
        # escaping issues with quotes / backslashes inside onclick attributes.
        _pre_style = ('position:relative;background:var(--glass);padding:.8rem 2.5rem .8rem .8rem;'
                      'border-radius:6px;overflow-x:auto;font-size:.82rem')
        _btn_style = ('position:absolute;top:.4rem;right:.4rem;background:var(--primary);color:#fff;'
                      'border:none;border-radius:4px;padding:2px 8px;font-size:.7rem;cursor:pointer')
        _copy_id_counter = [0]

        def _code_block(code_text):
            """Return HTML for a <pre> block with a copy-to-clipboard button."""
            import html as _html
            safe = _html.escape(code_text)
            cid = f'_cb{_copy_id_counter[0]}'
            _copy_id_counter[0] += 1
            return (f'<div style="position:relative">'
                    f'<pre id="{cid}" style="{_pre_style}">{safe}</pre>'
                    f'<button style="{_btn_style}" '
                    f'onclick="var t=document.getElementById(\'{cid}\').textContent;'
                    f'navigator.clipboard.writeText(t);'
                    f'this.textContent=\'Copied!\';'
                    f'setTimeout(function(){{this.textContent=\'Copy\'}}.bind(this),1500)">'
                    f'Copy</button></div>')

        # --- OS detection (early, so all subsequent code can use these) ---
        _os_platform = platform.system()   # 'Darwin', 'Linux', or 'Windows'
        _is_windows = (_os_platform == 'Windows')

        if _os_platform == 'Darwin':
            _os_label = 'macOS (Darwin)'
        elif _os_platform == 'Linux':
            _os_label = 'Linux'
        elif _os_platform == 'Windows':
            _os_label = 'Windows'
        else:
            _os_label = _os_platform

        # Python binary: Windows typically uses "python", Unix uses "python3"
        _py_bin = 'python' if _is_windows else 'python3'

        # cd command: Windows CMD needs /d flag for cross-drive navigation
        def _cd_cmd(path):
            return f'cd /d {path}' if _is_windows else f'cd {path}'

        # file:// URI helpers
        def _file_uri(path):
            """Build a file:// URI from an absolute path."""
            if _is_windows:
                return 'file:///' + path.replace('\\', '/')
            return 'file://' + path

        _file_uri_prefix = 'file:///' if _is_windows else 'file://'

        H.append('<div class="section visible" id="nextsteps"><h2>Next Steps</h2>')
        H.append('<div class="card primary">')
        H.append('<p>If you are happy with this configuration, follow these steps '
                 'to run the model using Docker. '
                 'Just copy-paste the code snippets below into your terminal!</p>')

        # OS badge indicating which platform the snippets target
        _badge_bg = '#0078d4' if _is_windows else '#2d7d46'
        H.append(
            f'<p style="margin-top:.5rem;font-size:.82rem;">'
            f'<span style="display:inline-block;background:{_badge_bg};color:#fff;'
            f'padding:2px 10px;border-radius:12px;font-size:.78rem;font-weight:600;">'
            f'&#x1F4BB; {_os_label}'
            f'</span> '
            f'These snippets are generated for <strong>{_os_label}</strong>. '
            f'If you are running on a different OS, some commands may need adjustment.</p>'
        )

        # Step 1: Start Docker container
        _script_dir = os.path.dirname(os.path.abspath(__file__))
        _containers_dir = os.path.normpath(
            os.path.join(_script_dir, '..', '..', '..', 'containers'))
        _step1_cmd = f'{_cd_cmd(_containers_dir)}\ndocker compose up -d'
        H.append('<h3 style="margin-top:1rem">1. Start the Docker container</h3>')
        H.append('<p>If the container is not already running:</p>')
        H.append(_code_block(_step1_cmd))

        # Step 2: Run the model
        H.append('<h3 style="margin-top:1rem">2. Run the model</h3>')

        # Build the docker exec command with resolved container paths
        _docker_cmd = None
        if executable_path and file_manager_path:
            try:
                import Gen_Input_Driver as _gJSON
                _dc_path = os.path.normpath(
                    os.path.join(_script_dir, '..', '..', '..', 'containers',
                                 'docker-compose.yml'))
                _host_root, _cont_root = _gJSON._parse_docker_volume_mount(_dc_path)
                if _host_root and _cont_root:
                    _cont_exec = _gJSON._correct_path_for_docker(
                        os.path.abspath(executable_path), _host_root, _cont_root)
                    _cont_fm = _gJSON._correct_path_for_docker(
                        os.path.abspath(file_manager_path), _host_root, _cont_root)
                    _abs_out = os.path.abspath(output_dir)
                    if not _abs_out.endswith('/'):
                        _abs_out += '/'
                    _cont_wd = _gJSON._correct_path_for_docker(
                        _abs_out, _host_root, _cont_root).rstrip('/')
                    # SUMMA requires '-m' flag before the file manager path;
                    # mizuRoute takes it as a positional argument.
                    _fm_arg = f'-m {_cont_fm}' if hostmodel.lower() == 'summa' else _cont_fm
                    # SUMMA-OpenWQ is serial (no MPI); force -np 1
                    _np = 1 if hostmodel.lower() == 'summa' else mpi_np
                    _docker_cmd = (
                        f'docker exec {docker_container_name} /bin/bash -c '
                        f'"export HDF5_USE_FILE_LOCKING=FALSE && '
                        f'cd {_cont_wd} && '
                        f'mpirun --allow-run-as-root -np {_np} '
                        f'{_cont_exec} '
                        f'{_fm_arg}"'
                    )
            except Exception as _e:
                print(f"  WARNING: Could not resolve Docker paths for report: {_e}")

        if not _docker_cmd and executable_path and file_manager_path:
            # Fallback: show host paths (user will need to convert manually)
            _fm_flag = '-m ' if hostmodel.lower() == 'summa' else ''
            _np = 1 if hostmodel.lower() == 'summa' else mpi_np
            _docker_cmd = (
                f'docker exec {docker_container_name} /bin/bash -c '
                f'"export HDF5_USE_FILE_LOCKING=FALSE && '
                f'cd <container_path_to:{os.path.abspath(output_dir)}> && '
                f'mpirun --allow-run-as-root -np {_np} '
                f'<container_path_to:{os.path.abspath(executable_path)}> '
                f'{_fm_flag}<container_path_to:{os.path.abspath(file_manager_path)}>"'
            )

        if _docker_cmd:
            H.append(_code_block(_docker_cmd))
        else:
            H.append('<p style="color:var(--muted);font-style:italic">'
                     'Could not build Docker command &mdash; executable_path '
                     'or file_manager_path not provided.</p>')

        # Step 3: Check outputs
        H.append('<h3 style="margin-top:1rem">3. Check outputs</h3>')
        _h5_dir = os.path.join(output_dir, 'openwq_out', 'HDF5')
        H.append('<p>After a successful run, HDF5 results will be saved to:</p>')
        H.append(_code_block(_cd_cmd(_h5_dir)))

        # Step 4: Read & visualize results
        _read_outputs_dir = os.path.normpath(
            os.path.join(_script_dir, '..', '..', '2_Read_Outputs'))
        _hdf5_lib_dir = os.path.join(_read_outputs_dir, 'hdf5_support_lib')
        _abs_output_dir = os.path.abspath(output_dir)
        _supporting_scripts_dir = os.path.normpath(
            os.path.join(_script_dir, '..', '..'))

        # Build species list string
        _species_list = chemical_species if isinstance(chemical_species, (list, tuple)) \
            else [chemical_species]
        _species_str_all = ', '.join(f'"{s}"' for s in _species_list)
        # Default snippet shows the first two species (matching default checkbox state)
        _default_n = min(2, len(_species_list))
        _species_str = ', '.join(f'"{s}"' for s in _species_list[:_default_n]) \
            if _species_list else _species_str_all

        # Compartment names — use hostmodel-specific list so users see all
        # valid compartments, not just the ones configured for output.
        _HOSTMODEL_COMPARTMENTS = {
            "mizuroute": ["RIVER_NETWORK_REACHES"],
            "summa": ["SCALARCANOPYWAT", "ILAYERVOLFRACWAT_SNOW", "RUNOFF",
                       "ILAYERVOLFRACWAT_SOIL", "SCALARAQUIFER"],
        }
        _cmp_names = _HOSTMODEL_COMPARTMENTS.get(
            hostmodel.lower(),
            list(compartments_and_cells.keys())
                if isinstance(compartments_and_cells, dict)
                else ['RIVER_NETWORK_REACHES']
        )
        # Default selection: compartments that were configured for output
        _cmp_configured = set(compartments_and_cells.keys()) \
            if isinstance(compartments_and_cells, dict) else set(_cmp_names)
        _cmp_default = [c for c in _cmp_names if c in _cmp_configured]
        if not _cmp_default:
            _cmp_default = _cmp_names[:1]
        _cmp_str = ', '.join(f'"{c}"' for c in _cmp_default)

        # Shapefile info for mapping — hostmodel-dependent defaults.
        # NOTE: _shp_key is the column in the RIVER-NETWORK shapefile used by
        # the 3D viewer (it must match the host-model IDs in the HDF5 output).
        # _basin_mapping_key is the column in the BASIN shapefile used for
        # SUMMA basin filtering. The two shapefiles can have different schemas
        # (e.g. TauDEM-derived river networks have LINKNO, not HRU_ID), so we
        # inspect each one independently below.
        _shp_path = river_network_shapefile or ''
        _basin_shp_path = basin_shapefile or ''

        def _detect_shp_mapping_key(shp_path, candidates, default_key, label):
            """Pick the first column from `candidates` that exists in the shapefile."""
            if not shp_path or not os.path.isfile(shp_path):
                return default_key
            try:
                import fiona
                with fiona.open(shp_path) as _src:
                    _props = list(_src.schema['properties'].keys())
                for _cand in candidates:
                    if _cand in _props:
                        print(f"  {label} mapping key: {_cand} "
                              f"(available: {_props})")
                        return _cand
                print(f"  WARNING: No standard mapping key found in {label} "
                      f"shapefile. Available: {_props}. Falling back to "
                      f"'{default_key}' — you may need to override "
                      f"shpfile_info['mapping_key'] manually.")
                return default_key
            except Exception as _e:
                print(f"  WARNING: Could not inspect {label} shapefile: {_e}")
                return default_key

        if hostmodel.lower() == 'summa':
            _h5_mapping_key = 'hruId'
            # Order matters: prefer exact SUMMA conventions before fallbacks
            _rn_candidates = ['hruId', 'HRU_ID', 'hru_id', 'HruId',
                              'gruId', 'GRU_ID', 'gru_id',
                              'LINKNO', 'COMID', 'reachID', 'SegId']
            _basin_candidates = ['HRU_ID', 'hruId', 'hru_id', 'HruId',
                                 'GRU_ID', 'gruId', 'gru_id']
            _shp_key = _detect_shp_mapping_key(
                river_network_shapefile, _rn_candidates, 'HRU_ID',
                'River network')
            _basin_mapping_key = _detect_shp_mapping_key(
                basin_shapefile, _basin_candidates, 'HRU_ID',
                'Basin')
        else:
            _h5_mapping_key = 'reachID'
            _rn_candidates = ['reachID', 'REACHID', 'reach_id', 'ReachID',
                              'SegId', 'segId', 'seg_id', 'COMID', 'LINKNO']
            _shp_key = _detect_shp_mapping_key(
                river_network_shapefile, _rn_candidates, 'SegId',
                'River network')
            _basin_mapping_key = ''
        _feature_label = openwq_h5_mapping_key or _shp_key

        # --- Venv activation (OS-aware, uses _is_windows from earlier block) ---
        if _is_windows:
            _venv_activate = (
                f'"{_supporting_scripts_dir}\\.venv\\Scripts\\activate.bat" 2>NUL || '
                f'"{_supporting_scripts_dir}\\venv\\Scripts\\activate.bat"'
            )
        else:
            _venv_activate = (
                f'source "{_supporting_scripts_dir}/.venv/bin/activate" 2>/dev/null || '
                f'source "{_supporting_scripts_dir}/venv/bin/activate" 2>/dev/null || true'
            )

        # Common Python preamble builder: imports + HDF5 reading
        _sed_flag = ts_module_name.upper() != "NONE"

        def _py_path(p):
            """Normalize a filesystem path for safe embedding in Python code.

            Python accepts forward slashes on all platforms, so converting
            backslashes avoids Windows paths like ``C:\\Users\\foo`` being
            misinterpreted as escape sequences (``\\U``, ``\\f``, etc.)
            when embedded inside string literals.
            """
            return p.replace('\\', '/')

        def _python_preamble(chem_spec_str=None):
            """Build the preamble with a specific chemSpec list."""
            cs = chem_spec_str if chem_spec_str is not None else _species_str
            _lib_dir_safe = _py_path(_hdf5_lib_dir)
            _results_dir_safe = _py_path(
                os.path.join(_abs_output_dir, "openwq_out"))
            return (
                f'import sys\n'
                f'sys.path.insert(0, "{_lib_dir_safe}")\n'
                f'import Read_h5_driver as h5_rlib\n'
                f'\n'
                f'openwq_results = h5_rlib.Read_h5_driver(\n'
                f'    openwq_info={{\n'
                f'        "path_to_results": "{_results_dir_safe}",\n'
                f'        "mapping_key": "{_h5_mapping_key}"\n'
                f'    }},\n'
                f'    output_format="HDF5",\n'
                f'    debugmode=True,\n'
                f'    cmp=[{_cmp_str}],\n'
                f'    space_elem="all",\n'
                f'    chemSpec=[{cs}],\n'
                f'    chemUnits="{units}",\n'
                f'    noDataFlag=-9999,\n'
                f'    sediment_as_well={_sed_flag}\n'
                f')'
            )

        def _terminal_snippet(python_body):
            """Wrap Python code into a self-contained terminal command.

            On macOS/Linux: uses a bash heredoc.
            On Windows: writes to a temporary .py file via PowerShell, then executes it.
            """
            if _is_windows:
                _tmp = '%TEMP%\\_openwq_snippet.py'
                return (
                    f'{_venv_activate}\n'
                    f"powershell -Command \"Set-Content -Path '{_tmp}' -Value @'\n"
                    f'{python_body}\n'
                    f"'@\"\n"
                    f'{_py_bin} {_tmp}\n'
                    f'del {_tmp}'
                )
            else:
                return (
                    f'{_venv_activate}\n'
                    f"{_py_bin} << 'PYEOF'\n"
                    f'{python_body}\n'
                    f'PYEOF'
                )

        H.append('<h3 style="margin-top:1rem">4. Read &amp; visualize results</h3>')
        H.append('<p style="font-size:.85rem;color:var(--muted)">'
                 'The snippets below are <strong>self-contained</strong> &mdash; '
                 'copy each into your terminal and it will run as-is '
                 '(activates the venv, reads HDF5, and produces output).</p>')
        H.append('<p style="font-size:.82rem;color:var(--accent);margin:.4rem 0 .5rem">'
                 '&#x26A0; Selecting multiple combinations of species, compartments, '
                 'and debug mode may slow down rendering and increase memory usage. '
                 'For best performance, keep only a few active at a time.</p>')

        # --- Species checkboxes ---
        import html as _html_mod
        H.append('<p style="font-weight:600;margin-top:1rem">Select chemical species to visualize:</p>')
        H.append('<div id="speciesCbRow" style="display:flex;flex-wrap:wrap;gap:.5rem .8rem;'
                 'margin:.5rem 0 1rem;align-items:center">')
        for _i_sp, _sp in enumerate(_species_list):
            _sp_esc = _html_mod.escape(_sp)
            _chk = " checked" if _i_sp < 2 else ""
            H.append(
                f'<label style="display:inline-flex;align-items:center;gap:.3rem;'
                f'font-size:.85rem;cursor:pointer;padding:.25rem .5rem;'
                f'border:1px solid var(--border);border-radius:6px;'
                f'background:var(--surface);transition:border-color .15s">'
                f'<input type="checkbox" class="species-cb" '
                f'data-species="{_sp_esc}"{_chk}> {_sp_esc}</label>')
        H.append('<a href="#" id="speciesToggleAll" '
                 'style="font-size:.78rem;margin-left:.5rem;color:var(--primary)">'
                 'Deselect all</a>')
        H.append('</div>')

        # --- Compartment checkboxes ---
        H.append('<p style="font-weight:600;margin-top:1rem">Select compartments to read:</p>')
        H.append('<p style="font-size:.82rem;color:var(--muted);margin:0 0 .5rem">'
                 f'Available compartments for <strong>{hostmodel}</strong>. '
                 'Only checked compartments are included in the code snippets.</p>')
        H.append('<div id="cmpCbRow" style="display:flex;flex-wrap:wrap;gap:.5rem .8rem;'
                 'margin:.5rem 0 1rem;align-items:center">')
        for _cmp in _cmp_names:
            _cmp_esc = _html_mod.escape(_cmp)
            _chk = " checked" if _cmp in _cmp_default else ""
            H.append(
                f'<label style="display:inline-flex;align-items:center;gap:.3rem;'
                f'font-size:.85rem;cursor:pointer;padding:.25rem .5rem;'
                f'border:1px solid var(--border);border-radius:6px;'
                f'background:var(--surface);transition:border-color .15s">'
                f'<input type="checkbox" class="cmp-cb" '
                f'data-cmp="{_cmp_esc}"{_chk}> {_cmp_esc}</label>')
        H.append('<a href="#" id="cmpToggleAll" '
                 'style="font-size:.78rem;margin-left:.5rem;color:var(--primary)">'
                 'Select all</a>')
        H.append('</div>')

        # --- Debug mode checkbox ---
        H.append('<p style="font-weight:600;margin-top:1rem">Debug mode:</p>')
        H.append('<p style="font-size:.82rem;color:var(--muted);margin:0 0 .5rem">'
                 'When enabled, the Read driver loads derivative / source outputs '
                 '(<code>dC/dt chem</code>, <code>dC/dt transp</code>, '
                 '<code>SS</code>, <code>EWF</code>, <code>IC</code>) and the '
                 'Plot driver overlays them as dashed lines on each graph. '
                 'A per-plot <em>Debug</em> toggle lets you show / hide them.</p>')
        H.append('<div id="debugCbRow" style="display:flex;flex-wrap:wrap;gap:.5rem .8rem;'
                 'margin:.5rem 0 1rem;align-items:center">')
        H.append(
            '<label style="display:inline-flex;align-items:center;gap:.3rem;'
            'font-size:.85rem;cursor:pointer;padding:.25rem .5rem;'
            'border:1px solid var(--border);border-radius:6px;'
            'background:var(--surface);transition:border-color .15s">'
            '<input type="checkbox" id="debugModeCb" checked> '
            'Include debug outputs</label>')
        H.append('</div>')

        # --- 4a: Interactive 3D Spatial Viewer (WebGL) ---
        # For SUMMA, results are HRU/GRU based, so we drive the viewer with
        # basin polygons (the catchment shapefile). For mizuRoute, results are
        # per river reach, so we keep using the river-network polylines.
        _webgl_out_dir = os.path.join(_abs_output_dir, "openwq_out",
                                       "openwq_webgl_viewer")
        _webgl_out_dir_safe = _py_path(_webgl_out_dir)
        _srv_dir_safe = _py_path(os.path.join(_abs_output_dir, "openwq_out"))
        _shp_path_safe = _py_path(river_network_shapefile) if river_network_shapefile else ''
        _basin_shp_path_safe = _py_path(basin_shapefile) if basin_shapefile else ''

        # Pick the right shapefile + mapping key for the spatial viewer.
        if hostmodel.lower() == 'summa':
            if _basin_shp_path_safe:
                _viewer_shp_path_safe = _basin_shp_path_safe
                _viewer_shp_key       = _basin_mapping_key
                _viewer_geom_label    = 'HRU/GRU polygons'
            else:
                # No basin shapefile — fall back to river network with a warning
                _viewer_shp_path_safe = _shp_path_safe
                _viewer_shp_key       = _shp_key
                _viewer_geom_label    = 'river network (no basin shapefile available)'
                print("  WARNING: SUMMA host but no basin_shapefile provided -- "
                      "falling back to river_network_shapefile for the spatial "
                      "viewer (HRU polygons would be preferred).")
            _viewer_title_html = 'Generate interactive HRU/GRU viewer:'
        else:
            _viewer_shp_path_safe = _shp_path_safe
            _viewer_shp_key       = _shp_key
            _viewer_geom_label    = 'river network'
            _viewer_title_html    = 'Generate interactive 3D river viewer:'

        # Try to find the host-model output directory for hydromodel_info
        _hm_out_dir_safe = ''
        _hm_prefix = ''
        if file_manager_path and executable_path:
            try:
                _hm_dir, _hm_pfx = _parse_hostmodel_output_info(
                    file_manager_path, hostmodel,
                    os.path.dirname(os.path.abspath(executable_path)))
                if _hm_dir:
                    _hm_out_dir_safe = _py_path(os.path.abspath(_hm_dir))
                    _hm_prefix = _hm_pfx or ''
            except Exception:
                pass

        _webgl_body = (
            f'{_python_preamble()}\n'
            f'\n'
            f'# --- Build 3D WebGL viewer ---\n'
            f'# Geometry: {_viewer_geom_label}\n'
            f'import glob, os\n'
            f'\n'
            f'shpfile_info = {{\n'
            f'    "path_to_shp": "{_viewer_shp_path_safe}",\n'
            f'    "mapping_key": "{_viewer_shp_key}"\n'
            f'}}\n'
            f'\n'
        )

        # Host-specific flow-variable candidates (used at SNIPPET runtime to
        # pick whichever variable actually exists in the netCDF).
        if hostmodel.lower() == 'summa':
            _flow_var_candidates = [
                'averageRoutedRunoff', 'basRunoff', 'scalarTotalRunoff',
                'scalarSurfaceRunoff', 'scalarSoilBaseflow',
                'scalarAquiferBaseflow', 'scalarSoilDrainage']
        else:  # mizuroute (default routing methods)
            _flow_var_candidates = [
                'DWroutedRunoff', 'IRFroutedRunoff', 'KWTroutedRunoff',
                'KWroutedRunoff', 'MCroutedRunoff', 'instRunoff']
        _flow_cands_str = '[' + ', '.join(f'"{v}"' for v in _flow_var_candidates) + ']'

        # Add hydromodel_info: auto-detect .nc file AND flow variable name.
        # If no usable flow variable is found, hydromodel_info is set to None
        # so the WebGL driver runs without flow animation (it handles None).
        if _hm_out_dir_safe:
            _webgl_body += (
                f'# Auto-detect host-model output (netCDF)\n'
                f'_nc_files = sorted(glob.glob("{_hm_out_dir_safe}/{_hm_prefix}*.nc"))\n'
                f'if not _nc_files:\n'
                f'    print("WARNING: No .nc files found in {_hm_out_dir_safe}")\n'
                f'\n'
                f'# Auto-detect a flow variable in the netCDF\n'
                f'_flow_var = None\n'
                f'_flow_candidates = {_flow_cands_str}\n'
                f'if _nc_files:\n'
                f'    try:\n'
                f'        from netCDF4 import Dataset as _NCDS\n'
                f'        with _NCDS(_nc_files[0]) as _f:\n'
                f'            _available = list(_f.variables.keys())\n'
                f'        for _c in _flow_candidates:\n'
                f'            if _c in _available:\n'
                f'                _flow_var = _c\n'
                f'                break\n'
                f'        if _flow_var:\n'
                f'            print(f"Flow variable detected: {{_flow_var}}")\n'
                f'        else:\n'
                f'            print(f"WARNING: None of {{_flow_candidates}} found in {{_nc_files[0]}}")\n'
                f'            print(f"  Available variables: {{_available}}")\n'
                f'            print("  3D viewer will run without flow animation.")\n'
                f'    except Exception as _e:\n'
                f'        print(f"WARNING: Could not inspect netCDF for flow variable: {{_e}}")\n'
                f'\n'
                f'if _nc_files and _flow_var is not None:\n'
                f'    hydromodel_info = {{\n'
                f'        "path_to_results": _nc_files[0],\n'
                f'        "mapping_key": "{_h5_mapping_key}"\n'
                f'    }}\n'
                f'else:\n'
                f'    hydromodel_info = None    # WebGL driver skips flow when None\n'
                f'\n'
            )
        else:
            _webgl_body += (
                f'# TODO: set the path to the host-model netCDF output below\n'
                f'# (or leave hydromodel_info = None to skip flow animation)\n'
                f'hydromodel_info = {{\n'
                f'    "path_to_results": "<path/to/hostmodel_output.nc>",\n'
                f'    "mapping_key": "{_h5_mapping_key}"\n'
                f'}}\n'
                f'_flow_var = "{_flow_var_candidates[0]}"   # adjust to match the variable in your netCDF\n'
                f'\n'
            )

        _webgl_body += (
            f'import WebGL_h5_driver as h5_wlib\n'
            f'\n'
            f'_result = h5_wlib.WebGL_h5_driver(\n'
            f'    what2map="openwq",\n'
            f'    hostmodel="{hostmodel}",\n'
            f'    shpfile_info=shpfile_info,\n'
            f'    openwq_results=openwq_results,\n'
            f'    chemSpec=[{_species_str}],\n'
            f'    sediment_as_well={_sed_flag},\n'
            f'    hydromodel_info=hydromodel_info,\n'
            f'    hydromodel_var2print=_flow_var,\n'
            f'    output_dir="{_webgl_out_dir_safe}",\n'
            f'    timeframes=50,\n'
            f'    n_particles=65536,\n'
            f'    river_width_cells=3,\n'
            f'    grid_resolution=None,\n'
            f'    satellite_resolution=2\n'
            f')\n'
            f'\n'
            f'if _result:\n'
            f'    # Start local HTTP server and open the 3D viewer.\n'
            f'    # Walk a small port range so the snippet works even if a\n'
            f'    # previous run is still holding 8080, then fall back to an\n'
            f'    # OS-assigned free port if all preferred ones are taken.\n'
            f'    import threading, http.server, functools, time, webbrowser\n'
            f'    _handler = functools.partial(\n'
            f'        http.server.SimpleHTTPRequestHandler,\n'
            f'        directory="{_srv_dir_safe}")\n'
            f'    _httpd = None\n'
            f'    _port = None\n'
            f'    for _try_port in range(8080, 8100):\n'
            f'        try:\n'
            f'            _httpd = http.server.HTTPServer(("localhost", _try_port), _handler)\n'
            f'            _port = _try_port\n'
            f'            break\n'
            f'        except OSError:\n'
            f'            continue\n'
            f'    if _httpd is None:\n'
            f'        _httpd = http.server.HTTPServer(("localhost", 0), _handler)\n'
            f'        _port = _httpd.server_address[1]\n'
            f'        print(f"All preferred ports 8080-8099 in use; using OS-assigned port {{_port}}.")\n'
            f'    _url = f"http://localhost:{{_port}}/openwq_webgl_viewer/index.html"\n'
            f'    threading.Thread(target=_httpd.serve_forever, daemon=True).start()\n'
            f'    print(f"Serving at http://localhost:{{_port}}")\n'
            f'    webbrowser.open(_url)\n'
            f'\n'
            f'    print("Press Ctrl+C to stop the server...")\n'
            f'    try:\n'
            f'        while True: time.sleep(1)\n'
            f'    except KeyboardInterrupt:\n'
            f'        _httpd.shutdown()\n'
            f'        print("\\nServer stopped.")\n'
            f'else:\n'
            f'    print("\\nNo 3D viewer was generated.")\n'
            f'    print("Check that the model has been run and HDF5 output files exist.")'
        )

        # --- 4b: Interactive time-series plots (all species) ---
        _report_name = "openwq_simulation_report.html"
        _out_all_html = os.path.join(_abs_output_dir, "openwq_out",
                                     _report_name)
        _out_all_html_safe = _py_path(_out_all_html)
        # Build optional observation parameter for the plotting snippet
        _obs_plot_param = ''
        if _obs_source == "grqa" and obs_stats and obs_stats.get('n_stations', 0) > 0:
            _grqa_clip_dir = os.path.join(_abs_output_dir, 'openwq_in',
                                          'grqa_clipped_data')
            _obs_plot_param = (
                f'    observation_dir="{_py_path(_grqa_clip_dir)}",\n')
        elif _obs_source == "user_csv" and user_observation_csv:
            _obs_plot_param = (
                f'    observation_csv="{_py_path(os.path.abspath(user_observation_csv))}",\n')
        # Add compartment filter for observations
        if observation_compartments:
            _obs_plot_param += (
                f'    observation_compartments={repr(observation_compartments)},\n')

        _plot_all_body = (
            f'{_python_preamble()}\n'
            f'\n'
            f'import Plot_h5_driver as h5_plib\n'
            f'\n'
            f'_result = h5_plib.Plot_h5_driver(\n'
            f'    what2map="openwq",\n'
            f'    hostmodel="{hostmodel}",\n'
            f'    mapping_key_values="all",\n'
            f'    openwq_results=openwq_results,\n'
            f'    chemSpec=[{_species_str}],\n'
            f'    debugmode=True,\n'
            f'    output_path="{_out_all_html_safe}",\n'
            f'    river_network_shp="{_shp_path_safe}",\n'
            f'    basin_shapefile="{_basin_shp_path_safe}",\n'
            f'    basin_mapping_key="{_basin_mapping_key}",\n'
            f'    mapping_key="{_shp_key}",\n'
            f'    feature_label="{_feature_label}",\n'
            f'    separator="{plot_separator}",\n'
            f'{_obs_plot_param}'
            f')\n'
            f'\n'
            f'if _result:\n'
            f'    import webbrowser\n'
            f'    webbrowser.open("{_file_uri(_out_all_html)}")\n'
            f'else:\n'
            f'    print("\\nNo plots were generated.")\n'
            f'    print("Check that the model has been run and HDF5 output files exist.")'
        )

        # Two-column layout: 4a (left) and 4b (right) side by side
        H.append('<div style="display:flex;gap:1.5rem;margin-top:1rem;align-items:flex-start">')
        H.append('<div style="flex:1;min-width:0">')
        H.append(f'<p style="font-weight:600">{_viewer_title_html}</p>')
        _cb_id_4a = f'_cb{_copy_id_counter[0]}'  # ID before _code_block increments
        H.append(_code_block(_terminal_snippet(_webgl_body)))
        H.append('</div>')
        H.append('<div style="flex:1;min-width:0">')
        H.append('<p style="font-weight:600">Plot interactive '
                 'time series (all species):</p>')
        _cb_id_4b = f'_cb{_copy_id_counter[0]}'  # ID before _code_block increments
        H.append(_code_block(_terminal_snippet(_plot_all_body)))
        H.append('</div>')
        H.append('</div>')

        # JavaScript: update code snippets when species, compartment, or debug checkboxes change
        H.append(f"""<script>
(function(){{
  var preIds = ['{_cb_id_4a}', '{_cb_id_4b}'];
  var origTexts = preIds.map(function(id){{ return document.getElementById(id).textContent; }});
  var specCbs = document.querySelectorAll('.species-cb');
  var cmpCbs = document.querySelectorAll('.cmp-cb');
  var debugCb = document.getElementById('debugModeCb');
  var specToggle = document.getElementById('speciesToggleAll');
  var cmpToggle = document.getElementById('cmpToggleAll');
  var reSpec = /chemSpec=\\[.*?\\]/g;
  var reCmp = /cmp=\\[.*?\\]/g;
  var reDebug = /debugmode=(True|False)/g;

  function updateSnippets(){{
    var checkedSpec = [];
    specCbs.forEach(function(cb){{ if(cb.checked) checkedSpec.push('"'+cb.dataset.species+'"'); }});
    var newSpec = 'chemSpec=[' + checkedSpec.join(', ') + ']';

    var checkedCmp = [];
    cmpCbs.forEach(function(cb){{ if(cb.checked) checkedCmp.push('"'+cb.dataset.cmp+'"'); }});
    var newCmp = 'cmp=[' + checkedCmp.join(', ') + ']';

    var newDebug = 'debugmode=' + (debugCb && debugCb.checked ? 'True' : 'False');

    preIds.forEach(function(id, idx){{
      var el = document.getElementById(id);
      var txt = origTexts[idx]
        .replace(reSpec, newSpec)
        .replace(reCmp, newCmp)
        .replace(reDebug, newDebug);
      el.textContent = txt;
    }});

    if(specToggle){{
      specToggle.textContent = checkedSpec.length === specCbs.length ? 'Deselect all' : 'Select all';
    }}
    if(cmpToggle){{
      cmpToggle.textContent = checkedCmp.length === cmpCbs.length ? 'Deselect all' : 'Select all';
    }}
  }}

  specCbs.forEach(function(cb){{ cb.addEventListener('change', updateSnippets); }});
  cmpCbs.forEach(function(cb){{ cb.addEventListener('change', updateSnippets); }});
  if(debugCb){{ debugCb.addEventListener('change', updateSnippets); }}

  if(specToggle){{
    specToggle.addEventListener('click', function(e){{
      e.preventDefault();
      var allChecked = Array.from(specCbs).every(function(cb){{ return cb.checked; }});
      specCbs.forEach(function(cb){{ cb.checked = !allChecked; }});
      updateSnippets();
    }});
  }}
  if(cmpToggle){{
    cmpToggle.addEventListener('click', function(e){{
      e.preventDefault();
      var allChecked = Array.from(cmpCbs).every(function(cb){{ return cb.checked; }});
      cmpCbs.forEach(function(cb){{ cb.checked = !allChecked; }});
      updateSnippets();
    }});
  }}
}})();
</script>""")

        H.append('</div></div>')

        H.append('</div>')  # container

    # Footer
    H.append(f"""<div class="footer">
<div class="logo">Open<span>WQ</span></div>
<p>Simulation report generated on {now}</p>
</div>""")

    H.append('</div>')  # main
    H.append('</div>')  # layout

    # =========================================================================
    # JAVASCRIPT
    # =========================================================================
    try:
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

        # Basin map — river network + basin + GRQA stations
        if (map_layers and map_center) or obs_geojson_str:
            layer_js_blocks = []
            overlay_entries = []
            fit_var = None

            # Shapefile layers (basin polygons + river network)
            # Interactive layers get popups; background layers are lighter
            # and togglable via the Leaflet overlay checkbox control.
            # Sort so background layers render first (underneath).
            _sorted_layers = sorted(map_layers, key=lambda l: l.get('interactive', True))
            bg_layer_vars = []  # track background layer vars for toggle
            for idx, layer in enumerate(_sorted_layers):
                var_name = f"lyr_{idx}"
                geojson_var = f"gj_{idx}"
                style = layer['style']
                label = layer['label'].replace("'", "\\'")
                layer_type = layer['type']
                is_interactive = layer.get('interactive', True)

                layer_js_blocks.append(f"  var {geojson_var} = {layer['geojson_str']};")

                color = style.get('color', '#0066cc')
                weight = style.get('weight', 2.5)
                opacity = style.get('opacity', 0.85)
                fill_color = style.get('fillColor', color)
                fill_opacity = style.get('fillOpacity', 0.2)

                if is_interactive:
                    # Interactive layer: full styling + click popups
                    layer_js_blocks.append(f"""  var {var_name} = L.geoJSON({geojson_var}, {{
    style: function(f){{ return {{color:'{color}',weight:{weight},opacity:{opacity},fillColor:'{fill_color}',fillOpacity:{fill_opacity}}}; }},
    onEachFeature: function(f,l){{
      var props = f.properties;
      var html = '<b>{label}</b><br>';
      for(var k in props){{ if(props[k]!==null && props[k]!=='') html += '<b>'+k+'</b>: '+props[k]+'<br>'; }}
      l.bindPopup(html);
    }}
  }}).addTo(map);""")
                else:
                    # Background reference layer: lighter, no popups, non-interactive
                    layer_js_blocks.append(f"""  var {var_name} = L.geoJSON({geojson_var}, {{
    style: function(f){{ return {{color:'{color}',weight:{weight},opacity:{opacity},fillColor:'{fill_color}',fillOpacity:{fill_opacity}}}; }},
    interactive: false
  }}).addTo(map);""")
                    bg_layer_vars.append((var_name, label))

                overlay_label = f'{label} (ref.)' if not is_interactive else label
                overlay_entries.append(f"'{overlay_label}': {var_name}")
                if fit_var is None:
                    fit_var = var_name

            # Observation stations layer (circle markers)
            if obs_geojson_str:
                _stn_popup_title = ('Observation Station' if _obs_source == 'user_csv'
                                    else 'GRQA Station')
                _stn_layer_label = (f'{_obs_label} Stations').replace("'", "\\'")
                layer_js_blocks.append(f"  var gj_obs = {obs_geojson_str};")
                layer_js_blocks.append(f"""  var lyr_obs = L.geoJSON(gj_obs, {{
    pointToLayer: function(f, latlng){{
      return L.circleMarker(latlng, {{
        radius: 7, fillColor: '#e63946', color: '#333',
        weight: 1.5, opacity: 0.9, fillOpacity: 0.85
      }});
    }},
    onEachFeature: function(f,l){{
      var props = f.properties;
      var html = '<b>{_stn_popup_title}</b><br>';
      for(var k in props){{ if(props[k]!==null && props[k]!=='') html += '<b>'+k+'</b>: '+props[k]+'<br>'; }}
      l.bindPopup(html);
    }}
  }}).addTo(map);""")
                overlay_entries.append(f"'{_stn_layer_label}': lyr_obs")
                if fit_var is None:
                    fit_var = "lyr_obs"

            overlay_obj = "{" + ", ".join(overlay_entries) + "}"
            center_js = _js(map_center) if map_center else "[0, 0]"
            zoom = 9 if map_center else 2

            H.append(f"""
(function(){{
  var map = L.map('basinMap',{{
    zoomControl:false,dragging:false,scrollWheelZoom:false,
    doubleClickZoom:false,touchZoom:false,boxZoom:false,keyboard:false
  }}).setView({center_js}, {zoom});
  var satTile = L.tileLayer('https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{{z}}/{{y}}/{{x}}',{{
    attribution:'Esri, Maxar, Earthstar Geographics', maxZoom:19}}).addTo(map);
  var lightTile = L.tileLayer('https://{{s}}.basemaps.cartocdn.com/light_all/{{z}}/{{x}}/{{y}}{{r}}.png',{{
    attribution:'&copy; OpenStreetMap &copy; CARTO', maxZoom:19
  }});
  var topoTile = L.tileLayer('https://{{s}}.tile.opentopomap.org/{{z}}/{{x}}/{{y}}.png',{{
    attribution:'OpenTopoMap', maxZoom:17}});
  // Lock / unlock control (top-left, locked by default)
  var mapLocked=true;
  var LockCtrl=L.Control.extend({{
    options:{{position:'topleft'}},
    onAdd:function(){{
      var btn=L.DomUtil.create('button','leaflet-bar');
      btn.style.cssText='width:34px;height:34px;background:#fff;border:none;'
        +'border-radius:4px;cursor:pointer;font-size:16px;line-height:34px;text-align:center;'
        +'box-shadow:0 1px 5px rgba(0,0,0,.3);';
      btn.title='Unlock map';
      btn.innerHTML='&#x1F512;';
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
  L.control.zoom({{position:'topleft'}}).addTo(map);
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
    except Exception as _e:
        print(f"  WARNING: JavaScript generation failed: {_e}")
        H.append('<script>console.error("OpenWQ report: JS generation failed");</script>')

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
        observation_data_source="grqa",
        grqa_local_data_path=None,
        grqa_buffer_km=10,
        user_observation_csv=None,
        observation_compartments=None,
        container_runtime="docker",
        mpi_np=2,
        docker_container_name="docker_openwq",
        executable_path=None,
        file_manager_path=None,
        config_errors=None,
        bgc_template_path=None,
        openwq_h5_mapping_key=None,
        seasonal_loads_method=None,
        plot_separator=' | ',
):
    """Generate an HTML simulation report (entry point for template).

    Wraps generate_simulation_report() with logging and error handling.

    Parameters:
        river_network_shapefile: path to river network shapefile
        basin_shapefile: path to basin/catchment shapefile
        observation_data_source: "grqa", "user_csv", or "skip"
        grqa_local_data_path: "auto"/path for GRQA
        grqa_buffer_km: buffer distance in km for GRQA station search
        user_observation_csv: path to user-provided CSV file

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
            observation_data_source=observation_data_source,
            grqa_local_data_path=grqa_local_data_path,
            grqa_buffer_km=grqa_buffer_km,
            user_observation_csv=user_observation_csv,
            observation_compartments=observation_compartments,
            container_runtime=container_runtime,
            mpi_np=mpi_np,
            docker_container_name=docker_container_name,
            executable_path=executable_path,
            file_manager_path=file_manager_path,
            config_errors=config_errors,
            bgc_template_path=bgc_template_path,
            openwq_h5_mapping_key=openwq_h5_mapping_key,
            seasonal_loads_method=seasonal_loads_method,
            plot_separator=plot_separator,
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
