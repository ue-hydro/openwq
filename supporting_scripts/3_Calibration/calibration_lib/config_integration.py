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
Config Integration Module
==========================

Loads model_config_template.py as a Python module and exposes its variables
for use by the calibration framework.  Also provides functions to generate
config files for individual calibration evaluations using Gen_Input_Driver.

Manifest handoff
----------------
When the user runs the *config* template (1_Model_Config/...), Gen_Report
writes ``openwq_baseline_manifest.json`` alongside the HTML report.  The
calibration side prefers values from that manifest over re-detecting them,
so the calibration runs use exactly the same spatial mapping keys, host-
model conventions, and shapefile paths the user saw in the config viewer.
The manifest is OPTIONAL — if it's not there, calibration falls back to
re-exec'ing the template the way it always has.
"""

import os
import sys
import copy
import json
import logging
from pathlib import Path
from typing import Dict, Any, Optional, List

logger = logging.getLogger(__name__)

# Manifest filename written by Gen_Report.generate_simulation_report.
_MANIFEST_FILENAME = "openwq_baseline_manifest.json"
# Where load_model_config stores the loaded manifest dict in the returned
# config so downstream consumers (get_observation_config, get_spatial_mapping)
# can find it without re-reading the disk.
_MANIFEST_CFG_KEY = "_baseline_manifest"

# Variables that should NOT be forwarded to Gen_Input_Driver
_EXCLUDED_VARS = {
    # Python builtins / imports
    'sys', 'os', 'io', 'webbrowser', '__builtins__', '__name__',
    '__doc__', '__file__', '__loader__', '__spec__', '__cached__',
    # model_config_template internal imports
    'gJSON_lib', 'uniform_param', '_generate_report',
    '_has_errors', '_config_error_msg', '_report_path', '_tb',
    '_excluded', '_config', '_basin_shp', '_h5_mapping_key',
    '_seasonal_method',
    # Bookkeeping entries stored on the returned config dict — never
    # forwarded to Gen_Input_Driver.
    _MANIFEST_CFG_KEY,
    '_model_config_path',
}


def _try_load_manifest(model_config: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    """Locate and load ``openwq_baseline_manifest.json`` for this config.

    The manifest is written next to the config HTML report under
    ``{dir2save_input_files}/openwq_baseline_manifest.json``.  We probe a
    couple of directories so this still works when the user has moved the
    output around.  Returns the parsed dict or ``None`` if no manifest is
    found / parseable.
    """
    candidates = []
    d = model_config.get('dir2save_input_files')
    if d:
        candidates.append(os.path.join(d, _MANIFEST_FILENAME))
    ex = model_config.get('executable_path')
    if ex:
        candidates.append(
            os.path.join(os.path.dirname(os.path.abspath(ex)),
                         _MANIFEST_FILENAME))
    # Deduplicate while preserving order
    seen = set()
    uniq = []
    for c in candidates:
        if c and c not in seen:
            seen.add(c)
            uniq.append(c)
    for path in uniq:
        if not os.path.isfile(path):
            continue
        try:
            with open(path, 'r', encoding='utf-8') as f:
                manifest = json.load(f)
            logger.info(f"Loaded baseline manifest: {path}")
            return manifest
        except Exception as e:
            logger.warning(f"Found manifest at {path} but failed to parse: {e}")
            return None
    logger.debug(
        f"No baseline manifest found (looked in: {uniq}). "
        f"Calibration will re-detect mapping keys."
    )
    return None


def load_model_config(model_config_path: str) -> Dict[str, Any]:
    """
    Import model_config_template.py and extract user-defined variables.

    Reads the file up to the "END OF USER CONFIGURATION" marker and
    executes only that portion in an isolated namespace.  This gives us
    all user-defined variables (paths, module settings, species, etc.)
    without triggering config generation or report opening.

    Parameters
    ----------
    model_config_path : str
        Absolute path to the user's copy of model_config_template.py.

    Returns
    -------
    Dict[str, Any]
        Dictionary of all user-defined configuration variables.

    Raises
    ------
    FileNotFoundError
        If model_config_path does not exist.
    RuntimeError
        If the file cannot be parsed or executed.
    """
    model_config_path = os.path.abspath(model_config_path)

    if not os.path.isfile(model_config_path):
        raise FileNotFoundError(
            f"Model config file not found: {model_config_path}"
        )

    logger.info(f"Loading model config from: {model_config_path}")

    with open(model_config_path, 'r') as f:
        content = f.read()

    # Find the "END OF USER CONFIGURATION" marker
    marker = "END OF USER CONFIGURATION"
    marker_pos = content.find(marker)
    if marker_pos > 0:
        # Find the beginning of the line containing the marker
        line_start = content.rfind('\n', 0, marker_pos) + 1
        config_content = content[:line_start]
    else:
        logger.warning(
            f"Could not find '{marker}' marker in {model_config_path}. "
            f"Using entire file content."
        )
        config_content = content

    # Set up sys.path so that imports in the config file work
    config_dir = os.path.dirname(model_config_path)
    old_path = sys.path[:]
    old_cwd = os.getcwd()

    namespace = {'__builtins__': __builtins__}

    try:
        sys.path.insert(0, config_dir)
        os.chdir(config_dir)
        exec(compile(config_content, model_config_path, 'exec'), namespace)
    except Exception as e:
        raise RuntimeError(
            f"Failed to load model config from {model_config_path}: {e}"
        ) from e
    finally:
        sys.path = old_path
        os.chdir(old_cwd)

    # Filter to user-defined variables
    config = {}
    for k, v in namespace.items():
        if k.startswith('_'):
            continue
        if k in _EXCLUDED_VARS:
            continue
        if callable(v) and k not in ('chemical_species',):
            # Skip functions/classes but keep list/dict values
            if not isinstance(v, (list, dict, tuple, str, int, float, bool, type(None))):
                continue
        config[k] = v

    # Remember the source .py path so downstream consumers (e.g. the
    # calibration results report's "Run with Best Parameters" section,
    # which dumps the template content) can find the file again
    # without it having to be passed through every helper signature.
    config['_model_config_path'] = os.path.abspath(model_config_path)

    # Derive dir2save_input_files (where config files are generated)
    if 'executable_path' in config and 'dir2save_input_files' not in config:
        config['dir2save_input_files'] = os.path.dirname(
            os.path.abspath(config['executable_path'])
        )

    # Overlay the baseline manifest if Gen_Report.generate_simulation_report
    # wrote one.  We attach the parsed dict under _baseline_manifest so
    # callers (get_observation_config / get_spatial_mapping) can prefer
    # manifest values over local re-detection.  We also use it to backfill
    # river_network_mapping_key into the config dict itself when the user
    # left it as the default None in the template — that way the very
    # first time calibration runs after the user runs the config template,
    # the resolved key is available without any extra wiring.
    manifest = _try_load_manifest(config)
    if manifest:
        config[_MANIFEST_CFG_KEY] = manifest
        # Backfill river_network_mapping_key from the manifest when the user
        # left it unset in the template.
        if not config.get('river_network_mapping_key') and manifest.get(
                'river_network_mapping_key'):
            config['river_network_mapping_key'] = manifest[
                'river_network_mapping_key']

    logger.info(f"Loaded {len(config)} config variables"
                + (" + baseline manifest" if manifest else ""))
    logger.debug(f"Config keys: {sorted(config.keys())}")

    return config


def get_gen_input_driver_module(model_config: Dict[str, Any]):
    """
    Import Gen_Input_Driver from the Model_Config support library.

    Parameters
    ----------
    model_config : Dict[str, Any]
        Loaded model configuration (from load_model_config).

    Returns
    -------
    module
        The Gen_Input_Driver module.
    """
    # The Gen_Input_Driver module is in Model_Config/config_support_lib/
    # relative to the model config file.
    # We can locate it from executable_path or from a known relative path.

    # ── Preferred: load the engine from the openWQ CLONE that ships THIS
    # calibration_lib, so a `git pull` of the clone updates it — never a stale
    # project copy.  Layout (mirrors how calibration_lib itself is used):
    #   <clone>/supporting_scripts/3_Calibration/calibration_lib/   (this file)
    #   <clone>/supporting_scripts/1_Model_Config/config_support_lib/
    import importlib, types
    this_dir = os.path.dirname(os.path.abspath(__file__))
    clone_csl = os.path.normpath(
        os.path.join(this_dir, '..', '..', '1_Model_Config', 'config_support_lib'))
    if os.path.isfile(os.path.join(clone_csl, 'Gen_Input_Driver.py')):
        # Evict any config_support_lib modules already imported from a PROJECT
        # copy (e.g. cached when the model config was exec-loaded), so the
        # clone's versions — and their sibling imports — load fresh.
        _marker = os.path.join('1_Model_Config', 'config_support_lib')
        _clone_prefix = os.path.abspath(clone_csl) + os.sep
        for _name in list(sys.modules):
            _m = sys.modules.get(_name)
            if not isinstance(_m, types.ModuleType):
                continue
            _f = getattr(_m, '__file__', None) or ''
            if _marker in _f and not os.path.abspath(_f).startswith(_clone_prefix):
                del sys.modules[_name]
        # Make the clone's config_support_lib win on sys.path.
        while clone_csl in sys.path:
            sys.path.remove(clone_csl)
        sys.path.insert(0, clone_csl)
        return importlib.import_module('Gen_Input_Driver')

    # ── Legacy fallback (no clone layout found) ──────────────────────────────
    # Try multiple strategies to find it
    search_paths = []

    # Strategy 1: It was imported during config loading — check if it's
    # already in sys.modules or in the config namespace
    if 'gJSON_lib' in sys.modules:
        return sys.modules['gJSON_lib']

    # Strategy 2: Look relative to common known locations
    for base_key in ('executable_path', 'dir2save_input_files'):
        if base_key not in model_config:
            continue
        base = os.path.dirname(os.path.abspath(model_config[base_key]))
        # Walk up to find supporting_scripts
        for _ in range(10):
            candidate = os.path.join(
                base, 'supporting_scripts', '1_Model_Config', 'config_support_lib'
            )
            if os.path.isdir(candidate):
                search_paths.append(candidate)
                break
            base = os.path.dirname(base)

    # Strategy 3: Look relative to this file
    this_dir = os.path.dirname(os.path.abspath(__file__))
    # We're in: supporting_scripts/3_Calibration/calibration_lib/
    # Target:   supporting_scripts/1_Model_Config/config_support_lib/
    candidate = os.path.join(
        this_dir, '..', '..', '1_Model_Config', 'config_support_lib'
    )
    candidate = os.path.normpath(candidate)
    if os.path.isdir(candidate):
        search_paths.append(candidate)

    for path in search_paths:
        if path not in sys.path:
            sys.path.insert(0, path)
        try:
            import Gen_Input_Driver as gJSON_lib
            return gJSON_lib
        except ImportError:
            continue

    raise ImportError(
        "Could not find Gen_Input_Driver module. "
        "Searched in: " + ", ".join(search_paths)
    )


def generate_config_for_eval(
    model_config: Dict[str, Any],
    eval_dir: str,
    suppress_report: bool = True
) -> None:
    """
    Generate OpenWQ config files for a calibration evaluation directory.

    Calls Gen_Input_Driver with the model config, but redirects output
    to the evaluation directory.

    Parameters
    ----------
    model_config : Dict[str, Any]
        Loaded model configuration (from load_model_config).
    eval_dir : str
        Path to the evaluation working directory.
    suppress_report : bool
        If True, suppress HTML report generation during config generation.
    """
    gJSON_lib = get_gen_input_driver_module(model_config)

    # Clone config and redirect output directory
    eval_config = copy.deepcopy(model_config)
    eval_config['dir2save_input_files'] = str(eval_dir)

    # Reuse the heavy Copernicus LULC processing: seed this eval's
    # ss_copernicus_files with the baseline's already-computed LULC areas, so the
    # per-eval (dynamic-coefficient) regeneration only re-applies the cheap
    # climate adjustment instead of re-clipping the multi-GB ESA-CCI rasters.
    # The climate params still affect the loading every eval; only the slow
    # raster clipping is skipped.
    try:
        import shutil as _shutil
        _exe = model_config.get('executable_path', '')
        if _exe:
            _base_areas = os.path.join(os.path.dirname(_exe), 'openwq_in',
                                       'ss_copernicus_files', 'lulc_areas_all.csv')
            if os.path.isfile(_base_areas):
                _eval_ss = os.path.join(str(eval_dir), 'openwq_in', 'ss_copernicus_files')
                os.makedirs(_eval_ss, exist_ok=True)
                _dst = os.path.join(_eval_ss, 'lulc_areas_all.csv')
                if not os.path.isfile(_dst):
                    _shutil.copy2(_base_areas, _dst)
    except Exception:
        pass  # fall back to full regeneration if seeding fails

    # Suppress report generation during calibration
    if suppress_report:
        eval_config.pop('generate_report', None)
        eval_config.pop('open_report', None)

    # Remove variables that Gen_Input_Driver doesn't accept
    keys_to_remove = []
    for k in eval_config:
        if k in _EXCLUDED_VARS:
            keys_to_remove.append(k)
    for k in keys_to_remove:
        eval_config.pop(k, None)

    # Remove report-only keys that Gen_Input_Driver does not consume.
    #
    # NOTE: do NOT drop ``file_manager_path`` — Gen_Input_Driver reads it
    # (via **kwargs) to auto-detect the simulation period for the Copernicus
    # source/sink loads when ``ss_method_copernicus_period`` is None.  It is
    # one of the model-config settings that must be carried over, exactly as
    # the model-config script passes it when run directly.  (``executable_path``
    # is kept too — harmless via **kwargs and lets any executable-relative
    # derivation behave as in the user's own run.)
    for k in ('generate_report', 'open_report',
              'river_network_shapefile',
              # Report-side configuration not consumed by Gen_Input_Driver.
              'river_network_mapping_key',
              'observation_compartments',
              'observation_data_source', 'grqa_local_data_path',
              'grqa_buffer_km', 'user_observation_csv',
              'mpi_np', 'ss_method_copernicus_use_proxy_if_outside_range'):
        eval_config.pop(k, None)

    # Run Gen_Input_Driver with the cwd set to the model-config's own
    # directory, so RELATIVE paths inside the model config (e.g.
    # path2selected_NATIVE_BGC_FLEX_framework =
    # "config_support_lib/BGC_templates/…/SWAT_full_nutrients.json") resolve
    # exactly as they do when the user runs the model-config script from
    # that directory.  This mirrors load_model_config(), which exec's the
    # config with cwd set to its directory.  Output still goes to the
    # (absolute) eval_dir via dir2save_input_files, so cwd doesn't affect it.
    _mcp = model_config.get('_model_config_path')
    _cfg_dir = os.path.dirname(os.path.abspath(_mcp)) if _mcp else None
    _old_cwd = os.getcwd()
    try:
        if _cfg_dir and os.path.isdir(_cfg_dir):
            os.chdir(_cfg_dir)

        # Calibration deliberately (re)generates each evaluation's config every
        # time - reusing the template's already-processed GRQA/Copernicus via the
        # baseline manifest, NOT re-running them.  So the interactive
        # "inputs already exist" guard must never fire here.
        eval_config["force_regenerate"] = True
        # Call Gen_Input_Driver
        try:
            gJSON_lib.Gen_Input_Driver(**eval_config)
        except TypeError as e:
            # If there are unexpected kwargs, filter them out and retry
            logger.warning(f"Gen_Input_Driver got unexpected kwargs: {e}")
            # Get the function signature to determine valid kwargs
            import inspect
            sig = inspect.signature(gJSON_lib.Gen_Input_Driver)
            valid_params = set(sig.parameters.keys())
            has_var_keyword = any(
                p.kind == inspect.Parameter.VAR_KEYWORD
                for p in sig.parameters.values()
            )
            if not has_var_keyword:
                filtered = {k: v for k, v in eval_config.items()
                            if k in valid_params}
                gJSON_lib.Gen_Input_Driver(**filtered)
            else:
                raise
    finally:
        os.chdir(_old_cwd)


def get_observation_config(model_config: Dict[str, Any]) -> Dict[str, Any]:
    """
    Extract observation data configuration from the model config.

    Parameters
    ----------
    model_config : Dict[str, Any]
        Loaded model configuration.

    Returns
    -------
    Dict[str, Any]
        Observation data config with keys:
        - source: "grqa", "user_csv", or "skip"
        - grqa_local_data_path: path to GRQA data (if source == "grqa")
        - grqa_buffer_km: search radius (if source == "grqa")
        - user_observation_csv: path to user CSV (if source == "user_csv")
        - basin_shapefile: path to basin shapefile (if available)
        - basin_mapping_key: mapping key in basin shapefile
        - river_network_shapefile: path to river network shapefile
    """
    obs_config = {
        "source": model_config.get("observation_data_source", "skip"),
    }

    # GRQA settings
    obs_config["grqa_local_data_path"] = model_config.get(
        "grqa_local_data_path", None
    )
    obs_config["grqa_buffer_km"] = model_config.get("grqa_buffer_km", 100)

    # User CSV settings
    obs_config["user_observation_csv"] = model_config.get(
        "user_observation_csv", None
    )

    # Basin shapefile
    basin_info = model_config.get("ss_method_copernicus_basin_info", {})
    if isinstance(basin_info, dict):
        obs_config["basin_shapefile"] = basin_info.get("path_to_shp", None)
        obs_config["basin_mapping_key"] = basin_info.get("mapping_key", None)
    else:
        obs_config["basin_shapefile"] = None
        obs_config["basin_mapping_key"] = None

    # River network
    obs_config["river_network_shapefile"] = model_config.get(
        "river_network_shapefile", None
    )
    # User-supplied (or manifest-resolved) override for the column in the
    # river-network shapefile whose values match the reach IDs OpenWQ
    # writes to HDF5.  ``None`` means "let downstream auto-detect".
    obs_config["river_network_mapping_key"] = model_config.get(
        "river_network_mapping_key", None
    )

    # Species
    obs_config["chemical_species"] = model_config.get(
        "chemical_species", []
    )

    return obs_config


def get_observation_period(model_config: Dict[str, Any],
                           work_dir: Optional[str] = None,
                           log=None) -> Optional[Dict[str, Any]]:
    """Return the time span covered by the available observation data.

    Reads the SAME already-prepared data as the obs-prep functions (no
    re-download): the GRQA clipped observations CSV or the user observation
    CSV.  Used by the calibration setup report to *propose* a calibration /
    validation split (first half = calibration, second half = validation).

    Returns
    -------
    Optional[Dict[str, Any]]
        ``{"obs_start": "YYYY-MM-DD HH:MM", "obs_end": ..., "n_obs": int,
        "dates_by_species": {species: [epoch_ms, ...]}}`` or ``None`` when no
        observation data / dates can be found.  ``dates_by_species`` lets the
        setup-report slider show a live "observations in window" count that
        updates with the selected target species.
    """
    import pandas as pd
    _log = log or (lambda *a, **k: None)
    source = model_config.get("observation_data_source", "skip")
    dd = None      # primary-only obs → slider range + in-window count
    dd_all = None  # ALL obs (primary + secondary) → per-species coverage lanes

    # 1) Prefer the already-prepared objective CSV in the calibration work
    #    dir — it carries the ``is_primary`` flag, so the slider's range and
    #    in-window count reflect exactly the observations the metric uses by
    #    default (primary / pouring-point stations).  This makes the
    #    calibratable range START at the first primary obs, structurally
    #    preventing the user from choosing a window the metric can't score.
    #    Only present after a prior run; first-time setup falls back to (2).
    try:
        if work_dir:
            _prep = os.path.join(work_dir, 'calibration_observations.csv')
            if os.path.isfile(_prep):
                pf = pd.read_csv(_prep)
                if 'datetime' in pf.columns and 'species' in pf.columns:
                    _base = pd.DataFrame({
                        'datetime': pd.to_datetime(pf['datetime'],
                                                   errors='coerce'),
                        'species': pf['species'].astype(str),
                        'is_primary': (
                            pf['is_primary'].fillna(True).astype(bool)
                            if 'is_primary' in pf.columns else True),
                    }).dropna(subset=['datetime'])
                    if len(_base):
                        # ALL obs → coverage lanes (so secondary-only species
                        # like NH4 still show up); primary subset → range/count.
                        dd_all = _base[['datetime', 'species']].copy()
                        _prim = _base[_base['is_primary']]
                        dd = _prim[['datetime', 'species']].copy()
                        if len(dd) == 0:
                            dd = None
    except Exception as exc:
        _log(f"Could not read prepared observation CSV: {exc}")
        dd = None

    # 2) Fallback: the raw clipped GRQA / user observation CSV.
    try:
        if dd is None and source == "grqa":
            dir2save = model_config.get('dir2save_input_files')
            if not dir2save:
                exe = model_config.get('executable_path', '')
                dir2save = os.path.dirname(os.path.abspath(exe)) if exe else None
            if dir2save:
                obs_csv = os.path.join(dir2save, 'openwq_in',
                                       'grqa_clipped_data',
                                       'grqa_clipped_observations.csv')
                if os.path.isfile(obs_csv):
                    df = pd.read_csv(obs_csv)
                    col = ('obs_date' if 'obs_date' in df.columns
                           else ('datetime' if 'datetime' in df.columns
                                 else None))
                    spc = ('model_species' if 'model_species' in df.columns
                           else ('species' if 'species' in df.columns else None))
                    if col:
                        dd = pd.DataFrame({
                            'datetime': pd.to_datetime(df[col], errors='coerce'),
                            'species': (df[spc].astype(str) if spc else 'all'),
                        }).dropna(subset=['datetime'])
        elif dd is None and source == "user_csv":
            csv_path = model_config.get('user_observation_csv')
            if csv_path and os.path.isfile(csv_path):
                df = pd.read_csv(csv_path)
                cols = {c.lower(): c for c in df.columns}
                spc_c = cols.get('parameter') or cols.get('species')
                _spc = (df[spc_c].astype(str) if spc_c else 'all')
                if 'datetime' in cols:
                    dt = pd.to_datetime(df[cols['datetime']], errors='coerce')
                    dd = pd.DataFrame({'datetime': dt, 'species': _spc}
                                      ).dropna(subset=['datetime'])
                elif all(k in cols for k in ('year', 'month', 'day')):
                    parts = pd.DataFrame({
                        'year': pd.to_numeric(df[cols['year']], errors='coerce'),
                        'month': pd.to_numeric(df[cols['month']], errors='coerce'),
                        'day': pd.to_numeric(df[cols['day']], errors='coerce'),
                    })
                    dt = pd.to_datetime(parts, errors='coerce')
                    dd = pd.DataFrame({'datetime': dt, 'species': _spc}
                                      ).dropna(subset=['datetime'])
    except Exception as exc:
        _log(f"Could not determine observation period: {exc}")
        return None

    if dd is None or len(dd) == 0:
        return None

    # Raw fallback sources (GRQA / user CSV) carry no is_primary flag, so the
    # "all obs" coverage is the same data as the primary set there.
    if dd_all is None:
        dd_all = dd

    dates = dd['datetime']
    # Epoch-ms timestamps grouped by species. Cap per species to keep the
    # embedded payload small (counts stay representative for the window).
    _MAX_PER_SP = 6000

    def _by_species(frame):
        out: Dict[str, list] = {}
        for sp, grp in frame.groupby('species'):
            ms = (grp['datetime'].astype('int64') // 1_000_000).tolist()
            if len(ms) > _MAX_PER_SP:
                step = len(ms) // _MAX_PER_SP + 1
                ms = ms[::step]
            out[str(sp)] = ms
        return out

    # Primary-only → drives the in-window obs count (what the metric scores).
    dates_by_species = _by_species(dd)
    # All obs (primary + secondary) → per-species coverage lanes, so a species
    # with only secondary obs (e.g. NH4 here) still appears under the slider.
    dates_by_species_all = _by_species(dd_all)

    # TRUE per-species totals + full date range (computed on the UNCAPPED
    # frames). The tick arrays above are subsampled for very large datasets,
    # so the coverage-lane LABELS must read their counts/ranges from here to
    # stay correct for ANY number of observations.
    _prim_n: Dict[str, int] = {
        str(sp): int(len(grp)) for sp, grp in dd.groupby('species')
    }
    species_summary: Dict[str, Dict[str, Any]] = {}
    for sp, grp in dd_all.groupby('species'):
        d = grp['datetime']
        species_summary[str(sp)] = {
            "n": int(len(grp)),
            "n_primary": _prim_n.get(str(sp), 0),
            "start": int(d.min().value // 1_000_000),
            "end": int(d.max().value // 1_000_000),
        }

    return {
        "obs_start": dates.min().strftime('%Y-%m-%d %H:%M'),
        "obs_end": dates.max().strftime('%Y-%m-%d %H:%M'),
        "n_obs": int(len(dates)),
        "dates_by_species": dates_by_species,
        "dates_by_species_all": dates_by_species_all,
        "species_summary": species_summary,
    }


def get_model_sim_period(model_config: Dict[str, Any],
                         log=None) -> Optional[Dict[str, Any]]:
    """Return the model's full simulation period from its control file.

    SUMMA   → ``simStartTime`` / ``simEndTime`` in the fileManager.
    mizuRoute → ``<sim_start>`` / ``<sim_end>`` in the control file.

    Used as the *outer* bound for the calibration/validation slider (the
    span outside the observation period is greyed out).

    Returns
    -------
    Optional[Dict[str, Any]]
        ``{"sim_start": "...", "sim_end": "...", "hostmodel": "..."}`` or
        ``None`` when the control file / dates can't be read.
    """
    import re
    _log = log or (lambda *a, **k: None)
    fm = (model_config.get('file_manager_path')
          or model_config.get('control_file_path') or '')
    if not fm or not os.path.isfile(fm):
        return None
    hostmodel = (model_config.get('hostmodel') or '').lower()
    try:
        with open(fm, 'r') as fh:
            text = fh.read()
    except Exception as exc:
        _log(f"Could not read control file {fm}: {exc}")
        return None

    # SUMMA fileManager: simStartTime '1990-10-01 01:00'
    m1 = re.search(r"simStartTime\s+'([^']*)'", text)
    m2 = re.search(r"simEndTime\s+'([^']*)'", text)
    if m1 and m2:
        return {"sim_start": m1.group(1).strip(),
                "sim_end": m2.group(1).strip(),
                "hostmodel": hostmodel or "summa"}

    # mizuRoute control file: <sim_start>  1990-10-01 00:00:00  ! comment
    m1 = re.search(r"<sim_start>\s*([0-9][^!<\n]*)", text)
    m2 = re.search(r"<sim_end>\s*([0-9][^!<\n]*)", text)
    if m1 and m2:
        return {"sim_start": m1.group(1).strip(),
                "sim_end": m2.group(1).strip(),
                "hostmodel": hostmodel or "mizuroute"}

    return None


def get_model_forcing_period(model_config: Dict[str, Any],
                             log=None) -> Optional[Dict[str, Any]]:
    """Return the ACTUAL forcing-data availability window (host-readable).

    Unlike :func:`get_model_sim_period` (the *declared* period in the control
    file), this reads the real time axis of the forcing/runoff file.  A
    control file can request a period the forcing doesn't cover; calibrating
    outside the forcing yields empty output (no matched obs-sim pairs), so the
    report draws this as its own bar so the mismatch is visible.

      mizuRoute → time span of the runoff file  <input_dir>/<fname_qsim>

    Control-file paths are *container* paths (e.g. ``/code/...``); they are
    translated back to host paths using the known host location of the control
    file itself (longest-common-suffix prefix swap).  Returns
    ``{"forcing_start", "forcing_end", "source"}`` (``'YYYY-MM-DD HH:MM'``
    strings) or ``None`` when it can't be determined (missing file, no
    netCDF4, unreadable) — the report then simply omits the forcing bar.
    """
    import re
    _log = log or (lambda *a, **k: None)
    fm = (model_config.get('file_manager_path')
          or model_config.get('control_file_path') or '')
    if not fm or not os.path.isfile(fm):
        return None
    try:
        with open(fm, 'r') as fh:
            text = fh.read()
    except Exception:
        return None

    def _c2h(cpath, cref, href):
        """Translate a container path to host using a reference dir known on
        both sides (the control-file dir == its <ancil_dir>)."""
        if not cpath:
            return cpath
        cp = cref.rstrip('/').split('/')
        hp = href.rstrip('/').split('/')
        i = 0
        while i < len(cp) and i < len(hp) and cp[-1 - i] == hp[-1 - i]:
            i += 1
        if i == 0:
            return cpath
        c_pre = '/'.join(cp[:len(cp) - i])
        h_pre = '/'.join(hp[:len(hp) - i])
        return (h_pre + cpath[len(c_pre):]) if (c_pre and cpath.startswith(c_pre)) else cpath

    def _nc_span(p):
        """First/last 'time' value of a NetCDF as 'YYYY-MM-DD HH:MM' strings."""
        try:
            from netCDF4 import Dataset, num2date  # type: ignore
        except Exception:
            return None
        if not p or not os.path.isfile(p):
            return None
        try:
            ds = Dataset(p)
            try:
                tv = ds.variables.get('time')
                if tv is None or len(tv) == 0:
                    return None
                u = getattr(tv, 'units', None)
                if not u:
                    return None
                cal = getattr(tv, 'calendar', 'standard')
                d0 = num2date(tv[0], u, cal)
                d1 = num2date(tv[-1], u, cal)
                f = lambda d: ('%04d-%02d-%02d %02d:%02d'
                               % (d.year, d.month, d.day,
                                  getattr(d, 'hour', 0), getattr(d, 'minute', 0)))
                return (f(d0), f(d1))
            finally:
                ds.close()
        except Exception as exc:
            _log("forcing nc read failed (%s): %s" % (p, exc))
            return None

    host_dir = os.path.dirname(os.path.abspath(fm))

    # mizuRoute: runoff file = <input_dir>/<fname_qsim>.  <ancil_dir> (container)
    # corresponds to the control-file dir (host), giving the prefix swap.
    m_anc = re.search(r"<ancil_dir>\s*([^!<\n]+)", text)
    m_in = re.search(r"<input_dir>\s*([^!<\n]+)", text)
    m_q = re.search(r"<fname_qsim>\s*([^!<\n]+)", text)
    if m_in and m_q:
        in_dir = m_in.group(1).strip()
        cref = m_anc.group(1).strip() if m_anc else os.path.dirname(in_dir)
        host_in = _c2h(in_dir, cref, host_dir)
        span = _nc_span(os.path.join(host_in, m_q.group(1).strip()))
        if span:
            return {"forcing_start": span[0], "forcing_end": span[1],
                    "source": m_q.group(1).strip(),
                    "hostmodel": (model_config.get('hostmodel') or 'mizuroute').lower()}

    return None


def validate_ss_reach_mapping(model_config: Dict[str, Any], log=None) -> Optional[str]:
    """Preflight for mizuRoute + Copernicus SS: do the basin shapefile's
    ``mapping_key`` values intersect the river-network reach ``segId``s?

    Returns a warning STRING when they're DISJOINT (the loads would silently
    apply to no reach → all-zero SS output — the lumped-vs-delineated-catchment
    trap), else ``None``.  Best-effort: returns ``None`` on any read failure so
    it never blocks a run.
    """
    import re
    _log = log or (lambda *a, **k: None)
    try:
        hostmodel = (model_config.get('hostmodel') or '').lower()
        ss_method = (model_config.get('ss_method') or '').lower()
        if 'mizuroute' not in hostmodel or 'copernicus' not in ss_method:
            return None
        binfo = model_config.get('ss_method_copernicus_basin_info') or {}
        shp = binfo.get('path_to_shp')
        key = binfo.get('mapping_key')
        fm = (model_config.get('file_manager_path')
              or model_config.get('control_file_path') or '')
        if not (shp and key and os.path.isfile(shp) and fm and os.path.isfile(fm)):
            return None
        text = open(fm, encoding='utf-8', errors='ignore').read()
        _m = re.search(r"<fname_ntopOld>\s*([^\n!<]+)", text)
        topo_name = _m.group(1).strip() if _m else 'topology.nc'
        _ms = re.search(r"<varname_segId>\s*([^\n!<]+)", text)
        seg_var = _ms.group(1).strip() if _ms else 'segId'
        # topology lives in the control file's own dir (= <ancil_dir>).
        topo = os.path.join(os.path.dirname(os.path.abspath(fm)), topo_name)
        if not os.path.isfile(topo):
            return None

        def _norm(v):
            try:
                return str(int(float(v)))
            except Exception:
                return str(v).strip()

        try:
            import geopandas as gpd
            _vals = list(gpd.read_file(shp)[key])
        except Exception:
            import fiona
            with fiona.open(shp) as src:
                _vals = [f['properties'].get(key) for f in src]
        basin_ids = set(_norm(v) for v in _vals if v is not None)

        from netCDF4 import Dataset
        ds = Dataset(topo)
        try:
            seg = ds.variables.get(seg_var)
            if seg is None:
                return None
            reach_ids = set(str(int(v)) for v in seg[:])
        finally:
            ds.close()

        if not basin_ids or not reach_ids or (basin_ids & reach_ids):
            return None
        return (
            "Copernicus SS reach-mapping mismatch: NONE of the basin shapefile's "
            f"'{key}' values match the river-network reach segIds — openWQ would "
            "apply the loads to no reach (ALL-ZERO source/sink output).\n"
            f"    basin shapefile : {shp}\n"
            f"    basin '{key}'    : {sorted(basin_ids)[:5]}{' ...' if len(basin_ids) > 5 else ''}\n"
            f"    reach segIds     : {sorted(reach_ids)[:5]}{' ...' if len(reach_ids) > 5 else ''}\n"
            "  Fix: point ss_method_copernicus_basin_info at the DELINEATED "
            "per-reach catchment whose mapping_key == segId (not the lumped "
            "HRUs_GRUs catchment).")
    except Exception as exc:
        _log(f"validate_ss_reach_mapping skipped: {exc}")
        return None


def get_spatial_mapping(model_config: Dict[str, Any]) -> Dict[str, Any]:
    """Return the resolved spatial-mapping convention for this model run.

    Mirrors the host-specific branches in
    ``Gen_Report.generate_simulation_report`` (around lines 2640-2668) so
    that calibration code uses the same conventions as the config viewer
    without duplicating the detection logic.

    Resolution order for each field:
      1. The baseline manifest written by Gen_Report (most authoritative —
         it captures what the user actually saw on the report).
      2. Explicit override on the model config (e.g.
         ``river_network_mapping_key`` in the template).
      3. Host-defaulted value (``reachID``/``hruId`` for the H5 key,
         ``HRU_ID``/``SegId`` for the shapefile key, etc.).

    Returns a dict with:
      - ``hostmodel`` — lower-cased host model name
      - ``h5_mapping_key`` — column name OpenWQ writes to HDF5 (``reachID``
        for mizuRoute, ``hruId`` for SUMMA)
      - ``river_network_mapping_key`` — column in river-network shapefile
      - ``basin_mapping_key`` — column in basin shapefile
      - ``feature_label`` — label used in viewer trace names
    """
    manifest = model_config.get(_MANIFEST_CFG_KEY) or {}
    hostmodel = (model_config.get('hostmodel') or
                 manifest.get('hostmodel') or '').lower()

    # Sensible host-aware defaults (mirrors Gen_Report.py:2645/2660)
    if hostmodel == 'summa':
        default_h5_key = 'hruId'
        default_shp_key = 'HRU_ID'
    else:
        default_h5_key = 'reachID'
        default_shp_key = 'SegId'

    # h5_mapping_key: manifest > user override > host default
    h5_key = (manifest.get('h5_mapping_key')
              or model_config.get('openwq_h5_mapping_key')
              or default_h5_key)

    # river_network_mapping_key: manifest > model config > default
    rn_key = (manifest.get('river_network_mapping_key')
              or model_config.get('river_network_mapping_key')
              or default_shp_key)

    # basin_mapping_key: manifest > model config (via ss_method_copernicus_basin_info)
    basin_info = model_config.get('ss_method_copernicus_basin_info', {})
    basin_user_key = (basin_info.get('mapping_key')
                      if isinstance(basin_info, dict) else None)
    basin_key = (manifest.get('basin_mapping_key')
                 or basin_user_key
                 or ('HRU_ID' if hostmodel == 'summa' else ''))

    feature_label = (manifest.get('feature_label')
                     or model_config.get('openwq_h5_mapping_key')
                     or rn_key)

    return {
        'hostmodel': hostmodel,
        'h5_mapping_key': h5_key,
        'river_network_mapping_key': rn_key,
        'basin_mapping_key': basin_key,
        'feature_label': feature_label,
        'from_manifest': bool(manifest),
    }


def prepare_calibration_observations_csv(model_config: Dict[str, Any],
                                         output_dir: str,
                                         log=None) -> Optional[str]:
    """Resolve a calibration-ready observations CSV straight from the model
    config the user set up, honouring ``observation_data_source``:

      * ``"grqa"``    → reshape the GRQA data the model-config setup already
                        clipped to the basin (folder *or* "auto" download —
                        both end up in ``openwq_in/grqa_clipped_data/``).
      * ``"user_csv"``→ reshape the user-provided CSV
                        (``user_observation_csv``).
      * ``"skip"``    → no observations.

    All paths come from ``model_config`` (no re-download).  Returns the path
    to ``<output_dir>/calibration_observations.csv`` or ``None`` when nothing
    is configured / preparable.
    """
    src = (model_config.get('observation_data_source') or 'skip').strip().lower()
    if src == 'grqa':
        return prepare_grqa_calibration_csv(model_config, output_dir, log)
    if src == 'user_csv':
        return prepare_user_csv_calibration_csv(model_config, output_dir, log)
    return None


def _match_and_write_obs(model_config: Dict[str, Any],
                         station_locations: Dict[str, tuple],
                         records: List[Dict[str, Any]],
                         output_dir: str,
                         source_prefix: str,
                         log=None) -> Optional[str]:
    """Match station locations → model feature (mizuRoute reach / SUMMA HRU)
    with the SAME shared ``spatial_matching`` helper the results-report map
    uses, then write the objective-format CSV
    (``datetime, reach_id, species, value, units, source, is_primary``).

    ``records`` is a list of dicts with keys
    ``site_id, datetime, species, value, units``.  Returns the CSV path or
    ``None`` on failure.
    """
    import pandas as pd
    _log = log or (lambda *a, **k: None)
    if not station_locations or not records:
        _log("No stations / records to prepare observations from.")
        return None

    spatial = get_spatial_mapping(model_config)
    hostmodel = spatial.get('hostmodel', '')
    obs_cfg = get_observation_config(model_config)
    river_shp = obs_cfg.get('river_network_shapefile')
    basin_shp = obs_cfg.get('basin_shapefile')

    # Lazy imports (avoid a circular import with the results-report module,
    # and keep config_integration importable without geopandas/fiona).
    try:
        from .Gen_Calibration_Results_Report import (
            _import_spatial_matching, _load_shapefile_as_geojson)
        sm = _import_spatial_matching()
    except Exception as e:  # pragma: no cover - defensive
        _log(f"spatial_matching unavailable for obs preparation: {e}")
        return None

    river_gj = _load_shapefile_as_geojson(river_shp)[0] if river_shp else None
    basin_gj = _load_shapefile_as_geojson(basin_shp)[0] if basin_shp else None

    s2f, primary = sm.match_stations(
        station_locations,
        hostmodel=hostmodel,
        river_geojson=river_gj,
        basin_geojson=basin_gj,
        river_mapping_key=spatial.get('river_network_mapping_key') or 'SegId',
        basin_mapping_key=spatial.get('basin_mapping_key') or 'HRU_ID',
        log=_log,
    )
    if not s2f:
        _log("No observation stations matched a model feature.")
        return None

    rows = []
    for rec in records:
        sid = str(rec.get('site_id', ''))
        feat = s2f.get(sid)
        if not feat:
            continue
        rows.append({
            'datetime': rec.get('datetime'),
            'reach_id': feat,
            'species': rec.get('species'),
            'value': rec.get('value'),
            'units': rec.get('units', ''),
            'source': f'{source_prefix}:{sid}',
            'is_primary': sid in primary,
        })
    if not rows:
        _log("No observations fell on a matched feature.")
        return None

    out_df = pd.DataFrame(rows).dropna(
        subset=['datetime', 'value', 'species', 'reach_id'])
    if out_df.empty:
        return None

    os.makedirs(output_dir, exist_ok=True)
    out_path = os.path.join(output_dir, 'calibration_observations.csv')
    out_df.to_csv(out_path, index=False)
    _log(f"Prepared {len(out_df)} calibration observations "
         f"({out_df['reach_id'].nunique()} features, "
         f"{out_df['species'].nunique()} species) → {out_path}")
    return out_path


def prepare_grqa_calibration_csv(model_config: Dict[str, Any],
                                 output_dir: str,
                                 log=None) -> Optional[str]:
    """Reshape the GRQA observations already extracted by the model-config
    setup into the objective-function CSV format — *no re-download*.

    When the user runs their model-config script, ``Gen_Report`` clips the
    GRQA database (folder *or* "auto" Zenodo download) to the basin and
    writes::

        <dir2save_input_files>/openwq_in/grqa_clipped_data/
            grqa_clipped_observations.csv   (station-based: obs_date, site_id,
                                             model_species, obs_value, unit …)
            grqa_clipped_stations.csv       (site_id, lat_wgs84, lon_wgs84 …)

    Returns the prepared CSV path, or ``None`` when the clipped data is
    absent (caller surfaces a clear "run the model-config setup first"
    message).
    """
    import pandas as pd
    _log = log or (lambda *a, **k: None)

    dir2save = model_config.get('dir2save_input_files')
    if not dir2save:
        exe = model_config.get('executable_path', '')
        dir2save = os.path.dirname(os.path.abspath(exe)) if exe else None
    if not dir2save:
        return None

    clipped_dir = os.path.join(dir2save, 'openwq_in', 'grqa_clipped_data')
    obs_csv = os.path.join(clipped_dir, 'grqa_clipped_observations.csv')
    stn_csv = os.path.join(clipped_dir, 'grqa_clipped_stations.csv')
    if not os.path.isfile(obs_csv):
        _log(f"No clipped GRQA observations at {obs_csv}")
        return None

    obs = pd.read_csv(obs_csv)
    if obs.empty:
        _log("Clipped GRQA observations CSV is empty.")
        return None

    # Column aliases (clipped/GRQA names → objective names)
    date_col = 'obs_date' if 'obs_date' in obs.columns else 'datetime'
    spc_col = 'model_species' if 'model_species' in obs.columns else 'species'
    val_col = 'obs_value' if 'obs_value' in obs.columns else 'value'
    unit_col = next((c for c in ('unit', 'units') if c in obs.columns), None)
    lat_col = 'lat_wgs84' if 'lat_wgs84' in obs.columns else 'lat'
    lon_col = 'lon_wgs84' if 'lon_wgs84' in obs.columns else 'lon'

    # Station locations {site_id: (lat, lon)} — from the stations CSV, or
    # derived from the observations as a fallback.
    station_locations: Dict[str, tuple] = {}
    if os.path.isfile(stn_csv):
        st = pd.read_csv(stn_csv)
        s_lat = 'lat_wgs84' if 'lat_wgs84' in st.columns else 'lat'
        s_lon = 'lon_wgs84' if 'lon_wgs84' in st.columns else 'lon'
        for _, r in st.iterrows():
            sid = str(r.get('site_id', ''))
            if sid:
                try:
                    station_locations[sid] = (float(r[s_lat]), float(r[s_lon]))
                except (TypeError, ValueError, KeyError):
                    pass
    if not station_locations and lat_col in obs.columns and lon_col in obs.columns:
        for sid, grp in obs.groupby('site_id'):
            r = grp.iloc[0]
            try:
                station_locations[str(sid)] = (float(r[lat_col]), float(r[lon_col]))
            except (TypeError, ValueError):
                pass

    records = [{
        'site_id': str(r.get('site_id', '')),
        'datetime': r.get(date_col),
        'species': r.get(spc_col),
        'value': r.get(val_col),
        'units': (r.get(unit_col) if unit_col else ''),
    } for _, r in obs.iterrows()]

    return _match_and_write_obs(
        model_config, station_locations, records, output_dir, 'GRQA', _log)


def prepare_user_csv_calibration_csv(model_config: Dict[str, Any],
                                     output_dir: str,
                                     log=None) -> Optional[str]:
    """Reshape the user-provided observation CSV into the objective format.

    Expected user columns (case-insensitive, any order, as documented in the
    model-config template)::

        station_id, lat, lon, parameter, year, month, day, minute, value, units

    Each station is matched to its model feature with the shared
    ``spatial_matching`` helper.  If the CSV is ALREADY in objective format
    (has both ``datetime`` and ``reach_id`` columns) it is returned as-is.
    Returns the prepared CSV path or ``None``.
    """
    import pandas as pd
    _log = log or (lambda *a, **k: None)

    csv_path = model_config.get('user_observation_csv')
    if not csv_path or not os.path.isfile(csv_path):
        _log(f"User observation CSV not found: {csv_path}")
        return None

    df = pd.read_csv(csv_path)
    if df.empty:
        _log("User observation CSV is empty.")
        return None

    cols = {c.lower(): c for c in df.columns}

    # Already in objective format → use as-is (no reshape needed).
    if 'datetime' in cols and 'reach_id' in cols:
        _log("User CSV already in objective format; using as-is.")
        return csv_path

    required = ['station_id', 'lat', 'lon', 'parameter', 'year', 'month',
                'day', 'value']
    missing = [n for n in required if n not in cols]
    if missing:
        _log(f"User CSV missing required columns {missing}; cannot reshape.")
        return None

    def C(name):
        return cols.get(name)

    station_locations: Dict[str, tuple] = {}
    for sid, grp in df.groupby(C('station_id')):
        r = grp.iloc[0]
        try:
            station_locations[str(sid)] = (float(r[C('lat')]), float(r[C('lon')]))
        except (TypeError, ValueError):
            pass

    minute_c = cols.get('minute')
    unit_c = cols.get('units') or cols.get('unit')
    records = []
    for _, r in df.iterrows():
        try:
            y = int(r[C('year')]); mo = int(r[C('month')]); d = int(r[C('day')])
            mins = int(r[minute_c]) if (minute_c and pd.notna(r[minute_c])) else 0
            dt = pd.Timestamp(year=y, month=mo, day=d) + pd.Timedelta(minutes=mins)
        except (TypeError, ValueError):
            continue
        records.append({
            'site_id': str(r[C('station_id')]),
            'datetime': dt.isoformat(sep=' '),
            'species': r[C('parameter')],
            'value': r[C('value')],
            'units': (r[unit_c] if unit_c else ''),
        })

    return _match_and_write_obs(
        model_config, station_locations, records, output_dir, 'USER', _log)


def get_container_config(model_config: Dict[str, Any]) -> Dict[str, Any]:
    """
    Extract container runtime configuration from the model config.

    Values not explicitly set in model_config fall back to sensible defaults
    (Docker runtime, ``docker_openwq`` container name, docker-compose.yml
    auto-derived from the ``containers/`` directory relative to the
    executable path).

    Parameters
    ----------
    model_config : Dict[str, Any]
        Loaded model configuration.

    Returns
    -------
    Dict[str, Any]
        Container config with keys matching calibration_driver expectations.
    """
    # Resolve docker_compose_path.  Needed so ModelRunner parses the SAME
    # host→/code volume mount the container actually uses (otherwise it falls
    # back to a wrong ~→/code guess and the per-eval `cd` into the container
    # fails).  Try, in order:
    #   1. an explicit model-config value;
    #   2. .../openwq/bin/exe → .../openwq/containers/docker-compose.yml
    #      (only works when the executable lives inside the openWQ repo);
    #   3. the openWQ repo's containers/ dir relative to THIS file
    #      (config_integration lives at
    #       <openwq>/supporting_scripts/3_Calibration/calibration_lib/), which
    #      is the robust path when the executable is a host-model binary in a
    #      domain folder (e.g. SUMMA/mizuRoute).
    docker_compose = model_config.get("docker_compose_path", "")
    if not docker_compose:
        exe = model_config.get("executable_path", "")
        if exe:
            _openwq_root = os.path.dirname(os.path.dirname(os.path.abspath(exe)))
            _candidate = os.path.join(_openwq_root, "containers", "docker-compose.yml")
            if os.path.isfile(_candidate):
                docker_compose = _candidate
    if not docker_compose:
        _here = os.path.dirname(os.path.abspath(__file__))
        # calibration_lib → 3_Calibration → supporting_scripts → <openwq>
        _openwq = os.path.dirname(os.path.dirname(os.path.dirname(_here)))
        _candidate = os.path.join(_openwq, "containers", "docker-compose.yml")
        if os.path.isfile(_candidate):
            docker_compose = _candidate

    return {
        "executable_path": model_config.get("executable_path", ""),
        "file_manager_path": model_config.get("file_manager_path", ""),
        "control_file_path": model_config.get("control_file_path", ""),
        "mpi_np": model_config.get("mpi_np", 2),
        "hostmodel": model_config.get("hostmodel", "mizuroute"),
        "container_runtime": model_config.get("container_runtime", "docker"),
        "docker_container_name": model_config.get(
            "docker_container_name", "docker_openwq"
        ),
        "docker_compose_path": docker_compose,
    }


def get_bgc_template_path(model_config: Dict[str, Any]) -> Optional[str]:
    """
    Get the absolute path to the BGC template file.

    Parameters
    ----------
    model_config : Dict[str, Any]
        Loaded model configuration.

    Returns
    -------
    Optional[str]
        Absolute path to the BGC template, or None if not applicable.
    """
    bgc_module = model_config.get("bgc_module_name", "")
    if bgc_module != "NATIVE_BGC_FLEX":
        return None

    template_path = model_config.get(
        "path2selected_NATIVE_BGC_FLEX_framework", ""
    )
    if not template_path:
        return None

    # The path may be relative to the model_config_template.py directory
    if not os.path.isabs(template_path):
        # Try relative to executable dir
        exe_dir = os.path.dirname(
            model_config.get("executable_path", "")
        )
        if exe_dir:
            candidate = os.path.join(exe_dir, template_path)
            if os.path.isfile(candidate):
                return os.path.abspath(candidate)

        # Try relative to Model_Config dir
        this_dir = os.path.dirname(os.path.abspath(__file__))
        model_config_dir = os.path.normpath(
            os.path.join(this_dir, '..', '..', '1_Model_Config')
        )
        candidate = os.path.join(model_config_dir, template_path)
        if os.path.isfile(candidate):
            return os.path.abspath(candidate)

    if os.path.isfile(template_path):
        return os.path.abspath(template_path)

    logger.warning(f"BGC template not found: {template_path}")
    return None


def get_module_selections(model_config: Dict[str, Any]) -> Dict[str, Any]:
    """
    Summarize which modules are active and their key settings.

    Parameters
    ----------
    model_config : Dict[str, Any]
        Loaded model configuration.

    Returns
    -------
    Dict[str, Any]
        Summary of active modules and their settings.
    """
    bgc_module = model_config.get("bgc_module_name", "NONE")
    td_module = model_config.get("td_module_name", "NONE")
    ts_module = model_config.get("ts_module_name", "NONE")
    le_module = model_config.get("le_module_name", "NONE")
    si_module = model_config.get("si_module_name", "NONE")
    ss_method = model_config.get("ss_method", "none")

    # BGC template basename
    bgc_template = ""
    bgc_path = model_config.get("path2selected_NATIVE_BGC_FLEX_framework", "")
    if bgc_path:
        bgc_template = os.path.basename(bgc_path)

    # Lateral exchange interface count
    le_config = model_config.get("le_module_config", [])
    le_count = len(le_config) if isinstance(le_config, list) else 0

    # Sorption species
    si_species_params = model_config.get("si_species_params", None)
    si_species = []
    if isinstance(si_species_params, dict):
        si_species = list(si_species_params.keys())

    return {
        "bgc_module": bgc_module,
        "bgc_active": bgc_module not in ("NONE", "", None),
        "bgc_template": bgc_template,
        "td_module": td_module,
        "td_active": td_module == "OPENWQ_NATIVE_TD_ADVDISP",
        "ts_module": ts_module,
        "ts_active": ts_module in ("HYPE_MMF", "HYPE_HBVSED"),
        "le_module": le_module,
        "le_active": le_module == "NATIVE_LE_BOUNDMIX",
        "le_interface_count": le_count,
        "si_module": si_module,
        "si_active": si_module in ("FREUNDLICH", "LANGMUIR"),
        "si_species": si_species,
        "ss_method": ss_method,
        "ss_active": ss_method not in ("none", "", None),
    }


def get_ss_species_with_loads(model_config: Dict[str, Any]) -> set:
    """
    Read the generated SS JSON config to find which model species
    actually have source/sink loads.

    The config template already resolved stoichiometric conversions
    (e.g. TN → NO3-N, NH4 → NH4-N) when generating the SS JSON.
    This function simply reads those resolved species names.

    Returns
    -------
    set
        Set of model species names that have SS load entries
        (e.g. {"NH4-N", "NO3-N"}).
    """
    import json as _json

    output_dir = model_config.get("dir2save_input_files", "")
    ss_method = model_config.get("ss_method", "none")
    if not output_dir or ss_method in ("none", "", None):
        return set()

    # Construct SS JSON path (same naming convention as Gen_SS_Driver)
    ss_json_path = os.path.join(
        output_dir, "openwq_in", f"openWQ_SS_{ss_method}.json"
    )
    if not os.path.isfile(ss_json_path):
        logger.debug(f"SS JSON not found: {ss_json_path}")
        return set()

    try:
        with open(ss_json_path, 'r') as f:
            content = f.read()
        # Strip comment lines (// ...)
        json_start = content.find('{')
        if json_start < 0:
            return set()
        data = _json.loads(content[json_start:])
    except Exception as e:
        logger.warning(f"Could not parse SS JSON: {e}")
        return set()

    species = set()
    for key, entry in data.items():
        if isinstance(entry, dict):
            chem = entry.get("CHEMICAL_NAME", "")
            if chem:
                species.add(chem)
    logger.info(f"SS species with loads: {species}")
    return species


def get_species_observation_availability(
    model_config: Dict[str, Any],
) -> Dict[str, Dict[str, Any]]:
    """
    Check which model species have observation data available.

    For GRQA source, scans the local data directory for CSV files matching
    each species.  For user CSV source, reads the file to identify unique
    species.  Returns a dict keyed by model species name with availability
    metadata.

    Parameters
    ----------
    model_config : Dict[str, Any]
        Loaded model configuration (must include ``chemical_species``,
        ``observation_data_source``, and related fields).

    Returns
    -------
    Dict[str, Dict[str, Any]]
        Mapping of model species → ``{"has_obs": bool, "source": str,
        "details": str}``.  Species not in the model are omitted.
    """
    species_list = model_config.get("chemical_species", [])

    # Resolve "all" → actual species list from BGC template
    _is_all = (
        (isinstance(species_list, str) and species_list.lower() == "all")
        or (isinstance(species_list, list) and len(species_list) == 1
            and isinstance(species_list[0], str)
            and species_list[0].lower() == "all")
    )
    if _is_all:
        try:
            # Use the same resolver that model_config_template uses,
            # but pass the fully-resolved BGC template path
            this_dir = os.path.dirname(os.path.abspath(__file__))
            model_config_lib = os.path.normpath(
                os.path.join(this_dir, '..', '..', '1_Model_Config',
                             'config_support_lib')
            )
            if model_config_lib not in sys.path:
                sys.path.insert(0, model_config_lib)
            from Gen_Input_Driver import extract_species_from_bgc

            resolved_bgc_path = get_bgc_template_path(model_config) or ""
            species_list = extract_species_from_bgc(
                bgc_module_name=model_config.get("bgc_module_name", ""),
                path2selected_NATIVE_BGC_FLEX_framework=resolved_bgc_path,
                phreeqc_mobile_species=model_config.get(
                    "phreeqc_mobile_species", None
                ),
            )
        except Exception as exc:
            logger.debug(f"Could not resolve 'all' species: {exc}")
            species_list = []

    if isinstance(species_list, str):
        species_list = [species_list]
    species_list = list(species_list)

    source = model_config.get("observation_data_source", "skip")
    result: Dict[str, Dict[str, Any]] = {}

    if source == "skip":
        for sp in species_list:
            result[sp] = {
                "has_obs": False,
                "source": "skip",
                "details": "Observations disabled",
            }
        return result

    if source == "grqa":
        grqa_path = model_config.get("grqa_local_data_path", "")

        # Import GRQA mapping tables
        try:
            from .observation_data.grqa_extract_stations import (
                GRQA_V14_FILENAME_MAP,
                GRQA_PARAMETERS,
            )
        except ImportError:
            GRQA_V14_FILENAME_MAP = {}
            GRQA_PARAMETERS = {}

        # Build reverse map: model species → GRQA parameter code
        # The species_mapping in model_config (if any) maps
        # GRQA param → model species.
        grqa_to_model = {}
        species_mapping = model_config.get("grqa_species_mapping", {})
        if isinstance(species_mapping, dict) and species_mapping:
            grqa_to_model = species_mapping
        else:
            # Use identity mapping (species name == GRQA param)
            grqa_to_model = {sp: sp for sp in species_list}

        model_to_grqa = {v: k for k, v in grqa_to_model.items()}

        def _resolve_grqa_param(model_species: str) -> str:
            """Try to map a model species name to a GRQA parameter code.

            Handles common naming conventions like NH4-N → NH4,
            PO4-P → PO4, etc.
            """
            # Direct match first
            if model_species in GRQA_PARAMETERS:
                return model_species
            # Strip common suffixes: -N, -P, _N, _P
            for suffix in ('-N', '-P', '_N', '_P'):
                if model_species.endswith(suffix):
                    base = model_species[:-len(suffix)]
                    if base in GRQA_PARAMETERS:
                        return base
            # Case-insensitive search
            lower = model_species.lower()
            for code in GRQA_PARAMETERS:
                if code.lower() == lower:
                    return code
            return model_species

        for sp in species_list:
            grqa_param = model_to_grqa.get(sp)
            # Fall back to fuzzy resolution if no explicit mapping
            # or if the mapped name is not a known GRQA parameter
            if not grqa_param or grqa_param not in GRQA_PARAMETERS:
                grqa_param = _resolve_grqa_param(sp)
            has_obs = False
            details = ""

            # Check if GRQA parameter is known
            if grqa_param in GRQA_PARAMETERS:
                info = GRQA_PARAMETERS[grqa_param]
                details = info.get("name", grqa_param)

                # Check if file exists locally
                if grqa_path and os.path.isdir(grqa_path):
                    candidates = [grqa_param]
                    if grqa_param in GRQA_V14_FILENAME_MAP:
                        candidates.append(GRQA_V14_FILENAME_MAP[grqa_param])
                    for c in candidates:
                        fpath = os.path.join(grqa_path, f"{c}_GRQA.csv")
                        if os.path.isfile(fpath):
                            has_obs = True
                            break
                    if not has_obs:
                        details += " (file not found locally)"
                else:
                    # GRQA path not set but parameter is valid
                    has_obs = True
                    details += " (available in GRQA)"
            else:
                details = "Not available in GRQA"

            result[sp] = {
                "has_obs": has_obs,
                "source": "grqa",
                "details": details,
            }
        return result

    if source in ("user_csv", "csv"):
        csv_path = model_config.get("user_observation_csv", "")
        csv_species: set = set()
        if csv_path and os.path.isfile(csv_path):
            try:
                import csv as csv_mod
                with open(csv_path, "r", encoding="utf-8") as f:
                    reader = csv_mod.DictReader(f)
                    # Find the species/parameter column
                    fieldnames = [n.lower().strip() for n in (reader.fieldnames or [])]
                    sp_col = None
                    for candidate_col in ("parameter", "species", "param"):
                        if candidate_col in fieldnames:
                            # Get original-case column name
                            idx = fieldnames.index(candidate_col)
                            sp_col = (reader.fieldnames or [])[idx]
                            break
                    if sp_col:
                        for row in reader:
                            val = row.get(sp_col, "").strip()
                            if val:
                                csv_species.add(val)
            except Exception as exc:
                logger.warning(f"Could not read user CSV: {exc}")

        for sp in species_list:
            has_obs = sp in csv_species
            details = "Found in observation CSV" if has_obs else "Not in CSV"
            result[sp] = {
                "has_obs": has_obs,
                "source": "user_csv",
                "details": details,
            }
        return result

    # Unknown source
    for sp in species_list:
        result[sp] = {
            "has_obs": False,
            "source": source,
            "details": f"Unknown source: {source}",
        }
    return result
