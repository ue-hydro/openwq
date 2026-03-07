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
"""

import os
import sys
import copy
import logging
from pathlib import Path
from typing import Dict, Any, Optional, List

logger = logging.getLogger(__name__)

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
}


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

    # Derive dir2save_input_files (where config files are generated)
    if 'executable_path' in config and 'dir2save_input_files' not in config:
        config['dir2save_input_files'] = os.path.dirname(
            os.path.abspath(config['executable_path'])
        )

    logger.info(f"Loaded {len(config)} config variables")
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

    # Remove non-Gen_Input_Driver keys that may be in the config
    # (report-related, executable path, etc.)
    for k in ('executable_path', 'file_manager_path',
              'generate_report', 'open_report',
              'river_network_shapefile',
              'observation_data_source', 'grqa_local_data_path',
              'grqa_buffer_km', 'user_observation_csv',
              'mpi_np', 'ss_method_copernicus_use_proxy_if_outside_range'):
        eval_config.pop(k, None)

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

    # Species
    obs_config["chemical_species"] = model_config.get(
        "chemical_species", []
    )

    return obs_config


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
    # Try to derive docker_compose_path from the executable location.
    # Layout: .../openwq/bin/exe → .../openwq/containers/docker-compose.yml
    docker_compose = model_config.get("docker_compose_path", "")
    if not docker_compose:
        exe = model_config.get("executable_path", "")
        if exe:
            _openwq_root = os.path.dirname(os.path.dirname(os.path.abspath(exe)))
            _candidate = os.path.join(_openwq_root, "containers", "docker-compose.yml")
            if os.path.isfile(_candidate):
                docker_compose = _candidate

    return {
        "executable_path": model_config.get("executable_path", ""),
        "file_manager_path": model_config.get("file_manager_path", ""),
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
