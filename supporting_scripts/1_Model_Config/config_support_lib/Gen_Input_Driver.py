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
# !/usr/bin/env python3

from typing import Dict, List, Union, Any, Optional
import datetime
import json
import os
import re

import Gen_Master_file as mJSON_lib
import Gen_Config_file as cJSON_lib
import Gen_TDmodule_file as tdmJSON_lib
import Gen_LEmodule_file as lemJSON_lib
import Load_BGQmodule_file as bgqmJSON_lib
import Gen_PHREEQCmodule_file as phreeqcJSON_lib
import Gen_TSmodule_file as tsmJSON_lib
import Gen_SImodule_file as simJSON_lib
import Gen_SS_Driver as ssJSON_lib
import Gen_EWF_Driver as ewfJSON_lib

# Upper bound on the number of source/sink JSON entries the climate-adjusted
# (time_series) path may generate before it is considered too large for openWQ
# to load (nlohmann::json balloons in RAM; a >1M-entry file OOM-kills the run).
# The driver estimates the count BEFORE generating and refuses/prompts above it.
_SS_ENTRY_GUARD_MAX = 200_000

# SS entries written per year per spatial-unit per species, by resolution.
# ("native" is detected from the climate file's own time step at run time.)
_SS_PERIODS_PER_YEAR = {"monthly": 12, "daily": 365, "hourly": 8760}


def _count_basin_spatial_units(basin_info):
    """Number of features (HRUs for SUMMA / reaches for mizuRoute) in the basin
    shapefile — the spatial dimension of the SS file.  Returns 1 on any failure
    (keeps the estimate conservative-but-finite rather than crashing the guard).
    """
    try:
        shp = (basin_info or {}).get("path_to_shp") or (basin_info or {}).get("path")
        if not shp or not os.path.isfile(shp):
            return 1
        try:
            import fiona
            with fiona.open(shp) as src:
                return max(int(len(src)), 1)
        except Exception:
            import geopandas as gpd
            return max(int(len(gpd.read_file(shp))), 1)
    except Exception:
        return 1


def _detect_native_periods_per_year(climate_source):
    """Estimate SS periods/year for resolution='native' from the climate file's
    own median time step (e.g. hourly → 8760, 3-hourly → 2920, daily → 365)."""
    try:
        from netCDF4 import Dataset, num2date
        spec = climate_source["precip"]
        ds = Dataset(spec["path"])
        try:
            tv = ds.variables[spec.get("time_key_or_column", "time")]
            if tv is None or len(tv) < 2:
                return 8760
            import numpy as np
            d0 = num2date(tv[0], tv.units, getattr(tv, "calendar", "standard"))
            d1 = num2date(tv[1], tv.units, getattr(tv, "calendar", "standard"))
            dt_h = abs((d1 - d0).total_seconds()) / 3600.0
            return int(round(8760.0 / dt_h)) if dt_h > 0 else 8760
        finally:
            ds.close()
    except Exception:
        return 8760


# [calibration fix] Cache resolved SS resolutions so the size-guard question is
# asked at most ONCE per process. During calibration the SS is regenerated in
# THIS SAME process for every eval (parameter_handler -> Gen_Input_Driver /
# Gen_SS_Driver), so without this the prompt would reappear on every eval and
# block the run. Keyed by the inputs that determine the entry-count estimate.
_SS_RESOLUTION_GUARD_CACHE = {}


def _ss_resolution_size_guard(resolution, n_years, n_units, n_species,
                              climate_source,
                              scope_label="full simulation period"):
    """Estimate the SS entry count for *resolution* and enforce
    ``_SS_ENTRY_GUARD_MAX``.  Returns the (possibly user-downscaled) resolution.

    *n_years* is the span the estimate is taken over and *scope_label*
    describes it in the messages.  In the calibration step the caller passes
    the (shorter) calibration-window span instead of the full simulation
    period — the per-eval SS is trimmed to that window before openWQ reads it,
    so a short window can make 'daily'/'hourly' perfectly affordable.

    Interactive terminal → explains the size and offers coarser options.
    Non-interactive (HPC/batch) → raises with the exact line to set, instead of
    silently generating a multi-GB file or hanging on input().
    """
    import sys

    def _ppy(res):
        if res == "native":
            return _detect_native_periods_per_year(climate_source)
        return _SS_PERIODS_PER_YEAR.get(res, 365)

    def _est(res):
        return int(max(n_years, 1) * _ppy(res)
                   * max(n_units, 1) * max(n_species, 1))

    res = str(resolution).lower()

    # Asked-once cache (see _SS_RESOLUTION_GUARD_CACHE note above): reuse a prior
    # resolution for identical sizing instead of re-prompting on every eval.
    _cache_key = (res, int(max(n_years, 1)), int(max(n_units, 1)),
                  int(max(n_species, 1)), str(scope_label))
    if _cache_key in _SS_RESOLUTION_GUARD_CACHE:
        return _SS_RESOLUTION_GUARD_CACHE[_cache_key]

    est = _est(res)
    if est <= _SS_ENTRY_GUARD_MAX:
        _SS_RESOLUTION_GUARD_CACHE[_cache_key] = res
        return res

    # Coarser options that would fit, in increasing fineness.
    order = ["monthly", "daily", "hourly", "native"]
    fits = [r for r in order if _est(r) <= _SS_ENTRY_GUARD_MAX]
    approx_mb = est * 124 / 1e6   # ~124 bytes/entry observed
    msg = (
        f"\n  ⚠  Source/sink file too large at '{res}' resolution:\n"
        f"     ~{est:,} entries (~{approx_mb:.0f} MB) over the {scope_label}\n"
        f"     = {n_years} yr x {_ppy(res):,}/yr x {n_units} units x {n_species} species.\n"
        f"     openWQ will likely run out of memory loading it.\n"
        f"     Guard limit: {_SS_ENTRY_GUARD_MAX:,} entries "
        f"(set Gen_Input_Driver._SS_ENTRY_GUARD_MAX to change).")
    if fits:
        opts = ", ".join(f"'{r}' (~{_est(r):,})" for r in fits)
        msg += f"\n     Resolutions that fit: {opts}."
    print(msg)

    if not sys.stdin or not sys.stdin.isatty():
        raise ValueError(
            "ss_climate_data_source_adjusting_resolution='" + res + "' produces ~"
            + f"{est:,}" + " SS entries over the " + scope_label
            + " (> " + f"{_SS_ENTRY_GUARD_MAX:,}" + "). "
            "Coarsen it — e.g. set "
            "ss_climate_data_source_adjusting_resolution = "
            + ("'" + fits[0] + "'" if fits else "'monthly'")
            + " — then re-run. (Non-interactive run: not prompting.)")

    # Interactive: let the user pick a coarser option.
    choices = fits + [res]
    print("\n  Choose a resolution:")
    for i, r in enumerate(choices, 1):
        tag = "  (keep — may OOM)" if r == res else ""
        print(f"     {i}) {r} (~{_est(r):,} entries){tag}")
    try:
        sel = input("  Selection [1]: ").strip() or "1"
        idx = int(sel) - 1
        chosen = choices[idx] if 0 <= idx < len(choices) else choices[0]
    except Exception:
        chosen = choices[0]
    print(f"  → using '{chosen}' resolution.")
    _SS_RESOLUTION_GUARD_CACHE[_cache_key] = chosen
    return chosen


def uniform_param(value):
    """Build a spatially-varying parameter table that applies *value* to all cells.

    Returns the dict format expected by the TS module parameters:
        {"0": ["IX", "IY", "IZ", "VALUE"], "1": ["ALL", 1, 1, value]}
    """
    return {"0": ["IX", "IY", "IZ", "VALUE"], "1": ["ALL", 1, 1, value]}


# Copernicus ESA CCI Land Cover availability range
COPERNICUS_LC_YEAR_MIN = 1992
COPERNICUS_LC_YEAR_MAX = 2022


def parse_sim_period_from_control_file(
        file_manager_path: str,
        hostmodel: str,
) -> List[int]:
    """Extract simulation start/end years from a host-model control file.

    Supports two formats:
      - **mizuRoute**: ``<sim_start> YYYY-MM-DD HH:MM:SS``
      - **SUMMA**:     ``simStartTime  'YYYY-MM-DD HH:MM'``

    Parameters
    ----------
    file_manager_path : str
        Path to the host-model control file on the HOST machine.
    hostmodel : str
        ``"mizuroute"`` or ``"summa"``.

    Returns
    -------
    list of int
        ``[start_year, end_year]`` clamped to the Copernicus CCI LC range
        (1992–2022).

    Raises
    ------
    FileNotFoundError
        If *file_manager_path* does not exist.
    ValueError
        If the expected tags cannot be found in the file.
    """
    if not os.path.isfile(file_manager_path):
        raise FileNotFoundError(
            f"Cannot auto-detect simulation period: "
            f"control file not found: {file_manager_path}"
        )

    with open(file_manager_path, 'r') as f:
        content = f.read()

    hm = hostmodel.strip().lower()
    start_year = end_year = None

    if hm == "mizuroute":
        # Format:  <sim_start>  YYYY-MM-DD HH:MM:SS  ! comment
        m_start = re.search(r'<sim_start>\s+(\d{4})-', content)
        m_end = re.search(r'<sim_end>\s+(\d{4})-', content)
        if m_start:
            start_year = int(m_start.group(1))
        if m_end:
            end_year = int(m_end.group(1))

    elif hm == "summa":
        # Format:  simStartTime  'YYYY-MM-DD HH:MM'  ! comment
        m_start = re.search(r"simStartTime\s+'(\d{4})-", content)
        m_end = re.search(r"simEndTime\s+'(\d{4})-", content)
        if m_start:
            start_year = int(m_start.group(1))
        if m_end:
            end_year = int(m_end.group(1))

    else:
        raise ValueError(
            f"Cannot auto-detect simulation period: "
            f"unknown hostmodel '{hostmodel}'. Expected 'mizuroute' or 'summa'."
        )

    if start_year is None or end_year is None:
        raise ValueError(
            f"Cannot auto-detect simulation period from '{file_manager_path}'.\n"
            f"  Expected tags: "
            + ("<sim_start>/<sim_end>" if hm == "mizuroute"
               else "simStartTime/simEndTime")
            + " with YYYY-MM-DD format.\n"
            f"  Set ss_method_copernicus_period = [start_year, end_year] manually."
        )

    return [start_year, end_year]


def extract_species_from_bgc(
        bgc_module_name: str,
        path2selected_NATIVE_BGC_FLEX_framework: str = "",
        phreeqc_mobile_species: Optional[List[str]] = None,
) -> List[str]:
    """Extract all chemical species names from the BGC module configuration.

    Parameters
    ----------
    bgc_module_name : str
        ``"NATIVE_BGC_FLEX"`` or ``"PHREEQC"``.
    path2selected_NATIVE_BGC_FLEX_framework : str
        Path to the NATIVE_BGC_FLEX JSON file (only used when
        bgc_module_name = ``"NATIVE_BGC_FLEX"``).
    phreeqc_mobile_species : list of str, optional
        User-provided list of PHREEQC mobile species (only used when
        bgc_module_name = ``"PHREEQC"``).

    Returns
    -------
    list of str
        Species names extracted from the BGC configuration.
    """
    if bgc_module_name == "NATIVE_BGC_FLEX":
        if not path2selected_NATIVE_BGC_FLEX_framework:
            raise ValueError(
                "chemical_species = 'all' with NATIVE_BGC_FLEX requires "
                "path2selected_NATIVE_BGC_FLEX_framework to be set."
            )
        if not os.path.isfile(path2selected_NATIVE_BGC_FLEX_framework):
            raise FileNotFoundError(
                f"BGC framework file not found: "
                f"{path2selected_NATIVE_BGC_FLEX_framework}"
            )

        with open(path2selected_NATIVE_BGC_FLEX_framework, 'r') as f:
            bgc_data = json.load(f)

        # Runtime format: CHEMICAL_SPECIES.LIST at top level
        # Legacy format: BIOGEOCHEMISTRY_CONFIGURATION.CHEMICAL_SPECIES
        if "CHEMICAL_SPECIES" in bgc_data:
            chem_section = bgc_data["CHEMICAL_SPECIES"]
            if "LIST" in chem_section:
                # Runtime format: LIST is a numbered dict {"1": "NO3-N", ...}
                species = list(chem_section["LIST"].values())
                return species
            else:
                # Old descriptive format at top level
                species = [k for k in chem_section.keys()
                           if not k.startswith("_")]
                return species
        elif "BIOGEOCHEMISTRY_CONFIGURATION" in bgc_data:
            chem_species_dict = (bgc_data
                                 .get("BIOGEOCHEMISTRY_CONFIGURATION", {})
                                 .get("CHEMICAL_SPECIES", {}))
            species = [k for k in chem_species_dict.keys()
                       if not k.startswith("_")]
            return species
        else:
            return []

    elif bgc_module_name == "PHREEQC":
        if phreeqc_mobile_species and len(phreeqc_mobile_species) > 0:
            return list(phreeqc_mobile_species)
        else:
            raise ValueError(
                "chemical_species = 'all' with PHREEQC requires "
                "phreeqc_mobile_species to be set (list of mobile species "
                "from the PHREEQC database)."
            )

    else:
        raise ValueError(
            f"Cannot extract species: unknown bgc_module_name '{bgc_module_name}'. "
            f"Expected 'NATIVE_BGC_FLEX' or 'PHREEQC'."
        )


def extract_cycling_frameworks_from_bgc(
        path2selected_NATIVE_BGC_FLEX_framework: str,
) -> List[str]:
    """Extract cycling framework names from a NATIVE_BGC_FLEX JSON file.

    Parameters
    ----------
    path2selected_NATIVE_BGC_FLEX_framework : str
        Path to the NATIVE_BGC_FLEX JSON file.

    Returns
    -------
    list of str
        Framework names (keys of ``CYCLING_FRAMEWORKS``, excluding
        underscore-prefixed metadata keys).
    """
    if not path2selected_NATIVE_BGC_FLEX_framework:
        raise ValueError(
            "extract_cycling_frameworks_from_bgc requires "
            "path2selected_NATIVE_BGC_FLEX_framework to be set."
        )
    if not os.path.isfile(path2selected_NATIVE_BGC_FLEX_framework):
        raise FileNotFoundError(
            f"BGC framework file not found: "
            f"{path2selected_NATIVE_BGC_FLEX_framework}"
        )

    with open(path2selected_NATIVE_BGC_FLEX_framework, 'r') as f:
        bgc_data = json.load(f)

    frameworks_section = bgc_data.get("CYCLING_FRAMEWORKS", {})
    # Filter out metadata keys starting with underscore
    framework_names = [k for k in frameworks_section.keys()
                       if not k.startswith("_")]
    return framework_names


def _get_docker_container_name(docker_compose_path: str) -> str:
    """
    Extract the container_name from docker-compose.yml.
    Returns the name string or 'docker_openwq' as fallback.
    """
    try:
        with open(docker_compose_path, 'r') as f:
            content = f.read()
        match = re.search(r'container_name:\s*(\S+)', content)
        if match:
            return match.group(1).strip()
    except Exception:
        pass
    return 'docker_openwq'


def _query_running_container_mount(container_name: str) -> tuple:
    """
    Query a running Docker container for its bind-mount volume mapping.
    Returns (host_root, container_root) or (None, None) if unavailable.

    This is the most reliable method because it returns the ACTUAL resolved
    mount point, regardless of which docker-compose.yml was used to start
    the container or how many '../' levels were in the volume definition.
    """
    import subprocess
    try:
        result = subprocess.run(
            ['docker', 'inspect', container_name,
             '--format', '{{range .Mounts}}{{.Source}}||{{.Destination}}{{end}}'],
            capture_output=True, text=True, timeout=5
        )
        if result.returncode == 0 and result.stdout.strip():
            # Parse "source||destination" — we use '||' as separator to avoid
            # confusion with ':' in Windows paths
            mount_str = result.stdout.strip()
            parts = mount_str.split('||')
            if len(parts) >= 2:
                host_root = parts[0].strip()
                container_root = parts[1].strip()
                if not host_root.endswith('/'):
                    host_root += '/'
                if not container_root.endswith('/'):
                    container_root += '/'
                return host_root, container_root
    except FileNotFoundError:
        pass  # Docker CLI not installed
    except Exception:
        pass
    return None, None


def _parse_docker_volume_mount(docker_compose_path: str) -> tuple:
    """
    Determine the Docker volume mount mapping (host_root → container_root).
    Returns (host_root, container_root) or (None, None) if resolution fails.

    Strategy:
      1. Query the running container via 'docker inspect' — this returns the
         ACTUAL resolved mount, which is always correct regardless of which
         docker-compose.yml started the container.
      2. Fall back to parsing docker-compose.yml and resolving relative paths
         from its directory (works when the container is not yet running).
    """
    # --- Strategy 1: query the running container ---
    container_name = _get_docker_container_name(docker_compose_path)
    host_root, container_root = _query_running_container_mount(container_name)
    if host_root and container_root:
        return host_root, container_root

    # --- Strategy 2: parse docker-compose.yml ---
    try:
        with open(docker_compose_path, 'r') as f:
            content = f.read()

        # Find volume mount lines matching pattern like '- host_path:container_path:options'
        volume_pattern = re.compile(r'^\s*-\s+([^:]+):([^:]+)(?::.*)?$', re.MULTILINE)
        matches = volume_pattern.findall(content)

        if not matches:
            print("WARNING: No volume mounts found in docker-compose.yml")
            return None, None

        # Use the first volume mount
        host_path_raw, container_path = matches[0]
        host_path_raw = host_path_raw.strip()
        container_path = container_path.strip()

        # Resolve relative host path from docker-compose.yml directory
        compose_dir = os.path.dirname(os.path.abspath(docker_compose_path))
        host_root = os.path.normpath(os.path.join(compose_dir, host_path_raw))

        # Ensure paths end with '/'
        if not host_root.endswith('/'):
            host_root += '/'
        if not container_path.endswith('/'):
            container_path += '/'

        return host_root, container_path

    except FileNotFoundError:
        print(f"WARNING: docker-compose.yml not found at: {docker_compose_path}")
        return None, None
    except Exception as e:
        print(f"WARNING: Failed to parse docker-compose.yml: {e}")
        return None, None


def _correct_path_for_docker(path: str, host_root: str, container_root: str) -> str:
    """
    Replace host root prefix with container root in a file path.
    If the path doesn't start with host_root, return it unchanged.
    """
    if path and path.startswith(host_root):
        return container_root + path[len(host_root):]
    return path


def _correct_dict_paths_for_docker(d: dict, host_root: str, container_root: str) -> dict:
    """
    Recursively correct file paths in a dictionary.
    Looks for string values that contain the host root path and replaces them.
    """
    if d is None:
        return d
    corrected = {}
    for key, value in d.items():
        if isinstance(value, str) and host_root in value:
            corrected[key] = _correct_path_for_docker(value, host_root, container_root)
        elif isinstance(value, dict):
            corrected[key] = _correct_dict_paths_for_docker(value, host_root, container_root)
        elif isinstance(value, list):
            corrected[key] = _correct_list_paths_for_docker(value, host_root, container_root)
        else:
            corrected[key] = value
    return corrected


def _correct_list_paths_for_docker(lst: list, host_root: str, container_root: str) -> list:
    """
    Recursively correct file paths in a list.
    """
    if lst is None:
        return lst
    corrected = []
    for item in lst:
        if isinstance(item, str) and host_root in item:
            corrected.append(_correct_path_for_docker(item, host_root, container_root))
        elif isinstance(item, dict):
            corrected.append(_correct_dict_paths_for_docker(item, host_root, container_root))
        elif isinstance(item, list):
            corrected.append(_correct_list_paths_for_docker(item, host_root, container_root))
        else:
            corrected.append(item)
    return corrected


def Gen_Input_Driver(
        # where to save input files
        dir2save_input_files: str,

        # Top-level settings
        project_name: str,
        geographical_location: str,
        authors: str,
        date: str,
        comment: str,
        hostmodel: str,

        # Computational settings
        solver: str,
        run_mode_debug: bool,
        use_num_threads: Union[int, str],

        # Initial conditions: set equal to all
        ic_all_value: float,
        ic_all_units: str,

        # Biogeochemistry module
        bgc_module_name: str,
        path2selected_NATIVE_BGC_FLEX_framework: str = "",

        # PHREEQC-specific settings (used when bgc_module_name = "PHREEQC")
        phreeqc_input_filepath: str = "",
        phreeqc_database_filepath: str = "",
        phreeqc_mobile_species: Optional[List[str]] = None,
        phreeqc_component_h2o: bool = True,
        phreeqc_temperature_mapping: Optional[Dict[str, str]] = None,
        phreeqc_pressure_mapping: Optional[Dict[str, str]] = None,
        phreeqc_chemical_species_names: Optional[List[str]] = None,

        # Transport Dissolved module
        td_module_name: str = "NATIVE_TD_ADV",
        td_module_dispersion_xyz: List[Union[int, float]] = None,
        td_module_characteristic_length_m: Union[int, float] = 1.0,

        # Lateral Exchange module
        le_module_name: str = "NONE",
        le_module_config: List[Dict[str, Union[str, int, float]]] = None,

        # Transport Sediments module
        ts_module_name: str = "NONE",
        ts_sediment_compartment: str = "RIVER_NETWORK_REACHES",
        ts_transport_compartment: str = "RUNOFF",
        ts_direction: str = "z",
        ts_erosion_inhibit_compartment: str = "ILAYERVOLFRACWAT_SNOW",
        ts_data_format: str = "JSON",
        ts_mmf_defaults: Optional[Dict[str, float]] = None,
        ts_mmf_parameters: Optional[Dict[str, Dict[str, list]]] = None,
        ts_hbvsed_defaults: Optional[Dict[str, float]] = None,
        ts_hbvsed_parameters: Optional[Dict[str, Dict[str, list]]] = None,
        ts_hbvsed_monthly_erosion_factor: Optional[List[float]] = None,

        # Sorption Isotherm module
        si_module_name: str = "NONE",
        si_sediment_compartment: str = "RIVER_NETWORK_REACHES",
        si_bulk_density_kg_m3: float = 1500.0,
        si_layer_thickness_m: float = 1.0,
        si_species_params: Optional[Dict[str, Dict[str, float]]] = None,

        # Sink and Source settings
        ss_method: str = "load_from_csv",
        ss_metadata_source: str = "",
        ss_metadata_comment: str = "",
        ss_method_csv_config: List[Dict[str, Union[str, int]]] = None,
        ss_method_copernicus_basins_hrus: Dict[str, str] = None,
        ss_method_copernicus_compartment_name_for_load: str = "",
        ss_method_copernicus_nc_lc_dir = None,
        ss_method_copernicus_period: List[Union[int, float]] = None,
        ss_method_copernicus_default_loads_bool: bool = True,
        ss_method_copernicus_annual_to_seasonal_loads_method: str = "uniform",

        # External water fluxes
        ewf_method: str = "fixed_value",
        ewf_method_fixedval_comment: str = "",
        ewf_method_fixedval_source: str = "",
        ewf_method_fixedval_chem_name: str = "",
        ewf_method_fixedval_value: str = "",
        ewf_method_fixedval_units: str = "",
        ewf_method_fixedval_external_inputflux_name: str = "",

        # Output settings
        output_format: str = "HDF5",
        chemical_species: List[str] = None,
        units: str = "MG/L",
        no_water_conc_flag: int = -9999,
        export_sediment: bool = False,
        compartments_and_cells: Dict[str, Dict[str, List]] = None,
        timestep: List[Union[int, str]] = None,

        # Climate-adjusted export coefficient parameters
        # (used for based_on_lulc dynamic loads with climate_dependency=True)
        # How the monthly climate is sourced:
        #   "fixed_parameters" → use the hand-entered ss_climate_data dict below
        #   "time_series"      → read precip + temperature from a NetCDF/CSV file
        #                        (ss_climate_data_source), so openWQ stays
        #                        independent of any host model.
        ss_climate_data_type: str = "fixed_parameters",
        ss_climate_data: Optional[Dict[int, Dict[str, List[float]]]] = None,
        # Used only when ss_climate_data_type = "time_series":
        #   {'precip': {file_type:'nc'|'csv', path:..., nc_key_or_column:...},
        #    'temp':   {file_type:'nc'|'csv', path:..., nc_key_or_column:...}}
        ss_climate_data_source: Optional[Dict[str, Dict[str, str]]] = None,
        # Temporal resolution at which the time_series climate adjusts the SS
        # load (and thus the number of SS entries written): one of
        #   "native"  → the climate file's own step (auto-detected)
        #   "hourly"  → 1 entry / hour
        #   "daily"   → 1 entry / day   (default)
        #   "monthly" → 1 entry / month (smallest file; constant within a month)
        # A finer step = more realistic dynamics but a MUCH larger SS file; the
        # driver estimates the size first and refuses/prompts if it's too big
        # (see _SS_ENTRY_GUARD_MAX).
        ss_climate_data_source_adjusting_resolution: str = "daily",
        # Calibration-only: (start, end) of the calibration window.  When set,
        # the SS size guard is evaluated over this window (not the full sim
        # period), because the per-eval SS is later trimmed to it — so a short
        # window keeps 'daily'/'hourly' affordable.  None for a normal run.
        calibration_period=None,
        ss_climate_precip_scaling_power: float = 1.0,
        ss_climate_temp_q10: float = 2.0,
        ss_climate_temp_reference_c: float = 15.0,

        # ML model source/sink parameters
        # (used when ss_method = "ml_model")
        ss_ml_training_data_csv: Optional[str] = None,
        ss_ml_model_type: str = "xgboost",
        ss_ml_target_species: Optional[List[str]] = None,
        ss_ml_feature_columns: Optional[List[str]] = None,
        ss_ml_n_estimators: int = 200,
        ss_ml_max_depth: int = 6,

        # Optional parameters (MUST be at the end)
        ss_method_copernicus_optional_custom_annual_load_coeffs_per_lulc_class: Optional[Dict[int, Dict[str, float]]] = None,

        # Skip the "inputs already exist" prompt and regenerate unconditionally
        # (set force_regenerate = True in your config for non-interactive reruns).
        force_regenerate: bool = False,

        # Accept extra kwargs so model_config_template.py can pass locals() cleanly
        **kwargs

) -> None:

    print(f"Project: {project_name}")
    print(f"Authors: {authors}")

    # Calibration / non-interactive reruns pass force_regenerate=True.  Use that
    # to suppress interactive "reuse vs rebuild?" prompts (LULC cache in
    # Gen_SS_Driver, GRQA cache in Gen_Report) so they are NOT asked once per
    # evaluation — only a standalone template run (force_regenerate=False) prompts.
    os.environ['OPENWQ_SUPPRESS_PROMPTS'] = '1' if force_regenerate else '0'

    # ── Guard: don't silently overwrite existing input files ────────────────
    # Regenerating re-runs the (slow) GRQA / Copernicus / BGC processing and
    # overwrites the openWQ inputs in dir2save_input_files.  If they already
    # exist, ask before proceeding.  Skip with force_regenerate=True, or in a
    # non-interactive shell (where we default to NOT touching existing inputs).
    import sys as _sys
    # Key only on INPUT artifacts (openwq_out is model OUTPUT, not input).
    _existing = [n for n in ("openwq_in", "openWQ_master.json",
                             "openwq_baseline_manifest.json")
                 if os.path.exists(os.path.join(dir2save_input_files, n))]
    if _existing and not force_regenerate:
        print(f"\n⚠  openWQ input files already exist in:\n     {dir2save_input_files}")
        print(f"     (found: {', '.join(_existing)})")
        print("   Regenerating will OVERWRITE them and re-run the slow "
              "GRQA / Copernicus processing.")
        if _sys.stdin.isatty():
            _ans = input("   Delete & regenerate the input files? [y/N]: ").strip().lower()
        else:
            _ans = ""   # non-interactive: default to keeping existing inputs
            print("   (non-interactive shell — set force_regenerate=True to regenerate)")
        if _ans not in ("y", "yes"):
            print("   ↳ Keeping existing inputs; skipping generation.\n")
            return
        print("   ↳ Regenerating…\n")

    # Docker path correction setup
    # This script runs on the HOST but generates JSON files with paths that will be
    # used INSIDE the Docker container.  Always enabled — model_config_template runs
    # with Docker only (Apptainer is handled separately by the calibration workflow).
    _docker_host_root = None
    _docker_container_root = None

    # Locate docker-compose.yml relative to this script
    # This script is in: .../supporting_scripts/Model_Config/config_support_lib/
    # docker-compose.yml is in: .../containers/
    script_dir = os.path.dirname(os.path.abspath(__file__))
    docker_compose_path = os.path.normpath(
        os.path.join(script_dir, '..', '..', '..', 'containers', 'docker-compose.yml'))

    _docker_host_root, _docker_container_root = _parse_docker_volume_mount(docker_compose_path)

    if _docker_host_root and _docker_container_root:
        print(f"  Docker mode: mapping host '{_docker_host_root}' -> container '{_docker_container_root}'")

        # Correct paths that are ONLY embedded in JSON content (NOT used for file I/O)
        # These are paths that OpenWQ reads at runtime inside the container.
        #
        # NOTE: The following paths are NOT corrected here because they are used
        # for file I/O during config generation (the script reads/copies these files):
        #   - path2selected_NATIVE_BGC_FLEX_framework (read by Load_BGQmodule_file.py)
        #   - phreeqc_input_filepath / phreeqc_database_filepath (copied by Gen_PHREEQCmodule_file.py)
        #   - ss_method_copernicus_basins_hrus (shapefile read by Gen_SS_Driver.py)
        #   - ss_method_copernicus_nc_lc_dir (NetCDF files read by Gen_SS_Driver.py)
        #
        # Only CSV source/sink Filepath entries are purely embedded in JSON (the CSV files
        # themselves are read by OpenWQ at runtime, not by the config generator).

        # Correct CSV source/sink file paths (embedded in SS JSON, read at runtime by OpenWQ)
        if ss_method_csv_config is not None:
            ss_method_csv_config = _correct_list_paths_for_docker(
                ss_method_csv_config, _docker_host_root, _docker_container_root)
    else:
        # The docker-compose mapping is ONLY used to rewrite CSV source/sink
        # Filepaths (the only absolute paths embedded in the generated JSON).
        # openWQ's config otherwise uses relative paths, and Apptainer/Singularity
        # runs (e.g. the HPC calibration) bind host->container at run time — so
        # there is nothing to correct unless CSV source/sink entries are used.
        # Only warn when that mapping is actually needed.
        if ss_method_csv_config:
            print("WARNING: Could not determine Docker path mapping from "
                  "docker-compose.yml. CSV source/sink paths will NOT be corrected.")

    # =============================================
    # Cross-module compatibility validation
    # =============================================
    # PHREEQC natively handles sorption via SURFACE/EXCHANGE keywords.
    # If the SI module (Freundlich/Langmuir) is also active, both will
    # write sorption fluxes to d_chemass_dt_chem at runtime, causing
    # double-counting of sorption processes.
    if (bgc_module_name.upper() == "PHREEQC"
            and si_module_name.upper() not in ("NONE", "")):
        raise ValueError(
            f"\n\n"
            f"  ERROR: Incompatible module combination detected!\n"
            f"  ================================================\n"
            f"  BIOGEOCHEMISTRY = '{bgc_module_name}'\n"
            f"  SORPTION_ISOTHERM = '{si_module_name}'\n\n"
            f"  PHREEQC natively handles sorption processes via SURFACE and\n"
            f"  EXCHANGE blocks in the .pqi input file. Running both PHREEQC\n"
            f"  and the Freundlich/Langmuir SI module will double-count\n"
            f"  sorption on the chemistry derivative array (d_chemass_dt_chem).\n\n"
            f"  FIX: Set si_module_name = \"NONE\" when using PHREEQC,\n"
            f"  or use bgc_module_name = \"NATIVE_BGC_FLEX\" if you need the\n"
            f"  simplified Freundlich/Langmuir isotherms.\n"
        )

    # Set defaults for mutable args
    if td_module_dispersion_xyz is None:
        td_module_dispersion_xyz = [0.0, 0.0, 0.0]
    if le_module_config is None:
        le_module_config = []
    if ss_method_csv_config is None:
        ss_method_csv_config = []
    if chemical_species is None:
        chemical_species = ["NO3-N"]

    # ── Resolve chemical_species = "all" ───────────────────────────────
    if (isinstance(chemical_species, str) and chemical_species.lower() == "all") or \
       (isinstance(chemical_species, list) and len(chemical_species) == 1
        and isinstance(chemical_species[0], str) and chemical_species[0].lower() == "all"):

        chemical_species = extract_species_from_bgc(
            bgc_module_name=bgc_module_name,
            path2selected_NATIVE_BGC_FLEX_framework=path2selected_NATIVE_BGC_FLEX_framework,
            phreeqc_mobile_species=phreeqc_mobile_species,
        )

        print(f"\n  ── Output Species (auto-detected from {bgc_module_name}) ──")
        print(f"     chemical_species = 'all' → resolved to {len(chemical_species)} species:")
        for i, sp in enumerate(chemical_species, 1):
            print(f"       {i}. {sp}")
        print()
    if timestep is None:
        timestep = [1, "hour"]

    # Add comment lines at the top
    json_header_comment = [
        "",
        "// #####################",
        "// JSON file automatically generated using OpenWQ's supporting scripts",
        f"// {datetime.datetime.now().replace(microsecond=0)}\n"
        "// #####################",
        ""
    ]

    # Generate paths for file I/O (host filesystem — where files are actually written)
    master_file_fullpath = os.path.join(f'{dir2save_input_files}', "openWQ_master.json")
    general_json_input_dir = os.path.join(f'{dir2save_input_files}', "openwq_in")
    config_file_fullpath = os.path.join(f'{general_json_input_dir}', "openWQ_config.json")
    bgc_config_filepath = os.path.join(f'{general_json_input_dir}', f"openWQ_MODULE_{bgc_module_name}.json")
    td_config_filepath = os.path.join(f'{general_json_input_dir}', f"openWQ_MODULE_{td_module_name}.json") if td_module_name != "NONE" else ""
    le_config_filepath = os.path.join(f'{general_json_input_dir}', f"openWQ_MODULE_{le_module_name}.json") if le_module_name != "NONE" else ""
    ts_config_filepath = os.path.join(f'{general_json_input_dir}', f"openWQ_MODULE_{ts_module_name}.json") if ts_module_name != "NONE" else ""
    si_config_filepath = os.path.join(f'{general_json_input_dir}', f"openWQ_MODULE_{si_module_name}.json") if si_module_name != "NONE" else ""
    ss_config_filepath = os.path.join(f'{general_json_input_dir}', f"openWQ_SS_{ss_method}.json")
    ewf_config_filepath = os.path.join(f'{general_json_input_dir}', f"openWQ_EWF_fixed_value.json")
    output_file_fullpath = os.path.join(f'{dir2save_input_files}', "openwq_out/")

    # =============================================
    # RELATIVE PATHS for master JSON content
    # =============================================
    # The master JSON uses RELATIVE paths (openwq_in/..., openwq_out/) so it is
    # portable across Docker, Apptainer, and native execution. The coupled model
    # (OpenWQ_hydrolink.cpp) reads "openWQ_master.json" from the CWD, so all paths
    # inside it resolve relative to where the master JSON lives.
    # This replaces the old Docker absolute-path correction for master JSON content.
    rel_config_file = "openwq_in/openWQ_config.json"
    rel_bgc_config  = f"openwq_in/openWQ_MODULE_{bgc_module_name}.json"
    rel_td_config   = f"openwq_in/openWQ_MODULE_{td_module_name}.json" if td_module_name != "NONE" else ""
    rel_le_config   = f"openwq_in/openWQ_MODULE_{le_module_name}.json" if le_module_name != "NONE" else ""
    rel_ts_config   = f"openwq_in/openWQ_MODULE_{ts_module_name}.json" if ts_module_name != "NONE" else ""
    rel_si_config   = f"openwq_in/openWQ_MODULE_{si_module_name}.json" if si_module_name != "NONE" else ""
    rel_ss_config   = f"openwq_in/openWQ_SS_{ss_method}.json"
    rel_ewf_config  = f"openwq_in/openWQ_EWF_fixed_value.json"
    rel_output      = "openwq_out/"

    ###############
    # Call create_master_json
    ###############

    mJSON_lib.create_master_json(
        json_header_comment=json_header_comment,
        master_file_fullpath=master_file_fullpath,   # Host path: where to write the file
        config_file_fullpath=rel_config_file,        # Relative path: embedded in JSON content
        bgc_config_filepath=rel_bgc_config,
        td_config_filepath=rel_td_config,
        le_config_filepath=rel_le_config,
        ts_config_filepath=rel_ts_config,
        si_config_filepath=rel_si_config,
        ss_config_filepath=rel_ss_config,
        ewf_config_filepath=rel_ewf_config,
        output_file_fullpath=rel_output,
        project_name=project_name,
        geographical_location=geographical_location,
        authors=authors,
        date=date,
        comment=comment,
        solver=solver,
        run_mode_debug=run_mode_debug,
        use_num_threads=use_num_threads,
        bgc_module_name=bgc_module_name,
        td_module_name=td_module_name,
        le_module_name=le_module_name,
        ts_module_name=ts_module_name,
        ts_sediment_compartment=ts_sediment_compartment,
        ts_transport_compartment=ts_transport_compartment,
        si_module_name=si_module_name,
        si_sediment_compartment=si_sediment_compartment,
        ss_method=ss_method,
        ss_metadata_source=ss_metadata_source,
        ewf_method=ewf_method,
        ewf_method_fixedval_source=ewf_method_fixedval_source,
        output_format=output_format,
        chemical_species=chemical_species,
        units=units,
        no_water_conc_flag=no_water_conc_flag,
        export_sediment=export_sediment,
        compartments_and_cells=compartments_and_cells,
        timestep=timestep
    )

    ###############
    # Determine compartment names from host model
    ###############

    if (hostmodel == "mizuroute"):
        compartment_names = ["RIVER_NETWORK_REACHES"]
    elif (hostmodel == "summa"):
        compartment_names = ["SCALARCANOPYWAT",
                             "ILAYERVOLFRACWAT_SNOW",
                             "RUNOFF",
                             "ILAYERVOLFRACWAT_SOIL",
                             "SCALARAQUIFER"]

    ###############
    # Call create_config_json
    # For PHREEQC: species come from the user-provided list (auto-discovered at runtime);
    #              no CYCLING_FRAMEWORK is needed.
    # For NATIVE_BGC_FLEX: species and cycling frameworks are hardcoded from the framework file.
    ###############

    if bgc_module_name == "PHREEQC":

        # PHREEQC species are discovered at runtime from the database/input file.
        # The user provides expected species names for initial conditions.
        if phreeqc_chemical_species_names is None or len(phreeqc_chemical_species_names) == 0:
            print("WARNING: phreeqc_chemical_species_names is empty. "
                  "Config file will have no initial conditions. "
                  "PHREEQC will use default initial conditions from the .pqi file.")
            chem_names = []
        else:
            chem_names = phreeqc_chemical_species_names

        cJSON_lib.create_config_json(
            config_file_fullpath=config_file_fullpath,
            json_header_comment=json_header_comment,
            compartment_names=compartment_names,
            cycling_framework=None,  # PHREEQC does not use cycling frameworks
            chemical_species_names=chem_names,
            ic_all_value=ic_all_value,
            ic_all_units=ic_all_units,
        )
    else:
        # NATIVE_BGC_FLEX path
        # Extract cycling framework names dynamically from the BGC template
        cycling_fw_names = extract_cycling_frameworks_from_bgc(
            path2selected_NATIVE_BGC_FLEX_framework=path2selected_NATIVE_BGC_FLEX_framework
        )
        cJSON_lib.create_config_json(
            config_file_fullpath=config_file_fullpath,
            json_header_comment=json_header_comment,
            compartment_names=compartment_names,
            cycling_framework=cycling_fw_names,
            chemical_species_names=chemical_species,
            ic_all_value=ic_all_value,
            ic_all_units=ic_all_units,
        )

    ###############
    # Call create_td_module_json
    ###############

    if td_module_name != "NONE":
        tdmJSON_lib.create_td_module_json(
            td_config_filepath=td_config_filepath,
            json_header_comment=json_header_comment,
            td_module_name=td_module_name,
            td_module_dispersion_xyz=td_module_dispersion_xyz,
            td_module_characteristic_length_m=td_module_characteristic_length_m
        )

    ###############
    # Call create_le_module_json
    ###############

    if le_module_name != "NONE":
        lemJSON_lib.create_le_module_json(
            le_config_filepath=le_config_filepath,
            json_header_comment=json_header_comment,
            le_module_name=le_module_name,
            le_module_config=le_module_config
        )

    ###############
    # Call create_ts_module_json
    ###############

    if ts_module_name != "NONE":
        tsmJSON_lib.create_ts_module_json(
            ts_config_filepath=ts_config_filepath,
            json_header_comment=json_header_comment,
            ts_module_name=ts_module_name,
            ts_direction=ts_direction,
            ts_erosion_inhibit_compartment=ts_erosion_inhibit_compartment,
            ts_data_format=ts_data_format,
            ts_mmf_defaults=ts_mmf_defaults,
            ts_mmf_parameters=ts_mmf_parameters,
            ts_hbvsed_defaults=ts_hbvsed_defaults,
            ts_hbvsed_parameters=ts_hbvsed_parameters,
            ts_hbvsed_monthly_erosion_factor=ts_hbvsed_monthly_erosion_factor
        )

    ###############
    # Call create_si_module_json
    ###############

    if si_module_name != "NONE":
        simJSON_lib.create_si_module_json(
            si_config_filepath=si_config_filepath,
            json_header_comment=json_header_comment,
            si_module_name=si_module_name,
            si_bulk_density_kg_m3=si_bulk_density_kg_m3,
            si_layer_thickness_m=si_layer_thickness_m
        )

    ###############
    # Call BGC module setup (NATIVE_BGC_FLEX or PHREEQC)
    ###############

    if bgc_module_name == "PHREEQC":
        # Generate PHREEQC module JSON + copy database and input files
        phreeqcJSON_lib.create_phreeqc_module_json(
            bgc_config_filepath=bgc_config_filepath,
            json_header_comment=json_header_comment,
            phreeqc_input_filepath=phreeqc_input_filepath,
            phreeqc_database_filepath=phreeqc_database_filepath,
            phreeqc_mobile_species=phreeqc_mobile_species if phreeqc_mobile_species else ["H"],
            phreeqc_component_h2o=phreeqc_component_h2o,
            phreeqc_temperature_mapping=phreeqc_temperature_mapping,
            phreeqc_pressure_mapping=phreeqc_pressure_mapping,
        )
    elif bgc_module_name == "NATIVE_BGC_FLEX":
        # Copy the selected BGC framework file
        bgqmJSON_lib.load_bgq_module_json(
            json_header_comment=json_header_comment,
            path2selected_NATIVE_BGC_FLEX_framework=path2selected_NATIVE_BGC_FLEX_framework,
            bgc_config_filepath=bgc_config_filepath
        )
    else:
        print(f"WARNING: Unknown bgc_module_name '{bgc_module_name}'. "
              f"Valid options are: NATIVE_BGC_FLEX, PHREEQC")

    ###############
    # Call gen_ss_driver
    ###############

    # based_on_lulc sub-schema (passed from the template via kwargs):
    #   lulc_loads          'static' | 'dynamic'
    #   climate_dependency  True/False — multiply the loads by the climate signal
    # climate weighting is active only for dynamic loads with climate_dependency
    # on; when its source is a time series it uses the sub-annual generation path.
    _lulc_loads = str(kwargs.get("lulc_loads", "static")).lower()
    _climate_dep = bool(kwargs.get("climate_dependency", False))
    _clim_active = (_lulc_loads == "dynamic" and _climate_dep)
    _clim_timeseries = _clim_active and str(ss_climate_data_type).lower() == "time_series"

    # ── Auto-detect Copernicus period from host-model control file ─────
    if ss_method == "based_on_lulc":

        if ss_method_copernicus_period is None:
            # Auto-detect from the host-model control file
            _fm_path = kwargs.get('file_manager_path')
            if _fm_path is None:
                raise ValueError(
                    "ss_method_copernicus_period is None (auto-detect), but "
                    "file_manager_path was not provided. Either set "
                    "ss_method_copernicus_period = [start_year, end_year] "
                    "or provide file_manager_path."
                )
            ss_method_copernicus_period = parse_sim_period_from_control_file(
                file_manager_path=_fm_path,
                hostmodel=hostmodel,
            )
            print(f"\n  ℹ  Copernicus period auto-detected from control file:")
            print(f"     Simulation period: {ss_method_copernicus_period[0]}"
                  f"–{ss_method_copernicus_period[1]}")

        # Report Copernicus availability vs requested period
        req_start = int(ss_method_copernicus_period[0])
        req_end = int(ss_method_copernicus_period[1])
        cop_start = max(req_start, COPERNICUS_LC_YEAR_MIN)
        cop_end = min(req_end, COPERNICUS_LC_YEAR_MAX)

        # User preference: use proxy year when outside Copernicus range?
        _use_proxy = kwargs.get('ss_method_copernicus_use_proxy_if_outside_range', True)

        _is_fully_outside = (req_start > COPERNICUS_LC_YEAR_MAX
                             or req_end < COPERNICUS_LC_YEAR_MIN)
        _is_partially_outside = (req_start < COPERNICUS_LC_YEAR_MIN
                                 or req_end > COPERNICUS_LC_YEAR_MAX)

        print(f"\n  ── Copernicus ESA CCI Land Cover ──")
        print(f"     Available range:  {COPERNICUS_LC_YEAR_MIN}–{COPERNICUS_LC_YEAR_MAX}")
        print(f"     Requested period: {req_start}–{req_end}")
        print(f"     Use proxy year if outside range: {_use_proxy}")

        if _is_fully_outside:
            cop_start = min(max(req_start, COPERNICUS_LC_YEAR_MIN), COPERNICUS_LC_YEAR_MAX)
            cop_end = cop_start
            if _use_proxy:
                print(f"     ⚠  WARNING: Simulation period ({req_start}–{req_end}) "
                      f"is entirely outside Copernicus range "
                      f"({COPERNICUS_LC_YEAR_MIN}–{COPERNICUS_LC_YEAR_MAX}).")
                print(f"        → Using nearest available year ({cop_start}) as proxy.")
                print(f"          (ss_method_copernicus_use_proxy_if_outside_range = True)")
            else:
                raise ValueError(
                    f"\n  ✗  ABORTED: Simulation period ({req_start}–{req_end}) "
                    f"is entirely outside Copernicus range "
                    f"({COPERNICUS_LC_YEAR_MIN}–{COPERNICUS_LC_YEAR_MAX}).\n"
                    f"     Proxy years are disabled "
                    f"(ss_method_copernicus_use_proxy_if_outside_range = False).\n"
                    f"     Options:\n"
                    f"       1. Set ss_method_copernicus_use_proxy_if_outside_range = True "
                    f"to use the nearest available year as proxy.\n"
                    f"       2. Set ss_method_copernicus_period = [{COPERNICUS_LC_YEAR_MIN}, "
                    f"{COPERNICUS_LC_YEAR_MAX}] to use the full Copernicus range.\n"
                    f"       3. Use ss_method = 'load_from_csv' to provide your own load data."
                )
        elif _is_partially_outside:
            if _use_proxy:
                print(f"     ⚠  Note: Simulation extends beyond Copernicus range.")
                print(f"        → Using Copernicus data for: {cop_start}–{cop_end}")
                print(f"          Years outside {COPERNICUS_LC_YEAR_MIN}–{COPERNICUS_LC_YEAR_MAX} "
                      f"will use the nearest available year as proxy.")
                print(f"          (ss_method_copernicus_use_proxy_if_outside_range = True)")
            else:
                raise ValueError(
                    f"\n  ✗  ABORTED: Simulation period ({req_start}–{req_end}) "
                    f"extends beyond Copernicus range "
                    f"({COPERNICUS_LC_YEAR_MIN}–{COPERNICUS_LC_YEAR_MAX}).\n"
                    f"     Proxy years are disabled "
                    f"(ss_method_copernicus_use_proxy_if_outside_range = False).\n"
                    f"     Options:\n"
                    f"       1. Set ss_method_copernicus_use_proxy_if_outside_range = True "
                    f"to use the nearest available year for out-of-range periods.\n"
                    f"       2. Set ss_method_copernicus_period = [{cop_start}, {cop_end}] "
                    f"to restrict to the available range.\n"
                    f"       3. Use ss_method = 'load_from_csv' to provide your own load data."
                )
        else:
            print(f"     ✓  Full coverage: Copernicus data available for "
                  f"{cop_start}–{cop_end}")

        # Clamp the period to what Copernicus actually has
        ss_method_copernicus_period = [cop_start, cop_end]
        print()

    # Build BGC engine label for SS METADATA traceability
    if bgc_module_name == "NATIVE_BGC_FLEX":
        _bgc_engine_label = f"NATIVE_BGC_FLEX: {path2selected_NATIVE_BGC_FLEX_framework}"
    elif bgc_module_name == "PHREEQC":
        _bgc_engine_label = "PHREEQC"
    else:
        _bgc_engine_label = bgc_module_name or ""
    _chem_species_list = list(chemical_species) if isinstance(chemical_species, (list, tuple)) else [chemical_species]

    if (ss_method == "load_from_csv"):

        ssJSON_lib.set_ss_from_csv(
            ss_config_filepath=ss_config_filepath,
            json_header_comment=json_header_comment,
            ss_metadata_source=ss_metadata_source,
            ss_metadata_comment=ss_metadata_comment,
            ss_method_csv_config=ss_method_csv_config,
            bgc_engine_label=_bgc_engine_label,
            chemical_species_list=_chem_species_list,
        )
    elif ss_method == "based_on_lulc" and not _clim_timeseries:

        ssJSON_lib.set_ss_from_copernicus_lulc_with_loads(
            ss_config_filepath=ss_config_filepath,
            json_header_comment=json_header_comment,
            ss_metadata_source=ss_metadata_source,
            ss_metadata_comment=ss_metadata_comment,
            ss_method_copernicus_basins_hrus=ss_method_copernicus_basins_hrus,
            ss_method_copernicus_nc_lc_dir=ss_method_copernicus_nc_lc_dir,
            ss_method_copernicus_compartment_name_for_load=ss_method_copernicus_compartment_name_for_load,
            ss_method_copernicus_period=ss_method_copernicus_period,
            ss_method_copernicus_default_loads_bool=ss_method_copernicus_default_loads_bool,
            optional_load_coefficients=ss_method_copernicus_optional_custom_annual_load_coeffs_per_lulc_class,
            ss_method_copernicus_annual_to_seasonal_loads_method=("seasonal" if _clim_active else "uniform"),
            simulation_period=[req_start, req_end],
            bgc_engine_label=_bgc_engine_label,
            chemical_species_list=_chem_species_list,
            climate_data=ss_climate_data,
            precip_scaling_power=ss_climate_precip_scaling_power if ss_climate_data else 1.0,
            temp_q10=ss_climate_temp_q10 if ss_climate_data else 2.0,
            temp_reference_c=ss_climate_temp_reference_c if ss_climate_data else 15.0,
        )
    elif ss_method == "based_on_lulc" and _clim_timeseries:

        # Resolve the monthly climate (precip + temperature) per the selected
        # source: a hand-entered fixed dict, or a NetCDF/CSV time series read
        # straight from file (host-model independent).
        _daily_climate = None
        _subannual_climate = None
        _subannual_unit = "hour"
        if str(ss_climate_data_type).lower() == "time_series":
            if not ss_climate_data_source:
                raise ValueError(
                    "ss_climate_data_type='time_series' requires ss_climate_data_source "
                    "= {'precip': {file_type, path, nc_key_or_column}, "
                    "'temp': {file_type, path, nc_key_or_column}}."
                )
            _years = list(range(int(req_start), int(req_end) + 1))

            # Pick the SS temporal resolution, but FIRST estimate the resulting
            # file size and refuse/downscale if it would be too big for openWQ
            # (a per-step SS file over many units x species x years explodes —
            # a 1M-entry / 100+ MB JSON OOM-kills the run).
            # Size the guard to the calibration window when one is set: the
            # per-eval SS is trimmed to that window before openWQ reads it, so
            # a short window can make a fine resolution affordable even when the
            # full period would not.  Window year-span is clamped to [1, full].
            _guard_years = len(_years)
            _scope = "full simulation period"
            if calibration_period:
                try:
                    _ws, _we = calibration_period
                    _wy0, _wy1 = int(str(_ws)[:4]), int(str(_we)[:4])
                    _guard_years = min(len(_years), max(1, _wy1 - _wy0 + 1))
                    _scope = f"calibration window {_ws}..{_we}"
                except (TypeError, ValueError):
                    pass
            _res = _ss_resolution_size_guard(
                resolution=ss_climate_data_source_adjusting_resolution,
                n_years=_guard_years,
                n_units=_count_basin_spatial_units(ss_method_copernicus_basins_hrus),
                n_species=len(_chem_species_list or []) or 1,
                climate_source=ss_climate_data_source,
                scope_label=_scope)

            # Monthly climate is always built — it's both the proxy-year fallback
            # and the 'monthly' resolution itself.
            _resolved_climate = ssJSON_lib.build_climate_data_from_timeseries(
                ss_climate_data_source, years=_years)
            if _res == "monthly":
                print(f"  Climate (dynamic SS): time series, MONTHLY resolution "
                      f"({len(_resolved_climate)} year(s) "
                      f"[{min(_resolved_climate)}-{max(_resolved_climate)}])")
            elif _res == "daily":
                _daily_climate = ssJSON_lib.build_daily_climate_from_timeseries(
                    ss_climate_data_source, years=_years)
                _ndays = sum(len(v.get('day', [])) for v in _daily_climate.values())
                print(f"  Climate (dynamic SS): time series, DAILY resolution "
                      f"({_ndays} days; SS adjusted every day)")
            else:  # 'hourly' or 'native'
                _subannual_climate = ssJSON_lib.build_subannual_climate_from_timeseries(
                    ss_climate_data_source, years=_years, resolution=_res)
                _subannual_unit = "hour"
                _nper = sum(len(v.get('hour', [])) for v in _subannual_climate.values())
                print(f"  Climate (dynamic SS): time series, {_res.upper()} resolution "
                      f"({_nper} periods; SS adjusted every step)")
        else:
            if ss_climate_data is None:
                raise ValueError(
                    "based_on_lulc dynamic loads with climate_dependency=True and "
                    "ss_climate_data_type='fixed_parameters' requires ss_climate_data.\n"
                    "Provide {year: {'precip_mm': [12], 'temp_c': [12]}}, or set "
                    "ss_climate_data_type='time_series' + ss_climate_data_source."
                )
            _resolved_climate = ss_climate_data

        ssJSON_lib.set_ss_climate_adjusted_export_coefficients(
            ss_config_filepath=ss_config_filepath,
            json_header_comment=json_header_comment,
            ss_method_copernicus_basins_hrus=ss_method_copernicus_basins_hrus,
            ss_method_copernicus_nc_lc_dir=ss_method_copernicus_nc_lc_dir,
            ss_method_copernicus_period=ss_method_copernicus_period,
            ss_method_copernicus_default_loads_bool=ss_method_copernicus_default_loads_bool,
            ss_method_copernicus_compartment_name_for_load=ss_method_copernicus_compartment_name_for_load,
            climate_data=_resolved_climate,
            precip_scaling_power=ss_climate_precip_scaling_power,
            temp_q10=ss_climate_temp_q10,
            temp_reference_c=ss_climate_temp_reference_c,
            optional_load_coefficients=ss_method_copernicus_optional_custom_annual_load_coeffs_per_lulc_class,
            simulation_period=[req_start, req_end],
            bgc_engine_label=_bgc_engine_label,
            chemical_species_list=_chem_species_list,
            climate_data_type=ss_climate_data_type,
            climate_data_source=ss_climate_data_source,
            daily_climate_data=_daily_climate,
            subannual_climate_data=_subannual_climate,
            subannual_unit=_subannual_unit,
        )

    elif (ss_method == "ml_model"):

        import Gen_MLmodel_SS as mlJSON_lib

        if ss_ml_training_data_csv is None:
            raise ValueError(
                "ss_method='ml_model' requires ss_ml_training_data_csv.\n"
                "Provide path to a CSV with columns: date, discharge_m3s, precip_mm, temp_c, <species_concentration>"
            )

        mlJSON_lib.train_and_generate_ss_json(
            ss_config_filepath=ss_config_filepath,
            json_header_comment=json_header_comment,
            training_data_csv=ss_ml_training_data_csv,
            model_type=ss_ml_model_type,
            target_species=ss_ml_target_species,
            feature_columns=ss_ml_feature_columns,
            n_estimators=ss_ml_n_estimators,
            max_depth=ss_ml_max_depth,
            ss_metadata_comment=ss_metadata_comment,
            ss_metadata_source=ss_metadata_source,
            compartment_name=ss_method_copernicus_compartment_name_for_load,
        )

    elif (ss_method == "none"):
        print("  Skipping sink/source generation (ss_method='none')")

    else:
        print(
            f"WARNING: The SS method '{ss_method}' is unknown or not available for automatic generation. "
            f"Available methods: 'load_from_csv', 'based_on_lulc', 'ml_model', 'none'")

    ###############
    # Call gen_ewf_driver
    ###############
    if (ewf_method == "fixed_value"):
        ewfJSON_lib.set_ewf_fixed_value(
            ewf_config_filepath=ewf_config_filepath,
            json_header_comment=json_header_comment,
            ewf_method_fixedval_comment=ewf_method_fixedval_comment,
            ewf_method_fixedval_source=ewf_method_fixedval_source,
            ewf_method_fixedval_chem_name=ewf_method_fixedval_chem_name,
            ewf_method_fixedval_value=ewf_method_fixedval_value,
            ewf_method_fixedval_units=ewf_method_fixedval_units,
            ewf_method_fixedval_external_inputflux_name=ewf_method_fixedval_external_inputflux_name
        )
    elif (ewf_method == "none"):
        print("  Skipping external water flux generation (ewf_method='none')")
    else:
        print(
            f"WARNING: The EWF method '{ewf_method}' is unknown or not available for automatic generation. "
            f"Available methods: 'fixed_value', 'none'")
