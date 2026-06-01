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
Calibration Setup Report
=========================

Provides two report generators:

1. ``generate_setup_report()`` — static HTML report (used by the generated
   calibration script after the user has finalised settings).

2. ``generate_interactive_setup()`` — **interactive** HTML report where all
   calibration options are presented as form elements (dropdowns, editable
   tables, checkboxes).  The user adjusts settings in the browser and then
   downloads a ready-to-run Python calibration script.
"""

import os
import json
import inspect
import html as html_lib
import logging
from datetime import datetime
from typing import List, Dict, Any, Optional

from . import report_helpers as rh
from . import config_integration as _ci


def _caller_template_stem(default="calibration"):
    """Basename (no extension) of the *calibration template* that called
    into this module — used to name the setup + results reports after it.

    Walks the stack outward and returns the first frame whose file lives
    OUTSIDE this calibration_lib package (i.e. the user's
    ``calibration_config_template*.py``).  Falls back to ``default``."""
    try:
        _here = os.path.dirname(os.path.abspath(__file__))
        for fr in inspect.stack()[1:]:
            fn = fr.filename or ""
            if not fn or fn.startswith("<"):
                continue
            if os.path.dirname(os.path.abspath(fn)) == _here:
                continue
            stem = os.path.splitext(os.path.basename(fn))[0]
            if stem:
                return stem
    except Exception:
        pass
    return default

logger = logging.getLogger(__name__)


def generate_setup_report(
    output_dir: str,
    model_config: Dict[str, Any],
    calibration_parameters: List[Dict],
    calibration_settings: Dict[str, Any],
    observation_config: Dict[str, Any],
    container_config: Dict[str, Any],
) -> Optional[str]:
    """
    Generate the calibration setup HTML report.

    Parameters
    ----------
    output_dir : str
        Directory where the report HTML file will be saved.
    model_config : Dict
        Loaded model configuration variables.
    calibration_parameters : List[Dict]
        Final list of calibration parameters (auto-extracted + overrides + additional).
    calibration_settings : Dict
        Calibration settings (algorithm, max_evaluations, objective_function, etc.).
    observation_config : Dict
        Observation data configuration (source, paths, species, etc.).
    container_config : Dict
        Container runtime config (executable path, file manager, mpi_np, etc.).

    Returns
    -------
    Optional[str]
        Path to the generated HTML report, or None on failure.
    """
    try:
        os.makedirs(output_dir, exist_ok=True)
        report_path = os.path.join(output_dir, "calibration_setup_report.html")

        H = []  # HTML content accumulator

        # ── Head ──
        H.append(rh.build_html_head("OpenWQ Calibration Setup Report"))
        H.append("<body>")
        H.append('<div class="layout">')

        # ── Sidebar ──
        nav_items = [
            {"id": "summary", "label": "Summary"},
            {"id": "parameters", "label": "Calibration Parameters"},
            {"id": "settings", "label": "Settings"},
            {"id": "model-config", "label": "Model Configuration"},
            {"id": "observations", "label": "Observations"},
            {"id": "run-commands", "label": "Run Commands"},
        ]
        H.append(rh.build_sidebar(nav_items, logo_text="OpenWQ Calibration"))

        H.append('<div class="main">')

        # ── Header ──
        project_name = model_config.get("project_name", "Calibration")
        authors = model_config.get("authors", "")
        date_str = model_config.get("date", datetime.now().strftime("%B %Y"))

        meta_items = [
            f"Authors: {authors}" if authors else "",
            f"Date: {date_str}",
            f"Algorithm: {calibration_settings.get('algorithm', 'DDS')}",
            f"Parameters: {len(calibration_parameters)}",
        ]
        meta_items = [m for m in meta_items if m]

        H.append(rh.build_header(
            title=project_name,
            subtitle="Calibration Setup Report",
            meta_items=meta_items,
            badge_text="SETUP",
            badge_class="badge-accent"
        ))

        H.append('<div class="container">')

        # ── Section: Summary ──
        H.append(_build_summary_section(
            calibration_parameters, calibration_settings, model_config
        ))

        # ── Section: Calibration Parameters ──
        H.append(_build_parameters_section(calibration_parameters))

        # ── Section: Calibration Settings ──
        H.append(_build_settings_section(calibration_settings))

        # ── Section: Model Configuration ──
        H.append(_build_model_config_section(model_config))

        # ── Section: Observation Data ──
        H.append(_build_observations_section(observation_config))

        # ── Section: Run Commands ──
        H.append(_build_run_commands_section(
            model_config, container_config, calibration_settings
        ))

        H.append('</div>')  # container
        H.append(rh.build_footer())
        H.append('</div>')  # main
        H.append('</div>')  # layout

        H.append(rh.build_theme_toggle_js())
        H.append("</body></html>")

        # Write report
        with open(report_path, 'w', encoding='utf-8') as f:
            f.write("\n".join(H))

        logger.info(f"Setup report saved to: {report_path}")
        return report_path

    except Exception as e:
        logger.error(f"Failed to generate setup report: {e}")
        import traceback
        logger.debug(traceback.format_exc())
        return None


# =========================================================================
# Section Builders
# =========================================================================

def _build_summary_section(
    parameters: List[Dict],
    settings: Dict[str, Any],
    model_config: Dict[str, Any]
) -> str:
    """Build the summary KPI section."""

    # Count parameter sources
    auto_count = sum(1 for p in parameters if p.get("source", "").startswith("auto"))
    user_count = len(parameters) - auto_count

    # Count BGC vs other parameters
    bgc_count = sum(1 for p in parameters if p.get("file_type") == "bgc_json")
    other_count = len(parameters) - bgc_count

    kpis = [
        {"icon": "\U0001f4ca", "value": str(len(parameters)), "label": "Parameters"},
        {"icon": "\U0001f3af", "value": settings.get("algorithm", "DDS"), "label": "Algorithm"},
        {"icon": "\U0001f504", "value": str(settings.get("max_evaluations", 500)), "label": "Max Evaluations"},
        {"icon": "\U0001f4cf", "value": settings.get("objective_function", "KGE"), "label": "Objective"},
    ]

    kpi_html = rh.build_kpi_grid(kpis)

    # Info boxes
    info_parts = []
    if bgc_count > 0:
        info_parts.append(f"{bgc_count} BGC reaction rate parameters")
    if other_count > 0:
        info_parts.append(f"{other_count} non-BGC parameters (transport, SS, etc.)")
    if auto_count > 0:
        info_parts.append(f"{auto_count} auto-extracted from template")
    if user_count > 0:
        info_parts.append(f"{user_count} user-defined")

    info_html = rh.build_highlight_box(
        "<strong>Parameter breakdown:</strong> " + " &bull; ".join(info_parts),
        "info"
    )

    # Sensitivity analysis info
    sa_html = ""
    if settings.get("run_sensitivity_first", False):
        sa_method = settings.get("sensitivity_method", "morris").upper()
        sa_html = rh.build_highlight_box(
            f"<strong>Sensitivity analysis</strong> ({sa_method}) will run "
            f"before calibration to identify influential parameters.",
            "warning"
        )

    bgc_name = model_config.get("bgc_module_name", "")
    bgc_template = os.path.basename(
        model_config.get("path2selected_NATIVE_BGC_FLEX_framework", "")
    )

    model_info = rh.build_highlight_box(
        f"<strong>Model:</strong> {model_config.get('hostmodel', '')} "
        f"&bull; <strong>BGC:</strong> {bgc_name} "
        + (f"(<code>{bgc_template}</code>)" if bgc_template else ""),
        "success"
    )

    return f"""
<div class="section" id="summary">
    <h2>Summary</h2>
    {kpi_html}
    {info_html}
    {model_info}
    {sa_html}
</div>
"""


def _build_parameters_section(parameters: List[Dict]) -> str:
    """Build the calibration parameters table section."""

    if not parameters:
        return """
<div class="section" id="parameters">
    <h2>Calibration Parameters</h2>
    <div class="card"><p>No calibration parameters defined.</p></div>
</div>
"""

    # Format functions
    def fmt_num(v):
        if isinstance(v, (int, float)):
            if abs(v) < 0.01 or abs(v) >= 10000:
                return f"{v:.4g}"
            return f"{v:.4f}"
        return str(v)

    def fmt_bounds(v):
        if isinstance(v, (list, tuple)) and len(v) == 2:
            return f"[{fmt_num(v[0])}, {fmt_num(v[1])}]"
        return str(v)

    def fmt_source(v):
        v = str(v)
        if "auto" in v.lower():
            return '<span class="badge badge-secondary">auto</span>'
        return '<span class="badge badge-primary">user</span>'

    columns = [
        {"key": "name", "label": "Parameter", "align": "left"},
        {"key": "initial", "label": "Initial", "align": "right",
         "format": fmt_num},
        {"key": "bounds", "label": "Bounds", "align": "right",
         "format": fmt_bounds},
        {"key": "transform", "label": "Transform", "align": "left"},
        {"key": "units", "label": "Units", "align": "left"},
        {"key": "description", "label": "Description", "align": "left"},
        {"key": "source", "label": "Source", "align": "left",
         "format": fmt_source},
    ]

    table_html = rh.build_param_table(parameters, columns, table_id="param-table")

    # Group by file_type for summary
    type_counts = {}
    for p in parameters:
        ft = p.get("file_type", "unknown")
        type_counts[ft] = type_counts.get(ft, 0) + 1

    type_summary = " &bull; ".join(
        f"{ft}: {count}" for ft, count in sorted(type_counts.items())
    )

    return f"""
<div class="section" id="parameters">
    <h2>Calibration Parameters ({len(parameters)})</h2>
    <p style="color:var(--text2);margin-bottom:1rem;">
        Parameter types: {type_summary}
    </p>
    {table_html}
</div>
"""


def _build_settings_section(settings: Dict[str, Any]) -> str:
    """Build the calibration settings section."""

    # Main settings card
    rows = [
        ("Algorithm", settings.get("algorithm", "DDS")),
        ("Max Evaluations", str(settings.get("max_evaluations", 500))),
        ("Objective Function", settings.get("objective_function", "KGE")),
        ("Temporal Resolution", settings.get("temporal_resolution", "native")),
        ("Aggregation Method", settings.get("aggregation_method", "mean")),
        ("Random Seed", str(settings.get("random_seed", "None"))),
    ]

    settings_rows = "".join(
        f"<tr><td><strong>{k}</strong></td><td>{v}</td></tr>"
        for k, v in rows
    )

    # Target species
    targets = settings.get("calibration_targets", {})
    species = targets.get("species", [])
    reaches = targets.get("reach_ids", "all")
    compartments = targets.get("compartments", [])

    species_str = ", ".join(species) if isinstance(species, list) else str(species)
    comp_str = ", ".join(compartments) if isinstance(compartments, list) else str(compartments)

    # Weights
    weights = settings.get("objective_weights", {})
    weights_str = ", ".join(
        f"{k}: {v}" for k, v in weights.items()
    ) if weights else "Equal"

    target_rows = f"""
    <tr><td><strong>Target Species</strong></td><td>{html_lib.escape(species_str)}</td></tr>
    <tr><td><strong>Target Reaches</strong></td><td>{reaches}</td></tr>
    <tr><td><strong>Compartments</strong></td><td>{html_lib.escape(comp_str)}</td></tr>
    <tr><td><strong>Objective Weights</strong></td><td>{html_lib.escape(weights_str)}</td></tr>
    """

    # Sensitivity analysis settings
    sa_html = ""
    if settings.get("run_sensitivity_first", False):
        sa_method = settings.get("sensitivity_method", "morris")
        sa_rows = [
            ("Method", sa_method.upper()),
        ]
        if sa_method == "morris":
            sa_rows.append(("Trajectories", str(settings.get("sensitivity_morris_trajectories", 10))))
            sa_rows.append(("Levels", str(settings.get("sensitivity_morris_levels", 4))))
        elif sa_method == "sobol":
            sa_rows.append(("Samples", str(settings.get("sensitivity_sobol_samples", 1024))))
        sa_rows.append(("Threshold", str(settings.get("sensitivity_threshold", 0.1))))

        sa_table_rows = "".join(
            f"<tr><td><strong>{k}</strong></td><td>{v}</td></tr>"
            for k, v in sa_rows
        )

        sa_html = f"""
        <div class="card secondary">
            <h3>Sensitivity Analysis</h3>
            <table class="param-table" style="max-width:400px;">
                <tbody>{sa_table_rows}</tbody>
            </table>
        </div>
        """

    return f"""
<div class="section" id="settings">
    <h2>Calibration Settings</h2>
    <div class="card primary">
        <h3>Optimization</h3>
        <table class="param-table" style="max-width:500px;">
            <tbody>{settings_rows}</tbody>
        </table>
    </div>
    <div class="card primary">
        <h3>Calibration Targets</h3>
        <table class="param-table" style="max-width:500px;">
            <tbody>{target_rows}</tbody>
        </table>
    </div>
    {sa_html}
</div>
"""


def _build_model_config_section(
    model_config: Dict[str, Any],
    species_obs_availability: Optional[Dict[str, Dict[str, Any]]] = None,
) -> str:
    """Build the model configuration summary section."""

    bgc = model_config.get("bgc_module_name", "N/A")
    td = model_config.get("td_module_name", "N/A")
    le = model_config.get("le_module_name", "N/A")
    ts = model_config.get("ts_module_name", "N/A")
    si = model_config.get("si_module_name", "N/A")
    ss = model_config.get("ss_method", "N/A")
    solver = model_config.get("solver", "N/A")
    hostmodel = model_config.get("hostmodel", "N/A")

    species = model_config.get("chemical_species", [])
    if isinstance(species, list):
        species_str = ", ".join(species)
    else:
        species_str = str(species)

    def _module_badge(name):
        if name in ("NONE", "none", "N/A", ""):
            return f'<span class="badge badge-none">{name}</span>'
        return f'<span class="badge badge-secondary">{name}</span>'

    bgc_template = os.path.basename(
        model_config.get("path2selected_NATIVE_BGC_FLEX_framework", "")
    )

    rows = [
        ("Host Model", hostmodel),
        ("Solver", solver),
        ("BGC Module", f"{_module_badge(bgc)}" + (f" &mdash; <code>{bgc_template}</code>" if bgc_template else "")),
        ("Transport Dissolved", _module_badge(td)),
        ("Lateral Exchange", _module_badge(le)),
        ("Sediment Transport", _module_badge(ts)),
        ("Sorption Isotherm", _module_badge(si)),
        ("Source/Sink Method", _module_badge(ss)),
        ("Chemical Species", species_str),
    ]

    table_rows = "".join(
        f"<tr><td><strong>{k}</strong></td><td>{v}</td></tr>"
        for k, v in rows
    )

    ic_value = model_config.get("ic_all_value", "N/A")
    ic_units = model_config.get("ic_all_units", "N/A")

    return f"""
<div class="section" id="model-config">
    <h2>Model Configuration</h2>
    <div class="card">
        <h3>Module Summary</h3>
        <table class="param-table" style="max-width:700px;">
            <tbody>{table_rows}</tbody>
        </table>
    </div>
    <div class="card">
        <h3>Initial Conditions</h3>
        <p>Uniform IC: <strong>{ic_value} {ic_units}</strong> (all species, all compartments)</p>
    </div>
</div>
"""


def _candidate_bgc_data_sources(
    model_config: Dict[str, Any]
) -> List[str]:
    """Return BGC source-file candidates in priority order.

    Handles BOTH ``NATIVE_BGC_FLEX`` (JSON) and ``PHREEQC`` (``.pqi``)
    modules so the calibratability check and the diagram work for any
    template the user selects.

    Order:
      1. The *runtime* combined BGC JSON at
         ``<exe_dir>/openwq_in/openWQ_MODULE_NATIVE_BGC_FLEX.json`` —
         present after the user has run ``Gen_Input_Driver``.
      2. The user's individual NATIVE_BGC_FLEX template returned by
         ``config_integration.get_bgc_template_path`` — available
         before the model has been run.
      3. The PHREEQC ``.pqi`` file pointed to by
         ``model_config["phreeqc_input_filepath"]`` (when the user
         selected ``bgc_module_name == "PHREEQC"``).
    """
    out: List[str] = []
    exe_path = model_config.get("executable_path", "")
    if exe_path:
        out.append(os.path.join(
            os.path.dirname(os.path.abspath(exe_path)),
            "openwq_in", "openWQ_MODULE_NATIVE_BGC_FLEX.json",
        ))
    try:
        tpl = _ci.get_bgc_template_path(model_config)
    except Exception:
        tpl = None
    if tpl:
        out.append(tpl)
    # PHREEQC: look up the .pqi file from the model config.  Resolve
    # relative paths the same way get_bgc_template_path does (relative
    # to the model_config dir).
    pqi_path = model_config.get("phreeqc_input_filepath", "") or \
               model_config.get("phreeqc_pqi_path", "")
    if pqi_path:
        if not os.path.isabs(pqi_path):
            # Try relative to executable_path's grandparent, mirroring
            # the layout used by Gen_Input_Driver.
            if exe_path:
                _candidate = os.path.join(
                    os.path.dirname(os.path.abspath(exe_path)),
                    pqi_path,
                )
                if os.path.isfile(_candidate):
                    out.append(os.path.abspath(_candidate))
                    pqi_path = ""   # mark as handled
            if pqi_path:
                # Try relative to 1_Model_Config (where pqi templates ship)
                _here = os.path.dirname(os.path.abspath(__file__))
                _config_dir = os.path.normpath(
                    os.path.join(_here, "..", "..", "1_Model_Config")
                )
                _candidate = os.path.join(_config_dir, pqi_path)
                if os.path.isfile(_candidate):
                    out.append(os.path.abspath(_candidate))
                    pqi_path = ""
        if pqi_path and os.path.isfile(pqi_path):
            out.append(os.path.abspath(pqi_path))
    # Deduplicate while preserving order.
    seen: set = set()
    uniq: List[str] = []
    for p in out:
        if p and p not in seen:
            seen.add(p)
            uniq.append(p)
    return uniq


def _parse_bgc_source(path: str) -> Optional[dict]:
    """Dispatch a BGC source path to the right parser based on extension.

    * ``.json`` → :func:`report_helpers.safe_load_bgc_json` →
      :func:`extract_bgc_network`.
    * ``.pqi``  → :func:`report_helpers.extract_bgc_network_from_pqi`
      (PHREEQC).

    Returns the network dict (``{"species": [...], "frameworks": ...}``)
    or ``None`` when the file is unreadable / produces an empty network.
    Centralising the dispatch here means ``_resolve_bgc_data_source``,
    ``_load_bgc_network`` and ``_build_bgc_diagram`` cannot disagree
    about which path is "usable".
    """
    if not path:
        return None
    ext = os.path.splitext(path)[1].lower()
    try:
        if ext == ".pqi":
            net = rh.extract_bgc_network_from_pqi(path)
        else:
            data = rh.safe_load_bgc_json(path)
            net = rh.extract_bgc_network(data) if data is not None else None
    except Exception:
        return None
    if not net or not net.get("frameworks"):
        return None
    return net


def _resolve_bgc_data_source(
    model_config: Dict[str, Any]
) -> Optional[str]:
    """Return the first candidate path whose contents actually PARSE.

    Before this used to return the first candidate that merely existed
    on disk — which broke on the runtime JSON because it begins with a
    ``// ...`` comment header that vanilla ``json.load`` rejects.  Now
    each candidate is round-tripped through :func:`_parse_bgc_source`
    so PHREEQC ``.pqi`` files are also accepted and an unparseable
    runtime JSON falls through to the user's template.
    """
    for path in _candidate_bgc_data_sources(model_config):
        if _parse_bgc_source(path) is not None:
            return path
    return None


def _load_bgc_network(
    model_config: Dict[str, Any],
) -> Optional[dict]:
    """Parse the best available BGC source (JSON or PHREEQC) into a
    network dict.

    Walks the candidate list and returns the first one that yields a
    non-empty network.  Format dispatch lives in
    :func:`_parse_bgc_source`.
    """
    for path in _candidate_bgc_data_sources(model_config):
        net = _parse_bgc_source(path)
        if net:
            return net
    return None


def _build_bgc_diagram(
    model_config: Dict[str, Any],
    species_obs_availability: Optional[Dict[str, Dict[str, Any]]] = None,
) -> str:
    """Build the BGC reaction network diagram HTML.

    Uses :func:`_resolve_bgc_data_source` to find the best source file
    (runtime combined JSON, NATIVE_BGC_FLEX template, or PHREEQC
    ``.pqi``) and routes it to the right
    :func:`report_helpers.generate_bgc_reaction_diagram` argument so
    BOTH file formats render.
    """
    path = _resolve_bgc_data_source(model_config)
    if not path:
        return ""
    try:
        kwargs: Dict[str, Any] = dict(
            module_name=model_config.get("bgc_module_name", ""),
            species_obs_availability=species_obs_availability,
        )
        if path.lower().endswith(".pqi"):
            kwargs["pqi_filepath"] = path
        else:
            kwargs["bgc_template_path"] = path
        return rh.generate_bgc_reaction_diagram(**kwargs)
    except Exception:
        return ""


def _build_observations_section(obs_config: Dict[str, Any]) -> str:
    """Build the observation data section."""

    source = obs_config.get("source", "skip")

    if source == "skip":
        return f"""
<div class="section" id="observations">
    <h2>Observation Data</h2>
    {rh.build_highlight_box(
        "Observation data is set to <strong>skip</strong>. "
        "No model-observation comparison will be performed.",
        "warning"
    )}
</div>
"""

    rows = [("Source", source.upper())]

    if source == "grqa":
        grqa_path = obs_config.get("grqa_local_data_path", "N/A")
        buffer_km = obs_config.get("grqa_buffer_km", 100)
        rows.append(("GRQA Data Path", str(grqa_path)))
        rows.append(("Buffer Distance", f"{buffer_km} km"))
    elif source == "user_csv":
        csv_path = obs_config.get("user_observation_csv", "N/A")
        rows.append(("CSV Path", str(csv_path)))

    species = obs_config.get("chemical_species", [])
    if isinstance(species, list) and species:
        rows.append(("Species", ", ".join(species)))

    table_rows = "".join(
        f"<tr><td><strong>{k}</strong></td><td>{v}</td></tr>"
        for k, v in rows
    )

    return f"""
<div class="section" id="observations">
    <h2>Observation Data</h2>
    <div class="card">
        <table class="param-table" style="max-width:600px;">
            <tbody>{table_rows}</tbody>
        </table>
    </div>
</div>
"""


def _build_run_commands_section(
    model_config: Dict[str, Any],
    container_config: Dict[str, Any],
    calibration_settings: Dict[str, Any]
) -> str:
    """Build the run commands section with code snippets."""

    # Determine the calibration command
    run_cmd = "python calibration_config_template.py"
    resume_cmd = "python calibration_config_template.py --resume"
    sa_cmd = "python calibration_config_template.py --sensitivity-only"
    dry_run_cmd = "python calibration_config_template.py --dry-run"
    show_params_cmd = "python calibration_config_template.py --show-parameters"

    # Docker model execution info
    exe_path = container_config.get("executable_path", "")
    file_mgr = container_config.get("file_manager_path", "")
    mpi_np = container_config.get("mpi_np", 2)

    run_code = f"""# Run calibration
{run_cmd}

# Resume from checkpoint (if interrupted)
{resume_cmd}

# Run sensitivity analysis only
{sa_cmd}

# Dry run (validate configuration)
{dry_run_cmd}

# Show auto-extracted parameters
{show_params_cmd}"""

    # Best-params application snippet
    best_params_code = """# After calibration completes, the best parameters are saved in:
#   <calibration_work_dir>/results/best_parameters.json
#
# To apply best parameters to your model configuration:
import json
with open("results/best_parameters.json") as f:
    best_params = json.load(f)
print("Best parameters:", best_params)"""

    return f"""
<div class="section" id="run-commands">
    <h2>Run Commands</h2>
    <div class="card primary">
        <h3>Calibration Commands</h3>
        {rh.build_code_block(run_code)}
    </div>
    <div class="card secondary">
        <h3>After Calibration</h3>
        {rh.build_code_block(best_params_code, "python")}
    </div>
</div>
"""


# =========================================================================
# Interactive Setup Report
# =========================================================================

def generate_interactive_setup(
    output_dir: str,
    model_config: Dict[str, Any],
    auto_extracted_parameters: List[Dict],
    observation_config: Dict[str, Any],
    container_config: Dict[str, Any],
    model_config_path: str,
    calibration_work_dir: str,
    module_parameters: Optional[Dict[str, List[Dict]]] = None,
    module_selections: Optional[Dict[str, Any]] = None,
    species_obs_availability: Optional[Dict[str, Dict[str, Any]]] = None,
    ss_species_with_loads: Optional[set] = None,
) -> Optional[str]:
    """
    Generate an **interactive** calibration setup HTML report.

    The report contains form elements (dropdowns, editable tables,
    checkboxes) that the user adjusts in the browser.  A live-updating
    code block shows the generated Python calibration script.  The user
    clicks "Download Script" to get a ready-to-run ``.py`` file.

    Parameters
    ----------
    output_dir : str
        Directory where the report HTML file will be saved.
    model_config : Dict
        Loaded model configuration variables.
    auto_extracted_parameters : List[Dict]
        Auto-extracted calibration parameters from BGC template.
    observation_config : Dict
        Observation data configuration.
    container_config : Dict
        Container runtime config.
    model_config_path : str
        Absolute path to the user's model_config_template.py.
    calibration_work_dir : str
        Path to the calibration workspace directory.
    module_parameters : Dict[str, List[Dict]], optional
        Parameters grouped by module. If provided, renders collapsible
        module groups instead of a flat parameter table.
    module_selections : Dict[str, Any], optional
        Summary of active modules (from config_integration.get_module_selections).
    species_obs_availability : Dict[str, Dict[str, Any]], optional
        Per-species observation data availability from
        ``config_integration.get_species_observation_availability()``.

    Returns
    -------
    Optional[str]
        Path to the generated HTML report, or None on failure.
    """
    try:
        os.makedirs(output_dir, exist_ok=True)
        # Name the setup report after the calibration template that
        # invoked us, e.g. calibration_config_template_X.py ->
        #   <stem>_config_report.html
        # The matching calibration RESULTS report (produced later by the
        # generated run script) becomes <stem>_results_report.html.
        _calib_stem = _caller_template_stem(default="calibration")
        _setup_report_name = f"{_calib_stem}_config_report.html"
        report_path = os.path.join(output_dir, _setup_report_name)

        H = []  # HTML accumulator

        # ── Head (with interactive layout CSS) ──
        H.append(rh.build_html_head(
            "OpenWQ Calibration — Interactive Setup",
            extra_css=rh.get_interactive_layout_css(),
        ))
        H.append("<body>")

        # ── Collect metadata ──
        project_name = model_config.get("project_name", "Calibration Setup")
        total_params = (sum(len(v) for v in module_parameters.values())
                        if module_parameters
                        else len(auto_extracted_parameters))

        # ── Split Layout ──
        H.append('<div class="split-layout">')

        # ── Top Bar ──
        H.append('<div class="top-bar">')
        H.append(rh.build_compact_header(
            title='Open<span style="color:var(--secondary);">WQ</span>'
                  ' &mdash; Interactive Calibration Setup',
            badge_text="INTERACTIVE",
            badge_class="badge-accent",
        ))
        H.append("""
<div class="top-bar-right">
    <span class="progress-heading">Setup checklist</span>
    <div class="progress-segments" id="progressBar"
         title="Configuration completeness: hover each segment for details">
        <span class="progress-seg" data-seg="1" title="Species: at least 1 selected"></span>
        <span class="progress-seg" data-seg="2" title="Algorithm & objective function set"></span>
        <span class="progress-seg" data-seg="3" title="Parameters: at least 1 enabled"></span>
        <span class="progress-seg" data-seg="4" title="Reach IDs specified"></span>
        <span class="progress-seg" data-seg="5" title="Bounds valid (min < max)"></span>
    </div>
    <span class="progress-label" id="progressLabel"
          title="Configuration completeness check"></span>
</div>
""")
        H.append('</div>')  # top-bar

        # ── Panes Row (flex container for left + resizer + right) ──
        H.append('<div class="panes-row">')

        # ── Config Pane (Left) ──
        H.append('<div class="config-pane" id="configPane">')

        # Tab bar
        H.append("""
<div class="tab-bar">
    <button class="tab-btn active" data-tab="overview" onclick="switchTab('overview')">Overview</button>
    <button class="tab-btn" data-tab="settings" onclick="switchTab('settings')">Settings</button>
    <button class="tab-btn" data-tab="targets" onclick="switchTab('targets')">Targets</button>
    <button class="tab-btn" data-tab="parameters" onclick="switchTab('parameters')">Parameters</button>
</div>
""")

        # ── Tab: Overview ──
        H.append('<div class="tab-panel active" data-tab="overview">')
        H.append(_build_interactive_summary(
            auto_extracted_parameters, model_config, observation_config,
            total_params=total_params,
            module_selections=module_selections,
        ))
        H.append(_build_model_config_section(
            model_config,
            species_obs_availability=species_obs_availability,
        ))
        H.append(_build_observations_section(observation_config))
        H.append('</div>')

        # ── Tab: Settings ──
        H.append('<div class="tab-panel" data-tab="settings">')
        H.append(_build_interactive_settings_section(
            container_config.get("container_runtime", "docker")))
        H.append(_build_interactive_sensitivity_section())
        H.append('</div>')

        # ── Tab: Targets ──
        species_list = model_config.get("chemical_species", [])
        if isinstance(species_list, (list, tuple)):
            species_list = list(species_list)
        else:
            species_list = [str(species_list)]
        obs_source = observation_config.get("source", "skip")

        # Host-model-aware reach/HRU + compartment options
        _spatial = _ci.get_spatial_mapping(model_config)
        _hostmodel = (_spatial.get("hostmodel") or "").lower()
        if _hostmodel == "summa":
            _feat_label = "HRU"
            _feat_shp = observation_config.get("basin_shapefile")
            _feat_key = _spatial.get("basin_mapping_key") or "HRU_ID"
        else:
            _feat_label = "Reach"
            _feat_shp = observation_config.get("river_network_shapefile")
            _feat_key = _spatial.get("river_network_mapping_key") or "SegId"
        _available_reaches = _read_feature_ids(_feat_shp, _feat_key)

        # Compartments: keys defined in the model config are authoritative;
        # fall back to host-aware defaults when not present.
        _cc = model_config.get("compartments_and_cells", {})
        if isinstance(_cc, dict) and _cc:
            _available_compartments = list(_cc.keys())
        elif _hostmodel == "summa":
            _available_compartments = [
                "SCALARCANOPYWAT", "ILAYERVOLFRACWAT_SNOW", "RUNOFF",
                "ILAYERVOLFRACWAT_SOIL", "SCALARAQUIFER",
            ]
        else:
            _available_compartments = ["RIVER_NETWORK_REACHES"]

        # Pre-tick the compartments the user already restricts obs to;
        # otherwise tick all available.
        _obs_comp = model_config.get("observation_compartments") or []
        if isinstance(_obs_comp, str):
            _obs_comp = [_obs_comp]
        _obs_comp_set = set(_obs_comp)
        if _obs_comp_set:
            _selected_compartments = [c for c in _available_compartments
                                      if c in _obs_comp_set] or list(_available_compartments)
        else:
            _selected_compartments = list(_available_compartments)

        H.append('<div class="tab-panel" data-tab="targets">')
        H.append(_build_interactive_targets_section(
            species_list, species_obs_availability or {}, obs_source,
            hostmodel=_hostmodel,
            feature_label=_feat_label,
            available_reaches=_available_reaches,
            available_compartments=_available_compartments,
            selected_compartments=_selected_compartments,
        ))
        H.append('</div>')

        # ── Tab: Parameters ──
        H.append('<div class="tab-panel" data-tab="parameters">')
        # Reaction network diagram (above parameter list)
        _param_diagram = _build_bgc_diagram(
            model_config,
            species_obs_availability=species_obs_availability,
        )
        if _param_diagram:
            H.append('<div class="card" style="margin-bottom:1rem;">')
            H.append('<h3>Reaction Network</h3>')
            H.append(_param_diagram)
            H.append('</div>')
        H.append('<input type="text" class="param-search" id="paramSearch" '
                 'placeholder="Search parameters by name, module, or type..." '
                 'oninput="filterParams(this.value)"/>')
        if module_parameters:
            # Load BGC network for the calibratability check.  Uses the
            # SAME path resolver as the BGC diagram above so the two
            # sides can never disagree about which frameworks are
            # "active" (have at least one observed species).  Before
            # the helper existed the diagram had a template-fallback
            # but the parameter loader didn't, so on a first-time run
            # every BGC parameter ended up grayed out even though the
            # diagram correctly showed reactions with observed species.
            _bgc_net = (_load_bgc_network(model_config)
                        if species_obs_availability else None)
            H.append(_build_interactive_parameters_section_grouped(
                module_parameters, module_selections or {},
                species_obs_availability=species_obs_availability,
                bgc_network=_bgc_net,
                ss_species_with_loads=ss_species_with_loads,
            ))
        else:
            H.append(_build_interactive_parameters_section(
                auto_extracted_parameters
            ))
        H.append('</div>')

        H.append('</div>')  # config-pane

        # ── Resizer Handle ──
        H.append('<div class="pane-resizer" id="paneResizer"></div>')

        # ── Script Pane (Right) ──
        # Recommended save path: the openWQ "3_Calibration" folder, because
        # the generated script does `from calibration_lib...` and that
        # package lives there.  (The script also bakes this folder into its
        # sys.path, so it still runs if saved elsewhere.)  Results are
        # written to calibration_work_dir as set in the calibration config
        # template.  Filename follows the originating template stem
        # (e.g. "model_config_template_mizuRoute_run.py").
        cal_script_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        _run_script_name = f"{_calib_stem}_run.py"
        save_hint = html_lib.escape(os.path.join(cal_script_dir, _run_script_name))

        H.append('<div class="script-pane">')
        H.append(f"""
<div class="script-pane-header">
    <h3>Generated Script</h3>
    <div style="display:flex;gap:.4rem;align-items:center;">
        <button class="copy-btn" id="wrapToggleBtn" style="position:static;font-size:.75rem;padding:.25rem .6rem;
            background:rgba(0,102,204,.1);border:1px solid var(--border);color:var(--primary);
            border-radius:6px;cursor:pointer;"
            onclick="toggleWrap()" title="Toggle line wrapping / horizontal scroll">Wrap: on</button>
        <button class="copy-btn" style="position:static;font-size:.75rem;padding:.25rem .6rem;
            background:rgba(0,102,204,.1);border:1px solid var(--border);color:var(--primary);
            border-radius:6px;cursor:pointer;"
            onclick="copyScript()">Copy</button>
        <button id="downloadBtnHeader" style="position:static;font-size:.75rem;padding:.25rem .6rem;
            background:linear-gradient(135deg,rgba(0,102,204,.15),rgba(0,168,107,.15));
            border:1px solid var(--border);color:var(--secondary);
            border-radius:6px;cursor:pointer;font-weight:600;"
            onclick="downloadScript()">Save the script</button>
    </div>
</div>
<div style="padding:.2rem 1.2rem .1rem;font-size:.7rem;color:var(--text3);word-break:break-all;line-height:1.45;"
     title="{save_hint}">
    Save to: <code id="saveHint" style="font-size:.68rem;">{save_hint}</code>
    <button id="copyHintBtn" onclick="copySaveHint()" title="Copy path to clipboard"
        style="margin-left:.4rem;font-size:.62rem;padding:.05rem .4rem;
        background:rgba(0,102,204,.1);border:1px solid var(--border);
        color:var(--primary);border-radius:4px;cursor:pointer;
        font-family:inherit;white-space:nowrap;vertical-align:baseline;">Copy</button>
</div>
""")
        # Build the how-to section with OS-aware commands
        import platform
        os_name = platform.system()  # "Darwin", "Linux", "Windows"

        containers_dir = os.path.normpath(
            os.path.join(cal_script_dir, '..', '..', '..', '..', 'containers'))
        # The run-script is recommended to live in cal_script_dir (the
        # "3_Calibration" folder that contains calibration_lib).  The how-to
        # therefore cd's into that folder first, then runs the script by name
        # so the `from calibration_lib...` import always resolves.
        results_path = os.path.join(calibration_work_dir,
                                    f'{_calib_stem}_results_report.html')

        if os_name == "Windows":
            docker_cmd = html_lib.escape(
                f'cd /d "{containers_dir}" && docker compose up -d')
            cd_cmd = html_lib.escape(f'cd /d "{cal_script_dir}"')
            run_cmd = html_lib.escape(f'python "{_run_script_name}"')
            resume_cmd = html_lib.escape(f'python "{_run_script_name}" --resume')
            dryrun_cmd = html_lib.escape(f'python "{_run_script_name}" --dry-run')
            open_cmd = "start"
        else:
            docker_cmd = html_lib.escape(
                f"cd {containers_dir} && docker compose up -d")
            cd_cmd = html_lib.escape(f"cd {cal_script_dir}")
            run_cmd = html_lib.escape(f"python {_run_script_name}")
            resume_cmd = html_lib.escape(f"python {_run_script_name} --resume")
            dryrun_cmd = html_lib.escape(f"python {_run_script_name} --dry-run")
            open_cmd = "open" if os_name == "Darwin" else "xdg-open"

        # Inline CSS for the mini code snippets with copy buttons.
        # align-items:flex-start + wrapping code so the full command is
        # visible (paths can be long) with the copy button pinned top-right.
        snippet_css = (
            "display:flex;align-items:flex-start;"
            "background:var(--dark);color:#e2e8f0;border-radius:6px;"
            "padding:.35rem .6rem;margin:.3rem 0 .5rem;font-family:'JetBrains Mono',monospace;"
            "font-size:.75rem;line-height:1.4;gap:.5rem;"
        )
        snippet_code_css = (
            "min-width:0;flex:1;white-space:pre-wrap;"
            "overflow-wrap:anywhere;word-break:break-all;"
        )
        copy_btn_css = (
            "flex-shrink:0;background:rgba(255,255,255,.1);"
            "border:1px solid rgba(255,255,255,.2);color:#e2e8f0;"
            "padding:.15rem .4rem;border-radius:4px;font-size:.65rem;"
            "cursor:pointer;font-family:inherit;white-space:nowrap;"
        )

        H.append(f"""
<details style="margin:.6rem 1rem 0;border:1px solid var(--border);border-radius:8px;padding:.2rem .8rem;font-size:.82rem;">
  <summary style="cursor:pointer;font-weight:600;padding:.4rem 0;color:var(--primary);user-select:none;">
    How to use this script
  </summary>
  <div style="padding:.3rem 0 .8rem;line-height:1.65;color:var(--text2);">

    <div style="margin-bottom:.6rem;">
      <strong>1.</strong> Configure calibration settings, species, and parameters
      using the tabs on the left.
    </div>

    <div style="margin-bottom:.6rem;">
      <strong>2.</strong> Save the script using the <em>Save the script</em> button above.
    </div>

    <div style="margin-bottom:.6rem;">
      <strong>3.</strong> <span id="step3Text">Start the Docker container:</span>
      <div id="step3DockerCmd" style="{snippet_css}">
        <code id="cmdDocker" style="{snippet_code_css}">{docker_cmd}</code>
        <button style="{copy_btn_css}"
          onclick="navigator.clipboard.writeText(document.getElementById('cmdDocker').textContent);this.textContent='Copied!';setTimeout(()=>this.textContent='Copy',1500)">Copy</button>
      </div>
    </div>

    <div style="margin-bottom:.6rem;">
      <strong>4.</strong> Go to the calibration folder (the one that contains
      <code>calibration_lib</code>) so the script's imports resolve:
      <div style="{snippet_css}">
        <code id="cmdCd" style="{snippet_code_css}">{cd_cmd}</code>
        <button style="{copy_btn_css}"
          onclick="navigator.clipboard.writeText(document.getElementById('cmdCd').textContent);this.textContent='Copied!';setTimeout(()=>this.textContent='Copy',1500)">Copy</button>
      </div>
    </div>

    <div style="margin-bottom:.6rem;">
      <strong>5.</strong> Run the calibration:
      <div style="{snippet_css}">
        <code id="cmdRun" style="{snippet_code_css}">{run_cmd}</code>
        <button style="{copy_btn_css}"
          onclick="navigator.clipboard.writeText(document.getElementById('cmdRun').textContent);this.textContent='Copied!';setTimeout(()=>this.textContent='Copy',1500)">Copy</button>
      </div>
    </div>

    <div style="margin-bottom:.6rem;">
      <strong>6.</strong> Resume if interrupted:
      <div style="{snippet_css}">
        <code id="cmdResume" style="{snippet_code_css}">{resume_cmd}</code>
        <button style="{copy_btn_css}"
          onclick="navigator.clipboard.writeText(document.getElementById('cmdResume').textContent);this.textContent='Copied!';setTimeout(()=>this.textContent='Copy',1500)">Copy</button>
      </div>
    </div>

    <div style="margin-bottom:.6rem;">
      <strong>7.</strong> Validate config (optional dry run):
      <div style="{snippet_css}">
        <code id="cmdDryrun" style="{snippet_code_css}">{dryrun_cmd}</code>
        <button style="{copy_btn_css}"
          onclick="navigator.clipboard.writeText(document.getElementById('cmdDryrun').textContent);this.textContent='Copied!';setTimeout(()=>this.textContent='Copy',1500)">Copy</button>
      </div>
    </div>

    <div style="margin-bottom:.3rem;">
      <strong>8.</strong> View results &mdash; the script auto-generates an interactive
      results report and opens it in your browser. To reopen it:
      <div style="{snippet_css}">
        <code id="cmdResults" style="{snippet_code_css}">{open_cmd} {html_lib.escape(results_path)}</code>
        <button style="{copy_btn_css}"
          onclick="navigator.clipboard.writeText(document.getElementById('cmdResults').textContent);this.textContent='Copied!';setTimeout(()=>this.textContent='Copy',1500)">Copy</button>
      </div>
    </div>

  </div>
</details>
""")

        # ── HPC (Apptainer / Singularity) folder checklist ──
        # Derive the concrete folders the user must copy to the cluster from
        # the paths already in their model config.  The heavy GRQA /
        # Copernicus processing is reused from the model-config setup, so the
        # raw GRQA database is NOT needed on the HPC.
        _supp_dir = os.path.dirname(cal_script_dir)  # supporting_scripts (1_/2_/3_)
        _exe = (container_config.get("executable_path", "")
                or model_config.get("executable_path", ""))
        _model_run_dir = (os.path.dirname(os.path.abspath(_exe)) if _exe
                          else model_config.get("dir2save_input_files", "") or "")
        _river_shp = observation_config.get("river_network_shapefile") or ""
        _basin_shp = observation_config.get("basin_shapefile") or ""
        _shp_dirs = [os.path.dirname(p) for p in (_river_shp, _basin_shp) if p]
        _dom_parts = [d for d in ([_model_run_dir] + _shp_dirs) if d]
        try:
            _domain_dir = (os.path.commonpath(_dom_parts)
                           if len(_dom_parts) > 1
                           else (_dom_parts[0] if _dom_parts else ""))
        except ValueError:
            _domain_dir = _model_run_dir

        def _hpc_row(path_html, why_html):
            return (
                '<tr>'
                '<td style="vertical-align:top;padding:.25rem .5rem;">'
                f'<code style="font-size:.72rem;word-break:break-all;">{path_html}</code></td>'
                '<td style="padding:.25rem .5rem;font-size:.82rem;color:var(--text2);">'
                f'{why_html}</td></tr>')

        _hpc_rows = [
            _hpc_row(html_lib.escape(_supp_dir or "(supporting_scripts)"),
                     "<strong>Supporting scripts</strong> &mdash; the openWQ "
                     "<code>1_Model_Config</code>, <code>2_Read_Outputs</code> and "
                     "<code>3_Calibration</code> tree (calibration_lib, the model-config "
                     "template, BGC templates, spatial_matching)."),
            _hpc_row(html_lib.escape(_domain_dir or "(domain folder)"),
                     "<strong>Domain folder</strong> &mdash; the model run dir (the "
                     "<em>already-processed</em> openWQ config in <code>openwq_in/</code> + the "
                     "clipped GRQA in <code>grqa_clipped_data/</code>), the "
                     "<code>shapefiles/</code>, and the host-model settings / forcing / "
                     "control file the executable reads."),
            _hpc_row("&lt;your openwq.sif&gt;",
                     "<strong>Apptainer / Singularity image</strong> &mdash; build/transfer it "
                     "and set its path in the <em>Apptainer SIF path</em> field (Execution settings)."),
            _hpc_row(html_lib.escape(calibration_work_dir or "(calibration_work_dir)"),
                     "<strong>Calibration work dir</strong> &mdash; created on the HPC; "
                     "evaluations + results are written here (no need to pre-copy)."),
        ]

        H.append(f"""
<details style="margin:.2rem 1rem .6rem;border:1px solid var(--border);border-radius:8px;padding:.2rem .8rem;font-size:.82rem;">
  <summary style="cursor:pointer;font-weight:600;padding:.4rem 0;color:var(--primary);user-select:none;">
    Running on HPC (Apptainer / Singularity) &mdash; folders to copy
  </summary>
  <div style="padding:.3rem 0 .8rem;line-height:1.6;color:var(--text2);">
    <p style="margin:.2rem 0 .5rem;">Copy these to the cluster. The heavy GRQA / Copernicus
    processing is <strong>reused</strong> from your model-config setup &mdash; nothing is
    re-downloaded or re-processed per evaluation:</p>
    <table style="width:100%;border-collapse:collapse;font-size:.8rem;">
      <tr style="border-bottom:1px solid var(--border);">
        <th style="text-align:left;padding:.25rem .5rem;">Folder / file</th>
        <th style="text-align:left;padding:.25rem .5rem;">Why</th></tr>
      {''.join(_hpc_rows)}
    </table>
    <div class="hint" style="margin-top:.6rem;">
      <strong>Paths:</strong> keep the same directory layout on the HPC, <em>or</em> edit the
      absolute paths in your model config (<code>executable_path</code>, <code>file_manager_path</code>,
      <code>river_network_shapefile</code>, the basin shapefile, <code>dir2save_input_files</code>) to the
      HPC locations and set the <strong>Apptainer bind path</strong> to the shared root that contains them.
    </div>
    <div class="hint" style="margin-top:.4rem;">
      <strong>Not needed on the HPC:</strong> the raw GRQA database
      (<code>grqa_local_data_path</code>, often &gt;1&nbsp;GB) &mdash; the calibration reuses the
      already-clipped GRQA inside the model run dir.
    </div>
    <div class="hint" style="margin-top:.4rem;">
      Finally, set <strong>Container runtime = apptainer</strong> and fill the
      <strong>Apptainer SIF / bind path</strong> fields above.
    </div>
  </div>
</details>
""")

        H.append("""
<div class="script-pane-body">
    <pre class="code-block" id="scriptBlock">
<code id="scriptCode">Loading...</code>
    </pre>
</div>
""")
        H.append('</div>')  # script-pane
        H.append('</div>')  # panes-row

        # ── Action Bar (Bottom) ──
        H.append('<div class="action-bar">')
        H.append('<span class="action-status" id="scriptStatus"></span>')
        H.append("""
<button class="theme-btn-compact" onclick="toggleTheme()" title="Toggle theme">
    <svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
        <circle cx="12" cy="12" r="5"/><line x1="12" y1="1" x2="12" y2="3"/>
        <line x1="12" y1="21" x2="12" y2="23"/><line x1="4.22" y1="4.22" x2="5.64" y2="5.64"/>
        <line x1="18.36" y1="18.36" x2="19.78" y2="19.78"/><line x1="1" y1="12" x2="3" y2="12"/>
        <line x1="21" y1="12" x2="23" y2="12"/><line x1="4.22" y1="19.78" x2="5.64" y2="18.36"/>
        <line x1="18.36" y1="5.64" x2="19.78" y2="4.22"/>
    </svg>
</button>
""")
        H.append('</div>')  # action-bar

        H.append('</div>')  # split-layout

        # ── Theme toggle JS ──
        H.append(rh.build_theme_toggle_js())

        # ── Main interactive JS ──
        H.append(_build_interactive_js(
            model_config_path=model_config_path,
            calibration_work_dir=calibration_work_dir,
            auto_extracted_parameters=auto_extracted_parameters,
            module_parameters=module_parameters,
            container_config=container_config,
            report_stem=_calib_stem,
        ))

        H.append("</body></html>")

        with open(report_path, 'w', encoding='utf-8') as f:
            f.write("\n".join(H))

        logger.info(f"Interactive setup report saved to: {report_path}")
        return report_path

    except Exception as e:
        logger.error(f"Failed to generate interactive setup report: {e}")
        import traceback
        logger.debug(traceback.format_exc())
        return None


# =========================================================================
# Interactive Section Builders
# =========================================================================

def _build_interactive_summary(
    parameters: List[Dict],
    model_config: Dict[str, Any],
    observation_config: Dict[str, Any],
    total_params: int = 0,
    module_selections: Optional[Dict[str, Any]] = None,
) -> str:
    """Summary KPIs for the interactive report."""
    if not total_params:
        total_params = len(parameters)
    bgc_template = os.path.basename(
        model_config.get("path2selected_NATIVE_BGC_FLEX_framework", "")
    )
    hostmodel = model_config.get("hostmodel", "N/A")
    obs_source = observation_config.get("source", "skip").upper()

    # Count active modules
    n_modules = 0
    if module_selections:
        for key in ("bgc_active", "td_active", "ts_active",
                    "le_active", "si_active", "ss_active"):
            if module_selections.get(key):
                n_modules += 1

    kpis = [
        {"icon": "\U0001f4ca", "value": str(total_params), "label": "Calibration Params"},
        {"icon": "\U0001f9ea", "value": str(n_modules) if n_modules else "N/A",
         "label": "Active Modules"},
        {"icon": "\U0001f30a", "value": hostmodel, "label": "Host Model"},
        {"icon": "\U0001f4c8", "value": obs_source, "label": "Observation Source"},
    ]

    return f"""
<div class="section" id="summary">
    <h2>Summary</h2>
    {rh.build_kpi_grid(kpis)}
    {rh.build_highlight_box(
        "<strong>How to use:</strong> Configure calibration settings below, "
        "then scroll to <em>Generated Script</em> at the bottom to download "
        "a ready-to-run Python script.",
        "info"
    )}
</div>
"""


def _build_interactive_settings_section(container_runtime_default: str = "docker") -> str:
    """Calibration settings with interactive form elements."""
    return f"""
<div class="section" id="settings">
    <h2>Calibration Settings</h2>
    <div class="card primary">
        <h3>Optimization</h3>
        <div class="form-row">
            {rh.build_form_select("algorithm", "Algorithm",
                ["DDS", "RANDOM"], "DDS",
                "DDS (Dynamically Dimensioned Search) is recommended")}
            {rh.build_form_number("max_evaluations", "Max Evaluations",
                500, min_val=1, step=1,
                hint="Total number of model evaluations")}
        </div>
        <div class="form-row">
            {rh.build_form_select("objective_function", "Objective Function",
                ["KGE", "NSE", "RMSE", "PBIAS"], "KGE",
                "KGE (Kling-Gupta Efficiency) is recommended")}
            {rh.build_form_number("random_seed", "Random Seed",
                42, min_val=0, step=1)}
        </div>
        <div class="form-row">
            {rh.build_form_select("temporal_resolution", "Temporal Resolution",
                ["native", "daily", "weekly", "monthly"], "monthly",
                "Aggregate obs/model to this scale before computing metrics")}
            {rh.build_form_select("aggregation_method", "Aggregation Method",
                ["mean", "median"], "mean")}
        </div>
        <div class="form-row">
            {rh.build_form_checkbox(
                "use_primary_only",
                "Use only pouring-point (primary) observations in the metric",
                checked=True,
                hint="SUMMA case: only the obs closest to each HRU's outlet "
                     "drives the objective. Secondary (gray) stations still "
                     "appear on the map. For mizuRoute every matched station "
                     "is already primary, so this has no effect on the metric.")}
        </div>
        <div class="form-row">
            {rh.build_form_select(
                "container_runtime", "Container runtime",
                ["docker", "apptainer"], container_runtime_default,
                "How to run the model: 'docker' (local Docker Desktop / "
                "docker compose) or 'apptainer' (Singularity, typical on HPC "
                "clusters). The Apptainer SIF / bind paths below apply only "
                "when 'apptainer' is selected.")}
            {rh.build_form_number(
                "n_parallel", "Parallel model runs", 1, min_val=1, step=1,
                hint="How many model evaluations to run concurrently during "
                     "the sensitivity analysis (the runs are independent). "
                     "Works for both Docker and Apptainer. DDS calibration is "
                     "sequential by design, so this only speeds up the "
                     "sensitivity stage. Set to roughly your available CPU "
                     "cores / container capacity.")}
        </div>
        <div class="form-row">
            {rh.build_form_input(
                "apptainer_sif_path",
                "Apptainer SIF path (HPC)", "",
                hint="Optional. Used by the results report's Apptainer / "
                     "SLURM snippets and forwarded to ModelRunner when "
                     "running on HPC. Leave blank to fall back to 'openwq.sif'.",
                placeholder="e.g. /scratch/$USER/openwq.sif")}
            {rh.build_form_input(
                "apptainer_bind_path",
                "Apptainer bind path (HPC)", "",
                hint="Optional. Bind-mount root containing the config and "
                     "outputs. Leave blank to derive from executable_path.",
                placeholder="e.g. /scratch/$USER/openwq_root")}
        </div>
    </div>
</div>
"""


def _norm_feature_id(val) -> str:
    """Normalise a shapefile feature ID to a clean string.

    Float-typed IDs such as ``740457190.0`` become ``"740457190"`` so they
    match the integer-cast IDs OpenWQ writes to HDF5 / the objective uses.
    """
    try:
        f = float(val)
        if f == int(f):
            return str(int(f))
        return str(f)
    except (TypeError, ValueError):
        return str(val).strip()


def _read_feature_ids(shapefile_path: Optional[str],
                      mapping_key: Optional[str],
                      max_features: int = 5000) -> List[str]:
    """Return sorted unique feature IDs from *shapefile_path*'s *mapping_key*
    column (river reaches for mizuRoute, HRUs for SUMMA).

    Best-effort and dependency-tolerant: returns ``[]`` on any problem
    (missing fiona, missing file, missing column), in which case the UI
    falls back to a free-text reach/HRU input.
    """
    if not shapefile_path or not mapping_key:
        return []
    try:
        import fiona  # noqa: F401
    except Exception:
        return []
    try:
        if not os.path.isfile(shapefile_path):
            return []
        ids = set()
        with fiona.open(shapefile_path) as src:
            for feat in src:
                props = (feat.get("properties", {}) or {})
                if mapping_key in props:
                    val = props[mapping_key]
                else:
                    # case-insensitive fallback
                    lk = {k.lower(): k for k in props}
                    real = lk.get(mapping_key.lower())
                    if real is None:
                        continue
                    val = props[real]
                if val is None:
                    continue
                ids.add(_norm_feature_id(val))
                if len(ids) >= max_features:
                    break
        # numeric sort when every ID is numeric, else lexical
        def _key(x):
            try:
                return (0, float(x))
            except (TypeError, ValueError):
                return (1, str(x))
        return sorted(ids, key=_key)
    except Exception as e:  # pragma: no cover - defensive
        logger.warning(f"Could not read feature IDs from {shapefile_path}: {e}")
        return []


def _build_interactive_targets_section(
    species_list: List[str],
    species_obs_availability: Dict[str, Dict[str, Any]],
    obs_source: str,
    hostmodel: str = "",
    feature_label: str = "Reach",
    available_reaches: Optional[List[str]] = None,
    available_compartments: Optional[List[str]] = None,
    selected_compartments: Optional[List[str]] = None,
) -> str:
    """Calibration targets with species checkboxes and observation info.

    The reach/compartment selectors are host-model-aware:
      * mizuRoute → river reaches + ``RIVER_NETWORK_REACHES``-style compartments
      * SUMMA     → HRUs + the SUMMA land/soil compartments

    *available_reaches* populates a multi-select scroll list (falls back to a
    free-text field if empty); *available_compartments* renders one checkbox
    per compartment, pre-ticked from *selected_compartments*.
    """

    obs_count = sum(
        1 for sp in species_list
        if species_obs_availability.get(sp, {}).get("has_obs", False)
    )
    total = len(species_list)

    no_obs_count = total - obs_count

    if obs_source == "skip":
        obs_note = rh.build_highlight_box(
            "&#9888; <strong>Observations disabled.</strong> Set "
            "<code>observation_data_source</code> in your model config to "
            "<code>\"grqa\"</code> or <code>\"user_csv\"</code> to enable "
            "calibration against observations.",
            "warning"
        )
    elif obs_count == 0:
        obs_note = rh.build_highlight_box(
            "&#9888; <strong>No observation data found</strong> for any of "
            f"the {total} model species. Calibration requires observations "
            "to compare against. Check your "
            "<code>observation_data_source</code> and GRQA data path "
            "settings.",
            "warning"
        )
    elif obs_count < total:
        obs_note = rh.build_highlight_box(
            f"&#9888; <strong>Observation data available for {obs_count} of "
            f"{total} species</strong> (source: {obs_source.upper()}). "
            f"The remaining {no_obs_count} species are greyed out because "
            f"no observation data was found for them &mdash; they cannot be "
            f"included in calibration.",
            "warning"
        )
    else:
        obs_note = rh.build_highlight_box(
            f"&#10003; <strong>Observation data available for all {total} "
            f"species.</strong> Source: {obs_source.upper()}.",
            "success"
        )

    species_html = rh.build_species_checkboxes(
        species_list, species_obs_availability
    )

    # ── Host-aware reach/HRU + compartment selectors ──
    available_reaches = available_reaches or []
    available_compartments = available_compartments or []
    _sel_comp = set(selected_compartments
                    if selected_compartments is not None
                    else available_compartments)
    feat_lower = (feature_label or "Reach").lower()
    feat_plural = "Reaches" if feature_label == "Reach" else f"{feature_label}s"

    if available_reaches:
        _opts = [f'<option value="all" selected>all '
                 f'({len(available_reaches)} {feat_plural})</option>']
        for _rid in available_reaches:
            _r = html_lib.escape(str(_rid))
            _opts.append(f'<option value="{_r}">{_r}</option>')
        reach_field_html = (
            f'<label for="reach_ids">Target {feature_label} IDs</label>'
            f'<select class="form-input" id="reach_ids" multiple size="8" '
            f'style="height:auto;min-height:9rem;font-family:monospace;">'
            f'{"".join(_opts)}</select>'
            f'<div class="hint">Ctrl / Cmd-click to select multiple '
            f'{feat_lower}s; choose <strong>all</strong> to use every '
            f'{feat_lower}.</div>'
        )
    else:
        reach_field_html = (
            f'<label for="reach_ids">Target {feature_label} IDs</label>'
            f'<input class="form-input" type="text" id="reach_ids" '
            f'value="all" placeholder="all or comma-separated IDs"/>'
            f'<div class="hint">"all" or comma-separated IDs '
            f'(no shapefile found to list {feat_lower}s).</div>'
        )

    if available_compartments:
        _checks = []
        for _comp in available_compartments:
            _c = html_lib.escape(_comp)
            _chk = "checked" if _comp in _sel_comp else ""
            _checks.append(
                f'<label class="comp-check" style="display:flex;'
                f'align-items:center;gap:.45rem;padding:.25rem 0;'
                f'cursor:pointer;">'
                f'<input type="checkbox" class="comp-cb" value="{_c}" {_chk}/>'
                f'<code style="font-size:.78rem;">{_c}</code></label>'
            )
        compartment_field_html = (
            f'<label>Compartments</label>'
            f'<div class="compartment-checks" style="border:1px solid '
            f'var(--border);border-radius:8px;padding:.45rem .75rem;'
            f'max-height:9rem;overflow:auto;">{"".join(_checks)}</div>'
            f'<div class="hint">Tick the compartments to include '
            f'({hostmodel or "host"} model).</div>'
        )
    else:
        compartment_field_html = (
            '<label for="compartments">Compartments</label>'
            '<input class="form-input" type="text" id="compartments" '
            'value="RIVER_NETWORK_REACHES" '
            'placeholder="Comma-separated compartment names"/>'
            '<div class="hint">Comma-separated (e.g. RIVER_NETWORK_REACHES)</div>'
        )

    return f"""
<div class="section" id="targets">
    <h2>Calibration Targets</h2>
    {obs_note}
    <div class="card primary">
        <h3>Target Species &amp; Weights</h3>
        <p style="font-size:.85rem;color:var(--text2);margin-bottom:.8rem;">
            Select which species to include in the objective function and set
            their relative weights. Only species with observation data can be
            used for calibration.
        </p>
        {species_html}
    </div>
    <div class="card">
        <h3>{feature_label} &amp; Compartment Selection</h3>
        <div class="form-row">
            <div class="form-group">
                {reach_field_html}
            </div>
            <div class="form-group">
                {compartment_field_html}
            </div>
        </div>
    </div>
</div>
"""


def _build_interactive_parameters_section(parameters: List[Dict]) -> str:
    """Editable parameter table."""
    count = len(parameters)
    if count == 0:
        return f"""
<div class="section" id="parameters">
    <h2>Calibration Parameters</h2>
    {rh.build_highlight_box(
        "<strong>No parameters auto-extracted.</strong> "
        "Your BGC template has no _PARAMETERS_INFO entries.",
        "warning"
    )}
</div>
"""

    return f"""
<div class="section" id="parameters">
    <h2>Calibration Parameters ({count})</h2>
    <p style="color:var(--text2);margin-bottom:1rem;font-size:.9rem;">
        Edit bounds and transform below. Uncheck parameters to exclude them
        from calibration.
    </p>
    {rh.build_editable_param_table(parameters)}
</div>
"""


# Group ordering and display labels
_GROUP_ORDER = [
    "bgc", "transport_dissolved", "sediment_transport",
    "lateral_exchange", "sorption_isotherm", "source_sink",
]
_GROUP_LABELS = {
    "bgc": "BGC Parameters",
    "transport_dissolved": "Transport Dissolved",
    "sediment_transport": "Sediment Transport",
    "lateral_exchange": "Lateral Exchange",
    "sorption_isotherm": "Sorption Isotherm",
    "source_sink": "Source/Sink Loads",
}
_GROUP_MODULE_KEY = {
    "bgc": "bgc_module",
    "transport_dissolved": "td_module",
    "sediment_transport": "ts_module",
    "lateral_exchange": "le_module",
    "sorption_isotherm": "si_module",
    "source_sink": "ss_method",
}


def _build_interactive_parameters_section_grouped(
    module_parameters: Dict[str, List[Dict]],
    module_selections: Dict[str, Any],
    species_obs_availability: Optional[Dict[str, Dict[str, Any]]] = None,
    bgc_network: Optional[dict] = None,
    ss_species_with_loads: Optional[set] = None,
) -> str:
    """Build the parameters section with collapsible module groups.

    *species_obs_availability*: species -> {"has_obs": bool} dict.  When
    provided, parameters that affect only species without observations are
    grayed out.  *bgc_network*: parsed network dict (from
    ``extract_bgc_network``) used to resolve which species each BGC
    reaction touches.  *ss_species_with_loads*: set of model species
    names that have SS load entries (from the generated SS JSON).
    """
    import html as html_lib

    # Pre-compute calibratability per parameter.
    #
    # A reaction is calibratable when tuning its rate can perturb a
    # species that is observed — either directly (its consumed/produced
    # species has obs) OR transitively (the perturbation propagates
    # through one or more downstream/upstream reactions that share a
    # species with the current frontier).
    #
    # Concrete example: a user with NO3 + NH4 observations should be
    # able to calibrate not only the nitrification reaction NH4 -> NO3,
    # but also any upstream reaction A -> NH4 (because the rate of
    # A -> NH4 sets the NH4 pool that drives nitrification, which sets
    # NO3) and any downstream NO3 -> X reaction (because it changes
    # the NO3 concentration the obs are compared against).
    #
    # `compute_obs_reachable` performs the fixed-point species-reaction
    # closure once for the entire network, so per-parameter checks
    # collapse to O(1) set lookups below.
    _reach = rh.compute_obs_reachable(bgc_network, species_obs_availability)
    _active_fws: set = _reach["reachable_frameworks"]
    _active_reactions: set = _reach["reachable_reactions"]   # (fw, rxn_name)

    # Set of model species that have SS loads (already resolved by the
    # config template, which applied stoichiometric conversions).
    _ss_load_species = ss_species_with_loads or set()

    # Observed model species (species with has_obs=True)
    _observed_species: set = set()
    if species_obs_availability:
        _observed_species = {
            sp for sp, info in species_obs_availability.items()
            if info.get("has_obs", False)
        }

    # SS species that can be calibrated = model species that have both
    # loads in the generated SS JSON AND observation data.
    _ss_calibratable_species = _ss_load_species & _observed_species

    def _is_calibratable(p, group_key):
        """Return True if this parameter belongs to an active sub-cycle."""
        if species_obs_availability is None:
            return True  # no obs info → all calibratable
        if not _observed_species:
            return True  # no obs at all → don't gray everything
        if group_key == "bgc":
            # Active = framework has at least one observed species
            # Per-REACTION gating: only enable the parameter if THIS
            # specific reaction is in the obs-reachable closure.
            # `compute_obs_reachable` stores BOTH the (framework, name)
            # tuple AND the (framework, idx) tuple — we try both
            # because `extract_parameters.py` and `extract_bgc_network`
            # disagree on the `_NAME` fallback when the template
            # doesn't explicitly set one (one uses ``rxn{idx}``, the
            # other uses ``LIST_TRANSFORMATIONS[idx]``).  Matching by
            # idx as a backup makes the check robust across templates.
            fw = p.get("_framework", "")
            rxn = p.get("_reaction", "")
            rxn_num = str(p.get("_reaction_num", "") or "")
            if fw and (rxn or rxn_num):
                if (fw, rxn) in _active_reactions:
                    return True
                if (fw, rxn_num) in _active_reactions:
                    return True
                return False
            if fw:
                # Last-resort framework-level fallback for legacy /
                # hand-written parameters that have no reaction tag at
                # all: treat as calibratable iff the framework has any
                # reachable reaction.
                return fw in _active_fws
            return True  # can't determine → keep
        elif group_key == "source_sink":
            # SS params now use the same model species names as the
            # generated SS JSON (e.g. NH4-N, NO3-N).  A param is
            # calibratable if its species is in _ss_calibratable_species
            # (species with both SS loads AND observation data).
            if not _ss_calibratable_species:
                return False
            path_info = p.get("path", {})
            if isinstance(path_info, dict):
                sp = path_info.get("species", "")
                if sp == "all":
                    return bool(_ss_calibratable_species)
                if sp:
                    return sp in _ss_calibratable_species
            # Climate response params (no species) are calibratable
            # if any SS species is calibratable
            return bool(_ss_calibratable_species)
        else:
            # Transport, LE, sediment → affects all species
            return True

    # Build flat list and track group boundaries
    all_params = []
    group_meta = []  # (group_key, label, start_idx, count)

    for group_key in _GROUP_ORDER:
        params = module_parameters.get(group_key, [])
        if not params:
            continue
        # Build label with model name
        label = _GROUP_LABELS.get(group_key, group_key)
        mod_key = _GROUP_MODULE_KEY.get(group_key, "")
        mod_name = module_selections.get(mod_key, "")
        if mod_name:
            label += f" &mdash; {html_lib.escape(str(mod_name))}"

        start_idx = len(all_params)
        # Annotate each param with calibratability
        annotated = []
        for p in params:
            pp = dict(p)  # shallow copy
            pp["_calibratable"] = _is_calibratable(p, group_key)
            annotated.append(pp)
        all_params.extend(annotated)
        group_meta.append((group_key, label, start_idx, len(annotated)))

    total = len(all_params)

    if total == 0:
        return f"""
<div class="section" id="parameters">
    <h2>Calibration Parameters</h2>
    {rh.build_highlight_box(
        "<strong>No calibration parameters found.</strong> "
        "All modules are set to NONE or have no extractable parameters.",
        "warning"
    )}
</div>
"""

    H = [
        f'<div class="section" id="parameters">',
        f'<h2>Calibration Parameters ({total})</h2>',
        '<p style="color:var(--text2);margin-bottom:1rem;font-size:.9rem;">'
        'Edit bounds and transform below. Uncheck parameters to exclude '
        'them from calibration. Parameters are grouped by module.</p>',
    ]

    for group_key, label, start_idx, count in group_meta:
        group_params = all_params[start_idx:start_idx + count]
        has_disabled = any(not p.get("_calibratable", True)
                          for p in group_params)
        H.append(f'<details class="module-group" open>')
        H.append(f'<summary>{label} '
                 f'<span class="group-badge">{count}</span></summary>')
        H.append('<div class="module-content">')
        if has_disabled:
            H.append(rh.build_highlight_box(
                'Grayed-out parameters belong to sub-cycles for which '
                'no observation data is available, so they cannot be '
                'calibrated.',
                'warning'))
        H.append(rh.build_editable_param_table(
            group_params, idx_offset=start_idx))
        H.append('</div></details>')

    H.append('</div>')
    return '\n'.join(H)


def _build_interactive_sensitivity_section() -> str:
    """Sensitivity analysis section with collapsible fields."""
    return f"""
<div class="section" id="sensitivity">
    <h2>Sensitivity Analysis</h2>
    <div class="card secondary">
        {rh.build_form_checkbox("run_sensitivity_first",
            "Run sensitivity analysis before calibration",
            checked=False,
            hint="Identifies influential parameters; less-sensitive ones are fixed")}
        <div class="sa-fields hidden" id="saFields">
            <div class="form-row" style="margin-top:.8rem;">
                {rh.build_form_select("sensitivity_method", "Method",
                    ["morris", "sobol"], "morris")}
                {rh.build_form_number("sensitivity_threshold", "Threshold",
                    0.1, min_val=0, max_val=1, step=0.01,
                    hint="Parameters below this sensitivity index are fixed")}
            </div>
            <div id="morrisFields">
                <div class="form-row">
                    {rh.build_form_number("sensitivity_morris_trajectories",
                        "Morris Trajectories", 10, min_val=2, step=1)}
                    {rh.build_form_number("sensitivity_morris_levels",
                        "Morris Levels", 4, min_val=2, step=1)}
                </div>
            </div>
            <div id="sobolFields" style="display:none;">
                <div class="form-row">
                    {rh.build_form_number("sensitivity_sobol_samples",
                        "Sobol Samples", 1024, min_val=64, step=64)}
                </div>
            </div>
        </div>
    </div>
</div>
"""


def _build_interactive_script_section() -> str:
    """Generated script code block with download button."""
    return """
<div class="section" id="script">
    <h2>Generated Script</h2>
    <div class="highlight-box info">
        <strong>Live preview:</strong> The script below updates as you change
        settings above. Click <em>Download Script</em> when ready.
    </div>
    <details style="margin-bottom:1rem;border:1px solid var(--border);border-radius:10px;padding:.2rem 1rem;">
      <summary style="cursor:pointer;font-weight:600;padding:.6rem 0;color:var(--primary);user-select:none;">
        How to use this script (step-by-step)
      </summary>
      <div style="padding:.4rem 0 .8rem;font-size:.88rem;line-height:1.7;color:var(--text2);">
        <ol style="margin:0;padding-left:1.4rem;">
          <li><strong>Configure settings above</strong> &mdash; Select calibration
              parameters, adjust bounds, choose the algorithm, objective function,
              and sensitivity analysis options. The script preview below updates
              automatically.</li>
          <li><strong>Download the script</strong> &mdash; Click the
              <em>Download Script</em> button. Save the file (e.g.
              <code>my_calibration_run.py</code>) into the
              <code>supporting_scripts/3_Calibration/</code> directory of your
              OpenWQ build, alongside <code>calibration_config_template.py</code>.</li>
          <li><strong>Ensure the Docker container is running</strong> &mdash; Start
              the container before launching calibration:
              <code style="display:block;margin:.3rem 0 .3rem 0;padding:.3rem .6rem;background:var(--dark);color:#e2e8f0;border-radius:6px;font-size:.82rem;">cd containers &amp;&amp; docker compose up -d</code></li>
          <li><strong>Run the script</strong> &mdash; From the
              <code>3_Calibration/</code> directory:
              <code style="display:block;margin:.3rem 0 .3rem 0;padding:.3rem .6rem;background:var(--dark);color:#e2e8f0;border-radius:6px;font-size:.82rem;">python my_calibration_run.py</code></li>
          <li><strong>Monitor progress</strong> &mdash; The calibration driver
              prints objective-function values after each evaluation. Progress is
              saved to checkpoints automatically.</li>
          <li><strong>Resume if interrupted</strong> &mdash; If calibration stops
              (timeout, crash, etc.), resume from the last checkpoint:
              <code style="display:block;margin:.3rem 0 .3rem 0;padding:.3rem .6rem;background:var(--dark);color:#e2e8f0;border-radius:6px;font-size:.82rem;">python my_calibration_run.py --resume</code></li>
          <li><strong>Review results</strong> &mdash; On completion the script
              generates an interactive HTML results report in your
              <code>calibration_work_dir</code> with convergence plots, parameter
              traces, and best-fit values. It opens automatically in your browser.</li>
        </ol>
      </div>
    </details>
    <pre class="code-block" id="scriptBlock" style="max-height:600px;overflow-y:auto;">
<button class="copy-btn">Copy</button>
<code id="scriptCode">Loading...</code>
    </pre>
    <div style="margin-top:1rem;display:flex;align-items:center;gap:1rem;">
""" + rh.build_download_button() + """
        <span style="font-size:.82rem;color:var(--text3);" id="scriptStatus"></span>
    </div>
</div>
"""


# =========================================================================
# Interactive JavaScript
# =========================================================================

def _build_interactive_js(model_config_path, calibration_work_dir,
                          auto_extracted_parameters,
                          module_parameters=None,
                          container_config=None,
                          report_stem="calibration",
                          calib_lib_dir=""):
    """Build the JavaScript for the interactive setup report."""
    import json as json_mod
    # Absolute path to the openWQ "3_Calibration" folder (the parent of
    # calibration_lib).  Baked into the generated run-script's sys.path so
    # the import works regardless of where the user saves the script.
    if not calib_lib_dir:
        calib_lib_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    # Build flat PARAMS array: use module_parameters if available,
    # otherwise fall back to auto_extracted_parameters only
    if module_parameters:
        flat_params = []
        for gk in _GROUP_ORDER:
            flat_params.extend(module_parameters.get(gk, []))
        params_json = rh._js(flat_params)
    else:
        params_json = rh._js(auto_extracted_parameters)

    # Escape paths for JS string embedding
    mcp = json_mod.dumps(model_config_path)
    cwd = json_mod.dumps(calibration_work_dir)
    rstem = json_mod.dumps(report_stem or "calibration")
    clibdir = json_mod.dumps(calib_lib_dir)

    # Build the JS as a plain string (not f-string) to avoid issues
    # with JS curly braces and Python triple-quote conflicts.
    # Use {0}, {1}, {2} placeholders for the 3 injected values.
    js_template = r'''<script>
(function() {
  'use strict';

  // Data injected from Python
  var MODEL_CONFIG_PATH = ''' + mcp + r''';
  var CALIBRATION_WORK_DIR = ''' + cwd + r''';
  var REPORT_STEM = ''' + rstem + r''';
  var CALIB_LIB_DIR = ''' + clibdir + r''';
  var PARAMS = ''' + params_json + r''';

  // Helper: Python repr
  function pyRepr(v, indent) {
    indent = indent || 0;
    var pad = '    '.repeat(indent);
    var pad1 = '    '.repeat(indent + 1);
    if (v === true) return 'True';
    if (v === false) return 'False';
    if (v === null || v === undefined) return 'None';
    if (typeof v === 'number') return String(v);
    if (typeof v === 'string') return JSON.stringify(v);
    if (v && v.__tuple__) {
      return '(' + v.items.map(function(x){ return pyRepr(x); }).join(', ') + ')';
    }
    if (Array.isArray(v)) {
      if (v.length === 0) return '[]';
      if (v.length <= 5 && v.every(function(x){ return typeof x !== 'object' || x === null; })) {
        return '[' + v.map(function(x){ return pyRepr(x); }).join(', ') + ']';
      }
      var lines = v.map(function(x){ return pad1 + pyRepr(x, indent+1) + ','; });
      return '[\n' + lines.join('\n') + '\n' + pad + ']';
    }
    if (typeof v === 'object') {
      var keys = Object.keys(v);
      if (keys.length === 0) return '{}';
      if (keys.length <= 3 && keys.every(function(k){ var val = v[k]; return val === null || typeof val !== 'object'; })) {
        return '{' + keys.map(function(k){ return pyRepr(k) + ': ' + pyRepr(v[k]); }).join(', ') + '}';
      }
      var entries = keys.map(function(k) {
        return pad1 + pyRepr(k) + ': ' + pyRepr(v[k], indent+1) + ',';
      });
      return '{\n' + entries.join('\n') + '\n' + pad + '}';
    }
    return String(v);
  }

  // Collect form state
  function collectFormState() {
    var s = {};
    s.model_config_path = MODEL_CONFIG_PATH;
    s.calibration_work_dir = CALIBRATION_WORK_DIR;

    s.algorithm = document.getElementById('algorithm').value;
    s.max_evaluations = parseInt(document.getElementById('max_evaluations').value) || 500;
    s.objective_function = document.getElementById('objective_function').value;
    s.temporal_resolution = document.getElementById('temporal_resolution').value;
    s.aggregation_method = document.getElementById('aggregation_method').value;
    s.random_seed = parseInt(document.getElementById('random_seed').value) || 42;
    // Container runtime: 'docker' or 'apptainer' (Singularity).
    var _crEl = document.getElementById('container_runtime');
    s.container_runtime = _crEl ? _crEl.value : 'docker';
    // Number of parallel model runs (sensitivity stage).
    var _npEl = document.getElementById('n_parallel');
    s.n_parallel = _npEl ? (parseInt(_npEl.value) || 1) : 1;
    if (s.n_parallel < 1) s.n_parallel = 1;
    // Spatial-matching + HPC options (these are optional — the form
    // fields may not exist on older / customised reports so guard them).
    var _upoEl = document.getElementById('use_primary_only');
    s.use_primary_only = _upoEl ? !!_upoEl.checked : true;
    var _sifEl = document.getElementById('apptainer_sif_path');
    s.apptainer_sif_path = _sifEl ? (_sifEl.value.trim() || null) : null;
    var _bindEl = document.getElementById('apptainer_bind_path');
    s.apptainer_bind_path = _bindEl ? (_bindEl.value.trim() || null) : null;

    var speciesCbs = document.querySelectorAll('.species-cb');
    var species = [];
    var weights = {};
    speciesCbs.forEach(function(cb) {
      if (cb.checked) {
        var sp = cb.dataset.species;
        species.push(sp);
        var wEl = document.getElementById('w_' + sp.replace(/-/g,'_').replace(/ /g,'_'));
        weights[sp] = parseFloat(wEl ? wEl.value : 1) || 1.0;
      }
    });
    s.species = species;
    s.objective_weights = weights;

    // Target reach/HRU IDs: multi-select scroll list when a shapefile was
    // available, else a free-text field.  Normalised to "all" or a
    // comma-separated string so generateScript() stays format-agnostic.
    var reachEl = document.getElementById('reach_ids');
    if (reachEl && reachEl.tagName === 'SELECT') {
      var picked = Array.prototype.slice.call(reachEl.selectedOptions)
        .map(function(o){ return o.value; });
      if (picked.length === 0 || picked.indexOf('all') >= 0) {
        s.reach_ids = 'all';
      } else {
        s.reach_ids = picked.join(',');
      }
    } else {
      s.reach_ids = reachEl ? reachEl.value.trim() : 'all';
    }

    // Compartments: checkbox list when available, else free-text field.
    var compCbs = document.querySelectorAll('.comp-cb');
    if (compCbs.length) {
      var comps = [];
      compCbs.forEach(function(cb){ if (cb.checked) comps.push(cb.value); });
      s.compartments = comps;
    } else {
      var compEl = document.getElementById('compartments');
      var compVal = compEl ? compEl.value.trim() : '';
      s.compartments = compVal.split(',').map(function(c){ return c.trim(); }).filter(Boolean);
    }

    s.run_sensitivity_first = document.getElementById('run_sensitivity_first').checked;
    s.sensitivity_method = document.getElementById('sensitivity_method').value;
    s.sensitivity_morris_trajectories = parseInt(document.getElementById('sensitivity_morris_trajectories').value) || 10;
    s.sensitivity_morris_levels = parseInt(document.getElementById('sensitivity_morris_levels').value) || 4;
    s.sensitivity_sobol_samples = parseInt(document.getElementById('sensitivity_sobol_samples').value) || 1024;
    s.sensitivity_threshold = parseFloat(document.getElementById('sensitivity_threshold').value) || 0.1;

    var paramRows = document.querySelectorAll('.param-row');
    var params = [];
    paramRows.forEach(function(row) {
      var idx = parseInt(row.dataset.paramIdx);
      var cb = row.querySelector('.param-cb');
      if (!cb || !cb.checked) return;
      var orig = PARAMS[idx];
      if (!orig) return;
      var minVal = parseFloat(row.querySelector('.param-min').value);
      var maxVal = parseFloat(row.querySelector('.param-max').value);
      var transform = row.querySelector('.param-transform').value;
      var p = {};
      p.name = orig.name;
      p.file_type = orig.file_type;
      p.path = orig.path;
      p.initial = orig.initial;
      p.bounds = {__tuple__: true, items: [minVal, maxVal]};
      p.transform = transform;
      if (orig.units) p.units = orig.units;
      if (orig.description) p.description = orig.description;
      params.push(p);
    });
    s.parameters = params;

    // Excluded frameworks (set by obs-only toggle)
    s.excluded_frameworks = window._excludedFrameworks || [];

    return s;
  }

  // Generate script text
  function generateScript(s) {
    var lines = [];
    lines.push('#!/usr/bin/env python3');
    lines.push('# Auto-generated calibration script from OpenWQ Interactive Setup Report.');
    lines.push('import sys, os');
    lines.push('# calibration_lib lives in the openWQ "3_Calibration" folder.  Add');
    lines.push('# that folder to sys.path (absolute, baked at generation time) so');
    lines.push('# this script imports calibration_lib no matter where it is saved.');
    lines.push('sys.path.insert(0, ' + pyRepr(CALIB_LIB_DIR) + ')');
    lines.push('sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))');
    lines.push('');
    lines.push('from calibration_lib.calibration_driver import run_calibration');
    lines.push('from calibration_lib import config_integration');
    lines.push('from calibration_lib import Gen_Calibration_Results_Report');
    lines.push('');
    lines.push('');
    lines.push('# ' + '='.repeat(68));
    lines.push('# Configuration (generated from interactive report)');
    lines.push('# ' + '='.repeat(68));
    lines.push('');
    lines.push('model_config_path = ' + pyRepr(s.model_config_path));
    lines.push('calibration_work_dir = ' + pyRepr(s.calibration_work_dir));
    lines.push('report_stem = ' + pyRepr(REPORT_STEM));
    lines.push('');
    lines.push('algorithm = ' + pyRepr(s.algorithm));
    lines.push('max_evaluations = ' + pyRepr(s.max_evaluations));
    lines.push('objective_function = ' + pyRepr(s.objective_function));
    lines.push('temporal_resolution = ' + pyRepr(s.temporal_resolution));
    lines.push('aggregation_method = ' + pyRepr(s.aggregation_method));
    lines.push('random_seed = ' + pyRepr(s.random_seed));
    lines.push('');
    lines.push('# Container runtime selected in the setup report:');
    lines.push('#   "docker"    — local Docker (docker compose up -d first)');
    lines.push('#   "apptainer" — Singularity/Apptainer (set the SIF + bind paths below)');
    lines.push('container_runtime = ' + pyRepr(s.container_runtime || 'docker'));
    lines.push('');
    lines.push('# Number of model evaluations to run concurrently during the');
    lines.push('# sensitivity analysis (independent runs; works for Docker and');
    lines.push('# Apptainer). DDS calibration is sequential, so this only speeds');
    lines.push('# up the sensitivity stage.');
    lines.push('n_parallel = ' + pyRepr(s.n_parallel || 1));
    lines.push('');
    lines.push('# Spatial-matching: when True the objective uses only the');
    lines.push('# pouring-point observation per HRU (SUMMA case). For');
    lines.push('# mizuRoute every matched station is already primary so this');
    lines.push('# has no effect on the metric — but it still drives the');
    lines.push('# primary/secondary visualisation in the results-report map.');
    lines.push('use_primary_only = ' + pyRepr(s.use_primary_only !== false));
    lines.push('');
    lines.push('# HPC defaults (used by the results-report Apptainer / SLURM');
    lines.push('# snippets and forwarded to ModelRunner when running on HPC).');
    lines.push('# Leave as None to fall back to "openwq.sif" + a bind path');
    lines.push('# derived from executable_path.');
    lines.push('apptainer_sif_path = ' + pyRepr(s.apptainer_sif_path || null));
    lines.push('apptainer_bind_path = ' + pyRepr(s.apptainer_bind_path || null));
    lines.push('');

    var reachRepr = s.reach_ids === 'all' ? pyRepr('all') : pyRepr(
      s.reach_ids.split(',').map(function(x){ return parseInt(x.trim()); }).filter(function(x){ return !isNaN(x); })
    );
    lines.push('calibration_targets = {');
    lines.push('    "species": ' + pyRepr(s.species) + ',');
    lines.push('    "reach_ids": ' + reachRepr + ',');
    lines.push('    "compartments": ' + pyRepr(s.compartments) + ',');
    lines.push('}');
    lines.push('');
    lines.push('objective_weights = ' + pyRepr(s.objective_weights));
    lines.push('');

    lines.push('# Sensitivity analysis');
    lines.push('run_sensitivity_first = ' + pyRepr(s.run_sensitivity_first));
    lines.push('sensitivity_method = ' + pyRepr(s.sensitivity_method));
    lines.push('sensitivity_morris_trajectories = ' + pyRepr(s.sensitivity_morris_trajectories));
    lines.push('sensitivity_morris_levels = ' + pyRepr(s.sensitivity_morris_levels));
    lines.push('sensitivity_sobol_samples = ' + pyRepr(s.sensitivity_sobol_samples));
    lines.push('sensitivity_threshold = ' + pyRepr(s.sensitivity_threshold));
    lines.push('');

    // Excluded frameworks (obs-only filter)
    if (s.excluded_frameworks && s.excluded_frameworks.length > 0) {
      lines.push('# Sub-cycles excluded from calibration (no observation data)');
      lines.push('excluded_frameworks = ' + pyRepr(s.excluded_frameworks));
      lines.push('');
    } else {
      lines.push('excluded_frameworks = []');
      lines.push('');
    }

    lines.push('# Calibration parameters (' + s.parameters.length + ' selected)');
    lines.push('calibration_parameters = [');
    s.parameters.forEach(function(p) {
      var entries = [];
      entries.push('"name": ' + pyRepr(p.name));
      entries.push('"file_type": ' + pyRepr(p.file_type));
      entries.push('"path": ' + pyRepr(p.path));
      entries.push('"initial": ' + pyRepr(p.initial));
      entries.push('"bounds": ' + pyRepr(p.bounds));
      entries.push('"transform": ' + pyRepr(p.transform));
      if (p.units) entries.push('"units": ' + pyRepr(p.units));
      if (p.description) entries.push('"description": ' + pyRepr(p.description));
      lines.push('    {' + entries.join(', ') + '},');
    });
    lines.push(']');
    lines.push('');
    lines.push('');

    lines.push('# ' + '='.repeat(68));
    lines.push('# Execution');
    lines.push('# ' + '='.repeat(68));
    lines.push('');
    lines.push('if __name__ == "__main__":');
    lines.push('    import argparse');
    lines.push('    import webbrowser');
    lines.push('');
    lines.push('    parser = argparse.ArgumentParser(description="OpenWQ Calibration")');
    lines.push('    parser.add_argument("--resume", action="store_true",');
    lines.push('                        help="Resume from checkpoint")');
    lines.push('    args = parser.parse_args()');
    lines.push('');
    lines.push('    print("\\n" + "=" * 60)');
    lines.push('    print("OPENWQ CALIBRATION")');
    lines.push('    print("=" * 60)');
    lines.push('');
    lines.push('    print("\\nLoading model config...")');
    lines.push('    model_cfg = config_integration.load_model_config(model_config_path)');
    lines.push('    container_config = config_integration.get_container_config(model_cfg)');
    lines.push('    obs_config = config_integration.get_observation_config(model_cfg)');
    lines.push('');
    lines.push('    print(f"Algorithm: {algorithm}")');
    lines.push('    print(f"Max evaluations: {max_evaluations}")');
    lines.push('    print(f"Objective: {objective_function}")');
    lines.push('    print(f"Parameters: {len(calibration_parameters)}")');
    lines.push('    print()');
    lines.push('');
    lines.push('    results = run_calibration(');
    lines.push('        calibration_work_dir=calibration_work_dir,');
    lines.push('        model_config=model_cfg,');
    lines.push('        calibration_parameters=calibration_parameters,');
    lines.push('        algorithm=algorithm,');
    lines.push('        max_evaluations=max_evaluations,');
    lines.push('        n_parallel=n_parallel,');
    lines.push('        objective_function=objective_function,');
    lines.push('        objective_weights=objective_weights,');
    lines.push('        calibration_targets=calibration_targets,');
    lines.push('        random_seed=random_seed,');
    lines.push('        run_sensitivity_first=run_sensitivity_first,');
    lines.push('        sensitivity_method=sensitivity_method,');
    lines.push('        sensitivity_morris_trajectories=sensitivity_morris_trajectories,');
    lines.push('        sensitivity_morris_levels=sensitivity_morris_levels,');
    lines.push('        sensitivity_sobol_samples=sensitivity_sobol_samples,');
    lines.push('        sensitivity_threshold=sensitivity_threshold,');
    lines.push('        resume=args.resume,');
    lines.push('        temporal_resolution=temporal_resolution,');
    lines.push('        aggregation_method=aggregation_method,');
    // Runtime chosen in the setup report (docker / apptainer), overriding
    // whatever the model config defaulted to.
    lines.push('        container_runtime=container_runtime,');
    lines.push('        docker_container_name=container_config.get("docker_container_name", "docker_openwq"),');
    lines.push('        docker_compose_path=container_config.get("docker_compose_path", ""),');
    lines.push('        executable_full_path=container_config.get("executable_path", ""),');
    lines.push('        file_manager_path=container_config.get("file_manager_path", ""),');
    lines.push('        excluded_frameworks=excluded_frameworks,');
    lines.push('        # Spatial-matching options consumed by ObjectiveFunction.');
    lines.push('        # use_primary_only=True restricts the metric to pouring-');
    lines.push('        # point observations (SUMMA case).  For mizuRoute every');
    lines.push('        # matched station is already primary so it has no effect.');
    lines.push('        use_primary_only=use_primary_only,');
    lines.push('        # HPC paths (forwarded to ModelRunner + surfaced in the');
    lines.push('        # results report\'s Apptainer / SLURM snippets).  Defaults');
    lines.push('        # are None — set them in the setup-report HTML form.');
    lines.push('        apptainer_sif_path=apptainer_sif_path,');
    lines.push('        apptainer_bind_path=apptainer_bind_path,');
    lines.push('    )');
    lines.push('');
    lines.push('    # Generate results report');
    lines.push('    cal_settings = {');
    lines.push('        "algorithm": algorithm, "max_evaluations": max_evaluations,');
    lines.push('        "objective_function": objective_function,');
    lines.push('        "temporal_resolution": temporal_resolution,');
    lines.push('        "aggregation_method": aggregation_method,');
    lines.push('        "calibration_targets": calibration_targets,');
    lines.push('        "objective_weights": objective_weights,');
    lines.push('        "random_seed": random_seed,');
    lines.push('        "run_sensitivity_first": run_sensitivity_first,');
    lines.push('        "sensitivity_method": sensitivity_method,');
    lines.push('        # Echo HPC + primary-only choices so the results report\'s');
    lines.push('        # Apptainer / SLURM snippets show the actual user paths.');
    lines.push('        "use_primary_only": use_primary_only,');
    lines.push('        "apptainer_sif_path": apptainer_sif_path,');
    lines.push('        "apptainer_bind_path": apptainer_bind_path,');
    lines.push('    }');
    lines.push('    report_path = Gen_Calibration_Results_Report.generate_results_report(');
    lines.push('        output_dir=calibration_work_dir,');
    lines.push('        model_config=model_cfg,');
    lines.push('        calibration_parameters=calibration_parameters,');
    lines.push('        calibration_settings=cal_settings,');
    lines.push('        calibration_results=results,');
    lines.push('        sensitivity_results=results.get("sensitivity_results"),');
    lines.push('        performance_metrics=results.get("performance_metrics"),');
    lines.push('        matched_data=results.get("matched_data"),');
    lines.push('        report_stem=report_stem,');
    lines.push('    )');
    lines.push('');
    lines.push('    print("\\n" + "=" * 60)');
    lines.push('    print("CALIBRATION COMPLETE")');
    lines.push('    print("=" * 60)');
    lines.push('    best = results.get("best_objective", "N/A")');
    lines.push('    nevals = results.get("n_evaluations", 0)');
    lines.push('    print(f"Best {objective_function}: {best}")');
    lines.push('    print(f"Total evaluations: {nevals}")');
    lines.push('    if report_path:');
    lines.push('        print(f"Results report: {report_path}")');
    lines.push('        webbrowser.open("file://" + os.path.abspath(report_path))');
    lines.push('');

    return lines.join('\n');
  }

  // ── Syntax highlighting ──
  function highlightPython(code) {
    var escaped = code
      .replace(/&/g, '&amp;')
      .replace(/</g, '&lt;')
      .replace(/>/g, '&gt;');
    // Order matters: comments first, then strings, then keywords/numbers
    // Comments
    escaped = escaped.replace(/(#[^\n]*)/g, '<span class="syn-comment">$1</span>');
    // Single/double quoted strings (simple, non-nested)
    escaped = escaped.replace(/("(?:[^"\\]|\\.)*"|'(?:[^'\\]|\\.)*')/g, function(m) {
      if (m.indexOf('syn-') >= 0) return m;
      return '<span class="syn-str">' + m + '</span>';
    });
    // Constants: True, False, None
    escaped = escaped.replace(/\b(True|False|None)\b/g, '<span class="syn-const">$1</span>');
    // Keywords
    escaped = escaped.replace(/\b(import|from|def|if|elif|else|return|as|for|in|not|and|or|class|try|except|finally|with|raise|pass|break|continue|while|yield|lambda|global|nonlocal|assert|del|is)\b/g, function(m) {
      // Avoid re-highlighting inside already-highlighted spans
      return '<span class="syn-kw">' + m + '</span>';
    });
    // Numbers (integers and floats)
    escaped = escaped.replace(/\b(\d+\.?\d*(?:[eE][+-]?\d+)?)\b/g, function(m, num, offset, str) {
      // Avoid highlighting numbers inside already-highlighted spans
      var before = str.substring(Math.max(0, offset - 20), offset);
      if (before.indexOf('syn-') >= 0 && before.lastIndexOf('>') < before.lastIndexOf('syn-')) return m;
      return '<span class="syn-num">' + num + '</span>';
    });
    // Function calls
    escaped = escaped.replace(/\b([a-zA-Z_]\w*)\s*\(/g, function(m, fn) {
      if (['if','for','while','return','and','or','not','in','is','import','from','def','class','print'].indexOf(fn) >= 0) return m;
      if (m.indexOf('syn-') >= 0) return m;
      return '<span class="syn-fn">' + fn + '</span>(';
    });
    return escaped;
  }

  // ── Tab switching ──
  window.switchTab = function(tabName) {
    document.querySelectorAll('.tab-panel').forEach(function(p) {
      p.classList.toggle('active', p.dataset.tab === tabName);
    });
    document.querySelectorAll('.tab-btn').forEach(function(b) {
      b.classList.toggle('active', b.dataset.tab === tabName);
    });
  };

  // ── Progress indicator ──
  function updateProgress() {
    var state = collectFormState();
    var segments = document.querySelectorAll('.progress-seg');
    var count = 0;

    // Seg 1: At least 1 species selected
    var s1 = state.species && state.species.length > 0;
    if (s1) count++;
    if (segments[0]) segments[0].classList.toggle('active', s1);

    // Seg 2: Algorithm and objective set (always true with defaults)
    var s2 = !!state.algorithm && !!state.objective_function;
    if (s2) count++;
    if (segments[1]) segments[1].classList.toggle('active', s2);

    // Seg 3: At least 1 parameter enabled
    var s3 = state.parameters && state.parameters.length > 0;
    if (s3) count++;
    if (segments[2]) segments[2].classList.toggle('active', s3);

    // Seg 4: Reach IDs specified
    var s4 = !!state.reach_ids && state.reach_ids.trim() !== '';
    if (s4) count++;
    if (segments[3]) segments[3].classList.toggle('active', s4);

    // Seg 5: All bounds valid (min < max for each param)
    var s5 = true;
    if (state.parameters && state.parameters.length > 0) {
      state.parameters.forEach(function(p) {
        if (p.bounds && p.bounds.items) {
          if (p.bounds.items[0] >= p.bounds.items[1]) s5 = false;
        }
      });
    } else {
      s5 = false;
    }
    if (s5) count++;
    if (segments[4]) segments[4].classList.toggle('active', s5);

    var label = document.getElementById('progressLabel');
    if (label) {
      if (count === 5) {
        label.textContent = 'Ready \u2714';
        label.style.opacity = '1';
        label.style.color = '#4ade80';
      } else {
        var remaining = 5 - count;
        label.textContent = remaining + (remaining === 1 ? ' step left' : ' steps left');
        label.style.opacity = '.7';
        label.style.color = '';
      }
    }
  }

  // ── Parameter search ──
  window.filterParams = function(query) {
    var q = query.toLowerCase();
    document.querySelectorAll('.param-row').forEach(function(row) {
      var cells = row.querySelectorAll('td');
      var text = '';
      cells.forEach(function(c) { text += ' ' + c.textContent.toLowerCase(); });
      row.style.display = (!q || text.indexOf(q) >= 0) ? '' : 'none';
    });
    // Also filter module group headers
    document.querySelectorAll('.module-group').forEach(function(group) {
      var visibleRows = group.querySelectorAll('.param-row:not([style*="display: none"])');
      if (visibleRows.length === 0 && q) {
        // Check if group title matches
        var title = group.querySelector('summary');
        if (title && title.textContent.toLowerCase().indexOf(q) >= 0) {
          group.style.display = '';
          group.querySelectorAll('.param-row').forEach(function(r) { r.style.display = ''; });
        } else {
          group.style.display = 'none';
        }
      } else {
        group.style.display = '';
      }
    });
  };

  // ── Copy script ──
  window.copyScript = function() {
    var codeEl = document.getElementById('scriptCode');
    navigator.clipboard.writeText(codeEl.textContent.trim());
    // Target the Copy button specifically (not the Wrap toggle, which
    // also carries the .copy-btn class).
    var btn = document.querySelector('.script-pane-header .copy-btn:not(#wrapToggleBtn)');
    if (btn) {
      btn.textContent = 'Copied!';
      setTimeout(function(){ btn.textContent = 'Copy'; }, 2000);
    }
  };

  // Toggle line wrapping vs single-line + horizontal scroll for the
  // generated-script viewer.  Wrapping is ON by default so the whole
  // script is visible without horizontal scrolling.
  window.toggleWrap = function() {
    var pre = document.getElementById('scriptBlock');
    var btn = document.getElementById('wrapToggleBtn');
    if (!pre) return;
    var nowrap = pre.classList.toggle('nowrap');
    if (btn) btn.textContent = nowrap ? 'Wrap: off' : 'Wrap: on';
  };

  // Reflect the selected container runtime in the "How to use this script"
  // steps: Docker shows the `docker compose up -d` command; Apptainer
  // (Singularity) drops it and points the user at the SIF / bind paths.
  function updateRuntimeUI() {
    var crEl = document.getElementById('container_runtime');
    var rt = crEl ? crEl.value : 'docker';
    var txt = document.getElementById('step3Text');
    var cmd = document.getElementById('step3DockerCmd');
    if (rt === 'apptainer') {
      if (txt) txt.textContent = 'Apptainer / Singularity: build your openwq ' +
        '.sif image and set the SIF & bind paths in the Execution settings ' +
        'tab (no "docker compose" needed).';
      if (cmd) cmd.style.display = 'none';
    } else {
      if (txt) txt.textContent = 'Start the Docker container:';
      if (cmd) cmd.style.display = '';
    }
  }

  // Copy the recommended save path (shown after "Save to:") to the
  // clipboard.  Uses the live text, which the load-time hook refines to
  // wherever the report HTML is actually opened from.
  window.copySaveHint = function() {
    var el = document.getElementById('saveHint');
    var btn = document.getElementById('copyHintBtn');
    if (!el) return;
    navigator.clipboard.writeText(el.textContent.trim());
    if (btn) {
      var prev = btn.textContent;
      btn.textContent = 'Copied!';
      setTimeout(function(){ btn.textContent = prev; }, 1500);
    }
  };


  // Update script display with syntax highlighting
  function updateScript() {
    var state = collectFormState();
    var script = generateScript(state);
    var codeEl = document.getElementById('scriptCode');
    codeEl.innerHTML = highlightPython(script);
    // Update status bar
    var statusEl = document.getElementById('scriptStatus');
    if (statusEl) {
      statusEl.textContent = state.parameters.length + ' params \u2022 ' +
        state.species.length + ' species \u2022 ' +
        state.algorithm + ' \u2022 ' + state.objective_function;
    }
    updateProgress();
  }

  // Save script — uses Save-As dialog when available.
  // The suggested filename follows the originating template stem
  // (<REPORT_STEM>_run.py).  The recommended location (shown in the status
  // hint below the header) is the openWQ "3_Calibration" folder, where
  // calibration_lib lives — though the script bakes that folder into its
  // sys.path so it runs from anywhere.  Browsers cannot pre-set an
  // arbitrary save directory for security reasons, so the hint tells the
  // user exactly where to save it.
  window.downloadScript = function() {
    var state = collectFormState();
    var script = generateScript(state);
    var blob = new Blob([script], {type: 'text/x-python'});
    var btn = document.getElementById('downloadBtnHeader');
    var sugName = (REPORT_STEM ? REPORT_STEM : 'my_calibration') + '_run.py';

    function onSaved() {
      if (btn) {
        btn.textContent = '\u2713 Saved!';
        setTimeout(function(){ btn.textContent = 'Save the script'; }, 2000);
      }
    }

    // Use File System Access API (Save As dialog) when supported
    if (window.showSaveFilePicker) {
      window.showSaveFilePicker({
        suggestedName: sugName,
        types: [{
          description: 'Python Script',
          accept: {'text/x-python': ['.py']},
        }],
      }).then(function(handle) {
        return handle.createWritable().then(function(writable) {
          return writable.write(blob).then(function() {
            return writable.close();
          });
        });
      }).then(onSaved).catch(function(err) {
        // User cancelled the dialog — do nothing
        if (err.name !== 'AbortError') console.error(err);
      });
    } else {
      // Fallback: standard download
      var url = URL.createObjectURL(blob);
      var a = document.createElement('a');
      a.href = url;
      a.download = sugName;
      document.body.appendChild(a);
      a.click();
      document.body.removeChild(a);
      URL.revokeObjectURL(url);
      onSaved();
    }
  };

  // Event binding
  var formInputs = document.querySelectorAll(
    '.form-select, .form-input, .species-cb, .weight-input, ' +
    '.param-cb, .inline-input, .inline-select, .comp-cb'
  );
  formInputs.forEach(function(el) {
    el.addEventListener('change', updateScript);
    el.addEventListener('input', updateScript);
  });

  ['reach_ids', 'compartments'].forEach(function(id) {
    var el = document.getElementById(id);
    if (el) {
      // 'change' covers <select multiple>; 'input' covers text fields.
      el.addEventListener('change', updateScript);
      el.addEventListener('input', updateScript);
    }
  });

  // SA enable checkbox
  var saCheckbox = document.getElementById('run_sensitivity_first');
  var saFields = document.getElementById('saFields');
  if (saCheckbox && saFields) {
    saCheckbox.addEventListener('change', function() {
      saFields.classList.toggle('hidden', !this.checked);
      updateScript();
    });
  }

  // SA method toggle
  var saMethodEl = document.getElementById('sensitivity_method');
  if (saMethodEl) {
    saMethodEl.addEventListener('change', function() {
      var morrisF = document.getElementById('morrisFields');
      var sobolF = document.getElementById('sobolFields');
      if (morrisF) morrisF.style.display = this.value === 'morris' ? '' : 'none';
      if (sobolF) sobolF.style.display = this.value === 'sobol' ? '' : 'none';
      updateScript();
    });
  }

  // Per-group toggle-all checkboxes
  document.querySelectorAll('.groupToggleAll').forEach(function(toggle) {
    toggle.addEventListener('change', function() {
      var table = this.closest('table');
      if (!table) return;
      var cbs = table.querySelectorAll('.param-cb');
      var checked = this.checked;
      cbs.forEach(function(cb) {
        cb.checked = checked;
        var row = cb.closest('tr');
        if (row) row.classList.toggle('param-disabled', !checked);
      });
      updateScript();
    });
  });

  // Individual param checkboxes
  document.querySelectorAll('.param-cb').forEach(function(cb) {
    cb.addEventListener('change', function() {
      var row = this.closest('tr');
      if (row) row.classList.toggle('param-disabled', !this.checked);
      updateScript();
    });
  });

  // Species toggle-all (only enabled species)
  var speciesToggle = document.getElementById('speciesToggleAll');
  if (speciesToggle) {
    speciesToggle.addEventListener('click', function(e) {
      e.preventDefault();
      var cbs = document.querySelectorAll('.species-cb:not(:disabled)');
      var allChecked = Array.from(cbs).every(function(cb){ return cb.checked; });
      cbs.forEach(function(cb){ cb.checked = !allChecked; });
      this.textContent = allChecked ? 'Select all' : 'Deselect all';
      updateScript();
    });
  }

  // ── Pane resizer drag ──
  (function() {
    var resizer = document.getElementById('paneResizer');
    var configPane = document.getElementById('configPane');
    if (!resizer || !configPane) return;

    var dragging = false;

    resizer.addEventListener('mousedown', function(e) {
      e.preventDefault();
      dragging = true;
      resizer.classList.add('dragging');
      document.body.style.cursor = 'col-resize';
      document.body.style.userSelect = 'none';
    });

    document.addEventListener('mousemove', function(e) {
      if (!dragging) return;
      var container = configPane.parentElement;
      var rect = container.getBoundingClientRect();
      var x = e.clientX - rect.left;
      var pct = (x / rect.width) * 100;
      pct = Math.max(20, Math.min(80, pct));
      configPane.style.width = pct + '%';
    });

    document.addEventListener('mouseup', function() {
      if (!dragging) return;
      dragging = false;
      resizer.classList.remove('dragging');
      document.body.style.cursor = '';
      document.body.style.userSelect = '';
    });
  })();

  // ── "Obs data only" diagram toggle → update script only ──
  // The parameter list is independent of the diagram filter.
  // Non-obs sub-cycle params are always visible but disabled/unchecked
  // (set at render time via _is_calibratable). The diagram toggle only
  // affects the visual diagram, not the parameter table.
  document.addEventListener('obs-filter-changed', function(evt) {
    updateScript();
  });

  // Collect excluded frameworks from the static _is_calibratable flags
  // (these are baked in at report generation time, not driven by the
  // diagram toggle).
  window._excludedFrameworks = (function() {
    var excluded = [];
    document.querySelectorAll('.param-row.param-no-obs').forEach(function(row) {
      var fw = row.getAttribute('data-fw');
      if (fw && excluded.indexOf(fw) === -1) excluded.push(fw);
    });
    return excluded;
  })();

  // Initial render
  updateScript();

  // Container-runtime UI: reflect the current choice and update on change.
  updateRuntimeUI();
  var _crEl = document.getElementById('container_runtime');
  if (_crEl) _crEl.addEventListener('change', updateRuntimeUI);

  // Auto-activate "Obs data only" in the diagram on page load
  // (calibration should default to only calibratable sub-cycles)
  setTimeout(function() {
    var obsCb = document.querySelector('.obs-only-cb');
    if (obsCb && !obsCb.checked) {
      obsCb.checked = true;
      obsCb.dispatchEvent(new Event('change', {bubbles: true}));
    }
  }, 100);

})();
</script>
'''
    return js_template
