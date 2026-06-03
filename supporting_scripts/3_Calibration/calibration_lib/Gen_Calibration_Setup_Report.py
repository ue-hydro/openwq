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
    import platform

    # Prefix the long-running commands (calibration, resume, and the
    # influential-parameters / sensitivity run) with the OS-native keep-awake
    # wrapper so the machine doesn't sleep / throttle Docker mid-run.  Fast
    # validation/inspection commands are left un-prefixed.
    #   macOS  : caffeinate -i -s        Linux: systemd-inhibit …
    #   Windows: none (no per-command wrapper) — see report note
    _sys = platform.system()
    if _sys == "Darwin":
        _wake = "caffeinate -i -s "
    elif _sys == "Linux":
        _wake = 'systemd-inhibit --what=idle:sleep --why="OpenWQ calibration run" '
    else:  # Windows
        _wake = ""

    # Determine the calibration command
    run_cmd = f"{_wake}python calibration_config_template.py"
    resume_cmd = f"{_wake}python calibration_config_template.py --resume"
    sa_cmd = f"{_wake}python calibration_config_template.py --sensitivity-only"
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
    hpc_settings_path: Optional[str] = None,
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
    <button class="tab-btn" data-tab="execution" onclick="switchTab('execution')">Execution</button>
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
        # Order follows the workflow decision above: the influential-parameters
        # (sensitivity) settings come first, then the calibration settings.
        H.append('<div class="tab-panel" data-tab="settings">')
        H.append(_build_workflow_mode_section())
        H.append(_build_interactive_sensitivity_section())
        H.append(_build_interactive_settings_section(
            container_config.get("container_runtime", "docker")))
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

        # ── Tab: Execution ──
        # Model-execution backend + HPC/Apptainer settings live here (moved
        # out of Settings).  The HPC card also auto-fills a SLURM sbatch file.
        # HPC fields are pre-filled from hpc_settings.json — the user-provided
        # path (their own copy) if given, else the shipped template.
        _hpc_defaults = _load_hpc_settings(hpc_settings_path)
        H.append('<div class="tab-panel" data-tab="execution">')
        H.append(_build_interactive_execution_section(
            container_config.get("container_runtime", "docker"),
            hpc_defaults=_hpc_defaults))
        H.append('</div>')

        H.append('</div>')  # config-pane

        # ── Resizer Handle ──
        H.append('<div class="pane-resizer" id="paneResizer"></div>')

        # ── Script Pane (Right) ──
        # Where to save the run-script + SLURM job file: the *calibration working
        # directory* (calibration_work_dir) - the host-model-agnostic folder where
        # evaluations + results also live - NOT inside a specific host-model clone
        # (the shared "3_Calibration" folder can resolve into another model's
        # checkout via a symlink, e.g. mizuRoute's, which is wrong for a SUMMA run).
        # The generated script bakes CALIB_LIB_DIR into sys.path, so
        # `from calibration_lib...` works wherever the script is saved.  Falls back
        # to the 3_Calibration folder only when no work dir is given.  Filename
        # follows the originating template stem.
        cal_script_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        _save_dir = (os.path.abspath(calibration_work_dir)
                     if calibration_work_dir else cal_script_dir)
        _run_script_name = f"{_calib_stem}_run.py"
        save_hint = html_lib.escape(os.path.join(_save_dir, _run_script_name))

        def _step_header(n, title, collapsible=False, body_id=None):
            """Numbered step heading for the flat 'Next Steps'-style guide
            (matches the model-configuration report's step presentation).

            When ``collapsible`` is set (with a ``body_id``) the heading becomes
            a clickable toggle and an opening ``<div id=body_id>`` (collapsed by
            default) is returned with it; that div must be closed (``</div>``)
            before the next step header."""
            circle = (
                '<span style="flex-shrink:0;display:inline-flex;width:1.5rem;'
                'height:1.5rem;border-radius:50%;background:var(--primary);'
                'color:#fff;align-items:center;justify-content:center;'
                f'font-size:.8rem;font-weight:700;">{n}</span>')
            label = ('<span style="font-weight:700;color:var(--text);'
                     f'font-size:.95rem;">{title}</span>')
            if collapsible and body_id:
                return (
                    f'<div onclick="toggleStep(\'{body_id}\', this)" '
                    'title="Click to expand / collapse" '
                    'style="display:flex;align-items:center;gap:.55rem;'
                    'margin:1.15rem 1.2rem .25rem;cursor:pointer;user-select:none;">'
                    + circle + label +
                    '<span class="step-caret" style="margin-left:auto;'
                    'margin-right:1.2rem;color:var(--text3);font-size:.7rem;'
                    'transition:transform .15s;">&#9656;</span></div>'
                    f'<div id="{body_id}" class="step-body" style="display:none;">')
            return (
                '<div style="display:flex;align-items:center;gap:.55rem;'
                'margin:1.15rem 1.2rem .25rem;">'
                + circle + label + '</div>')

        H.append('<div class="script-pane">')
        H.append(f"""
<div class="script-pane-header">
    <h3 title="Save the calibration script, then run it locally (Docker) or on HPC (Apptainer). The script updates live as you change the settings on the left.">
        Let's run the calibration&#8230;</h3>
    <div style="display:flex;gap:.4rem;align-items:center;">
        <button class="copy-btn" id="wrapToggleBtn" style="position:static;font-size:.75rem;padding:.25rem .6rem;
            background:rgba(0,102,204,.1);border:1px solid var(--border);color:var(--primary);
            border-radius:6px;cursor:pointer;"
            onclick="toggleWrap()" title="Toggle line wrapping / horizontal scroll">Wrap: on</button>
        <button class="copy-btn" style="position:static;font-size:.75rem;padding:.25rem .6rem;
            background:rgba(0,102,204,.1);border:1px solid var(--border);color:var(--primary);
            border-radius:6px;cursor:pointer;"
            onclick="copyScript()">Copy</button>
    </div>
</div>
<div style="padding:.5rem 1.2rem .1rem;font-size:.78rem;color:var(--text3);line-height:1.5;">
    You've set things up with the tabs on the left. Now follow these steps &mdash; run the
    calibration <strong>either locally (step&nbsp;2) or on HPC (step&nbsp;3)</strong>.
</div>
{_step_header(1, "Save the calibration script")}
<div style="padding:0 1.2rem;">
    <div style="font-size:.72rem;color:var(--text3);word-break:break-all;line-height:1.45;margin-bottom:.45rem;"
         title="{save_hint}">
        Save to: <code id="saveHint" style="font-size:.68rem;">{save_hint}</code>
        <button id="copyHintBtn" onclick="copySaveHint()" title="Copy path to clipboard"
            style="margin-left:.4rem;font-size:.62rem;padding:.05rem .4rem;
            background:rgba(0,102,204,.1);border:1px solid var(--border);
            color:var(--primary);border-radius:4px;cursor:pointer;
            font-family:inherit;white-space:nowrap;vertical-align:baseline;">Copy</button>
    </div>
    <button id="downloadBtnHeader" onclick="downloadScript()"
        style="font-size:.8rem;padding:.4rem .95rem;font-weight:600;
        background:linear-gradient(135deg,rgba(0,102,204,.18),rgba(0,168,107,.18));
        border:1px solid var(--border);color:var(--secondary);border-radius:7px;
        cursor:pointer;">Save the script</button>
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

        # Keep the machine awake for the (long) calibration / sensitivity run
        # so it doesn't sleep or throttle Docker mid-run — important for long
        # overnight runs.  The mechanism is OS-specific; the same prefix covers
        # both the influential-parameters and calibration phases (one run
        # command launches whichever workflow mode is selected).  The fast
        # dry-run is left un-prefixed.
        #   macOS  : caffeinate -i -s CMD       (built-in, releases on exit)
        #   Linux  : systemd-inhibit … CMD      (built-in on systemd distros)
        #   Windows: no per-command wrapper     -> guidance note instead
        if os_name == "Darwin":
            _wake = "caffeinate -i -s "
        elif os_name == "Linux":
            _wake = ('systemd-inhibit --what=idle:sleep '
                     '--why="OpenWQ calibration run" ')
        else:  # Windows
            _wake = ""

        if os_name == "Windows":
            docker_cmd = html_lib.escape(
                f'cd /d "{containers_dir}" && docker compose up -d')
            cd_cmd = html_lib.escape(f'cd /d "{_save_dir}"')
            run_cmd = html_lib.escape(f'python "{_run_script_name}"')
            resume_cmd = html_lib.escape(f'python "{_run_script_name}" --resume')
            dryrun_cmd = html_lib.escape(f'python "{_run_script_name}" --dry-run')
            report_cmd = html_lib.escape(
                f'cd /d "{_save_dir}" && python "{_run_script_name}" --report')
            open_cmd = "start"
        else:
            docker_cmd = html_lib.escape(
                f"cd {containers_dir} && docker compose up -d")
            cd_cmd = html_lib.escape(f"cd {_save_dir}")
            run_cmd = html_lib.escape(f"{_wake}python {_run_script_name}")
            resume_cmd = html_lib.escape(f"{_wake}python {_run_script_name} --resume")
            dryrun_cmd = html_lib.escape(f"python {_run_script_name} --dry-run")
            report_cmd = html_lib.escape(
                f"cd {_save_dir} && python {_run_script_name} --report")
            open_cmd = "open" if os_name == "Darwin" else "xdg-open"

        # OS-specific explanatory note appended to the run / resume steps so the
        # reader understands the keep-awake prefix and how to adapt it.
        _note_css = 'color:var(--text3);font-size:.78rem;'
        if os_name == "Darwin":
            _caf_note = (
                f' <span style="{_note_css}">'
                '(the <code>caffeinate&nbsp;-i&nbsp;-s</code> prefix keeps macOS '
                'awake for the whole run &mdash; covers both the influential-'
                'parameters and calibration phases; releases automatically when '
                'it finishes)</span>')
        elif os_name == "Linux":
            _caf_note = (
                f' <span style="{_note_css}">'
                '(the <code>systemd-inhibit&nbsp;&hellip;</code> prefix stops the '
                'machine idle-sleeping for the whole run; needs systemd &mdash; '
                'drop it on non-systemd hosts or HPC schedulers, which don\'t '
                'idle-sleep)</span>')
        else:  # Windows
            _caf_note = (
                f' <span style="{_note_css}">'
                '(Windows has no per-command keep-awake; before a long run '
                'disable sleep with <code>powercfg /change standby-timeout-ac 0'
                '</code> and re-enable it afterwards, or use Settings &rarr; '
                'System &rarr; Power)</span>')

        # Inline CSS for the mini code snippets with copy buttons.
        # align-items:flex-start + wrapping code so the full command is
        # visible (paths can be long) with the copy button pinned top-right.
        snippet_css = (
            "display:flex;align-items:flex-start;"
            "background:var(--code-bg);color:#e2e8f0;"
            "border:1px solid var(--code-border);border-radius:6px;"
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
{_step_header(2, "Run it locally (Docker)", collapsible=True, body_id="step2body")}
<div style="padding:0 1.2rem;font-size:.82rem;color:var(--text2);line-height:1.6;">
  <p style="margin:.1rem 0 .35rem;"><span id="step3Text">Start the Docker
  container</span>, go to the calibration folder (the one that contains
  <code>calibration_lib</code>), then run the saved script:</p>
  <div id="step3DockerCmd" style="{snippet_css}">
    <code id="cmdDocker" style="{snippet_code_css}">{docker_cmd}</code>
    <button style="{copy_btn_css}"
      onclick="navigator.clipboard.writeText(document.getElementById('cmdDocker').textContent);this.textContent='Copied!';setTimeout(()=>this.textContent='Copy',1500)">Copy</button>
  </div>
  <div style="{snippet_css}">
    <code id="cmdCd" style="{snippet_code_css}">{cd_cmd}</code>
    <button style="{copy_btn_css}"
      onclick="navigator.clipboard.writeText(document.getElementById('cmdCd').textContent);this.textContent='Copied!';setTimeout(()=>this.textContent='Copy',1500)">Copy</button>
  </div>
  <div style="{snippet_css}">
    <code id="cmdRun" style="{snippet_code_css}">{run_cmd}</code>
    <button style="{copy_btn_css}"
      onclick="navigator.clipboard.writeText(document.getElementById('cmdRun').textContent);this.textContent='Copied!';setTimeout(()=>this.textContent='Copy',1500)">Copy</button>
  </div>
  <div style="margin:-.2rem 0 .3rem;">{_caf_note}</div>
  <p style="margin:.5rem 0 .2rem;">Resume if interrupted, or validate the config
  without running (optional dry run):</p>
  <div style="{snippet_css}">
    <code id="cmdResume" style="{snippet_code_css}">{resume_cmd}</code>
    <button style="{copy_btn_css}"
      onclick="navigator.clipboard.writeText(document.getElementById('cmdResume').textContent);this.textContent='Copied!';setTimeout(()=>this.textContent='Copy',1500)">Copy</button>
  </div>
  <div style="{snippet_css}">
    <code id="cmdDryrun" style="{snippet_code_css}">{dryrun_cmd}</code>
    <button style="{copy_btn_css}"
      onclick="navigator.clipboard.writeText(document.getElementById('cmdDryrun').textContent);this.textContent='Copied!';setTimeout(()=>this.textContent='Copy',1500)">Copy</button>
  </div>
</div>
""")

        # ── HPC (Apptainer / Singularity) — copy-to-cluster code snippets ──
        # Replace the old "folders to copy" checklist with ready-to-run shell:
        # an auto-generated SLURM sbatch, an rsync copy that preserves the
        # local layout under ONE HPC base dir, a one-time path-remap (so every
        # baked-in absolute path re-points to its HPC copy), and a submit
        # command.  All the local source paths are baked in; the user only
        # edits HPC_USER / HPC_HOST / HPC_BASE / SIF and the SLURM account.
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
        _mc_dir = (os.path.dirname(os.path.abspath(model_config_path))
                   if model_config_path else "")
        _run_local = os.path.join(_save_dir, _run_script_name)
        # Source trees to copy + the common local root for the path-remap.
        # (rsync -R with the "/./" marker copies each tree RELATIVE to the common
        # root, so re-pointing is just: <common_root>  ->  $HPC_BASE.)
        # Only the calibration_lib PACKAGE is needed at run time (the run-script
        # does `from calibration_lib import ...`) - NOT the whole supporting_scripts
        # tree (other numbered folders, templates, HTML guides...).
        #
        # Anchor it on the ACTIVE PROJECT folder - the common parent of the model
        # config + domain + work dir (e.g. ".../century_basins_tested") - via the
        # project's OWN "3_Calibration" entry (typically a symlink to the openWQ
        # clone).  This way nothing is pulled from an unrelated host-model clone
        # and every copied path stays under the project (=> clean
        # $HPC_BASE/3_Calibration/calibration_lib).  Falls back to calibration_lib's
        # resolved location only when the project has no 3_Calibration of its own.
        _proj_parts = [os.path.abspath(d) for d in (_mc_dir, _domain_dir, _save_dir) if d]
        try:
            _project_root = (os.path.commonpath(_proj_parts) if len(_proj_parts) > 1
                             else (_proj_parts[0] if _proj_parts else cal_script_dir))
        except ValueError:
            _project_root = cal_script_dir
        _proj_calib_lib = os.path.join(_project_root, "3_Calibration", "calibration_lib")
        if os.path.isdir(_proj_calib_lib):
            _calib_lib_src = _proj_calib_lib                       # project path (often a symlink)
            _calib_lib_dir = os.path.join(_project_root, "3_Calibration")
        else:
            _calib_lib_src = os.path.join(cal_script_dir, "calibration_lib")
            _calib_lib_dir = cal_script_dir
        # calibration_lib is NOT copied: on the HPC it is used from the user's
        # already-cloned openWQ repo (the "openWQ clone path on the HPC" field),
        # via the OWQ_CALIB_LIB_DIR env var the sbatch exports.  We only ship the
        # project's model config + domain inputs + the calibration work dir.
        _src_roots = [os.path.abspath(d) for d in
                      (_mc_dir, _domain_dir, _save_dir) if d]
        # de-dup + drop any root that's a subpath of another (rsync -R of the
        # parent already carries it).
        _src_roots = sorted(set(_src_roots))
        _src_roots = [d for d in _src_roots
                      if not any(d != o and d.startswith(o + os.sep)
                                 for o in _src_roots)]
        try:
            _common_root = (os.path.commonpath(
                _src_roots + [os.path.abspath(_run_local),
                              os.path.abspath(model_config_path or _supp_dir),
                              os.path.abspath(calibration_work_dir or _supp_dir)])
                if _src_roots else "")
        except ValueError:
            _common_root = ""
        # If everything to copy lives inside the project folder, the common root
        # collapses to that folder - which would STRIP it from the HPC layout
        # (giving $HPC_BASE/1_Model_Config...).  Anchor one level UP so the project
        # folder is preserved as a clean top level: $HPC_BASE/<project>/...  (not
        # the full local path, not its bare contents).
        if (_common_root
                and os.path.normpath(_common_root) == os.path.normpath(_project_root)):
            _parent = os.path.dirname(_common_root)
            if _parent and _parent not in ("/", ""):
                _common_root = _parent
        _remap_ok = bool(_common_root) and _common_root not in ("/", "")

        # ── 0. SLURM job file ──
        # The sbatch is built in the browser (Execution tab) from the SLURM
        # form fields + this run-script path, and saved there with the
        # "Save .sbatch" button.  We only need its name + local save path
        # here so the copy snippet can rsync it to the cluster.
        _sbatch_name = f"{_calib_stem}.sbatch"
        _sbatch_local = os.path.join(_save_dir, _sbatch_name)

        # Baked local paths handed to the browser so the "copy & run" block in
        # the HPC section can be assembled (pre-filled) from the Execution-tab
        # fields — the user edits nothing, just pastes & runs.
        # Host-model setup dirs in the domain to KEEP, by basename: the active
        # model-run dir (e.g. 0_SUMMA_OPENWQ) + the calibration work dir.  The copy
        # then drops sibling host-model setups (e.g. 0_MIZUROUTE_OPENWQ) that are
        # irrelevant to THIS run.
        _keep_dirs = []
        for _kd in (_model_run_dir,
                    os.path.abspath(calibration_work_dir) if calibration_work_dir else ""):
            _b = os.path.basename(_kd.rstrip(os.sep)) if _kd else ""
            if _b and _b not in _keep_dirs:
                _keep_dirs.append(_b)
        # CONTAINER prefix the model's config files use for the project (e.g.
        # "/code/openwq_code"), parsed from the fileManager/control file.  The HPC
        # Apptainer bind maps $HPC_BASE -> this prefix so every baked container
        # path (SUMMA fileManager, mizuRoute control, ...) resolves to the copied
        # data - identically for any host-model coupling.
        _container_prefix = ""
        _proj_folder = (os.path.basename(_project_root.rstrip(os.sep))
                        if _project_root else "")
        # SUMMA -> fileManager, mizuRoute -> .control; current templates put
        # either under file_manager_path, but accept control_file_path too.
        _fm = ((container_config.get("file_manager_path")
                or container_config.get("control_file_path") or "")
               if container_config else "")
        if _proj_folder and _fm and os.path.isfile(_fm):
            try:
                import re as _re
                _fm_txt = open(_fm, encoding="utf-8", errors="ignore").read()
                _pm = _re.search(r"(/\S*?)/" + _re.escape(_proj_folder) + r"/", _fm_txt)
                if _pm:
                    _container_prefix = _pm.group(1)   # e.g. /code/openwq_code
            except Exception:
                _container_prefix = ""
        _hpc_baked = {
            "src_roots": _src_roots,
            "common_root": _common_root,
            "run_local": _run_local,
            "mc_path": model_config_path or "",
            "sbatch_local": _sbatch_local,
            "sbatch_name": _sbatch_name,
            "work_dir": calibration_work_dir or "",
            "remap_ok": bool(_remap_ok),
            "keep_dirs": _keep_dirs,
            "container_prefix": _container_prefix,
        }

        H.append(f"""
</div>
{_step_header(3, "Run on HPC (Apptainer / Singularity)", collapsible=True, body_id="step3body")}
<div style="padding:0 1.2rem;font-size:.82rem;color:var(--text2);line-height:1.6;">
  <p style="margin:.1rem 0 .5rem;">In the <strong>Execution</strong> tab set
  <strong>Container runtime&nbsp;=&nbsp;apptainer</strong> and fill in the
  <strong>HPC details</strong> + <strong>SLURM</strong> fields. The job file and the
  copy&amp;run commands below fill in automatically from those fields &mdash; no editing
  required.</p>

  <p style="margin:.5rem 0 .15rem;font-weight:600;color:var(--text);">
    Save the SLURM job file (<code>{html_lib.escape(_sbatch_name)}</code>) next to your
    run script:</p>
  <div style="font-size:.72rem;color:var(--text3);word-break:break-all;line-height:1.45;margin:.1rem 0 .5rem;"
       title="{html_lib.escape(_sbatch_local)}">
    Save to: <code id="sbatchSaveHint" style="font-size:.68rem;">{html_lib.escape(_sbatch_local)}</code>
    <button id="copySbatchHintBtn" onclick="copySbatchHint()" title="Copy path to clipboard"
        style="margin-left:.4rem;font-size:.62rem;padding:.05rem .4rem;
        background:rgba(0,102,204,.1);border:1px solid var(--border);
        color:var(--primary);border-radius:4px;cursor:pointer;
        font-family:inherit;white-space:nowrap;vertical-align:baseline;">Copy</button>
  </div>
  <div style="position:relative;margin:.35rem 0 .9rem;">
    <button id="saveSbatchBtn" onclick="downloadSbatch()"
      style="position:absolute;top:.4rem;right:.4rem;background:rgba(0,168,107,.25);
      border:1px solid rgba(0,168,107,.55);color:#d1fae5;padding:.18rem .6rem;
      border-radius:5px;font-size:.66rem;cursor:pointer;font-family:inherit;
      font-weight:600;">Save .sbatch</button>
    <pre id="sbatchPreview" style="background:var(--code-bg,#1e293b);color:#e2e8f0;
      border:1px solid var(--code-border,#334155);
      border-radius:8px;padding:.7rem .9rem;overflow-x:auto;margin:0;white-space:pre;
      font-family:'JetBrains Mono',monospace;font-size:.73rem;line-height:1.55;">Loading&#8230;</pre>
  </div>

  <p style="margin:.5rem 0 .15rem;font-weight:600;color:var(--text);">
    Run these blocks in order. The simplest way is from <strong>your LOCAL terminal</strong>
    &mdash; each one connects to the HPC for you (<code>rsync</code>/<code>ssh</code>), so
    you'll be asked for your HPC password and you don't need to log in first.
    Blocks <strong>b</strong> and <strong>c</strong> <em>also</em> work if you've already
    <code>ssh</code>-ed into the cluster &mdash; they auto-detect that and run there
    directly. (Block <strong>a</strong> copies <em>from</em> your laptop, so it always
    runs locally.)</p>

  <p style="margin:.4rem 0 .1rem;font-size:.78rem;">
    a. Copy code &amp; inputs to the HPC</p>
  <div style="margin:.15rem 0 .4rem;padding:.4rem .6rem;font-size:.74rem;
       background:rgba(217,119,6,.10);border:1px solid rgba(217,119,6,.45);
       border-left:3px solid #d97706;border-radius:6px;color:var(--text);line-height:1.5;">
    <strong style="color:#b45309;">&#9888; Note:</strong> the openWQ
    <code>.sif</code> is <strong>not</strong> copied &mdash; you must
    <strong>build it on the HPC</strong> (the copy block below shows the command).
    It is expected at:
    <code id="sifExpectedPath" style="word-break:break-all;">$HPC_BASE/openwq.sif</code>
  </div>
  <div style="position:relative;margin:.2rem 0 .7rem;">
    <button onclick="copyHpcRun(this,'hpcRunCopy')"
      style="position:absolute;top:.4rem;right:.4rem;background:rgba(255,255,255,.12);
      border:1px solid rgba(255,255,255,.25);color:#e2e8f0;padding:.18rem .6rem;
      border-radius:5px;font-size:.66rem;cursor:pointer;font-family:inherit;">Copy</button>
    <pre id="hpcRunCopy" style="background:var(--code-bg,#1e293b);color:#e2e8f0;
      border:1px solid var(--code-border,#334155);
      border-radius:8px;padding:.7rem .9rem;overflow-x:auto;margin:0;white-space:pre;
      font-family:'JetBrains Mono',monospace;font-size:.73rem;line-height:1.55;">Loading&#8230;</pre>
  </div>

  <p style="margin:.4rem 0 .1rem;font-size:.78rem;">
    b. Re-point the absolute paths to the HPC (run once)</p>
  <div style="position:relative;margin:.2rem 0 .7rem;">
    <button onclick="copyHpcRun(this,'hpcRunRemap')"
      style="position:absolute;top:.4rem;right:.4rem;background:rgba(255,255,255,.12);
      border:1px solid rgba(255,255,255,.25);color:#e2e8f0;padding:.18rem .6rem;
      border-radius:5px;font-size:.66rem;cursor:pointer;font-family:inherit;">Copy</button>
    <pre id="hpcRunRemap" style="background:var(--code-bg,#1e293b);color:#e2e8f0;
      border:1px solid var(--code-border,#334155);
      border-radius:8px;padding:.7rem .9rem;overflow-x:auto;margin:0;white-space:pre;
      font-family:'JetBrains Mono',monospace;font-size:.73rem;line-height:1.55;">Loading&#8230;</pre>
  </div>

  <p style="margin:.4rem 0 .1rem;font-size:.78rem;">
    c. Submit the SLURM job</p>
  <div style="margin:.15rem 0 .4rem;padding:.4rem .6rem;font-size:.74rem;
       background:rgba(37,99,235,.10);border:1px solid rgba(37,99,235,.40);
       border-left:3px solid #2563eb;border-radius:6px;color:var(--text);line-height:1.5;">
    <strong style="color:#1d4ed8;">&#8505; This runs on the cluster.</strong>
    Submitting hands the job to the HPC's <strong>SLURM scheduler</strong>, which only
    exists on the cluster &mdash; so the actual <code>sbatch</code> always executes there.
    The block figures out how to reach it: run it <strong>from your laptop</strong> and it
    <code>ssh</code>-es in to submit for you; run it <strong>after logging into the
    HPC</strong> (<code>ssh&nbsp;$HPC_USER@$HPC_HOST</code>) and it calls <code>sbatch</code>
    directly. Either way works &mdash; nothing to edit.
  </div>
  <div style="position:relative;margin:.2rem 0 .7rem;">
    <button onclick="copyHpcRun(this,'hpcRunSubmit')"
      style="position:absolute;top:.4rem;right:.4rem;background:rgba(255,255,255,.12);
      border:1px solid rgba(255,255,255,.25);color:#e2e8f0;padding:.18rem .6rem;
      border-radius:5px;font-size:.66rem;cursor:pointer;font-family:inherit;">Copy</button>
    <pre id="hpcRunSubmit" style="background:var(--code-bg,#1e293b);color:#e2e8f0;
      border:1px solid var(--code-border,#334155);
      border-radius:8px;padding:.7rem .9rem;overflow-x:auto;margin:0;white-space:pre;
      font-family:'JetBrains Mono',monospace;font-size:.73rem;line-height:1.55;">Loading&#8230;</pre>
  </div>

  <div class="hint" style="margin-top:.3rem;">
    Together these copy your code + inputs + the <code>.sif</code> under one
    <code>$HPC_BASE</code>, re-point every absolute path to the cluster, and submit the
    job. The heavy GRQA / Copernicus processing is <strong>reused</strong> and the
    &gt;1&nbsp;GB raw GRQA database is <strong>not</strong> copied.
  </div>
</div>

</div>
{_step_header(4, "View the results", collapsible=True, body_id="step4body")}
<div style="padding:0 1.2rem .6rem;font-size:.82rem;color:var(--text2);line-height:1.6;">
  <p style="margin:.1rem 0 .35rem;">The script auto-generates the interactive results
  report and opens it when the run finishes. To watch progress <strong>while it's still
  running</strong> (e.g. during the long influential-parameters screening), open a
  <em>second</em> terminal and regenerate the report from the current on-disk state
  &mdash; it rebuilds from the newest checkpoint with an &ldquo;in&nbsp;progress&rdquo;
  note until the run completes:</p>
  <div style="{snippet_css}">
    <code id="cmdReport" style="{snippet_code_css}">{report_cmd}</code>
    <button style="{copy_btn_css}"
      onclick="navigator.clipboard.writeText(document.getElementById('cmdReport').textContent);this.textContent='Copied!';setTimeout(()=>this.textContent='Copy',1500)">Copy</button>
  </div>
  <div style="font-size:.72rem;color:var(--text3);margin:.1rem 0 .3rem;">
    Or just open the last-rendered report file directly (no regeneration):</div>
  <div style="{snippet_css}">
    <code id="cmdResults" style="{snippet_code_css}">{open_cmd} {html_lib.escape(results_path)}</code>
    <button style="{copy_btn_css}"
      onclick="navigator.clipboard.writeText(document.getElementById('cmdResults').textContent);this.textContent='Copied!';setTimeout(()=>this.textContent='Copy',1500)">Copy</button>
  </div>
</div>
</div>
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
        # Observation + model simulation periods drive the calibration /
        # validation split slider.  Either may be None (no obs data / control
        # file unreadable); the slider JS hides itself gracefully.
        try:
            _obs_period = _ci.get_observation_period(
                model_config, work_dir=calibration_work_dir)
        except Exception:
            _obs_period = None
        try:
            _sim_period = _ci.get_model_sim_period(model_config)
        except Exception:
            _sim_period = None

        H.append(_build_interactive_js(
            model_config_path=model_config_path,
            calibration_work_dir=calibration_work_dir,
            auto_extracted_parameters=auto_extracted_parameters,
            module_parameters=module_parameters,
            container_config=container_config,
            calib_lib_dir=_calib_lib_dir,
            report_stem=_calib_stem,
            observation_period=_obs_period,
            model_sim_period=_sim_period,
            hpc_baked=_hpc_baked,
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
        "<strong>How to use:</strong> Configure calibration settings using the "
        "tabs, then use the <em>Let's run the calibration&#8230;</em> panel on the "
        "right to save the script and run it locally or on HPC.",
        "info"
    )}
</div>
"""


def _build_interactive_settings_section(container_runtime_default: str = "docker") -> str:
    """Calibration settings with interactive form elements."""
    return f"""
<div class="section" id="settings">
    <h2>Calibration Settings</h2>
    <!-- Theme: run mode / output (top — biggest lever on per-eval cost) -->
    <div class="card primary">
        <h3>Run mode &amp; output</h3>
        <div class="form-row">
            {rh.build_form_checkbox(
                "run_mode_debug",
                "Enable openWQ debug mode (RUN_MODE_DEBUG) during calibration",
                checked=False,
                hint="Overrides the model config's run_mode_debug. Debug mode "
                     "writes 5 extra diagnostic HDF5 files per "
                     "species/compartment (chemistry, transport, ewf, ic, ss) "
                     "- handy when configuring/testing the model, but a large "
                     "speed / memory / output drag over many evaluations. "
                     "Leave OFF (recommended): each eval then writes only the "
                     "main concentration output. Tick only if you need the "
                     "diagnostics during calibration.")}
        </div>
        <div style="margin-top:.6rem;padding:.55rem .75rem;border-radius:8px;
             background:rgba(245,158,11,.12);border:1px solid rgba(245,158,11,.45);
             border-left:4px solid #f59e0b;font-size:.8rem;color:var(--text);line-height:1.5;">
            <strong style="color:#d97706;">&#9888;&nbsp;WARNING&nbsp;(!)</strong>
            &mdash; enabling debug mode <strong>significantly slows down the
            calibration</strong> (5 extra diagnostic HDF5 files written per
            species/compartment on <em>every</em> evaluation) and is
            <strong>NOT recommended</strong>. Keep it <strong>OFF</strong> unless
            you specifically need the diagnostics.
        </div>
    </div>

    <!-- Theme: optimization algorithm -->
    <div class="card" id="optCard">
        <h3>Optimization algorithm</h3>
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
    </div>

    <!-- Theme: objective metric (how model vs obs is scored) -->
    <div class="card">
        <h3>Objective metric</h3>
        <div class="form-row">
            {rh.build_form_select("temporal_resolution", "Temporal Resolution",
                ["native", "daily", "weekly", "monthly"], "daily",
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
    </div>

    <div class="card" id="calibPeriodCard">
        <h3>Calibration / validation period</h3>
        <p class="hint" style="margin-top:0;">
            Each evaluation simulates <strong>only the calibration window</strong>
            instead of the model's full period &mdash; this keeps runtime and
            memory low (a full multi-year run can exhaust container RAM and be
            OOM-killed). By default the <strong>first third</strong> of your
            observation period is used to <strong>calibrate</strong> and the
            remaining <strong>two thirds</strong> are reserved to
            <strong>validate</strong> later. Drag the <strong>two handles</strong>
            to set the calibration window's <strong>start</strong> (left handle) and
            the calibration&nbsp;/&nbsp;validation <strong>split</strong> (right
            handle); any data before the start handle is <strong>excluded</strong>.
            The greyed span has no observations.
        </p>
        <div class="calib-slider-wrap" id="calibSliderWrap">
            <!-- Per-species observation coverage, ABOVE the slider on the same
                 time axis as the track below (one lane per species with data). -->
            <div id="calibSpeciesLanes" style="margin-bottom:.6rem;"></div>
            <div class="calib-track" id="calibTrack">
                <div class="calib-seg calib-gray"  id="segGrayLeft"></div>
                <div class="calib-seg calib-excl"  id="segExclLeft">
                    <span class="calib-seg-label">Excluded</span></div>
                <div class="calib-seg calib-calib" id="segCalib">
                    <span class="calib-seg-label">Calibration</span></div>
                <div class="calib-seg calib-valid" id="segValid">
                    <span class="calib-seg-label">Validation</span></div>
                <div class="calib-seg calib-gray"  id="segGrayRight"></div>
                <div class="calib-handle" id="calibHandleStart"
                     title="Drag to move the calibration window START"></div>
                <div class="calib-handle" id="calibHandle"
                     title="Drag to move the calibration / validation split"></div>
            </div>
            <div class="calib-axis">
                <span id="calibAxisStart"></span>
                <span id="calibAxisEnd"></span>
            </div>
            <div class="form-row" style="margin-top:.6rem;">
                <div class="form-group">
                    <label>Calibration window <span class="calib-swatch calib-sw-c"></span></label>
                    <div class="calib-window-text" id="calibWindowText">&mdash;</div>
                </div>
                <div class="form-group">
                    <label>Validation window <span class="calib-swatch calib-sw-v"></span></label>
                    <div class="calib-window-text" id="valWindowText">&mdash;</div>
                </div>
            </div>
            <div class="form-row">
                {rh.build_form_checkbox(
                    "use_calibration_period",
                    "Restrict each evaluation to the calibration window "
                    "(recommended)",
                    checked=True,
                    hint="When on, the per-evaluation control file is rewritten "
                         "to run only the calibration window. Turn off to run "
                         "the model's full period per evaluation (slow, "
                         "memory-heavy).")}
            </div>
        </div>
        <p class="hint" id="calibNoDataMsg" style="display:none;">
            &#9888; No observation period or model simulation period could be
            detected, so the period cannot be restricted automatically. Each
            evaluation will run the model's full period. To enable this, make
            sure the model config points to valid observation data and a
            readable control file (fileManager / mizuRoute control).
        </p>
        <!-- Hidden fields populated by the slider JS, read by collectFormState. -->
        <input type="hidden" id="calibration_period_start" value="">
        <input type="hidden" id="calibration_period_end" value="">
        <input type="hidden" id="validation_period_start" value="">
        <input type="hidden" id="validation_period_end" value="">
    </div>
</div>
"""


# hpc_settings.json uses explicit, self-describing key names; map them to the
# (shorter) internal form-field ids the report/JS use. Unknown keys pass
# through unchanged, so older files using the internal names still work.
_HPC_KEY_ALIASES = {
    "hpc_username":              "hpc_user",
    "hpc_login_hostname":        "hpc_host",
    "hpc_working_dir":           "hpc_base",
    "apptainer_sif_path_on_hpc": "sif_hpc",
    "openwq_repo_path_on_hpc":   "hpc_openwq_dir",
    "executable_path_on_hpc":    "hpc_exe",
    "slurm_walltime":            "slurm_time",
    "slurm_cpus_per_task":       "slurm_cpus",
    "slurm_mem_per_node":        "slurm_mem",
    "slurm_modules_to_load":     "slurm_modules",
}


def _load_hpc_settings(path: Optional[str] = None) -> Dict[str, Any]:
    """Load HPC settings JSON so the HPC fields come pre-filled.

    Tries, in order: the explicit ``path`` (a user's own copy, set in the
    calibration template), then the shipped ``hpc_settings.json`` next to the
    run-script templates.  Returns ``{}`` if none readable.  Keys starting with
    ``_`` (e.g. ``_README``) are ignored, and the explicit JSON key names are
    translated to internal field ids via ``_HPC_KEY_ALIASES`` (unknown keys
    pass through unchanged for backward compatibility).
    """
    _cal_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    _candidates = []
    if path:
        _candidates.append(path)
    _candidates.append(os.path.join(_cal_dir, "hpc_settings.json"))
    for _p in _candidates:
        try:
            if _p and os.path.isfile(_p):
                with open(_p, encoding="utf-8") as _f:
                    _d = json.load(_f)
                if isinstance(_d, dict):
                    return {_HPC_KEY_ALIASES.get(k, k): v
                            for k, v in _d.items()
                            if not str(k).startswith("_")}
        except Exception:
            continue
    return {}


def _build_interactive_execution_section(
        container_runtime_default: str = "docker",
        hpc_defaults: Optional[Dict[str, Any]] = None) -> str:
    """Model-execution backend + HPC / Apptainer settings (own tab).

    Moved out of the Settings tab so the execution choices live together.
    The HPC card additionally exposes the SLURM directives that auto-fill a
    SLURM job (sbatch) file: the file is assembled in JavaScript from these
    fields plus the baked run-script path, shown as a live preview, and
    written to disk with a Save button (no Copy button - it's a file).

    The HPC fields are pre-filled from ``hpc_settings.json`` (shipped next to
    the run-script templates) when present; the user can still edit them.
    """
    if hpc_defaults is None:
        hpc_defaults = _load_hpc_settings()

    def _hd(key, fallback):
        """String default for an HPC field: JSON value if set, else fallback."""
        v = hpc_defaults.get(key)
        return fallback if v is None or v == "" else str(v)

    try:
        _slurm_cpus_default = int(hpc_defaults.get("slurm_cpus") or 4)
    except (TypeError, ValueError):
        _slurm_cpus_default = 4

    # Container-runtime selector: two mutually-exclusive buttons (not a dropdown).
    _rt = (container_runtime_default
           if container_runtime_default in ("docker", "apptainer") else "docker")
    _rt_docker_cls = "rt-btn rt-active" if _rt == "docker" else "rt-btn"
    _rt_appt_cls   = "rt-btn rt-active" if _rt == "apptainer" else "rt-btn"
    # Stale-folder handling defaults to the runtime: docker -> prompt (interactive),
    # apptainer -> clean (HPC/batch jobs can't answer an interactive prompt).
    _stale_default = "clean" if _rt == "apptainer" else "prompt"

    return f"""
<div class="section" id="execution">
    <h2>Model execution</h2>

    <!-- Theme: how the models are executed -->
    <div class="card">
        <h3>Execution backend</h3>
        <div class="form-row">
            <div class="form-group">
                <label>Container runtime</label>
                <div id="runtimeToggle" style="display:flex;gap:.4rem;margin-top:.25rem;">
                    <button type="button" class="{_rt_docker_cls}" data-rt="docker"
                        onclick="setRuntime('docker')">&#128051;&nbsp;Docker</button>
                    <button type="button" class="{_rt_appt_cls}" data-rt="apptainer"
                        onclick="setRuntime('apptainer')">&#128230;&nbsp;Apptainer</button>
                </div>
                <input type="hidden" id="container_runtime" value="{_rt}">
                <div class="hint" style="margin-top:.3rem;">How to run the model:
                    <strong>Docker</strong> (local Docker&nbsp;Desktop) or
                    <strong>Apptainer</strong> (Singularity, typical on HPC clusters).
                    The Apptainer SIF&nbsp;/&nbsp;bind paths and the SLURM job below
                    apply only when Apptainer is selected. Selecting a runtime also
                    sets the stale&#8209;folder default below
                    (<em>prompt</em> for Docker, <em>clean</em> for Apptainer).</div>
            </div>
            {rh.build_form_number(
                "n_parallel", "Parallel model runs", 1, min_val=1, step=1,
                hint="How many model evaluations to run concurrently during "
                     "the sensitivity analysis (the runs are independent). "
                     "Works for both Docker and Apptainer. DDS calibration is "
                     "sequential by design, so this only speeds up the "
                     "sensitivity stage. Auto-capped to fit container memory.")}
        </div>
        <div class="form-row">
            {rh.build_form_select(
                "clean_work_dir", "Stale evaluation folders",
                ["prompt", "clean", "keep"], _stale_default,
                "What to do with leftover eval_* folders from a previous run "
                "in the calibration work dir. 'prompt' asks at the terminal "
                "(interactive only). 'clean' always deletes them first - use "
                "this on HPC / batch jobs where you can't answer a prompt. "
                "'keep' leaves them in place. (Resume runs never clean.)")}
        </div>
    </div>

    <!-- Theme: HPC / Apptainer (optional) -->
    <div class="card" id="hpcCard">
        <h3>HPC/Apptainer Settings</h3>
        <p class="hint" style="margin-top:0;">
            Fill these in to run the calibration on an HPC cluster with
            Apptainer / Singularity. You won't edit any code: the values below
            are dropped straight into the ready-to-run <strong>SLURM job file</strong>
            and <strong>copy&nbsp;&amp;&nbsp;run</strong> commands shown under
            <em>step&nbsp;3 (Run it on HPC)</em> in the script pane &mdash; just
            save the <code>.sbatch</code> and paste the commands into a terminal.
        </p>
        <p class="hint" style="margin-top:.45rem;background:rgba(217,119,6,.07);
            border:1px solid var(--border);border-left:3px solid #d97706;
            border-radius:6px;padding:.5rem .65rem;">
            <strong>&#9881;&nbsp;Build the <code>.sif</code> on the HPC:</strong> the openWQ
            Apptainer image must be built <strong>on the cluster itself</strong> (images are
            architecture-specific, and Deucalion-class systems have ARM and x86 nodes) &mdash; it
            is <strong>not</strong> uploaded. Build it once with
            <code>apptainer build openwq.sif &lt;openwq&gt;/containers/openwq_apptainer.def</code>,
            then put its absolute cluster path in <em>Apptainer image (.sif) on the HPC</em> below.
        </p>
        <p class="hint" style="margin-top:.45rem;background:rgba(0,102,204,.06);
            border:1px solid var(--border);border-left:3px solid var(--primary);
            border-radius:6px;padding:.5rem .65rem;">
            <strong>&#128190;&nbsp;Reuse across runs (recommended):</strong> these fields are
            pre-filled from a <code>hpc_settings.json</code> file. To avoid re-typing your HPC
            details every time, <strong>make a copy</strong> of <code>hpc_settings.json</code>,
            <strong>complete it</strong> with your values, and <strong>point your calibration
            template at it</strong> &mdash; set
            <code>hpc_settings_json&nbsp;=&nbsp;"/path/to/your_copy.json"</code> near the top of
            <code>calibration_config_template.py</code>. Every regenerated report then opens with
            these fields already filled in (you can still tweak them here).
        </p>

        <h4 style="margin:.6rem 0 .3rem;font-size:.92rem;color:var(--text);">
            Your HPC connection</h4>
        <div class="form-row">
            {rh.build_form_input(
                "hpc_user", "HPC username", _hd("hpc_user", ""),
                hint="Your login name on the cluster (the part before @).",
                placeholder="your_username")}
            {rh.build_form_input(
                "hpc_host", "HPC hostname", _hd("hpc_host", ""),
                hint="The cluster you ssh into.",
                placeholder="hpc.your-institution.edu")}
        </div>
        <div class="form-row">
            {rh.build_form_input(
                "hpc_base", "HPC working directory",
                _hd("hpc_base", "/scratch/$USER/openwq_cal"),
                hint="A writable dir on the cluster. Everything is copied "
                     "under here and all paths re-point to it automatically.",
                placeholder="/scratch/$USER/openwq_cal")}
            {rh.build_form_input(
                "sif_hpc", "Apptainer image (.sif) on the HPC", _hd("sif_hpc", ""),
                hint="Absolute path to the openwq .sif ON the cluster. Build the "
                     ".sif on the HPC itself (Apptainer images must be built where "
                     "they run, so the CPU architecture matches the compute nodes) "
                     "- it is NOT uploaded. The model is pointed at this path "
                     "automatically.",
                placeholder="$HPC_BASE/openwq.sif")}
        </div>
        <div class="form-row">
            {rh.build_form_input(
                "hpc_openwq_dir", "openWQ clone path on the HPC",
                _hd("hpc_openwq_dir", ""),
                hint="Absolute path to the openWQ checkout INSIDE your host-model "
                     "clone on the cluster - the folder that contains "
                     "supporting_scripts/ (the inner '.../openwq/openwq').  Point at "
                     "that folder (or its parent); the sbatch auto-locates "
                     "supporting_scripts/3_Calibration/calibration_lib under it.  The "
                     "code is used FROM there (not copied) - exported as "
                     "OWQ_CALIB_LIB_DIR.",
                placeholder="/home/$USER/.../openwq/openwq")}
        </div>
        <div class="form-row">
            {rh.build_form_input(
                "hpc_exe", "Model executable on the HPC",
                _hd("hpc_exe", ""),
                hint="Absolute path to the host-model binary you COMPILED on the "
                     "cluster (summa / mizuroute / mizuroute_cslm - whichever you "
                     "built).  The .sif is dependency-only, so the binary is run from "
                     "here.  This is the only per-coupling path the tool can't guess; "
                     "everything else (inputs, .sif, calibration_lib) is automatic.  "
                     "Exported as OWQ_EXEC_PATH.",
                placeholder="$HPC_BASE/.../bin/summa_openwq_Release")}
        </div>

        <h4 style="margin:1rem 0 .3rem;font-size:.92rem;color:var(--text);">
            SLURM job file (sbatch)</h4>
        <p class="hint" style="margin-top:0;">
            These populate the <code>#SBATCH</code> directives in the job file.
        </p>
        <div class="form-row">
            {rh.build_form_input(
                "slurm_job_name", "Job name", _hd("slurm_job_name", "openwq_calib"),
                hint="SLURM --job-name (also names the .out log).")}
            {rh.build_form_input(
                "slurm_account", "Account / allocation", _hd("slurm_account", ""),
                hint="SLURM --account. Your compute allocation (required on "
                     "most clusters, e.g. Compute Canada / Slurm fairshare).",
                placeholder="e.g. def-yourpi")}
        </div>
        <div class="form-row">
            {rh.build_form_input(
                "slurm_partition", "Partition / queue", _hd("slurm_partition", ""),
                hint="SLURM --partition. Leave blank to use the cluster "
                     "default queue.",
                placeholder="e.g. batch / cpu2023")}
            {rh.build_form_input(
                "slurm_time", "Wall time (HH:MM:SS)", _hd("slurm_time", "48:00:00"),
                hint="SLURM --time. Max run time before the job is killed.",
                placeholder="48:00:00")}
        </div>
        <div class="form-row">
            {rh.build_form_number(
                "slurm_cpus", "CPUs per task", _slurm_cpus_default,
                min_val=1, step=1,
                hint="SLURM --cpus-per-task. Should be >= parallel model "
                     "runs.")}
            {rh.build_form_input(
                "slurm_mem", "Memory", _hd("slurm_mem", "16G"),
                hint="SLURM --mem (per node). Include the unit, e.g. 16G.",
                placeholder="16G")}
        </div>
        <div class="form-row">
            {rh.build_form_input(
                "slurm_modules",
                "Module load / env activation (optional)", _hd("slurm_modules", ""),
                hint="Command(s) run on the compute node to make Python + "
                     "Apptainer available. Leave blank to insert a commented "
                     "placeholder you can edit later.",
                placeholder="module load apptainer python")}
        </div>
        <p class="hint" style="margin:.6rem 0 0;">
            &#10142;&nbsp;The filled-in <code>.sbatch</code> file and the
            copy&nbsp;&amp;&nbsp;run terminal block are in the
            <em>Run the calibration on HPC</em> section of the script pane.
        </p>
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
        _n = len(available_reaches)
        # Count-aware word so a single feature reads "1 HRU", not "1 HRUs".
        _word = feature_label if _n == 1 else feat_plural
        # Size the list to its contents (+1 for the "all" row), clamped so it's
        # neither a cramped 1-liner nor a tall box full of empty rows.
        _size = max(3, min(10, _n + 1))
        _opts = [f'<option value="all" selected '
                 f'style="font-weight:700;color:var(--primary);">'
                 f'&#10003; all ({_n} {_word})</option>']
        for _rid in available_reaches:
            _r = html_lib.escape(str(_rid))
            _opts.append(f'<option value="{_r}">{_r}</option>')
        reach_field_html = (
            f'<label for="reach_ids">Target {feature_label} IDs'
            f'<span style="font-weight:500;font-size:.68rem;color:var(--text3);'
            f'margin-left:.45rem;background:var(--bg);border:1px solid var(--border);'
            f'border-radius:10px;padding:.05rem .45rem;vertical-align:middle;">'
            f'{_n} available</span></label>'
            f'<select class="form-input" id="reach_ids" multiple size="{_size}" '
            f'style="height:auto;font-family:\'JetBrains Mono\',monospace;'
            f'padding:.3rem;line-height:1.65;">'
            f'{"".join(_opts)}</select>'
            f'<div class="hint">The highlighted <strong>all</strong> option '
            f'uses every {feat_lower}. To target specific {feat_plural.lower()} '
            f'instead, Ctrl / Cmd-click to pick one or more.</div>'
        )
    else:
        reach_field_html = (
            f'<label for="reach_ids">Target {feature_label} IDs</label>'
            f'<input class="form-input" type="text" id="reach_ids" '
            f'value="all" placeholder="all or comma-separated IDs"/>'
            f'<div class="hint">"all" or comma-separated IDs '
            f'(no shapefile found to list {feat_plural.lower()}).</div>'
        )

    if available_compartments:
        # Single-select: exactly ONE compartment can be targeted. Pre-select
        # the first previously-selected compartment, else the first available.
        _default_comp = next((c for c in available_compartments
                              if c in _sel_comp), available_compartments[0])
        _checks = []
        for _comp in available_compartments:
            _c = html_lib.escape(_comp)
            _chk = "checked" if _comp == _default_comp else ""
            _checks.append(
                f'<label class="comp-check" style="display:flex;'
                f'align-items:center;gap:.45rem;padding:.25rem 0;'
                f'cursor:pointer;">'
                f'<input type="radio" name="compartment" class="comp-cb" '
                f'value="{_c}" {_chk}/>'
                f'<code style="font-size:.78rem;">{_c}</code></label>'
            )
        compartment_field_html = (
            f'<label>Compartment</label>'
            f'<div class="compartment-checks" style="border:1px solid '
            f'var(--border);border-radius:8px;padding:.45rem .75rem;'
            f'max-height:9rem;overflow:auto;">{"".join(_checks)}</div>'
            f'<div class="hint">Select <strong>one</strong> compartment to '
            f'target ({hostmodel or "host"} model).</div>'
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

    # Map each reaction to the species it consumes/produces so we can show
    # "A -> B" under every BGC parameter.  Keyed by BOTH (framework, 1-based
    # index) and (framework, reaction name) — mirroring compute_obs_reachable
    # — so it matches whichever tag a parameter carries (_reaction_num /
    # _reaction).  Fully data-driven: works for any NATIVE_BGC_FLEX template.
    _rxn_io_lookup: Dict[tuple, Dict[str, list]] = {}
    if bgc_network:
        for _fw_name, _fw_info in (bgc_network.get("frameworks") or {}).items():
            for _i, _r in enumerate(_fw_info.get("reactions") or []):
                _io = {
                    "consumed": list(rh._yield_species(_r.get("consumed"))),
                    "produced": list(rh._yield_species(_r.get("produced"))),
                }
                _rxn_io_lookup[(_fw_name, str(_i + 1))] = _io
                _nm = _r.get("name") or ""
                if _nm:
                    _rxn_io_lookup[(_fw_name, _nm)] = _io

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
            # Attach the reaction's consumed/produced species (BGC params only)
            # so the table can render "A -> B" under the parameter name.
            if group_key == "bgc":
                _fw = pp.get("_framework", "")
                _rnum = str(pp.get("_reaction_num", "") or "")
                _rname = str(pp.get("_reaction", "") or "")
                _io = (_rxn_io_lookup.get((_fw, _rnum))
                       or _rxn_io_lookup.get((_fw, _rname)))
                if _io and (_io["consumed"] or _io["produced"]):
                    pp["reaction_io"] = _io
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


def _build_workflow_mode_section() -> str:
    """Top-of-settings 3-way workflow selector rendered as a segmented control
    (three distinct, clickable options — not a continuous slider):
    Influential parameters (Morris) · Both · Calibration."""
    return """
<div class="section" id="workflow">
    <h2>Workflow</h2>
    <div class="card primary">
        <h3>What to run</h3>
        <p class="hint" style="margin-top:0;">
            Choose <strong>one</strong> of the three options below. The selected
            one is highlighted, and the settings panels update to match it.
        </p>
        <div class="mode-slider-wrap">
            <div class="mode-seg" role="radiogroup" aria-label="Workflow to run">
                <span class="mode-seg-thumb" id="modeThumb" aria-hidden="true"></span>
                <button type="button" class="mode-seg-btn" data-m="0"
                        role="radio" aria-checked="false">
                    <span class="mode-seg-title">Influential parameters</span>
                    <span class="mode-seg-sub">Screening only &mdash; no calibration</span>
                </button>
                <button type="button" class="mode-seg-btn" data-m="1"
                        role="radio" aria-checked="true">
                    <span class="mode-seg-title">Both</span>
                    <span class="mode-seg-sub">Screening first, then calibrate</span>
                </button>
                <button type="button" class="mode-seg-btn" data-m="2"
                        role="radio" aria-checked="false">
                    <span class="mode-seg-title">Calibration</span>
                    <span class="mode-seg-sub">Calibration only</span>
                </button>
            </div>
            <input type="hidden" id="calibration_mode" value="1"/>
            <div class="mode-desc" id="modeDesc"></div>
        </div>
    </div>
</div>
"""


def _build_interactive_sensitivity_section() -> str:
    """Influential-parameters (sensitivity) settings.  Hidden by the workflow
    slider when the user selects 'Calibration' only."""
    return f"""
<div class="section" id="sensitivity">
    <h2>Influential parameters (sensitivity)</h2>
    <div class="card secondary">
        <p class="hint" style="margin-top:0;">
            Morris (or Sobol) screening ranks parameters by how much they move
            the objective. In <strong>Both</strong> mode the least-sensitive
            ones can be fixed during calibration.
        </p>
        <div class="sa-fields" id="saFields">
            <div class="form-row">
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
                          calib_lib_dir="",
                          observation_period=None,
                          model_sim_period=None,
                          hpc_baked=None):
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

    # Observation + model simulation periods → power the calibration /
    # validation split slider in the Settings tab.  Computed by the caller
    # (where model_config is available) and passed in; either may be None
    # (no obs data, or control file unreadable) and the slider JS copes.
    _obs_period = observation_period
    _sim_period = model_sim_period

    # Escape paths for JS string embedding
    mcp = json_mod.dumps(model_config_path)
    cwd = json_mod.dumps(calibration_work_dir)
    rstem = json_mod.dumps(report_stem or "calibration")
    clibdir = json_mod.dumps(calib_lib_dir)
    hpc_json = json_mod.dumps(hpc_baked or {})
    obs_period_json = json_mod.dumps(_obs_period)
    sim_period_json = json_mod.dumps(_sim_period)

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
  var HPC = ''' + hpc_json + r''';
  var PARAMS = ''' + params_json + r''';
  var OBS_PERIOD = ''' + obs_period_json + r''';
  var SIM_PERIOD = ''' + sim_period_json + r''';

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
    // Stale-folder handling: 'prompt' | 'clean' | 'keep'.
    var _cwdEl = document.getElementById('clean_work_dir');
    s.clean_work_dir_mode = _cwdEl ? _cwdEl.value : 'prompt';
    // openWQ debug-mode override (default OFF — diagnostics are a drag here).
    var _rmdEl = document.getElementById('run_mode_debug');
    s.run_mode_debug = _rmdEl ? !!_rmdEl.checked : false;
    // Calibration / validation period split (from the slider).  The hidden
    // inputs are populated by the slider JS in 'YYYY-MM-DD HH:MM' format.
    var _ucpEl = document.getElementById('use_calibration_period');
    s.use_calibration_period = _ucpEl ? !!_ucpEl.checked : false;
    var _cps = document.getElementById('calibration_period_start');
    var _cpe = document.getElementById('calibration_period_end');
    var _vps = document.getElementById('validation_period_start');
    var _vpe = document.getElementById('validation_period_end');
    s.calibration_period_start = (_cps && _cps.value) ? _cps.value : null;
    s.calibration_period_end   = (_cpe && _cpe.value) ? _cpe.value : null;
    s.validation_period_start  = (_vps && _vps.value) ? _vps.value : null;
    s.validation_period_end    = (_vpe && _vpe.value) ? _vpe.value : null;
    // Spatial-matching + HPC options (these are optional — the form
    // fields may not exist on older / customised reports so guard them).
    var _upoEl = document.getElementById('use_primary_only');
    s.use_primary_only = _upoEl ? !!_upoEl.checked : true;
    var _sifEl = document.getElementById('apptainer_sif_path');
    s.apptainer_sif_path = _sifEl ? (_sifEl.value.trim() || null) : null;
    var _bindEl = document.getElementById('apptainer_bind_path');
    s.apptainer_bind_path = _bindEl ? (_bindEl.value.trim() || null) : null;
    // SLURM directives -> auto-filled sbatch (Execution tab). All optional;
    // guarded so older/customised reports without these fields still work.
    function _gv(id, dflt) {
      var el = document.getElementById(id);
      return (el && el.value != null && String(el.value).trim() !== '')
             ? String(el.value).trim() : dflt;
    }
    s.slurm_job_name  = _gv('slurm_job_name', 'openwq_calib');
    s.slurm_account   = _gv('slurm_account', '');
    s.slurm_partition = _gv('slurm_partition', '');
    s.slurm_time      = _gv('slurm_time', '48:00:00');
    s.slurm_cpus      = parseInt(_gv('slurm_cpus', '4')) || 4;
    s.slurm_mem       = _gv('slurm_mem', '16G');
    s.slurm_modules   = _gv('slurm_modules', '');
    // HPC connection (drives the pre-filled copy & run block).
    s.hpc_user  = _gv('hpc_user', 'your_username');
    s.hpc_host  = _gv('hpc_host', 'hpc.your-institution.edu');
    s.hpc_base  = _gv('hpc_base', '/scratch/$USER/openwq_cal');
    s.sif_hpc = _gv('sif_hpc', '$HPC_BASE/openwq.sif');
    s.hpc_openwq_dir = _gv('hpc_openwq_dir', '');
    s.hpc_exe = _gv('hpc_exe', '');

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

    // Workflow mode slider: 0=influential params only, 1=both, 2=calibration only.
    var _modeEl = document.getElementById('calibration_mode');
    var _modeVal = _modeEl ? (parseInt(_modeEl.value)) : 1;
    s.calibration_mode = (_modeVal === 0) ? 'sensitivity'
                       : (_modeVal === 2) ? 'calibration' : 'both';
    // Legacy flag derived from the mode (sensitivity runs unless 'calibration').
    s.run_sensitivity_first = (s.calibration_mode !== 'calibration');
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
    lines.push('# calibration_lib location: on HPC the sbatch exports');
    lines.push('# OWQ_CALIB_LIB_DIR (the cloned openWQ 3_Calibration); otherwise fall');
    lines.push('# back to the path baked here at generation time (local runs).');
    lines.push('_CALIB_LIB_DIR = os.environ.get("OWQ_CALIB_LIB_DIR") or ' + pyRepr(CALIB_LIB_DIR));
    lines.push('sys.path.insert(0, _CALIB_LIB_DIR)');
    lines.push('sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))');
    lines.push('');
    lines.push('from calibration_lib.calibration_driver import run_calibration, regenerate_results_report');
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
    lines.push('# Handling of leftover eval_* folders from a previous run:');
    lines.push('#   None  -> prompt at the terminal (interactive only)');
    lines.push('#   True  -> always delete first (use on HPC / batch jobs)');
    lines.push('#   False -> keep them');
    lines.push('# Overridable at runtime with --clean / --keep.');
    var _cwdMode = s.clean_work_dir_mode || 'prompt';
    var _cwdRepr = (_cwdMode === 'clean') ? 'True'
                 : (_cwdMode === 'keep') ? 'False' : 'None';
    lines.push('clean_work_dir = ' + _cwdRepr);
    lines.push('');
    lines.push('# openWQ debug-mode override for calibration. The model config may');
    lines.push('# set run_mode_debug=True (useful when configuring), but it makes');
    lines.push('# openWQ write 5 extra diagnostic HDF5 files per species/compartment');
    lines.push('# - a big speed/memory/output drag. Forced OFF here unless ticked.');
    lines.push('run_mode_debug = ' + pyRepr(!!s.run_mode_debug));
    lines.push('');
    lines.push('# Calibration window (start, end) that each evaluation simulates.');
    lines.push('# Restricting the per-eval period keeps runtime + memory low and');
    lines.push('# avoids out-of-memory kills from full multi-year runs (openWQ');
    lines.push('# holds its output in RAM). The complementary validation window');
    lines.push('# is recorded for a later validation run.');
    lines.push('#   (start, end) -> rewrite simStart/simEnd in each eval control file');
    lines.push('#   None         -> run the model\'s full period per evaluation');
    if (s.use_calibration_period && s.calibration_period_start && s.calibration_period_end) {
      lines.push('calibration_period = (' + pyRepr(s.calibration_period_start) +
                 ', ' + pyRepr(s.calibration_period_end) + ')');
      lines.push('validation_period = (' + pyRepr(s.validation_period_start) +
                 ', ' + pyRepr(s.validation_period_end) + ')');
    } else {
      lines.push('calibration_period = None');
      lines.push('validation_period = None');
    }
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

    lines.push('# Workflow mode (from the report slider):');
    lines.push('#   "sensitivity" -> identify influential parameters only (no calibration)');
    lines.push('#   "both"        -> sensitivity first, then calibration');
    lines.push('#   "calibration" -> calibration only (no sensitivity)');
    lines.push('calibration_mode = ' + pyRepr(s.calibration_mode || 'both'));
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
    lines.push('    parser.add_argument("--clean", action="store_true",');
    lines.push('                        help="Delete stale eval_* folders before running "');
    lines.push('                             "(no prompt; ideal for HPC/batch jobs)")');
    lines.push('    parser.add_argument("--keep", action="store_true",');
    lines.push('                        help="Keep stale eval_* folders (no prompt)")');
    lines.push('    parser.add_argument("--report", action="store_true",');
    lines.push('                        help="Regenerate the results report from the latest "');
    lines.push('                             "on-disk state and open it, then exit. Works "');
    lines.push('                             "while a run is still going (run in a 2nd terminal).")');
    lines.push('    args = parser.parse_args()');
    lines.push('    # CLI flags override the report default (clean_work_dir).');
    lines.push('    if args.clean:');
    lines.push('        clean_work_dir = True');
    lines.push('    elif args.keep:');
    lines.push('        clean_work_dir = False');
    lines.push('');
    lines.push('    print("\\n" + "=" * 60)');
    lines.push('    print("OPENWQ CALIBRATION")');
    lines.push('    print("=" * 60)');
    lines.push('');
    lines.push('    print("\\nLoading model config...")');
    lines.push('    model_cfg = config_integration.load_model_config(model_config_path)');
    lines.push('');
    lines.push('    # Settings echoed into the results report (shared by the final');
    lines.push('    # report and the --report regeneration below).');
    lines.push('    cal_settings = {');
    lines.push('        "algorithm": algorithm, "max_evaluations": max_evaluations,');
    lines.push('        "objective_function": objective_function,');
    lines.push('        "temporal_resolution": temporal_resolution,');
    lines.push('        "aggregation_method": aggregation_method,');
    lines.push('        "calibration_targets": calibration_targets,');
    lines.push('        "objective_weights": objective_weights,');
    lines.push('        "random_seed": random_seed,');
    lines.push('        "run_sensitivity_first": run_sensitivity_first,');
    lines.push('        "calibration_mode": calibration_mode,');
    lines.push('        "sensitivity_method": sensitivity_method,');
    lines.push('        # Lets --report estimate Morris/Sobol screening progress.');
    lines.push('        "sensitivity_morris_trajectories": sensitivity_morris_trajectories,');
    lines.push('        # Echo HPC + primary-only choices so the results report\'s');
    lines.push('        # Apptainer / SLURM snippets show the actual user paths.');
    lines.push('        "use_primary_only": use_primary_only,');
    lines.push('        "apptainer_sif_path": apptainer_sif_path,');
    lines.push('        "apptainer_bind_path": apptainer_bind_path,');
    lines.push('    }');
    lines.push('');
    lines.push('    # --report: rebuild the results report from the latest on-disk state');
    lines.push('    #   (works while a run is still going — launch in a 2nd terminal),');
    lines.push('    #   open it, and stop without running any evaluations.');
    lines.push('    if args.report:');
    lines.push('        rp = regenerate_results_report(');
    lines.push('            calibration_work_dir=calibration_work_dir,');
    lines.push('            model_config=model_cfg,');
    lines.push('            calibration_parameters=calibration_parameters,');
    lines.push('            calibration_settings=cal_settings,');
    lines.push('            report_stem=report_stem,');
    lines.push('        )');
    lines.push('        if rp:');
    lines.push('            print(f"Results report (regenerated from latest state): {rp}")');
    lines.push('            webbrowser.open("file://" + os.path.abspath(rp))');
    lines.push('        else:');
    lines.push('            print("No results available yet — nothing to report.")');
    lines.push('        sys.exit(0)');
    lines.push('');
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
    lines.push('        # Workflow: "sensitivity" | "both" | "calibration".');
    lines.push('        calibration_mode=calibration_mode,');
    lines.push('        # Lets the driver refresh the results report mid-run.');
    lines.push('        report_stem=report_stem,');
    lines.push('        sensitivity_method=sensitivity_method,');
    lines.push('        sensitivity_morris_trajectories=sensitivity_morris_trajectories,');
    lines.push('        sensitivity_morris_levels=sensitivity_morris_levels,');
    lines.push('        sensitivity_sobol_samples=sensitivity_sobol_samples,');
    lines.push('        sensitivity_threshold=sensitivity_threshold,');
    lines.push('        resume=args.resume,');
    lines.push('        clean_work_dir=clean_work_dir,');
    lines.push('        # Force openWQ debug mode off (default) so evals skip the');
    lines.push('        # 5 diagnostic HDF5 files per species/compartment.');
    lines.push('        run_mode_debug=run_mode_debug,');
    lines.push('        # Per-eval calibration window (None = full model period).');
    lines.push('        calibration_period=calibration_period,');
    lines.push('        temporal_resolution=temporal_resolution,');
    lines.push('        aggregation_method=aggregation_method,');
    // Runtime chosen in the setup report (docker / apptainer), overriding
    // whatever the model config defaulted to.
    lines.push('        container_runtime=container_runtime,');
    lines.push('        docker_container_name=container_config.get("docker_container_name", "docker_openwq"),');
    lines.push('        docker_compose_path=container_config.get("docker_compose_path", ""),');
    lines.push('        # On HPC the sbatch exports OWQ_EXEC_PATH (the binary you');
    lines.push('        # compiled there); locally fall back to the model config path.');
    lines.push('        executable_full_path=(os.environ.get("OWQ_EXEC_PATH")');
    lines.push('                              or container_config.get("executable_path", "")),');
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
    lines.push('    # Generate the final results report (cal_settings built above).');
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
    // Each highlighting pass below wraps tokens in <span class="syn-...">.
    // The trouble: the *next* pass would then re-scan that markup and, e.g.,
    // the keyword pass would match the Python keyword "class" *inside* the
    // class="syn-..." attribute and mangle the tag.  To avoid this entirely
    // we replace every highlighted token with a placeholder made of
    // Private-Use-Area characters (which no later regex can match: they are
    // not \w, not digits, not [a-zA-Z_]) and restore the real spans only at
    // the very end.
    var stash = [];
    // SENT = Private-Use delimiter; built via fromCharCode so this
    // source file stays pure ASCII (no invisible characters).
    var SENT = String.fromCharCode(0xE000);
    function stub(cls, text) {
      var i = stash.length;
      stash.push('<span class="' + cls + '">' + text + '</span>');
      var hi = 0xE001 + Math.floor(i / 6399);
      var lo = 0xE001 + (i % 6399);
      return SENT + String.fromCharCode(hi) + String.fromCharCode(lo) + SENT;
    }
    // Comments and strings together, left-to-right, so a '#' inside a string
    // is not mistaken for a comment (and a quote inside a comment is ignored).
    escaped = escaped.replace(
      /(#[^\n]*)|("(?:[^"\\]|\\.)*"|'(?:[^'\\]|\\.)*')/g,
      function(m, comment, str) {
        return comment ? stub('syn-comment', comment) : stub('syn-str', str);
      });
    // Constants: True, False, None
    escaped = escaped.replace(/\b(True|False|None)\b/g,
      function(m) { return stub('syn-const', m); });
    // Keywords
    escaped = escaped.replace(/\b(import|from|def|if|elif|else|return|as|for|in|not|and|or|class|try|except|finally|with|raise|pass|break|continue|while|yield|lambda|global|nonlocal|assert|del|is)\b/g,
      function(m) { return stub('syn-kw', m); });
    // Numbers (integers and floats)
    escaped = escaped.replace(/\b(\d+\.?\d*(?:[eE][+-]?\d+)?)\b/g,
      function(m) { return stub('syn-num', m); });
    // Function calls
    escaped = escaped.replace(/\b([a-zA-Z_]\w*)\s*\(/g, function(m, fn) {
      if (['if','for','while','return','and','or','not','in','is','import','from','def','class','print'].indexOf(fn) >= 0) return m;
      return stub('syn-fn', fn) + '(';
    });
    // Restore the real <span> markup.
    var RANGE = '[' + String.fromCharCode(0xE001) + '-' +
                String.fromCharCode(0xF8FF) + ']';
    var restoreRe = new RegExp(SENT + '(' + RANGE + ')(' + RANGE + ')' + SENT, 'g');
    escaped = escaped.replace(restoreRe,
      function(m, a, b) {
        return stash[(a.charCodeAt(0) - 0xE001) * 6399 + (b.charCodeAt(0) - 0xE001)];
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
  // Container-runtime toggle (two buttons, single-select).  Mirrors the choice
  // into the hidden #container_runtime input, flips the active button, and sets
  // the stale-folder default (docker -> prompt, apptainer -> clean).
  window.setRuntime = function(val) {
    if (val !== 'docker' && val !== 'apptainer') val = 'docker';
    var hid = document.getElementById('container_runtime');
    if (hid) hid.value = val;
    document.querySelectorAll('#runtimeToggle .rt-btn').forEach(function(b) {
      var on = (b.getAttribute('data-rt') === val);
      b.classList.toggle('rt-active', on);
      b.setAttribute('aria-pressed', on ? 'true' : 'false');
    });
    var cwd = document.getElementById('clean_work_dir');
    if (cwd) cwd.value = (val === 'apptainer') ? 'clean' : 'prompt';
    if (typeof updateRuntimeUI === 'function') updateRuntimeUI();
    if (typeof updateScript === 'function') updateScript();
  };

  function updateRuntimeUI() {
    var crEl = document.getElementById('container_runtime');
    var rt = crEl ? crEl.value : 'docker';
    var txt = document.getElementById('step3Text');
    var cmd = document.getElementById('step3DockerCmd');
    if (rt === 'apptainer') {
      if (txt) txt.textContent = 'Apptainer / Singularity: build your openwq ' +
        '.sif image ON the HPC (so it matches the compute-node architecture), ' +
        'fill in the HPC fields in the Execution tab, then use the ' +
        '"Run the calibration on HPC" section below (no "docker compose" needed).';
      if (cmd) cmd.style.display = 'none';
    } else {
      if (txt) txt.textContent = 'Start the Docker container:';
      if (cmd) cmd.style.display = '';
    }
  }

  // Expand / collapse a numbered run step (steps 2 & 3 start collapsed).
  window.toggleStep = function(id, hdr) {
    var body = document.getElementById(id);
    if (!body) return;
    var hidden = (body.style.display === 'none' || body.style.display === '');
    body.style.display = hidden ? 'block' : 'none';
    var caret = hdr ? hdr.querySelector('.step-caret') : null;
    if (caret) caret.style.transform = hidden ? 'rotate(90deg)' : 'rotate(0deg)';
  };

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

  // Copy the .sbatch save location (next to the run script) to the clipboard.
  window.copySbatchHint = function() {
    var el = document.getElementById('sbatchSaveHint');
    var btn = document.getElementById('copySbatchHintBtn');
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
    updateSbatch(state);
    updateHpcRun(state);
    updateProgress();
  }

  // \u2500\u2500 SLURM sbatch (Execution tab) \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500
  // Assemble the SLURM job file from the Execution-tab fields + the baked
  // run-script path.  Kept consistent with the copy & submit snippets:
  // it is parameterised by HPC_BASE (passed at submit time) and runs the
  // generated calibration script with --clean.
  function sbatchName() {
    return (REPORT_STEM ? REPORT_STEM : 'calibration') + '.sbatch';
  }
  // Map a baked LOCAL absolute path to its clean location under $HPC_BASE.
  // Everything is rsync'd RELATIVE to the common local root, so on the cluster
  // it lives at $HPC_BASE/<path-relative-to-common-root> (no /Users/... prefix).
  function _owqRelToRoot(p) {
    var r = (HPC && HPC.common_root) || '';
    return (r && p && p.indexOf(r) === 0) ? p.slice(r.length) : p;
  }
  function buildSbatch(s) {
    var runLocalAbs = (HPC && HPC.run_local) ? HPC.run_local
                 : (CALIB_LIB_DIR + '/' +
                    (REPORT_STEM ? REPORT_STEM : 'calibration') + '_run.py');
    var runLocal = _owqRelToRoot(runLocalAbs);
    var name = sbatchName();
    var L = [];
    L.push('#!/bin/bash');
    L.push('#SBATCH --job-name=' + (s.slurm_job_name || 'openwq_calib'));
    if (s.slurm_account)   L.push('#SBATCH --account=' + s.slurm_account);
    else L.push('#SBATCH --account=def-YOUR_ACCOUNT      # <-- EDIT: your SLURM allocation');
    if (s.slurm_partition) L.push('#SBATCH --partition=' + s.slurm_partition);
    L.push('#SBATCH --time=' + (s.slurm_time || '48:00:00'));
    L.push('#SBATCH --cpus-per-task=' + (s.slurm_cpus || 4) + '               # >= n_parallel');
    L.push('#SBATCH --mem=' + (s.slurm_mem || '16G'));
    L.push('#SBATCH --output=' + (s.slurm_job_name || 'openwq_calib') + '_%j.out');
    L.push('set -euo pipefail');
    L.push(': "${HPC_BASE:?Pass HPC_BASE, e.g. sbatch --export=ALL,HPC_BASE=/scratch/$USER/openwq_cal ' + name + '}"');
    L.push('');
    if (s.slurm_modules) {
      L.push('# Make Python + Apptainer available on the compute node:');
      L.push(s.slurm_modules);
    } else {
      L.push('# --- EDIT: make Python + Apptainer available on the compute node ---');
      L.push('# module load apptainer python      # or: source ~/miniconda3/bin/activate openwq');
    }
    L.push('');
    L.push('# calibration_lib is used from your cloned openWQ on the HPC (not copied).');
    var owqDir = (s.hpc_openwq_dir || '').replace(/\/+$/, '') || '/path/to/openwq';
    L.push('# OWQ_REPO = your cloned openWQ (the folder that contains supporting_scripts/).');
    L.push('OWQ_REPO="' + owqDir + '"' + (s.hpc_openwq_dir ? '' : '   # <-- EDIT: your cloned openWQ'));
    L.push('# Locate calibration_lib (accepts the openWQ repo folder OR its parent).');
    L.push('for _c in "$OWQ_REPO/supporting_scripts/3_Calibration" "$OWQ_REPO/openwq/supporting_scripts/3_Calibration"; do');
    L.push('  [ -d "$_c/calibration_lib" ] && export OWQ_CALIB_LIB_DIR="$_c" && break');
    L.push('done');
    L.push(': "${OWQ_CALIB_LIB_DIR:?No supporting_scripts/3_Calibration/calibration_lib under $OWQ_REPO - check the openWQ clone path}"');
    L.push('');
    L.push('# Model binary you compiled ON the HPC (the .sif is dependency-only).');
    var _exe = (s.hpc_exe || '').trim();
    if (_exe)
      L.push('export OWQ_EXEC_PATH="' + _exe + '"');
    else
      L.push('export OWQ_EXEC_PATH="/path/to/your/bin/<model>_openwq_Release"   # <-- EDIT: your HPC-compiled binary');
    L.push('');
    L.push('RUN="$HPC_BASE' + runLocal + '"');
    L.push('cd "$(dirname "$RUN")"');
    L.push('python "$(basename "$RUN")" --clean');
    L.push('');
    return L.join('\n');
  }
  function updateSbatch(state) {
    var pre = document.getElementById('sbatchPreview');
    if (!pre) return;            // Execution tab not present
    var s = state || collectFormState();
    pre.textContent = buildSbatch(s);
    var nm = document.getElementById('sbatchName');
    if (nm) nm.textContent = '\u2014 ' + sbatchName();
  }

  // Save the sbatch to disk (Save-As dialog when available).  Deliberately
  // a SAVE (not Copy): it's a file the copy snippet then rsyncs to the HPC.
  window.downloadSbatch = function() {
    var s = collectFormState();
    var text = buildSbatch(s);
    var blob = new Blob([text], {type: 'text/x-sh'});
    var btn = document.getElementById('saveSbatchBtn');
    var sugName = sbatchName();
    function onSaved() {
      if (btn) {
        btn.textContent = '\u2713 Saved!';
        setTimeout(function(){ btn.textContent = 'Save .sbatch'; }, 2000);
      }
    }
    if (window.showSaveFilePicker) {
      window.showSaveFilePicker({
        suggestedName: sugName,
        id: 'owqCalibSave',            // shared dir memory with the run-script save
        startIn: window._owqSaveHandle || 'documents',  // reopen the last-used folder
        types: [{description: 'SLURM job file',
                 accept: {'text/x-sh': ['.sbatch', '.sh']}}],
      }).then(function(handle) {
        window._owqSaveHandle = handle;   // remember this folder for the next save
        return handle.createWritable().then(function(writable) {
          return writable.write(blob).then(function() { return writable.close(); });
        });
      }).then(onSaved).catch(function(err) {
        if (err.name !== 'AbortError') console.error(err);
      });
    } else {
      var url = URL.createObjectURL(blob);
      var a = document.createElement('a');
      a.href = url; a.download = sugName;
      document.body.appendChild(a); a.click(); document.body.removeChild(a);
      URL.revokeObjectURL(url);
      onSaved();
    }
  };

  // ── HPC copy & run blocks (Run-on-HPC section) ─────────────────────────
  // Three independently-pasteable blocks, each pre-filled from the
  // Execution-tab HPC fields + the paths baked into HPC.  Each block restates
  // the HPC_USER/HOST/BASE vars so it works on its own (a fresh terminal),
  // and they run in order: (a) copy, (b) re-point paths, (c) submit.
  function _hpcVars(s) {
    return [
      // zsh (macOS default) does NOT treat '#' as a comment when you paste into
      // an interactive shell, so the comment lines below would error.  This line
      // enables it; in bash/sh `setopt` does not exist, so it is a silent no-op.
      'setopt interactive_comments 2>/dev/null || true',
      'HPC_USER=' + (s.hpc_user || 'your_username'),
      'HPC_HOST=' + (s.hpc_host || 'hpc.your-institution.edu'),
      'HPC_BASE=' + (s.hpc_base || '/scratch/$USER/openwq_cal'),
    ];
  }
  // Parent directory of a path (no trailing slash, never returns "").
  function _owqDirOf(p) {
    var q = (p || '').replace(/\/+$/, '');
    var i = q.lastIndexOf('/');
    return i > 0 ? q.slice(0, i) : q;
  }
  function buildHpcCopy(s) {
    var roots = (HPC && HPC.src_roots) || [];
    var L = _hpcVars(s);
    L.push('SIF_HPC=' + (s.sif_hpc || '$HPC_BASE/openwq.sif'));
    L.push('');
    L.push('# The openwq .sif is NOT uploaded - build it ONCE on the HPC itself,');
    L.push('# so its CPU architecture matches the compute nodes, e.g.:');
    L.push('#   ssh "$HPC_USER@$HPC_HOST"');
    L.push('#   mkdir -p "$(dirname "$SIF_HPC")" && cd "$(dirname "$SIF_HPC")"');
    L.push('#   apptainer build "$SIF_HPC" <openwq>/containers/openwq_apptainer.def');
    L.push('');
    L.push('# Copy your model config + inputs to the HPC.  The calibration_lib CODE');
    L.push('# is NOT copied - it is used from your cloned openWQ on the HPC (see the');
    L.push('# OWQ_CALIB_LIB_DIR export in the .sbatch).');
    L.push('# Each tree is dropped straight into $HPC_BASE/<project>/... : only the');
    L.push('# path BELOW your project root is kept, never the local machine prefix.');
    L.push('# -L dereferences any input symlinks to real files.  eval_*/results +');
    L.push('# old HDF5 + __pycache__ are skipped.');
    // Shared rsync filter lines, applied to every source.  None are "/"-anchored,
    // so they match at any depth regardless of which source tree is being copied.
    var FILT = [
      "  --exclude '0_*_calibration/evaluations' \\",
      "  --exclude '0_*_calibration/results' \\",
      "  --exclude 'openwq_out/HDF5' \\",
      "  --exclude '__pycache__' \\",
      "  --exclude '.DS_Store' \\",
    ];
    // Keep this run's host-model dirs (e.g. 0_SUMMA_OPENWQ + its calibration dir);
    // drop sibling setups for OTHER host models (e.g. 0_MIZUROUTE_OPENWQ).
    ((HPC && HPC.keep_dirs) || []).forEach(function(d){
      FILT.push("  --include '" + d + "/***' \\");
    });
    FILT.push("  --exclude '0_*_OPENWQ' \\");
    FILT.push("  --exclude '0_*_OPENWQ_calibration' \\");
    var _croot = (HPC && HPC.common_root) || '';
    var _clean = !!(HPC && HPC.remap_ok) && !!_croot;
    if (_clean) {
      // Group sources by the parent dir they share UNDER the common root, then
      // rsync each group (NO -R) straight into $HPC_BASE/<that-parent>/.  macOS
      // ships openrsync, which ignores the "/./" dot-anchor and would otherwise
      // recreate the whole local /Users/... path on the cluster - so we set the
      // destination explicitly and copy each source by its basename instead.
      var groups = {};      // relParent -> [src, ...]
      var order = [];
      roots.forEach(function(r){
        var rel = _owqRelToRoot(_owqDirOf(r));   // e.g. "/century_basins_tested" or ""
        if (!(rel in groups)) { groups[rel] = []; order.push(rel); }
        groups[rel].push(r);
      });
      order.forEach(function(rel){
        var dest = '$HPC_BASE' + rel;
        L.push('ssh "$HPC_USER@$HPC_HOST" "mkdir -p ' + dest + '"');
        L.push('rsync -avzL \\');
        FILT.forEach(function(f){ L.push(f); });
        groups[rel].forEach(function(src){ L.push('  "' + src + '" \\'); });
        L.push('  "$HPC_USER@$HPC_HOST:' + dest + '/"');
      });
    } else {
      // Fallback: no clean common root (inputs span drives).  Copy each tree by
      // its basename into $HPC_BASE and re-point the absolute paths by hand.
      L.push('ssh "$HPC_USER@$HPC_HOST" "mkdir -p $HPC_BASE"');
      L.push('rsync -avzL \\');
      FILT.forEach(function(f){ L.push(f); });
      roots.forEach(function(r){ L.push('  "' + r + '" \\'); });
      L.push('  "$HPC_USER@$HPC_HOST:$HPC_BASE/"');
    }
    return L.join('\n');
  }
  function buildHpcRemap(s) {
    var L = _hpcVars(s);
    L.push('SIF_HPC=' + (s.sif_hpc || '$HPC_BASE/openwq.sif'));
    L.push('');
    if (HPC && HPC.remap_ok) {
      L.push('# Re-point the absolute paths to their HPC copies (runs once).');
      L.push('# Run this on your LAPTOP or after logging into the HPC - it auto-detects:');
      L.push('# if sbatch is on PATH it edits the files here, else it ssh-es into the HPC.');
      L.push('if command -v sbatch >/dev/null 2>&1; then');
      L.push('  owq_sh() { HPC_BASE="$HPC_BASE" SIF_HPC="$SIF_HPC" bash; }');
      L.push('else');
      L.push('  owq_sh() { ssh "$HPC_USER@$HPC_HOST" "HPC_BASE=\'$HPC_BASE\' SIF_HPC=\'$SIF_HPC\' bash -s"; }');
      L.push('fi');
      L.push('owq_sh <<\'EOF\'');
      L.push('set -euo pipefail');
      L.push('ROOT="' + (HPC.common_root || '') + '"');
      L.push('RUN="$HPC_BASE' + _owqRelToRoot(HPC.run_local || '') + '"');
      L.push('MC="$HPC_BASE' + _owqRelToRoot(HPC.mc_path || '') + '"');
      var _bindVal = (HPC && HPC.container_prefix) ? ('$HPC_BASE:' + HPC.container_prefix) : '$HPC_BASE';
      L.push('# Idempotent + self-healing: re-point only while local $ROOT paths remain,');
      L.push('# so re-copying fresh inputs later just re-points them again (no stale flag).');
      L.push('if grep -qF "$ROOT" "$RUN" "$MC"; then');
      L.push('  # Re-point every baked-in local root ($ROOT) to its clean $HPC_BASE copy.');
      L.push('  sed -i "s#$ROOT#$HPC_BASE#g" "$RUN" "$MC"');
      L.push('  sed -i "s#^apptainer_sif_path = .*#apptainer_sif_path = \\"$SIF_HPC\\"#" "$RUN"');
      L.push('  sed -i "s#^apptainer_bind_path = .*#apptainer_bind_path = \\"' + _bindVal + '\\"#" "$RUN"');
      L.push('  sed -i "s#^container_runtime = .*#container_runtime = \\"apptainer\\"#" "$RUN"');
      L.push('  echo "Re-pointed absolute paths to $HPC_BASE."');
      L.push('else');
      L.push('  echo "Already re-pointed - no $ROOT paths remain."');
      L.push('fi');
      L.push('EOF');
    } else {
      L.push('# NOTE: your inputs span different drives, so paths cannot be re-pointed');
      L.push('# automatically. After copying, edit the absolute paths in your run script');
      L.push('# + model config to their $HPC_BASE locations by hand.');
    }
    return L.join('\n');
  }
  function buildHpcSubmit(s) {
    var L = _hpcVars(s);
    L.push('');
    var _sbatchHpc = '$HPC_BASE' + ((HPC && HPC.sbatch_local)
                       ? _owqRelToRoot(HPC.sbatch_local)
                       : ('/' + ((HPC && HPC.sbatch_name) || 'calibration.sbatch')));
    L.push('# Submit the SLURM job (it was copied inside the work dir by step a).');
    L.push('# Run this on your LAPTOP or after logging into the HPC - it auto-detects:');
    L.push('# if sbatch is on PATH it submits here, otherwise it ssh-es into the HPC.');
    L.push('SBATCH_FILE="' + _sbatchHpc + '"');
    L.push('if command -v sbatch >/dev/null 2>&1; then');
    L.push('  cd "$HPC_BASE" && sbatch --export=ALL,HPC_BASE="$HPC_BASE" "$SBATCH_FILE"');
    L.push('else');
    L.push('  ssh "$HPC_USER@$HPC_HOST" "cd $HPC_BASE && sbatch --export=ALL,HPC_BASE=$HPC_BASE $SBATCH_FILE"');
    L.push('fi');
    L.push('');
    L.push('# Monitor (on the HPC):  squeue -u $HPC_USER');
    L.push('#       (from laptop):   ssh "$HPC_USER@$HPC_HOST" "squeue -u $HPC_USER"');
    if (HPC && HPC.work_dir)
      L.push('# Results (from laptop): rsync -avz "$HPC_USER@$HPC_HOST:$HPC_BASE' + _owqRelToRoot(HPC.work_dir) + '/" ./hpc_results/');
    return L.join('\n');
  }
  function updateHpcRun(state) {
    var s = state || collectFormState();
    var c = document.getElementById('hpcRunCopy');   if (c) c.textContent = buildHpcCopy(s);
    var r = document.getElementById('hpcRunRemap');  if (r) r.textContent = buildHpcRemap(s);
    var b = document.getElementById('hpcRunSubmit'); if (b) b.textContent = buildHpcSubmit(s);
    // Keep the "expected .sif location" warning in sync with the HPC sif field.
    var sp = document.getElementById('sifExpectedPath');
    if (sp) sp.textContent = (s.sif_hpc && String(s.sif_hpc).trim())
                           ? String(s.sif_hpc).trim() : '$HPC_BASE/openwq.sif';
  }
  window.copyHpcRun = function(btn, targetId) {
    var pre = document.getElementById(targetId);
    if (!pre) return;
    navigator.clipboard.writeText(pre.textContent);
    if (btn) {
      var t = btn.textContent;
      btn.textContent = 'Copied!';
      setTimeout(function(){ btn.textContent = t; }, 1500);
    }
  };

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
        id: 'owqCalibSave',            // shared dir memory with the .sbatch save
        startIn: window._owqSaveHandle || 'documents',  // reopen the last-used folder
        types: [{
          description: 'Python Script',
          accept: {'text/x-python': ['.py']},
        }],
      }).then(function(handle) {
        window._owqSaveHandle = handle;   // remember this folder for the next save
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

  // Workflow mode slider (0=influential params, 1=both, 2=calibration).
  // Shows/hides the Sensitivity card and the Optimization card accordingly.
  function updateWorkflowUI() {
    var el = document.getElementById('calibration_mode');
    var m = el ? parseInt(el.value) : 1;
    if (isNaN(m)) m = 1;
    var sensCard = document.getElementById('sensitivity');
    var optCard = document.getElementById('optCard');
    var desc = document.getElementById('modeDesc');
    if (sensCard) sensCard.style.display = (m === 2) ? 'none' : '';
    if (optCard)  optCard.style.display  = (m === 0) ? 'none' : '';
    // The screening method (Morris or Sobol) is chosen later, in the
    // Influential-parameters settings below — so name both here.
    var txt = (m === 0)
        ? 'Only influential parameters will be identified (Morris or Sobol screening). No calibration will run.'
      : (m === 2)
        ? 'Only calibration will run, using the settings you provide. No sensitivity screening.'
        : 'Both, in sequence: influential parameters are identified first (Morris or Sobol screening), then calibration runs.';
    if (desc) desc.textContent = txt;
    document.querySelectorAll('.mode-seg-btn').forEach(function(b) {
      var on = parseInt(b.dataset.m) === m;
      b.classList.toggle('active', on);
      b.setAttribute('aria-checked', on ? 'true' : 'false');
    });
    var thumb = document.getElementById('modeThumb');
    if (thumb) thumb.style.transform = 'translateX(' + (m * 100) + '%)';
  }
  // Each option is a discrete button; clicking sets the hidden value.
  document.querySelectorAll('.mode-seg-btn').forEach(function(b) {
    b.addEventListener('click', function() {
      var el = document.getElementById('calibration_mode');
      if (el) el.value = this.dataset.m;
      updateWorkflowUI();
      updateScript();
    });
  });

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

  // ── Calibration / validation period slider ──
  // Greys out the span without observations.  Two draggable handles set the
  // calibration window START (left) and the Calibration / Validation split
  // (right); any data before the start handle is excluded.  Defaults: start at
  // the first observation, split at 1/3 (calibration = first third, validation
  // = last two-thirds).  Populates the hidden period
  // inputs that collectFormState() / generateScript() read.
  (function() {
    var wrap   = document.getElementById('calibSliderWrap');
    var track  = document.getElementById('calibTrack');
    var handle = document.getElementById('calibHandle');
    var handleStart = document.getElementById('calibHandleStart');
    if (!wrap || !track || !handle) return;

    function parseDate(s) {
      if (!s) return null;
      // Accept 'YYYY-MM-DD HH:MM[:SS]' or ISO; normalise space->T (Safari).
      var d = new Date(String(s).trim().replace(' ', 'T'));
      return isNaN(d.getTime()) ? null : d;
    }
    function p2(n) { return (n < 10 ? '0' : '') + n; }
    function fmt(d) {
      if (!d) return '—';
      return d.getFullYear() + '-' + p2(d.getMonth() + 1) + '-' + p2(d.getDate())
           + ' ' + p2(d.getHours()) + ':' + p2(d.getMinutes());
    }

    var simStart = SIM_PERIOD ? parseDate(SIM_PERIOD.sim_start) : null;
    var simEnd   = SIM_PERIOD ? parseDate(SIM_PERIOD.sim_end)   : null;
    var obsStart = OBS_PERIOD ? parseDate(OBS_PERIOD.obs_start) : null;
    var obsEnd   = OBS_PERIOD ? parseDate(OBS_PERIOD.obs_end)   : null;

    // Outer track range [lo, hi] = the model's simulation period (the hard
    // bound — you can't calibrate where the model doesn't run).  The active
    // (calibratable) sub-range [aStart, aEnd] = observations ∩ model period;
    // the rest of the track (model period without obs) is greyed out.
    var haveObs = obsStart && obsEnd && obsEnd > obsStart;
    var haveSim = simStart && simEnd && simEnd > simStart;
    var lo, hi, aStart, aEnd;
    if (haveSim) {
      lo = simStart; hi = simEnd;
      if (haveObs) {
        aStart = (obsStart > simStart) ? obsStart : simStart;
        aEnd   = (obsEnd   < simEnd)   ? obsEnd   : simEnd;
        if (aEnd <= aStart) {        // obs disjoint from model period
          aStart = simStart; aEnd = simEnd;
        }
      } else {
        aStart = simStart; aEnd = simEnd;
      }
    } else if (haveObs) {
      lo = obsStart; hi = obsEnd; aStart = obsStart; aEnd = obsEnd;
    } else {
      wrap.style.display = 'none';
      var msg = document.getElementById('calibNoDataMsg');
      if (msg) msg.style.display = '';
      var ucp0 = document.getElementById('use_calibration_period');
      if (ucp0) ucp0.checked = false;
      return;
    }
    var span = hi.getTime() - lo.getTime();
    if (span <= 0) { wrap.style.display = 'none'; return; }

    function frac(d) { return (d.getTime() - lo.getTime()) / span; }
    function dateAtFrac(f) { return new Date(lo.getTime() + f * span); }
    var aStartF = frac(aStart), aEndF = frac(aEnd);
    // default: calibration = first 1/3 of the obs period, validation = last 2/3
    var splitF = aStartF + (aEndF - aStartF) / 3;
    var calStartF = aStartF;              // calibration-window START (left handle)

    // Observation timestamps (epoch ms) per species → live in-window count.
    // Counts only the SELECTED target species when any are ticked, else all
    // species present in the obs data.  (Optimistic: spatial / primary-only
    // filtering happens at run time, so this is an upper bound, but it
    // cleanly separates an empty window from a populated one.)
    var OBS_MS = (OBS_PERIOD && OBS_PERIOD.dates_by_species) || {};
    function selectedSpecies() {
      var sel = [];
      document.querySelectorAll('.species-cb:checked').forEach(function(cb) {
        if (cb.dataset && cb.dataset.species) sel.push(cb.dataset.species);
      });
      return sel;
    }
    function countObs(loMs, hiMs) {
      var keys = Object.keys(OBS_MS);
      if (keys.length === 0) return null;   // no per-species obs dates available
      var sp = selectedSpecies();
      var use = sp.length ? keys.filter(function(k){ return sp.indexOf(k) >= 0; }) : keys;
      if (use.length === 0) use = keys;     // selected species absent in obs → count all
      var n = 0;
      use.forEach(function(k) {
        var arr = OBS_MS[k];
        for (var i = 0; i < arr.length; i++) {
          if (arr[i] >= loMs && arr[i] <= hiMs) n++;
        }
      });
      return n;
    }

    function render() {
      var pct = function(f) { return (f * 100) + '%'; };
      var gl = document.getElementById('segGrayLeft');
      var ex = document.getElementById('segExclLeft');
      var sc = document.getElementById('segCalib');
      var sv = document.getElementById('segValid');
      var gr = document.getElementById('segGrayRight');
      gl.style.left = '0%';            gl.style.width = pct(aStartF);
      ex.style.left = pct(aStartF);    ex.style.width = pct(calStartF - aStartF);
      sc.style.left = pct(calStartF);  sc.style.width = pct(splitF - calStartF);
      sv.style.left = pct(splitF);     sv.style.width = pct(aEndF - splitF);
      gr.style.left = pct(aEndF);      gr.style.width = pct(1 - aEndF);
      handle.style.left = pct(splitF);
      if (handleStart) handleStart.style.left = pct(calStartF);

      var cStart = dateAtFrac(calStartF);
      var cEnd = dateAtFrac(splitF);
      var cN = countObs(cStart.getTime(), cEnd.getTime());
      var vN = countObs(cEnd.getTime(), aEnd.getTime());
      var cSuf = (cN === null) ? '' : ('  ·  ' + cN + ' obs');
      var vSuf = (vN === null) ? '' : ('  ·  ' + vN + ' obs');
      var cEl = document.getElementById('calibWindowText');
      cEl.textContent = fmt(cStart) + '  →  ' + fmt(cEnd) + cSuf;
      // Flag an empty calibration window (no eligible obs) in red.
      cEl.style.color = (cN === 0) ? 'var(--danger, #dc2626)' : '';
      cEl.style.fontWeight = (cN === 0) ? '700' : '';
      document.getElementById('valWindowText').textContent   = fmt(cEnd)   + '  →  ' + fmt(aEnd) + vSuf;
      document.getElementById('calibAxisStart').textContent  = fmt(lo);
      document.getElementById('calibAxisEnd').textContent    = fmt(hi);
      document.getElementById('calibration_period_start').value = fmt(cStart);
      document.getElementById('calibration_period_end').value   = fmt(cEnd);
      document.getElementById('validation_period_start').value  = fmt(cEnd);
      document.getElementById('validation_period_end').value    = fmt(aEnd);
    }

    var active = 'split';   // which handle is being dragged: 'start' | 'split'
    function setFromClientX(clientX) {
      var rect = track.getBoundingClientRect();
      var f = (clientX - rect.left) / rect.width;
      if (active === 'start') {
        calStartF = Math.max(aStartF, Math.min(splitF, f));   // can't pass the split
      } else {
        splitF = Math.max(calStartF, Math.min(aEndF, f));     // can't pass the start
      }
      render();
      if (typeof updateScript === 'function') updateScript();
    }

    var dragging = false;
    handle.addEventListener('mousedown', function(e) {
      e.preventDefault(); active = 'split'; dragging = true;
      document.body.style.userSelect = 'none';
    });
    if (handleStart) handleStart.addEventListener('mousedown', function(e) {
      e.preventDefault(); active = 'start'; dragging = true;
      document.body.style.userSelect = 'none';
    });
    track.addEventListener('mousedown', function(e) {
      // Clicking the track grabs whichever handle is nearer the click point.
      var rect = track.getBoundingClientRect();
      var f = (e.clientX - rect.left) / rect.width;
      active = (Math.abs(f - calStartF) < Math.abs(f - splitF)) ? 'start' : 'split';
      dragging = true; document.body.style.userSelect = 'none';
      setFromClientX(e.clientX);
    });
    document.addEventListener('mousemove', function(e) {
      if (dragging) setFromClientX(e.clientX);
    });
    document.addEventListener('mouseup', function() {
      if (dragging) { dragging = false; document.body.style.userSelect = ''; }
    });

    // Enable/disable the slider with the "use calibration window" checkbox.
    var ucp = document.getElementById('use_calibration_period');
    if (ucp) {
      ucp.addEventListener('change', function() {
        track.style.opacity = this.checked ? '1' : '.4';
        track.style.pointerEvents = this.checked ? '' : 'none';
        if (typeof updateScript === 'function') updateScript();
      });
    }

    // Per-species observation coverage lanes, drawn ABOVE the track on the
    // SAME time axis [lo, hi].  One lane per species present in the obs data
    // (ALL observations, not just the metric's primary ones) so a species
    // with only secondary obs (e.g. NH4) still shows up — making it obvious
    // why a configured species can contribute nothing to the fit.  Primary
    // (scored) obs are drawn as solid ticks; secondary as faint ticks.
    // Colours deliberately avoid the calibration (blue) / validation (green)
    // hues used by the track below.
    //
    // Robust for ANY number of observations: the COUNTS + date range come from
    // species_summary (true totals computed in Python on the uncapped data),
    // while the ticks are subsampled for display so the DOM stays light even
    // for hundreds of thousands of points.
    var OBS_ALL = (OBS_PERIOD && OBS_PERIOD.dates_by_species_all) || OBS_MS;
    var OBS_SUM = (OBS_PERIOD && OBS_PERIOD.species_summary) || {};
    var LANE_MAX_TICKS = 600;   // max rendered ticks per layer per species
    function renderSpeciesLanes() {
      var host = document.getElementById('calibSpeciesLanes');
      if (!host) return;
      var keys = Object.keys(OBS_ALL).sort();
      if (!keys.length) { host.innerHTML = ''; return; }
      var loMs = lo.getTime(), hiMs = hi.getTime(), rng = hiMs - loMs;
      // Warm / purple palette — distinct from track blue (#2563eb) & green (#10b981).
      var palette = ['#ea580c', '#9333ea', '#db2777', '#ca8a04', '#be123c', '#7c3aed'];
      // The metric scores primary obs only while "use_primary_only" is ON;
      // when it's OFF every observation is scored. Reflect that live.
      var _upo = document.getElementById('use_primary_only');
      var primaryOnly = _upo ? !!_upo.checked : true;
      // Some host models (mizuRoute) have no "secondary" stations — every
      // matched obs is primary. Detect that from the data so the legend +
      // per-species notes don't reference a category that doesn't exist.
      var hasSecondary = keys.some(function(k) {
        var sm = OBS_SUM[k]; return sm ? (sm.n > sm.n_primary) : false;
      });
      var legend;
      if (!hasSecondary) {
        legend = 'all observations are at matched target station(s), scored';
      } else if (primaryOnly) {
        legend = '<span style="opacity:.95;">solid&nbsp;=&nbsp;primary (scored)</span>, '
               + '<span style="opacity:.55;">faint&nbsp;=&nbsp;secondary (not scored)</span>';
      } else {
        legend = '<span style="opacity:.95;">solid&nbsp;=&nbsp;primary</span>, '
               + '<span style="opacity:.55;">faint&nbsp;=&nbsp;secondary</span> &mdash; all scored';
      }
      var html = '<div style="font-size:.62rem;color:var(--text3);margin:.1rem 0 .25rem;">'
               + 'Observations by species (same time axis as the slider below) &mdash; '
               + legend + ':</div>';
      // Build ticks for one layer: filter to the visible axis, subsample to
      // LANE_MAX_TICKS for display, and count the in-window points (fallback).
      function tickLayer(arr, color, opacity, width) {
        var inr = [];
        for (var i = 0; i < arr.length; i++) {
          var t = arr[i];
          if (t >= loMs && t <= hiMs) inr.push(t);
        }
        var draw = inr;
        if (inr.length > LANE_MAX_TICKS) {
          var step = Math.ceil(inr.length / LANE_MAX_TICKS);
          draw = [];
          for (var j = 0; j < inr.length; j += step) draw.push(inr[j]);
        }
        var s = '';
        for (var m = 0; m < draw.length; m++) {
          var f = rng > 0 ? (draw[m] - loMs) / rng : 0;
          s += '<span style="position:absolute;left:' + (f * 100) + '%;top:0;bottom:0;width:'
             + width + 'px;margin-left:-' + (width / 2) + 'px;background:' + color
             + ';opacity:' + opacity + ';"></span>';
        }
        return { html: s, n: inr.length };
      }
      keys.forEach(function(k, idx) {
        var color = palette[idx % palette.length];
        var allT  = tickLayer(OBS_ALL[k] || [], color, 0.28, 2);   // all (faint)
        var primT = tickLayer(OBS_MS[k]  || [], color, 0.95, 2);   // primary (solid)
        // True counts + full date range (accurate for any number of obs);
        // fall back to the in-window tick counts if the summary is missing.
        var sm = OBS_SUM[k] || null;
        var nAll  = sm ? sm.n         : allT.n;
        var nPrim = sm ? sm.n_primary : primT.n;
        var startMs = sm ? sm.start : null, endMs = sm ? sm.end : null;
        var span = '';
        if (startMs !== null && rng > 0) {
          var cs = Math.max(startMs, loMs), ce = Math.min(endMs, hiMs);
          if (ce >= cs) {
            var f0 = (cs - loMs) / rng, f1 = (ce - loMs) / rng;
            span = '<span style="position:absolute;left:' + (f0 * 100) + '%;width:'
                 + ((f1 - f0) * 100) + '%;top:50%;height:2px;transform:translateY(-50%);'
                 + 'background:' + color + ';opacity:.22;"></span>';
          }
        }
        var txt;
        if (nAll > 0) {
          var rangeTxt = (startMs !== null)
            ? (fmt(new Date(startMs)).slice(0, 7) + ' → '
               + fmt(new Date(endMs)).slice(0, 7) + '  ·  ')
            : '';
          if (hasSecondary) {
            txt = rangeTxt + nAll
                + ' obs (' + nPrim
                + ' observations in primary target observation station(s))';
            if (nPrim === 0) {
              txt += primaryOnly
                ? ' <span style="color:var(--danger,#dc2626);font-weight:700;">'
                  + '— secondary only, not scored by the metric</span>'
                : ' <span style="color:var(--text3);">'
                  + '— secondary only (scored, since primary-only is off)</span>';
            }
          } else {
            txt = rangeTxt + nAll + ' obs';
          }
        } else {
          txt = 'no observations';
        }
        html += '<div style="margin:.2rem 0;">'
              + '<div style="font-size:.62rem;color:var(--text2);margin-bottom:2px;">'
              + '<span style="display:inline-block;width:8px;height:8px;border-radius:2px;'
              + 'background:' + color + ';margin-right:.35rem;"></span>'
              + '<strong>' + k + '</strong> '
              + '<span style="color:var(--text3);">' + txt + '</span></div>'
              + '<div style="position:relative;height:9px;background:rgba(127,127,127,.12);'
              + 'border-radius:3px;">' + span
              + (hasSecondary ? (allT.html + primT.html) : (primT.html || allT.html))
              + '</div></div>';
      });
      host.innerHTML = html;
    }

    // Re-count when the target-species selection changes (Targets tab).
    document.querySelectorAll('.species-cb').forEach(function(cb) {
      cb.addEventListener('change', render);
    });

    // Re-draw the species lanes when "use only primary observations" toggles,
    // so the scored / not-scored wording tracks the actual metric setting.
    var _upoLane = document.getElementById('use_primary_only');
    if (_upoLane) _upoLane.addEventListener('change', renderSpeciesLanes);

    renderSpeciesLanes();
    render();
  })();

  // Initial render
  updateScript();

  // Workflow mode: reflect the default (Both) — show/hide sensitivity & opt cards.
  updateWorkflowUI();

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
