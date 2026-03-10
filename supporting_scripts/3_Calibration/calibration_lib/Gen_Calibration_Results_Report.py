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
Calibration Results Report
============================

Generates an HTML report AFTER calibration completes. Shows convergence,
best parameters, sensitivity analysis results, performance plots, and
a code snippet to run a simulation with the best parameters.

This is the companion to Gen_Calibration_Setup_Report.py (which is
generated BEFORE calibration).
"""

import os
import io
import json
import html as html_lib
import logging
from datetime import datetime
from typing import List, Dict, Any, Optional, Tuple

from . import report_helpers as rh
from . import config_integration as _ci

logger = logging.getLogger(__name__)

# Try to import plotting libraries (graceful degradation if missing)
try:
    import matplotlib
    matplotlib.use('Agg')  # Non-interactive backend
    import matplotlib.pyplot as plt
    import matplotlib.ticker as mticker
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    logger.warning("matplotlib not available — plots will be skipped in results report")

try:
    import numpy as np
    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False


def generate_results_report(
    output_dir: str,
    model_config: Dict[str, Any],
    calibration_parameters: List[Dict],
    calibration_settings: Dict[str, Any],
    calibration_results: Dict[str, Any],
    sensitivity_results: Optional[Dict[str, Any]] = None,
    performance_metrics: Optional[List[Dict[str, Any]]] = None,
    matched_data: Optional[Any] = None,
) -> Optional[str]:
    """
    Generate the calibration results HTML report.

    Parameters
    ----------
    output_dir : str
        Directory where the report HTML file will be saved.
    model_config : Dict
        Loaded model configuration variables.
    calibration_parameters : List[Dict]
        Parameter definitions (name, initial, bounds, transform, etc.).
    calibration_settings : Dict
        Calibration settings (algorithm, max_evaluations, objective_function, etc.).
    calibration_results : Dict
        Calibration results with keys:
        - best_params : Dict[str, float] — best parameter values
        - best_objective : float
        - n_evaluations : int
        - converged : bool
        - convergence_reason : str
        - history : List[Dict] — [{eval_id, objective, parameters, ...}]
        - runtime_hours : float (optional)
        - convergence_eval : int (optional) — evaluation that found the best
    sensitivity_results : Dict, optional
        Sensitivity analysis results with keys:
        - method : str — "morris" or "sobol"
        - parameter_names : List[str]
        - mu_star / sigma (morris) or S1 / ST (sobol)
    performance_metrics : List[Dict], optional
        Per-species performance metrics. Each dict has:
        - species : str
        - reach_id : str/int
        - KGE, NSE, RMSE, PBIAS, R2 : float
        - obs_mean, sim_mean : float (optional)
    matched_data : DataFrame, optional
        Matched observations vs simulations DataFrame.

    Returns
    -------
    Optional[str]
        Path to the generated HTML report, or None on failure.
    """
    try:
        os.makedirs(output_dir, exist_ok=True)
        report_path = os.path.join(output_dir, "calibration_results_report.html")

        H = []  # HTML content accumulator

        # ── Extra CSS for results-specific elements ──
        extra_css = _get_results_css()

        # ── Head ──
        H.append(rh.build_html_head(
            "OpenWQ Calibration Results Report",
            extra_css=extra_css
        ))
        H.append("<body>")
        H.append('<div class="layout">')

        # ── Sidebar ──
        nav_items = [
            {"id": "summary", "label": "Summary"},
            {"id": "best-params", "label": "Best Parameters"},
            {"id": "convergence", "label": "Convergence"},
            {"id": "param-evolution", "label": "Parameter Evolution"},
            {"id": "param-correlations", "label": "Correlations"},
        ]
        if sensitivity_results:
            nav_items.append({"id": "sensitivity", "label": "Sensitivity"})
        try:
            _bgc_path = _ci.get_bgc_template_path(model_config)
        except Exception:
            _bgc_path = None
        if _bgc_path and os.path.isfile(_bgc_path):
            nav_items.append({"id": "bgc-network", "label": "BGC Network"})
        if performance_metrics:
            nav_items.append({"id": "performance", "label": "Performance"})
        nav_items.append({"id": "run-best", "label": "Run Best Params"})

        H.append(rh.build_sidebar(nav_items, logo_text="OpenWQ Calibration"))

        H.append('<div class="main">')

        # ── Header ──
        project_name = model_config.get("project_name", "Calibration")
        authors = model_config.get("authors", "")
        date_str = datetime.now().strftime("%Y-%m-%d %H:%M")

        best_obj = calibration_results.get("best_objective", float('nan'))
        n_evals = calibration_results.get("n_evaluations", 0)

        meta_items = [
            f"Authors: {authors}" if authors else "",
            f"Completed: {date_str}",
            f"Best {calibration_settings.get('objective_function', 'KGE')}: {best_obj:.4f}",
            f"Evaluations: {n_evals}",
        ]
        meta_items = [m for m in meta_items if m]

        H.append(rh.build_header(
            title=project_name,
            subtitle="Calibration Results Report",
            meta_items=meta_items,
            badge_text="RESULTS",
            badge_class="badge-secondary"
        ))

        H.append('<div class="container">')

        # ── Section: Summary ──
        H.append(_build_summary_section(
            calibration_results, calibration_settings, calibration_parameters
        ))

        # ── Section: Best Parameters ──
        H.append(_build_best_params_section(
            calibration_parameters, calibration_results
        ))

        # ── Section: Convergence ──
        H.append(_build_convergence_section(
            calibration_results, calibration_settings, output_dir
        ))

        # ── Section: Parameter Evolution ──
        H.append(_build_param_evolution_section(
            calibration_results, calibration_parameters, output_dir
        ))

        # ── Section: Parameter Correlations ──
        H.append(_build_correlations_section(
            calibration_results, calibration_parameters, output_dir
        ))

        # ── Section: Sensitivity Analysis (if available) ──
        if sensitivity_results:
            H.append(_build_sensitivity_section(
                sensitivity_results, output_dir
            ))

        # ── Section: BGC Reaction Network (if template available) ──
        if _bgc_path and os.path.isfile(_bgc_path):
            try:
                _diag_html = rh.generate_bgc_reaction_diagram(
                    bgc_template_path=_bgc_path,
                    module_name=model_config.get("bgc_module_name", ""),
                )
                if _diag_html:
                    H.append('<div class="section" id="bgc-network">')
                    H.append('<h2>BGC Reaction Network</h2>')
                    H.append(_diag_html)
                    H.append('</div>')
            except Exception:
                pass

        # ── Section: Performance Metrics (if available) ──
        if performance_metrics:
            H.append(_build_performance_section(
                performance_metrics, matched_data,
                calibration_settings, output_dir
            ))

        # ── Section: Run Best Parameters ──
        H.append(_build_run_best_section(
            calibration_results, model_config, calibration_settings
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

        logger.info(f"Results report saved to: {report_path}")
        return report_path

    except Exception as e:
        logger.error(f"Failed to generate results report: {e}")
        import traceback
        logger.debug(traceback.format_exc())
        return None


# =========================================================================
# Results-specific CSS
# =========================================================================

def _get_results_css() -> str:
    """Extra CSS specific to the results report."""
    return """
    /* ── Result Indicators ───────────────────────────────────── */
    .result-good { color: var(--secondary); }
    .result-ok   { color: var(--primary); }
    .result-poor { color: var(--accent); }
    .result-bad  { color: #e74c3c; }

    .metric-cell {
        font-family: 'JetBrains Mono', monospace;
        font-weight: 600;
    }

    /* ── Change indicator ────────────────────────────────────── */
    .change-up   { color: var(--secondary); }
    .change-down { color: var(--accent); }
    .change-zero { color: var(--text3); }

    /* ── Sensitivity bars ────────────────────────────────────── */
    .sa-bar-container {
        width: 100%; background: var(--border); border-radius: 4px;
        overflow: hidden; height: 20px; margin: 2px 0;
    }
    .sa-bar {
        height: 100%; border-radius: 4px; min-width: 2px;
        background: linear-gradient(90deg, var(--primary), var(--secondary));
        transition: width .3s;
    }

    /* ── Performance card ─────────────────────────────────────── */
    .perf-card {
        display: grid; grid-template-columns: 1fr 1fr;
        gap: .5rem; margin-bottom: .5rem;
    }
    .perf-metric {
        text-align: center; padding: .6rem; border-radius: 8px;
        background: rgba(0,0,0,.02); border: 1px solid var(--border);
    }
    [data-theme="dark"] .perf-metric { background: rgba(255,255,255,.02); }
    .perf-metric .perf-value {
        font-size: 1.3rem; font-weight: 700;
        font-family: 'JetBrains Mono', monospace;
    }
    .perf-metric .perf-label {
        font-size: .72rem; color: var(--text2);
        text-transform: uppercase; letter-spacing: .8px;
    }
    """


# =========================================================================
# Section Builders
# =========================================================================

def _build_summary_section(
    results: Dict[str, Any],
    settings: Dict[str, Any],
    parameters: List[Dict]
) -> str:
    """Build the results summary section with KPIs."""

    best_obj = results.get("best_objective", float('nan'))
    n_evals = results.get("n_evaluations", 0)
    converged = results.get("converged", False)
    conv_reason = results.get("convergence_reason", "N/A")
    conv_eval = results.get("convergence_eval", "N/A")
    runtime_h = results.get("runtime_hours", 0)
    obj_fn = settings.get("objective_function", "KGE")

    # Format runtime
    if runtime_h >= 1:
        runtime_str = f"{runtime_h:.1f} h"
    elif runtime_h > 0:
        runtime_str = f"{runtime_h * 60:.0f} min"
    else:
        runtime_str = "N/A"

    # Objective value quality assessment
    obj_quality = _assess_objective(best_obj, obj_fn)

    kpis = [
        {"icon": "\U0001f3c6", "value": f"{best_obj:.4f}",
         "label": f"Best {obj_fn}"},
        {"icon": "\U0001f504", "value": str(n_evals),
         "label": "Evaluations"},
        {"icon": "\u23f1\ufe0f", "value": runtime_str,
         "label": "Runtime"},
        {"icon": "\U0001f4ca", "value": str(len(parameters)),
         "label": "Parameters"},
    ]

    kpi_html = rh.build_kpi_grid(kpis)

    # Convergence info box
    conv_icon = "\u2705" if converged else "\u26a0\ufe0f"
    conv_style = "success" if converged else "warning"
    conv_html = rh.build_highlight_box(
        f"<strong>{conv_icon} Convergence:</strong> {html_lib.escape(str(conv_reason))} "
        f"&bull; Best found at evaluation <strong>{conv_eval}</strong> "
        f"of <strong>{n_evals}</strong>",
        conv_style
    )

    # Quality assessment box
    quality_html = rh.build_highlight_box(
        f"<strong>{obj_quality['icon']} {obj_fn} Quality:</strong> "
        f"{obj_quality['label']} ({obj_quality['description']})",
        obj_quality['style']
    )

    return f"""
<div class="section" id="summary">
    <h2>Results Summary</h2>
    {kpi_html}
    {conv_html}
    {quality_html}
</div>
"""


def _build_best_params_section(
    parameters: List[Dict],
    results: Dict[str, Any]
) -> str:
    """Build the best parameter values table."""

    best_params = results.get("best_params", {})

    if not parameters or not best_params:
        return """
<div class="section" id="best-params">
    <h2>Best Parameter Values</h2>
    <div class="card"><p>No parameter results available.</p></div>
</div>
"""

    # Build table data
    table_data = []
    for param in parameters:
        name = param.get("name", "")
        initial = param.get("initial", 0)
        bounds = param.get("bounds", (0, 1))
        units = param.get("units", "")
        transform = param.get("transform", "linear")

        # Best value
        best_val = best_params.get(name, initial)

        # Percent change
        if initial != 0:
            pct_change = ((best_val - initial) / abs(initial)) * 100
        else:
            pct_change = 0.0 if best_val == 0 else float('inf')

        table_data.append({
            "name": name,
            "initial": initial,
            "best": best_val,
            "pct_change": pct_change,
            "bounds": bounds,
            "transform": transform,
            "units": units,
        })

    # Format functions
    def fmt_num(v):
        if isinstance(v, float):
            if abs(v) < 0.001 or abs(v) >= 10000:
                return f"{v:.4g}"
            return f"{v:.4f}"
        return str(v)

    def fmt_bounds(v):
        if isinstance(v, (list, tuple)) and len(v) == 2:
            return f"[{fmt_num(v[0])}, {fmt_num(v[1])}]"
        return str(v)

    def fmt_change(v):
        if not isinstance(v, (int, float)):
            return str(v)
        if abs(v) == float('inf'):
            return '\u221e'
        sign = "+" if v > 0 else ""
        cls = "change-up" if v > 0 else ("change-down" if v < 0 else "change-zero")
        return f'<span class="{cls}">{sign}{v:.1f}%</span>'

    columns = [
        {"key": "name", "label": "Parameter", "align": "left"},
        {"key": "initial", "label": "Initial", "align": "right",
         "format": fmt_num},
        {"key": "best", "label": "Best", "align": "right",
         "format": fmt_num},
        {"key": "pct_change", "label": "Change", "align": "right",
         "format": fmt_change},
        {"key": "bounds", "label": "Bounds", "align": "right",
         "format": fmt_bounds},
        {"key": "transform", "label": "Transform", "align": "left"},
        {"key": "units", "label": "Units", "align": "left"},
    ]

    table_html = rh.build_param_table(table_data, columns, table_id="best-params-table")

    # Arrow indicator summary
    improved = sum(1 for p in table_data if abs(p["pct_change"]) > 1)
    unchanged = len(table_data) - improved

    return f"""
<div class="section" id="best-params">
    <h2>Best Parameter Values</h2>
    <p style="color:var(--text2);margin-bottom:1rem;">
        {improved} parameters changed significantly (&gt;1%),
        {unchanged} remained near initial values.
    </p>
    {table_html}
</div>
"""


def _build_convergence_section(
    results: Dict[str, Any],
    settings: Dict[str, Any],
    output_dir: str
) -> str:
    """Build the convergence plot section."""

    history = results.get("history", [])
    obj_fn = settings.get("objective_function", "KGE")

    if not history:
        return """
<div class="section" id="convergence">
    <h2>Convergence</h2>
    <div class="card"><p>No evaluation history available.</p></div>
</div>
"""

    # Generate convergence plot
    plot_html = ""
    if HAS_MATPLOTLIB and HAS_NUMPY:
        plot_html = _generate_convergence_plot(history, obj_fn, output_dir)

    # Also build an evaluation history summary table
    n_evals = len(history)
    if n_evals > 0:
        # Get objectives
        objectives = [h.get("objective", h.get("best_objective", 0))
                      for h in history]

        best_obj = min(objectives) if objectives else 0
        worst_obj = max(objectives) if objectives else 0
        mean_obj = sum(objectives) / len(objectives) if objectives else 0

        stats_html = f"""
    <div class="card">
        <h3>Objective Function Statistics</h3>
        <div class="perf-card" style="grid-template-columns: repeat(4, 1fr);">
            <div class="perf-metric">
                <div class="perf-value result-good">{best_obj:.4f}</div>
                <div class="perf-label">Best</div>
            </div>
            <div class="perf-metric">
                <div class="perf-value">{mean_obj:.4f}</div>
                <div class="perf-label">Mean</div>
            </div>
            <div class="perf-metric">
                <div class="perf-value result-poor">{worst_obj:.4f}</div>
                <div class="perf-label">Worst</div>
            </div>
            <div class="perf-metric">
                <div class="perf-value">{n_evals}</div>
                <div class="perf-label">Total Evals</div>
            </div>
        </div>
    </div>
"""
    else:
        stats_html = ""

    return f"""
<div class="section" id="convergence">
    <h2>Convergence</h2>
    {plot_html}
    {stats_html}
</div>
"""


def _build_param_evolution_section(
    results: Dict[str, Any],
    parameters: List[Dict],
    output_dir: str
) -> str:
    """Build the parameter evolution plots section."""

    history = results.get("history", [])

    if not history or not parameters:
        return """
<div class="section" id="param-evolution">
    <h2>Parameter Evolution</h2>
    <div class="card"><p>No parameter evolution data available.</p></div>
</div>
"""

    # Generate parameter evolution plots
    plot_html = ""
    if HAS_MATPLOTLIB and HAS_NUMPY:
        plot_html = _generate_param_evolution_plots(
            history, parameters, output_dir
        )

    return f"""
<div class="section" id="param-evolution">
    <h2>Parameter Evolution</h2>
    <p style="color:var(--text2);margin-bottom:1rem;">
        How each parameter changed across calibration evaluations.
    </p>
    {plot_html}
</div>
"""


def _build_correlations_section(
    results: Dict[str, Any],
    parameters: List[Dict],
    output_dir: str
) -> str:
    """Build the parameter correlations heatmap section."""

    history = results.get("history", [])

    if not history or len(parameters) < 2:
        return """
<div class="section" id="param-correlations">
    <h2>Parameter Correlations</h2>
    <div class="card"><p>Insufficient data for correlation analysis.</p></div>
</div>
"""

    plot_html = ""
    if HAS_MATPLOTLIB and HAS_NUMPY:
        plot_html = _generate_correlation_plot(
            history, parameters, output_dir
        )

    return f"""
<div class="section" id="param-correlations">
    <h2>Parameter Correlations</h2>
    <p style="color:var(--text2);margin-bottom:1rem;">
        Pearson correlation between parameter values across all evaluations.
        High correlations may indicate parameter interactions or redundancy.
    </p>
    {plot_html}
</div>
"""


def _build_sensitivity_section(
    sa_results: Dict[str, Any],
    output_dir: str
) -> str:
    """Build the sensitivity analysis results section."""

    method = sa_results.get("method", "unknown").upper()
    param_names = sa_results.get("parameter_names", [])

    if not param_names:
        return """
<div class="section" id="sensitivity">
    <h2>Sensitivity Analysis</h2>
    <div class="card"><p>No sensitivity results available.</p></div>
</div>
"""

    # Build sensitivity table and bar chart
    if method == "MORRIS":
        table_html, plot_html = _build_morris_results(
            sa_results, param_names, output_dir
        )
    elif method == "SOBOL":
        table_html, plot_html = _build_sobol_results(
            sa_results, param_names, output_dir
        )
    else:
        table_html = f'<div class="card"><p>Unknown sensitivity method: {method}</p></div>'
        plot_html = ""

    return f"""
<div class="section" id="sensitivity">
    <h2>Sensitivity Analysis ({method})</h2>
    {plot_html}
    {table_html}
</div>
"""


def _build_performance_section(
    perf_metrics: List[Dict[str, Any]],
    matched_data: Optional[Any],
    settings: Dict[str, Any],
    output_dir: str
) -> str:
    """Build the model performance section."""

    if not perf_metrics:
        return """
<div class="section" id="performance">
    <h2>Model Performance</h2>
    <div class="card"><p>No performance metrics available.</p></div>
</div>
"""

    # Build performance metrics table
    metrics_html = _build_perf_metrics_table(perf_metrics)

    # Build per-species performance cards
    species_cards = _build_species_perf_cards(perf_metrics)

    # Generate performance plots if possible
    plot_html = ""
    if HAS_MATPLOTLIB and HAS_NUMPY and matched_data is not None:
        plot_html = _generate_performance_plots(
            matched_data, perf_metrics, output_dir
        )

    temporal_res = settings.get("temporal_resolution", "native")
    agg_method = settings.get("aggregation_method", "mean")

    return f"""
<div class="section" id="performance">
    <h2>Model Performance</h2>
    <p style="color:var(--text2);margin-bottom:1rem;">
        Comparison at <strong>{temporal_res}</strong> resolution
        (aggregation: {agg_method}).
    </p>
    {species_cards}
    {metrics_html}
    {plot_html}
</div>
"""


def _build_run_best_section(
    results: Dict[str, Any],
    model_config: Dict[str, Any],
    settings: Dict[str, Any]
) -> str:
    """Build the section with code to run the best parameters."""

    best_params = results.get("best_params", {})

    # Format best params as JSON
    best_params_json = json.dumps(best_params, indent=2, default=str)

    apply_code = f"""# Best parameters found by calibration
import json

best_parameters = {best_params_json}

# Save to file
with open("best_parameters.json", "w") as f:
    json.dump(best_parameters, f, indent=2)

# To apply these parameters and run the model:
# 1. These parameters have already been saved to:
#    <calibration_work_dir>/results/best_parameters.json
#
# 2. To re-run with best parameters:
python calibration_config_template.py --apply-best

# 3. Or manually apply to your model config and run:
python model_config_template.py"""

    # Resume / rerun snippets
    rerun_code = """# Continue calibration from checkpoint (more evaluations)
python calibration_config_template.py --resume

# Run sensitivity analysis on final best parameters
python calibration_config_template.py --sensitivity-only

# Generate this report again from saved results
python calibration_config_template.py --report-only"""

    return f"""
<div class="section" id="run-best">
    <h2>Run Best Parameters</h2>
    <div class="card primary">
        <h3>Best Parameter Values</h3>
        {rh.build_code_block(apply_code, "python")}
    </div>
    <div class="card secondary">
        <h3>Additional Commands</h3>
        {rh.build_code_block(rerun_code)}
    </div>
</div>
"""


# =========================================================================
# Plot Generators (require matplotlib)
# =========================================================================

def _fig_to_base64(fig) -> str:
    """Convert a matplotlib figure to base64-encoded PNG data URI."""
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=150, bbox_inches='tight',
                facecolor='#1a1b2e', edgecolor='none')
    buf.seek(0)
    data = rh.embed_image_base64.__wrapped__(buf) if hasattr(rh.embed_image_base64, '__wrapped__') else None

    # Direct encoding
    buf.seek(0)
    import base64
    encoded = base64.b64encode(buf.read()).decode('ascii')
    plt.close(fig)
    return f"data:image/png;base64,{encoded}"


def _setup_plot_style():
    """Configure matplotlib for the report's dark theme."""
    plt.rcParams.update({
        'figure.facecolor': '#1a1b2e',
        'axes.facecolor': '#1e1f2e',
        'axes.edgecolor': '#2a2b3d',
        'axes.labelcolor': '#a0aec0',
        'text.color': '#e2e8f0',
        'xtick.color': '#a0aec0',
        'ytick.color': '#a0aec0',
        'grid.color': '#2d2e42',
        'grid.alpha': 0.5,
        'font.size': 10,
        'font.family': 'sans-serif',
    })


def _generate_convergence_plot(
    history: List[Dict],
    obj_fn: str,
    output_dir: str
) -> str:
    """Generate convergence plot and return HTML with embedded image."""
    _setup_plot_style()

    evals = []
    objectives = []
    best_so_far = []

    current_best = float('inf')
    for h in history:
        eval_id = h.get("eval_id", len(evals))
        obj = h.get("objective", h.get("best_objective", 0))
        evals.append(eval_id)
        objectives.append(obj)

        if obj < current_best:
            current_best = obj
        best_so_far.append(current_best)

    evals = np.array(evals)
    objectives = np.array(objectives)
    best_so_far = np.array(best_so_far)

    fig, ax = plt.subplots(figsize=(10, 5))

    # All evaluations as scatter
    ax.scatter(evals, objectives, s=12, alpha=0.4, color='#4d9ee8',
               label='All evaluations', zorder=2)

    # Best-so-far line
    ax.plot(evals, best_so_far, color='#34d399', linewidth=2.5,
            label='Best so far', zorder=3)

    # Mark the best point
    best_idx = np.argmin(objectives)
    ax.scatter([evals[best_idx]], [objectives[best_idx]],
               s=100, color='#fb923c', marker='*', zorder=4,
               label=f'Best: {objectives[best_idx]:.4f}')

    ax.set_xlabel(f'Evaluation Number')
    ax.set_ylabel(f'Objective ({obj_fn})')
    ax.set_title('Calibration Convergence', fontsize=14, fontweight='bold',
                 color='#e2e8f0')
    ax.legend(facecolor='#1e1f2e', edgecolor='#2a2b3d',
              fontsize=9, loc='upper right')
    ax.grid(True, alpha=0.3)

    data_uri = _fig_to_base64(fig)

    return f"""
    <div class="card">
        <img src="{data_uri}" class="plot-img" alt="Convergence plot"/>
    </div>
    """


def _generate_param_evolution_plots(
    history: List[Dict],
    parameters: List[Dict],
    output_dir: str
) -> str:
    """Generate parameter evolution subplots."""
    _setup_plot_style()

    n_params = len(parameters)
    if n_params == 0:
        return ""

    # Determine grid layout
    n_cols = min(3, n_params)
    n_rows = (n_params + n_cols - 1) // n_cols
    fig_height = max(4, n_rows * 3.2)

    fig, axes = plt.subplots(n_rows, n_cols,
                             figsize=(4 * n_cols, fig_height),
                             squeeze=False)

    for i, param in enumerate(parameters):
        row, col = divmod(i, n_cols)
        ax = axes[row][col]

        pname = param.get("name", f"param_{i}")
        initial = param.get("initial", 0)
        bounds = param.get("bounds", (0, 1))

        # Extract parameter values from history
        eval_ids = []
        values = []
        for h in history:
            params_dict = h.get("parameters", {})
            if pname in params_dict:
                eval_ids.append(h.get("eval_id", len(eval_ids)))
                values.append(params_dict[pname])

        if not values:
            ax.text(0.5, 0.5, 'No data', ha='center', va='center',
                    color='#636e72', transform=ax.transAxes)
            ax.set_title(pname, fontsize=8, color='#a0aec0')
            continue

        eval_ids = np.array(eval_ids)
        values = np.array(values)

        # Plot scatter
        ax.scatter(eval_ids, values, s=8, alpha=0.5, color='#4d9ee8')

        # Mark initial and best
        ax.axhline(y=initial, color='#636e72', linestyle='--',
                   alpha=0.5, linewidth=1, label='Initial')

        # Bounds
        if isinstance(bounds, (list, tuple)) and len(bounds) == 2:
            ax.axhline(y=bounds[0], color='#e74c3c', linestyle=':',
                       alpha=0.3, linewidth=1)
            ax.axhline(y=bounds[1], color='#e74c3c', linestyle=':',
                       alpha=0.3, linewidth=1)

        # Shorten long names
        display_name = pname
        if len(display_name) > 25:
            display_name = display_name[:12] + '..' + display_name[-11:]

        ax.set_title(display_name, fontsize=8, color='#a0aec0')
        ax.grid(True, alpha=0.2)
        ax.tick_params(labelsize=7)

    # Remove unused subplots
    for i in range(n_params, n_rows * n_cols):
        row, col = divmod(i, n_cols)
        axes[row][col].set_visible(False)

    fig.suptitle('Parameter Evolution', fontsize=14, fontweight='bold',
                 color='#e2e8f0', y=1.02)
    fig.tight_layout()

    data_uri = _fig_to_base64(fig)

    return f"""
    <div class="card">
        <img src="{data_uri}" class="plot-img" alt="Parameter evolution plots"/>
    </div>
    """


def _generate_correlation_plot(
    history: List[Dict],
    parameters: List[Dict],
    output_dir: str
) -> str:
    """Generate parameter correlation heatmap."""
    _setup_plot_style()

    param_names = [p.get("name", f"p{i}") for i, p in enumerate(parameters)]

    # Build parameter matrix
    n_evals = len(history)
    n_params = len(param_names)

    if n_evals < 3 or n_params < 2:
        return '<div class="card"><p>Insufficient data for correlation analysis (need &ge;3 evaluations and &ge;2 parameters).</p></div>'

    matrix = np.full((n_evals, n_params), np.nan)
    for i, h in enumerate(history):
        params_dict = h.get("parameters", {})
        for j, pname in enumerate(param_names):
            if pname in params_dict:
                matrix[i, j] = params_dict[pname]

    # Remove columns/rows with all NaN
    valid_cols = ~np.all(np.isnan(matrix), axis=0)
    if np.sum(valid_cols) < 2:
        return '<div class="card"><p>Insufficient valid parameter data for correlation analysis.</p></div>'

    matrix = matrix[:, valid_cols]
    valid_names = [n for n, v in zip(param_names, valid_cols) if v]

    # Compute correlations (handle NaN)
    # Use pairwise complete observations
    n_valid = matrix.shape[1]
    corr = np.eye(n_valid)
    for i in range(n_valid):
        for j in range(i + 1, n_valid):
            mask = ~(np.isnan(matrix[:, i]) | np.isnan(matrix[:, j]))
            if np.sum(mask) > 2:
                r = np.corrcoef(matrix[mask, i], matrix[mask, j])[0, 1]
                corr[i, j] = r
                corr[j, i] = r

    # Shorten names for display
    short_names = []
    for name in valid_names:
        if len(name) > 18:
            short_names.append(name[:8] + '..' + name[-8:])
        else:
            short_names.append(name)

    fig, ax = plt.subplots(figsize=(max(6, n_valid * 0.8),
                                    max(5, n_valid * 0.7)))

    im = ax.imshow(corr, cmap='RdBu_r', vmin=-1, vmax=1, aspect='auto')

    ax.set_xticks(range(n_valid))
    ax.set_yticks(range(n_valid))
    ax.set_xticklabels(short_names, rotation=45, ha='right', fontsize=7)
    ax.set_yticklabels(short_names, fontsize=7)

    # Add correlation values in cells
    for i in range(n_valid):
        for j in range(n_valid):
            val = corr[i, j]
            color = 'white' if abs(val) > 0.5 else '#a0aec0'
            ax.text(j, i, f'{val:.2f}', ha='center', va='center',
                    fontsize=7, color=color)

    cbar = fig.colorbar(im, ax=ax, shrink=0.8)
    cbar.ax.tick_params(labelsize=8)
    cbar.set_label('Pearson Correlation', fontsize=9, color='#a0aec0')

    ax.set_title('Parameter Correlations', fontsize=14, fontweight='bold',
                 color='#e2e8f0')
    fig.tight_layout()

    data_uri = _fig_to_base64(fig)

    return f"""
    <div class="card">
        <img src="{data_uri}" class="plot-img" alt="Parameter correlation heatmap"/>
    </div>
    """


def _generate_performance_plots(
    matched_data,
    perf_metrics: List[Dict],
    output_dir: str
) -> str:
    """Generate observed vs modeled performance plots."""
    _setup_plot_style()

    try:
        import pandas as pd
    except ImportError:
        return '<div class="card"><p>pandas not available — performance plots skipped.</p></div>'

    if not isinstance(matched_data, pd.DataFrame) or matched_data.empty:
        return ""

    plots_html = []

    # Get unique species
    species_list = matched_data['species'].unique() if 'species' in matched_data.columns else []

    for species in species_list:
        sp_data = matched_data[matched_data['species'] == species]

        if sp_data.empty:
            continue

        # Create 2-panel figure: time series + scatter
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

        # Panel 1: Time series
        if 'datetime' in sp_data.columns:
            dates = pd.to_datetime(sp_data['datetime'])
            ax1.plot(dates, sp_data['observed'], 'o', color='#4d9ee8',
                     markersize=4, alpha=0.7, label='Observed')
            ax1.plot(dates, sp_data['simulated'], '-', color='#fb923c',
                     linewidth=1.5, alpha=0.8, label='Simulated')
        else:
            idx = range(len(sp_data))
            ax1.plot(idx, sp_data['observed'], 'o', color='#4d9ee8',
                     markersize=4, alpha=0.7, label='Observed')
            ax1.plot(idx, sp_data['simulated'], '-', color='#fb923c',
                     linewidth=1.5, alpha=0.8, label='Simulated')

        ax1.set_title(f'{species} — Time Series', fontsize=11, color='#e2e8f0')
        ax1.set_ylabel('Concentration')
        ax1.legend(facecolor='#1e1f2e', edgecolor='#2a2b3d', fontsize=8)
        ax1.grid(True, alpha=0.3)

        # Panel 2: Scatter (obs vs sim)
        obs = sp_data['observed'].values
        sim = sp_data['simulated'].values

        ax2.scatter(obs, sim, s=20, alpha=0.5, color='#4d9ee8')

        # 1:1 line
        all_vals = np.concatenate([obs, sim])
        vmin, vmax = np.nanmin(all_vals), np.nanmax(all_vals)
        margin = (vmax - vmin) * 0.05
        ax2.plot([vmin - margin, vmax + margin],
                 [vmin - margin, vmax + margin],
                 '--', color='#636e72', linewidth=1, label='1:1 line')

        ax2.set_xlabel('Observed')
        ax2.set_ylabel('Simulated')
        ax2.set_title(f'{species} — Obs vs Sim', fontsize=11, color='#e2e8f0')
        ax2.legend(facecolor='#1e1f2e', edgecolor='#2a2b3d', fontsize=8)
        ax2.grid(True, alpha=0.3)
        ax2.set_aspect('equal', adjustable='box')

        fig.suptitle(f'Performance: {species}', fontsize=14,
                     fontweight='bold', color='#e2e8f0', y=1.02)
        fig.tight_layout()

        data_uri = _fig_to_base64(fig)

        plots_html.append(f"""
    <details class="module-details" open>
        <summary>{html_lib.escape(species)}</summary>
        <div class="module-content">
            <img src="{data_uri}" class="plot-img"
                 alt="Performance plot for {html_lib.escape(species)}"/>
        </div>
    </details>
        """)

    return "\n".join(plots_html)


# =========================================================================
# Sensitivity Results Helpers
# =========================================================================

def _build_morris_results(
    sa_results: Dict[str, Any],
    param_names: List[str],
    output_dir: str
) -> Tuple[str, str]:
    """Build Morris sensitivity analysis table and plot."""

    mu_star = sa_results.get("mu_star", [])
    sigma = sa_results.get("sigma", [])

    if not mu_star or len(mu_star) != len(param_names):
        return '<div class="card"><p>Incomplete Morris results.</p></div>', ""

    # Sort by mu_star (most influential first)
    indices = list(range(len(param_names)))
    indices.sort(key=lambda i: mu_star[i], reverse=True)

    max_mu = max(mu_star) if mu_star else 1

    # Build table with bar indicators
    rows = []
    for rank, i in enumerate(indices, 1):
        pname = param_names[i]
        mu = mu_star[i]
        sig = sigma[i] if i < len(sigma) else 0

        # Bar width as percentage of max
        bar_pct = (mu / max_mu * 100) if max_mu > 0 else 0

        bar_html = f"""
        <div class="sa-bar-container">
            <div class="sa-bar" style="width:{bar_pct:.1f}%"></div>
        </div>
        """

        # Shorten name
        display_name = pname
        if len(display_name) > 30:
            display_name = display_name[:14] + '..' + display_name[-14:]

        rows.append(f"""
        <tr>
            <td>{rank}</td>
            <td>{html_lib.escape(display_name)}</td>
            <td class="num">{mu:.4g}</td>
            <td class="num">{sig:.4g}</td>
            <td style="min-width:120px;">{bar_html}</td>
        </tr>
        """)

    table_html = f"""
    <div class="card">
        <h3>Morris Screening Results</h3>
        <p style="color:var(--text2);font-size:.85rem;margin-bottom:.8rem;">
            \u03bc* (mu_star) measures overall importance.
            \u03c3 (sigma) measures non-linearity and interactions.
        </p>
        <div class="table-wrap">
        <table class="param-table">
        <thead>
            <tr>
                <th>Rank</th>
                <th>Parameter</th>
                <th>μ*</th>
                <th>σ</th>
                <th>Importance</th>
            </tr>
        </thead>
        <tbody>
        {"".join(rows)}
        </tbody>
        </table>
        </div>
    </div>
    """

    # Generate plot
    plot_html = ""
    if HAS_MATPLOTLIB and HAS_NUMPY:
        plot_html = _generate_morris_plot(
            param_names, mu_star, sigma, output_dir
        )

    return table_html, plot_html


def _build_sobol_results(
    sa_results: Dict[str, Any],
    param_names: List[str],
    output_dir: str
) -> Tuple[str, str]:
    """Build Sobol sensitivity analysis table and plot."""

    S1 = sa_results.get("S1", [])
    ST = sa_results.get("ST", [])

    if not S1 or len(S1) != len(param_names):
        return '<div class="card"><p>Incomplete Sobol results.</p></div>', ""

    # Sort by total-order index (most influential first)
    indices = list(range(len(param_names)))
    indices.sort(key=lambda i: ST[i] if i < len(ST) else S1[i], reverse=True)

    max_st = max(ST) if ST else max(S1) if S1 else 1

    rows = []
    for rank, i in enumerate(indices, 1):
        pname = param_names[i]
        s1 = S1[i]
        st = ST[i] if i < len(ST) else s1

        bar_pct = (st / max_st * 100) if max_st > 0 else 0

        bar_html = f"""
        <div class="sa-bar-container">
            <div class="sa-bar" style="width:{bar_pct:.1f}%"></div>
        </div>
        """

        display_name = pname
        if len(display_name) > 30:
            display_name = display_name[:14] + '..' + display_name[-14:]

        rows.append(f"""
        <tr>
            <td>{rank}</td>
            <td>{html_lib.escape(display_name)}</td>
            <td class="num">{s1:.4f}</td>
            <td class="num">{st:.4f}</td>
            <td style="min-width:120px;">{bar_html}</td>
        </tr>
        """)

    table_html = f"""
    <div class="card">
        <h3>Sobol Sensitivity Indices</h3>
        <p style="color:var(--text2);font-size:.85rem;margin-bottom:.8rem;">
            S1 = first-order (main effect). ST = total-order (including interactions).
        </p>
        <div class="table-wrap">
        <table class="param-table">
        <thead>
            <tr>
                <th>Rank</th>
                <th>Parameter</th>
                <th>S1</th>
                <th>ST</th>
                <th>Importance</th>
            </tr>
        </thead>
        <tbody>
        {"".join(rows)}
        </tbody>
        </table>
        </div>
    </div>
    """

    # Generate plot
    plot_html = ""
    if HAS_MATPLOTLIB and HAS_NUMPY:
        plot_html = _generate_sobol_plot(
            param_names, S1, ST, output_dir
        )

    return table_html, plot_html


def _generate_morris_plot(
    param_names: List[str],
    mu_star: List[float],
    sigma: List[float],
    output_dir: str
) -> str:
    """Generate Morris mu* vs sigma scatter plot."""
    _setup_plot_style()

    fig, ax = plt.subplots(figsize=(8, 6))

    mu_arr = np.array(mu_star)
    sig_arr = np.array(sigma)

    ax.scatter(mu_arr, sig_arr, s=60, color='#4d9ee8', alpha=0.7, zorder=3)

    # Label each point
    for i, name in enumerate(param_names):
        short = name if len(name) <= 15 else name[:7] + '..' + name[-6:]
        ax.annotate(short, (mu_arr[i], sig_arr[i]),
                    fontsize=7, color='#a0aec0',
                    xytext=(5, 5), textcoords='offset points')

    ax.set_xlabel('μ* (Mean Absolute Elementary Effect)')
    ax.set_ylabel('σ (Std of Elementary Effect)')
    ax.set_title('Morris Screening: μ* vs σ', fontsize=14,
                 fontweight='bold', color='#e2e8f0')
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    data_uri = _fig_to_base64(fig)

    return f"""
    <div class="card">
        <img src="{data_uri}" class="plot-img" alt="Morris sensitivity plot"/>
    </div>
    """


def _generate_sobol_plot(
    param_names: List[str],
    S1: List[float],
    ST: List[float],
    output_dir: str
) -> str:
    """Generate Sobol indices bar chart."""
    _setup_plot_style()

    n = len(param_names)
    x = np.arange(n)
    width = 0.35

    # Shorten names
    short_names = []
    for name in param_names:
        if len(name) > 15:
            short_names.append(name[:7] + '..' + name[-6:])
        else:
            short_names.append(name)

    fig, ax = plt.subplots(figsize=(max(8, n * 0.8), 5))

    ax.bar(x - width / 2, S1, width, label='S1 (First-order)',
           color='#4d9ee8', alpha=0.8)
    ax.bar(x + width / 2, ST, width, label='ST (Total-order)',
           color='#34d399', alpha=0.8)

    ax.set_xticks(x)
    ax.set_xticklabels(short_names, rotation=45, ha='right', fontsize=8)
    ax.set_ylabel('Sensitivity Index')
    ax.set_title('Sobol Sensitivity Indices', fontsize=14,
                 fontweight='bold', color='#e2e8f0')
    ax.legend(facecolor='#1e1f2e', edgecolor='#2a2b3d', fontsize=9)
    ax.grid(True, alpha=0.3, axis='y')

    fig.tight_layout()
    data_uri = _fig_to_base64(fig)

    return f"""
    <div class="card">
        <img src="{data_uri}" class="plot-img" alt="Sobol sensitivity indices"/>
    </div>
    """


# =========================================================================
# Performance Helpers
# =========================================================================

def _build_perf_metrics_table(perf_metrics: List[Dict]) -> str:
    """Build performance metrics table for all species."""

    def fmt_metric(v):
        if isinstance(v, (int, float)):
            return f"{v:.4f}"
        return str(v) if v is not None else "N/A"

    def color_kge(v):
        if not isinstance(v, (int, float)):
            return f'<span class="metric-cell">{v}</span>'
        cls = "result-good" if v >= 0.7 else (
            "result-ok" if v >= 0.4 else (
                "result-poor" if v >= 0 else "result-bad"
            )
        )
        return f'<span class="metric-cell {cls}">{v:.4f}</span>'

    def color_nse(v):
        return color_kge(v)  # Same thresholds

    def color_pbias(v):
        if not isinstance(v, (int, float)):
            return f'<span class="metric-cell">{v}</span>'
        absv = abs(v)
        cls = "result-good" if absv <= 10 else (
            "result-ok" if absv <= 25 else (
                "result-poor" if absv <= 50 else "result-bad"
            )
        )
        sign = "+" if v > 0 else ""
        return f'<span class="metric-cell {cls}">{sign}{v:.1f}%</span>'

    columns = [
        {"key": "species", "label": "Species", "align": "left"},
    ]

    # Check if reach_id is present
    if any("reach_id" in m for m in perf_metrics):
        columns.append({"key": "reach_id", "label": "Reach", "align": "left"})

    columns.extend([
        {"key": "KGE", "label": "KGE", "align": "right", "format": color_kge},
        {"key": "NSE", "label": "NSE", "align": "right", "format": color_nse},
        {"key": "RMSE", "label": "RMSE", "align": "right", "format": fmt_metric},
        {"key": "PBIAS", "label": "PBIAS", "align": "right", "format": color_pbias},
    ])

    # Add R2 if available
    if any("R2" in m for m in perf_metrics):
        columns.append({"key": "R2", "label": "R²", "align": "right",
                        "format": fmt_metric})

    return rh.build_param_table(perf_metrics, columns,
                                table_id="perf-metrics-table")


def _build_species_perf_cards(perf_metrics: List[Dict]) -> str:
    """Build per-species performance summary cards."""

    # Group by species
    species_data = {}
    for m in perf_metrics:
        sp = m.get("species", "unknown")
        if sp not in species_data:
            species_data[sp] = []
        species_data[sp].append(m)

    cards = []
    for species, metrics in species_data.items():
        # Average metrics across reaches
        kge_vals = [m.get("KGE") for m in metrics
                    if m.get("KGE") is not None]
        nse_vals = [m.get("NSE") for m in metrics
                    if m.get("NSE") is not None]
        rmse_vals = [m.get("RMSE") for m in metrics
                     if m.get("RMSE") is not None]
        pbias_vals = [m.get("PBIAS") for m in metrics
                      if m.get("PBIAS") is not None]

        avg_kge = sum(kge_vals) / len(kge_vals) if kge_vals else None
        avg_nse = sum(nse_vals) / len(nse_vals) if nse_vals else None
        avg_rmse = sum(rmse_vals) / len(rmse_vals) if rmse_vals else None
        avg_pbias = sum(pbias_vals) / len(pbias_vals) if pbias_vals else None

        # Direct formatting for metric values
        kge_str = f"{avg_kge:.4f}" if avg_kge is not None else "N/A"
        nse_str = f"{avg_nse:.4f}" if avg_nse is not None else "N/A"
        rmse_str = f"{avg_rmse:.4f}" if avg_rmse is not None else "N/A"
        pbias_str = f"{avg_pbias:.1f}%" if avg_pbias is not None else "N/A"

        kge_cls = ""
        if avg_kge is not None:
            kge_cls = ("result-good" if avg_kge >= 0.7 else
                       "result-ok" if avg_kge >= 0.4 else
                       "result-poor" if avg_kge >= 0 else "result-bad")

        nse_cls = ""
        if avg_nse is not None:
            nse_cls = ("result-good" if avg_nse >= 0.7 else
                       "result-ok" if avg_nse >= 0.4 else
                       "result-poor" if avg_nse >= 0 else "result-bad")

        n_reaches = len(metrics)
        reaches_str = f"{n_reaches} reach{'es' if n_reaches > 1 else ''}"

        cards.append(f"""
        <div class="card primary">
            <h3>{html_lib.escape(species)}
                <span class="badge badge-primary" style="margin-left:.5rem;">{reaches_str}</span>
            </h3>
            <div class="perf-card">
                <div class="perf-metric">
                    <div class="perf-value {kge_cls}">{kge_str}</div>
                    <div class="perf-label">KGE</div>
                </div>
                <div class="perf-metric">
                    <div class="perf-value {nse_cls}">{nse_str}</div>
                    <div class="perf-label">NSE</div>
                </div>
                <div class="perf-metric">
                    <div class="perf-value">{rmse_str}</div>
                    <div class="perf-label">RMSE</div>
                </div>
                <div class="perf-metric">
                    <div class="perf-value">{pbias_str}</div>
                    <div class="perf-label">PBIAS</div>
                </div>
            </div>
        </div>
        """)

    return "\n".join(cards)


# =========================================================================
# Objective Quality Assessment
# =========================================================================

def _assess_objective(value: float, obj_fn: str) -> Dict[str, str]:
    """Assess the quality of the best objective function value."""
    import math
    if math.isnan(value) or math.isinf(value):
        return {
            "icon": "\u2753",
            "label": "Unknown",
            "description": "Could not assess — invalid value",
            "style": "warning"
        }

    if obj_fn in ("KGE", "NSE"):
        if value >= 0.75:
            return {"icon": "\U0001f31f", "label": "Very Good",
                    "description": f"{obj_fn} ≥ 0.75 indicates excellent model performance",
                    "style": "success"}
        elif value >= 0.5:
            return {"icon": "\u2705", "label": "Good",
                    "description": f"{obj_fn} ≥ 0.5 indicates satisfactory model performance",
                    "style": "success"}
        elif value >= 0.0:
            return {"icon": "\u26a0\ufe0f", "label": "Acceptable",
                    "description": f"{obj_fn} ≥ 0 but model has room for improvement",
                    "style": "warning"}
        else:
            return {"icon": "\u274c", "label": "Poor",
                    "description": f"Negative {obj_fn} — model performs worse than mean of observations",
                    "style": "warning"}

    elif obj_fn == "RMSE":
        # RMSE is context-dependent, so give generic assessment
        return {"icon": "\U0001f4ca", "label": f"RMSE = {value:.4f}",
                "description": "Lower is better — compare with observation range",
                "style": "info"}

    elif obj_fn == "PBIAS":
        abs_val = abs(value)
        if abs_val <= 10:
            return {"icon": "\U0001f31f", "label": "Very Good",
                    "description": f"|PBIAS| ≤ 10% indicates excellent volume balance",
                    "style": "success"}
        elif abs_val <= 25:
            return {"icon": "\u2705", "label": "Good",
                    "description": f"|PBIAS| ≤ 25% indicates satisfactory volume balance",
                    "style": "success"}
        else:
            return {"icon": "\u26a0\ufe0f", "label": "Needs Improvement",
                    "description": f"|PBIAS| > 25% indicates significant volume bias",
                    "style": "warning"}

    return {"icon": "\U0001f4ca", "label": f"{obj_fn} = {value:.4f}",
            "description": "Custom objective function",
            "style": "info"}
