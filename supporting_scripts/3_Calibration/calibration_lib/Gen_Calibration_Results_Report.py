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
import re
import sys
import json
import html as html_lib
import logging
from pathlib import Path
from datetime import datetime
from typing import List, Dict, Any, Optional, Tuple

import pandas as pd

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
    report_stem: Optional[str] = None,
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
    report_stem : str, optional
        Stem of the originating calibration template (e.g.
        "model_config_template_mizuRoute"). When provided the report is
        named "<report_stem>_results_report.html"; otherwise it falls
        back to "calibration_results_report.html".

    Returns
    -------
    Optional[str]
        Path to the generated HTML report, or None on failure.
    """
    try:
        os.makedirs(output_dir, exist_ok=True)
        _results_report_name = (
            f"{report_stem}_results_report.html"
            if report_stem else "calibration_results_report.html"
        )
        report_path = os.path.join(output_dir, _results_report_name)

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
        ]
        # Observation map nav item — placed right after Summary so it
        # mirrors the section order on the page.
        if _has_observation_map_inputs(model_config):
            nav_items.append({"id": "obs-map", "label": "Observation Map"})
        nav_items.extend([
            {"id": "best-params", "label": "Best Parameters"},
            {"id": "convergence", "label": "Convergence"},
            {"id": "param-evolution", "label": "Parameter Evolution"},
            {"id": "param-correlations", "label": "Correlations"},
        ])
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
            calibration_results, calibration_settings, calibration_parameters,
            performance_metrics=performance_metrics,
        ))

        # ── Section: Observation Map (now right after Summary) ──
        # Same matching as the config-side viewer: stations are coloured
        # by per-station performance (when available); marker size scales
        # with the number of observations.  Skipped silently if no obs
        # CSV or no shapefile is loaded.
        try:
            _map_section = _build_observation_map_section(
                model_config, calibration_settings,
                performance_metrics=performance_metrics,
                matched_data=matched_data,
            )
        except Exception as _e:
            logger.warning(f"Observation map skipped: {_e}")
            _map_section = ""
        if _map_section:
            H.append(_map_section)

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
                calibration_settings, output_dir,
                hostmodel=(model_config.get("hostmodel") or "mizuroute"),
            ))

        # ── Section: Run Best Parameters ──
        H.append(_build_run_best_section(
            calibration_results, model_config, calibration_settings,
            calibration_parameters=calibration_parameters,
        ))

        H.append('</div>')  # container
        H.append(rh.build_footer())
        H.append('</div>')  # main
        H.append('</div>')  # layout

        H.append(rh.build_theme_toggle_js())
        H.append(_build_pdf_export_block(project_name))
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
    parameters: List[Dict],
    performance_metrics: Optional[List[Dict[str, Any]]] = None,
) -> str:
    """Build the results summary section with KPIs.

    Surfaces three views of the objective metric:
      * **Best**   \u2014 the single best score found during optimisation
                     (across all evaluations; what the optimizer
                     converged on).
      * **Mean**   \u2014 arithmetic average of the metric across every
                     ``(species, reach)`` row in ``performance_metrics``
                     (a population-level summary).
      * **Median** \u2014 robust central tendency across the same rows
                     (less sensitive to a single poorly-fitting reach).

    The mean/median tiles only appear when at least one
    ``performance_metrics`` row is available, so the summary stays
    compact for sensitivity-only or no-obs runs.
    """

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

    # Objective value quality assessment (best)
    obj_quality = _assess_objective(best_obj, obj_fn)

    # \u2500\u2500 Mean + median across stations (when performance_metrics avail.) \u2500\u2500
    # For PBIAS the magnitude convention is |PBIAS|, so we average
    # absolute values; for everything else we use the raw metric.
    metric_key = obj_fn.upper()
    USE_ABS = {"PBIAS"}
    vals: List[float] = []
    n_rows = 0
    if performance_metrics:
        for row in performance_metrics:
            v = row.get(metric_key)
            if isinstance(v, (int, float)) and not (isinstance(v, float)
                                                    and v != v):  # NaN check
                vals.append(abs(v) if metric_key in USE_ABS else float(v))
                n_rows += 1
    mean_obj = sum(vals) / len(vals) if vals else None
    median_obj: Optional[float] = None
    if vals:
        s = sorted(vals)
        m = len(s) // 2
        median_obj = s[m] if len(s) % 2 == 1 else 0.5 * (s[m - 1] + s[m])

    # Build KPI grid \u2014 best always, then mean+median when available.
    kpis = [
        {"icon": "\U0001f3c6", "value": f"{best_obj:.4f}",
         "label": f"Best {obj_fn}"},
    ]
    if mean_obj is not None:
        label_prefix = f"|{obj_fn}|" if metric_key in USE_ABS else obj_fn
        kpis.append({
            "icon": "\u2236",   # ratio symbol (\u2236) \u2014 "mean" feel
            "value": f"{mean_obj:.4f}",
            "label": f"Mean {label_prefix} ({n_rows} rows)",
        })
        kpis.append({
            "icon": "\u00f7",   # division sign \u2014 "median" feel
            "value": f"{median_obj:.4f}",
            "label": f"Median {label_prefix}",
        })
    kpis.extend([
        {"icon": "\U0001f504", "value": str(n_evals),
         "label": "Evaluations"},
        {"icon": "\u23f1\ufe0f", "value": runtime_str,
         "label": "Runtime"},
        {"icon": "\U0001f4ca", "value": str(len(parameters)),
         "label": "Parameters"},
    ])

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

    # Quality assessment box (based on the BEST score)
    quality_html = rh.build_highlight_box(
        f"<strong>{obj_quality['icon']} {obj_fn} Quality (best):</strong> "
        f"{obj_quality['label']} ({obj_quality['description']})",
        obj_quality['style']
    )

    # When we have per-station metrics, also evaluate the MEAN against
    # the same quality scale so the user can compare the optimizer's
    # peak fit to the population behaviour at a glance.
    pop_quality_html = ""
    if mean_obj is not None:
        pop = _assess_objective(mean_obj, obj_fn)
        pop_quality_html = rh.build_highlight_box(
            f"<strong>{pop['icon']} {obj_fn} Quality (mean across "
            f"{n_rows} {('|' + obj_fn + '|') if metric_key in USE_ABS else obj_fn}"
            f" rows):</strong> {pop['label']} "
            f"({pop['description']})",
            pop['style']
        )

    return f"""
<div class="section" id="summary">
    <h2>Results Summary</h2>
    {kpi_html}
    {conv_html}
    {quality_html}
    {pop_quality_html}
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
    output_dir: str,
    hostmodel: str = "mizuroute",
) -> str:
    """Build the model performance section.

    ``hostmodel`` drives the per-feature terminology and the explanatory
    note above the plots: ``"summa"`` shows HRU outlet observations
    (the pouring-point data that drives ``ObjectiveFunction`` when
    ``use_primary_only=True``); anything else shows river-reach
    observations.
    """

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
            matched_data, perf_metrics, output_dir,
            hostmodel=hostmodel,
        )

    temporal_res = settings.get("temporal_resolution", "native")
    agg_method = settings.get("aggregation_method", "mean")

    # Host-aware note: SUMMA reports outlet-only observations per HRU
    # because that's what the objective function uses; mizuRoute reports
    # every observation matched to a river reach.
    _is_summa = (hostmodel or "").lower() == "summa"
    if _is_summa:
        host_note = (
            '<div class="highlight-box info" style="margin-bottom:1rem">'
            '<strong>SUMMA:</strong> performance is reported at the '
            '<strong>outlet of each HRU</strong> &mdash; the pouring-point '
            'observation that drives the calibration objective function '
            '(<code>use_primary_only=True</code>).  Other observations that '
            'fall inside an HRU but away from its outlet are not shown here.'
            '</div>'
        )
    else:
        host_note = (
            '<div class="highlight-box info" style="margin-bottom:1rem">'
            '<strong>mizuRoute:</strong> performance is reported per '
            '<strong>river reach</strong>, showing every observation matched '
            'to that reach.'
            '</div>'
        )

    return f"""
<div class="section" id="performance">
    <h2>Model Performance</h2>
    <p style="color:var(--text2);margin-bottom:1rem;">
        Comparison at <strong>{temporal_res}</strong> resolution
        (aggregation: {agg_method}).
    </p>
    {host_note}
    {species_cards}
    {metrics_html}
    {plot_html}
</div>
"""


def _build_docker_run_snippet(container_cfg: Dict[str, Any],
                              spatial: Dict[str, Any]) -> str:
    """Build the same kind of ``docker exec ... mpirun ...`` line the
    config-side report shows, but with the post-``--apply-best`` config in
    mind.  Falls back to a placeholder with ``<container_path_to:...>``
    when we cannot resolve the docker volume mount."""
    exe = container_cfg.get("executable_path") or ""
    fm = container_cfg.get("file_manager_path") or ""
    name = container_cfg.get("docker_container_name") or "docker_openwq"
    dc_path = container_cfg.get("docker_compose_path") or ""
    hostmodel = (spatial.get("hostmodel") or "mizuroute").lower()
    mpi_np = container_cfg.get("mpi_np", 2)
    _np = 1 if hostmodel == "summa" else mpi_np
    _fm_flag = "-m " if hostmodel == "summa" else ""

    # Try to resolve container paths from docker-compose.yml so the user
    # gets a fully literal command they can paste.  Use the same helper the
    # config report uses (Gen_Input_Driver._parse_docker_volume_mount /
    # _correct_path_for_docker) — but tolerate it not being importable
    # (the calibration runs may not have sys.path set the same way).
    _cont_exe = _cont_fm = _cont_wd = None
    if exe and fm and dc_path and os.path.isfile(dc_path):
        try:
            # Locate the config-side support library
            here = Path(__file__).resolve().parent
            candidates = [
                here.parent.parent / "1_Model_Config" / "config_support_lib",
                here.parent.parent.parent / "1_Model_Config" / "config_support_lib",
            ]
            for c in candidates:
                if (c / "Gen_Input_Driver.py").is_file():
                    if str(c) not in sys.path:
                        sys.path.insert(0, str(c))
                    break
            import Gen_Input_Driver as _gid  # type: ignore
            _host_root, _cont_root = _gid._parse_docker_volume_mount(dc_path)
            if _host_root and _cont_root:
                _cont_exe = _gid._correct_path_for_docker(
                    os.path.abspath(exe), _host_root, _cont_root)
                _cont_fm = _gid._correct_path_for_docker(
                    os.path.abspath(fm), _host_root, _cont_root)
                _wd = os.path.dirname(os.path.abspath(exe))
                _cont_wd = _gid._correct_path_for_docker(
                    _wd + "/", _host_root, _cont_root).rstrip("/")
        except Exception:
            pass

    if _cont_exe and _cont_fm and _cont_wd:
        return (f'docker exec {name} /bin/bash -c '
                f'"export HDF5_USE_FILE_LOCKING=FALSE && '
                f'cd {_cont_wd} && '
                f'mpirun --allow-run-as-root -np {_np} '
                f'{_cont_exe} {_fm_flag}{_cont_fm}"')

    return (f'docker exec {name} /bin/bash -c '
            f'"export HDF5_USE_FILE_LOCKING=FALSE && '
            f'cd <container_path_to:{os.path.dirname(os.path.abspath(exe)) if exe else "<exec_dir>"}> && '
            f'mpirun --allow-run-as-root -np {_np} '
            f'<container_path_to:{os.path.abspath(exe) if exe else "<executable>"}> '
            f'{_fm_flag}<container_path_to:{os.path.abspath(fm) if fm else "<file_manager>"}>"')


def _resolve_hpc_settings(container_cfg: Dict[str, Any],
                          settings: Dict[str, Any],
                          model_config: Dict[str, Any]
                          ) -> Tuple[str, str]:
    """Resolve the Apptainer SIF + bind path from any of the places the
    user might have set them, falling back to sensible placeholders.

    Lookup order (most-specific wins):
      1. ``calibration_settings`` — produced by the interactive setup
         report and forwarded into the generated run script.
      2. ``container_cfg`` — passed through from ``model_config`` via
         ``config_integration.get_container_config``.
      3. ``model_config`` — direct read for clones / older runs.
      4. A literal ``openwq.sif`` + the parent-of-executable directory.
    """
    sif = (settings.get("apptainer_sif_path")
           or container_cfg.get("apptainer_sif_path")
           or model_config.get("apptainer_sif_path")
           or "openwq.sif")
    bind = (settings.get("apptainer_bind_path")
            or container_cfg.get("apptainer_bind_path")
            or model_config.get("apptainer_bind_path"))
    if not bind:
        exe = container_cfg.get("executable_path") or ""
        if exe and "<" not in exe:
            bind = os.path.dirname(os.path.dirname(os.path.abspath(exe)))
        else:
            bind = "/scratch/$USER"
    return sif, bind


def _build_apptainer_snippet(container_cfg: Dict[str, Any],
                             spatial: Dict[str, Any],
                             settings: Optional[Dict[str, Any]] = None,
                             model_config: Optional[Dict[str, Any]] = None
                             ) -> str:
    """Build a single-run Apptainer command (HPC).  Format mirrors the one
    in ``hpc_templates/slurm_array_template.sh`` so users who copy this
    from the report can paste it into their own batch script.

    SIF + bind path are resolved through ``_resolve_hpc_settings`` so the
    snippet shows the user's configured paths when they exist (avoids the
    user having to hand-edit the ``openwq.sif`` placeholder every time).
    """
    settings = settings or {}
    model_config = model_config or {}
    exe_path = container_cfg.get("executable_path") or "<path/to/executable>"
    fm_path = container_cfg.get("file_manager_path") or "<path/to/fileManager.txt>"
    hostmodel = (spatial.get("hostmodel") or "mizuroute").lower()
    sif, bind = _resolve_hpc_settings(container_cfg, settings, model_config)
    fm_flag = "-m " if hostmodel == "summa" else ""
    return (
        '# 1) Module load (adjust for your cluster)\n'
        '# module load apptainer\n'
        '\n'
        '# 2) Set paths\n'
        f'SIF_PATH={sif}                       # your openwq Apptainer image\n'
        f'BIND_PATH={bind}        # bind-mount root (must contain config + outputs)\n'
        f'EXECUTABLE={os.path.basename(exe_path)}\n'
        f'FILE_MANAGER={fm_path}\n'
        '\n'
        '# 3) Run\n'
        'cd "$(dirname "$EXECUTABLE")"\n'
        'apptainer exec \\\n'
        '    --bind "$BIND_PATH" \\\n'
        '    --env master_json="$PWD/openWQ_master.json" \\\n'
        '    "$SIF_PATH" \\\n'
        f'    ./"$EXECUTABLE" {fm_flag}"$FILE_MANAGER"'
    )


def _build_slurm_snippet(container_cfg: Dict[str, Any],
                         spatial: Dict[str, Any],
                         settings: Optional[Dict[str, Any]] = None,
                         model_config: Optional[Dict[str, Any]] = None
                         ) -> str:
    """Minimal ready-to-submit SLURM batch script wrapping the single-run
    Apptainer call.  Useful for users who just want to verify the
    best-parameter run on a worker node before launching a calibration
    array on top."""
    settings = settings or {}
    model_config = model_config or {}
    hostmodel = (spatial.get("hostmodel") or "mizuroute").lower()
    fm_flag = "-m " if hostmodel == "summa" else ""
    fm_path = container_cfg.get("file_manager_path") or "<path/to/fileManager.txt>"
    exe_path = container_cfg.get("executable_path") or "<path/to/executable>"
    sif, bind = _resolve_hpc_settings(container_cfg, settings, model_config)
    # SLURM resource hints — use the calibration's mpi_np / hostmodel
    # rather than hard-coding 4 tasks.  SUMMA-OpenWQ is serial.
    _np = 1 if hostmodel == "summa" else container_cfg.get("mpi_np", 2)
    return (
        '#!/bin/bash\n'
        '#SBATCH --job-name=openwq_best\n'
        '#SBATCH --time=02:00:00\n'
        '#SBATCH --nodes=1\n'
        f'#SBATCH --ntasks-per-node={_np}\n'
        '#SBATCH --mem=16G\n'
        '#SBATCH --output=openwq_best_%j.out\n'
        '\n'
        '# module load apptainer   # adjust for your cluster\n'
        f'SIF_PATH={sif}\n'
        f'BIND_PATH={bind}\n'
        f'EXECUTABLE={os.path.basename(exe_path)}\n'
        f'FILE_MANAGER={fm_path}\n'
        '\n'
        'cd "$(dirname "$EXECUTABLE")"\n'
        'apptainer exec \\\n'
        '    --bind "$BIND_PATH" \\\n'
        '    --env master_json="$PWD/openWQ_master.json" \\\n'
        '    "$SIF_PATH" \\\n'
        f'    ./"$EXECUTABLE" {fm_flag}"$FILE_MANAGER"'
    )


def _build_visualization_snippet(model_config: Dict[str, Any],
                                 spatial: Dict[str, Any]) -> str:
    """Generate a self-contained Python snippet that reads the best-params
    HDF5 output and re-uses the SAME ``Plot_h5_driver`` the config viewer
    uses — so the plots match what the user already validated."""
    hostmodel = (spatial.get("hostmodel") or "mizuroute").lower()
    feature_label = spatial.get("feature_label") or spatial.get(
        "river_network_mapping_key") or "SegId"
    river_shp = model_config.get("river_network_shapefile") or ""
    # Resolve basin shapefile (lives in ss_method_copernicus_basin_info)
    basin_info = model_config.get("ss_method_copernicus_basin_info") or {}
    basin_shp = basin_info.get("path_to_shp") if isinstance(basin_info, dict) else ""
    species = model_config.get("chemical_species") or []
    # User's openwq_out directory — derived from executable_path
    exe = model_config.get("executable_path") or ""
    out_dir = (os.path.dirname(os.path.abspath(exe))
               if exe else "<dir2save_input_files>")
    species_repr = json.dumps(list(species))

    body = (
        "import sys, os\n"
        "# Adjust if your 2_Read_Outputs path differs:\n"
        f"sys.path.insert(0, os.path.abspath('../2_Read_Outputs/hdf5_support_lib'))\n"
        "import Read_h5_driver as h5_reader\n"
        "import Plot_h5_driver as h5_plib\n"
        "\n"
        "openwq_info = {\n"
        f"    'path_to_results': r'{out_dir}',\n"
        f"    'mapping_key': '{spatial.get('h5_mapping_key') or ('hruId' if hostmodel == 'summa' else 'reachID')}',\n"
        "}\n"
        "\n"
        "openwq_results = h5_reader.Read_h5_driver(\n"
        "    openwq_info=openwq_info,\n"
        "    output_format='HDF5',\n"
        "    debugmode=True,\n"
        f"    chemSpec={species_repr},\n"
        "    space_elem='all',\n"
        "    chemUnits='MG/L',\n"
        ")\n"
        "\n"
        "# Reproduce the same interactive viewer the config report links to\n"
        "h5_plib.Plot_h5_driver(\n"
        "    what2map='openwq',\n"
        f"    hostmodel='{hostmodel}',\n"
        "    mapping_key_values='all',\n"
        "    openwq_results=openwq_results,\n"
        f"    chemSpec={species_repr},\n"
        "    debugmode=True,\n"
        f"    output_path=r'{os.path.join(out_dir, 'openwq_best_params_report.html')}',\n"
        f"    river_network_shp=r'{river_shp}',\n"
        f"    basin_shapefile=r'{basin_shp}',\n"
        f"    mapping_key='{feature_label}',\n"
        f"    feature_label='{feature_label}',\n"
        ")\n"
    )
    return body


def _build_run_best_section(
    results: Dict[str, Any],
    model_config: Dict[str, Any],
    settings: Dict[str, Any],
    calibration_parameters: Optional[List[Dict[str, Any]]] = None,
) -> str:
    """One copy-paste block: the user's ``model_config_template.py``
    script with the calibrated setup.

    After the user runs ``python calibration_config_template.py
    --apply-best``, the BGC JSON file referenced by the model template
    contains the best parameter values.  So the same Python template
    re-run produces the calibrated model.  We surface it here as a
    single code block the user can copy verbatim into a new file.
    """
    template_path = model_config.get("_model_config_path", "") or ""

    if not template_path or not os.path.isfile(template_path):
        # Defensive fallback — couldn't find the template on disk.
        best_params = results.get("best_params") or {}
        return f"""
<div class="section" id="run-best">
    <h2>Run with Best Parameters</h2>
    <div class="card">
        <p>Could not locate <code>model_config_template.py</code> on
           disk.  Best parameter values are shown below — re-apply
           them manually to your model setup.</p>
        {rh.build_code_block(json.dumps(best_params, indent=2, default=str), "json")}
    </div>
</div>
"""

    try:
        with open(template_path, "r", encoding="utf-8", errors="replace") as f:
            template_text = f.read()
    except Exception as _e:
        logger.warning(f"Run-best section: failed to read {template_path}: {_e}")
        template_text = ""

    template_filename = os.path.basename(template_path)
    sub = (
        '<p style="margin-top:.2rem">'
        'Run-ready model setup. After applying the calibrated parameters '
        'with <code>python calibration_config_template.py --apply-best</code>, '
        'the BGC JSON pointed to by this template contains the best '
        'values — so re-running this script produces the calibrated '
        f'model.  Source: <code>{html_lib.escape(template_filename)}</code>.'
        '</p>'
    )

    return f"""
<div class="section" id="run-best">
    <h2>Run with Best Parameters
        <span style="font-weight:400;font-size:.7em;color:var(--muted,#888)">
        (configuration template python script with the best parameters)
        </span>
    </h2>
    <div class="card primary">
        {sub}
        {rh.build_code_block(template_text, "python")}
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
    output_dir: str,
    max_individual_plots: int = 60,
    hostmodel: str = "mizuroute",
) -> str:
    """Per-species performance figures.

    Host-aware labelling:
      * ``mizuroute`` — each spatial element is a **river reach**.
        The per-reach plot shows every observation matched to that
        reach (the standard "show results points within the river
        reach" view).
      * ``summa`` — each spatial element is an **HRU**.  Because
        ``ObjectiveFunction`` runs with ``use_primary_only=True`` by
        default, ``matched_data`` only carries the pouring-point
        observation per HRU — so what's drawn is literally "the
        outlet observation that drives the calibration objective".
        The plot titles + a banner above the species block call this
        out explicitly so the reader doesn't mistake the single
        outlet trace for "all observations inside the HRU".

    For each chemical species we emit:

      1. **One overview figure** — every spatial element drawn as a
         distinct colour on the same axes (time series + obs-vs-sim
         scatter), so the user can see the network-level distribution
         at a glance.

      2. **One titled figure per spatial element** — collapsed into
         ``<details>`` elements (closed by default) so the page stays
         scannable; each carries its own ``<summary>`` heading
         (``"HRU 740457190 — outlet observation"`` or
         ``"Reach 740457190"``) and a two-panel chart focused on that
         single element.

    The individual plots are skipped when there are more than
    ``max_individual_plots`` elements per species (default 60).
    """
    # Host-aware terminology used everywhere below.
    _is_summa = (hostmodel or "").lower() == "summa"
    feat_label = "HRU" if _is_summa else "Reach"
    feat_label_plural = "HRUs" if _is_summa else "reaches"
    # For SUMMA we make it explicit that each plot is the OUTLET point
    # (per pouring-point matching from prepare_calibration_observations).
    outlet_suffix = " — outlet observation" if _is_summa else ""
    _setup_plot_style()

    try:
        import pandas as pd
    except ImportError:
        return '<div class="card"><p>pandas not available — performance plots skipped.</p></div>'

    if not isinstance(matched_data, pd.DataFrame) or matched_data.empty:
        return ""

    plots_html = []

    species_list = (matched_data['species'].unique()
                    if 'species' in matched_data.columns else [])
    has_reach = 'reach_id' in matched_data.columns

    for species in species_list:
        sp_data = matched_data[matched_data['species'] == species]
        if sp_data.empty:
            continue

        # Collect reach IDs for this species (stable order so a reach
        # gets the same colour in every figure).
        if has_reach:
            reach_ids = sorted(sp_data['reach_id'].dropna().unique(),
                               key=lambda x: (str(type(x).__name__), x))
        else:
            reach_ids = [None]
        n_reaches = len(reach_ids)

        # Pick a colour map for the reaches.  ``tab20`` gives 20 visually
        # distinct hues; for larger networks we fall through to
        # ``viridis`` which stays readable at 100+ traces.
        import matplotlib.cm as _cm
        if n_reaches <= 20:
            cmap = _cm.get_cmap('tab20', max(n_reaches, 2))
            show_legend = True
        else:
            cmap = _cm.get_cmap('viridis', n_reaches)
            show_legend = False

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

        # Track the global min/max so the 1:1 line spans the data.
        all_obs, all_sim = [], []

        for i, rid in enumerate(reach_ids):
            if rid is None:
                rd = sp_data
            else:
                rd = sp_data[sp_data['reach_id'] == rid]
            if rd.empty:
                continue
            colour = cmap(i)
            label = (f'{feat_label} {rid}' if rid is not None else 'all')

            # ── Time series (panel 1) ──
            if 'datetime' in rd.columns:
                _x = pd.to_datetime(rd['datetime'])
            else:
                _x = range(len(rd))
            # Observed: small open circles
            ax1.plot(_x, rd['observed'], 'o', color=colour, markersize=3.5,
                     markerfacecolor='none', markeredgewidth=0.8,
                     alpha=0.75, label=f'obs · {label}')
            # Simulated: thin line in the same colour
            ax1.plot(_x, rd['simulated'], '-', color=colour, linewidth=1.0,
                     alpha=0.85, label=f'sim · {label}')

            # ── Scatter (panel 2) ──
            ax2.scatter(rd['observed'], rd['simulated'],
                        s=18, alpha=0.55, color=colour,
                        edgecolors='none', label=label)

            all_obs.extend(rd['observed'].dropna().tolist())
            all_sim.extend(rd['simulated'].dropna().tolist())

        ax1.set_title(
            f'{species} — time series ({n_reaches} {feat_label_plural}'
            + (', outlet observations' if _is_summa else '') + ')',
            fontsize=11, color='#e2e8f0')
        ax1.set_ylabel('Concentration')
        ax1.grid(True, alpha=0.3)

        # ── 1:1 line on the scatter panel ──
        if all_obs and all_sim:
            v = np.array(all_obs + all_sim)
            vmin = float(np.nanmin(v))
            vmax = float(np.nanmax(v))
            margin = (vmax - vmin) * 0.05 if vmax > vmin else 0.05
            ax2.plot([vmin - margin, vmax + margin],
                     [vmin - margin, vmax + margin],
                     '--', color='#9aa5b1', linewidth=1, label='1:1 line')

        ax2.set_xlabel('Observed')
        ax2.set_ylabel('Simulated')
        ax2.set_title(f'{species} — observed vs simulated',
                      fontsize=11, color='#e2e8f0')
        ax2.grid(True, alpha=0.3)
        ax2.set_aspect('equal', adjustable='box')

        # ── Legend vs colorbar ──
        if show_legend and n_reaches > 1:
            # Compact per-reach legend, deduplicated (one entry per reach).
            handles, labels = ax2.get_legend_handles_labels()
            # tab20 distinct labels per reach + the 1:1 entry; cap at 20.
            ax2.legend(handles, labels, loc='upper left',
                       facecolor='#1e1f2e', edgecolor='#2a2b3d',
                       fontsize=7, ncol=1,
                       bbox_to_anchor=(1.02, 1), borderaxespad=0)
        elif n_reaches > 1:
            # Too many reaches for an inline legend — add a colorbar
            # bridging the viridis ramp and the index range.
            import matplotlib as mpl
            sm = mpl.cm.ScalarMappable(
                cmap=cmap, norm=mpl.colors.Normalize(vmin=0, vmax=n_reaches - 1))
            sm.set_array([])
            cbar = fig.colorbar(sm, ax=ax2, fraction=0.04, pad=0.04)
            cbar.set_label(f'{feat_label} index (0 … {n_reaches - 1})',
                           color='#a0aec0', fontsize=8)
            cbar.ax.tick_params(labelsize=7)
        else:
            # Single reach — just show a basic legend.
            ax1.legend(facecolor='#1e1f2e', edgecolor='#2a2b3d', fontsize=8)
            ax2.legend(facecolor='#1e1f2e', edgecolor='#2a2b3d', fontsize=8)

        fig.suptitle(f'Performance: {species}', fontsize=14,
                     fontweight='bold', color='#e2e8f0', y=1.02)
        fig.tight_layout()

        overview_uri = _fig_to_base64(fig)

        # ── Per-reach individual figures ───────────────────────────
        # Each gets its own 2-panel chart + its own clearly titled
        # <details> block so the user can pin down behaviour at a
        # specific reach without squinting at the overview.  Skipped
        # past the cap to keep the report a manageable size.
        per_reach_blocks: List[str] = []
        if has_reach and 0 < n_reaches <= max_individual_plots:
            for i, rid in enumerate(reach_ids):
                rd = sp_data[sp_data['reach_id'] == rid].sort_values('datetime')
                if rd.empty:
                    continue
                colour = cmap(i)
                rfig, (rax1, rax2) = plt.subplots(1, 2, figsize=(13, 4.2))
                # Panel 1: time series for THIS reach
                if 'datetime' in rd.columns:
                    rx = pd.to_datetime(rd['datetime'])
                else:
                    rx = range(len(rd))
                rax1.plot(rx, rd['observed'], 'o', color=colour, markersize=5,
                          markerfacecolor='none', markeredgewidth=1.0,
                          alpha=0.85, label='Observed')
                rax1.plot(rx, rd['simulated'], '-', color=colour, linewidth=1.6,
                          alpha=0.95, label='Simulated')
                rax1.set_title(
                    f'{species} — {feat_label} {rid}{outlet_suffix} — time series',
                    fontsize=11, color='#e2e8f0')
                rax1.set_ylabel('Concentration')
                rax1.grid(True, alpha=0.3)
                rax1.legend(facecolor='#1e1f2e', edgecolor='#2a2b3d', fontsize=8)

                # Panel 2: obs-vs-sim scatter for THIS reach
                ro = rd['observed'].values
                rs = rd['simulated'].values
                rax2.scatter(ro, rs, s=22, alpha=0.7, color=colour,
                             edgecolors='none')
                rv = np.array(list(ro) + list(rs))
                if rv.size:
                    rmin = float(np.nanmin(rv)); rmax = float(np.nanmax(rv))
                    rmar = (rmax - rmin) * 0.05 if rmax > rmin else 0.05
                    rax2.plot([rmin - rmar, rmax + rmar],
                              [rmin - rmar, rmax + rmar],
                              '--', color='#9aa5b1', linewidth=1, label='1:1 line')
                rax2.set_xlabel('Observed')
                rax2.set_ylabel('Simulated')
                rax2.set_title(
                    f'{species} — {feat_label} {rid}{outlet_suffix} — observed vs simulated',
                    fontsize=11, color='#e2e8f0')
                rax2.set_aspect('equal', adjustable='box')
                rax2.grid(True, alpha=0.3)
                rax2.legend(facecolor='#1e1f2e', edgecolor='#2a2b3d', fontsize=8)

                # Header: include the per-element KGE etc. if known
                _hdr = f'Performance: {species} — {feat_label} {rid}{outlet_suffix}'
                if perf_metrics:
                    for r in perf_metrics:
                        if (str(r.get('reach_id', '')) == str(rid)
                                and r.get('species') == species):
                            _bits = []
                            for k in ('KGE', 'NSE', 'RMSE', 'PBIAS', 'R2'):
                                v = r.get(k)
                                if isinstance(v, (int, float)):
                                    _bits.append(f'{k}={v:.3f}')
                            if _bits:
                                _hdr += '  (' + ', '.join(_bits) + ')'
                            break
                rfig.suptitle(_hdr, fontsize=13, fontweight='bold',
                              color='#e2e8f0', y=1.02)
                rfig.tight_layout()
                per_reach_uri = _fig_to_base64(rfig)
                per_reach_blocks.append(f"""
        <details class="module-details" style="margin:8px 0">
            <summary>{html_lib.escape(feat_label)} {html_lib.escape(str(rid))}{html_lib.escape(outlet_suffix)}</summary>
            <div class="module-content">
                <img src="{per_reach_uri}" class="plot-img"
                     alt="Performance plot for {html_lib.escape(species)} {html_lib.escape(feat_label)} {html_lib.escape(str(rid))}"/>
            </div>
        </details>
                """)

        # If too many elements for individual plots, surface that fact:
        skipped_notice = ''
        if has_reach and n_reaches > max_individual_plots:
            skipped_notice = (
                f'<p style="margin:.5rem 0;font-size:.85rem;color:var(--muted,#888)">'
                f'<i>Per-{feat_label.lower()} individual plots suppressed for '
                f'{n_reaches} {feat_label_plural} (cap: {max_individual_plots}).  '
                f'Overview above shows all {feat_label_plural} colour-coded.</i></p>'
            )

        _count_label = f"{n_reaches} {feat_label if n_reaches == 1 else feat_label_plural}"
        _overview_hdr = (f"All {feat_label_plural} overview"
                         + (" (outlet observations)" if _is_summa else ""))
        _detail_hdr = f"Per-{feat_label.lower()} detail"
        plots_html.append(f"""
    <details class="module-details" open>
        <summary>{html_lib.escape(species)} &mdash; {_count_label}</summary>
        <div class="module-content">
            <h3 style="margin:.4rem 0 .3rem;font-size:1.05rem">{_overview_hdr}</h3>
            <img src="{overview_uri}" class="plot-img"
                 alt="Overview performance plot for {html_lib.escape(species)}"/>
            {skipped_notice}
            {f'<h3 style="margin:1rem 0 .3rem;font-size:1.05rem">{_detail_hdr}</h3>' if per_reach_blocks else ''}
            {''.join(per_reach_blocks)}
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


# =========================================================================
# Observation Map Section — shows where each obs station was matched to
# (reach for mizuRoute, basin/HRU for SUMMA), using the SAME helpers as
# the config-side viewer (`Plot_h5_driver._match_stations_*`).
# =========================================================================

def _has_observation_map_inputs(model_config: Dict[str, Any]) -> bool:
    """Cheap sanity check used by the sidebar — do we have *any* of the
    pieces needed to render the obs map?  Avoids adding a dangling nav
    item when the section is going to come back empty anyway."""
    try:
        obs_cfg = _ci.get_observation_config(model_config)
    except Exception:
        return False
    has_obs = (obs_cfg.get("source") in ("grqa", "user_csv")
               or bool(obs_cfg.get("user_observation_csv"))
               or bool(obs_cfg.get("grqa_local_data_path")))
    has_shp = bool(obs_cfg.get("river_network_shapefile")) \
        or bool(obs_cfg.get("basin_shapefile"))
    return has_obs and has_shp


def _load_shapefile_as_geojson(shp_path: str):
    """Read a shapefile, reproject to WGS84, return ``(geojson_dict, bounds)``
    or ``(None, None)`` on failure.  Keeps the calibration report
    independent of geopandas-only callers — uses fiona+shapely directly."""
    try:
        import fiona
        from shapely.geometry import shape, mapping as shp_mapping
    except ImportError:
        return None, None
    try:
        feats = []
        src_proj4 = None
        with fiona.open(shp_path) as src:
            if src.crs:
                if isinstance(src.crs, dict):
                    src_proj4 = ' '.join(
                        f'+{k}' if v is True else f'+{k}={v}'
                        for k, v in src.crs.items())
                else:
                    try:
                        src_proj4 = src.crs.to_proj4()
                    except Exception:
                        src_proj4 = None
            for f in src:
                geom = shape(f['geometry'])
                feats.append({
                    'type': 'Feature',
                    'geometry': geom,
                    'properties': dict(f['properties']),
                })
        if not feats:
            return None, None
        sb = feats[0]['geometry'].bounds
        if (abs(sb[0]) > 360 or abs(sb[1]) > 360) and src_proj4:
            try:
                from pyproj import Proj
                from shapely.ops import transform as shp_transform
                p = Proj(src_proj4)
                _to_wgs = lambda x, y: p(x, y, inverse=True)  # noqa: E731
                for f in feats:
                    f['geometry'] = shp_transform(_to_wgs, f['geometry'])
            except Exception:
                pass
        gj_feats = [{
            'type': 'Feature',
            'geometry': shp_mapping(f['geometry']),
            'properties': {k: (str(v) if v is not None else '')
                           for k, v in f['properties'].items()},
        } for f in feats]
        all_b = [f['geometry'].bounds for f in feats]
        bounds = (min(b[0] for b in all_b), min(b[1] for b in all_b),
                  max(b[2] for b in all_b), max(b[3] for b in all_b))
        return ({'type': 'FeatureCollection', 'features': gj_feats}, bounds)
    except Exception:
        return None, None


def _import_spatial_matching():
    """Locate ``spatial_matching.py`` in 2_Read_Outputs/hdf5_support_lib and
    return the module.  Tries the conventional repo layout, then falls
    back to the env-var-based ``OPENWQ_H5_SUPPORT_LIB`` if set."""
    import importlib
    here = Path(__file__).resolve().parent
    candidates = [
        # century_basins_tested layout
        here.parent.parent / "2_Read_Outputs" / "hdf5_support_lib",
        # other clones' layout (just in case)
        here.parent.parent.parent / "2_Read_Outputs" / "hdf5_support_lib",
    ]
    env = os.environ.get("OPENWQ_H5_SUPPORT_LIB")
    if env:
        candidates.insert(0, Path(env))
    for c in candidates:
        if c and (c / "spatial_matching.py").is_file():
            if str(c) not in sys.path:
                sys.path.insert(0, str(c))
            return importlib.import_module("spatial_matching")
    raise ImportError(
        "Could not locate spatial_matching.py in any of: "
        + ", ".join(str(c) for c in candidates)
    )


# Colour palette identical to the config-side viewer for visual consistency.
_RESULTS_MAP_COLORWAY = [
    '#0066cc', '#00a86b', '#ff6b35', '#004499', '#34d399', '#fb923c',
    '#667eea', '#764ba2', '#e63946', '#2ec4b6', '#e9c46a', '#264653',
]


def _build_observation_map_section(
        model_config: Dict[str, Any],
        calibration_settings: Dict[str, Any],
        performance_metrics: Optional[List[Dict[str, Any]]] = None,
        matched_data: Optional[Any] = None,
) -> str:
    """Render a Leaflet map of basins / reaches + observation stations.

    Markers carry CALIBRATION-AWARE styling:
      * radius  ∝ √n_obs   (so a station with 10× more data renders
                            ~3× larger; capped to keep the map readable)
      * fill    = colormap(objective_metric)
                  RdYlGn for higher-is-better metrics (KGE / NSE / R2),
                  reversed for lower-is-better (RMSE / |PBIAS|).
      * hover   = popup with the station ID, the matched feature, n_obs,
                  and the full metric table (KGE / NSE / RMSE / PBIAS / R²)
                  per species — same columns as the Performance section.

    The map starts LOCKED (no pan/zoom, no scroll-wheel) with the Esri
    Satellite basemap as the default — both at the user's request.  A
    lock button at the top right toggles interactivity.  An OSM/Satellite
    layer switcher and Basins/River/Stations/Gray-obs toggles live in
    the same control panel as before.
    """
    obs_cfg = _ci.get_observation_config(model_config)
    spatial = _ci.get_spatial_mapping(model_config)
    hostmodel = spatial.get("hostmodel") or "mizuroute"

    river_shp = obs_cfg.get("river_network_shapefile")
    basin_shp = obs_cfg.get("basin_shapefile")
    river_gj, river_bounds = (_load_shapefile_as_geojson(river_shp)
                              if river_shp and os.path.isfile(river_shp)
                              else (None, None))
    basin_gj, basin_bounds = (_load_shapefile_as_geojson(basin_shp)
                              if basin_shp and os.path.isfile(basin_shp)
                              else (None, None))
    if not river_gj and not basin_gj:
        return ""

    # Pick the primary interactive layer the same way Plot_h5_driver does.
    if hostmodel.lower() == "summa":
        primary_gj = basin_gj or river_gj
        primary_key = spatial.get("basin_mapping_key") or "HRU_ID"
        is_polygon = bool(basin_gj)
    else:
        primary_gj = river_gj or basin_gj
        primary_key = spatial.get("river_network_mapping_key") or "SegId"
        is_polygon = not bool(river_gj)
    if primary_gj is None:
        return ""

    # Load observations and pull station locations.
    obs_csv = obs_cfg.get("user_observation_csv")
    grqa_dir = obs_cfg.get("grqa_local_data_path")
    station_locations: Dict[str, tuple] = {}
    if obs_csv and os.path.isfile(obs_csv):
        try:
            df = pd.read_csv(obs_csv)
            for _, row in df.drop_duplicates(
                    'station_id' if 'station_id' in df.columns else df.columns[0]
            ).iterrows():
                if 'lat' in df.columns and 'lon' in df.columns:
                    sid = str(row.get('station_id', row.get('source', '')))
                    if sid:
                        station_locations[sid] = (float(row['lat']),
                                                  float(row['lon']))
        except Exception:
            pass
    if not station_locations and grqa_dir:
        stn_csv = os.path.join(grqa_dir, "grqa_clipped_stations.csv")
        if os.path.isfile(stn_csv):
            try:
                df = pd.read_csv(stn_csv)
                for _, row in df.iterrows():
                    sid = str(row.get('site_id', ''))
                    if sid:
                        station_locations[sid] = (float(row['lat_wgs84']),
                                                  float(row['lon_wgs84']))
            except Exception:
                pass
    if not station_locations:
        return ""

    # Match using the SAME shared helpers Plot_h5_driver uses.
    try:
        sm = _import_spatial_matching()
    except ImportError as e:
        logger.warning(f"Observation map: shared spatial_matching missing — {e}")
        return ""
    s2f, primary = sm.match_stations(
        station_locations,
        hostmodel=hostmodel,
        river_geojson=river_gj,
        basin_geojson=basin_gj,
        river_mapping_key=spatial.get("river_network_mapping_key") or "SegId",
        basin_mapping_key=spatial.get("basin_mapping_key") or "HRU_ID",
        log=lambda *a, **kw: None,
    )

    # ── Per-station performance summary ─────────────────────────────
    # Build station_stats[sid] = {
    #     "n_obs":       int,
    #     "metrics":     {metric_name: aggregate_value across species},
    #     "by_species":  {species: {metric_name: value, ...}},
    # }
    #
    # n_obs comes from the obs CSV (count of rows per station).
    # Metrics come from performance_metrics + station_to_feature: a
    # station's reach is the feature it was matched to, and per
    # (species, reach) we have one performance row.  When a station
    # corresponds to multiple species we average them per metric for
    # the marker fill, but keep the per-species breakdown for the
    # hover tooltip.
    metric_keys = ("KGE", "NSE", "RMSE", "PBIAS", "R2")
    station_stats: Dict[str, Dict[str, Any]] = {
        sid: {"n_obs": 0, "metrics": {}, "by_species": {}}
        for sid in station_locations
    }
    # 1) n_obs from obs CSV
    if obs_csv and os.path.isfile(obs_csv):
        try:
            _df = pd.read_csv(obs_csv)
            if "station_id" in _df.columns:
                _counts = _df.groupby("station_id").size().to_dict()
            elif "source" in _df.columns:
                _counts = _df.groupby("source").size().to_dict()
            else:
                _counts = {}
            for sid, n in _counts.items():
                if str(sid) in station_stats:
                    station_stats[str(sid)]["n_obs"] = int(n)
        except Exception:
            pass

    # 2) Per-station, per-species metrics by reverse-lookup through s2f
    #    + performance_metrics (which is keyed on reach_id + species).
    #
    # ID normalisation is critical here: shapefiles store integer-like
    # IDs as FLOATS (``GRU_ID=740457190.0``), which the GeoJSON loader
    # then stringifies to ``"740457190.0"``.  Meanwhile ObjectiveFunction
    # converts H5 reach IDs through ``int()`` so ``performance_metrics``
    # carries ``"740457190"``.  Those two strings never match and every
    # station ends up unlinked even when the underlying IDs are the
    # same numeric value.  We canonicalise to ``str(int(x))`` whenever
    # possible — drops the ``.0``, drops leading zeros, drops whitespace
    # — so both sides converge on a single key form.
    def _norm_id(v: Any) -> str:
        if v is None:
            return ""
        s = str(v).strip()
        if not s:
            return ""
        if s.endswith(".0"):
            s = s[:-2]
        try:
            return str(int(float(s)))
        except (ValueError, TypeError):
            return s

    # SUMMA bridge: the H5 output (and therefore performance_metrics)
    # may also use sequential reach IDs ("1", "2", ...) instead of the
    # real shapefile GRU_ID.  We register every performance row under
    # the canonical normalised ID AND the sequential-position remap, so
    # the lookup succeeds for either ID space — mizuRoute (already
    # aligned), SUMMA-with-real-IDs (needs only normalisation), and
    # SUMMA-with-sequential-IDs (needs the position remap).
    if performance_metrics:
        # Position remap: "1" → real_id_of_feat_0, "2" → real_id_of_feat_1, ...
        # All real_id values pass through _norm_id so the keys match
        # whatever normalised form the performance rows produce.
        seq_to_real: Dict[str, str] = {}
        real_set: set = set()
        for i, feat in enumerate(primary_gj.get('features', [])):
            real_id = _norm_id(feat['properties'].get(primary_key, ''))
            if real_id:
                seq_to_real[str(i + 1)] = real_id
                real_set.add(real_id)
        _pm_by_reach: Dict[str, List[Dict[str, Any]]] = {}
        for row in performance_metrics:
            _rid = _norm_id(row.get("reach_id", ""))
            if not _rid:
                continue
            keys_to_register = {_rid}
            # If the row's reach_id looks like a sequential index and we
            # have a remap, also register the real shapefile ID.
            if _rid in seq_to_real:
                keys_to_register.add(seq_to_real[_rid])
            # If the row's reach_id IS a real shapefile ID, also
            # register the sequential form.
            for seq, real in seq_to_real.items():
                if real == _rid:
                    keys_to_register.add(seq)
                    break
            for k in keys_to_register:
                _pm_by_reach.setdefault(k, []).append(row)
        for sid, fid in s2f.items():
            for row in _pm_by_reach.get(_norm_id(fid), []):
                sp = str(row.get("species", ""))
                m = {k: row.get(k) for k in metric_keys
                     if row.get(k) is not None}
                if sp:
                    station_stats[sid]["by_species"][sp] = m
        # Diagnostic — surfaces in the calibration run log so the user
        # can verify the lookup actually flowed metrics into stations.
        _n_with_metrics = sum(
            1 for st in station_stats.values() if st["by_species"])
        logger.info(
            f"[obs-map] linked metrics to {_n_with_metrics}/{len(station_stats)} "
            f"stations  (performance rows: {len(performance_metrics)}, "
            f"seq-real remap entries: {len(seq_to_real)})"
        )
        # Aggregate per-metric mean across the species for that station.
        for sid, st in station_stats.items():
            agg: Dict[str, float] = {}
            for k in metric_keys:
                vals = [m[k] for m in st["by_species"].values()
                        if isinstance(m.get(k), (int, float))]
                if vals:
                    agg[k] = float(sum(vals) / len(vals))
            st["metrics"] = agg

    # 3) The objective metric drives the marker colour.
    #
    # Each metric has its own "good vs bad" semantics — for some the
    # range is fixed by definition (KGE/NSE/R² are bounded above at 1,
    # below 0 is "no better than the mean"), for others (RMSE, MAE,
    # MAPE) the magnitude depends entirely on the data units and there
    # is no theoretical max.  The colormap should stay calibratable
    # between runs for the bounded metrics so a station's colour means
    # the same thing every time the report is generated — i.e. KGE=0.7
    # is the same shade no matter what the rest of the dataset looks
    # like.  Unbounded metrics fall back to the observed-range
    # behaviour from before, but with a note in the legend so the
    # user knows the scale is data-dependent.
    #
    # The literature standard ranges I use here:
    #   KGE   ∈ (-∞, 1]   — bounded above, 1 = perfect.  Visualise [0,1];
    #                       values <0 (worse than the mean of obs) clamp
    #                       to "worst" colour rather than blowing up the
    #                       scale.
    #   NSE   ∈ (-∞, 1]   — same.  0 = "as good as predicting the mean".
    #   R²    ∈ [0, 1]    — by definition non-negative for non-trivial
    #                       linear fits; clamp negative values to 0.
    #   r     ∈ [-1, 1]   — Pearson correlation, 1 = perfect positive.
    #   PBIAS ∈ (-∞, ∞)   — sign-aware bias; we use |PBIAS| for the
    #                       colormap with Moriasi et al. (2007) "good"
    #                       cutoffs around 25 % for streamflow.
    #   RMSE  ∈ [0, ∞)    — no theoretical max; observed-range.
    #   MAE   ∈ [0, ∞)    — same.
    #   MAPE  ∈ [0, ∞) %  — same; almost always clamped to [0, 100].
    #   LogNSE ∈ (-∞, 1]  — like NSE but on log-transformed flows.
    #
    # Sources: Gupta et al. (2009), Nash & Sutcliffe (1970), Moriasi
    # et al. (2007).
    METRIC_RANGES: Dict[str, Optional[Dict[str, Any]]] = {
        # higher = better, fixed range
        "KGE":   {"vmin": 0.0, "vmax": 1.0, "reversed": False, "use_abs": False},
        "NSE":   {"vmin": 0.0, "vmax": 1.0, "reversed": False, "use_abs": False},
        "R2":    {"vmin": 0.0, "vmax": 1.0, "reversed": False, "use_abs": False},
        "LOGNSE": {"vmin": 0.0, "vmax": 1.0, "reversed": False, "use_abs": False},
        "R":     {"vmin": -1.0, "vmax": 1.0, "reversed": False, "use_abs": False},
        # symmetric around zero — visualise |.|, low is good
        "PBIAS": {"vmin": 0.0, "vmax": 25.0, "reversed": True,  "use_abs": True},
        "MAPE":  {"vmin": 0.0, "vmax": 50.0, "reversed": True,  "use_abs": False},
        # unbounded — use observed range
        "RMSE":  None,
        "MAE":   None,
    }
    objective_metric = str(
        (calibration_settings or {}).get("objective_function") or "KGE"
    ).upper()
    _preset = METRIC_RANGES.get(objective_metric)

    metric_vals = []
    for st in station_stats.values():
        v = st["metrics"].get(objective_metric)
        if isinstance(v, (int, float)):
            metric_vals.append(abs(v) if (_preset and _preset.get("use_abs"))
                               else v)

    if _preset is not None:
        # Bounded metric — keep the literature range so the colour scale
        # means the same thing across runs.  Tag the legend so the user
        # knows the bounds are conceptual, not derived from this run.
        colorbar = {
            "metric":     objective_metric,
            "vmin":       _preset["vmin"],
            "vmax":       _preset["vmax"],
            "reversed":   _preset["reversed"],
            "use_abs":    _preset["use_abs"],
            "range_kind": "theoretical",
        }
    elif metric_vals:
        # Unbounded metric (RMSE / MAE / ...) — derive from the data.
        _vmin, _vmax = min(metric_vals), max(metric_vals)
        if abs(_vmax - _vmin) < 1e-9:
            # All stations identical → expand by 5% so the colorbar
            # gradient still renders and one value doesn't get the
            # "worst" colour purely by accident.
            _pad = max(abs(_vmin) * 0.05, 1e-6)
            _vmin = _vmin - _pad
            _vmax = _vmax + _pad
        colorbar = {
            "metric":     objective_metric,
            "vmin":       _vmin,
            "vmax":       _vmax,
            "reversed":   True,   # RMSE / MAE: lower is better
            "use_abs":    False,
            "range_kind": "observed",
        }
    else:
        # Unbounded metric AND no values at all — surface neutral range
        # so the legend doesn't display garbage; markers will all gray.
        colorbar = {
            "metric":     objective_metric,
            "vmin":       0.0,
            "vmax":       1.0,
            "reversed":   False,
            "use_abs":    False,
            "range_kind": "fallback",
        }

    # Stable colour by sorted feature ID (matches the config viewer's
    # _smart_sort_key behaviour for numeric-string IDs).
    def _sk(s):
        s = str(s)
        try:
            return (0, int(s), s)
        except ValueError:
            return (1, s.lower(), s)
    fids_sorted = sorted(
        {str(f['properties'].get(primary_key, ''))
         for f in primary_gj.get('features', [])
         if f['properties'].get(primary_key)},
        key=_sk,
    )
    map_fid_color = {fid: _RESULTS_MAP_COLORWAY[i % len(_RESULTS_MAP_COLORWAY)]
                     for i, fid in enumerate(fids_sorted)}

    # Compute map centre + bounds.
    bb = primary_gj.get('features') and basin_bounds or river_bounds
    if not bb:
        return ""
    centre = [(bb[1] + bb[3]) / 2, (bb[0] + bb[2]) / 2]

    # Build the HTML block.
    div_id = "obs-map-leaflet"
    n_primary = sum(1 for s in s2f if s in primary)
    n_secondary = len(s2f) - n_primary
    n_unmatched = len(station_locations) - len(s2f)
    sub = (f"{len(s2f)} matched ({n_primary} primary"
           + (f", {n_secondary} secondary" if n_secondary else "")
           + (f", {n_unmatched} unmatched" if n_unmatched else "")
           + f") — hostmodel: {hostmodel}")
    html = ['<div class="section" id="obs-map">']
    html.append('<h2>Observation Map</h2>')
    html.append(f'<p class="muted" style="margin-top:-8px">{sub}</p>')
    # CSS + Leaflet (inline so the section is self-contained).
    html.append(
        '<link rel="stylesheet" href="https://unpkg.com/leaflet@1.9.4/dist/leaflet.css"/>'
        '<script src="https://unpkg.com/leaflet@1.9.4/dist/leaflet.js"></script>'
    )
    html.append(f'<div id="{div_id}" style="height:520px;border-radius:8px;'
                f'border:1px solid var(--border,#ddd);"></div>')

    has_river = bool(river_gj and river_gj is not primary_gj or
                     (river_gj and not is_polygon))
    has_basin_overlay = bool(basin_gj and basin_gj is not primary_gj)
    n_secondary = sum(1 for sid in s2f if sid not in primary)
    has_any_gray = bool(
        n_secondary or
        any(sid not in s2f for sid in station_locations)
    )
    payload = {
        'primary': primary_gj,
        'river': river_gj if river_gj is not primary_gj else None,
        'basin': basin_gj if basin_gj is not primary_gj else None,
        'primary_key': primary_key,
        'is_polygon': is_polygon,
        'station_data': {sid: [lat, lon]
                         for sid, (lat, lon) in station_locations.items()},
        'station_to_feature': s2f,
        'pouring_point': {sid: True for sid in primary},
        'map_fid_color': map_fid_color,
        'centre': centre,
        'bounds': [[bb[1], bb[0]], [bb[3], bb[2]]],
        'feature_label': spatial.get("feature_label") or primary_key,
        'fids_sorted': fids_sorted,
        'has_river': bool(river_gj is not None),
        'has_basin_overlay': has_basin_overlay,
        'has_secondary': has_any_gray,
        'station_stats': station_stats,
        'colorbar': colorbar,
    }
    # Rich Leaflet block:
    #  - Satellite default basemap (Esri WorldImagery)
    #  - Lock button (default LOCKED, no pan/zoom) at top-right
    #  - Marker radius ∝ sqrt(n_obs), fill = colormap(metric)
    #  - Hover tooltip + popup with the full per-species metric table
    #  - Colorbar legend (bottom-left) showing the objective metric scale
    html.append(
        '<script>(function(){\n'
        # ── Wait for Leaflet to load before building the map ───────
        # When the report is opened on a slow connection, or if the
        # unpkg CDN is briefly unreachable, the external Leaflet
        # <script> tag above may not have finished executing by the
        # time this inline script runs.  Without this poll the map
        # silently fails ("L is not defined") and the container stays
        # blank.  We poll up to 10s, then surface a clear fallback
        # message so the user knows it's a network problem, not a
        # calibration bug.
        f'var div = document.getElementById("{div_id}");\n'
        'if(div){div.innerHTML='
        '\'<div style="display:flex;align-items:center;justify-content:center;'
        'height:100%;color:#666;font:13px/1.4 -apple-system,Segoe UI,sans-serif">'
        'Loading map\\u2026</div>\';}\n'
        'var _tries = 0;\n'
        'function _go(){\n'
        '  if(typeof L === "undefined"){\n'
        '    if(++_tries > 200){\n'
        '      if(div){div.innerHTML='
        '\'<div style="padding:20px;color:#a33;font:13px/1.4 -apple-system,Segoe UI,sans-serif">'
        '<b>Map failed to load.</b><br/>The Leaflet library could not be fetched from unpkg.com '
        '\\u2014 check the machine\\u2019s internet access or mirror Leaflet locally '
        '(see <code>supporting_scripts/requirements.txt</code>).</div>\';}\n'
        '      return;\n'
        '    }\n'
        '    return setTimeout(_go, 50);\n'
        '  }\n'
        '  if(div) div.innerHTML = "";\n'
        '  _build();\n'
        '}\n'
        'function _build(){\n'
        f'var P = {json.dumps(payload)};\n'
        f'var map = L.map("{div_id}",{{zoomControl:true,'
        'dragging:false,scrollWheelZoom:false,doubleClickZoom:false,'
        'touchZoom:false,boxZoom:false,keyboard:false});\n'
        'map.fitBounds(P.bounds,{padding:[10,10]});\n'
        '// ── Base layers (Satellite is DEFAULT now) ──\n'
        'var osm = L.tileLayer("https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png",\n'
        '  {maxZoom:19, attribution:"&copy; OpenStreetMap contributors"});\n'
        'var sat = L.tileLayer(\n'
        '  "https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}",\n'
        '  {maxZoom:19, attribution:"Tiles &copy; Esri \\u2014 Source: Esri, Maxar, Earthstar Geographics, USDA, USGS"});\n'
        'sat.addTo(map);\n'
        'L.control.layers({"Satellite":sat,"OpenStreetMap":osm},null,{position:"topright"}).addTo(map);\n'
        '// ── Feature colour + reachable styling ──\n'
        'function colOf(fid){return P.map_fid_color[fid]||"#888";}\n'
        'function styleFn(feat){\n'
        '  var fid=String(feat.properties[P.primary_key]||"");\n'
        '  var c=colOf(fid);\n'
        '  if(P.is_polygon){return {color:"#333",weight:1,fillColor:c,fillOpacity:0.35};}\n'
        '  return {color:c,weight:3,opacity:0.9};\n'
        '}\n'
        'var primaryLayer = L.geoJSON(P.primary,{style:styleFn,\n'
        '  onEachFeature:function(f,l){\n'
        '    var fid=String(f.properties[P.primary_key]||"");\n'
        '    l.bindTooltip(P.feature_label+" "+fid,{sticky:true,opacity:0.9});\n'
        '  }}).addTo(map);\n'
        'var riverRefLayer = null;\n'
        'if(P.river && P.is_polygon){\n'
        '  riverRefLayer = L.geoJSON(P.river,{\n'
        '    style:{color:"#fff",weight:2.5,opacity:0.9},interactive:false}).addTo(map);\n'
        '}\n'
        'var basinRefLayer = null;\n'
        'if(P.has_basin_overlay && !P.is_polygon){\n'
        '  basinRefLayer = L.geoJSON(P.basin,{\n'
        '    style:{color:"#bbb",weight:1,fillOpacity:0.06},interactive:false}).addTo(map);\n'
        '}\n'
        '// ── Colormap: RdYlGn for higher-is-better, reversed otherwise ──\n'
        'var CB = P.colorbar;\n'
        'var PALETTE = ["#d73027","#f46d43","#fdae61","#fee08b","#ffffbf",\n'
        '               "#d9ef8b","#a6d96a","#66bd63","#1a9850"];\n'
        'if(CB.reversed) PALETTE = PALETTE.slice().reverse();\n'
        'function cmap(v){\n'
        '  if(v===null||v===undefined||isNaN(v)) return "#888";\n'
        '  var x = (v - CB.vmin) / (CB.vmax - CB.vmin);\n'
        '  if(!isFinite(x)) x = 0.5;\n'
        '  x = Math.max(0, Math.min(1, x));\n'
        '  var i = Math.floor(x * (PALETTE.length - 1));\n'
        '  return PALETTE[i];\n'
        '}\n'
        '// ── Marker radius from n_obs ──\n'
        'var ALL_N = Object.keys(P.station_stats).map(function(k){return P.station_stats[k].n_obs||0;});\n'
        'var MAX_N = Math.max.apply(null, ALL_N.concat([1]));\n'
        'function radiusOf(n){\n'
        '  if(!n) return 4;\n'
        '  return 4 + 10 * Math.sqrt(n / MAX_N);\n'
        '}\n'
        '// ── Hover/popup HTML builder with per-species metric table ──\n'
        'function popupHtml(sid){\n'
        '  var st = P.station_stats[sid] || {};\n'
        '  var fid = P.station_to_feature[sid] || null;\n'
        '  var lines = ["<b>Station "+sid+"</b>"];\n'
        '  if(fid){\n'
        '    var role = P.pouring_point[sid] ? "pouring-point" : "secondary";\n'
        '    lines.push("Matched "+P.feature_label+": "+fid+" <i>("+role+")</i>");\n'
        '  } else lines.push("<i>unmatched</i>");\n'
        '  lines.push("Observations: <b>"+(st.n_obs||0)+"</b>");\n'
        '  if(st.by_species && Object.keys(st.by_species).length){\n'
        '    var rows = "";\n'
        '    var allMetrics = {};\n'
        '    Object.values(st.by_species).forEach(function(m){\n'
        '      Object.keys(m||{}).forEach(function(k){allMetrics[k]=true;});\n'
        '    });\n'
        '    var metrics = Object.keys(allMetrics);\n'
        '    var header = "<th style=\\"text-align:left;padding:2px 6px\\">Species</th>";\n'
        '    metrics.forEach(function(k){\n'
        '      header += "<th style=\\"text-align:right;padding:2px 6px\\">"+k+"</th>";\n'
        '    });\n'
        '    Object.keys(st.by_species).forEach(function(sp){\n'
        '      var m = st.by_species[sp];\n'
        '      var row = "<td style=\\"padding:2px 6px\\">"+sp+"</td>";\n'
        '      metrics.forEach(function(k){\n'
        '        var v = m[k];\n'
        '        row += "<td style=\\"text-align:right;padding:2px 6px\\">"+\n'
        '               (typeof v === "number" ? v.toFixed(3) : "&mdash;")+"</td>";\n'
        '      });\n'
        '      rows += "<tr>"+row+"</tr>";\n'
        '    });\n'
        '    lines.push("<table style=\\"border-collapse:collapse;font-size:11px;margin-top:4px\\">"+\n'
        '               "<tr>"+header+"</tr>"+rows+"</table>");\n'
        '  }\n'
        '  return lines.join("<br>");\n'
        '}\n'
        '// ── Station markers ──\n'
        'var ZERO = {radius:0.001,opacity:0,fillOpacity:0,weight:0};\n'
        'var stationLayer = L.layerGroup();\n'
        'var stationStyles = {};\n'
        'Object.keys(P.station_data).forEach(function(sid){\n'
        '  var ll = P.station_data[sid];\n'
        '  var st = P.station_stats[sid] || {n_obs:0,metrics:{}};\n'
        '  var fid = P.station_to_feature[sid] || null;\n'
        '  var isPrim = !!P.pouring_point[sid];\n'
        '  var v = st.metrics ? st.metrics[CB.metric] : undefined;\n'
        '  if(typeof v === "number" && CB.use_abs) v = Math.abs(v);\n'
        '  // Marker fill comes from the colormap when we have a metric\n'
        '  // value for this station; otherwise GRAY so it\\u2019s obvious the\n'
        '  // metric lookup didn\\u2019t hit (rather than silently falling\n'
        '  // through to the basin/reach polygon colour, which made every\n'
        '  // marker look blue when the link failed).\n'
        '  var fill = (typeof v === "number") ? cmap(v) : "#888";\n'
        '  var r = radiusOf(st.n_obs);\n'
        '  var prim = {radius:r,fillColor:fill,color:"#fff",weight:1.6,opacity:1,fillOpacity:0.95};\n'
        '  var sec  = {radius:Math.max(3,r-1),fillColor:"#888",color:"#fff",weight:1,opacity:0.9,fillOpacity:0.7};\n'
        '  stationStyles[sid] = {primary:prim,secondary:sec,isPrimary:isPrim};\n'
        '  var m = L.circleMarker([ll[0],ll[1]], isPrim?prim:sec);\n'
        '  var tip = "Station "+sid + (fid?" \\u2192 "+fid:" (unmatched)") +\n'
        '            "  \\u2014  "+CB.metric+": "+\n'
        '            (typeof v === "number" ? v.toFixed(3) : "n/a") +\n'
        '            "  ("+st.n_obs+" obs)";\n'
        '  m.bindTooltip(tip,{direction:"top",sticky:true,opacity:0.95});\n'
        '  m.bindPopup(popupHtml(sid));\n'
        '  m.__sid = sid;\n'
        '  stationLayer.addLayer(m);\n'
        '});\n'
        'stationLayer.addTo(map);\n'
        '// ── Legend + toggle control (top-left) ──\n'
        'var CtrlCfg = L.Control.extend({\n'
        '  options:{position:"topleft"},\n'
        '  onAdd:function(){\n'
        '    var div = L.DomUtil.create("div","leaflet-bar");\n'
        '    div.style.cssText="background:rgba(255,255,255,.94);color:#222;padding:8px 10px;'
        'font:12px/1.4 -apple-system,Segoe UI,sans-serif;border-radius:6px;'
        'min-width:200px;max-width:260px;max-height:60vh;overflow:auto;'
        'box-shadow:0 1px 5px rgba(0,0,0,.3);";\n'
        '    L.DomEvent.disableClickPropagation(div);\n'
        '    L.DomEvent.disableScrollPropagation(div);\n'
        '    var html_ = "<div style=\\"font-weight:600;margin-bottom:6px\\">Map controls</div>";\n'
        '    function tog(id,label,checked){\n'
        '      return "<label style=\\"display:flex;align-items:center;gap:6px;margin:3px 0\\">"+\n'
        '             "<input type=\\"checkbox\\" id=\\""+id+"\\" "+(checked?"checked":"")+\n'
        '             " style=\\"accent-color:#0066cc\\"/>"+label+"</label>";\n'
        '    }\n'
        '    if(P.has_basin_overlay && !P.is_polygon) html_ += tog("tg_basin","Basins",true);\n'
        '    if(P.has_river && P.is_polygon)         html_ += tog("tg_river","River network",true);\n'
        '    html_ += tog("tg_stations","Observation stations",true);\n'
        '    if(P.has_secondary)                     html_ += tog("tg_gray","Gray obs (secondary)",true);\n'
        '    div.innerHTML = html_;\n'
        '    return div;\n'
        '  }\n'
        '});\n'
        'new CtrlCfg().addTo(map);\n'
        '// ── Colorbar legend (bottom-left) ──\n'
        'var Colorbar = L.Control.extend({\n'
        '  options:{position:"bottomleft"},\n'
        '  onAdd:function(){\n'
        '    var d = L.DomUtil.create("div","leaflet-bar");\n'
        '    d.style.cssText="background:rgba(255,255,255,.94);color:#222;padding:8px 10px;'
        'font:12px/1.3 -apple-system,Segoe UI,sans-serif;border-radius:6px;'
        'box-shadow:0 1px 5px rgba(0,0,0,.3);min-width:170px;";\n'
        '    var stops = PALETTE.map(function(c,i){\n'
        '      return c+" "+((i/(PALETTE.length-1))*100).toFixed(1)+"%";\n'
        '    }).join(",");\n'
        '    var title = (CB.use_abs?"|":"")+CB.metric+(CB.use_abs?"|":"");\n'
        '    var kindHint = "";\n'
        '    if(CB.range_kind === "theoretical"){\n'
        '      kindHint = "<span style=\\"font-size:10px;color:#666;margin-left:6px\\">'
        '(theoretical range)</span>";\n'
        '    } else if(CB.range_kind === "observed"){\n'
        '      kindHint = "<span style=\\"font-size:10px;color:#666;margin-left:6px\\">'
        '(observed range)</span>";\n'
        '    }\n'
        '    d.innerHTML = "<div style=\\"font-weight:600;margin-bottom:4px\\">"+title+'
        '" (marker fill)"+kindHint+"</div>"+\n'
        '      "<div style=\\"height:12px;border-radius:3px;background:linear-gradient(to right,"+stops+\n'
        '      ");border:1px solid #888\\"></div>"+\n'
        '      "<div style=\\"display:flex;justify-content:space-between;'
        'font-size:11px;margin-top:2px\\"><span>"+CB.vmin.toFixed(2)+'
        '"</span><span>"+CB.vmax.toFixed(2)+"</span></div>"+\n'
        '      "<div style=\\"font-size:10px;color:#555;margin-top:6px\\">'
        'Marker <b>size</b> &prop; observations</div>";\n'
        '    return d;\n'
        '  }\n'
        '});\n'
        'new Colorbar().addTo(map);\n'
        '// ── Lock button (top-right, default LOCKED) ──\n'
        'var locked = true;\n'
        'function applyLock(){\n'
        '  if(locked){\n'
        '    map.dragging.disable(); map.scrollWheelZoom.disable();\n'
        '    map.doubleClickZoom.disable(); map.boxZoom.disable();\n'
        '    map.touchZoom && map.touchZoom.disable();\n'
        '    map.keyboard && map.keyboard.disable();\n'
        '    if(map.zoomControl){ map.removeControl(map.zoomControl); }\n'
        '  } else {\n'
        '    map.dragging.enable(); map.scrollWheelZoom.enable();\n'
        '    map.doubleClickZoom.enable(); map.boxZoom.enable();\n'
        '    map.touchZoom && map.touchZoom.enable();\n'
        '    map.keyboard && map.keyboard.enable();\n'
        '    if(!map._zoomControlPresent){\n'
        '      map.zoomControl = L.control.zoom({position:"topleft"}).addTo(map);\n'
        '      map._zoomControlPresent = true;\n'
        '    }\n'
        '  }\n'
        '}\n'
        'applyLock();\n'
        'var LockBtn = L.Control.extend({\n'
        '  options:{position:"topright"},\n'
        '  onAdd:function(){\n'
        '    var b = L.DomUtil.create("button","leaflet-bar");\n'
        '    function refresh(){ b.innerHTML = locked ? "&#x1F512;" : "&#x1F513;";\n'
        '      b.title = locked ? "Map locked \\u2014 click to unlock" : "Map unlocked \\u2014 click to lock"; }\n'
        '    refresh();\n'
        '    b.style.cssText="width:34px;height:34px;background:#fff;color:#333;border:none;'
        'border-radius:4px;cursor:pointer;font-size:16px;line-height:34px;text-align:center;'
        'box-shadow:0 1px 5px rgba(0,0,0,.3);margin-top:4px";\n'
        '    L.DomEvent.disableClickPropagation(b);\n'
        '    b.addEventListener("click",function(){ locked = !locked; applyLock(); refresh(); });\n'
        '    return b;\n'
        '  }\n'
        '});\n'
        'new LockBtn().addTo(map);\n'
        '// ── Recenter button (top-right) ──\n'
        'var Recenter = L.Control.extend({\n'
        '  options:{position:"topright"},\n'
        '  onAdd:function(){\n'
        '    var b = L.DomUtil.create("button","leaflet-bar");\n'
        '    b.title = "Recenter map";\n'
        '    b.innerHTML = "&#x2316;";\n'
        '    b.style.cssText="width:34px;height:34px;background:#fff;color:#333;border:none;'
        'border-radius:4px;cursor:pointer;font-size:18px;line-height:34px;text-align:center;'
        'box-shadow:0 1px 5px rgba(0,0,0,.3);margin-top:4px";\n'
        '    L.DomEvent.disableClickPropagation(b);\n'
        '    b.addEventListener("click",function(){ map.fitBounds(P.bounds,{padding:[10,10]}); });\n'
        '    return b;\n'
        '  }\n'
        '});\n'
        'new Recenter().addTo(map);\n'
        '// ── Toggle handlers ──\n'
        'setTimeout(function(){\n'
        '  var c;\n'
        '  c = document.getElementById("tg_basin");\n'
        '  if(c && basinRefLayer){ c.onchange=function(){ '
        'if(this.checked) basinRefLayer.addTo(map); else map.removeLayer(basinRefLayer); }; }\n'
        '  c = document.getElementById("tg_river");\n'
        '  if(c && riverRefLayer){ c.onchange=function(){ '
        'if(this.checked) riverRefLayer.addTo(map); else map.removeLayer(riverRefLayer); }; }\n'
        '  c = document.getElementById("tg_stations");\n'
        '  if(c){ c.onchange=function(){ '
        'if(this.checked) stationLayer.addTo(map); else map.removeLayer(stationLayer); }; }\n'
        '  c = document.getElementById("tg_gray");\n'
        '  if(c){ c.onchange=function(){\n'
        '    var show = this.checked;\n'
        '    stationLayer.eachLayer(function(mk){\n'
        '      var st = stationStyles[mk.__sid];\n'
        '      if(!st) return;\n'
        '      if(st.isPrimary){ mk.setStyle(st.primary); return; }\n'
        '      mk.setStyle(show ? st.secondary : ZERO);\n'
        '    });\n'
        '  }; }\n'
        '},0);\n'
        '}\n'  # close _build()
        '_go();\n'  # kick off the readiness poll
        '})();</script>\n'
    )
    html.append('</div>')
    return "\n".join(html)


# =========================================================================
# PDF Export — floating "Download PDF" button + html2pdf.js bootstrap.
# Captures every section with a proper heading + white background, so the
# resulting PDF is a clean printable summary of the entire report.
# =========================================================================

def _build_pdf_export_block(project_name: str) -> str:
    """Floating button that exports every section as a multi-page PDF.

    Uses ``html2pdf.js`` (= ``html2canvas`` + ``jsPDF``) loaded from the
    same CDN convention everything else in the report uses.  At export
    time we:

      1. Build a temporary off-screen wrapper with a forced white
         background and dark text so the PDF doesn't carry the report's
         dark theme.
      2. Add a title page (project name, "OpenWQ Calibration Results
         Report" subtitle, generated-at timestamp).
      3. Walk a fixed list of section IDs in render order; for each
         present section we clone its DOM into the wrapper, strip
         interactive-only elements (Copy buttons, hover styles), and
         force every nested element to a printable colour palette via
         the ``.pdf-wrap`` CSS scope.  Section headings (``<h2>``)
         carry over from the source so each map / figure / table has
         its own "context heading" in the PDF.
      4. Call ``html2pdf().from(wrapper).save()`` with ``backgroundColor:
         '#ffffff'`` and ``useCORS: true`` so the Leaflet satellite
         tiles get captured too.

    The button itself is fixed-positioned bottom-right and excluded
    from the export via the ``.no-pdf`` class.
    """
    project_name_js = json.dumps(project_name or "OpenWQ Calibration")
    # CSS: scope every printable element to known light colours; hide
    # everything we don't want printed (the Theme toggle, the sidebar,
    # the floating button itself).
    css = """
.pdf-wrap, .pdf-wrap * {
    background:#ffffff !important;
    color:#222 !important;
    border-color:#ccc !important;
    box-shadow:none !important;
    text-shadow:none !important;
}
.pdf-wrap pre, .pdf-wrap code {
    background:#f5f5f5 !important;
    color:#222 !important;
    border:1px solid #ddd !important;
}
.pdf-wrap a { color:#0066cc !important; }
.pdf-wrap .copy-btn, .pdf-wrap .no-pdf, .pdf-wrap button { display:none !important; }
.pdf-wrap table { border-collapse:collapse; }
.pdf-wrap th, .pdf-wrap td { border:1px solid #ddd; padding:4px 8px; }
.pdf-wrap h1, .pdf-wrap h2, .pdf-wrap h3 { color:#111 !important; }
.pdf-page-break { page-break-before:always; }

#pdfExportBtn {
    position:fixed; right:18px; bottom:18px;
    z-index:9999; padding:10px 18px;
    background:#0066cc; color:#fff;
    border:none; border-radius:30px; font-weight:600; font-size:14px;
    cursor:pointer; box-shadow:0 3px 12px rgba(0,0,0,.3);
    display:inline-flex; align-items:center; gap:8px;
}
#pdfExportBtn:hover  { background:#004499; }
#pdfExportBtn:disabled{ background:#888; cursor:wait; }
#pdfExportBtn .spin { display:none; width:14px; height:14px;
    border:2px solid #fff; border-top-color:transparent;
    border-radius:50%; animation:pdfspin .9s linear infinite; }
#pdfExportBtn.busy .spin { display:inline-block; }
@keyframes pdfspin { to { transform:rotate(360deg); } }
@media print { #pdfExportBtn { display:none !important; } }
"""
    return f"""
<style>{css}</style>
<button id="pdfExportBtn" class="no-pdf" title="Export the entire report as a PDF" onclick="_owqExportPDF()">
  <span class="spin"></span>
  <span class="label">&#x1F4C4; Download PDF</span>
</button>
<script src="https://cdnjs.cloudflare.com/ajax/libs/html2pdf.js/0.10.1/html2pdf.bundle.min.js"></script>
<script>
function _owqExportPDF() {{
    var btn = document.getElementById('pdfExportBtn');
    if(!btn) return;
    if(typeof html2pdf === 'undefined') {{
        alert('PDF library is still loading — please try again in a moment.');
        return;
    }}
    btn.disabled = true; btn.classList.add('busy');
    btn.querySelector('.label').textContent = ' Building PDF…';

    // Sections to include, in render order.  Skipped silently if absent.
    var SECTIONS = ['summary','obs-map','best-params','convergence',
                    'param-evolution','param-correlations','sensitivity',
                    'performance','run-best'];

    // Off-screen wrapper — forced white background, dark text.
    var wrap = document.createElement('div');
    wrap.className = 'pdf-wrap';
    wrap.style.cssText =
        'position:absolute;left:-10000px;top:0;width:794px;'+   // 210mm @ 96dpi
        'background:#ffffff;color:#222;'+
        'font-family:-apple-system,Segoe UI,Helvetica,Arial,sans-serif;'+
        'padding:24px;line-height:1.45;';

    // ── Title page ─────────────────────────────────────────────
    var title = document.createElement('div');
    title.style.cssText = 'text-align:center;padding:40px 0 20px';
    title.innerHTML =
        '<h1 style="font-size:30px;margin:0 0 10px;color:#111">'+ {project_name_js} +'</h1>'+
        '<p style="font-size:15px;color:#444;margin:0">OpenWQ Calibration Results Report</p>'+
        '<p style="font-size:12px;color:#666;margin-top:14px">Generated: '+
        new Date().toLocaleString() +'</p>'+
        '<hr style="margin-top:30px;border:0;border-top:1px solid #ccc"/>';
    wrap.appendChild(title);

    // ── Clone each section, force its background, prepend page-break ─
    SECTIONS.forEach(function(sid, idx) {{
        var src = document.getElementById(sid);
        if(!src) return;
        var clone = src.cloneNode(true);
        clone.classList.add('pdf-page-break');
        // Strip elements that don't render in the PDF (copy buttons,
        // floating tooltips, code-block expand/collapse triangles).
        clone.querySelectorAll('.copy-btn, .no-pdf, button').forEach(function(el){{
            el.parentNode && el.parentNode.removeChild(el);
        }});
        // Some plot images are emitted inside <details><summary>; force-open
        // them so they appear in the captured DOM.
        clone.querySelectorAll('details').forEach(function(el){{ el.open = true; }});
        wrap.appendChild(clone);
    }});

    document.body.appendChild(wrap);

    // Give the cloned Leaflet tiles a moment to load before capture.
    // Map tiles in the original are already on screen; the clone shares
    // their <img> srcs so they should render almost immediately, but a
    // 250ms grace covers slow CDN responses on the satellite tiles.
    setTimeout(function() {{
        var opt = {{
            margin: [10, 10, 12, 10],
            filename: ('openwq_calibration_results_'+
                       new Date().toISOString().slice(0,10) +'.pdf'),
            image:    {{ type: 'jpeg', quality: 0.92 }},
            html2canvas: {{
                scale: 2,
                useCORS: true,
                backgroundColor: '#ffffff',
                logging: false,
                allowTaint: false,
                ignoreElements: function(el) {{
                    return el.id === 'pdfExportBtn' ||
                           (el.classList && el.classList.contains('no-pdf'));
                }}
            }},
            jsPDF: {{ unit: 'mm', format: 'a4', orientation: 'portrait' }},
            pagebreak: {{ mode: ['css'], before: '.pdf-page-break' }}
        }};
        html2pdf().set(opt).from(wrap).save().then(function() {{
            wrap.remove();
            btn.disabled = false; btn.classList.remove('busy');
            btn.querySelector('.label').textContent = ' \\u2705 Saved \\u2014 Download PDF';
            setTimeout(function() {{
                btn.querySelector('.label').textContent = ' \\u1F4C4 Download PDF';
            }}, 2500);
        }}).catch(function(err) {{
            console.error('PDF export failed:', err);
            wrap.remove();
            btn.disabled = false; btn.classList.remove('busy');
            btn.querySelector('.label').textContent = ' \\u26A0 Export failed \\u2014 retry';
            alert('PDF export failed: '+(err && err.message ? err.message : err));
        }});
    }}, 250);
}}
</script>
"""
