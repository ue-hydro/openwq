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


# Failure penalty returned by the objective function for any evaluation that
# could not be scored (no model output, no obs-sim overlap, etc.).  Must match
# the ``1e10`` used in objective_functions.py.  When the optimiser's BEST score
# is still this penalty, NO evaluation ever succeeded — converting it back to a
# metric would read as "Best KGE = 1 - 1e10 = -9999999999", which is meaningless.
_FAIL_PENALTY = 1e10


def _no_successful_eval(best_obj) -> bool:
    """True when the best objective is missing/NaN or still the failure
    penalty — i.e. the calibration produced no scorable evaluation."""
    if not isinstance(best_obj, (int, float)):
        return True
    if best_obj != best_obj:  # NaN
        return True
    return best_obj >= _FAIL_PENALTY * 0.9


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
    in_progress: bool = False,
    defer_heavy_plots: bool = False,
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
        # Plotly.js + shared theme-aware chart helpers + per-figure Lock-zoom /
        # branded PDF export.  Pass the project name so exported PDFs are stamped.
        H.append(_plotly_bootstrap(model_config.get("project_name", "OpenWQ Calibration")))
        H.append('<div class="layout">')

        # ── Determine which of the two parts actually ran ──
        _did_sens = bool(sensitivity_results)
        _did_calib = (bool(calibration_results.get("calibration_ran", True))
                      and bool(calibration_results.get("history")))
        # Hybrid physics-ML LAYER 2: closure diagnostics JSON (self-contained,
        # written by openWQ). Its presence == ML was active on the best fit.
        try:
            _ml_diag = _ml_diag_path(output_dir) if _did_calib else None
        except Exception:
            _ml_diag = None
        _mode = calibration_settings.get("calibration_mode") or (
            "both" if (_did_sens and _did_calib)
            else "sensitivity" if _did_sens else "calibration")

        # ── Sidebar (Part 1: influential params, Part 2: calibration) ──
        try:
            _bgc_path = _ci.get_bgc_template_path(model_config)
        except Exception:
            _bgc_path = None
        _has_bgc = bool(_bgc_path) and os.path.isfile(_bgc_path)
        nav_items = [{"id": "model-config", "label": "Model Setup"}]
        # Observation Map then BGC Reaction Network sit between Model Setup and
        # the rest (body order: Model Config → Obs Map → BGC → Output Config).
        if _did_calib and _has_observation_map_inputs(model_config):
            nav_items.append({"id": "obs-map", "label": "Observation Map"})
        if _has_bgc:
            nav_items.append({"id": "bgc-network", "label": "BGC Network"})
        nav_items.append({"id": "summary", "label": "Summary"})
        # Performance + Convergence sit right after the Summary (body order).
        if performance_metrics:
            nav_items.append({"id": "performance", "label": "Performance"})
        if _did_calib:
            nav_items.append({"id": "convergence", "label": "Convergence"})
        nav_items.append({"id": "sensitivity", "label": "Influential params"})
        if _did_calib:
            nav_items.append({"id": "best-params", "label": "Best Parameters"})
            if _ml_diag:
                nav_items.append({"id": "ml-closures", "label": "ML Closures"})
            nav_items.extend([
                {"id": "param-evolution", "label": "Parameter Evolution"},
                {"id": "param-correlations", "label": "Correlations"},
            ])
        else:
            nav_items.append({"id": "best-params", "label": "Calibration"})
        if _did_calib:
            # The time-series section is always rendered for a calibration run
            # (when no obs-sim pairs exist it shows a placeholder explaining
            # why), so always provide its nav link.
            nav_items.append({"id": "timeseries",
                               "label": "Observed vs Simulated: calibration"})
            nav_items.append({"id": "valid-ts",
                               "label": "Observed vs Simulated: Best Sim "
                                        "(Calibration &amp; Validation)"})
            nav_items.append({"id": "next-steps", "label": "Next Steps"})
            # "Run Best Params" snippet lives at the very bottom of the report.
            nav_items.append({"id": "run-best", "label": "Run Best Params"})

        H.append(rh.build_sidebar(nav_items, logo_text="OpenWQ Calibration"))

        H.append('<div class="main">')

        # ── Header ──
        project_name = model_config.get("project_name", "Calibration")
        authors = model_config.get("authors", "")
        date_str = datetime.now().strftime("%Y-%m-%d %H:%M")

        best_obj = calibration_results.get("best_objective", float('nan'))
        n_evals = calibration_results.get("n_evaluations", 0)

        _mode_label = {"sensitivity": "Influential parameters only",
                       "both": "Influential parameters + calibration",
                       "calibration": "Calibration only"}.get(_mode, "")
        meta_items = [
            f"Authors: {authors}" if authors else "",
            f"Completed: {date_str}",
            f"Workflow: {_mode_label}" if _mode_label else "",
        ]
        if _did_calib:
            _ofn = calibration_settings.get('objective_function', 'KGE')
            # Header "Best {metric}" = the BEST per-species score (matching the
            # Summary tile), so a poorly-observed species can't drag the headline
            # below the well-fit target species.  Falls back to the optimiser's
            # best objective (metric = 1 - objective for KGE/NSE) when no
            # per-species metrics exist.  "n/a" when no evaluation succeeded.
            if _no_successful_eval(best_obj):
                _bobj = "n/a (no successful evaluations)"
            else:
                _mk = _ofn.upper()
                _sv = []
                for _r in (performance_metrics or []):
                    _v = _r.get(_mk)
                    if isinstance(_v, (int, float)) and not (isinstance(_v, float) and _v != _v):
                        _sv.append(abs(_v) if _mk == "PBIAS" else float(_v))
                if _sv:
                    _hdr_metric = min(_sv) if _mk in ("RMSE", "PBIAS") else max(_sv)
                elif _ofn in ("KGE", "NSE"):
                    _hdr_metric = 1.0 - best_obj
                else:
                    _hdr_metric = best_obj
                _bobj = (f"{_hdr_metric:.4f}"
                         if isinstance(_hdr_metric, (int, float)) and _hdr_metric == _hdr_metric
                         else "N/A")
            meta_items.append(f"Best {_ofn}: {_bobj}")
            meta_items.append(f"Evaluations: {n_evals}")
        meta_items = [m for m in meta_items if m]

        H.append(rh.build_header(
            title="Calibration Results Report",
            subtitle=project_name,
            meta_items=meta_items,
            badge_text="RESULTS",
            badge_class="badge-secondary"
        ))

        H.append('<div class="container">')

        # ── In-progress banner (partial report opened mid-run) ──
        if in_progress:
            _n_so_far = calibration_results.get("n_evaluations", 0)
            _sens_prog = calibration_results.get("sensitivity_progress")
            _sens_line = (
                f" Influential-parameter screening: "
                f"<strong>{_sens_prog}</strong> model runs completed."
                if _sens_prog else "")
            H.append(f"""
<div class="section" id="in-progress">
    <div class="card" style="border-left:5px solid #3b82f6;
         background:rgba(59,130,246,.08);">
        <p style="margin:0;font-size:.95rem;">
            &#9203; <strong>In progress.</strong> This is a partial report
            generated from the latest results available right now
            ({_n_so_far} calibration evaluation{'s' if _n_so_far != 1 else ''} so far).{_sens_line}
            Re-run the <code>--report</code> command (or reopen this file) any
            time for a fresher snapshot, and once more after the run finishes
            for the complete report. (Sections still being computed may be
            empty or show a placeholder.)
        </p>
    </div>
</div>
""")

        # ── Warning: a targeted species has NO primary observations ──
        # The objective function scores only PRIMARY observations (the basin
        # outlet / pouring-point station).  If a species set as a calibration
        # target has none, it cannot be matched or calibrated (e.g. SUMMA: the
        # outlet station simply has no data for that species).  Surface it so the
        # user isn't left wondering why a targeted species is absent from the
        # metric and the observed-vs-simulated plots.
        try:
            _targets = [str(s) for s in
                        ((calibration_settings or {}).get("calibration_targets") or {})
                        .get("species", [])]
        except Exception:
            _targets = []
        if _did_calib and _targets:
            _species_with_primary = set()
            _obs_csv = os.path.join(output_dir, "calibration_observations.csv")
            if os.path.isfile(_obs_csv):
                try:
                    import pandas as _pd
                    _odf = _pd.read_csv(_obs_csv)
                    if "species" in _odf.columns:
                        if "is_primary" in _odf.columns:
                            _odf = _odf[_odf["is_primary"].fillna(False).astype(bool)]
                        _species_with_primary = {str(s) for s in _odf["species"].unique()}
                except Exception:
                    _species_with_primary = set()
            _no_primary = [s for s in _targets if s not in _species_with_primary]
            # Only warn when we could actually read the obs AND some species DO
            # have primary obs — i.e. this is "these species specifically", not a
            # totally empty/failed observation set (covered elsewhere).
            if _no_primary and _species_with_primary:
                _n = len(_no_primary)
                _ws = ", ".join(f"<strong>{html_lib.escape(s)}</strong>"
                                for s in _no_primary)
                _subj = "These species were" if _n > 1 else "This species was"
                _tgt = "calibration targets" if _n > 1 else "a calibration target"
                _it = "them" if _n > 1 else "it"
                _itwas = "they were" if _n > 1 else "it was"
                _itdoes = "they do" if _n > 1 else "it does"
                _have = ", ".join(html_lib.escape(s)
                                  for s in sorted(_species_with_primary)) or "none"
                _objfn = str((calibration_settings or {}).get(
                    "objective_function", "objective-function"))
                H.append(f"""
<div class="section" id="no-primary-obs">
    <div class="card" style="border-left:5px solid #ff6b35;
         background:rgba(255,107,53,.08);">
        <p style="margin:0;font-size:.95rem;line-height:1.6;">
            &#9888;&#65039; <strong>No primary observations for {_ws}.</strong>
            {_subj} set as {_tgt} for the objective function, but the metric only
            uses <em>primary</em> observations (the basin outlet / pouring-point
            station) and that station has no data for {_it}. So {_itwas}
            <strong>not matched or calibrated</strong>, and {_itdoes} not appear
            in the {html_lib.escape(_objfn)} metric or the
            observed-vs-simulated plots. Species that DO have primary
            observations: <strong>{_have}</strong>.
        </p>
    </div>
</div>
""")

        # ── Section: Model Configuration (general model setup) ──
        # Reuse the calibration setup report's "Module Summary" builder so the
        # results report stands on its own — the reader sees the host model,
        # solver, modules and Source/Sink climate detail up front, without
        # needing to cross-reference the config report.  Defensive: never let a
        # missing key break the results report.
        try:
            from .Gen_Calibration_Setup_Report import (
                _build_model_config_section, _build_output_config_section)
            H.append(_build_model_config_section(
                model_config,
                config_template_path=(model_config.get("_model_config_path")
                                      if isinstance(model_config, dict) else None)))
            # ── Observation Map (between Model Configuration and BGC) ──
            if _did_calib:
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
            # ── BGC Reaction Network (between Model & Output configuration) ──
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
            H.append(_build_output_config_section(model_config))
        except Exception as _e:
            logger.warning(f"Model configuration summary skipped: {_e}")

        # ── Section: Summary ──
        H.append(_build_summary_section(
            calibration_results, calibration_settings, calibration_parameters,
            performance_metrics=performance_metrics,
        ))

        # ── Section: Model Performance (placed right after the Summary) ──
        if performance_metrics:
            H.append(_build_performance_section(
                performance_metrics, matched_data,
                calibration_settings, output_dir,
                hostmodel=(model_config.get("hostmodel") or "mizuroute"),
            ))

        # ── Section: Convergence (placed right after Model Performance) ──
        if _did_calib:
            H.append(_build_convergence_section(
                calibration_results, calibration_settings, output_dir))

        # ═══════════════════════════════════════════════════════════════
        # PART 1 — Influential parameters (sensitivity screening)
        # ═══════════════════════════════════════════════════════════════
        H.append(_build_part_divider(
            "Part 1 &mdash; Influential parameters",
            "Morris / Sobol screening: which parameters most move the objective."))
        if _did_sens:
            H.append(_build_sensitivity_section(sensitivity_results, output_dir))
        elif in_progress and _mode in ("sensitivity", "both"):
            _sp = calibration_results.get("sensitivity_progress")
            _msg = ("Influential-parameter screening is running"
                    + (f" &mdash; {_sp} model runs completed so far."
                       if _sp else
                       " &mdash; results will appear here once it finishes.")
                    + " Morris/Sobol indices are only computed after all "
                    "screening runs complete; reopen this report then.")
            H.append(_build_not_run_notice(
                "sensitivity", "Influential parameters", _msg,
                status="&#9203; In progress.", colour="#3b82f6"))
        else:
            H.append(_build_not_run_notice(
                "sensitivity", "Influential parameters",
                f"Sensitivity screening was not run for this workflow "
                f"(mode = &lsquo;{_mode}&rsquo;). Set the workflow slider to "
                "&lsquo;Influential parameters&rsquo; or &lsquo;Both&rsquo; to compute it."))

        # ═══════════════════════════════════════════════════════════════
        # PART 2 — Calibration
        # ═══════════════════════════════════════════════════════════════
        H.append(_build_part_divider(
            "Part 2 &mdash; Calibration",
            "Optimised parameters, convergence, and model performance."))
        if _did_calib:
            H.append(_build_best_params_section(
                calibration_parameters, calibration_results))
            # Hybrid physics-ML LAYER 2 — closure correction ("how much the ML
            # pulled"), rendered only when openWQ wrote the diagnostics JSON.
            if _ml_diag:
                H.append(_build_ml_closure_section(_ml_diag))
            H.append(_build_param_evolution_section(
                calibration_results, calibration_parameters, output_dir))
            H.append(_build_correlations_section(
                calibration_results, calibration_parameters, output_dir))
            # Observed vs simulated time series at the best parameters, per
            # calibrated species (driven by matched_data; independent of the
            # numeric performance-metrics table above).  The "Run with Best
            # Parameters" snippet is appended at the very END of the report.
            if defer_heavy_plots:
                # Light / interim mode (per-evaluation refresh during a long
                # calibration).  The observed-vs-simulated section reads and
                # processes the FULL per-reach simulated series (the report's
                # heaviest work, and it grows with stations x evaluations), so we
                # skip it here and point the user at the on-demand `--report`
                # (or the automatic end-of-run report), which builds it once.
                H.append(_build_not_run_notice(
                    "best-params", "Observed vs simulated",
                    "Full observed-vs-simulated time series, scatter and "
                    "validation plots are generated on demand to keep the "
                    "calibration fast &mdash; run the calibration script with "
                    "<code>--report</code>, or wait for the automatic "
                    "end-of-run report, to see them.",
                    status="&#9889; Deferred for speed.", colour="#3b82f6"))
            else:
                try:
                    _ts_section = _build_calibrated_timeseries_section(
                        matched_data, calibration_settings, model_config,
                        output_dir, n_evals=n_evals)
                except Exception as _e:
                    logger.warning(f"Time-series section skipped: {_e}")
                    _ts_section = ""
                if _ts_section:
                    H.append(_ts_section)
                # Combined calibration + validation best-fit (best params re-run
                # over the validation period; produced automatically when a run
                # completes).  Shows a placeholder until that data exists.
                try:
                    H.append(_build_validation_section(
                        output_dir, calibration_settings, model_config))
                except Exception as _e:
                    logger.warning(f"Validation section skipped: {_e}")
        elif in_progress and _mode in ("both", "calibration"):
            H.append(_build_not_run_notice(
                "best-params", "Calibration",
                "Calibration is in progress &mdash; results will appear here as "
                "evaluations complete. Reopen this report shortly.",
                status="&#9203; In progress.", colour="#3b82f6"))
        else:
            H.append(_build_not_run_notice(
                "best-params", "Calibration",
                f"Calibration was not run for this workflow "
                f"(mode = &lsquo;{_mode}&rsquo; &mdash; influential parameters only). "
                "Set the workflow slider to &lsquo;Calibration&rsquo; or "
                "&lsquo;Both&rsquo; to calibrate."))

        # ── Section: Run with Best Parameters (kept LAST in the report, per
        # request — it's the actionable "what to run next" snippet) ──
        if _did_calib:
            H.append(_build_next_steps_section(
                model_config, calibration_settings))
            H.append(_build_run_best_section(
                calibration_results, model_config, calibration_settings,
                calibration_parameters=calibration_parameters,
                output_dir=output_dir,
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


def _build_part_divider(title: str, subtitle: str = "") -> str:
    """A bold banner that separates the two report parts (influential
    parameters vs calibration)."""
    sub = (f'<p style="margin:.3rem 0 0;color:var(--text2);font-size:.92rem;">'
           f'{subtitle}</p>' if subtitle else "")
    return f"""
<div class="section part-divider">
    <h2 style="margin-bottom:0;border-bottom:3px solid var(--primary);
               padding-bottom:.45rem;font-size:1.5rem;">{title}</h2>
    {sub}
</div>
"""


def _build_not_run_notice(section_id: str, title: str, message: str,
                          status: str = "&#9888; Not run.",
                          colour: str = "#f59e0b") -> str:
    """A notice card shown in a part that was skipped (or is still running)."""
    return f"""
<div class="section" id="{section_id}">
    <h2>{title}</h2>
    <div class="card" style="border-left:5px solid {colour};
         background:{colour}14;">
        <p style="margin:0;font-size:.95rem;">
            <strong>{status}</strong> {message}
        </p>
    </div>
</div>
"""


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

    /* ── ML closure diagnostics ("how much the ML pulled") ───────── */
    .mlc-grid { display: grid; gap: 1rem;
        grid-template-columns: repeat(auto-fit, minmax(280px, 1fr)); }
    .mlc-card { border: 1px solid var(--border); border-radius: 12px;
        padding: 1rem 1.1rem 1.1rem; background: var(--surface); }
    .mlc-head { display: flex; align-items: baseline; justify-content: space-between;
        gap: .5rem; margin-bottom: .5rem; }
    .mlc-title { font-weight: 700; font-size: 1.05rem; }
    .mlc-badge { font-family: 'JetBrains Mono', monospace; font-size: .68rem;
        text-transform: uppercase; letter-spacing: .6px; color: var(--text2);
        border: 1px solid var(--border2); border-radius: 999px; padding: .12rem .55rem;
        white-space: nowrap; }
    .mlc-big { font-size: 2rem; font-weight: 800; line-height: 1;
        font-family: 'JetBrains Mono', monospace; }
    .mlc-amp  { color: var(--secondary); }
    .mlc-damp { color: var(--accent); }
    .mlc-caption { font-size: .82rem; color: var(--text2); margin: .2rem 0 .1rem; }
    /* correction gauge: track from (1-max) .. (1+max), band = observed range */
    .mlc-gauge { position: relative; height: 26px; margin: .9rem 0 .25rem;
        border-radius: 6px; border: 1px solid var(--border);
        background: linear-gradient(90deg, rgba(255,107,53,.14), rgba(128,128,128,.10) 50%, rgba(0,168,107,.16)); }
    .mlc-band { position: absolute; top: 3px; bottom: 3px; border-radius: 3px;
        background: var(--primary); opacity: .30; }
    .mlc-one  { position: absolute; top: -3px; bottom: -3px; width: 2px;
        background: var(--text); opacity: .75; }
    .mlc-mean { position: absolute; top: 50%; width: 12px; height: 12px;
        border-radius: 50%; background: var(--primary);
        border: 2px solid var(--surface); transform: translate(-50%, -50%); }
    .mlc-scale { display: flex; justify-content: space-between; font-size: .66rem;
        color: var(--text3); font-family: 'JetBrains Mono', monospace; margin-bottom: .45rem; }
    .mlc-scale .lo { color: var(--accent); } .mlc-scale .hi { color: var(--secondary); }
    .mlc-legend { display: flex; flex-wrap: wrap; gap: .1rem .9rem; font-size: .68rem;
        color: var(--text2); margin-bottom: .3rem; align-items: center; }
    .mlc-legend span.item { display: inline-flex; align-items: center; gap: .3rem; }
    .mlc-legend .sw-line { width: 2px; height: 12px; background: var(--text); opacity: .75; }
    .mlc-legend .sw-band { width: 14px; height: 10px; border-radius: 2px;
        background: var(--primary); opacity: .30; }
    .mlc-legend .sw-mean { width: 11px; height: 11px; border-radius: 50%;
        background: var(--primary); border: 2px solid var(--surface); }
    .mlc-stats { display: grid; grid-template-columns: auto 1fr; gap: .32rem .8rem;
        font-size: .82rem; align-items: baseline; }
    .mlc-stats .k { color: var(--text2); }
    .mlc-stats .v { text-align: right; font-family: 'JetBrains Mono', monospace;
        font-weight: 600; }
    .mlc-note { font-size: .82rem; color: var(--text2); line-height: 1.5;
        margin-top: .8rem; padding-top: .7rem; border-top: 1px dashed var(--border); }

    /* ── Run-with-best "save template" toolbar ─────────────────── */
    .rb-toolbar { display: flex; flex-wrap: wrap; align-items: center;
        gap: .55rem .9rem; margin: .3rem 0 .9rem; }
    .rb-save-btn { display: inline-flex; align-items: center; gap: .45rem;
        cursor: pointer; font: inherit; font-size: .85rem; font-weight: 600;
        color: #fff; background: var(--primary); border: none; border-radius: 8px;
        padding: .5rem .95rem; }
    .rb-save-btn:hover { filter: brightness(1.08); }
    .rb-path { display: flex; flex-wrap: wrap; align-items: center; gap: .4rem;
        font-size: .8rem; color: var(--text2); }
    .rb-path .p { font-family: 'JetBrains Mono', monospace; background: rgba(0,0,0,.05);
        padding: .25rem .55rem; border-radius: 5px; word-break: break-all; color: var(--text); }
    [data-theme="dark"] .rb-path .p { background: rgba(255,255,255,.06); }
    .rb-copy-path { cursor: pointer; font: inherit; font-size: .78rem; color: var(--text2);
        border: 1px solid var(--border); background: var(--surface); border-radius: 6px;
        padding: .28rem .6rem; white-space: nowrap; }
    .rb-copy-path:hover { border-color: var(--primary); color: var(--primary); }
    .rb-hint { font-size: .78rem; color: var(--text3); margin: -.3rem 0 .9rem; line-height: 1.5; }
    .rb-hint kbd { font-family: 'JetBrains Mono', monospace; font-size: .72rem;
        background: rgba(0,0,0,.06); border: 1px solid var(--border); border-radius: 4px;
        padding: .05rem .35rem; }
    [data-theme="dark"] .rb-hint kbd { background: rgba(255,255,255,.08); }
    /* collapsible code box: top visible, scroll for the rest, expandable */
    .rb-code-controls { display: flex; align-items: center; gap: .6rem; margin: .2rem 0 .3rem; }
    .rb-expand { cursor: pointer; font: inherit; font-size: .8rem; font-weight: 600;
        color: var(--text2); border: 1px solid var(--border); background: var(--surface);
        border-radius: 6px; padding: .3rem .7rem; }
    .rb-expand:hover { border-color: var(--primary); color: var(--primary); }
    .rb-code-meta { font-size: .74rem; color: var(--text3); }
    .rb-codewrap { position: relative; max-height: 11rem; overflow: auto;
        overscroll-behavior: contain; border: 1px solid var(--border); border-radius: 8px; }
    .rb-codewrap.expanded { max-height: none; }
    /* locked: box stops capturing the wheel so the PAGE scrolls past it */
    .rb-codewrap.locked { overflow: hidden; }
    .rb-codewrap .code-block { margin: 0; border-radius: 0; overflow-x: auto; }
    /* fade hint at the bottom while collapsed */
    .rb-codewrap:not(.expanded)::after { content: ""; position: sticky; left: 0; bottom: 0;
        display: block; height: 2.4rem; margin-top: -2.4rem; pointer-events: none;
        background: linear-gradient(to bottom, transparent, var(--surface)); }
    .rb-expand.locked-on { border-color: var(--primary); color: var(--primary); }
    .rb-warn { font-size: .85rem; line-height: 1.5; border-left: 4px solid var(--accent);
        background: rgba(255,107,53,.08); border-radius: 6px; padding: .55rem .8rem; }
    /* output-report interactive controls (tickboxes) */
    .viz-ctrl { display: flex; flex-direction: column; gap: .5rem; margin: .2rem 0 .6rem; }
    .viz-group { display: flex; flex-wrap: wrap; align-items: center; gap: .4rem .7rem; }
    .viz-label { font-size: .74rem; font-weight: 700; text-transform: uppercase;
        letter-spacing: .5px; color: var(--text2); min-width: 6.5rem; }
    .viz-cb { display: inline-flex; align-items: center; gap: .35rem; cursor: pointer;
        font-size: .84rem; font-family: 'JetBrains Mono', monospace;
        border: 1px solid var(--border); border-radius: 6px; padding: .2rem .5rem;
        background: var(--surface); }
    .viz-cb input { accent-color: var(--primary); cursor: pointer; margin: 0; }
    .viz-in { font: inherit; font-size: .82rem; font-family: 'JetBrains Mono', monospace;
        border: 1px solid var(--border); border-radius: 6px; padding: .25rem .5rem;
        background: var(--surface); color: var(--text); width: 10rem; }
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

    # best_obj is the value the optimiser MINIMISES.  For KGE/NSE that is
    # (1 - metric), so convert it back to the ACTUAL metric for both the KPI
    # value and the quality label.  Without this, a poor KGE of e.g. -5.7 is
    # stored as objective 6.7 and was being mislabelled "Very Good" (and shown
    # as "Best KGE = 6.7340", which is impossible — KGE ≤ 1).  For RMSE the
    # minimised objective already IS the metric, so no conversion.  This also
    # puts "Best" on the same raw-metric scale as the Mean/Median KPIs below.
    # When NO evaluation succeeded the optimiser's "best" is still the failure
    # penalty (1e10).  Converting that back would read as "Best KGE = 1 - 1e10
    # = -9999999999", which is meaningless, so flag it and show a clear
    # diagnostic instead of a fake number.
    no_success = _no_successful_eval(best_obj)
    # Headline "Best {metric}" = the BEST score achieved for ANY single species
    # (max KGE/NSE/R2, or min RMSE / |PBIAS|), taken from the per-species
    # performance metrics — NOT the multi-species objective the optimiser
    # minimises.  The objective is an (equal-weighted) mean across species, so a
    # species with very few / unfittable observations can drag it far below the
    # well-fit target species; the headline should instead show the best species'
    # fit.  Falls back to the optimiser's best objective when no per-species
    # metrics are available (e.g. sensitivity-only runs).
    _mk = obj_fn.upper()
    _spec_vals = []
    if performance_metrics:
        for _row in performance_metrics:
            _v = _row.get(_mk)
            if isinstance(_v, (int, float)) and not (isinstance(_v, float) and _v != _v):
                _spec_vals.append(abs(_v) if _mk == "PBIAS" else float(_v))
    _lower_is_better = _mk in ("RMSE", "PBIAS")
    if no_success:
        best_metric = None
        best_value_str = "n/a"
        obj_quality = None
    elif _spec_vals:
        # best across species (best-fit species, not the mean objective)
        best_metric = min(_spec_vals) if _lower_is_better else max(_spec_vals)
        best_value_str = f"{best_metric:.4f}"
        obj_quality = _assess_objective(best_metric, obj_fn)
    elif obj_fn in ("KGE", "NSE"):
        best_metric = 1.0 - best_obj
        best_value_str = f"{best_metric:.4f}"
        # Objective value quality assessment (best), on the real metric scale.
        obj_quality = _assess_objective(best_metric, obj_fn)
    else:
        best_metric = best_obj
        best_value_str = f"{best_metric:.4f}"
        obj_quality = _assess_objective(best_metric, obj_fn)

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
        {"icon": "\U0001f3c6", "value": best_value_str,
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

    # Quality assessment box (based on the BEST score) — or, when NO evaluation
    # ever succeeded, a prominent diagnostic explaining why there is no score
    # and no observed-vs-simulated time series (instead of a fake metric).
    if no_success:
        quality_html = rh.build_highlight_box(
            f"<strong>&#9888;&#65039; No successful evaluations.</strong> "
            f"All <strong>{n_evals}</strong> model evaluation(s) returned the "
            f"failure penalty, so no <strong>{obj_fn}</strong> score could be "
            f"computed and there are <strong>no observed-vs-simulated pairs to "
            f"plot</strong>. This means the model did not produce output that "
            f"overlaps the observations &mdash; <em>not</em> that the fit was "
            f"poor. Common causes: the model / container did not run to "
            f"completion (check the SLURM <code>.out</code> log and "
            f"<code>calibration.log</code>), the simulated output file or "
            f"variable was not found, or the simulated period / reach does not "
            f"overlap the observation dates. Fix the underlying model run, then "
            f"re-run the calibration.",
            "warning"
        )
    else:
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
    <details class="module-details">
        <summary>Show / hide the full parameter table</summary>
        {table_html}
    </details>
</div>
"""


def _ml_diag_path(output_dir):
    """Path to the Layer-2 closure-diagnostics JSON written by openWQ for the
    GLOBAL-best evaluation (``<eval>/openwq_out/**/ml_closure_diagnostics.json``),
    or ``None`` when no closures ran / the engine predates the feature.  The
    file is fully self-contained, so the report needs nothing from the Python
    calibration side to render the ML section."""
    import os as _os, re as _re, glob as _glob
    n = _best_eval_number(output_dir)
    evdir = None
    if n is not None:
        for d in _glob.glob(_os.path.join(str(output_dir), "evaluations", "eval_*")):
            m = _re.search(r"eval_(\d+)", _os.path.basename(d))
            if m and int(m.group(1)) == n:
                evdir = d
                break
    search_roots = ([evdir] if evdir else []) + [
        d for d in _glob.glob(_os.path.join(str(output_dir), "evaluations", "eval_*"))
    ]
    for d in search_roots:
        hits = _glob.glob(_os.path.join(d, "openwq_out", "*", "ml_closure_diagnostics.json")) \
             + _glob.glob(_os.path.join(d, "openwq_out", "ml_closure_diagnostics.json"))
        if hits:
            return hits[0]
    return None


def _build_ml_closure_section(diag_path):
    """Hybrid physics-ML LAYER 2 — "how much the ML pulled".  Renders the exact,
    solver-side closure diagnostics openWQ accumulated during the best-fit run:
    the multiplicative factor each closure applied to a process tendency
    (factor = 1 -> pure physics, |factor-1| = the "pull"), the flux-weighted
    correction, the net vs gross mass it moved, and whether it saturated. Only
    called when the JSON exists (i.e. ML was actually active)."""
    import json as _json, os as _os
    if not diag_path or not _os.path.isfile(diag_path):
        return ""
    try:
        with open(diag_path) as _f:
            diag = _json.load(_f)
    except Exception as _e:
        logger.warning(f"ML closure diagnostics unreadable: {_e}")
        return ""
    closures = [c for c in (diag.get("closures") or [])
                if int(c.get("n_samples", 0) or 0) > 0]
    if not closures:
        return ""

    units = diag.get("mass_units", "g") or "g"
    solver = diag.get("solver", "")

    def _pct(x):
        try: return f"{float(x) * 100:.1f}%"
        except Exception: return "n/a"

    def _mass(v):
        # v is already in the model's output mass unit; auto-scale for legibility
        try: v = float(v)
        except Exception: return "n/a"
        a = abs(v)
        if a >= 1e6:  return f"{v/1e6:.2f}×10⁶ {units}"
        if a >= 1e3:  return f"{v/1e3:.2f}×10³ {units}"
        if a < 1e-3 and a > 0: return f"{v:.2e} {units}"
        return f"{v:.3g} {units}"

    cards = []
    tot_gross = 0.0
    worst_pull = 0.0
    for c in closures:
        sp    = str(c.get("species", "?"))
        term  = str(c.get("term", "?"))
        comp  = str(c.get("compartment", "?"))
        alpha = float(c.get("alpha", 0.0) or 0.0)
        maxc  = float(c.get("max_correction", 0.0) or 0.0)
        fmean = float(c.get("factor_mean", 1.0) or 1.0)
        fmin  = float(c.get("factor_min", 1.0) or 1.0)
        fmax  = float(c.get("factor_max", 1.0) or 1.0)
        pull_m  = float(c.get("pull_mean", 0.0) or 0.0)
        pull_x  = float(c.get("pull_max", 0.0) or 0.0)
        pull_fw = float(c.get("pull_flux_weighted", 0.0) or 0.0)
        net_fc  = float(c.get("net_flux_change", 0.0) or 0.0)
        clampf  = float(c.get("clamp_fraction", 0.0) or 0.0)
        nsamp   = int(c.get("n_samples", 0) or 0)
        gross   = float(c.get("gross_mass_moved", 0.0) or 0.0)
        net     = float(c.get("net_mass_moved", 0.0) or 0.0)
        tot_gross += abs(gross)
        worst_pull = max(worst_pull, pull_fw)

        amp = fmean >= 1.0
        dir_word = "amplified" if fmean > 1.0005 else ("damped" if fmean < 0.9995 else "left unchanged")
        dir_cls = "mlc-amp" if amp else "mlc-damp"

        # gauge geometry: track [1-max .. 1+max]
        lo = 1.0 - maxc if maxc > 0 else min(0.5, fmin)
        hi = 1.0 + maxc if maxc > 0 else max(1.5, fmax)
        span = (hi - lo) or 1.0
        def _p(x): return max(0.0, min(100.0, (x - lo) / span * 100.0))
        b_l, b_r = _p(min(fmin, fmax)), _p(max(fmin, fmax))
        b_w = max(1.5, b_r - b_l)

        conserved = abs(net_fc) < 5e-3
        net_disp = "≈0 (conserved)" if conserved else _pct(net_fc)
        netmass_disp = "≈0 (conserved)" if abs(net) < abs(gross) * 1e-3 else _mass(net)

        saturated = (fmax - fmin) < 0.01
        if saturated:
            shape = (f"The closure <strong>saturated</strong> — it applied a near-constant "
                     f"&times;{fmean:.2f} to the {term} tendency (a state-independent "
                     f"{_pct(pull_m)} rescale, i.e. equivalent to a constant rate multiplier).")
        else:
            shape = (f"<strong>State-dependent</strong>: the correction varied from "
                     f"&times;{fmin:.2f} to &times;{fmax:.2f} with the {sp} state.")
        cons_txt = ""
        if conserved and gross > 0:
            cons_txt = (f" Net mass moved &asymp; 0 &mdash; the closure redistributed mass "
                        f"<em>within</em> the reaction group without adding or removing any "
                        f"(mass-conserving); {_mass(gross)} was rescaled gross.")

        clamp_row = (f'<span class="k">Hit &plusmn;max bound</span>'
                     f'<span class="v">{_pct(clampf)}</span>') if clampf > 1e-9 else (
                     f'<span class="k">Hit &plusmn;max bound</span><span class="v">never</span>')

        cards.append(f"""
    <div class="mlc-card">
      <div class="mlc-head">
        <span class="mlc-title">{sp}</span>
        <span class="mlc-badge">{term} &middot; {comp}</span>
      </div>
      <div class="mlc-big {dir_cls}">&times;{fmean:.2f}</div>
      <div class="mlc-caption">mean factor on the {term} tendency &mdash; {dir_word} the physics</div>
      <div class="mlc-caption" style="margin-top:.7rem;">Correction applied vs pure physics
        (bounded to &plusmn;{maxc*100:.0f}% by <em>max</em>):</div>
      <div class="mlc-gauge" title="Left = strongest damping the closure may apply; centre = pure physics (no correction); right = strongest amplification. Shaded bar = range of factors actually applied; dot = average; vertical line = physics.">
        <div class="mlc-band" style="left:{b_l:.1f}%;width:{b_w:.1f}%;"></div>
        <div class="mlc-one"  style="left:{_p(1.0):.1f}%;"></div>
        <div class="mlc-mean" style="left:{_p(fmean):.1f}%;"></div>
      </div>
      <div class="mlc-scale">
        <span class="lo">&times;{lo:.2f}<br>max damping</span>
        <span>&times;1.00<br>physics</span>
        <span class="hi" style="text-align:right;">&times;{hi:.2f}<br>max amplifying</span>
      </div>
      <div class="mlc-legend">
        <span class="item"><span class="sw-line"></span>physics (&times;1.00)</span>
        <span class="item"><span class="sw-band"></span>range applied</span>
        <span class="item"><span class="sw-mean"></span>average (&times;{fmean:.2f})</span>
      </div>
      <div class="mlc-stats">
        <span class="k">Mean pull |f&minus;1|</span><span class="v">{_pct(pull_m)}</span>
        <span class="k">Peak pull</span><span class="v">{_pct(pull_x)}</span>
        <span class="k">Flux-weighted pull</span><span class="v">{_pct(pull_fw)}</span>
        <span class="k">Net flux change</span><span class="v">{net_disp}</span>
        {clamp_row}
        <span class="k">Gross mass moved</span><span class="v">{_mass(gross)}</span>
        <span class="k">Net mass moved</span><span class="v">{netmass_disp}</span>
        <span class="k">&alpha; / max</span><span class="v">{alpha:.3g} / {maxc:.3g}</span>
        <span class="k">Samples (cell&middot;step)</span><span class="v">{nsamp:,}</span>
      </div>
      <div class="mlc-note">{shape}{cons_txt}</div>
    </div>""")

    _solver_note = ("factor sampled once per step at the step's final state"
                    if solver == "sundials" else
                    "factor tracked exactly for every cell and step"
                    if solver == "forward_euler" else "")
    lead = (f"During the best-fit run the {len(closures)} active closure"
            f"{'s' if len(closures) != 1 else ''} corrected the physics by up to "
            f"<strong>{_pct(worst_pull)}</strong> (flux-weighted), moving a gross "
            f"<strong>{_mass(tot_gross)}</strong> in total.")

    return f"""
<div class="section" id="ml-closures">
    <h2>ML Closure Correction &mdash; how much the ML pulled</h2>
    <p style="color:var(--text2);margin-bottom:.4rem;max-width:70ch;">
        Each Layer-2 closure multiplies a process tendency by a bounded factor
        <code>1 + clamp(&alpha;&middot;g(state), &plusmn;max)</code>. A factor of
        <strong>1.0 is exact physics</strong>, so its departure from 1 &mdash;
        <code>|factor&minus;1|</code> &mdash; is exactly how hard the ML had to pull
        to improve the fit. These numbers are measured <em>inside the solver</em>
        during the best evaluation{(' (' + _solver_note + ')') if _solver_note else ''}.
    </p>
    <p style="color:var(--text2);margin-bottom:1rem;max-width:70ch;">{lead}</p>
    <div class="mlc-grid">{''.join(cards)}</div>
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
    # KGE/NSE are reported as (1 - minimised objective); error metrics pass
    # through.  Match the Summary so the stats + plot show the actual metric.
    _maximise = str(obj_fn).upper() in ("KGE", "NSE")

    def _obj_to_metric(o):
        if not isinstance(o, (int, float)) or o != o or o >= _FAIL_PENALTY * 0.9:
            return None
        return (1.0 - o) if _maximise else o

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
        # Convert each minimised objective back to the reported metric so
        # Best/Mean/Worst are on the SAME scale as the Summary KPIs and the
        # convergence plot (e.g. KGE = 1 - objective).  Failed evaluations
        # (penalty) become None and are excluded from the statistics.
        _metrics = [_obj_to_metric(h.get("objective", h.get("best_objective", 0)))
                    for h in history]
        _valid = [m for m in _metrics if m is not None]
        if _valid:
            best_m = max(_valid) if _maximise else min(_valid)
            worst_m = min(_valid) if _maximise else max(_valid)
            mean_m = sum(_valid) / len(_valid)
        else:
            best_m = worst_m = mean_m = None

        def _fmt_m(v):
            return (f"{v:.4f}" if isinstance(v, (int, float)) and v == v
                    else "n/a")

        _dir = "higher is better" if _maximise else "lower is better"
        stats_html = f"""
    <div class="card">
        <h3>{html_lib.escape(obj_fn)} statistics <span style="font-weight:400;font-size:.8rem;color:var(--text3);">({html_lib.escape(obj_fn)}; {_dir})</span></h3>
        <div class="perf-card" style="grid-template-columns: repeat(4, 1fr);">
            <div class="perf-metric">
                <div class="perf-value result-good">{_fmt_m(best_m)}</div>
                <div class="perf-label">Best</div>
            </div>
            <div class="perf-metric">
                <div class="perf-value">{_fmt_m(mean_m)}</div>
                <div class="perf-label">Mean</div>
            </div>
            <div class="perf-metric">
                <div class="perf-value result-poor">{_fmt_m(worst_m)}</div>
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

    # When the run used DDS parallel chains, note how the budget was split and
    # that each coloured path is one chain.
    _chain_ids = sorted({h["chain"] for h in history
                         if isinstance(h.get("chain"), int)})
    chains_note = ""
    if len(_chain_ids) > 1:
        _per = {}
        for h in history:
            c = h.get("chain")
            if isinstance(c, int):
                _per[c] = _per.get(c, 0) + 1
        _ev_each = min(_per.values()) if _per else 0
        chains_note = (
            f' <strong style="color:var(--text);">DDS parallel chains:</strong> '
            f'{len(_chain_ids)} independent chains (~{_ev_each} evaluations each)'
            f' &mdash; each coloured path is one chain; the best parameter set'
            f' is taken across all of them.')

    return f"""
<div class="section" id="param-evolution">
    <h2>Parameter Evolution</h2>
    <p style="color:var(--text2);margin-bottom:1rem;">
        How each parameter changed across calibration evaluations.{chains_note}
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

    # ── Normalise the saved format ────────────────────────────────────
    # The driver writes mu_star / sigma / S1 / ST as DICTs keyed by parameter
    # name, with no "method" / "parameter_names".  Re-shape into the
    # name-list + value-list form the table/plot builders expect.
    sa = dict(sa_results or {})
    method = (sa.get("method")
              or ("sobol" if ("S1" in sa or "ST" in sa) else "morris")).upper()
    sa["method"] = method.lower()

    param_names = sa.get("parameter_names")
    if not param_names:
        for k in ("mu_star", "ST", "S1", "mu"):
            if isinstance(sa.get(k), dict):
                param_names = list(sa[k].keys())
                break
        if not param_names and sa.get("rankings"):
            param_names = [r[0] for r in sa["rankings"]]
    param_names = param_names or []
    sa["parameter_names"] = param_names

    # Dict -> list aligned with param_names
    for k in ("mu", "mu_star", "sigma", "S1", "ST"):
        if isinstance(sa.get(k), dict):
            sa[k] = [sa[k].get(n, 0.0) for n in param_names]
    sa_results = sa

    if not param_names:
        return """
<div class="section" id="sensitivity">
    <h2>Influential parameters</h2>
    <div class="card"><p>No sensitivity results available.</p></div>
</div>
"""

    # Detect a DEGENERATE screening: every sensitivity index is ~0, so no
    # parameter shows any influence.  This happens when all screening
    # evaluations produced the SAME objective (typically the failure penalty -
    # i.e. the model runs failed / produced no comparable output, so changing a
    # parameter changed nothing).  It does NOT mean the parameters are
    # unimportant.  Surface a clear warning instead of a table of zeros that
    # reads as "nothing matters".
    _all_idx = []
    for _k in ("mu_star", "ST", "S1", "mu"):
        _v = sa.get(_k)
        if isinstance(_v, list):
            _all_idx.extend(x for x in _v if isinstance(x, (int, float)))
    _degenerate = bool(_all_idx) and all(abs(x) < 1e-12 for x in _all_idx)
    _n_sens = sa.get("n_evaluations") or 0
    _degenerate_html = ""
    if _degenerate:
        _degenerate_html = rh.build_highlight_box(
            "<strong>&#9888;&#65039; Screening could not rank any parameter.</strong> "
            f"Every {method} sensitivity index is zero across all "
            f"{len(param_names)} parameters &mdash; meaning the "
            f"{(str(_n_sens) + ' ') if _n_sens else ''}screening evaluation(s) all "
            "produced the <strong>same</strong> objective value. That is the "
            "signature of a <strong>failed screening run</strong> (e.g. every model "
            "run returned the failure penalty, so varying a parameter changed "
            "nothing). It does <em>not</em> mean the parameters are uninfluential. "
            "Re-run the influential-parameter screening on a working model, then "
            "regenerate this report.",
            "warning")

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
    <h2>Influential parameters &mdash; {method} screening</h2>
    {_degenerate_html}
    {plot_html}
    {table_html}
</div>
"""


def _collapse_subunit_dups(df, settings=None, value_col="simulated",
                           extra_keys=()):
    """Collapse multiple simulated values that share the same timestamp / reach
    / species into one.

    SUMMA writes several zoned HRU columns (``hruId_1_z1``, ``_z2`` ...) that
    all normalise to the same feature id, so a raw simulated series carries
    several values per timestamp for one reach/HRU.  Plotting them as one line
    makes it zig-zag between zones (looks "unstable").  The objective function
    averages those sub-units before scoring, so we do the same here — one value
    per timestamp, matching what the metric actually sees.  No-op when there
    are no duplicates (e.g. mizuRoute, one cell per reach)."""
    if df is None or getattr(df, "empty", True):
        return df
    keys = [k for k in (["datetime", "reach_id", "species"] + list(extra_keys))
            if k in getattr(df, "columns", [])]
    if value_col not in getattr(df, "columns", []) or not keys:
        return df
    try:
        if not df.duplicated(subset=keys).any():
            return df
        how = "median" if str((settings or {}).get(
            "aggregation_method", "mean")).lower() == "median" else "mean"
        agg = {value_col: how}
        for c in df.columns:
            if c not in keys and c != value_col:
                agg[c] = "first"
        return df.groupby(keys, as_index=False, sort=False).agg(agg)
    except Exception:
        return df


def _best_eval_number(output_dir):
    """Integer simulation number of the GLOBAL-best evaluation — the lowest
    successful, finite objective with HDF5 output across all ``eval_*`` folders
    — or ``None``.  Mirrors the driver's ``_best_eval_dir`` so the "best fit"
    legend can name the simulation it corresponds to, and that number matches
    the gray-overlay "sim N" labels."""
    import os as _os, re as _re, glob as _glob
    best_n, best_o = None, None
    for d in _glob.glob(_os.path.join(str(output_dir), "evaluations", "eval_*")):
        f = _os.path.join(d, "objective.txt")
        if not _os.path.isfile(f):
            continue
        try:
            t = open(f).read()
        except Exception:
            continue
        if "success: True" not in t:
            continue
        m = _re.search(r"objective:\s*([-0-9eE.+]+)", t)
        if not m:
            continue
        try:
            o = float(m.group(1))
        except ValueError:
            continue
        if o >= 1e10:                       # _REPORT_PENALTY → failed eval
            continue
        if not _glob.glob(_os.path.join(d, "openwq_out", "HDF5", "*.h5")):
            continue
        mn = _re.search(r"eval_(\d+)", _os.path.basename(d))
        if mn and (best_o is None or o < best_o):
            best_o, best_n = o, int(mn.group(1))
    return best_n


def _build_calibrated_timeseries_section(
    matched_data: Optional[Any],
    settings: Dict[str, Any],
    model_config: Dict[str, Any],
    output_dir: str,
    n_evals: int = 0,
) -> str:
    """Observed-vs-simulated time series (+ obs-vs-sim scatter) for each
    calibrated species, evaluated at the best-fit parameters.

    Driven purely by ``matched_data`` (the obs/sim pairs the objective
    function matched) so it renders even when the numeric performance-metrics
    table is unavailable.  Placed after "Run with Best Parameters".

    When there are no obs-sim pairs (e.g. every evaluation failed), it renders
    an explanatory placeholder rather than silently disappearing.
    """
    # Do we have any obs-sim pairs to plot?  (Check before the matplotlib guard
    # so the placeholder still shows when plotting libs are missing.)
    _has_pairs = False
    try:
        import pandas as pd
        _has_pairs = isinstance(matched_data, pd.DataFrame) and not matched_data.empty
    except ImportError:
        _has_pairs = (matched_data is not None
                      and getattr(matched_data, "empty", True) is False)

    # Load the persisted full simulated series so objective species with NO
    # observations at the target reach still get a simulated-only plot.
    _sim_all = None
    _all_sim = None
    _of_species = []
    _simonly = []
    try:
        import os as _os, pandas as _pd
        _of_species = [str(s) for s in
                       ((settings or {}).get("calibration_targets") or {}).get(
                           "species", [])]
        _sim_csv = _os.path.join(str(output_dir), "results",
                                 "simulated_data.csv")
        if _os.path.isfile(_sim_csv):
            _sim_all = _pd.read_csv(_sim_csv)
        # All evaluations' simulated series (for the light-gray overlay).
        _allsim_csv = _os.path.join(str(output_dir), "results",
                                    "all_simulated.csv")
        if _os.path.isfile(_allsim_csv):
            _all_sim = _pd.read_csv(_allsim_csv)
        # Collapse multi-zone HRU duplicates (SUMMA) so the simulated lines
        # don't zig-zag between sub-units at each timestamp; matches how the
        # objective averages sub-units.  No-op for one-cell-per-reach data.
        _sim_all = _collapse_subunit_dups(_sim_all, settings)
        _all_sim = _collapse_subunit_dups(
            _all_sim, settings, extra_keys=["eval_id", "chain"])
        _matched_sp = set()
        if _has_pairs and "species" in matched_data.columns:
            _matched_sp = {str(s) for s in matched_data["species"].unique()}
        if (_sim_all is not None and not _sim_all.empty
                and "species" in _sim_all.columns):
            _sim_sp = {str(s) for s in _sim_all["species"].unique()}
            _simonly = [s for s in (_of_species or sorted(_sim_sp))
                        if s not in _matched_sp and s in _sim_sp]
    except Exception:
        _sim_all = None

    # NOTE: the authoritative "best fit" is rebuilt upstream in
    # regenerate_results_report()/the driver from the global-best eval folder
    # (ground truth across all runs/resumes) and handed in via matched_data /
    # simulated_data.csv — so here we simply TRUST matched_data + _sim_all and
    # do not re-derive a (possibly worse) best from the sampled all_simulated.

    if not _has_pairs and not _simonly:
        _n = n_evals if isinstance(n_evals, int) and n_evals > 0 else 0
        _eval_txt = (f"across the {_n} evaluation(s)" if _n
                     else "in this run")
        return f"""
<div class="section" id="timeseries">
    <h2>Observed vs Simulated: calibration</h2>
    <div class="card" style="border-left:4px solid #ff6b35;">
        <p style="margin:.2rem 0;color:var(--text2);font-size:.92rem;line-height:1.6;">
            &#9888;&#65039; <strong>No observed-vs-simulated pairs are available,
            so there is nothing to plot here.</strong> The objective function
            matched zero observations to model output {_eval_txt} &mdash; usually
            because the model did not produce output overlapping the observation
            dates/locations (see the diagnostic in the <a href="#summary"
            style="color:var(--primary);font-weight:600;">Results&nbsp;Summary</a>
            above). Once at least one evaluation produces comparable output, this
            section will show the observed-vs-simulated time series and scatter for
            each calibrated species.
        </p>
    </div>
</div>
"""

    if not (HAS_MATPLOTLIB and HAS_NUMPY):
        return ""
    import pandas as pd

    hostmodel = (model_config.get("hostmodel") or "mizuroute")
    # Interactive (Plotly) obs-vs-sim time series + scatter, per species; plus
    # a simulated-only series for objective species with no obs at the reach.
    plots_html = _generate_timeseries_charts(
        matched_data, output_dir, hostmodel=hostmodel,
        simulated_data=_sim_all, of_species=_of_species, all_sim=_all_sim,
        best_eval=_best_eval_number(output_dir))
    if not plots_html:
        return ""

    _species = []
    try:
        if "species" in matched_data.columns:
            _species = [str(s) for s in matched_data["species"].unique()]
    except Exception:
        _species = []
    _sp_txt = (", ".join(f"<code>{s}</code>" for s in _species)
               if _species else "each calibrated species")

    return f"""
<div class="section" id="timeseries">
    <h2>Observed vs Simulated: calibration</h2>
    <div class="card" style="border-left:4px solid var(--secondary);">
        <p style="margin:.2rem 0;color:var(--text2);font-size:.92rem;line-height:1.6;">
            Time series of <strong>observed</strong> vs <strong>simulated</strong>
            concentrations at the best-fit parameters, for {_sp_txt}. Each block
            pairs the <strong>time series</strong> (left) with an
            <strong>observed-vs-simulated scatter</strong> (right) &mdash; points
            on the dashed 1:1 line indicate perfect agreement; points above it
            mean the model over-predicts, below it under-predicts. Per-element
            traces are collapsed into expandable panels below each overview.
        </p>
    </div>
    {plots_html}
</div>
"""


def _build_validation_section(output_dir, settings, model_config):
    """Combined calibration + validation best-fit time series.

    Reads ``results/validation_matched_data.csv`` (period-tagged obs/sim pairs
    produced when a calibration run completes by re-running the best parameters
    over the validation window).  Renders the full record per species with the
    calibration window shaded green on the left, the validation window shaded
    orange on the right, a dashed boundary line, and the objective metric for
    each period.  Shows a placeholder until that data exists."""
    import os
    _head = ('Observed vs Simulated: Best Sim '
             '(Calibration &amp; Validation)')
    vpath = os.path.join(str(output_dir), "results",
                         "validation_matched_data.csv")
    vd = None
    have = os.path.isfile(vpath)
    if have:
        try:
            import pandas as pd
            vd = pd.read_csv(vpath)
            have = (not vd.empty and "period" in vd.columns
                    and {"observed", "simulated", "species", "datetime"}
                    <= set(vd.columns))
        except Exception:
            have = False
    if not have:
        return (
            '\n<div class="section" id="valid-ts">\n'
            f'    <h2>{_head}</h2>\n'
            '    <div class="card"><p style="color:var(--text2);'
            'font-size:.92rem;line-height:1.6;">This section is produced '
            'automatically when a calibration run <strong>completes</strong>: '
            'the best parameters are re-run over the <strong>validation</strong> '
            'period (held-out data), and the full record is shown here &mdash; '
            'calibration on the left, validation on the right &mdash; with the '
            'objective metric for each period.</p></div>\n</div>\n')
    if not (HAS_MATPLOTLIB and HAS_NUMPY):
        return ""
    import pandas as pd
    import numpy as np
    from .objective_functions import ObjectiveFunction as _OF
    # Collapse any multi-zone HRU duplicates per period so the line is smooth
    # (defensive — validation pairs come from the aggregated matcher already).
    vd = _collapse_subunit_dups(vd, settings, extra_keys=["period"])
    # Full validation-run simulated series (spans the ENTIRE cal+val record at
    # native resolution) — drives a CONTINUOUS best-fit line like the
    # calibration chart, instead of a line through only the sparse obs-matched
    # points.  Falls back to the sparse pairs when this file is absent (so
    # already-generated reports without it don't break).
    vsim = None
    try:
        _vspath = os.path.join(str(output_dir), "results",
                               "validation_simulated_data.csv")
        if os.path.isfile(_vspath):
            vsim = pd.read_csv(_vspath)
            vsim = _collapse_subunit_dups(vsim, settings)
    except Exception:
        vsim = None
    metric = str((settings or {}).get("objective_function", "KGE")).upper()
    _mfn = {"KGE": _OF.kge, "NSE": _OF.nse, "RMSE": _OF.rmse,
            "PBIAS": _OF.pbias}.get(metric, _OF.kge)

    def _metric(sub):
        if sub is None or sub.empty:
            return None
        o = pd.to_numeric(sub["observed"], errors="coerce").to_numpy(float)
        s = pd.to_numeric(sub["simulated"], errors="coerce").to_numpy(float)
        m = ~(np.isnan(o) | np.isnan(s))
        if int(m.sum()) < 2:
            return None
        try:
            return float(_mfn(o[m], s[m]))
        except Exception:
            return None

    def _fmt(v):
        return "n/a" if v is None else f"{v:.3f}"

    # One consistent calibration->validation boundary for ALL species charts:
    # the first validation observation (else the last calibration observation).
    # A single global boundary keeps the NO3 and NH4 charts aligned even when
    # their sparse obs happen to fall in different periods.
    _global_bnd = None
    try:
        _valr = vd[vd["period"] == "validation"]
        _calr = vd[vd["period"] == "calibration"]
        if not _valr.empty:
            _global_bnd = str(
                pd.to_datetime(_valr["datetime"], errors="coerce").min())
        elif not _calr.empty:
            _global_bnd = str(
                pd.to_datetime(_calr["datetime"], errors="coerce").max())
    except Exception:
        _global_bnd = None

    # Render every species with obs pairs PLUS any that only have a full
    # simulated series (so an objective species with no obs still gets a line).
    _species = list(vd["species"].unique())
    if vsim is not None and "species" in vsim.columns:
        for _s in vsim["species"].unique():
            if _s not in _species:
                _species.append(_s)

    blocks = []
    for sp in _species:
        sd = vd[vd["species"] == sp].copy()
        sd["_dt"] = pd.to_datetime(sd["datetime"], errors="coerce")
        sd = sd.dropna(subset=["_dt"]).sort_values("_dt")
        cal = sd[sd["period"] == "calibration"]
        val = sd[sd["period"] == "validation"]
        bnd = _global_bnd
        # Observed markers (sparse, at observation times only)
        x_obs = [str(v) for v in sd["_dt"]]
        obs = [None if v != v else float(v)
               for v in pd.to_numeric(sd["observed"], errors="coerce")]
        # Best-fit LINE: prefer this species' FULL validation series (continuous,
        # native resolution, thinned like the calibration chart); fall back to
        # the sparse matched sim when the full series is unavailable.
        sim_x = x_obs
        sim_y = [None if v != v else float(v)
                 for v in pd.to_numeric(sd["simulated"], errors="coerce")]
        if vsim is not None and "species" in vsim.columns:
            fs = vsim[vsim["species"].astype(str) == str(sp)].copy()
            if not fs.empty:
                fs["_dt"] = pd.to_datetime(fs["datetime"], errors="coerce")
                fs = fs.dropna(subset=["_dt"]).sort_values("_dt")
                _xl = [str(v) for v in fs["_dt"]]
                _yl = [None if v != v else float(v)
                       for v in pd.to_numeric(fs["simulated"], errors="coerce")]
                if len(_xl) > 4000:                  # thin a long line
                    _step = (len(_xl) // 4000) + 1
                    _xl, _yl = _xl[::_step], _yl[::_step]
                sim_x, sim_y = _xl, _yl
        if not any(v is not None for v in obs) and \
           not any(v is not None for v in sim_y):
            continue
        traces = [
            {"type": "scatter", "mode": "markers", "x": x_obs, "y": obs,
             "name": "observed",
             "marker": {"size": 6, "symbol": "circle", "color": "#ef4444",
                        "line": {"width": 0.5, "color": "#7f1d1d"}},
             "hovertemplate": "%{x|%Y-%m-%d}<br>obs %{y:.4g}<extra></extra>"},
            {"type": "scatter", "mode": "lines", "x": sim_x, "y": sim_y,
             "name": "best fit", "line": {"width": 1.5, "color": "#10b981"},
             "hovertemplate": "%{x|%Y-%m-%d}<br>best %{y:.4g}<extra></extra>"},
        ]
        shapes, annos = [], []
        _x0 = sim_x[0] if sim_x else (x_obs[0] if x_obs else None)
        _x1 = sim_x[-1] if sim_x else (x_obs[-1] if x_obs else None)
        if bnd and _x0 is not None and _x1 is not None:
            shapes.append({"type": "rect", "xref": "x", "x0": _x0, "x1": bnd,
                           "yref": "paper", "y0": 0, "y1": 1,
                           "fillcolor": "rgba(16,185,129,.06)",
                           "line": {"width": 0}, "layer": "below"})
            shapes.append({"type": "rect", "xref": "x", "x0": bnd, "x1": _x1,
                           "yref": "paper", "y0": 0, "y1": 1,
                           "fillcolor": "rgba(251,146,60,.08)",
                           "line": {"width": 0}, "layer": "below"})
            shapes.append({"type": "line", "xref": "x", "x0": bnd, "x1": bnd,
                           "yref": "paper", "y0": 0, "y1": 1,
                           "line": {"color": "#888", "width": 1.2,
                                    "dash": "dash"}})
        annos.append({"xref": "paper", "x": 0.01, "yref": "paper", "y": 1.0,
                      "yanchor": "bottom", "showarrow": False, "align": "left",
                      "font": {"size": 12, "color": "#0e9f6e"},
                      "text": f"<b>Calibration</b>  {metric}={_fmt(_metric(cal))}"})
        annos.append({"xref": "paper", "x": 0.99, "yref": "paper", "y": 1.0,
                      "yanchor": "bottom", "xanchor": "right", "showarrow": False,
                      "align": "right", "font": {"size": 12, "color": "#d97706"},
                      "text": f"<b>Validation</b>  {metric}={_fmt(_metric(val))}"})
        layout = {
            "title": {"text": f"{sp} — best fit: calibration + validation"},
            "xaxis": {"title": "Date"}, "yaxis": {"title": "Concentration"},
            "hovermode": "closest", "showlegend": True,
            "legend": {"orientation": "h", "yanchor": "bottom", "y": 1.06,
                       "xanchor": "left", "x": 0},
            "margin": {"l": 64, "r": 24, "t": 100, "b": 52},
            "shapes": shapes, "annotations": annos}
        c = _plotly_chart(_plot_id("valts"), traces, layout, height=470,
                          card=False, log_axes='y')
        blocks.append(
            f'<div class="card"><h3 style="margin:.2rem 0 .6rem;">'
            f'{html_lib.escape(str(sp))}</h3>{c}</div>')
    if not blocks:
        return ""
    intro = ('<p style="color:var(--text2);margin-bottom:1rem;">Best-fit '
             'parameters re-run over the <strong>validation</strong> period. '
             'Green = calibration window, orange = validation (held-out); the '
             'objective metric is shown for each.</p>')
    return ('\n<div class="section" id="valid-ts">\n'
            f'    <h2>{_head}</h2>\n    {intro}\n    '
            + "\n    ".join(blocks) + '\n</div>\n')


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

    # Cumulative / non-cumulative distribution of the metrics across the
    # reaches/HRUs that carry observations (needs >= 2 stations; otherwise a
    # note explains why the plots are not shown).
    dist_html = _generate_perf_distribution_plots(perf_metrics,
                                                  hostmodel=hostmodel)

    # The observed-vs-simulated time-series / scatter plots now live in their
    # own "Observed vs Simulated — calibrated time series" section (placed
    # after "Run with Best Parameters"), so they aren't duplicated here.  This
    # section keeps the numeric metrics; a pointer links to the plots.
    _ts_pointer = (
        '<div class="highlight-box info" style="margin-top:1rem">'
        'See the <a href="#timeseries" style="color:var(--primary);'
        'font-weight:600;">Observed vs Simulated time series</a> section '
        '(below &ldquo;Run with Best Parameters&rdquo;) for the per-species '
        'time-series and observed-vs-simulated plots.</div>'
        if matched_data is not None else "")

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
    {dist_html}
    {_ts_pointer}
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
    # Resolve basin shapefile (lives in ss_method_copernicus_basins_hrus)
    basin_info = model_config.get("ss_method_copernicus_basins_hrus") or {}
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


def _build_next_steps_section(
    model_config: Dict[str, Any],
    settings: Dict[str, Any],
) -> str:
    """Flat, numbered 'what to do now' guide (mirrors the model-configuration
    report's step-by-step 'Next Steps' presentation)."""
    # Validation window, if it's in the settings, to make step 3 concrete.
    vp = (settings.get("validation_period")
          or settings.get("validation_window") or None)
    vtxt = ""
    if isinstance(vp, (list, tuple)) and len(vp) == 2 and vp[0] and vp[1]:
        vtxt = (f" (your validation window: "
                f"<code>{html_lib.escape(str(vp[0]))}</code> &rarr; "
                f"<code>{html_lib.escape(str(vp[1]))}</code>)")

    def _step(n, title, body):
        return (
            '<h3 style="margin-top:1.1rem;display:flex;align-items:center;'
            'gap:.5rem;">'
            '<span style="flex-shrink:0;display:inline-flex;width:1.5rem;'
            'height:1.5rem;border-radius:50%;background:var(--primary);'
            'color:#fff;align-items:center;justify-content:center;'
            f'font-size:.8rem;font-weight:700;">{n}</span>{title}</h3>'
            f'<div style="margin:.1rem 0 .2rem;line-height:1.6;">{body}</div>')

    return f"""
<div class="section" id="next-steps">
    <h2>Next Steps</h2>
    <div class="card primary">
      <p style="margin-top:0;">Now that the calibration has finished, here's what
      to do with the results:</p>
      {_step(1, "Review the fit",
             "Check the headline metric (KGE / NSE / RMSE / PBIAS) and the "
             "<em>Convergence</em>, <em>parameter</em> and "
             "<em>observed&nbsp;vs&nbsp;simulated time&nbsp;series</em> plots above. "
             "The observation map colours each station by its performance.")}
      {_step(2, "Apply the best parameters",
             "Write the calibrated values into your model's BGC configuration "
             "(<code>python &lt;your_calibration_template&gt;.py --apply-best</code>), "
             "then re-run your model with them. The complete ready-to-run setup "
             "is in <a href='#run-best'><em>Run with Best Parameters</em></a> "
             "below.")}
      {_step(3, "Validate",
             "Run the model over the <strong>validation period</strong>" + vtxt +
             " &mdash; the window held out after your calibration split &mdash; to "
             "confirm the fit generalises to data the calibration never saw.")}
      {_step(4, "Refine if needed",
             "If the fit is weak, reopen the calibration setup report, widen or "
             "adjust the parameter bounds (or pick different parameters / a "
             "different calibration window), and re-calibrate.")}
    </div>
</div>
"""


def _minify_python(src: str) -> str:
    """Strip comments (full-line and inline) and collapse blank lines from a
    Python source, keeping it byte-for-byte RUNNABLE. Uses ``tokenize`` so it
    never touches ``#`` inside strings, and re-parses the result as a guard —
    on any doubt it returns the original source unchanged."""
    try:
        import tokenize as _tk, io as _io
        lines = src.splitlines(keepends=True)
        drop_line = set()          # 1-based full-line comment lines
        inline_at = {}             # 1-based line -> col where an inline comment starts
        for t in _tk.generate_tokens(_io.StringIO(src).readline):
            if t.type == _tk.COMMENT:
                r, c = t.start
                if lines[r - 1][:c].strip() == "":
                    drop_line.add(r)
                else:
                    inline_at[r] = min(inline_at.get(r, 1 << 30), c)
        out, prev_blank = [], False
        for i, ln in enumerate(lines, 1):
            if i in drop_line:
                continue
            if i in inline_at:
                ln = ln[:inline_at[i]].rstrip() + "\n"
            if ln.strip() == "":
                if prev_blank:
                    continue
                prev_blank = True
            else:
                prev_blank = False
            out.append(ln)
        result = "".join(out).strip() + "\n"
        compile(result, "<minified>", "exec")   # must still parse -> still runnable
        return result
    except Exception:
        return src


def _build_code_toolbox(code_text: str, prefix: str, lang: str = "python",
                        save_name: str = "", save_path: str = "",
                        meta: str = "") -> str:
    """A reusable, ID-scoped code box: the code in a collapsible/lockable scroll
    box (Expand/Collapse + Lock-scroll) plus optional Save-As (.py) and Copy-path
    controls. Every element id is namespaced by `prefix`, so multiple boxes can
    coexist on one page. Self-contained (its own IIFE reads the code from the DOM)."""
    n_lines = code_text.count("\n") + 1
    save_btn = (f'<button type="button" id="{prefix}-save" class="rb-save-btn" '
                f'title="Open a Save-As dialog to choose where to save this .py file">'
                f'&#128190; Save (.py)&hellip;</button>') if save_name else ""
    path_row = (f'<span class="rb-path"><span>Path:</span>'
                f'<span class="p" id="{prefix}-path">{html_lib.escape(save_path)}</span>'
                f'<button type="button" id="{prefix}-copypath" class="rb-copy-path" '
                f'title="Copy this path">&#128203; Copy path</button></span>') if save_path else ""
    toolbar = (f'<div class="rb-toolbar">{save_btn}{path_row}</div>'
               if (save_btn or path_row) else "")
    controls = (
        '<div class="rb-code-controls">'
        f'<button type="button" id="{prefix}-expand" class="rb-expand">&#9660; Expand</button>'
        f'<button type="button" id="{prefix}-lock" class="rb-expand" '
        f'title="Lock the box so the page scrolls past it">&#128274; Lock scroll</button>'
        f'<span class="rb-code-meta">{n_lines} lines'
        f'{(" &mdash; " + meta) if meta else ""} &mdash; scroll inside, lock, or expand</span>'
        '</div>')
    body = (f'{toolbar}{controls}'
            f'<div id="{prefix}-wrap" class="rb-codewrap">'
            f'<div id="{prefix}-code">{rh.build_code_block(code_text, lang)}</div></div>')
    fname_j = json.dumps(save_name or "script.py")
    fpath_j = json.dumps(save_path or "")
    script = f"""
    <script>(function() {{
      var P = "{prefix}";
      var wrap = document.getElementById(P + "-wrap");
      var expBtn = document.getElementById(P + "-expand");
      var lockBtn = document.getElementById(P + "-lock");
      var saveBtn = document.getElementById(P + "-save");
      var copyBtn = document.getElementById(P + "-copypath");
      function codeText() {{ var c = document.querySelector("#" + P + "-code code");
                            return c ? c.textContent : ""; }}
      function flash(b) {{ var o = b.innerHTML; b.innerHTML = "&#10003; Copied";
                          setTimeout(function() {{ b.innerHTML = o; }}, 1500); }}
      function fbCopy(t, b) {{ var ta = document.createElement("textarea"); ta.value = t;
        ta.style.position = "fixed"; ta.style.opacity = "0"; document.body.appendChild(ta);
        ta.focus(); ta.select(); try {{ document.execCommand("copy"); }} catch (e) {{}}
        ta.remove(); flash(b); }}
      if (expBtn && wrap) expBtn.addEventListener("click", function() {{
        var e = wrap.classList.toggle("expanded");
        expBtn.innerHTML = e ? "&#9650; Collapse" : "&#9660; Expand";
        if (!e) wrap.scrollTop = 0; }});
      if (lockBtn && wrap) lockBtn.addEventListener("click", function() {{
        var l = wrap.classList.toggle("locked");
        lockBtn.classList.toggle("locked-on", l);
        lockBtn.innerHTML = l ? "&#128275; Unlock scroll" : "&#128274; Lock scroll"; }});
      if (copyBtn) copyBtn.addEventListener("click", function() {{
        var fp = {fpath_j};
        if (navigator.clipboard && navigator.clipboard.writeText)
          navigator.clipboard.writeText(fp).then(function() {{ flash(copyBtn); }},
                                                  function() {{ fbCopy(fp, copyBtn); }});
        else fbCopy(fp, copyBtn); }});
      function blobDL(t) {{ var b = new Blob([t], {{type: "text/x-python"}});
        var a = document.createElement("a"); a.href = URL.createObjectURL(b);
        a.download = {fname_j}; document.body.appendChild(a); a.click();
        setTimeout(function() {{ URL.revokeObjectURL(a.href); a.remove(); }}, 0); }}
      if (saveBtn) saveBtn.addEventListener("click", async function() {{
        var t = codeText();
        if (window.showSaveFilePicker) {{
          try {{
            var h = await window.showSaveFilePicker({{ suggestedName: {fname_j},
              types: [{{description: "Python script", accept: {{"text/x-python": [".py"]}}}}] }});
            var w = await h.createWritable(); await w.write(t); await w.close(); return;
          }} catch (e) {{ if (e && e.name === "AbortError") return; }}
        }}
        blobDL(t); }});
    }})();</script>"""
    return body + script


def _build_run_best_section(
    results: Dict[str, Any],
    model_config: Dict[str, Any],
    settings: Dict[str, Any],
    calibration_parameters: Optional[List[Dict[str, Any]]] = None,
    output_dir: str = "",
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
    template_abspath = os.path.abspath(template_path)
    template_dir = os.path.dirname(template_abspath)

    # ── PRIMARY: run the calibrated model = re-run the GLOBAL-BEST evaluation's
    #    config, which ALREADY has the best parameters applied (that is how its
    #    result was produced). Re-running the model-config TEMPLATE instead only
    #    regenerates the template's DEFAULT parameter values — so it must not be
    #    mistaken for the calibrated model. ──
    import glob as _glob2, re as _re2
    _best_n = None
    try:
        _best_n = _best_eval_number(output_dir) if output_dir else None
    except Exception:
        _best_n = None
    _best_dir = None
    if _best_n is not None and output_dir:
        for _d in _glob2.glob(os.path.join(str(output_dir), "evaluations", "eval_*")):
            _m = _re2.search(r"eval_(\d+)", os.path.basename(_d))
            if _m and int(_m.group(1)) == _best_n:
                _best_dir = os.path.abspath(_d)
                break
    _exe = model_config.get("executable_path", "") or ""
    _host = (model_config.get("hostmodel") or "").lower()
    _bp_path = os.path.join(os.path.abspath(str(output_dir)), "results",
                            "best_parameters.json") if output_dir else ""

    # Container config: the executable is a LINUX binary built in the openWQ
    # container — it cannot run natively (e.g. on macOS). So the reproduce
    # command must `docker exec` into the SAME container the calibration used,
    # with host paths translated to the container mount (…/Documents -> /code).
    _container, _compose, _runtime, _cc = "docker_openwq", "", "docker", {}
    try:
        from calibration_lib import config_integration as _ci2
        _cc = _ci2.get_container_config(model_config) or {}
        _container = _cc.get("docker_container_name") or _container
        _compose = _cc.get("docker_compose_path") or model_config.get("docker_compose_path") or ""
        _runtime = (_cc.get("container_runtime") or "docker").lower()
    except Exception:
        _compose = model_config.get("docker_compose_path") or ""
    _hroot, _croot = None, "/code"
    if _compose and os.path.isfile(_compose):
        try:
            _vm = _re2.search(r'volumes:\s*\n\s*-\s*([^:\n]+):([^:\s]+)', open(_compose).read())
            if _vm:
                _croot = _vm.group(2).strip()
                _hroot = os.path.abspath(os.path.join(os.path.dirname(_compose),
                                                       _vm.group(1).strip()))
        except Exception:
            pass
    def _c(p):
        return (_croot + p[len(_hroot):]) if (_hroot and p and p.startswith(_hroot)) else p

    calibrated_block = ""
    if _best_dir and _exe and os.path.isdir(_best_dir):
        _master = os.path.join(_best_dir, "openWQ_master.json")
        if _host == "summa":
            _ctrl, _flag, _np = os.path.join(_best_dir, "fileManager_eval.txt"), "-m ", "1"
        else:
            _cands = (_glob2.glob(os.path.join(_best_dir, "*.control"))
                      + _glob2.glob(os.path.join(_best_dir, "*.control.txt")))
            _ctrl, _flag, _np = (_cands[0] if _cands else os.path.join(_best_dir, "<control_file>")), "", "2"
        _out_dir = os.path.join(_best_dir, "openwq_out")
        _sif, _bind = _resolve_hpc_settings(_cc, settings or {}, model_config)

        # Debug export = re-run with RUN_MODE_DEBUG=true (openWQ then also writes
        # the per-process terms: dt_chem/sorpt/transp, ss, ewf, ic). Done
        # NON-destructively: a `sed` on the host writes an openWQ_master_debug.json
        # copy with the flag flipped, and master_json points at that copy (the
        # eval's original master is untouched; all other paths in it are relative).
        # NOTE: every host link (SUMMA + mizuRoute OpenWQ_hydrolink.cpp) HARDCODES
        # `set_OpenWQ_masterjson("openWQ_master.json")` -> openWQ always reads
        # ./openWQ_master.json from the run cwd and IGNORES the master_json env
        # var. So a separate debug copy is never read; debug export must flip
        # RUN_MODE_DEBUG in the cwd master ITSELF: back it up, write the flipped
        # version in place, run, then restore the original (net non-destructive).
        _orig_bak = os.path.join(_best_dir, ".openWQ_master.orig.json")
        _swap_in = (
            f'cd "{_best_dir}" && cp openWQ_master.json "{_orig_bak}" && '
            "sed -E 's/(\"RUN_MODE_DEBUG\"[[:space:]]*:[[:space:]]*)false/\\1true/' "
            f'"{_orig_bak}" > openWQ_master.json')
        _swap_back = f'mv "{_orig_bak}" "{_best_dir}/openWQ_master.json"'

        def _mk_cmd(pre="", post=""):
            if _runtime == "docker":
                body = (
                    f'docker exec -e OMP_NUM_THREADS=1 -e master_json={_c(_master)} {_container} \\\n'
                    f'  /bin/bash -lc "cd {_c(_best_dir)} && \\\n'
                    f'    mpirun --allow-run-as-root -np {_np} -x master_json -x OMP_NUM_THREADS \\\n'
                    f'    {_c(_exe)} {_flag}{_c(_ctrl)}"')
            else:
                # cd into the eval dir so ./openWQ_master.json is the one read.
                body = (
                    f'cd "{_best_dir}" && OMP_NUM_THREADS=1 apptainer exec \\\n'
                    f'  --bind "{_bind}" \\\n'
                    f'  --env master_json="{_master}" \\\n'
                    f'  "{_sif}" \\\n'
                    f'  mpirun -np {_np} "{_exe}" {_flag}"{_ctrl}"')
            return (pre + "\n" if pre else "") + body + ("\n" + post if post else "")

        _cal_cmd = _mk_cmd()                                     # normal run
        _cal_cmd_dbg = _mk_cmd(pre=_swap_in, post=_swap_back)    # debug-export run
        if _runtime == "docker":
            _how = ('run it <strong>inside the openWQ Docker container</strong> '
                    f'(<code>{html_lib.escape(_container)}</code> &mdash; started by the '
                    'calibration; if stopped, <code>docker compose up -d</code> from the '
                    '<code>containers/</code> folder):')
        else:
            _how = ('run it with <strong>Apptainer/Singularity</strong> on your HPC '
                    '(bind the folder holding the config &amp; outputs):')

        # Debug-export tickbox: swaps the command between the two variants live.
        _cal_ctrl = (
            '<div style="display:flex;flex-wrap:wrap;align-items:center;gap:.5rem .7rem;'
            'margin:.1rem 0 .5rem;">'
            '<label class="viz-cb"><input type="checkbox" id="rbcal-debug"> Debug export</label>'
            '<span style="font-size:.8rem;color:var(--text2);line-height:1.5;">'
            '<code>RUN_MODE_DEBUG = true</code> &mdash; also writes the per-process terms '
            '(larger output) and enables the <em>Debug</em> option in the output report below.'
            '</span></div>'
            '<div id="rbcal-codewrap">' + rh.build_code_block(_cal_cmd, "bash") + '</div>'
            '<script>(function(){'
            f'var N={json.dumps(_cal_cmd)},D={json.dumps(_cal_cmd_dbg)};'
            'var cb=document.getElementById("rbcal-debug");'
            'var code=document.querySelector("#rbcal-codewrap code");'
            'if(cb&&code)cb.addEventListener("change",function(){code.textContent=cb.checked?D:N;});'
            '})();</script>')

        # ── Output-report controls: INTERACTIVE tickboxes (species/compartments)
        #    + space-elements + debug that rebuild a copy-paste heredoc snippet
        #    LIVE (like the config report). Options come from THIS run's actual
        #    output files. Mirrors 2_Read_Outputs/reading_plotting_results_template.py. ──
        _h5lib = os.path.abspath(os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            "..", "..", "2_Read_Outputs", "hdf5_support_lib"))
        _mapkey = (model_config.get("h5_mapping_key")
                   or ("hruId" if _host == "summa" else "reachID"))
        _feat = (model_config.get("river_network_mapping_key")
                 or ("hruId" if _host == "summa" else "SegId"))
        _river_shp = model_config.get("river_network_shapefile") or ""
        _basin = model_config.get("ss_method_copernicus_basins_hrus") or {}
        _basin_shp = (_basin.get("path_to_shp") if isinstance(_basin, dict) else "") or ""
        _rep_out = os.path.join(_out_dir, "openwq_calibrated_report.html")
        # Observations for the plot: exactly the ones the calibration was run
        # against (`<calibration_dir>/calibration_observations.csv`, written by
        # the driver with the feature id each record was compared to), plus the
        # compartment(s) they belong to (from the model config).
        _cal_obs_csv = os.path.join(str(output_dir), "calibration_observations.csv")
        _obs_comps = list(model_config.get("observation_compartments") or [])
        _basin_key = (_basin.get("mapping_key") if isinstance(_basin, dict) else "") or ""
        _obs_lines = ""
        if os.path.isfile(_cal_obs_csv):
            _obs_lines = (f"    observation_calibration_csv=r'{os.path.abspath(_cal_obs_csv)}',\n"
                          f"    observation_compartments={json.dumps(_obs_comps)},\n")
        if _basin_key:
            _obs_lines += f"    basin_mapping_key='{_basin_key}',\n"
        # discover available compartments + species from the output filenames
        _comps, _specs = set(), set()
        for _f in _glob2.glob(os.path.join(_out_dir, "HDF5", "*.h5")):
            _b = os.path.basename(_f)
            if "@" in _b:
                _cn, _rest = _b.split("@", 1)
                _comps.add(_cn)
                if "#" in _rest:
                    _specs.add(_rest.split("#")[0])
        if not _specs:                               # CSV output (no species in names)
            for _f in _glob2.glob(os.path.join(_out_dir, "CSV", "*.txt")):
                _b = os.path.basename(_f)
                if "@" in _b:
                    _comps.add(_b.split("@", 1)[0])
            _specs = set(model_config.get("chemical_species") or [])
        if not _specs:
            _specs = set(model_config.get("chemical_species") or ["NO3-N"])
        _def_cmp_name = "ILAYERVOLFRACWAT_SOIL" if _host == "summa" else "RIVER_NETWORK_REACHES"
        if not _comps:
            _comps = {_def_cmp_name}
        _specs, _comps = sorted(_specs), sorted(_comps)
        _cmp_default = _def_cmp_name if _def_cmp_name in _comps else (_comps[0] if _comps else "")
        # output format actually produced by this run (the driver appends /<FMT>/)
        _out_fmt = ("HDF5" if _glob2.glob(os.path.join(_out_dir, "HDF5", "*.h5"))
                    else "CSV" if _glob2.glob(os.path.join(_out_dir, "CSV", "*"))
                    else "HDF5")

        def _cbs(items, cls, checked):
            return "".join(
                f'<label class="viz-cb"><input type="checkbox" class="{cls}" '
                f'value="{html_lib.escape(it)}"' + (" checked" if it in checked else "")
                + f'> {html_lib.escape(it)}</label>' for it in items)
        _species_cbs = _cbs(_specs, "rbviz-sp", set(_specs))
        _comp_cbs = _cbs(_comps, "rbviz-cmp", {_cmp_default})

        _viz_template = (
            "python3 - <<'PY'\n"
            "# ---- Output report for the calibrated best simulation (runs on host) ----\n"
            "SPECIES      = __SPECIES__\n"
            "COMPARTMENTS = __COMPARTMENTS__\n"
            "SPACE_ELEM   = __SPACE_ELEM__\n"
            "DEBUG        = __DEBUG__\n"
            f"H5LIB        = r'{_h5lib}'   # openWQ supporting_scripts/2_Read_Outputs/hdf5_support_lib (edit if relocated)\n"
            "\n"
            "import sys, os\n"
            "if not os.path.isdir(H5LIB):\n"
            "    raise SystemExit('H5LIB not found - set it to your openWQ supporting_scripts/2_Read_Outputs/hdf5_support_lib')\n"
            "sys.path.insert(0, H5LIB)\n"
            "import Read_h5_driver as h5_rlib\n"
            "import Plot_h5_driver as h5_plib\n"
            "\n"
            f"openwq_info = {{'path_to_results': r'{_out_dir}', 'mapping_key': '{_mapkey}'}}\n"
            "openwq_results = h5_rlib.Read_h5_driver(\n"
            f"    openwq_info=openwq_info, output_format='{_out_fmt}',\n"
            "    debugmode=DEBUG, cmp=COMPARTMENTS, space_elem=SPACE_ELEM,\n"
            "    chemSpec=SPECIES, chemUnits='MG/L')\n"
            "\n"
            "h5_plib.Plot_h5_driver(\n"
            f"    what2map='openwq', hostmodel='{_host}', mapping_key_values='all',\n"
            "    openwq_results=openwq_results, chemSpec=SPECIES, debugmode=DEBUG,\n"
            f"    output_path=r'{_rep_out}',\n"
            f"    mapping_key='{_feat}', feature_label='{_feat}',\n"
            f"    river_network_shp=r'{_river_shp}', basin_shapefile=r'{_basin_shp}',\n"
            + _obs_lines +
            "    )\n"
            f"print('Report:', r'{_rep_out}')\n"
            "import webbrowser\n"
            f"webbrowser.open('file://' + os.path.abspath(r'{_rep_out}'))   # open it automatically\n"
            "PY\n"
        )
        _tpl_json = json.dumps(_viz_template)
        _viz_js = (
            "<script>(function(){\n"
            "var TPL=" + _tpl_json + ";\n"
            'var code=document.querySelector("#rbviz-codewrap code");\n'
            'var sp=document.getElementById("rbviz-space");\n'
            'var dbg=document.getElementById("rbviz-debug");\n'
            'function L(sel){var a=[];document.querySelectorAll(sel).forEach('
            'function(n){if(n.checked)a.push(JSON.stringify(n.value));});return "["+a.join(", ")+"]";}\n'
            "function build(){if(!code)return;\n"
            'var sv=(sp&&sp.value.trim())?sp.value.trim():"\'all\'";\n'
            'var d=(dbg&&dbg.checked)?"True":"False";\n'
            'code.textContent=TPL.replace("__SPECIES__",L(".rbviz-sp"))'
            '.replace("__COMPARTMENTS__",L(".rbviz-cmp"))'
            '.replace("__SPACE_ELEM__",sv).replace("__DEBUG__",d);}\n'
            'document.querySelectorAll(".rbviz-sp,.rbviz-cmp").forEach('
            'function(n){n.addEventListener("change",build);});\n'
            "if(sp)sp.addEventListener(\"input\",build);"
            "if(dbg)dbg.addEventListener(\"change\",build);\n"
            "build();\n"
            "})();</script>"
        )
        viz_block = (
            '<p style="margin:1rem 0 .2rem;font-weight:600;">&#9654; Generate the output '
            'report for this simulation</p>'
            '<p style="margin:.1rem 0 .5rem;font-size:.85rem;color:var(--text2);">'
            'Tick the species / compartments (from this run&rsquo;s output), set the options, '
            'then copy &amp; paste the <strong>whole block</strong> into your terminal:</p>'
            '<div class="viz-ctrl">'
            f'<div class="viz-group"><span class="viz-label">Species</span>{_species_cbs}</div>'
            f'<div class="viz-group"><span class="viz-label">Compartments</span>{_comp_cbs}</div>'
            '<div class="viz-group"><span class="viz-label">Space elements</span>'
            '<input type="text" id="rbviz-space" class="viz-in" value="\'all\'"'
            ' title="\'all\' or a Python list of ids, e.g. [1, 2]"></div>'
            '<div class="viz-group"><label class="viz-cb"><input type="checkbox" id="rbviz-debug"> '
            'Debug (per-process terms)</label></div>'
            '</div>'
            '<div id="rbviz-codewrap">' + rh.build_code_block("", "bash") + '</div>'
            + _viz_js)

        calibrated_block = (
            '<div class="card primary" style="margin-bottom:1rem;">'
            '<h3 style="margin-top:0;">&#9654; Run the calibrated model</h3>'
            '<p style="font-size:.9rem;margin:.2rem 0 .5rem;">'
            f'The best evaluation (<code>{html_lib.escape(os.path.basename(_best_dir))}</code>) '
            'already ran with the <strong>calibrated best parameters</strong>, so '
            '<strong>its output already IS the calibrated best simulation</strong> '
            f'(no re-run needed to see results). To reproduce it, {_how}</p>'
            + _cal_ctrl +
            '<p style="font-size:.8rem;color:var(--text2);margin:.5rem 0 0;line-height:1.6;">'
            f'&#128193; Calibrated output (already on disk): <code>{html_lib.escape(_out_dir)}</code><br>'
            f'&#9881;&#65039; Calibrated config folder: <code>{html_lib.escape(_best_dir)}</code>'
            + (f'<br>&#128202; Best parameter values: <code>{html_lib.escape(_bp_path)}</code>'
               if _bp_path else '')
            + '</p>'
            + viz_block
            + '</div>')
    else:
        # Best-eval config not locatable — fall back to the raw best values.
        _bp = results.get("best_params") or {}
        calibrated_block = (
            '<div class="card primary" style="margin-bottom:1rem;">'
            '<h3 style="margin-top:0;">&#9654; Run the calibrated model</h3>'
            '<p style="font-size:.9rem;">The calibrated best simulation is the '
            '<strong>best evaluation</strong> under <code>evaluations/</code> '
            '(its config already has the best parameters applied). The best values '
            + (f'are saved in <code>{html_lib.escape(_bp_path)}</code> and '
               if _bp_path else '')
            + 'are listed below &mdash; apply them to your model config, then run.</p>'
            + (rh.build_code_block(json.dumps(_bp, indent=2, default=str), "json") if _bp else "")
            + '</div>')

    # Minified, still-runnable version (comments + blank lines stripped) so the
    # config shown here is compact; the full template stays on disk untouched.
    _orig_lines = template_text.count("\n") + 1
    template_text = _minify_python(template_text)
    _min_lines = template_text.count("\n") + 1
    _trim = (f" (minimized {_orig_lines}&rarr;{_min_lines} lines)"
             if _min_lines < _orig_lines else "")

    sub = (
        '<div class="rb-warn"><strong>&#9888; This is your model setup template with its '
        'DEFAULT parameter values &mdash; not the calibrated ones.</strong> Re-running it '
        'regenerates the openWQ config from these defaults (handy to inspect or edit your '
        'setup). To run the <em>calibrated</em> model, use '
        '<em>&#9654; Run the calibrated model</em> above.</div>'
        f'<p style="margin-top:.4rem">Model config template{_trim}. Source: '
        f'<code>{html_lib.escape(template_filename)}</code>.</p>'
    )

    # "Code to run it": execute the saved template. Its paths are absolute, so
    # this works from any working directory once the file is saved.
    run_cmd = f'python "{template_abspath}"'
    run_block = (
        '<p style="margin:.9rem 0 .2rem;font-weight:600;">&#9654; Run the template '
        '(regenerates DEFAULT config)</p>'
        '<p style="margin:.1rem 0 .3rem;font-size:.85rem;color:var(--text2);">'
        'Regenerates the openWQ input config (master + <code>openwq_in/</code>) from these '
        '<strong>default</strong> settings &mdash; not the calibrated parameters:</p>'
        + rh.build_code_block(run_cmd, "bash")
    )

    # Save-template button + the exact path to save it at (so its baked paths
    # resolve when run). The button downloads the script text verbatim; the
    # path is the template's canonical location on disk.
    fname_json = json.dumps(template_filename)
    path_json = json.dumps(template_abspath)
    toolbar = f"""
        <div class="rb-toolbar">
            <button type="button" id="rb-save-btn" class="rb-save-btn"
                    title="Open a Save-As dialog to choose where to save this .py file">
                &#128190; Save config template (.py)&hellip;
            </button>
            <span class="rb-path">
                <span>Save it here so its paths resolve:</span>
                <span class="p" id="rb-save-path">{html_lib.escape(template_abspath)}</span>
                <button type="button" id="rb-copy-path-btn" class="rb-copy-path"
                        title="Copy this path">&#128203; Copy path</button>
            </span>
        </div>
        <p class="rb-hint">Copy the path, click <em>Save config template</em>, then paste it
           in the Save dialog (on macOS press <kbd>&#8679;&#8984;G</kbd> to go to the folder)
           so the file lands where its baked paths resolve.</p>"""
    script = f"""
        <script>
        (function() {{
            var FNAME = {fname_json}, FPATH = {path_json};
            var saveBtn = document.getElementById('rb-save-btn');
            var copyBtn = document.getElementById('rb-copy-path-btn');
            function blobDownload(txt) {{
                var blob = new Blob([txt], {{type: 'text/x-python'}});
                var a = document.createElement('a');
                a.href = URL.createObjectURL(blob); a.download = FNAME;
                document.body.appendChild(a); a.click();
                setTimeout(function() {{ URL.revokeObjectURL(a.href); a.remove(); }}, 0);
            }}
            if (saveBtn) saveBtn.addEventListener('click', async function() {{
                var code = document.querySelector('#run-best-code code');
                var txt = code ? code.textContent : '';
                // Preferred: native Save-As dialog (File System Access API) so the
                // user can choose the location / paste the path. Falls back to a
                // plain download on browsers without it (Firefox/Safari) or file://.
                if (window.showSaveFilePicker) {{
                    try {{
                        var handle = await window.showSaveFilePicker({{
                            suggestedName: FNAME,
                            types: [{{description: 'Python script',
                                      accept: {{'text/x-python': ['.py']}}}}]
                        }});
                        var w = await handle.createWritable();
                        await w.write(txt); await w.close();
                        return;
                    }} catch (e) {{
                        if (e && e.name === 'AbortError') return;   // user cancelled
                        // any other error (e.g. blocked on file://) -> fall back
                    }}
                }}
                blobDownload(txt);
            }});
            function flash(btn) {{ var o = btn.innerHTML; btn.innerHTML = '&#10003; Copied';
                setTimeout(function() {{ btn.innerHTML = o; }}, 1500); }}
            function fallbackCopy(t, btn) {{ var ta = document.createElement('textarea');
                ta.value = t; ta.style.position = 'fixed'; ta.style.opacity = '0';
                document.body.appendChild(ta); ta.focus(); ta.select();
                try {{ document.execCommand('copy'); }} catch (e) {{}} ta.remove(); flash(btn); }}
            if (copyBtn) copyBtn.addEventListener('click', function() {{
                if (navigator.clipboard && navigator.clipboard.writeText) {{
                    navigator.clipboard.writeText(FPATH).then(
                        function() {{ flash(copyBtn); }},
                        function() {{ fallbackCopy(FPATH, copyBtn); }});
                }} else {{ fallbackCopy(FPATH, copyBtn); }}
            }});
            var wrap = document.getElementById('rb-codewrap');
            var expBtn = document.getElementById('rb-expand-btn');
            if (expBtn && wrap) expBtn.addEventListener('click', function() {{
                var expanded = wrap.classList.toggle('expanded');
                expBtn.innerHTML = expanded ? '&#9650; Collapse' : '&#9660; Expand';
                if (!expanded) wrap.scrollTop = 0;
            }});
            var lockBtn = document.getElementById('rb-lock-btn');
            if (lockBtn && wrap) lockBtn.addEventListener('click', function() {{
                var locked = wrap.classList.toggle('locked');
                lockBtn.classList.toggle('locked-on', locked);
                lockBtn.innerHTML = locked ? '&#128275; Unlock scroll' : '&#128274; Lock scroll';
            }});
        }})();
        </script>"""

    # ── Calibration config + run script (_run.py) — the script that configured &
    #    ran THIS calibration. Shown BEFORE the run-the-calibrated-model snippet.
    #    Uses the ID-scoped code toolbox (prefix "rbrun") so it coexists with the
    #    template box (prefix "rb-…"). ──
    runpy_card = ""
    _runpy_cands = (_glob2.glob(os.path.join(str(output_dir), "*_run.py"))
                    if output_dir else [])
    if _runpy_cands:
        _runpy_path = os.path.abspath(_runpy_cands[0])
        _runpy_name = os.path.basename(_runpy_path)
        try:
            with open(_runpy_path, "r", encoding="utf-8", errors="replace") as _rf:
                _runpy_text = _rf.read()
        except Exception:
            _runpy_text = ""
        _ro = _runpy_text.count("\n") + 1
        _runpy_text = _minify_python(_runpy_text)
        _rm = _runpy_text.count("\n") + 1
        _rmeta = f"minimized {_ro}&rarr;{_rm}" if _rm < _ro else ""
        _runpy_box = _build_code_toolbox(
            _runpy_text, "rbrun", save_name=_runpy_name,
            save_path=_runpy_path, meta=_rmeta)
        runpy_card = (
            '<div class="card">'
            '<h3 style="margin-top:0;">Best parameter simulation run script</h3>'
            '<p style="margin:.2rem 0 .3rem;font-size:.9rem;">The script that configured and '
            f'ran this calibration: <code>{html_lib.escape(_runpy_name)}</code>.</p>'
            + _runpy_box
            + '</div>')

    return f"""
<div class="section" id="run-best">
    <h2>Run with Best Parameters</h2>
    {runpy_card}
    {calibrated_block}
</div>
"""




# =========================================================================
# Plot Generators (require matplotlib)
# =========================================================================

def _fig_to_base64(fig) -> str:
    """Convert a matplotlib figure to base64-encoded PNG data URI."""
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=150, bbox_inches='tight',
                facecolor='#ffffff', edgecolor='none')
    buf.seek(0)
    data = rh.embed_image_base64.__wrapped__(buf) if hasattr(rh.embed_image_base64, '__wrapped__') else None

    # Direct encoding
    buf.seek(0)
    import base64
    encoded = base64.b64encode(buf.read()).decode('ascii')
    plt.close(fig)
    return f"data:image/png;base64,{encoded}"


def _setup_plot_style():
    """Configure matplotlib for clean, light, print-friendly figures.

    Matches the light/vector aesthetic of the 1_Model_Config and
    2_Read_Outputs workflows (white background, dark text) so the figures
    look professional both in the report and when exported to PDF.
    """
    plt.rcParams.update({
        'figure.facecolor': '#ffffff',
        'axes.facecolor': '#ffffff',
        'axes.edgecolor': '#cbd5e1',
        'axes.labelcolor': '#374151',
        'text.color': '#1f2937',
        'xtick.color': '#4b5563',
        'ytick.color': '#4b5563',
        'grid.color': '#e5e7eb',
        'grid.alpha': 0.9,
        'font.size': 10,
        'font.family': 'sans-serif',
    })


# =========================================================================
# Interactive Plotly charts (built in JS via Plotly.newPlot, loaded from the
# same CDN the 1_Model_Config / 2_Read_Outputs reports use).  Charts are
# theme-aware (light/dark follow the report toggle) and carry the mode-bar
# SVG download — the native, vector "export" those workflows use.
# =========================================================================

_PLOT_SEQ = [0]


def _plot_id(prefix: str) -> str:
    """Unique DOM id for a chart container."""
    _PLOT_SEQ[0] += 1
    return f"owq-{prefix}-{_PLOT_SEQ[0]}"


def _plotly_bootstrap(project_name: str = "") -> str:
    """Plotly.js (CDN) + shared theme-aware layout/config helpers, sortable
    tables, and the per-figure "Lock zoom" / branded "PDF" export (replicated
    from the FLUXOS results report).  Emitted once, near the top of <body>."""
    _head = (
        "<style>\n"
        ".plot-controls{display:flex;justify-content:flex-end;align-items:center;"
        "gap:.4rem;margin:.5rem 0 -.1rem;font-size:.82rem;color:var(--text2);flex-wrap:wrap;}\n"
        ".plot-controls .hint{font-size:.72rem;color:var(--text3);margin-right:auto;font-style:italic;}\n"
        ".plot-controls label{display:inline-flex;align-items:center;gap:.4rem;cursor:pointer;"
        "user-select:none;padding:.2rem .55rem;border-radius:6px;border:1px solid var(--border);"
        "background:var(--surface);}\n"
        ".plot-controls label:hover{border-color:var(--primary);}\n"
        ".plot-controls input[type=checkbox]{cursor:pointer;accent-color:var(--primary);margin:0;}\n"
        ".plot-controls .pdf-btn{display:inline-flex;align-items:center;gap:.35rem;cursor:pointer;"
        "font:inherit;font-size:.82rem;color:var(--text2);padding:.2rem .6rem;border-radius:6px;"
        "border:1px solid var(--border);background:var(--surface);}\n"
        ".plot-controls .pdf-btn:hover{border-color:var(--primary);color:var(--primary);}\n"
        ".plot-controls .pdf-btn:disabled{opacity:.55;cursor:wait;}\n"
        "</style>\n"
        "<script>\n"
        "window.__OWQ_PROJECT_NAME = "
        + json.dumps(project_name or "OpenWQ Calibration") + ";\n"
        "window.__OWQ_GENERATED_AT = "
        + json.dumps(datetime.now().isoformat()) + ";\n"
        "</script>\n"
    )
    return _head + """
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
<script>
window._owqPlots = window._owqPlots || [];
function _owqDark(){
  var t=(document.documentElement.getAttribute('data-theme')||
         document.body.getAttribute('data-theme')||'');
  return t==='dark';
}
function _owqLayout(extra){
  var t=_owqDark();
  var base={ autosize:true, margin:{l:64,r:24,t:48,b:52},
    plot_bgcolor:t?'#1e1f2e':'#ffffff', paper_bgcolor:t?'#1e1f2e':'#ffffff',
    font:{color:t?'#cbd5e0':'#374151', family:'Inter,system-ui,sans-serif', size:12},
    title:{font:{size:15, color:t?'#e2e8f0':'#1f2937'}},
    colorway:['#0066cc','#00a86b','#ff6b35','#004499','#34d399','#fb923c',
              '#667eea','#764ba2','#e63946','#2ec4b6','#e9c46a','#264653'],
    xaxis:{gridcolor:t?'#2a2b3d':'#eef2f7', zerolinecolor:t?'#3a3b4d':'#dfe6ee', linecolor:t?'#3a3b4d':'#dfe6ee'},
    yaxis:{gridcolor:t?'#2a2b3d':'#eef2f7', zerolinecolor:t?'#3a3b4d':'#dfe6ee', linecolor:t?'#3a3b4d':'#dfe6ee'},
    legend:{font:{size:11}, bgcolor:'rgba(0,0,0,0)'} };
  extra=extra||{};
  for(var k in extra){
    if((k==='xaxis'||k==='yaxis'||k==='title'||k==='legend') &&
       extra[k] && typeof extra[k]==='object'){
      for(var kk in extra[k]) base[k][kk]=extra[k][kk];
    } else base[k]=extra[k];
  }
  return base;
}
var _owqPCFG={responsive:true, displayModeBar:true, displaylogo:false, scrollZoom:true,
  modeBarButtonsToRemove:['lasso2d','select2d'],
  toImageButtonOptions:{format:'svg', filename:'openwq_calibration_plot', scale:2}};
function _owqDraw(id, traces, extra){
  var el=document.getElementById(id);
  if(typeof Plotly==='undefined'){
    if(el) el.innerHTML='<div style="padding:1rem;color:#888;font-size:.85rem;">'+
      'Interactive chart could not load Plotly.js (needs internet).</div>';
    return;
  }
  Plotly.newPlot(id, traces, _owqLayout(extra), _owqPCFG);
  window._owqPlots.push({id:id, ex:extra});
}
function _owqRetheme(){ (window._owqPlots||[]).forEach(function(p){
  try{ Plotly.relayout(p.id, _owqLayout(p.ex)); }catch(e){} }); }
try{
  new MutationObserver(_owqRetheme).observe(document.documentElement,{attributes:true,attributeFilter:['data-theme']});
  new MutationObserver(_owqRetheme).observe(document.body,{attributes:true,attributeFilter:['data-theme']});
}catch(e){}

// ── Sortable result tables (click a column header to sort) ──────────────
function _owqMakeSortable(tbl){
  var head=tbl.tHead; if(!head) return;
  var ths=head.rows[0] ? [].slice.call(head.rows[0].cells) : [];
  ths.forEach(function(th, idx){
    if(th.classList.contains('no-sort')) return;
    th.style.cursor='pointer';
    th.title='Click to sort';
    if(!/[\\u21C5\\u25B2\\u25BC]/.test(th.textContent))
      th.insertAdjacentHTML('beforeend',' <span class="owq-sort-ind" style="opacity:.45;font-size:.85em;">\\u21C5</span>');
    th.addEventListener('click', function(){
      var tb=tbl.tBodies[0]; if(!tb) return;
      var rows=[].slice.call(tb.rows);
      var asc=th.getAttribute('data-asc')!=='1';
      rows.sort(function(a,b){
        var x=(a.cells[idx]?a.cells[idx].textContent:'').trim();
        var y=(b.cells[idx]?b.cells[idx].textContent:'').trim();
        var nx=parseFloat(x.replace(/[^0-9eE.+\\-]/g,''));
        var ny=parseFloat(y.replace(/[^0-9eE.+\\-]/g,''));
        var cmp=(!isNaN(nx)&&!isNaN(ny)) ? (nx-ny) : x.localeCompare(y);
        return asc?cmp:-cmp;
      });
      rows.forEach(function(r){ tb.appendChild(r); });
      ths.forEach(function(o){ o.removeAttribute('data-asc');
        var ind=o.querySelector('.owq-sort-ind'); if(ind) ind.textContent='\\u21C5'; });
      th.setAttribute('data-asc', asc?'1':'0');
      var myInd=th.querySelector('.owq-sort-ind'); if(myInd) myInd.textContent=asc?'\\u25B2':'\\u25BC';
    });
  });
}
document.addEventListener('DOMContentLoaded', function(){
  document.querySelectorAll('table.sortable').forEach(_owqMakeSortable);
});

// ── Per-figure "Lock zoom" + branded PDF export (FLUXOS-style) ──────────
function _owqAxisKeys(layout){
  var keys=[]; Object.keys(layout||{}).forEach(function(k){ if(/^(x|y)axis(\\d*)$/.test(k)) keys.push(k); });
  if(!keys.length) keys=['xaxis','yaxis']; return keys;
}
function _owqApplyLock(plot, locked){
  plot.dataset.locked = locked?'true':'false';
  var u={}; _owqAxisKeys(plot.layout||{}).forEach(function(k){ u[k+'.fixedrange']=!!locked; });
  try{ Plotly.relayout(plot, u); }catch(e){}
}
function _owqInstallWheel(plot){
  plot.addEventListener('wheel', function(e){
    if(plot.dataset.locked==='true') e.stopImmediatePropagation();
  }, {capture:true, passive:true});
}
var _OWQ_JSPDF_CDN='https://cdn.jsdelivr.net/npm/jspdf@2.5.1/dist/jspdf.umd.min.js';
function _owqLoadJsPDF(){
  if(window.jspdf&&window.jspdf.jsPDF) return Promise.resolve();
  return new Promise(function(res,rej){
    var s=document.createElement('script'); s.src=_OWQ_JSPDF_CDN; s.async=true;
    s.onload=function(){ (window.jspdf&&window.jspdf.jsPDF)?res():rej(new Error('jsPDF global missing')); };
    s.onerror=function(){ rej(new Error('Failed to load jsPDF from CDN')); };
    document.head.appendChild(s);
  });
}
function _owqFigTitle(plot){
  var n=plot.closest('.section')||plot.parentElement;
  while(n&&n!==document.body){ var h=n.querySelector('h2,h3'); if(h&&h.textContent.trim()) return h.textContent.trim(); n=n.parentElement; }
  return plot.id;
}
function _owqSlug(s){ return (s||'').replace(/[^a-z0-9_-]+/gi,'_').replace(/^_+|_+$/g,''); }
var _OPF={BG:'#ffffff',FG:'#1a1a1a',GRID:'#d4d4d4',FONT:'Inter,Helvetica,Arial,sans-serif',BASE:20,TICK:24,AXIS:28,HEAD:30};
function _owqExportLayout(src){
  var L=JSON.parse(JSON.stringify(src||{}));
  L.paper_bgcolor=_OPF.BG; L.plot_bgcolor=_OPF.BG;
  L.font=Object.assign({},L.font||{},{color:_OPF.FG,size:_OPF.BASE,family:_OPF.FONT});
  if(L.title){ if(typeof L.title==='string') L.title={text:L.title}; L.title.font=Object.assign({},L.title.font||{},{color:_OPF.FG,size:_OPF.HEAD,family:_OPF.FONT}); }
  L.margin=Object.assign({},L.margin||{}); L.margin.l=Math.max(L.margin.l||0,95); L.margin.r=Math.max(L.margin.r||0,50); L.margin.t=Math.max(L.margin.t||0,60); L.margin.b=Math.max(L.margin.b||0,95);
  Object.keys(L).forEach(function(k){ if(!/^(x|y)axis\\d*$/.test(k)) return; var ax=L[k]=Object.assign({},L[k]||{});
    ax.color=_OPF.FG; ax.linecolor=_OPF.FG; ax.gridcolor=_OPF.GRID; ax.zerolinecolor=_OPF.GRID; ax.tickcolor=_OPF.FG; ax.automargin=true;
    ax.tickfont=Object.assign({},ax.tickfont||{},{color:_OPF.FG,size:_OPF.TICK,family:_OPF.FONT});
    if(ax.title){ if(typeof ax.title==='string') ax.title={text:ax.title}; ax.title.font=Object.assign({},ax.title.font||{},{color:_OPF.FG,size:_OPF.AXIS,family:_OPF.FONT}); }
  });
  if(L.legend){ L.legend=Object.assign({},L.legend,{bgcolor:_OPF.BG,bordercolor:_OPF.GRID,font:Object.assign({},(L.legend.font||{}),{color:_OPF.FG,size:_OPF.BASE,family:_OPF.FONT})}); }
  return L;
}
async function _owqDownloadPDF(plot, btn){
  if(!window.Plotly){ alert('Plotly not loaded yet — try again.'); return; }
  var orig=btn?btn.innerHTML:null;
  try{
    if(btn){ btn.disabled=true; btn.innerHTML='\\u231B Building\\u2026'; }
    await _owqLoadJsPDF();
    var cssW=Math.max(plot.clientWidth||1000,800), cssH=Math.max(plot.clientHeight||500,460);
    var png=await Plotly.toImage({data:plot.data||[], layout:_owqExportLayout(plot.layout||{}), config:{displayModeBar:false,staticPlot:true}}, {format:'png', width:cssW*2, height:cssH*2});
    var jsPDF=window.jspdf.jsPDF; var landscape=cssW>=cssH;
    var pdf=new jsPDF({unit:'mm',format:'a4',orientation:landscape?'landscape':'portrait',compress:true});
    var pageW=pdf.internal.pageSize.getWidth(), pageH=pdf.internal.pageSize.getHeight(), margin=12;
    var proj=window.__OWQ_PROJECT_NAME||'OpenWQ', fig=_owqFigTitle(plot), gen=window.__OWQ_GENERATED_AT||new Date().toISOString();
    var BR=[0,102,204], INK=[26,26,26], MUT=[110,110,110]; var ly=margin+2;
    pdf.setFont('helvetica','bold'); pdf.setFontSize(20);
    pdf.setTextColor(INK[0],INK[1],INK[2]); pdf.text('open',margin,ly);
    var ow=pdf.getTextWidth('open'); pdf.setTextColor(BR[0],BR[1],BR[2]); pdf.text('WQ',margin+ow,ly);
    pdf.setFont('helvetica','bold'); pdf.setFontSize(12); pdf.setTextColor(INK[0],INK[1],INK[2]); pdf.text(proj,pageW-margin,ly,{align:'right'});
    pdf.setFont('helvetica','normal'); pdf.setFontSize(10); pdf.setTextColor(MUT[0],MUT[1],MUT[2]); pdf.text(fig,pageW-margin,ly+5,{align:'right'});
    pdf.setDrawColor(BR[0],BR[1],BR[2]); pdf.setLineWidth(0.4); var ruleY=margin+9; pdf.line(margin,ruleY,pageW-margin,ruleY);
    var headerH=12, footerH=8, availW=pageW-2*margin, availH=pageH-2*margin-headerH-footerH, ratio=cssW/cssH;
    var iw=availW, ih=availW/ratio; if(ih>availH){ ih=availH; iw=availH*ratio; }
    var ox=margin+(availW-iw)/2, oy=margin+headerH+(availH-ih)/2;
    pdf.addImage(png,'PNG',ox,oy,iw,ih,undefined,'FAST');
    var fY=pageH-margin/2; pdf.setFont('helvetica','bold'); pdf.setFontSize(9); pdf.setTextColor(INK[0],INK[1],INK[2]); pdf.text('open',margin,fY);
    var owf=pdf.getTextWidth('open'); pdf.setTextColor(BR[0],BR[1],BR[2]); pdf.text('WQ',margin+owf,fY);
    pdf.setFont('helvetica','normal'); pdf.setFontSize(8); pdf.setTextColor(MUT[0],MUT[1],MUT[2]);
    pdf.text('  \\u00b7  generated '+gen.slice(0,19).replace('T',' '), margin+owf+pdf.getTextWidth('WQ'), fY);
    pdf.text(plot.id, pageW-margin, fY, {align:'right'});
    pdf.setProperties({title:fig+' \\u2014 '+proj, subject:'OpenWQ figure export', author:proj, creator:'OpenWQ results report'});
    pdf.save((_owqSlug(proj)?_owqSlug(proj)+'_':'')+_owqSlug(plot.id||fig)+'.pdf');
  }catch(e){ console.error('[owq-pdf]',e); alert('Could not export PDF: '+(e.message||e)); }
  finally{ if(btn){ btn.disabled=false; btn.innerHTML=orig; } }
}
function _owqWireFigTools(){
  document.querySelectorAll('.lock-zoom').forEach(function(cb){
    var plot=document.getElementById(cb.getAttribute('data-target')); if(!plot) return;
    _owqInstallWheel(plot); _owqApplyLock(plot, cb.checked);
    cb.addEventListener('change', function(){ _owqApplyLock(plot, cb.checked); });
  });
  document.querySelectorAll('.pdf-btn').forEach(function(btn){
    if(btn.__owqWired) return; btn.__owqWired=true;
    btn.addEventListener('click', function(){ var plot=document.getElementById(btn.getAttribute('data-target')); if(plot) _owqDownloadPDF(plot, btn); });
  });
  document.querySelectorAll('.log-scale').forEach(function(cb){
    if(cb.__owqWired) return; cb.__owqWired=true;
    var plot=document.getElementById(cb.getAttribute('data-target')); if(!plot) return;
    cb.addEventListener('change', function(){
      var ax=cb.getAttribute('data-logaxes')||'y', t=cb.checked?'log':'linear', up={};
      if(ax.indexOf('y')>=0) up['yaxis.type']=t;
      if(ax.indexOf('x')>=0) up['xaxis.type']=t;
      if(window.Plotly) try{ Plotly.relayout(plot, up); }catch(e){}
    });
  });
}
// Charts are newPlot'd inline as the body parses; wire after a short tick so
// each plot's layout exists before we lock the axes.
if(document.readyState==='loading') document.addEventListener('DOMContentLoaded', function(){ setTimeout(_owqWireFigTools,250); });
else setTimeout(_owqWireFigTools,250);
</script>
"""


def _plotly_chart(div_id: str, traces: list, layout_extra: dict = None,
                  height: int = 440, card: bool = True, log_axes: str = None) -> str:
    """Emit a chart container + the Plotly.newPlot call (theme-aware), with a
    per-figure toolbar carrying a "Lock zoom" toggle and a branded "PDF"
    export button (same controls/aesthetics as the FLUXOS results report).

    ``log_axes`` (``'y'``, ``'x'`` or ``'xy'``) adds a "Log scale" checkbox that
    switches those axes between linear and logarithmic."""
    payload_t = json.dumps(traces)
    payload_l = json.dumps(layout_extra or {})
    _log_cb = (
        f'<label title="Toggle logarithmic axis">'
        f'<input type="checkbox" class="log-scale" data-target="{div_id}" '
        f'data-logaxes="{log_axes}"> Log scale</label>' if log_axes else '')
    toolbar = (
        '<div class="plot-controls">'
        '<span class="hint">Uncheck to scroll / drag-zoom &middot; double-click resets</span>'
        + _log_cb
        + f'<label title="Disable / enable scroll &amp; box zoom">'
        f'<input type="checkbox" class="lock-zoom" data-target="{div_id}" checked> '
        'Lock zoom</label>'
        f'<button type="button" class="pdf-btn" data-target="{div_id}" '
        f'title="Download this figure as a branded PDF (A4). First click '
        f'fetches jsPDF from the CDN, then it is cached.">&#128196; PDF</button>'
        '</div>')
    inner = (toolbar
             + f'<div id="{div_id}" class="plotly-chart" '
               f'style="width:100%;height:{height}px;"></div>'
               f'<script>_owqDraw({json.dumps(div_id)}, {payload_t}, {payload_l});</script>')
    return f'<div class="card">{inner}</div>' if card else inner


def _generate_convergence_plot(
    history: List[Dict],
    obj_fn: str,
    output_dir: str
) -> str:
    """Interactive convergence chart on the SAME scale the Summary reports.

    The optimiser minimises an internal objective; for KGE/NSE that is
    ``1 - metric``.  Plotting that raw objective made the chart's best read
    e.g. ``1.1489`` while the Summary reported the actual ``KGE = -0.1489``.
    Here we convert back (``metric = 1 - objective`` for KGE/NSE; error metrics
    pass through), so the plotted best matches the headline and "higher is
    better" for skill metrics, "lower is better" for error metrics."""
    _maximise = str(obj_fn).upper() in ("KGE", "NSE")

    def _to_metric(o):
        if not isinstance(o, (int, float)) or o != o or o >= _FAIL_PENALTY * 0.9:
            return None
        return (1.0 - o) if _maximise else o

    evals, raw_objs = [], []
    for h in history:
        eval_id = h.get("eval_id", len(evals))
        evals.append(eval_id)
        raw_objs.append(h.get("objective", h.get("best_objective", 0)))

    if not evals:
        return '<div class="card"><p>No evaluation history available.</p></div>'

    # If every evaluation hit the failure penalty there is nothing meaningful
    # to plot, so show a short note instead of a chart full of penalty values.
    if all((not isinstance(o, (int, float))) or o != o or o >= _FAIL_PENALTY * 0.9
           for o in raw_objs):
        return ('<div class="card"><p>&#9888;&#65039; Convergence plot omitted '
                '&mdash; every evaluation returned the failure penalty, so there '
                'is no successful model run to plot. See the diagnostic in the '
                '<a href="#summary" style="color:var(--primary);font-weight:600;">'
                'Results&nbsp;Summary</a>.</p></div>')

    metric_vals = [_to_metric(o) for o in raw_objs]

    # Overall best = the minimised-objective winner, expressed as the metric.
    _valid = [(i, o) for i, o in enumerate(raw_objs)
              if isinstance(o, (int, float)) and o == o and o < _FAIL_PENALTY * 0.9]
    best_idx = min(_valid, key=lambda t: t[1])[0]
    metric_best = _to_metric(raw_objs[best_idx])

    def _running_best(seq):
        """Running best of a metric sequence (None carries the previous best
        forward), honouring the maximise/minimise direction."""
        out, cur = [], None
        for v in seq:
            if v is not None:
                cur = v if cur is None else (max(cur, v) if _maximise
                                             else min(cur, v))
            out.append(cur)
        return out

    # Parallel-chain runs: one running-best line PER chain (so you can see each
    # chain converge and which one won), over a faint cloud of all evaluations.
    chain_ids = sorted({h["chain"] for h in history
                        if isinstance(h.get("chain"), int)})
    has_chains = len(chain_ids) > 1
    _PAL = ["#4d9ee8", "#fb923c", "#10b981", "#a78bfa", "#ef4444", "#14b8a6",
            "#f59e0b", "#ec4899", "#6366f1", "#84cc16", "#06b6d4", "#f43f5e"]
    _hov = f"eval %{{x}}<br>{obj_fn} %{{y:.4f}}<extra></extra>"

    if has_chains:
        traces = [{
            "type": "scatter", "mode": "markers", "x": evals, "y": metric_vals,
            "name": "All evaluations", "showlegend": False,
            "marker": {"size": 5, "color": "#9aa5b1", "opacity": 0.30},
            "hovertemplate": _hov}]
        for ci, cid in enumerate(chain_ids):
            pts = sorted((h.get("eval_id", 0),
                          _to_metric(h.get("objective", h.get("best_objective", 0))))
                         for h in history if h.get("chain") == cid)
            xs = [e for e, _ in pts]
            ys = _running_best([m for _, m in pts])
            traces.append({
                "type": "scatter", "mode": "lines", "x": xs, "y": ys,
                "name": f"Chain {cid + 1}",
                "line": {"color": _PAL[ci % len(_PAL)], "width": 2},
                "hovertemplate": (f"Chain {cid + 1}<br>eval %{{x}}<br>"
                                  f"best %{{y:.4f}}<extra></extra>")})
    else:
        traces = [
            {"type": "scatter", "mode": "markers", "x": evals, "y": metric_vals,
             "name": "All evaluations",
             "marker": {"size": 6, "color": "#4d9ee8", "opacity": 0.55},
             "hovertemplate": _hov},
            {"type": "scatter", "mode": "lines", "x": evals,
             "y": _running_best(metric_vals),
             "name": "Best so far", "line": {"color": "#10b981", "width": 2.5}},
        ]
    # Global best star (both modes).
    traces.append({
        "type": "scatter", "mode": "markers", "x": [evals[best_idx]],
        "y": [metric_best],
        "name": f"Best: {metric_best:.4f}",
        "marker": {"size": 16, "color": "#fb923c", "symbol": "star",
                   "line": {"color": "#fff", "width": 1}}})

    _better = "higher is better" if _maximise else "lower is better"
    layout = {"title": {"text": "Calibration Convergence"},
              "xaxis": {"title": "Evaluation Number"},
              "yaxis": {"title": f"{obj_fn} ({_better})"},
              "showlegend": True,
              "hovermode": "closest"}
    # A log axis can only show strictly-positive values.  KGE/NSE routinely go
    # negative (poor exploratory evaluations), and on a log scale those points
    # silently vanish — making the chart look empty.  So only offer the
    # log-scale toggle when every plotted value is > 0 (e.g. RMSE / |PBIAS|, or
    # an all-positive skill run); otherwise omit it so it can't blank the plot.
    _vm = [v for v in metric_vals if v is not None]
    _use_log = bool(_vm) and all(v > 0 for v in _vm)
    return _plotly_chart(_plot_id("conv"), traces, layout, height=430,
                         log_axes=('y' if _use_log else None))


def _generate_param_evolution_plots(
    history: List[Dict],
    parameters: List[Dict],
    output_dir: str
) -> str:
    """Interactive parameter-evolution chart with a dropdown to pick the
    parameter; shows its sampled value per evaluation plus its initial value
    and bounds as reference lines.  For DDS parallel-chain runs each parameter
    is drawn as ONE coloured line per chain (Chain 1, Chain 2, …)."""
    if not parameters:
        return ""

    # Detect a parallel-chain run: each evaluation may carry a 'chain' id.
    chain_ids = sorted({h["chain"] for h in history
                        if isinstance(h.get("chain"), int)})
    has_chains = len(chain_ids) > 1
    _PAL = ["#4d9ee8", "#fb923c", "#10b981", "#a78bfa", "#ef4444", "#14b8a6",
            "#f59e0b", "#ec4899", "#6366f1", "#84cc16", "#06b6d4", "#f43f5e"]

    # One "group" of traces per parameter; a group is N chain-lines (parallel
    # run) or a single line (sequential).  group_sizes drives the dropdown.
    traces, names, metas, group_sizes = [], [], [], []
    for i, param in enumerate(parameters):
        pname = param.get("name", f"param_{i}")
        recs = []  # (eval_id, value, chain)
        for h in history:
            pd_ = h.get("parameters", {})
            if pname in pd_:
                recs.append((h.get("eval_id", len(recs)), pd_[pname],
                             h.get("chain")))
        if not recs:
            continue
        first_param = (len(group_sizes) == 0)
        added = 0
        if has_chains:
            for ci, cid in enumerate(chain_ids):
                ev = [r[0] for r in recs if r[2] == cid]
                vals = [r[1] for r in recs if r[2] == cid]
                if not vals:
                    continue
                col = _PAL[ci % len(_PAL)]
                traces.append({
                    "type": "scatter", "mode": "markers+lines",
                    "x": ev, "y": vals, "name": f"Chain {cid + 1}",
                    "visible": first_param, "legendgroup": f"chain{cid}",
                    "marker": {"size": 5, "color": col},
                    "line": {"width": 1, "color": col},
                    "hovertemplate": (f"Chain {cid + 1}<br>eval %{{x}}<br>"
                                      f"%{{y:.5g}}<extra></extra>")})
                added += 1
        else:
            ev = [r[0] for r in recs]
            vals = [r[1] for r in recs]
            traces.append({
                "type": "scatter", "mode": "markers+lines", "x": ev, "y": vals,
                "name": pname, "visible": first_param,
                "marker": {"size": 6, "color": "#4d9ee8"},
                "line": {"width": 1, "color": "#4d9ee8"},
                "hovertemplate": "eval %{x}<br>%{y:.5g}<extra></extra>"})
            added += 1
        if added == 0:
            continue
        names.append(pname)
        metas.append({"initial": param.get("initial"),
                      "bounds": param.get("bounds")})
        group_sizes.append(added)
    if not traces:
        return ""

    n = len(traces)
    # First trace index of each parameter group (for dropdown toggling).
    starts, _acc = [], 0
    for g in group_sizes:
        starts.append(_acc)
        _acc += g

    def _shapes_for(i):
        sh = []
        meta = metas[i]
        ini = meta.get("initial")
        if isinstance(ini, (int, float)):
            sh.append({"type": "line", "xref": "paper", "x0": 0, "x1": 1,
                       "y0": ini, "y1": ini,
                       "line": {"color": "#9aa5b1", "width": 1, "dash": "dash"}})
        b = meta.get("bounds")
        if isinstance(b, (list, tuple)) and len(b) == 2:
            for bv in b:
                if isinstance(bv, (int, float)):
                    sh.append({"type": "line", "xref": "paper", "x0": 0, "x1": 1,
                               "y0": bv, "y1": bv,
                               "line": {"color": "#e74c3c", "width": 1,
                                        "dash": "dot"}})
        return sh

    # Per-parameter trace-visibility groups + reference-line shapes, fed to a
    # NATIVE <select> picker.  (Plotly's own updatemenus dropdown does not
    # scroll, so with 30+ parameters the lower ones get clipped/unreachable —
    # a native <select> scrolls like any other dropdown.)
    groups = [[starts[i], group_sizes[i]] for i in range(len(names))]
    shapes = [_shapes_for(i) for i in range(len(names))]

    layout = {
        "title": {"text": "Parameter Evolution"},
        "xaxis": {"title": "Evaluation"},
        "yaxis": {"title": names[0]},
        "shapes": _shapes_for(0),
        "showlegend": has_chains,
        "hovermode": "closest",
    }
    pid = _plot_id("evo")
    chart = _plotly_chart(pid, traces, layout, height=470, log_axes='y')

    # Native (scrollable) parameter picker wired to Plotly.restyle/relayout.
    opts = "\n".join(
        f'<option value="{i}">{html_lib.escape(str(pn))}</option>'
        for i, pn in enumerate(names))
    sel_id = f"{pid}-sel"
    selector = (
        f'<div style="margin:.4rem 0 .2rem;display:flex;align-items:center;'
        f'gap:.5rem;flex-wrap:wrap;">'
        f'<label for="{sel_id}" style="font-size:.85rem;color:var(--text2);'
        f'font-weight:600;">Parameter:</label>'
        f'<select id="{sel_id}" style="max-width:520px;padding:.35rem .5rem;'
        f'border:1px solid var(--border);border-radius:6px;'
        f'background:var(--surface);color:var(--text);font-size:.85rem;">'
        f'{opts}</select></div>')
    script = (
        '<script>(function(){'
        f'var sel=document.getElementById({json.dumps(sel_id)});if(!sel)return;'
        f'var pid={json.dumps(pid)},groups={json.dumps(groups)},'
        f'shapes={json.dumps(shapes)},titles={json.dumps(names)},n={n};'
        'sel.addEventListener("change",function(){'
        'var i=parseInt(this.value,10)||0,vis=[],t;'
        'for(t=0;t<n;t++){vis.push(false);}'
        'for(t=groups[i][0];t<groups[i][0]+groups[i][1];t++){vis[t]=true;}'
        'if(window.Plotly){'
        'Plotly.restyle(pid,{visible:vis});'
        'Plotly.relayout(pid,{"yaxis.title.text":titles[i],"shapes":shapes[i]});'
        '}});'
        '})();</script>')
    return selector + chart + script


def _generate_correlation_plot(
    history: List[Dict],
    parameters: List[Dict],
    output_dir: str
) -> str:
    """Interactive parameter-correlation heatmap (Pearson r)."""
    param_names = [p.get("name", f"p{i}") for i, p in enumerate(parameters)]
    n_evals = len(history)
    n_params = len(param_names)
    if n_evals < 3 or n_params < 2 or not HAS_NUMPY:
        return ('<div class="card"><p>Insufficient data for correlation '
                'analysis (need &ge;3 evaluations and &ge;2 parameters).</p></div>')

    matrix = np.full((n_evals, n_params), np.nan)
    for i, h in enumerate(history):
        pd_ = h.get("parameters", {})
        for j, pname in enumerate(param_names):
            if pname in pd_:
                matrix[i, j] = pd_[pname]

    valid_cols = ~np.all(np.isnan(matrix), axis=0)
    if int(np.sum(valid_cols)) < 2:
        return ('<div class="card"><p>Insufficient valid parameter data for '
                'correlation analysis.</p></div>')
    matrix = matrix[:, valid_cols]
    valid_names = [n for n, v in zip(param_names, valid_cols) if v]
    n_valid = matrix.shape[1]

    corr = np.eye(n_valid)
    for i in range(n_valid):
        for j in range(i + 1, n_valid):
            mask = ~(np.isnan(matrix[:, i]) | np.isnan(matrix[:, j]))
            if int(np.sum(mask)) > 2:
                r = np.corrcoef(matrix[mask, i], matrix[mask, j])[0, 1]
                corr[i, j] = r
                corr[j, i] = r

    z = [[float(corr[i][j]) for j in range(n_valid)] for i in range(n_valid)]
    # No in-cell value labels — just the colour scale (exact r is in the hover).
    traces = [{
        "type": "heatmap", "z": z, "x": valid_names, "y": valid_names,
        "colorscale": "RdBu", "reversescale": True,
        "zmid": 0, "zmin": -1, "zmax": 1,
        "colorbar": {"title": "Pearson r"},
        "hovertemplate": "%{x}<br>%{y}<br>r = %{z:.3f}<extra></extra>"}]
    layout = {"title": {"text": "Parameter Correlations"},
              "xaxis": {"tickangle": -45, "automargin": True},
              "yaxis": {"automargin": True, "autorange": "reversed"}}
    return _plotly_chart(_plot_id("corr"), traces, layout,
                         height=max(420, n_valid * 36 + 140))


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
            fontsize=11, color='#1f2937')
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
                      fontsize=11, color='#1f2937')
        ax2.grid(True, alpha=0.3)
        ax2.set_aspect('equal', adjustable='box')

        # ── Legend vs colorbar ──
        if show_legend and n_reaches > 1:
            # Compact per-reach legend, deduplicated (one entry per reach).
            handles, labels = ax2.get_legend_handles_labels()
            # tab20 distinct labels per reach + the 1:1 entry; cap at 20.
            ax2.legend(handles, labels, loc='upper left',
                       facecolor='#ffffff', edgecolor='#d1d5db',
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
                           color='#4b5563', fontsize=8)
            cbar.ax.tick_params(labelsize=7)
        else:
            # Single reach — just show a basic legend.
            ax1.legend(facecolor='#ffffff', edgecolor='#d1d5db', fontsize=8)
            ax2.legend(facecolor='#ffffff', edgecolor='#d1d5db', fontsize=8)

        fig.suptitle(f'Performance: {species}', fontsize=14,
                     fontweight='bold', color='#1f2937', y=1.02)
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
                    fontsize=11, color='#1f2937')
                rax1.set_ylabel('Concentration')
                rax1.grid(True, alpha=0.3)
                rax1.legend(facecolor='#ffffff', edgecolor='#d1d5db', fontsize=8)

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
                    fontsize=11, color='#1f2937')
                rax2.set_aspect('equal', adjustable='box')
                rax2.grid(True, alpha=0.3)
                rax2.legend(facecolor='#ffffff', edgecolor='#d1d5db', fontsize=8)

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
                              color='#1f2937', y=1.02)
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


def _gray_overlay(all_sim, species, pd, max_g=300, max_pts=200):
    """ONE light-gray WebGL trace overlaying every evaluation's simulated
    series for ``species`` (None-separated per eval/reach).  Per-point
    ``customdata`` carries the evaluation number, so hovering shows ``sim N``
    and a click can read which run it is.  A single trace regardless of eval
    count → efficient and conflict-free.  Returns ``(trace_or_None, eval_ids)``.
    """
    if not (isinstance(all_sim, pd.DataFrame) and not all_sim.empty
            and 'eval_id' in all_sim.columns):
        return None, []
    # Split all_sim by species ONCE per frame (it can be 10M+ rows) and cache it
    # on the frame's id: re-running `.astype(str) == species` on every species
    # call was a second multi-minute cost (each full pass is ~30s under Rosetta).
    _cache = _gray_overlay.__dict__.setdefault('_by_sp', {})
    _by_sp = _cache.get(id(all_sim))
    if _by_sp is None:
        _by_sp = {str(k): g for k, g in
                  all_sim.groupby(all_sim['species'].astype(str), sort=False)}
        _cache.clear()                 # keep only the most recent frame's split
        _cache[id(all_sim)] = _by_sp
    asp = _by_sp.get(str(species))
    if asp is None or asp.empty:
        return None, []
    eids = sorted(asp['eval_id'].unique())
    if len(eids) > max_g:                           # sample to keep it light
        st = len(eids) / max_g
        eids = [eids[int(k * st)] for k in range(max_g)]
    asp = asp[asp['eval_id'].isin(eids)]
    has_reach = 'reach_id' in asp.columns
    gx, gy, gcd = [], [], []
    for eid in eids:
        ie = int(eid)
        g0 = asp[asp['eval_id'] == eid]
        rgroups = ([grp for _, grp in g0.groupby('reach_id')]
                   if has_reach else [g0])
        for g in rgroups:
            if 'datetime' in g.columns:
                g = g.sort_values('datetime')
                xs = [str(v) for v in pd.to_datetime(g['datetime'])]
            else:
                xs = list(range(len(g)))
            ys = [None if v != v else float(v) for v in g['simulated']]
            if len(xs) > max_pts:                   # thin long series
                step = len(xs) // max_pts + 1
                xs = xs[::step]
                ys = ys[::step]
            gx += xs + [None]
            gy += ys + [None]
            gcd += [ie] * len(xs) + [None]
    if not gx:
        return None, []
    trace = {"type": "scattergl", "mode": "lines", "x": gx, "y": gy,
             "customdata": gcd, "name": f"all simulations ({len(eids)})",
             "_role": "graysim", "legendgroup": "allsim", "showlegend": True,
             "line": {"color": "#c9ced6", "width": 0.6}, "opacity": 0.5,
             "hovertemplate":
                 "sim %{customdata}<br>%{x}<br>%{y:.4g}<extra></extra>"}
    return trace, [int(e) for e in eids]


def _isolated_trace():
    """Empty placeholder trace (role 'isolated') that the JS fills in when the
    user isolates a single simulation — kept on top, in a darker gray (it is
    just one of the gray runs, not the best fit, so it stays gray)."""
    return {"type": "scatter", "mode": "lines", "x": [], "y": [],
            "name": "isolated simulation", "_role": "isolated",
            "visible": False, "showlegend": True,
            "line": {"color": "#6b7280", "width": 1.8},
            "hovertemplate": "%{x}<br>%{y:.4g}<extra>isolated</extra>"}


def _sim_controls(chart_id, y_values, eval_ids):
    """``(yaxis_dict, controls_html)`` for a chart that carries a gray overlay.

    Controls: a "Zoom to best fit" button (the chart also OPENS at that range),
    an "Isolate" dropdown (highlights one run, hides the gray spaghetti), and a
    "Show all simulations" button.  Traces are found by a ``_role`` tag, never
    by index, so the observed/best-fit traces are never touched.  Clicking a
    gray line also isolates that run (read from its per-point customdata)."""
    yaxis = {"title": "Concentration"}
    idj = json.dumps(chart_id)
    bs = ("padding:.25rem .55rem;font-size:.78rem;border:1px solid "
          "var(--border);border-radius:6px;background:var(--surface);"
          "color:var(--text);cursor:pointer;margin-right:.4rem;")
    vals = [v for v in (y_values or []) if v is not None]
    zoom_btn = ""
    if vals:
        lo, hi = min(vals), max(vals)
        pad = (hi - lo) * 0.08 if hi > lo else (abs(hi) * 0.08 or 0.1)
        y0, y1 = round(lo - pad, 6), round(hi + pad, 6)
        yaxis["range"] = [y0, y1]
        zr = json.dumps({"yaxis.range": [y0, y1], "xaxis.autorange": True})
        zoom_btn = ("<button type='button' style='" + bs + "' onclick='"
                    "if(window.Plotly)Plotly.relayout(" + idj + "," + zr
                    + ")'>&#128269; Zoom to best fit</button>")
    # "Zoom to active simulations": fit the y-axis to whatever traces are
    # CURRENTLY visible (best fit, an isolated run, all of them, obs — whatever
    # the user has toggled on), so the view tracks the live state of the graph.
    act_id = chart_id + "-zact"
    active_btn = ("<button id='" + act_id + "' type='button' style='" + bs
                  + "'>&#128269; Zoom to active simulations</button>")
    sel_id = chart_id + "-iso"
    all_id = chart_id + "-all"
    opts = "<option value=''>— all simulations —</option>" + "".join(
        "<option value='%d'>sim %d</option>" % (e, e) for e in eval_ids)
    dropdown = ("<label style='font-size:.78rem;color:var(--text2);'>Isolate: "
                "<select id='" + sel_id + "' style='" + bs
                + "max-width:150px;margin-left:.3rem;'>" + opts
                + "</select></label>")
    all_btn = ("<button id='" + all_id + "' type='button' style='" + bs
               + "'>&#8617; Show all simulations</button>")
    script = (
        "<script>(function(){var id=" + idj + ";"
        "function gd(){return document.getElementById(id);}"
        "function ti(g,r){var d=g.data||[];for(var i=0;i<d.length;i++)"
        "if(d[i]._role===r)return i;return -1;}"
        "function iso(v){var g=gd();if(!window.Plotly||!g)return;"
        "var gi=ti(g,'graysim'),hi=ti(g,'isolated');if(gi<0)return;"
        "if(!v){Plotly.restyle(id,{visible:true},[gi]);"
        "if(hi>=0)Plotly.restyle(id,{visible:false},[hi]);return;}"
        "if(hi<0)return;var t=g.data[gi],cd=t.customdata||[],xs=[],ys=[];"
        "for(var i=0;i<cd.length;i++){if(String(cd[i])===String(v)){"
        "xs.push(t.x[i]);ys.push(t.y[i]);}}"
        "Plotly.restyle(id,{x:[xs],y:[ys],name:'sim '+v,visible:true},[hi]);"
        "Plotly.restyle(id,{visible:false},[gi]);}"
        "var sel=document.getElementById(" + json.dumps(sel_id) + ");"
        "if(sel)sel.onchange=function(){iso(this.value);};"
        "var ab=document.getElementById(" + json.dumps(all_id) + ");"
        "if(ab)ab.onclick=function(){if(sel)sel.value='';iso('');"
        "if(window.Plotly)Plotly.relayout(id,{'yaxis.autorange':true,"
        "'xaxis.autorange':true});};"
        # Zoom to active sims: union y-extent of every currently-visible trace.
        "var azb=document.getElementById(" + json.dumps(act_id) + ");"
        "if(azb)azb.onclick=function(){var g=gd();if(!window.Plotly||!g)return;"
        "var d=g.data||[],lo=Infinity,hi=-Infinity;"
        "for(var i=0;i<d.length;i++){var t=d[i];"
        "if(t.visible===false||t.visible==='legendonly')continue;"
        "var ys=t.y||[];for(var j=0;j<ys.length;j++){var v=ys[j];"
        "if(v==null||v!==v)continue;if(v<lo)lo=v;if(v>hi)hi=v;}}"
        "if(lo===Infinity)return;"
        "var pad=(hi>lo)?(hi-lo)*0.08:(Math.abs(hi)*0.08||0.1);"
        "Plotly.relayout(id,{'yaxis.range':[lo-pad,hi+pad],"
        "'xaxis.autorange':true});};"
        "function at(){var g=gd();if(!g||!g.on){setTimeout(at,300);return;}"
        "g.on('plotly_click',function(ev){if(!ev||!ev.points||!ev.points.length)"
        "return;var p=ev.points[0],g2=gd(),d=g2.data||[],tr=d[p.curveNumber];"
        "if(!tr)return;"
        "if(tr._role==='isolated'){if(sel)sel.value='';iso('');return;}"
        "if(tr._role!=='graysim')return;"
        "var e=p.customdata;if(e==null)return;"
        "if(sel)sel.value=String(e);iso(String(e));});}at();"
        "})();</script>")
    controls = ("<div style='margin:.1rem 0 .3rem;display:flex;"
                "align-items:center;gap:.5rem;flex-wrap:wrap;'>"
                + zoom_btn + active_btn + dropdown + all_btn + "</div>" + script)
    return yaxis, controls


def _generate_timeseries_charts(matched_data, output_dir,
                                hostmodel: str = "mizuroute",
                                simulated_data=None, of_species=None,
                                all_sim=None, best_eval=None) -> str:
    """Interactive (Plotly) observed-vs-simulated time series + obs-vs-sim
    scatter, one block per calibrated species.  Each reach/HRU is a legend
    entry the user can toggle; the mode-bar provides SVG export.

    Objective species that have simulated output but NO observations at the
    target reach(es) get a simulated-only time series with a clear note (so a
    selected species is never silently missing from the report)."""
    try:
        import pandas as pd
    except ImportError:
        return ""
    _has_matched = (isinstance(matched_data, pd.DataFrame)
                    and not matched_data.empty)
    _has_sim = (isinstance(simulated_data, pd.DataFrame)
                and not simulated_data.empty)
    if not _has_matched and not _has_sim:
        return ""

    _is_summa = (hostmodel or "").lower() == "summa"
    feat_label = "HRU" if _is_summa else "Reach"
    feat_label_plural = "HRUs" if _is_summa else "reaches"
    has_reach = _has_matched and 'reach_id' in matched_data.columns
    species_list = (list(matched_data['species'].unique())
                    if _has_matched and 'species' in matched_data.columns
                    else [])
    matched_species = {str(s) for s in species_list}
    # "best fit" legend label, named with its simulation number when known so
    # it ties to the gray-overlay "sim N" lines.
    _best_lbl = (f"best fit (sim {int(best_eval)})"
                 if best_eval is not None else "best fit")
    # Per-station/reach colour palette.  When more than one reach/HRU is shown,
    # each station gets its own colour, used for BOTH its observations and its
    # best-fit line on the left time-series AND its markers on the right
    # obs-vs-sim scatter, so a station is the same colour across both panels.
    _STN_PAL = ["#636EFA", "#EF553B", "#00CC96", "#AB63FA", "#FFA15A",
                "#19D3F3", "#FF6692", "#B6E880", "#FF97FF", "#FECB52",
                "#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD",
                "#8C564B", "#E377C2", "#BCBD22", "#17BECF", "#393B79",
                "#637939", "#8C6D31", "#843C39", "#7B4173"]
    blocks = []

    for species in species_list:
        sp = matched_data[matched_data['species'] == species]
        if sp.empty:
            continue
        if has_reach:
            reach_ids = sorted(sp['reach_id'].dropna().unique(),
                               key=lambda x: (str(type(x).__name__), x))
        else:
            reach_ids = [None]
        n_reaches = len(reach_ids)

        ts_traces, sc_traces, all_v = [], [], []
        _multi = n_reaches > 1
        # Pre-slice the (potentially 10M+ row) simulated frame to THIS species
        # ONCE and group it by reach ONCE.  The old code re-filtered the FULL
        # frame with `.astype(str) == species` (and again per reach) INSIDE the
        # reach loop -> O(reaches x rows) string conversions, the report's main
        # cost (minutes under Rosetta on long multi-reach runs).
        _sp_groups, _sp_sim_all = {}, None
        if _has_sim and simulated_data is not None and not simulated_data.empty:
            _sp_sim_all = simulated_data[
                simulated_data['species'].astype(str) == str(species)]
            if 'reach_id' in _sp_sim_all.columns:
                _sp_groups = {str(k): g for k, g in _sp_sim_all.groupby(
                    _sp_sim_all['reach_id'].astype(str), sort=False)}
        for _i, rid in enumerate(reach_ids):
            rd = sp if rid is None else sp[sp['reach_id'] == rid]
            if rd.empty:
                continue
            # One colour per station/reach (only when several are shown).
            # A single reach keeps the clear observed=red / best-fit=green scheme.
            _c = _STN_PAL[_i % len(_STN_PAL)]
            _obs_c = _c if _multi else "#ef4444"
            _sim_c = _c if _multi else "#10b981"
            _obs_edge = "rgba(0,0,0,.45)" if _multi else "#7f1d1d"
            if 'datetime' in rd.columns:
                rd = rd.sort_values('datetime')
                x = [str(v) for v in pd.to_datetime(rd['datetime'])]
            else:
                x = list(range(len(rd)))
            obs = [None if v != v else float(v) for v in rd['observed']]
            sim = [None if v != v else float(v) for v in rd['simulated']]
            # Best-fit LINE: draw the best eval's FULL simulated series (the
            # same representation as the gray overlay) so the highlighted best
            # is directly comparable to the gray runs — NOT a sparse zig-zag
            # through observation dates only (which always looks worse than a
            # full hydrograph).  Observations stay as markers and the scatter
            # still uses the matched obs-date pairs.  Falls back to the matched
            # simulated values when the full series isn't available.
            x_line, sim_line = x, sim
            if _has_sim and _sp_sim_all is not None:
                _sdr = (_sp_groups.get(str(rid)) if rid is not None
                        else _sp_sim_all)
                if _sdr is not None and not _sdr.empty and 'datetime' in _sdr.columns:
                    _sdr = _sdr.sort_values('datetime')
                    _xl = [str(v) for v in pd.to_datetime(_sdr['datetime'])]
                    _yl = [None if v != v else float(v)
                           for v in _sdr['simulated']]
                    if len(_xl) > 4000:                  # thin a long line
                        _st = len(_xl) // 4000 + 1
                        _xl, _yl = _xl[::_st], _yl[::_st]
                    x_line, sim_line = _xl, _yl
            label = f'{feat_label} {rid}' if rid is not None else 'all'
            # Short legend labels when there's a single reach (the reach id is
            # redundant then) so the horizontal top legend fits on one row.
            _sfx = f" · {label}" if n_reaches > 1 else ""
            ts_traces.append({
                "type": "scatter", "mode": "markers", "x": x, "y": obs,
                "name": f"observed{_sfx}", "legendgroup": f"obs-{label}",
                "marker": {"size": 7, "symbol": "circle", "color": _obs_c,
                           "line": {"width": 0.6, "color": _obs_edge}},
                "hovertemplate": "%{x|%Y-%m-%d}<br>obs %{y:.4g}<extra></extra>"})
            ts_traces.append({
                "type": "scatter", "mode": "lines", "x": x_line, "y": sim_line,
                "name": f"{_best_lbl}{_sfx}", "legendgroup": f"best-{label}",
                "line": {"width": 2, "color": _sim_c},
                "hovertemplate":
                    "%{x|%Y-%m-%d}<br>best %{y:.4g}<extra></extra>"})
            sc_traces.append({
                "type": "scatter", "mode": "markers", "x": obs, "y": sim,
                "name": label,
                "marker": {"size": 7, "opacity": 0.6, "color": _c},
                "hovertemplate": "obs %{x:.4g}<br>sim %{y:.4g}<extra></extra>"})
            all_v += [v for v in obs if v is not None]
            all_v += [v for v in sim_line if v is not None]

        if not ts_traces:
            continue
        if all_v:
            vmin, vmax = min(all_v), max(all_v)
            m = (vmax - vmin) * 0.05 if vmax > vmin else (abs(vmax) * 0.05 or 0.05)
            sc_traces.append({
                "type": "scatter", "mode": "lines",
                "x": [vmin - m, vmax + m], "y": [vmin - m, vmax + m],
                "line": {"dash": "dash", "color": "#9aa5b1", "width": 1},
                "name": "1:1 line", "hoverinfo": "skip"})

        # Light-gray overlay: ONE combined trace (hover shows the sim number),
        # behind the obs/best lines; an empty 'isolated' trace sits on top for
        # the isolate dropdown / click-to-isolate.
        _gray, _evals = _gray_overlay(all_sim, species, pd)
        ts_id = _plot_id("ts")
        if _gray is not None:
            ts_traces = [_gray] + ts_traces + [_isolated_trace()]
            _yaxis, _zoom_btn = _sim_controls(ts_id, all_v, _evals)
        else:
            _yaxis, _zoom_btn = {"title": "Concentration"}, ""

        _leg = n_reaches <= 25
        _outlet = " · outlet obs" if _is_summa else ""
        ts_layout = {
            "title": {"text": f"{species} — time series "
                              f"({n_reaches} {feat_label_plural}{_outlet})"},
            "xaxis": {"title": "Date"}, "yaxis": _yaxis,
            "hovermode": "closest", "showlegend": True,
            "legend": {"orientation": "h", "yanchor": "bottom", "y": 1.02,
                       "xanchor": "left", "x": 0},
            "margin": {"l": 64, "r": 24, "t": 92, "b": 52}}
        sc_layout = {
            "title": {"text": f"{species} — observed vs simulated"},
            "xaxis": {"title": "Observed"}, "yaxis": {"title": "Simulated"},
            "hovermode": "closest", "showlegend": _leg,
            "legend": {"orientation": "h", "yanchor": "bottom", "y": 1.02,
                       "xanchor": "left", "x": 0},
            "margin": {"l": 64, "r": 24, "t": 92, "b": 52}}
        c1 = _plotly_chart(ts_id, ts_traces, ts_layout, height=440, card=False,
                           log_axes='y')
        c2 = _plotly_chart(_plot_id("sc"), sc_traces, sc_layout,
                           height=440, card=False, log_axes='xy')
        blocks.append(
            f'<div class="card">'
            f'<h3 style="margin:.2rem 0 .7rem;">{html_lib.escape(str(species))} '
            f'&mdash; {n_reaches} '
            f'{feat_label if n_reaches == 1 else feat_label_plural}</h3>'
            f'{_zoom_btn}'
            f'<div style="display:flex;flex-wrap:wrap;gap:1rem;">'
            f'<div style="flex:1 1 520px;min-width:300px;">{c1}</div>'
            f'<div style="flex:1 1 360px;min-width:280px;">{c2}</div>'
            f'</div></div>')

    # ── Simulated-only blocks for objective species lacking observations ──
    # (e.g. NH4-N selected for the OF but no obs at the target reach).
    if _has_sim and of_species:
        _sim_has_reach = 'reach_id' in simulated_data.columns
        for species in of_species:
            if str(species) in matched_species:
                continue  # already shown with observations above
            sp = simulated_data[
                simulated_data['species'].astype(str) == str(species)]
            if sp.empty:
                continue
            if _sim_has_reach:
                reach_ids = sorted(sp['reach_id'].dropna().unique(),
                                   key=lambda x: (str(type(x).__name__), x))
            else:
                reach_ids = [None]
            ts_traces = []
            sim_v = []
            for rid in reach_ids:
                rd = sp if rid is None else sp[sp['reach_id'] == rid]
                if rd.empty:
                    continue
                if 'datetime' in rd.columns:
                    rd = rd.sort_values('datetime')
                    x = [str(v) for v in pd.to_datetime(rd['datetime'])]
                else:
                    x = list(range(len(rd)))
                sim = [None if v != v else float(v) for v in rd['simulated']]
                sim_v += [v for v in sim if v is not None]
                label = f'{feat_label} {rid}' if rid is not None else 'all'
                ts_traces.append({
                    "type": "scatter", "mode": "lines", "x": x, "y": sim,
                    "name": f"{_best_lbl} · {label}",
                    "line": {"width": 2, "color": "#10b981"},
                    "hovertemplate":
                        "%{x|%Y-%m-%d}<br>best %{y:.4g}<extra></extra>"})
            if not ts_traces:
                continue
            # Same treatment as the obs species: combined gray overlay (hover
            # shows sim number) + isolate dropdown/click + show-all/zoom.
            _gray, _evals = _gray_overlay(all_sim, species, pd)
            _reach_txt = ", ".join(str(r) for r in reach_ids if r is not None)
            ts_id = _plot_id("tsonly")
            if _gray is not None:
                ts_traces = [_gray] + ts_traces + [_isolated_trace()]
                _yaxis, _zoom_btn = _sim_controls(ts_id, sim_v, _evals)
            else:
                _yaxis, _zoom_btn = {"title": "Concentration"}, ""
            ts_layout = {
                "title": {"text": f"{species} — simulated (no observations)"},
                "xaxis": {"title": "Date"},
                "yaxis": _yaxis,
                "hovermode": "closest", "showlegend": True,
                "legend": {"orientation": "h", "yanchor": "bottom", "y": 1.02,
                           "xanchor": "left", "x": 0},
                "margin": {"l": 64, "r": 24, "t": 92, "b": 52}}
            c1 = _plotly_chart(ts_id, ts_traces, ts_layout,
                               height=440, card=False, log_axes='y')
            note = (
                f'<div style="margin:.2rem 0 .7rem;padding:.55rem .75rem;'
                f'border-radius:8px;background:rgba(245,158,11,.12);'
                f'border-left:4px solid #f59e0b;font-size:.85rem;'
                f'color:var(--text);line-height:1.5;">'
                f'<strong style="color:#d97706;">&#9888; No observation '
                f'data</strong> for <code>{html_lib.escape(str(species))}</code>'
                f' at the target {feat_label_plural} '
                f'({html_lib.escape(_reach_txt) or "—"}). Showing the '
                f'<strong>simulated</strong> series only &mdash; this species '
                f'did <strong>not</strong> contribute to the calibration '
                f'objective (nothing to match against here).</div>')
            blocks.append(
                f'<div class="card">'
                f'<h3 style="margin:.2rem 0 .7rem;">'
                f'{html_lib.escape(str(species))} &mdash; simulated only</h3>'
                f'{note}{_zoom_btn}<div>{c1}</div></div>')

    return "\n".join(blocks)


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
        <table class="param-table sortable">
        <thead>
            <tr>
                <th>Rank</th>
                <th>Parameter</th>
                <th>μ*</th>
                <th>σ</th>
                <th class="no-sort">Importance</th>
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
        <table class="param-table sortable">
        <thead>
            <tr>
                <th>Rank</th>
                <th>Parameter</th>
                <th>S1</th>
                <th>ST</th>
                <th class="no-sort">Importance</th>
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
    """Interactive 2-panel Morris summary: μ* ranking bar chart +
    μ* vs σ scatter (above the 1:1 line ⇒ interactions / non-linearity)."""
    names = list(param_names)
    mu = [float(x) for x in mu_star]
    sg = [float(x) for x in sigma]
    n = len(names)
    order = sorted(range(n), key=lambda i: mu[i])   # ascending → biggest on top

    bar = [{"type": "bar", "orientation": "h",
            "x": [mu[i] for i in order], "y": [names[i] for i in order],
            "marker": {"color": "#4d9ee8"},
            "hovertemplate": "%{y}<br>μ* = %{x:.4g}<extra></extra>"}]
    bar_layout = {"title": {"text": "Importance ranking (μ*)"},
                  "xaxis": {"title": "μ* (overall importance)"},
                  "yaxis": {"automargin": True}, "showlegend": False}

    mx = max(mu) if mu else 1.0
    scat = [{"type": "scatter", "mode": "markers", "x": mu, "y": sg,
             "text": names, "marker": {"size": 11, "color": "#10b981"},
             "name": "parameters",
             "hovertemplate": "%{text}<br>μ* = %{x:.4g}<br>σ = %{y:.4g}<extra></extra>"}]
    if mx > 0:
        scat.append({"type": "scatter", "mode": "lines", "x": [0, mx],
                     "y": [0, mx], "name": "1:1",
                     "line": {"dash": "dash", "color": "#9aa5b1", "width": 1},
                     "hoverinfo": "skip"})
    scat_layout = {"title": {"text": "μ* vs σ  (above 1:1 ⇒ interactions)"},
                   "xaxis": {"title": "μ* (mean |elementary effect|)"},
                   "yaxis": {"title": "σ (std of elementary effect)"},
                   "showlegend": False, "hovermode": "closest"}

    _h = max(360, n * 28 + 130)
    c1 = _plotly_chart(_plot_id("morris-bar"), bar, bar_layout, height=_h, card=False)
    c2 = _plotly_chart(_plot_id("morris-sc"), scat, scat_layout,
                       height=max(360, _h), card=False)
    return (f'<div class="card"><div style="display:flex;flex-wrap:wrap;gap:1rem;">'
            f'<div style="flex:1 1 400px;min-width:300px;">{c1}</div>'
            f'<div style="flex:1 1 340px;min-width:280px;">{c2}</div>'
            f'</div></div>')


def _generate_sobol_plot(
    param_names: List[str],
    S1: List[float],
    ST: List[float],
    output_dir: str
) -> str:
    """Interactive grouped-bar Sobol indices chart (S1 + ST)."""
    names = list(param_names)
    traces = [
        {"type": "bar", "x": names, "y": [float(v) for v in S1],
         "name": "S1 (First-order)", "marker": {"color": "#4d9ee8"},
         "hovertemplate": "%{x}<br>S1 = %{y:.4f}<extra></extra>"},
        {"type": "bar", "x": names, "y": [float(v) for v in ST],
         "name": "ST (Total-order)", "marker": {"color": "#34d399"},
         "hovertemplate": "%{x}<br>ST = %{y:.4f}<extra></extra>"},
    ]
    layout = {"title": {"text": "Sobol Sensitivity Indices"},
              "barmode": "group",
              "xaxis": {"tickangle": -45, "automargin": True},
              "yaxis": {"title": "Sensitivity index"},
              "hovermode": "closest"}
    return _plotly_chart(_plot_id("sobol"), traces, layout, height=470)


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


def _generate_perf_distribution_plots(
    perf_metrics: List[Dict],
    hostmodel: str = "mizuroute",
) -> str:
    """Cumulative + non-cumulative distribution of each metric across the
    reaches/HRUs that carry observations.

    For every performance metric (KGE, NSE, R2, RMSE, PBIAS) we draw a
    two-panel figure:

      * left  — **non-cumulative** distribution (histogram: how many
        reaches/HRUs fall in each value band; dashed line marks the median);
      * right — **cumulative** distribution (empirical CDF: the fraction of
        reaches/HRUs at or below a given value).

    Multiple calibrated species are overlaid (one colour each).  Drawing a
    distribution needs at least **two** stations (reaches/HRUs) with
    observations — with a single station the section is kept but replaced by a
    short note explaining that two or more are required.
    """
    import math

    _is_summa = (hostmodel or "").lower() == "summa"
    feat_singular = "HRU" if _is_summa else "reach"
    feat_plural = "HRUs" if _is_summa else "reaches"

    heading = (
        '<h3 style="margin-top:1.5rem;">Spatial distribution of performance '
        f'metrics across {feat_plural}</h3>'
    )

    if not perf_metrics:
        return ""

    # Number of distinct stations (reaches/HRUs) with observations.  When no
    # ``reach_id`` is present the setup is a single aggregated station.
    has_reach = any(m.get("reach_id") is not None for m in perf_metrics)
    if has_reach:
        n_stations = len({m.get("reach_id") for m in perf_metrics
                          if m.get("reach_id") is not None})
    else:
        n_stations = 1

    # A distribution is only meaningful with >= 2 stations.  The charts are drawn
    # client-side with Plotly (interactive, like the rest of the report), so no
    # server-side plotting library is needed.
    if n_stations < 2:
        msg = (
            'These <strong>cumulative</strong> and '
            '<strong>non-cumulative</strong> distribution plots summarise '
            f'how each metric varies across the {feat_plural} that have '
            'observations. <strong>At least two '
            f'{feat_plural} with observations are required</strong> to draw '
            f'them &mdash; this calibration uses a single {feat_singular}, so '
            'the distributions are not shown. Add observations at one or more '
            f'additional {feat_plural} to enable this view.'
        )
        return heading + (
            f'<div class="highlight-box info" style="margin-top:.5rem;">{msg}</div>'
        )

    from collections import defaultdict

    by_species: Dict[str, List[Dict]] = defaultdict(list)
    for m in perf_metrics:
        by_species[m.get("species", "unknown")].append(m)

    species_list = list(by_species.keys())
    # D3 "category10" — one stable colour per species, overlaid in both panels.
    _SP_PAL = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
               "#8c564b", "#e377c2", "#17becf", "#bcbd22", "#7f7f7f"]
    sp_colour = {sp: _SP_PAL[i % len(_SP_PAL)]
                 for i, sp in enumerate(species_list)}

    # (dict key, axis label, interpretation hint)
    metric_meta = [
        ("KGE", "KGE", "higher is better"),
        ("NSE", "NSE", "higher is better"),
        ("R2", "R²", "higher is better"),
        ("RMSE", "RMSE", "lower is better"),
        ("PBIAS", "PBIAS (%)", "0 is best"),
    ]

    def _clean(vals):
        out = []
        for v in vals:
            if isinstance(v, (int, float)) and not (
                    isinstance(v, float) and (math.isnan(v) or math.isinf(v))):
                out.append(float(v))
        return out

    def _median(xs):
        s = sorted(xs)
        k = len(s)
        return s[k // 2] if k % 2 else 0.5 * (s[k // 2 - 1] + s[k // 2])

    # Shared chart furniture (matches the obs/sim time-series charts above).
    _leg = {"orientation": "h", "yanchor": "bottom", "y": 1.02,
            "xanchor": "left", "x": 0}
    _mgn = {"l": 60, "r": 20, "t": 86, "b": 52}
    figs_html = []
    for key, label, hint in metric_meta:
        # One value per reach/HRU for this species; need >= 2 to draw.
        series = []
        for sp in species_list:
            vals = _clean([mm.get(key) for mm in by_species[sp]])
            if len(vals) >= 2:
                series.append((sp, sorted(vals)))
        if not series:
            continue

        _multi_sp = len(series) > 1
        hist_traces, cdf_traces, med_shapes = [], [], []
        for sp, arr in series:
            c = sp_colour.get(sp, "#1f77b4")
            n = len(arr)
            nbins = int(min(12, max(4, math.ceil(math.sqrt(n)))))
            _nm = f"{sp} (n={n})"
            # Left panel: overlaid histogram (non-cumulative).
            hist_traces.append({
                "type": "histogram", "x": arr, "name": _nm,
                "legendgroup": sp, "opacity": 0.55, "nbinsx": nbins,
                "marker": {"color": c, "line": {"color": c, "width": 1.3}},
                "hovertemplate": (f"{label} %{{x}}<br>%{{y}} "
                                  f"{feat_plural}<extra>{_nm}</extra>")})
            # Median marker as a dashed vertical line spanning the plot.
            med = _median(arr)
            med_shapes.append({
                "type": "line", "x0": med, "x1": med, "xref": "x",
                "yref": "paper", "y0": 0, "y1": 1,
                "line": {"color": c, "dash": "dash", "width": 1.2}})
            # Right panel: empirical CDF (step, post) + markers.
            ecdf_y = [(i + 1) / float(n) for i in range(n)]
            cdf_traces.append({
                "type": "scatter", "mode": "lines+markers",
                "x": arr, "y": ecdf_y, "name": _nm, "legendgroup": sp,
                "line": {"color": c, "width": 2, "shape": "hv"},
                "marker": {"color": c, "size": 5},
                "hovertemplate": (f"{label} %{{x}}<br>fraction ≤ %{{y:.2f}}"
                                  f"<extra>{_nm}</extra>")})

        hist_layout = {
            "title": {"text": f"{label} — non-cumulative (histogram)"},
            "xaxis": {"title": f"{label}  ({hint})"},
            "yaxis": {"title": f"number of {feat_plural}"},
            "barmode": "overlay", "shapes": med_shapes,
            "showlegend": _multi_sp, "legend": _leg,
            "hovermode": "closest", "margin": _mgn}
        cdf_layout = {
            "title": {"text": f"{label} — cumulative (CDF)"},
            "xaxis": {"title": f"{label}  ({hint})"},
            "yaxis": {"title": f"fraction of {feat_plural} ≤ value",
                      "range": [0, 1.02]},
            "showlegend": _multi_sp, "legend": _leg,
            "hovermode": "closest", "margin": _mgn}
        cL = _plotly_chart(_plot_id("dist-h"), hist_traces, hist_layout,
                           height=380, card=False)
        cR = _plotly_chart(_plot_id("dist-c"), cdf_traces, cdf_layout,
                           height=380, card=False)
        figs_html.append(
            '<div class="card" style="margin-top:1rem;">'
            '<div style="display:flex;flex-wrap:wrap;gap:1rem;">'
            f'<div style="flex:1 1 420px;min-width:300px;">{cL}</div>'
            f'<div style="flex:1 1 420px;min-width:300px;">{cR}</div>'
            '</div></div>')

    if not figs_html:
        return heading + (
            '<div class="highlight-box info" style="margin-top:.5rem;">'
            f'Not enough valid metric values across {feat_plural} to draw '
            'distributions.</div>'
        )

    intro = (
        '<p style="color:var(--text2);margin:.4rem 0 .2rem;">'
        f'Each metric is summarised across the <strong>{n_stations}</strong> '
        f'{feat_plural} with observations, two complementary ways: a '
        f'<strong>non-cumulative</strong> histogram (how many {feat_plural} fall '
        'in each value band; dashed line = median) and a '
        f'<strong>cumulative</strong> distribution (the fraction of {feat_plural} '
        'at or below a given value). A tight distribution means performance is '
        f'consistent across the network; a wide one flags {feat_plural} that lag '
        'behind.</p>'
    )
    return heading + intro + "\n".join(figs_html)


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
    # The model-config setup writes the *clipped* GRQA station/observation CSVs
    # into <run_dir>/openwq_in/grqa_clipped_data/ — NOT into the raw GRQA
    # database dir (grqa_local_data_path).  Resolve that folder so the map can
    # find station coordinates.  (This was why the map came up empty for GRQA
    # runs: it only looked inside the raw DB dir, where the clipped CSV never
    # exists.)
    def _resolve_clipped_dir():
        cands = []
        _rd = model_config.get("dir2save_input_files")
        if not _rd:
            try:
                _exe = (_ci.get_container_config(model_config) or {}).get(
                    "executable_path", "")
            except Exception:
                _exe = ""
            _rd = os.path.dirname(os.path.abspath(_exe)) if _exe else ""
        if _rd:
            cands.append(os.path.join(_rd, "openwq_in", "grqa_clipped_data"))
        if grqa_dir:
            cands.append(os.path.join(grqa_dir, "grqa_clipped_data"))
            cands.append(grqa_dir)
        for _c in cands:
            if _c and os.path.isdir(_c):
                return _c
        return None
    _clipped_dir = _resolve_clipped_dir()
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
    if not station_locations and _clipped_dir:
        stn_csv = os.path.join(_clipped_dir, "grqa_clipped_stations.csv")
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
    elif _clipped_dir:
        # GRQA case: count rows per station from the clipped observations CSV
        # (keyed by site_id) so the marker radius reflects data density.
        _obs_clip = os.path.join(_clipped_dir, "grqa_clipped_observations.csv")
        if os.path.isfile(_obs_clip):
            try:
                _df = pd.read_csv(
                    _obs_clip, usecols=lambda c: c == "site_id")
                for sid, n in _df.groupby("site_id").size().to_dict().items():
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

    # 2b) Colour stations directly from matched_data (decoupled from the
    # performance-metrics table, which may be unavailable).  Also handles the
    # common case where match_stations() returns 0 matches — e.g. a SUMMA
    # domain with only a lumped river network and no HRU polygons: then every
    # marker would otherwise be gray.  We compute the objective metric from
    # the obs/sim pairs and:
    #   • per matched reach → colour that reach's stations
    #   • if NOTHING matched but obs/sim pairs exist (lumped) → colour every
    #     station by the single overall calibrated metric and treat them as
    #     primary so the colormap (not gray) is used.
    _obj_name = str((calibration_settings or {}).get("objective_function")
                    or "KGE").upper()
    _lumped_colour = False
    try:
        import pandas as _pd
        import numpy as _np
        from .objective_functions import ObjectiveFunction as _OF

        def _calc_metric(o, s):
            o = _np.asarray(o, dtype=float); s = _np.asarray(s, dtype=float)
            mask = ~(_np.isnan(o) | _np.isnan(s))
            o, s = o[mask], s[mask]
            if o.size < 2:
                return None
            if _obj_name == "NSE":
                return float(_OF.nse(o, s))
            if _obj_name == "RMSE":
                return float(_OF.rmse(o, s))
            if _obj_name == "PBIAS":
                return float(_OF.pbias(o, s))
            return float(_OF.kge(o, s))

        if isinstance(matched_data, _pd.DataFrame) and not matched_data.empty \
                and {'observed', 'simulated'} <= set(matched_data.columns):
            _by_reach: Dict[str, float] = {}
            if 'reach_id' in matched_data.columns:
                for _rid, _g in matched_data.groupby('reach_id'):
                    _mv = _calc_metric(_g['observed'], _g['simulated'])
                    if _mv is not None:
                        _by_reach[_norm_id(_rid)] = _mv
            _overall = _calc_metric(matched_data['observed'],
                                    matched_data['simulated'])
            # Fill matched stations (don't overwrite perf-metrics values).
            for _sid, _fid in s2f.items():
                if station_stats.get(_sid) and \
                        station_stats[_sid]["metrics"].get(_obj_name) is None:
                    _mv = _by_reach.get(_norm_id(_fid))
                    if _mv is not None:
                        station_stats[_sid]["metrics"][_obj_name] = _mv
            _have_any = any(st["metrics"].get(_obj_name) is not None
                            for st in station_stats.values())
            if (not _have_any) and (_overall is not None):
                primary = set(primary)
                _lumped_reach = (next(iter(_by_reach)) if _by_reach else
                                 (str(matched_data['reach_id'].iloc[0])
                                  if 'reach_id' in matched_data.columns else ""))
                for _sid in station_stats:
                    station_stats[_sid]["metrics"][_obj_name] = _overall
                    primary.add(_sid)
                    if _sid not in s2f and _lumped_reach:
                        s2f[_sid] = _lumped_reach
                _lumped_colour = True
    except Exception as _e:
        logger.debug(f"obs-map matched_data colouring skipped: {_e}")

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
        # KGE / NSE / logNSE are bounded ABOVE at 1 (perfect) but unbounded
        # below, so we visualise the 0–1 SKILL scale (0 = no better than the
        # mean of obs) and clamp negative values to the "worst" colour.
        "KGE":    {"vmin": 0.0,  "vmax": 1.0, "reversed": False, "use_abs": False, "kind": "skill"},
        "NSE":    {"vmin": 0.0,  "vmax": 1.0, "reversed": False, "use_abs": False, "kind": "skill"},
        "LOGNSE": {"vmin": 0.0,  "vmax": 1.0, "reversed": False, "use_abs": False, "kind": "skill"},
        # Genuinely bounded by definition.
        "R2":     {"vmin": 0.0,  "vmax": 1.0, "reversed": False, "use_abs": False, "kind": "theoretical"},
        "R":      {"vmin": -1.0, "vmax": 1.0, "reversed": False, "use_abs": False, "kind": "theoretical"},
        # No theoretical max — 0 is best; coloured by Moriasi et al. (2007)
        # performance cutoffs (|PBIAS| ≤ 25 % = satisfactory).  Uses |.|.
        "PBIAS":  {"vmin": 0.0,  "vmax": 25.0, "reversed": True, "use_abs": True,  "kind": "performance"},
        "MAPE":   {"vmin": 0.0,  "vmax": 50.0, "reversed": True, "use_abs": False, "kind": "performance"},
        # Unbounded and units-dependent → derive the scale from the data.
        "RMSE":   None,
        "MAE":    None,
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
            "range_kind": _preset.get("kind", "theoretical"),
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
    if _lumped_colour:
        sub += (f" &middot; stations coloured by the overall calibrated "
                f"{_obj_name} (lumped network — no per-feature spatial match)")
    html = ['<div class="section" id="obs-map">']
    html.append('<h2>Observation Map</h2>')
    html.append(f'<p class="muted" style="margin-top:-8px">{sub}</p>')
    # Map export button (uses html2pdf on the map div; crossOrigin tiles let
    # html2canvas capture the basemap).  Mirrors the per-figure PDF buttons.
    html.append(
        '<div style="display:flex;justify-content:flex-end;margin:0 0 6px;">'
        '<button type="button" class="no-pdf" onclick="_owqExportMap()" '
        'title="Download the map as PDF" '
        'style="display:inline-flex;align-items:center;gap:5px;padding:5px 11px;'
        'font-size:11px;font-weight:600;background:#0066cc;color:#fff;border:none;'
        'border-radius:7px;cursor:pointer;box-shadow:0 2px 6px rgba(0,0,0,.25);">'
        '⬇ Download map (PDF)</button></div>'
        '<script>function _owqExportMap(){'
        'var el=document.getElementById("obs-map-leaflet"); if(!el) return;'
        'if(typeof html2pdf==="undefined"){ window.print(); return; }'
        'var r=el.getBoundingClientRect(), mm=25.4/96;'
        'var wmm=r.width*mm, hmm=r.height*mm;'
        'html2pdf().set({margin:5, filename:"observation_map.pdf",'
        'image:{type:"jpeg",quality:0.95},'
        'html2canvas:{useCORS:true, scale:2, backgroundColor:"#ffffff", logging:false},'
        'jsPDF:{unit:"mm", orientation:"portrait", format:[wmm+10,hmm+10]}'
        '}).from(el).save().catch(function(e){console.error(e);'
        'alert("Map export failed \\u2014 the basemap tiles may block capture. '
        'Switch to the OpenStreetMap basemap and retry, or take a screenshot.");});'
        '}</script>')
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
        '  {maxZoom:19, crossOrigin:true, attribution:"&copy; OpenStreetMap contributors"});\n'
        'var sat = L.tileLayer(\n'
        '  "https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}",\n'
        '  {maxZoom:19, crossOrigin:true, attribution:"Tiles &copy; Esri \\u2014 Source: Esri, Maxar, Earthstar Geographics, USDA, USGS"});\n'
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
        '    } else if(CB.range_kind === "skill"){\n'
        '      kindHint = "<span style=\\"font-size:10px;color:#666;margin-left:6px\\">'
        '(skill scale \\u00b7 1 = perfect, 0 = mean)</span>";\n'
        '    } else if(CB.range_kind === "performance"){\n'
        '      kindHint = "<span style=\\"font-size:10px;color:#666;margin-left:6px\\">'
        '(0 = best \\u00b7 performance cutoff)</span>";\n'
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
/* Per-figure "Download PDF" buttons, pinned top-right of each plot card. */
.fig-pdf-btn {
    position:absolute; top:10px; right:10px; z-index:5;
    display:inline-flex; align-items:center; gap:5px;
    padding:5px 10px; font-size:11px; font-weight:600; line-height:1;
    background:rgba(0,102,204,.92); color:#fff;
    border:none; border-radius:7px; cursor:pointer;
    box-shadow:0 2px 6px rgba(0,0,0,.25);
    font-family:inherit; opacity:.85; transition:opacity .15s, background .15s;
}
.fig-pdf-btn:hover { opacity:1; background:#004499; }
@media print { .fig-pdf-btn { display:none !important; } }
"""
    return f"""
<style>{css}</style>
<script src="https://cdnjs.cloudflare.com/ajax/libs/html2pdf.js/0.10.1/html2pdf.bundle.min.js"></script>
<script src="https://cdnjs.cloudflare.com/ajax/libs/jspdf/2.5.1/jspdf.umd.min.js"></script>
<script>
// ── Per-figure "Download PDF" buttons ──────────────────────────────────
// Every matplotlib figure is a base64 PNG data-URI (same-origin → no canvas
// taint), so a single-image jsPDF export is fast and reliable.  We inject a
// button into each plot card on load.  jsPDF ships inside the html2pdf
// bundle (window.jspdf); if that didn't load we fall back to a PNG download
// so the button always does something.
function _owqFigName(img) {{
    var n = (img.getAttribute('alt') || 'figure');
    n = n.replace(/[^A-Za-z0-9_\\- ]+/g, '').trim().replace(/\\s+/g, '_');
    return n || 'figure';
}}
function _owqFigPdfImg(img) {{
    if(!img) return;
    var name = _owqFigName(img);
    if(typeof html2pdf === 'undefined') {{   // CDN missing — download the PNG
        var a = document.createElement('a');
        a.href = img.src; a.download = name + '.png'; a.click();
        return;
    }}
    // Build a clean WHITE export page with a header, then let html2pdf size
    // and paginate it.  Scaling the figure to the page width guarantees the
    // whole graph fits (no more cut-off), regardless of the figure's pixel
    // dimensions — the same idea as the 1_config / 2_read_output exports.
    var title = (img.getAttribute('alt') || 'Figure');
    var wrap = document.createElement('div');
    wrap.style.cssText = 'position:absolute;left:-10000px;top:0;width:1000px;'+
        'background:#ffffff;padding:26px;box-sizing:border-box;'+
        'font-family:-apple-system,Segoe UI,Helvetica,Arial,sans-serif;';
    var head = document.createElement('div');
    head.style.cssText = 'border-bottom:2px solid #0066cc;padding-bottom:10px;margin-bottom:16px;';
    var h1 = document.createElement('div');
    h1.style.cssText = 'font-size:20px;font-weight:700;color:#111;';
    h1.textContent = title;
    var sub = document.createElement('div');
    sub.style.cssText = 'font-size:12px;color:#555;margin-top:3px;';
    sub.textContent = 'OpenWQ \\u2014 Calibration Results  \\u00b7  ' + new Date().toLocaleString();
    head.appendChild(h1); head.appendChild(sub);
    var im = document.createElement('img');
    im.style.cssText = 'display:block;width:100%;height:auto;';
    wrap.appendChild(head); wrap.appendChild(im);
    document.body.appendChild(wrap);
    var opt = {{
        margin:[10, 10, 12, 10], filename: name + '.pdf',
        image:{{ type:'jpeg', quality:0.95 }},
        html2canvas:{{ scale:2, backgroundColor:'#ffffff', useCORS:true, logging:false }},
        jsPDF:{{ unit:'mm', format:'a4', orientation:'portrait' }},
        pagebreak:{{ mode:['css','legacy'] }}
    }};
    var _done = false;
    function _run() {{
        if(_done) return; _done = true;
        html2pdf().set(opt).from(wrap).save().then(function(){{ wrap.remove(); }})
          .catch(function(e){{ console.error('Figure PDF failed:', e); wrap.remove();
            var a2 = document.createElement('a');
            a2.href = img.src; a2.download = name + '.png'; a2.click(); }});
    }}
    // Wait for the (cloned) image to decode before capturing, otherwise
    // html2canvas grabs an empty box and the figure is missing from the PDF.
    im.onload = _run;
    im.onerror = function(){{ wrap.remove();
        var a3 = document.createElement('a');
        a3.href = img.src; a3.download = name + '.png'; a3.click(); }};
    im.src = img.src;
    if(im.complete && im.naturalWidth) _run();   // already cached
}}
function _owqInjectFigButtons() {{
    document.querySelectorAll('img.plot-img').forEach(function(img) {{
        // Pin the button to the figure's container — a plot card, or (for the
        // time-series figures that aren't wrapped in a card) the image's
        // immediate parent.
        var host = img.closest('.card') || img.parentElement;
        if(!host || host.querySelector(':scope > .fig-pdf-btn')) return;
        if(getComputedStyle(host).position === 'static')
            host.style.position = 'relative';
        var b = document.createElement('button');
        b.type = 'button';
        b.className = 'fig-pdf-btn no-pdf';
        b.title = 'Download this figure as PDF';
        b.innerHTML = '\\u2B07 PDF';
        b.addEventListener('click', function() {{ _owqFigPdfImg(img); }});
        host.appendChild(b);
    }});
}}
if(document.readyState === 'loading')
    document.addEventListener('DOMContentLoaded', _owqInjectFigButtons);
else
    _owqInjectFigButtons();
</script>
"""
