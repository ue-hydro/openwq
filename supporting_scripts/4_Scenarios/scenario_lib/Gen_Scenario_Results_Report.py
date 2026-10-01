# Copyright 2026, Diogo Costa, diogo.costa@uevora.pt
# This file is part of OpenWQ model.
#
# This program, openWQ, is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""Scenario comparison report (HTML, Plotly) — written by ``scenarios.run_scenarios``.

Sections: run summary · scenario recipes (levers + parameters + what was
actually applied) · change versus the baseline per species (basin aggregate +
per unit) · threshold exceedance · time-series overlays per species for the
focus units · reproducibility."""
from __future__ import annotations

import os
import re
import json
import html as html_lib
import logging
from pathlib import Path
from typing import Any, Dict, List, Optional

logger = logging.getLogger(__name__)

try:
    from calibration_lib import report_helpers as rh
    from calibration_lib import Gen_Calibration_Results_Report as _GRR
    from . import scenarios as _S
except ImportError:                          # pragma: no cover
    import report_helpers as rh
    import Gen_Calibration_Results_Report as _GRR
    import scenarios as _S

_PAL = ["#2563eb", "#dc2626", "#059669", "#d97706", "#7c3aed", "#0891b2", "#be185d",
        "#65a30d", "#ea580c", "#4f46e5", "#0d9488", "#b91c1c"]


def _fmt(v, nd=3):
    try:
        if v is None or v != v:
            return "&mdash;"
        v = float(v)
        if abs(v) >= 100:
            return f"{v:,.0f}"
        return f"{v:.{nd}g}"
    except (TypeError, ValueError):
        return html_lib.escape(str(v))


def _pct(v):
    try:
        if v is None or v != v:
            return "&mdash;"
        v = float(v)
        col = "#dc2626" if v > 0.5 else ("#059669" if v < -0.5 else "var(--text2)")
        return f'<span style="color:{col};font-weight:600;">{v:+.1f}%</span>'
    except (TypeError, ValueError):
        return "&mdash;"


def _focus_units(work_dir: Path, summary, max_units: int = 6) -> List[str]:
    """Units to chart: the calibration's observation units first, then the
    units with the largest baseline mean."""
    units: List[str] = []
    obs = work_dir / "calibration_observations.csv"
    if obs.is_file():
        try:
            import pandas as pd
            o = pd.read_csv(obs, usecols=lambda c: c in ("reach_id", "is_primary"))
            if "reach_id" in o.columns:
                prim = o[o["is_primary"].fillna(True).astype(bool)] if "is_primary" in o.columns else o
                for u in prim["reach_id"].astype(str).unique().tolist():
                    units.append(_S._norm_id(u))
        except Exception:
            pass
    if summary is not None and not summary.empty:
        b = summary[(summary["unit"] != "ALL")]
        base = b[b["scenario"] == "baseline"] if "baseline" in set(b["scenario"]) else b
        for u in base.sort_values("mean", ascending=False)["unit"].astype(str).tolist():
            if u not in units:
                units.append(u)
    return units[:max_units]


def _fmt_units(u, max_show: int = 8) -> str:
    """'all units' | '3 units: a, b, c' | '42 units: a, b, ... (+34 more)'."""
    if u is None or (isinstance(u, str) and u.strip().lower() in ("", "all")):
        return "all units"
    ids = [str(x).strip() for x in (u if isinstance(u, (list, tuple, set)) else str(u).replace(";", ",").split(","))
           if str(x).strip()]
    if not ids or any(x.lower() == "all" for x in ids):
        return "all units"
    shown = ", ".join(ids[:max_show])
    more = f" (+{len(ids) - max_show} more)" if len(ids) > max_show else ""
    return f"{len(ids)} unit{'' if len(ids) == 1 else 's'}: {shown}{more}"


def generate_scenario_report(*, model_config: Dict[str, Any], results: Dict[str, Any],
                             scenarios: List[Dict[str, Any]], work_dir: Optional[str] = None,
                             calibration_work_dir: Optional[str] = None,
                             report_stem: str = "scenarios") -> Optional[str]:
    import pandas as pd
    work_dir = Path(work_dir or calibration_work_dir)
    scen_root = work_dir / "scenarios"
    report_path = work_dir / f"{report_stem}_scenarios_report.html"
    summary = results.get("summary")
    if summary is None:
        p = scen_root / "scenarios_summary.csv"
        summary = pd.read_csv(p) if p.is_file() else pd.DataFrame()
    order = [n for n in results.get("order", []) if n in results.get("scenarios", {})]
    ok_names = [n for n in order if results["scenarios"][n].get("success")]
    species = results.get("species") or []
    thresholds = results.get("thresholds") or {}
    period = results.get("period")

    H: List[str] = []
    H.append(rh.build_html_head("OpenWQ Scenario Analysis", extra_css=_GRR._get_results_css() + _extra_css()))
    H.append("<body>")
    H.append(_GRR._plotly_bootstrap(model_config.get("project_name", "OpenWQ Scenarios")))
    H.append('<div class="layout">')
    H.append(rh.build_sidebar([
        {"id": "summary", "label": "Summary"}, {"id": "recipes", "label": "Scenario recipes"},
        {"id": "change", "label": "Change vs baseline"}, {"id": "exceed", "label": "Threshold exceedance"},
        {"id": "series", "label": "Time series"}, {"id": "repro", "label": "Reproducibility"}]))
    H.append('<div class="main">')
    from datetime import datetime as _dt
    H.append(rh.build_header(
        'Open<span style="color:var(--secondary);">WQ</span> &mdash; Scenario Analysis',
        html_lib.escape(str(model_config.get("project_name", ""))),
        [f"Generated {_dt.now():%Y-%m-%d %H:%M}", f"{len(order)} scenario run(s)"],
        badge_text="SCENARIOS", badge_class="badge-accent"))
    H.append('<div class="container">')

    # ── Summary ──
    n_fail = len(order) - len(ok_names)
    H.append('<div class="section" id="summary"><h2>Summary</h2>')
    H.append(rh.build_kpi_grid([
        {"icon": "&#128203;", "value": str(len(order)), "label": "scenarios run (incl. baseline)"},
        {"icon": "&#10003;", "value": str(len(ok_names)), "label": "completed"},
        {"icon": "&#9888;", "value": str(n_fail), "label": "failed"},
        {"icon": "&#128197;", "value": (f"{period[0][:10]} &rarr; {period[1][:10]}" if period else "model period"),
         "label": "scored period"},
        {"icon": "&#9881;", "value": ("calibrated" if results.get("calibrated_setup") else "as configured"),
         "label": "base model"},
    ]))
    prov = results.get("calibrated_setup") or {}
    _base_txt = (f"Base: the calibrated model config (best evaluation "
                 f"<code>{html_lib.escape(str(prov.get('best_eval') or '?'))}</code>"
                 + (f", {html_lib.escape(str(prov.get('objective_function') or 'objective'))} = "
                    f"{_fmt(prov.get('best_objective'))}" if prov.get('best_objective') is not None else "")
                 + f"; inputs in <code>{html_lib.escape(str(results.get('base_dir') or ''))}</code>)."
                 if prov else
                 f"Base: the model config as is (<code>{html_lib.escape(str(results.get('base_dir') or ''))}</code>; "
                 "no calibrated setup found).")
    H.append(rh.build_highlight_box(
        _base_txt + " Every scenario is that SAME model with its input files changed by the "
        "levers listed below (loads, point sources, initial stores, in-stream rates, "
        "host forcing). <b>baseline</b> = the base model with no lever; it is drawn on every "
        "graph next to the scenarios, and every change is reported relative to it. "
        "Concentrations are the simulated values in the "
        f"<code>{html_lib.escape(', '.join(results.get('compartments') or []))}</code> "
        "compartment(s) over the scored period.", "info"))
    if n_fail:
        for n in order:
            e = results["scenarios"][n]
            if not e.get("success"):
                H.append(rh.build_error_card(f"Scenario '{html_lib.escape(n)}' failed",
                                             html_lib.escape(str(e.get("error", ""))[:800])))
    H.append("</div>")

    # ── Recipes ──
    H.append('<div class="section" id="recipes"><h2>Scenario recipes</h2>')
    H.append('<p class="hint">What each scenario changes, and what the runner actually applied '
             '(an option that could not act on this run says so).</p>')
    for n in order:
        e = results["scenarios"][n]
        sc = next((s for s in scenarios if s.get("name") == n), {})
        applied = e.get("applied") or {}
        _open = "open" if n != "baseline" else ""
        _fail = "" if e.get("success") else ' <span style="color:#dc2626;font-weight:700;">FAILED</span>'
        H.append(f'<details class="module-group" {_open}>'
                 f'<summary>{html_lib.escape(n)} '
                 f'<span class="group-badge">{len(_S.scenario_levers(sc))} option(s)</span>'
                 f'{_fail}</summary><div class="module-content">')
        if sc.get("description"):
            H.append(f'<p style="margin-top:0;">{html_lib.escape(str(sc["description"]))}</p>')
        _zones = [z for z in (sc.get("management") or []) if isinstance(z, dict) and z.get("options")]
        if _zones:
            _zt = []
            for z in _zones:
                _zu = _fmt_units(z.get("where"))
                _zn = "All units" if _zu == "all units" else (
                    (f"Group '{z.get('name')}'" if z.get("name") else "Group") + f" ({_zu})")
                _zt.append(f"<b>{html_lib.escape(_zn)}</b>: {len(z.get('options') or [])} option(s)")
            H.append('<p style="margin:.1rem 0 .4rem;font-size:.84rem;">Management &amp; policy options by '
                     'where they apply &mdash; ' + " &middot; ".join(_zt) + '. In a group, an option replaces '
                     'the same option set for all units; climate options apply to the whole domain.</p>')
        levers = applied.get("levers") or _S.scenario_levers(sc)
        if not levers:
            H.append('<p class="hint">No option &mdash; the base model as is.</p>')
        else:
            H.append('<table class="data-table" style="table-layout:fixed;"><thead><tr>'
                     '<th style="width:9%;">Group</th><th style="width:22%;">Option</th>'
                     '<th style="width:28%;">Parameters</th><th style="width:16%;">Where</th>'
                     '<th>Applied</th></tr></thead><tbody>')
            for lv in levers:
                lev = _S.LEVER_BY_ID.get(lv.get("id"), {})
                _prm = lv.get("params") or {}
                prm = ", ".join(f"{html_lib.escape(str(k))}={html_lib.escape(str(v))}"
                                for k, v in _prm.items() if k not in ("units", "unit", "exclude_units"))
                _grp = "Climate" if _S.lever_group(lv.get("id")) == "climate" else "Management &amp; policy"
                _sp = _S.lever_spatial(lev) if lev else ""
                _where = lv.get("where") or (_fmt_units(_prm.get("units")) if _sp == "units"
                                             else (f"unit {_prm.get('unit')}" if _sp == "unit" else "whole domain"))
                res = lv.get("result") or lv.get("ts") or {}
                if lv.get("forcing"):
                    fo = lv["forcing"]
                    _files = fo.get("files") or (applied.get("forcing_override") or {}).get("files") or []
                    res = {"forcing": "perturbed copy in forcing_scenario/ (%d file%s)" % (
                        len(_files), "" if len(_files) == 1 else "s")}
                    if fo.get("note"):
                        res["note"] = fo["note"]
                err = lv.get("error")
                _app = (('<span style="color:#dc2626">' + html_lib.escape(str(err)) + '</span>') if err
                        else html_lib.escape(json.dumps(res, default=str)[:300]) if res else "&mdash;")
                H.append(f'<tr><td style="font-size:.78rem;">{_grp}</td>'
                         f'<td title="{html_lib.escape(lev.get("desc", ""))}"><b>'
                         f'{html_lib.escape(lev.get("label", lv.get("id", "?")))}</b>'
                         f'<br><span class="hint">{html_lib.escape(lev.get("ref", ""))[:120]}</span></td>'
                         f'<td style="font-size:.8rem;overflow-wrap:anywhere;">{prm}</td>'
                         f'<td style="font-size:.8rem;overflow-wrap:anywhere;">{html_lib.escape(str(_where))}</td>'
                         f'<td style="font-size:.8rem;overflow-wrap:anywhere;">{_app}</td></tr>')
            H.append("</tbody></table>")
        ssf = applied.get("ss_factor")
        if ssf:
            H.append(f'<p class="hint">Diffuse-load rows rescaled: {ssf.get("rows", 0)} in {ssf.get("files", 0)} file(s).</p>')
        for note in applied.get("notes") or []:
            H.append(f'<p style="font-size:.82rem;color:#d97706;margin:.15rem 0;">&#9888; {html_lib.escape(str(note))}</p>')
        H.append("</div></details>")
    H.append("</div>")

    # ── Change vs baseline ──
    H.append('<div class="section" id="change"><h2>Change versus the baseline</h2>')
    if summary is None or summary.empty:
        H.append('<p class="hint">No simulated series available.</p>')
    else:
        allrows = summary[summary["unit"] == "ALL"]
        # bar chart: % change of the mean per scenario × species (basin aggregate)
        traces = []
        for i, sp in enumerate(species):
            sub = allrows[allrows["species"] == sp]
            xs = [n for n in order if n != "baseline" and n in set(sub["scenario"])]
            ys = [float(sub[sub["scenario"] == n]["mean_pct_change"].iloc[0]) if "mean_pct_change" in sub.columns
                  and len(sub[sub["scenario"] == n]) else None for n in xs]
            traces.append({"type": "bar", "name": sp, "x": xs, "y": ys,
                           "marker": {"color": _PAL[i % len(_PAL)]},
                           "hovertemplate": "%{x}<br>" + sp + ": %{y:+.1f}%<extra></extra>"})
        # absolute means, baseline INCLUDED, one group per species
        abs_traces = []
        for i, sp in enumerate(species):
            sub = allrows[allrows["species"] == sp]
            xs = [n for n in order if n in set(sub["scenario"])]
            abs_traces.append({"type": "bar", "name": sp, "x": xs,
                               "y": [float(sub[sub["scenario"] == n]["mean"].iloc[0]) for n in xs],
                               "marker": {"color": _PAL[i % len(_PAL)]},
                               "hovertemplate": "%{x}<br>" + sp + ": %{y:.4g} mg/L<extra></extra>"})
        if abs_traces and any(t["x"] for t in abs_traces):
            H.append('<h3>Mean concentration, basin aggregate &mdash; baseline and scenarios</h3>')
            H.append(_GRR._plotly_chart("scn-bar-abs", abs_traces, {"barmode": "group",
                                                                     "yaxis": {"title": "mean (mg/L)"},
                                                                     "xaxis": {"title": ""}}, height=380))
        if traces and any(t["x"] for t in traces):
            H.append('<h3>Mean concentration, basin aggregate &mdash; % change vs baseline</h3>')
            H.append(_GRR._plotly_chart("scn-bar", traces, {"barmode": "group",
                                                             "yaxis": {"title": "% change of the mean"},
                                                             "xaxis": {"title": ""}}, height=380))
        # table per species: units × scenarios (mean, %Δ)
        for sp in species:
            sub = summary[summary["species"] == sp]
            if sub.empty:
                continue
            units = ["ALL"] + [u for u in sub["unit"].astype(str).unique().tolist() if u != "ALL"]
            H.append(f'<h3>{html_lib.escape(sp)} &mdash; mean over the period (and % change)</h3>')
            H.append('<div class="table-wrap"><table class="data-table sortable"><thead><tr><th>Unit</th>'
                     + "".join(f"<th>{html_lib.escape(n)}</th>" for n in order) + "</tr></thead><tbody>")
            for u in units[:60]:
                cells = []
                for n in order:
                    r = sub[(sub["unit"].astype(str) == u) & (sub["scenario"] == n)]
                    if r.empty:
                        cells.append("<td>&mdash;</td>")
                    else:
                        r = r.iloc[0]
                        pc = r.get("mean_pct_change") if "mean_pct_change" in sub.columns else None
                        cells.append(f'<td>{_fmt(r["mean"])}'
                                     + (f' <span style="font-size:.75rem;">({_pct(pc)})</span>' if n != "baseline" else "")
                                     + "</td>")
                lab = "<b>basin (all units)</b>" if u == "ALL" else html_lib.escape(u)
                H.append(f"<tr><td>{lab}</td>{''.join(cells)}</tr>")
            H.append("</tbody></table></div>")
            if len(units) > 60:
                H.append(f'<p class="hint">{len(units) - 60} more units in <code>scenarios/scenarios_summary.csv</code>.</p>')
    H.append("</div>")

    # ── Exceedance ──
    H.append('<div class="section" id="exceed"><h2>Threshold exceedance</h2>')
    if thresholds and summary is not None and not summary.empty and "exceed_frac" in summary.columns:
        H.append('<p class="hint">Share of the scored time steps above the threshold (basin aggregate and per unit).</p>')
        H.append('<div class="table-wrap"><table class="data-table"><thead><tr><th>Species</th><th>Threshold</th><th>Unit</th>'
                 + "".join(f"<th>{html_lib.escape(n)}</th>" for n in order) + "</tr></thead><tbody>")
        for sp, t in thresholds.items():
            sub = summary[(summary["species"].map(_S._canon) == _S._canon(sp))]
            if sub.empty:
                continue
            for u in ["ALL"] + [x for x in sub["unit"].astype(str).unique().tolist() if x != "ALL"][:30]:
                cells = []
                for n in order:
                    r = sub[(sub["unit"].astype(str) == u) & (sub["scenario"] == n)]
                    cells.append("<td>&mdash;</td>" if r.empty or r.iloc[0]["exceed_frac"] != r.iloc[0]["exceed_frac"]
                                 else f'<td>{100 * float(r.iloc[0]["exceed_frac"]):.1f}%</td>')
                H.append(f'<tr><td>{html_lib.escape(sp)}</td><td>{_fmt(t)}</td>'
                         f'<td>{"<b>basin</b>" if u == "ALL" else html_lib.escape(u)}</td>{"".join(cells)}</tr>')
        H.append("</tbody></table></div>")
    else:
        H.append('<p class="hint">No thresholds were set (Scenarios tab &rarr; thresholds), so nothing to report here.</p>')
    H.append("</div>")

    # ── Time series ──
    H.append('<div class="section" id="series"><h2>Time series</h2>')
    big_p = scen_root / "scenarios_simulated.csv"
    if big_p.is_file():
        big = pd.read_csv(big_p, parse_dates=["datetime"])
        focus = _focus_units(work_dir, summary)
        if not focus and not big.empty:
            focus = big["unit"].astype(str).unique().tolist()[:4]
        H.append(f'<p class="hint">Units shown: {", ".join(html_lib.escape(u) for u in focus)} '
                 '(observation units first, then the largest baseline means). Daily means. '
                 'Every unit is in <code>scenarios/scenarios_simulated.csv</code>.</p>')
        k = 0
        for sp in species:
            for u in focus:
                sub = big[(big["species"] == sp) & (big["unit"].astype(str) == str(u))]
                if sub.empty:
                    continue
                sub = sub.set_index("datetime").groupby("scenario").resample("D")["value"].mean().reset_index()
                traces = []
                for i, n in enumerate(order):
                    g = sub[sub["scenario"] == n]
                    if g.empty:
                        continue
                    traces.append({"type": "scatter", "mode": "lines", "name": n,
                                   "x": [d.strftime("%Y-%m-%d") for d in g["datetime"]],
                                   "y": [None if v != v else float(v) for v in g["value"]],
                                   # baseline: a mid grey that reads on the light AND the dark theme
                                   "line": {"width": 3.0 if n == "baseline" else 1.6,
                                            "color": "#6b7280" if n == "baseline" else _PAL[i % len(_PAL)],
                                            "dash": "solid"}})
                H.append(f'<h3>{html_lib.escape(sp)} &mdash; unit {html_lib.escape(str(u))}</h3>')
                H.append(_GRR._plotly_chart(f"scn-ts-{k}", traces,
                                            {"yaxis": {"title": f"{sp} (mg/L)"}, "xaxis": {"title": ""},
                                             "legend": {"orientation": "h"}}, height=380, log_axes="y"))
                k += 1
    else:
        H.append('<p class="hint">No simulated series file found.</p>')
    H.append("</div>")

    # ── Reproducibility ──
    H.append('<div class="section" id="repro"><h2>Reproducibility</h2>')
    H.append('<table class="data-table"><thead><tr><th>Scenario</th><th>Run folder</th><th>Runtime</th></tr></thead><tbody>')
    for n in order:
        e = results["scenarios"][n]
        H.append(f'<tr><td>{html_lib.escape(n)}</td><td><code style="font-size:.75rem;">'
                 f'{html_lib.escape(str(e.get("eval_dir") or e.get("dir", "")))}</code></td>'
                 f'<td>{_fmt(e.get("runtime_s"), 3)} s</td></tr>')
    H.append("</tbody></table>")
    H.append('<p class="hint">Each run folder holds the modified <code>openwq_in/</code>, '
             '<code>scenario_applied.json</code> (every lever and what it changed), '
             '<code>simulated.csv</code> and, for climate levers, <code>forcing_scenario/</code>.</p>')
    H.append("</div>")

    H.append("</div>")   # container
    H.append(rh.build_footer())
    H.append("</div>")   # main
    H.append("</div>")   # layout
    H.append(rh.build_theme_toggle_js())
    H.append("</body></html>")
    report_path.write_text("\n".join(H))
    logger.info(f"Scenario report saved to: {report_path}")
    return str(report_path)


def _extra_css() -> str:
    return """
.table-wrap{overflow-x:auto;}
.data-table{width:100%;border-collapse:collapse;font-size:.84rem;margin:.4rem 0 1rem;}
.data-table th{text-align:left;padding:.4rem .5rem;border-bottom:2px solid var(--border);white-space:nowrap;}
.data-table td{padding:.35rem .5rem;border-bottom:1px solid var(--border);vertical-align:top;}
"""
