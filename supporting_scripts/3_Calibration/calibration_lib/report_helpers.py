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
Report Helpers
===============

Shared CSS styles and HTML scaffolding for calibration reports.
Uses the same visual style as the model config report (Gen_Report.py).
"""

import json
import base64
import math
import os
import re
from datetime import datetime
from typing import List, Dict, Any, Optional
import html as html_lib


def _js(obj):
    """Convert a Python object to a safe JSON string for HTML embedding."""
    try:
        import numpy as np
        if isinstance(obj, (np.integer,)):
            obj = int(obj)
        elif isinstance(obj, (np.floating,)):
            obj = float(obj)
        elif isinstance(obj, np.ndarray):
            obj = obj.tolist()
    except ImportError:
        pass
    return json.dumps(obj, default=str)


def get_css_styles():
    """Return the CSS styles string matching Gen_Report.py visual style."""
    return """
    /* ── Reset & Typography ─────────────────────────────────── */
    *, *::before, *::after { box-sizing: border-box; margin: 0; padding: 0; }
    html { scroll-behavior: smooth; }
    body {
        font-family: 'Inter', -apple-system, BlinkMacSystemFont, sans-serif;
        background: var(--bg); color: var(--text);
        line-height: 1.65; font-size: 15px;
    }

    /* ── Theme Variables ─────────────────────────────────────── */
    :root, [data-theme="light"] {
        --primary: #0066cc; --primary-dark: #004499;
        --secondary: #00a86b; --accent: #ff6b35;
        --dark: #1a1a2e; --light: #f8f9fa;
        --bg: #f0f2f5; --surface: #fff;
        --text: #2d3436; --text2: #636e72; --text3: #a0aec0;
        --border: #e2e8f0; --border2: #d0d4da;
        --shadow: 0 4px 20px rgba(0,0,0,.06);
        --shadow-lg: 0 10px 40px rgba(0,0,0,.08);
        --plot-bg: #fff; --plot-grid: #f0f0f0; --plot-font: #333;
        --glass: rgba(255,255,255,.85);
    }
    [data-theme="dark"] {
        --primary: #4d9ee8; --primary-dark: #3a8bd4;
        --secondary: #34d399; --accent: #fb923c;
        --bg: #0f1117; --surface: #1a1b2e;
        --text: #e2e8f0; --text2: #a0aec0; --text3: #636e72;
        --border: #2a2b3d; --border2: #3a3b4d;
        --shadow: 0 4px 20px rgba(0,0,0,.3);
        --shadow-lg: 0 10px 40px rgba(0,0,0,.4);
        --plot-bg: #1e1f2e; --plot-grid: #2d2e42; --plot-font: #ccc;
        --glass: rgba(15,17,23,.92);
    }

    /* ── Layout ──────────────────────────────────────────────── */
    .layout { display: flex; min-height: 100vh; }
    .sidebar {
        position: sticky; top: 0; height: 100vh; width: 220px;
        background: var(--glass); backdrop-filter: blur(12px);
        border-right: 1px solid var(--border); padding: 1.2rem .8rem;
        display: flex; flex-direction: column; z-index: 100;
        overflow-y: auto;
    }
    .sidebar .logo { font-weight: 700; font-size: .95rem; color: var(--primary);
        margin-bottom: 1.5rem; padding: 0 .7rem; }
    .sidebar nav { display: flex; flex-direction: column; gap: 2px; }
    .sidebar nav a {
        display: block; padding: .4rem .7rem; font-size: .82rem;
        color: var(--text2); text-decoration: none; border-radius: 6px;
        transition: all .15s;
    }
    .sidebar nav a:hover { background: rgba(0,102,204,.06); }
    .sidebar nav a.active { background: rgba(0,102,204,.1); color: var(--primary); font-weight: 600; }

    .main { flex: 1; overflow-x: hidden; }

    .header {
        background: linear-gradient(135deg, var(--dark) 0%, #16213e 60%, var(--primary-dark) 100%);
        color: #fff; padding: 2.5rem 2rem 2rem;
    }
    .header h1 { font-size: clamp(1.4rem, 4vw, 2rem); font-weight: 700; margin-bottom: .3rem; }
    .header .subtitle { font-size: .95rem; opacity: .85; margin-bottom: 1rem; }
    .header .meta { display: flex; flex-wrap: wrap; gap: .6rem 1.5rem; font-size: .82rem; opacity: .7; }

    .container { max-width: 1100px; margin: 2rem auto; padding: 0 1.5rem; }
    .section { margin-bottom: 2.5rem; }
    .section h2 {
        font-size: 1.25rem; font-weight: 700; margin-bottom: 1.2rem;
        padding-bottom: .5rem; position: relative;
    }
    .section h2::before {
        content: ''; position: absolute; bottom: 0; left: 0;
        width: 60px; height: 3px; border-radius: 3px;
        background: linear-gradient(90deg, var(--primary), var(--secondary));
    }

    .footer {
        text-align: center; padding: 2rem; font-size: .8rem;
        color: var(--text3); border-top: 1px solid var(--border);
    }

    /* ── Cards ────────────────────────────────────────────────── */
    .card {
        background: var(--surface); border: 1px solid var(--border);
        border-radius: 12px; padding: 1.5rem; margin-bottom: 1rem;
        box-shadow: var(--shadow); transition: transform .2s, box-shadow .2s;
    }
    .card:hover { box-shadow: var(--shadow-lg); }
    .card.primary { border-left: 5px solid var(--primary); }
    .card.secondary { border-left: 5px solid var(--secondary); }
    .card.accent { border-left: 5px solid var(--accent); }
    .card h3 { font-size: 1rem; font-weight: 600; margin-bottom: .8rem; }

    /* ── KPI Grid ─────────────────────────────────────────────── */
    .kpi-grid {
        display: grid; grid-template-columns: repeat(auto-fit, minmax(160px, 1fr));
        gap: 1rem; margin-bottom: 1.5rem;
    }
    .kpi {
        background: var(--surface); border: 1px solid var(--border);
        border-radius: 12px; padding: 1.2rem .8rem; text-align: center;
        box-shadow: var(--shadow); transition: transform .2s;
        position: relative; overflow: hidden;
    }
    .kpi:hover { transform: translateY(-3px); box-shadow: var(--shadow-lg); }
    .kpi::before {
        content: ''; position: absolute; top: 0; left: 0; right: 0; height: 3px;
        background: linear-gradient(90deg, var(--primary), var(--secondary));
    }
    .kpi .icon { font-size: 1.6rem; margin-bottom: .4rem; }
    .kpi .value {
        font-size: clamp(.85rem, 3.5vw, 1.5rem); font-weight: 700;
        color: var(--primary); font-family: 'JetBrains Mono', monospace;
    }
    .kpi .label {
        font-size: .72rem; color: var(--text2); text-transform: uppercase;
        letter-spacing: .8px; margin-top: .2rem;
    }

    /* ── Tables ───────────────────────────────────────────────── */
    .table-wrap {
        overflow-x: auto; border-radius: 12px; border: 1px solid var(--border);
        margin-bottom: 1rem;
    }
    table.param-table {
        width: 100%; border-collapse: collapse; font-size: .85rem;
    }
    table.param-table th {
        background: var(--dark); color: #fff; padding: .6rem .8rem;
        text-align: left; font-weight: 600; position: sticky; top: 0;
    }
    table.param-table td {
        padding: .5rem .8rem; border-bottom: 1px solid var(--border);
    }
    table.param-table tr:nth-child(even) { background: rgba(0,0,0,.02); }
    [data-theme="dark"] table.param-table tr:nth-child(even) { background: rgba(255,255,255,.02); }
    table.param-table td.num {
        text-align: right; font-family: 'JetBrains Mono', monospace;
    }

    /* ── Badges ───────────────────────────────────────────────── */
    .badge {
        display: inline-flex; align-items: center; padding: .2rem .6rem;
        border-radius: 20px; font-size: .75rem; font-weight: 500;
    }
    .badge-primary { background: rgba(0,102,204,.1); color: var(--primary); }
    .badge-secondary { background: rgba(0,168,107,.1); color: var(--secondary); }
    .badge-accent { background: rgba(255,107,53,.1); color: var(--accent); }
    .badge-none { background: rgba(0,0,0,.04); color: var(--text3); }

    /* ── Highlight Boxes ──────────────────────────────────────── */
    .highlight-box {
        padding: 1.2rem; border-radius: 10px; margin-bottom: 1rem;
    }
    .highlight-box.info {
        background: linear-gradient(135deg, #eff6ff, #dbeafe);
        border-left: 5px solid var(--primary);
    }
    .highlight-box.success {
        background: linear-gradient(135deg, #ecfdf5, #d1fae5);
        border-left: 5px solid var(--secondary);
    }
    .highlight-box.warning {
        background: linear-gradient(135deg, #fff7ed, #ffedd5);
        border-left: 5px solid var(--accent);
    }
    [data-theme="dark"] .highlight-box.info { background: rgba(0,102,204,.08); }
    [data-theme="dark"] .highlight-box.success { background: rgba(0,168,107,.08); }
    [data-theme="dark"] .highlight-box.warning { background: rgba(255,107,53,.08); }

    /* ── Code Blocks ──────────────────────────────────────────── */
    pre.code-block {
        background: var(--dark); color: #e2e8f0; padding: 1.2rem;
        border-radius: 10px; overflow-x: auto; font-family: 'JetBrains Mono', monospace;
        font-size: .82rem; line-height: 1.6; margin: .8rem 0;
        position: relative;
    }
    pre.code-block .copy-btn {
        position: absolute; top: .5rem; right: .5rem; background: rgba(255,255,255,.1);
        border: 1px solid rgba(255,255,255,.2); color: #e2e8f0; padding: .3rem .6rem;
        border-radius: 6px; font-size: .72rem; cursor: pointer;
    }
    pre.code-block .copy-btn:hover { background: rgba(255,255,255,.2); }

    /* ── Collapsible Details ──────────────────────────────────── */
    details.module-details { margin-top: .8rem; }
    details.module-details > summary {
        padding: .6rem 1rem; font-weight: 600; cursor: pointer;
        background: var(--surface); border: 1px solid var(--border);
        border-radius: 8px; font-size: .9rem;
    }
    details.module-details > summary::before {
        content: '\\25B6'; margin-right: .5rem; display: inline-block;
        transition: transform .2s;
    }
    details.module-details[open] > summary::before { transform: rotate(90deg); }
    .module-content {
        padding: 1rem; border: 1px solid var(--border);
        border-top: none; border-radius: 0 0 8px 8px;
    }

    /* ── Theme Toggle ─────────────────────────────────────────── */
    .theme-toggle { margin-top: auto; padding-top: 1rem; }
    .theme-btn {
        display: flex; align-items: center; gap: .5rem;
        width: 100%; padding: .5rem .7rem; border: 1px solid var(--border);
        border-radius: 8px; background: transparent; color: var(--text2);
        cursor: pointer; font-size: .82rem;
    }
    .theme-btn:hover { border-color: var(--primary); }

    /* ── Form Elements ────────────────────────────────────────── */
    .form-group { margin-bottom: 1rem; }
    .form-group label {
        display: block; font-weight: 600; font-size: .85rem;
        margin-bottom: .3rem; color: var(--text);
    }
    .form-group .hint { font-size: .75rem; color: var(--text3); margin-top: .2rem; }
    .form-row { display: flex; gap: 1.5rem; flex-wrap: wrap; }
    .form-row .form-group { flex: 1; min-width: 180px; }

    select.form-select, input.form-input {
        width: 100%; padding: .45rem .7rem; font-size: .85rem;
        border: 1px solid var(--border); border-radius: 8px;
        background: var(--surface); color: var(--text);
        font-family: inherit; transition: border-color .15s;
    }
    select.form-select:focus, input.form-input:focus {
        outline: none; border-color: var(--primary);
        box-shadow: 0 0 0 3px rgba(0,102,204,.1);
    }
    input.form-input[type="number"] { font-family: 'JetBrains Mono', monospace; }

    /* Editable parameter table */
    table.param-table input.inline-input {
        width: 80px; padding: .2rem .4rem; font-size: .8rem;
        border: 1px solid var(--border); border-radius: 4px;
        background: var(--surface); color: var(--text);
        font-family: 'JetBrains Mono', monospace; text-align: right;
    }
    table.param-table input.inline-input:focus {
        outline: none; border-color: var(--primary);
    }
    table.param-table select.inline-select {
        padding: .2rem .3rem; font-size: .8rem;
        border: 1px solid var(--border); border-radius: 4px;
        background: var(--surface); color: var(--text);
    }
    table.param-table tr.param-disabled { opacity: .4; }

    /* Species checkboxes */
    .species-grid { display: flex; flex-wrap: wrap; gap: .8rem; margin: .5rem 0; }
    .species-item {
        display: flex; align-items: center; gap: .5rem;
        padding: .4rem .8rem; border: 1px solid var(--border); border-radius: 8px;
        background: var(--surface); transition: all .15s;
    }
    .species-item:has(input:checked) {
        border-color: var(--primary); background: rgba(0,102,204,.05);
    }
    .species-item input[type="checkbox"] { accent-color: var(--primary); }
    .species-item .weight-input {
        width: 55px; padding: .15rem .3rem; font-size: .8rem;
        border: 1px solid var(--border); border-radius: 4px;
        background: var(--surface); color: var(--text);
        font-family: 'JetBrains Mono', monospace; text-align: right;
    }

    /* Download button */
    .download-btn {
        display: inline-flex; align-items: center; gap: .5rem;
        padding: .7rem 1.5rem;
        background: linear-gradient(135deg, var(--primary), var(--secondary));
        color: #fff; border: none; border-radius: 10px; font-size: .95rem;
        font-weight: 600; cursor: pointer;
        transition: transform .15s, box-shadow .15s;
        margin-right: 1rem;
    }
    .download-btn:hover { transform: translateY(-2px); box-shadow: var(--shadow-lg); }
    .download-btn svg { width: 18px; height: 18px; }

    /* SA conditional fields */
    .sa-fields { margin-top: .8rem; }
    .sa-fields.hidden { display: none; }

    /* Module group collapsibles */
    details.module-group {
        margin-bottom: 1rem; border: 1px solid var(--border);
        border-radius: 10px; overflow: hidden;
    }
    details.module-group summary {
        padding: .7rem 1rem; background: var(--surface); font-weight: 600;
        font-size: .9rem; cursor: pointer; border-bottom: 1px solid transparent;
        list-style: none; display: flex; align-items: center; gap: .5rem;
    }
    details.module-group summary::before {
        content: '\25B6'; font-size: .7rem; transition: transform .2s;
    }
    details.module-group[open] summary::before { transform: rotate(90deg); }
    details.module-group[open] summary { border-bottom-color: var(--border); }
    details.module-group .module-content { padding: .5rem; }
    details.module-group summary .group-badge {
        background: var(--primary); color: #fff; font-size: .7rem;
        padding: .1rem .5rem; border-radius: 10px; font-weight: 500;
    }

    /* ── Responsive ───────────────────────────────────────────── */
    @media (max-width: 900px) {
        .sidebar { display: none; }
        .kpi-grid { grid-template-columns: repeat(2, 1fr); }
    }

    /* ── Plot Container ──────────────────────────────────────── */
    .plot-container { width: 100%; min-height: 400px; margin: 1rem 0; }
    .plot-img { width: 100%; border-radius: 10px; border: 1px solid var(--border); }

    /* ── Progress Bar ─────────────────────────────────────────── */
    .progress-bar-container {
        width: 100%; background: var(--border); border-radius: 10px;
        overflow: hidden; height: 8px; margin: .5rem 0;
    }
    .progress-bar {
        height: 100%; border-radius: 10px;
        background: linear-gradient(90deg, var(--primary), var(--secondary));
        transition: width .3s;
    }

    /* ── BGC Reaction Network Diagram ─────────────────────── */
    .bgc-diagram-container {
        overflow-x: auto; margin: 1rem 0;
        border: 1px solid var(--border); border-radius: 12px;
        padding: 1rem 1rem .6rem; background: var(--surface);
    }
    .bgc-diagram-container svg {
        display: block; margin: 0 auto; max-width: 100%;
        font-family: 'Inter', -apple-system, BlinkMacSystemFont, sans-serif;
    }
    .bgc-diagram-header {
        display: flex; justify-content: space-between; align-items: center;
        margin-bottom: .4rem;
    }
    .bgc-diagram-title {
        font-size: .9rem; font-weight: 600;
        color: var(--text);
    }
    .bgc-zoom-bar {
        display: flex; align-items: center; gap: 4px;
    }
    .zoom-btn {
        width: 26px; height: 26px; border-radius: 6px;
        border: 1px solid var(--border); background: var(--surface);
        color: var(--text); font-size: .85rem; cursor: pointer;
        display: inline-flex; align-items: center; justify-content: center;
        transition: background .15s; font-family: inherit; line-height: 1;
        padding: 0;
    }
    .zoom-btn:hover { background: var(--border); }
    .zoom-btn.zoom-lock { font-size: .7rem; }
    .zoom-btn.zoom-lock.locked { background: var(--border); }
    .zoom-btn[disabled] { opacity: .35; cursor: default; }
    .zoom-btn[disabled]:hover { background: var(--surface); }
    .zoom-hint {
        font-size: .62rem; color: var(--text3); margin-left: 4px;
    }
    /* ── Framework toggle bar ── */
    .bgc-fw-toggle-bar {
        display: flex; flex-wrap: wrap; gap: 6px;
        margin-bottom: .6rem; align-items: center;
    }
    .fw-toggle-chip {
        display: inline-flex; align-items: center; gap: 4px;
        padding: 3px 10px; font-size: .72rem; font-weight: 500;
        border: 1.5px solid var(--chip-color, var(--primary));
        border-radius: 14px; background: transparent;
        color: var(--chip-color, var(--primary));
        cursor: pointer; transition: background .2s, opacity .2s;
        font-family: inherit; line-height: 1.4;
    }
    .fw-toggle-chip.active {
        background: var(--chip-color, var(--primary));
        color: var(--surface, #fff);
    }
    .fw-toggle-chip:not(.active) { opacity: .45; }
    .fw-toggle-chip[data-fw="__all__"] { --chip-color: var(--text2); }
    .chip-swatch {
        width: 10px; height: 10px; border-radius: 50%;
        display: inline-block;
    }
    .fw-toggle-chip.active .chip-swatch {
        border: 1.5px solid var(--surface, #fff);
    }
    .chip-tick {
        font-size: .7rem; line-height: 1; margin-right: 1px;
        display: inline-block; transition: opacity .15s;
    }
    .fw-toggle-chip:not(.active) .chip-tick { opacity: 0; width: 0; margin: 0; overflow: hidden; }
    /* ── Hover tooltip ── */
    .bgc-tooltip {
        position: fixed; z-index: 9999; pointer-events: none;
        max-width: 420px; padding: 8px 12px;
        background: var(--surface, #1a1b2e); color: var(--text, #e2e8f0);
        border: 1px solid var(--border, #2a2b3d);
        border-radius: 8px; font-size: .78rem; line-height: 1.55;
        box-shadow: 0 6px 24px rgba(0,0,0,.25);
        font-family: 'Inter', -apple-system, BlinkMacSystemFont, sans-serif;
    }
    .bgc-tooltip strong { display: block; margin-bottom: 3px; font-weight: 600; }
    .bgc-tip-reaction {
        display: block; color: var(--text2, #a0aec0); font-size: .72rem;
        margin-bottom: 4px;
    }
    .bgc-tip-expr {
        display: block; font-family: 'JetBrains Mono', monospace;
        font-size: .7rem; color: var(--secondary, #34d399);
        background: rgba(0,0,0,.15); padding: 3px 6px;
        border-radius: 4px; margin: 3px 0; word-break: break-all;
    }
    .bgc-tip-params {
        display: block; font-size: .68rem; color: var(--text3, #636e72);
        margin-top: 2px;
    }
    /* ── Diagram interactivity ── */
    .bgc-diagram-container .reaction-arrow {
        transition: opacity .25s ease, filter .25s ease; cursor: pointer;
    }
    .bgc-diagram-container .species-node {
        transition: opacity .25s ease, filter .25s ease;
        cursor: pointer;
    }
    .bgc-diagram-container svg { cursor: crosshair; }
    .bgc-diagram-container .reaction-arrow path.hit-area { pointer-events: stroke; }
    .bgc-diagram-container .reaction-arrow path.vis-path { pointer-events: none; }
    .bgc-diagram-container .reaction-arrow text { pointer-events: none; }
    .bgc-diagram-container .label-bg { pointer-events: none; }
    """


def build_html_head(title: str, extra_css: str = "") -> str:
    """Build the <head> section of the HTML report."""
    return f"""<!DOCTYPE html>
<html lang="en" data-theme="dark">
<head>
<meta charset="utf-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1"/>
<title>{title}</title>
<link rel="preconnect" href="https://fonts.googleapis.com"/>
<link href="https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&family=JetBrains+Mono:wght@400;500&display=swap" rel="stylesheet"/>
<style>
{get_css_styles()}
{extra_css}
</style>
</head>
"""


def build_header(
    title: str,
    subtitle: str,
    meta_items: List[str],
    badge_text: str = "",
    badge_class: str = "badge-primary"
) -> str:
    """Build the page header with title, subtitle, and metadata."""
    badge_html = ""
    if badge_text:
        badge_html = f'<span class="badge {badge_class}" style="font-size:.85rem;padding:.3rem .8rem;margin-left:.8rem;">{badge_text}</span>'

    meta_html = "".join(f"<span>{m}</span>" for m in meta_items)

    return f"""
<div class="header">
    <h1>{title}{badge_html}</h1>
    <div class="subtitle">{subtitle}</div>
    <div class="meta">{meta_html}</div>
</div>
"""


def build_sidebar(
    nav_items: List[Dict[str, str]],
    logo_text: str = "OpenWQ"
) -> str:
    """Build the sidebar navigation.

    Parameters
    ----------
    nav_items : List[Dict]
        Each dict has keys: 'id' (section id), 'label' (display text)
    logo_text : str
        Logo text shown at top of sidebar
    """
    nav_links = "\n".join(
        f'<a href="#{item["id"]}">{item["label"]}</a>'
        for item in nav_items
    )

    return f"""
<aside class="sidebar">
    <div class="logo">{logo_text}</div>
    <nav>
{nav_links}
    </nav>
    <div class="theme-toggle">
        <button class="theme-btn" onclick="toggleTheme()">
            <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                <circle cx="12" cy="12" r="5"/><line x1="12" y1="1" x2="12" y2="3"/>
                <line x1="12" y1="21" x2="12" y2="23"/><line x1="4.22" y1="4.22" x2="5.64" y2="5.64"/>
                <line x1="18.36" y1="18.36" x2="19.78" y2="19.78"/><line x1="1" y1="12" x2="3" y2="12"/>
                <line x1="21" y1="12" x2="23" y2="12"/><line x1="4.22" y1="19.78" x2="5.64" y2="18.36"/>
                <line x1="18.36" y1="5.64" x2="19.78" y2="4.22"/>
            </svg>
            Toggle Theme
        </button>
    </div>
</aside>
"""


def build_footer() -> str:
    """Build the page footer."""
    now = datetime.now().strftime("%Y-%m-%d %H:%M")
    return f"""
<div class="footer">
    OpenWQ &mdash; Open Water Quality Model &bull; Report generated {now}
    &bull; <a href="https://openwq.readthedocs.io" style="color:var(--primary);">Documentation</a>
</div>
"""


def build_theme_toggle_js() -> str:
    """Build the theme toggle JavaScript."""
    return """
<script>
function toggleTheme() {
    const html = document.documentElement;
    const current = html.getAttribute('data-theme');
    const next = current === 'dark' ? 'light' : 'dark';
    html.setAttribute('data-theme', next);
    localStorage.setItem('owq_theme', next);
}
// Restore saved theme
(function() {
    const saved = localStorage.getItem('owq_theme');
    if (saved) document.documentElement.setAttribute('data-theme', saved);
})();

// Copy button handler
document.addEventListener('click', function(e) {
    if (e.target.classList.contains('copy-btn')) {
        const pre = e.target.closest('pre');
        const code = pre.querySelector('code') || pre;
        navigator.clipboard.writeText(code.textContent.trim());
        e.target.textContent = 'Copied!';
        setTimeout(() => { e.target.textContent = 'Copy'; }, 2000);
    }
});
</script>
"""


def build_kpi(icon: str, value: str, label: str) -> str:
    """Build a single KPI card."""
    return f"""
<div class="kpi">
    <div class="icon">{icon}</div>
    <div class="value">{value}</div>
    <div class="label">{label}</div>
</div>
"""


def build_kpi_grid(kpis: List[Dict[str, str]]) -> str:
    """Build a KPI grid from a list of KPI dicts.

    Each dict has keys: 'icon', 'value', 'label'
    """
    items = "".join(
        build_kpi(k["icon"], k["value"], k["label"])
        for k in kpis
    )
    return f'<div class="kpi-grid">{items}</div>'


def build_param_table(
    parameters: List[Dict],
    columns: List[Dict[str, str]],
    table_id: str = ""
) -> str:
    """Build an HTML parameter table.

    Parameters
    ----------
    parameters : List[Dict]
        Parameter data (list of dicts)
    columns : List[Dict]
        Each dict has 'key' (data key), 'label' (header text),
        'align' ('left' or 'right'), and optional 'format' (callable)
    table_id : str
        Optional table ID attribute
    """
    id_attr = f' id="{table_id}"' if table_id else ""

    headers = "".join(
        f'<th style="text-align:{c.get("align", "left")}">{c["label"]}</th>'
        for c in columns
    )

    rows = []
    for param in parameters:
        cells = []
        for col in columns:
            val = param.get(col["key"], "")
            fmt = col.get("format")
            if fmt and val != "":
                val = fmt(val)
            align = col.get("align", "left")
            cls = ' class="num"' if align == "right" else ""
            cells.append(f"<td{cls}>{val}</td>")
        rows.append("<tr>" + "".join(cells) + "</tr>")

    return f"""
<div class="table-wrap">
<table class="param-table"{id_attr}>
<thead><tr>{headers}</tr></thead>
<tbody>
{"".join(rows)}
</tbody>
</table>
</div>
"""


def build_code_block(code: str, language: str = "") -> str:
    """Build a code block with copy button."""
    import html as html_lib
    escaped = html_lib.escape(code)
    return f"""
<pre class="code-block">
<button class="copy-btn">Copy</button>
<code>{escaped}</code>
</pre>
"""


def build_highlight_box(content: str, style: str = "info") -> str:
    """Build a highlight box (info, success, warning)."""
    return f'<div class="highlight-box {style}">{content}</div>'


def build_error_card(title: str, message: str) -> str:
    """Build an error display card."""
    import html as html_lib
    return f"""
<div class="card" style="border-left:5px solid #e74c3c;">
    <h3 style="color:#e74c3c;">\u26a0 {html_lib.escape(title)}</h3>
    <p style="color:var(--text2);font-size:.9rem;">{html_lib.escape(message)}</p>
</div>
"""


def embed_image_base64(image_path: str) -> str:
    """Read an image file and return a base64-encoded data URI."""
    if not os.path.isfile(image_path):
        return ""

    ext = os.path.splitext(image_path)[1].lower()
    mime_map = {'.png': 'image/png', '.jpg': 'image/jpeg',
                '.jpeg': 'image/jpeg', '.svg': 'image/svg+xml',
                '.gif': 'image/gif'}
    mime = mime_map.get(ext, 'image/png')

    with open(image_path, 'rb') as f:
        data = base64.b64encode(f.read()).decode('ascii')

    return f"data:{mime};base64,{data}"


# =========================================================================
# Interactive Form Helpers
# =========================================================================

def build_form_select(
    field_id: str,
    label: str,
    options: List[str],
    default: str = "",
    hint: str = ""
) -> str:
    """Build a form group with a <select> dropdown.

    Parameters
    ----------
    field_id : str
        The HTML id for the <select> element.
    label : str
        Label text displayed above the dropdown.
    options : List[str]
        List of option values (used as both value and display text).
    default : str
        Which option is selected by default.
    hint : str
        Optional hint text shown below the dropdown.
    """
    import html as html_lib
    opts_html = "\n".join(
        f'<option value="{html_lib.escape(o)}"'
        f'{" selected" if o == default else ""}>'
        f'{html_lib.escape(o)}</option>'
        for o in options
    )
    hint_html = f'<div class="hint">{html_lib.escape(hint)}</div>' if hint else ""
    return f"""
<div class="form-group">
    <label for="{field_id}">{html_lib.escape(label)}</label>
    <select class="form-select" id="{field_id}" name="{field_id}">
        {opts_html}
    </select>
    {hint_html}
</div>
"""


def build_form_number(
    field_id: str,
    label: str,
    default: Any = 0,
    min_val: Optional[Any] = None,
    max_val: Optional[Any] = None,
    step: Optional[Any] = None,
    hint: str = ""
) -> str:
    """Build a form group with a number <input>.

    Parameters
    ----------
    field_id : str
        The HTML id for the <input> element.
    label : str
        Label text.
    default : Any
        Default value.
    min_val, max_val, step : Optional[Any]
        HTML input attributes for min, max, step.
    hint : str
        Optional hint text.
    """
    import html as html_lib
    attrs = f'value="{default}"'
    if min_val is not None:
        attrs += f' min="{min_val}"'
    if max_val is not None:
        attrs += f' max="{max_val}"'
    if step is not None:
        attrs += f' step="{step}"'
    hint_html = f'<div class="hint">{html_lib.escape(hint)}</div>' if hint else ""
    return f"""
<div class="form-group">
    <label for="{field_id}">{html_lib.escape(label)}</label>
    <input class="form-input" type="number" id="{field_id}" name="{field_id}" {attrs}/>
    {hint_html}
</div>
"""


def build_form_checkbox(
    field_id: str,
    label: str,
    checked: bool = False,
    hint: str = ""
) -> str:
    """Build a form group with a checkbox.

    Parameters
    ----------
    field_id : str
        The HTML id for the <input> element.
    label : str
        Label text shown next to the checkbox.
    checked : bool
        Whether the checkbox is checked by default.
    hint : str
        Optional hint text.
    """
    import html as html_lib
    checked_attr = " checked" if checked else ""
    hint_html = f'<div class="hint">{html_lib.escape(hint)}</div>' if hint else ""
    return f"""
<div class="form-group" style="display:flex;align-items:center;gap:.5rem;">
    <input type="checkbox" id="{field_id}" name="{field_id}"{checked_attr}
           style="accent-color:var(--primary);width:16px;height:16px;"/>
    <label for="{field_id}" style="margin-bottom:0;">{html_lib.escape(label)}</label>
    {hint_html}
</div>
"""


def build_species_checkboxes(
    species_list: List[str],
    species_obs_availability: Optional[Dict[str, Dict]] = None,
) -> str:
    """Build a grid of species checkboxes with weight inputs.

    Parameters
    ----------
    species_list : List[str]
        List of chemical species names from the model config.
    species_obs_availability : Dict[str, Dict], optional
        Per-species observation availability from
        ``config_integration.get_species_observation_availability()``.
        Each entry has ``has_obs`` (bool) and ``details`` (str).
    """
    import html as html_lib
    obs_avail = species_obs_availability or {}
    items = []
    for sp in species_list:
        safe_id = sp.replace("-", "_").replace(" ", "_")
        info = obs_avail.get(sp, {})
        has_obs = info.get("has_obs", True)  # default True if no info
        details = html_lib.escape(info.get("details", ""))

        # Visual indicator, checked/disabled state
        if has_obs:
            item_style = ""
            cb_attrs = "checked"
            weight_disabled = ""
            badge = (
                f'<span style="font-size:.68rem;color:#00a86b;'
                f'margin-left:.4rem;font-weight:600;" '
                f'title="{details}">&#x2713; obs</span>'
            )
        else:
            item_style = (
                "opacity:.45;pointer-events:none;"
                "filter:grayscale(1);cursor:not-allowed;"
            )
            cb_attrs = "disabled"
            weight_disabled = "disabled"
            badge = (
                f'<span style="font-size:.68rem;color:var(--text3);'
                f'margin-left:.4rem;" '
                f'title="{details}">no obs</span>'
            )

        items.append(f"""
<div class="species-item" style="{item_style}">
    <input type="checkbox" class="species-cb" id="sp_{safe_id}"
           data-species="{html_lib.escape(sp)}" {cb_attrs}/>
    <label for="sp_{safe_id}" style="margin:0;font-size:.85rem;">
        {html_lib.escape(sp)}</label>
    {badge}
    <span style="font-size:.75rem;color:var(--text3);margin-left:.3rem;">w:</span>
    <input type="number" class="weight-input" id="w_{safe_id}"
           data-species="{html_lib.escape(sp)}"
           value="1.0" min="0" step="0.1" {weight_disabled}
           style="width:55px;padding:.15rem .3rem;font-size:.8rem;
                  border:1px solid var(--border);border-radius:4px;
                  background:var(--surface);color:var(--text);
                  font-family:'JetBrains Mono',monospace;text-align:right;"/>
</div>""")

    toggle = ('<a href="#" id="speciesToggleAll" '
              'style="font-size:.8rem;color:var(--primary);cursor:pointer;">'
              'Deselect all</a>')

    return f"""
<div style="margin-bottom:.5rem;">{toggle}</div>
<div class="species-grid">
{"".join(items)}
</div>
"""


def build_editable_param_table(parameters: List[Dict],
                               idx_offset: int = 0) -> str:
    """Build an editable parameter table with inline inputs.

    Each row has: checkbox (include), name, initial, min bound, max bound,
    transform select, units, description.

    Parameters
    ----------
    parameters : List[Dict]
        Auto-extracted calibration parameters.
    idx_offset : int
        Offset added to each row index so that data-param-idx values
        are globally unique when multiple tables are rendered.
    """
    import html as html_lib

    headers = """
<tr>
    <th style="width:30px;text-align:center;">
        <input type="checkbox" class="groupToggleAll" checked
               style="accent-color:var(--primary);"/>
    </th>
    <th>Parameter</th>
    <th style="text-align:right;">Initial</th>
    <th style="text-align:right;">Min</th>
    <th style="text-align:right;">Max</th>
    <th>Transform</th>
    <th>Units</th>
    <th>Description</th>
</tr>"""

    rows = []
    for i, p in enumerate(parameters):
        gi = i + idx_offset
        name = p.get("name", "")
        initial = p.get("initial", 0)
        bounds = p.get("bounds", (0, 0))
        transform = p.get("transform", "linear")
        units = p.get("units", "")
        desc = p.get("description", "")

        # Format initial value
        if isinstance(initial, float):
            if abs(initial) < 0.01 or abs(initial) >= 10000:
                initial_str = f"{initial:.4g}"
            else:
                initial_str = f"{initial:.4f}"
        else:
            initial_str = str(initial)

        # Log/linear select
        log_sel = " selected" if transform == "log" else ""
        lin_sel = " selected" if transform != "log" else ""

        # Serialise the full param dict as a data attribute for JS
        param_json = _js(p)

        # Non-calibratable params: disabled checkbox + note
        calibratable = p.get("_calibratable", True)
        has_explicit_range = p.get("has_explicit_range", True)
        row_cls = "param-row" if calibratable else "param-row param-no-obs"
        row_style = "" if calibratable else ' style="opacity:0.45;"'
        notes = []
        if not calibratable:
            notes.append(
                '<span style="color:#c47800;font-size:.7rem;'
                'font-style:italic;">'
                '(no observation data to calibrate this sub-cycle)'
                '</span>')
        if not has_explicit_range:
            notes.append(
                '<span style="color:#888;font-size:.7rem;'
                'font-style:italic;">'
                '(auto-generated bounds &mdash; review before calibrating)'
                '</span>')
        no_obs_note = " ".join(notes)
        if no_obs_note:
            no_obs_note = " " + no_obs_note
        cb_checked = " checked" if calibratable else ""
        cb_disabled = "" if calibratable else " disabled"
        fw_attr = ""
        if p.get("_framework"):
            fw_attr = f' data-fw="{html_lib.escape(p["_framework"])}"'

        rows.append(f"""
<tr class="{row_cls}" data-param-idx="{gi}" data-param='{html_lib.escape(param_json)}'{fw_attr}{row_style}>
    <td style="text-align:center;">
        <input type="checkbox" class="param-cb" data-idx="{gi}"{cb_checked}{cb_disabled}
               style="accent-color:var(--primary);"/>
    </td>
    <td><code style="font-size:.8rem;">{html_lib.escape(name)}</code>{no_obs_note}</td>
    <td class="num">{initial_str}</td>
    <td><input type="number" class="inline-input param-min" data-idx="{gi}"
               value="{bounds[0]}" step="any"/></td>
    <td><input type="number" class="inline-input param-max" data-idx="{gi}"
               value="{bounds[1]}" step="any"/></td>
    <td>
        <select class="inline-select param-transform" data-idx="{gi}">
            <option value="log"{log_sel}>log</option>
            <option value="linear"{lin_sel}>linear</option>
        </select>
    </td>
    <td style="font-size:.8rem;">{html_lib.escape(units)}</td>
    <td style="font-size:.8rem;color:var(--text2);max-width:200px;">
        {html_lib.escape(desc)}</td>
</tr>""")

    return f"""
<div class="table-wrap">
<table class="param-table" id="paramTable">
<thead>{headers}</thead>
<tbody>
{"".join(rows)}
</tbody>
</table>
</div>
"""


def build_download_button(button_id: str = "downloadBtn") -> str:
    """Build a download button with an icon."""
    return f"""
<button class="download-btn" id="{button_id}" onclick="downloadScript()">
    <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"
         stroke-linecap="round" stroke-linejoin="round">
        <path d="M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4"/>
        <polyline points="7 10 12 15 17 10"/>
        <line x1="12" y1="15" x2="12" y2="3"/>
    </svg>
    Download Script
</button>
"""


def get_interactive_layout_css() -> str:
    """Return extra CSS for the split-pane interactive report layout.

    This is passed via ``extra_css`` to ``build_html_head()`` and is
    additive — the base ``get_css_styles()`` remains unchanged.
    """
    return """
    /* ── Split-Pane Layout ──────────────────────────────────── */
    .split-layout {
        display: flex;
        flex-direction: column;
        height: 100vh;
        overflow: hidden;
    }
    .panes-row {
        display: flex;
        flex: 1;
        min-height: 0;
        overflow: hidden;
    }

    /* ── Top Bar ────────────────────────────────────────────── */
    .top-bar {
        background: linear-gradient(135deg, var(--dark) 0%, #16213e 60%, var(--primary-dark) 100%);
        color: #fff;
        padding: .8rem 1.5rem;
        display: flex;
        align-items: center;
        justify-content: space-between;
        gap: 1rem;
        min-height: 56px;
        flex-shrink: 0;
    }
    .top-bar-left {
        display: flex;
        align-items: center;
        gap: .8rem;
        min-width: 0;
    }
    .top-bar h1 {
        font-size: 1.15rem;
        font-weight: 700;
        white-space: nowrap;
        margin: 0;
    }
    .top-bar .subtitle {
        font-size: .85rem;
        opacity: .7;
        white-space: nowrap;
        overflow: hidden;
        text-overflow: ellipsis;
    }
    .top-bar-right {
        display: flex;
        align-items: center;
        gap: 1rem;
        flex-shrink: 0;
    }

    /* ── Progress Segments ──────────────────────────────────── */
    .progress-heading {
        font-size: .68rem;
        text-transform: uppercase;
        letter-spacing: .04em;
        opacity: .55;
        margin-right: .4rem;
        white-space: nowrap;
    }
    .progress-segments {
        display: inline-flex;
        align-items: center;
        gap: 4px;
    }
    .progress-seg {
        width: 36px;
        height: 6px;
        border-radius: 3px;
        background: rgba(255,255,255,.2);
        transition: background .3s;
    }
    .progress-seg.active {
        background: linear-gradient(90deg, var(--primary), var(--secondary));
    }
    .progress-label {
        font-size: .72rem;
        opacity: .6;
        margin-left: .4rem;
        white-space: nowrap;
    }

    /* ── Config Pane (Left) ─────────────────────────────────── */
    .config-pane {
        overflow-y: auto;
        padding: 0;
        display: flex;
        flex-direction: column;
        width: 55%;
        min-width: 300px;
        flex-shrink: 0;
    }

    /* ── Tab Bar ────────────────────────────────────────────── */
    .tab-bar {
        display: flex;
        gap: 0;
        border-bottom: 2px solid var(--border);
        background: var(--surface);
        position: sticky;
        top: 0;
        z-index: 50;
        flex-shrink: 0;
    }
    .tab-btn {
        flex: 1;
        padding: .65rem .5rem;
        font-size: .82rem;
        font-weight: 500;
        text-align: center;
        cursor: pointer;
        background: none;
        border: none;
        border-bottom: 3px solid transparent;
        color: var(--text2);
        transition: all .15s;
        font-family: inherit;
        white-space: nowrap;
    }
    .tab-btn:hover {
        color: var(--primary);
        background: rgba(0,102,204,.04);
    }
    .tab-btn.active {
        color: var(--primary);
        border-bottom-color: var(--primary);
        font-weight: 600;
    }

    /* ── Tab Panel ──────────────────────────────────────────── */
    .tab-panel {
        display: none;
        padding: 1.2rem 1.5rem;
        flex: 1;
        overflow-y: auto;
    }
    .tab-panel.active {
        display: block;
    }

    /* ── Resize Handle ──────────────────────────────────────── */
    .pane-resizer {
        width: 5px;
        cursor: col-resize;
        background: var(--border);
        transition: background .15s;
        position: relative;
        z-index: 20;
        flex-shrink: 0;
    }
    .pane-resizer:hover,
    .pane-resizer.dragging {
        background: var(--primary);
    }
    .pane-resizer::after {
        content: '';
        position: absolute;
        top: 50%;
        left: 50%;
        transform: translate(-50%, -50%);
        width: 3px;
        height: 32px;
        border-radius: 2px;
        background: var(--text3);
        opacity: 0;
        transition: opacity .15s;
    }
    .pane-resizer:hover::after,
    .pane-resizer.dragging::after {
        opacity: .6;
    }

    /* ── Script Pane (Right) ────────────────────────────────── */
    .script-pane {
        overflow-y: auto;
        display: flex;
        flex-direction: column;
        background: var(--surface);
        flex: 1;
        min-width: 250px;
    }
    .script-pane-header {
        padding: .8rem 1.2rem;
        border-bottom: 1px solid var(--border);
        display: flex;
        align-items: center;
        justify-content: space-between;
        flex-shrink: 0;
        background: var(--surface);
        position: sticky;
        top: 0;
        z-index: 10;
    }
    .script-pane-header h3 {
        font-size: .92rem;
        font-weight: 600;
        margin: 0;
    }
    .script-pane-body {
        flex: 1;
        overflow-y: auto;
        padding: 0;
    }
    .script-pane-body pre.code-block {
        margin: 0;
        border-radius: 0;
        min-height: 100%;
        font-size: .78rem;
        line-height: 1.55;
    }

    /* ── Action Bar (Bottom) ────────────────────────────────── */
    .action-bar {
        display: flex;
        align-items: center;
        gap: 1rem;
        padding: .6rem 1.5rem;
        background: var(--surface);
        border-top: 1px solid var(--border);
        box-shadow: 0 -2px 10px rgba(0,0,0,.04);
        flex-shrink: 0;
    }
    .action-bar .download-btn {
        margin: 0;
        padding: .55rem 1.2rem;
        font-size: .88rem;
    }
    .action-status {
        font-size: .8rem;
        color: var(--text2);
        margin-left: auto;
    }
    .action-bar .theme-btn-compact {
        background: none;
        border: 1px solid var(--border);
        border-radius: 6px;
        padding: .3rem .5rem;
        cursor: pointer;
        color: var(--text2);
        display: flex;
        align-items: center;
    }
    .action-bar .theme-btn-compact:hover {
        border-color: var(--primary);
    }

    /* ── Parameter Search ───────────────────────────────────── */
    .param-search {
        width: 100%;
        padding: .5rem .8rem;
        font-size: .85rem;
        border: 1px solid var(--border);
        border-radius: 8px;
        background: var(--surface);
        color: var(--text);
        font-family: inherit;
        margin-bottom: 1rem;
        transition: border-color .15s;
    }
    .param-search:focus {
        outline: none;
        border-color: var(--primary);
        box-shadow: 0 0 0 3px rgba(0,102,204,.1);
    }
    .param-search::placeholder {
        color: var(--text3);
    }

    /* ── Syntax Highlighting ────────────────────────────────── */
    .syn-kw { color: #c792ea; }
    .syn-str { color: #c3e88d; }
    .syn-comment { color: #546e7a; font-style: italic; }
    .syn-num { color: #f78c6c; }
    .syn-fn { color: #82aaff; }
    .syn-const { color: #ff5370; }
    [data-theme="light"] .syn-kw { color: #7c3aed; }
    [data-theme="light"] .syn-str { color: #16a34a; }
    [data-theme="light"] .syn-comment { color: #94a3b8; font-style: italic; }
    [data-theme="light"] .syn-num { color: #ea580c; }
    [data-theme="light"] .syn-fn { color: #2563eb; }
    [data-theme="light"] .syn-const { color: #dc2626; }

    /* ── Responsive: single-column below 1000px ─────────────── */
    @media (max-width: 1000px) {
        .split-layout {
            height: auto;
            min-height: 100vh;
        }
        .panes-row {
            flex-direction: column;
        }
        .config-pane {
            width: 100% !important;
            border-bottom: 1px solid var(--border);
            max-height: 60vh;
        }
        .pane-resizer {
            display: none;
        }
        .script-pane {
            max-height: 50vh;
        }
    }
    """


def build_compact_header(
    title: str,
    subtitle: str = "",
    badge_text: str = "",
    badge_class: str = "badge-primary",
) -> str:
    """Build a compact single-line header for the split-pane layout."""
    badge_html = ""
    if badge_text:
        badge_html = (
            f'<span class="badge {badge_class}" '
            f'style="font-size:.72rem;padding:.15rem .5rem;margin-left:.5rem;">'
            f'{badge_text}</span>'
        )

    subtitle_html = (
        f'\n    <span class="subtitle">{subtitle}</span>'
        if subtitle else ""
    )

    return f"""
<div class="top-bar-left">
    <h1>{title}{badge_html}</h1>{subtitle_html}
</div>
"""


# =====================================================================
#  BGC Reaction Network Diagram Generator
# =====================================================================

# PHREEQC element notation → friendly label
_PHREEQC_LABELS = {
    "N(-3)": "NH\u2084\u207a", "N(5)": "NO\u2083\u207b",
    "N(3)": "NO\u2082\u207b", "N(0)": "N\u2082",
    "O(0)": "O\u2082", "C(4)": "CO\u2083/HCO\u2083",
    "C(-4)": "CH\u2084", "S(6)": "SO\u2084\u00b2\u207b",
    "S(-2)": "H\u2082S", "Fe(2)": "Fe\u00b2\u207a",
    "Fe(3)": "Fe\u00b3\u207a", "Mn(2)": "Mn\u00b2\u207a",
    "P(5)": "PO\u2084\u00b3\u207b",
}

# Color palette for multi-framework diagrams (arrow colour, bg fill)
_FW_COLORS = [
    ("#4d9ee8", "rgba(77,158,232,.10)"),   # blue
    ("#34d399", "rgba(52,211,153,.10)"),    # green
    ("#fb923c", "rgba(251,146,60,.10)"),    # orange
    ("#a78bfa", "rgba(167,139,250,.10)"),   # purple
    ("#f472b6", "rgba(244,114,182,.10)"),   # pink
    ("#2dd4bf", "rgba(45,212,191,.10)"),    # teal
]


def extract_bgc_network(bgc_data: dict) -> dict:
    """Extract a reaction network from a NATIVE_BGC_FLEX JSON dict.

    Accepts both top-level format (CHEMICAL_SPECIES at root) and nested
    format (under BIOGEOCHEMISTRY_CONFIGURATION).

    Returns ``{"species": [...], "frameworks": {...}}`` or empty dict.
    """
    if not bgc_data or not isinstance(bgc_data, dict):
        return {}

    # Unwrap nested format
    if "BIOGEOCHEMISTRY_CONFIGURATION" in bgc_data:
        bgc = bgc_data["BIOGEOCHEMISTRY_CONFIGURATION"]
    else:
        bgc = bgc_data

    # --- Species ---
    species = []
    chem = bgc.get("CHEMICAL_SPECIES", {})
    descs = chem.get("_LIST_DESCRIPTIONS", {})
    sp_list = chem.get("LIST", {})
    if sp_list:
        for idx in sorted(sp_list.keys(), key=lambda x: int(x)):
            name = sp_list[idx]
            species.append({
                "name": name,
                "description": descs.get(name, ""),
            })
    if not species:
        return {}

    # --- Frameworks ---
    frameworks = {}
    fw_data = bgc.get("CYCLING_FRAMEWORKS", {})
    for fw_name, fw in fw_data.items():
        if not isinstance(fw, dict):
            continue
        desc = fw.get("_DESCRIPTION", fw_name)
        reactions = []
        # Runtime format with LIST_TRANSFORMATIONS
        lt = fw.get("LIST_TRANSFORMATIONS", {})
        if lt:
            for idx in sorted(lt.keys(), key=lambda x: int(x)):
                t = fw.get(idx, {})
                rxn: dict = {
                    "name": t.get("_NAME", lt[idx]),
                    "consumed": t.get("CONSUMED", "NONE"),
                    "produced": t.get("PRODUCED", "NONE"),
                }
                # Capture kinetics expression & parameter info
                kin = t.get("KINETICS", [])
                if isinstance(kin, list) and kin:
                    rxn["kinetics_expr"] = str(kin[0])
                    if len(kin) > 1:
                        rxn["kinetics_units"] = str(kin[1])
                pinfo = t.get("_PARAMETERS_INFO", {})
                if pinfo:
                    params = []
                    for pname, pd in pinfo.items():
                        if isinstance(pd, dict):
                            pstr = pname
                            val = pd.get("VALUE", "")
                            units = pd.get("UNITS", "")
                            if val != "":
                                pstr += f" = {val}"
                            if units:
                                pstr += f" {units}"
                            rng = pd.get("RANGE")
                            if isinstance(rng, list) and len(rng) == 2:
                                pstr += f" [{rng[0]}\u2013{rng[1]}]"
                            params.append(pstr)
                    if params:
                        rxn["parameters"] = params
                reactions.append(rxn)
        else:
            # Legacy dict format
            for key, val in fw.items():
                if key.startswith("_") or not isinstance(val, dict):
                    continue
                if "KINETICS" not in val and "FORMULA" not in val:
                    continue
                rxn = {
                    "name": key,
                    "consumed": val.get("CONSUMED", "NONE"),
                    "produced": val.get("PRODUCED", "NONE"),
                }
                kin = val.get("KINETICS", [])
                if isinstance(kin, list) and kin:
                    rxn["kinetics_expr"] = str(kin[0])
                    if len(kin) > 1:
                        rxn["kinetics_units"] = str(kin[1])
                reactions.append(rxn)
        if reactions:
            frameworks[fw_name] = {"description": desc, "reactions": reactions}

    if not frameworks:
        return {}

    return {"species": species, "frameworks": frameworks}


def extract_bgc_network_from_pqi(pqi_path: str) -> dict:
    """Extract a reaction network from a PHREEQC .pqi file.

    Parses KINETICS blocks for reaction names and ``-formula`` lines.
    Returns the same structure as :func:`extract_bgc_network`.
    """
    if not pqi_path or not os.path.isfile(pqi_path):
        return {}

    with open(pqi_path, "r", errors="replace") as fh:
        lines = fh.readlines()

    reactions = []
    species_set = set()
    pqi_parms: Dict[str, list] = {}
    in_kinetics = False
    current_rxn = None

    for raw in lines:
        line = raw.strip()
        # Detect KINETICS block start
        if re.match(r"^KINETICS\b", line, re.IGNORECASE):
            in_kinetics = True
            current_rxn = None
            continue
        # Detect block end
        if in_kinetics and (line.upper() == "END"
                           or re.match(r"^(RATES|SOLUTION|SELECTED_OUTPUT|"
                                       r"EQUILIBRIUM_PHASES|USER_PUNCH)\b",
                                       line, re.IGNORECASE)):
            in_kinetics = False
            current_rxn = None
            continue
        if not in_kinetics:
            continue
        # Reaction name: non-indented line that is not a keyword
        if not raw[0:1].isspace() and line and not line.startswith("-"):
            current_rxn = line.split()[0]
            continue
        # -parms line (capture parameter values)
        if current_rxn and line.lower().startswith("-parms"):
            pqi_parms[current_rxn] = line.split()[1:]
            continue
        # -formula line
        if current_rxn and line.lower().startswith("-formula"):
            formula_raw = line[len("-formula"):].strip()
            tokens = line.split()[1:]  # e.g. ["N(-3)", "-1", "N(5)", "1"]
            consumed, produced = [], []
            i = 0
            while i < len(tokens) - 1:
                elem = tokens[i]
                try:
                    coeff = float(tokens[i + 1])
                except (ValueError, IndexError):
                    i += 1
                    continue
                label = _PHREEQC_LABELS.get(elem, elem)
                species_set.add(label)
                if coeff < 0:
                    consumed.append(label)
                elif coeff > 0:
                    produced.append(label)
                i += 2
            for c in (consumed or ["NONE"]):
                for p in (produced or ["NONE"]):
                    rxn: dict = {
                        "name": current_rxn,
                        "consumed": c,
                        "produced": p,
                    }
                    if formula_raw:
                        rxn["kinetics_expr"] = formula_raw
                    reactions.append(rxn)

    if not reactions:
        return {}

    # Attach -parms values to reactions
    for rxn in reactions:
        pvals = pqi_parms.get(rxn["name"])
        if pvals:
            rxn["parameters"] = [f"PARM({i + 1}) = {v}"
                                 for i, v in enumerate(pvals)]

    species = [{"name": s, "description": ""} for s in sorted(species_set)]
    return {
        "species": species,
        "frameworks": {"PHREEQC_Reactions": {
            "description": "Reactions parsed from PHREEQC input",
            "reactions": reactions,
        }},
    }


def _compute_diagram_layout(
    species_names: List[str],
    all_reactions: List[dict],
    width: int = 860,
    height: int = 500,
) -> Dict[str, tuple]:
    """Compute (x, y) positions for species nodes using a circular layout.

    Nodes involved in reactions are placed on the main ellipse.
    Unreferenced species are placed in a secondary row below.
    """
    involved = set()
    for r in all_reactions:
        if r["consumed"] != "NONE":
            involved.add(r["consumed"])
        if r["produced"] != "NONE":
            involved.add(r["produced"])

    main_nodes = [s for s in species_names if s in involved]
    extra_nodes = [s for s in species_names if s not in involved]

    positions = {}
    cx, cy = width / 2, height / 2 - 10
    n = max(len(main_nodes), 1)
    rx = min(width * 0.38, 340)
    ry = min(height * 0.34, 200)

    for i, sp in enumerate(main_nodes):
        angle = 2 * math.pi * i / n - math.pi / 2
        x = cx + rx * math.cos(angle)
        y = cy + ry * math.sin(angle)
        positions[sp] = (x, y)

    if extra_nodes:
        y_extra = height - 35
        spacing = min(width / (len(extra_nodes) + 1), 140)
        x_start = (width - spacing * (len(extra_nodes) - 1)) / 2
        for i, sp in enumerate(extra_nodes):
            positions[sp] = (x_start + i * spacing, y_extra)

    return positions


def _svg_escape(text: str) -> str:
    """Escape text for safe embedding in SVG/XML."""
    return html_lib.escape(str(text), quote=True)


# ── Label collision avoidance helpers ────────────────────────────────

def _bezier_point(t: float, a: float, q: float, b: float) -> float:
    """Evaluate a quadratic Bézier component at parameter *t*."""
    return (1 - t) ** 2 * a + 2 * (1 - t) * t * q + t ** 2 * b


def _boxes_overlap(a: tuple, b: tuple, pad: float = 2.0) -> bool:
    """Return True if two axis-aligned boxes (x, y, w, h) overlap."""
    ax, ay, aw, ah = a
    bx, by, bw, bh = b
    return not (ax + aw + pad < bx or bx + bw + pad < ax or
                ay + ah + pad < by or by + bh + pad < ay)


def _resolve_label_positions(labels: List[dict],
                             placed_boxes: List[tuple]) -> None:
    """Greedy label placement with collision avoidance.

    Each *label* dict must contain keys:
        x, y          – initial position (centre)
        w, h          – estimated bounding-box size
        ax, ay, qx, qy, bx, by  – Bézier control points
        nx, ny        – perpendicular unit vector for offsets

    The function mutates ``x`` and ``y`` in-place to resolved positions.
    *placed_boxes* is a list of (x, y, w, h) already occupied (e.g. nodes)
    and is appended to as labels are placed.
    """
    _CANDIDATES = [
        (0.50, 0), (0.35, 0), (0.65, 0),
        (0.50, 16), (0.50, -16),
        (0.30, 0), (0.70, 0),
        (0.35, 14), (0.65, 14),
        (0.35, -14), (0.65, -14),
        (0.20, 0), (0.80, 0),
        (0.50, 24), (0.50, -24),
    ]

    for lab in labels:
        w, h = lab["w"], lab["h"]
        best = None
        for t_val, perp_off in _CANDIDATES:
            cx = _bezier_point(t_val, lab["ax"], lab["qx"], lab["bx"])
            cy = _bezier_point(t_val, lab["ay"], lab["qy"], lab["by"])
            cx += lab["nx"] * perp_off
            cy += lab["ny"] * perp_off
            box = (cx - w / 2, cy - h / 2, w, h)
            if not any(_boxes_overlap(box, pb) for pb in placed_boxes):
                best = (cx, cy, box)
                break
        if best:
            lab["x"], lab["y"] = best[0], best[1]
            placed_boxes.append(best[2])
        else:
            # Fallback: use first candidate (midpoint) even if it overlaps
            cx = _bezier_point(0.5, lab["ax"], lab["qx"], lab["bx"])
            cy = _bezier_point(0.5, lab["ay"], lab["qy"], lab["by"])
            box = (cx - w / 2, cy - h / 2, w, h)
            lab["x"], lab["y"] = cx, cy
            placed_boxes.append(box)


# ── SVG renderer ─────────────────────────────────────────────────────

def _build_kinetics_tip(r: dict) -> str:
    """Build a tooltip string from reaction kinetics data."""
    parts: List[str] = []
    expr = r.get("kinetics_expr", "")
    if expr:
        units = r.get("kinetics_units", "")
        parts.append(expr + (f"  ({units})" if units else ""))
    params = r.get("parameters", [])
    if params:
        parts.append("Parameters: " + "; ".join(params))
    return " &#10; ".join(parts) if parts else ""


def _render_bgc_svg(
    network: dict,
    max_width: int = 860,
    species_obs_availability: Optional[dict] = None,
) -> str:
    """Render a BGC reaction network as an inline SVG string.

    Includes data-* attributes for JS interactivity, framework group
    wrappers, hit-area paths, and collision-free label placement.
    """

    all_species = [s["name"] for s in network["species"]]
    species_desc = {s["name"]: s.get("description", "") for s in network["species"]}
    # Collect reactions, sorted by framework for grouping
    all_reactions = []
    fw_names = list(network["frameworks"].keys())
    for fw_name in fw_names:
        for r in network["frameworks"][fw_name]["reactions"]:
            all_reactions.append({**r, "_fw": fw_name})

    if not all_reactions:
        return ""

    n_species = len(all_species)
    svg_w = max_width
    svg_h = max(380, min(250 + n_species * 22, 650))

    positions = _compute_diagram_layout(all_species, all_reactions, svg_w, svg_h)

    node_h = 28
    char_w = 8
    _LABEL_FONT = 8.5
    _LABEL_CHAR_W = 5.2

    # Pre-compute species node boxes for label collision avoidance
    placed_boxes: List[tuple] = []
    node_boxes: Dict[str, tuple] = {}
    for sp in all_species:
        if sp not in positions:
            continue
        x, y = positions[sp]
        nw = max(len(sp) * char_w + 24, 64)
        box = (x - nw / 2, y - node_h / 2, nw, node_h)
        placed_boxes.append(box)
        node_boxes[sp] = box

    # ── Build SVG ──
    S: List[str] = []
    S.append(f'<svg xmlns="http://www.w3.org/2000/svg" '
             f'viewBox="0 0 {svg_w} {svg_h}" '
             f'width="{svg_w}" height="{svg_h}" '
             f'style="max-width:100%;height:auto;">')

    # Defs: arrowhead markers per framework colour
    S.append("<defs>")
    for fi, fw_name in enumerate(fw_names):
        colour = _FW_COLORS[fi % len(_FW_COLORS)][0]
        mid = f"arrow-{fi}"
        S.append(
            f'<marker id="{mid}" markerWidth="8" markerHeight="6" '
            f'refX="8" refY="3" orient="auto" markerUnits="strokeWidth">'
            f'<path d="M0,0 L8,3 L0,6 Z" fill="{colour}"/></marker>')
    S.append("</defs>")

    # Pannable / zoomable canvas group
    S.append('<g class="diagram-canvas">')

    # ── Reaction arrows (behind nodes) ──
    # Pre-count arrows between each pair for offset calculation
    pair_counts: Dict[tuple, int] = {}
    for r in all_reactions:
        c, p = r["consumed"], r["produced"]
        if c == "NONE" or p == "NONE":
            continue
        key = (min(c, p), max(c, p))
        pair_counts[key] = pair_counts.get(key, 0) + 1

    pair_idx: Dict[tuple, int] = {}
    label_candidates: List[dict] = []  # collect for collision resolution
    current_fw = None

    for r in all_reactions:
        fw = r["_fw"]
        fi = fw_names.index(fw)
        colour = _FW_COLORS[fi % len(_FW_COLORS)][0]
        c, p = r["consumed"], r["produced"]

        if c == "NONE" and p == "NONE":
            continue

        # Open / switch framework group wrapper
        if fw != current_fw:
            if current_fw is not None:
                S.append("</g>")
            S.append(f'<g class="fw-group" data-fw="{_svg_escape(fw)}">')
            current_fw = fw

        # ── Source arrows (NONE → species) ──
        if c == "NONE" and p in positions:
            px, py = positions[p]
            sx, sy = px - 60, py - 30
            qsx, qsy = (sx + px) / 2, sy
            d_str = (f"M{sx:.0f},{sy:.0f} Q{qsx:.0f},{qsy:.0f} "
                     f"{px:.0f},{py - node_h / 2:.0f}")
            kin_tip = _build_kinetics_tip(r)
            kin_attr = f' data-kinetics="{_svg_escape(kin_tip)}"' if kin_tip else ""
            _src_obs_attr = ""
            if species_obs_availability is not None:
                p_obs = species_obs_availability.get(p, {}).get("has_obs", False)
                _src_obs_attr = f' data-has-obs="{"true" if p_obs else "false"}"'
            S.append(
                f'<g class="reaction-arrow source-sink-arrow" '
                f'data-fw="{_svg_escape(fw)}" '
                f'data-rxn="{_svg_escape(r["name"])}" '
                f'data-consumed="NONE" data-produced="{_svg_escape(p)}"'
                f'{kin_attr}{_src_obs_attr}>'
                f'<path d="{d_str}" fill="none" stroke="transparent" '
                f'stroke-width="14" class="hit-area"/>'
                f'<path d="{d_str}" fill="none" stroke="{colour}" '
                f'stroke-width="1.5" stroke-dasharray="5 3" '
                f'marker-end="url(#arrow-{fi})" class="vis-path"/>')
            # Collect label
            lab_text = r["name"].replace("_", " ")
            if len(lab_text) > 22:
                lab_text = lab_text[:20] + "\u2026"
            dx_s = px - sx
            dy_s = (py - node_h / 2) - sy
            ln = math.sqrt(dx_s * dx_s + dy_s * dy_s) or 1
            label_candidates.append({
                "text": lab_text, "colour": colour, "italic": True,
                "fw": fw, "rxn": r["name"],
                "x": qsx - 10, "y": sy - 4,
                "w": len(lab_text) * _LABEL_CHAR_W + 8, "h": 12,
                "ax": sx, "ay": sy,
                "qx": qsx, "qy": qsy,
                "bx": px, "by": py - node_h / 2,
                "nx": -dy_s / ln, "ny": dx_s / ln,
            })
            S.append("</g>")
            continue

        # ── Sink arrows (species → NONE) ──
        if p == "NONE" and c in positions:
            cx_p, cy_p = positions[c]
            ex, ey = cx_p + 60, cy_p + 30
            qex, qey = (cx_p + ex) / 2, ey
            d_str = (f"M{cx_p:.0f},{cy_p + node_h / 2:.0f} "
                     f"Q{qex:.0f},{qey:.0f} {ex:.0f},{ey:.0f}")
            kin_tip = _build_kinetics_tip(r)
            kin_attr = f' data-kinetics="{_svg_escape(kin_tip)}"' if kin_tip else ""
            _snk_obs_attr = ""
            if species_obs_availability is not None:
                c_obs = species_obs_availability.get(c, {}).get("has_obs", False)
                _snk_obs_attr = f' data-has-obs="{"true" if c_obs else "false"}"'
            S.append(
                f'<g class="reaction-arrow source-sink-arrow" '
                f'data-fw="{_svg_escape(fw)}" '
                f'data-rxn="{_svg_escape(r["name"])}" '
                f'data-consumed="{_svg_escape(c)}" data-produced="NONE"'
                f'{kin_attr}{_snk_obs_attr}>'
                f'<path d="{d_str}" fill="none" stroke="transparent" '
                f'stroke-width="14" class="hit-area"/>'
                f'<path d="{d_str}" fill="none" stroke="{colour}" '
                f'stroke-width="1.5" stroke-dasharray="5 3" '
                f'marker-end="url(#arrow-{fi})" class="vis-path"/>')
            lab_text = r["name"].replace("_", " ")
            if len(lab_text) > 22:
                lab_text = lab_text[:20] + "\u2026"
            dx_s = ex - cx_p
            dy_s = ey - (cy_p + node_h / 2)
            ln = math.sqrt(dx_s * dx_s + dy_s * dy_s) or 1
            label_candidates.append({
                "text": lab_text, "colour": colour, "italic": True,
                "fw": fw, "rxn": r["name"],
                "x": qex + 4, "y": ey + 12,
                "w": len(lab_text) * _LABEL_CHAR_W + 8, "h": 12,
                "ax": cx_p, "ay": cy_p + node_h / 2,
                "qx": qex, "qy": qey,
                "bx": ex, "by": ey,
                "nx": -dy_s / ln, "ny": dx_s / ln,
            })
            S.append("</g>")
            continue

        # ── Regular reaction arrows (species → species) ──
        if c not in positions or p not in positions:
            continue

        x1, y1 = positions[c]
        x2, y2 = positions[p]

        key = (min(c, p), max(c, p))
        idx = pair_idx.get(key, 0)
        pair_idx[key] = idx + 1
        total = pair_counts.get(key, 1)
        offset = (idx - (total - 1) / 2) * 18

        mx, my = (x1 + x2) / 2, (y1 + y2) / 2
        dx, dy = x2 - x1, y2 - y1
        length = math.sqrt(dx * dx + dy * dy) or 1
        nx, ny = -dy / length, dx / length
        curve = length * 0.18 + offset

        # Shorten to avoid node overlap
        node_w_c = max(len(c) * char_w + 20, 60)
        node_w_p = max(len(p) * char_w + 20, 60)
        t_start = max((node_w_c / 2 + 4) / length, 0.05)
        t_end = max((node_w_p / 2 + 4) / length, 0.05)
        ax_ = x1 + dx * t_start
        ay_ = y1 + dy * t_start
        bx_ = x2 - dx * t_end
        by_ = y2 - dy * t_end

        smx, smy = (ax_ + bx_) / 2, (ay_ + by_) / 2
        qx2 = smx + nx * curve
        qy2 = smy + ny * curve

        d_str = f"M{ax_:.1f},{ay_:.1f} Q{qx2:.1f},{qy2:.1f} {bx_:.1f},{by_:.1f}"

        kin_tip = _build_kinetics_tip(r)
        kin_attr = f' data-kinetics="{_svg_escape(kin_tip)}"' if kin_tip else ""
        # Obs availability for this arrow
        _arrow_obs_attr = ""
        _arrow_dash = ""
        if species_obs_availability is not None:
            c_obs = species_obs_availability.get(c, {}).get("has_obs", False)
            p_obs = species_obs_availability.get(p, {}).get("has_obs", False)
            _arrow_has = c_obs or p_obs
            _arrow_obs_attr = f' data-has-obs="{"true" if _arrow_has else "false"}"'
            if not _arrow_has:
                _arrow_dash = " stroke-dasharray='6 3'"
        S.append(
            f'<g class="reaction-arrow" '
            f'data-fw="{_svg_escape(fw)}" '
            f'data-rxn="{_svg_escape(r["name"])}" '
            f'data-consumed="{_svg_escape(c)}" '
            f'data-produced="{_svg_escape(p)}"'
            f'{kin_attr}{_arrow_obs_attr}>'
            f'<path d="{d_str}" fill="none" stroke="transparent" '
            f'stroke-width="14" class="hit-area"/>'
            f'<path d="{d_str}" fill="none" stroke="{colour}" '
            f'stroke-width="1.8" marker-end="url(#arrow-{fi})"{_arrow_dash} '
            f'class="vis-path"/>')

        # Collect label candidate
        lab_text = r["name"].replace("_", " ")
        if len(lab_text) > 22:
            lab_text = lab_text[:20] + "\u2026"
        lx_init = (ax_ + 2 * qx2 + bx_) / 4
        ly_init = (ay_ + 2 * qy2 + by_) / 4
        label_candidates.append({
            "text": lab_text, "colour": colour, "italic": False,
            "fw": fw, "rxn": r["name"],
            "x": lx_init, "y": ly_init - 4,
            "w": len(lab_text) * _LABEL_CHAR_W + 8, "h": 12,
            "ax": ax_, "ay": ay_,
            "qx": qx2, "qy": qy2,
            "bx": bx_, "by": by_,
            "nx": nx, "ny": ny,
        })
        S.append("</g>")

    # Close last framework group
    if current_fw is not None:
        S.append("</g>")

    # ── Resolve label overlaps and emit labels ──
    _resolve_label_positions(label_candidates, placed_boxes)

    S.append('<g class="label-layer">')
    for lab in label_candidates:
        lx, ly = lab["x"], lab["y"]
        w = lab["w"]
        h = lab["h"]
        style_extra = "font-style:italic;" if lab.get("italic") else ""
        fw_attr = f' data-fw="{_svg_escape(lab["fw"])}"' if lab.get("fw") else ""
        rxn_attr = f' data-rxn="{_svg_escape(lab["rxn"])}"' if lab.get("rxn") else ""
        # Background rect for readability
        S.append(
            f'<rect class="label-bg"{fw_attr}{rxn_attr} x="{lx - w / 2 - 1:.1f}" '
            f'y="{ly - h / 2 - 1:.1f}" '
            f'width="{w + 2:.0f}" height="{h + 2:.0f}" rx="3" ry="3" '
            f'fill="var(--surface,#fff)" fill-opacity="0.82" stroke="none"/>')
        S.append(
            f'<text{fw_attr}{rxn_attr} x="{lx:.1f}" y="{ly + 1:.1f}" '
            f'text-anchor="middle" dominant-baseline="central" '
            f'font-size="{_LABEL_FONT}" fill="{lab["colour"]}" '
            f'style="font-weight:500;pointer-events:none;{style_extra}">'
            f'{_svg_escape(lab["text"])}</text>')
    S.append("</g>")

    # ── Build species → framework-index mapping ──
    species_to_fws: dict[str, list[int]] = {}
    for fi, fw_name in enumerate(fw_names):
        fw_info = network["frameworks"].get(fw_name, {})
        for r in fw_info.get("reactions", []):
            for sp_ref in (r.get("consumed", "NONE"), r.get("produced", "NONE")):
                if sp_ref != "NONE" and sp_ref in positions:
                    species_to_fws.setdefault(sp_ref, [])
                    if fi not in species_to_fws[sp_ref]:
                        species_to_fws[sp_ref].append(fi)

    # ── Species nodes (on top) ──
    for sp in all_species:
        if sp not in positions:
            continue
        x, y = positions[sp]
        label = sp
        desc = species_desc.get(sp, "")
        desc_attr = f' data-desc="{_svg_escape(desc)}"' if desc else ""
        # Observation availability (calibration report only)
        obs_attr = ""
        if species_obs_availability is not None:
            has_obs = species_obs_availability.get(sp, {}).get("has_obs", False)
            obs_attr = f' data-has-obs="{"true" if has_obs else "false"}"'
        node_w = max(len(label) * char_w + 24, 64)
        is_loss = "loss" in sp.lower() or "atm" in sp.lower()
        cls = "species-node sink" if is_loss else "species-node"
        # Dashed border for species without observation data (preserves colour)
        no_obs_style = ""
        if species_obs_availability is not None:
            has_obs = species_obs_availability.get(sp, {}).get("has_obs", False)
            if not has_obs:
                no_obs_style = "stroke-dasharray:5 3;"

        # Framework colours for this species
        sp_fw_indices = species_to_fws.get(sp, [])
        if sp_fw_indices:
            sp_colours = [
                _FW_COLORS[fi % len(_FW_COLORS)][0] for fi in sp_fw_indices
            ]
        else:
            sp_colours = ["var(--primary,#4d9ee8)"]
        n_layers = len(sp_colours)

        S.append(
            f'<g class="{cls}" data-species="{_svg_escape(sp)}"{desc_attr}{obs_attr}>')
        # Emit nested rects: outermost first, innermost last (gets fill)
        _inset = 3  # px per nesting level
        for lvl, col in enumerate(sp_colours):
            ix = _inset * lvl
            rw = node_w - 2 * ix
            rh = node_h - 2 * ix
            rx_val = max(10 - lvl * 2, 4)
            is_innermost = (lvl == n_layers - 1)
            fill = "var(--surface,#fff)" if is_innermost else "none"
            dash = ""
            if is_loss:
                dash += " stroke-dasharray:4 2;"
            if lvl == 0 and no_obs_style:
                dash += no_obs_style
            S.append(
                f'<rect x="{x - rw / 2:.1f}" y="{y - rh / 2:.1f}" '
                f'width="{rw:.0f}" height="{rh}" rx="{rx_val}" ry="{rx_val}" '
                f'style="fill:{fill};stroke:{col};stroke-width:2;{dash}"/>')
        S.append(
            f'<text x="{x:.1f}" y="{y:.1f}" text-anchor="middle" '
            f'dominant-baseline="central" font-size="11.5" font-weight="600" '
            f'style="fill:var(--text,#222);">'
            f'{_svg_escape(label)}</text></g>')

    S.append("</g>")  # close diagram-canvas

    # Selection rectangle for area-zoom (hidden by default)
    S.append('<rect class="zoom-sel" x="0" y="0" width="0" height="0" '
             'fill="rgba(77,158,232,.12)" stroke="var(--primary,#4d9ee8)" '
             'stroke-width="1" stroke-dasharray="4 2" '
             'style="display:none;pointer-events:none;"/>')

    S.append("</svg>")
    return "\n".join(S)


# ── Interactivity JavaScript ─────────────────────────────────────────

def _bgc_diagram_js() -> str:
    """Return a <script> block for click-to-highlight, framework toggles, hover tooltips, and zoom."""
    return """<script>
(function(){
  'use strict';
  var activeEl = null;   /* currently highlighted arrow or species node */
  var activeKind = null; /* 'arrow' | 'species' | null */

  /* ── Tooltip element (shared, created once) ───── */
  var tip = document.createElement('div');
  tip.className = 'bgc-tooltip';
  tip.style.display = 'none';
  document.body.appendChild(tip);

  function showTip(el, evt){
    var rxn = el.getAttribute('data-rxn') || '';
    var consumed = el.getAttribute('data-consumed') || '';
    var produced = el.getAttribute('data-produced') || '';
    var kin = el.getAttribute('data-kinetics') || '';

    var lines = [];
    lines.push('<strong>' + rxn.replace(/_/g,' ') + '</strong>');
    if(consumed!=='NONE' || produced!=='NONE'){
      lines.push('<span class="bgc-tip-reaction">' +
        (consumed!=='NONE'?consumed:'(source)') + ' \\u2192 ' +
        (produced!=='NONE'?produced:'(sink)') + '</span>');
    }
    if(kin){
      var raw = el.getAttribute('data-kinetics');
      var segs = raw.split(' &#10; ');
      for(var i=0;i<segs.length;i++){
        var s=segs[i].trim();
        if(!s) continue;
        if(s.indexOf('Parameters:')===0){
          lines.push('<span class="bgc-tip-params">' + s + '</span>');
        } else {
          lines.push('<span class="bgc-tip-expr">' + s + '</span>');
        }
      }
    }
    if(!kin && consumed==='NONE' && produced==='NONE') return;

    tip.innerHTML = lines.join('');
    tip.style.display = 'block';
    positionTip(evt);
  }

  function positionTip(evt){
    var x = evt.clientX + 14;
    var y = evt.clientY + 14;
    var tw = tip.offsetWidth;
    var th = tip.offsetHeight;
    if(x + tw > window.innerWidth - 8) x = evt.clientX - tw - 10;
    if(y + th > window.innerHeight - 8) y = evt.clientY - th - 10;
    tip.style.left = x + 'px';
    tip.style.top  = y + 'px';
  }

  function hideTip(){ tip.style.display='none'; }

  function showSpeciesTip(el, evt){
    var sp = el.getAttribute('data-species') || '';
    var desc = el.getAttribute('data-desc') || '';
    var hasObs = el.getAttribute('data-has-obs');
    if(!desc && hasObs === null) return;
    var html = '<strong>' + sp + '</strong>';
    if(desc) html += '<span class="bgc-tip-reaction">' + desc + '</span>';
    if(hasObs === 'false')
      html += '<span class="bgc-tip-reaction" style="color:#c47800;font-style:italic;">No observation data for calibration</span>';
    tip.innerHTML = html;
    tip.style.display = 'block';
    positionTip(evt);
  }

  /* Hover events */
  document.addEventListener('mouseover', function(e){
    var arrow = e.target.closest('.reaction-arrow');
    if(arrow){ showTip(arrow, e); return; }
    var species = e.target.closest('.species-node');
    if(species) showSpeciesTip(species, e);
  });
  document.addEventListener('mousemove', function(e){
    if(tip.style.display==='block') positionTip(e);
  });
  document.addEventListener('mouseout', function(e){
    var arrow = e.target.closest('.reaction-arrow');
    var species = e.target.closest('.species-node');
    if(arrow || species) hideTip();
  });

  /* ── Dim / undim helpers ──────────────────────── */
  function dimAll(svg){
    svg.querySelectorAll('.reaction-arrow').forEach(function(g){
      g.style.opacity='0.12'; g.style.filter='grayscale(100%)';
    });
    svg.querySelectorAll('.species-node').forEach(function(g){
      g.style.opacity='0.18'; g.style.filter='grayscale(100%)';
    });
    svg.querySelectorAll('.label-layer text').forEach(function(t){
      t.style.opacity='0.12';
    });
    svg.querySelectorAll('.label-layer .label-bg').forEach(function(r){
      r.style.opacity='0.12';
    });
  }

  function undimEl(el){ el.style.opacity='1'; el.style.filter=''; }

  function undimLabel(svg, rxnName){
    if(!rxnName) return;
    svg.querySelectorAll('.label-layer [data-rxn="'+rxnName+'"]').forEach(function(el){
      el.style.opacity='1';
    });
  }

  /* ── Click-to-highlight ────────────────────────── */
  function resetHighlight(svg){
    svg.querySelectorAll('.reaction-arrow').forEach(function(g){
      g.style.opacity=''; g.style.filter='';
    });
    svg.querySelectorAll('.species-node').forEach(function(g){
      g.style.opacity=''; g.style.filter='';
    });
    svg.querySelectorAll('.label-layer text').forEach(function(t){
      t.style.opacity='';
    });
    svg.querySelectorAll('.label-layer .label-bg').forEach(function(r){
      r.style.opacity='';
    });
    activeEl = null; activeKind = null;
    /* Reset zoom to full view */
    var container = svg.closest('.bgc-diagram-container');
    if(container && container._zoomState){
      var z = container._zoomState;
      z.curVB.x = z.origVB.x; z.curVB.y = z.origVB.y;
      z.curVB.w = z.origVB.w; z.curVB.h = z.origVB.h;
      z.applyVB();
    }
  }

  function highlightReaction(svg, el){
    var consumed = el.getAttribute('data-consumed');
    var produced = el.getAttribute('data-produced');
    var rxnName  = el.getAttribute('data-rxn');

    dimAll(svg);
    undimEl(el);
    undimLabel(svg, rxnName);

    /* Highlight consumed + produced species */
    var keep = {};
    if(consumed && consumed!=='NONE') keep[consumed]=true;
    if(produced && produced!=='NONE') keep[produced]=true;
    svg.querySelectorAll('.species-node').forEach(function(g){
      if(keep[g.getAttribute('data-species')]) undimEl(g);
    });
    activeEl = el; activeKind = 'arrow';
    zoomToHighlighted(svg);
  }

  /* ── NEW: click a species node ─────────────────── */
  function highlightSpecies(svg, speciesNode){
    var sp = speciesNode.getAttribute('data-species');
    dimAll(svg);
    undimEl(speciesNode);

    /* Find all arrows that consume or produce this species */
    var connectedSpecies = {};
    connectedSpecies[sp] = true;
    svg.querySelectorAll('.reaction-arrow').forEach(function(g){
      var c = g.getAttribute('data-consumed');
      var p = g.getAttribute('data-produced');
      if(c === sp || p === sp){
        undimEl(g);
        undimLabel(svg, g.getAttribute('data-rxn'));
        /* Also highlight the other endpoint species */
        if(c && c!=='NONE') connectedSpecies[c] = true;
        if(p && p!=='NONE') connectedSpecies[p] = true;
      }
    });
    /* Highlight all connected species nodes */
    svg.querySelectorAll('.species-node').forEach(function(g){
      if(connectedSpecies[g.getAttribute('data-species')]) undimEl(g);
    });
    activeEl = speciesNode; activeKind = 'species';
    zoomToHighlighted(svg);
  }

  /* ── Auto zoom-fit to highlighted elements ─────── */
  function zoomToHighlighted(svg){
    var container = svg.closest('.bgc-diagram-container');
    if(!container || !container._zoomState) return;
    var z = container._zoomState;

    /* Collect bounding boxes of all visible (non-dimmed) elements */
    var minX = Infinity, minY = Infinity, maxX = -Infinity, maxY = -Infinity;
    var found = false;

    function expandBB(el){
      try {
        var bb = el.getBBox();
        if(bb.width === 0 && bb.height === 0) return;
        /* Transform to SVG root coordinates */
        var ctm = el.getCTM();
        var svgCtm = svg.getCTM();
        if(!ctm || !svgCtm) return;
        /* Get inverse of SVG CTM to go from screen to SVG coords */
        var inv = svgCtm.inverse();
        var m = inv.multiply(ctm);
        /* Transform all four corners */
        var pts = [
          {x: bb.x, y: bb.y},
          {x: bb.x + bb.width, y: bb.y},
          {x: bb.x, y: bb.y + bb.height},
          {x: bb.x + bb.width, y: bb.y + bb.height}
        ];
        for(var i=0;i<4;i++){
          var px = m.a * pts[i].x + m.c * pts[i].y + m.e;
          var py = m.b * pts[i].x + m.d * pts[i].y + m.f;
          if(px < minX) minX = px;
          if(py < minY) minY = py;
          if(px > maxX) maxX = px;
          if(py > maxY) maxY = py;
        }
        found = true;
      } catch(e){}
    }

    svg.querySelectorAll('.reaction-arrow').forEach(function(g){
      if(g.style.opacity === '' || g.style.opacity === '1') expandBB(g);
    });
    svg.querySelectorAll('.species-node').forEach(function(g){
      if(g.style.opacity === '' || g.style.opacity === '1') expandBB(g);
    });

    if(!found) return;

    /* Add padding (15% on each side) */
    var padX = (maxX - minX) * 0.15;
    var padY = (maxY - minY) * 0.15;
    minX -= padX; minY -= padY; maxX += padX; maxY += padY;
    var w = maxX - minX;
    var h = maxY - minY;
    if(w < 10 || h < 10) return;

    /* Preserve aspect ratio */
    var aspect = z.origVB.w / z.origVB.h;
    if(w / h > aspect){ h = w / aspect; }
    else { w = h * aspect; }
    /* Center the box */
    var cx = (minX + maxX) / 2;
    var cy = (minY + maxY) / 2;

    /* Don't zoom out beyond original */
    if(w > z.origVB.w){ w = z.origVB.w; h = z.origVB.h; cx = z.origVB.x + z.origVB.w/2; cy = z.origVB.y + z.origVB.h/2; }

    z.curVB.x = cx - w/2;
    z.curVB.y = cy - h/2;
    z.curVB.w = w;
    z.curVB.h = h;
    z.applyVB();
  }

  /* ── Click handler (arrows, species, background) ── */
  document.addEventListener('click', function(e){
    /* Ignore clicks on toggle chips, zoom buttons etc */
    if(e.target.closest('.fw-toggle-chip') || e.target.closest('.zoom-btn')) return;
    var svg = e.target.closest('.bgc-diagram-container svg');
    if(!svg) return;

    /* Check if this was a drag (rect-zoom) — ignore click after drag */
    var container = svg.closest('.bgc-diagram-container');
    if(container && container._wasDrag){ container._wasDrag = false; return; }

    var arrow = e.target.closest('.reaction-arrow');
    var species = e.target.closest('.species-node');

    if(arrow){
      if(arrow===activeEl && activeKind==='arrow') resetHighlight(svg);
      else highlightReaction(svg, arrow);
    } else if(species){
      if(species===activeEl && activeKind==='species') resetHighlight(svg);
      else highlightSpecies(svg, species);
    } else {
      if(activeEl) resetHighlight(svg);
    }
  });

  /* ── Update species visibility based on active frameworks ── */
  function updateSpeciesVisibility(svg, container){
    /* Collect active framework names */
    var activeFW = {};
    container.querySelectorAll('.fw-toggle-chip.active').forEach(function(c){
      var f = c.getAttribute('data-fw');
      if(f && f!=='__all__') activeFW[f] = true;
    });
    /* Collect species used by active frameworks */
    var usedSpecies = {};
    svg.querySelectorAll('.reaction-arrow').forEach(function(g){
      var f = g.getAttribute('data-fw');
      if(!activeFW[f]) return;
      var c = g.getAttribute('data-consumed');
      var p = g.getAttribute('data-produced');
      if(c && c!=='NONE') usedSpecies[c] = true;
      if(p && p!=='NONE') usedSpecies[p] = true;
    });
    /* Gray out species not used by any active framework */
    svg.querySelectorAll('.species-node').forEach(function(g){
      var sp = g.getAttribute('data-species');
      if(usedSpecies[sp]){
        g.style.opacity=''; g.style.filter='';
      } else {
        g.style.opacity='0.18'; g.style.filter='grayscale(100%)';
      }
    });
  }

  /* ── Framework toggle chips (checkbox style) ───── */
  document.addEventListener('click', function(e){
    var chip = e.target.closest('.fw-toggle-chip');
    if(!chip) return;
    var container = chip.closest('.bgc-diagram-container');
    if(!container) return;
    var svg = container.querySelector('svg');
    if(!svg) return;
    var fw = chip.getAttribute('data-fw');

    /* "All" button */
    if(fw==='__all__'){
      container.querySelectorAll('.fw-toggle-chip').forEach(function(c){
        c.classList.add('active');
      });
      svg.querySelectorAll('.fw-group').forEach(function(g){
        g.style.display='';
      });
      svg.querySelectorAll('.label-layer text, .label-layer .label-bg').forEach(function(el){
        el.style.display='';
      });
      svg.querySelectorAll('.species-node').forEach(function(g){
        g.style.opacity=''; g.style.filter='';
      });
      if(activeEl) resetHighlight(svg);
      return;
    }

    /* Toggle individual framework */
    chip.classList.toggle('active');
    var isActive = chip.classList.contains('active');
    svg.querySelectorAll('.fw-group[data-fw="'+fw+'"]').forEach(function(g){
      g.style.display = isActive ? '' : 'none';
    });
    /* Also hide/show corresponding labels */
    svg.querySelectorAll('.label-layer [data-fw="'+fw+'"]').forEach(function(el){
      el.style.display = isActive ? '' : 'none';
    });
    /* Gray out species no longer used by any active framework */
    updateSpeciesVisibility(svg, container);

    if(activeEl && activeKind==='arrow' && activeEl.getAttribute('data-fw')===fw && !isActive){
      resetHighlight(svg);
    }
  });

  /* ── Zoom & Pan ────────────────────────────────── */
  document.querySelectorAll('.bgc-diagram-container').forEach(function(container){
    var svg = container.querySelector('svg');
    if(!svg) return;
    var canvas = svg.querySelector('.diagram-canvas');
    var selRect = svg.querySelector('.zoom-sel');
    if(!canvas) return;

    var vb = svg.viewBox.baseVal;
    var origVB = {x:vb.x, y:vb.y, w:vb.width, h:vb.height};
    var curVB  = {x:vb.x, y:vb.y, w:vb.width, h:vb.height};
    var zoomLocked = true;  /* +/- and scroll zoom locked by default */

    function applyVB(){
      svg.setAttribute('viewBox',
        curVB.x+' '+curVB.y+' '+curVB.w+' '+curVB.h);
    }

    /* Expose zoom state on container so highlight functions can access it */
    container._zoomState = { origVB: origVB, curVB: curVB, applyVB: applyVB };

    /* Mouse-to-SVG coordinate conversion */
    function svgPt(evt){
      var r = svg.getBoundingClientRect();
      return {
        x: curVB.x + (evt.clientX - r.left) / r.width  * curVB.w,
        y: curVB.y + (evt.clientY - r.top)  / r.height * curVB.h
      };
    }

    /* Wheel zoom (disabled when locked) */
    svg.addEventListener('wheel', function(e){
      e.preventDefault();
      if(zoomLocked) return;
      var pt = svgPt(e);
      var factor = e.deltaY > 0 ? 1.15 : 1/1.15;
      var nw = curVB.w * factor;
      var nh = curVB.h * factor;
      if(nw > origVB.w * 1.05){ nw=origVB.w; nh=origVB.h; }
      if(nw < origVB.w * 0.1){ return; }
      curVB.x = pt.x - (pt.x - curVB.x) * (nw / curVB.w);
      curVB.y = pt.y - (pt.y - curVB.y) * (nh / curVB.h);
      curVB.w = nw;
      curVB.h = nh;
      applyVB();
    }, {passive:false});

    /* Rectangle-select zoom — primary drag action (left-click drag) */
    var dragStart = null;
    var dragStartScreen = null;
    var dragging = false;

    svg.addEventListener('mousedown', function(e){
      if(e.button !== 0) return; /* left button only */
      e.preventDefault();
      dragStart = svgPt(e);
      dragStartScreen = {x: e.clientX, y: e.clientY};
      dragging = false; /* becomes true once movement threshold is met */
    });

    svg.addEventListener('mousemove', function(e){
      if(!dragStart) return;
      /* Threshold: start rectangle only after 5px of movement */
      if(!dragging){
        var dx = e.clientX - dragStartScreen.x;
        var dy = e.clientY - dragStartScreen.y;
        if(Math.sqrt(dx*dx + dy*dy) < 5) return;
        dragging = true;
        selRect.style.display = '';
      }
      var cur = svgPt(e);
      var x = Math.min(dragStart.x, cur.x);
      var y = Math.min(dragStart.y, cur.y);
      var w = Math.abs(cur.x - dragStart.x);
      var h = Math.abs(cur.y - dragStart.y);
      selRect.setAttribute('x', x);
      selRect.setAttribute('y', y);
      selRect.setAttribute('width', w);
      selRect.setAttribute('height', h);
    });

    svg.addEventListener('mouseup', function(e){
      if(!dragStart) return;
      var wasDrag = dragging;
      selRect.style.display = 'none';
      selRect.setAttribute('width', 0);
      selRect.setAttribute('height', 0);

      if(wasDrag){
        container._wasDrag = true; /* signal click handler to ignore this */
        var end = svgPt(e);
        var x = Math.min(dragStart.x, end.x);
        var y = Math.min(dragStart.y, end.y);
        var w = Math.abs(end.x - dragStart.x);
        var h = Math.abs(end.y - dragStart.y);
        dragStart = null; dragStartScreen = null; dragging = false;
        if(w < origVB.w * 0.02 || h < origVB.h * 0.02) return;
        /* Preserve aspect ratio */
        var aspect = origVB.w / origVB.h;
        if(w / h > aspect){ h = w / aspect; }
        else { w = h * aspect; }
        curVB.x = x; curVB.y = y; curVB.w = w; curVB.h = h;
        applyVB();
      }
      dragStart = null; dragStartScreen = null; dragging = false;
    });

    /* Cancel drag if mouse leaves SVG */
    svg.addEventListener('mouseleave', function(){
      if(dragStart){
        selRect.style.display = 'none';
        selRect.setAttribute('width', 0);
        selRect.setAttribute('height', 0);
        dragStart = null; dragStartScreen = null; dragging = false;
      }
    });

    /* Zoom control buttons */
    container.querySelectorAll('.zoom-btn').forEach(function(btn){
      btn.addEventListener('click', function(e){
        e.stopPropagation();
        var action = btn.getAttribute('data-zoom');
        if(action === 'lock'){
          zoomLocked = !zoomLocked;
          btn.classList.toggle('locked', zoomLocked);
          btn.innerHTML = zoomLocked ? '&#x1f512;' : '&#x1f513;';
          btn.title = zoomLocked ? 'Unlock +/\u2212 zoom' : 'Lock +/\u2212 zoom';
          /* Enable/disable +/- buttons */
          container.querySelectorAll('.zoom-btn[data-zoom="in"],.zoom-btn[data-zoom="out"]').forEach(function(b){
            if(zoomLocked) b.setAttribute('disabled','');
            else b.removeAttribute('disabled');
          });
          /* Update hint text */
          var hint = container.querySelector('.zoom-hint');
          if(hint) hint.textContent = zoomLocked
            ? 'Drag to select area'
            : 'Scroll to zoom \u00b7 Drag to select area';
          return;
        }
        if(action === 'in'){
          var cx = curVB.x + curVB.w/2;
          var cy = curVB.y + curVB.h/2;
          curVB.w *= 0.7; curVB.h *= 0.7;
          curVB.x = cx - curVB.w/2;
          curVB.y = cy - curVB.h/2;
        } else if(action === 'out'){
          var cx = curVB.x + curVB.w/2;
          var cy = curVB.y + curVB.h/2;
          curVB.w = Math.min(curVB.w * 1.4, origVB.w);
          curVB.h = Math.min(curVB.h * 1.4, origVB.h);
          curVB.x = cx - curVB.w/2;
          curVB.y = cy - curVB.h/2;
          if(curVB.w >= origVB.w * 0.99){
            curVB.x = origVB.x; curVB.y = origVB.y;
            curVB.w = origVB.w; curVB.h = origVB.h;
          }
        } else if(action === 'reset'){
          curVB.x = origVB.x; curVB.y = origVB.y;
          curVB.w = origVB.w; curVB.h = origVB.h;
        }
        applyVB();
      });
    });

    /* ── "Only species with observation data" toggle ── */
    var obsCb = container.querySelector('.obs-only-cb');
    if(obsCb){
      obsCb.addEventListener('change', function(){
        var obsOnly = obsCb.checked;

        /* Gray out species without obs (keep visible & interactive) */
        svg.querySelectorAll('.species-node').forEach(function(g){
          var hasObs = g.getAttribute('data-has-obs');
          if(obsOnly && hasObs === 'false'){
            g.style.opacity = '0.25';
            g.style.filter = 'grayscale(100%)';
          } else {
            g.style.opacity = '';
            g.style.filter = '';
          }
        });

        /* Gray out arrows without obs (keep visible & interactive) */
        svg.querySelectorAll('.reaction-arrow').forEach(function(g){
          var hasObs = g.getAttribute('data-has-obs');
          if(obsOnly && hasObs === 'false'){
            g.style.opacity = '0.2';
            g.style.filter = 'grayscale(100%)';
          } else {
            g.style.opacity = '';
            g.style.filter = '';
          }
        });

        /* Gray out labels of no-obs arrows */
        svg.querySelectorAll('.label-layer [data-rxn]').forEach(function(el){
          var rxn = el.getAttribute('data-rxn');
          var arrow = svg.querySelector('.reaction-arrow[data-rxn="'+rxn+'"]');
          if(arrow){
            var hasObs = arrow.getAttribute('data-has-obs');
            if(obsOnly && hasObs === 'false'){
              el.style.opacity = '0.2';
            } else {
              el.style.opacity = '';
            }
          }
        });

        /* Gray out framework chips that have no obs species */
        container.querySelectorAll('.fw-toggle-chip').forEach(function(chip){
          var fw = chip.getAttribute('data-fw');
          if(fw === '__all__') return;
          if(obsOnly){
            var hasObsArrow = false;
            svg.querySelectorAll('.reaction-arrow[data-fw="'+fw+'"]').forEach(function(a){
              if(a.getAttribute('data-has-obs') === 'true') hasObsArrow = true;
            });
            if(!hasObsArrow){
              chip.style.opacity = '0.35';
            }
          } else {
            chip.style.opacity = '';
          }
        });

        /* Dispatch custom event for calibration report integration */
        container.dispatchEvent(new CustomEvent('obs-filter-changed', {
          detail: { obsOnly: obsOnly },
          bubbles: true
        }));
      });
    }
  });
})();
</script>"""


# ── Main entry point ─────────────────────────────────────────────────

def generate_bgc_reaction_diagram(
    bgc_data: Optional[dict] = None,
    bgc_template_path: Optional[str] = None,
    pqi_filepath: Optional[str] = None,
    module_name: str = "NATIVE_BGC_FLEX",
    max_width: int = 860,
    species_obs_availability: Optional[dict] = None,
) -> str:
    """Generate an HTML/SVG reaction network diagram for BGC modules.

    Accepts either pre-parsed *bgc_data*, a path to a JSON template, or
    a path to a PHREEQC .pqi file.  Returns an HTML ``<div>`` containing
    the SVG diagram with interactive framework toggles and click-to-
    highlight, or ``""`` if no data is available.

    *species_obs_availability* (calibration report only): dict mapping
    species name -> {"has_obs": bool, ...}.  When provided, species and
    frameworks without observation data get a visual annotation.
    """
    network = {}

    # 1. Load data
    if bgc_data is not None:
        network = extract_bgc_network(bgc_data)
    elif bgc_template_path and os.path.isfile(bgc_template_path):
        try:
            with open(bgc_template_path, "r", errors="replace") as fh:
                data = json.load(fh)
            network = extract_bgc_network(data)
        except Exception:
            return ""
    elif pqi_filepath and os.path.isfile(pqi_filepath):
        network = extract_bgc_network_from_pqi(pqi_filepath)
    elif module_name.upper() == "PHREEQC" and pqi_filepath:
        network = extract_bgc_network_from_pqi(pqi_filepath)

    if not network or not network.get("frameworks"):
        return ""

    # 2. Render SVG
    svg_html = _render_bgc_svg(
        network, max_width=max_width,
        species_obs_availability=species_obs_availability,
    )
    if not svg_html:
        return ""

    # 3. Framework toggle bar (above diagram) — checkbox-style chips
    fw_names = list(network["frameworks"].keys())
    toggle_items = [
        '<button class="fw-toggle-chip active" data-fw="__all__">'
        '<span class="chip-tick">\u2714</span> All</button>'
    ]
    for fi, fw_name in enumerate(fw_names):
        colour = _FW_COLORS[fi % len(_FW_COLORS)][0]
        fw_info = network["frameworks"][fw_name]
        n_rxn = len(fw_info["reactions"])
        # Determine if this framework has any species with observation data
        no_obs_tag = ""
        if species_obs_availability is not None:
            fw_species = set()
            for rxn in fw_info["reactions"]:
                c = rxn.get("consumed", "NONE")
                p = rxn.get("produced", "NONE")
                if c and c != "NONE":
                    fw_species.add(c)
                if p and p != "NONE":
                    fw_species.add(p)
            has_any_obs = any(
                species_obs_availability.get(sp, {}).get("has_obs", False)
                for sp in fw_species
            )
            if not has_any_obs:
                no_obs_tag = (
                    ' <span style="font-size:0.75em;opacity:0.7;'
                    'font-style:italic;">(no obs data)</span>'
                )
        toggle_items.append(
            f'<button class="fw-toggle-chip active" '
            f'data-fw="{_svg_escape(fw_name)}" '
            f'style="--chip-color:{colour};">'
            f'<span class="chip-tick">\u2714</span>'
            f'<span class="chip-swatch" style="background:{colour};"></span>'
            f'{_svg_escape(fw_name)} ({n_rxn}){no_obs_tag}</button>')
    # "Only species with observation data" toggle
    obs_toggle_html = ""
    if species_obs_availability is not None:
        has_any_obs = any(v.get("has_obs") for v in species_obs_availability.values())
        has_any_no_obs = any(not v.get("has_obs") for v in species_obs_availability.values())
        if has_any_obs and has_any_no_obs:
            obs_toggle_html = (
                '<div style="margin-bottom:.5rem;">'
                '<label class="obs-only-toggle" style="display:inline-flex;'
                'align-items:center;gap:6px;cursor:pointer;user-select:none;'
                'font-size:.82rem;font-weight:600;color:var(--text,#222);'
                'padding:5px 12px;border:1.5px solid var(--primary,#4d9ee8);'
                'border-radius:8px;background:rgba(77,158,232,.06);'
                'transition:background .2s;">'
                '<input type="checkbox" class="obs-only-cb" '
                'style="width:16px;height:16px;accent-color:var(--primary,#4d9ee8);'
                'cursor:pointer;"/>'
                '<span>Only species with observation data</span></label>'
                '</div>'
            )

    toggle_bar = (
        '<div class="bgc-fw-toggle-bar">'
        + "".join(toggle_items)
        + "</div>"
    )

    # 4. Zoom controls
    zoom_controls = (
        '<div class="bgc-zoom-bar">'
        '<button class="zoom-btn zoom-lock locked" data-zoom="lock" '
        'title="Unlock +/\u2212 zoom">&#x1f512;</button>'
        '<button class="zoom-btn" data-zoom="in" title="Zoom in" disabled>+</button>'
        '<button class="zoom-btn" data-zoom="out" title="Zoom out" disabled>\u2212</button>'
        '<button class="zoom-btn" data-zoom="reset" title="Reset zoom">&#8634;</button>'
        '<span class="zoom-hint">Drag to select area</span>'
        '</div>'
    )

    # 5. Warning banner for sub-cycles without observation data
    obs_warning = ""
    if species_obs_availability is not None:
        no_obs_fws = []
        for fw_name, fw_info in network["frameworks"].items():
            fw_species = set()
            for rxn in fw_info["reactions"]:
                c = rxn.get("consumed", "NONE")
                p = rxn.get("produced", "NONE")
                if c and c != "NONE":
                    fw_species.add(c)
                if p and p != "NONE":
                    fw_species.add(p)
            has_any = any(
                species_obs_availability.get(sp, {}).get("has_obs", False)
                for sp in fw_species
            )
            if not has_any:
                no_obs_fws.append(fw_name)
        no_obs_species = [
            sp["name"] for sp in network["species"]
            if not species_obs_availability.get(
                sp["name"], {}
            ).get("has_obs", False)
        ]
        if no_obs_fws or no_obs_species:
            parts = []
            if no_obs_fws:
                fw_list = ", ".join(f"<strong>{f}</strong>" for f in no_obs_fws)
                parts.append(
                    f"Sub-cycles without observation data: {fw_list}"
                )
            if no_obs_species:
                sp_list = ", ".join(
                    f"<code>{s}</code>" for s in no_obs_species
                )
                parts.append(
                    f"Species without observations: {sp_list}"
                )
            obs_warning = (
                '<div style="margin-top:0.5rem;padding:0.5rem 0.8rem;'
                'border-left:4px solid #e0a040;background:rgba(224,160,64,.08);'
                'border-radius:4px;font-size:0.85em;">'
                '<span style="color:#c47800;font-weight:600;">'
                '&#x26a0; Calibration Note</span><br>'
                + "<br>".join(parts)
                + '<br><em style="opacity:0.7;">Species boxes and arrows with dashed '
                'lines lack observation data for calibration.</em>'
                '</div>'
            )

    # 6. Wrap in container: title → toggle bar → zoom → SVG → JS → warning
    n_sp = len(network["species"])
    title = f"Reaction Network \u2014 {n_sp} species"
    return (
        f'<div class="bgc-diagram-container">'
        f'<div class="bgc-diagram-header">'
        f'<div class="bgc-diagram-title">{title}</div>'
        f'{zoom_controls}'
        f'</div>'
        f'{obs_toggle_html}'
        f'{toggle_bar}'
        f'{svg_html}'
        f'{_bgc_diagram_js()}'
        f'{obs_warning}'
        f'</div>'
    )
