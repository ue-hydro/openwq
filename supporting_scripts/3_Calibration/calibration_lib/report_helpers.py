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
    .bgc-diagram-title {
        font-size: .9rem; font-weight: 600; margin-bottom: .6rem;
        color: var(--text);
    }
    .bgc-diagram-legend {
        display: flex; flex-wrap: wrap; gap: .8rem;
        margin-top: .6rem; font-size: .72rem; color: var(--text2);
    }
    .bgc-diagram-legend .legend-item {
        display: flex; align-items: center; gap: .3rem;
    }
    .bgc-diagram-legend .legend-swatch {
        width: 16px; height: 3px; border-radius: 2px;
    }
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

        rows.append(f"""
<tr class="param-row" data-param-idx="{gi}" data-param='{html_lib.escape(param_json)}'>
    <td style="text-align:center;">
        <input type="checkbox" class="param-cb" data-idx="{gi}" checked
               style="accent-color:var(--primary);"/>
    </td>
    <td><code style="font-size:.8rem;">{html_lib.escape(name)}</code></td>
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
                reactions.append({
                    "name": t.get("_NAME", lt[idx]),
                    "consumed": t.get("CONSUMED", "NONE"),
                    "produced": t.get("PRODUCED", "NONE"),
                })
        else:
            # Legacy dict format
            for key, val in fw.items():
                if key.startswith("_") or not isinstance(val, dict):
                    continue
                if "KINETICS" not in val and "FORMULA" not in val:
                    continue
                reactions.append({
                    "name": key,
                    "consumed": val.get("CONSUMED", "NONE"),
                    "produced": val.get("PRODUCED", "NONE"),
                })
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
        # -formula line
        if current_rxn and line.lower().startswith("-formula"):
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
                    reactions.append({
                        "name": current_rxn,
                        "consumed": c,
                        "produced": p,
                    })

    if not reactions:
        return {}

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
    # Determine which species are actually involved in reactions
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
    # Ellipse radii (wider than tall)
    rx = min(width * 0.38, 340)
    ry = min(height * 0.34, 200)

    for i, sp in enumerate(main_nodes):
        angle = 2 * math.pi * i / n - math.pi / 2  # start from top
        x = cx + rx * math.cos(angle)
        y = cy + ry * math.sin(angle)
        positions[sp] = (x, y)

    # Place extra nodes in a row at the bottom
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


def _render_bgc_svg(
    network: dict,
    max_width: int = 860,
) -> str:
    """Render a BGC reaction network as an inline SVG string."""

    all_species = [s["name"] for s in network["species"]]
    # Collect all reactions across frameworks
    all_reactions = []
    fw_names = list(network["frameworks"].keys())
    for fw_name in fw_names:
        for r in network["frameworks"][fw_name]["reactions"]:
            all_reactions.append({**r, "_fw": fw_name})

    if not all_reactions:
        return ""

    # Determine SVG dimensions
    n_species = len(all_species)
    svg_w = max_width
    svg_h = max(380, min(250 + n_species * 22, 650))

    positions = _compute_diagram_layout(all_species, all_reactions, svg_w, svg_h)

    # Node dimensions
    node_h = 28
    char_w = 8  # approximate width per character

    # Build SVG
    S = []
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
    # Dashed arrow marker for source/sink
    S.append(
        '<marker id="arrow-ss" markerWidth="7" markerHeight="5" '
        'refX="7" refY="2.5" orient="auto" markerUnits="strokeWidth">'
        '<path d="M0,0 L7,2.5 L0,5 Z" fill="currentColor" '
        'class="source-sink-arrow"/></marker>')
    S.append("</defs>")

    # --- Draw reaction arrows first (behind nodes) ---
    # Track how many arrows exist between each pair to offset them
    pair_counts: Dict[tuple, int] = {}
    for r in all_reactions:
        c, p = r["consumed"], r["produced"]
        if c == "NONE" or p == "NONE":
            continue
        key = (min(c, p), max(c, p))
        pair_counts[key] = pair_counts.get(key, 0) + 1

    pair_idx: Dict[tuple, int] = {}
    for r in all_reactions:
        fi = fw_names.index(r["_fw"])
        colour = _FW_COLORS[fi % len(_FW_COLORS)][0]
        c, p = r["consumed"], r["produced"]

        if c == "NONE" and p == "NONE":
            continue

        # Source/sink arrows (dashed)
        if c == "NONE" and p in positions:
            px, py = positions[p]
            sx = px - 60
            sy = py - 30
            S.append(
                f'<g class="source-sink-arrow">'
                f'<path d="M{sx:.0f},{sy:.0f} Q{(sx+px)/2:.0f},{sy:.0f} '
                f'{px:.0f},{py - node_h/2:.0f}" '
                f'fill="none" stroke="{colour}" stroke-width="1.5" '
                f'stroke-dasharray="5 3" marker-end="url(#arrow-{fi})"/>'
                f'<text x="{(sx+px)/2 - 10:.0f}" y="{sy - 4:.0f}" '
                f'font-size="8" fill="{colour}" '
                f'style="font-style:italic;">'
                f'{_svg_escape(r["name"])}</text></g>')
            continue

        if p == "NONE" and c in positions:
            cx_p, cy_p = positions[c]
            ex = cx_p + 60
            ey = cy_p + 30
            S.append(
                f'<g class="source-sink-arrow">'
                f'<path d="M{cx_p:.0f},{cy_p + node_h/2:.0f} '
                f'Q{(cx_p+ex)/2:.0f},{ey:.0f} {ex:.0f},{ey:.0f}" '
                f'fill="none" stroke="{colour}" stroke-width="1.5" '
                f'stroke-dasharray="5 3" marker-end="url(#arrow-{fi})"/>'
                f'<text x="{(cx_p+ex)/2 + 4:.0f}" y="{ey + 12:.0f}" '
                f'font-size="8" fill="{colour}" '
                f'style="font-style:italic;">'
                f'{_svg_escape(r["name"])}</text></g>')
            continue

        if c not in positions or p not in positions:
            continue

        x1, y1 = positions[c]
        x2, y2 = positions[p]

        # Offset for multiple arrows between same pair
        key = (min(c, p), max(c, p))
        idx = pair_idx.get(key, 0)
        pair_idx[key] = idx + 1
        total = pair_counts.get(key, 1)
        offset = (idx - (total - 1) / 2) * 18

        # Compute bezier control point perpendicular to midpoint
        mx, my = (x1 + x2) / 2, (y1 + y2) / 2
        dx, dy = x2 - x1, y2 - y1
        length = math.sqrt(dx * dx + dy * dy) or 1
        # Perpendicular unit vector
        nx, ny = -dy / length, dx / length
        # Curve amount proportional to distance + offset
        curve = length * 0.18 + offset
        qx = mx + nx * curve
        qy = my + ny * curve

        # Shorten start/end to not overlap nodes
        node_w_c = max(len(c) * char_w + 20, 60)
        node_w_p = max(len(p) * char_w + 20, 60)
        # Move start point away from centre of node c
        t_start = max((node_w_c / 2 + 4) / length, 0.05)
        t_end = max((node_w_p / 2 + 4) / length, 0.05)
        ax = x1 + dx * t_start
        ay = y1 + dy * t_start
        bx = x2 - dx * t_end
        by = y2 - dy * t_end

        # Recalculate control relative to shortened segment
        smx, smy = (ax + bx) / 2, (ay + by) / 2
        qx2 = smx + nx * curve
        qy2 = smy + ny * curve

        S.append(
            f'<g class="reaction-arrow">'
            f'<path d="M{ax:.1f},{ay:.1f} Q{qx2:.1f},{qy2:.1f} '
            f'{bx:.1f},{by:.1f}" '
            f'fill="none" stroke="{colour}" stroke-width="1.8" '
            f'marker-end="url(#arrow-{fi})"/>')

        # Label at curve midpoint
        lx = (ax + 2 * qx2 + bx) / 4
        ly = (ay + 2 * qy2 + by) / 4
        label = r["name"].replace("_", " ")
        if len(label) > 22:
            label = label[:20] + "\u2026"
        S.append(
            f'<text x="{lx:.1f}" y="{ly - 4:.1f}" '
            f'text-anchor="middle" font-size="8.5" fill="{colour}" '
            f'style="font-weight:500;pointer-events:none;">'
            f'{_svg_escape(label)}</text>')
        S.append("</g>")

    # --- Draw species nodes (on top of arrows) ---
    for sp in all_species:
        if sp not in positions:
            continue
        x, y = positions[sp]
        label = sp
        node_w = max(len(label) * char_w + 24, 64)
        is_loss = "loss" in sp.lower() or "atm" in sp.lower()
        cls = "species-node sink" if is_loss else "species-node"
        S.append(
            f'<g class="{cls}">'
            f'<rect x="{x - node_w/2:.1f}" y="{y - node_h/2:.1f}" '
            f'width="{node_w:.0f}" height="{node_h}" rx="10" ry="10" '
            f'style="fill:var(--surface,#fff);stroke:var(--primary,#4d9ee8);'
            f'stroke-width:2;{" stroke-dasharray:4 2;" if is_loss else ""}"/>'
            f'<text x="{x:.1f}" y="{y:.1f}" text-anchor="middle" '
            f'dominant-baseline="central" font-size="11.5" font-weight="600" '
            f'style="fill:var(--text,#222);">'
            f'{_svg_escape(label)}</text></g>')

    S.append("</svg>")
    return "\n".join(S)


def generate_bgc_reaction_diagram(
    bgc_data: Optional[dict] = None,
    bgc_template_path: Optional[str] = None,
    pqi_filepath: Optional[str] = None,
    module_name: str = "NATIVE_BGC_FLEX",
    max_width: int = 860,
) -> str:
    """Generate an HTML/SVG reaction network diagram for BGC modules.

    Accepts either pre-parsed *bgc_data*, a path to a JSON template, or
    a path to a PHREEQC .pqi file.  Returns an HTML ``<div>`` containing
    the SVG diagram, or ``""`` if no data is available.
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
    svg_html = _render_bgc_svg(network, max_width=max_width)
    if not svg_html:
        return ""

    # 3. Framework legend
    fw_names = list(network["frameworks"].keys())
    legend_items = []
    for fi, fw_name in enumerate(fw_names):
        colour = _FW_COLORS[fi % len(_FW_COLORS)][0]
        n_rxn = len(network["frameworks"][fw_name]["reactions"])
        legend_items.append(
            f'<span class="legend-item">'
            f'<span class="legend-swatch" '
            f'style="background:{colour};"></span>'
            f'{_svg_escape(fw_name)} ({n_rxn} reactions)</span>')

    legend_html = ""
    if len(fw_names) > 1 or True:  # always show legend
        legend_html = (
            '<div class="bgc-diagram-legend">'
            + "".join(legend_items)
            + "</div>"
        )

    # 4. Wrap in container
    n_sp = len(network["species"])
    title = f"Reaction Network \u2014 {n_sp} species"
    return (
        f'<div class="bgc-diagram-container">'
        f'<div class="bgc-diagram-title">{title}</div>'
        f'{svg_html}'
        f'{legend_html}'
        f'</div>'
    )
