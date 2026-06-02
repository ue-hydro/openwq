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
Gen_BGC_Diagram.py — Standalone BGC reaction network SVG diagram generator.

Generates interactive SVG diagrams showing biogeochemical cycling frameworks,
species nodes, and reaction arrows. Supports both NATIVE_BGC_FLEX JSON
templates and PHREEQC .pqi files.

This module is self-contained — it does NOT depend on calibration_lib or
any external packages beyond the Python standard library.
"""

import json
import math
import os
import re
from typing import List, Dict, Any, Optional
import html as html_lib


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
