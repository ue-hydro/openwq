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
"""
The CALIBRATED MODEL CONFIG (``<model_template>_config_run.py``).

A calibration leaves its best parameter values inside the best evaluation's
input files, while the model-config template keeps generating the DEFAULT
values. The calibration results report therefore offers to save a
*calibrated model config*: the user's own model-config template with a short
block appended that (1) writes the inputs to a new run folder, (2) bakes the
generation-time calibrated values (climate-response parameters) as template
variables, and (3) tells the input generator to apply the calibrated best
parameters and the ML layers to the freshly generated, FULL-PERIOD inputs
exactly as the best evaluation was configured.

That file is a normal model-config template: running it produces the
calibrated model's inputs + config report (with the usual run / output-report
snippets), and the scenario template uses it as the base of what-if runs.
"""
from __future__ import annotations

import os
import re
import json
import logging
import importlib.util
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)

CALIBRATED_SETUP_FILE = "calibrated_setup.json"
_MARKER = "END OF USER CONFIGURATION"

# generation-time calibration parameters -> model-config template variables
_GEN_TIME_PARAMS = {
    "precip_scaling_power": "ss_climate_precip_scaling_power",
    "Q10_biological": "ss_climate_temp_q10",
    "T_reference": "ss_climate_temp_reference_c",
}


def config_run_path_for(template_path: str) -> str:
    """``<dir>/<template stem>_config_run.py`` next to the model template."""
    p = Path(template_path)
    return str(p.with_name(p.stem + "_config_run.py"))


def calibrated_dir_for(model_config: Dict[str, Any]) -> str:
    """The calibrated run folder: ``<dir2save_input_files>_calibrated``."""
    base = str(model_config.get("dir2save_input_files") or "").rstrip("/\\")
    if not base:
        exe = model_config.get("executable_path", "")
        base = os.path.join(os.path.dirname(os.path.abspath(exe)) if exe else os.getcwd(), "openwq_run")
    return base + "_calibrated"


def build_config_run_text(template_path: str, calibration_dir: str, run_script_path: str,
                          model_config: Dict[str, Any], best_params: Optional[Dict[str, float]] = None,
                          parameter_definitions: Optional[List[Dict[str, Any]]] = None,
                          best_eval: Optional[str] = None, best_objective=None,
                          objective_name: str = "") -> Tuple[str, str]:
    """The calibrated model config as text (+ its recommended path). The
    user's template is copied verbatim; the calibrated block goes right
    before the END-OF-USER-CONFIGURATION marker (later assignments override
    the template's own values, so the original lines stay untouched)."""
    text = Path(template_path).read_text(encoding="utf-8", errors="replace")
    out_dir = calibrated_dir_for(model_config)
    best_params = best_params or {}
    gen_lines = []
    for p in (parameter_definitions or []):
        if str(p.get("file_type", "")).startswith("ss_climate"):
            key = ((p.get("path") or {}).get("param") if isinstance(p.get("path"), dict) else None)
            var = _GEN_TIME_PARAMS.get(str(key))
            if var and p.get("name") in best_params:
                gen_lines.append(f"{var} = {float(best_params[p['name']])!r}")
    stamp = datetime.now().strftime("%Y-%m-%d %H:%M")
    block = [
        "",
        "# ╔" + "═" * 72 + "╗",
        "# ║  CALIBRATED MODEL CONFIG — written by the calibration results report" + " " * 2 + "║",
        "# ╚" + "═" * 72 + "╝",
        f"#  Saved {stamp}. Base template: {os.path.basename(template_path)}",
        f"#  Calibration: {calibration_dir}",
        (f"#  Best evaluation: {best_eval}" + (f"  ({objective_name} = {best_objective})" if best_objective is not None else ""))
        if best_eval else "#  Best evaluation: (see results/best_parameters.json)",
        "#",
        "#  Running this file generates the inputs for the FULL model period in the",
        "#  folder below, then applies the calibrated best parameter values and the",
        "#  ML layers to them (exactly as the best evaluation was configured), and",
        "#  writes the usual config report with the run / output-report snippets.",
        "#  It is the base of scenario analyses (supporting_scripts/4_Scenarios).",
        f'dir2save_input_files = {out_dir!r}',
        f'calibrated_from_calibration_dir = {os.path.abspath(calibration_dir)!r}',
        f'calibrated_from_run_script = {os.path.abspath(run_script_path)!r}',
        "force_regenerate = True   # a calibrated run folder is always rebuilt from scratch",
    ]
    if gen_lines:
        block.append("# generation-time calibrated values (climate response of the loads)")
        block.extend(gen_lines)
    block.append("")
    block_txt = "\n".join(block) + "\n"
    lines = text.split("\n")
    idx = next((i for i, l in enumerate(lines) if _MARKER in l), None)
    if idx is None:
        # no marker: insert before the first import of the generator call block
        idx = next((i for i, l in enumerate(lines) if l.startswith("from Gen_Report import") or
                    l.startswith("import Gen_Report")), len(lines))
    # step back over the box-drawing lines that frame the marker
    while idx > 0 and lines[idx - 1].startswith("# ╔"):
        idx -= 1
    new = "\n".join(lines[:idx]) + "\n" + block_txt + "\n".join(lines[idx:])
    return new, config_run_path_for(template_path)


def _load_run_script(run_script: str):
    spec = importlib.util.spec_from_file_location("owq_calibration_run_for_apply", run_script)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def apply_calibrated_best(dir2save: str, calibration_dir: str, run_script: str,
                          model_config: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    """Apply the calibration's best parameter values (+ ML layers) to the
    inputs in ``dir2save`` (freshly generated, full period). Returns the
    provenance dict also written to ``<dir2save>/calibrated_setup.json``."""
    import numpy as np
    from .parameter_handler import ParameterHandler
    calibration_dir = os.path.abspath(calibration_dir)
    dir2save = os.path.abspath(dir2save)
    bp = os.path.join(calibration_dir, "results", "best_parameters.json")
    if not os.path.isfile(bp):
        raise FileNotFoundError(f"no calibrated values found: {bp} (has the calibration finished?)")
    best = json.load(open(bp))
    if not os.path.isfile(run_script):
        raise FileNotFoundError(f"calibration run script not found: {run_script}")
    cal = _load_run_script(run_script)
    params = list(getattr(cal, "calibration_parameters", []) or [])
    values = np.array([float(best.get(p["name"], p.get("initial", 0.0))) for p in params], dtype=float)
    n_hit = sum(1 for p in params if p["name"] in best)
    cfg = dict(model_config or {})
    cfg.setdefault("dir2save_input_files", dir2save)
    handler = ParameterHandler(calibration_work_dir=calibration_dir, model_config=cfg,
                               running_on_docker=True, calibration_period=None,
                               ml_closures=getattr(cal, "ml_closures", None),
                               ml_runtime=getattr(cal, "ml_runtime", None))
    if params or handler.ml_closures or handler.ml_runtime:
        handler.apply_parameters(Path(dir2save), params, values)
    # provenance
    best_eval, best_obj = None, None
    try:
        ores = json.load(open(os.path.join(calibration_dir, "results", "optimization_results.json")))
        best_obj = ores.get("best_objective")
    except Exception:
        pass
    try:
        from .calibration_driver import _best_eval_dir
        bd = _best_eval_dir(Path(calibration_dir))
        best_eval = os.path.basename(bd) if bd else None
    except Exception:
        pass
    prov = {"schema": "openwq_calibrated_setup/1",
            "applied": datetime.now().isoformat(timespec="seconds"),
            "calibration_dir": calibration_dir, "run_script": os.path.abspath(run_script),
            "best_eval": best_eval, "best_objective": best_obj,
            "objective_function": getattr(cal, "objective_function", None),
            "n_parameters": len(params), "n_from_best": n_hit,
            "parameters": {p["name"]: float(v) for p, v in zip(params, values)},
            "ml_closures": bool(handler.ml_closures), "ml_runtime": bool(handler.ml_runtime),
            "mapping_json": (os.path.join(calibration_dir, "ml_attributes", "mapping.json")
                             if os.path.isfile(os.path.join(calibration_dir, "ml_attributes", "mapping.json")) else None)}
    json.dump(prov, open(os.path.join(dir2save, CALIBRATED_SETUP_FILE), "w"), indent=2)
    logger.info(f"Calibrated best applied to {dir2save}: {n_hit}/{len(params)} parameter(s) "
                f"from {bp}" + (f", best eval {best_eval}" if best_eval else ""))
    return prov


def read_calibrated_setup(dir2save: str) -> Optional[Dict[str, Any]]:
    p = os.path.join(str(dir2save or ""), CALIBRATED_SETUP_FILE)
    if os.path.isfile(p):
        try:
            return json.load(open(p))
        except Exception:
            return None
    return None
