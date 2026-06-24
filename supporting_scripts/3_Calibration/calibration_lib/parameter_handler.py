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
Parameter Handler Module
========================

Maps calibration parameters to their target configuration files and applies values.
Supports all OpenWQ module types: BGC, PHREEQC, Sorption, Sediment, Transport,
Lateral Exchange, and Source/Sink methods.

Two modes of operation:
  1. **Legacy mode** (test_case_dir): copies config files from a test case directory.
  2. **Config-integrated mode** (model_config): generates configs via Gen_Input_Driver.
"""

import json
import os
import re
import copy
import shutil
from pathlib import Path
from typing import Dict, List, Any, Tuple, Optional, Union
import numpy as np
import pandas as pd
import logging

logger = logging.getLogger(__name__)


class ParameterHandler:
    """
    Handles reading, modifying, and writing parameter values to OpenWQ config files.
    """

    def __init__(self,
                 calibration_work_dir: str,
                 model_config: Optional[Dict[str, Any]] = None,
                 base_model_config_dir: Optional[str] = None,
                 test_case_dir: Optional[str] = None,
                 running_on_docker: bool = True,
                 calibration_period: Optional[Tuple[str, str]] = None):
        """
        Initialize the parameter handler.

        Use **model_config** for the integrated workflow (recommended).
        Use **base_model_config_dir** + **test_case_dir** for legacy mode.

        Parameters
        ----------
        calibration_work_dir : str
            Path to working directory for calibration evaluations.
        model_config : Dict[str, Any], optional
            Loaded model configuration (from config_integration.load_model_config).
            When provided, eval directories are populated via Gen_Input_Driver.
        base_model_config_dir : str, optional
            [LEGACY] Path to the Model_Config directory.
        test_case_dir : str, optional
            [LEGACY] Path to the test case directory (openwq_in, mizuroute_in, …).
        running_on_docker : bool
            Whether the model runs in Docker (affects path handling).
        """
        self.calibration_work_dir = Path(calibration_work_dir)
        self.running_on_docker = running_on_docker
        # Calibration window — forwarded to the per-eval SS generation so the
        # size guard is evaluated over the window (the SS is trimmed to it).
        self.calibration_period = calibration_period

        # Integrated-mode config
        self.model_config = model_config

        # Legacy-mode paths
        self.base_model_config_dir = Path(base_model_config_dir) if base_model_config_dir else None
        self.test_case_dir = Path(test_case_dir) if test_case_dir else None

        # Ensure directories exist
        self.calibration_work_dir.mkdir(parents=True, exist_ok=True)
        (self.calibration_work_dir / "evaluations").mkdir(exist_ok=True)

    # =====================================================================
    # Evaluation directory setup
    # =====================================================================

    # Generation-time SS parameters: these reshape/retrain the *generated* SS
    # JSON rather than editing it post-hoc, so they must be injected into the
    # model config BEFORE Gen_Input_Driver runs.  Keyed by file_type; each
    # entry maps the parameter's path["param"] (set in extract_parameters) to
    # a (model_config key, value-type) pair.  The value-type coerces the
    # optimiser's real-scale float to what the config key expects — integer
    # hyperparameters (e.g. ML tree count/depth) are rounded.
    _GENERATION_TIME_SS_PARAMS = {
        "ss_climate": {
            "precip_scaling_power": ("ss_climate_precip_scaling_power", "float"),
            "Q10_biological":       ("ss_climate_temp_q10", "float"),
            "T_reference":          ("ss_climate_temp_reference_c", "float"),
        },
        "ss_ml": {
            "n_estimators": ("ss_ml_n_estimators", "int"),
            "max_depth":    ("ss_ml_max_depth", "int"),
        },
    }

    def setup_working_directory(self,
                                eval_id: int,
                                parameters: Optional[List[Dict]] = None,
                                values: Optional[np.ndarray] = None) -> Path:
        """
        Create isolated working directory for a single model evaluation.

        If *model_config* was provided at init, the directory is populated
        by calling Gen_Input_Driver.  Otherwise it falls back to copying
        files from test_case_dir (legacy behaviour).

        Parameters
        ----------
        eval_id : int
            Unique evaluation identifier.
        parameters, values : optional
            The calibration parameter definitions and their real-scale
            values for this evaluation.  When provided, *generation-time*
            parameters (currently the dynamic SS climate-response params,
            ``file_type == "ss_climate"``) are baked into the per-eval model
            config before the SS JSON is generated, because they change how
            the annual load is distributed across months rather than scaling
            an already-written value.  All other parameters are applied
            afterwards via :meth:`apply_parameters`.

        Returns
        -------
        Path
            Path to the evaluation working directory.
        """
        if self.model_config is not None:
            eval_config = self._build_eval_config(parameters, values)
            # Fast path — reuse the config the model-config script already
            # produced.  When there are NO generation-time parameter overrides
            # (climate / ML), the per-eval config is byte-for-byte what
            # Gen_Input_Driver would regenerate, so we copy the baseline
            # (including the already-processed Copernicus SS JSON, BGC, etc.)
            # instead of re-running the whole generator — much faster and it
            # skips re-processing GRQA/Copernicus on every evaluation.
            # apply_parameters() then patches the calibrated values on the
            # copy.  (`_build_eval_config` returns the SAME object when no
            # generation-time overrides were injected.)
            if eval_config is self.model_config:
                baseline_dir = self._setup_from_baseline_config(eval_id)
                if baseline_dir is not None:
                    return baseline_dir
            return self._setup_from_config(eval_id, eval_config)
        elif self.test_case_dir is not None:
            return self._setup_from_testcase(eval_id)
        else:
            raise RuntimeError(
                "ParameterHandler: either model_config or test_case_dir "
                "must be provided."
            )

    def _build_eval_config(self,
                           parameters: Optional[List[Dict]],
                           values: Optional[np.ndarray]) -> Dict[str, Any]:
        """
        Return the model config to use for generating this evaluation's
        input files, with any generation-time parameters overridden by their
        calibrated values.

        Returns ``self.model_config`` unchanged when there are no
        generation-time parameters to inject (the common case), so the hot
        path allocates nothing extra.
        """
        if not parameters or values is None:
            return self.model_config

        overrides: Dict[str, Any] = {}
        for i, param in enumerate(parameters):
            fmap = self._GENERATION_TIME_SS_PARAMS.get(param.get("file_type"))
            if not fmap:
                continue
            path = param.get("path", {})
            pkey = path.get("param") if isinstance(path, dict) else None
            spec = fmap.get(pkey)
            if spec is None:
                continue
            cfg_key, vtype = spec
            try:
                raw = float(values[i])
            except (TypeError, ValueError, IndexError):
                continue
            overrides[cfg_key] = int(round(raw)) if vtype == "int" else raw

        if not overrides:
            return self.model_config

        eval_config = copy.deepcopy(self.model_config)
        eval_config.update(overrides)
        logger.debug(
            f"Baked generation-time SS params into eval config: {overrides}")
        return eval_config

    def _setup_from_config(self,
                           eval_id: int,
                           model_config: Optional[Dict[str, Any]] = None) -> Path:
        """
        Generate config files via Gen_Input_Driver for an evaluation.

        This is the recommended workflow when the calibration template
        imports model_config_template.py.

        ``model_config`` may be an override (e.g. with generation-time
        climate params baked in); it defaults to ``self.model_config``.
        """
        from . import config_integration

        cfg = model_config if model_config is not None else self.model_config

        eval_dir = self.calibration_work_dir / "evaluations" / f"eval_{eval_id:04d}"

        # Clean up if exists
        if eval_dir.exists():
            shutil.rmtree(eval_dir)

        eval_dir.mkdir(parents=True)

        # Create output directory (model will write here)
        (eval_dir / "openwq_out" / "HDF5").mkdir(parents=True, exist_ok=True)

        # Generate config files using Gen_Input_Driver
        try:
            config_integration.generate_config_for_eval(
                cfg, str(eval_dir), suppress_report=True,
                calibration_period=self.calibration_period
            )
        except Exception as e:
            logger.error(
                f"Failed to generate config for eval_{eval_id:04d}: {e}"
            )
            raise

        logger.debug(
            f"Created evaluation directory (Gen_Input_Driver): {eval_dir}"
        )
        return eval_dir

    def _setup_from_baseline_config(self, eval_id: int) -> Optional[Path]:
        """Fast eval setup: copy the config the model-config script already
        generated (under ``dir2save_input_files``) into the eval directory,
        reusing the processed Copernicus SS JSON, BGC, transport, etc.
        verbatim — instead of re-running Gen_Input_Driver.

        Safe because, with no generation-time parameter overrides, the
        regenerated config would be identical to this baseline; the per-eval
        differences are applied afterwards by :meth:`apply_parameters`, which
        edits the copied JSONs in place.

        The OpenWQ master file uses paths relative to the run directory
        (``openwq_in/…``, ``openwq_out/``) and the model is launched with its
        CWD set to the eval dir, so the copied config is location-independent
        — no path rewriting needed.

        Returns the eval dir, or ``None`` when no baseline config is available
        (caller then falls back to full regeneration).
        """
        mc = self.model_config or {}
        base = mc.get('dir2save_input_files')
        if not base:
            exe = mc.get('executable_path', '')
            base = os.path.dirname(os.path.abspath(exe)) if exe else None
        if not base:
            return None

        base = Path(base)
        src_in = base / "openwq_in"
        src_master = base / "openWQ_master.json"
        if not src_in.is_dir() or not src_master.is_file():
            logger.debug(
                f"No baseline config at {base}; will regenerate per eval.")
            return None

        eval_dir = self.calibration_work_dir / "evaluations" / f"eval_{eval_id:04d}"
        if eval_dir.exists():
            shutil.rmtree(eval_dir)
        eval_dir.mkdir(parents=True)

        # Copy the processed input config, skipping report-only / bulky
        # artefacts (clipped GRQA data, any previously-written HDF5 output).
        shutil.copytree(
            src_in, eval_dir / "openwq_in",
            ignore=shutil.ignore_patterns(
                'grqa_clipped_data', 'HDF5', '*.h5', '*.nc'))
        shutil.copy2(src_master, eval_dir / "openWQ_master.json")
        (eval_dir / "openwq_out" / "HDF5").mkdir(parents=True, exist_ok=True)

        # Apply the calibration's run_mode_debug override to the copied master.
        # The baseline master keeps whatever RUN_MODE_DEBUG the model-config
        # script wrote, and this fast path copies it VERBATIM — so without this
        # patch the setup-report "debug OFF" choice is silently ignored and
        # openWQ keeps writing the 5 extra diagnostic HDF5 files per
        # species/compartment (chemistry/transport/ewf/ic/ss).  The override is
        # set on self.model_config by run_calibration; mirror it into the
        # COMPUTATIONAL_SETTINGS.RUN_MODE_DEBUG flag of the per-eval master.
        rmd = mc.get("run_mode_debug")
        if rmd is not None:
            try:
                dst_master = eval_dir / "openWQ_master.json"
                data, header = self._read_json_with_header(dst_master)
                cs = data.get("COMPUTATIONAL_SETTINGS")
                if isinstance(cs, dict):
                    if bool(cs.get("RUN_MODE_DEBUG")) != bool(rmd):
                        cs["RUN_MODE_DEBUG"] = bool(rmd)
                        self._write_json_with_header(dst_master, data, header)
                        logger.debug(
                            f"eval_{eval_id:04d}: set master RUN_MODE_DEBUG="
                            f"{bool(rmd)} (calibration override)")
            except Exception as e:
                logger.warning(
                    f"Could not apply run_mode_debug override to baseline "
                    f"master (eval_{eval_id:04d}): {e}")

        logger.debug(
            f"Created evaluation directory (reused baseline config from "
            f"{base}): {eval_dir}")
        return eval_dir

    def _setup_from_testcase(self, eval_id: int) -> Path:
        """
        Legacy: copy files from test_case_dir to create eval directory.
        """
        eval_dir = self.calibration_work_dir / "evaluations" / f"eval_{eval_id:04d}"

        # Clean up if exists
        if eval_dir.exists():
            shutil.rmtree(eval_dir)

        eval_dir.mkdir(parents=True)

        # Copy test case structure
        src_openwq_in = self.test_case_dir / "openwq_in"
        dst_openwq_in = eval_dir / "openwq_in"
        if src_openwq_in.exists():
            shutil.copytree(src_openwq_in, dst_openwq_in)

        # Copy openWQ_master.json
        src_master = self.test_case_dir / "openWQ_master.json"
        if src_master.exists():
            shutil.copy2(src_master, eval_dir / "openWQ_master.json")

        # Create output directory
        (eval_dir / "openwq_out" / "HDF5").mkdir(parents=True, exist_ok=True)

        logger.debug(f"Created evaluation directory (copy): {eval_dir}")
        return eval_dir

    def apply_parameters(self,
                        eval_dir: Path,
                        parameters: List[Dict],
                        values: np.ndarray) -> None:
        """
        Apply parameter values to configuration files in evaluation directory.

        Parameters
        ----------
        eval_dir : Path
            Path to evaluation working directory
        parameters : List[Dict]
            Parameter definitions from calibration config
        values : np.ndarray
            Parameter values to apply (in real scale, after transform)
        """
        ss_items = []   # source/sink load params, applied together after the
                        # loop (seasonal modes group sibling params by species)
        _SS_LOAD_TYPES = ("ss_csv_scale", "ss_seasonal_amp",
                          "ss_seasonal_phase", "ss_seasonal_month")
        for i, param in enumerate(parameters):
            value = values[i]
            file_type = param["file_type"]
            path = param["path"]
            dtype = param.get("dtype", "float")

            # Convert to appropriate type
            if dtype == "int":
                value = int(round(value))

            logger.debug(f"Applying {param['name']} = {value} (type={file_type})")

            if file_type in _SS_LOAD_TYPES:
                ss_items.append((param, value))
                continue

            # Route to appropriate handler
            if file_type == "bgc_json":
                self._apply_bgc_json_param(eval_dir, path, value)
            elif file_type == "phreeqc_pqi":
                self._apply_phreeqc_pqi_param(eval_dir, path, value)
            elif file_type == "sorption_json":
                self._apply_sorption_param(eval_dir, path, value)
            elif file_type == "sediment_json":
                self._apply_sediment_param(eval_dir, path, value)
            elif file_type == "transport_json":
                self._apply_transport_param(eval_dir, path, value)
            elif file_type == "lateral_exchange_json":
                self._apply_lateral_exchange_param(eval_dir, path, value)
            elif file_type in ("ss_copernicus_static", "ss_copernicus_dynamic"):
                initial = param.get("initial", None)
                self._apply_ss_copernicus_coeff(eval_dir, path, value, initial_value=initial)
            elif file_type in ("ss_climate", "ss_climate_param",
                                "ss_ml", "ss_ml_param"):
                # Generation-time parameter: the dynamic SS climate-response
                # params reshape the monthly load distribution, and the ML SS
                # hyperparameters (tree count/depth) retrain the model — both
                # change how the SS JSON is *generated*, so they are baked
                # into the eval config at setup_working_directory() time (see
                # _build_eval_config).  Nothing to edit post-hoc here.
                logger.debug(
                    f"{file_type} param '{param['name']}' is applied at "
                    f"generation time; skipping post-hoc edit.")
            else:
                logger.warning(f"Unknown file_type: {file_type} for parameter {param['name']}")

        # Apply all source/sink load params together (grouped by species) so
        # seasonal modes build a per-month factor from sibling params.
        if ss_items:
            self._apply_ss_loads(eval_dir, ss_items)

    # =========================================================================
    # BGC JSON (NATIVE_BGC_FLEX)
    # =========================================================================

    def _apply_bgc_json_param(self,
                              eval_dir: Path,
                              json_path: List[str],
                              value: float) -> None:
        """
        Modify a parameter in the BGC JSON file.

        Parameters
        ----------
        eval_dir : Path
            Evaluation directory
        json_path : List[str]
            Path within JSON: e.g., ["CYCLING_FRAMEWORKS", "N_cycle", "3", "parameter_values", "k"]
        value : float
            New parameter value
        """
        bgc_file = eval_dir / "openwq_in" / "openWQ_MODULE_NATIVE_BGC_FLEX.json"

        if not bgc_file.exists():
            logger.warning(f"BGC file not found: {bgc_file}")
            return

        data, header = self._read_json_with_header(bgc_file)

        # Navigate to the parameter
        obj = data
        for key in json_path[:-1]:
            obj = obj[key]
        obj[json_path[-1]] = value

        self._write_json_with_header(bgc_file, data, header)

    # =========================================================================
    # PHREEQC PQI File
    # =========================================================================

    def _apply_phreeqc_pqi_param(self,
                                  eval_dir: Path,
                                  path: Dict,
                                  value: float) -> None:
        """
        Modify a parameter in the PHREEQC .pqi input file.

        Parameters
        ----------
        eval_dir : Path
            Evaluation directory
        path : Dict
            Parameter location: {"block": "SOLUTION", "species": "N(5)"} or
            {"block": "EQUILIBRIUM_PHASES", "phase": "CO2(g)", "field": "si"}
        value : float
            New parameter value
        """
        # Find .pqi file in openwq_in
        pqi_files = list((eval_dir / "openwq_in").glob("*.pqi"))
        if not pqi_files:
            logger.warning(f"No .pqi file found in {eval_dir / 'openwq_in'}")
            return

        pqi_file = pqi_files[0]

        with open(pqi_file, 'r') as f:
            lines = f.readlines()

        block = path.get("block", "")
        modified = False

        if block == "SOLUTION":
            species = path.get("species", "")
            lines, modified = self._modify_pqi_solution(lines, species, value)
        elif block == "EQUILIBRIUM_PHASES":
            phase = path.get("phase", "")
            field = path.get("field", "si")
            lines, modified = self._modify_pqi_equilibrium(lines, phase, field, value)
        elif block == "KINETICS":
            reaction = path.get("reaction", "")
            param = path.get("param", "-parms")
            index = path.get("index", 0)
            lines, modified = self._modify_pqi_kinetics(lines, reaction, param, index, value)
        elif block == "SURFACE":
            surface = path.get("surface", "")
            field = path.get("field", "sites")
            lines, modified = self._modify_pqi_surface(lines, surface, field, value)

        if modified:
            with open(pqi_file, 'w') as f:
                f.writelines(lines)
            logger.debug(f"Modified PHREEQC parameter in {block}")
        else:
            logger.warning(f"Could not find PHREEQC parameter: {path}")

    def _is_phreeqc_block_start(self, line: str) -> bool:
        """Check if line starts a new PHREEQC block."""
        block_keywords = [
            "SOLUTION", "EQUILIBRIUM_PHASES", "KINETICS", "SURFACE",
            "EXCHANGE", "GAS_PHASE", "SOLID_SOLUTIONS", "REACTION",
            "MIX", "USE", "SAVE", "END", "TITLE", "PRINT", "SELECTED_OUTPUT"
        ]
        stripped = line.strip().upper()
        for keyword in block_keywords:
            if stripped.startswith(keyword):
                return True
        return False

    def _modify_pqi_solution(self, lines: List[str], species: str, value: float) -> Tuple[List[str], bool]:
        """Modify a species concentration in SOLUTION block."""
        in_solution = False
        for i, line in enumerate(lines):
            stripped = line.strip().upper()
            if stripped.startswith("SOLUTION"):
                in_solution = True
            elif in_solution and stripped and not stripped.startswith("#"):
                # Check if we've left the SOLUTION block
                if self._is_phreeqc_block_start(line) and not stripped.startswith("SOLUTION"):
                    in_solution = False
                # Match species name (handle cases like "N(5)" or "Ca")
                elif stripped:
                    first_word = stripped.split()[0]
                    if species.upper() == first_word or species.upper() in first_word:
                        # Found the species line - replace value
                        parts = line.split()
                        if len(parts) >= 2:
                            # Handle "N(5)  2.0 as N" format
                            parts[1] = str(value)
                            lines[i] = "    " + "  ".join(parts) + "\n"
                            return lines, True
        return lines, False

    def _modify_pqi_equilibrium(self, lines: List[str], phase: str, field: str, value: float) -> Tuple[List[str], bool]:
        """
        Modify equilibrium phase parameters.

        PHREEQC format: "Phase_name  SI  amount"
        - field="si": modify saturation index (2nd value)
        - field="amount": modify initial moles (3rd value)
        """
        in_equil = False
        for i, line in enumerate(lines):
            stripped = line.strip().upper()
            if stripped.startswith("EQUILIBRIUM_PHASES"):
                in_equil = True
            elif in_equil and stripped and not stripped.startswith("#"):
                # Check if we've left the block
                if self._is_phreeqc_block_start(line) and not stripped.startswith("EQUILIBRIUM_PHASES"):
                    in_equil = False
                elif phase.upper() in stripped:
                    # Found the phase line
                    parts = line.split()
                    if len(parts) >= 2:
                        if field.lower() == "si":
                            parts[1] = str(value)  # SI is second value
                        elif field.lower() == "amount" and len(parts) >= 3:
                            parts[2] = str(value)  # Amount is third value
                        lines[i] = "    " + "    ".join(parts) + "\n"
                        return lines, True
        return lines, False

    def _modify_pqi_kinetics(self, lines: List[str], reaction: str, param: str,
                             index: int, value: float) -> Tuple[List[str], bool]:
        """Modify kinetics parameters."""
        # Simplified implementation - would need expansion for complex KINETICS blocks
        in_kinetics = False
        in_reaction = False
        for i, line in enumerate(lines):
            stripped = line.strip().upper()
            if stripped.startswith("KINETICS"):
                in_kinetics = True
            elif in_kinetics:
                if stripped.startswith("END"):
                    in_kinetics = False
                elif reaction.upper() in stripped:
                    in_reaction = True
                elif in_reaction and param.upper() in stripped:
                    parts = line.split()
                    if len(parts) > index + 1:
                        parts[index + 1] = str(value)
                        lines[i] = "    " + " ".join(parts) + "\n"
                        return lines, True
        return lines, False

    def _modify_pqi_surface(self, lines: List[str], surface: str, field: str,
                            value: float) -> Tuple[List[str], bool]:
        """Modify surface parameters."""
        # Simplified implementation
        in_surface = False
        for i, line in enumerate(lines):
            stripped = line.strip().upper()
            if stripped.startswith("SURFACE"):
                in_surface = True
            elif in_surface:
                if stripped.startswith("END"):
                    in_surface = False
                elif surface.upper() in stripped:
                    parts = line.split()
                    if len(parts) >= 2:
                        parts[1] = str(value)
                        lines[i] = "    " + " ".join(parts) + "\n"
                        return lines, True
        return lines, False

    # =========================================================================
    # Sorption Isotherm JSON
    # =========================================================================

    def _apply_sorption_param(self,
                              eval_dir: Path,
                              path: Dict,
                              value: float) -> None:
        """
        Modify a parameter in the sorption isotherm JSON file.

        Parameters
        ----------
        eval_dir : Path
            Evaluation directory
        path : Dict
            Parameter location: {"module": "FREUNDLICH", "species": "NH4-N", "param": "Kfr"}
            or {"module": "FREUNDLICH", "param": "BULK_DENSITY_KG_M3"}
        value : float
            New parameter value
        """
        module = path.get("module", "")
        si_file = eval_dir / "openwq_in" / f"openWQ_MODULE_{module}.json"

        if not si_file.exists():
            logger.warning(f"Sorption file not found: {si_file}")
            return

        data, header = self._read_json_with_header(si_file)

        param_name = path.get("param", "")
        species = path.get("species")

        if species:
            # Per-species parameter
            if "SPECIES" in data and species in data["SPECIES"]:
                data["SPECIES"][species][param_name] = value
        else:
            # Global parameter
            data[param_name] = value

        self._write_json_with_header(si_file, data, header)

    # =========================================================================
    # Sediment Transport JSON
    # =========================================================================

    def _apply_sediment_param(self,
                              eval_dir: Path,
                              path: Dict,
                              value: float) -> None:
        """
        Modify a parameter in the sediment transport JSON file.

        Parameters
        ----------
        eval_dir : Path
            Evaluation directory
        path : Dict
            Parameter location: {"module": "HYPE_HBVSED", "param": "EROSION_INDEX"}
            or {"module": "HYPE_HBVSED", "param": "MONTHLY_EROSION_FACTOR", "index": 0}
        value : float
            New parameter value
        """
        module = path.get("module", "")
        ts_file = eval_dir / "openwq_in" / f"openWQ_MODULE_{module}.json"

        if not ts_file.exists():
            logger.warning(f"Sediment file not found: {ts_file}")
            return

        data, header = self._read_json_with_header(ts_file)

        param_name = path.get("param", "")
        index = path.get("index")

        if index is not None:
            # Monthly erosion factor or similar array
            if param_name in data:
                data[param_name][index] = value
        else:
            # Scalar parameter - may be in PARAMETER_DEFAULTS or directly
            if "PARAMETER_DEFAULTS" in data and param_name in data["PARAMETER_DEFAULTS"]:
                data["PARAMETER_DEFAULTS"][param_name] = value
            elif param_name in data:
                data[param_name] = value
            else:
                # Try spatially-varying format {"0": [...], "1": [...]}
                for key in data:
                    if isinstance(data[key], dict) and "1" in data[key]:
                        # Spatially-varying parameter table
                        # Update the default value (usually in row "1", last column)
                        if isinstance(data[key]["1"], list):
                            data[key]["1"][-1] = value
                        break

        self._write_json_with_header(ts_file, data, header)

    # =========================================================================
    # Transport Dissolved JSON
    # =========================================================================

    def _apply_transport_param(self,
                               eval_dir: Path,
                               path: Dict,
                               value: float) -> None:
        """
        Modify a parameter in the transport dissolved JSON file.
        """
        # Try to find the transport module file
        td_files = list((eval_dir / "openwq_in").glob("openWQ_MODULE_*TD*.json"))
        if not td_files:
            td_files = list((eval_dir / "openwq_in").glob("openWQ_MODULE_OPENWQ_NATIVE_TD_ADVDISP.json"))

        if not td_files:
            logger.warning("Transport module file not found")
            return

        td_file = td_files[0]
        data, header = self._read_json_with_header(td_file)

        param_name = path.get("param", "")
        if param_name in data:
            data[param_name] = value
        elif "TRANSPORT_CONFIGURATION" in data:
            data["TRANSPORT_CONFIGURATION"][param_name] = value

        self._write_json_with_header(td_file, data, header)

    # =========================================================================
    # Lateral Exchange JSON
    # =========================================================================

    def _apply_lateral_exchange_param(self,
                                      eval_dir: Path,
                                      path: Dict,
                                      value: float) -> None:
        """
        Modify K_val in the lateral exchange JSON file.
        """
        le_file = eval_dir / "openwq_in" / "openWQ_MODULE_NATIVE_LE_BOUNDMIX.json"

        if not le_file.exists():
            logger.warning(f"Lateral exchange file not found: {le_file}")
            return

        data, header = self._read_json_with_header(le_file)

        exchange_id = str(path.get("exchange_id", 1))
        param_name = path.get("param", "K_val")

        if exchange_id in data:
            data[exchange_id][param_name] = value

        self._write_json_with_header(le_file, data, header)

    # =========================================================================
    # Source/Sink: CSV Scaling
    # =========================================================================

    def _apply_ss_csv_scale(self,
                            eval_dir: Path,
                            path: Dict,
                            scale_factor: float) -> None:
        """
        Scale source/sink load values, covering both load formats:

        * ``DATA_FORMAT == "ASCII"`` — external CSV loads
          (``ss_method = "load_from_csv"``); the LOAD column of the
          referenced CSV is scaled.
        * ``DATA_FORMAT == "JSON"`` — inline loads such as the Copernicus
          LULC source/sink entries (static or dynamic coefficient, uniform
          or seasonal monthly distribution); the load value (index 9) of
          every DATA row in the matching entry is scaled in place.

        Parameters
        ----------
        eval_dir : Path
            Evaluation directory
        path : Dict
            {"species": "all"} or {"species": "N_ORG_active"}
        scale_factor : float
            Multiplier for load values
        """
        species_filter = path.get("species", "all")

        # Find the SS JSON to get CSV paths
        ss_files = list((eval_dir / "openwq_in").glob("openWQ_SS_*.json"))
        if not ss_files:
            logger.warning("No SS JSON file found")
            return

        ss_file = ss_files[0]
        data, header = self._read_json_with_header(ss_file)

        json_modified = False  # track whether any inline-JSON load was scaled

        # Process each SS entry
        for key, entry in data.items():
            if key == "METADATA":
                continue
            if not isinstance(entry, dict):
                continue

            chem_name = entry.get("CHEMICAL_NAME", "")
            data_format = entry.get("DATA_FORMAT", "")

            # Check species filter
            if species_filter != "all" and chem_name != species_filter:
                continue

            # ``scale_factor`` is either a constant multiplier OR a callable
            # month(1-12) -> multiplier (seasonal calibration). For JSON loads
            # the month is row[1]; for ASCII we apply the annual-mean factor
            # (seasonal scaling targets the Copernicus inline-JSON loads).
            _is_seasonal = callable(scale_factor)
            if data_format == "ASCII":
                # External CSV loads — scale the LOAD column in the file.
                filepath = entry.get("DATA", {}).get("FILEPATH", "")
                if filepath:
                    _sf = ((sum(scale_factor(m) for m in range(1, 13)) / 12.0)
                           if _is_seasonal else scale_factor)
                    self._scale_csv_loads(eval_dir, filepath, _sf)
            elif data_format == "JSON":
                # Inline JSON loads (e.g. Copernicus LULC source/sink entries,
                # both static- and dynamic-coefficient, and uniform or seasonal
                # monthly distributions).  DATA is a dict (or list) of rows;
                # each row is a list of the form
                #   [yyyy, mm, ..., load, "continuous"/"discrete", ...]
                # load value at index 9, month (MM) at index 1.  Scale in place
                # so the SS_COP_scale_<species> multiplier (constant or seasonal
                # per-month) actually takes effect.
                json_data = entry.get("DATA", {})
                if isinstance(json_data, dict):
                    rows = json_data.values()
                elif isinstance(json_data, list):
                    rows = json_data
                else:
                    rows = []
                n_scaled = 0
                for row in rows:
                    if isinstance(row, list) and len(row) > 9:
                        try:
                            if _is_seasonal:
                                try:
                                    _m = int(float(row[1]))
                                except (ValueError, TypeError, IndexError):
                                    _m = 1
                                _f = scale_factor(_m)
                            else:
                                _f = scale_factor
                            row[9] = float(row[9]) * _f
                            n_scaled += 1
                        except (ValueError, TypeError):
                            continue
                if n_scaled:
                    json_modified = True
                    logger.debug(
                        f"Scaled {n_scaled} inline-JSON SS load rows for "
                        f"'{chem_name}' ({'seasonal' if _is_seasonal else scale_factor})")

        # Persist inline-JSON edits back to the SS file.  (ASCII edits are
        # written to their own CSV files by _scale_csv_loads directly.)
        if json_modified:
            self._write_json_with_header(ss_file, data, header)

    def _apply_ss_loads(self, eval_dir: Path, ss_items: list) -> None:
        """Apply ALL source/sink load-scaling parameters together.

        Groups the params by species and builds one month(1-12)->multiplier
        function per species (constant / harmonic / per-month), then scales the
        loads. Grouping is required because the harmonic mode needs three
        sibling params (S0, amplitude, phase) combined, and per-month mode needs
        all twelve. Seasonal modes reshape the within-year load distribution, so
        they can shift PHASE - which a single constant scale never can."""
        by_sp = {}
        for param, value in ss_items:
            sp = (param.get("path") or {}).get("species", "all")
            by_sp.setdefault(sp, []).append((param, value))
        for sp, items in by_sp.items():
            factor = self._build_month_factor(items)
            self._apply_ss_csv_scale(eval_dir, {"species": sp}, factor)

    @staticmethod
    def _build_month_factor(items: list):
        """month(1-12) -> load multiplier for one species' SS params.

        Per-month params -> each month its own multiplier; else a harmonic
        S0*(1+A*cos(2pi*(m-phi)/12)) when amplitude+phase are present; else the
        constant S0 (returned as a plain float so the scalar fast-path is used
        and behaviour is unchanged when no seasonal params exist)."""
        import math
        S0 = 1.0
        amp = phase = None
        monthly = {}
        for param, value in items:
            try:
                v = float(value)
            except (TypeError, ValueError):
                continue
            ft = param.get("file_type")
            if ft == "ss_csv_scale":
                S0 = v
            elif ft == "ss_seasonal_amp":
                amp = v
            elif ft == "ss_seasonal_phase":
                phase = v
            elif ft == "ss_seasonal_month":
                mo = (param.get("path") or {}).get("month")
                if mo is not None:
                    monthly[int(mo)] = v
        if monthly:
            return lambda m: monthly.get(int(m), 1.0)
        if amp is not None and phase is not None:
            return lambda m: S0 * (1.0 + amp * math.cos(
                2.0 * math.pi * (m - phase) / 12.0))
        return S0   # constant -> scalar path (unchanged)

    def _scale_csv_loads(self, eval_dir: Path, csv_path: str, scale_factor: float) -> None:
        """Scale load values in a CSV file."""
        # The CSV path may be absolute or relative
        if csv_path.startswith("/code"):
            # Docker path - find corresponding local file
            # Extract the relative part after the mount point
            csv_path = csv_path.replace("/code/", "")

        # Try to find the CSV file
        csv_file = Path(csv_path)
        if not csv_file.exists():
            # Try relative to eval_dir
            csv_file = eval_dir / csv_path
        if not csv_file.exists():
            logger.warning(f"CSV file not found: {csv_path}")
            return

        # Read CSV with header rows
        with open(csv_file, 'r') as f:
            lines = f.readlines()

        # Find header row (contains YYYY)
        header_idx = None
        for i, line in enumerate(lines):
            if "YYYY" in line.upper():
                header_idx = i
                break

        if header_idx is None:
            logger.warning(f"Could not find header in CSV: {csv_file}")
            return

        # Read data portion
        header_lines = lines[:header_idx + 1]
        data_lines = lines[header_idx + 1:]

        # Parse header to find load column
        header = lines[header_idx].strip().split(',')
        load_idx = None
        for i, col in enumerate(header):
            if col.strip().lower() == 'load':
                load_idx = i
                break

        if load_idx is None:
            logger.warning(f"Could not find 'load' column in CSV: {csv_file}")
            return

        # Scale load values
        scaled_lines = []
        for line in data_lines:
            if not line.strip():
                scaled_lines.append(line)
                continue
            parts = line.strip().split(',')
            if len(parts) > load_idx:
                try:
                    original_load = float(parts[load_idx])
                    parts[load_idx] = str(original_load * scale_factor)
                except ValueError:
                    pass
            scaled_lines.append(','.join(parts) + '\n')

        # Write back
        with open(csv_file, 'w') as f:
            f.writelines(header_lines)
            f.writelines(scaled_lines)

        logger.debug(f"Scaled loads in {csv_file} by factor {scale_factor}")

    # =========================================================================
    # Source/Sink: Copernicus LULC Coefficients
    # =========================================================================

    def _apply_ss_copernicus_coeff(self,
                                   eval_dir: Path,
                                   path: Dict,
                                   value: float,
                                   initial_value: float = None) -> None:
        """
        Modify Copernicus LULC export coefficients in the SS JSON.

        The SS JSON has numbered entries, each with:
          - CHEMICAL_NAME: species name (e.g., "NH4_N", "TN", "TP")
          - DATA: dict of rows, each row is a list:
            [year, month, day, hour, "all", "all", lulc_idx, seg_id, subreach, load, "discrete"]

        The 'value' parameter is the NEW export coefficient. We compute a scale
        factor relative to the initial coefficient and apply it to all matching
        load entries.

        Parameters
        ----------
        eval_dir : Path
            Evaluation directory (with freshly-copied SS JSON)
        path : Dict
            {"lulc_class": int, "species": str}
            - lulc_class: LULC index used in the DATA arrays (position 6)
            - species: species name to match against CHEMICAL_NAME
              Use "TN"/"TP" to match aggregate species, or specific names
              like "NH4_N", "NO3_N" to match individual species.
        value : float
            New export coefficient value
        initial_value : float, optional
            Initial (baseline) export coefficient. If provided, the scale
            factor = value / initial_value is applied to all matching loads.
            If None, 'value' is treated as a direct scale factor.
        """
        lulc_class = path.get("lulc_class")
        species = path.get("species")

        # Compute scale factor
        if initial_value is not None and initial_value > 0:
            scale_factor = value / initial_value
        else:
            scale_factor = value

        # Find SS JSON
        ss_files = list((eval_dir / "openwq_in").glob("openWQ_SS_*copernicus*.json"))
        if not ss_files:
            ss_files = list((eval_dir / "openwq_in").glob("openWQ_SS_*.json"))
        if not ss_files:
            logger.warning("SS JSON not found for Copernicus coefficient modification")
            return

        ss_file = ss_files[0]
        data, header = self._read_json_with_header(ss_file)

        # Also load original (unmodified) SS from test_case_dir for base loads
        orig_data = None
        if self.test_case_dir is not None:
            original_ss = self.test_case_dir / "openwq_in" / ss_file.name
            if original_ss.exists() and original_ss != ss_file:
                orig_data, _ = self._read_json_with_header(original_ss)

        modified = False
        n_modified = 0

        for key, entry in data.items():
            if key == "METADATA" or not isinstance(entry, dict):
                continue

            chem_name = entry.get("CHEMICAL_NAME", "")

            # Match species: exact match, or "TN" matches N-species, "TP" matches P-species
            if not self._ss_species_match(chem_name, species):
                continue

            if "DATA" not in entry:
                continue

            # Get original entry's DATA for base load values
            orig_entry_data = None
            if orig_data and key in orig_data and "DATA" in orig_data[key]:
                orig_entry_data = orig_data[key]["DATA"]

            for row_key, row in entry["DATA"].items():
                if not isinstance(row, list) or len(row) < 10:
                    continue

                # Check LULC class at position 6
                row_lulc = row[6]
                if row_lulc != lulc_class and str(row_lulc) != str(lulc_class):
                    continue

                # Get original base load value
                if orig_entry_data and row_key in orig_entry_data:
                    orig_row = orig_entry_data[row_key]
                    base_load = orig_row[9] if isinstance(orig_row, list) and len(orig_row) > 9 else row[9]
                else:
                    base_load = row[9]

                # Apply scale factor to base load
                row[9] = base_load * scale_factor
                n_modified += 1
                modified = True

        if modified:
            self._write_json_with_header(ss_file, data, header)
            logger.debug(f"Modified {n_modified} SS load entries for LULC={lulc_class}, species={species} "
                        f"(scale_factor={scale_factor:.4f})")
        else:
            logger.warning(f"Could not find SS entries for LULC={lulc_class}, species={species}")

    @staticmethod
    def _ss_species_match(chem_name: str, target_species: str) -> bool:
        """
        Check if a SS entry's CHEMICAL_NAME matches the target species.

        Supports exact match and aggregate matching:
          - "TN" matches TN, NH4_N, NO3_N, NO2_N, DON, PON (all N species)
          - "TP" matches TP, PO4_P (all P species)
          - Exact names match directly (e.g., "NH4_N" matches "NH4_N")
        """
        if not chem_name or not target_species:
            return False

        chem_upper = chem_name.upper()
        target_upper = target_species.upper()

        # Exact match
        if chem_upper == target_upper:
            return True

        # Aggregate TN matches all nitrogen species
        if target_upper == "TN":
            n_species = {"TN", "NH4_N", "NO3_N", "NO2_N", "DON", "PON",
                        "NH4-N", "NO3-N", "NO2-N", "N_ORG_ACTIVE", "N_ORG_PASSIVE"}
            return chem_upper in n_species

        # Aggregate TP matches all phosphorus species
        if target_upper == "TP":
            p_species = {"TP", "PO4_P", "PO4-P", "P_ORG_ACTIVE", "P_ORG_PASSIVE",
                        "P_INORG", "DIP", "DOP"}
            return chem_upper in p_species

        return False

    # =========================================================================
    # Source/Sink: Climate Parameters
    # =========================================================================

    def _apply_ss_climate_param(self,
                                eval_dir: Path,
                                path: Union[str, Dict],
                                value: float) -> None:
        """
        [SUPERSEDED] Dynamic SS climate-response parameters
        (precip_scaling_power, Q10_biological, T_reference) are now applied
        at *generation time*: setup_working_directory() bakes their
        calibrated values into the per-evaluation model config via
        _build_eval_config(), so Gen_Input_Driver/Gen_SS_Driver regenerate
        the SS JSON with the correct climate-weighted monthly distribution.

        This post-hoc hook is therefore no longer wired into
        apply_parameters() and is kept only for backward compatibility.
        """
        logger.debug(
            f"_apply_ss_climate_param({path}) is a no-op; climate params are "
            f"baked into the SS JSON at generation time.")

    # =========================================================================
    # Source/Sink: ML Model Parameters
    # =========================================================================

    def _apply_ss_ml_param(self,
                           eval_dir: Path,
                           path: str,
                           value: float) -> None:
        """
        [SUPERSEDED] ML source/sink hyperparameters (n_estimators, max_depth)
        are now applied at *generation time*: setup_working_directory() bakes
        their calibrated values into the per-evaluation model config via
        _build_eval_config(), so Gen_Input_Driver/Gen_MLmodel_SS retrains the
        model and regenerates the SS JSON with the chosen hyperparameters.

        This post-hoc hook is therefore no longer wired into
        apply_parameters() and is kept only for backward compatibility.
        """
        logger.debug(
            f"_apply_ss_ml_param({path}) is a no-op; ML hyperparameters are "
            f"baked into the SS JSON at generation time.")

    # =========================================================================
    # Utility Methods
    # =========================================================================

    def _read_json_with_header(self, filepath: Path) -> Tuple[Dict, str]:
        """Read a JSON file that may have comment header lines."""
        with open(filepath, 'r') as f:
            content = f.read()

        # Find where JSON starts (first '{')
        json_start = content.find('{')
        if json_start > 0:
            header = content[:json_start]
            json_content = content[json_start:]
        else:
            header = ""
            json_content = content

        # Parse JSON
        data = json.loads(json_content)
        return data, header

    def _write_json_with_header(self, filepath: Path, data: Dict, header: str) -> None:
        """Write a JSON file with optional comment header."""
        json_string = json.dumps(data, indent=4)

        # Compact arrays onto single lines
        def compact_list(match):
            list_content = match.group(0)
            compact = re.sub(r'\[\s+', '[', list_content)
            compact = re.sub(r'\s+\]', ']', compact)
            compact = re.sub(r',\s+', ', ', compact)
            return compact

        pattern = r'\[[\s\n]*(?:[^\[\]]*?)[\s\n]*\]'
        prev_string = ""
        while prev_string != json_string:
            prev_string = json_string
            json_string = re.sub(pattern, compact_list, json_string)

        with open(filepath, 'w') as f:
            f.write(header)
            f.write(json_string)
            if not json_string.endswith('\n'):
                f.write('\n')

    @staticmethod
    def transform_to_real(value: float, transform: str) -> float:
        """Transform from optimization space to real parameter space."""
        if transform == "log":
            return 10 ** value
        return value

    @staticmethod
    def transform_to_opt(value: float, transform: str) -> float:
        """Transform from real parameter space to optimization space."""
        if transform == "log":
            return np.log10(value)
        return value

    @staticmethod
    def get_bounds_in_opt_space(bounds: Tuple[float, float],
                                transform: str) -> Tuple[float, float]:
        """Convert bounds to optimization space."""
        if transform == "log":
            return (np.log10(bounds[0]), np.log10(bounds[1]))
        return bounds
