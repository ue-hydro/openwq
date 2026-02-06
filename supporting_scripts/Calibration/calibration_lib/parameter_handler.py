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
"""

import json
import os
import re
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
                 base_model_config_dir: str,
                 test_case_dir: str,
                 calibration_work_dir: str,
                 running_on_docker: bool = True):
        """
        Initialize with paths to base configuration and working directory.

        Parameters
        ----------
        base_model_config_dir : str
            Path to the Model_Config directory containing template_model_config.py
        test_case_dir : str
            Path to the test case directory (containing openwq_in, mizuroute_in, etc.)
        calibration_work_dir : str
            Path to working directory for calibration evaluations
        running_on_docker : bool
            Whether the model runs in Docker (affects path handling)
        """
        self.base_model_config_dir = Path(base_model_config_dir)
        self.test_case_dir = Path(test_case_dir)
        self.calibration_work_dir = Path(calibration_work_dir)
        self.running_on_docker = running_on_docker

        # Ensure directories exist
        self.calibration_work_dir.mkdir(parents=True, exist_ok=True)
        (self.calibration_work_dir / "evaluations").mkdir(exist_ok=True)

    def setup_working_directory(self, eval_id: int) -> Path:
        """
        Create isolated working directory for a single model evaluation.
        Copies all necessary files from test case directory.

        Parameters
        ----------
        eval_id : int
            Unique evaluation identifier

        Returns
        -------
        Path
            Path to the evaluation working directory
        """
        eval_dir = self.calibration_work_dir / "evaluations" / f"eval_{eval_id:04d}"

        # Clean up if exists
        if eval_dir.exists():
            shutil.rmtree(eval_dir)

        eval_dir.mkdir(parents=True)

        # Copy test case structure
        # Copy openwq_in directory
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

        logger.debug(f"Created evaluation directory: {eval_dir}")
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
        for i, param in enumerate(parameters):
            value = values[i]
            file_type = param["file_type"]
            path = param["path"]
            dtype = param.get("dtype", "float")

            # Convert to appropriate type
            if dtype == "int":
                value = int(round(value))

            logger.debug(f"Applying {param['name']} = {value} (type={file_type})")

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
            elif file_type == "ss_csv_scale":
                self._apply_ss_csv_scale(eval_dir, path, value)
            elif file_type in ("ss_copernicus_static", "ss_copernicus_dynamic"):
                self._apply_ss_copernicus_coeff(eval_dir, path, value)
            elif file_type == "ss_climate_param":
                self._apply_ss_climate_param(eval_dir, path, value)
            elif file_type == "ss_ml_param":
                self._apply_ss_ml_param(eval_dir, path, value)
            else:
                logger.warning(f"Unknown file_type: {file_type} for parameter {param['name']}")

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
        Scale source/sink load values in CSV files.

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

            if data_format == "ASCII":
                filepath = entry.get("DATA", {}).get("FILEPATH", "")
                if filepath:
                    self._scale_csv_loads(eval_dir, filepath, scale_factor)

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
                                   value: float) -> None:
        """
        Modify Copernicus LULC export coefficients.

        This requires modifying the template config and regenerating the SS JSON.
        For now, we modify the generated JSON directly if possible.
        """
        lulc_class = path.get("lulc_class")
        species = path.get("species")

        # Find SS JSON
        ss_files = list((eval_dir / "openwq_in").glob("openWQ_SS_*copernicus*.json"))
        if not ss_files:
            logger.warning("Copernicus SS JSON not found - may need template regeneration")
            return

        ss_file = ss_files[0]
        data, header = self._read_json_with_header(ss_file)

        # The Copernicus SS JSON structure may vary based on the method
        # Try to find and modify the load coefficient
        modified = False
        for key, entry in data.items():
            if key == "METADATA":
                continue
            if isinstance(entry, dict):
                # Check if this entry matches our LULC class and species
                if entry.get("LULC_CLASS") == lulc_class or \
                   str(entry.get("LULC_CLASS")) == str(lulc_class):
                    if "LOAD_COEFFICIENTS" in entry:
                        if species in entry["LOAD_COEFFICIENTS"]:
                            entry["LOAD_COEFFICIENTS"][species] = value
                            modified = True

        if modified:
            self._write_json_with_header(ss_file, data, header)
        else:
            logger.warning(f"Could not find Copernicus coefficient for LULC={lulc_class}, species={species}")

    # =========================================================================
    # Source/Sink: Climate Parameters
    # =========================================================================

    def _apply_ss_climate_param(self,
                                eval_dir: Path,
                                path: Union[str, Dict],
                                value: float) -> None:
        """
        Modify climate response parameters for dynamic SS method.

        These affect the template config - for now we note this limitation.
        """
        # Climate parameters like ss_climate_precip_scaling_power are in the Python template
        # To properly handle these, we would need to:
        # 1. Copy template_model_config.py to eval_dir
        # 2. Modify the parameter value
        # 3. Run the generator to create new JSONs
        #
        # For efficiency, we log a warning and skip for now
        # This should be implemented via the template regeneration approach
        logger.warning(f"Climate parameter modification ({path}) requires template regeneration - not yet implemented")

    # =========================================================================
    # Source/Sink: ML Model Parameters
    # =========================================================================

    def _apply_ss_ml_param(self,
                           eval_dir: Path,
                           path: str,
                           value: float) -> None:
        """
        Modify ML model hyperparameters.

        Similar to climate params, these are in the Python template.
        """
        logger.warning(f"ML parameter modification ({path}) requires template regeneration - not yet implemented")

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
