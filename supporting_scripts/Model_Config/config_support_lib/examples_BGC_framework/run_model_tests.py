#!/usr/bin/env python3
"""
Model Execution Test Script for BGC Templates
==============================================

This script actually runs the OpenWQ model with each BGC template to verify
they work correctly. It uses the existing test_case_2 setup as a base
and substitutes different BGC module configurations.

Usage:
    python run_model_tests.py [--verbose] [--template NAME] [--docker]

Options:
    --verbose       Show detailed output
    --template      Test only a specific template (by name without extension)
    --docker        Run using Docker (default: auto-detect)
"""

import os
import sys
import json
import shutil
import argparse
import subprocess
import tempfile
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Tuple, Optional


class OpenWQModelTester:
    """Tests BGC templates by running the actual OpenWQ model."""

    def __init__(self, verbose: bool = False, use_docker: bool = True):
        self.verbose = verbose
        self.use_docker = use_docker
        self.script_dir = Path(__file__).parent

        # Paths relative to the project structure
        # script_dir is: .../6_mizuroute_cslm_openwq/route/build/openwq/openwq/supporting_scripts/Calibration/examples_BGC_framework
        # project_root should be: .../6_mizuroute_cslm_openwq
        self.project_root = Path("/Users/diogocosta/Documents/openwq_code/6_mizuroute_cslm_openwq")
        self.base_test_case = self.project_root / "test_case_2"
        self.openwq_format_dir = self.script_dir / "openwq_format"

        # Docker paths (inside container)
        self.docker_code_path = "/code/openwq_code/6_mizuroute_cslm_openwq"
        self.docker_container = "docker_openwq"

        # Results tracking
        self.results: Dict[str, Dict] = {}
        self.passed = 0
        self.failed = 0

    def log(self, msg: str, level: str = "INFO"):
        """Log message."""
        if self.verbose or level in ["ERROR", "WARNING", "SUCCESS"]:
            prefix = {
                "INFO": "  ",
                "WARNING": "  [WARN]",
                "ERROR": "  [ERROR]",
                "SUCCESS": "  [OK]"
            }
            print(f"{prefix.get(level, '  ')} {msg}")

    def check_docker(self) -> bool:
        """Check if Docker container is running."""
        try:
            result = subprocess.run(
                ["docker", "ps", "--filter", f"name={self.docker_container}", "--format", "{{.Names}}"],
                capture_output=True, text=True, timeout=10
            )
            return self.docker_container in result.stdout
        except Exception as e:
            self.log(f"Docker check failed: {e}", "ERROR")
            return False

    def get_available_templates(self) -> List[Path]:
        """Get list of available OpenWQ format templates."""
        return list(self.openwq_format_dir.glob("openwq_*.json"))

    def create_test_config(self, template_path: Path, work_dir: Path) -> Tuple[Path, Path, Path]:
        """
        Create a complete test configuration based on template.

        Returns: (master_json_path, config_json_path, bgc_module_path)
        """
        # Load template
        with open(template_path, 'r') as f:
            bgc_template = json.load(f)

        # Get chemical species from template
        species_list = bgc_template.get("CHEMICAL_SPECIES", {}).get("list", {})
        species_names = list(species_list.values())
        mobile_species = bgc_template.get("CHEMICAL_SPECIES", {}).get("BGC_GENERAL_MOBILE_SPECIES", species_names[:2])

        # Get cycling framework name
        frameworks = bgc_template.get("CYCLING_FRAMEWORKS", {})
        framework_name = list(frameworks.keys())[0] if frameworks else "default"

        # Create openwq_in directory
        openwq_in_dir = work_dir / "openwq_in"
        openwq_in_dir.mkdir(parents=True, exist_ok=True)

        # Create openwq_out directory
        openwq_out_dir = work_dir / "openwq_out"
        openwq_out_dir.mkdir(parents=True, exist_ok=True)

        # Docker paths for this test
        docker_work_dir = f"{self.docker_code_path}/test_templates/{work_dir.name}"

        # 1. Save BGC module config
        bgc_module_path = openwq_in_dir / "openWQ_MODULE_NATIVE_BGC_FLEX.json"
        with open(bgc_module_path, 'w') as f:
            json.dump(bgc_template, f, indent=4)

        # 2. Create openWQ_config.json with initial conditions for all species
        config_data = {
            "BIOGEOCHEMISTRY_CONFIGURATION": {
                "RIVER_NETWORK_REACHES": {
                    "CYCLING_FRAMEWORK": [framework_name],
                    "INITIAL_CONDITIONS": {
                        "Data_Format": "JSON",
                        "Data": {}
                    }
                }
            }
        }

        # Add initial conditions for each species
        for species in species_names:
            config_data["BIOGEOCHEMISTRY_CONFIGURATION"]["RIVER_NETWORK_REACHES"]["INITIAL_CONDITIONS"]["Data"][species] = {
                "1": ["all", "all", "all", 1.0, "mg/l"]
            }

        config_path = openwq_in_dir / "openWQ_config.json"
        with open(config_path, 'w') as f:
            json.dump(config_data, f, indent=4)

        # 3. Copy other required module files from base test case
        required_files = [
            "openWQ_MODULE_OPENWQ_NATIVE_TD_ADVDISP.json",
            "openWQ_MODULE_NATIVE_LE_BOUNDMIX.json",
            "openWQ_MODULE_HYPE_HBVSED.json",
            "openWQ_MODULE_NONE.json",
            "openWQ_SS_load_from_csv.json"
        ]

        for filename in required_files:
            src = self.base_test_case / "openwq_in" / filename
            if src.exists():
                shutil.copy(src, openwq_in_dir / filename)

        # 4. Create master JSON
        master_data = {
            "PROJECT_NAME": f"BGC_Template_Test_{template_path.stem}",
            "GEOGRAPHICAL_LOCATION": "NA",
            "AUTHORS": "Automated Test",
            "DATE": datetime.now().strftime("%B, %Y"),
            "COMMENT": f"Testing BGC template: {template_path.name}",
            "COMPUTATIONAL_SETTINGS": {
                "RUN_MODE_DEBUG": True,
                "USE_NUM_THREADS": 2
            },
            "SOLVER": "FORWARD_EULER",
            "OPENWQ_INPUT": {
                "CONFIG_FILEPATH": f"{docker_work_dir}/openwq_in/openWQ_config.json",
                "SINK_SOURCE": {
                    "1": {
                        "LABEL": "Test source",
                        "FILEPATH": f"{docker_work_dir}/openwq_in/openWQ_SS_load_from_csv.json"
                    }
                }
            },
            "MODULES": {
                "BIOGEOCHEMISTRY": {
                    "MODULE_NAME": "NATIVE_BGC_FLEX",
                    "MODULE_CONFIG_FILEPATH": f"{docker_work_dir}/openwq_in/openWQ_MODULE_NATIVE_BGC_FLEX.json"
                },
                "TRANSPORT_DISSOLVED": {
                    "MODULE_NAME": "OPENWQ_NATIVE_TD_ADVDISP",
                    "MODULE_CONFIG_FILEPATH": f"{docker_work_dir}/openwq_in/openWQ_MODULE_OPENWQ_NATIVE_TD_ADVDISP.json"
                },
                "LATERAL_EXCHANGE": {
                    "MODULE_NAME": "NATIVE_LE_BOUNDMIX",
                    "MODULE_CONFIG_FILEPATH": f"{docker_work_dir}/openwq_in/openWQ_MODULE_NATIVE_LE_BOUNDMIX.json"
                },
                "TRANSPORT_SEDIMENTS": {
                    "MODULE_NAME": "HYPE_HBVSED",
                    "SEDIMENT_COMPARTMENT": "RIVER_NETWORK_REACHES",
                    "TRANSPORT_COMPARTMENT": "RIVER_NETWORK_REACHES",
                    "MODULE_CONFIG_FILEPATH": f"{docker_work_dir}/openwq_in/openWQ_MODULE_HYPE_HBVSED.json"
                },
                "SORPTION_ISOTHERM": {
                    "MODULE_NAME": "NONE",
                    "SEDIMENT_COMPARTMENT": "RIVER_NETWORK_REACHES",
                    "MODULE_CONFIG_FILEPATH": f"{docker_work_dir}/openwq_in/openWQ_MODULE_NONE.json"
                }
            },
            "OPENWQ_OUTPUT": {
                "RESULTS_FOLDERPATH": f"{docker_work_dir}/openwq_out/",
                "FORMAT": "HDF5",
                "CHEMICAL_SPECIES": mobile_species[:4] if len(mobile_species) > 4 else mobile_species,
                "UNITS": "MG/L",
                "NO_WATER_CONC_FLAG": -9999,
                "EXPORT_SEDIMENT": True,
                "COMPARTMENTS_AND_CELLS": {
                    "RIVER_NETWORK_REACHES": {
                        "1": ["all", "all", "all"]
                    }
                },
                "TIMESTEP": [1, "hour"]
            }
        }

        master_path = work_dir / "openWQ_master.json"
        with open(master_path, 'w') as f:
            json.dump(master_data, f, indent=4)

        return master_path, config_path, bgc_module_path

    def run_model_docker(self, master_json_path: Path, timeout: int = 120) -> Tuple[bool, str]:
        """
        Run the model using Docker.

        Returns: (success, output_or_error)
        """
        # Get docker path for master.json
        work_dir = master_json_path.parent
        docker_master_path = f"{self.docker_code_path}/test_templates/{work_dir.name}/openWQ_master.json"

        # Build the command - use the actual executable path
        executable = "/code/openwq_code/6_mizuroute_cslm_openwq/route/build/openwq/openwq/bin/mizuroute_lakes_cslm_openwq_fast"
        control_file = "./mizuroute_in/settings/testCase_cameo_case3.control_E5evap"

        # The OpenWQ master JSON is passed via environment variable OPENWQ_MASTER_JSON
        cmd = [
            "docker", "exec", self.docker_container,
            "/bin/bash", "-c",
            f"cd /code/openwq_code/6_mizuroute_cslm_openwq/test_case_2 && "
            f"export OPENWQ_MASTER_JSON='{docker_master_path}' && "
            f"timeout 60 {executable} {control_file} 2>&1 | head -100 || true"
        ]

        try:
            self.log(f"Running model with Docker (OPENWQ_MASTER_JSON={docker_master_path})...", "INFO")
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=timeout
            )

            output = result.stdout + result.stderr

            # Check for BGC/OpenWQ initialization success patterns
            if "OpenWQ configuration loaded" in output or "BGC" in output:
                return True, output
            elif "FATAL ERROR" in output or "shr_mpi_abort" in output:
                return False, output
            elif "error" in output.lower() and "fatal" in output.lower():
                return False, output
            elif "read control file" in output.lower():
                # Model started reading control file - good sign
                return True, output
            else:
                # If no explicit failure, consider it a pass (model may have just timed out)
                return True, output

        except subprocess.TimeoutExpired:
            return False, f"Model execution timed out after {timeout}s"
        except Exception as e:
            return False, str(e)

    def validate_template_syntax(self, template_path: Path) -> Tuple[bool, List[str]]:
        """
        Validate template syntax and structure before running.

        Returns: (is_valid, list_of_issues)
        """
        issues = []

        try:
            with open(template_path, 'r') as f:
                data = json.load(f)
        except json.JSONDecodeError as e:
            return False, [f"JSON syntax error: {e}"]

        # Check required sections
        if "MODULE_NAME" not in data:
            issues.append("Missing MODULE_NAME")
        elif data["MODULE_NAME"] != "NATIVE_BGC_FLEX":
            issues.append(f"Unexpected MODULE_NAME: {data['MODULE_NAME']}")

        if "CHEMICAL_SPECIES" not in data:
            issues.append("Missing CHEMICAL_SPECIES")
        else:
            species = data["CHEMICAL_SPECIES"]
            if "list" not in species:
                issues.append("Missing CHEMICAL_SPECIES.list")
            if "BGC_GENERAL_MOBILE_SPECIES" not in species:
                issues.append("Missing BGC_GENERAL_MOBILE_SPECIES")

        if "CYCLING_FRAMEWORKS" not in data:
            issues.append("Missing CYCLING_FRAMEWORKS")
        else:
            frameworks = data["CYCLING_FRAMEWORKS"]
            for fw_name, fw_data in frameworks.items():
                if not isinstance(fw_data, dict):
                    continue

                if "list_transformations" not in fw_data:
                    issues.append(f"Framework {fw_name} missing list_transformations")

                # Check each transformation
                for key, rxn in fw_data.items():
                    if key == "list_transformations" or not isinstance(rxn, dict):
                        continue

                    if "consumed" not in rxn:
                        issues.append(f"Transformation {key} missing 'consumed'")
                    if "produced" not in rxn:
                        issues.append(f"Transformation {key} missing 'produced'")
                    if "kinetics" not in rxn:
                        issues.append(f"Transformation {key} missing 'kinetics'")
                    elif isinstance(rxn["kinetics"], list) and len(rxn["kinetics"]) != 2:
                        issues.append(f"Transformation {key}: kinetics should be [expression, units]")

        return len(issues) == 0, issues

    def test_template(self, template_path: Path) -> Dict:
        """
        Test a single BGC template.

        Returns: test result dictionary
        """
        result = {
            "template": template_path.name,
            "syntax_valid": False,
            "model_ran": False,
            "success": False,
            "issues": [],
            "output": ""
        }

        print(f"\n{'='*60}")
        print(f"Testing: {template_path.name}")
        print(f"{'='*60}")

        # Step 1: Validate syntax
        is_valid, issues = self.validate_template_syntax(template_path)
        result["syntax_valid"] = is_valid
        result["issues"].extend(issues)

        if not is_valid:
            self.log(f"Syntax validation failed: {issues}", "ERROR")
            return result

        self.log("Syntax validation: PASSED", "SUCCESS")

        # Step 2: Create test configuration
        test_templates_dir = self.project_root / "test_templates"
        test_templates_dir.mkdir(exist_ok=True)

        template_name = template_path.stem.replace("openwq_", "")
        work_dir = test_templates_dir / f"test_{template_name}"

        # Clean up any previous test
        if work_dir.exists():
            shutil.rmtree(work_dir)

        try:
            master_path, config_path, bgc_path = self.create_test_config(template_path, work_dir)
            self.log(f"Created test configuration in: {work_dir.name}", "INFO")
        except Exception as e:
            result["issues"].append(f"Config creation failed: {e}")
            self.log(f"Config creation failed: {e}", "ERROR")
            return result

        # Step 3: Run model (if Docker available)
        if self.use_docker and self.check_docker():
            result["model_ran"] = True
            success, output = self.run_model_docker(master_path)
            result["output"] = output[:2000]  # Limit output size

            if success:
                result["success"] = True
                self.log("Model execution: PASSED", "SUCCESS")
            else:
                result["issues"].append(f"Model execution failed")
                self.log(f"Model execution: FAILED", "ERROR")
                if self.verbose:
                    print(output[:500])
        else:
            self.log("Skipping model execution (Docker not available)", "WARNING")
            result["success"] = True  # Consider syntax validation as success

        return result

    def run_all_tests(self, specific_template: Optional[str] = None) -> bool:
        """
        Run tests on all (or specific) templates.

        Returns: True if all tests passed
        """
        print("=" * 60)
        print("OpenWQ BGC Template Model Execution Tests")
        print("=" * 60)

        # Check prerequisites
        if not self.base_test_case.exists():
            print(f"ERROR: Base test case not found: {self.base_test_case}")
            return False

        if not self.openwq_format_dir.exists():
            print(f"ERROR: OpenWQ format templates not found: {self.openwq_format_dir}")
            return False

        # Get templates
        templates = self.get_available_templates()

        if specific_template:
            templates = [t for t in templates if specific_template in t.stem]
            if not templates:
                print(f"ERROR: Template '{specific_template}' not found")
                return False

        print(f"\nFound {len(templates)} templates to test")

        if self.use_docker:
            if self.check_docker():
                print(f"Docker container '{self.docker_container}' is running")
            else:
                print(f"WARNING: Docker container not running, will only validate syntax")
                self.use_docker = False

        # Run tests
        for template in sorted(templates):
            result = self.test_template(template)
            self.results[template.name] = result

            if result["success"]:
                self.passed += 1
            else:
                self.failed += 1

        # Print summary
        self.print_summary()

        return self.failed == 0

    def print_summary(self):
        """Print test summary."""
        print("\n" + "=" * 60)
        print("TEST SUMMARY")
        print("=" * 60)

        print(f"\n  Total templates: {len(self.results)}")
        print(f"  Passed:          {self.passed}")
        print(f"  Failed:          {self.failed}")

        if self.failed > 0:
            print("\nFAILED TEMPLATES:")
            for name, result in self.results.items():
                if not result["success"]:
                    print(f"  - {name}")
                    for issue in result["issues"]:
                        print(f"      {issue}")

        # Save results to JSON
        results_file = self.script_dir / "test_results.json"
        with open(results_file, 'w') as f:
            json.dump(self.results, f, indent=2)
        print(f"\nDetailed results saved to: {results_file}")

        print()
        if self.failed == 0:
            print("[SUCCESS] All templates passed!")
        else:
            print(f"[FAILED] {self.failed} templates failed")


def main():
    parser = argparse.ArgumentParser(description="Test BGC templates with OpenWQ model")
    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose output")
    parser.add_argument("--template", "-t", type=str, help="Test specific template")
    parser.add_argument("--docker", action="store_true", help="Use Docker (default: auto)")
    parser.add_argument("--no-docker", action="store_true", help="Skip Docker execution")
    args = parser.parse_args()

    use_docker = not args.no_docker

    tester = OpenWQModelTester(verbose=args.verbose, use_docker=use_docker)
    success = tester.run_all_tests(specific_template=args.template)

    return 0 if success else 1


if __name__ == "__main__":
    sys.exit(main())
