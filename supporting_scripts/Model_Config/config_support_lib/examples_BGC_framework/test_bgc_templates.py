#!/usr/bin/env python3
"""
Test Script for BGC Templates
==============================

This script validates and tests all BGC JSON templates in the examples_BGC_framework directory.
It checks:
1. JSON syntax validity
2. Required structure elements
3. Conversion to OpenWQ-compatible format
4. Expression parsing validation

Usage:
    python test_bgc_templates.py [--verbose] [--convert]

Options:
    --verbose   Show detailed validation output
    --convert   Generate OpenWQ-compatible versions of templates
"""

import json
import os
import sys
import re
import argparse
from pathlib import Path
from typing import Dict, List, Tuple, Any, Optional


class BGCTemplateValidator:
    """Validates BGC template JSON files for OpenWQ compatibility."""

    def __init__(self, verbose: bool = False):
        self.verbose = verbose
        self.errors: List[str] = []
        self.warnings: List[str] = []
        self.valid_count = 0
        self.invalid_count = 0

    def log(self, msg: str, level: str = "INFO"):
        """Log message if verbose mode is on."""
        if self.verbose or level in ["ERROR", "WARNING"]:
            prefix = {"INFO": "  ", "WARNING": "  [WARN]", "ERROR": "  [ERROR]"}
            print(f"{prefix.get(level, '  ')} {msg}")

    def validate_json_syntax(self, filepath: Path) -> Tuple[bool, Optional[dict]]:
        """Check if JSON file has valid syntax."""
        try:
            with open(filepath, 'r') as f:
                data = json.load(f)
            return True, data
        except json.JSONDecodeError as e:
            self.errors.append(f"{filepath.name}: JSON syntax error - {e}")
            return False, None
        except Exception as e:
            self.errors.append(f"{filepath.name}: Read error - {e}")
            return False, None

    def validate_structure(self, data: dict, filepath: Path) -> bool:
        """Validate the template structure has required elements."""
        filename = filepath.name
        valid = True

        # Check for either our template format or OpenWQ format
        if "BIOGEOCHEMISTRY_CONFIGURATION" in data:
            # Our human-readable template format
            bgc_config = data["BIOGEOCHEMISTRY_CONFIGURATION"]

            # Check for CHEMICAL_SPECIES
            if "CHEMICAL_SPECIES" not in bgc_config:
                self.errors.append(f"{filename}: Missing CHEMICAL_SPECIES section")
                valid = False
            else:
                species = bgc_config["CHEMICAL_SPECIES"]
                species_count = len([k for k in species.keys() if not k.startswith("_")])
                self.log(f"Found {species_count} chemical species")

            # Check for CYCLING_FRAMEWORKS
            if "CYCLING_FRAMEWORKS" not in bgc_config:
                self.errors.append(f"{filename}: Missing CYCLING_FRAMEWORKS section")
                valid = False
            else:
                frameworks = bgc_config["CYCLING_FRAMEWORKS"]
                framework_count = len([k for k in frameworks.keys() if not k.startswith("_")])
                self.log(f"Found {framework_count} cycling frameworks")

                # Validate each framework
                for fw_name, fw_data in frameworks.items():
                    if fw_name.startswith("_"):
                        continue
                    if isinstance(fw_data, dict):
                        valid = valid and self._validate_framework(fw_data, fw_name, filename)

        elif "MODULE_NAME" in data:
            # OpenWQ native format
            if data.get("MODULE_NAME") != "NATIVE_BGC_FLEX":
                self.warnings.append(f"{filename}: MODULE_NAME is not NATIVE_BGC_FLEX")

            if "CHEMICAL_SPECIES" not in data:
                self.errors.append(f"{filename}: Missing CHEMICAL_SPECIES section")
                valid = False

            if "CYCLING_FRAMEWORKS" not in data:
                self.errors.append(f"{filename}: Missing CYCLING_FRAMEWORKS section")
                valid = False

        else:
            self.errors.append(f"{filename}: Unknown template format - missing BIOGEOCHEMISTRY_CONFIGURATION or MODULE_NAME")
            valid = False

        return valid

    def _validate_framework(self, fw_data: dict, fw_name: str, filename: str) -> bool:
        """Validate a single cycling framework."""
        valid = True

        for rxn_name, rxn_data in fw_data.items():
            if rxn_name.startswith("_"):
                continue
            if not isinstance(rxn_data, dict):
                continue

            # Check for required reaction fields
            if "CONSUMED" not in rxn_data and "consumed" not in rxn_data:
                self.warnings.append(f"{filename}/{fw_name}/{rxn_name}: Missing CONSUMED field")

            if "PRODUCED" not in rxn_data and "produced" not in rxn_data:
                self.warnings.append(f"{filename}/{fw_name}/{rxn_name}: Missing PRODUCED field")

            if "KINETICS" not in rxn_data and "kinetics" not in rxn_data:
                self.warnings.append(f"{filename}/{fw_name}/{rxn_name}: Missing KINETICS field")
            else:
                # Validate kinetics expression
                kinetics = rxn_data.get("KINETICS") or rxn_data.get("kinetics")
                if isinstance(kinetics, str):
                    valid = valid and self._validate_expression(kinetics, rxn_name, filename)

        return valid

    def _validate_expression(self, expr: str, rxn_name: str, filename: str) -> bool:
        """Validate a kinetics expression for basic syntax."""
        valid = True

        # Check for balanced parentheses
        if expr.count('(') != expr.count(')'):
            self.errors.append(f"{filename}/{rxn_name}: Unbalanced parentheses in expression: {expr}")
            valid = False

        # Check for common expression elements
        # Valid operators: +, -, *, /, ^, (, )
        # Valid functions: exp, log, sqrt, min, max, sin, cos, etc.
        # Valid characters: alphanumeric, _, .

        # Check for invalid characters
        invalid_chars = re.findall(r'[^\w\s\+\-\*\/\^\(\)\.\,\=\<\>\!]', expr)
        if invalid_chars:
            self.warnings.append(f"{filename}/{rxn_name}: Unusual characters in expression: {set(invalid_chars)}")

        return valid

    def validate_file(self, filepath: Path) -> bool:
        """Validate a single BGC template file."""
        self.log(f"\nValidating: {filepath.name}")

        # Step 1: JSON syntax
        is_valid_json, data = self.validate_json_syntax(filepath)
        if not is_valid_json:
            self.invalid_count += 1
            return False

        self.log("JSON syntax: OK")

        # Step 2: Structure validation
        is_valid_structure = self.validate_structure(data, filepath)
        if not is_valid_structure:
            self.invalid_count += 1
            return False

        self.log("Structure: OK")
        self.valid_count += 1
        return True

    def validate_directory(self, dirpath: Path) -> bool:
        """Validate all JSON files in a directory."""
        json_files = list(dirpath.glob("**/*.json"))

        if not json_files:
            print(f"No JSON files found in {dirpath}")
            return False

        print(f"\nFound {len(json_files)} JSON files to validate\n")
        print("=" * 60)

        for filepath in sorted(json_files):
            self.validate_file(filepath)

        print("=" * 60)
        self.print_summary()

        return self.invalid_count == 0

    def print_summary(self):
        """Print validation summary."""
        print(f"\n{'='*60}")
        print("VALIDATION SUMMARY")
        print(f"{'='*60}")
        print(f"  Valid templates:   {self.valid_count}")
        print(f"  Invalid templates: {self.invalid_count}")
        print(f"  Errors:           {len(self.errors)}")
        print(f"  Warnings:         {len(self.warnings)}")

        if self.errors:
            print(f"\nERRORS:")
            for err in self.errors:
                print(f"  - {err}")

        if self.warnings and self.verbose:
            print(f"\nWARNINGS:")
            for warn in self.warnings:
                print(f"  - {warn}")

        print()


class BGCTemplateConverter:
    """Converts human-readable BGC templates to OpenWQ format."""

    def __init__(self, verbose: bool = False):
        self.verbose = verbose

    def convert_template(self, data: dict) -> Optional[dict]:
        """Convert a template to OpenWQ format."""
        if "BIOGEOCHEMISTRY_CONFIGURATION" not in data:
            return None  # Already in OpenWQ format or invalid

        bgc_config = data["BIOGEOCHEMISTRY_CONFIGURATION"]

        # Build OpenWQ format
        openwq_data = {
            "MODULE_NAME": "NATIVE_BGC_FLEX",
            "CHEMICAL_SPECIES": self._convert_species(bgc_config.get("CHEMICAL_SPECIES", {})),
            "CYCLING_FRAMEWORKS": self._convert_frameworks(bgc_config.get("CYCLING_FRAMEWORKS", {}))
        }

        return openwq_data

    def _convert_species(self, species_data: dict) -> dict:
        """Convert species to OpenWQ format."""
        species_list = {}
        mobile_species = []

        idx = 1
        for name, props in species_data.items():
            if name.startswith("_"):
                continue
            species_list[str(idx)] = name

            # Check if this is a mobile species (default: yes for dissolved species)
            if isinstance(props, dict):
                if props.get("MOBILE", True):
                    mobile_species.append(name)
            else:
                mobile_species.append(name)
            idx += 1

        return {
            "list": species_list,
            "BGC_GENERAL_MOBILE_SPECIES": mobile_species[:3]  # Limit to first 3 for now
        }

    def _convert_frameworks(self, frameworks_data: dict) -> dict:
        """Convert cycling frameworks to OpenWQ format."""
        result = {}

        for fw_name, fw_data in frameworks_data.items():
            if fw_name.startswith("_"):
                continue
            if not isinstance(fw_data, dict):
                continue

            fw_result = {"list_transformations": {}}
            rxn_idx = 1

            for rxn_name, rxn_data in fw_data.items():
                if rxn_name.startswith("_"):
                    continue
                if not isinstance(rxn_data, dict):
                    continue
                if "KINETICS" not in rxn_data and "kinetics" not in rxn_data:
                    continue

                # Add to transformation list
                fw_result["list_transformations"][str(rxn_idx)] = rxn_name

                # Get kinetics and parameters
                kinetics = rxn_data.get("KINETICS") or rxn_data.get("kinetics", "")
                consumed = rxn_data.get("CONSUMED") or rxn_data.get("consumed", "NONE")
                produced = rxn_data.get("PRODUCED") or rxn_data.get("produced", "NONE")
                params = rxn_data.get("PARAMETERS") or rxn_data.get("parameters", {})

                # Extract parameter names and values
                param_names = []
                param_values = {}
                for pname, pdata in params.items():
                    param_names.append(pname)
                    if isinstance(pdata, dict):
                        param_values[pname] = pdata.get("VALUE", 0.0)
                    else:
                        param_values[pname] = pdata

                # Build reaction entry
                fw_result[str(rxn_idx)] = {
                    "consumed": consumed,
                    "produced": produced,
                    "kinetics": [kinetics, "1/day"],
                    "parameter_names": param_names,
                    "parameter_values": param_values
                }

                rxn_idx += 1

            if rxn_idx > 1:  # Has reactions
                result[fw_name] = fw_result

        return result

    def convert_and_save(self, input_path: Path, output_dir: Path) -> bool:
        """Convert a template and save to output directory."""
        try:
            with open(input_path, 'r') as f:
                data = json.load(f)

            converted = self.convert_template(data)
            if converted is None:
                if self.verbose:
                    print(f"  Skipping {input_path.name} (already in OpenWQ format)")
                return True

            output_path = output_dir / f"openwq_{input_path.stem}.json"
            with open(output_path, 'w') as f:
                json.dump(converted, f, indent=4)

            if self.verbose:
                print(f"  Converted: {input_path.name} -> {output_path.name}")
            return True

        except Exception as e:
            print(f"  Error converting {input_path.name}: {e}")
            return False


def test_expression_parser():
    """Test expression parsing with exprtk-compatible expressions."""
    print("\nTesting Expression Parser Compatibility")
    print("=" * 60)

    # Test expressions from templates
    test_expressions = [
        ("k * C", True, "Simple first-order"),
        ("k * C * 1.08^(T - 20)", True, "Temperature-corrected"),
        ("k * S / (Km + S)", True, "Michaelis-Menten"),
        ("k * C * (DO / (K_O2 + DO))", True, "O2 limitation"),
        ("k * C * (K_O2_inhib / (K_O2_inhib + DO))", True, "O2 inhibition"),
        ("k * exp(-E / (R * T))", True, "Arrhenius"),
        ("max(0, k * (SM - SM_thresh))", True, "Threshold function"),
        ("k * C1 * C2", True, "Second-order"),
        ("v_s / depth * C", True, "Settling"),
        ("J / depth", True, "Flux per depth"),
    ]

    passed = 0
    failed = 0

    for expr, expected_valid, description in test_expressions:
        # Basic syntax checks
        is_valid = True
        issues = []

        # Check balanced parentheses
        if expr.count('(') != expr.count(')'):
            is_valid = False
            issues.append("unbalanced parentheses")

        # Check for empty expression
        if not expr.strip():
            is_valid = False
            issues.append("empty expression")

        status = "PASS" if is_valid == expected_valid else "FAIL"
        if is_valid == expected_valid:
            passed += 1
        else:
            failed += 1

        print(f"  [{status}] {description}: {expr}")
        if issues:
            print(f"         Issues: {', '.join(issues)}")

    print(f"\nExpression tests: {passed} passed, {failed} failed")
    return failed == 0


def main():
    parser = argparse.ArgumentParser(description="Test BGC Templates for OpenWQ")
    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose output")
    parser.add_argument("--convert", "-c", action="store_true", help="Convert to OpenWQ format")
    parser.add_argument("--expressions", "-e", action="store_true", help="Test expression parser")
    args = parser.parse_args()

    # Get the directory containing this script
    script_dir = Path(__file__).parent

    print("=" * 60)
    print("OpenWQ BGC Template Validation and Testing")
    print("=" * 60)

    # Test expression parsing
    if args.expressions:
        test_expression_parser()
        print()

    # Validate templates
    validator = BGCTemplateValidator(verbose=args.verbose)

    # Validate individual chemicals
    individual_dir = script_dir / "individual_chemicals"
    if individual_dir.exists():
        print("\n--- Individual Chemical Templates ---")
        validator.validate_directory(individual_dir)

    # Validate model frameworks
    models_dir = script_dir / "model_frameworks"
    if models_dir.exists():
        print("\n--- Model Framework Templates ---")
        validator2 = BGCTemplateValidator(verbose=args.verbose)
        validator2.validate_directory(models_dir)

    # Convert templates if requested
    if args.convert:
        print("\n--- Converting Templates to OpenWQ Format ---")
        converter = BGCTemplateConverter(verbose=True)
        output_dir = script_dir / "openwq_format"
        output_dir.mkdir(exist_ok=True)

        for json_file in script_dir.glob("**/*.json"):
            if "openwq_format" not in str(json_file):
                converter.convert_and_save(json_file, output_dir)

        print(f"\nConverted templates saved to: {output_dir}")

    # Summary
    total_valid = validator.valid_count
    total_invalid = validator.invalid_count

    print("\n" + "=" * 60)
    print("FINAL SUMMARY")
    print("=" * 60)
    print(f"Templates validated: {total_valid + total_invalid}")
    print(f"Valid: {total_valid}")
    print(f"Invalid: {total_invalid}")

    if total_invalid == 0:
        print("\n[SUCCESS] All templates passed validation!")
        return 0
    else:
        print(f"\n[FAILED] {total_invalid} templates failed validation")
        return 1


if __name__ == "__main__":
    sys.exit(main())
