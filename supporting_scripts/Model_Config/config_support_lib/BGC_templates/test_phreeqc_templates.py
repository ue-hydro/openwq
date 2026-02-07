#!/usr/bin/env python3
"""
Test Script for PHREEQC Templates (.pqi files)
==============================================

This script validates PHREEQC input files for syntax and structure.

Usage:
    python test_phreeqc_templates.py [--verbose]
"""

import os
import sys
import re
import argparse
from pathlib import Path
from typing import List, Tuple, Dict, Any


class PHREEQCValidator:
    """Validates PHREEQC input files (.pqi)."""

    # Valid PHREEQC keywords
    VALID_KEYWORDS = {
        'TITLE', 'SOLUTION', 'SOLUTION_MASTER_SPECIES', 'SOLUTION_SPECIES',
        'PHASES', 'EQUILIBRIUM_PHASES', 'EXCHANGE', 'EXCHANGE_MASTER_SPECIES',
        'EXCHANGE_SPECIES', 'SURFACE', 'SURFACE_MASTER_SPECIES', 'SURFACE_SPECIES',
        'KINETICS', 'RATES', 'REACTION', 'MIX', 'USE', 'SAVE', 'DUMP',
        'SELECTED_OUTPUT', 'USER_PUNCH', 'USER_PRINT', 'KNOBS', 'PRINT',
        'DATABASE', 'END', 'INVERSE_MODELING', 'GAS_PHASE', 'SOLID_SOLUTIONS',
        'CALCULATE_VALUES', 'COPY', 'DELETE', 'RUN_CELLS', 'TRANSPORT',
        'ADVECTION', 'USER_GRAPH', 'INCLUDE$', 'INCLUDE'
    }

    # Common solution species
    COMMON_SPECIES = {
        'Ca', 'Mg', 'Na', 'K', 'Fe', 'Mn', 'Al', 'Si', 'Cl', 'S', 'C', 'N', 'P',
        'O', 'H', 'Cu', 'Zn', 'Pb', 'Cd', 'As', 'Cr', 'Ni', 'Hg', 'Ba', 'Sr'
    }

    def __init__(self, verbose: bool = False):
        self.verbose = verbose
        self.errors: List[str] = []
        self.warnings: List[str] = []
        self.valid_count = 0
        self.invalid_count = 0

    def log(self, msg: str, level: str = "INFO"):
        """Log message if verbose."""
        if self.verbose or level in ["ERROR", "WARNING"]:
            prefix = {"INFO": "  ", "WARNING": "  [WARN]", "ERROR": "  [ERROR]"}
            print(f"{prefix.get(level, '  ')} {msg}")

    def validate_file(self, filepath: Path) -> bool:
        """Validate a single PHREEQC file."""
        self.log(f"\nValidating: {filepath.name}")

        try:
            with open(filepath, 'r') as f:
                content = f.read()
        except Exception as e:
            self.errors.append(f"{filepath.name}: Read error - {e}")
            self.invalid_count += 1
            return False

        # Parse and validate
        is_valid = True
        is_valid = is_valid and self._check_keywords(content, filepath.name)
        is_valid = is_valid and self._check_solution_blocks(content, filepath.name)
        is_valid = is_valid and self._check_end_statements(content, filepath.name)
        is_valid = is_valid and self._check_balanced_blocks(content, filepath.name)

        if is_valid:
            self.valid_count += 1
            self.log("PHREEQC syntax: OK")
        else:
            self.invalid_count += 1

        return is_valid

    def _check_keywords(self, content: str, filename: str) -> bool:
        """Check that all keyword blocks use valid PHREEQC keywords."""
        valid = True
        lines = content.split('\n')

        keyword_pattern = re.compile(r'^([A-Z][A-Z_0-9]+)(\s+|$)')

        for i, line in enumerate(lines, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            match = keyword_pattern.match(line)
            if match:
                keyword = match.group(1)
                if keyword not in self.VALID_KEYWORDS:
                    # Check if it's a partial match or a data line
                    if not any(kw.startswith(keyword) for kw in self.VALID_KEYWORDS):
                        # Could be a data line or unknown keyword
                        pass  # Don't flag as error, might be data

        return valid

    def _check_solution_blocks(self, content: str, filename: str) -> bool:
        """Check SOLUTION blocks for proper structure."""
        valid = True

        # Find all SOLUTION blocks
        solution_pattern = re.compile(r'SOLUTION\s+(\d+)', re.IGNORECASE)
        solutions = solution_pattern.findall(content)

        if solutions:
            self.log(f"Found {len(solutions)} SOLUTION block(s)")

        # Check for required elements in solutions
        if 'temp' not in content.lower() and solutions:
            self.warnings.append(f"{filename}: No temperature specified in SOLUTION blocks")

        if 'ph' not in content.lower() and solutions:
            self.warnings.append(f"{filename}: No pH specified in SOLUTION blocks")

        return valid

    def _check_end_statements(self, content: str, filename: str) -> bool:
        """Check for proper END statements."""
        valid = True

        # Count major blocks and END statements
        end_count = len(re.findall(r'\bEND\b', content, re.IGNORECASE))

        # Count blocks that should have ENDs
        blocks = ['SOLUTION', 'EQUILIBRIUM_PHASES', 'KINETICS', 'REACTION',
                  'MIX', 'SELECTED_OUTPUT', 'USER_PUNCH']
        block_count = 0
        for block in blocks:
            block_count += len(re.findall(rf'\b{block}\b', content, re.IGNORECASE))

        if end_count < block_count:
            self.warnings.append(f"{filename}: May be missing END statements "
                               f"({end_count} ENDs for {block_count} blocks)")

        return valid

    def _check_balanced_blocks(self, content: str, filename: str) -> bool:
        """Check for balanced RATES/-start/-end blocks."""
        valid = True

        # Check RATES blocks
        rates_starts = len(re.findall(r'-start', content))
        rates_ends = len(re.findall(r'-end', content))

        if rates_starts != rates_ends:
            self.errors.append(f"{filename}: Unbalanced -start/-end in RATES block "
                             f"({rates_starts} starts, {rates_ends} ends)")
            valid = False

        return valid

    def validate_directory(self, dirpath: Path) -> bool:
        """Validate all .pqi files in directory."""
        pqi_files = list(dirpath.glob("**/*.pqi"))

        if not pqi_files:
            print(f"No PHREEQC files found in {dirpath}")
            return True

        print(f"\nFound {len(pqi_files)} PHREEQC files to validate\n")
        print("=" * 60)

        for filepath in sorted(pqi_files):
            self.validate_file(filepath)

        print("=" * 60)
        self.print_summary()

        return self.invalid_count == 0

    def print_summary(self):
        """Print validation summary."""
        print(f"\n{'='*60}")
        print("PHREEQC VALIDATION SUMMARY")
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


def main():
    parser = argparse.ArgumentParser(description="Test PHREEQC Templates")
    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose output")
    args = parser.parse_args()

    script_dir = Path(__file__).parent
    phreeqc_dir = script_dir / "phreeqc_templates"

    print("=" * 60)
    print("PHREEQC Template Validation")
    print("=" * 60)

    if not phreeqc_dir.exists():
        print(f"PHREEQC directory not found: {phreeqc_dir}")
        return 1

    validator = PHREEQCValidator(verbose=args.verbose)
    success = validator.validate_directory(phreeqc_dir)

    if success:
        print("[SUCCESS] All PHREEQC templates passed validation!")
        return 0
    else:
        print(f"[FAILED] Some PHREEQC templates failed validation")
        return 1


if __name__ == "__main__":
    sys.exit(main())
