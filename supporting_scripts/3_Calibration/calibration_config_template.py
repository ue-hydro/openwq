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
OpenWQ Calibration Configuration Template
==========================================

PREREQUISITE: model_config_template.py must be set up and generating correct
config files FIRST. All model settings (paths, modules, species, observation
data, container settings, etc.) are read from there automatically.

QUICK START:
  1. cp calibration_config_template.py my_calibration.py
  2. Edit the TWO paths below
  3. Run:  python my_calibration.py
  4. An interactive HTML report opens in your browser
  5. Configure settings, parameters, species in the report
  6. Click "Download Script" to get my_calibration_run.py
  7. Run:  python my_calibration_run.py

COMMAND-LINE OPTIONS:
  python my_calibration.py                  # Generate interactive setup report
  python my_calibration.py --show-parameters  # Show auto-extracted parameters
  python my_calibration.py --dry-run          # Validate config without report

Full documentation: https://openwq.readthedocs.io
"""

import sys
import os

# Add calibration library to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


# ╔═══════════════════════════════════════════════════════════════════════╗
# ║                                                                       ║
# ║   USER CONFIGURATION — Only 2 paths to set                           ║
# ║                                                                       ║
# ╚═══════════════════════════════════════════════════════════════════════╝

# Absolute path to your model_config_template.py
model_config_path = "/Users/diogocosta/Documents/openwq_code/diogo_test/mizuRoute-OpenWQ/route/build/openwq/openwq/supporting_scripts/1_Model_Config/model_config_template_test.py"

# Directory where calibration evaluations are stored
calibration_work_dir = "/Users/diogocosta/Documents/openwq_code/calibration_workflow_test"


# ╔═══════════════════════════════════════════════════════════════════════╗
# ║                                                                       ║
# ║   END OF USER CONFIGURATION                                          ║
# ║                                                                       ║
# ╚═══════════════════════════════════════════════════════════════════════╝
# ── Execution logic below — do not modify ────────────────────────────

import argparse
import webbrowser
import traceback

from calibration_lib import config_integration
from calibration_lib import extract_parameters
from calibration_lib import Gen_Calibration_Setup_Report


def _main():
    """Main execution logic."""

    # ── Parse CLI arguments ──
    parser = argparse.ArgumentParser(
        description="OpenWQ Calibration — Interactive Setup"
    )
    parser.add_argument("--dry-run", action="store_true",
                        help="Validate configuration without generating report")
    parser.add_argument("--show-parameters", action="store_true",
                        help="Show auto-extracted parameters and exit")
    args = parser.parse_args()

    # ── Step 1: Load model configuration ──
    print("\n" + "=" * 60)
    print("OPENWQ CALIBRATION — INTERACTIVE SETUP")
    print("=" * 60)

    print(f"\n[1/4] Loading model config: {model_config_path}")
    try:
        model_cfg = config_integration.load_model_config(model_config_path)
        print(f"      Loaded {len(model_cfg)} configuration variables")
    except Exception as e:
        print(f"      ERROR: {e}")
        traceback.print_exc()
        sys.exit(1)

    # ── Step 2: Auto-extract BGC parameters ──
    print(f"\n[2/4] Extracting calibration parameters...")

    bgc_template_path = config_integration.get_bgc_template_path(model_cfg)

    auto_params = []
    if bgc_template_path:
        print(f"      BGC template: {os.path.basename(bgc_template_path)}")
        try:
            auto_params = extract_parameters.extract_calibration_parameters(
                bgc_template_path
            )
            print(f"      Auto-extracted {len(auto_params)} parameters "
                  f"with RANGE fields")
        except Exception as e:
            print(f"      WARNING: Could not extract parameters: {e}")
    else:
        print("      No BGC template found (NATIVE_BGC_FLEX not selected)")

    # ── Step 2b: Extract parameters for ALL active modules ──
    print(f"\n[2b/5] Extracting parameters for all active modules...")

    module_selections = config_integration.get_module_selections(model_cfg)
    ss_load_species = config_integration.get_ss_species_with_loads(model_cfg)
    module_parameters = extract_parameters.extract_all_module_parameters(
        model_cfg, bgc_params=auto_params, ss_load_species=ss_load_species
    )

    total_module_params = sum(len(v) for v in module_parameters.values())
    for group_key, params in module_parameters.items():
        if params:
            print(f"      {group_key}: {len(params)} parameters")
    print(f"      Total across all modules: {total_module_params}")

    # ── --show-parameters: display and exit ──
    if args.show_parameters:
        print("\n" + "=" * 60)
        print("AUTO-EXTRACTED PARAMETERS (ALL MODULES)")
        print("=" * 60)
        for group_key, params in module_parameters.items():
            if params:
                print(f"\n--- {group_key.upper()} ---")
                extract_parameters.print_parameter_table(params)
        sys.exit(0)

    # ── Step 3: Extract observation & container config ──
    print(f"\n[3/5] Extracting observation and container settings...")

    obs_config = config_integration.get_observation_config(model_cfg)
    container_config = config_integration.get_container_config(model_cfg)
    species_obs_availability = (
        config_integration.get_species_observation_availability(model_cfg)
    )

    # If chemical_species was "all", replace it with the resolved list
    # so the report generator gets actual species names
    if species_obs_availability:
        model_cfg["chemical_species"] = list(species_obs_availability.keys())

    print(f"      Observation source: {obs_config.get('source', 'skip')}")
    obs_count = sum(1 for v in species_obs_availability.values() if v["has_obs"])
    print(f"      Species with observations: {obs_count}/{len(species_obs_availability)}")
    if ss_load_species:
        print(f"      SS species with loads: {', '.join(sorted(ss_load_species))}")
    print(f"      Host model: {container_config.get('hostmodel', 'unknown')}")

    # ── Validate ──
    if args.dry_run:
        print(f"\n[DRY RUN] Configuration validated.")
        print(f"  Model config: OK ({len(model_cfg)} variables)")
        print(f"  Total calibration parameters: {total_module_params}")
        for gk, params in module_parameters.items():
            if params:
                print(f"    {gk}: {len(params)}")
        print(f"  Observation source: {obs_config.get('source')}")
        print(f"\n  All checks passed. Ready to generate interactive report.")
        sys.exit(0)

    # ── Step 4: Generate interactive setup report ──
    print(f"\n[4/5] Generating interactive setup report...")

    os.makedirs(calibration_work_dir, exist_ok=True)

    report_path = Gen_Calibration_Setup_Report.generate_interactive_setup(
        output_dir=calibration_work_dir,
        model_config=model_cfg,
        auto_extracted_parameters=auto_params,
        observation_config=obs_config,
        container_config=container_config,
        model_config_path=os.path.abspath(model_config_path),
        calibration_work_dir=os.path.abspath(calibration_work_dir),
        module_parameters=module_parameters,
        module_selections=module_selections,
        species_obs_availability=species_obs_availability,
        ss_species_with_loads=ss_load_species,
    )

    if report_path:
        print(f"      Report: {report_path}")
        print(f"\n{'=' * 60}")
        print("NEXT STEPS:")
        print("  1. Open the report in your browser")
        print("  2. Configure calibration settings, species, parameters")
        print("  3. Click 'Download Script' to get my_calibration_run.py")
        print("  4. Run: python my_calibration_run.py")
        print(f"{'=' * 60}\n")
        webbrowser.open(f"file://{os.path.abspath(report_path)}")
    else:
        print("      ERROR: Report generation failed")
        sys.exit(1)


# Entry point
if __name__ == "__main__":
    _main()
