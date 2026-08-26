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
  6. Click "Save the script" to get <template>_run.py
  7. Run:  python <template>_run.py

COMMAND-LINE OPTIONS (this template — generates the report):
  python my_calibration.py                     # Generate interactive setup report
  python my_calibration.py --show-parameters   # Show auto-extracted parameters
  python my_calibration.py --dry-run           # Validate config without report

COMMAND-LINE OPTIONS (the generated <template>_run.py — runs the calibration):
  python <template>_run.py                     # Run the calibration
  python <template>_run.py --resume            # Resume from the last checkpoint
  python <template>_run.py --clean             # Delete stale eval_* folders first
                                               #   (no prompt — ideal for HPC/batch jobs)
  python <template>_run.py --keep              # Keep stale eval_* folders (no prompt)
  # Without --clean/--keep, the "Stale evaluation folders" choice from the
  # report applies (prompt when interactive, keep otherwise).

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

# ─────────────────────────────────────────────────────────────────────────────
#  MODEL CHAIN  (single model = list of one; chained models = list of many)
# ─────────────────────────────────────────────────────────────────────────────
#  Ordered list of model config templates to run for EVERY calibration
#  evaluation. THE ORDER IS SEQUENTIAL AND SIGNIFICANT — it is the exact order
#  in which the models are executed, from most UPSTREAM to most DOWNSTREAM:
#
#    • Index 0 runs FIRST  (the most upstream model, e.g. SUMMA-openWQ).
#    • Each later model runs AFTER the ones before it, and typically INGESTS
#      the previous model's openWQ flux output via its EWF. That coupling is
#      declared INSIDE the downstream template itself
#      (ewf_h5_source_folder = "../<upstream>/openwq_out/HDF5"), NOT here — so
#      this list only fixes the RUN ORDER (the calibration rewires the EWF to
#      each evaluation's per-model folders automatically).
#    • A model must therefore appear AFTER every model it depends on; getting
#      the order wrong means a downstream model reads a stale/missing flux file.
#
#  VALIDATION: the LAST entry is the model scored against observations (its
#  openWQ output drives the calibration objective). Build the chain up to the
#  model you want to validate — do not add models beyond it. (There is no
#  per-model "validate" flag: validating anything but the terminal model after
#  running the whole chain would be meaningless.)
#
#  Backward compatible: a single-model calibration is just a one-element list
#  (that one model is both the chain and the validation target). A bare string
#  is also accepted and treated as a one-element chain.
model_chain = [
    # 1st — UPSTREAM: SUMMA-openWQ (produces the flux h5 that mizuRoute ingests)
    "/Users/diogocosta/Documents/openwq_code/diogo_test/mizuRoute-OpenWQ/route/build/openwq/openwq/supporting_scripts/1_Model_Config/model_config_template_test.py",
    # 2nd — DOWNSTREAM + VALIDATED: mizuRoute-openWQ (scored against observations)
    # "/path/to/model_config_template_MIZUROUTE_....py",
]

# Directory where calibration evaluations are stored
calibration_work_dir = "/Users/diogocosta/Documents/openwq_code/calibration_workflow_test"

# Optional: path to an HPC settings JSON that pre-fills the HPC / Apptainer
# fields in the interactive report (Execution tab) — HPC username, host,
# working dir, .sif path, and the SLURM directives.  Defaults to the shipped
# "hpc_settings.json" next to this template.  Copy that file anywhere, fill in
# your HPC details once, and point this at your copy to reuse it across runs.
hpc_settings_json = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                 "hpc_settings.json")

# ── Normalize the chain + resolve the validated model ────────────────────────
# A bare string is treated as a one-element chain. The VALIDATED model (scored
# against observations) is the LAST entry; the setup report and single-model
# parameter extraction below use it, while the calibration runtime receives the
# full ordered chain. Backward compatible: a legacy `model_config_path = "..."`
# still works if `model_chain` is left undefined.
try:
    model_chain
except NameError:
    model_chain = [model_config_path]  # legacy single-model configs
if isinstance(model_chain, str):
    model_chain = [model_chain]
model_chain = [str(p) for p in model_chain if p and str(p).strip()]
if not model_chain:
    raise ValueError("model_chain is empty — set at least one model config template path")
# The last model in the chain is the validation target (obs/objective + report).
model_config_path = model_chain[-1]


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

    # ── Chained calibration: extract + tag parameters from EVERY model ──
    # obs/container/species already come from the VALIDATED (last) model above.
    # Here we also pull the UPSTREAM models' parameters, tag each with its chain
    # position (model_index) + prefix its name (m{i}:) for uniqueness, and merge
    # them into module_parameters so the one report lists all models' params.
    # For a single-model chain this branch is skipped (byte-for-byte unchanged).
    _chain_models = None
    _model_chain_paths = None
    if len(model_chain) > 1:
        _model_chain_paths = list(model_chain)
        _combined = {}
        _chain_models = []
        for _mi, _mpath in enumerate(model_chain):
            _is_last = (_mi == len(model_chain) - 1)
            if _is_last:
                _mcfg, _mp = model_cfg, module_parameters   # already extracted above
            else:
                _mcfg = config_integration.load_model_config(_mpath)
                _bgc = config_integration.get_bgc_template_path(_mcfg)
                _auto = (extract_parameters.extract_calibration_parameters(_bgc)
                         if _bgc else [])
                _ssl = config_integration.get_ss_species_with_loads(_mcfg)
                _mp = extract_parameters.extract_all_module_parameters(
                    _mcfg, bgc_params=_auto, ss_load_species=_ssl)
            for _gk, _plist in _mp.items():
                for _p in _plist:
                    _p['model_index'] = _mi
                    _p['model_label'] = _mcfg.get('hostmodel', f'model{_mi}')
                    _p['name'] = f"m{_mi}:{_p['name']}"
                _combined.setdefault(_gk, []).extend(_plist)
            _chain_models.append({
                'label': f"Model {_mi + 1} — {_mcfg.get('hostmodel', '?')}",
                'hostmodel': _mcfg.get('hostmodel', ''),
                'is_target': _is_last})
        module_parameters = _combined
        print(f"      Chain: {len(model_chain)} models; "
              f"{sum(len(v) for v in module_parameters.values())} params total "
              f"(validating '{model_cfg.get('hostmodel','?')}')")

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
        calibration_template_path=os.path.abspath(__file__),
        calibration_work_dir=os.path.abspath(calibration_work_dir),
        module_parameters=module_parameters,
        module_selections=module_selections,
        species_obs_availability=species_obs_availability,
        ss_species_with_loads=ss_load_species,
        hpc_settings_path=hpc_settings_json,
        model_chain_paths=_model_chain_paths,
        chain_models=_chain_models,
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
