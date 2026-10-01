"""
OpenWQ Scenario Analysis Configuration Template
================================================

What-if runs on a CALIBRATED model: management practices (BMPs), land-use
change, point sources & wastewater, policy targets, climate change, in-stream
measures, legacy stores — every lever applied through the model's input files.

PREREQUISITE: a calibrated model config. When a calibration finishes, its
results report ("Run the calibrated model") saves
``<model_template>_config_run.py`` next to your model-config template. Run that
file once (like any model-config template): it generates the calibrated inputs
for the full model period and its config report gives the snippets to run the
model and to produce its output report. Point ``config_run_path`` below at it.
(Any model-config template works as the base; a calibrated one is the point.)

QUICK START:
  1. cp scenario_config_template.py my_scenarios.py
  2. Edit config_run_path (and, optionally, scenario_work_dir) below
  3. Run:  python my_scenarios.py
  4. The interactive scenario report opens in your browser
  5. Build the scenarios (presets or your own levers), set the run settings
  6. Click "Save the scenario script" to get <template>_scenarios_run.py
  7. Run:  python <template>_scenarios_run.py
     -> one model run per scenario (+ the baseline) and a comparison report
        with the baseline and every scenario on the same graphs

COMMAND-LINE OPTIONS (this template — generates the report):
  python my_scenarios.py                 # Generate the interactive scenario report
  python my_scenarios.py --dry-run       # Validate the base config without a report

COMMAND-LINE OPTIONS (the generated <template>_scenarios_run.py):
  python <template>_scenarios_run.py             # Run baseline + scenarios, write the report
  python <template>_scenarios_run.py --keep      # Keep existing scenario run folders
  python <template>_scenarios_run.py --report    # Rebuild the comparison report only

Full documentation: https://openwq.readthedocs.io  (Scenario Analysis)
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))          # scenario_lib (this folder)
# calibration_lib (model runner, parameter handler, report helpers) lives in ../3_Calibration
sys.path.insert(0, os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                 "..", "3_Calibration")))

# ═══════════════════════════════════════════════════════════════════════════
#  USER SETTINGS — edit these
# ═══════════════════════════════════════════════════════════════════════════

# The CALIBRATED model config saved by the calibration results report
# (<model_template>_config_run.py), already run once so its inputs exist.
config_run_path = "/path/to/your/model_config_template_X_config_run.py"

# Where the scenario runs + reports go. None -> "<dir2save_input_files>_scenarios"
# next to the base run folder of config_run_path.
scenario_work_dir = None

# Optional: HPC / Apptainer settings JSON (same format as the calibration's
# hpc_settings.json). None -> local Docker runs only.
hpc_settings_json = None

# ═══════════════════════════════════════════════════════════════════════════
#  END OF USER SETTINGS — do not edit below this line
# ═══════════════════════════════════════════════════════════════════════════

import argparse
import webbrowser
import traceback


def _locate_libs():
    """This file may be copied anywhere (next to the case, like the other
    templates). If ``scenario_lib`` / ``calibration_lib`` are not importable
    from the paths above, derive the ``supporting_scripts`` folder from the base
    config's own ``config_support_lib`` path
    (``.../supporting_scripts/1_Model_Config/config_support_lib`` ->
    ``.../supporting_scripts/4_Scenarios`` and ``.../3_Calibration``)."""
    try:
        import scenario_lib  # noqa: F401
        import calibration_lib  # noqa: F401
        return
    except ImportError:
        pass
    import re
    try:
        txt = open(config_run_path, encoding="utf-8", errors="replace").read()
    except Exception:
        return
    for m in re.finditer(r"""sys\.path\.insert\(\s*0\s*,\s*['"]([^'"]+)['"]""", txt):
        cand = m.group(1)
        i = cand.replace("\\", "/").find("/1_Model_Config/")
        if i >= 0:
            for sub, pkg in (("4_Scenarios", "scenario_lib"), ("3_Calibration", "calibration_lib")):
                d = os.path.join(cand[:i], sub)
                if os.path.isdir(os.path.join(d, pkg)) and d not in sys.path:
                    sys.path.insert(0, d)
            return


_locate_libs()
from scenario_lib import Gen_Scenario_Setup_Report
from calibration_lib import config_integration


def _main():
    parser = argparse.ArgumentParser(description="OpenWQ Scenario Analysis — Interactive Setup")
    parser.add_argument("--dry-run", action="store_true",
                        help="Validate the base configuration without generating the report")
    args = parser.parse_args()

    print("\n" + "=" * 60)
    print("OPENWQ SCENARIO ANALYSIS — INTERACTIVE SETUP")
    print("=" * 60)

    if not os.path.isfile(config_run_path):
        print(f"\nERROR: config_run_path not found:\n  {config_run_path}\n"
              "Save the calibrated model config from the calibration results report "
              "(section 'Run the calibrated model'), run it once, then point config_run_path at it.")
        sys.exit(1)

    print(f"\n[1/3] Loading the base (calibrated) model config: {config_run_path}")
    try:
        model_cfg = config_integration.load_model_config(config_run_path)
        print(f"      Loaded {len(model_cfg)} configuration variables")
    except Exception as e:
        print(f"      ERROR: {e}")
        traceback.print_exc()
        sys.exit(1)

    base_dir = str(model_cfg.get("dir2save_input_files") or "")
    work_dir = scenario_work_dir or (base_dir.rstrip("/\\") + "_scenarios")
    print(f"      Base run folder: {base_dir}")
    print(f"      Scenario work dir: {work_dir}")
    if not os.path.isdir(os.path.join(base_dir, "openwq_in")):
        print("      WARNING: the base run folder has no openwq_in/ yet — run the config_run "
              "template once so the calibrated inputs exist (the scenarios copy them).")

    if args.dry_run:
        print("\n[DRY RUN] Base configuration loaded. Ready to generate the scenario report.")
        sys.exit(0)

    print("\n[2/3] Generating the interactive scenario report...")
    os.makedirs(work_dir, exist_ok=True)
    report_path = Gen_Scenario_Setup_Report.generate_interactive_scenario_setup(
        output_dir=work_dir,
        model_config=model_cfg,
        config_run_path=os.path.abspath(config_run_path),
        hpc_settings_path=hpc_settings_json,
        template_path=os.path.abspath(__file__),     # names the report/script after this file
    )
    if not report_path:
        print("      ERROR: report generation failed")
        sys.exit(1)
    print(f"      Report: {report_path}")

    print("\n[3/3] Opening the report in your browser...")
    if os.environ.get("BROWSER", "").lower() not in ("true", "none", "0", "false"):
        webbrowser.open("file://" + os.path.abspath(report_path))
    print("\nNext: build the scenarios in the report, save the scenario script (right pane) and run it.")


if __name__ == "__main__":
    _main()
