# OpenWQ Scenario Analysis - Usage Guide

What-if runs on a CALIBRATED openWQ model: climate change, land-use change,
agricultural practices, urban and point sources, policy targets and in-stream
measures, every option applied through the model's input files.

## Folder structure

```
4_Scenarios/
  scenario_config_template.py   # copy next to your case, set config_run_path, run
  scenario_lib/
    scenarios.py                # option catalogue (GROUPS, CATEGORIES, LEVERS) + runner
    Gen_Scenario_Setup_Report.py   # the interactive scenario setup report
    Gen_Scenario_Results_Report.py # the comparison report (baseline + every scenario)
  USAGE_GUIDE.md
```

`scenario_lib` builds on `../3_Calibration/calibration_lib` (model runner,
parameter handler, reach mapping, report helpers); it puts that folder on the
Python path itself.

## Workflow

1. Model config template -> config report -> run the model -> output report.
2. Calibration template -> calibration report -> calibration -> results report.
3. In the calibration results report, save the calibrated model config
   `<model_template>_config_run.py` and run it once (it writes the calibrated
   inputs for the full period plus `calibrated_setup.json`).
4. Copy `scenario_config_template.py` next to your case, set `config_run_path`
   to that file and run it: the scenario setup report opens.
5. In the report: `+ New scenario`, tick the climate options, then the
   management & policy options for all units and, where needed, for groups of
   units picked on the map. Options with a published relation ask for
   measurable inputs (buffer width, slope, rates, treatment class, ...) and
   compute the N and P change themselves.
6. Save the scenario script `<template>_scenarios_run.py` and run it:

```
python <template>_scenarios_run.py            # baseline + every scenario, then the report
python <template>_scenarios_run.py --report   # rebuild the comparison report only
python <template>_scenarios_run.py --keep     # keep existing scenario run folders
```

On HPC set the container runtime to apptainer in the report before saving the
script; the script accepts `OWQ_SCENARIO_LIB_DIR` and `OWQ_CALIB_LIB_DIR` to
point at copies of the two library folders.

Full documentation: `wikipage/source/4_5_Scenarios.rst` (readthedocs,
"Scenario Analysis").
