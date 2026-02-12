# OpenWQ Assistant Instructions

This file provides context for Claude when working with the OpenWQ water quality modeling framework.

## Project Overview

**OpenWQ** is a flexible, open-source water quality modeling module designed to couple with hydrological models (mizuRoute, SUMMA, CRHM, etc.). It supports multiple biogeochemistry engines, sediment transport, sorption processes, and comprehensive calibration workflows.

## Repository Structure

```
openwq/
├── src/                          # C++ source code
│   ├── OpenWQ_couplercalls.cpp   # Host model coupling interface
│   ├── OpenWQ_solver.cpp         # Numerical solvers
│   ├── OpenWQ_chem_bgc.cpp       # Biogeochemistry engine
│   └── ...
├── supporting_scripts/           # Python utilities
│   ├── Model_Config/             # Configuration generators
│   │   ├── config_support_lib/   # Configuration libraries
│   │   │   ├── BGC_templates/    # BGC JSON templates
│   │   │   │   ├── NATIVE_BGC_FLEX/  # Flexible BGC framework
│   │   │   │   └── PHREEQC/          # Geochemistry templates
│   │   │   ├── sorption_module/      # Sorption isotherms
│   │   │   └── sediment_transport/   # Sediment modules
│   │   └── source_sink_module/       # Load generation
│   ├── Calibration/              # Calibration framework
│   │   ├── calibration_config_template.py  # User template (copy this)
│   │   └── calibration_lib/      # Core calibration modules
│   │       ├── calibration_driver.py
│   │       ├── objective_functions.py
│   │       ├── dds_optimizer.py
│   │       ├── sensitivity_analysis.py
│   │       └── observation_data/
│   └── Read_Outputs/             # HDF5 output readers
├── in_a_nutshell/                # HTML slide presentations
├── wikipage/source/              # Sphinx RST documentation
├── containers/                   # Docker/Apptainer configs
└── test_cases/                   # Example configurations
```

## Key Configuration Files

### Master Configuration (JSON)
- `openwq_config.json` - Main OpenWQ configuration
- `bgc_config.json` - Biogeochemistry reactions
- `sorption_config.json` - Sorption isotherms
- `sediment_config.json` - Sediment transport
- `source_sink_config.json` - External loads

### Supporting Scripts Entry Points
- `Model_Config/generate_all_configs.py` - Generate all config files
- `Calibration/calibration_config_template.py` - Calibration setup (copy and edit)
- `Read_Outputs/read_openwq_outputs.py` - Read HDF5 results

## Common Tasks

### 1. Setting Up BGC Reactions (NATIVE_BGC_FLEX)
```python
# Use the flexible BGC framework for custom reactions
# Templates in: supporting_scripts/Model_Config/config_support_lib/BGC_templates/NATIVE_BGC_FLEX/

# Key reaction types:
# - TRANSFORMATIONS: A -> B (e.g., nitrification)
# - CONSUMPTION: A consumed by process
# - PRODUCTION: A produced by process
```

### 2. Running Calibration
```bash
# 1. Copy template
cp calibration_config_template.py my_calibration.py

# 2. Edit paths and parameters in my_calibration.py

# 3. Run
python my_calibration.py                    # Full calibration
python my_calibration.py --sensitivity-only # SA only
python my_calibration.py --dry-run          # Validate config
python my_calibration.py --resume           # Resume from checkpoint
```

### 3. Calibration Temporal Resolution
```python
# Aggregate obs/model to common scale before computing metrics
temporal_resolution = "monthly"  # native, daily, weekly, monthly, yearly
aggregation_method = "mean"      # mean, sum, median, min, max
```

### 4. Reading Model Outputs
```python
from Read_Outputs.hdf5_support_lib import read_openwq_hdf5
data = read_openwq_hdf5.load_results("output.h5")
```

## Objective Functions (Calibration)
- **KGE** (recommended): Kling-Gupta Efficiency
- **NSE**: Nash-Sutcliffe Efficiency
- **RMSE**: Root Mean Square Error
- **PBIAS**: Percent Bias

## Parameter Types for Calibration
- `bgc_json`: BGC rate constants
- `phreeqc_pqi`: PHREEQC geochemistry
- `sorption_json`: Sorption isotherms
- `sediment_json`: Sediment transport
- `ss_csv_scale`: Source/sink scaling
- `transport_json`: Dispersion coefficients

## Docker/Apptainer Execution
```bash
# Docker (local)
docker-compose up

# Apptainer (HPC)
apptainer exec openwq.sif ./mizuroute_openwq -g 1 1
```

## Code Conventions
- C++ source uses OpenMP for parallelization
- Python scripts use pandas/numpy for data processing
- Configuration files are JSON (except PHREEQC .pqi files)
- Output files are HDF5 format

## Documentation
- Wiki: `wikipage/source/*.rst` (Sphinx RST format)
- Presentations: `in_a_nutshell/*.html` (self-contained slides)
- API docs: Generated from source docstrings

## When Helping Users

1. **Configuration questions**: Check `supporting_scripts/Model_Config/` for templates
2. **Calibration questions**: Reference `Calibration/calibration_lib/` modules
3. **BGC setup**: Look at `BGC_templates/NATIVE_BGC_FLEX/` examples
4. **Output analysis**: Use `Read_Outputs/` utilities
5. **Compilation issues**: Check `CMakeLists.txt` and `containers/`

## Important Notes

- Always use the template pattern: copy template → edit → run
- Calibration supports checkpoint/resume for long runs
- GRQA database integration for automatic observation data extraction
- HPC support via SLURM/PBS job arrays
- Source/sink parameters are typically the most sensitive - start calibration there
