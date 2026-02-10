# OpenWQ Calibration Framework - Step-by-Step Usage Guide

This guide provides detailed instructions for setting up and running model calibration using the OpenWQ Calibration Framework.

## Table of Contents

1. [Prerequisites](#1-prerequisites)
2. [Installation Verification](#2-installation-verification)
3. [Directory Structure Setup](#3-directory-structure-setup)
4. [Preparing Observation Data](#4-preparing-observation-data)
5. [Creating Your Calibration Configuration](#5-creating-your-calibration-configuration)
6. [Running Calibration (Docker - Local)](#6-running-calibration-docker---local)
7. [Running Calibration (Apptainer - HPC)](#7-running-calibration-apptainer---hpc)
8. [Monitoring Progress](#8-monitoring-progress)
9. [Analyzing Results](#9-analyzing-results)
10. [Troubleshooting](#10-troubleshooting)

---

## 1. Prerequisites

### Required Software

```bash
# Python 3.8+ with required packages
pip install numpy pandas h5py matplotlib

# Optional but recommended for sensitivity analysis
pip install SALib
```

### Required Files

1. **Working OpenWQ model setup** - A test case that runs successfully
2. **Observation data** - CSV file with measured concentrations
3. **Either:**
   - **Docker** (for local runs) - Container `docker_openwq` running
   - **Apptainer** (for HPC runs) - `.sif` image file available

### Verify Docker Setup (Local)

```bash
# Check Docker container is running
docker ps | grep openwq

# If not running, start it:
cd /path/to/openwq/containers
docker-compose up -d
```

### Verify Apptainer Setup (HPC)

```bash
# Check Apptainer/Singularity is available
apptainer --version
# or
singularity --version

# Verify your .sif file exists
ls -la /path/to/openwq.sif
```

---

## 2. Installation Verification

Run the test suite to verify the framework is working:

```bash
cd openwq/supporting_scripts/Calibration
python test_calibration_framework.py
```

**Expected output:**
```
============================================================
 Test Summary
============================================================
  [PASS] Module Imports
  [PASS] DDS Optimizer
  [PASS] Morris Screening
  [PASS] Checkpoint Manager
  [PASS] Parameter Handler
  [PASS] Model Runner
  [PASS] Results Analyzer

  Total: 7/7 tests passed

  All tests PASSED!
  The calibration framework is ready to use.
```

If any tests fail, check the error messages and verify your Python environment.

---

## 3. Directory Structure Setup

Create a calibration workspace with this structure:

```
calibration_workspace/
├── my_calibration.py          # Your calibration config (from template)
├── observations.csv           # Your observation data
├── evaluations/               # Created automatically - model runs
├── results/                   # Created automatically - output files
└── checkpoints/               # Created automatically - restart files
```

**Commands:**
```bash
# Create workspace
mkdir -p /path/to/calibration_workspace
cd /path/to/calibration_workspace

# Copy configuration template
cp /path/to/openwq/supporting_scripts/Calibration/calibration_config_template.py my_calibration.py
```

---

## 4. Preparing Observation Data

You have two options for preparing observation data:

### Option A: Manual CSV File

Create a CSV file with your observation data. Required columns:

| Column | Description | Example |
|--------|-------------|---------|
| `datetime` | Timestamp (ISO format) | `2010-06-01 00:00:00` |
| `reach_id` | River reach ID (must match model) | `1200014181` |
| `species` | Chemical species name | `NO3-N` |
| `value` | Measured concentration | `2.50` |
| `units` | Concentration units | `mg/l` |

**Example `observations.csv`:**
```csv
datetime,reach_id,species,value,units,source,uncertainty,quality_flag
2010-06-01 00:00:00,1200014181,NO3-N,2.50,mg/l,USGS_station_A,0.25,GOOD
2010-06-01 00:00:00,1200014181,NH4-N,0.15,mg/l,USGS_station_A,0.02,GOOD
2010-06-15 00:00:00,1200014181,NO3-N,2.80,mg/l,USGS_station_A,0.28,GOOD
2010-07-01 00:00:00,1200014181,NO3-N,3.20,mg/l,USGS_station_A,0.32,GOOD
```

**Important:**
- `reach_id` must match the reach IDs in your model output HDF5 files
- `species` must match exactly the species names in your BGC configuration
- Use consistent units throughout (typically `mg/l` or `MG/L`)

### Option B: GRQA Database Extraction (Recommended for Large-Scale Calibration)

The framework can automatically download and process data from the **Global River Water Quality Archive (GRQA)** database. This is ideal for large-scale calibrations.

**Step 4B.1: Enable GRQA in your configuration**

```python
# Enable automatic GRQA extraction
use_grqa = True
```

**Step 4B.2: Configure GRQA settings**

```python
grqa_config = {
    # River network shapefile for spatial matching
    "river_network_shapefile": "/path/to/river_network.shp",

    # Column in shapefile containing reach IDs
    "reach_id_column": "seg_id",

    # Bounding box [min_lon, min_lat, max_lon, max_lat]
    "bounding_box": [-125, 24, -66, 50],  # CONUS example

    # Time period
    "start_date": "2000-01-01",
    "end_date": "2020-12-31",

    # Maximum distance to match stations to reaches
    "max_station_distance_m": 500,

    # Map GRQA parameters to model species
    "species_mapping": {
        "NO3": "NO3-N",
        "NH4": "NH4-N",
        "TN": "TN",
        "PO4": "PO4-P",
        "TP": "TP",
    },
}
```

**Step 4B.3: Run GRQA extraction**

```bash
# Extract only (without calibration)
python my_calibration.py --extract-grqa-only

# Or run calibration (extraction happens automatically)
python my_calibration.py
```

**GRQA Species Mapping Reference:**

| GRQA Parameter | Description | Typical Model Species |
|----------------|-------------|----------------------|
| NO3 | Nitrate | NO3-N |
| NH4 | Ammonium | NH4-N |
| TN | Total Nitrogen | TN |
| NO2 | Nitrite | NO2-N |
| PO4 | Orthophosphate | PO4-P |
| TP | Total Phosphorus | TP |
| DO | Dissolved Oxygen | DO |
| BOD | Biochemical Oxygen Demand | BOD |
| DOC | Dissolved Organic Carbon | DOC |
| TSS | Total Suspended Solids | TSS |

**GRQA Dependencies:**

```bash
pip install geopandas requests shapely
```

---

## 5. Creating Your Calibration Configuration

Edit `my_calibration.py` with your specific settings:

### Step 5.1: Set Paths

```python
# =============================================================================
# PATH CONFIGURATION
# =============================================================================
base_model_config_dir = "/path/to/openwq/supporting_scripts/Model_Config"
test_case_dir = "/path/to/your/test_case"  # Contains openwq_in, mizuroute_in, etc.
calibration_work_dir = "/path/to/calibration_workspace"
observation_data_path = "/path/to/calibration_workspace/observations.csv"
```

### Step 5.2: Configure Container Runtime

**For Docker (local):**
```python
container_runtime = "docker"
docker_container_name = "docker_openwq"
docker_compose_path = "/path/to/openwq/containers/docker-compose.yml"
```

**For Apptainer (HPC):**
```python
container_runtime = "apptainer"
apptainer_sif_path = "/path/to/openwq.sif"
apptainer_bind_path = "/scratch/user/openwq_code:/code"
```

### Step 5.3: Define Parameters to Calibrate

```python
calibration_parameters = [
    # BGC rate constant (example: nitrification rate)
    {
        "name": "k_nitrification",
        "file_type": "bgc_json",
        "path": ["CYCLING_FRAMEWORKS", "N_cycle", "3", "parameter_values", "k"],
        "initial": 0.03,
        "bounds": (0.001, 0.5),
        "transform": "log"  # Use log-space for rate constants
    },

    # Source/sink scaling factor
    {
        "name": "fertilizer_scale",
        "file_type": "ss_csv_scale",
        "path": {"species": "all"},
        "initial": 1.0,
        "bounds": (0.1, 5.0),
        "transform": "linear"
    },

    # Add more parameters as needed...
]
```

**Finding the correct `path` for BGC parameters:**
1. Open your `openWQ_MODULE_NATIVE_BGC_FLEX.json`
2. Navigate to the parameter you want to calibrate
3. The path is the sequence of keys: `["CYCLING_FRAMEWORKS", "cycle_name", "reaction_index", "parameter_values", "param_name"]`

### Step 5.4: Configure Optimization Settings

```python
# Algorithm settings
algorithm = "DDS"           # Recommended for expensive models
max_evaluations = 300       # Start with 200-500

# Objective function
objective_function = "KGE"  # Kling-Gupta Efficiency (recommended)

# Target species and locations
calibration_targets = {
    "species": ["NO3-N", "NH4-N"],      # Species to calibrate against
    "reach_ids": "all",                  # Or list specific: [1200014181, 200014181]
    "compartments": ["RIVER_NETWORK_REACHES"]
}

# Random seed for reproducibility
random_seed = 42
```

### Step 5.5: (Optional) Configure Sensitivity Analysis

For calibrations with many parameters (10+), run sensitivity analysis first:

```python
run_sensitivity_first = True
sensitivity_method = "morris"
sensitivity_morris_trajectories = 10
sensitivity_morris_levels = 4
sensitivity_threshold = 0.1  # Parameters with mu_star < 10% of max are non-influential
```

---

## 6. Running Calibration (Docker - Local)

### Step 6.1: Verify Docker is Running

```bash
docker ps | grep openwq
# Should show: docker_openwq   Up ...
```

### Step 6.2: Run Calibration

```bash
cd /path/to/calibration_workspace

# Normal run
python my_calibration.py

# Dry run (validate configuration without running model)
python my_calibration.py --dry-run

# Sensitivity analysis only
python my_calibration.py --sensitivity-only
```

### Step 6.3: Resume Interrupted Calibration

If calibration is interrupted, resume from the last checkpoint:

```bash
python my_calibration.py --resume
```

---

## 7. Running Calibration (Apptainer - HPC)

### Step 7.1: Update Configuration for HPC

```python
# Container runtime
container_runtime = "apptainer"
apptainer_sif_path = "/scratch/user/openwq/openwq.sif"
apptainer_bind_path = "/scratch/user/openwq_code:/code"

# HPC settings
hpc_enabled = True
hpc_scheduler = "slurm"       # or "pbs"
hpc_partition = "standard"    # Your cluster partition
hpc_walltime = "24:00:00"     # Per evaluation
hpc_nodes = 1
hpc_tasks_per_node = 4
hpc_memory = "8G"
hpc_max_concurrent_jobs = 50  # Limit concurrent jobs
```

### Step 7.2: Transfer Files to HPC

```bash
# Copy calibration workspace to HPC
rsync -avz /local/calibration_workspace/ user@hpc:/scratch/user/calibration/

# Copy observations
scp observations.csv user@hpc:/scratch/user/calibration/
```

### Step 7.3: Submit Job

**Option A: Interactive mode (for testing)**
```bash
# SSH to HPC
ssh user@hpc

# Load modules
module load python/3.9
module load apptainer  # or singularity

# Run calibration
cd /scratch/user/calibration
python my_calibration.py
```

**Option B: Batch submission**

Create a wrapper script `run_calibration.sh`:
```bash
#!/bin/bash
#SBATCH --job-name=openwq_calib
#SBATCH --partition=standard
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH --output=calibration_%j.out
#SBATCH --error=calibration_%j.err

module load python/3.9
module load apptainer

cd /scratch/user/calibration
python my_calibration.py
```

Submit:
```bash
sbatch run_calibration.sh
```

### Step 7.4: Monitor HPC Jobs

```bash
# Check job status
squeue -u $USER

# View SLURM logs
tail -f calibration_*.out

# Check evaluation progress
ls -la evaluations/ | wc -l
```

---

## 8. Monitoring Progress

### Console Output

During calibration, you'll see progress updates:

```
[INFO] Starting calibration with DDS algorithm
[INFO] Max evaluations: 300
[INFO] Number of parameters: 5

[INFO] Eval 1/300: objective = 0.4523
[INFO] Eval 2/300: objective = 0.3891 (new best!)
[INFO] Eval 3/300: objective = 0.4102
...
[INFO] Eval 50/300: objective = 0.1234 (new best!)
```

### Checkpoint Files

The framework automatically saves checkpoints after each evaluation:

```
calibration_workspace/
├── calibration_state.json     # Current state (can view with text editor)
├── calibration_history.pkl    # Full history (Python pickle)
└── rng_state.pkl              # Random state (for exact reproducibility)
```

View current best parameters:
```bash
cat calibration_state.json | python -m json.tool
```

---

## 9. Analyzing Results

### Output Files

After calibration completes, results are in the `results/` directory:

```
results/
├── best_parameters.json       # Optimal parameter values
├── calibration_history.json   # All evaluations
├── convergence.png            # Convergence plot
├── parameter_evolution.png    # Parameter values over iterations
├── sensitivity_results.json   # Sensitivity analysis (if run)
└── calibration_report.txt     # Summary report
```

### View Best Parameters

```bash
cat results/best_parameters.json
```

```json
{
  "k_nitrification": 0.0234,
  "fertilizer_scale": 1.45,
  "best_objective": 0.1234
}
```

### Generate Additional Plots

```python
from calibration_lib.postprocessing import ResultsAnalyzer

analyzer = ResultsAnalyzer(results_dir="results/")

# Convergence plot
analyzer.plot_convergence(output_path="convergence.png", show=True)

# Parameter evolution
analyzer.plot_parameter_evolution(output_path="param_evolution.png", show=True)

# Full report
analyzer.generate_report(output_dir="results/")
```

### Apply Best Parameters to Model

Copy the best parameters to your model configuration:

1. Open `results/best_parameters.json`
2. Manually update your BGC/sorption/sediment JSON files with the optimal values
3. Or use the parameter handler to apply them programmatically

---

## 10. Troubleshooting

### Common Issues

#### "Container not found" error

**Docker:**
```bash
# Start the container
cd /path/to/openwq/containers
docker-compose up -d

# Verify it's running
docker ps
```

**Apptainer:**
```bash
# Verify SIF file exists
ls -la /path/to/openwq.sif

# Test it runs
apptainer exec /path/to/openwq.sif echo "Container works!"
```

#### "No observations loaded" error

- Check `observation_data_path` is correct
- Verify CSV format matches expected columns
- Ensure `species` names match exactly (case-sensitive)
- Verify `reach_id` values exist in model output

#### Model runs but objective is NaN or Inf

- Check model output HDF5 files are being created
- Verify observation datetime range overlaps with model output
- Check units match between observations and model

#### HPC jobs fail immediately

```bash
# Check SLURM logs
cat slurm_logs/array_*.err

# Common issues:
# - Module not loaded (add module load commands)
# - Path not found (check bind paths)
# - Memory exceeded (increase hpc_memory)
```

#### Calibration stuck at same objective

- Increase `max_evaluations`
- Widen parameter `bounds`
- Try different `initial` values
- Increase perturbation factor `r` in DDS (default 0.2)

### Getting Help

1. Check the test script: `python test_calibration_framework.py`
2. Enable verbose logging in your config: `logging_level = "DEBUG"`
3. Review the RST documentation in `wikipage/source/4_2_2_Calibration.rst`

---

## Quick Reference

### Command Line Options

```bash
python my_calibration.py                          # Normal calibration run
python my_calibration.py --resume                 # Resume from checkpoint
python my_calibration.py --dry-run                # Validate config only
python my_calibration.py --sensitivity-only       # Run sensitivity analysis only
python my_calibration.py --extract-grqa-only      # Extract GRQA data only
python my_calibration.py --generate-copernicus-only  # Generate Copernicus observations only
```

### Recommended Workflow

1. **Start small**: 50-100 evaluations with 2-3 parameters
2. **Verify setup**: Check outputs are created correctly
3. **Run sensitivity**: Identify influential parameters
4. **Full calibration**: 300-500 evaluations on influential parameters
5. **Validate**: Run model with optimized parameters, compare to observations

### Typical Runtime Estimates

| Model Runtime | Evaluations | Total Time |
|--------------|-------------|------------|
| 5 min | 300 | ~25 hours |
| 15 min | 300 | ~75 hours |
| 30 min | 300 | ~150 hours |

For long calibrations, use HPC with job arrays to run evaluations in parallel.
