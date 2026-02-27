## OpenWQ Supporting Scripts

This folder contains supporting scripts for OpenWQ model configuration, calibration, post-processing, and visualization.

## Quick Start

```bash
# 1. Create virtual environment
./setup_venv.sh

# 2. Activate environment
source .venv/bin/activate

# 3. Run scripts
cd Model_Config
python model_config_template.py
```

## Installation

### Automated Setup (Recommended)

```bash
cd supporting_scripts
./setup_venv.sh
```

This creates a virtual environment in `.venv/` with all dependencies.

### Manual Setup

```bash
# Create virtual environment
python3 -m venv .venv

# Activate
source .venv/bin/activate  # Linux/macOS
# .venv\Scripts\activate   # Windows

# Install dependencies
pip install -r requirements.txt
```

### Conda Environment

```bash
conda create -n openwq python=3.10
conda activate openwq
pip install -r requirements.txt
```

## Contents

#### Calibration Framework
Location: `Calibration/`

A comprehensive Python calibration framework for OpenWQ that includes:
- **DDS Optimizer**: Dynamically Dimensioned Search for efficient parameter optimization
- **Sensitivity Analysis**: Morris screening and Sobol variance-based methods
- **HPC Support**: SLURM/PBS job templates for high-performance computing
- **Parameter Defaults**: Scientifically-grounded parameter ranges from peer-reviewed literature
- **Post-Calibration Reporting**: Auto-generated interactive HTML reports with Plotly.js charts, Leaflet.js basin maps, and per-basin multi-variant comparisons

Key documentation:
- [Calibration Usage Guide](Calibration/USAGE_GUIDE.md) - Step-by-step calibration guide
- [Sensitivity Analysis Guide](Calibration/SENSITIVITY_ANALYSIS_GUIDE.md) - Comprehensive parameter sensitivity documentation
- [Calibration Config Template](Calibration/calibration_config_template.py) - User-facing configuration

#### Model Configuration Generator
Location: `Model_Config/`

Python-based configuration generator that creates all OpenWQ input JSON files from a single template script (`model_config_template.py`). Edit one file, run it, and get:

- **All JSON configuration files** -- master, config, BGC, transport, source/sink, sorption, sediment, lateral exchange
- **Interactive HTML report** (`openwq_report.html`) -- a self-contained configuration review tool featuring:
  - Project summary with KPI grid (species, compartments, solver, timestep)
  - Interactive basin map (Leaflet.js) with river network, basin polygons, and GRQA monitoring stations
  - GRQA observation data -- nearby stations from the Global River Water Quality Archive with stoichiometric conversion tables
  - Full model configuration tables with collapsible module parameter details
  - Source/sink load statistics
  - Time series charts (Plotly.js) and spatial statistics (if HDF5 outputs exist)
  - **Ready-to-copy Docker commands** to start the container and run the model
  - **Ready-to-copy Python code snippets** to read HDF5 outputs, generate WebGL 3D maps, and plot time series
  - Dark/light theme toggle, error-resilient rendering, copy-to-clipboard buttons

**Workflow:** Edit the template -> run it -> open the HTML report to review your configuration -> copy the Docker commands to run the model -> copy the Python snippets to visualize results.

#### Output Reading
Location: `Read_Outputs/`

Scripts for reading and processing OpenWQ HDF5 output files:
- HDF5 support library
- Visualization utilities

### Quick Reference: Parameter Sensitivity

| Category | Most Sensitive Parameters | Typical Impact on Results |
|----------|--------------------------|---------------------------|
| **Source/Sink** | Export coefficients, CSV scale factors | 🔴 Very High (sets mass balance) |
| **BGC Rates** | k_nitrification, k_denitrification | 🔴 High (N cycling) |
| **Sorption** | Kfr_PO4, qmax_PO4, bulk_density | 🔴 High (P retention) |
| **Sediment** | erosion_index, slope_exponent | 🟡 Medium-High (erosion) |
| **Transport** | dispersion_x | 🟡 Medium (mixing) |
| **Lateral Exchange** | K_val coefficients | 🟢 Low-Medium |

For detailed parameter-by-parameter analysis, see the [Sensitivity Analysis Guide](Calibration/SENSITIVITY_ANALYSIS_GUIDE.md).

### Usage

```python
# Example: Run sensitivity analysis
from Calibration.calibration_lib.sensitivity.morris_screening import MorrisScreening
from Calibration.calibration_lib.parameter_defaults import BGC_PARAMETERS

# List all available parameters
from Calibration.calibration_lib.parameter_defaults import list_all_parameters
print(list_all_parameters())
```

## Dependencies

All dependencies are listed in `requirements.txt`:

| Package | Purpose |
|---------|---------|
| numpy | Numerical operations |
| pandas | Data manipulation |
| scipy | Optimization, statistics |
| h5py | HDF5 file reading |
| netCDF4 | NetCDF file reading |
| geopandas | Geospatial operations |
| shapely | Geometry operations |
| fiona | Shapefile I/O |
| matplotlib | Plotting |
| requests | GRQA data download |
| SALib | Sensitivity analysis (optional) |

## Troubleshooting

### CRS/Projection Errors

If you encounter PROJ database errors when reading shapefiles:

```bash
pip install --upgrade pyproj
```

### SALib Not Found

SALib is optional but recommended for sensitivity analysis:

```bash
pip install SALib
```

### HDF5 Issues

If h5py fails to install:

```bash
# macOS
brew install hdf5
pip install h5py

# Ubuntu/Debian
sudo apt-get install libhdf5-dev
pip install h5py
```
