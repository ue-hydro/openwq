Model Calibration Framework
===========================

OpenWQ provides a comprehensive Python calibration framework for parameter optimization and sensitivity analysis. The framework is designed for computationally expensive models, supporting both local (Docker) and HPC (Apptainer/SLURM) execution environments.

Key features:

* **DDS Optimization**: Dynamically Dimensioned Search algorithm, efficient for expensive models (200-1000 evaluations)
* **Sensitivity Analysis**: Morris screening and Sobol variance-based methods for parameter importance ranking
* **All Parameter Types**: Supports BGC rates, PHREEQC geochemistry, sorption isotherms, sediment transport, source/sink loads, and more
* **HPC Support**: SLURM/PBS job arrays for parallel evaluation on clusters
* **Checkpointing**: Save and resume long-running calibrations

Location
~~~~~~~~

The calibration framework is located at::

    openwq/supporting_scripts/Calibration/

Key files:

* ``calibration_config_template.py`` - User configuration template (copy and edit)
* ``calibration_lib/calibration_driver.py`` - Main entry point
* ``calibration_lib/`` - Core library modules (including observation_data/)


Quick Start
~~~~~~~~~~~

1. Copy the configuration template::

    cd openwq/supporting_scripts/Calibration/
    cp calibration_config_template.py my_calibration.py

2. Edit paths and parameters in ``my_calibration.py``

3. Run calibration::

    python my_calibration.py

4. Resume if interrupted::

    python my_calibration.py --resume

5. Run sensitivity analysis only::

    python my_calibration.py --sensitivity-only

6. Validate configuration (dry run)::

    python my_calibration.py --dry-run


Configuration Overview
~~~~~~~~~~~~~~~~~~~~~~

The configuration file has seven main sections:

1. **Path Configuration** - Model config, observations, working directory
2. **Container Runtime** - Docker (local) or Apptainer (HPC) settings
3. **Calibration Parameters** - Parameters to optimize with bounds
4. **Calibration Settings** - Algorithm, max evaluations, objective function
5. **Sensitivity Analysis** - Morris/Sobol settings
6. **HPC Settings** - SLURM/PBS job configuration
7. **Run Calibration** - Entry point


Supported Parameter Types
~~~~~~~~~~~~~~~~~~~~~~~~~

The framework supports calibration of all major OpenWQ parameter types:

**BGC Parameters (NATIVE_BGC_FLEX)**

.. code-block:: python

    {
        "name": "k_nitrification",
        "file_type": "bgc_json",
        "path": ["CYCLING_FRAMEWORKS", "N_cycle", "3", "parameter_values", "k"],
        "initial": 0.03,
        "bounds": (0.001, 0.5),
        "transform": "log"
    }

**PHREEQC Geochemistry Parameters**

.. code-block:: python

    # Initial solution concentrations
    {
        "name": "initial_NO3",
        "file_type": "phreeqc_pqi",
        "path": {"block": "SOLUTION", "species": "N(5)"},
        "initial": 2.0,
        "bounds": (0.1, 10.0)
    }

    # Equilibrium phase partial pressures
    {
        "name": "pCO2_log",
        "file_type": "phreeqc_pqi",
        "path": {"block": "EQUILIBRIUM_PHASES", "phase": "CO2(g)", "field": "si"},
        "initial": -3.5,
        "bounds": (-4.5, -2.5)
    }

**Sorption Isotherm Parameters**

.. code-block:: python

    # Freundlich isotherm
    {
        "name": "NH4_Kfr",
        "file_type": "sorption_json",
        "path": {"module": "FREUNDLICH", "species": "NH4-N", "param": "Kfr"},
        "initial": 0.5,
        "bounds": (0.01, 10.0),
        "transform": "log"
    }

    # Langmuir isotherm
    {
        "name": "NH4_qmax",
        "file_type": "sorption_json",
        "path": {"module": "LANGMUIR", "species": "NH4-N", "param": "qmax_mg_per_kg"},
        "initial": 100.0,
        "bounds": (10.0, 500.0)
    }

**Sediment Transport Parameters**

.. code-block:: python

    # HYPE_HBVSED module
    {
        "name": "erosion_index",
        "file_type": "sediment_json",
        "path": {"module": "HYPE_HBVSED", "param": "EROSION_INDEX"},
        "initial": 0.4,
        "bounds": (0.1, 1.0)
    }

    # HYPE_MMF module
    {
        "name": "cohesion",
        "file_type": "sediment_json",
        "path": {"module": "HYPE_MMF", "param": "COHESION"},
        "initial": 7.5,
        "bounds": (1.0, 20.0)
    }

**Source/Sink Parameters**

.. code-block:: python

    # CSV load scaling
    {
        "name": "fertilizer_N_scale",
        "file_type": "ss_csv_scale",
        "path": {"species": "all"},
        "initial": 1.0,
        "bounds": (0.1, 5.0)
    }

    # Copernicus LULC export coefficients
    {
        "name": "cropland_TN_coeff",
        "file_type": "ss_copernicus_static",
        "path": {"lulc_class": 10, "species": "TN"},
        "initial": 20.0,
        "bounds": (5.0, 50.0)
    }

    # Climate response parameters
    {
        "name": "precip_scaling_power",
        "file_type": "ss_climate_param",
        "path": "ss_climate_precip_scaling_power",
        "initial": 1.0,
        "bounds": (0.5, 2.0)
    }

**Transport Parameters**

.. code-block:: python

    {
        "name": "dispersion_x",
        "file_type": "transport_json",
        "path": {"param": "DISPERSION_X_M2_S"},
        "initial": 0.3,
        "bounds": (0.01, 10.0),
        "transform": "log"
    }


Parameter Options
~~~~~~~~~~~~~~~~~

Each parameter definition supports these options:

+---------------+-------------------------------------------------------------------+
| Option        | Description                                                       |
+===============+===================================================================+
| ``name``      | Unique identifier (used in results)                               |
+---------------+-------------------------------------------------------------------+
| ``file_type`` | Type of configuration file (see list above)                       |
+---------------+-------------------------------------------------------------------+
| ``path``      | Location of parameter within the file                             |
+---------------+-------------------------------------------------------------------+
| ``initial``   | Starting value for optimization                                   |
+---------------+-------------------------------------------------------------------+
| ``bounds``    | (min, max) tuple defining search space                            |
+---------------+-------------------------------------------------------------------+
| ``transform`` | ``"linear"`` (default) or ``"log"`` for log-space optimization    |
+---------------+-------------------------------------------------------------------+
| ``dtype``     | ``"float"`` (default) or ``"int"`` for integer parameters         |
+---------------+-------------------------------------------------------------------+


Objective Functions
~~~~~~~~~~~~~~~~~~~

The framework supports four objective functions:

* **RMSE** - Root Mean Square Error
* **NSE** - Nash-Sutcliffe Efficiency (1 - NSE for minimization)
* **KGE** - Kling-Gupta Efficiency (1 - KGE for minimization) - *Recommended*
* **PBIAS** - Percent Bias (absolute value for minimization)

Configure in your calibration file:

.. code-block:: python

    objective_function = "KGE"

    # Multi-species weighting
    objective_weights = {
        "NO3-N": 1.0,
        "NH4-N": 0.5,
    }

    # Target species and locations
    calibration_targets = {
        "species": ["NO3-N", "NH4-N"],
        "reach_ids": "all",  # or [1200014181, 200014181]
        "compartments": ["RIVER_NETWORK_REACHES"]
    }


Temporal Resolution
~~~~~~~~~~~~~~~~~~~

The calibration framework supports flexible temporal resolution for objective function calculation and performance evaluation. Both model outputs and observations are aggregated to the specified resolution before computing metrics.

**Available Resolutions:**

* ``native`` - Use original model/observation timestamps (no aggregation)
* ``daily`` - Aggregate to daily averages
* ``weekly`` - Aggregate to weekly averages
* ``monthly`` - Aggregate to monthly averages
* ``yearly`` - Aggregate to yearly averages

**Aggregation Methods:**

* ``mean`` - Average value (default, recommended for concentrations)
* ``sum`` - Total sum (useful for loads/fluxes)
* ``median`` - Median value (robust to outliers)
* ``min`` / ``max`` - Minimum/maximum values

**Configuration:**

.. code-block:: python

    # Temporal resolution for calibration
    temporal_resolution = "monthly"    # Options: native, daily, weekly, monthly, yearly
    aggregation_method = "mean"        # Options: mean, sum, median, min, max

**When to Use Different Resolutions:**

+---------------+---------------------------------------------------------------+
| Resolution    | Recommended Use Case                                          |
+===============+===============================================================+
| ``native``    | High-frequency data, sub-daily dynamics important             |
+---------------+---------------------------------------------------------------+
| ``daily``     | Standard choice when daily patterns matter                    |
+---------------+---------------------------------------------------------------+
| ``weekly``    | Smooth out day-to-day variability                             |
+---------------+---------------------------------------------------------------+
| ``monthly``   | Seasonal patterns, reduce noise from sparse observations      |
+---------------+---------------------------------------------------------------+
| ``yearly``    | Long-term trends, annual budgets                              |
+---------------+---------------------------------------------------------------+

**How Aggregation Works:**

1. Observations are grouped by reach_id, species, and temporal period
2. Model outputs are extracted and matched to observation timestamps
3. Both are aggregated to the specified resolution using the chosen method
4. Objective functions are computed on the aggregated data
5. Performance plots show both native and aggregated comparisons

**Example Use Case:**

If you have sparse monthly grab samples but an hourly model, using ``temporal_resolution = "monthly"`` ensures fair comparison:

.. code-block:: python

    # Monthly calibration for sparse observations
    temporal_resolution = "monthly"
    aggregation_method = "mean"

    # This aggregates both model outputs and observations to monthly means
    # before computing KGE, NSE, RMSE, etc.

**Performance Plots:**

When temporal aggregation is enabled, the calibration framework generates additional diagnostic plots:

* Time series comparison (observed vs simulated at specified resolution)
* Scatter plot with 1:1 line and performance metrics
* Residual analysis (temporal patterns in errors)

These are saved in ``results/performance_plots/`` after calibration completes.


Observation Data Format
~~~~~~~~~~~~~~~~~~~~~~~

Observation data must be in CSV format with these columns:

.. code-block:: text

    datetime,reach_id,species,value,units,source,uncertainty,quality_flag
    2010-06-01 00:00:00,1200014181,NO3-N,2.50,mg/l,USGS_station_A,0.25,GOOD
    2010-06-01 00:00:00,1200014181,NH4-N,0.15,mg/l,USGS_station_A,0.02,GOOD

Required columns:

* ``datetime`` - Observation timestamp (YYYY-MM-DD HH:MM:SS)
* ``reach_id`` - River network reach ID matching model output
* ``species`` - Chemical species name (must match model species)
* ``value`` - Measured concentration value
* ``units`` - Concentration units (e.g., mg/l)

Optional columns:

* ``source`` - Data source identifier
* ``uncertainty`` - Measurement uncertainty
* ``quality_flag`` - Quality flag (GOOD, SUSPECT, BAD)

An example file is provided at ``observation_data/example_observations.csv``.


GRQA Database Integration
~~~~~~~~~~~~~~~~~~~~~~~~~

The calibration framework can automatically extract observation data from the **Global River Water Quality Archive (GRQA)** database (https://zenodo.org/records/15335450). This is recommended for large-scale calibrations where manually preparing observation data would be impractical.

**Enabling GRQA Extraction**

In your calibration configuration file:

.. code-block:: python

    # Enable GRQA data extraction
    use_grqa = True

    # Configure GRQA extraction
    grqa_config = {
        # River network shapefile for spatial matching
        "river_network_shapefile": "/path/to/river_network.shp",

        # Column in shapefile containing reach IDs
        "reach_id_column": "seg_id",

        # Bounding box to filter stations [min_lon, min_lat, max_lon, max_lat]
        "bounding_box": [-125, 24, -66, 50],  # e.g., CONUS

        # Time period for observations
        "start_date": "2000-01-01",
        "end_date": "2020-12-31",

        # Maximum distance (meters) to match stations to river reaches
        "max_station_distance_m": 500,

        # Species mapping: GRQA parameter names -> Model species names
        "species_mapping": {
            "NO3": "NO3-N",      # Nitrate
            "NH4": "NH4-N",      # Ammonium
            "TN": "TN",          # Total Nitrogen
            "PO4": "PO4-P",      # Orthophosphate
            "TP": "TP",          # Total Phosphorus
            "DO": "DO",          # Dissolved Oxygen
        },
    }

**How GRQA Extraction Works**

1. **Download**: Only required GRQA parameter files are downloaded based on your species mapping
2. **Filter**: Stations are filtered by bounding box and time period
3. **Match**: Stations are spatially matched to model river reaches using buffer distance
4. **Convert**: Units are converted and data is formatted for calibration
5. **Output**: Calibration-ready CSV file is generated automatically

**Command-Line Options**

.. code-block:: bash

    # Extract GRQA data only (without running calibration)
    python my_calibration.py --extract-grqa-only

    # Normal run with automatic GRQA extraction
    python my_calibration.py

**Species Mapping**

The species mapping dictionary maps GRQA parameter names to your model species names. Only parameters listed in this mapping will be downloaded and processed. Available GRQA parameters include:

* **Nitrogen**: NO3, NH4, TN, NO2, DON
* **Phosphorus**: PO4, TP, DP
* **Carbon**: DOC, TOC, BOD, COD
* **Other**: DO, TSS, Chl-a, EC, pH, TDS, Temp

**Dependencies**

GRQA extraction requires additional Python packages:

.. code-block:: bash

    pip install geopandas requests shapely

**Output Files**

GRQA extraction creates the following files in your calibration workspace:

* ``grqa_stations.csv`` - Extracted station metadata
* ``grqa_observations.csv`` - Raw observations from GRQA
* ``calibration_observations.csv`` - Formatted for calibration (auto-used)
* ``station_reach_matching.csv`` - Station-to-reach spatial matches


DDS Optimization Algorithm
~~~~~~~~~~~~~~~~~~~~~~~~~~

The Dynamically Dimensioned Search (DDS) algorithm is a global optimization method designed for computationally expensive models. Key features:

* **Efficient convergence**: Typically requires 200-500 evaluations
* **Adaptive search**: Starts with global exploration, progressively focuses on local refinement
* **No gradient needed**: Uses perturbation-based search

Configure DDS settings:

.. code-block:: python

    algorithm = "DDS"
    max_evaluations = 500
    random_seed = 42  # For reproducibility


Sensitivity Analysis
~~~~~~~~~~~~~~~~~~~~

For calibrations with many parameters (10+), run sensitivity analysis first to identify influential parameters.

**Morris Screening** (fast, qualitative)

.. code-block:: python

    run_sensitivity_first = True
    sensitivity_method = "morris"
    sensitivity_morris_trajectories = 10
    sensitivity_morris_levels = 4
    sensitivity_threshold = 0.1  # Relative to max mu_star

**Sobol Analysis** (slower, quantitative)

.. code-block:: python

    sensitivity_method = "sobol"
    sensitivity_sobol_samples = 1024


Recommended HPC Workflow
~~~~~~~~~~~~~~~~~~~~~~~~

For large-scale calibrations with dozens of parameters:

1. **Morris Screening** (~200 evaluations)

   - Identifies non-influential parameters
   - Typically eliminates 50-70% of parameters

2. **Sobol Sensitivity** (~1000 evaluations on reduced set)

   - Quantifies parameter importance and interactions
   - Further reduces to 10-15 most influential parameters

3. **DDS Optimization** (~300-500 evaluations on final set)

   - Calibrates only influential parameters
   - Fixed non-influential parameters at defaults

**Total**: ~1500-2000 evaluations (vs 10,000+ without screening)


HPC Configuration (SLURM)
~~~~~~~~~~~~~~~~~~~~~~~~~

For running on HPC clusters with Apptainer:

.. code-block:: python

    # Container runtime
    container_runtime = "apptainer"
    apptainer_sif_path = "/path/to/openwq.sif"
    apptainer_bind_path = "/scratch/user/openwq_code:/code"

    # HPC settings
    hpc_enabled = True
    hpc_scheduler = "slurm"
    hpc_partition = "standard"
    hpc_walltime = "24:00:00"
    hpc_nodes = 1
    hpc_tasks_per_node = 4
    hpc_memory = "8G"
    hpc_max_concurrent_jobs = 50

SLURM job templates are provided in ``hpc_templates/``:

* ``slurm_array_template.sh`` - For sensitivity analysis (parallel array)
* ``slurm_worker_template.sh`` - For DDS optimization (worker pool)
* ``pbs_array_template.sh`` - For PBS/Torque schedulers


Checkpointing and Resume
~~~~~~~~~~~~~~~~~~~~~~~~

The framework automatically saves checkpoints after each evaluation. To resume an interrupted calibration::

    python my_calibration.py --resume

Checkpoint files are saved in your working directory:

* ``calibration_state.json`` - Optimization state
* ``calibration_history.pkl`` - Full evaluation history
* ``rng_state.pkl`` - Random number generator state


Results and Outputs
~~~~~~~~~~~~~~~~~~~

After calibration completes, results are saved to ``results/``:

* ``best_parameters.json`` - Optimal parameter values
* ``calibration_history.json`` - All evaluated parameter sets and objectives
* ``convergence.png`` - Convergence plot
* ``parameter_evolution.png`` - Parameter values over evaluations
* ``sensitivity_results.json`` - Sensitivity analysis results (if run)
* ``calibration_report.txt`` - Summary report


Complete Configuration Example
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    # =============================================================================
    # PATH CONFIGURATION
    # =============================================================================
    base_model_config_dir = "/path/to/Model_Config"
    calibration_work_dir = "/path/to/calibration_workspace"
    observation_data_path = "/path/to/observations.csv"
    test_case_dir = "/path/to/test_case"

    # =============================================================================
    # CONTAINER RUNTIME
    # =============================================================================
    container_runtime = "docker"  # or "apptainer"
    docker_container_name = "docker_openwq"
    docker_compose_path = "/path/to/docker-compose.yml"
    executable_name = "mizuroute_lakes_cslm_openwq_fast"
    executable_args = "-g 1 1"
    file_manager_path = "/code/.../fileManager.txt"

    # =============================================================================
    # CALIBRATION PARAMETERS
    # =============================================================================
    calibration_parameters = [
        {
            "name": "k_nitrification",
            "file_type": "bgc_json",
            "path": ["CYCLING_FRAMEWORKS", "N_cycle", "3", "parameter_values", "k"],
            "initial": 0.03,
            "bounds": (0.001, 0.5),
            "transform": "log"
        },
        {
            "name": "fertilizer_N_scale",
            "file_type": "ss_csv_scale",
            "path": {"species": "all"},
            "initial": 1.0,
            "bounds": (0.1, 5.0)
        },
    ]

    # =============================================================================
    # CALIBRATION SETTINGS
    # =============================================================================
    algorithm = "DDS"
    max_evaluations = 500
    n_parallel = 1
    objective_function = "KGE"
    objective_weights = {"NO3-N": 1.0, "NH4-N": 0.5}
    calibration_targets = {
        "species": ["NO3-N", "NH4-N"],
        "reach_ids": "all",
        "compartments": ["RIVER_NETWORK_REACHES"]
    }
    random_seed = 42

    # =============================================================================
    # TEMPORAL RESOLUTION
    # =============================================================================
    temporal_resolution = "monthly"  # native, daily, weekly, monthly, yearly
    aggregation_method = "mean"      # mean, sum, median, min, max


Dependencies
~~~~~~~~~~~~

**Required:**

* numpy
* pandas
* h5py
* matplotlib (for plotting)

**Optional:**

* SALib (for Morris/Sobol sensitivity analysis)::

    pip install SALib

**Runtime:**

* Docker (local) or Apptainer/Singularity (HPC)


Troubleshooting
~~~~~~~~~~~~~~~

**"No observations loaded" error:**

Check that your observation file path is correct and the file contains data matching your target species and reach IDs.

**"Container not found" error:**

Ensure Docker is running (local) or the Apptainer SIF path is correct (HPC).

**Calibration converges to bounds:**

Consider expanding parameter bounds or using log-transform for rate constants.

**HPC jobs fail immediately:**

Check SLURM logs in ``slurm_logs/`` for error messages. Verify bind paths and module loads in the job script.

**Checkpoint resume fails:**

Ensure you're running from the same working directory. Check that ``calibration_state.json`` exists and is not corrupted.
