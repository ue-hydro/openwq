Automatic Model Configuration
========================================

OpenWQ provides a Python-based configuration generator that automates the creation of all input JSON files from a single template script. This approach ensures consistency across configuration files and simplifies the setup of complex simulations.

Location
~~~~~~~~

The configuration scripts are located at::

    openwq/supporting_scripts/Model_Config/

Key files:

* ``model_config_template.py`` - Main template script (edit this file)
* ``config_support_lib/Gen_Input_Driver.py`` - Core generation engine
* ``config_support_lib/Gen_Master_file.py`` - Master JSON generator
* ``config_support_lib/Gen_TSmodule_file.py`` - Transport Sediments module JSON generator
* ``config_support_lib/Gen_LEmodule_file.py`` - Lateral Exchange module JSON generator
* ``config_support_lib/examples_BGC_frameworks/`` - Example biogeochemistry configurations
* ``config_support_lib/examples_PHREEQC/`` - Example PHREEQC configurations
* ``config_support_lib/examples_SS_in_cvs/`` - Example source/sink CSV files


Quick start
~~~~~~~~~~~

1. Copy the template script::

    cp model_config_template.py my_project_config.py

2. Edit the configuration parameters in ``my_project_config.py``

3. Run the script::

    python my_project_config.py

4. The script generates all OpenWQ input files in the specified output directory


Configuration parameters
~~~~~~~~~~~~~~~~~~~~~~~~~

**General information:**

.. code-block:: python

    project_name = "My Project"
    geographical_location = "Study Area"
    authors = "Your Name"
    date = "January, 2026"
    comment = "Project description"

    hostmodel = "mizuroute"  # "mizuroute" or "summa"
    dir2save_input_files = "/path/to/output/"

    # Docker support (see Docker section below)
    running_on_docker = False

**Computational settings:**

.. code-block:: python

    solver = "FORWARD_EULER"           # "FORWARD_EULER" (Forward Euler) or "SUNDIALS" (CVode)
    run_mode_debug = True   # Enable debug output
    use_num_threads = 4     # Number of threads (or "all")

**Initial conditions:**

.. code-block:: python

    ic_all_value = 2.0        # Initial concentration value
    ic_all_units = "mg/l"     # Units for initial conditions


Module configuration
~~~~~~~~~~~~~~~~~~~~~

**Biogeochemistry module:**

.. code-block:: python

    # Option 1: User-defined kinetic expressions
    bgc_module_name = "NATIVE_BGC_FLEX"
    path2selected_NATIVE_BGC_FLEX_framework = "config_support_lib/examples_BGC_frameworks/openWQ_Ncycling_example.json"

    # Option 2: PHREEQC geochemical engine
    bgc_module_name = "PHREEQC"
    phreeqc_input_filepath = "config_support_lib/examples_PHREEQC/phreeqc_river.pqi"
    phreeqc_database_filepath = "config_support_lib/examples_PHREEQC/phreeqc.dat"
    phreeqc_mobile_species = ["H", "O", "Charge", "Ca"]
    phreeqc_component_h2o = True
    phreeqc_temperature_mapping = {
        "RIVER_NETWORK_REACHES": "air_temperature"
    }

**Transport module:**

.. code-block:: python

    # Advection only
    td_module_name = "NATIVE_TD_ADV"

    # Advection-dispersion
    td_module_name = "OPENWQ_NATIVE_TD_ADVDISP"
    td_module_dispersion_xyz = [0.3, 0.3, 0.3]  # Dispersion coefficients [Dx, Dy, Dz] in m2/s
    td_module_characteristic_length_m = 100.0     # Distance between cell centers [m]

    # No transport
    td_module_name = "NONE"

**Lateral exchange module:**

.. code-block:: python

    le_module_name = "NATIVE_LE_BOUNDMIX"
    le_module_config = [
        {
            "direction": "z",
            "upper_compartment": "RUNOFF",
            "lower_compartment": "ILAYERVOLFRACWAT_SOIL",
            "K_val": 0.000000001
        }
    ]

**Sediment transport module:**

Three options are available: ``HYPE_MMF`` (Morgan-Morgan-Finney), ``HYPE_HBVSED`` (HBV-sed), or ``NONE``.

.. code-block:: python

    ts_module_name = "HYPE_HBVSED"  # or "HYPE_MMF", "NONE"
    ts_sediment_compartment = "RIVER_NETWORK_REACHES"
    ts_transport_compartment = "RIVER_NETWORK_REACHES"

    # Common configuration
    ts_direction = "z"                                        # Exchange direction: "x", "y", or "z"
    ts_erosion_inhibit_compartment = "ILAYERVOLFRACWAT_SNOW"  # Compartment inhibiting erosion (e.g., snow)
    ts_data_format = "JSON"                                   # "JSON" or "ASCII"

HYPE_MMF settings (Morgan-Morgan-Finney):

.. code-block:: python

    # Default parameter values (applied to all cells unless overridden)
    # IMPORTANT: CROPCOVER + GROUNDCOVER must sum to 1.0
    ts_mmf_defaults = {
        "COHESION": 7.5,              # Soil cohesion (kPa)
        "ERODIBILITY": 2.0,           # Soil erodibility (g/J)
        "SREROEXP": 1.2,              # Surface runoff erosion exponent (-)
        "CROPCOVER": 0.5,             # Crop cover fraction (0-1)
        "GROUNDCOVER": 0.5,           # Ground cover fraction (0-1)
        "SLOPE": 0.4,                 # Basin slope
        "TRANSPORT_FACTOR_1": 0.2,    # Transport capacity factor 1
        "TRANSPORT_FACTOR_2": 0.5     # Transport capacity factor 2
    }

    # Spatially-varying parameters (cell-specific overrides)
    # Row "0" is the header; rows "1","2",... are data entries
    # Use "ALL" for ix/iy/iz to apply to all cells in that dimension
    ts_mmf_parameters = {
        "COHESION": {
            "0": ["IX", "IY", "IZ", "VALUE"],
            "1": ["ALL", 1, 1, 7.5]
        },
        "ERODIBILITY": {
            "0": ["IX", "IY", "IZ", "VALUE"],
            "1": ["ALL", 1, 1, 2.0]
        },
        # ... repeat for all parameters
    }

HYPE_HBVSED settings (HBV-sed):

.. code-block:: python

    # Default parameter values
    ts_hbvsed_defaults = {
        "SLOPE": 0.4,                                   # Basin slope
        "EROSION_INDEX": 0.4,                           # Erosion index
        "SOIL_EROSION_FACTOR_LAND_DEPENDENCE": 0.4,     # Land use factor (lusepar)
        "SOIL_EROSION_FACTOR_SOIL_DEPENDENCE": 0.4,     # Soil factor (soilpar)
        "SLOPE_EROSION_FACTOR_EXPONENT": 1.5,           # Slope exponent (slopepar)
        "PRECIP_EROSION_FACTOR_EXPONENT": 1.5,          # Precipitation exponent (precexppar)
        "PARAM_SCALING_EROSION_INDEX": 0.5              # Erosion index scaling (eroindexpar)
    }

    # Spatially-varying parameters (same format as HYPE_MMF)
    ts_hbvsed_parameters = {
        "SLOPE": {
            "0": ["IX", "IY", "IZ", "VALUE"],
            "1": ["ALL", 1, 1, 0.4]
        },
        # ... repeat for all parameters
    }

    # Monthly erosion factor (12 values, Jan-Dec)
    # Used as: erodmonth = 1.0 + MONTHLY_EROSION_FACTOR[month]
    ts_hbvsed_monthly_erosion_factor = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                        0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

.. note::

    The ``PARAMETERS`` section provides spatially-varying overrides. If a parameter is not specified
    in ``PARAMETERS``, the value from ``PARAMETER_DEFAULTS`` is used for all cells.
    The pseudo day-of-year (MMF) and monthly erosion factor (HBVSED) are dynamically
    computed from the simulation time at runtime.

**Sorption isotherm module:**

.. code-block:: python

    si_module_name = "NONE"  # or "FREUNDLICH", "LANGMUIR"
    si_sediment_compartment = "SUMMA_RUNOFF"

    # Soil/medium properties (used by both Freundlich and Langmuir)
    si_bulk_density_kg_m3 = 1500.0   # Bulk density [kg/m3]
    si_layer_thickness_m = 1.0       # Representative layer thickness [m]

    # Per-species isotherm parameters
    # For FREUNDLICH:
    si_species_params = {
        "NH4-N": {"Kfr": 1.2, "Nfr": 0.8, "Kadsdes_1_per_s": 0.001}
    }
    # For LANGMUIR:
    si_species_params = {
        "NH4-N": {"qmax_mg_per_kg": 200.0, "KL_L_per_mg": 0.05, "Kadsdes_1_per_s": 0.002}
    }
    # Set to None when si_module_name = "NONE"
    si_species_params = None

See :doc:`Sorption Isotherms <4_1_3Modules>` for detailed parameter descriptions and JSON format.


Source/sink configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~

**Method 1: Load from CSV files:**

.. code-block:: python

    ss_method = "load_from_csv"
    ss_method_csv_config = [
        {
            "Chemical_name": "NO3-N",
            "Compartment_name": "RIVER_NETWORK_REACHES",
            "Type": "source",
            "Units": "kg",
            "Filepath": "/path/to/source_data.csv",
            "Delimiter": ","
        }
    ]

The header row is auto-detected by scanning for the ``YYYY`` column.
The CSV must contain a row with columns: ``YYYY, MM, DD, HH, MIN, SEC, ix, iy, iz, load, load_type, time_units``.

**Method 2: Load from Copernicus land use data:**

.. code-block:: python

    ss_method = "using_copernicus_lulc_with_static_coeff"
    ss_method_copernicus_nc_lc_dir = '/path/to/ESACCI-LC/'
    ss_method_copernicus_period = [1993, 2020]
    ss_method_copernicus_default_loads_bool = True


Output configuration
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    output_format = "HDF5"              # "HDF5" or "CSV"
    chemical_species = ["NO3-N", "NH4-N", "N_ORG_fresh", "N_ORG_stable"]  # Species names to export
    units = "MG/L"                      # Output units
    no_water_conc_flag = -9999          # No-data flag
    export_sediment = False             # Export sediment results
    timestep = [1, "hour"]              # Output timestep
    compartments_and_cells = {
        "RIVER_NETWORK_REACHES": {
            "1": ["all", "all", "all"]
        }
    }


Generated files
~~~~~~~~~~~~~~~~

The configuration generator creates the following files in your output directory:

* ``openWQ_master.json`` -- Master configuration file (references all module files)
* ``openwq_in/openWQ_config.json`` -- Initial conditions and compartment mapping
* ``openwq_in/openWQ_MODULE_<BGC>.json`` -- Biogeochemistry module configuration (NATIVE_BGC_FLEX or PHREEQC)
* ``openwq_in/openWQ_MODULE_<TD>.json`` -- Transport dissolved module configuration
* ``openwq_in/openWQ_MODULE_<LE>.json`` -- Lateral exchange module configuration
* ``openwq_in/openWQ_MODULE_<TS>.json`` -- Transport sediments module configuration (HYPE_MMF or HYPE_HBVSED)
* ``openwq_in/openWQ_MODULE_<SI>.json`` -- Sorption isotherm module configuration
* ``openwq_in/openWQ_SS_<method>.json`` -- Source/sink configuration (if applicable)
* ``openwq_in/openWQ_EWF_fixed_value.json`` -- External water flux configuration (if applicable)


Docker support
~~~~~~~~~~~~~~~

When running OpenWQ inside a Docker container, file paths in the generated JSON files must reference the container filesystem rather than the host filesystem. The ``running_on_docker`` flag automates this path correction.

.. code-block:: python

    running_on_docker = True   # Set to True when model runs inside Docker

When enabled, the generator:

1. Reads ``docker-compose.yml`` (located at the openwq root) to determine the volume mount mapping
2. Resolves the relative host path from the compose file location
3. Replaces host path prefixes with container path prefixes in all JSON content

For example, with the default volume mount ``../../../../../../:/code:Z``:

* **Host path**: ``/Users/username/Documents/project/test_case/``
* **Container path**: ``/code/project/test_case/``

.. note::

    The path correction only applies to paths **embedded inside generated JSON files**
    (which OpenWQ reads at runtime inside the container). The files themselves are still
    written to the host filesystem so the script can create them on disk.
    Paths used for file I/O during generation (e.g., BGC framework files, PHREEQC database,
    shapefiles) are **not** corrected since the script runs on the host.


Tips
~~~~

1. **PHREEQC species discovery**: Run once with ``phreeqc_chemical_species_names = []`` to see the component list in the log, then update accordingly

2. **Debug mode**: Set ``run_mode_debug = True`` during development to get detailed solver output

3. **Compartment names**: Enable debug mode to see the exact compartment names exported by the host model

4. **Units consistency**: Ensure initial condition units match the expected units in the biogeochemistry module

5. **Sediment parameters**: When ``PARAMETERS`` is left empty (``{}``), warnings are shown but the model uses ``PARAMETER_DEFAULTS`` for all cells. Populate ``PARAMETERS`` with at least one entry per parameter to suppress warnings

6. **Docker paths**: When ``running_on_docker = True``, verify the generated master JSON contains container paths (e.g., ``/code/...``) rather than host paths

7. **Spatially-varying overrides**: Use ``"ALL"`` in the IX/IY/IZ fields to apply a value to all cells in that dimension. Add multiple rows (``"2"``, ``"3"``, ...) for cell-specific overrides
