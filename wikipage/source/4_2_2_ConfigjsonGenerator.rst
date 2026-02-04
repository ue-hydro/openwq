Automatic Model Configuration
========================================

OpenWQ provides a Python-based configuration generator that automates the creation of all input JSON files from a single template script. This approach ensures consistency across configuration files and simplifies the setup of complex simulations.

Location
~~~~~~~~

The configuration scripts are located at::

    openwq/supporting_scripts/Model_Config/

Key files:

* ``template_model_config.py`` - Main template script (edit this file)
* ``config_support_lib/Gen_Input_Driver.py`` - Core generation engine
* ``config_support_lib/examples_BGC_frameworks/`` - Example biogeochemistry configurations
* ``config_support_lib/examples_PHREEQC/`` - Example PHREEQC configurations


Quick start
~~~~~~~~~~~

1. Copy the template script::

    cp template_model_config.py my_project_config.py

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

**Computational settings:**

.. code-block:: python

    solver = "BE"           # "BE" (Backward Euler) or "SUNDIALS" (CVode)
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
    phreeqc_mobile_species = [1, 2, 3, 4]
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
    td_module_dispersion_xyz = [0.3, 0.3, 0.3]  # Dispersion coefficients

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

**Sediment and sorption modules:**

.. code-block:: python

    # Sediment transport
    ts_module_name = "NONE"  # or "HYPE_MMF", "HYPE_HBVSED"
    ts_sediment_compartment = "RIVER_NETWORK_REACHES"

    # Sorption isotherms
    si_module_name = "NONE"  # or "FREUNDLICH", "LANGMUIR"
    si_sediment_compartment = "SUMMA_RUNOFF"


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
            "Delimiter": ",",
            "Number_of_header_rows": 3,
            "Header_key_row": 3
        }
    ]

**Method 2: Load from Copernicus land use data:**

.. code-block:: python

    ss_method = "using_copernicus_lulc"
    ss_method_copernicus_nc_lc_dir = '/path/to/ESACCI-LC/'
    ss_method_copernicus_period = [1993, 2020]
    ss_method_copernicus_default_loads_bool = True


Output configuration
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    output_format = "HDF5"              # "HDF5" or "CSV"
    chemical_species = [1, 2, 3, 4]     # Species indices to export
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

* ``openwq_master.json`` - Master configuration file
* ``openwq_in/openwq_config.json`` - Initial conditions and compartment mapping
* ``openwq_in/openwq_bgc.json`` - Biogeochemistry configuration
* ``openwq_in/openwq_te.json`` - Transport configuration
* ``openwq_in/openwq_ss.json`` - Source/sink configuration (if applicable)
* ``openwq_in/openwq_ewf.json`` - External water flux configuration (if applicable)
* ``openwq_in/openwq_output.json`` - Output configuration


Tips
~~~~

1. **PHREEQC species discovery**: Run once with ``phreeqc_chemical_species_names = []`` to see the component list in the log, then update accordingly

2. **Debug mode**: Set ``run_mode_debug = True`` during development to get detailed solver output

3. **Compartment names**: Enable debug mode to see the exact compartment names exported by the host model

4. **Units consistency**: Ensure initial condition units match the expected units in the biogeochemistry module
