Repository Structure
==================================

This page describes the organization of the OpenWQ GitHub repository to help you navigate the codebase.

.. code-block:: text

    openwq/
    ├── containers/              # Container definitions for Docker and Apptainer
    │   ├── Dockerfile           # Docker image build instructions
    │   ├── docker-compose.yml   # Docker Compose configuration
    │   └── openwq_apptainer.def # Apptainer/Singularity definition file
    │
    ├── src/                     # C++ source code
    │   ├── compute/             # Core computational routines
    │   ├── extwatflux_ss/       # External water flux and source/sink handling
    │   ├── models_CH/           # Chemistry/BGC module implementations
    │   ├── models_LE/           # Lateral exchange module
    │   ├── models_SI/           # Sorption isotherm module
    │   ├── models_TD/           # Transport dissolved module
    │   ├── models_TS/           # Transport sediment module
    │   └── output/              # Output handling routines
    │
    ├── supporting_scripts/      # Python tools for configuration and analysis
    │   ├── Model_Config/        # Configuration file generators
    │   │   ├── model_config_template.py    # Main configuration template
    │   │   └── config_support_lib/         # Generator libraries
    │   ├── Read_Outputs/        # Post-processing and visualization
    │   │   ├── reading_plotting_results_template.py
    │   │   └── hdf5_support_lib/           # HDF5 reader, plotter, mapper
    │   └── Calibration/         # Calibration and sensitivity analysis
    │       ├── calibration_config_template.py  # Calibration configuration template
    │       └── calibration_lib/            # Optimization algorithms and drivers
    │
    ├── openwq_in_a_nutshell/    # Interactive HTML presentations
    │   ├── 1_OpenWQ_Overview.html
    │   ├── 2_Biogeochemistry_Engines.html
    │   ├── 3_Sources_Sinks_Loads.html
    │   ├── 4_Model_Setup_Guide.html
    │   ├── 5_Calibration_Guide.html
    │   └── 6_Supporting_Scripts.html
    │
    ├── json_files_repos/        # Example JSON configuration files
    │
    ├── testing/                 # Test cases and validation
    │
    ├── utils/                   # Utility scripts (compilation helpers)
    │
    ├── wikipage/                # Documentation source (Sphinx/RST)
    │   └── source/              # RST files for ReadTheDocs
    │
    ├── CMakeLists.txt           # Main CMake build configuration
    ├── main.cpp                 # Entry point (standalone mode)
    ├── README.md                # Quick start guide
    └── LICENSE                  # GNU GPL v3 license


Folder Descriptions
~~~~~~~~~~~~~~~~~~~~

**containers/**

Contains all containerization files for running OpenWQ in isolated environments:

- ``Dockerfile``: Builds a Docker image with Ubuntu 24.04, all dependencies pre-installed (PHREEQC-RM, SUNDIALS, Armadillo, ParallelIO)
- ``docker-compose.yml``: Simplifies container management with volume mounts
- ``openwq_apptainer.def``: Definition file for Apptainer/Singularity (used on HPC clusters)

To build with Docker::

    cd containers
    docker-compose build
    docker-compose up -d

To build with Apptainer::

    cd containers
    apptainer build openwq.sif openwq_apptainer.def


**src/**

The C++ source code organized by module:

- ``compute/``: Core solver routines, time stepping, parallel execution
- ``extwatflux_ss/``: Parsing and application of external water fluxes and source/sink terms
- ``models_CH/``: Biogeochemistry implementations (NATIVE_BGC_FLEX, PHREEQC)
- ``models_LE/``: Lateral exchange between compartments
- ``models_SI/``: Sorption isotherm calculations (Freundlich, Langmuir)
- ``models_TD/``: Dissolved transport (advection, dispersion)
- ``models_TS/``: Sediment transport (HYPE-based erosion models)
- ``output/``: HDF5/CSV output generation


**supporting_scripts/**

Python tools organized into three main categories. All scripts require Python 3.8+ with dependencies.

**Environment Setup**::

    cd supporting_scripts
    ./setup_venv.sh           # Create virtual environment
    source .venv/bin/activate # Activate environment

Or with conda::

    conda create -n openwq python=3.10
    conda activate openwq
    pip install -r requirements.txt

**Available Tools:**

1. **Model_Config/**: Generate all JSON configuration files from a single Python template

   - ``model_config_template.py``: Edit this file to configure your model
   - Run it to auto-generate all required JSON files

2. **Read_Outputs/**: Post-processing and visualization

   - ``Read_h5_driver.py``: Load HDF5 output files into Python
   - ``Plot_h5_driver.py``: Create time series plots
   - ``Map_h5_driver.py``: Generate spatial maps and animated GIFs

3. **Calibration/**: Parameter optimization and sensitivity analysis

   - Morris screening for parameter importance
   - Sobol variance-based sensitivity analysis
   - DDS (Dynamically Dimensioned Search) optimization
   - HPC templates for SLURM/PBS clusters


**openwq_in_a_nutshell/**

Interactive HTML presentations that can be viewed in any web browser. Use arrow keys to navigate. Topics include:

1. Overview of OpenWQ architecture
2. Biogeochemistry engines (NATIVE_BGC_FLEX, PHREEQC)
3. Source/sink and load configuration options
4. Model setup step-by-step guide
5. Calibration and sensitivity analysis guide
6. Supporting scripts usage guide


**wikipage/**

Documentation source files in reStructuredText (RST) format. These are built by Sphinx and hosted on ReadTheDocs. To build locally::

    cd wikipage
    pip install -r requirements.txt
    make html


Key Files
~~~~~~~~~~

- ``CMakeLists.txt``: CMake build configuration. Set ``COMPILE_TARGET`` and ``SOLVER_TYPE`` here.
- ``main.cpp``: Standalone entry point (used when OpenWQ runs independently, not coupled to a host model)
- ``README.md``: Quick start guide with compilation and Docker instructions
- ``.readthedocs.yaml``: Configuration for automatic documentation builds
