Software Documentation
==================================

Architecture Overview
~~~~~~~~~~~~~~~~~~~~~~

OpenWQ is written in C++17 and uses a modular architecture organized into the following major components:

**Computational Engine** (``src/compute/``)

* ``SOLVER_run_driver.cpp`` -- Main solver orchestration, selects between Backward Euler and SUNDIALS CVode
* ``solver_backeuler/`` -- Backward Euler implicit solver with OpenMP parallelization
* ``solver_sundials/`` -- SUNDIALS CVode adaptive ODE solver integration

**Chemistry Modules** (``src/models_CH/``)

* ``bgc_flex/`` -- Flexible biogeochemical reaction networks using ExprTk expressions
* ``PHREEQC/`` -- PhreeqcRM integration for geochemical modeling (setup, runtime, configuration)

**Transport Modules** (``src/models_TD/``)

* ``native_adv/`` -- Pure advection transport (concentration-based mixing)
* ``native_advdisp/`` -- Advection-dispersion transport with configurable dispersion coefficients

**Sorption Modules** (``src/models_SI/``)

* ``freundlich/`` -- Freundlich isotherm model
* ``langmuir/`` -- Langmuir isotherm model

**Sediment Transport** (``src/models_TS/``)

* ``hype_hvb/`` -- HYPE_HBVSED: HBV-based sediment erosion model with land-use, soil-type, slope, precipitation, and monthly erosion factors
* ``hype_mmf/`` -- HYPE_MMF: Morgan-Morgan-Finney erosion model with rainfall kinetic energy, soil cohesion, erodibility, and transport capacity
* ``json_driver.cpp`` -- TS module JSON configuration loader (reads MODULE_CONFIG_FILEPATH and populates TS_module)
* ``TS_config_driver.cpp`` -- TS configuration driver
* ``TS_runspace_driver.cpp`` -- TS runtime driver

**Lagrangian Elements** (``src/models_LE/``)

* ``native_boundmix/`` -- Boundary mixing model for inter-compartment exchange

**I/O and Configuration** (``src/readjson/``, ``src/output/``)

* JSON parsing using a custom reader with C/C++ comment support
* Output in HDF5, CSV, or VTU formats

**Integration** (``src/couplercalls/``, ``src/global/``, ``src/initiate/``)

* Coupler interface for host model communication
* Global state management
* Initialization routines


Build System
~~~~~~~~~~~~~~

OpenWQ uses CMake (minimum version 3.10) with support for multiple build targets:

* ``openwq`` -- Standalone
* ``summa_openwq`` -- Coupled to SUMMA
* ``crhm_openwq`` -- Coupled to CRHM
* ``mizuroute_openwq`` -- Coupled to mizuRoute
* ``mizuroute_lakes_openwq`` -- mizuRoute with lake routing
* ``mizuroute_lakes_cslm_openwq`` -- Full coupling with CSLM lake model
* ``fluxos_openwq`` -- Coupled to FLUXOS-OVERLAND

Two build modes are available:

* ``debug`` -- Debug symbols, warnings enabled (``-g -Wall -pedantic``)
* ``fast`` -- Optimized (``-O3 -march=native -funroll-loops``, Armadillo debug disabled)


Dependencies
~~~~~~~~~~~~~

* **Armadillo** -- C++ linear algebra (matrix/cube operations)
* **HDF5** -- Binary I/O for output files
* **SUNDIALS** -- CVode adaptive ODE solver (optional, required for SUNDIALS solver)
* **PhreeqcRM** -- PHREEQC geochemical engine (optional, required for PHREEQC module)
* **OpenMP** -- Shared-memory parallelization
* **LAPACK/BLAS** -- Linear algebra backend
* **NetCDF** -- Climate data I/O (for host model coupling)
* **MPI** -- Distributed-memory parallelization (for mizuRoute coupling)
* **VTK** -- Visualization output (optional)

See :doc:`Installation <3_0_Installation>` for detailed build instructions.


API Reference
~~~~~~~~~~~~~~

OpenWQ exposes four main API functions for host model coupling:

1. ``InitialConfig`` -- Initialize OpenWQ with compartment structure and threading
2. ``RunTimeLoopStart`` -- Begin a time step (pass volumes and simulation time)
3. ``RunSpaceStep`` -- Process inter-compartment water fluxes and transport
4. ``RunTimeLoopEnd`` -- Finalize time step (run solver, write output)

See :doc:`API documentation <5_3_00_APIs>` for function signatures and usage.
