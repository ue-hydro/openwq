Troubleshooting
==================================

This section covers common issues encountered when building and running OpenWQ.


Build Issues
~~~~~~~~~~~~~

**CMake cannot find Armadillo**

Ensure Armadillo is installed and that CMake can locate it. On Linux::

    sudo apt-get install libarmadillo-dev

On macOS::

    brew install armadillo

If installed in a non-standard location, set the ``ARMADILLO_ROOT`` environment variable.

**CMake cannot find SUNDIALS**

SUNDIALS is only required if you plan to use the SUNDIALS CVode solver. Install it via your package manager::

    # Linux
    sudo apt-get install libsundials-dev

    # macOS
    brew install sundials

**CMake cannot find PhreeqcRM**

PhreeqcRM is only required if you plan to use the PHREEQC geochemistry module. Follow the installation instructions from the `USGS PhreeqcRM repository <https://www.usgs.gov/software/phreeqcrm>`_.

**Linker errors with HDF5**

Ensure HDF5 is compiled with C++ support. On some systems, you may need to specify the HDF5 path explicitly::

    cmake -DHDF5_ROOT=/path/to/hdf5 ...

**Fortran compiler mismatch (mizuRoute coupling)**

When building the mizuRoute-coupled version, ensure that the Fortran and C++ compilers are compatible. Mixing compiler vendors (e.g., Intel Fortran with GNU C++) can cause linking errors.


Runtime Issues
~~~~~~~~~~~~~~~

**"Number of components does not match" error**

When using PHREEQC, this error means the number of chemical species defined in the OpenWQ configuration does not match the number of components discovered by PhreeqcRM's ``FindComponents()``. Enable ``RUN_MODE_DEBUG`` in the master config to see the full component list and adjust your configuration accordingly.

**NaN or Inf values in output**

This typically indicates numerical instability. Try:

* Reducing the simulation time step
* Switching to the SUNDIALS solver (which uses adaptive time stepping)
* Checking that initial conditions are physically reasonable
* Verifying that source/sink magnitudes are not excessively large

**Segmentation fault at startup**

Common causes:

* Missing or incorrectly formatted JSON configuration files
* File paths in JSON that do not exist
* Mismatched compartment names between the host model and OpenWQ configuration

**Output files are empty**

Verify that:

* The ``RESULTS_FOLDERPATH`` directory exists and is writable
* The ``CHEMICAL_SPECIES`` list in the master config contains valid species names
* The ``COMPARTMENTS_AND_CELLS`` section references valid compartment names
* The ``TIMESTEP`` output interval is not larger than the simulation duration

**OpenMP thread count warning**

If ``USE_NUM_THREADS`` is set higher than the number of available cores, OpenMP may produce warnings. Set it to ``"all"`` to automatically use all available cores, or specify a reasonable integer value.


JSON Configuration Issues
~~~~~~~~~~~~~~~~~~~~~~~~~~~

**JSON parse error**

OpenWQ's JSON parser supports C/C++ style comments (``//`` and ``/* */``), but the JSON structure must be valid. Common issues:

* Missing commas between key-value pairs
* Trailing commas after the last element in an object or array
* Unquoted string values
* Mismatched brackets or braces

Use a JSON validator to check your configuration files before running.

**"Unknown compartment name" error**

Compartment names in OpenWQ configuration files must exactly match the names exported by the host hydrological model. These are case-sensitive. Enable ``RUN_MODE_DEBUG`` to see the list of available compartment names.
