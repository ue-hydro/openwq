Apptainer/Singularity (HPC)
==================================

Apptainer (formerly Singularity) is the recommended containerization solution for High-Performance Computing (HPC) environments.
Unlike Docker, Apptainer is designed for multi-user HPC systems where root privileges are not available.

The OpenWQ repository includes a pre-configured definition file (``containers/openwq_apptainer.def``) that builds a complete environment with all dependencies.


Building the container
~~~~~~~~~~~~~~~~~~~~~~~

1. Ensure Apptainer is installed on your system (most HPC clusters have it available as a module)::

    module load apptainer    # or: module load singularity

2. Navigate to the containers directory inside the OpenWQ repository::

    cd /path/to/openwq/containers/

3. Build the container image (.sif file)::

    apptainer build openwq.sif openwq_apptainer.def

   This process may take 15-30 minutes as it downloads and compiles all dependencies.

4. The resulting ``openwq.sif`` file can be copied to any HPC system with Apptainer installed.


What's included in the container
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The definition file (``openwq_apptainer.def``) builds from Ubuntu 24.04 and installs:

**Development tools:**

* GCC/G++ and GFortran compilers
* CMake build system
* GDB debugger
* Git version control

**Scientific libraries:**

* **LAPACK/BLAS/OpenBLAS** -- Linear algebra
* **HDF5** -- Hierarchical Data Format I/O
* **NetCDF/NetCDF-Fortran** -- Network Common Data Form
* **Boost** -- C++ libraries

**OpenWQ-specific dependencies:**

* **PhreeqcRM** -- PHREEQC geochemical reaction module (built from usgs-coupled/phreeqcrm_use_cmake at commit cf03df1)
* **SUNDIALS 7.1.1** -- Suite of nonlinear and differential/algebraic equation solvers with Fortran interface
* **Armadillo 10.6.0** -- C++ linear algebra library with HDF5 support enabled
* **ParallelIO (PIO)** -- Parallel I/O library for NetCDF/PnetCDF (required for mizuRoute lake coupling)


Running with Apptainer
~~~~~~~~~~~~~~~~~~~~~~~

**Interactive shell:**

.. code-block:: bash

    apptainer shell openwq.sif

**Execute a command:**

.. code-block:: bash

    apptainer exec openwq.sif ./mizuroute_openwq control_file.txt

**Bind directories** (to access files outside the container):

.. code-block:: bash

    apptainer exec --bind /path/to/data:/data openwq.sif ./mizuroute_openwq /data/control_file.txt


Compiling inside the container
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To compile OpenWQ or a coupled model using the container:

.. code-block:: bash

    # Enter the container with your source code mounted
    apptainer shell --bind /path/to/source:/code openwq.sif

    # Inside the container, set environment variables
    export PhreeqcRM_DIR=/usr/local/phreeqcrm
    export SUNDIALS_DIR=/usr/local/sundials

    # Navigate to code and build
    cd /code/route/build/openwq/openwq/
    cmake -DHOST_MODEL_TARGET=mizuroute_lakes_cslm_openwq -DCMAKE_BUILD_TYPE=fast .
    make -j 4

A helper script (``utils/compile_summa_apptainer.sh``) is provided for SUMMA compilation.


HPC job submission
~~~~~~~~~~~~~~~~~~~

Example SLURM batch script:

.. code-block:: bash

    #!/bin/bash
    #SBATCH --job-name=openwq_run
    #SBATCH --nodes=1
    #SBATCH --ntasks=4
    #SBATCH --time=02:00:00
    #SBATCH --mem=16G

    module load apptainer

    # Run the model
    apptainer exec --bind $PWD:/run openwq.sif \
        /run/mizuroute_openwq /run/control_file.txt

    # Or with MPI
    srun apptainer exec --bind $PWD:/run openwq.sif \
        /run/mizuroute_openwq /run/control_file.txt


Troubleshooting
~~~~~~~~~~~~~~~~

**"FATAL: could not open image" error:**

Ensure the .sif file has read permissions and the path is correct.

**Library not found at runtime:**

The container includes all dependencies. If you see missing library errors, you may be accidentally using host libraries. Check your ``LD_LIBRARY_PATH`` environment variable.

**Out of memory during build:**

The container build requires substantial memory. On systems with limited resources, build on a high-memory node or workstation, then transfer the .sif file.
