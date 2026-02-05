Numerical Solvers
==================================

OpenWQ provides two numerical solvers for integrating the water quality ordinary differential equations (ODEs) at each time step. The solver is selected in the :doc:`Master configuration <4_1_1Master>`.


Solver selection
~~~~~~~~~~~~~~~~~

The solver is specified in the master configuration file under ``COMPUTATIONAL_SETTINGS``:

.. code-block:: json

    {
        "COMPUTATIONAL_SETTINGS": {
            "SOLVER": "FORWARD_EULER",
            "USE_NUM_THREADS": 4
        }
    }

Available options:

+------------------------+--------------------------------------+--------------------------------------------------+
| Value                  | Solver                               | Description                                      |
+========================+======================================+==================================================+
| ``"FORWARD_EULER"``    | Forward Euler                        | Simple explicit solver; fast, first-order         |
+------------------------+--------------------------------------+--------------------------------------------------+
| ``"SUNDIALS"``         | SUNDIALS CVode                       | Adaptive multi-step solver; higher accuracy       |
+------------------------+--------------------------------------+--------------------------------------------------+


Forward Euler (FORWARD_EULER)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Forward Euler solver is a first-order explicit method that updates chemical mass at each time step by summing all source/sink terms:

.. math::

    M^{n+1} = M^{n} + \Delta t \left( \dot{M}_{ic} + \dot{M}_{ss} + \dot{M}_{ewf} + \dot{M}_{chem} + \dot{M}_{transp} \right)

where:

* :math:`M^{n}` is the chemical mass at the current time step
* :math:`\dot{M}_{ic}` is the initial condition contribution (first time step only)
* :math:`\dot{M}_{ss}` is the source/sink rate
* :math:`\dot{M}_{ewf}` is the external water flux contribution
* :math:`\dot{M}_{chem}` is the biogeochemical reaction rate
* :math:`\dot{M}_{transp}` is the transport rate (advection and dispersion)

The Forward Euler solver processes all compartments and chemical species in parallel using OpenMP with dynamic scheduling.
Inner loops are further optimized with SIMD vectorization.

**Advantages**:

* Simple implementation with low computational overhead per step
* Good performance for problems with moderate dynamics and large spatial domains

**Limitations**:

* First-order accuracy; may require small time steps for high accuracy
* No adaptive time stepping
* Conditionally stable; time step must be small enough for the problem dynamics


SUNDIALS CVode
~~~~~~~~~~~~~~~

The `SUNDIALS <https://computing.llnl.gov/projects/sundials>`_ CVode solver is an advanced adaptive multi-step ODE solver. It provides higher-order accuracy and automatic time step control.

CVode supports both Adams methods (for non-stiff problems) and BDF methods (for stiff problems) with variable-order, variable-step-size control (up to order 5 for BDF, order 12 for Adams).

The right-hand side (RHS) function computes the total flux rate for each chemical species:

1. Chemistry derivatives are reset to zero
2. The current SUNDIALS state vector is copied to the OpenWQ chemical mass arrays
3. The chemistry driver computes reaction rates for all species
4. The total flux derivative is assembled from chemistry, transport, and source/sink contributions

**Advantages**:

* Higher-order accuracy (adaptive order selection)
* Adaptive time stepping automatically adjusts to problem stiffness
* Better accuracy for problems with rapidly changing concentrations

**Limitations**:

* Higher computational cost per step compared to Forward Euler
* Requires the SUNDIALS library to be installed (see :doc:`Installation <3_0_Installation>`)


Performance optimizations
~~~~~~~~~~~~~~~~~~~~~~~~~~

Both solvers leverage several optimizations for high performance:

* **OpenMP parallelization**: Compartment-chemical pairs are distributed across threads with dynamic scheduling for load balancing
* **SIMD vectorization**: Innermost loops use ``#pragma omp simd`` for automatic vectorization
* **Memory optimization**: Bulk data transfers use ``std::memcpy`` instead of element-wise copies
* **Cached indices**: Pre-computed compartment dimensions and species offsets avoid repeated lookups

The number of threads is controlled by the ``USE_NUM_THREADS`` setting in the master configuration.
Set to ``"all"`` to use all available CPU cores, or specify an integer value.


Solver choice guidelines
~~~~~~~~~~~~~~~~~~~~~~~~~~

+------------------------------------------+----------------------+--------------------------------------+
| Scenario                                 | Recommended solver   | Reason                               |
+==========================================+======================+======================================+
| Simple kinetic reactions (BGC-Flex)      | FORWARD_EULER        | Fast and sufficient                  |
+------------------------------------------+----------------------+--------------------------------------+
| Complex PHREEQC geochemistry             | FORWARD_EULER or     | Depends on problem stiffness         |
|                                          | SUNDIALS             |                                      |
+------------------------------------------+----------------------+--------------------------------------+
| Rapidly varying concentrations           | SUNDIALS             | Adaptive stepping captures dynamics  |
+------------------------------------------+----------------------+--------------------------------------+
| Large networks with many species         | FORWARD_EULER        | Lower per-step overhead              |
+------------------------------------------+----------------------+--------------------------------------+
| High-accuracy requirements               | SUNDIALS             | Higher-order methods available       |
+------------------------------------------+----------------------+--------------------------------------+
| Quick exploratory runs                   | FORWARD_EULER        | Simpler, faster setup                |
+------------------------------------------+----------------------+--------------------------------------+
