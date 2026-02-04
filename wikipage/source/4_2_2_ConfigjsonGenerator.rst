Configuration JSON Generator
========================================

The Configuration JSON file defines the initial conditions and compartment mappings for OpenWQ simulations.
While this file can be created manually (see :doc:`Configuration file <4_1_2Config>`), the following guidelines help generate it correctly.


Manual creation
~~~~~~~~~~~~~~~~

1. Identify the compartments exported by your host hydrological model (e.g., ``RUNOFF``, ``SOIL_LAYER``, ``SNOW``)
2. Define the chemical species to simulate (e.g., ``NO3``, ``NH4``, ``DOC``)
3. Set initial conditions for each species in each compartment
4. Specify the units for initial concentrations

The JSON structure follows the format described in :doc:`Configuration file <4_1_2Config>`.

.. tip::

    Enable ``RUN_MODE_DEBUG`` in the :doc:`Master configuration <4_1_1Master>` to see the list of compartment names and dimensions exported by the host model. This helps ensure that compartment names in the configuration file match exactly.
