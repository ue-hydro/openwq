Sink-Source JSON Generator
========================================

The Sink-Source JSON file defines external chemical loads (point sources, non-point sources, atmospheric deposition) and removal processes (water treatment, abstraction) for OpenWQ simulations.
While this file can be created manually (see :doc:`Source and Sink configuration <4_1_4SS>`), the following guidelines help generate it correctly.


Manual creation
~~~~~~~~~~~~~~~~

1. Identify the type of source or sink:
   - **Source**: Adds chemical mass to a compartment cell (e.g., wastewater discharge, fertilizer application)
   - **Sink**: Removes chemical mass from a compartment cell (e.g., water treatment withdrawal)

2. Define the spatial location (compartment name, cell indices ix, iy, iz)

3. Define the temporal pattern:
   - Continuous load at a constant rate
   - Time-varying load with a schedule
   - Instantaneous pulse load

4. Specify the chemical species and mass loading rates with units

The JSON structure follows the format described in :doc:`Source and Sink configuration <4_1_4SS>`.

.. tip::

    Use the :doc:`Units reference <4_3_units>` to ensure consistent units between source/sink inputs and the rest of the simulation configuration.
