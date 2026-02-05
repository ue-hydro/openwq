Example of Case Studies
==================================

Synthetic Test Cases
~~~~~~~~~~~~~~~~~~~~~

A total of 10 synthetic test cases have been developed to verify OpenWQ's outputs against analytical solutions.
These tests cover batch reactors with different reaction orders and advection-dispersion transport scenarios.
See :doc:`Synthetic tests <synthetic_tests>` for the full list.

The JSON input files for each synthetic test can be obtained from: https://github.com/ue-hydro/synthetic_tests


PHREEQC River Geochemistry
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A test case demonstrating PHREEQC geochemical modeling coupled with mizuRoute river routing is available in the ``test_case_phreeqc/`` directory of the mizuRoute-OpenWQ repository.

This example simulates the transport and geochemical evolution of Ca, Mg, and Na through a river network, including:

* Equilibrium speciation of dissolved species
* Advective transport between river reaches
* Temperature-dependent reaction calculations
* HDF5 output of species concentrations

Configuration files include:

* ``openWQ_master.json`` -- Master configuration with Forward Euler solver
* ``openwq_in/openWQ_config.json`` -- Initial conditions (Ca, Mg, Na in mol/kgw)
* ``openwq_in/openWQ_MODULE_PHREEQC.json`` -- PHREEQC module settings
* ``openwq_in/openWQ_MODULE_TD.json`` -- Pure advection transport
* ``openwq_in/phreeqc_river.pqi`` -- PHREEQC input file
* ``openwq_in/phreeqc.dat`` -- Thermodynamic database


Nutrient Cycling (BGC-Flex)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

OpenWQ's flexible reaction network approach is well suited for nutrient cycling simulations.
Using the BGC-Flex module, users can define nitrogen and phosphorus cycle reactions with temperature-dependent kinetics.

Example applications include:

* **Nitrogen cycle**: Nitrification, denitrification, mineralization, and plant uptake with temperature and moisture dependencies
* **Phosphorus cycle**: Sorption-desorption, mineralization, and particulate transport
* **Dissolved oxygen**: BOD decay, reaeration, and sediment oxygen demand

These reaction networks are defined entirely in JSON configuration files. See :doc:`Modules <4_1_3Modules>` for the file format and the :doc:`Theoretical Foundation <2_1_Theoretical>` for the underlying equations.
