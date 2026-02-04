Existing couplings
==============================

OpenWQ has been coupled to the following hydrological and hydrodynamic models.
Each coupling follows the same general architecture: the host model provides water volumes and fluxes, and OpenWQ handles constituent transport and biogeochemical reactions.

To couple OpenWQ to a new model, see the :doc:`Coupler Guide <5_3_Coupler_guide>`.


**SUMMA**
~~~~~~~~~
Structure for Unifying Multiple Modeling Alternatives

SUMMA is a flexible hydrological modeling framework that allows users to select from multiple options for each model component (snow, soil, runoff generation, etc.). The SUMMA-OpenWQ coupling enables water quality simulations within SUMMA's multi-layer soil and snow compartments.

- summa: https://summa.readthedocs.io/en/latest/
- summa-openwq: https://github.com/ue-hydro/Summa-openWQ

**CRHM**
~~~~~~~~
Cold Regions Hydrological Model

CRHM is a platform for modeling hydrological processes in cold regions, including snow redistribution, sublimation, frozen soil dynamics, and prairie hydrology. The CRHM-OpenWQ coupling was the first integration of OpenWQ and enabled nutrient transport simulations in cold-region catchments.

- crhm: https://research-groups.usask.ca/hydrology/modelling/crhm.php
- crhm-openwq: https://github.com/ue-hydro/CRHM

**MizuRoute**
~~~~~~~~~~~~~
A river network routing tool for continental domain water resources applications

MizuRoute provides five numerical routing methods (IRF, KWT, KWE, Muskingum-Cunge, Diffusive Wave) for computing streamflow through river networks. The mizuRoute-OpenWQ coupling enables reactive transport of dissolved and particulate constituents through continental-scale river networks, with optional PHREEQC geochemistry and CSLM lake energy balance.

- mizuroute: https://mizuroute.readthedocs.io/en/main/
- mizuroute-openwq: https://github.com/ue-hydro/mizuRoute

**FLUXOS-OVERLAND**
~~~~~~~~~~~~~~~~~~~
A hydrodynamic-solute transport modelling tool suitable for basin-scale, event-based simulations

FLUXOS solves the 2D shallow water equations coupled with solute transport for event-based simulations of runoff, flooding, and contaminant transport in Prairie landscapes.

- fluxos: https://fluxos-cpp.readthedocs.io/en/latest/
- fluxos-openwq: https://github.com/ue-hydro/FLUXOS_cpp
