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

**SUMMA with internally coupled mizuRoute**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Land and river network in one executable

The current SUMMA development branch (CH-Earth/summa, ``develop``) calls mizuRoute inside its own time loop: each GRU hands its routed runoff to the river network every step, and only the kinematic wave, Muskingum-Cunge and diffusive wave methods are available. OpenWQ is coupled to this combination as ONE instance that holds the SUMMA land compartments (canopy, snow, surface runoff, soil, aquifer, runoff on its way to the stream, optional lake layers) and the river compartment ``RIVER_NETWORK_REACHES`` of mizuRoute, so species are balanced in parallel over land and river and solute reaches the river in the same step as the water. The GRU to reach mapping is built from mizuRoute's own remapping and hydrofabric, and the coupler verifies it on every run (area of each GRU seen by the river network against its SUMMA area at initialization, and a water check of the land to river transfer and of the reach budgets printed at the end of the run).

The executable is ``summa-mizuroute_openwq_<build type>`` (CMake target ``summa-mizuroute_openwq``) and is run as a single process with ``-m <SUMMA fileManager> -c <mizuRoute TOML>`` from the folder that holds ``openWQ_master.json``. The TOML must set ``[simulation] use_mizuroute = true``; in the OpenWQ configuration template the host model is ``hostmodel = "summa-mizuroute"``. The flux-concentration export ``Qlocal_out`` gives the reach outflow concentrations. The coupling files live in ``build/source/openwq`` and ``build/source/mizuroute`` of the SUMMA tree (``summa_openWQ.f90``, ``OpenWQ_hydrolink.cpp``, ``wq_exchange.f90``).

- summa (develop, internal mizuRoute): https://github.com/CH-Earth/summa

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
- mizuroute-openwq: https://github.com/ue-hydro/mizuRoute-OpenWQ

**FLUXOS**
~~~~~~~~~~
A hydrodynamic-solute transport modelling tool suitable for basin-scale, event-based simulations

FLUXOS solves the 2D shallow water equations coupled with solute transport for event-based simulations of runoff, flooding, and contaminant transport in Prairie landscapes.

- fluxos: https://fluxos.readthedocs.io
- fluxos-openwq: https://github.com/ue-hydro/FLUXOS_cpp
