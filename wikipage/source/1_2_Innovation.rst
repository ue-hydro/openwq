Motivation & Innovation
==================================

OpenWQ aims to enable a more flexible representation of local and global nuances in terrestrial biogeochemical cycling across regions, climate zones, and aquatic ecosystems through an adaptive reaction-network design. This concept was developed based on 10 years of experience with the development of hydro-biogeochemical models at ECCC's Research Branch and the University of Saskatchewan.

OpenWQ addresses three main challenges with existing biogeochemical and water quality models:

* ``Structural rigidity`` that hinders the effective use of models across landscapes and in complex, diverse environments (permafrost, peatlands, variable contributing areas) that require more investigative, open-ended, and interactive simulation approaches;
* ``Hydro-flux calculations`` embedded in these modelling tools are often outdated or limited when compared to dedicated, disciplinary hydro-models (e.g., hydrological or hydrodynamic modelling tools); and
* ``Uncertainty`` in model structure, process representation, and future scenarios (e.g., climate change) cannot be properly addressed because the current models are not suitable for adequate testing of modelling hypotheses.


Structural rigidity
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Water quality models have improved significantly over the last decades. There are now more water quality models than ever before, many of them quite well established.
However, a common limitation across most of these models is their rigid internal structure: the biogeochemical process representations are hard-coded, making it difficult to test alternative hypotheses or adapt the model to new environments without modifying the source code.

OpenWQ overcomes this with two complementary approaches:

1. **JSON-based flexible reaction networks (BGC-Flex)**: Users define arbitrary kinetic rate expressions using the ExprTk mathematical expression library. Reactions can depend on hydrological variables (soil moisture, temperature, etc.) provided by the host model. No recompilation is needed to change the chemistry.

2. **PHREEQC geochemical engine**: For advanced geochemistry, OpenWQ integrates the PhreeqcRM library, enabling equilibrium speciation, mineral dissolution/precipitation, surface complexation, ion exchange, and user-defined kinetic reactions. This is particularly useful for problems involving acid mine drainage, nutrient cycling with pH effects, or contaminant fate in complex aqueous systems.

.. image:: inca.png
    :width: 250 px

.. image:: SWAT.png
    :width: 250 px

.. image:: myLake.jpg
    :width: 250 px


Hydro-flux calculations
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Many water quality models implement their own simplified hydrological routines, which often lag behind the state of the art in dedicated hydro-models.
OpenWQ takes a fundamentally different approach: it does not compute its own hydrology. Instead, it **couples to existing, well-tested hydrological models** that provide water volumes, fluxes, and state variables.

This means that users benefit from the latest advances in hydrological science (e.g., advanced snow processes in CRHM, flexible model structures in SUMMA, continental routing in mizuRoute, 2D hydrodynamics in FLUXOS) without sacrificing water quality modeling capabilities.

The coupling architecture is designed to be host-model agnostic, requiring only four API calls to integrate OpenWQ with any hydrological model. See the :doc:`Coupler Guide <5_3_Coupler_guide>` for details.

.. image:: review_paper_process_representation.jpg
    :width: 350 px

.. image:: 1-3D.png
    :width: 350 px


Uncertainty and modelling hypothesis testing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The structural flexibility of OpenWQ enables systematic exploration of modeling hypotheses and uncertainty quantification:

* **Reaction network uncertainty**: Users can test different biogeochemical formulations (e.g., first-order vs. Michaelis-Menten kinetics for nutrient cycling) by simply modifying the JSON configuration files.
* **Parameter uncertainty**: The JSON-based configuration makes it straightforward to script ensemble simulations with varying parameter values for sensitivity analysis.
* **Structural uncertainty**: By coupling to different host models (e.g., SUMMA vs. CRHM vs. mizuRoute), users can assess how the choice of hydrological model affects water quality predictions.
* **Chemistry module comparison**: Users can compare flexible BGC-Flex reactions against PHREEQC-based geochemical calculations to evaluate the importance of detailed speciation.
* **Solver uncertainty**: Users can compare results between the Backward Euler and SUNDIALS CVode solvers to assess numerical accuracy.

.. image:: OpenWQ_structure.png
    :width: 350 px
