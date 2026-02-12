.. FLUXOS documentation master file, created by
   sphinx-quickstart on Thu Mar  4 19:33:42 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to OpenWQ's documentation!
==================================

OpenWQ is a flexible biogeochemical modeling framework written in C++ that takes a fundamentally different approach from traditional water quality models. It supports both kinetic reaction networks (via the native BGC-Flex engine) and equilibrium geochemistry (via PHREEQC integration).

**What makes OpenWQ different?**

Unlike conventional water quality models that embed rigid, hard-coded biogeochemical processes alongside simplified hydrology, OpenWQ is designed as a **coupler** that integrates with existing, well-tested hydrological and hydrodynamic models. This architecture solves three major challenges:

1. **No hard-coded chemistry**: Most water quality models require source code modifications to change reaction formulations. OpenWQ uses JSON-configurable reaction networks — users define arbitrary kinetic expressions without recompilation, enabling rapid hypothesis testing and model adaptation to diverse environments (permafrost, peatlands, estuaries, etc.).

2. **Leverages state-of-the-art hydrology**: Rather than implementing its own simplified hydrology (which often lags behind dedicated hydro-models), OpenWQ inherits water volumes, fluxes, and state variables from the host model. This means users benefit from advanced process representations (e.g., multilayer snowpack in CRHM, flexible hydrological structures in SUMMA, continental-scale routing in mizuRoute) while gaining full biogeochemical capabilities.

3. **Enables structural uncertainty analysis**: The flexible architecture supports systematic comparison of reaction formulations, parameter sensitivity studies, and cross-model evaluations — essential for quantifying uncertainty in water quality predictions under climate change.

OpenWQ adapts to the host model's spatial and temporal discretization, supporting 1D, 2D, or 3D domains on structured or unstructured meshes. It is open-source, HPC-ready, and designed for extensibility.

Key capabilities:

* **Flexible biogeochemistry**: User-defined reaction networks via JSON configuration, powered by the ExprTk expression engine
* **PHREEQC geochemistry**: Full geochemical modeling (speciation, mineral reactions, ion exchange, kinetics) via the PhreeqcRM library
* **Multiple solvers**: Forward Euler (explicit) and SUNDIALS CVode (adaptive multi-step) ODE solvers
* **Transport modules**: Advection-only or advection-dispersion transport of dissolved constituents
* **Sorption isotherms**: Freundlich and Langmuir isotherm models for dissolved-sorbed partitioning
* **Sediment transport**: HBV-based erosion and mobile-immobile fractionation models
* **Multi-model coupling**: Already coupled to SUMMA, CRHM, mizuRoute, and FLUXOS
* **High performance**: OpenMP parallelization with SIMD vectorization, MPI support for distributed computing

.. toctree::
   :maxdepth: 4

   1_0_Home
   2_0_Documentation
   3_0_Installation
   5_3_base_Hydro_models
   4_1_OpenWQ_IO
   4_2_2support_scripts
   4_2_2_Calibration
   2_4_AI_Assistant
   5_0_Developer
