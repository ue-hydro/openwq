.. FLUXOS documentation master file, created by
   sphinx-quickstart on Thu Mar  4 19:33:42 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to OpenWQ's documentation!
==================================

OpenWQ is a generalised biogeochemical reaction-network model written in C++.
It was designed as a coupler for integration with existing hydro-models to enable flexible, multi-scale, multi-chemistry simulations.
OpenWQ adapts to the host hydro-model spatial and temporal configurations, which can include 1D, 2D or 3D spatial discretizations based on structured or unstructured meshes.
It aims to be flexible, fast, free, and extendable.

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
   5_0_Developer
   6_0_Contact
