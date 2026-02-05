Theoretical Foundation
==================================


Conservation Equations
~~~~~~~~~~~~~~~~~~~~~~

There are two options available for computing the conservative transport of solutes (dissolved in water) or fully suspended sediments.

Since OpenWQ operates as a plug-in water quality module coupled to external host hydrological models, it does not have direct access to the spatial grid topology or neighbouring cells. Instead, the host model provides pairwise water flux exchanges between source and recipient cells. The transport equations are therefore formulated in terms of **discrete mass fluxes between cell pairs** rather than continuous PDEs.


Advection
"""""""""

The advective mass flux of a dissolved species from a source cell to a recipient cell is:

.. math::

    F_{adv} = \frac{Q_{s \rightarrow r}}{W_s} \cdot m_s

where :math:`Q_{s \rightarrow r}` is the water flux from source to recipient :math:`[MT^{-1}]`, :math:`W_s` is the water mass in the source cell :math:`[M]`, and :math:`m_s` is the chemical mass in the source cell :math:`[M]`. This is equivalent to transporting mass at the source cell concentration: :math:`F_{adv} = Q_{s \rightarrow r} \cdot C_s`.


Dispersion (Fickian approximation)
"""""""""""""""""""""""""""""""""""

Since only pairwise source-recipient cell connections are available (without access to neighbouring cells for computing second-order spatial derivatives), the dispersion term is approximated using a **Fickian dispersive flux between cell pairs**:

.. math::

    F_{disp} = D_{eff} \cdot (C_s - C_r) \cdot W_s

where:

* :math:`F_{disp}` is the dispersive mass flux :math:`[MT^{-1}]`
* :math:`C_s = m_s / W_s` is the concentration in the source cell
* :math:`C_r = m_r / W_r` is the concentration in the recipient cell
* :math:`W_s` is the water mass in the source cell :math:`[M]`
* :math:`D_{eff}` is the effective dispersion rate :math:`[T^{-1}]`

The effective dispersion rate is pre-computed from user-provided parameters:

.. math::

    D_{eff} = \frac{D_{avg}}{L^2}

where :math:`D_{avg} = (D_x + D_y + D_z) / 3` is the average of the directional dispersion coefficients :math:`[L^2T^{-1}]` (averaged because the orientation of each source-recipient connection is not known), and :math:`L` is the user-provided characteristic length :math:`[L]` representing the distance between cell centers.

This formulation is **bidirectional**: the dispersive flux can move mass in either direction depending on the concentration gradient, smoothing concentration differences regardless of flow direction.


Total transport flux
""""""""""""""""""""

The first option (``OPENWQ_NATIVE_TD_ADVDISP``) computes the total transport as the sum of advective and dispersive fluxes:

.. math::

    F_{total} = F_{adv} + F_{disp}

The second option (``NATIVE_TD_ADV``) accounts for transport through advection only:

.. math::

    F_{total} = F_{adv}


Biogeochemistry: Reaction networks
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These expressions are flexible and defined by the user.
They are powered by the comprehensive C++ Mathematical Expression Toolkit Library (``ExprTk``) developed by `Arash Partow (1999-2020) <http://www.partow.net/programming/exprtk/index.html>`_.
The implementation of ``ExprTk`` in OpenWQ is simple to use and provides an extremely efficient run-time mathematical expression parser and evaluation engine.

``ExprTk`` supports numerous forms of functional, logical and vector processing semantics and is very easily extendible. In computational biogeochemistry, reaction kinetics usually take the general form of :math:`1^{st}`, :math:`2^{nd}` or :math:`3^{rd}`-multispecies `reaction kinetics <https://en.wikipedia.org/wiki/Chemical_kinetics>`_.

.. math::
    :nowrap:

    \begin{equation}
    \frac{\partial{c_a}}{\partial t} = - k \cdot B \cdot c_a
    \label{1storder_eq}
    \end{equation}

.. math::
    :nowrap:

    \begin{equation}
    \frac{\partial{c_a}}{\partial t} = - k \cdot B \cdot c_a^2
    \label{2ndorder_eq}
    \end{equation}

.. math::
    :nowrap:

    \begin{equation}
    \frac{\partial{c_a}}{\partial t} = - k \cdot B \cdot c_a^2 * c_b
    \label{3rdorder_eq}
    \end{equation}

where :math:`c_a` and :math:`c_b` :math:`[ML^{-3}]` are the concentrations of chemical species :math:`a` and :math:`b`, parameter/variable :math:`B` represents weather/hydrological dependencies (such as soil moisture and temperature) and :math:`k` is the reaction rate :math:`[ML^{-3}T^{-1}]`. The reaction rate :math:`k` can be provided as the reaction rate (usually using standard maximum at a reference temperature, often :math:`20^oC`) or using expressions that can include relationships with the hydrological/weather dependency variables/parameters.


PHREEQC Geochemistry
~~~~~~~~~~~~~~~~~~~~~~

As an alternative to the flexible BGC reaction networks, OpenWQ integrates the PHREEQC geochemical engine via PhreeqcRM.
PHREEQC computes equilibrium speciation and kinetic reactions based on thermodynamic databases, solving the full set of mass-action equations for aqueous species, minerals, surfaces, and exchangers.

For a given solution, PHREEQC solves:

.. math::

    K_i = \prod_j a_j^{\nu_{ij}}

where :math:`K_i` is the equilibrium constant for reaction :math:`i`, :math:`a_j` is the activity of species :math:`j`, and :math:`\nu_{ij}` is the stoichiometric coefficient.

This approach is more computationally expensive than the BGC-Flex approach, but provides thermodynamically consistent speciation calculations, which are essential for problems involving pH-dependent reactions, mineral interactions, or complex aqueous chemistry. See :doc:`PHREEQC configuration <4_1_3bPHREEQC>` for setup details.


Sorption Isotherms
~~~~~~~~~~~~~~~~~~~~~~

OpenWQ includes two sorption isotherm models for simulating partitioning between dissolved and sorbed phases. Both use a **kinetic approach**: the equilibrium sorbed concentration is computed from the isotherm equation, and mass is transferred at a user-defined kinetic rate.

**Freundlich isotherm** -- Models nonlinear sorption:

.. math::

    q = K_{fr} \cdot C^{N_{fr}}

where :math:`q` is the sorbed concentration :math:`[ML^{-3}_{soil}]`, :math:`C` is the dissolved concentration :math:`[ML^{-3}]`, :math:`K_{fr}` is the Freundlich coefficient, and :math:`N_{fr}` is the Freundlich exponent. The equilibrium dissolved concentration is found via Newton-Raphson iteration on the mass conservation equation :math:`C_{eq} \cdot V + K_{fr} \cdot C_{eq}^{N_{fr}} \cdot \rho \cdot L = m_{total}`.

**Langmuir isotherm** -- Models sorption with a finite number of binding sites:

.. math::

    q = \frac{q_{max} \cdot K_L \cdot C}{1 + K_L \cdot C}

where :math:`q_{max}` is the maximum adsorption capacity :math:`[MM^{-1}_{soil}]` and :math:`K_L` is the Langmuir equilibrium constant :math:`[L^3M^{-1}]`. The equilibrium concentration is solved analytically via the quadratic formula.

**Kinetic adsorption/desorption** -- For both isotherms, the mass flux from dissolved to sorbed phase is:

.. math::

    \Delta q = (q_{eq} - q_{current}) \cdot (1 - e^{-K_{adsdes} \cdot \Delta t})

.. math::

    F_{sorption} = \Delta q \cdot \rho \cdot L

where :math:`K_{adsdes}` is the kinetic adsorption/desorption rate :math:`[T^{-1}]`, :math:`\rho` is bulk density :math:`[ML^{-3}]`, and :math:`L` is layer thickness :math:`[L]`. This flux is applied as a sink/source on the dissolved mass. See :doc:`Sorption Isotherm configuration <4_1_7SI>` for setup details.


Sediment Transport
~~~~~~~~~~~~~~~~~~~~~~

OpenWQ includes two sediment transport erosion models based on the HYPE hydrological model framework. Both models compute erosion rates and transport sediments between compartments along a user-defined direction.

**HYPE_HBVSED** -- HBV-based sediment erosion model

This model computes soil erosion as a function of precipitation, slope, soil properties, and land use. The erosion rate is computed as:

.. math::

    \text{erosion} = \text{eroindex} \cdot \text{eroindexpar} \cdot \text{lusepar} \cdot \text{soilpar} \cdot \text{slope}^{\text{slopepar}} \cdot \text{precip}^{\text{precexppar}} \cdot \text{erodmonth}

where:

* :math:`\text{eroindex}` -- erosion index (spatially varying)
* :math:`\text{eroindexpar}` -- scaling parameter for the erosion index
* :math:`\text{lusepar}` -- land-use dependent soil erosion factor
* :math:`\text{soilpar}` -- soil-type dependent erosion factor
* :math:`\text{slope}` -- basin slope
* :math:`\text{slopepar}` -- slope erosion exponent
* :math:`\text{precip}` -- precipitation intensity
* :math:`\text{precexppar}` -- precipitation erosion exponent
* :math:`\text{erodmonth} = 1 + \text{monthpar}(m)` -- monthly erosion factor, where :math:`\text{monthpar}` is a 12-value array (Jan--Dec) allowing seasonal variation

Erosion only occurs when precipitation exceeds zero and the inhibiting compartment (e.g., snow layer) has no water. Sediments are removed from the source compartment and added to the recipient compartment.

**HYPE_MMF** -- Morgan-Morgan-Finney erosion model

This model uses the MMF approach to compute erosion based on kinetic energy from rainfall, soil cohesion, erodibility, and ground cover fractions. The erosion rate depends on:

.. math::

    F = K \cdot KE \cdot 10^{-3}

.. math::

    H = \text{cohesion}^{0.5} \cdot \text{cropcover} \cdot e^{-\text{erodibility} \cdot SR^{\text{sreroexp}}}

.. math::

    TC = \text{trans}_1 \cdot (SR)^{\text{trans}_2} \cdot \sin(\text{slope}) \cdot (1 - \text{groundcover})

where:

* :math:`F` -- soil particle detachment by raindrop impact
* :math:`H` -- resistance to erosion (soil cohesion and crop cover)
* :math:`TC` -- transport capacity of overland flow
* :math:`KE` -- total kinetic energy of rainfall (direct throughfall + leaf drainage)
* :math:`SR` -- surface runoff depth
* :math:`\text{cohesion}` -- soil cohesion (kPa)
* :math:`\text{erodibility}` -- soil erodibility (g/J)
* :math:`\text{sreroexp}` -- surface runoff erosion exponent
* :math:`\text{cropcover}` + :math:`\text{groundcover}` -- must sum to 1.0
* :math:`\text{trans}_1`, :math:`\text{trans}_2` -- transport capacity factors
* :math:`\text{slope}` -- basin slope

The pseudo day-of-year is dynamically computed from the simulation time and used to estimate kinetic energy from rainfall. Erosion only occurs when the inhibiting compartment has no water and crop/ground cover fractions are valid.

Both models support spatially-varying parameters via JSON or ASCII input formats, allowing cell-by-cell parameter overrides over the default values.


Numerical Solvers
~~~~~~~~~~~~~~~~~~~~~~

OpenWQ provides two numerical solvers for integrating the water quality ODEs:

* **Forward Euler**: A first-order explicit method. Simple and fast. See :doc:`Solvers <4_1_9Solvers>`.
* **SUNDIALS CVode**: An adaptive multi-step method (Adams or BDF) with variable order and automatic step-size control. Provides higher-order accuracy for stiff and non-stiff problems. See :doc:`Solvers <4_1_9Solvers>`.
