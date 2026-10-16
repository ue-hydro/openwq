SUMMA with internally coupled mizuRoute
=======================================

Since October 2026 the ``develop`` branch of `CH-Earth/summa <https://github.com/CH-Earth/summa>`_ runs mizuRoute inside
SUMMA's own time loop and carries the OpenWQ coupling for that pair (pull request #655; the ue-hydro fork is
`summa-mizuroute-openwq <https://github.com/ue-hydro/summa-mizuroute-openwq>`_). One OpenWQ instance holds the land
compartments of SUMMA and the river reaches of mizuRoute, so species are balanced in parallel over land and river and solute
reaches the river in the same step as the water.

What OpenWQ sees
~~~~~~~~~~~~~~~~

+---------------------------+-----------------------------+----------------------------------------------------------+
| compartment               | cells                       | content                                                  |
+===========================+=============================+==========================================================+
| ``SCALARCANOPYWAT``       | (HRU, domain) columns x 1   | canopy storage                                           |
+---------------------------+-----------------------------+----------------------------------------------------------+
| ``ILAYERVOLFRACWAT_SNOW`` | columns x max snow layers   | snow layers                                              |
+---------------------------+-----------------------------+----------------------------------------------------------+
| ``RUNOFF``                | columns x 1                 | liquid water reaching the surface in the step            |
+---------------------------+-----------------------------+----------------------------------------------------------+
| ``ILAYERVOLFRACWAT_SOIL`` | columns x max soil layers   | soil layers                                              |
+---------------------------+-----------------------------+----------------------------------------------------------+
| ``SCALARAQUIFER``         | columns x 1                 | aquifer (basin-wide aquifer: first column of the GRU)    |
+---------------------------+-----------------------------+----------------------------------------------------------+
| ``RUNOFF_TO_STREAM``      | columns x 1                 | runoff held while SUMMA routes it within the GRU         |
+---------------------------+-----------------------------+----------------------------------------------------------+
| ``ILAYERVOLFRACWAT_LAKE`` | columns x max lake layers   | lake layers, only if a domain has them                   |
+---------------------------+-----------------------------+----------------------------------------------------------+
| ``RIVER_NETWORK_REACHES`` | reaches x 1                 | mizuRoute reaches                                        |
+---------------------------+-----------------------------+----------------------------------------------------------+

Cell ids are ``<hruId>_z<layer>`` (``<hruId>_d<domain>_z<layer>`` when an HRU has several domains) and ``<segId>`` for the
reaches. External water fluxes: ``PRECIP`` and ``GLACIER_ICE_MELT``. Flux-concentration exports: ``scalarRunoffVol_m3``,
``averageRoutedRunoff``, ``scalarTotalRunoff`` and ``Qlocal_out`` (reach outflow). Dependencies for the kinetic expressions:
``SM``, ``Tair_K``, ``Tsoil_K``, ``SWrad_Wm2``, ``cellArea_m2``, ``T``.

Land to river mapping and its checks
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Each GRU delivers its routed runoff (the same value SUMMA hands to mizuRoute) to the reaches in proportion to the inflow each
reach receives from that GRU per unit runoff; the map is built once at initialization through mizuRoute's own ``remap_runoff``
and ``basin2reach``, so it follows whatever hydrofabric, remapping file and ``hw_drain_point`` the run uses. Each reach then
mixes its start-of-step storage, upstream inflow, lateral inflow and any water-management injection, and passes the outflow on,
using the volumes mizuRoute accumulated over its sub-steps. The coupler verifies this on every run:

* at initialization: runoff elements that reach no reach, the river-network area of each GRU against its SUMMA area
  (tolerance 1 %), several active routing methods (water quality follows the first);
* at the end of the run: runoff routed by SUMMA versus lateral inflow received by the reaches, outlet outflow, end storage,
  and the summed reach budget residual, which should be at round-off;
* ``openwq_compartments.json`` is written next to the outputs with the compartment names and indices.

Verified on the Bow test case (kinematic wave, Muskingum-Cunge, diffusive wave, remapping), a two-HRU variant with lateral
flow between HRUs, and the Athabasca River (69 GRUs, 69 reaches, 17 years): tracer closure to 1e-8 kg and water check at
round-off. Lake, wetland and glacier domains have no test case yet; water management is not exposed by the coupled TOML.

Build and run
~~~~~~~~~~~~~

.. code-block:: bash

    git clone -b develop --recurse-submodules https://github.com/CH-Earth/summa.git
    cd summa/build/source/openwq
    git clone -b develop https://github.com/ue-hydro/openwq.git
    cd openwq
    cmake -DHOST_MODEL_TARGET=summa-mizuroute_openwq -DCMAKE_BUILD_TYPE=Release .
    make -j 2
    # bin/summa-mizuroute_openwq_Release
    cd <folder with openWQ_master.json>
    summa-mizuroute_openwq_Release -m fileManager.txt -c mizuroute.toml

The target builds the SUMMA tree as an external project (``-DUSE_OPENWQ=ON -DUSE_MIZUROUTE=ON``, SUNDIALS according to
``SOLVER_TYPE``) and needs a shared PhreeqcRM (``build/source/openwq/_deps`` or ``-DSUMMA_PHREEQCRM_DIR``). The TOML must set
``[simulation] use_mizuroute = true``; the internal coupling offers the kinematic wave, Muskingum-Cunge and diffusive wave
methods. The run is a single process. ``build/source/openwq/README.md`` in the SUMMA tree has the details.

Supporting scripts
~~~~~~~~~~~~~~~~~~

* Configuration: ``hostmodel = "summa-mizuroute"`` with ``mizuroute_config_path`` in the model config template
  (:doc:`4_2_2_ConfigjsonGenerator`); the generator adds ``RIVER_NETWORK_REACHES`` and ``-c <toml>`` to the run command.
* Results: land compartments join the basin shapefile and river outputs the river shapefile, both layers are clickable and
  stations are matched per layer (:doc:`4_2_2_OutputReadScript`).
* Calibration and scenarios: one model in ``model_chain``, river targets matched to reaches, sub-basin cascade and ML
  regionalization with explicit compartment indices (:doc:`4_2_2_Calibration`, :doc:`4_5_Scenarios`).
