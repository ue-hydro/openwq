PHREEQC Geochemistry Module
==================================

OpenWQ integrates the `PHREEQC <https://www.usgs.gov/software/phreeqc-version-3>`_ geochemical modeling engine through the PhreeqcRM library.
This enables advanced equilibrium and kinetic geochemical calculations as an alternative to the flexible BGC reaction networks.

PHREEQC is selected by setting the ``BIOGEOCHEMISTRY`` module to ``PHREEQC`` in the :doc:`Master configuration <4_1_1Master>`.


Overview
~~~~~~~~

PHREEQC provides capabilities for:

* **Aqueous speciation**: Computing the distribution of dissolved species at chemical equilibrium
* **Mineral saturation**: Determining saturation indices and driving dissolution/precipitation
* **Surface complexation**: Modeling sorption onto mineral surfaces (e.g., iron oxides, clays)
* **Ion exchange**: Simulating cation exchange processes on soil surfaces
* **Kinetic reactions**: User-defined kinetic rate expressions (BASIC-language scripts in the .pqi file)
* **pH and pe calculations**: Computing solution pH and redox potential
* **Gas-phase equilibria**: Modeling gas dissolution and exsolution
* **Temperature-dependent reactions**: Equilibrium constants vary with temperature from the host model


Configuration
~~~~~~~~~~~~~~

The PHREEQC module is configured through a JSON file referenced in the master configuration:

.. code-block:: json

    {
        "FILEPATH": "openwq_in/phreeqc_river.pqi",
        "DATABASE": "openwq_in/phreeqc.dat",
        "COMPONENT_H2O": true,
        "BGC_GENERAL_MOBILE_SPECIES": ["H", "O", "Charge", "Ca", "Mg", "Na"],
        "TEMPERATURE": {
            "RIVER_NETWORK_REACHES": "air_temperature"
        }
    }

Configuration fields:

+-----------------------------------+-----------+--------------------------------------------------------------------+
| Field                             | Type      | Description                                                        |
+===================================+===========+====================================================================+
| ``FILEPATH``                      | string    | Path to the PHREEQC input file (``.pqi``)                          |
+-----------------------------------+-----------+--------------------------------------------------------------------+
| ``DATABASE``                      | string    | Path to the PHREEQC thermodynamic database file                    |
+-----------------------------------+-----------+--------------------------------------------------------------------+
| ``COMPONENT_H2O``                 | bool      | Include water as a component (default: ``true``)                   |
+-----------------------------------+-----------+--------------------------------------------------------------------+
| ``BGC_GENERAL_MOBILE_SPECIES``    | array     | Names of chemical species subject to advective transport           |
+-----------------------------------+-----------+--------------------------------------------------------------------+
| ``TEMPERATURE``                   | object    | Temperature source variable name for each compartment              |
+-----------------------------------+-----------+--------------------------------------------------------------------+

The JSON file supports C/C++ syntax for comments: single-line comment (``//``) or comment blocks (``/*`` and ``*/``).


PHREEQC input file (.pqi)
~~~~~~~~~~~~~~~~~~~~~~~~~~

The PHREEQC input file defines the initial solution chemistry, mineral assemblages, exchange compositions,
and any kinetic reactions using standard PHREEQC syntax. Example::

    SOLUTION 1
        units     mol/kgw
        temp      10.0
        pH        7.0
        Ca        0.001
        Mg        0.0005
        Na        0.002
    END

Refer to the `PHREEQC documentation <https://www.usgs.gov/software/phreeqc-version-3>`_ for the complete input file syntax and keyword reference.


Thermodynamic database
~~~~~~~~~~~~~~~~~~~~~~~

The thermodynamic database contains equilibrium constants and thermodynamic data for aqueous species, minerals, gases, and surface species. Several databases are available from USGS:

* ``phreeqc.dat`` -- General-purpose database (recommended starting point)
* ``wateq4f.dat`` -- Extended database with additional minerals and trace metals
* ``minteq.v4.dat`` -- MINTEQA2-compatible database
* ``llnl.dat`` -- Lawrence Livermore National Laboratory database (comprehensive)


Mobile species
~~~~~~~~~~~~~~~

The ``BGC_GENERAL_MOBILE_SPECIES`` array specifies which PHREEQC components are subject to advective transport between compartment cells. Only mobile species are transported by the OpenWQ transport module; immobile species remain in their original grid cell.

Species are identified by their component name (string), which must match the names returned by PHREEQC's ``FindComponents()`` function. The component list can be inspected in the OpenWQ debug output when ``RUN_MODE_DEBUG`` is enabled in the master configuration.


Temperature coupling
~~~~~~~~~~~~~~~~~~~~~

PHREEQC reaction rates and equilibrium constants are temperature-dependent. The ``TEMPERATURE`` field maps each host-model compartment to a temperature data source variable name.

Temperature values are passed to PhreeqcRM at each time step for accurate geochemical calculations.
The variable name must match a variable exported by the host hydrological model (e.g., ``"air_temperature"`` or ``"soil_temperature"``).


Initialization process
~~~~~~~~~~~~~~~~~~~~~~~~

During model initialization, the PHREEQC module performs the following steps:

1. **Grid calculation**: Computes the total number of grid cells across all compartments (``nxyz = nx * ny * nz``)
2. **PhreeqcRM creation**: Creates a PhreeqcRM instance with the computed grid size and the configured number of threads
3. **Database loading**: Loads the thermodynamic database file
4. **Input file execution**: Runs the PHREEQC input file to set up initial solution chemistry
5. **Component discovery**: Identifies all chemical components using ``FindComponents()``
6. **Verification**: Confirms that the number of discovered components matches the OpenWQ configuration


Example PHREEQC input file
~~~~~~~~~~~~~~~~~~~~~~~~~~~

A complete PHREEQC input file (``.pqi``) defines the geochemical system. Here is a more detailed example for river water quality modeling::

    # Example: River water geochemistry with carbonate equilibrium

    TITLE River carbonate system

    SOLUTION 1   River water - Initial composition
        units       mol/kgw
        temp        15.0
        pH          7.5
        pe          4.0      # Oxidizing conditions
        Ca          0.002    # Calcium
        Mg          0.001    # Magnesium
        Na          0.003    # Sodium
        K           0.0002   # Potassium
        Cl          0.002    # Chloride
        S(6)        0.001    # Sulfate as S
        C(4)        0.003    # Carbonate as C
        N(-3)       0.00001  # Ammonium as N
        N(5)        0.00005  # Nitrate as N
    END

    EQUILIBRIUM_PHASES 1
        Calcite     0.0      # Equilibrium with calcite (SI=0)
        CO2(g)     -3.5      # Atmospheric CO2
    END

    KINETICS 1
        # Example kinetic reaction: organic matter decomposition
        -formula  CH2O  -1  O2  -1  CO2  1  H2O  1
        -m0       0.001     # Initial moles
        -parms    0.01      # Rate constant (1/day)
    END

**Key PHREEQC blocks:**

* ``SOLUTION``: Defines initial water chemistry with concentrations and master species
* ``EQUILIBRIUM_PHASES``: Minerals and gases at chemical equilibrium
* ``KINETICS``: User-defined kinetic reactions with rate expressions
* ``EXCHANGE``: Cation exchange surfaces (e.g., clay minerals)
* ``SURFACE``: Surface complexation modeling (e.g., iron oxide adsorption)


Switching between NATIVE_BGC_FLEX and PHREEQC
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

OpenWQ provides two biogeochemistry engines:

**NATIVE_BGC_FLEX** (default):

* User-defined kinetic expressions using the ExprTk library
* Fast computation, suitable for simple reaction networks
* Direct control over reaction equations and parameters
* No thermodynamic database required

**PHREEQC**:

* Full thermodynamic equilibrium calculations
* Comprehensive mineral/gas equilibria
* Built-in aqueous speciation
* Requires thermodynamic database (.dat file)
* More computationally intensive

To switch between engines, modify the ``BIOGEOCHEMISTRY`` field in the :doc:`Master configuration <4_1_1Master>`:

.. code-block:: json

    // For NATIVE_BGC_FLEX:
    "BIOGEOCHEMISTRY": "NATIVE_BGC_FLEX"

    // For PHREEQC:
    "BIOGEOCHEMISTRY": "PHREEQC"

When using PHREEQC, the biogeochemistry JSON file must contain the PHREEQC-specific configuration fields (``FILEPATH``, ``DATABASE``, ``BGC_GENERAL_MOBILE_SPECIES``) instead of the ``CYCLING_FRAMEWORKS`` structure used by NATIVE_BGC_FLEX.


Discovering PHREEQC component names
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PHREEQC determines chemical components dynamically from the database and input file. To discover the available component names:

1. Set ``RUN_MODE_DEBUG: true`` in the master configuration
2. Run OpenWQ with an initial PHREEQC configuration
3. Check the console output for the component list from ``FindComponents()``

The component names used in ``BGC_GENERAL_MOBILE_SPECIES`` must match the names discovered by PHREEQC. Common components (with the standard ``phreeqc.dat`` database) include:

* ``H`` - Hydrogen
* ``O`` - Oxygen
* ``Charge`` - Electrical charge balance
* ``Ca`` - Calcium
* ``Mg`` - Magnesium
* ``Na`` - Sodium
* ``K`` - Potassium
* ``C`` - Carbon (total)
* ``N`` - Nitrogen (total)
* ``S`` - Sulfur
* ``Cl`` - Chloride

.. warning::

    The species names in ``BGC_GENERAL_MOBILE_SPECIES`` must match valid component names found by PHREEQC. If a name is not found in the component list, OpenWQ will report a warning during initialization and skip that species.


Example test case
~~~~~~~~~~~~~~~~~~

A complete PHREEQC test case is available in the ``test_case_phreeqc/`` directory of the mizuRoute-OpenWQ repository. It demonstrates:

* Configuration of a simple river network with geochemical reactions
* Initial conditions for Ca, Mg, and Na in ``mol/kgw``
* PHREEQC input file with solution chemistry definitions
* Transport using the pure advection module (``OPENWQ_NATIVE_TRANSP_DISS_ADV``)
* HDF5 output of chemical species concentrations


Best practices
~~~~~~~~~~~~~~~

1. **Start simple**: Begin with a few species and simple equilibrium reactions, then add complexity
2. **Validate with PHREEQC standalone**: Test your .pqi file in standalone PHREEQC before coupling with OpenWQ
3. **Check convergence**: Monitor PhreeqcRM warnings in the log output for convergence issues
4. **Balance charges**: Ensure your solutions are electrically balanced to avoid numerical issues
5. **Match units**: PHREEQC uses mol/kgw internally; ensure input concentrations are converted appropriately
