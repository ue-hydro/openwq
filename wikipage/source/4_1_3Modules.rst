Modules
==================================

The ``MODULES`` block in the :doc:`Master configuration <4_1_1Master>` specifies the module selections for biogeochemistry, dissolved transport, lateral exchange, sediment transport, and sorption isotherms. Each module is configured through its own JSON file, referenced from the master configuration. The sections below describe each module type and its available options.


BIOGEOCHEMISTRY
~~~~~~~~~~~~~~~~

OpenWQ provides two biogeochemistry engines. **NATIVE_BGC_FLEX** offers a flexible, user-defined kinetic reaction framework using the ExprTk library. **PHREEQC** integrates the USGS PHREEQC geochemical modeling engine for advanced equilibrium and kinetic geochemical calculations. The engine is selected by setting the ``BIOGEOCHEMISTRY`` field in the :doc:`Master configuration <4_1_1Master>`.


NATIVE_BGC_FLEX
^^^^^^^^^^^^^^^^

The NATIVE_BGC_FLEX biogeochemistry module uses a flexible reaction framework. The configuration JSON file defines the chemical species and the transformation kinetics for each cycling framework.

For the PHREEQC geochemistry engine alternative, see the `PHREEQC`_ section below.

The biogeochemical cycling file is a JSON file that defines both the chemical species to include in the model configuration and the transformations that define each cycling framework.

**Principal Key 1**: ``MODULE_NAME``

+--------------------+-------------------------------------------------------+
| ``MODULE_NAME``    | Select the module name desired                        |
+--------------------+-------------------------------------------------------+

**Principal Key 2**: ``CHEMICAL_SPECIES``

+-----------------------------------+-------------------------------------------------------+
| ``LIST`` -> ``(i#)``              | - Numbered list of chemical species                   |
|                                   | - Example: ``"1":"NO3"``                              |
+-----------------------------------+-------------------------------------------------------+
| ``BGC_GENERAL_MOBILE_SPECIES``    | List of mobile species names (strings)                |
+-----------------------------------+-------------------------------------------------------+

**Principal Key 3**: ``CYCLING_FRAMEWORK``

+---------------------------------------------------------------+---------------------------------------------------------+
| ``<BGQ_cycle_name>``                                          | Define cycling framework names                          |
+---------------------------------------------------------------+---------------------------------------------------------+
| ``<BGQ_cycle_name>`` -> ``LIST_TRANSFORMATIONS`` -> ``(i#)``  | - List transformation names for each cycling framework  |
|                                                               | - Content: ``"index"``: ``"transformation_name"``       |
|                                                               | - Format: ``(i#): (s#)``                                |
|                                                               | - Example: ``"1": "mineralization"``                    |
+---------------------------------------------------------------+---------------------------------------------------------+
| ``(i#)`` -> ``CONSUMED``                                      | Species consumed. e.g, ``"NH4"``                        |
+---------------------------------------------------------------+---------------------------------------------------------+
| ``(i#)`` -> ``PRODUCED``                                      | Species produced, e.g., ``"NO3"``                       |
+---------------------------------------------------------------+---------------------------------------------------------+
| ``(i#)`` -> ``KINETICS``                                      | - Kinetics equation and time units                      |
|                                                               | - Content: [``equation``, ``units``]                    |
|                                                               | - Format: ``[(s#), (s#)]``                              |
|                                                               | - Example: ``["NO3 * k", "1/day"]``                     |
+---------------------------------------------------------------+---------------------------------------------------------+
| ``(i#)`` -> ``PARAMETER_NAMES``                               | List of parameter names in transformation equation      |
| ``(i#)`` -> ``PARAMETER_VALUES``                              | - Define parameter values for equation                  |
|                                                               | - Content: ``parameter_name``: ``value``                |
|                                                               | - Format: ``(s#): (f#)``                                |
|                                                               | - Example: ``"k": 0.001``                               |
+---------------------------------------------------------------+---------------------------------------------------------+

The symbol ``(i#)`` refers to a integer number sequence.. The symbol ``(s#)`` refers to a string input. The symbol ``<f#>`` refers to a float input value.

The Kinetics Equation is supported by the C++ Mathematical Expression Toolkit Library
`C++ Mathematical Expression Toolkit Library <http://www.partow.net/programming/exprtk/index.html>`_. The equations can be written with (1) multiple chemical species, (2) user-defined parameters (``PARAMETER_NAMES`` and ``PARAMETER_VALUES``), and (3) in-built hydrological model dependencies.

These model dependencies are tailored to each hydrological model. The available dependency variables are listed below for each supported host model.


Available Dependency Variables by Host Model
"""""""""""""""""""""""""""""""""""""""""""""

**MIZUROUTE** (river routing model):

+-------------------+-------------------------------------------------------+---------------+
| Variable Name     | Description                                           | Units         |
+===================+=======================================================+===============+
| ``Treach_K``      | Reach water temperature (from air temperature proxy)  | Kelvin [K]    |
+-------------------+-------------------------------------------------------+---------------+
| ``SWrad_Wm2``     | Incoming shortwave radiation at the reach             | W/m²          |
+-------------------+-------------------------------------------------------+---------------+

**SUMMA** (land surface / snow / soil model):

+-------------------+-------------------------------------------------------+---------------+
| Variable Name     | Description                                           | Units         |
+===================+=======================================================+===============+
| ``SM``            | Soil moisture (volumetric water content fraction)     | [-]           |
+-------------------+-------------------------------------------------------+---------------+
| ``Tair_K``        | Air temperature                                       | Kelvin [K]    |
+-------------------+-------------------------------------------------------+---------------+
| ``Tsoil_K``       | Soil temperature                                      | Kelvin [K]    |
+-------------------+-------------------------------------------------------+---------------+
| ``SWrad_Wm2``     | Incoming shortwave radiation at the surface           | W/m²          |
+-------------------+-------------------------------------------------------+---------------+

.. note::

   Dependency variables can be used directly in kinetic expressions. For example, a temperature-dependent reaction rate could be written as:
   ``"KINETICS": ["NO3 * k * exp(-Ea / (R * Tsoil_K))", "1/day"]``
   where ``Ea`` and ``R`` are user-defined parameters.



Example:

.. code-block:: json

    {
        "CHEMICAL_SPECIES": {
            "LIST": {
                "1": "NO3",
                "2": "NH4",
                "3": "SRP",
                "4": "partP",
            },
            "MOBILE_SPECIES": ["NO3", "NH4", "partP"]
        },
        "CYCLING_FRAMEWORKS": {
            "N_inorg": {
                "LIST_TRANSFORMATIONS":{
                    "1": "nitrification",
                    "2": "denitrification"
                },
                "1":{
                    "CONSUMED": "NH4",
                    "PRODUCED": "NO3",
                    "KINETICS": ["NH4 * k", "1/day"],
                    "PARAMETER_NAMES": ["k"],
                    "PARAMETER_VALUES":{
                        "k": 0.01
                    }
                },
                "2":{
                    "CONSUMED": "NO3",
                    "PRODUCED": "N2",
                    "KINETICS": ["NO3 * k / (p^2)", "1/day"],
                    "PARAMETER_NAMES": ["k","p"],
                    "PARAMETER_VALUES":{
                        "k": 0.01,
                        "p": 10
                    }
                }
            },
            "P_inorg": {
                "LIST_TRANSFORMATIONS":{
                    "1": "dynamic_equilibrium"
                },
                "1":{
                    "CONSUMED": "SRP",
                    "PRODUCED": "partP",
                    "KINETICS": ["SRP * k * Tsoil_K / 273.15", "1/day"],
                    "PARAMETER_NAMES": ["k"],
                    "PARAMETER_VALUES":{
                        "k": 0.01
                    }
                }
            }
        }
    }


PHREEQC
^^^^^^^^^^^^^^^^

OpenWQ integrates the `PHREEQC <https://www.usgs.gov/software/phreeqc-version-3>`_ geochemical modeling engine through the PhreeqcRM library.
This enables advanced equilibrium and kinetic geochemical calculations as an alternative to the flexible BGC reaction networks.

PHREEQC is selected by setting the ``BIOGEOCHEMISTRY`` module to ``PHREEQC`` in the :doc:`Master configuration <4_1_1Master>`.


Overview
"""""""""

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
""""""""""""""

The PHREEQC module is configured through a JSON file referenced in the master configuration:

.. code-block:: json

    {
        "MODULE_NAME": "PHREEQC",
        "FILEPATH": "openwq_in/phreeqc_river.pqi",
        "DATABASE": "openwq_in/phreeqc.dat",
        "COMPONENT_H2O": true,
        "BGC_GENERAL_MOBILE_SPECIES": ["H", "O", "Charge", "Ca", "Mg", "Na"],
        "TEMPERATURE": {
            "RIVER_NETWORK_REACHES": "Treach_K"
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


PHREEQC input file (.pqi)
""""""""""""""""""""""""""

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
"""""""""""""""""""""""

The thermodynamic database contains equilibrium constants and thermodynamic data for aqueous species, minerals, gases, and surface species. Several databases are available from USGS:

* ``phreeqc.dat`` -- General-purpose database (recommended starting point)
* ``wateq4f.dat`` -- Extended database with additional minerals and trace metals
* ``minteq.v4.dat`` -- MINTEQA2-compatible database
* ``llnl.dat`` -- Lawrence Livermore National Laboratory database (comprehensive)


Mobile species
"""""""""""""""

The ``BGC_GENERAL_MOBILE_SPECIES`` array specifies which PHREEQC components are subject to advective transport between compartment cells. Only mobile species are transported by the OpenWQ transport module; immobile species remain in their original grid cell.

Species are identified by their component name (string), which must match the names returned by PHREEQC's ``FindComponents()`` function. The component list can be inspected in the OpenWQ debug output when ``RUN_MODE_DEBUG`` is enabled in the master configuration.


Temperature coupling
"""""""""""""""""""""

PHREEQC reaction rates and equilibrium constants are temperature-dependent. The ``TEMPERATURE`` field maps each host-model compartment to a temperature data source variable name.

Temperature values are passed to PhreeqcRM at each time step for accurate geochemical calculations.
The variable name must match a variable exported by the host hydrological model.

.. note::

   **Automatic temperature unit conversion:** OpenWQ auto-detects if temperatures are in Kelvin (values > 100) and automatically converts to Celsius for PHREEQC. This allows direct use of host model temperature variables without manual conversion.

**Available temperature variables by host model:**

+----------------+---------------------------+------------------------------------------+
| Host Model     | Variable Name             | Description                              |
+================+===========================+==========================================+
| MIZUROUTE      | ``Treach_K``              | Reach water temperature (Kelvin)         |
+----------------+---------------------------+------------------------------------------+
| MIZUROUTE      | ``Tlake_K``               | Lake temperature (Kelvin)                |
+----------------+---------------------------+------------------------------------------+
| SUMMA          | ``Tair_K``                | Air temperature (Kelvin)                 |
+----------------+---------------------------+------------------------------------------+
| SUMMA          | ``Tsoil_K``               | Soil temperature (Kelvin)                |
+----------------+---------------------------+------------------------------------------+

Example configuration:

.. code-block:: json

    {
        "TEMPERATURE": {
            "RIVER_NETWORK_REACHES": "Treach_K",
            "LAKES": "Tlake_K"
        }
    }


Initialization process
"""""""""""""""""""""""

During model initialization, the PHREEQC module performs the following steps:

1. **Grid calculation**: Computes the total number of grid cells across all compartments (``nxyz = nx * ny * nz``)
2. **PhreeqcRM creation**: Creates a PhreeqcRM instance with the computed grid size and the configured number of threads
3. **Database loading**: Loads the thermodynamic database file
4. **Input file execution**: Runs the PHREEQC input file to set up initial solution chemistry
5. **Component discovery**: Identifies all chemical components using ``FindComponents()``
6. **Verification**: Confirms that the number of discovered components matches the OpenWQ configuration


Example PHREEQC input file
""""""""""""""""""""""""""""

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
"""""""""""""""""""""""""""""""""""""""""""""""

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
"""""""""""""""""""""""""""""""""""""

PHREEQC determines chemical components dynamically from the database and input file. To discover the available component names:

1. Set ``RUN_MODE_DEBUG: true`` in the master configuration
2. Run OpenWQ with an initial PHREEQC configuration
3. Check the console output for the component list from ``FindComponents()``

The component names used in ``BGC_GENERAL_MOBILE_SPECIES`` must match the names discovered by PHREEQC.

.. note::

   Species names are **case-insensitive** in OpenWQ. Both ``"Ca"`` and ``"CA"`` will match the PHREEQC component ``Ca``.


Complete PHREEQC Species Reference (phreeqc.dat)
""""""""""""""""""""""""""""""""""""""""""""""""""

The following tables list all chemical species supported by the standard ``phreeqc.dat`` thermodynamic database.

**Master Species (Elements)**

These are the primary elements that PHREEQC tracks as components. Use these element names in ``BGC_GENERAL_MOBILE_SPECIES``.

+----------------+------------------------------+--------------------+------------------+
| Element        | Description                  | Default Species    | GFW (g/mol)      |
+================+==============================+====================+==================+
| ``H``          | Hydrogen                     | H+                 | 1.008            |
+----------------+------------------------------+--------------------+------------------+
| ``O``          | Oxygen                       | H2O                | 16.0             |
+----------------+------------------------------+--------------------+------------------+
| ``Ca``         | Calcium                      | Ca+2               | 40.08            |
+----------------+------------------------------+--------------------+------------------+
| ``Mg``         | Magnesium                    | Mg+2               | 24.312           |
+----------------+------------------------------+--------------------+------------------+
| ``Na``         | Sodium                       | Na+                | 22.9898          |
+----------------+------------------------------+--------------------+------------------+
| ``K``          | Potassium                    | K+                 | 39.102           |
+----------------+------------------------------+--------------------+------------------+
| ``Fe``         | Iron                         | Fe+2               | 55.847           |
+----------------+------------------------------+--------------------+------------------+
| ``Mn``         | Manganese                    | Mn+2               | 54.938           |
+----------------+------------------------------+--------------------+------------------+
| ``Al``         | Aluminum                     | Al+3               | 26.9815          |
+----------------+------------------------------+--------------------+------------------+
| ``Ba``         | Barium                       | Ba+2               | 137.34           |
+----------------+------------------------------+--------------------+------------------+
| ``Sr``         | Strontium                    | Sr+2               | 87.62            |
+----------------+------------------------------+--------------------+------------------+
| ``Si``         | Silicon (silica)             | H4SiO4             | 28.0843          |
+----------------+------------------------------+--------------------+------------------+
| ``Cl``         | Chloride                     | Cl-                | 35.453           |
+----------------+------------------------------+--------------------+------------------+
| ``C``          | Carbon (inorganic)           | CO3-2              | 12.0111          |
+----------------+------------------------------+--------------------+------------------+
| ``S``          | Sulfur                       | SO4-2              | 32.064           |
+----------------+------------------------------+--------------------+------------------+
| ``N``          | Nitrogen                     | NO3-               | 14.0067          |
+----------------+------------------------------+--------------------+------------------+
| ``B``          | Boron                        | H3BO3              | 10.81            |
+----------------+------------------------------+--------------------+------------------+
| ``P``          | Phosphorus                   | PO4-3              | 30.9738          |
+----------------+------------------------------+--------------------+------------------+
| ``F``          | Fluoride                     | F-                 | 18.9984          |
+----------------+------------------------------+--------------------+------------------+
| ``Li``         | Lithium                      | Li+                | 6.939            |
+----------------+------------------------------+--------------------+------------------+
| ``Br``         | Bromide                      | Br-                | 79.904           |
+----------------+------------------------------+--------------------+------------------+
| ``Zn``         | Zinc                         | Zn+2               | 65.37            |
+----------------+------------------------------+--------------------+------------------+
| ``Cd``         | Cadmium                      | Cd+2               | 112.4            |
+----------------+------------------------------+--------------------+------------------+
| ``Pb``         | Lead                         | Pb+2               | 207.19           |
+----------------+------------------------------+--------------------+------------------+
| ``Cu``         | Copper                       | Cu+2               | 63.546           |
+----------------+------------------------------+--------------------+------------------+
| ``Charge``     | Electrical charge balance    | --                 | --               |
+----------------+------------------------------+--------------------+------------------+
| ``Alkalinity`` | Alkalinity (as CaCO3)        | CO3-2              | 50.05            |
+----------------+------------------------------+--------------------+------------------+

**Redox States**

For redox-sensitive elements, PHREEQC tracks multiple oxidation states:

+----------------+--------------------------------------------------+
| Element        | Available Redox States                           |
+================+==================================================+
| ``H``          | H(0) = H2, H(1) = H+                             |
+----------------+--------------------------------------------------+
| ``O``          | O(0) = O2, O(-2) = H2O                           |
+----------------+--------------------------------------------------+
| ``Fe``         | Fe(+2) = Fe+2, Fe(+3) = Fe+3                     |
+----------------+--------------------------------------------------+
| ``Mn``         | Mn(+2) = Mn+2, Mn(+3) = Mn+3                     |
+----------------+--------------------------------------------------+
| ``Cu``         | Cu(+1) = Cu+, Cu(+2) = Cu+2                      |
+----------------+--------------------------------------------------+
| ``C``          | C(+4) = CO3-2/HCO3-, C(-4) = CH4                 |
+----------------+--------------------------------------------------+
| ``S``          | S(6) = SO4-2, S(-2) = HS-/H2S                    |
+----------------+--------------------------------------------------+
| ``N``          | N(+5) = NO3-, N(+3) = NO2-, N(0) = N2, N(-3) = NH4+ |
+----------------+--------------------------------------------------+

**Mineral Phases (for EQUILIBRIUM_PHASES)**

Carbonate minerals:

* Calcite, Aragonite, Dolomite, Siderite, Rhodochrosite
* Strontianite, Witherite, Cerussite, Smithsonite, Otavite

Sulfate minerals:

* Gypsum, Anhydrite, Celestite, Barite, Anglesite
* Arcanite, Mirabilite, Thenardite, Epsomite, Melanterite

Silicate minerals:

* SiO2(a), Chalcedony, Quartz, Gibbsite, Kaolinite
* Albite, Anorthite, K-feldspar, K-mica, Illite
* Ca-Montmorillonite, Chlorite(14A), Talc, Chrysotile, Sepiolite

Iron/Manganese minerals:

* Hematite, Goethite, Fe(OH)3(a), Pyrite, FeS(ppt), Mackinawite
* Hausmannite, Manganite, Pyrochroite

Other minerals:

* Hydroxyapatite, Vivianite, Fluorite, Halite, Sylvite
* Sulfur, Sphalerite, Willemite, Alunite, Jarosite-K

**Gas Phases**

Available for equilibrium with aqueous solutions:

* CO2(g), O2(g), H2(g), N2(g), H2S(g), CH4(g), NH3(g), H2O(g)

Redox-uncoupled gas species (for independent gas equilibria):

* Hdg (H2), Oxg (O2), Mtg (CH4), Sg (H2S), Ntg (N2)


**Recommended Species for Common Applications**

+---------------------------+--------------------------------------------------+
| Application               | Recommended Mobile Species                       |
+===========================+==================================================+
| Major ion chemistry       | Ca, Mg, Na, K, Cl, S, C                          |
+---------------------------+--------------------------------------------------+
| Nutrient modeling         | N, P, C                                          |
+---------------------------+--------------------------------------------------+
| Trace metal transport     | Fe, Mn, Zn, Cu, Pb, Cd                           |
+---------------------------+--------------------------------------------------+
| Complete water quality    | Ca, Mg, Na, K, Cl, S, C, N, P, Si, Fe, Mn        |
+---------------------------+--------------------------------------------------+

.. warning::

    The species names in ``BGC_GENERAL_MOBILE_SPECIES`` must match valid component names found by PHREEQC. If a name is not found in the component list, OpenWQ will report a warning during initialization and skip that species.


Example test case
""""""""""""""""""

A complete PHREEQC test case is available in the ``test_case_phreeqc/`` directory of the mizuRoute-OpenWQ repository. It demonstrates:

* Configuration of a simple river network with geochemical reactions
* Initial conditions for Ca, Mg, and Na in ``mol/kgw``
* PHREEQC input file with solution chemistry definitions
* Transport using the pure advection module (``OPENWQ_NATIVE_TRANSP_DISS_ADV``)
* HDF5 output of chemical species concentrations


Best practices
"""""""""""""""

1. **Start simple**: Begin with a few species and simple equilibrium reactions, then add complexity
2. **Validate with PHREEQC standalone**: Test your .pqi file in standalone PHREEQC before coupling with OpenWQ
3. **Check convergence**: Monitor PhreeqcRM warnings in the log output for convergence issues
4. **Balance charges**: Ensure your solutions are electrically balanced to avoid numerical issues
5. **Match units**: PHREEQC uses mol/kgw internally; ensure input concentrations are converted appropriately


.. _phreeqc-process-overlap:

Process overlap with other OpenWQ modules
""""""""""""""""""""""""""""""""""""""""""

PHREEQC is a comprehensive geochemistry engine that natively handles several processes also available as standalone OpenWQ modules. Understanding these overlaps is critical to avoid double-counting.

**Sorption (PHREEQC vs. SORPTION_ISOTHERM module)**

.. warning::

   When PHREEQC is the active biogeochemistry engine, the ``SORPTION_ISOTHERM`` module **must** be set to ``NONE``. OpenWQ will abort with an error if both are active simultaneously.

PHREEQC natively handles sorption processes through:

* **Surface complexation** (``SURFACE`` keyword) — Dzombak & Morel diffuse double-layer model, CD-MUSIC
* **Ion exchange** (``EXCHANGE`` keyword) — Gaines-Thomas, Gapon, Vanselow conventions
* **Langmuir isotherms** — directly expressible as mass-action reactions with ``SURFACE_SPECIES``
* **Freundlich isotherms** — achievable via workaround approaches with ``-no_check`` / ``-mole_balance``

These are mechanistically more rigorous than the simplified Freundlich/Langmuir isotherms in the ``SORPTION_ISOTHERM`` module. To configure sorption when using PHREEQC, define ``SURFACE`` or ``EXCHANGE`` blocks in the ``.pqi`` input file rather than enabling the SI module.

**Transport (PHREEQC vs. TRANSPORT_DISSOLVED module)**

PHREEQC also supports 1D advection-dispersion via the ``TRANSPORT`` and ``ADVECTION`` keywords. However, the OpenWQ integration uses PHREEQC strictly as a **batch reaction engine** (cell-by-cell equilibration via ``RunCells()``). Spatial transport is handled by the OpenWQ ``TRANSPORT_DISSOLVED`` module, which operates on a separate derivative array.

.. note::

   Do **not** include ``TRANSPORT`` or ``ADVECTION`` keywords in the ``.pqi`` input file when using PHREEQC within OpenWQ. These keywords are processed during the ``RunFile()`` initialization call and may alter the initial solution state unpredictably. All spatial transport should be handled by the ``TRANSPORT_DISSOLVED`` module.

**Kinetic reactions (PHREEQC vs. NATIVE_BGC_FLEX)**

PHREEQC and NATIVE_BGC_FLEX are mutually exclusive — only one biogeochemistry engine runs per simulation. This is enforced in the code and poses no risk of double-counting.


TRANSPORT_DISSOLVED
~~~~~~~~~~~~~~~~~~~~

The Transport Dissolved file is a JSON file that configures how dissolved chemical species are transported with the water flux between model cells.

Module Options
"""""""""""""""

+----------------------------------+----------------------------------------------------------+
| ``MODULE_NAME``                  | Description                                              |
+==================================+==========================================================+
| ``OPENWQ_NATIVE_TD_ADVDISP``    | Advection + Dispersion transport                         |
+----------------------------------+----------------------------------------------------------+
| ``NATIVE_TD_ADV``               | Advection-only transport                                 |
+----------------------------------+----------------------------------------------------------+
| ``NONE``                         | No dissolved transport                                   |
+----------------------------------+----------------------------------------------------------+


Transport Configuration
""""""""""""""""""""""""

The ``TRANSPORT_CONFIGURATION`` key contains the dispersion parameters. These parameters are used only when the module is ``OPENWQ_NATIVE_TD_ADVDISP``.

+-----------------------------------+--------------------------------------------------------------+
| Key                               | Description                                                  |
+===================================+==============================================================+
| ``dispersion_x_m2/s``            | Dispersion coefficient in x direction [m2/s]                 |
+-----------------------------------+--------------------------------------------------------------+
| ``dispersion_y_m2/s``            | Dispersion coefficient in y direction [m2/s]                 |
+-----------------------------------+--------------------------------------------------------------+
| ``dispersion_z_m2/s``            | Dispersion coefficient in z direction [m2/s]                 |
+-----------------------------------+--------------------------------------------------------------+
| ``characteristic_length_m``      | Characteristic distance between cell centers [m]             |
+-----------------------------------+--------------------------------------------------------------+

The effective dispersion rate is computed as:

.. math::

    D_{eff} = \frac{D_{avg}}{L^2}

where :math:`D_{avg} = (D_x + D_y + D_z) / 3` and :math:`L` is the ``characteristic_length_m``.

.. note::

    The dispersion uses a Fickian approximation between cell pairs. The dispersive flux is proportional to the concentration difference between source and recipient cells, and is bidirectional (it can smooth concentration gradients in either direction regardless of flow).


Example (Advection + Dispersion)
"""""""""""""""""""""""""""""""""

.. code-block:: json

    {
        "MODULE_NAME": "OPENWQ_NATIVE_TD_ADVDISP",

        "TRANSPORT_CONFIGURATION": {
            "dispersion_x_m2/s": 0.5,
            "dispersion_y_m2/s": 0.5,
            "dispersion_z_m2/s": 0.1,
            "characteristic_length_m": 100.0
        }
    }


Example (Advection-only)
""""""""""""""""""""""""""

.. code-block:: json

    {
        "MODULE_NAME": "NATIVE_TD_ADV",

        "TRANSPORT_CONFIGURATION": {
            "dispersion_x_m2/s": 0.0,
            "dispersion_y_m2/s": 0.0,
            "dispersion_z_m2/s": 0.0,
            "characteristic_length_m": 100.0
        }
    }


LATERAL_EXCHANGE
~~~~~~~~~~~~~~~~~~~~

The Lateral Exchange module controls boundary mixing between adjacent compartments.
It simulates the diffusive exchange of dissolved chemical species across compartment interfaces (e.g., between surface water and soil layers).


Module Options
"""""""""""""""

+----------------------------------+----------------------------------------------------------+
| ``MODULE_NAME``                  | Description                                              |
+==================================+==========================================================+
| ``NATIVE_LE_BOUNDMIX``          | Boundary mixing between upper and lower compartments     |
+----------------------------------+----------------------------------------------------------+
| ``NONE``                         | No lateral exchange                                      |
+----------------------------------+----------------------------------------------------------+


Configuration
""""""""""""""

The ``CONFIGURATION`` key contains numbered entries, each defining an exchange interface between two compartments.

+-------------------------------+--------------------------------------------------------------+
| Key                           | Description                                                  |
+===============================+==============================================================+
| ``direction``                | Exchange direction: ``"x"``, ``"y"``, or ``"z"``             |
+-------------------------------+--------------------------------------------------------------+
| ``upper_compartment``        | Name of the upper compartment                                |
+-------------------------------+--------------------------------------------------------------+
| ``lower_compartment``        | Name of the lower compartment                                |
+-------------------------------+--------------------------------------------------------------+
| ``K_val``                    | Exchange coefficient [1/s]                                   |
+-------------------------------+--------------------------------------------------------------+

Multiple exchange interfaces can be defined by adding numbered entries (``"1"``, ``"2"``, etc.) under ``CONFIGURATION``.


Example
""""""""

.. code-block:: json

    {
        "MODULE_NAME": "NATIVE_LE_BOUNDMIX",
        "CONFIGURATION": {
            "1": {
                "direction": "z",
                "upper_compartment": "RUNOFF",
                "lower_compartment": "ILAYERVOLFRACWAT_SOIL",
                "K_val": 0.000000001
            }
        }
    }


TRANSPORT_SEDIMENTS
~~~~~~~~~~~~~~~~~~~~

The Transport Sediments module simulates soil erosion and sediment transport processes.
Two erosion models are available, both based on the HYPE (HYdrological Predictions for the Environment) framework.


Module Options
"""""""""""""""

+----------------------------------+----------------------------------------------------------+
| ``MODULE_NAME``                  | Description                                              |
+==================================+==========================================================+
| ``HYPE_MMF``                    | Morgan-Morgan-Finney erosion model                       |
+----------------------------------+----------------------------------------------------------+
| ``HYPE_HBVSED``                 | HBV-SED erosion model                                    |
+----------------------------------+----------------------------------------------------------+
| ``NONE``                         | No sediment transport                                    |
+----------------------------------+----------------------------------------------------------+

In the :doc:`Master configuration <4_1_1Master>`, the ``TRANSPORT_SEDIMENTS`` block also requires a ``SEDIMENT_COMPARTMENT`` key specifying which compartment the erosion process acts on.


Common Configuration
"""""""""""""""""""""

Both models share a common JSON structure with the following blocks:

+-------------------------------+--------------------------------------------------------------+
| Key                           | Description                                                  |
+===============================+==============================================================+
| ``MODULE_NAME``              | ``"HYPE_MMF"`` or ``"HYPE_HBVSED"``                         |
+-------------------------------+--------------------------------------------------------------+
| ``CONFIGURATION``            | General settings (direction, inhibit compartment, format)    |
+-------------------------------+--------------------------------------------------------------+
| ``PARAMETER_DEFAULTS``       | Default values for erosion parameters                        |
+-------------------------------+--------------------------------------------------------------+
| ``PARAMETERS``               | Spatially distributed parameter overrides per cell           |
+-------------------------------+--------------------------------------------------------------+

CONFIGURATION
..............

+--------------------------------------+--------------------------------------------------------------+
| Key                                  | Description                                                  |
+======================================+==============================================================+
| ``direction``                       | Erosion direction: ``"x"``, ``"y"``, or ``"z"``              |
+--------------------------------------+--------------------------------------------------------------+
| ``EROSION_INHIBIT_COMPARTMENT``     | Compartment that inhibits erosion (e.g., snow cover)         |
+--------------------------------------+--------------------------------------------------------------+
| ``Data_Format``                     | ``"JSON"`` or ``"ASCII"`` for parameter input                |
+--------------------------------------+--------------------------------------------------------------+


HYPE_MMF
^^^^^^^^^

The Morgan-Morgan-Finney (MMF) model computes soil detachment and sediment transport using soil erodibility, rainfall kinetic energy, and ground cover.

+------------------------------+--------------------------------------------------------------+
| Parameter                    | Description                                                  |
+==============================+==============================================================+
| ``COHESION``                | Soil cohesion                                                |
+------------------------------+--------------------------------------------------------------+
| ``ERODIBILITY``             | Soil erodibility factor                                      |
+------------------------------+--------------------------------------------------------------+
| ``SREROEXP``                | Surface runoff erosion exponent                              |
+------------------------------+--------------------------------------------------------------+
| ``CROPCOVER``               | Crop cover factor [0-1]                                      |
+------------------------------+--------------------------------------------------------------+
| ``GROUNDCOVER``             | Ground cover factor [0-1]                                    |
+------------------------------+--------------------------------------------------------------+
| ``SLOPE``                   | Terrain slope                                                |
+------------------------------+--------------------------------------------------------------+
| ``TRANSPORT_FACTOR_1``      | Sediment transport capacity factor 1                         |
+------------------------------+--------------------------------------------------------------+
| ``TRANSPORT_FACTOR_2``      | Sediment transport capacity factor 2                         |
+------------------------------+--------------------------------------------------------------+


HYPE_HBVSED
^^^^^^^^^^^^^

The HBV-SED model uses a simpler erosion formulation based on precipitation intensity, soil erodibility, and slope characteristics.

+---------------------------------------------+--------------------------------------------------------------+
| Parameter                                   | Description                                                  |
+=============================================+==============================================================+
| ``SLOPE``                                  | Terrain slope                                                |
+---------------------------------------------+--------------------------------------------------------------+
| ``EROSION_INDEX``                          | Erosion index                                                |
+---------------------------------------------+--------------------------------------------------------------+
| ``SOIL_EROSION_FACTOR_LAND_DEPENDENCE``    | Soil erosion factor (land use dependence)                    |
+---------------------------------------------+--------------------------------------------------------------+
| ``SOIL_EROSION_FACTOR_SOIL_DEPENDENCE``    | Soil erosion factor (soil type dependence)                   |
+---------------------------------------------+--------------------------------------------------------------+
| ``SLOPE_EROSION_FACTOR_EXPONENT``          | Slope erosion factor exponent                                |
+---------------------------------------------+--------------------------------------------------------------+
| ``PRECIP_EROSION_FACTOR_EXPONENT``         | Precipitation erosion factor exponent                        |
+---------------------------------------------+--------------------------------------------------------------+
| ``PARAM_SCALING_EROSION_INDEX``            | Scaling parameter for erosion index                          |
+---------------------------------------------+--------------------------------------------------------------+

The HYPE_HBVSED model also accepts a ``MONTHLY_EROSION_FACTOR`` array (12 values, January--December) that adjusts the erosion rate seasonally. The factor is applied as ``erodmonth = 1.0 + MONTHLY_EROSION_FACTOR[month]``.


Spatially Distributed Parameters
""""""""""""""""""""""""""""""""""

Under the ``PARAMETERS`` key, each parameter can be overridden for specific cells.
Default values from ``PARAMETER_DEFAULTS`` are used for any cells not listed.

The format for each parameter block is:

.. code-block:: json

    "PARAMETER_NAME": {
        "0": ["IX", "IY", "IZ", "VALUE"],
        "1": [13, 1, 1, 0.094],
        "2": [19, 1, 1, 0.418]
    }

Row ``"0"`` is the header, and subsequent numbered rows provide cell-specific values.


Example (HYPE_HBVSED)
"""""""""""""""""""""""

.. code-block:: json

    {
        "MODULE_NAME": "HYPE_HBVSED",
        "CONFIGURATION": {
            "direction": "z",
            "EROSION_INHIBIT_COMPARTMENT": "ILAYERVOLFRACWAT_SNOW",
            "DATA_FORMAT": "JSON"
        },
        "PARAMETER_DEFAULTS": {
            "SLOPE": 0.4,
            "EROSION_INDEX": 0.4,
            "SOIL_EROSION_FACTOR_LAND_DEPENDENCE": 0.4,
            "SOIL_EROSION_FACTOR_SOIL_DEPENDENCE": 0.4,
            "SLOPE_EROSION_FACTOR_EXPONENT": 1.5,
            "PRECIP_EROSION_FACTOR_EXPONENT": 1.5,
            "PARAM_SCALING_EROSION_INDEX": 0.5
        },
        "MONTHLY_EROSION_FACTOR": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        "PARAMETERS": {
            "SLOPE": {
                "0": ["IX", "IY", "IZ", "VALUE"],
                "1": [13, 1, 1, 0.094],
                "2": [19, 1, 1, 0.418]
            }
        }
    }


Example (HYPE_MMF)
"""""""""""""""""""

.. code-block:: json

    {
        "MODULE_NAME": "HYPE_MMF",
        "CONFIGURATION": {
            "direction": "z",
            "EROSION_INHIBIT_COMPARTMENT": "ILAYERVOLFRACWAT_SNOW",
            "DATA_FORMAT": "JSON"
        },
        "PARAMETER_DEFAULTS": {
            "COHESION": 7.5,
            "ERODIBILITY": 2.0,
            "SREROEXP": 1.2,
            "CROPCOVER": 0.5,
            "GROUNDCOVER": 0.5,
            "SLOPE": 0.4,
            "TRANSPORT_FACTOR_1": 0.2,
            "TRANSPORT_FACTOR_2": 0.5
        },
        "PARAMETERS": {
            "COHESION": {
                "0": ["IX", "IY", "IZ", "VALUE"],
                "1": ["ALL", 1, 1, 7.5]
            }
        }
    }


SORPTION_ISOTHERM
~~~~~~~~~~~~~~~~~~~~

The Sorption Isotherm module configures equilibrium partitioning between dissolved and sorbed phases. It applies a kinetic approach: the equilibrium sorbed concentration is computed from the isotherm equation, and mass is transferred at a user-defined kinetic rate.

.. warning::

   This module is designed for use with **NATIVE_BGC_FLEX** as the biogeochemistry engine. When using **PHREEQC**, set ``SORPTION_ISOTHERM`` to ``NONE`` and configure sorption via PHREEQC's ``SURFACE`` or ``EXCHANGE`` blocks in the ``.pqi`` file instead. See :ref:`phreeqc-process-overlap` for details.

Module Options
"""""""""""""""

+----------------------------------+----------------------------------------------------------+
| ``MODULE_NAME``                  | Description                                              |
+==================================+==========================================================+
| ``FREUNDLICH``                  | Freundlich nonlinear sorption isotherm                   |
+----------------------------------+----------------------------------------------------------+
| ``LANGMUIR``                    | Langmuir sorption isotherm (finite binding sites)        |
+----------------------------------+----------------------------------------------------------+
| ``NONE``                         | No sorption                                              |
+----------------------------------+----------------------------------------------------------+


Soil Properties
""""""""""""""""

Both isotherm models require soil/medium properties, specified under the ``SOIL_PROPERTIES`` key:

+-----------------------------------+--------------------------------------------------------------+
| Key                               | Description                                                  |
+===================================+==============================================================+
| ``bulk_density_kg/m3``           | Bulk density of the soil/medium [kg/m3]                      |
+-----------------------------------+--------------------------------------------------------------+
| ``layer_thickness_m``            | Representative layer thickness [m]                           |
+-----------------------------------+--------------------------------------------------------------+


Per-Species Parameters
"""""""""""""""""""""""

Species are listed under the ``SPECIES`` key. Each species name must match a chemical species defined in the biogeochemistry module. The isotherm is applied independently to each listed species, across all compartments and cells.


FREUNDLICH
^^^^^^^^^^^

The Freundlich isotherm models nonlinear sorption:

.. math::

    q = K_{fr} \cdot C^{N_{fr}}

where :math:`q` is the sorbed concentration [mg/kg], :math:`C` is the dissolved concentration [mg/L], :math:`K_{fr}` is the Freundlich coefficient, and :math:`N_{fr}` is the Freundlich exponent.

+-----------------------------------+--------------------------------------------------------------+
| Key                               | Description                                                  |
+===================================+==============================================================+
| ``Kfr``                          | Freundlich coefficient [mg/kg / (mg/L)^Nfr]                 |
+-----------------------------------+--------------------------------------------------------------+
| ``Nfr``                          | Freundlich exponent [-] (1.0 = linear)                       |
+-----------------------------------+--------------------------------------------------------------+
| ``Kadsdes_1/s``                  | Kinetic adsorption/desorption rate [1/s]                     |
+-----------------------------------+--------------------------------------------------------------+

The equilibrium dissolved concentration is found via Newton-Raphson iteration on the mass conservation equation:

.. math::

    C_{eq} \cdot V + K_{fr} \cdot C_{eq}^{N_{fr}} \cdot \rho \cdot L = m_{total}


LANGMUIR
^^^^^^^^^

The Langmuir isotherm models sorption with a finite number of binding sites:

.. math::

    q = \frac{q_{max} \cdot K_L \cdot C}{1 + K_L \cdot C}

where :math:`q_{max}` is the maximum adsorption capacity [mg/kg] and :math:`K_L` is the Langmuir equilibrium constant [L/mg].

+-----------------------------------+--------------------------------------------------------------+
| Key                               | Description                                                  |
+===================================+==============================================================+
| ``qmax_mg/kg``                   | Maximum adsorption capacity [mg/kg_soil]                     |
+-----------------------------------+--------------------------------------------------------------+
| ``KL_L/mg``                      | Langmuir equilibrium constant [L/mg]                         |
+-----------------------------------+--------------------------------------------------------------+
| ``Kadsdes_1/s``                  | Kinetic adsorption/desorption rate [1/s]                     |
+-----------------------------------+--------------------------------------------------------------+

The equilibrium dissolved concentration is found analytically via the quadratic formula from mass conservation.


Kinetic Approach
"""""""""""""""""

Both models use a kinetic adsorption/desorption approach. The mass flux from dissolved to sorbed phase is:

.. math::

    \Delta q = (q_{eq} - q_{current}) \cdot (1 - e^{-K_{adsdes} \cdot \Delta t})

.. math::

    F_{sorption} = \Delta q \cdot \rho \cdot L

where :math:`q_{eq}` is the equilibrium sorbed concentration from the isotherm, :math:`q_{current}` is the current sorbed concentration, :math:`K_{adsdes}` is the kinetic rate, and :math:`\Delta t` is the timestep. The flux is applied as a sink on the dissolved mass in ``d_chemass_dt_chem``.


Example (Freundlich)
"""""""""""""""""""""

.. code-block:: json

    {
        "MODULE_NAME": "FREUNDLICH",
        "SOIL_PROPERTIES": {
            "bulk_density_kg/m3": 1500.0,
            "layer_thickness_m": 1.0
        },
        "SPECIES": {
            "NH4-N": {
                "Kfr": 1.2,
                "Nfr": 0.8,
                "Kadsdes_1/s": 0.001
            }
        }
    }


Example (Langmuir)
"""""""""""""""""""

.. code-block:: json

    {
        "MODULE_NAME": "LANGMUIR",
        "SOIL_PROPERTIES": {
            "bulk_density_kg/m3": 1500.0,
            "layer_thickness_m": 1.0
        },
        "SPECIES": {
            "NH4-N": {
                "qmax_mg/kg": 200.0,
                "KL_L/mg": 0.05,
                "Kadsdes_1/s": 0.002
            }
        }
    }


Example (No sorption)
""""""""""""""""""""""

.. code-block:: json

    {
        "MODULE_NAME": "NONE"
    }


|

The JSON file supports C/C++ syntax for comments: single-line comment (``//``) or comment blocks (``/*`` and ``*/``).
