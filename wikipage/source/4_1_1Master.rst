Master Configuration
==================================

The master configuration is a JSON file that provides OpenWQ with the information and instructions to run a simulation.
It references all other input files and controls which modules are enabled.
The major blocks are described below.


General Project Information
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Top-level keys that identify the simulation.

+-----------------------------+--------------------------------------+
| ``PROJECT_NAME``            | Project name string                  |
+-----------------------------+--------------------------------------+
| ``GEOGRAPHICAL_LOCATION``   | Project location string              |
+-----------------------------+--------------------------------------+
| ``AUTHORS``                 | Model run authors                    |
+-----------------------------+--------------------------------------+
| ``DATE``                    | Model run date                       |
+-----------------------------+--------------------------------------+
| ``COMMENT``                 | Additional documentation information |
+-----------------------------+--------------------------------------+


``COMPUTATIONAL_SETTINGS``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

+----------------------+-------------------------------------------------------------------------------------+
| ``RUN_MODE_DEBUG``   | Run in debug mode with extra console output (true/false)                            |
+----------------------+-------------------------------------------------------------------------------------+
| ``USE_NUM_THREADS``  | ``"all"`` or integer number of threads to use (e.g., ``4``)                         |
+----------------------+-------------------------------------------------------------------------------------+


``SOLVER``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Top-level key that selects the numerical solver.

+-------------------+-----------------------------------------------------------------------------------+
| ``"FORWARD_EULER"`` | Forward Euler explicit solver (first-order)                                     |
+-------------------+-----------------------------------------------------------------------------------+
| ``"SUNDIALS"``      | SUNDIALS CVode adaptive multi-step solver                                       |
+-------------------+-----------------------------------------------------------------------------------+

See :doc:`Numerical Solvers <4_1_9Solvers>` for details.


``OPENWQ_INPUT``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

References the main configuration file, external water fluxes, and sink/source forcing files.

CONFIG_FILEPATH
^^^^^^^^^^^^^^^^

+----------------------+-------------------------------------------------------------+
| ``CONFIG_FILEPATH``  | Path to :doc:`Configuration file <4_1_2Config>`             |
+----------------------+-------------------------------------------------------------+

EXTERNAL_WATER_FLUXES
^^^^^^^^^^^^^^^^^^^^^^

Numbered pairs of ``LABEL`` and ``FILEPATH`` for external water flux input files.
Each entry references a :doc:`Source/Sink & External Fluxes file <4_1_8SS>`. *Optional*.

+--------------------------+----------------------------------------------+
| ``(#)`` -> ``LABEL``     | Label identifying this flux entry            |
+--------------------------+----------------------------------------------+
| ``(#)`` -> ``FILEPATH``  | Path to the external water fluxes JSON file  |
+--------------------------+----------------------------------------------+

SINK_SOURCE
^^^^^^^^^^^^^

Numbered pairs of ``LABEL`` and ``FILEPATH`` for chemical source/sink input files.
Each entry references a :doc:`Source/Sink file <4_1_8SS>`. *Optional*.

+--------------------------+----------------------------------------------+
| ``(#)`` -> ``LABEL``     | Label identifying this source/sink entry     |
+--------------------------+----------------------------------------------+
| ``(#)`` -> ``FILEPATH``  | Path to the source/sink JSON file            |
+--------------------------+----------------------------------------------+


``MODULES``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Selects which process modules to enable and points to their configuration files.
Each module sub-block contains a ``MODULE_NAME`` and a ``MODULE_CONFIG_FILEPATH``.

BIOGEOCHEMISTRY
^^^^^^^^^^^^^^^^

+-------------------------------+---------------------------------------------------------------+
| ``MODULE_NAME``               | ``NATIVE_BGC_FLEX`` or ``PHREEQC``                            |
+-------------------------------+---------------------------------------------------------------+
| ``MODULE_CONFIG_FILEPATH``    | Path to biogeochemistry config file                           |
+-------------------------------+---------------------------------------------------------------+

See :doc:`NATIVE_BGC_FLEX <4_1_3BGC>` or :doc:`PHREEQC <4_1_3bPHREEQC>`.

TRANSPORT_DISSOLVED
^^^^^^^^^^^^^^^^^^^^^

+-------------------------------+---------------------------------------------------------------+
| ``MODULE_NAME``               | ``OPENWQ_NATIVE_TD_ADVDISP``, ``NATIVE_TD_ADV``, or ``NONE`` |
+-------------------------------+---------------------------------------------------------------+
| ``MODULE_CONFIG_FILEPATH``    | Path to transport dissolved config file                       |
+-------------------------------+---------------------------------------------------------------+

See :doc:`Transport Dissolved <4_1_4TD>`.

LATERAL_EXCHANGE
^^^^^^^^^^^^^^^^^^

+-------------------------------+---------------------------------------------------------------+
| ``MODULE_NAME``               | ``NATIVE_LE_BOUNDMIX`` or ``NONE``                            |
+-------------------------------+---------------------------------------------------------------+
| ``MODULE_CONFIG_FILEPATH``    | Path to lateral exchange config file                          |
+-------------------------------+---------------------------------------------------------------+

See :doc:`Lateral Exchange <4_1_5LE>`.

TRANSPORT_SEDIMENTS
^^^^^^^^^^^^^^^^^^^^^

+-------------------------------+---------------------------------------------------------------+
| ``MODULE_NAME``               | ``HYPE_MMF``, ``HYPE_HBVSED``, or ``NONE``                   |
+-------------------------------+---------------------------------------------------------------+
| ``SEDIMENT_COMPARTMENT``      | Compartment name where sediment is modeled                    |
+-------------------------------+---------------------------------------------------------------+
| ``MODULE_CONFIG_FILEPATH``    | Path to transport sediments config file                       |
+-------------------------------+---------------------------------------------------------------+

See :doc:`Transport Sediments <4_1_6TS>`.

SORPTION_ISOTHERM
^^^^^^^^^^^^^^^^^^^

+-------------------------------+---------------------------------------------------------------+
| ``MODULE_NAME``               | ``FREUNDLICH``, ``LANGMUIR``, or ``NONE``                     |
+-------------------------------+---------------------------------------------------------------+
| ``SEDIMENT_COMPARTMENT``      | Compartment name for sorption processes                       |
+-------------------------------+---------------------------------------------------------------+
| ``MODULE_CONFIG_FILEPATH``    | Path to sorption isotherm config file                         |
+-------------------------------+---------------------------------------------------------------+

See :doc:`Sorption Isotherms <4_1_7SI>`.


``OPENWQ_OUTPUT``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Controls the output format, species, compartments, and temporal resolution.

+----------------------------+-----------------------------------------------------------------------------------------------------+
| ``RESULTS_FOLDERPATH``     | Output folder path                                                                                  |
+----------------------------+-----------------------------------------------------------------------------------------------------+
| ``FORMAT``                 | Output file format: ``"HDF5"`` or ``"CSV"``                                                         |
+----------------------------+-----------------------------------------------------------------------------------------------------+
| ``CHEMICAL_SPECIES``       | List of chemical species names to include in output, e.g., ``["NO3-N","NH4-N"]``                    |
+----------------------------+-----------------------------------------------------------------------------------------------------+
| ``UNITS``                  | Units of chemical species output, e.g., ``"MG/L"``                                                  |
+----------------------------+-----------------------------------------------------------------------------------------------------+
| ``NO_WATER_CONC_FLAG``     | Integer flag value for concentration when there is no water, e.g., ``-9999``                        |
+----------------------------+-----------------------------------------------------------------------------------------------------+
| ``EXPORT_SEDIMENT``        | Export sediment transport results (true/false)                                                       |
+----------------------------+-----------------------------------------------------------------------------------------------------+
| ``COMPARTMENTS_AND_CELLS`` | Names of compartments and cell indices to output (see below)                                        |
+----------------------------+-----------------------------------------------------------------------------------------------------+
| ``TIMESTEP``               | Temporal resolution of output, e.g., ``[1, "hour"]``                                                |
+----------------------------+-----------------------------------------------------------------------------------------------------+

COMPARTMENTS_AND_CELLS
^^^^^^^^^^^^^^^^^^^^^^^

Each compartment is a sub-key under ``COMPARTMENTS_AND_CELLS``.
Cell selection can use ``"all"`` for all cells or numbered entries with ``[ix, iy, iz]`` indices.

+------------------------------------------------------------+---------------------------------------------------------------------+
| ``<Compt_Name>`` -> ``"1"`` -> ``["all","all","all"]``     | Output all cells of the compartment                                 |
+------------------------------------------------------------+---------------------------------------------------------------------+
| ``<Compt_Name>`` -> ``"1"`` -> ``[ix, iy, iz]``           | Output specific cell location                                       |
+------------------------------------------------------------+---------------------------------------------------------------------+


Example
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: json

    {
        "PROJECT_NAME": "steppler_N",
        "GEOGRAPHICAL_LOCATION": "steppler",
        "AUTHORS": "Author Name",
        "DATE": "March_2024",
        "COMMENT": "Example simulation",
        "COMPUTATIONAL_SETTINGS": {
            "RUN_MODE_DEBUG": true,
            "USE_NUM_THREADS": 4
        },
        "SOLVER": "FORWARD_EULER",
        "OPENWQ_INPUT": {
            "CONFIG_FILEPATH": "openwq_in/openWQ_config.json",
            "EXTERNAL_WATER_FLUXES": {
                "1": {
                    "LABEL": "SUMMA",
                    "FILEPATH": "openwq_in/openwq_EWF_mizu.json"
                }
            },
            "SINK_SOURCE": {
                "1": {
                    "LABEL": "fertilizer_N",
                    "FILEPATH": "openwq_in/openWQ_fertilizer_N.json"
                }
            }
        },
        "MODULES": {
            "BIOGEOCHEMISTRY": {
                "MODULE_NAME": "NATIVE_BGC_FLEX",
                "MODULE_CONFIG_FILEPATH": "openwq_in/openWQ_MODULE_BGC.json"
            },
            "TRANSPORT_DISSOLVED": {
                "MODULE_NAME": "NATIVE_TD_ADV",
                "MODULE_CONFIG_FILEPATH": "openwq_in/openWQ_MODULE_TD.json"
            },
            "LATERAL_EXCHANGE": {
                "MODULE_NAME": "NATIVE_LE_BOUNDMIX",
                "MODULE_CONFIG_FILEPATH": "openwq_in/openWQ_MODULE_LE.json"
            },
            "TRANSPORT_SEDIMENTS": {
                "MODULE_NAME": "NONE",
                "SEDIMENT_COMPARTMENT": "RIVER_NETWORK_REACHES",
                "MODULE_CONFIG_FILEPATH": "openwq_in/openWQ_MODULE_TS.json"
            },
            "SORPTION_ISOTHERM": {
                "MODULE_NAME": "NONE",
                "SEDIMENT_COMPARTMENT": "SUMMA_RUNOFF",
                "MODULE_CONFIG_FILEPATH": "openwq_in/openWQ_MODULE_SI.json"
            }
        },
        "OPENWQ_OUTPUT": {
            "RESULTS_FOLDERPATH": "openwq_out",
            "FORMAT": "HDF5",
            "CHEMICAL_SPECIES": ["NO3-N", "NH4-N"],
            "UNITS": "MG/L",
            "NO_WATER_CONC_FLAG": -9999,
            "EXPORT_SEDIMENT": false,
            "COMPARTMENTS_AND_CELLS": {
                "RIVER_NETWORK_REACHES": {
                    "1": ["all", "all", "all"]
                }
            },
            "TIMESTEP": [1, "hour"]
        }
    }

The JSON file supports C/C++ syntax for comments: single-line comment (``//``) or comment blocks (``/*`` and ``*/``).

The symbol ``(#)`` refers to an integer number sequence.
