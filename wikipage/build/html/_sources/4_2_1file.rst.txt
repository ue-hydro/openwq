Output formats
==================================

OpenWQ supports the following output data formats:

* `HDF5 <https://en.wikipedia.org/wiki/Hierarchical_Data_Format>`_ (json key/value: ``"FORMAT"``: ``"HDF5"``) -- **Preferred format**. Full support including the post-processing supporting scripts for reading, plotting, and mapping results.
* `CSV <https://en.wikipedia.org/wiki/Comma-separated_values>`_ (json key/value: ``"FORMAT"``: ``"CSV"``)

.. note::

    HDF5 is the recommended output format as it is the only format with full support, including the supporting scripts for reading, plotting, and mapping results.

The results will be printed in a folder with the name ``Output_OpenWQ``.

To select the data type go to the ``Master`` file and change the value of the sub-key ```"FORMAT"`` inside key ``"OPENWQ_OUTPUT"``

Example:

.. code-block:: json

    {
        "...": "...",
        "OPENWQ_OUTPUT": {
            "...": "...",
            "FORMAT": "HDF5",
            "...": "...",
            }
        "...": "...",
        }
    }



