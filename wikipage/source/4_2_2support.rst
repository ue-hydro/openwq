Reading Outputs
========================================

OpenWQ outputs can be saved in three standard scientific formats: `CSV <https://en.wikipedia.org/wiki/Comma-separated_values>`_, `HDF5 <https://en.wikipedia.org/wiki/Hierarchical_Data_Format>`_, and `VTU <https://people.math.sc.edu/Burkardt/data/vtu/vtu.html>`_.
These formats are supported by widely used programming languages including C++, R, Python, and MATLAB.

**HDF5 is recommended** for most applications due to its efficient handling of large datasets and support for structured data.

Supporting scripts
~~~~~~~~~~~~~~~~~~

OpenWQ provides a comprehensive Python toolkit for reading, analyzing, and visualizing model outputs. The toolkit includes:

* **Read_h5_driver**: Reads HDF5 files into pandas DataFrames
* **Plot_h5_driver**: Creates time-series plots for selected features
* **Map_h5_driver**: Creates spatial maps and animated GIFs

For detailed documentation on using these scripts, see the :doc:`Supporting scripts <4_2_2support_scripts>` section.


Quick read with Python
~~~~~~~~~~~~~~~~~~~~~~~

For simple data extraction, you can read HDF5 files directly with ``h5py``:

.. code-block:: python

    import h5py
    import numpy as np

    # Open the output file
    f = h5py.File('openwq_out/HDF5/COMPARTMENT@SPECIES#UNITS-main.h5', 'r')

    # List available timestamps
    print(list(f.keys()))

    # Read data for a specific timestamp
    data = f['2020Jan01-12:00:00'][:]
    print(data.shape)

    f.close()

For more advanced analysis including time-series extraction, spatial mapping, and animated visualizations, use the provided :doc:`supporting scripts <4_2_2support_scripts>`.

