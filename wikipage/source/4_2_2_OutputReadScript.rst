Reading model outputs
========================================

OpenWQ supports three output formats: HDF5, CSV, and VTU (see :doc:`Output formats <4_2_1file>`).
Below are examples for reading each format.


Reading HDF5 output (Python)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

HDF5 is the recommended output format for large simulations. Use the ``h5py`` library in Python::

    import h5py
    import numpy as np

    # Open the output file
    f = h5py.File('output/openwq_output.h5', 'r')

    # List available groups (compartments)
    print(list(f.keys()))

    # Read chemical species data for a compartment
    data = f['COMPARTMENT_NAME']['species_A'][:]

    # Data shape: (time_steps, ix, iy, iz)
    print(data.shape)

    f.close()

Or using ``pandas`` for tabular analysis::

    import pandas as pd
    import h5py

    f = h5py.File('output/openwq_output.h5', 'r')
    # Extract a time series at a specific cell
    ts = f['COMPARTMENT_NAME']['species_A'][:, 0, 0, 0]
    f.close()


Reading CSV output (Python)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

CSV output files can be read directly with ``pandas``::

    import pandas as pd

    df = pd.read_csv('output/openwq_output.csv')
    print(df.head())
    print(df.columns)


Reading VTU output (ParaView)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

VTU (VTK Unstructured Grid) output files can be visualized using `ParaView <https://www.paraview.org/>`_:

1. Open ParaView
2. File > Open, navigate to the output directory
3. Select the ``.vtu`` files
4. Click Apply to render
5. Use the color map selector to visualize different chemical species

VTU output is particularly useful for spatially distributed simulations where 2D or 3D visualization is needed.
