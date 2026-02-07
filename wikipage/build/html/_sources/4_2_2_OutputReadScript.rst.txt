Reading and Visualizing Model Outputs
======================================

OpenWQ provides a comprehensive Python toolkit for reading, plotting, and mapping HDF5 output files. The toolkit includes three main drivers:

* **Read_h5_driver**: Reads HDF5 files into pandas DataFrames
* **Plot_h5_driver**: Creates time-series plots
* **Map_h5_driver**: Creates spatial maps and animated GIFs

Location
~~~~~~~~

The output reading scripts are located at::

    openwq/supporting_scripts/Read_Outputs/

Key files:

* ``example_reading_plotting_results.py`` - Complete usage example
* ``template_reading_plotting_results.py`` - Template for your projects
* ``hdf5_support_lib/Read_h5_driver.py`` - HDF5 reading engine
* ``hdf5_support_lib/Plot_h5_driver.py`` - Time-series plotting
* ``hdf5_support_lib/Map_h5_driver.py`` - Spatial mapping


Quick start
~~~~~~~~~~~

1. Navigate to the Read_Outputs directory::

    cd openwq/supporting_scripts/Read_Outputs/

2. Edit the example script with your paths::

    cp example_reading_plotting_results.py my_analysis.py

3. Run the script::

    python my_analysis.py


Reading HDF5 outputs
~~~~~~~~~~~~~~~~~~~~~

The ``Read_h5_driver`` function reads OpenWQ HDF5 outputs into a structured dictionary:

.. code-block:: python

    import sys
    sys.path.insert(0, 'hdf5_support_lib')
    import Read_h5_driver as h5_rlib

    # Configuration
    openwq_resDir_mapKey = {
        "path_to_results": '/path/to/openwq_out',
        "mapping_key": "reachID"  # or "hruId" for SUMMA
    }

    # Read data
    openwq_results = h5_rlib.Read_h5_driver(
        resDir_mapKey_dic=openwq_resDir_mapKey,
        output_format='HDF5',
        debugmode=False,
        cmp=['RIVER_NETWORK_REACHES'],
        space_elem='all',
        chemSpec=['NO3-N', 'NH4-N'],
        chemUnits='MG/L',
        noDataFlag=-9999
    )

**Parameters:**

+-------------------------+------------------------------------------------------------------+
| Parameter               | Description                                                      |
+=========================+==================================================================+
| ``resDir_mapKey_dic``   | Dictionary with ``path_to_results`` and ``mapping_key``          |
+-------------------------+------------------------------------------------------------------+
| ``output_format``       | Output format (``'HDF5'`` only)                                  |
+-------------------------+------------------------------------------------------------------+
| ``debugmode``           | If ``True``, also reads debug output files                       |
+-------------------------+------------------------------------------------------------------+
| ``cmp``                 | List of compartment names to read                                |
+-------------------------+------------------------------------------------------------------+
| ``space_elem``          | ``'all'`` or list of ``(ix, iy, iz)`` coordinates                |
+-------------------------+------------------------------------------------------------------+
| ``chemSpec``            | List of chemical species names                                   |
+-------------------------+------------------------------------------------------------------+
| ``chemUnits``           | Units string (must match output file naming)                     |
+-------------------------+------------------------------------------------------------------+
| ``noDataFlag``          | Value used for no-data cells (replaced with NaN)                 |
+-------------------------+------------------------------------------------------------------+

**Return structure:**

The function returns a dictionary keyed by ``"Compartment@Chemical#Units"``:

.. code-block:: python

    openwq_results = {
        'RIVER_NETWORK_REACHES@NO3-N#MG/L': [
            ('main', [(filename, DataFrame, xyz_coords)]),
            ('d_output_dt_chemistry', [...]),
            ...
        ],
        ...
    }


Creating time-series plots
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``Plot_h5_driver`` function creates time-series plots for specific features:

.. code-block:: python

    import Plot_h5_driver as h5_plib

    # Plot OpenWQ results
    figs = h5_plib.Plot_h5_driver(
        what2map='openwq',
        hostmodel='mizuroute',
        mapping_key_values=[200005899, 200002273],  # Feature IDs to plot
        openwq_results=openwq_results,
        chemSpec=['NO3-N', 'NH4-N'],
        hydromodel_info={
            'path_to_results': '/path/to/mizuroute.nc',
            'mapping_key': 'reachID'
        },
        output_path='plots/timeseries.png',
        debugmode=False,
        figsize=(12, 6)
    )

**Parameters:**

+------------------------+-------------------------------------------------------------------+
| Parameter              | Description                                                       |
+========================+===================================================================+
| ``what2map``           | ``'openwq'`` or ``'hostmodel'``                                   |
+------------------------+-------------------------------------------------------------------+
| ``hostmodel``          | ``'mizuroute'`` or ``'summa'``                                    |
+------------------------+-------------------------------------------------------------------+
| ``mapping_key_values`` | List of feature IDs to plot                                       |
+------------------------+-------------------------------------------------------------------+
| ``openwq_results``     | Results from ``Read_h5_driver`` (for ``what2map='openwq'``)       |
+------------------------+-------------------------------------------------------------------+
| ``chemSpec``           | Chemical species to plot (creates one plot per species)           |
+------------------------+-------------------------------------------------------------------+
| ``hydromodel_info``    | Dictionary with NetCDF path and mapping key                       |
+------------------------+-------------------------------------------------------------------+
| ``output_path``        | Base path for output plots                                        |
+------------------------+-------------------------------------------------------------------+
| ``debugmode``          | If ``True``, plot all debug extensions                            |
+------------------------+-------------------------------------------------------------------+

**Output files:**

For each species and extension, the function creates a separate plot:

* ``timeseries_NO3-N_main.png``
* ``timeseries_NH4-N_main.png``
* (If debugmode=True) ``timeseries_NO3-N_d_output_dt_chemistry.png``, etc.


Creating spatial maps
~~~~~~~~~~~~~~~~~~~~~~

The ``Map_h5_driver`` function creates static PNG maps and animated GIFs:

.. code-block:: python

    import Map_h5_driver as h5_mplib

    # Shapefile configuration
    shpfile_info = {
        "path_to_shp": "/path/to/rivers.shp",
        "mapping_key": "SegId"
    }

    # Create maps for OpenWQ results
    result = h5_mplib.Map_h5_driver(
        shpfile_info=shpfile_info,
        openwq_results=openwq_results,
        hydromodel_info='/path/to/mizuroute.nc',
        output_html_path='output/maps/',
        chemSpec=['NO3-N'],
        file_extension='main',
        timeframes=50,
        hostmodel='mizuroute',
        what2map='openwq',
        create_gif=True,
        gif_duration=500
    )

    # Create maps for host model outputs
    result = h5_mplib.Map_h5_driver(
        shpfile_info=shpfile_info,
        hydromodel_info={
            'path_to_results': '/path/to/mizuroute.nc',
            'mapping_key': 'SegId'
        },
        hydromodel_var2print='runoff',
        output_html_path='output/maps/',
        timeframes=50,
        hostmodel='mizuroute',
        what2map='hostmodel',
        create_gif=True,
        gif_duration=500
    )

**Parameters:**

+----------------------+--------------------------------------------------------------------+
| Parameter            | Description                                                        |
+======================+====================================================================+
| ``shpfile_info``     | Dictionary with shapefile path and mapping key                     |
+----------------------+--------------------------------------------------------------------+
| ``openwq_results``   | Results from ``Read_h5_driver`` (for ``what2map='openwq'``)        |
+----------------------+--------------------------------------------------------------------+
| ``hydromodel_info``  | Path to NetCDF (openwq) or dict with path/key (hostmodel)          |
+----------------------+--------------------------------------------------------------------+
| ``output_html_path`` | Output directory for maps                                          |
+----------------------+--------------------------------------------------------------------+
| ``chemSpec``         | Chemical species to map (list or string)                           |
+----------------------+--------------------------------------------------------------------+
| ``timeframes``       | Number of equally-spaced frames to generate                        |
+----------------------+--------------------------------------------------------------------+
| ``hostmodel``        | ``'mizuroute'`` (polylines) or ``'summa'`` (polygons)              |
+----------------------+--------------------------------------------------------------------+
| ``what2map``         | ``'openwq'`` or ``'hostmodel'``                                    |
+----------------------+--------------------------------------------------------------------+
| ``create_gif``       | If ``True``, create animated GIF from PNG frames                   |
+----------------------+--------------------------------------------------------------------+
| ``gif_duration``     | Duration per frame in milliseconds                                 |
+----------------------+--------------------------------------------------------------------+

**Output files:**

* ``maps_NO3-N/map_0000.png``, ``map_0001.png``, ... (individual frames)
* ``openwq_NO3-N_animation.gif`` (animated GIF)


GIF creation requirements
~~~~~~~~~~~~~~~~~~~~~~~~~~

The mapping tool attempts to create GIFs using multiple methods. Install one of:

* **imageio** (recommended)::

    pip install imageio

* **Pillow**::

    pip install Pillow

* **ImageMagick** (system-level)::

    # macOS
    brew install imagemagick

    # Ubuntu/Debian
    sudo apt install imagemagick


Complete example
~~~~~~~~~~~~~~~~~

.. code-block:: python

    import sys
    sys.path.insert(0, 'hdf5_support_lib')

    import Read_h5_driver as h5_rlib
    import Plot_h5_driver as h5_plib
    import Map_h5_driver as h5_mplib

    # Configuration
    openwq_resDir_mapKey = {
        "path_to_results": '/path/to/openwq_out',
        "mapping_key": "reachID"
    }

    chemical_species = ["NO3-N", "NH4-N", "N_ORG_fresh"]
    compartments = ['RIVER_NETWORK_REACHES']

    # Step 1: Read data
    openwq_results = h5_rlib.Read_h5_driver(
        resDir_mapKey_dic=openwq_resDir_mapKey,
        output_format='HDF5',
        debugmode=False,
        cmp=compartments,
        space_elem='all',
        chemSpec=chemical_species,
        chemUnits='MG/L',
        noDataFlag=-9999
    )

    # Step 2: Create time-series plots
    h5_plib.Plot_h5_driver(
        what2map='openwq',
        hostmodel='mizuroute',
        mapping_key_values=[200005899, 200002273],
        openwq_results=openwq_results,
        chemSpec=chemical_species,
        hydromodel_info={'path_to_results': '/path/to/nc', 'mapping_key': 'reachID'},
        output_path='plots/timeseries.png'
    )

    # Step 3: Create spatial maps
    h5_mplib.Map_h5_driver(
        shpfile_info={
            "path_to_shp": "/path/to/shapefile.shp",
            "mapping_key": "SegId"
        },
        openwq_results=openwq_results,
        hydromodel_info='/path/to/mizuroute.nc',
        output_html_path='maps/',
        chemSpec=chemical_species[0],
        timeframes=50,
        hostmodel='mizuroute',
        what2map='openwq',
        create_gif=True,
        gif_duration=500
    )


Troubleshooting
~~~~~~~~~~~~~~~~

**"Only supports HDF5 outputs" error:**

The reading tools currently only support HDF5 format. Ensure your OpenWQ output is configured with ``output_format = "HDF5"``.

**No matching cells found:**

Check that the ``mapping_key`` in your configuration matches the ID column in the HDF5 files. Enable debug mode to see available keys.

**GIF creation fails:**

Install one of the GIF creation packages (imageio, Pillow, or ImageMagick). The tool will try each method in order.

**Feature ID not found:**

Ensure the feature IDs in ``mapping_key_values`` exist in your shapefile. Check the shapefile's attribute table for the correct ID column.
