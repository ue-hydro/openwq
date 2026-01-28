"""
OpenWQ Simplified Time-Series Plotting Tool
Clean interface for plotting time-series
"""

import numpy as np
import pandas as pd
import netCDF4
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import os


def Plot_h5_driver(what2map=None,
                   hostmodel=None,
                   mapping_key_values=None,
                   openwq_results=None,
                   chemSpec=None,
                   hydromodel_info=None,
                   hydromodel_var2print=None,
                   output_path=None,
                   file_extension='main',
                   combine_features=True,
                   figsize=(12, 6)):
    """
    Plot time-series for specific features

    Parameters:
    -----------
    what2map : str
        'hostmodel' or 'openwq'
    hostmodel : str
        'mizuroute' or 'summa'
    mapping_key_values : list
        List of feature IDs to plot (e.g., [200005899, 200002273])
    openwq_results : dict
        OpenWQ results dictionary from Read_h5_driver
        Used if what2map='openwq'
    chemSpec : list or str
        Chemical species to plot (e.g., ['NO3-N'] or 'NO3-N')
        Used if what2map='openwq'
    hydromodel_info : dict
        Dictionary with:
            - 'path_to_shp': Path to NetCDF file
            - 'mapping_key': Variable for IDs (e.g., 'reachID')
        Used if what2map='hostmodel'
    hydromodel_var2print : str
        Variable to plot (e.g., 'basRunoff', 'DWroutedRunoff')
        Used if what2map='hostmodel'
    output_path : str
        Path where plot will be saved
    file_extension : str
        File extension for OpenWQ data (default: 'main')
    combine_features : bool
        If True, plot all features on one graph (default: True)
    figsize : tuple
        Figure size (width, height) in inches (default: (12, 6))

    Returns:
    --------
    str : Path to created plot

    Examples:
    ---------
    # 1. Load OpenWQ results
    from Read_h5_driver import Read_h5_driver
    openwq_results = Read_h5_driver(
        h5files_folderPath='/path/to/openwq_out',
        hostmodel='mizuroute'
    )

    # 2. Plot mizuRoute runoff
    h5_mplib.Plot_h5_driver(
        # 1) What results to map?
        what2map='hostmodel',
        hostmodel='mizuroute',
        mapping_key_values=[200005899, 200002273],
        # 2) openwq info (👉🏼 used if what2map=openwq)
        openwq_results=openwq_results,
        chemSpec=['NO3-N'],
        # 3) hostmodel info (👉🏼 used if what2map=hostmodel)
        hydromodel_info=hydromodel_info,
        hydromodel_var2print='DWroutedRunoff',
        # 4) output config
        output_path='plots/runoff.png'
    )

    # 3. Plot OpenWQ NO3
    h5_mplib.Plot_h5_driver(
        # 1) What results to map?
        what2map='openwq',
        hostmodel='mizuroute',
        mapping_key_values=[200005899, 200002273],
        # 2) openwq info (👉🏼 used if what2map=openwq)
        openwq_results=openwq_results,
        chemSpec=['NO3-N'],
        # 3) hostmodel info (👉🏼 used if what2map=hostmodel)
        hydromodel_info=hydromodel_info,
        hydromodel_var2print='DWroutedRunoff',
        # 4) output config
        output_path='plots/NO3.png'
    )
    """

    print("=" * 70)
    print("TIME-SERIES PLOTTING")
    print(f"Mode: {what2map.upper()}")
    print("=" * 70)

    # Validate inputs
    if what2map not in ['hostmodel', 'openwq']:
        print(f"✗ Error: what2map must be 'hostmodel' or 'openwq', got '{what2map}'")
        return None

    if mapping_key_values is None or len(mapping_key_values) == 0:
        print("✗ Error: mapping_key_values must be provided")
        return None

    # Convert to list if single value
    if not isinstance(mapping_key_values, list):
        mapping_key_values = [mapping_key_values]

    print(f"Features to plot: {mapping_key_values}")

    # Create output directory if needed
    if output_path:
        output_dir = os.path.dirname(output_path)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir, exist_ok=True)

    # BRANCH: Load data based on what2map
    if what2map == 'hostmodel':
        print("\n" + "=" * 70)
        print("Loading hostmodel data...")
        print("=" * 70)

        if hydromodel_info is None:
            print("✗ Error: hydromodel_info must be provided for hostmodel mode")
            return None

        if hydromodel_var2print is None:
            print("✗ Error: hydromodel_var2print must be specified for hostmodel mode")
            return None

        nc_path = hydromodel_info['path_to_shp']
        mapping_key = hydromodel_info['mapping_key']

        print(f"  NetCDF: {nc_path}")
        print(f"  Variable: {hydromodel_var2print}")
        print(f"  Mapping key: {mapping_key}")

        # Open NetCDF
        nc = netCDF4.Dataset(nc_path, 'r')

        try:
            # Get IDs and find indices for requested features
            all_ids = np.array(nc.variables[mapping_key][:])

            print(f"  Total features in file: {len(all_ids)}")

            # Get data variable
            data_var = nc.variables[hydromodel_var2print]

            # Get time
            if 'time' in nc.variables:
                time_var = nc.variables['time']
                time_data = time_var[:]

                if hasattr(time_var, 'units'):
                    from netCDF4 import num2date
                    dates = num2date(time_data, time_var.units, only_use_cftime_datetimes=False)
                    time_index = pd.DatetimeIndex(dates)
                else:
                    time_index = pd.RangeIndex(len(time_data))
            else:
                time_index = pd.RangeIndex(data_var.shape[0])

            print(f"  Time range: {time_index[0]} to {time_index[-1]}")
            print(f"  Time steps: {len(time_index)}")

            # Get units
            units = nc.variables[hydromodel_var2print].units if hasattr(nc.variables[hydromodel_var2print], 'units') else 'units'

            # Extract data for requested features
            feature_data = {}
            missing_features = []

            for fid in mapping_key_values:
                idx = np.where(all_ids == fid)[0]

                if len(idx) > 0:
                    data = data_var[:, idx[0]]
                    feature_data[fid] = pd.Series(data, index=time_index)
                    print(f"  ✓ Found data for feature {fid}")
                else:
                    missing_features.append(fid)
                    print(f"  ✗ Feature {fid} not found")

            nc.close()

            if len(feature_data) == 0:
                print("✗ Error: No data found for any requested features")
                return None

            if missing_features:
                print(f"\n⚠ Warning: Missing features: {missing_features}")

            plot_title = f'{hostmodel.upper()} - {hydromodel_var2print} Time Series'
            ylabel = f'{hydromodel_var2print} ({units})'

        except Exception as e:
            print(f"✗ Error loading hostmodel data: {e}")
            nc.close()
            return None

    else:  # what2map == 'openwq'
        print("\n" + "=" * 70)
        print("Loading OpenWQ data...")
        print("=" * 70)

        if openwq_results is None:
            print("✗ Error: openwq_results must be provided for openwq mode")
            return None

        if chemSpec is None:
            print("✗ Error: chemSpec must be specified for openwq mode")
            return None

        # Convert chemSpec to list if string
        if isinstance(chemSpec, str):
            chemSpec = [chemSpec]

        print(f"  Chemical species: {chemSpec}")

        try:
            # Find matching species in openwq_results
            species_key = None
            for key in openwq_results.keys():
                # Check if any of the requested species is in the key
                for spec in chemSpec:
                    if spec in key:
                        species_key = key
                        print(f"  Found species: {species_key}")
                        break
                if species_key:
                    break

            if species_key is None:
                print(f"✗ Error: No data found for species: {chemSpec}")
                print(f"  Available species:")
                for key in list(openwq_results.keys())[:5]:
                    print(f"    - {key}")
                return None

            # Extract data from openwq_results structure
            data_found = False
            for ext, data_list in openwq_results[species_key]:
                if ext == file_extension:
                    if data_list and len(data_list) > 0:
                        filename, ttdata, xyz_coords = data_list[0]
                        data_found = True
                        print(f"  Data shape: {ttdata.shape}")
                        print(f"  Time range: {ttdata.index[0]} to {ttdata.index[-1]}")
                        break

            if not data_found:
                print(f"✗ Error: No data found for extension '{file_extension}'")
                return None

            # Parse units from species_key
            # Format: WATERSTATE@NO3#mg-N/L
            parts = species_key.split('#')
            if len(parts) > 1:
                units = parts[1]
            else:
                units = 'concentration'

            # Get chemical name
            chem_parts = species_key.split('@')
            if len(chem_parts) > 1:
                chem_name = chem_parts[1].split('#')[0]
            else:
                chem_name = chemSpec[0]

            # Extract data for requested features
            feature_data = {}
            missing_features = []

            # ttdata columns are like 'REACH_200005899'
            for fid in mapping_key_values:
                # Find column matching this feature ID
                matching_cols = [col for col in ttdata.columns if str(fid) in col.split('_')[-1]]

                if len(matching_cols) > 0:
                    data = ttdata[matching_cols[0]]
                    feature_data[fid] = data
                    print(f"  ✓ Found data for feature {fid}")
                else:
                    missing_features.append(fid)
                    print(f"  ✗ Feature {fid} not found")

            if len(feature_data) == 0:
                print("✗ Error: No data found for any requested features")
                return None

            if missing_features:
                print(f"\n⚠ Warning: Missing features: {missing_features}")

            plot_title = f'{chem_name} Concentrations Time Series'
            ylabel = f'{chem_name} ({units})'

        except Exception as e:
            print(f"✗ Error loading OpenWQ data: {e}")
            import traceback
            traceback.print_exc()
            return None

    # CREATE PLOT
    print("\n" + "=" * 70)
    print("Creating plot...")
    print("=" * 70)

    fig, ax = plt.subplots(figsize=figsize)

    # Plot each feature
    for fid, data in feature_data.items():
        ax.plot(data.index, data.values, label=f'Feature {fid}', linewidth=2)

    # Format plot
    ax.set_xlabel('Time', fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title(plot_title, fontsize=14, pad=20)

    if len(feature_data) > 1:
        ax.legend(loc='best', fontsize=10)

    ax.grid(True, alpha=0.3)

    # Format x-axis for dates
    if isinstance(data.index, pd.DatetimeIndex):
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
        fig.autofmt_xdate()

    plt.tight_layout()

    # Save plot
    if output_path is None:
        output_path = f'timeseries_plot.png'

    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()

    print(f"\n✓ Plot saved: {output_path}")

    # Summary
    print("\n" + "=" * 70)
    print("✓ COMPLETE!")
    print("=" * 70)
    print(f"  Features plotted: {len(feature_data)}")
    print(f"  Output: {output_path}")
    print("=" * 70)

    return output_path


if __name__ == '__main__':
    print("\nOpenWQ Simplified Time-Series Plotting")
    print("\nSimple interface for plotting time-series data")
    print("\nSupports:")
    print("  • Hostmodel outputs (NetCDF)")
    print("  • OpenWQ results (from Read_h5_driver)")