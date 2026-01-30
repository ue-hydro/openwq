# Copyright 2026, Diogo Costa, diogo.costa@uevora.pt
# This file is part of OpenWQ model.

# This program, openWQ, is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
OpenWQ Simplified Time-Series Plotting Tool
Clean interface for plotting time-series
Supports multiple chemical species (one plot per species per extension)
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
                   debugmode=False,
                   combine_features=True,
                   figsize=(12, 6)):
    """
    Plot time-series for specific features
    For OpenWQ: Creates one plot per chemical species per extension

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
        Chemical species to plot (e.g., ['NO3-N', 'NH4-N'] or 'NO3-N')
        Creates one plot per species per extension
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
        For multiple species/extensions, names are appended
    debugmode : bool
        If False: only process 'main' extension (default)
        If True: process 'main' + all debug extensions:
                 - d_output_dt_chemistry
                 - d_output_dt_transport
                 - d_output_dt_ss
                 - d_output_dt_ewf
                 - d_output_dt_ic
    combine_features : bool
        If True, plot all features on one graph (default: True)
    figsize : tuple
        Figure size (width, height) in inches (default: (12, 6))

    Returns:
    --------
    str or list : Path(s) to created plot(s)
        - Single path for hostmodel
        - List of paths for openwq (one per species per extension)

    Examples:
    ---------
    # Normal mode (only 'main')
    h5_mplib.Plot_h5_driver(
        what2map='openwq',
        hostmodel='mizuroute',
        mapping_key_values=[200005899, 200002273],
        openwq_results=openwq_results,
        chemSpec=['NO3-N', 'NH4-N'],
        hydromodel_info=hydromodel_info,
        hydromodel_var2print='DWroutedRunoff',
        output_path='plots/nutrients.png',
        debugmode=False
    )
    # Creates: plots/nutrients_NO3-N_main.png
    #          plots/nutrients_NH4-N_main.png

    # Debug mode (all extensions)
    h5_mplib.Plot_h5_driver(
        ...,
        debugmode=True
    )
    # Creates: plots/nutrients_NO3-N_main.png
    #          plots/nutrients_NO3-N_d_output_dt_chemistry.png
    #          plots/nutrients_NO3-N_d_output_dt_transport.png
    #          ... (and so on for each species)
    """

    print("=" * 70)
    print("TIME-SERIES PLOTTING")
    print(f"Mode: {what2map.upper()}")
    if what2map == 'openwq':
        print(f"Debug mode: {debugmode}")
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

        nc_path = hydromodel_info['path_to_results']
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
            units = nc.variables[hydromodel_var2print].units if hasattr(nc.variables[hydromodel_var2print],
                                                                        'units') else 'units'

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

            # CREATE PLOT FOR HOSTMODEL
            print("\n" + "=" * 70)
            print("Creating plot...")
            print("=" * 70)

            fig, ax = plt.subplots(figsize=figsize)

            # Plot each feature
            for fid, data in feature_data.items():
                ax.plot(data.index, data.values, label=f'Feature {fid}', linewidth=2)

            # Format plot
            plot_title = f'{hostmodel.upper()} - {hydromodel_var2print} Time Series'
            ylabel = f'{hydromodel_var2print} ({units})'

            ax.set_xlabel('Time', fontsize=12)
            ax.set_ylabel(ylabel, fontsize=12)
            ax.set_title(plot_title, fontsize=14, pad=20)

            if len(feature_data) > 1:
                ax.legend(loc='best', fontsize=10)

            ax.grid(True, alpha=0.3)

            # Format x-axis for dates
            if isinstance(time_index, pd.DatetimeIndex):
                ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
                fig.autofmt_xdate()

            plt.tight_layout()

            # Save plot
            if output_path is None:
                output_path = f'timeseries_{hydromodel_var2print}.png'
            else:
                # Insert hydromodel_var2print name and extension before file extension
                base_path = os.path.splitext(output_path)[0]
                extension = os.path.splitext(output_path)[1]
                hydro_output_path = f'{base_path}_{hydromodel_var2print}{extension}'


            plt.savefig(hydro_output_path, dpi=150, bbox_inches='tight')
            plt.close()

            print(f"\n✓ Plot saved: {hydro_output_path}")

            # Summary
            print("\n" + "=" * 70)
            print("✓ COMPLETE!")
            print("=" * 70)
            print(f"  Features plotted: {len(feature_data)}")
            print(f"  Output: {output_path}")
            print("=" * 70)

            return output_path

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

        # Define extensions to process based on debugmode
        if debugmode:
            file_extensions = [
                'main',
                'd_output_dt_chemistry',
                'd_output_dt_transport',
                'd_output_ss',
                'd_output_ewf',
                'd_output_ic'
            ]
        else:
            file_extensions = ['main']

        print(f"  Chemical species: {chemSpec}")
        print(f"  Extensions to process: {file_extensions}")
        print(f"  Will create {len(chemSpec) * len(file_extensions)} plot(s)")

        # Store all created plots
        all_output_paths = []

        # LOOP THROUGH EACH CHEMICAL SPECIES
        for spec_idx, spec in enumerate(chemSpec):
            print(f"\n{'=' * 70}")
            print(f"Processing species {spec_idx + 1}/{len(chemSpec)}: {spec}")
            print("=" * 70)

            try:
                # Find matching species in openwq_results
                species_key = None
                for key in openwq_results.keys():
                    if spec in key:
                        species_key = key
                        print(f"  Found species key: {species_key}")
                        break

                if species_key is None:
                    print(f"  ✗ No data found for species: {spec}")
                    print(f"  Available species:")
                    for key in list(openwq_results.keys())[:5]:
                        print(f"    - {key}")
                    continue

                # Parse units from species_key (do once per species)
                # Format: WATERSTATE@NO3-N#mg-N/L
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
                    chem_name = spec

                # LOOP THROUGH EACH FILE EXTENSION
                for ext_idx, file_extension in enumerate(file_extensions):
                    print(f"\n  {'- ' * 33}")
                    print(f"  Processing extension {ext_idx + 1}/{len(file_extensions)}: {file_extension}")
                    print(f"  {'- ' * 33}")

                    # Extract data from openwq_results structure
                    data_found = False
                    for ext, data_list in openwq_results[species_key]:
                        if ext == file_extension:
                            if data_list and len(data_list) > 0:
                                filename, ttdata, xyz_coords = data_list[0]
                                data_found = True
                                print(f"    Data shape: {ttdata.shape}")
                                print(f"    Time range: {ttdata.index[0]} to {ttdata.index[-1]}")
                                break

                    if not data_found:
                        print(f"    ✗ No data found for extension '{file_extension}'")
                        continue

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
                            print(f"    ✓ Found data for feature {fid}")
                        else:
                            missing_features.append(fid)
                            print(f"    ✗ Feature {fid} not found")

                    if len(feature_data) == 0:
                        print(f"    ✗ No data found for any requested features")
                        continue

                    if missing_features:
                        print(f"    ⚠ Missing features: {missing_features}")

                    # CREATE PLOT FOR THIS SPECIES AND EXTENSION
                    print(f"\n    Creating plot for {chem_name} ({file_extension})...")

                    fig, ax = plt.subplots(figsize=figsize)

                    # Plot each feature
                    for fid, data in feature_data.items():
                        ax.plot(data.index, data.values, label=f'Feature {fid}', linewidth=2)

                    # Format plot
                    plot_title = f'{chem_name} Concentrations Time Series ({file_extension})'
                    ylabel_text = f'{chem_name} ({units})'

                    ax.set_xlabel('Time', fontsize=12)
                    ax.set_ylabel(ylabel_text, fontsize=12)
                    ax.set_title(plot_title, fontsize=14, pad=20)

                    if len(feature_data) > 1:
                        ax.legend(loc='best', fontsize=10)

                    ax.grid(True, alpha=0.3)

                    # Format x-axis for dates
                    if isinstance(ttdata.index, pd.DatetimeIndex):
                        ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
                        fig.autofmt_xdate()

                    plt.tight_layout()

                    # Generate output path for this species and extension
                    if output_path is None:
                        species_output_path = f'timeseries_{spec}_{file_extension}.png'
                    else:
                        # Insert species name and extension before file extension
                        base_path = os.path.splitext(output_path)[0]
                        extension = os.path.splitext(output_path)[1]
                        species_output_path = f'{base_path}_{spec}_{file_extension}{extension}'

                    # Save plot
                    plt.savefig(species_output_path, dpi=150, bbox_inches='tight')
                    plt.close()

                    print(f"    ✓ Plot saved: {species_output_path}")
                    all_output_paths.append(species_output_path)

            except Exception as e:
                print(f"  ✗ Error processing species {spec}: {e}")
                import traceback
                traceback.print_exc()
                continue

        # Summary
        print("\n" + "=" * 70)
        print("✓ COMPLETE!")
        print("=" * 70)
        print(f"  Species requested: {len(chemSpec)}")
        print(f"  Extensions processed: {len(file_extensions)}")
        print(f"  Plots created: {len(all_output_paths)}")
        if len(all_output_paths) > 0:
            print(f"\n  Output files:")
            for path in all_output_paths:
                print(f"    - {os.path.basename(path)}")
        print("=" * 70)

        # Return list of created plots or None
        if len(all_output_paths) == 0:
            return None
        elif len(all_output_paths) == 1:
            return all_output_paths[0]
        else:
            return all_output_paths
