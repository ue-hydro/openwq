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
Read HDF5 data and save to timeseries collection - OPTIMIZED VERSION
"""

import os
import h5py
import numpy as np
import pandas as pd
from tqdm import tqdm


def _normalize_filename(file_name, file_extensions_i):
    """
    Normalize filename by removing spaces, converting to uppercase, and replacing slashes.

    Parameters
    ----------
    file_name : str
        Original filename
    file_extensions_i : str
        File extension identifier

    Returns
    -------
    str
        Normalized filename
    """
    return f"{file_name.replace(' ', '').upper().replace('/', '|')}-{file_extensions_i}"


def _find_matching_cells(xyz_elements_source, xyz_elements_requested):
    """
    Find indices of matching cells between source and requested coordinates.

    Parameters
    ----------
    xyz_elements_source : np.ndarray
        Source coordinates array (3, N)
    xyz_elements_requested : np.ndarray
        Requested coordinates array (3, M)

    Returns
    -------
    tuple
        (indices_validated, coordinates_validated)
    """
    indices_validated = []
    coordinates_validated = []

    num_requested = xyz_elements_requested.shape[1]

    for l in range(num_requested):
        # Vectorized matching for all three dimensions
        match_x = np.any(xyz_elements_source[0, :] == xyz_elements_requested[0, l])
        match_y = np.any(xyz_elements_source[1, :] == xyz_elements_requested[1, l])
        match_z = np.any(xyz_elements_source[2, :] == xyz_elements_requested[2, l])

        if match_x and match_y and match_z:
            indices_validated.append(l)
            coordinates_validated.append([
                xyz_elements_requested[0, l] - 1,
                xyz_elements_requested[1, l] - 1,
                xyz_elements_requested[2, l] - 1
            ])
        else:
            print(f"Requested cell not found: {xyz_elements_requested[:, l]} - entry skipped")

    return indices_validated, np.array(coordinates_validated)


def Read_h5_save_engine(
        openwq_info,
        file_extensions_i,
        file_name,
        space_elem,
        noDataFlag):
    """
    Read HDF5 data and save to timeseries collection (OPTIMIZED).

    Parameters
    ----------
    folderpath : str
        Path to folder containing HDF5 files
    file_extensions_i : str
        File extension identifier
    file_name : str
        Name of the file to read
    space_elem : array-like or str
        Spatial cell coordinates or 'all'
    noDataFlag: number

    Returns
    -------
    list
        List of tuples containing (filename, timetable, xyz_elements)
    """

    folderpath = openwq_info["path_to_results"]
    mappingKey = openwq_info["mapping_key"]

    # Normalize and find the file
    filename_fix = _normalize_filename(file_name, file_extensions_i)
    filepath_i = os.path.join(folderpath, f"{filename_fix}.h5")

    # Check file existence
    if not os.path.isfile(filepath_i):
        print(f"<Read_h5_save_tscollection> Warning: could not find \"{filename_fix}.h5\" file. Request skipped.")
        return []

    # Single file open for all operations to improve performance
    try:
        with h5py.File(filepath_i, 'r') as hf:
            # Get datasets (timestamps)
            all_keys = list(hf.keys())
            timestamps = all_keys[:-3]  # Exclude last 2 entries

            # Read xyz_elements
            if '/xyz_elements' not in hf:
                print(f"<main_hdf5> Warning: '/xyz_elements' not found in {filename_fix}.h5")
                return []
            # print(all_keys)
            xyz_elements_source = hf['/xyz_elements'][:]

            # Skip if no data
            if xyz_elements_source.size == 0:
                return []

            # Process requested cells
            if isinstance(space_elem, str) and space_elem.lower() == "all":
                xyz_elements_requested = xyz_elements_source
            else:
                xyz_elements_requested = np.array(space_elem)

            # Find matching cells (optimized)
            indices_validated, coordinates_validated = _find_matching_cells(
                xyz_elements_source, xyz_elements_requested
            )

            if len(indices_validated) == 0:
                print(f"<main_hdf5> Warning: No matching cells found for {filename_fix}.h5")
                return []

            num_validated = len(indices_validated)
            num_timesteps = len(timestamps)

            # Get corresponding hostmodel ids for the indices_validated
            mapKey_values = np.array([x.decode('utf-8') if isinstance(x, bytes) else x for x in hf[f'/{mappingKey}'][:]])

            # Pre-allocate arrays
            data_all = np.full((num_timesteps, num_validated), np.nan, dtype=np.float64)
            time_all = []

            # Read all timestep data in single file context
            for tstep, timestamp in enumerate(timestamps):
                dataset_path = f'/{timestamp}'
                if dataset_path not in hf:
                    print(f"Warning: Dataset {dataset_path} not found, skipping")
                    continue

                # Read data for this timestep
                data_i = hf[dataset_path][:]

                # Replace noDataFlag with NaN
                data_i[data_i==noDataFlag] = np.nan
                #data_i[data_i != noDataFlag] = 0.1 * tstep # TO REMOVE AFTER DEBUGGING

                data_all[tstep, :] = data_i[0, indices_validated]
                time_all.append(timestamp)

    except (FileNotFoundError, OSError) as e:
        print(f"<main_hdf5> Warning: Could not open file {filepath_i}: {e}")
        return []
    except Exception as e:
        print(f"<main_hdf5> Warning: Error reading {filepath_i}: {e}")
        return []

    # Convert and sort timestamps
    try:
        time_all_dt = pd.to_datetime(time_all, format='%Y%b%d-%H:%M:%S')
    except Exception:
        try:
            time_all_dt = pd.to_datetime(time_all)
        except Exception:
            print(f"Warning: Could not parse timestamps for {filename_fix}")
            return []

    # Sort by time (using argsort for efficiency)
    sorted_indices = np.argsort(time_all_dt)
    time_sorted = time_all_dt[sorted_indices]
    data_sorted = data_all[sorted_indices, :]

    # Create DataFrame
    ttdata = pd.DataFrame(
        data_sorted,
        index=time_sorted,
        columns=[f"{mappingKey}_{j}" for j in mapKey_values]
    )

    # Return as list with single tuple
    return [(filename_fix, ttdata, coordinates_validated)]


def Read_h5_driver(openwq_info=None,
                   output_format=None,
                   debugmode=None,
                   cmp=None,
                   space_elem=None,
                   chemSpec=None,
                   chemUnits=None,
                   noDataFlag=None):
    """
    Read and plot OpenWQ output data (HDF5 and CSV) - OPTIMIZED VERSION.

    Parameters
    ----------
    folderpath : str
        Fullpath to directory where the HDF5 files are located
    output_format : str
        Output format ('HDF5' or 'CSV')
    debugmode : bool
        Debug mode flag
    cmp : list
        List of cmp to extract
    space_elem : array-like or str
        Spatial cell coordinates or 'all'
    chemSpec : list
        List of chemical species to extract
    chemUnits : str
        Units for concentration

    Returns
    -------
    dict
        Dictionary with keys as "Compartment@Chemical#Units" and values as results
    """

    # Validate output format
    if output_format != 'HDF5':
        print(">> Only supports HDF5 outputs <<")
        return {}

    # Updating fullpath to outputs
    openwq_info["path_to_results"] = os.path.join(
        openwq_info["path_to_results"]
        , output_format, '')

    # Define file extensions
    file_extensions = [
        'main',
        'd_output_dt_chemistry',
        'd_output_dt_transport',
        'd_output_ss',
        'd_output_ewf',
        'd_output_ic'
    ]

    # Determine number of files to load
    load_files_num = len(file_extensions) if debugmode else 1

    # Pre-allocate results dictionary
    openwq_results = {}

    # Calculate total iterations for more accurate progress tracking
    total_iterations = len(chemSpec) * len(cmp) * load_files_num

    # Use single progress bar for all operations
    with tqdm(total=total_iterations, desc="Extracting OpenWQ data") as pbar:

        for chem in chemSpec:
            print(f"\n> Extracting results for: {chem}")

            for comp in cmp:
                print(f"  > Extracting results for: {comp}")

                # Initialize results list for this combination
                output_openwq_tscollect_all = []

                # Create composite key
                result_key = f"{comp}@{chem}#{chemUnits}"

                # Loop through file extensions
                for f in range(load_files_num):
                    file_extensions_i = file_extensions[f]

                    # Create file name
                    file_name = f"{comp}@{chem}#{chemUnits}"

                    # Read data
                    output_openwq_tscollect = Read_h5_save_engine(
                        openwq_info,
                        file_extensions_i,
                        file_name,
                        space_elem,
                        noDataFlag
                    )

                    # Save results
                    output_openwq_tscollect_all.append((
                        file_extensions_i,
                        output_openwq_tscollect
                    ))

                    # Update progress
                    pbar.update(1)

                # Store in results dictionary
                openwq_results[result_key] = output_openwq_tscollect_all

    return openwq_results

