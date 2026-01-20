# Copyright 2020, Diogo Costa, diogo.costa@uevora.pt
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
Interactive Plotly-based mapping for OpenWQ model results
Uses the results dictionary from Read_h5_driver
Creates interactive animated maps that can be viewed in browser
"""

import netCDF4
import h5py
import geopandas as gpd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from datetime import datetime
import os
import json
import pandas as pd


def Map_h5_driver(shpfile_fullpath,
                         openwq_results,
                         hydromodel_out_fullpath,
                         output_html_path,
                         species_name=None,
                         file_extension='main',
                         timestep='all',
                         vmin=None,
                         vmax=None,
                         colorscale='RdYlGn_r',
                         animation_frame_duration=500,
                         mapbox_style='open-street-map',
                         zoom=6,
                         center=None,
                         model='mizuroute'):
    """
    Create interactive Plotly map from OpenWQ results dictionary.

    Parameters
    ----------
    shpfile_fullpath : str
        Full path to shapefile
    openwq_results : dict
        Results dictionary from Read_h5_driver with structure:
        {
            'COMPARTMENT@CHEMICAL#UNITS': [
                ('main', [(filename, DataFrame, xyz_coords)]),
                ('d_output_dt_chemistry', [(filename, DataFrame, xyz_coords)]),
                ...
            ]
        }
    hydromodel_out_fullpath : str
        Full path to hydrological model (mizuRoute) NetCDF output file
    output_html_path : str
        Output HTML filename
    species_name : str, optional
        Which species to plot. If None, uses first species in dictionary.
        Example: 'SPECIES_A' or 'RIVER_NETWORK_REACHES@SPECIES_A#MG/L'
    file_extension : str, optional
        Which file extension to plot (default: 'main')
        Options: 'main', 'd_output_dt_chemistry', 'd_output_dt_transport',
                 'd_output_ss', 'd_output_ewf', 'd_output_ic'
    timestep : int, str, list, or 'all', optional
        Which timestep(s) to plot:
        - 'all': Create animated map with all timesteps (default)
        - int: Create static map for specific timestep index (e.g., 0, 10, 50)
        - str: Create static map for specific date (e.g., '2020-07-15')
        - list [min, max]: Create animated map for timestep range
          * [0, 100]: Indices 0 to 100
          * ['2020-01-01', '2020-12-31']: Date range
          * ['1961JAN01-08:00:00', '1961JAN04-03:00:00']: Date range with time
    vmin : float, optional
        Minimum value for color scale (default: auto from data)
    vmax : float, optional
        Maximum value for color scale (default: auto from data)
    colorscale : str, optional
        Plotly colorscale name (default: 'RdYlGn_r')
        Options: 'Viridis', 'Plasma', 'Inferno', 'Turbo', 'RdYlGn', 'RdYlGn_r'
    animation_frame_duration : int, optional
        Duration of each frame in milliseconds (default: 500)
        Only used when timestep='all'
    mapbox_style : str, optional
        Map style (default: 'open-street-map')
        Options: 'open-street-map', 'carto-positron', 'carto-darkmatter', 
                 'stamen-terrain', 'stamen-toner', 'stamen-watercolor'
    zoom : int, optional
        Initial map zoom level (default: 6)
    center : dict, optional
        Map center as {'lat': latitude, 'lon': longitude}
        If None, auto-centers on data
    model : str, optional
        Hydrological model type (default: 'mizuroute')

    Returns
    -------
    str
        Path to created HTML file

    Examples
    --------
    # Read data first
    >>> from Read_h5_optimized import Read_h5_driver
    >>> results = Read_h5_driver(
    ...     folderpath='/path/to/output/',
    ...     OpenWQ_output_format='HDF5',
    ...     Compartments=['RIVER_NETWORK_REACHES'],
    ...     Chemical_species=['SPECIES_A'],
    ...     Concentration_units='MG/L'
    ... )

    # Create animated map (all timesteps)
    >>> Map_h5_driver(
    ...     shpfile_fullpath='/path/to/river.shp',
    ...     openwq_results=results,
    ...     hydromodel_out_fullpath='/path/to/mizuroute.nc',
    ...     output_html_path='animation.html',
    ...     timestep='all'
    ... )

    # Create static map for timestep 0
    >>> Map_h5_driver(
    ...     shpfile_fullpath='/path/to/river.shp',
    ...     openwq_results=results,
    ...     hydromodel_out_fullpath='/path/to/mizuroute.nc',
    ...     output_html_path='snapshot.html',
    ...     timestep=0
    ... )

    # Create static map for specific date
    >>> Map_h5_driver(
    ...     shpfile_fullpath='/path/to/river.shp',
    ...     openwq_results=results,
    ...     hydromodel_out_fullpath='/path/to/mizuroute.nc',
    ...     output_html_path='snapshot_date.html',
    ...     timestep='2020-07-15'
    ... )

    # Create animated map for timestep range (by index)
    >>> Map_h5_driver(
    ...     shpfile_fullpath='/path/to/river.shp',
    ...     openwq_results=results,
    ...     hydromodel_out_fullpath='/path/to/mizuroute.nc',
    ...     output_html_path='animation_range.html',
    ...     timestep=[0, 100]  # First 101 timesteps
    ... )

    # Create animated map for date range
    >>> Map_h5_driver(
    ...     shpfile_fullpath='/path/to/river.shp',
    ...     openwq_results=results,
    ...     hydromodel_out_fullpath='/path/to/mizuroute.nc',
    ...     output_html_path='animation_dates.html',
    ...     timestep=['2020-01-01', '2020-03-31']  # Q1 2020
    ... )

    # With OpenWQ date format
    >>> Map_h5_driver(
    ...     shpfile_fullpath='/path/to/river.shp',
    ...     openwq_results=results,
    ...     hydromodel_out_fullpath='/path/to/mizuroute.nc',
    ...     output_html_path='animation_openwq_dates.html',
    ...     timestep=['1961JAN01-08:00:00', '1961JAN04-03:00:00']
    ... )
    """

    if model != "mizuroute":
        print("Only Mizuroute is supported for now!")
        return None

    print("=" * 70)
    print("Interactive Plotly Map Generator")
    print("=" * 70)

    # Step 1: Get data from results dictionary
    print("\nStep 1: Extracting data from results dictionary...")

    # Find the species to plot
    if species_name is None:
        # Use first species in dictionary
        species_key = list(openwq_results.keys())[0]
        print(f"  No species specified, using: {species_key}")
    else:
        # Find matching species
        species_key = None
        for key in openwq_results.keys():
            if species_name in key:
                species_key = key
                break

        if species_key is None:
            print(f"✗ Error: Species '{species_name}' not found in results!")
            print(f"  Available species: {list(openwq_results.keys())}")
            return None

        print(f"  Using species: {species_key}")

    # Parse the key to get compartment, chemical, and units
    parts = species_key.split('@')
    compartment = parts[0]
    chem_unit = parts[1].split('#')
    chemical = chem_unit[0]
    units = chem_unit[1]

    # Get data for specified file extension
    data_found = False
    for ext, data_list in openwq_results[species_key]:
        if ext == file_extension:
            if data_list and len(data_list) > 0:
                filename, ttdata, xyz_coords = data_list[0]
                data_found = True
                print(f"  File extension: {file_extension}")
                print(f"  Data shape: {ttdata.shape}")
                print(f"  Time range: {ttdata.index[0]} to {ttdata.index[-1]}")
                break

    if not data_found:
        print(f"✗ Error: No data found for file extension '{file_extension}'!")
        available_exts = [ext for ext, _ in openwq_results[species_key]]
        print(f"  Available extensions: {available_exts}")
        return None

    # Step 2: Load shapefile
    print("\nStep 2: Loading shapefile...")
    geodf = gpd.read_file(shpfile_fullpath, engine='fiona')
    geodf = geodf.to_crs("EPSG:4326")  # WGS84
    comid = np.array(geodf['ComID'])

    # Calculate center if not provided
    if center is None:
        bounds = geodf.total_bounds  # minx, miny, maxx, maxy
        center = {
            'lat': (bounds[1] + bounds[3]) / 2,
            'lon': (bounds[0] + bounds[2]) / 2
        }
        print(f"  Auto-centered at: lat={center['lat']:.4f}, lon={center['lon']:.4f}")

    # Step 3: Load mizuRoute output for mapping
    print("\nStep 3: Loading mizuRoute output...")
    mizuroute_out_file = netCDF4.Dataset(hydromodel_out_fullpath, 'r')
    mizuroute_out_idFeature = np.array(mizuroute_out_file.variables['reachID'])

    # Step 4: Determine which timesteps to plot
    print("\nStep 4: Determining timesteps...")

    if timestep == 'all':
        # Animated map with all timesteps
        print(f"  Creating animated map with {len(ttdata)} timesteps")
        is_animated = True
        timesteps_to_plot = list(range(len(ttdata)))

    elif isinstance(timestep, list):
        # Range of timesteps provided [min, max]
        if len(timestep) != 2:
            print(f"✗ Error: timestep range must have exactly 2 elements [min, max]")
            print(f"  Got: {timestep}")
            mizuroute_out_file.close()
            return None

        is_animated = True
        min_val, max_val = timestep[0], timestep[1]

        if isinstance(min_val, int) and isinstance(max_val, int):
            # Integer indices provided [0, 100]
            if min_val < 0 or max_val >= len(ttdata) or min_val > max_val:
                print(f"✗ Error: Invalid index range [{min_val}, {max_val}]")
                print(f"  Valid range: [0, {len(ttdata) - 1}]")
                mizuroute_out_file.close()
                return None

            timesteps_to_plot = list(range(min_val, max_val + 1))
            print(f"  Creating animated map for timesteps {min_val} to {max_val}")
            print(f"  Date range: {ttdata.index[min_val]} to {ttdata.index[max_val]}")
            print(f"  Total frames: {len(timesteps_to_plot)}")

        else:
            # Date strings provided ['2020-01-01', '2020-12-31']
            try:
                # Parse dates - handle OpenWQ format like '1961JAN01-08:00:00'
                def parse_openwq_date(date_str):
                    """Parse OpenWQ date format or standard datetime"""
                    try:
                        # Try pandas parsing first
                        return pd.to_datetime(date_str)
                    except:
                        # Try OpenWQ format: YYYYMMMDD-HH:MM:SS
                        try:
                            from datetime import datetime
                            return datetime.strptime(date_str, '%Y%b%d-%H:%M:%S')
                        except:
                            raise ValueError(f"Could not parse date: {date_str}")

                min_date = parse_openwq_date(str(min_val))
                max_date = parse_openwq_date(str(max_val))

                # Find timesteps in this date range
                timesteps_to_plot = []
                for i, ts in enumerate(ttdata.index):
                    ts_datetime = pd.to_datetime(ts)
                    if min_date <= ts_datetime <= max_date:
                        timesteps_to_plot.append(i)

                if len(timesteps_to_plot) == 0:
                    print(f"✗ Error: No timesteps found in date range")
                    print(f"  Requested: {min_val} to {max_val}")
                    print(f"  Available: {ttdata.index[0]} to {ttdata.index[-1]}")
                    mizuroute_out_file.close()
                    return None

                print(f"  Creating animated map for date range:")
                print(f"    From: {min_val}")
                print(f"    To: {max_val}")
                print(f"  Actual range: {ttdata.index[timesteps_to_plot[0]]} to {ttdata.index[timesteps_to_plot[-1]]}")
                print(f"  Total frames: {len(timesteps_to_plot)}")

            except Exception as e:
                print(f"✗ Error parsing date range: {e}")
                print(f"  Available dates: {ttdata.index[0]} to {ttdata.index[-1]}")
                mizuroute_out_file.close()
                return None

    else:
        # Static map for single timestep
        is_animated = False

        if isinstance(timestep, int):
            # Index provided
            if timestep < 0 or timestep >= len(ttdata):
                print(f"✗ Error: timestep {timestep} out of range [0, {len(ttdata) - 1}]")
                mizuroute_out_file.close()
                return None
            timestep_idx = timestep
            timestep_label = str(ttdata.index[timestep_idx])
            print(f"  Creating static map for timestep {timestep_idx}: {timestep_label}")
        else:
            # Date string provided
            try:
                timestep_idx = ttdata.index.get_loc(timestep)
                timestep_label = timestep
                print(f"  Creating static map for date: {timestep}")
            except KeyError:
                # Try finding nearest date
                try:
                    target_date = pd.to_datetime(timestep)
                    nearest_idx = ttdata.index.get_indexer([target_date], method='nearest')[0]
                    timestep_idx = nearest_idx
                    timestep_label = str(ttdata.index[timestep_idx])
                    print(f"  Date '{timestep}' not found. Using nearest: {timestep_label}")
                except:
                    print(f"✗ Error: Could not parse timestep '{timestep}'")
                    print(f"  Available dates: {ttdata.index[0]} to {ttdata.index[-1]}")
                    mizuroute_out_file.close()
                    return None

        timesteps_to_plot = [timestep_idx]

    # Step 5: Process data and create map
    print("\nStep 5: Processing data...")

    # Determine color scale range
    if vmin is None or vmax is None:
        # Sample data to get range
        sample_values = []
        for idx in timesteps_to_plot[::max(1, len(timesteps_to_plot) // 10)]:
            sample_values.extend(ttdata.iloc[idx, :].values)

        if vmin is None:
            vmin = np.nanmin(sample_values)
        if vmax is None:
            vmax = np.nanmax(sample_values)

        print(f"  Auto color range: {vmin:.3e} to {vmax:.3e}")

    if is_animated:
        # Create animated map
        print("\nStep 6: Creating animated map...")
        fig = _create_animated_plotly_map(
            ttdata, xyz_coords, geodf, comid, mizuroute_out_idFeature,
            compartment, chemical, units, file_extension,
            vmin, vmax, colorscale, mapbox_style, zoom, center,
            animation_frame_duration, timesteps_to_plot
        )
    else:
        # Create static map
        print("\nStep 6: Creating static map...")
        fig = _create_static_plotly_map(
            ttdata, xyz_coords, geodf, comid, mizuroute_out_idFeature,
            timestep_idx, timestep_label,
            compartment, chemical, units, file_extension,
            vmin, vmax, colorscale, mapbox_style, zoom, center
        )

    # Close mizuRoute file
    mizuroute_out_file.close()

    # Step 7: Save HTML
    print(f"\nStep 7: Saving to {output_html_path}...")
    fig.write_html(
        output_html_path,
        include_plotlyjs='cdn',
        config={
            'displayModeBar': True,
            'displaylogo': False,
            'modeBarButtonsToRemove': ['lasso2d', 'select2d']
        }
    )

    file_size = os.path.getsize(output_html_path) / 1024 / 1024  # MB

    print("=" * 70)
    print("✓ Interactive map created successfully!")
    print("=" * 70)
    print(f"  Output file: {output_html_path}")
    print(f"  File size: {file_size:.2f} MB")
    if is_animated:
        print(f"  Frames: {len(timesteps_to_plot)}")
        print(f"  Animation duration: {animation_frame_duration} ms per frame")
    else:
        print(f"  Timestep: {timestep_label}")
    print(f"\n  Open in browser to view interactive map!")
    print("=" * 70)

    return output_html_path


def _create_static_plotly_map(ttdata, xyz_coords, geodf, comid, mizuroute_out_idFeature,
                              timestep_idx, timestep_label,
                              compartment, chemical, units, file_extension,
                              vmin, vmax, colorscale, mapbox_style, zoom, center):
    """Create static Plotly map for single timestep"""

    # Get data for this timestep
    data_values = ttdata.iloc[timestep_idx, :].values

    # Map to geodataframe
    geodf_openwq_wq = _map_data_to_geodf(
        data_values, xyz_coords, geodf, comid, mizuroute_out_idFeature
    )

    # Convert to GeoJSON
    geodf_json = json.loads(geodf.to_json())

    # Create choropleth map
    fig = go.Figure(go.Choroplethmapbox(
        geojson=geodf_json,
        locations=geodf.index,
        z=geodf_openwq_wq,
        colorscale=colorscale,
        zmin=vmin,
        zmax=vmax,
        marker_opacity=0.7,
        marker_line_width=0.5,
        marker_line_color='white',
        colorbar=dict(
            title=f"{units}",
            thickness=15,
            len=0.7,
            x=1.02
        ),
        hovertemplate='<b>Reach ID</b>: %{location}<br>' +
                      f'<b>{chemical}</b>: ' + '%{z:.3e}<br>' +
                      '<extra></extra>'
    ))

    # Update layout
    fig.update_layout(
        mapbox_style=mapbox_style,
        mapbox=dict(
            center=center,
            zoom=zoom
        ),
        title={
            'text': f'{chemical} in {compartment}<br>{file_extension} - {timestep_label}',
            'x': 0.5,
            'xanchor': 'center',
            'font': {'size': 20}
        },
        height=800,
        margin={"r": 50, "t": 100, "l": 50, "b": 50}
    )

    return fig


def _create_animated_plotly_map(ttdata, xyz_coords, geodf, comid, mizuroute_out_idFeature,
                                compartment, chemical, units, file_extension,
                                vmin, vmax, colorscale, mapbox_style, zoom, center,
                                animation_frame_duration, timesteps_to_plot):
    """Create animated Plotly map with specified timesteps"""

    # Convert geodataframe to GeoJSON (once)
    geodf_json = json.loads(geodf.to_json())

    # Prepare frames for animation
    frames = []

    for i in timesteps_to_plot:
        # Get data for this timestep
        data_values = ttdata.iloc[i, :].values

        # Map to geodataframe
        geodf_openwq_wq = _map_data_to_geodf(
            data_values, xyz_coords, geodf, comid, mizuroute_out_idFeature
        )

        # Create frame
        frame = go.Frame(
            data=[go.Choroplethmapbox(
                geojson=geodf_json,
                locations=geodf.index,
                z=geodf_openwq_wq,
                colorscale=colorscale,
                zmin=vmin,
                zmax=vmax,
                marker_opacity=0.7,
                marker_line_width=0.5,
                marker_line_color='white',
                colorbar=dict(
                    title=f"{units}",
                    thickness=15,
                    len=0.7,
                    x=1.02
                ),
                hovertemplate='<b>Reach ID</b>: %{location}<br>' +
                              f'<b>{chemical}</b>: ' + '%{z:.3e}<br>' +
                              '<extra></extra>'
            )],
            name=str(ttdata.index[i]),
            layout=go.Layout(
                title_text=f'{chemical} in {compartment} - {file_extension}<br>Time: {ttdata.index[i]}'
            )
        )
        frames.append(frame)

        if (len(frames)) % 10 == 0:
            print(f"  Processed {len(frames)}/{len(timesteps_to_plot)} timesteps...")

    # Create initial figure with first timestep in range
    data_values_0 = ttdata.iloc[timesteps_to_plot[0], :].values
    geodf_openwq_wq_0 = _map_data_to_geodf(
        data_values_0, xyz_coords, geodf, comid, mizuroute_out_idFeature
    )

    fig = go.Figure(
        data=[go.Choroplethmapbox(
            geojson=geodf_json,
            locations=geodf.index,
            z=geodf_openwq_wq_0,
            colorscale=colorscale,
            zmin=vmin,
            zmax=vmax,
            marker_opacity=0.7,
            marker_line_width=0.5,
            marker_line_color='white',
            colorbar=dict(
                title=f"{units}",
                thickness=15,
                len=0.7,
                x=1.02
            ),
            hovertemplate='<b>Reach ID</b>: %{location}<br>' +
                          f'<b>{chemical}</b>: ' + '%{z:.3e}<br>' +
                          '<extra></extra>'
        )],
        frames=frames
    )

    # Add slider and play/pause buttons
    sliders = [dict(
        active=0,
        yanchor="top",
        y=0.99,
        xanchor="left",
        x=0.01,
        currentvalue=dict(
            prefix="Time: ",
            visible=True,
            xanchor="right"
        ),
        pad=dict(b=10, t=50),
        len=0.9,
        steps=[
            dict(
                args=[
                    [frame.name],
                    dict(
                        frame=dict(duration=animation_frame_duration, redraw=True),
                        mode="immediate",
                        transition=dict(duration=animation_frame_duration // 2)
                    )
                ],
                method="animate",
                label=frame.name
            )
            for frame in frames
        ]
    )]

    # Add play/pause buttons
    updatemenus = [dict(
        type="buttons",
        direction="left",
        x=0.01,
        y=0.01,
        xanchor="left",
        yanchor="bottom",
        showactive=False,
        buttons=[
            dict(
                label="▶ Play",
                method="animate",
                args=[
                    None,
                    dict(
                        frame=dict(duration=animation_frame_duration, redraw=True),
                        fromcurrent=True,
                        transition=dict(duration=animation_frame_duration // 2)
                    )
                ]
            ),
            dict(
                label="⏸ Pause",
                method="animate",
                args=[
                    [None],
                    dict(
                        frame=dict(duration=0, redraw=False),
                        mode="immediate",
                        transition=dict(duration=0)
                    )
                ]
            )
        ]
    )]

    # Update layout
    fig.update_layout(
        mapbox_style=mapbox_style,
        mapbox=dict(
            center=center,
            zoom=zoom
        ),
        title={
            'text': f'{chemical} in {compartment} - {file_extension}<br>Time: {ttdata.index[timesteps_to_plot[0]]}',
            'x': 0.5,
            'xanchor': 'center',
            'font': {'size': 20}
        },
        height=800,
        margin={"r": 50, "t": 120, "l": 50, "b": 80},
        sliders=sliders,
        updatemenus=updatemenus
    )

    return fig


def _map_data_to_geodf(data_values, xyz_coords, geodf, comid, mizuroute_out_idFeature):
    """Map OpenWQ data to geodataframe based on coordinates"""

    geodf_openwq_wq = np.zeros(len(comid))

    # Map data based on coordinates
    for j in range(len(comid)):
        comid_j = comid[j]
        index_j = np.where(mizuroute_out_idFeature == comid_j)[0]

        if len(index_j) > 0:
            # Find matching data point by coordinates
            # For now, just use first match
            # TODO: Could improve by matching xyz coordinates
            try:
                geodf_openwq_wq[j] = data_values[index_j[0]]
            except:
                geodf_openwq_wq[j] = np.nan
        else:
            geodf_openwq_wq[j] = np.nan

    return geodf_openwq_wq


if __name__ == "__main__":
    print("=" * 70)
    print("OpenWQ Plotly Interactive Mapping")
    print("Using results dictionary from Read_h5_driver")
    print("=" * 70)

    print("\nExample workflow:")
    print("""
# Step 1: Read OpenWQ data
from Read_h5_optimized import Read_h5_driver

results = Read_h5_driver(
    folderpath='/path/to/output/',
    OpenWQ_output_format='HDF5',
    Compartments=['RIVER_NETWORK_REACHES'],
    Chemical_species=['SPECIES_A', 'SPECIES_B'],
    Concentration_units='MG/L'
)

# Step 2: Create animated map (all timesteps)
from Map_h5_driver import Map_h5_driver

Map_h5_driver(
    shpfile_fullpath='/path/to/river.shp',
    openwq_results=results,
    hydromodel_out_fullpath='/path/to/mizuroute.nc',
    output_html_path='animation.html',
    timestep='all'  # Animated with play button
)

# Step 3: Create static map for specific timestep
Map_h5_driver(
    shpfile_fullpath='/path/to/river.shp',
    openwq_results=results,
    hydromodel_out_fullpath='/path/to/mizuroute.nc',
    output_html_path='snapshot.html',
    timestep=0  # First timestep
)

# Step 4: Create static map for specific date
Map_h5_driver(
    shpfile_fullpath='/path/to/river.shp',
    openwq_results=results,
    hydromodel_out_fullpath='/path/to/mizuroute.nc',
    output_html_path='snapshot_date.html',
    timestep='2020-07-15'  # Specific date
)

# Step 5: Create animated map for timestep range (indices)
Map_h5_driver(
    shpfile_fullpath='/path/to/river.shp',
    openwq_results=results,
    hydromodel_out_fullpath='/path/to/mizuroute.nc',
    output_html_path='animation_range.html',
    timestep=[0, 100]  # First 101 timesteps
)

# Step 6: Create animated map for date range
Map_h5_driver(
    shpfile_fullpath='/path/to/river.shp',
    openwq_results=results,
    hydromodel_out_fullpath='/path/to/mizuroute.nc',
    output_html_path='animation_dates.html',
    timestep=['2020-01-01', '2020-12-31']  # Entire year 2020
)

# Step 7: With OpenWQ date format
Map_h5_driver(
    shpfile_fullpath='/path/to/river.shp',
    openwq_results=results,
    hydromodel_out_fullpath='/path/to/mizuroute.nc',
    output_html_path='animation_openwq.html',
    timestep=['1961JAN01-08:00:00', '1961JAN04-03:00:00']
)
    """)