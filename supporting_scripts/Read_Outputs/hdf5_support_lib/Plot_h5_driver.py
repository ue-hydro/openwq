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
Plotly visualization function for OpenWQ model results dictionary
Supports both temporal (time series) and spatial (snapshot) plotting
"""

import plotly.graph_objects as go
import plotly.express as px
from pathlib import Path
import numpy as np


def Plot_h5_driver(results_dict,
                   plot_type='temporal',
                   output_dir='plots',
                   debug_mode=None,
                   # Temporal plot parameters
                   cells_to_plot=None,
                   # Spatial plot parameters
                   time_slice=0,
                   spatial_axis='x',
                   fixed_coords=None,
                   # Common parameters
                   colors=None,
                   line_width=2,
                   template="plotly_white",
                   height=600,
                   save_html=True,
                   show_plots=False):
    """
    Create Plotly plots for OpenWQ results - supports temporal and spatial visualization.

    Parameters
    ----------
    results_dict : dict
        Dictionary from Read_h5_driver

    plot_type : str, optional
        Type of plot to create:
        - 'temporal': Time series plots (default)
        - 'spatial': Spatial distribution at specific time

    output_dir : str, optional
        Directory to save plots (default: 'plots')

    debug_mode : bool or None, optional
        If True: plot all file extensions
        If False: plot only 'main'
        If None: plot all available extensions

    Temporal Plot Parameters
    ------------------------
    cells_to_plot : list of int, list of tuples, or 'all', optional
        Which cells to plot in time series:
        - None or 'all': plot all cells
        - list of int: plot cells by index, e.g., [0, 5, 10]
        - list of tuples: plot cells by (x,y,z) coords, e.g., [(10,0,0), (20,0,0)]

    Spatial Plot Parameters
    -----------------------
    time_slice : int or str, optional
        Which timestep to plot (default: 0 = first timestep)
        - int: index of timestep
        - str: date string matching index, e.g., '2020-01-15'

    spatial_axis : str, optional
        Which spatial dimension to plot on x-axis:
        - 'x': plot along x-axis (default)
        - 'y': plot along y-axis
        - 'z': plot along z-axis

    fixed_coords : dict, optional
        Fix certain coordinates for spatial plot. Examples:
        - {'y': 0, 'z': 0}: fix y and z, plot along x
        - {'x': 100, 'z': 0}: fix x and z, plot along y
        - {'x': 100, 'y': 50}: fix x and y, plot along z
        - None: plot all cells along specified axis

    Common Parameters
    -----------------
    colors : list of str, optional
        Custom colors for lines

    line_width : float, optional
        Width of lines (default: 2)

    template : str, optional
        Plotly template (default: "plotly_white")

    height : int, optional
        Plot height in pixels (default: 600)

    save_html : bool, optional
        Whether to save plots as HTML files (default: True)

    show_plots : bool, optional
        Whether to display plots (default: False)

    Returns
    -------
    dict
        Dictionary with keys as plot names and values as plotly figures

    Examples
    --------
    # Example 1: Temporal plot - time series for specific cells
    >>> figs = Plot_h5_driver(
    ...     results,
    ...     plot_type='temporal',
    ...     cells_to_plot=[0, 100, 500]
    ... )

    # Example 2: Temporal plot - time series for specific coordinates
    >>> figs = Plot_h5_driver(
    ...     results,
    ...     plot_type='temporal',
    ...     cells_to_plot=[(10, 0, 0), (50, 0, 0), (100, 0, 0)]
    ... )

    # Example 3: Spatial plot - concentration along x-axis at first timestep
    >>> figs = Plot_h5_driver(
    ...     results,
    ...     plot_type='spatial',
    ...     time_slice=0,
    ...     spatial_axis='x',
    ...     fixed_coords={'y': 0, 'z': 0}
    ... )

    # Example 4: Spatial plot - concentration along river at specific date
    >>> figs = Plot_h5_driver(
    ...     results,
    ...     plot_type='spatial',
    ...     time_slice='2020-07-15',
    ...     spatial_axis='x',
    ...     fixed_coords={'y': 0, 'z': 0}
    ... )

    # Example 5: Spatial plot - all y values at x=100, z=0
    >>> figs = Plot_h5_driver(
    ...     results,
    ...     plot_type='spatial',
    ...     time_slice=0,
    ...     spatial_axis='y',
    ...     fixed_coords={'x': 100, 'z': 0}
    ... )
    """

    # Validate plot_type
    if plot_type not in ['temporal', 'spatial']:
        raise ValueError("plot_type must be 'temporal' or 'spatial'")

    # Determine which file extensions to plot
    if debug_mode is True:
        file_extensions_to_plot = [
            'main',
            'd_output_dt_chemistry',
            'd_output_dt_transport',
            'd_output_ss',
            'd_output_ewf',
            'd_output_ic'
        ]
    elif debug_mode is False:
        file_extensions_to_plot = ['main']
    else:
        file_extensions_to_plot = None  # Plot all available

    # Create output directory
    if save_html:
        output_path = Path(output_dir)
        output_path.mkdir(exist_ok=True)

    # Store all figures
    all_figures = {}

    # Set default colors
    if colors is None:
        colors = px.colors.qualitative.Plotly

    # Counter for progress
    total_plots = 0
    created_plots = 0

    # Count total plots to create
    for key in results_dict.keys():
        for file_ext, data_list in results_dict[key]:
            if file_extensions_to_plot is None or file_ext in file_extensions_to_plot:
                total_plots += 1

    print(f"\nCreating {total_plots} {plot_type} plot(s)...")
    print("=" * 70)

    # Loop through each entry in the dictionary
    for key in results_dict.keys():

        # Parse the key to get compartment, chemical, and units
        parts = key.split('@')
        compartment = parts[0]
        chem_unit = parts[1].split('#')
        chemical = chem_unit[0]
        units = chem_unit[1]

        # Loop through each file extension
        for file_ext, data_list in results_dict[key]:

            # Skip if not in requested file extensions
            if file_extensions_to_plot is not None and file_ext not in file_extensions_to_plot:
                continue

            # Check if data exists
            if not data_list or len(data_list) == 0:
                print(f"⚠ Skipping {compartment}@{chemical} ({file_ext}): No data")
                continue

            # Extract data
            filename, ttdata, xyz_coords = data_list[0]

            # Create appropriate plot based on plot_type
            if plot_type == 'temporal':
                fig = _create_temporal_plot(
                    ttdata, xyz_coords, cells_to_plot,
                    compartment, chemical, units, file_ext,
                    colors, line_width, template, height
                )
                plot_suffix = "temporal"

            else:  # spatial
                fig = _create_spatial_plot(
                    ttdata, xyz_coords, time_slice, spatial_axis, fixed_coords,
                    compartment, chemical, units, file_ext,
                    colors, line_width, template, height
                )
                time_label = time_slice if isinstance(time_slice, str) else f"t{time_slice}"
                plot_suffix = f"spatial_{spatial_axis}_{time_label}"

            if fig is None:
                continue

            # Create filename for saving
            safe_compartment = compartment.replace('/', '_')
            safe_chemical = chemical.replace('/', '_')
            plot_name = f"{safe_compartment}_{safe_chemical}_{file_ext}_{plot_suffix}"

            # Save figure
            if save_html:
                output_file = output_path / f"{plot_name}.html"
                fig.write_html(str(output_file))
                print(f"✓ Created: {output_file.name}")

            # Show plot if requested
            if show_plots:
                fig.show()

            # Store figure
            all_figures[plot_name] = fig
            created_plots += 1

    print("=" * 70)
    print(f"✓ Successfully created {created_plots} plot(s)")
    if save_html:
        print(f"✓ Saved to directory: {output_path.absolute()}")

    return all_figures


def _create_temporal_plot(ttdata, xyz_coords, cells_to_plot,
                          compartment, chemical, units, file_ext,
                          colors, line_width, template, height):
    """Create temporal (time series) plot"""

    # Determine which cells to plot
    if cells_to_plot is None or cells_to_plot == 'all':
        selected_indices = list(range(len(ttdata.columns)))

    elif isinstance(cells_to_plot[0], (tuple, list)):
        # User provided (x,y,z) coordinates
        selected_indices = []
        for target_coord in cells_to_plot:
            # Find matching cell
            for i, coord in enumerate(xyz_coords):
                if (coord[0] == target_coord[0] and
                        coord[1] == target_coord[1] and
                        coord[2] == target_coord[2]):
                    selected_indices.append(i)
                    break
    else:
        # User provided indices
        selected_indices = cells_to_plot

    if len(selected_indices) == 0:
        print(f"⚠ Warning: No matching cells found")
        return None

    # Create figure
    fig = go.Figure()

    # Get x-axis (time)
    x_values = ttdata.index

    # Add traces (one per cell)
    for i, col_idx in enumerate(selected_indices):
        col_name = ttdata.columns[col_idx]
        y_values = ttdata[col_name].values

        # Create line name with coordinates
        if col_idx < len(xyz_coords):
            coords = xyz_coords[col_idx]
            line_name = f"Cell ({coords[0]:.0f}, {coords[1]:.0f}, {coords[2]:.0f})"
        else:
            line_name = col_name

        # Get color (cycle if needed)
        color = colors[i % len(colors)]

        fig.add_trace(go.Scatter(
            x=x_values,
            y=y_values,
            mode='lines',
            name=line_name,
            line=dict(color=color, width=line_width),
            hovertemplate='%{y:.3e}<extra></extra>'
        ))

    # Update layout
    title = f"{chemical} in {compartment} - {file_ext} (Temporal)"

    fig.update_layout(
        title=dict(text=title, x=0.5, xanchor='center'),
        xaxis_title="Time",
        yaxis_title=f"Concentration ({units})",
        template=template,
        showlegend=True,
        hovermode='x unified',
        height=height,
        legend=dict(
            yanchor="top",
            y=0.99,
            xanchor="left",
            x=0.01
        )
    )

    # Add grid
    fig.update_xaxes(showgrid=True, gridwidth=1, gridcolor='LightGray')
    fig.update_yaxes(showgrid=True, gridwidth=1, gridcolor='LightGray')

    return fig


def _create_spatial_plot(ttdata, xyz_coords, time_slice, spatial_axis, fixed_coords,
                         compartment, chemical, units, file_ext,
                         colors, line_width, template, height):
    """Create spatial (snapshot) plot"""

    # Get the time slice
    if isinstance(time_slice, str):
        # Time slice is a date string
        time_index = ttdata.index.get_loc(time_slice)
        time_label = time_slice
    else:
        # Time slice is an index
        time_index = time_slice
        time_label = str(ttdata.index[time_index])

    # Get data for this timestep
    data_snapshot = ttdata.iloc[time_index, :]

    # Determine which dimension is the spatial axis
    axis_map = {'x': 0, 'y': 1, 'z': 2}
    axis_idx = axis_map[spatial_axis.lower()]

    # Filter cells based on fixed_coords
    if fixed_coords is None:
        # Plot all cells
        selected_cells = list(range(len(xyz_coords)))
    else:
        # Filter cells that match fixed coordinates
        selected_cells = []
        for i, coord in enumerate(xyz_coords):
            match = True
            for axis, value in fixed_coords.items():
                axis_i = axis_map[axis.lower()]
                if coord[axis_i] != value:
                    match = False
                    break
            if match:
                selected_cells.append(i)

    if len(selected_cells) == 0:
        print(f"⚠ Warning: No cells match fixed_coords criteria")
        return None

    # Extract spatial coordinates and values
    spatial_coords = [xyz_coords[i][axis_idx] for i in selected_cells]
    values = [data_snapshot.iloc[i] for i in selected_cells]

    # Sort by spatial coordinate
    sorted_indices = np.argsort(spatial_coords)
    spatial_coords_sorted = [spatial_coords[i] for i in sorted_indices]
    values_sorted = [values[i] for i in sorted_indices]

    # Create figure
    fig = go.Figure()

    fig.add_trace(go.Scatter(
        x=spatial_coords_sorted,
        y=values_sorted,
        mode='lines+markers',
        name='Concentration',
        line=dict(color=colors[0], width=line_width),
        marker=dict(size=6, color=colors[0]),
        hovertemplate='%{x:.1f}: %{y:.3e}<extra></extra>'
    ))

    # Create title with fixed coordinates info
    fixed_info = ""
    if fixed_coords:
        fixed_parts = [f"{k}={v}" for k, v in fixed_coords.items()]
        fixed_info = f" (fixed: {', '.join(fixed_parts)})"

    title = f"{chemical} in {compartment} - {file_ext}<br>Spatial at {time_label}{fixed_info}"

    fig.update_layout(
        title=dict(text=title, x=0.5, xanchor='center'),
        xaxis_title=f"{spatial_axis.upper()}-coordinate",
        yaxis_title=f"Concentration ({units})",
        template=template,
        showlegend=False,
        height=height
    )

    # Add grid
    fig.update_xaxes(showgrid=True, gridwidth=1, gridcolor='LightGray')
    fig.update_yaxes(showgrid=True, gridwidth=1, gridcolor='LightGray')

    return fig

