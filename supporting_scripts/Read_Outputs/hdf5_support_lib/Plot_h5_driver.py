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
Plotly visualization function for OpenWQ model results dictionary
"""

import plotly.graph_objects as go
import plotly.express as px
from pathlib import Path


def Plot_h5_driver(results_dict,
                    output_dir='plots',
                    debug_mode=None,
                    cells_to_plot=None,
                    colors=None,
                    line_width=2,
                    template="plotly_white",
                    height=600,
                    save_html=True,
                    show_plots=False):
    """
    Create Plotly line plots for all entries in OpenWQ results dictionary.

    The function automatically creates one plot for each file extension
    (main, d_output_dt_chemistry, etc.) for each dictionary entry
    (e.g., SPECIES_A, SPECIES_B).

    Parameters
    ----------
    results_dict : dict
        Dictionary from Read_h5_driver with structure:
        {
            'COMPARTMENT@CHEMICAL#UNITS': [
                ('main', [(filename, DataFrame, xyz_coords)]),
                ('d_output_dt_chemistry', [(filename, DataFrame, xyz_coords)]),
                ...
            ]
        }
    output_dir : str, optional
        Directory to save plots (default: 'plots')
    file_extensions_to_plot : list of str, optional
        Which file extensions to plot. If None, plots all.
        Options: 'main', 'd_output_dt_chemistry', 'd_output_dt_transport',
                 'd_output_ss', 'd_output_ewf', 'd_output_ic'
    cells_to_plot : list of int or 'all', optional
        Which cell columns to plot. If None or 'all', plots all cells.
        Can be list of column indices, e.g., [0, 5, 10]
    colors : list of str, optional
        Custom colors for lines. If None, uses Plotly default palette
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
    >>> from Read_h5_optimized import Read_h5_driver
    >>> from plot_openwq_results import plot_openwq_results
    >>>
    >>> # Read data
    >>> results = Read_h5_driver(...)
    >>>
    >>> # Create all plots
    >>> figs = plot_openwq_results(results)
    >>>
    >>> # Create plots only for 'main' extension
    >>> figs = plot_openwq_results(
    ...     results,
    ...     file_extensions_to_plot=['main']
    ... )
    >>>
    >>> # Plot only specific cells
    >>> figs = plot_openwq_results(
    ...     results,
    ...     cells_to_plot=[0, 5, 10],
    ...     colors=['red', 'blue', 'green']
    ... )
    """

    if debug_mode:
        file_extensions_to_plot = [
            'main',
            'd_output_dt_chemistry',
            'd_output_dt_transport',
            'd_output_ss',
            'd_output_ewf',
            'd_output_ic'
        ]
    else:
        file_extensions_to_plot = ["main"]

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

    print(f"\nCreating {total_plots} plot(s)...")
    print("=" * 70)

    # Loop through each entry in the dictionary
    for key in results_dict.keys():

        # Parse the key to get compartment, chemical, and units
        # Format: 'COMPARTMENT@CHEMICAL#UNITS'
        parts = key.split('@')
        compartment = parts[0]
        chem_unit = parts[1].split('#')
        chemical = chem_unit[0]
        units = chem_unit[1]

        # Loop through each file extension (main, d_output_dt_chemistry, etc.)
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

            # Determine which cells to plot
            if cells_to_plot is None or cells_to_plot == 'all':
                data_to_plot = ttdata
                selected_indices = list(range(len(ttdata.columns)))
            else:
                data_to_plot = ttdata.iloc[:, cells_to_plot]
                selected_indices = cells_to_plot

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
            title = f"{chemical} in {compartment} - {file_ext}"

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

            # Create filename for saving
            safe_compartment = compartment.replace('/', '_')
            safe_chemical = chemical.replace('/', '_')
            plot_name = f"{safe_compartment}_{safe_chemical}_{file_ext}"

            # Save figure
            if save_html:
                output_file = output_path / f"{plot_name}.html"
                fig.write_html(str(output_file))
                print(f"✓ Created: {output_file.name} ({len(selected_indices)} cells)")

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
