"""
OpenWQ Simple Static Maps with Automatic GIF Creation
UPDATED: Supports both mizuRoute (polylines) and SUMMA (polygons)
Can map OpenWQ results AND hostmodel outputs
"""

import numpy as np
import pandas as pd
import geopandas as gpd
import netCDF4
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Polygon as MplPolygon
import os
import subprocess


def _map_data_to_geodf(data_values, xyz_coords, geodf, shpf_mapKey, data_hostmdlKeys):
    """Map OpenWQ data to geodataframe - works for both lines and polygons"""
    geodf_openwq_wq = np.zeros(len(shpf_mapKey))
    data_hostmdlKeys_mainInfo = [item.split('_')[1] for item in data_hostmdlKeys]

    for j in range(len(shpf_mapKey)):
        mapKey_j = shpf_mapKey[j]
        index_j = np.where(np.array(data_hostmdlKeys_mainInfo) == str(mapKey_j))[0]

        if len(index_j) > 0:
            try:
                geodf_openwq_wq[j] = data_values[index_j[0]]
            except:
                geodf_openwq_wq[j] = np.nan
        else:
            geodf_openwq_wq[j] = np.nan

    return geodf_openwq_wq


def _create_gif_from_pngs(image_dir, output_gif, duration_ms=500):
    """Try multiple methods to create GIF from PNG files"""

    print(f"\n  Attempting to create GIF: {output_gif}")
    print(f"  Frame duration: {duration_ms}ms")

    png_files = sorted([f for f in os.listdir(image_dir) if f.endswith('.png')])

    if len(png_files) == 0:
        print("  ✗ No PNG files found")
        return False

    print(f"  Found {len(png_files)} PNG files")

    # Method 1: Try imageio
    try:
        import imageio.v2 as imageio
        print("  → Trying imageio...")

        images = []
        for filename in png_files:
            filepath = os.path.join(image_dir, filename)
            images.append(imageio.imread(filepath))

        imageio.mimsave(output_gif, images, duration=duration_ms / 1000.0, loop=0)

        size_mb = os.path.getsize(output_gif) / (1024 * 1024)
        print(f"  ✓ SUCCESS with imageio!")
        print(f"    File: {output_gif}")
        print(f"    Size: {size_mb:.1f} MB")
        return True

    except ImportError:
        print("  ✗ imageio not available")
    except Exception as e:
        print(f"  ✗ imageio failed: {e}")

    # Method 2: Try Pillow/PIL
    try:
        from PIL import Image
        print("  → Trying Pillow/PIL...")

        image_files = [os.path.join(image_dir, f) for f in png_files]
        images = [Image.open(f) for f in image_files]

        images[0].save(
            output_gif,
            save_all=True,
            append_images=images[1:],
            duration=duration_ms,
            loop=0
        )

        size_mb = os.path.getsize(output_gif) / (1024 * 1024)
        print(f"  ✓ SUCCESS with Pillow!")
        print(f"    File: {output_gif}")
        print(f"    Size: {size_mb:.1f} MB")
        return True

    except ImportError:
        print("  ✗ Pillow not available")
    except Exception as e:
        print(f"  ✗ Pillow failed: {e}")

    # Method 3: Try ImageMagick
    try:
        print("  → Trying ImageMagick (convert)...")

        delay = int(duration_ms / 10)

        cmd = [
            'convert',
            '-delay', str(delay),
            '-loop', '0',
            os.path.join(image_dir, 'map_*.png'),
            output_gif
        ]

        result = subprocess.run(cmd, capture_output=True, text=True, check=True)

        size_mb = os.path.getsize(output_gif) / (1024 * 1024)
        print(f"  ✓ SUCCESS with ImageMagick!")
        print(f"    File: {output_gif}")
        print(f"    Size: {size_mb:.1f} MB")
        return True

    except FileNotFoundError:
        print("  ✗ ImageMagick not installed")
    except subprocess.CalledProcessError as e:
        print(f"  ✗ ImageMagick failed: {e}")
    except Exception as e:
        print(f"  ✗ ImageMagick error: {e}")

    # All methods failed
    print("\n  ✗ Could not create GIF - no suitable tool found")
    print("  Install one of:")
    print("    - imageio:      pip install imageio --break-system-packages")
    print("    - Pillow:       pip install Pillow --break-system-packages")
    print("    - ImageMagick:  brew install imagemagick (Mac)")

    return False


def Map_h5_driver(shpfile_info=None,
                  openwq_results=None,
                  hydromodel_info=None,
                  hydromodel_var2print=None,
                  output_html_path=None,
                  chemSpec=None,
                  file_extension='main',
                  timeframes=None,
                  hostmodel=None,
                  what2map=None,
                  create_gif=True,
                  gif_duration=500):
    """
    Create static PNG maps AND animated GIF automatically

    Supports:
    - OpenWQ results (water quality)
    - mizuRoute outputs (polylines - river segments)
    - SUMMA outputs (polygons - HRUs)

    Parameters:
    -----------
    shpfile_info : dict
        Dictionary with 'path_to_shp' and 'mapping_key' keys
    openwq_results : dict
        Results dictionary from Read_h5_driver (required if what2map='openwq')
    hydromodel_info : str or dict
        If what2map='openwq': Path to hydromodel NetCDF output file (string)
        If what2map='hostmodel': Dictionary with keys:
            - 'path_to_shp': Path to NetCDF file
            - 'mapping_key': Variable name for segment/HRU IDs (e.g., 'SegId' or 'hruId')
            - 'var2print': Variable to map (e.g., 'runoff', 'scalarSWE')
    output_html_path : str
        Base path for output
    chemSpec : list or None
        List of chemical species to map (only for what2map='openwq')
    file_extension : str
        File extension identifier (default: 'main')
    timeframes : int or None
        Number of equally spaced frames (default: 50)
    hostmodel : str
        Host model name: 'mizuroute' or 'summa' (default: None)
    what2map : str
        What to map: 'openwq' or 'hostmodel' (default: None)
    create_gif : bool
        Whether to create GIF automatically (default: True)
    gif_duration : int
        Duration per frame in milliseconds (default: 500)

    Returns:
    --------
    dict : {'image_dirs': [...], 'gifs': [...]}

    Examples:
    ---------
    # Map OpenWQ results
    Map_h5_driver(
        shpfile_info={'path_to_shp': 'rivers.shp', 'mapping_key': 'ID'},
        openwq_results=results,
        hydromodel_info='mizuroute.nc',
        output_html_path='output',
        chemSpec=['NO3'],
        what2map='openwq',
        timeframes=50
    )

    # Map mizuRoute runoff (polylines)
    Map_h5_driver(
        shpfile_info={'path_to_shp': 'rivers.shp', 'mapping_key': 'ID'},
        hydromodel_info={
            'path_to_shp': 'mizuroute_output.nc',
            'mapping_key': 'SegId',
            'var2print': 'runoff'
        },
        output_html_path='output',
        hostmodel='mizuroute',
        what2map='hostmodel',
        timeframes=50
    )

    # Map SUMMA SWE (polygons)
    Map_h5_driver(
        shpfile_info={'path_to_shp': 'hrus.shp', 'mapping_key': 'hruId'},
        hydromodel_info={
            'path_to_shp': 'summa_output.nc',
            'mapping_key': 'hruId',
            'var2print': 'scalarSWE'
        },
        output_html_path='output',
        hostmodel='summa',
        what2map='hostmodel',
        timeframes=50
    )
    """

    print("=" * 70)
    print("STATIC PNG MAPS + AUTOMATIC GIF CREATION")
    print(f"Mapping: {what2map.upper()}")
    if what2map == 'hostmodel':
        print(f"Model: {hostmodel.upper()}")
    print("=" * 70)

    # Load shapefile
    print("\n" + "=" * 70)
    print("STEP 1: Loading shapefile...")
    print("=" * 70)

    shpfile_fullpath = shpfile_info["path_to_shp"]
    geodf = gpd.read_file(shpfile_fullpath, engine='fiona')
    geodf = geodf.to_crs("EPSG:4326")

    shapefile_mappingKey = shpfile_info["mapping_key"]
    shpf_mapKey = np.array(geodf[shapefile_mappingKey])

    # Detect geometry type
    geom_types = geodf.geometry.geom_type.unique()
    print(f"  Loaded {len(geodf)} features")
    print(f"  Geometry types: {', '.join(geom_types)}")

    # Determine if we're working with lines or polygons
    is_polygon = any(gt in ['Polygon', 'MultiPolygon'] for gt in geom_types)
    is_line = any(gt in ['LineString', 'MultiLineString'] for gt in geom_types)

    if is_polygon and is_line:
        print("  ⚠ Warning: Mixed geometry types detected")
    elif is_polygon:
        print(f"  → Polygon shapefile (suitable for SUMMA HRUs)")
    elif is_line:
        print(f"  → Polyline shapefile (suitable for mizuRoute segments)")

    # Prepare output directory
    output_dir = os.path.dirname(output_html_path) if output_html_path else '.'

    created_dirs = []
    created_gifs = []

    # BRANCH: Handle openwq vs hostmodel
    if what2map == 'hostmodel':
        print("\n" + "=" * 70)
        print("STEP 2: Loading hostmodel output...")
        print("=" * 70)

        # Extract hostmodel NetCDF info
        if isinstance(hydromodel_info, dict):
            hostmodel_nc_path = hydromodel_info['path_to_shp']
            hostmodel_mapping_key = hydromodel_info['mapping_key']
            var2print = hydromodel_var2print
        else:
            print("✗ Error: For what2map='hostmodel', hydromodel_info must be a dictionary")
            print("  with keys: 'path_to_shp', 'mapping_key', 'var2print'")
            return None

        print(f"  NetCDF file: {hostmodel_nc_path}")
        print(f"  Mapping key: {hostmodel_mapping_key}")
        print(f"  Variable: {var2print}")
        print(f"  Model: {hostmodel}")

        # Open NetCDF file
        nc_file = netCDF4.Dataset(hostmodel_nc_path, 'r')

        # Get segment/HRU IDs and data
        try:
            segment_ids = np.array(nc_file.variables[hostmodel_mapping_key][:])
            data_var = nc_file.variables[var2print]

            print(f"  Data shape: {data_var.shape}")
            print(f"  Features: {len(segment_ids)}")

            # Get time variable
            if 'time' in nc_file.variables:
                time_var = nc_file.variables['time']
                time_data = time_var[:]

                # Try to convert to datetime
                if hasattr(time_var, 'units'):
                    from netCDF4 import num2date
                    dates = num2date(time_data, time_var.units, only_use_cftime_datetimes=False)
                    time_index = pd.DatetimeIndex(dates)
                else:
                    time_index = pd.RangeIndex(len(time_data))
            else:
                time_index = pd.RangeIndex(data_var.shape[0])

            print(f"  Time range: {time_index[0]} to {time_index[-1]}")

            # Create DataFrame for compatibility
            ttdata = pd.DataFrame(
                data_var[:, :],
                index=time_index,
                columns=[f'SEG_{sid}' for sid in segment_ids]
            )

            # Create species_list with single entry
            species_list = [f'{hostmodel}_{var2print}']
            units = nc_file.variables[var2print].units if hasattr(nc_file.variables[var2print], 'units') else 'units'

            print(f"  Units: {units}")

        except KeyError as e:
            print(f"✗ Error: Variable not found in NetCDF: {e}")
            nc_file.close()
            return None

    else:  # what2map == 'openwq'
        print("\n" + "=" * 70)
        print("STEP 2: Processing OpenWQ results...")
        print("=" * 70)

        # Determine which species to process
        if chemSpec is None:
            species_list = list(openwq_results.keys())
            print(f"\nProcessing ALL {len(species_list)} species")
        elif isinstance(chemSpec, list):
            species_list = []
            for chem in chemSpec:
                found = False
                for key in openwq_results.keys():
                    if chem in key:
                        species_list.append(key)
                        found = True
                        break
                if not found:
                    print(f"⚠ Warning: '{chem}' not found, skipping...")

            if len(species_list) == 0:
                print(f"✗ Error: No matching species!")
                return None

            print(f"\nProcessing {len(species_list)} species:")
            for sp in species_list:
                print(f"  - {sp}")
        else:
            species_list = []
            for key in openwq_results.keys():
                if chemSpec in key:
                    species_list.append(key)
                    break

            if len(species_list) == 0:
                print(f"✗ Error: '{chemSpec}' not found!")
                return None

        # Load mizuRoute output for reference (if needed)
        if isinstance(hydromodel_info, str):
            mizuroute_out_file = netCDF4.Dataset(hydromodel_info, 'r')
        else:
            mizuroute_out_file = None

    # Process each species/variable
    for idx, species_key in enumerate(species_list):
        print("\n" + "=" * 70)
        print(f"PROCESSING {idx + 1}/{len(species_list)}: {species_key}")
        print("=" * 70)

        if what2map == 'hostmodel':
            # Hostmodel: already have ttdata
            compartment = hostmodel
            var2print_i = var2print
            xyz_coords = None

        else:  # openwq
            # Parse species key
            parts = species_key.split('@')
            compartment = parts[0]
            chem_unit = parts[1].split('#')
            var2print_i = chem_unit[0]
            units = chem_unit[1]

            # Get data from results dictionary
            print(f"\nExtracting data for {var2print_i}...")
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
                print(f"✗ No data for extension '{file_extension}', skipping...")
                continue

        # Create output directory
        var2print_output_dir = os.path.join(output_dir, f'maps_{var2print_i}')
        os.makedirs(var2print_output_dir, exist_ok=True)

        # Determine timeframes
        print("\nDetermining timeframes...")
        total_timeframes = len(ttdata)

        if timeframes is None:
            timeframes_to_plot = list(range(total_timeframes))
            print(f"  Using ALL {total_timeframes} timeframes")
        else:
            if timeframes >= total_timeframes:
                timeframes_to_plot = list(range(total_timeframes))
                print(f"  Using ALL {total_timeframes} timeframes")
            else:
                timeframes_to_plot = list(range(0, total_timeframes, max(1, total_timeframes // timeframes)))
                timeframes_to_plot = timeframes_to_plot[:timeframes]
                print(f"  Using {len(timeframes_to_plot)} equally spaced timeframes")
                print(f"  Covering: {ttdata.index[timeframes_to_plot[0]]} to {ttdata.index[timeframes_to_plot[-1]]}")

        # Calculate color scale
        print("\nCalculating color scale...")
        all_values = []
        for t_idx in timeframes_to_plot[::max(1, len(timeframes_to_plot) // 10)]:
            data_values = ttdata.iloc[t_idx, :].values
            all_values.extend([v for v in data_values if not np.isnan(v)])

        vmin = np.percentile(all_values, 1)
        vmax = np.percentile(all_values, 99)
        print(f"  Color range: [{vmin:.3e}, {vmax:.3e}]")

        cmap = plt.cm.RdYlGn_r
        norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

        # CREATE MAPS
        print(f"\nCreating {len(timeframes_to_plot)} static maps...")

        for i, t_idx in enumerate(timeframes_to_plot):
            # Get data
            data_values = ttdata.iloc[t_idx, :].values
            data_hostmdlKeys = list(ttdata.iloc[t_idx, :].index)

            # Map to geodataframe
            geodf_openwq_wq = _map_data_to_geodf(
                data_values, xyz_coords, geodf, shpf_mapKey, data_hostmdlKeys
            )

            # Create figure
            fig, ax = plt.subplots(figsize=(12, 10))

            # Plot features - HANDLE BOTH LINES AND POLYGONS
            for j, geom in enumerate(geodf.geometry):
                value = geodf_openwq_wq[j]

                if np.isnan(value):
                    color = 'gray'
                    alpha = 0.3
                else:
                    color = cmap(norm(value))
                    alpha = 0.8

                # LINES (mizuRoute)
                if geom.geom_type == 'LineString':
                    coords = list(geom.coords)
                    xs, ys = zip(*coords)
                    ax.plot(xs, ys, color=color, linewidth=2, alpha=alpha)

                elif geom.geom_type == 'MultiLineString':
                    for line in geom.geoms:
                        coords = list(line.coords)
                        xs, ys = zip(*coords)
                        ax.plot(xs, ys, color=color, linewidth=2, alpha=alpha)

                # POLYGONS (SUMMA)
                elif geom.geom_type == 'Polygon':
                    # Get exterior coordinates
                    xs, ys = geom.exterior.xy
                    polygon = MplPolygon(
                        list(zip(xs, ys)),
                        facecolor=color,
                        edgecolor='black',
                        linewidth=0.5,
                        alpha=alpha
                    )
                    ax.add_patch(polygon)

                elif geom.geom_type == 'MultiPolygon':
                    for poly in geom.geoms:
                        xs, ys = poly.exterior.xy
                        polygon = MplPolygon(
                            list(zip(xs, ys)),
                            facecolor=color,
                            edgecolor='black',
                            linewidth=0.5,
                            alpha=alpha
                        )
                        ax.add_patch(polygon)

            # Add colorbar
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            cbar = plt.colorbar(sm, ax=ax, fraction=0.046, pad=0.04)
            cbar.set_label(f'{var2print_i} ({units})', rotation=270, labelpad=20)

            # Set title
            timestamp = str(ttdata.index[t_idx])
            if what2map == 'hostmodel':
                ax.set_title(f'{hostmodel.upper()} - {var2print_i} - {timestamp}', fontsize=14, pad=20)
            else:
                ax.set_title(f'{var2print_i} Concentrations - {timestamp}', fontsize=14, pad=20)

            ax.set_xlabel('Longitude', fontsize=12)
            ax.set_ylabel('Latitude', fontsize=12)
            ax.set_aspect('equal')
            ax.grid(True, alpha=0.3)

            # Set axis limits based on geodataframe bounds
            bounds = geodf.total_bounds
            ax.set_xlim(bounds[0], bounds[2])
            ax.set_ylim(bounds[1], bounds[3])

            # Save
            output_file = os.path.join(var2print_output_dir, f'map_{i:04d}.png')
            plt.savefig(output_file, dpi=150, bbox_inches='tight')
            plt.close()

            if (i + 1) % 10 == 0 or (i + 1) == len(timeframes_to_plot):
                print(f"    Created {i + 1}/{len(timeframes_to_plot)} maps...")

        created_dirs.append(var2print_output_dir)

        print("=" * 70)
        print(f"✓ PNG maps complete for {var2print_i}!")
        print(f"  Directory: {var2print_output_dir}")
        print(f"  Maps: {len(timeframes_to_plot)} PNG files")
        print("=" * 70)

        # CREATE GIF
        if create_gif and len(timeframes_to_plot) > 1:
            print("\n" + "=" * 70)
            print(f"Creating animated GIF for {var2print_i}...")
            print("=" * 70)

            if what2map == 'hostmodel':
                model_print = hostmodel
            else:
                model_print = 'openwq'

            gif_path = os.path.join(output_dir, f'{model_print}_{var2print_i}_animation.gif')

            if _create_gif_from_pngs(var2print_output_dir, gif_path, gif_duration):
                created_gifs.append(gif_path)

    # Close files
    if what2map == 'hostmodel':
        nc_file.close()
    elif what2map == 'openwq' and 'mizuroute_out_file' in locals() and mizuroute_out_file is not None:
        mizuroute_out_file.close()

    # SUMMARY
    print("\n" + "=" * 70)
    print("✓ ALL COMPLETE!")
    print("=" * 70)
    print(f"\nPNG Maps:")
    for d in created_dirs:
        png_count = len([f for f in os.listdir(d) if f.endswith('.png')])
        print(f"    - {d} ({png_count} images)")

    if len(created_gifs) > 0:
        print(f"\nAnimated GIFs:")
        for g in created_gifs:
            size_mb = os.path.getsize(g) / (1024 * 1024)
            print(f"    - {g} ({size_mb:.1f} MB)")
    else:
        print(f"\nAnimated GIFs: None created")

    print("=" * 70)

    return {
        'image_dirs': created_dirs,
        'gifs': created_gifs
    }


if __name__ == '__main__':
    print("\nOpenWQ Static PNG Maps with Automatic GIF Creation")
    print("\nSupports:")
    print("  • OpenWQ water quality results (what2map='openwq')")
    print("  • mizuRoute outputs - polylines (hostmodel='mizuroute')")
    print("  • SUMMA outputs - polygons (hostmodel='summa')")