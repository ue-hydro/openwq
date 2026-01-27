"""
OpenWQ Simple Static Maps with Automatic GIF Creation
Creates individual static maps for each timestep AND automatically creates GIF
No special installs needed - tries multiple methods to create GIF
"""

import numpy as np
import pandas as pd
import geopandas as gpd
import netCDF4
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import os
import subprocess


def _map_data_to_geodf(data_values, xyz_coords, geodf, shpf_mapKey, data_hostmdlKeys):
    """Map OpenWQ data to geodataframe - SAME AS WORKING VERSION"""
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
    """
    Try multiple methods to create GIF from PNG files
    Returns True if successful, False otherwise
    """

    print(f"\n  Attempting to create GIF: {output_gif}")
    print(f"  Frame duration: {duration_ms}ms")

    # Get all PNG files
    png_files = sorted([f for f in os.listdir(image_dir) if f.endswith('.png')])

    if len(png_files) == 0:
        print("  ✗ No PNG files found")
        return False

    print(f"  Found {len(png_files)} PNG files")

    # Method 1: Try imageio (most reliable Python method)
    try:
        import imageio.v2 as imageio
        print("  → Trying imageio...")

        images = []
        for filename in png_files:
            filepath = os.path.join(image_dir, filename)
            images.append(imageio.imread(filepath))

        # Duration in seconds for imageio
        imageio.mimsave(output_gif, images, duration=duration_ms/1000.0, loop=0)

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

    # Method 3: Try ImageMagick (command line)
    try:
        print("  → Trying ImageMagick (convert)...")

        # Calculate delay (ImageMagick uses hundredths of a second)
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


def Map_h5_driver(shpfile_fullpath_mapKey=None,
                  openwq_results=None,
                  hydromodel_out_fullpath=None,
                  output_html_path=None,
                  chemical_species=None,
                  file_extension='main',
                  timestep='all',
                  hostmodel='mizuroute',
                  max_timesteps=20,
                  create_gif=True,
                  gif_duration=None):
    """
    Create static PNG maps AND animated GIF automatically

    Parameters:
    -----------
    shpfile_fullpath_mapKey : dict
        Dictionary with 'path_to_shp' and 'mapping_key' keys
    openwq_results : dict
        Results dictionary from Read_h5_driver
    hydromodel_out_fullpath : str
        Path to hydromodel NetCDF output file
    output_html_path : str
        Base path for output (will be used for folder name)
    chemical_species : list or None
        List of chemical species to map
    file_extension : str
        File extension identifier (default: 'main')
    timestep : str, int, or list
        Timestep(s) to plot (default: 'all')
    max_timesteps : int, optional
        Maximum number of timesteps (default: 20)
    create_gif : bool, optional
        Whether to create GIF automatically (default: True)
    gif_duration : int, optional
        Duration per frame in milliseconds (default: 500)

    Returns:
    --------
    dict : {'image_dirs': [...], 'gifs': [...]}
    """

    print("=" * 70)
    print("STATIC PNG MAPS + AUTOMATIC GIF CREATION")
    print("=" * 70)

    # Determine which species to process
    if chemical_species is None:
        species_list = list(openwq_results.keys())
        print(f"\nProcessing ALL {len(species_list)} species")
    elif isinstance(chemical_species, list):
        species_list = []
        for chem in chemical_species:
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
            if chemical_species in key:
                species_list.append(key)
                break

        if len(species_list) == 0:
            print(f"✗ Error: '{chemical_species}' not found!")
            return None

    # Load shapefile ONCE
    print("\n" + "=" * 70)
    print("STEP 1: Loading shapefile...")
    print("=" * 70)

    shpfile_fullpath = shpfile_fullpath_mapKey["path_to_shp"]
    geodf = gpd.read_file(shpfile_fullpath, engine='fiona')
    geodf = geodf.to_crs("EPSG:4326")

    shapefile_mappingKey = shpfile_fullpath_mapKey["mapping_key"]
    shpf_mapKey = np.array(geodf[shapefile_mappingKey])

    print(f"  Loaded {len(geodf)} features")

    # Load mizuRoute output ONCE
    print("\nSTEP 2: Loading mizuRoute output...")
    mizuroute_out_file = netCDF4.Dataset(hydromodel_out_fullpath, 'r')

    # Prepare output directory
    output_dir = os.path.dirname(output_html_path) if output_html_path else '.'

    created_dirs = []
    created_gifs = []

    # Process each species
    for idx, species_key in enumerate(species_list):
        print("\n" + "=" * 70)
        print(f"SPECIES {idx + 1}/{len(species_list)}: {species_key}")
        print("=" * 70)

        # Parse species key
        parts = species_key.split('@')
        compartment = parts[0]
        chem_unit = parts[1].split('#')
        chemical = chem_unit[0]
        units = chem_unit[1]

        # Create output directory for this species
        species_output_dir = os.path.join(output_dir, f'maps_{chemical}')
        os.makedirs(species_output_dir, exist_ok=True)

        # Get data from results dictionary
        print(f"\nExtracting data for {chemical}...")
        data_found = False
        for ext, data_list in openwq_results[species_key]:
            if ext == file_extension:
                if data_list and len(data_list) > 0:
                    filename, ttdata, xyz_coords = data_list[0]
                    data_found = True
                    print(f"  Data shape: {ttdata.shape}")
                    print(f"  Time: {ttdata.index[0]} to {ttdata.index[-1]}")
                    break

        if not data_found:
            print(f"✗ No data for extension '{file_extension}', skipping...")
            continue

        # Determine timesteps
        print("\nDetermining timesteps...")

        if timestep == 'all':
            total_timesteps = len(ttdata)

            if max_timesteps and max_timesteps < total_timesteps:
                timesteps_to_plot = list(range(0, total_timesteps, max(1, total_timesteps // max_timesteps)))
                timesteps_to_plot = timesteps_to_plot[:max_timesteps]
                print(f"  Using {len(timesteps_to_plot)}/{total_timesteps} timesteps")
            else:
                timesteps_to_plot = list(range(total_timesteps))
                print(f"  Using all {total_timesteps} timesteps")
        elif isinstance(timestep, list):
            min_val, max_val = timestep[0], timestep[1]

            if isinstance(min_val, int) and isinstance(max_val, int):
                timesteps_to_plot = list(range(min_val, max_val + 1))
            else:
                try:
                    min_date = pd.to_datetime(str(min_val))
                    max_date = pd.to_datetime(str(max_val))

                    timesteps_to_plot = []
                    for i, ts in enumerate(ttdata.index):
                        ts_datetime = pd.to_datetime(ts)
                        if min_date <= ts_datetime <= max_date:
                            timesteps_to_plot.append(i)

                    if len(timesteps_to_plot) == 0:
                        print(f"✗ No timesteps in range, skipping...")
                        continue
                except Exception as e:
                    print(f"✗ Error parsing dates: {e}, skipping...")
                    continue

            print(f"  Using {len(timesteps_to_plot)} timesteps")
        else:
            timestep_idx = timestep if isinstance(timestep, int) else ttdata.index.get_loc(timestep)
            timesteps_to_plot = [timestep_idx]
            print(f"  Using timestep {timestep_idx}")

        # Calculate global color scale
        print("\nCalculating color scale...")
        all_values = []
        for t_idx in timesteps_to_plot[::max(1, len(timesteps_to_plot)//10)]:
            data_values = ttdata.iloc[t_idx, :].values
            all_values.extend([v for v in data_values if not np.isnan(v)])

        vmin = np.percentile(all_values, 1)
        vmax = np.percentile(all_values, 99)
        print(f"  Color range: [{vmin:.3e}, {vmax:.3e}]")

        # Create colormap
        cmap = plt.cm.RdYlGn_r
        norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

        # CREATE INDIVIDUAL MAPS
        print(f"\nCreating {len(timesteps_to_plot)} static maps...")

        for i, t_idx in enumerate(timesteps_to_plot):
            # Get data for this timestep
            data_values = ttdata.iloc[t_idx, :].values
            data_hostmdlKeys = list(ttdata.iloc[t_idx, :].index)

            # Map to geodataframe
            geodf_openwq_wq = _map_data_to_geodf(
                data_values, xyz_coords, geodf, shpf_mapKey, data_hostmdlKeys
            )

            # Create figure
            fig, ax = plt.subplots(figsize=(12, 10))

            # Plot each feature with color based on concentration
            for j, geom in enumerate(geodf.geometry):
                value = geodf_openwq_wq[j]

                if np.isnan(value):
                    color = 'gray'
                    alpha = 0.3
                else:
                    color = cmap(norm(value))
                    alpha = 0.8

                if geom.geom_type == 'LineString':
                    coords = list(geom.coords)
                    xs, ys = zip(*coords)
                    ax.plot(xs, ys, color=color, linewidth=2, alpha=alpha)
                elif geom.geom_type == 'MultiLineString':
                    for line in geom.geoms:
                        coords = list(line.coords)
                        xs, ys = zip(*coords)
                        ax.plot(xs, ys, color=color, linewidth=2, alpha=alpha)

            # Add colorbar
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            cbar = plt.colorbar(sm, ax=ax, fraction=0.046, pad=0.04)
            cbar.set_label(f'{chemical} ({units})', rotation=270, labelpad=20)

            # Set title and labels
            timestamp = str(ttdata.index[t_idx])
            ax.set_title(f'{chemical} Concentrations - {timestamp}', fontsize=14, pad=20)
            ax.set_xlabel('Longitude', fontsize=12)
            ax.set_ylabel('Latitude', fontsize=12)
            ax.set_aspect('equal')
            ax.grid(True, alpha=0.3)

            # Save figure
            output_file = os.path.join(species_output_dir, f'map_{i:04d}.png')
            plt.savefig(output_file, dpi=150, bbox_inches='tight')
            plt.close()

            if (i + 1) % 10 == 0 or (i + 1) == len(timesteps_to_plot):
                print(f"    Created {i + 1}/{len(timesteps_to_plot)} maps...")

        created_dirs.append(species_output_dir)

        print("=" * 70)
        print(f"✓ PNG maps complete for {chemical}!")
        print(f"  Directory: {species_output_dir}")
        print(f"  Maps: {len(timesteps_to_plot)} PNG files")
        print("=" * 70)

        # CREATE GIF
        if create_gif and len(timesteps_to_plot) > 1:
            print("\n" + "=" * 70)
            print(f"Creating animated GIF for {chemical}...")
            print("=" * 70)

            gif_path = os.path.join(output_dir, f'{chemical}_animation.gif')

            if _create_gif_from_pngs(species_output_dir, gif_path, gif_duration):
                created_gifs.append(gif_path)
            else:
                print("\n  Manual GIF creation commands:")
                print(f"  imageio:      python create_gif.py {species_output_dir}")
                print(f"  ImageMagick:  cd {species_output_dir} && convert -delay 50 -loop 0 *.png animation.gif")

    # Close mizuRoute file
    mizuroute_out_file.close()

    # SUMMARY
    print("\n" + "=" * 70)
    print("✓ ALL COMPLETE!")
    print("=" * 70)
    print(f"\nPNG Maps:")
    print(f"  Species: {len(created_dirs)}")
    for d in created_dirs:
        png_count = len([f for f in os.listdir(d) if f.endswith('.png')])
        print(f"    - {d} ({png_count} images)")

    if len(created_gifs) > 0:
        print(f"\nAnimated GIFs:")
        print(f"  Created: {len(created_gifs)}")
        for g in created_gifs:
            size_mb = os.path.getsize(g) / (1024 * 1024)
            print(f"    - {g} ({size_mb:.1f} MB)")
    else:
        print(f"\nAnimated GIFs: None created")
        print("  Install imageio, Pillow, or ImageMagick to enable")

    print("=" * 70)

    return {
        'image_dirs': created_dirs,
        'gifs': created_gifs
    }


if __name__ == '__main__':
    print("\nOpenWQ Static PNG Maps with Automatic GIF Creation")
    print("\nAdvantages:")
    print("  ✓ Always creates PNG maps")
    print("  ✓ Automatically creates GIF if possible")
    print("  ✓ Tries multiple GIF creation methods")
    print("  ✓ Works even if GIF creation fails")