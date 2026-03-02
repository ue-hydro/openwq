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
OpenWQ WebGL 3D Interactive River Viewer Export
================================================

Exports OpenWQ water quality and mizuRoute flow results as an interactive
WebGL 3D viewer. River network segments (polylines) are rasterized onto a
regular grid, with velocity vectors computed from polyline tangent direction
scaled by flow magnitude. The resulting PNG textures are used by a GPU-based
particle advection system (adapted from the FLUXOS viewer) to visualize
flow and water quality concentrations as animated particles.

Output:
    A self-contained directory with index.html + data files that can be served
    via any HTTP server for interactive visualization.

Usage:
    import WebGL_h5_driver as h5_wlib
    h5_wlib.WebGL_h5_driver(
        what2map='openwq', hostmodel='mizuroute',
        shpfile_info=shpfile_info,
        openwq_results=openwq_results,
        chemSpec=["NO3-N"],
        hydromodel_info=hydromodel_info,
        hydromodel_var2print='DWroutedRunoff',
        output_dir='openwq_webgl_viewer',
    )

Requirements:
    numpy, pandas, geopandas, netCDF4, shapely (all already used by Map_h5_driver)
"""

import numpy as np
import pandas as pd
import geopandas as gpd
import netCDF4
import os
import struct
import zlib
import json as _json
import pathlib
import math
import urllib.request
import io
import fiona
from shapely.geometry import shape, LineString, MultiLineString, Point, MultiPoint, Polygon, MultiPolygon


# ═══════════════════════════════════════════════════════════════════════════════
#  MINIMAL PNG WRITER  (no Pillow / PIL dependency)
#  Copied from FLUXOS fluxos_viewer.py
# ═══════════════════════════════════════════════════════════════════════════════

def _write_png(rgba_array):
    """
    Encode an (H, W, 4) uint8 RGBA array as a PNG file in memory.

    Returns bytes of the PNG file.
    """
    h, w = rgba_array.shape[:2]

    def _chunk(chunk_type, data):
        c = chunk_type + data
        crc = struct.pack(">I", zlib.crc32(c) & 0xFFFFFFFF)
        return struct.pack(">I", len(data)) + c + crc

    # Signature
    sig = b"\x89PNG\r\n\x1a\n"

    # IHDR
    ihdr_data = struct.pack(">IIBBBBB", w, h, 8, 6, 0, 0, 0)  # 8-bit RGBA
    ihdr = _chunk(b"IHDR", ihdr_data)

    # IDAT — raw image data with filter byte 0 (None) per row
    raw_rows = []
    for y in range(h):
        raw_rows.append(b"\x00")  # filter: None
        raw_rows.append(rgba_array[y].tobytes())
    raw = b"".join(raw_rows)
    compressed = zlib.compress(raw, 9)
    idat = _chunk(b"IDAT", compressed)

    # IEND
    iend = _chunk(b"IEND", b"")

    return sig + ihdr + idat + iend


# ═══════════════════════════════════════════════════════════════════════════════
#  HILLSHADE COMPUTATION
#  Adapted from FLUXOS fluxos_viewer.py
# ═══════════════════════════════════════════════════════════════════════════════

def _compute_hillshade(dem, cellsize, azimuth=315, altitude=45):
    """Compute a hillshade from a DEM.

    Returns an (H, W) uint8 grayscale image (0-255).
    NaN cells in the DEM are rendered as neutral gray (128).
    """
    filled = dem.copy()
    nan_mask = np.isnan(filled)
    if np.any(nan_mask):
        try:
            from scipy.ndimage import distance_transform_edt
            idx = distance_transform_edt(nan_mask, return_distances=False,
                                          return_indices=True)
            filled = filled[tuple(idx)]
        except Exception:
            filled[nan_mask] = np.nanmean(dem)

    dy, dx = np.gradient(filled, cellsize)
    slope = np.arctan(np.sqrt(dx**2 + dy**2))
    aspect = np.arctan2(-dy, dx)

    az_rad = np.radians(azimuth)
    alt_rad = np.radians(altitude)

    hs = (np.sin(alt_rad) * np.cos(slope) +
          np.cos(alt_rad) * np.sin(slope) * np.cos(az_rad - aspect))
    hs = np.clip(hs, 0, 1)
    hs_uint8 = (hs * 255).astype(np.uint8)
    hs_uint8[nan_mask] = 128
    return hs_uint8


# ═══════════════════════════════════════════════════════════════════════════════
#  ROBUST SHAPEFILE LOADER  (handles broken pyproj / CRS issues)
# ═══════════════════════════════════════════════════════════════════════════════

def _transform_geometry_to_wgs84(geom, src_proj):
    """Transform a shapely geometry from projected CRS to WGS84 lon/lat.

    Uses pyproj.Proj inverse projection, which avoids the broken SQLite
    database in pyproj 2.6.x that prevents EPSG lookups.

    Parameters
    ----------
    geom : shapely geometry
        Geometry in the source projected CRS.
    src_proj : pyproj.Proj
        Source projection object.

    Returns
    -------
    shapely geometry in WGS84 (lon/lat).
    """
    def _transform_coords(coords, proj):
        return [proj(x, y, inverse=True) for x, y in coords]

    if geom.geom_type == 'LineString':
        return LineString(_transform_coords(geom.coords, src_proj))
    elif geom.geom_type == 'MultiLineString':
        return MultiLineString([
            LineString(_transform_coords(line.coords, src_proj))
            for line in geom.geoms
        ])
    elif geom.geom_type == 'Point':
        lon, lat = src_proj(geom.x, geom.y, inverse=True)
        return Point(lon, lat)
    elif geom.geom_type == 'MultiPoint':
        return MultiPoint([
            Point(*src_proj(p.x, p.y, inverse=True))
            for p in geom.geoms
        ])
    elif geom.geom_type == 'Polygon':
        exterior = _transform_coords(geom.exterior.coords, src_proj)
        interiors = [_transform_coords(ring.coords, src_proj) for ring in geom.interiors]
        return Polygon(exterior, interiors)
    elif geom.geom_type == 'MultiPolygon':
        polys = []
        for poly in geom.geoms:
            exterior = _transform_coords(poly.exterior.coords, src_proj)
            interiors = [_transform_coords(ring.coords, src_proj) for ring in poly.interiors]
            polys.append(Polygon(exterior, interiors))
        return MultiPolygon(polys)
    else:
        return geom


def _load_shapefile_robust(shpfile_fullpath):
    """Load a shapefile and reproject to WGS84, with robust CRS handling.

    Handles broken pyproj installations (e.g., pyproj 2.6.x with outdated
    SQLite database) by falling back to fiona + manual inverse projection.

    Parameters
    ----------
    shpfile_fullpath : str
        Path to shapefile (.shp).

    Returns
    -------
    GeoDataFrame with geometries in WGS84 (lon/lat).
    """
    # ── Attempt 1: Standard geopandas loading ──
    try:
        geodf = gpd.read_file(shpfile_fullpath)
        if geodf.crs is not None and not geodf.crs.equals("EPSG:4326"):
            geodf = geodf.to_crs("EPSG:4326")
        print("  Loaded with geopandas (standard path)")
        return geodf
    except Exception as e:
        print(f"  ⚠ Standard loading failed: {e}")

    # ── Attempt 2: Fiona + manual reprojection ──
    print("  Falling back to fiona + manual reprojection...")
    try:
        from pyproj import Proj

        with fiona.open(shpfile_fullpath, 'r') as src:
            crs_dict = dict(src.crs)
            features = list(src)

        print(f"  Fiona read {len(features)} features, CRS: {crs_dict}")

        # Check if coordinates are already in lon/lat
        is_longlat = crs_dict.get('proj') in ('longlat', 'latlong', 'latlon')

        # Build geometries
        geoms_src = [shape(f['geometry']) for f in features]
        props = [dict(f['properties']) for f in features]

        if is_longlat:
            geoms_wgs84 = geoms_src
            print("  CRS is already lon/lat, no reprojection needed")
        else:
            # Reproject using inverse Proj (avoids broken SQLite lookups)
            src_proj = Proj(crs_dict)
            geoms_wgs84 = [_transform_geometry_to_wgs84(g, src_proj) for g in geoms_src]
            print("  Reprojected to WGS84 using pyproj inverse projection")

        # Build GeoDataFrame (without setting CRS to avoid pyproj issues)
        df = pd.DataFrame(props)
        geodf = gpd.GeoDataFrame(df, geometry=geoms_wgs84)

        # Verify coordinates are in valid WGS84 range
        bounds = geodf.total_bounds
        if (-200 < bounds[0] < 200 and -100 < bounds[1] < 100 and
                -200 < bounds[2] < 200 and -100 < bounds[3] < 100):
            print(f"  Coordinates look valid (bounds: lon [{bounds[0]:.4f}, {bounds[2]:.4f}], "
                  f"lat [{bounds[1]:.4f}, {bounds[3]:.4f}])")
        else:
            print(f"  ⚠ Warning: coordinates may not be in WGS84! Bounds: {bounds}")

        return geodf

    except Exception as e2:
        raise RuntimeError(
            f"Failed to load shapefile with both geopandas and fiona fallback.\n"
            f"  geopandas error: see above\n"
            f"  fiona error: {e2}\n"
            f"  File: {shpfile_fullpath}"
        )


# ═══════════════════════════════════════════════════════════════════════════════
#  GRID AND RASTERIZATION UTILITIES
# ═══════════════════════════════════════════════════════════════════════════════

def _compute_grid_params(geodf, grid_resolution=None, target_cells=600):
    """Compute raster grid parameters from shapefile bounds.

    Parameters
    ----------
    geodf : GeoDataFrame
        Loaded shapefile with geometry in EPSG:4326
    grid_resolution : float or None
        Cell size in degrees. If None, auto-compute targeting ~target_cells
        on the longest axis.
    target_cells : int
        Target number of cells on the longest axis (default: 600)

    Returns
    -------
    dict : {'xmin', 'ymin', 'xmax', 'ymax', 'cellsize', 'ncols', 'nrows'}
    """
    bounds = geodf.total_bounds  # (minx, miny, maxx, maxy)
    xmin, ymin, xmax, ymax = bounds

    # Add 5% padding
    dx = xmax - xmin
    dy = ymax - ymin
    pad = max(dx, dy) * 0.05
    xmin -= pad
    ymin -= pad
    xmax += pad
    ymax += pad
    dx = xmax - xmin
    dy = ymax - ymin

    if grid_resolution is not None:
        cellsize = grid_resolution
    else:
        cellsize = max(dx, dy) / target_cells

    ncols = max(1, int(np.ceil(dx / cellsize)))
    nrows = max(1, int(np.ceil(dy / cellsize)))

    # Cap at reasonable maximum
    if ncols > 1200 or nrows > 1200:
        print(f"  ⚠ Grid too large ({ncols}x{nrows}), capping at 1200 cells")
        cellsize = max(dx, dy) / 1200
        ncols = max(1, int(np.ceil(dx / cellsize)))
        nrows = max(1, int(np.ceil(dy / cellsize)))

    return {
        'xmin': xmin, 'ymin': ymin,
        'xmax': xmax, 'ymax': ymax,
        'cellsize': cellsize,
        'ncols': ncols, 'nrows': nrows
    }


def _rasterize_segments(geodf, shpf_mapKey, grid_params, river_width_cells=3):
    """Rasterize shapefile polylines onto a regular grid.

    For each polyline segment, walk along the geometry at sub-cell intervals
    and mark grid cells as 'river'. Also compute per-cell tangent vectors
    (normalized flow direction from the polyline geometry).

    Parameters
    ----------
    geodf : GeoDataFrame
        Shapefile with LineString/MultiLineString or Polygon/MultiPolygon geometries
    shpf_mapKey : ndarray
        Segment ID for each feature
    grid_params : dict
        From _compute_grid_params
    river_width_cells : int
        Width of rasterized river in grid cells (default 3).

    Returns
    -------
    tuple : (river_mask, tangent_x, tangent_y, segment_id_grid)
        river_mask : ndarray (nrows, ncols) bool
        tangent_x : ndarray (nrows, ncols) float — easting tangent component
        tangent_y : ndarray (nrows, ncols) float — northing tangent component
        segment_id_grid : ndarray (nrows, ncols) — segment ID at each river cell
    """
    nrows = grid_params['nrows']
    ncols = grid_params['ncols']
    cs = grid_params['cellsize']
    xmin = grid_params['xmin']
    ymax = grid_params['ymax']  # top edge (row 0 = north)

    river_mask = np.zeros((nrows, ncols), dtype=bool)
    tangent_x = np.zeros((nrows, ncols), dtype=np.float64)
    tangent_y = np.zeros((nrows, ncols), dtype=np.float64)
    segment_id_grid = np.full((nrows, ncols), -1, dtype=np.int64)

    # Dilation offsets for river width
    half_w = river_width_cells // 2
    offsets = []
    for dr in range(-half_w, half_w + 1):
        for dc in range(-half_w, half_w + 1):
            if abs(dr) + abs(dc) <= river_width_cells:  # diamond shape
                offsets.append((dr, dc))

    step_dist = cs * 0.4  # sub-cell step distance

    for feat_idx, geom in enumerate(geodf.geometry):
        seg_id = shpf_mapKey[feat_idx]

        # Collect all linestrings from the geometry
        lines = []
        if geom is None:
            continue
        if geom.geom_type == 'LineString':
            lines = [geom]
        elif geom.geom_type == 'MultiLineString':
            lines = list(geom.geoms)
        elif geom.geom_type == 'Polygon':
            lines = [geom.exterior]
        elif geom.geom_type == 'MultiPolygon':
            lines = [p.exterior for p in geom.geoms]
        else:
            continue

        for line in lines:
            coords = list(line.coords)
            if len(coords) < 2:
                continue

            # Walk along polyline at sub-cell intervals
            total_length = line.length
            if total_length < 1e-10:
                continue

            n_steps = max(2, int(np.ceil(total_length / step_dist)))

            for si in range(n_steps):
                frac = si / max(1, n_steps - 1)
                pt = line.interpolate(frac, normalized=True)
                x, y = pt.x, pt.y

                # Compute tangent from neighboring interpolation points
                frac_prev = max(0, frac - 0.01)
                frac_next = min(1, frac + 0.01)
                pt_prev = line.interpolate(frac_prev, normalized=True)
                pt_next = line.interpolate(frac_next, normalized=True)

                tx = pt_next.x - pt_prev.x
                ty = pt_next.y - pt_prev.y
                tmag = np.sqrt(tx * tx + ty * ty)
                if tmag > 1e-12:
                    tx /= tmag
                    ty /= tmag

                # Map to grid cell
                col = int((x - xmin) / cs)
                row = int((ymax - y) / cs)  # row 0 = north

                # Paint cell and neighbors
                for dr, dc in offsets:
                    r = row + dr
                    c = col + dc
                    if 0 <= r < nrows and 0 <= c < ncols:
                        river_mask[r, c] = True
                        tangent_x[r, c] = tx
                        tangent_y[r, c] = ty
                        segment_id_grid[r, c] = seg_id

    n_river_cells = np.sum(river_mask)
    print(f"  Rasterized {len(geodf)} segments → {n_river_cells} river cells "
          f"({100 * n_river_cells / (nrows * ncols):.1f}% coverage)")

    return river_mask, tangent_x, tangent_y, segment_id_grid


def _map_data_to_grid(data_values, data_hostmdlKeys, river_mask, segment_id_grid, shpf_mapKey):
    """Map segment-level data values to the rasterized grid.

    Parameters
    ----------
    data_values : ndarray
        Values per segment for one timestep
    data_hostmdlKeys : list
        Column names like 'reachID_123'
    river_mask : ndarray (nrows, ncols) bool
    segment_id_grid : ndarray (nrows, ncols)
    shpf_mapKey : ndarray
        Segment IDs from shapefile

    Returns
    -------
    ndarray (nrows, ncols) : mapped values (0 where no data)
    """
    nrows, ncols = river_mask.shape
    grid_values = np.zeros((nrows, ncols), dtype=np.float64)

    # Build lookup: segment_id → value
    data_hostmdlKeys_mainInfo = [item.split('_')[-1] for item in data_hostmdlKeys]
    seg_to_value = {}
    for j, key_str in enumerate(data_hostmdlKeys_mainInfo):
        try:
            seg_to_value[int(float(key_str))] = data_values[j]
        except (ValueError, IndexError):
            seg_to_value[key_str] = data_values[j]

    # Map to grid using segment_id_grid
    river_rows, river_cols = np.where(river_mask)
    for i in range(len(river_rows)):
        r, c = river_rows[i], river_cols[i]
        seg_id = segment_id_grid[r, c]
        if seg_id in seg_to_value:
            val = seg_to_value[seg_id]
            if not np.isnan(val):
                grid_values[r, c] = val

    return grid_values


def _build_velocity_field(flow_grid, river_mask, tangent_x, tangent_y):
    """Build velocity grids from flow magnitude and tangent direction.

    Parameters
    ----------
    flow_grid : ndarray (nrows, ncols)
        Flow magnitude mapped to grid
    river_mask : ndarray (nrows, ncols) bool
    tangent_x, tangent_y : ndarray (nrows, ncols) float

    Returns
    -------
    tuple : (grid_vx, grid_vy, wet_mask, grid_h)
    """
    grid_vx = tangent_x * flow_grid
    grid_vy = tangent_y * flow_grid
    wet_mask = river_mask & (flow_grid > 0)
    grid_h = flow_grid.copy()

    return grid_vx, grid_vy, wet_mask, grid_h


def _encode_velocity_png(grid_vx, grid_vy, wet_mask, grid_h,
                          vx_min, vx_max, vy_min, vy_max, h_max):
    """Encode velocity/depth/wet into RGBA PNG bytes.

    The WebGL shader (from FLUXOS) cross-reads the R/G channels:
      - R channel is decoded with u_vel_range_x → drives vertical (northing) motion
      - G channel is decoded with u_vel_range_y → drives horizontal (easting) motion

    Therefore we encode:
        R = vy (northing velocity), normalized with vy_min/vy_max
        G = vx (easting velocity), normalized with vx_min/vx_max
        B = h / h_max * 255
        A = 255 if wet, 0 if dry

    And the metadata must store vx_min/vx_max = vy range, vy_min/vx_max = vx range
    (done in the caller).
    """
    # R channel: northing velocity (vy) — drives vertical particle motion
    r_ch = ((grid_vy - vy_min) / max(vy_max - vy_min, 1e-12) * 255.0
            ).clip(0, 255).astype(np.uint8)
    # G channel: easting velocity (vx) — drives horizontal particle motion
    g_ch = ((grid_vx - vx_min) / max(vx_max - vx_min, 1e-12) * 255.0
            ).clip(0, 255).astype(np.uint8)
    b_ch = (grid_h / max(h_max, 1e-12) * 255.0).clip(0, 255).astype(np.uint8)
    a_ch = np.where(wet_mask, np.uint8(255), np.uint8(0))

    rgba = np.stack([r_ch, g_ch, b_ch, a_ch], axis=-1)
    return _write_png(rgba)


def _encode_concentration_png(grid_conc, conc_max, wet_mask):
    """Encode concentration into RGBA PNG bytes.

    Encoding matches FLUXOS format:
        R = conc / conc_max * 255
        G = 0
        B = 0
        A = 255 if wet and conc > 0, else 0
    """
    c_r = (grid_conc / max(conc_max, 1e-12) * 255.0).clip(0, 255).astype(np.uint8)
    c_g = np.zeros_like(c_r)
    c_b = np.zeros_like(c_r)
    c_a = np.where(wet_mask & (grid_conc > 0), np.uint8(255), np.uint8(0))

    c_rgba = np.stack([c_r, c_g, c_b, c_a], axis=-1)
    return _write_png(c_rgba)


def _generate_heightmap(grid_params, dem_path=None, river_mask=None):
    """Generate heightmap PNG and hillshade PNG.

    Parameters
    ----------
    grid_params : dict
        Grid parameters
    dem_path : str or None
        Path to DEM file (.asc or .tif). If None, generates pseudo-terrain
        with river channels carved into flat terrain.
    river_mask : ndarray or None
        Boolean mask of river cells (used for pseudo-terrain when no DEM).

    Returns
    -------
    tuple : (heightmap_bytes, hillshade_bytes, z_min, z_max)
    """
    nrows = grid_params['nrows']
    ncols = grid_params['ncols']
    cs = grid_params['cellsize']

    if dem_path is not None and os.path.isfile(dem_path):
        print(f"  Loading DEM: {dem_path}")
        try:
            import rasterio
            from rasterio.transform import from_bounds
            from rasterio.warp import reproject, Resampling

            with rasterio.open(dem_path) as src:
                dst_transform = from_bounds(
                    grid_params['xmin'], grid_params['ymin'],
                    grid_params['xmax'], grid_params['ymax'],
                    ncols, nrows
                )
                dem = np.empty((nrows, ncols), dtype=np.float64)
                reproject(
                    source=rasterio.band(src, 1),
                    destination=dem,
                    dst_transform=dst_transform,
                    dst_crs='EPSG:4326',
                    resampling=Resampling.bilinear
                )
                if src.nodata is not None:
                    dem[dem == src.nodata] = np.nan

            print(f"  DEM loaded: {nrows}x{ncols}, "
                  f"range=[{np.nanmin(dem):.1f}, {np.nanmax(dem):.1f}] m")

        except ImportError:
            print("  ⚠ rasterio not available, using pseudo-terrain")
            dem = None
        except Exception as e:
            print(f"  ⚠ DEM loading failed ({e}), using pseudo-terrain")
            dem = None
    else:
        if dem_path is not None:
            print(f"  ⚠ DEM file not found: {dem_path}, using pseudo-terrain")
        dem = None

    # Generate pseudo-terrain from river mask when no real DEM available
    if dem is None:
        dem = np.ones((nrows, ncols), dtype=np.float64) * 10.0  # base elevation
        if river_mask is not None and np.any(river_mask):
            # Carve river channels: river cells at elevation 0
            dem[river_mask] = 0.0
            # Smooth the terrain to create gentle slopes toward rivers
            try:
                from scipy.ndimage import gaussian_filter
                # Create distance-based terrain: smooth transition from 10 to 0
                dem = gaussian_filter(dem, sigma=5.0)
            except ImportError:
                # Manual simple smoothing: just set river and immediate neighbors
                pass
            print(f"  Generated pseudo-terrain with river channels")

    # Compute z range
    nan_mask = np.isnan(dem)
    valid_dem = dem.copy()
    valid_dem[nan_mask] = 0.0

    z_min = float(np.nanmin(dem)) if not np.all(nan_mask) else 0.0
    z_max = float(np.nanmax(dem)) if not np.all(nan_mask) else 1.0
    if z_max - z_min < 0.01:
        z_max = z_min + 1.0

    # Heightmap: 16-bit encoding in RGBA
    z_norm = np.clip((valid_dem - z_min) / (z_max - z_min), 0.0, 1.0)
    z_16bit = (z_norm * 65535.0).astype(np.uint16)
    r_ch = (z_16bit >> 8).astype(np.uint8)       # high byte
    g_ch = (z_16bit & 0xFF).astype(np.uint8)      # low byte
    b_ch = np.where(nan_mask, np.uint8(0), np.uint8(255))  # validity mask
    a_ch = np.full_like(r_ch, 255)

    hm_rgba = np.stack([r_ch, g_ch, b_ch, a_ch], axis=-1)
    heightmap_bytes = _write_png(hm_rgba)

    # Hillshade
    if np.all(dem == 0):
        # Flat terrain: uniform gray
        hs = np.full((nrows, ncols), 180, dtype=np.uint8)
    else:
        hs = _compute_hillshade(dem, cs)

    hs_rgba = np.stack([hs, hs, hs, np.full_like(hs, 255)], axis=-1)
    hillshade_bytes = _write_png(hs_rgba)

    return heightmap_bytes, hillshade_bytes, z_min, z_max


# ═══════════════════════════════════════════════════════════════════════════════
#  SATELLITE IMAGERY DOWNLOAD
#  Adapted from FLUXOS fluxos_viewer.py — downloads Esri World Imagery tiles
# ═══════════════════════════════════════════════════════════════════════════════

def _download_satellite(sw_lat, sw_lon, ne_lat, ne_lon,
                        img_w, img_h, output_path,
                        river_mask=None):
    """Download satellite imagery from Esri World Imagery tile service.

    Downloads map tiles, reprojects from Web Mercator to equirectangular,
    and saves as JPEG (or PNG fallback). Optionally darkens pixels outside
    the river mask.

    Parameters
    ----------
    sw_lat, sw_lon : float  — Southwest corner (lat/lon)
    ne_lat, ne_lon : float  — Northeast corner (lat/lon)
    img_w, img_h : int      — Output image dimensions
    output_path : str       — Output file path (.jpg)
    river_mask : ndarray or None — Boolean mask of river cells (True=river)
    """
    # Choose zoom level: aim for ~1 pixel per grid cell
    lat_mid = (sw_lat + ne_lat) / 2.0
    lon_span = ne_lon - sw_lon
    meters_per_pixel_z0 = 156543.03392 * math.cos(math.radians(lat_mid))
    target_mpp = (lon_span * 111320 * math.cos(math.radians(lat_mid))) / img_w
    zoom = max(1, min(18, int(round(math.log2(meters_per_pixel_z0 / target_mpp)))))
    print(f"  Satellite: zoom level {zoom}")

    def lat_lon_to_tile(lat, lon, z):
        n = 2 ** z
        x = int((lon + 180.0) / 360.0 * n)
        lat_r = math.radians(lat)
        y = int((1.0 - math.log(math.tan(lat_r) + 1.0 / math.cos(lat_r))
                 / math.pi) / 2.0 * n)
        return max(0, min(n - 1, x)), max(0, min(n - 1, y))

    def tile_bounds(tx, ty, z):
        n = 2 ** z
        west = tx / n * 360.0 - 180.0
        east = (tx + 1) / n * 360.0 - 180.0
        north = math.degrees(math.atan(math.sinh(
            math.pi * (1 - 2 * ty / n))))
        south = math.degrees(math.atan(math.sinh(
            math.pi * (1 - 2 * (ty + 1) / n))))
        return west, south, east, north

    tx_min, ty_min = lat_lon_to_tile(ne_lat, sw_lon, zoom)
    tx_max, ty_max = lat_lon_to_tile(sw_lat, ne_lon, zoom)

    n_tiles_x = tx_max - tx_min + 1
    n_tiles_y = ty_max - ty_min + 1
    print(f"  Satellite: {n_tiles_x}x{n_tiles_y} = "
          f"{n_tiles_x * n_tiles_y} tiles")

    # Download tiles into mosaic
    tile_px = 256
    mosaic_w = n_tiles_x * tile_px
    mosaic_h = n_tiles_y * tile_px
    mosaic = np.full((mosaic_h, mosaic_w, 3), 128, dtype=np.uint8)

    url_template = ("https://server.arcgisonline.com/ArcGIS/rest/services/"
                    "World_Imagery/MapServer/tile/{z}/{y}/{x}")

    for ty in range(ty_min, ty_max + 1):
        for tx in range(tx_min, tx_max + 1):
            url = url_template.format(z=zoom, y=ty, x=tx)
            try:
                req = urllib.request.Request(url, headers={
                    "User-Agent": "OpenWQ-viewer/1.0"})
                with urllib.request.urlopen(req, timeout=15) as resp:
                    tile_data = resp.read()
                try:
                    from PIL import Image as PILImage
                    tile_img = PILImage.open(io.BytesIO(tile_data))
                    tile_arr = np.array(tile_img.convert("RGB"))
                except ImportError:
                    continue
                oy = (ty - ty_min) * tile_px
                ox = (tx - tx_min) * tile_px
                th, tw = tile_arr.shape[:2]
                mosaic[oy:oy+th, ox:ox+tw] = tile_arr[:, :, :3]
            except Exception as e:
                print(f"    Tile {tx},{ty}: {e}")

    # Reproject: Mercator mosaic → equirectangular grid
    w0, _, _, n0 = tile_bounds(tx_min, ty_min, zoom)
    _, s1, e1, _ = tile_bounds(tx_max, ty_max, zoom)
    mosaic_west, mosaic_east = w0, e1
    mosaic_north, mosaic_south = n0, s1

    def lat_to_merc_y(lat):
        return math.log(math.tan(math.pi / 4 + math.radians(lat) / 2))

    merc_n = lat_to_merc_y(mosaic_north)
    merc_s = lat_to_merc_y(mosaic_south)
    merc_range = merc_n - merc_s
    lon_range = mosaic_east - mosaic_west

    out = np.full((img_h, img_w, 3), 40, dtype=np.uint8)

    # Pre-compute mosaic column for each output column
    col_lons = sw_lon + (np.arange(img_w) + 0.5) / img_w * (ne_lon - sw_lon)
    mosaic_cols = ((col_lons - mosaic_west) / lon_range * mosaic_w).astype(np.float64)

    for r in range(img_h):
        frac = (r + 0.5) / img_h
        lat = ne_lat - frac * (ne_lat - sw_lat)
        merc_y = lat_to_merc_y(lat)
        mosaic_row_f = (merc_n - merc_y) / merc_range * mosaic_h
        ry0 = max(0, min(mosaic_h - 2, int(mosaic_row_f)))
        ry1 = ry0 + 1
        fy = mosaic_row_f - ry0

        for c in range(img_w):
            mx = mosaic_cols[c]
            cx0 = max(0, min(mosaic_w - 2, int(mx)))
            cx1 = cx0 + 1
            fx = mx - cx0
            p00 = mosaic[ry0, cx0].astype(np.float64)
            p10 = mosaic[ry0, cx1].astype(np.float64)
            p01 = mosaic[ry1, cx0].astype(np.float64)
            p11 = mosaic[ry1, cx1].astype(np.float64)
            val = (p00 * (1 - fx) * (1 - fy) + p10 * fx * (1 - fy) +
                   p01 * (1 - fx) * fy + p11 * fx * fy)
            out[r, c] = np.clip(val, 0, 255).astype(np.uint8)

    # Brightness boost: auto-stretch histogram so the image is visible
    # (satellite imagery of boreal/arctic areas can be very dark)
    out_f = out.astype(np.float32)
    p2, p98 = np.percentile(out_f, 2), np.percentile(out_f, 98)
    if p98 > p2 + 5:
        out_f = (out_f - p2) / (p98 - p2) * 235 + 20  # map [p2,p98] → [20,255]
        out = np.clip(out_f, 0, 255).astype(np.uint8)

    # Gently darken pixels outside river area (keep them visible)
    if river_mask is not None:
        outside = ~river_mask
        out[outside] = (out[outside].astype(np.float32) * 0.55).astype(np.uint8)

    # Save as JPEG (PIL) or PNG (fallback)
    try:
        from PIL import Image as PILImage
        pil_img = PILImage.fromarray(out)
        pil_img.save(output_path, "JPEG", quality=85)
        fsize = os.path.getsize(output_path)
        print(f"  Satellite: {fsize / 1024:.0f} KB ({img_w}x{img_h} px)")
    except ImportError:
        sat_rgba = np.stack([out[:, :, 0], out[:, :, 1], out[:, :, 2],
                             np.full(out.shape[:2], 255, np.uint8)], axis=-1)
        sat_bytes = _write_png(sat_rgba)
        out_png = output_path.replace(".jpg", ".png")
        with open(out_png, "wb") as f:
            f.write(sat_bytes)
        print(f"  Satellite (PNG fallback): {len(sat_bytes) / 1024:.0f} KB")


def _generate_openwq_webgl_html(fluxos_viewer_path=None):
    """Generate the OpenWQ-adapted WebGL viewer HTML.

    Loads the pre-modified OpenWQ WebGL template (openwq_webgl_template.html)
    which includes all OpenWQ-specific features:
    1. Title: "OpenWQ Water Quality — Interactive River Viewer"
    2. Heading: "OpenWQ River Network"
    3. Time display: supports datetime strings
    4. Species dropdown for chemical species switching
    5. Colorbar labels adapt to variable type

    If the OpenWQ template is not found, falls back to the FLUXOS template
    and applies modifications on the fly.

    Parameters
    ----------
    fluxos_viewer_path : str or None
        Path to the FLUXOS fluxos_web/index.html template (fallback only).
        If None, searches common locations.

    Returns
    -------
    str : Complete HTML content
    """
    this_dir = pathlib.Path(__file__).parent

    # ── Primary: Load pre-modified OpenWQ template ──
    openwq_template = this_dir / 'openwq_webgl_template.html'
    if openwq_template.exists():
        html = openwq_template.read_text()
        print(f"  Loaded viewer template: {openwq_template}")
        return html

    # ── Fallback: Load FLUXOS template and modify on the fly ──
    search_paths = []
    if fluxos_viewer_path:
        search_paths.append(pathlib.Path(fluxos_viewer_path))
    search_paths.append(
        pathlib.Path.home() / 'Documents' / 'FLUXOS_cpp' / 'fluxos_web' / 'index.html'
    )

    html = None
    for path in search_paths:
        if path.exists():
            html = path.read_text()
            print(f"  Loaded FLUXOS viewer template: {path}")
            break

    if html is None:
        print("  ⚠ Could not find viewer template, using minimal fallback")
        return _generate_minimal_viewer_html()

    # Apply modifications
    # 1. Title
    html = html.replace(
        'FLUXOS Flood Simulation — 3D Flow Viewer',
        'OpenWQ Water Quality — Interactive River Viewer'
    )
    html = html.replace(
        'FLUXOS Flood Simulation',
        'OpenWQ River Network'
    )

    # 2. Add species selector after the concentration button div
    species_html = '''
    <div class="ctrl-row" id="species-row" style="display:none; margin-top:8px;">
        <label>Species</label>
        <select id="sel-species" class="btn" style="flex:1; margin:0 8px; padding:4px 8px; background:rgba(40,60,100,0.5); color:#b0c4de; border:1px solid rgba(80,120,180,0.4); border-radius:5px;">
        </select>
    </div>'''

    # Insert after the button row containing btn-conc
    insert_marker = '<!-- species-selector-placeholder -->'
    if insert_marker not in html:
        # Insert after the div containing btn-conc
        conc_btn_end = html.find('btn-conc')
        if conc_btn_end > 0:
            # Find the closing </div> of the button row
            next_div_close = html.find('</div>', conc_btn_end)
            if next_div_close > 0:
                insert_pos = next_div_close + len('</div>')
                html = html[:insert_pos] + species_html + html[insert_pos:]

    # 3. Add species switching JavaScript before the closing </script>
    species_js = '''

    // ── OpenWQ Species Switching ──
    initSpeciesSelector() {
        const m = this.meta;
        if (!m.chemical_species || m.chemical_species.length === 0) return;

        const sel = document.getElementById('sel-species');
        const row = document.getElementById('species-row');
        if (!sel || !row) return;

        // Populate dropdown
        m.chemical_species.forEach((sp, i) => {
            const opt = document.createElement('option');
            opt.value = sp;
            opt.textContent = sp;
            sel.appendChild(opt);
        });

        row.style.display = '';
        this.currentSpecies = m.chemical_species[0];

        sel.addEventListener('change', () => {
            this.currentSpecies = sel.value;
            // Clear concentration frame cache to force reload
            if (this.concFrameCache) this.concFrameCache.clear();
            console.log('Switched species to:', this.currentSpecies);
        });
    }
'''

    # 4. Modify time display to support datetime strings
    # Replace the formatTime function or add datetime support
    datetime_js = '''
    // Format time as datetime string if available
    formatTimeLabel(frameIdx) {
        const m = this.meta;
        if (m.time_format === 'datetime' && m.times && m.times[frameIdx]) {
            return m.times[frameIdx];
        }
        // Fallback: seconds
        const t = m.times ? m.times[frameIdx] : frameIdx;
        if (typeof t === 'number') {
            if (t < 60) return t.toFixed(0) + 's';
            if (t < 3600) return (t / 60).toFixed(1) + 'm';
            if (t < 86400) return (t / 3600).toFixed(1) + 'h';
            return (t / 86400).toFixed(1) + 'd';
        }
        return String(t);
    }
'''

    # 5. Modify concentration path to include species subdirectory
    # Replace 'data/conc/c_' with dynamic species path
    html = html.replace(
        "'data/conc/c_'",
        "'data/conc/' + (this.currentSpecies || '') + '/c_'"
    )
    # Also handle the fetch URL pattern
    html = html.replace(
        "`data/conc/c_${",
        "`data/conc/${this.currentSpecies || 'default'}/c_${"
    )

    # 6. Update colorbar labels
    html = html.replace(
        'Flow speed (m/s)',
        'Flow magnitude'
    )

    return html


def _generate_minimal_viewer_html():
    """Generate a minimal WebGL viewer HTML when FLUXOS template is not available.

    This is a simplified 2D particle viewer (no 3D terrain) that still uses
    the same PNG texture format.
    """
    return '''\
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>OpenWQ Water Quality — Interactive River Viewer</title>
<style>
* { margin: 0; padding: 0; box-sizing: border-box; }
html, body { width: 100%; height: 100%; overflow: hidden; background: #0a0a12; }
canvas#gl { width: 100%; height: 100%; display: block; }
#controls {
    position: fixed; bottom: 16px; left: 16px;
    background: rgba(5,8,20,0.82); color: #b0c4de;
    padding: 16px 20px; border-radius: 10px;
    font-family: 'SF Mono', 'Fira Code', 'Consolas', monospace;
    font-size: 12px; min-width: 300px;
    border: 1px solid rgba(60,90,140,0.3);
    backdrop-filter: blur(8px); user-select: none;
}
#controls h2 { font-size: 15px; color: #e0ecff; margin-bottom: 12px; font-weight: 600; }
.ctrl-row { margin-bottom: 8px; display: flex; align-items: center; justify-content: space-between; }
.ctrl-row label { flex: 0 0 100px; color: #8899bb; }
.ctrl-row input[type=range] { flex: 1; margin: 0 8px; }
.ctrl-row .val { flex: 0 0 80px; text-align: right; color: #ccddef; font-variant-numeric: tabular-nums; }
.btn { background: rgba(40,60,100,0.5); color: #b0c4de; border: 1px solid rgba(80,120,180,0.4);
    border-radius: 5px; padding: 6px 14px; cursor: pointer; font-family: inherit; font-size: 12px; }
.btn:hover { background: rgba(60,90,140,0.6); }
.btn.active { background: rgba(80,140,220,0.4); border-color: rgba(100,160,240,0.6); }
#info { position: fixed; top: 16px; right: 16px; background: rgba(5,8,20,0.7); color: #8899bb;
    padding: 10px 14px; border-radius: 8px; font-family: 'SF Mono', monospace; font-size: 11px;
    border: 1px solid rgba(60,90,140,0.2); pointer-events: none; }
#info .hl { color: #ccddef; }
#colorbar-wrap { position: fixed; bottom: 20px; right: 16px; background: rgba(5,8,20,0.7);
    padding: 8px 12px 6px; border-radius: 8px; border: 1px solid rgba(60,90,140,0.2); }
.cb-title { font-family: monospace; font-size: 11px; color: #8899bb; margin-bottom: 4px; }
.cb-labels { display: flex; justify-content: space-between; font-family: monospace; font-size: 10px; color: #8899bb; margin-top: 2px; }
canvas#colorbar { display: block; border-radius: 3px; }
#loading { position: fixed; top: 50%; left: 50%; transform: translate(-50%,-50%);
    color: #5577aa; font-family: monospace; font-size: 16px; }
input[type=range] { -webkit-appearance: none; height: 4px; background: rgba(60,90,140,0.4); border-radius: 2px; outline: none; }
input[type=range]::-webkit-slider-thumb { -webkit-appearance: none; width: 14px; height: 14px;
    border-radius: 50%; background: #5588cc; cursor: pointer; border: 2px solid #334466; }
</style>
</head>
<body>
<canvas id="gl"></canvas>
<div id="loading">Loading simulation data...</div>
<div id="controls" style="display:none;">
    <h2>OpenWQ River Network</h2>
    <div class="ctrl-row">
        <label>Time</label>
        <input type="range" id="sl-time" min="0" max="100" value="0">
        <span class="val" id="v-time">--</span>
    </div>
    <div class="ctrl-row" style="justify-content:center; gap:8px; margin:10px 0;">
        <button class="btn active" id="btn-play">&#9654; Play</button>
        <button class="btn" id="btn-reset">&#8634; Reset</button>
    </div>
    <div class="ctrl-row">
        <label>Speed</label>
        <input type="range" id="sl-speed" min="5" max="300" value="80">
        <span class="val" id="v-speed">0.8x</span>
    </div>
    <div class="ctrl-row">
        <label>Particles</label>
        <input type="range" id="sl-particles" min="10" max="18" value="14">
        <span class="val" id="v-particles">16K</span>
    </div>
    <div class="ctrl-row">
        <label>Trail length</label>
        <input type="range" id="sl-trail" min="880" max="998" value="970">
        <span class="val" id="v-trail">0.97</span>
    </div>
    <div class="ctrl-row" id="species-row" style="display:none; margin-top:8px;">
        <label>Species</label>
        <select id="sel-species" class="btn" style="flex:1; margin:0 8px;">
        </select>
    </div>
</div>
<div id="info" style="display:none;">
    <div>t = <span class="hl" id="i-time">--</span></div>
    <div><span class="hl" id="i-fps">--</span> fps</div>
    <div><span class="hl" id="i-particles">--</span> particles</div>
</div>
<div id="colorbar-wrap" style="display:none;">
    <div class="cb-title">Flow magnitude</div>
    <canvas id="colorbar" width="180" height="14"></canvas>
    <div class="cb-labels"><span id="cb-min">0</span><span id="cb-max">--</span></div>
</div>

<script>
// Minimal 2D particle viewer — for full 3D support, provide the FLUXOS viewer template
document.getElementById('loading').textContent =
    'Minimal viewer loaded. For full 3D support with terrain, provide the FLUXOS viewer template. ' +
    'See README for details.';
</script>
</body>
</html>
'''


# ═══════════════════════════════════════════════════════════════════════════════
#  MAIN DRIVER FUNCTION
# ═══════════════════════════════════════════════════════════════════════════════

def WebGL_h5_driver(shpfile_info=None,
                    openwq_results=None,
                    hydromodel_info=None,
                    hydromodel_var2print=None,
                    output_dir=None,
                    chemSpec=None,
                    file_extension='main',
                    timeframes=None,
                    hostmodel=None,
                    what2map=None,
                    sediment_as_well=False,
                    grid_resolution=None,
                    dem_path=None,
                    n_particles=65536,
                    river_width_cells=3,
                    step=1,
                    satellite_resolution=1,
                    fluxos_viewer_path=None):
    """
    Export OpenWQ/mizuRoute results as an interactive WebGL 3D particle viewer.

    Creates a directory with index.html + data files (velocity PNGs, concentration
    PNGs, heightmap, hillshade, metadata) that can be served via any HTTP server.

    River network segments are rasterized onto a regular grid with velocity vectors
    computed from polyline tangent direction × flow magnitude. GPU-based particle
    advection (adapted from the FLUXOS viewer) visualizes flow and concentrations.

    Parameters
    ----------
    shpfile_info : dict
        Dictionary with 'path_to_shp' and 'mapping_key' keys
    openwq_results : dict
        Results dictionary from Read_h5_driver (required if what2map='openwq')
    hydromodel_info : dict
        Dictionary with 'path_to_results' and 'mapping_key' keys
    hydromodel_var2print : str
        Variable name for flow data (e.g., 'DWroutedRunoff')
    output_dir : str
        Output directory for WebGL export
    chemSpec : list or None
        List of chemical species to export (only for what2map='openwq')
    file_extension : str
        File extension identifier (default: 'main')
    timeframes : int or None
        Number of equally spaced frames (default: all)
    hostmodel : str
        Host model name: 'mizuroute' or 'summa'
    what2map : str
        What to export: 'openwq' or 'hostmodel'
    sediment_as_well : bool
        If True, also export sediment transport results
    grid_resolution : float or None
        Grid cell size in degrees (auto-compute if None)
    dem_path : str or None
        Optional path to DEM .asc or .tif file for 3D terrain
    n_particles : int
        Default particle count for the viewer (default: 65536)
    river_width_cells : int
        Width of rasterized river in grid cells (default: 3)
    step : int
        Subsample every Nth timestep (default: 1 = all)
    satellite_resolution : int
        Satellite image resolution multiplier relative to the grid
        (1 = same as grid, 2 = double, 4 = high quality). Default: 1.
    fluxos_viewer_path : str or None
        Path to FLUXOS index.html viewer template. If None, searches common locations.

    Returns
    -------
    dict : {'output_dir', 'n_frames', 'metadata_path', 'html_path'}
    """

    print("=" * 70)
    print("WEBGL 3D INTERACTIVE RIVER VIEWER EXPORT")
    print(f"Mapping: {what2map.upper()}")
    if what2map == 'hostmodel':
        print(f"Model: {hostmodel.upper()}")
    print("=" * 70)

    # ── STEP 1: Load shapefile ───────────────────────────────
    print("\n" + "=" * 70)
    print("STEP 1: Loading shapefile...")
    print("=" * 70)

    shpfile_fullpath = shpfile_info["path_to_shp"]
    geodf = _load_shapefile_robust(shpfile_fullpath)

    shapefile_mappingKey = shpfile_info["mapping_key"]
    shpf_mapKey = np.array(geodf[shapefile_mappingKey])

    geom_types = geodf.geometry.geom_type.unique()
    print(f"  Loaded {len(geodf)} features")
    print(f"  Geometry types: {', '.join(geom_types)}")

    # ── STEP 2: Compute grid parameters ──────────────────────
    print("\n" + "=" * 70)
    print("STEP 2: Computing grid parameters...")
    print("=" * 70)

    grid_params = _compute_grid_params(geodf, grid_resolution)
    nrows = grid_params['nrows']
    ncols = grid_params['ncols']
    cs = grid_params['cellsize']

    print(f"  Grid: {ncols} x {nrows} cells")
    print(f"  Cell size: {cs:.6f} degrees")
    print(f"  Extent: ({grid_params['xmin']:.4f}, {grid_params['ymin']:.4f}) → "
          f"({grid_params['xmax']:.4f}, {grid_params['ymax']:.4f})")

    # ── STEP 3: Rasterize segments ───────────────────────────
    print("\n" + "=" * 70)
    print("STEP 3: Rasterizing river segments...")
    print("=" * 70)

    river_mask, tangent_x, tangent_y, segment_id_grid = _rasterize_segments(
        geodf, shpf_mapKey, grid_params, river_width_cells
    )

    # ── STEP 4: Load data ────────────────────────────────────
    print("\n" + "=" * 70)
    print("STEP 4: Loading simulation data...")
    print("=" * 70)

    # Flow data (always needed for velocity field)
    flow_ttdata = None
    flow_data_keys = None

    if hydromodel_info is not None:
        hostmodel_nc_path = hydromodel_info['path_to_results']
        hostmodel_mapping_key = hydromodel_info['mapping_key']

        print(f"  NetCDF file: {hostmodel_nc_path}")
        print(f"  Mapping key: {hostmodel_mapping_key}")
        print(f"  Flow variable: {hydromodel_var2print}")

        nc_file = netCDF4.Dataset(hostmodel_nc_path, 'r')
        try:
            segment_ids = np.array(nc_file.variables[hostmodel_mapping_key][:])
            data_var = nc_file.variables[hydromodel_var2print]

            # Get time
            if 'time' in nc_file.variables:
                time_var = nc_file.variables['time']
                time_data = time_var[:]
                if hasattr(time_var, 'units'):
                    from netCDF4 import num2date
                    dates = num2date(time_data, time_var.units,
                                     only_use_cftime_datetimes=False)
                    time_index = pd.DatetimeIndex(dates)
                else:
                    time_index = pd.RangeIndex(len(time_data))
            else:
                time_index = pd.RangeIndex(data_var.shape[0])

            flow_ttdata = pd.DataFrame(
                data_var[:, :],
                index=time_index,
                columns=[f'SEG_{sid}' for sid in segment_ids]
            )

            flow_data_keys = list(flow_ttdata.columns)
            units_flow = (nc_file.variables[hydromodel_var2print].units
                         if hasattr(nc_file.variables[hydromodel_var2print], 'units')
                         else 'units')

            print(f"  Data shape: {flow_ttdata.shape}")
            print(f"  Time range: {time_index[0]} to {time_index[-1]}")
            print(f"  Units: {units_flow}")
        finally:
            nc_file.close()

    # Concentration data (only for what2map='openwq')
    conc_data_per_species = {}  # species_name -> (ttdata, data_keys, units)

    if what2map == 'openwq' and openwq_results is not None:
        # Determine species to export
        if chemSpec is None:
            species_list = list(openwq_results.keys())
        elif isinstance(chemSpec, list):
            species_list = []
            for chem in chemSpec:
                for key in openwq_results.keys():
                    if chem in key:
                        species_list.append(key)
                        break
        else:
            species_list = []
            for key in openwq_results.keys():
                if chemSpec in key:
                    species_list.append(key)
                    break

        # Append sediment if requested
        if sediment_as_well:
            for key in openwq_results.keys():
                if 'Sediment' in key and key not in species_list:
                    species_list.append(key)

        print(f"\n  Chemical species to export: {len(species_list)}")
        for sp in species_list:
            print(f"    - {sp}")

        for species_key in species_list:
            # Parse species name
            parts = species_key.split('@')
            chem_unit = parts[1].split('#') if len(parts) > 1 else [species_key]
            species_name = chem_unit[0]
            units_conc = chem_unit[1] if len(chem_unit) > 1 else 'kg'

            # Extract data
            for ext, data_list in openwq_results[species_key]:
                if ext == file_extension:
                    if data_list and len(data_list) > 0:
                        filename, ttdata, xyz_coords = data_list[0]
                        conc_data_per_species[species_name] = (
                            ttdata, list(ttdata.columns), units_conc
                        )
                        print(f"    {species_name}: shape={ttdata.shape}, units={units_conc}")
                    break

    # ── Determine reference timesteps ────────────────────────
    # Use flow data timesteps as the reference
    if flow_ttdata is not None:
        ref_ttdata = flow_ttdata
    elif conc_data_per_species:
        first_species = list(conc_data_per_species.keys())[0]
        ref_ttdata = conc_data_per_species[first_species][0]
    else:
        print("✗ Error: No data available for export")
        return None

    total_timeframes = len(ref_ttdata)

    # Apply step and timeframes subsampling
    if step > 1:
        indices = list(range(0, total_timeframes, step))
    else:
        indices = list(range(total_timeframes))

    if timeframes is not None and timeframes < len(indices):
        indices = [indices[i] for i in
                   np.linspace(0, len(indices) - 1, timeframes, dtype=int)]

    n_frames = len(indices)
    print(f"\n  Exporting {n_frames} of {total_timeframes} timesteps")

    # ── Create output directories ────────────────────────────
    os.makedirs(output_dir, exist_ok=True)
    data_dir = os.path.join(output_dir, "data")
    vel_dir = os.path.join(data_dir, "velocity")
    os.makedirs(vel_dir, exist_ok=True)

    conc_dirs = {}
    for species_name in conc_data_per_species:
        species_dir = os.path.join(data_dir, "conc", species_name)
        os.makedirs(species_dir, exist_ok=True)
        conc_dirs[species_name] = species_dir

    # ── STEP 5: Pass 1 — compute global ranges ──────────────
    print("\n" + "=" * 70)
    print("STEP 5: Computing global velocity and concentration ranges...")
    print("=" * 70)

    all_vx = []
    all_vy = []
    all_h = []
    all_conc = {sp: [] for sp in conc_data_per_species}

    for idx in indices[::max(1, len(indices) // 20)]:  # Sample every 20th frame
        if flow_ttdata is not None:
            flow_grid = _map_data_to_grid(
                flow_ttdata.iloc[idx].values, flow_data_keys,
                river_mask, segment_id_grid, shpf_mapKey
            )
            grid_vx, grid_vy, wet, grid_h = _build_velocity_field(
                flow_grid, river_mask, tangent_x, tangent_y
            )
            if np.any(wet):
                all_vx.append(grid_vx[wet])
                all_vy.append(grid_vy[wet])
                all_h.append(grid_h[wet])

        for species_name, (ttdata, data_keys, _) in conc_data_per_species.items():
            if idx < len(ttdata):
                conc_grid = _map_data_to_grid(
                    ttdata.iloc[idx].values, data_keys,
                    river_mask, segment_id_grid, shpf_mapKey
                )
                conc_wet = conc_grid[river_mask & (conc_grid > 0)]
                if len(conc_wet) > 0:
                    all_conc[species_name].append(conc_wet)

    # Velocity ranges (99.5th percentile)
    if all_vx:
        cat_vx = np.concatenate(all_vx)
        cat_vy = np.concatenate(all_vy)
        vx_min = float(np.percentile(cat_vx, 0.5))
        vx_max = float(np.percentile(cat_vx, 99.5))
        vy_min = float(np.percentile(cat_vy, 0.5))
        vy_max = float(np.percentile(cat_vy, 99.5))
    else:
        vx_min, vx_max = -0.1, 0.1
        vy_min, vy_max = -0.1, 0.1

    if vx_max - vx_min < 1e-9:
        vx_min, vx_max = -0.1, 0.1
    if vy_max - vy_min < 1e-9:
        vy_min, vy_max = -0.1, 0.1

    if all_h:
        h_max = float(np.percentile(np.concatenate(all_h), 99.9))
    else:
        h_max = 1.0
    if h_max < 1e-6:
        h_max = 1.0
    h_max_raw = h_max  # Keep original for B-channel encoding

    # Concentration ranges per species (95th percentile)
    conc_max_per_species = {}
    for species_name, conc_list in all_conc.items():
        if conc_list:
            cat_conc = np.concatenate(conc_list)
            cat_pos = cat_conc[cat_conc > 0]
            if len(cat_pos) > 0:
                raw_max = float(cat_pos.max())
                threshold = raw_max * 0.01
                plume = cat_pos[cat_pos >= threshold]
                conc_max_per_species[species_name] = float(
                    np.percentile(plume, 95) if len(plume) > 0 else raw_max
                )
            else:
                conc_max_per_species[species_name] = 1.0
        else:
            conc_max_per_species[species_name] = 1.0

    print(f"  Velocity range: vx=[{vx_min:.4f}, {vx_max:.4f}], "
          f"vy=[{vy_min:.4f}, {vy_max:.4f}]")
    print(f"  Max flow magnitude: {h_max:.4f}")
    for sp, cm in conc_max_per_species.items():
        print(f"  Max concentration ({sp}): {cm:.4f}")

    # ── STEP 6: Pass 2 — export PNGs ────────────────────────
    print("\n" + "=" * 70)
    print(f"STEP 6: Exporting {n_frames} velocity + concentration PNGs...")
    print("=" * 70)

    times = []
    total_bytes = 0

    for fi, idx in enumerate(indices):
        # Timestamp
        ts = ref_ttdata.index[idx]
        if isinstance(ts, pd.Timestamp):
            times.append(str(ts))
        else:
            times.append(str(ts))

        # Build velocity field
        if flow_ttdata is not None:
            flow_grid = _map_data_to_grid(
                flow_ttdata.iloc[idx].values, flow_data_keys,
                river_mask, segment_id_grid, shpf_mapKey
            )
        else:
            flow_grid = np.ones_like(river_mask, dtype=np.float64)  # uniform flow

        grid_vx, grid_vy, wet_mask, grid_h = _build_velocity_field(
            flow_grid, river_mask, tangent_x, tangent_y
        )

        # Encode velocity PNG
        vel_bytes = _encode_velocity_png(
            grid_vx, grid_vy, wet_mask, grid_h,
            vx_min, vx_max, vy_min, vy_max, h_max
        )
        fname = f"v_{fi:04d}.png"
        with open(os.path.join(vel_dir, fname), "wb") as f:
            f.write(vel_bytes)
        total_bytes += len(vel_bytes)

        # Encode concentration PNGs per species
        for species_name, (ttdata, data_keys, _) in conc_data_per_species.items():
            if idx < len(ttdata):
                conc_grid = _map_data_to_grid(
                    ttdata.iloc[idx].values, data_keys,
                    river_mask, segment_id_grid, shpf_mapKey
                )
                conc_bytes = _encode_concentration_png(
                    conc_grid,
                    conc_max_per_species.get(species_name, 1.0),
                    wet_mask
                )
                with open(os.path.join(conc_dirs[species_name],
                                        f"c_{fi:04d}.png"), "wb") as f:
                    f.write(conc_bytes)
                total_bytes += len(conc_bytes)

        if fi % 50 == 0 or fi == n_frames - 1:
            print(f"    [{fi+1}/{n_frames}] t={times[-1]}")

    print(f"  Total data: {total_bytes / 1024 / 1024:.1f} MB")

    # ── STEP 7: Export heightmap + hillshade ─────────────────
    print("\n" + "=" * 70)
    print("STEP 7: Exporting heightmap and hillshade...")
    print("=" * 70)

    heightmap_bytes, hillshade_bytes, z_min, z_max = _generate_heightmap(
        grid_params, dem_path, river_mask=river_mask
    )

    with open(os.path.join(data_dir, "heightmap.png"), "wb") as f:
        f.write(heightmap_bytes)
    print(f"  Heightmap: {len(heightmap_bytes)/1024:.0f} KB "
          f"(z_min={z_min:.2f}, z_max={z_max:.2f})")

    # Cap h_max for shader so water offset stays proportional to terrain
    # Shader formula: water_offset = (B * h_max) / z_range * exaggeration
    # With h_max_visual, max offset = 0.15 * exaggeration (15% of terrain)
    z_range = max(z_max - z_min, 0.01)
    water_height_fraction = 0.15
    h_max_visual = z_range * water_height_fraction
    print(f"  Water height: raw h_max={h_max_raw:.2f}, "
          f"visual h_max={h_max_visual:.4f} (z_range={z_range:.2f})")

    with open(os.path.join(data_dir, "hillshade.png"), "wb") as f:
        f.write(hillshade_bytes)
    print(f"  Hillshade: {len(hillshade_bytes)/1024:.0f} KB")

    # ── STEP 7b: Download satellite imagery ───────────────────
    has_satellite = False
    print("\n  Downloading satellite imagery...")
    try:
        sat_path = os.path.join(data_dir, "satellite.jpg")
        sat_mult = max(1, int(satellite_resolution))
        sat_w, sat_h = ncols * sat_mult, nrows * sat_mult
        print(f"  Satellite resolution: {sat_w}x{sat_h} px "
              f"({sat_mult}x grid)")
        _download_satellite(
            sw_lat=grid_params['ymin'], sw_lon=grid_params['xmin'],
            ne_lat=grid_params['ymax'], ne_lon=grid_params['xmax'],
            img_w=sat_w, img_h=sat_h,
            output_path=sat_path,
            river_mask=river_mask if sat_mult == 1 else None
        )
        # Check if file was created (could be .jpg or .png fallback)
        if os.path.exists(sat_path):
            has_satellite = True
        elif os.path.exists(sat_path.replace(".jpg", ".png")):
            has_satellite = True
    except Exception as e:
        print(f"  ⚠ Satellite download failed: {e}")
        print("  Continuing without satellite imagery...")

    # ── STEP 8: Export metadata JSON ─────────────────────────
    print("\n" + "=" * 70)
    print("STEP 8: Exporting metadata...")
    print("=" * 70)

    # Compute bounding box in lat/lon
    bbox_ll = {
        "sw_lat": float(grid_params['ymin']),
        "sw_lon": float(grid_params['xmin']),
        "ne_lat": float(grid_params['ymax']),
        "ne_lon": float(grid_params['xmax']),
    }

    # Compute speed factor for particle advection
    # The shader uses: offset = vel * speed_factor * user_speed / grid_res
    # We want max velocity particles to move ~3 cells per frame at 1x speed
    max_vel_magnitude = max(abs(vx_max), abs(vx_min), abs(vy_max), abs(vy_min), 1e-6)
    target_cells_per_frame = 3.0
    computed_speed_factor = target_cells_per_frame / max_vel_magnitude

    metadata = {
        "width": int(ncols),
        "height": int(nrows),
        "n_timesteps": n_frames,
        "times": times,
        "time_format": "datetime",
        # NOTE: The FLUXOS shader cross-reads channels:
        #   R channel (northing/vy) decoded with u_vel_range_x → metadata.vx_min/max
        #   G channel (easting/vx) decoded with u_vel_range_y → metadata.vy_min/max
        # So we swap: vx_min/max = vy range (R channel), vy_min/max = vx range (G channel)
        "vx_min": float(vy_min),  # R channel range = northing velocity
        "vx_max": float(vy_max),
        "vy_min": float(vx_min),  # G channel range = easting velocity
        "vy_max": float(vx_max),
        "cellsize": float(cs),
        "speed_factor": float(computed_speed_factor),
        "default_particles": n_particles,
        "z_min": z_min,
        "z_max": z_max,
        "h_max": float(h_max_visual),
        "height_exaggeration": 0.1,
        "bbox": bbox_ll,
        "model": hostmodel or "unknown",
        "flow_variable": hydromodel_var2print or "flow",
        "has_satellite": has_satellite,
    }

    if conc_data_per_species:
        metadata["has_concentration"] = True
        species_names = list(conc_data_per_species.keys())
        metadata["chemical_species"] = species_names
        metadata["conc_max"] = {
            sp: conc_max_per_species.get(sp, 1.0) for sp in species_names
        }
        # Units per species
        metadata["conc_units"] = {
            sp: conc_data_per_species[sp][2] for sp in species_names
        }

    meta_path = os.path.join(data_dir, "metadata.json")
    with open(meta_path, "w") as f:
        _json.dump(metadata, f, indent=2)
    print(f"  Metadata: {meta_path}")

    # ── STEP 9: Generate HTML viewer ─────────────────────────
    print("\n" + "=" * 70)
    print("STEP 9: Generating WebGL viewer...")
    print("=" * 70)

    html = _generate_openwq_webgl_html(fluxos_viewer_path)

    # Inline metadata JSON into the HTML so the viewer works on file://
    # (fetch() is blocked by CORS on file:// protocol).
    # The template checks for window.__OPENWQ_METADATA__ before fetch().
    _meta_json_str = _json.dumps(metadata)
    _inline_script = (
        f'\n<script>window.__OPENWQ_METADATA__ = {_meta_json_str};</script>\n'
    )
    html = html.replace('</head>', f'{_inline_script}</head>', 1)

    html_path = os.path.join(output_dir, "index.html")
    with open(html_path, "w") as f:
        f.write(html)

    # ── SUMMARY ──────────────────────────────────────────────
    print("\n" + "=" * 70)
    print("✓ WEBGL VIEWER EXPORT COMPLETE!")
    print("=" * 70)
    print(f"\n  Output directory: {output_dir}/")
    print(f"  Velocity frames: {n_frames}")
    if conc_data_per_species:
        print(f"  Concentration species: {list(conc_data_per_species.keys())}")
    print(f"  Grid: {ncols}x{nrows} cells ({cs:.6f}° resolution)")
    print(f"  Total data: {total_bytes / 1024 / 1024:.1f} MB")
    print(f"\n  To view, start a local server:")
    print(f"    python3 -m http.server 8080 -d {output_dir}")
    print(f"  Then open:  http://localhost:8080")
    print("=" * 70)

    return {
        'output_dir': output_dir,
        'n_frames': n_frames,
        'metadata_path': meta_path,
        'html_path': html_path,
    }
