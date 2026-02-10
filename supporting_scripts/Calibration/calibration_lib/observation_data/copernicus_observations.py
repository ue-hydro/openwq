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
Copernicus Observation Generator for OpenWQ Calibration
========================================================

Generates synthetic water quality observations using Copernicus land cover data
and optional climate data. Uses export coefficients to estimate nutrient loads
based on land use within each catchment.

This is useful for:
- Calibration in data-sparse regions
- Sensitivity analysis with controlled observation uncertainty
- Testing calibration workflows before using real observations
- Regions where GRQA data is not available

Usage:
    from copernicus_observations import CopernicusObservationGenerator

    generator = CopernicusObservationGenerator(
        river_network_shapefile="river_network.shp",
        land_cover_path="copernicus_lc.tif",
        output_dir="output/"
    )

    result = generator.generate(
        target_species=["TN", "TP"],
        export_coefficients={10: {"TN": 25.0, "TP": 2.0}},
        start_date="2010-01-01",
        end_date="2020-12-31"
    )
"""

import os
import json
import numpy as np
import pandas as pd
from datetime import datetime, timedelta
from typing import Dict, List, Optional, Tuple, Any

try:
    import geopandas as gpd
    from shapely.geometry import Point
    HAS_GEOPANDAS = True
except ImportError:
    HAS_GEOPANDAS = False

try:
    import rasterio
    from rasterio.mask import mask as rasterio_mask
    HAS_RASTERIO = True
except ImportError:
    HAS_RASTERIO = False

try:
    import xarray as xr
    HAS_XARRAY = True
except ImportError:
    HAS_XARRAY = False


class CopernicusObservationGenerator:
    """
    Generates synthetic water quality observations from Copernicus data.

    Uses land cover-based export coefficients to estimate nutrient concentrations
    at each river reach, optionally adjusted by climate data.
    """

    # Default export coefficients (kg/ha/yr) by Copernicus land cover class
    DEFAULT_EXPORT_COEFFICIENTS = {
        # Cropland
        10: {"TN": 25.0, "TP": 2.0, "TSS": 500.0, "DOC": 15.0},
        # Forest (all types)
        20: {"TN": 2.0, "TP": 0.2, "TSS": 50.0, "DOC": 8.0},
        # Grassland
        30: {"TN": 10.0, "TP": 0.8, "TSS": 150.0, "DOC": 10.0},
        # Shrubland
        40: {"TN": 5.0, "TP": 0.4, "TSS": 100.0, "DOC": 6.0},
        # Wetland
        50: {"TN": 3.0, "TP": 0.3, "TSS": 30.0, "DOC": 20.0},
        # Water bodies
        60: {"TN": 0.5, "TP": 0.05, "TSS": 10.0, "DOC": 5.0},
        # Tundra
        70: {"TN": 1.0, "TP": 0.1, "TSS": 20.0, "DOC": 3.0},
        # Artificial/Urban
        80: {"TN": 15.0, "TP": 1.5, "TSS": 300.0, "DOC": 12.0},
        # Bare/sparse vegetation
        90: {"TN": 1.0, "TP": 0.1, "TSS": 200.0, "DOC": 2.0},
        # Snow/Ice
        100: {"TN": 0.1, "TP": 0.01, "TSS": 5.0, "DOC": 0.5},
    }

    def __init__(
        self,
        river_network_shapefile: str,
        reach_id_column: str = "seg_id",
        land_cover_path: Optional[str] = None,
        climate_data_path: Optional[str] = None,
        output_dir: str = "."
    ):
        """
        Initialize the Copernicus observation generator.

        Args:
            river_network_shapefile: Path to river network shapefile
            reach_id_column: Column name containing reach IDs
            land_cover_path: Path to Copernicus land cover raster (GeoTIFF or NetCDF)
                            Can be a single file or a directory containing yearly files
            climate_data_path: Optional path to climate data (NetCDF)
                              Can be a single file or directory
            output_dir: Output directory for generated files
        """
        if not HAS_GEOPANDAS:
            raise ImportError("geopandas is required: pip install geopandas")

        self.river_network_shapefile = river_network_shapefile
        self.reach_id_column = reach_id_column
        self.land_cover_path = land_cover_path
        self.climate_data_path = climate_data_path
        self.output_dir = output_dir

        os.makedirs(output_dir, exist_ok=True)

        # Load river network
        print(f"Loading river network from {river_network_shapefile}...")
        self.river_network = gpd.read_file(river_network_shapefile)
        print(f"  Loaded {len(self.river_network)} reaches")

        # Load land cover if provided
        self.land_cover = None
        self.land_cover_files = {}  # For multi-year data
        if land_cover_path:
            self._load_land_cover(land_cover_path)

        # Load climate data if provided
        self.climate_data = None
        self.climate_files = {}  # For multi-file data
        if climate_data_path:
            self._load_climate_data(climate_data_path)

    def _load_land_cover(self, path: str):
        """Load land cover data from file or directory."""
        path = Path(path)

        if not HAS_RASTERIO:
            print("  Warning: rasterio not available, using default land cover fractions")
            return

        if path.is_file():
            # Single file
            print(f"Loading land cover from {path}...")
            self.land_cover = rasterio.open(str(path))
            print(f"  Raster size: {self.land_cover.width} x {self.land_cover.height}")

        elif path.is_dir():
            # Directory with multiple files
            print(f"Scanning land cover directory: {path}")
            for ext in ['*.tif', '*.tiff', '*.nc']:
                for f in path.glob(ext):
                    # Try to extract year from filename
                    name = f.stem
                    for year in range(1990, 2030):
                        if str(year) in name:
                            self.land_cover_files[year] = f
                            break

            if self.land_cover_files:
                print(f"  Found {len(self.land_cover_files)} yearly land cover files")
                # Load first year as default
                first_year = min(self.land_cover_files.keys())
                self.land_cover = rasterio.open(str(self.land_cover_files[first_year]))
            else:
                print(f"  Warning: No land cover files found in {path}")

    def _load_climate_data(self, path: str):
        """Load climate data from file or directory."""
        path = Path(path)

        if not HAS_XARRAY:
            print("  Warning: xarray not available, climate adjustment disabled")
            return

        if path.is_file():
            # Single file
            print(f"Loading climate data from {path}...")
            self.climate_data = xr.open_dataset(str(path))

        elif path.is_dir():
            # Directory with multiple files
            print(f"Scanning climate data directory: {path}")
            nc_files = list(path.glob('*.nc'))
            if nc_files:
                print(f"  Found {len(nc_files)} NetCDF files")
                # Use xarray to open multiple files
                try:
                    self.climate_data = xr.open_mfdataset(
                        str(path / '*.nc'),
                        combine='by_coords'
                    )
                except Exception as e:
                    print(f"  Warning: Could not combine climate files: {e}")
                    # Load first file only
                    self.climate_data = xr.open_dataset(str(nc_files[0]))

    def _get_land_cover_fractions(self, geometry) -> Dict[int, float]:
        """
        Calculate land cover class fractions for a catchment geometry.

        Args:
            geometry: Shapely geometry of catchment

        Returns:
            Dictionary of {land_cover_class: fraction}
        """
        if self.land_cover is None:
            # Return default (assume mixed land use)
            return {10: 0.3, 20: 0.4, 30: 0.2, 80: 0.1}

        try:
            # Mask raster by geometry
            out_image, _ = rasterio_mask(
                self.land_cover, [geometry], crop=True, nodata=0
            )
            lc_data = out_image[0]

            # Count pixels per class
            unique, counts = np.unique(lc_data[lc_data > 0], return_counts=True)
            total = counts.sum()

            if total == 0:
                return {10: 0.3, 20: 0.4, 30: 0.2, 80: 0.1}

            return {int(u): c / total for u, c in zip(unique, counts)}

        except Exception:
            return {10: 0.3, 20: 0.4, 30: 0.2, 80: 0.1}

    def _get_climate_factor(
        self,
        geometry,
        date: datetime,
        climate_params: Dict[str, float]
    ) -> float:
        """
        Calculate climate adjustment factor for a given location and date.

        Args:
            geometry: Location geometry
            date: Date for climate data
            climate_params: Climate response parameters

        Returns:
            Climate adjustment factor (1.0 = no adjustment)
        """
        if self.climate_data is None:
            return 1.0

        try:
            # Get centroid
            centroid = geometry.centroid

            # Extract climate variables at location
            precip = self.climate_data.sel(
                time=date,
                lat=centroid.y,
                lon=centroid.x,
                method="nearest"
            ).get("tp", 1.0)  # Total precipitation

            temp = self.climate_data.sel(
                time=date,
                lat=centroid.y,
                lon=centroid.x,
                method="nearest"
            ).get("t2m", 288.0)  # 2m temperature (K)

            # Convert temperature to Celsius
            temp_c = float(temp) - 273.15

            # Calculate precipitation factor
            precip_power = climate_params.get("precip_scaling_power", 1.0)
            precip_factor = (float(precip) / 0.003) ** precip_power  # Normalize to ~3mm/day

            # Calculate temperature factor (Q10)
            q10 = climate_params.get("temp_q10", 2.0)
            t_ref = climate_params.get("temp_reference_c", 15.0)
            temp_factor = q10 ** ((temp_c - t_ref) / 10.0)

            # Combined factor
            return np.clip(precip_factor * temp_factor, 0.1, 10.0)

        except Exception:
            return 1.0

    def _calculate_concentration(
        self,
        reach_id: int,
        geometry,
        species: str,
        date: datetime,
        export_coefficients: Dict[int, Dict[str, float]],
        climate_params: Optional[Dict[str, float]] = None,
        catchment_area_km2: float = 100.0,
        flow_m3s: float = 10.0
    ) -> float:
        """
        Calculate estimated concentration for a reach.

        Args:
            reach_id: Reach identifier
            geometry: Reach geometry
            species: Chemical species name
            date: Date for the observation
            export_coefficients: Land cover export coefficients
            climate_params: Climate response parameters
            catchment_area_km2: Catchment area in km²
            flow_m3s: River flow in m³/s

        Returns:
            Estimated concentration in mg/L
        """
        # Get land cover fractions
        lc_fractions = self._get_land_cover_fractions(geometry)

        # Calculate weighted export coefficient
        total_export = 0.0  # kg/ha/yr
        for lc_class, fraction in lc_fractions.items():
            if lc_class in export_coefficients:
                coeff = export_coefficients[lc_class].get(species, 0.0)
            elif lc_class in self.DEFAULT_EXPORT_COEFFICIENTS:
                coeff = self.DEFAULT_EXPORT_COEFFICIENTS[lc_class].get(species, 0.0)
            else:
                coeff = 0.0
            total_export += fraction * coeff

        # Apply climate adjustment
        if climate_params:
            climate_factor = self._get_climate_factor(geometry, date, climate_params)
            total_export *= climate_factor

        # Convert to concentration
        # Export: kg/ha/yr, Area: km² = 100 ha, Flow: m³/s
        # Load = export * area * 100 ha/km² = kg/yr
        # Concentration = Load / (flow * seconds_per_year) * 1000 mg/kg / 1000 L/m³
        # Simplified: concentration (mg/L) = export * area * 100 / (flow * 31536000) * 1000

        load_kg_yr = total_export * catchment_area_km2 * 100  # kg/yr
        flow_m3_yr = flow_m3s * 31536000  # m³/yr
        concentration_mg_l = (load_kg_yr * 1e6) / flow_m3_yr  # mg/L

        return max(0.001, concentration_mg_l)  # Minimum detection limit

    def generate(
        self,
        target_species: List[str],
        export_coefficients: Optional[Dict[int, Dict[str, float]]] = None,
        start_date: str = "2010-01-01",
        end_date: str = "2020-12-31",
        temporal_resolution: str = "monthly",
        climate_params: Optional[Dict[str, float]] = None,
        uncertainty_fraction: float = 0.3,
        add_noise: bool = True,
        noise_cv: float = 0.2
    ) -> Dict[str, Any]:
        """
        Generate synthetic observations for calibration.

        Args:
            target_species: List of species to generate
            export_coefficients: Land cover export coefficients (kg/ha/yr)
            start_date: Start date (ISO format)
            end_date: End date (ISO format)
            temporal_resolution: "daily", "weekly", "monthly", or "annual"
            climate_params: Climate response parameters
            uncertainty_fraction: Uncertainty as fraction of value
            add_noise: Whether to add random noise
            noise_cv: Coefficient of variation for noise

        Returns:
            Dictionary with results and output paths
        """
        if export_coefficients is None:
            export_coefficients = self.DEFAULT_EXPORT_COEFFICIENTS

        print(f"\nGenerating observations for species: {target_species}")
        print(f"Period: {start_date} to {end_date} ({temporal_resolution})")

        # Generate date range
        start = pd.to_datetime(start_date)
        end = pd.to_datetime(end_date)

        if temporal_resolution == "daily":
            dates = pd.date_range(start, end, freq="D")
        elif temporal_resolution == "weekly":
            dates = pd.date_range(start, end, freq="W")
        elif temporal_resolution == "monthly":
            dates = pd.date_range(start, end, freq="MS")
        elif temporal_resolution == "annual":
            dates = pd.date_range(start, end, freq="YS")
        else:
            dates = pd.date_range(start, end, freq="MS")

        print(f"Generating {len(dates)} time points for {len(self.river_network)} reaches...")

        # Generate observations
        observations = []
        np.random.seed(42)  # Reproducibility

        for idx, row in self.river_network.iterrows():
            reach_id = row[self.reach_id_column]
            geometry = row.geometry

            # Get catchment properties (use defaults if not available)
            catchment_area = row.get("catchment_area_km2", 100.0)
            if pd.isna(catchment_area) or catchment_area <= 0:
                catchment_area = 100.0

            flow = row.get("mean_flow_m3s", 10.0)
            if pd.isna(flow) or flow <= 0:
                flow = 10.0

            for date in dates:
                for species in target_species:
                    # Calculate base concentration
                    conc = self._calculate_concentration(
                        reach_id=reach_id,
                        geometry=geometry,
                        species=species,
                        date=date.to_pydatetime(),
                        export_coefficients=export_coefficients,
                        climate_params=climate_params,
                        catchment_area_km2=catchment_area,
                        flow_m3s=flow
                    )

                    # Add noise if requested
                    if add_noise:
                        noise = np.random.normal(1.0, noise_cv)
                        conc *= max(0.1, noise)

                    # Calculate uncertainty
                    uncertainty = conc * uncertainty_fraction

                    observations.append({
                        "datetime": date.strftime("%Y-%m-%d %H:%M:%S"),
                        "reach_id": reach_id,
                        "species": species,
                        "value": round(conc, 4),
                        "units": "mg/l",
                        "source": "Copernicus_synthetic",
                        "uncertainty": round(uncertainty, 4),
                        "quality_flag": "SYNTHETIC"
                    })

        # Create DataFrame
        obs_df = pd.DataFrame(observations)
        print(f"Generated {len(obs_df)} observations")

        # Save to CSV
        output_path = os.path.join(self.output_dir, "calibration_observations.csv")
        obs_df.to_csv(output_path, index=False)
        print(f"Saved to: {output_path}")

        # Save metadata
        metadata = {
            "source": "Copernicus land cover export coefficient model",
            "start_date": start_date,
            "end_date": end_date,
            "temporal_resolution": temporal_resolution,
            "target_species": target_species,
            "reaches_covered": len(self.river_network),
            "total_observations": len(obs_df),
            "uncertainty_fraction": uncertainty_fraction,
            "noise_added": add_noise,
            "noise_cv": noise_cv if add_noise else None,
            "export_coefficients": {
                str(k): v for k, v in export_coefficients.items()
            },
            "climate_adjusted": climate_params is not None,
            "generated_at": datetime.now().isoformat()
        }

        metadata_path = os.path.join(self.output_dir, "copernicus_metadata.json")
        with open(metadata_path, "w") as f:
            json.dump(metadata, f, indent=2)

        return {
            "calibration_observations_path": output_path,
            "metadata_path": metadata_path,
            "total_observations": len(obs_df),
            "reaches_covered": len(self.river_network),
            "species_list": target_species,
            "start_date": start_date,
            "end_date": end_date
        }


def main():
    """Command-line interface for Copernicus observation generator."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Generate synthetic observations from Copernicus data"
    )
    parser.add_argument(
        "--river-network", "-r",
        required=True,
        help="Path to river network shapefile"
    )
    parser.add_argument(
        "--land-cover", "-l",
        help="Path to Copernicus land cover raster"
    )
    parser.add_argument(
        "--output-dir", "-o",
        default=".",
        help="Output directory"
    )
    parser.add_argument(
        "--species",
        nargs="+",
        default=["TN", "TP", "TSS"],
        help="Species to generate"
    )
    parser.add_argument(
        "--start-date",
        default="2010-01-01",
        help="Start date (ISO format)"
    )
    parser.add_argument(
        "--end-date",
        default="2020-12-31",
        help="End date (ISO format)"
    )
    parser.add_argument(
        "--resolution",
        choices=["daily", "weekly", "monthly", "annual"],
        default="monthly",
        help="Temporal resolution"
    )

    args = parser.parse_args()

    generator = CopernicusObservationGenerator(
        river_network_shapefile=args.river_network,
        land_cover_path=args.land_cover,
        output_dir=args.output_dir
    )

    result = generator.generate(
        target_species=args.species,
        start_date=args.start_date,
        end_date=args.end_date,
        temporal_resolution=args.resolution
    )

    print("\nGeneration complete!")
    print(f"Output: {result['calibration_observations_path']}")


if __name__ == "__main__":
    main()
