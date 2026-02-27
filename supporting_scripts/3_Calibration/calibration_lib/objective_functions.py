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
Objective Functions Module
==========================

Computes goodness-of-fit metrics by comparing model outputs with observations.

Supports temporal aggregation of both observations and simulations to different
resolutions (native, daily, weekly, monthly, yearly) before computing metrics.
"""

import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Union, Literal
import logging
import sys

logger = logging.getLogger(__name__)

# Valid temporal resolution options
TemporalResolution = Literal["native", "daily", "weekly", "monthly", "yearly"]


class ObjectiveFunction:
    """
    Computes objective function values for calibration.
    """

    def __init__(self,
                 observation_path: str,
                 target_species: List[str],
                 target_reaches: Union[List[int], str],
                 compartments: List[str],
                 metric: str = "KGE",
                 weights: Dict[str, float] = None,
                 no_data_flag: float = -9999,
                 h5_reader_path: str = None,
                 temporal_resolution: TemporalResolution = "native",
                 aggregation_method: str = "mean"):
        """
        Initialize objective function calculator.

        Parameters
        ----------
        observation_path : str
            Path to observation CSV file
        target_species : List[str]
            Chemical species to include in objective
        target_reaches : List[int] or "all"
            Reach IDs to include in objective
        compartments : List[str]
            Compartment names in model output
        metric : str
            Objective function type: "RMSE", "NSE", "KGE"
        weights : Dict[str, float]
            Species weights for multi-species objective
        no_data_flag : float
            Flag value for missing data
        h5_reader_path : str
            Path to the Read_h5_driver.py module
        temporal_resolution : str
            Temporal resolution for objective function calculation:
            - "native": Use original model timestep (no aggregation)
            - "daily": Aggregate to daily values
            - "weekly": Aggregate to weekly values
            - "monthly": Aggregate to monthly values
            - "yearly": Aggregate to yearly values
        aggregation_method : str
            Method for temporal aggregation: "mean", "sum", "median", "min", "max"
        """
        self.observation_path = observation_path
        self.target_species = target_species
        self.target_reaches = target_reaches
        self.compartments = compartments
        self.metric = metric
        self.weights = weights or {s: 1.0 for s in target_species}
        self.no_data_flag = no_data_flag
        self.temporal_resolution = temporal_resolution
        self.aggregation_method = aggregation_method

        # Validate temporal resolution
        valid_resolutions = ["native", "daily", "weekly", "monthly", "yearly"]
        if temporal_resolution not in valid_resolutions:
            raise ValueError(
                f"Invalid temporal_resolution '{temporal_resolution}'. "
                f"Valid options: {valid_resolutions}"
            )

        # Validate aggregation method
        valid_methods = ["mean", "sum", "median", "min", "max"]
        if aggregation_method not in valid_methods:
            raise ValueError(
                f"Invalid aggregation_method '{aggregation_method}'. "
                f"Valid options: {valid_methods}"
            )

        logger.info(f"Objective function configured with temporal_resolution='{temporal_resolution}', "
                    f"aggregation_method='{aggregation_method}'")

        # Setup H5 reader import path
        if h5_reader_path:
            sys.path.insert(0, str(Path(h5_reader_path)))

        self.observations = self._load_observations()

    def _load_observations(self) -> pd.DataFrame:
        """Load and preprocess observation data."""
        try:
            obs = pd.read_csv(self.observation_path, parse_dates=['datetime'])
        except FileNotFoundError:
            logger.warning(f"Observation file not found: {self.observation_path}")
            return pd.DataFrame()

        # Filter by species
        obs = obs[obs['species'].isin(self.target_species)]

        # Filter by reaches
        if self.target_reaches != "all" and isinstance(self.target_reaches, list):
            obs = obs[obs['reach_id'].isin(self.target_reaches)]

        logger.info(f"Loaded {len(obs)} observations for {obs['species'].nunique()} species")
        return obs

    def compute(self,
                output_dir: Path,
                units: str = "MG/L") -> float:
        """
        Compute objective function value.

        Parameters
        ----------
        output_dir : Path
            Path to openwq_out/ directory
        units : str
            Concentration units in model output

        Returns
        -------
        float
            Objective function value (lower is better for minimization)
        """
        if self.observations.empty:
            logger.warning("No observations loaded - returning penalty")
            return 1e10

        # Extract simulated values
        simulated = self._extract_simulated(output_dir, units)

        if simulated.empty:
            logger.warning("No simulated values extracted - returning penalty")
            return 1e10

        # Match observations with simulated
        matched = self._match_obs_sim(simulated)

        if matched.empty:
            logger.warning("No matched observation-simulation pairs - returning penalty")
            return 1e10

        # Apply temporal aggregation if specified
        if self.temporal_resolution != "native":
            matched = self._aggregate_matched_data(matched)

            if matched.empty:
                logger.warning("No data after temporal aggregation - returning penalty")
                return 1e10

        # Store matched data for later analysis/plotting
        self._last_matched_data = matched

        # Compute objective by species
        objectives = {}
        for species in self.target_species:
            species_data = matched[matched['species'] == species]
            if species_data.empty:
                continue

            obs_vals = species_data['observed'].values
            sim_vals = species_data['simulated'].values

            # Remove NaN pairs
            valid_mask = ~(np.isnan(obs_vals) | np.isnan(sim_vals))
            obs_vals = obs_vals[valid_mask]
            sim_vals = sim_vals[valid_mask]

            if len(obs_vals) < 2:
                continue

            if self.metric == "RMSE":
                obj = self.rmse(obs_vals, sim_vals)
            elif self.metric == "NSE":
                obj = self.nse_minimization(obs_vals, sim_vals)
            elif self.metric == "KGE":
                obj = self.kge_minimization(obs_vals, sim_vals)
            else:
                obj = self.rmse(obs_vals, sim_vals)

            objectives[species] = obj

        if not objectives:
            return 1e10

        # Weighted average
        total_weight = sum(self.weights.get(s, 1.0) for s in objectives.keys())
        weighted_obj = sum(
            self.weights.get(s, 1.0) * obj
            for s, obj in objectives.items()
        ) / total_weight

        logger.debug(f"Objective: {weighted_obj:.6f} (species objectives: {objectives})")
        return weighted_obj

    def _extract_simulated(self,
                           output_dir: Path,
                           units: str) -> pd.DataFrame:
        """
        Read model outputs and extract simulated values at observation points.
        """
        try:
            # Try to import the HDF5 reader
            from hdf5_support_lib import Read_h5_driver
            h5_lib = Read_h5_driver
        except ImportError:
            # Fallback: direct h5py reading
            return self._extract_simulated_direct(output_dir, units)

        openwq_info = {
            "path_to_results": str(output_dir),
            "mapping_key": "reachID"
        }

        try:
            results = h5_lib.Read_h5_driver(
                openwq_info=openwq_info,
                output_format='HDF5',
                debugmode=False,
                cmp=self.compartments,
                space_elem='all',
                chemSpec=self.target_species,
                chemUnits=units,
                noDataFlag=self.no_data_flag
            )

            # Convert results to DataFrame
            all_data = []
            for key, extensions in results.items():
                # Parse key: "COMPARTMENT@SPECIES#UNITS"
                parts = key.split('@')
                if len(parts) < 2:
                    continue
                compartment = parts[0]
                species_units = parts[1].split('#')
                species = species_units[0]

                # Get main results (first extension)
                for ext_name, data_list in extensions:
                    if ext_name != 'main':
                        continue
                    for filename, df, coords in data_list:
                        if df is None or df.empty:
                            continue
                        # df has time index and columns like "reachID_123"
                        for col in df.columns:
                            reach_id = col.replace("reachID_", "")
                            try:
                                reach_id = int(reach_id)
                            except ValueError:
                                continue

                            for dt, val in df[col].items():
                                all_data.append({
                                    'datetime': dt,
                                    'reach_id': reach_id,
                                    'species': species,
                                    'simulated': val
                                })

            return pd.DataFrame(all_data)

        except Exception as e:
            logger.warning(f"Error reading HDF5 with library: {e}")
            return self._extract_simulated_direct(output_dir, units)

    def _extract_simulated_direct(self,
                                  output_dir: Path,
                                  units: str) -> pd.DataFrame:
        """
        Direct HDF5 reading without the Read_h5_driver library.
        """
        import h5py

        all_data = []
        h5_dir = output_dir / "HDF5"

        for species in self.target_species:
            for compartment in self.compartments:
                # Build filename
                filename = f"{compartment}@{species}#{units.replace('/', '|')}-main.h5"
                filepath = h5_dir / filename

                if not filepath.exists():
                    # Try uppercase
                    filename_upper = filename.upper()
                    filepath = h5_dir / filename_upper
                    if not filepath.exists():
                        continue

                try:
                    with h5py.File(filepath, 'r') as hf:
                        # Get reach IDs
                        if '/reachID' in hf:
                            reach_ids = [x.decode() if isinstance(x, bytes) else x
                                        for x in hf['/reachID'][:]]
                        else:
                            reach_ids = list(range(hf[list(hf.keys())[0]].shape[1]))

                        # Get timestamps
                        timestamps = [k for k in hf.keys()
                                     if k not in ['xyz_elements', 'reachID', 'mapping_key']]

                        for ts in timestamps:
                            try:
                                dt = pd.to_datetime(ts, format='%Y%b%d-%H:%M:%S')
                            except:
                                try:
                                    dt = pd.to_datetime(ts)
                                except:
                                    continue

                            data = hf[f'/{ts}'][:]
                            if data.ndim > 1:
                                data = data[0, :]

                            for i, reach_id in enumerate(reach_ids):
                                try:
                                    rid = int(reach_id)
                                except:
                                    continue

                                val = data[i] if i < len(data) else np.nan
                                if val == self.no_data_flag:
                                    val = np.nan

                                all_data.append({
                                    'datetime': dt,
                                    'reach_id': rid,
                                    'species': species,
                                    'simulated': val
                                })

                except Exception as e:
                    logger.warning(f"Error reading {filepath}: {e}")
                    continue

        return pd.DataFrame(all_data)

    def _get_temporal_grouper(self, resolution: str) -> str:
        """
        Get pandas frequency string for temporal aggregation.

        Parameters
        ----------
        resolution : str
            Temporal resolution: "daily", "weekly", "monthly", "yearly"

        Returns
        -------
        str
            Pandas frequency string for resample/groupby
        """
        freq_map = {
            "daily": "D",
            "weekly": "W",
            "monthly": "MS",  # Month Start for consistent grouping
            "yearly": "YS"    # Year Start
        }
        return freq_map.get(resolution, "D")

    def _aggregate_temporal(self,
                            df: pd.DataFrame,
                            value_column: str = "value",
                            datetime_column: str = "datetime") -> pd.DataFrame:
        """
        Aggregate data to the specified temporal resolution.

        Parameters
        ----------
        df : pd.DataFrame
            DataFrame with datetime, reach_id, species, and value columns
        value_column : str
            Name of the column containing values to aggregate
        datetime_column : str
            Name of the datetime column

        Returns
        -------
        pd.DataFrame
            Aggregated DataFrame with same structure
        """
        if self.temporal_resolution == "native" or df.empty:
            return df

        # Ensure datetime column is datetime type
        df = df.copy()
        df[datetime_column] = pd.to_datetime(df[datetime_column])

        # Get the aggregation frequency
        freq = self._get_temporal_grouper(self.temporal_resolution)

        # Define aggregation function
        agg_func_map = {
            "mean": "mean",
            "sum": "sum",
            "median": "median",
            "min": "min",
            "max": "max"
        }
        agg_func = agg_func_map.get(self.aggregation_method, "mean")

        # Group by reach_id, species, and temporal period
        df['period'] = df[datetime_column].dt.to_period(freq[0] if freq != "MS" else "M")

        # Handle different frequencies
        if self.temporal_resolution == "weekly":
            # Use isocalendar for consistent week grouping
            df['period'] = df[datetime_column].dt.isocalendar().year.astype(str) + '-W' + \
                          df[datetime_column].dt.isocalendar().week.astype(str).str.zfill(2)

        # Aggregate
        grouped = df.groupby(['reach_id', 'species', 'period']).agg({
            value_column: agg_func,
            datetime_column: 'first'  # Keep first datetime as representative
        }).reset_index()

        # Convert period back to representative datetime
        if self.temporal_resolution == "daily":
            grouped[datetime_column] = pd.to_datetime(grouped['period'].astype(str))
        elif self.temporal_resolution == "weekly":
            # Parse year-week format
            grouped[datetime_column] = pd.to_datetime(
                grouped['period'].astype(str) + '-1',
                format='%Y-W%W-%w'
            )
        elif self.temporal_resolution == "monthly":
            grouped[datetime_column] = grouped['period'].dt.to_timestamp()
        elif self.temporal_resolution == "yearly":
            grouped[datetime_column] = grouped['period'].dt.to_timestamp()

        # Drop the temporary period column
        grouped = grouped.drop(columns=['period'])

        logger.debug(f"Aggregated {len(df)} records to {len(grouped)} records "
                     f"at {self.temporal_resolution} resolution")

        return grouped

    def _aggregate_matched_data(self, matched: pd.DataFrame) -> pd.DataFrame:
        """
        Aggregate matched observation-simulation pairs to the specified temporal resolution.

        Parameters
        ----------
        matched : pd.DataFrame
            DataFrame with datetime, reach_id, species, observed, simulated columns

        Returns
        -------
        pd.DataFrame
            Aggregated DataFrame
        """
        if self.temporal_resolution == "native" or matched.empty:
            return matched

        matched = matched.copy()
        matched['datetime'] = pd.to_datetime(matched['datetime'])

        # Get the aggregation frequency
        freq = self._get_temporal_grouper(self.temporal_resolution)

        # Define aggregation function
        agg_func_map = {
            "mean": "mean",
            "sum": "sum",
            "median": "median",
            "min": "min",
            "max": "max"
        }
        agg_func = agg_func_map.get(self.aggregation_method, "mean")

        # Create period column for grouping
        if self.temporal_resolution == "daily":
            matched['period'] = matched['datetime'].dt.date
        elif self.temporal_resolution == "weekly":
            matched['period'] = matched['datetime'].dt.isocalendar().year.astype(str) + '-W' + \
                               matched['datetime'].dt.isocalendar().week.astype(str).str.zfill(2)
        elif self.temporal_resolution == "monthly":
            matched['period'] = matched['datetime'].dt.to_period('M')
        elif self.temporal_resolution == "yearly":
            matched['period'] = matched['datetime'].dt.to_period('Y')

        # Aggregate both observed and simulated values
        aggregated = matched.groupby(['reach_id', 'species', 'period']).agg({
            'observed': agg_func,
            'simulated': agg_func,
            'datetime': 'first'
        }).reset_index()

        # Drop the period column
        aggregated = aggregated.drop(columns=['period'])

        logger.info(f"Temporal aggregation: {len(matched)} pairs → {len(aggregated)} pairs "
                    f"({self.temporal_resolution}, {self.aggregation_method})")

        return aggregated

    def _match_obs_sim(self, simulated: pd.DataFrame) -> pd.DataFrame:
        """
        Match observations with simulated values.
        """
        if simulated.empty:
            return pd.DataFrame()

        matched_data = []

        for _, obs_row in self.observations.iterrows():
            obs_dt = obs_row['datetime']
            obs_reach = obs_row['reach_id']
            obs_species = obs_row['species']
            obs_val = obs_row['value']

            # Find matching simulation
            mask = (
                (simulated['reach_id'] == obs_reach) &
                (simulated['species'] == obs_species)
            )
            sim_subset = simulated[mask]

            if sim_subset.empty:
                continue

            # Find closest time
            time_diffs = abs(sim_subset['datetime'] - obs_dt)
            closest_idx = time_diffs.idxmin()
            sim_val = sim_subset.loc[closest_idx, 'simulated']

            # Only match if within reasonable time window (e.g., 1 day)
            if time_diffs.loc[closest_idx] <= pd.Timedelta(days=1):
                matched_data.append({
                    'datetime': obs_dt,
                    'reach_id': obs_reach,
                    'species': obs_species,
                    'observed': obs_val,
                    'simulated': sim_val
                })

        return pd.DataFrame(matched_data)

    # =========================================================================
    # Objective Function Metrics
    # =========================================================================

    @staticmethod
    def rmse(obs: np.ndarray, sim: np.ndarray) -> float:
        """Root Mean Square Error."""
        return np.sqrt(np.nanmean((obs - sim) ** 2))

    @staticmethod
    def nse(obs: np.ndarray, sim: np.ndarray) -> float:
        """Nash-Sutcliffe Efficiency."""
        numerator = np.nansum((obs - sim) ** 2)
        denominator = np.nansum((obs - np.nanmean(obs)) ** 2)
        if denominator == 0:
            return -np.inf
        return 1.0 - numerator / denominator

    @staticmethod
    def nse_minimization(obs: np.ndarray, sim: np.ndarray) -> float:
        """Nash-Sutcliffe Efficiency transformed for minimization (1 - NSE)."""
        return 1.0 - ObjectiveFunction.nse(obs, sim)

    @staticmethod
    def kge(obs: np.ndarray, sim: np.ndarray) -> float:
        """Kling-Gupta Efficiency."""
        # Correlation coefficient
        if np.nanstd(obs) == 0 or np.nanstd(sim) == 0:
            r = 0
        else:
            r = np.corrcoef(obs[~np.isnan(obs) & ~np.isnan(sim)],
                           sim[~np.isnan(obs) & ~np.isnan(sim)])[0, 1]

        # Ratio of standard deviations
        alpha = np.nanstd(sim) / np.nanstd(obs) if np.nanstd(obs) != 0 else 0

        # Ratio of means
        beta = np.nanmean(sim) / np.nanmean(obs) if np.nanmean(obs) != 0 else 0

        # KGE
        kge_val = 1.0 - np.sqrt((r - 1) ** 2 + (alpha - 1) ** 2 + (beta - 1) ** 2)
        return kge_val

    @staticmethod
    def kge_minimization(obs: np.ndarray, sim: np.ndarray) -> float:
        """KGE transformed for minimization (1 - KGE)."""
        return 1.0 - ObjectiveFunction.kge(obs, sim)

    @staticmethod
    def pbias(obs: np.ndarray, sim: np.ndarray) -> float:
        """Percent Bias."""
        return 100 * np.nansum(sim - obs) / np.nansum(obs)

    def get_matched_data(self) -> pd.DataFrame:
        """
        Get the last matched observation-simulation data.

        Returns
        -------
        pd.DataFrame
            DataFrame with datetime, reach_id, species, observed, simulated columns
            (already aggregated to the specified temporal resolution)
        """
        if hasattr(self, '_last_matched_data'):
            return self._last_matched_data.copy()
        return pd.DataFrame()

    def get_temporal_resolution_info(self) -> Dict:
        """
        Get information about the temporal resolution settings.

        Returns
        -------
        Dict
            Dictionary with resolution and aggregation method info
        """
        return {
            "temporal_resolution": self.temporal_resolution,
            "aggregation_method": self.aggregation_method,
            "description": self._get_resolution_description()
        }

    def _get_resolution_description(self) -> str:
        """Get human-readable description of temporal resolution."""
        descriptions = {
            "native": "Native model timestep (no aggregation)",
            "daily": "Daily aggregated values",
            "weekly": "Weekly aggregated values",
            "monthly": "Monthly aggregated values",
            "yearly": "Yearly aggregated values"
        }
        return descriptions.get(self.temporal_resolution, "Unknown")
