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
                 aggregation_method: str = "mean",
                 hostmodel: str = "mizuroute",
                 h5_mapping_key: Optional[str] = None,
                 use_primary_only: bool = True,
                 zone_select=None,
                 metric_focus: str = "both"):
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
        # Host-aware HDF5 mapping_key: SUMMA writes `hruId`, mizuRoute
        # writes `reachID`.  Explicit `h5_mapping_key` wins; otherwise
        # we derive from `hostmodel` to match Gen_Report's convention.
        self.hostmodel = (hostmodel or "mizuroute").lower()
        if h5_mapping_key:
            self.h5_mapping_key = h5_mapping_key
        elif self.hostmodel == "summa":
            self.h5_mapping_key = "hruId"
        else:
            self.h5_mapping_key = "reachID"
        # When True, only observations from PRIMARY stations
        # (pouring-point obs in SUMMA; every matched station in
        # mizuRoute) contribute to the objective.  Secondary (gray)
        # SUMMA stations still appear on the calibration results-report
        # map but they're excluded from the metric so the calibration
        # doesn't average over multiple obs per basin.
        self.use_primary_only = bool(use_primary_only)
        # SUMMA writes one column per vertical zone / soil layer per HRU
        # (``hruId_1_z1``, ``_z2`` ...).  By default every zone of an HRU is
        # averaged into the HRU value.  Set ``zone_select`` to a layer index
        # (e.g. 1) to use ONLY that layer's output in the objective AND the
        # report; ``None``/``'mean'``/``'all'`` keeps the average.  No effect
        # on hosts without zoned columns (mizuRoute).
        self.zone_select = self._parse_zone(zone_select)
        if self.zone_select is not None:
            logger.info(f"Zone selection: using only soil layer / zone "
                        f"z{self.zone_select} for the objective.")
        # Objective focus — which error components matter most.  Works with ANY
        # base metric:
        #   "both"      → use the chosen metric as-is (KGE/NSE/RMSE/PBIAS)
        #   "phase"     → reward TIMING/SHAPE (correlation) — avoids the
        #                 "off-phase but right-magnitude wins" pathology
        #   "magnitude" → reward level + variability (bias + spread), not timing
        # "phase"/"magnitude" use the universal correlation vs (variability+bias)
        # decomposition (Gupta 2009 scaled KGE), because RMSE/PBIAS/NSE don't
        # separate timing from level/spread on their own.
        self.metric_focus = str(metric_focus or "both").strip().lower()
        if self.metric_focus not in self._FOCUS_WEIGHTS:
            self.metric_focus = "both"
        if self.metric_focus != "both":
            logger.info(f"Objective focus: '{self.metric_focus}' — using the "
                        f"correlation/variability/bias decomposition with "
                        f"weights {self._FOCUS_WEIGHTS[self.metric_focus]} "
                        f"(overrides the '{self.metric}' base metric).")

        # Validate temporal resolution
        valid_resolutions = ["native", "daily", "weekly", "monthly", "yearly"]
        if temporal_resolution not in valid_resolutions:
            raise ValueError(
                f"Invalid temporal_resolution '{temporal_resolution}'. "
                f"Valid options: {valid_resolutions}"
            )

        # Validate aggregation method.  Note: "sum" is intentionally
        # excluded — summing concentrations across a time window
        # produces a value with no physical meaning for water quality
        # metrics (you'd be adding mg/L at different timestamps).  Use
        # "mean" (default) or "median" for concentrations, or "min" /
        # "max" for extreme-value summaries.
        valid_methods = ["mean", "median", "min", "max"]
        if aggregation_method not in valid_methods:
            raise ValueError(
                f"Invalid aggregation_method '{aggregation_method}'. "
                f"Valid options: {valid_methods}"
            )

        logger.info(f"Objective function configured with temporal_resolution='{temporal_resolution}', "
                    f"aggregation_method='{aggregation_method}'")

        # Setup H5 reader import path.  `from hdf5_support_lib import
        # Read_h5_driver` needs the PARENT of the hdf5_support_lib package dir
        # on sys.path — so insert both the given path AND its parent, making
        # the import resolve whether h5_reader_path points at the package dir
        # (…/2_Read_Outputs/hdf5_support_lib) or at its parent (…/2_Read_Outputs).
        if h5_reader_path:
            _hp = Path(h5_reader_path)
            sys.path.insert(0, str(_hp))
            sys.path.insert(0, str(_hp.parent))

        self.observations = self._load_observations()

    def _load_observations(self) -> pd.DataFrame:
        """Load and preprocess observation data.

        Honours the ``is_primary`` column produced by
        ``prepare_calibration_observations.py``: when
        ``use_primary_only=True`` (default) only the pouring-point obs of
        each basin / matched reach contributes to the metric.  Secondary
        (gray) station rows are still kept in memory if present, exposed
        via ``self.observations_secondary`` for the results-report map.
        """
        try:
            obs = pd.read_csv(self.observation_path, parse_dates=['datetime'])
        except FileNotFoundError:
            logger.warning(f"Observation file not found: {self.observation_path}")
            self.observations_secondary = pd.DataFrame()
            return pd.DataFrame()

        # Filter by species
        obs = obs[obs['species'].isin(self.target_species)]

        # Filter by reaches
        if self.target_reaches != "all" and isinstance(self.target_reaches, list):
            obs = obs[obs['reach_id'].isin(self.target_reaches)]

        # Separate primary vs secondary so the results-report map can
        # still draw the gray ones even when they're excluded from the
        # objective.  Rows missing the column are treated as primary
        # (legacy CSVs without spatial-matching info).
        if 'is_primary' in obs.columns:
            _is_prim = obs['is_primary'].fillna(True).astype(bool)
            self.observations_secondary = obs[~_is_prim].copy()
            if self.use_primary_only:
                obs = obs[_is_prim].copy()
                logger.info(
                    f"use_primary_only=True → "
                    f"{len(obs)} primary obs kept for metric, "
                    f"{len(self.observations_secondary)} secondary "
                    f"obs reserved for map display only."
                )
        else:
            self.observations_secondary = pd.DataFrame()

        logger.info(f"Loaded {len(obs)} observations for "
                    f"{obs['species'].nunique()} species")
        return obs

    def check_observation_availability(self,
                                       sim_start=None,
                                       sim_end=None) -> Dict:
        """Pre-flight: are there observations the metric can actually score?

        A calibration only yields a usable objective when at least two
        observations fall within BOTH:

          * the target reaches/HRUs and species — already enforced upstream by
            :meth:`_load_observations` (so ``self.observations`` is the spatially
            and species-filtered set), AND
          * the simulated time window ``[sim_start, sim_end]``.

        Supplying ``sim_start`` / ``sim_end`` catches the classic temporal-
        mismatch trap: the gauge has data (e.g. 1999-2005) but the model only
        simulates 1990-1997, so there are no observed-vs-simulated pairs and
        EVERY evaluation silently returns the failure penalty (1e10).  This
        check lets the driver abort up front with a clear explanation instead.

        Parameters
        ----------
        sim_start, sim_end : str | datetime | None
            The simulated period.  Strings are parsed with ``pd.to_datetime``.
            When either is ``None`` the temporal filter is skipped (only
            spatial/species availability is checked).

        Returns
        -------
        Dict
            ``ok`` (bool), ``headline`` (one-line reason), ``message`` (multi-
            line, plain-text explanation suitable for the terminal/log),
            ``n_total``, ``n_in_window``, ``per_reach`` and the resolved
            ``sim_start`` / ``sim_end``.
        """
        obs = self.observations

        ws = pd.to_datetime(sim_start) if sim_start is not None else None
        we = pd.to_datetime(sim_end) if sim_end is not None else None
        win_txt = (f"{ws:%Y-%m-%d} -> {we:%Y-%m-%d}"
                   if (ws is not None and we is not None) else "unknown")

        if self.target_reaches == "all" or not isinstance(self.target_reaches,
                                                           list):
            reach_txt = "all reaches with observations"
        else:
            reach_txt = ", ".join(str(r) for r in self.target_reaches)
        species_txt = ", ".join(self.target_species)

        # Case 0 — nothing matches the target species/reaches at all.
        if obs is None or obs.empty:
            msg = (
                "No observations match the calibration targets.\n"
                f"  Target species  : {species_txt}\n"
                f"  Target reaches  : {reach_txt}\n"
                f"  Simulated window: {win_txt}\n\n"
                "No performance metric (KGE/NSE/RMSE) can be computed. Check "
                "that the observation CSV contains these species and that the "
                "target reach/HRU IDs match the observation 'reach_id' column "
                "(spatial matching)."
            )
            return {
                "ok": False,
                "headline": "no observations match the target species/reaches",
                "message": msg,
                "n_total": 0, "n_in_window": 0, "per_reach": {},
                "sim_start": str(sim_start) if sim_start else None,
                "sim_end": str(sim_end) if sim_end else None,
            }

        # Per-reach breakdown (total vs. inside the simulated window).
        dts = pd.to_datetime(obs['datetime'], errors='coerce')
        if ws is not None and we is not None:
            in_win_mask = (dts >= ws) & (dts <= we)
        else:
            in_win_mask = pd.Series(True, index=obs.index)

        per_reach = {}
        n_in_window = 0
        for rid, grp in obs.groupby('reach_id'):
            gmask = in_win_mask.loc[grp.index]
            n_win = int(gmask.sum())
            gdt = pd.to_datetime(grp['datetime'], errors='coerce').dropna()
            per_reach[str(rid)] = {
                "total": int(len(grp)),
                "in_window": n_win,
                "obs_start": gdt.min().strftime('%Y-%m-%d') if len(gdt) else "-",
                "obs_end": gdt.max().strftime('%Y-%m-%d') if len(gdt) else "-",
            }
            n_in_window += n_win

        n_total = int(len(obs))
        # KGE/NSE/RMSE need >=2 paired points to be meaningful.
        ok = n_in_window >= 2

        if ok:
            headline = f"{n_in_window} observation(s) inside the simulated window"
            msg = (
                f"OK: {n_in_window} of {n_total} target observation(s) fall "
                f"inside the simulated window ({win_txt}).\n"
                f"  Target species : {species_txt}\n"
                f"  Target reaches : {reach_txt}"
            )
            return {"ok": True, "headline": headline, "message": msg,
                    "n_total": n_total, "n_in_window": n_in_window,
                    "per_reach": per_reach,
                    "sim_start": (ws.strftime('%Y-%m-%d %H:%M')
                                  if ws is not None else None),
                    "sim_end": (we.strftime('%Y-%m-%d %H:%M')
                                if we is not None else None)}

        # Not OK — obs exist for the targets but (almost) none overlap the run.
        lines = [
            f"    reach {rid}: {info['in_window']} in-window / "
            f"{info['total']} total (obs span {info['obs_start']} ... "
            f"{info['obs_end']})"
            for rid, info in sorted(per_reach.items())
        ]
        if n_in_window == 0:
            why = ("NONE of them fall inside the simulated time window, so "
                   "there are no observed-vs-simulated pairs")
        else:
            why = (f"only {n_in_window} observation(s) fall inside the "
                   "simulated window (at least 2 are needed to compute a "
                   "metric)")
        headline = "observations exist but do not overlap the simulated period"
        msg = (
            "The target reach(es)/species have observations, but "
            f"{why}. No metric (KGE/NSE/RMSE) can be computed and every "
            "evaluation would return the failure penalty.\n\n"
            f"  Target species  : {species_txt}\n"
            f"  Target reaches  : {reach_txt}\n"
            f"  Simulated window: {win_txt}\n"
            f"  Obs on target reaches (total): {n_total}\n"
            f"  Obs inside the window        : {n_in_window}\n\n"
            "  Per-reach availability:\n" + "\n".join(lines) + "\n\n"
            "Fix one of the following, then re-run:\n"
            "  1. Extend the simulation period (control-file <sim_start> / "
            "<sim_end>, or the calibration window) to cover the observation "
            "dates above - make sure the model forcing also covers it; or\n"
            "  2. Target a reach/HRU whose observations fall inside the "
            "simulated window."
        )
        return {"ok": False, "headline": headline, "message": msg,
                "n_total": n_total, "n_in_window": n_in_window,
                "per_reach": per_reach,
                "sim_start": (ws.strftime('%Y-%m-%d %H:%M')
                              if ws is not None else None),
                "sim_end": (we.strftime('%Y-%m-%d %H:%M')
                            if we is not None else None)}

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

        # Keep the full simulated series (all target species at the target
        # reaches, including species that have NO observations) so the results
        # report can still plot a simulated-only time series for them.
        self._last_simulated = simulated

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

            if self.metric_focus in ("phase", "magnitude"):
                # PHASE / MAGNITUDE focus is defined via the universal KGE
                # decomposition (correlation r vs variability+bias) and so
                # applies to ANY base metric: RMSE conflates timing+level,
                # PBIAS is pure bias, and NSE doesn't split into a weightable
                # sum — only the r/alpha/beta decomposition cleanly isolates
                # "phase" from "magnitude".  ('both' keeps the chosen metric.)
                _sr, _sa, _sb = self._FOCUS_WEIGHTS[self.metric_focus]
                obj = 1.0 - self.kge_scaled(obs_vals, sim_vals, _sr, _sa, _sb)
            elif self.metric == "RMSE":
                obj = self.rmse(obs_vals, sim_vals)
            elif self.metric == "NSE":
                obj = self.nse_minimization(obs_vals, sim_vals)
            elif self.metric == "KGE":
                obj = self.kge_minimization(obs_vals, sim_vals)
            elif self.metric == "PBIAS":
                # Perfect bias = 0; minimise its absolute value.
                obj = abs(self.pbias(obs_vals, sim_vals))
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

        # Track the BEST evaluation's matched + simulated series so the report
        # shows the actual best fit — not merely the most recent evaluation
        # (which is what _last_* hold).  Lower objective = better.
        if (not hasattr(self, "_best_obj") or self._best_obj is None
                or weighted_obj < self._best_obj):
            self._best_obj = weighted_obj
            self._best_matched_data = self._last_matched_data
            self._best_simulated = getattr(self, "_last_simulated", None)
        return weighted_obj

    @staticmethod
    def _parse_zone(z):
        """Normalise a user zone selection to an int layer index or None.

        Accepts ``None``/``''``/``'all'``/``'mean'``/``'average'`` (→ None,
        meaning average all zones) or anything containing a number such as
        ``2``, ``'z2'``, ``'layer 2'`` (→ ``2``)."""
        if z is None:
            return None
        s = str(z).strip().lower()
        if s in ("", "none", "all", "mean", "average", "avg"):
            return None
        import re
        m = re.search(r"(\d+)", s)
        return int(m.group(1)) if m else None

    @staticmethod
    def _zone_of(raw):
        """Return the vertical-zone / soil-layer index of an openWQ column
        (``hruId_1_z2`` → ``2``), or ``None`` when the column carries no zone
        tag (e.g. mizuRoute ``reachID_123``)."""
        import re
        s = raw.decode() if isinstance(raw, (bytes, bytearray)) else str(raw)
        m = re.search(r"_z(\d+)$", s.strip(), re.IGNORECASE)
        return int(m.group(1)) if m else None

    @staticmethod
    def _normalize_feature_id(raw):
        """Normalise an openWQ feature/column ID so it matches the integer
        observation ``reach_id``.

        Handles:
          * ``bytes`` -> ``str``
          * ``hruId_`` / ``reachID_`` prefixes (library column names)
          * SUMMA vertical-zone / soil-layer suffixes, e.g. ``1_z1`` -> ``1``
            (so every layer of HRU 1 maps to HRU 1, matching the obs HRU id)

        Returns an ``int`` when possible (obs ``reach_id`` is integer),
        otherwise the cleaned string (so non-numeric ids still flow through
        rather than being silently dropped).
        """
        import re
        s = raw.decode() if isinstance(raw, (bytes, bytearray)) else str(raw)
        s = s.strip()
        for pre in ("hruId_", "reachID_", "hruid_", "reachid_"):
            if s.startswith(pre):
                s = s[len(pre):]
                break
        # Strip a trailing vertical zone / layer tag like '_z1', '_Z12'.
        s = re.sub(r"_z\d+$", "", s, flags=re.IGNORECASE)
        try:
            return int(s)
        except ValueError:
            try:
                return int(float(s))
            except ValueError:
                return s

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
            "mapping_key": self.h5_mapping_key,
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

            # Only the target reaches/HRUs are ever used (the rest are dropped
            # by _match_obs_sim).  Keep just those so we don't melt every reach
            # × every timestep into a multi-million-row frame per eval — slow
            # and a real OOM risk under n_parallel.  None => keep all.
            _target_set = (set(self.target_reaches)
                           if isinstance(self.target_reaches, list) else None)

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
                        # df has time index and columns like
                        # "reachID_123" (mizuRoute) or "hruId_123" (SUMMA).
                        for col in df.columns:
                            # SUMMA writes zoned HRU ids like 'hruId_1_z1';
                            # normalise to the base feature id (int 1) so it
                            # matches the integer observation 'reach_id'.
                            # When a specific zone is selected, keep only that
                            # layer (columns with no zone tag always pass).
                            if self.zone_select is not None:
                                _zc = self._zone_of(col)
                                if _zc is not None and _zc != self.zone_select:
                                    continue
                            reach_id = self._normalize_feature_id(col)
                            if _target_set is not None and reach_id not in _target_set:
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
        # Keep only target reaches/HRUs (see _extract_simulated) so the fallback
        # also avoids a multi-million-row frame per eval.  None => keep all.
        _target_set = (set(self.target_reaches)
                       if isinstance(self.target_reaches, list) else None)

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
                        # Get IDs — host-aware (mizuRoute writes
                        # /reachID, SUMMA writes /hruId).  We probe the
                        # configured key first, then fall back to the
                        # other so a misconfiguration still produces
                        # something usable.
                        _id_key = self.h5_mapping_key
                        _id_root = f"/{_id_key}"
                        if _id_root in hf:
                            reach_ids = [x.decode() if isinstance(x, bytes) else x
                                         for x in hf[_id_root][:]]
                        elif '/reachID' in hf:
                            _id_key = 'reachID'
                            reach_ids = [x.decode() if isinstance(x, bytes) else x
                                         for x in hf['/reachID'][:]]
                        elif '/hruId' in hf:
                            _id_key = 'hruId'
                            reach_ids = [x.decode() if isinstance(x, bytes) else x
                                         for x in hf['/hruId'][:]]
                        else:
                            reach_ids = list(range(hf[list(hf.keys())[0]].shape[1]))

                        # NEW extensible layout: a single /concentrations
                        # [time x cells] dataset + /timestamps (current openWQ
                        # writer).  Mirror hdf5_support_lib/Read_h5_driver so
                        # this fallback stays correct even when that library
                        # can't be imported.  (The block below handles the
                        # LEGACY one-dataset-per-timestep files.)
                        if 'concentrations' in hf and 'timestamps' in hf:
                            ts_raw = hf['timestamps'][:]
                            conc = np.asarray(hf['concentrations'][:],
                                              dtype=np.float64)  # (ntime,ncells)
                            conc[conc == self.no_data_flag] = np.nan
                            norm_ids = [self._normalize_feature_id(r)
                                        for r in reach_ids]
                            zone_ids = [self._zone_of(r) for r in reach_ids]
                            for ti in range(conc.shape[0]):
                                ts = ts_raw[ti]
                                ts = (ts.decode()
                                      if isinstance(ts, (bytes, bytearray))
                                      else str(ts))
                                try:
                                    dt = pd.to_datetime(
                                        ts, format='%Y%b%d-%H:%M:%S')
                                except Exception:
                                    try:
                                        dt = pd.to_datetime(ts)
                                    except Exception:
                                        continue
                                row = conc[ti, :]
                                for i, rid in enumerate(norm_ids):
                                    if (_target_set is not None
                                            and rid not in _target_set):
                                        continue
                                    if (self.zone_select is not None
                                            and zone_ids[i] is not None
                                            and zone_ids[i] != self.zone_select):
                                        continue
                                    val = (row[i] if i < row.shape[0]
                                           else np.nan)
                                    all_data.append({
                                        'datetime': dt, 'reach_id': rid,
                                        'species': species,
                                        'simulated': val})
                            continue  # file done; skip legacy per-timestep block

                        # Get timestamps (skip the well-known metadata keys)
                        _metadata_keys = {
                            'xyz_elements', 'reachID', 'hruId',
                            'mapping_key', _id_key,
                        }
                        timestamps = [k for k in hf.keys()
                                      if k not in _metadata_keys]

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
                                # Normalise SUMMA zoned HRU ids ('1_z1' -> 1)
                                # / strip 'reachID_' etc. so they match the
                                # integer observation 'reach_id'.
                                rid = self._normalize_feature_id(reach_id)
                                if (_target_set is not None
                                        and rid not in _target_set):
                                    continue
                                if self.zone_select is not None:
                                    _zc = self._zone_of(reach_id)
                                    if _zc is not None and _zc != self.zone_select:
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

    # Per-component KGE scaling for the objective "focus" (Gupta et al. 2009):
    # (s_r, s_alpha, s_beta) on (correlation, variability, bias).
    _FOCUS_WEIGHTS = {
        "both":      (1.0, 1.0, 1.0),    # standard KGE
        "phase":     (1.0, 0.25, 0.25),  # reward timing/shape (correlation)
        "magnitude": (0.25, 1.0, 1.0),   # reward level + spread, not timing
    }

    @staticmethod
    def kge_scaled(obs: np.ndarray, sim: np.ndarray,
                   sr: float = 1.0, sa: float = 1.0, sb: float = 1.0) -> float:
        """KGE with per-component scaling factors so the objective can favour
        PHASE (correlation r) or MAGNITUDE (variability alpha + bias beta).
        With sr=sa=sb=1 this is the standard KGE."""
        mask = ~(np.isnan(obs) | np.isnan(sim))
        o, s = obs[mask], sim[mask]
        if len(o) < 2:
            return -np.inf
        if np.nanstd(o) == 0 or np.nanstd(s) == 0:
            r = 0.0
        else:
            r = np.corrcoef(o, s)[0, 1]
        alpha = np.nanstd(s) / np.nanstd(o) if np.nanstd(o) != 0 else 0.0
        beta = np.nanmean(s) / np.nanmean(o) if np.nanmean(o) != 0 else 0.0
        return 1.0 - np.sqrt((sr * (r - 1)) ** 2
                             + (sa * (alpha - 1)) ** 2
                             + (sb * (beta - 1)) ** 2)

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
        # Prefer the BEST evaluation's matched data (so the report's "best fit"
        # is the actual best), falling back to the most recent one.
        _b = getattr(self, '_best_matched_data', None)
        if _b is not None:
            return _b.copy()
        if hasattr(self, '_last_matched_data'):
            return self._last_matched_data.copy()
        return pd.DataFrame()

    def get_simulated_data(self) -> pd.DataFrame:
        """Return the full simulated series (all target species at the target
        reaches, incl. species with NO observations) for the BEST evaluation —
        falling back to the most recent.  Lets the report draw the best-fit
        simulated-only series for objective species lacking observations."""
        _b = getattr(self, '_best_simulated', None)
        if _b is not None:
            return _b.copy()
        if hasattr(self, '_last_simulated') and self._last_simulated is not None:
            return self._last_simulated.copy()
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
