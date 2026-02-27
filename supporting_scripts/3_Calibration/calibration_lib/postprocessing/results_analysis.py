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
Results Analysis Module
=======================

Analyzes calibration results, generates convergence plots, parameter evolution,
and performance summaries.
"""

import json
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any
from dataclasses import dataclass
import logging

logger = logging.getLogger(__name__)

# Try to import matplotlib (optional for plotting)
try:
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    logger.warning("matplotlib not available - plotting disabled")

# Try to import fiona for spatial data reading (map generation)
try:
    import fiona
    HAS_FIONA = True
except ImportError:
    HAS_FIONA = False
    logger.info("fiona not available - map generation disabled")


@dataclass
class CalibrationSummary:
    """Summary of calibration results."""
    best_parameters: Dict[str, float]
    best_objective: float
    n_evaluations: int
    convergence_eval: int  # Evaluation where best was found
    runtime_hours: float
    parameter_bounds: Dict[str, Tuple[float, float]]
    objective_function: str

    def to_dict(self) -> Dict:
        """Convert to dictionary."""
        return {
            "best_parameters": self.best_parameters,
            "best_objective": self.best_objective,
            "n_evaluations": self.n_evaluations,
            "convergence_eval": self.convergence_eval,
            "runtime_hours": self.runtime_hours,
            "parameter_bounds": self.parameter_bounds,
            "objective_function": self.objective_function,
        }

    def save_json(self, filepath: Path):
        """Save summary to JSON file."""
        with open(filepath, 'w') as f:
            json.dump(self.to_dict(), f, indent=2)


class ResultsAnalyzer:
    """
    Analyzes and visualizes calibration results.

    Generates plots, statistics, and reports from calibration history.
    """

    def __init__(self,
                 results_dir: Path,
                 parameter_definitions: List[Dict] = None):
        """
        Initialize results analyzer.

        Parameters
        ----------
        results_dir : Path
            Directory containing calibration results
        parameter_definitions : List[Dict], optional
            Parameter definitions from calibration config
        """
        self.results_dir = Path(results_dir)
        self.parameter_definitions = parameter_definitions or []

        # Load results if available
        self.history: List[Dict] = []
        self.sensitivity_results: Optional[Dict] = None

        self._load_results()

    def _load_results(self):
        """Load calibration results from files."""
        # Load calibration history
        history_file = self.results_dir / "calibration_history.json"
        if history_file.exists():
            with open(history_file, 'r') as f:
                self.history = json.load(f)
            logger.info(f"Loaded {len(self.history)} evaluations from history")

        # Load sensitivity results if available
        sens_file = self.results_dir / "sensitivity_results.json"
        if sens_file.exists():
            with open(sens_file, 'r') as f:
                self.sensitivity_results = json.load(f)
            logger.info("Loaded sensitivity analysis results")

    def load_from_checkpoint(self, checkpoint_file: Path):
        """Load results from a checkpoint file."""
        with open(checkpoint_file, 'r') as f:
            checkpoint = json.load(f)

        self.history = checkpoint.get("history", [])
        logger.info(f"Loaded {len(self.history)} evaluations from checkpoint")

    def get_summary(self) -> CalibrationSummary:
        """
        Generate calibration summary.

        Returns
        -------
        CalibrationSummary
            Summary of calibration results
        """
        if not self.history:
            raise ValueError("No calibration history available")

        # Find best evaluation
        best_idx = 0
        best_obj = float('inf')

        for i, entry in enumerate(self.history):
            if entry.get("objective", float('inf')) < best_obj:
                best_obj = entry["objective"]
                best_idx = i

        best_params = self.history[best_idx].get("parameters", {})

        # Get parameter bounds
        bounds = {}
        for param_def in self.parameter_definitions:
            name = param_def["name"]
            bounds[name] = tuple(param_def.get("bounds", (0, 1)))

        # Calculate runtime
        if self.history:
            start_time = self.history[0].get("timestamp", 0)
            end_time = self.history[-1].get("timestamp", start_time)
            runtime_hours = (end_time - start_time) / 3600
        else:
            runtime_hours = 0

        return CalibrationSummary(
            best_parameters=best_params,
            best_objective=best_obj,
            n_evaluations=len(self.history),
            convergence_eval=best_idx,
            runtime_hours=runtime_hours,
            parameter_bounds=bounds,
            objective_function="KGE"  # Default
        )

    def get_convergence_data(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Get convergence data (evaluation number vs best objective).

        Returns
        -------
        Tuple[np.ndarray, np.ndarray]
            (evaluation_numbers, best_objectives)
        """
        if not self.history:
            return np.array([]), np.array([])

        evals = []
        best_objs = []
        best_so_far = float('inf')

        for i, entry in enumerate(self.history):
            obj = entry.get("objective", float('inf'))
            if obj < best_so_far:
                best_so_far = obj
            evals.append(i)
            best_objs.append(best_so_far)

        return np.array(evals), np.array(best_objs)

    def get_parameter_evolution(self,
                                param_name: str
                                ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Get parameter evolution over evaluations.

        Parameters
        ----------
        param_name : str
            Parameter name

        Returns
        -------
        Tuple[np.ndarray, np.ndarray]
            (evaluation_numbers, parameter_values)
        """
        evals = []
        values = []

        for i, entry in enumerate(self.history):
            params = entry.get("parameters", {})
            if param_name in params:
                evals.append(i)
                values.append(params[param_name])

        return np.array(evals), np.array(values)

    def get_parameter_statistics(self) -> pd.DataFrame:
        """
        Get statistics for all parameters across evaluations.

        Returns
        -------
        pd.DataFrame
            Parameter statistics
        """
        if not self.history:
            return pd.DataFrame()

        # Collect all parameters
        all_params = {}
        for entry in self.history:
            params = entry.get("parameters", {})
            for name, value in params.items():
                if name not in all_params:
                    all_params[name] = []
                all_params[name].append(value)

        # Compute statistics
        stats = []
        for name, values in all_params.items():
            values = np.array(values)
            stats.append({
                "parameter": name,
                "min": np.min(values),
                "max": np.max(values),
                "mean": np.mean(values),
                "std": np.std(values),
                "median": np.median(values),
            })

        return pd.DataFrame(stats)

    def get_objective_statistics(self) -> Dict[str, float]:
        """
        Get objective function statistics.

        Returns
        -------
        Dict[str, float]
            Objective statistics
        """
        if not self.history:
            return {}

        objectives = [entry.get("objective", float('inf'))
                     for entry in self.history
                     if entry.get("objective") is not None]

        if not objectives:
            return {}

        return {
            "best": min(objectives),
            "worst": max(objectives),
            "mean": np.mean(objectives),
            "std": np.std(objectives),
            "median": np.median(objectives),
            "q25": np.percentile(objectives, 25),
            "q75": np.percentile(objectives, 75),
        }

    def plot_convergence(self,
                         output_path: Optional[Path] = None,
                         show: bool = False) -> Optional[Any]:
        """
        Plot convergence curve.

        Parameters
        ----------
        output_path : Path, optional
            Path to save plot
        show : bool
            Whether to display plot

        Returns
        -------
        Figure or None
            Matplotlib figure if available
        """
        if not HAS_MATPLOTLIB:
            logger.warning("matplotlib not available - cannot create plot")
            return None

        evals, best_objs = self.get_convergence_data()
        if len(evals) == 0:
            logger.warning("No data for convergence plot")
            return None

        fig, ax = plt.subplots(figsize=(10, 6))

        # Plot all objectives
        all_objs = [entry.get("objective", np.nan) for entry in self.history]
        ax.scatter(range(len(all_objs)), all_objs, alpha=0.3, s=10,
                   label='All evaluations', color='gray')

        # Plot best-so-far
        ax.plot(evals, best_objs, 'b-', linewidth=2, label='Best objective')

        # Mark final best
        ax.scatter([evals[-1]], [best_objs[-1]], color='red', s=100,
                   zorder=5, marker='*', label=f'Best: {best_objs[-1]:.4f}')

        ax.set_xlabel('Evaluation Number', fontsize=12)
        ax.set_ylabel('Objective Function Value', fontsize=12)
        ax.set_title('Calibration Convergence', fontsize=14)
        ax.legend()
        ax.grid(True, alpha=0.3)

        plt.tight_layout()

        if output_path:
            plt.savefig(output_path, dpi=150, bbox_inches='tight')
            logger.info(f"Convergence plot saved to {output_path}")

        if show:
            plt.show()

        return fig

    def plot_parameter_evolution(self,
                                 param_names: List[str] = None,
                                 output_path: Optional[Path] = None,
                                 show: bool = False) -> Optional[Any]:
        """
        Plot parameter evolution over evaluations.

        Parameters
        ----------
        param_names : List[str], optional
            Parameters to plot (default: all)
        output_path : Path, optional
            Path to save plot
        show : bool
            Whether to display plot

        Returns
        -------
        Figure or None
            Matplotlib figure if available
        """
        if not HAS_MATPLOTLIB:
            logger.warning("matplotlib not available - cannot create plot")
            return None

        # Get parameter names if not specified
        if param_names is None:
            if self.history:
                param_names = list(self.history[0].get("parameters", {}).keys())
            else:
                return None

        n_params = len(param_names)
        if n_params == 0:
            return None

        # Create subplots
        n_cols = min(3, n_params)
        n_rows = int(np.ceil(n_params / n_cols))

        fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 4*n_rows))
        if n_params == 1:
            axes = np.array([[axes]])
        elif n_rows == 1:
            axes = axes.reshape(1, -1)
        elif n_cols == 1:
            axes = axes.reshape(-1, 1)

        for i, param_name in enumerate(param_names):
            row = i // n_cols
            col = i % n_cols
            ax = axes[row, col]

            evals, values = self.get_parameter_evolution(param_name)
            if len(evals) > 0:
                ax.scatter(evals, values, alpha=0.5, s=10)
                ax.set_xlabel('Evaluation')
                ax.set_ylabel('Value')
                ax.set_title(param_name)
                ax.grid(True, alpha=0.3)

                # Add bounds if available
                for param_def in self.parameter_definitions:
                    if param_def["name"] == param_name:
                        bounds = param_def.get("bounds", (None, None))
                        if bounds[0] is not None:
                            ax.axhline(bounds[0], color='r', linestyle='--',
                                       alpha=0.5, label='Bounds')
                        if bounds[1] is not None:
                            ax.axhline(bounds[1], color='r', linestyle='--',
                                       alpha=0.5)
                        break

        # Hide empty subplots
        for i in range(n_params, n_rows * n_cols):
            row = i // n_cols
            col = i % n_cols
            axes[row, col].set_visible(False)

        plt.tight_layout()

        if output_path:
            plt.savefig(output_path, dpi=150, bbox_inches='tight')
            logger.info(f"Parameter evolution plot saved to {output_path}")

        if show:
            plt.show()

        return fig

    def plot_sensitivity_results(self,
                                 output_path: Optional[Path] = None,
                                 show: bool = False) -> Optional[Any]:
        """
        Plot sensitivity analysis results.

        Parameters
        ----------
        output_path : Path, optional
            Path to save plot
        show : bool
            Whether to display plot

        Returns
        -------
        Figure or None
            Matplotlib figure if available
        """
        if not HAS_MATPLOTLIB:
            logger.warning("matplotlib not available - cannot create plot")
            return None

        if not self.sensitivity_results:
            logger.warning("No sensitivity results available")
            return None

        method = self.sensitivity_results.get("method", "morris")

        if method == "morris":
            return self._plot_morris_results(output_path, show)
        elif method == "sobol":
            return self._plot_sobol_results(output_path, show)
        else:
            logger.warning(f"Unknown sensitivity method: {method}")
            return None

    def _plot_morris_results(self,
                             output_path: Optional[Path] = None,
                             show: bool = False) -> Optional[Any]:
        """Plot Morris screening results."""
        results = self.sensitivity_results

        param_names = results.get("parameter_names", [])
        mu_star = np.array(results.get("mu_star", []))
        sigma = np.array(results.get("sigma", []))

        if len(param_names) == 0:
            return None

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

        # Bar chart of mu_star
        y_pos = np.arange(len(param_names))
        sorted_idx = np.argsort(mu_star)[::-1]

        ax1.barh(y_pos, mu_star[sorted_idx], color='steelblue')
        ax1.set_yticks(y_pos)
        ax1.set_yticklabels([param_names[i] for i in sorted_idx])
        ax1.set_xlabel(r'$\mu^*$ (Mean absolute elementary effect)')
        ax1.set_title('Morris Sensitivity: Parameter Importance')
        ax1.grid(True, alpha=0.3, axis='x')

        # Scatter plot of mu_star vs sigma
        ax2.scatter(mu_star, sigma, s=100, alpha=0.7)
        for i, name in enumerate(param_names):
            ax2.annotate(name, (mu_star[i], sigma[i]),
                        xytext=(5, 5), textcoords='offset points',
                        fontsize=9)

        ax2.set_xlabel(r'$\mu^*$ (Mean absolute elementary effect)')
        ax2.set_ylabel(r'$\sigma$ (Standard deviation)')
        ax2.set_title('Morris Sensitivity: Linear vs Non-linear Effects')
        ax2.grid(True, alpha=0.3)

        # Add diagonal reference line
        max_val = max(max(mu_star), max(sigma)) * 1.1
        ax2.plot([0, max_val], [0, max_val], 'k--', alpha=0.3,
                 label=r'$\sigma = \mu^*$')
        ax2.legend()

        plt.tight_layout()

        if output_path:
            plt.savefig(output_path, dpi=150, bbox_inches='tight')
            logger.info(f"Morris sensitivity plot saved to {output_path}")

        if show:
            plt.show()

        return fig

    def _plot_sobol_results(self,
                            output_path: Optional[Path] = None,
                            show: bool = False) -> Optional[Any]:
        """Plot Sobol sensitivity indices."""
        results = self.sensitivity_results

        param_names = results.get("parameter_names", [])
        S1 = np.array(results.get("S1", []))
        ST = np.array(results.get("ST", []))

        if len(param_names) == 0:
            return None

        fig, ax = plt.subplots(figsize=(12, 6))

        x = np.arange(len(param_names))
        width = 0.35

        # Sort by total effect
        sorted_idx = np.argsort(ST)[::-1]

        bars1 = ax.bar(x - width/2, S1[sorted_idx], width, label='First-order (S1)',
                       color='steelblue')
        bars2 = ax.bar(x + width/2, ST[sorted_idx], width, label='Total effect (ST)',
                       color='coral')

        ax.set_xlabel('Parameter')
        ax.set_ylabel('Sensitivity Index')
        ax.set_title('Sobol Sensitivity Indices')
        ax.set_xticks(x)
        ax.set_xticklabels([param_names[i] for i in sorted_idx], rotation=45, ha='right')
        ax.legend()
        ax.grid(True, alpha=0.3, axis='y')

        # Add interaction effect annotation
        interaction = ST[sorted_idx] - S1[sorted_idx]
        for i, (s1, st, inter) in enumerate(zip(S1[sorted_idx], ST[sorted_idx], interaction)):
            if inter > 0.05:  # Only annotate significant interactions
                ax.annotate(f'{inter:.2f}', (i, st), xytext=(0, 5),
                           textcoords='offset points', ha='center', fontsize=8,
                           color='darkred')

        plt.tight_layout()

        if output_path:
            plt.savefig(output_path, dpi=150, bbox_inches='tight')
            logger.info(f"Sobol sensitivity plot saved to {output_path}")

        if show:
            plt.show()

        return fig

    def plot_parameter_correlations(self,
                                    output_path: Optional[Path] = None,
                                    show: bool = False) -> Optional[Any]:
        """
        Plot parameter correlation matrix.

        Parameters
        ----------
        output_path : Path, optional
            Path to save plot
        show : bool
            Whether to display plot

        Returns
        -------
        Figure or None
            Matplotlib figure if available
        """
        if not HAS_MATPLOTLIB:
            logger.warning("matplotlib not available - cannot create plot")
            return None

        if not self.history:
            return None

        # Build parameter DataFrame
        param_data = {}
        for entry in self.history:
            params = entry.get("parameters", {})
            for name, value in params.items():
                if name not in param_data:
                    param_data[name] = []
                param_data[name].append(value)

        df = pd.DataFrame(param_data)

        # Compute correlation matrix
        corr = df.corr()

        fig, ax = plt.subplots(figsize=(10, 8))

        im = ax.imshow(corr, cmap='RdBu_r', vmin=-1, vmax=1)

        # Add colorbar
        plt.colorbar(im, ax=ax, label='Correlation')

        # Add labels
        ax.set_xticks(range(len(corr.columns)))
        ax.set_yticks(range(len(corr.index)))
        ax.set_xticklabels(corr.columns, rotation=45, ha='right')
        ax.set_yticklabels(corr.index)

        # Add correlation values
        for i in range(len(corr.index)):
            for j in range(len(corr.columns)):
                text = ax.text(j, i, f'{corr.iloc[i, j]:.2f}',
                              ha='center', va='center', fontsize=8)

        ax.set_title('Parameter Correlation Matrix')

        plt.tight_layout()

        if output_path:
            plt.savefig(output_path, dpi=150, bbox_inches='tight')
            logger.info(f"Correlation plot saved to {output_path}")

        if show:
            plt.show()

        return fig

    def generate_report(self, output_dir: Path) -> Path:
        """
        Generate comprehensive calibration report.

        Parameters
        ----------
        output_dir : Path
            Output directory for report

        Returns
        -------
        Path
            Path to report file
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # Generate summary
        try:
            summary = self.get_summary()
            summary.save_json(output_dir / "calibration_summary.json")
        except ValueError as e:
            logger.warning(f"Could not generate summary: {e}")
            summary = None

        # Generate plots
        if HAS_MATPLOTLIB:
            self.plot_convergence(output_dir / "convergence.png")
            self.plot_parameter_evolution(output_path=output_dir / "parameter_evolution.png")
            self.plot_parameter_correlations(output_dir / "parameter_correlations.png")

            if self.sensitivity_results:
                self.plot_sensitivity_results(output_dir / "sensitivity_results.png")

        # Generate text report
        report_path = output_dir / "calibration_report.txt"
        with open(report_path, 'w') as f:
            f.write("=" * 70 + "\n")
            f.write("OpenWQ Calibration Report\n")
            f.write("=" * 70 + "\n\n")

            if summary:
                f.write("SUMMARY\n")
                f.write("-" * 40 + "\n")
                f.write(f"Total evaluations: {summary.n_evaluations}\n")
                f.write(f"Best objective: {summary.best_objective:.6f}\n")
                f.write(f"Convergence at evaluation: {summary.convergence_eval}\n")
                f.write(f"Total runtime: {summary.runtime_hours:.2f} hours\n\n")

                f.write("BEST PARAMETERS\n")
                f.write("-" * 40 + "\n")
                for name, value in sorted(summary.best_parameters.items()):
                    bounds = summary.parameter_bounds.get(name, (None, None))
                    f.write(f"  {name}: {value:.6g}")
                    if bounds[0] is not None:
                        f.write(f"  (bounds: [{bounds[0]:.4g}, {bounds[1]:.4g}])")
                    f.write("\n")

            # Add objective statistics
            obj_stats = self.get_objective_statistics()
            if obj_stats:
                f.write("\nOBJECTIVE FUNCTION STATISTICS\n")
                f.write("-" * 40 + "\n")
                for key, value in obj_stats.items():
                    f.write(f"  {key}: {value:.6f}\n")

            # Add parameter statistics
            param_stats = self.get_parameter_statistics()
            if not param_stats.empty:
                f.write("\nPARAMETER STATISTICS\n")
                f.write("-" * 40 + "\n")
                f.write(param_stats.to_string(index=False))
                f.write("\n")

            f.write("\n" + "=" * 70 + "\n")
            f.write("End of Report\n")

        logger.info(f"Report generated at {output_dir}")
        return report_path

    def export_to_csv(self, output_path: Path):
        """
        Export calibration history to CSV.

        Parameters
        ----------
        output_path : Path
            Output CSV file path
        """
        if not self.history:
            logger.warning("No history to export")
            return

        # Flatten history to DataFrame
        rows = []
        for entry in self.history:
            row = {"evaluation": entry.get("eval_id", 0),
                   "objective": entry.get("objective")}
            row.update(entry.get("parameters", {}))
            rows.append(row)

        df = pd.DataFrame(rows)
        df.to_csv(output_path, index=False)
        logger.info(f"History exported to {output_path}")

    def plot_timeseries_comparison(self,
                                   matched_data: pd.DataFrame,
                                   species: str = None,
                                   reach_ids: list = None,
                                   temporal_resolution: str = "native",
                                   output_path: Optional[Path] = None,
                                   show: bool = False) -> Optional[Any]:
        """
        Plot time series comparison of observed vs simulated values.

        Parameters
        ----------
        matched_data : pd.DataFrame
            DataFrame with datetime, reach_id, species, observed, simulated columns
        species : str, optional
            Filter to specific species (default: all species)
        reach_ids : list, optional
            Filter to specific reach IDs (default: all reaches)
        temporal_resolution : str
            Temporal resolution label for title
        output_path : Path, optional
            Path to save plot
        show : bool
            Whether to display plot

        Returns
        -------
        Figure or None
            Matplotlib figure if available
        """
        if not HAS_MATPLOTLIB:
            logger.warning("matplotlib not available - cannot create plot")
            return None

        if matched_data.empty:
            logger.warning("No matched data for time series plot")
            return None

        data = matched_data.copy()
        data['datetime'] = pd.to_datetime(data['datetime'])

        # Filter by species
        if species:
            data = data[data['species'] == species]

        # Filter by reach_ids
        if reach_ids:
            data = data[data['reach_id'].isin(reach_ids)]

        if data.empty:
            logger.warning("No data after filtering")
            return None

        # Get unique species for subplots
        unique_species = data['species'].unique()
        n_species = len(unique_species)

        fig, axes = plt.subplots(n_species, 1, figsize=(14, 4 * n_species), squeeze=False)

        for i, sp in enumerate(unique_species):
            ax = axes[i, 0]
            sp_data = data[data['species'] == sp].sort_values('datetime')

            # Plot observed
            ax.scatter(sp_data['datetime'], sp_data['observed'],
                      color='blue', alpha=0.6, s=30, label='Observed', marker='o')

            # Plot simulated
            ax.scatter(sp_data['datetime'], sp_data['simulated'],
                      color='red', alpha=0.6, s=30, label='Simulated', marker='x')

            # If multiple reaches, show as different colors
            unique_reaches = sp_data['reach_id'].unique()
            if len(unique_reaches) <= 10:
                colors = plt.cm.tab10(np.linspace(0, 1, len(unique_reaches)))
                for j, reach in enumerate(unique_reaches):
                    reach_data = sp_data[sp_data['reach_id'] == reach].sort_values('datetime')
                    if len(reach_data) > 1:
                        ax.plot(reach_data['datetime'], reach_data['simulated'],
                               color=colors[j], alpha=0.3, linewidth=1)

            # Compute metrics for this species
            obs = sp_data['observed'].values
            sim = sp_data['simulated'].values
            valid_mask = ~(np.isnan(obs) | np.isnan(sim))
            if valid_mask.sum() >= 2:
                obs_valid = obs[valid_mask]
                sim_valid = sim[valid_mask]

                # Import metrics from objective_functions
                from ..objective_functions import ObjectiveFunction
                rmse = ObjectiveFunction.rmse(obs_valid, sim_valid)
                nse = ObjectiveFunction.nse(obs_valid, sim_valid)
                kge = ObjectiveFunction.kge(obs_valid, sim_valid)
                pbias = ObjectiveFunction.pbias(obs_valid, sim_valid)

                metrics_text = f'RMSE={rmse:.3f}  NSE={nse:.3f}  KGE={kge:.3f}  PBIAS={pbias:.1f}%'
                ax.text(0.02, 0.98, metrics_text, transform=ax.transAxes,
                       fontsize=10, verticalalignment='top',
                       bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

            ax.set_xlabel('Date', fontsize=11)
            ax.set_ylabel(f'{sp} Concentration', fontsize=11)
            ax.set_title(f'{sp} - Time Series Comparison ({temporal_resolution})', fontsize=12)
            ax.legend(loc='upper right')
            ax.grid(True, alpha=0.3)

            # Format x-axis dates
            fig.autofmt_xdate()

        plt.tight_layout()

        if output_path:
            plt.savefig(output_path, dpi=150, bbox_inches='tight')
            logger.info(f"Time series plot saved to {output_path}")

        if show:
            plt.show()

        return fig

    def plot_scatter_obs_sim(self,
                             matched_data: pd.DataFrame,
                             species: str = None,
                             temporal_resolution: str = "native",
                             output_path: Optional[Path] = None,
                             show: bool = False) -> Optional[Any]:
        """
        Plot scatter plot of observed vs simulated values (1:1 plot).

        Parameters
        ----------
        matched_data : pd.DataFrame
            DataFrame with observed and simulated columns
        species : str, optional
            Filter to specific species
        temporal_resolution : str
            Temporal resolution label for title
        output_path : Path, optional
            Path to save plot
        show : bool
            Whether to display plot

        Returns
        -------
        Figure or None
            Matplotlib figure if available
        """
        if not HAS_MATPLOTLIB:
            logger.warning("matplotlib not available - cannot create plot")
            return None

        if matched_data.empty:
            logger.warning("No matched data for scatter plot")
            return None

        data = matched_data.copy()

        # Filter by species
        if species:
            data = data[data['species'] == species]

        if data.empty:
            return None

        # Get unique species for subplots
        unique_species = data['species'].unique()
        n_species = len(unique_species)
        n_cols = min(3, n_species)
        n_rows = int(np.ceil(n_species / n_cols))

        fig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 5 * n_rows), squeeze=False)

        for i, sp in enumerate(unique_species):
            row = i // n_cols
            col = i % n_cols
            ax = axes[row, col]

            sp_data = data[data['species'] == sp]
            obs = sp_data['observed'].values
            sim = sp_data['simulated'].values

            # Remove NaN pairs
            valid_mask = ~(np.isnan(obs) | np.isnan(sim))
            obs = obs[valid_mask]
            sim = sim[valid_mask]

            if len(obs) < 2:
                ax.set_visible(False)
                continue

            # Scatter plot
            ax.scatter(obs, sim, alpha=0.5, s=30, c='steelblue')

            # 1:1 line
            min_val = min(obs.min(), sim.min())
            max_val = max(obs.max(), sim.max())
            ax.plot([min_val, max_val], [min_val, max_val], 'k--', linewidth=1.5, label='1:1 line')

            # Linear regression
            if len(obs) >= 3:
                z = np.polyfit(obs, sim, 1)
                p = np.poly1d(z)
                ax.plot([min_val, max_val], [p(min_val), p(max_val)],
                       'r-', linewidth=1, alpha=0.7, label=f'Fit: y={z[0]:.2f}x+{z[1]:.2f}')

            # Compute metrics
            from ..objective_functions import ObjectiveFunction
            rmse = ObjectiveFunction.rmse(obs, sim)
            nse = ObjectiveFunction.nse(obs, sim)
            kge = ObjectiveFunction.kge(obs, sim)
            r2 = np.corrcoef(obs, sim)[0, 1] ** 2

            metrics_text = f'n={len(obs)}\nRMSE={rmse:.3f}\nNSE={nse:.3f}\nKGE={kge:.3f}\nR²={r2:.3f}'
            ax.text(0.05, 0.95, metrics_text, transform=ax.transAxes,
                   fontsize=9, verticalalignment='top',
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

            ax.set_xlabel('Observed', fontsize=11)
            ax.set_ylabel('Simulated', fontsize=11)
            ax.set_title(f'{sp} ({temporal_resolution})', fontsize=12)
            ax.legend(loc='lower right', fontsize=9)
            ax.grid(True, alpha=0.3)
            ax.set_aspect('equal', adjustable='box')

        # Hide empty subplots
        for i in range(n_species, n_rows * n_cols):
            row = i // n_cols
            col = i % n_cols
            axes[row, col].set_visible(False)

        plt.suptitle('Observed vs Simulated Comparison', fontsize=14, y=1.02)
        plt.tight_layout()

        if output_path:
            plt.savefig(output_path, dpi=150, bbox_inches='tight')
            logger.info(f"Scatter plot saved to {output_path}")

        if show:
            plt.show()

        return fig

    def plot_residuals(self,
                       matched_data: pd.DataFrame,
                       species: str = None,
                       temporal_resolution: str = "native",
                       output_path: Optional[Path] = None,
                       show: bool = False) -> Optional[Any]:
        """
        Plot residuals (simulated - observed) analysis.

        Parameters
        ----------
        matched_data : pd.DataFrame
            DataFrame with observed, simulated, datetime columns
        species : str, optional
            Filter to specific species
        temporal_resolution : str
            Temporal resolution label for title
        output_path : Path, optional
            Path to save plot
        show : bool
            Whether to display plot

        Returns
        -------
        Figure or None
            Matplotlib figure if available
        """
        if not HAS_MATPLOTLIB:
            logger.warning("matplotlib not available - cannot create plot")
            return None

        if matched_data.empty:
            return None

        data = matched_data.copy()
        data['datetime'] = pd.to_datetime(data['datetime'])

        if species:
            data = data[data['species'] == species]

        if data.empty:
            return None

        # Calculate residuals
        data['residual'] = data['simulated'] - data['observed']

        unique_species = data['species'].unique()
        n_species = len(unique_species)

        fig, axes = plt.subplots(n_species, 2, figsize=(12, 4 * n_species), squeeze=False)

        for i, sp in enumerate(unique_species):
            sp_data = data[data['species'] == sp].dropna(subset=['residual'])

            if sp_data.empty:
                continue

            residuals = sp_data['residual'].values

            # Left plot: Residuals over time
            ax1 = axes[i, 0]
            ax1.scatter(sp_data['datetime'], residuals, alpha=0.5, s=20, c='steelblue')
            ax1.axhline(0, color='red', linestyle='--', linewidth=1)
            ax1.set_xlabel('Date', fontsize=11)
            ax1.set_ylabel('Residual (Sim - Obs)', fontsize=11)
            ax1.set_title(f'{sp} - Residuals Over Time ({temporal_resolution})', fontsize=12)
            ax1.grid(True, alpha=0.3)

            # Right plot: Residuals histogram
            ax2 = axes[i, 1]
            ax2.hist(residuals, bins=30, color='steelblue', alpha=0.7, edgecolor='white')
            ax2.axvline(0, color='red', linestyle='--', linewidth=1.5)
            ax2.axvline(np.mean(residuals), color='green', linestyle='-', linewidth=1.5,
                       label=f'Mean={np.mean(residuals):.3f}')
            ax2.set_xlabel('Residual (Sim - Obs)', fontsize=11)
            ax2.set_ylabel('Frequency', fontsize=11)
            ax2.set_title(f'{sp} - Residual Distribution', fontsize=12)
            ax2.legend()
            ax2.grid(True, alpha=0.3)

        fig.autofmt_xdate()
        plt.tight_layout()

        if output_path:
            plt.savefig(output_path, dpi=150, bbox_inches='tight')
            logger.info(f"Residuals plot saved to {output_path}")

        if show:
            plt.show()

        return fig

    def generate_performance_plots(self,
                                   matched_data: pd.DataFrame,
                                   output_dir: Path,
                                   temporal_resolution: str = "native"):
        """
        Generate all performance analysis plots.

        Parameters
        ----------
        matched_data : pd.DataFrame
            DataFrame with matched observation-simulation pairs
        output_dir : Path
            Directory to save plots
        temporal_resolution : str
            Temporal resolution label
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        if HAS_MATPLOTLIB and not matched_data.empty:
            self.plot_timeseries_comparison(
                matched_data,
                temporal_resolution=temporal_resolution,
                output_path=output_dir / "timeseries_comparison.png"
            )

            self.plot_scatter_obs_sim(
                matched_data,
                temporal_resolution=temporal_resolution,
                output_path=output_dir / "scatter_obs_sim.png"
            )

            self.plot_residuals(
                matched_data,
                temporal_resolution=temporal_resolution,
                output_path=output_dir / "residuals_analysis.png"
            )

            logger.info(f"Performance plots generated in {output_dir}")

    def generate_html_report(self,
                              output_dir: Path,
                              matched_data: pd.DataFrame = None,
                              config_summary: Dict = None,
                              temporal_resolution: str = "native",
                              shapefile_dir: str = None,
                              observation_data_path: str = None) -> Path:
        """
        Generate a modern, self-contained HTML report with interactive
        Plotly.js charts, an interactive Leaflet.js basin map, dark / light
        theme toggle, sticky sidebar navigation with scroll-spy, and
        scroll-triggered animations.

        Parameters
        ----------
        output_dir : Path
            Directory to save report (and where plots already exist)
        matched_data : pd.DataFrame, optional
            Matched observation-simulation pairs from best evaluation
        config_summary : Dict, optional
            Dictionary of calibration configuration settings
        temporal_resolution : str
            Temporal resolution label
        shapefile_dir : str, optional
            Path to directory containing basin GeoPackages
            (*_basinHru.gpkg, *_riverNetwork.gpkg) and optionally
            *_lulc_areas.csv for the map section
        observation_data_path : str, optional
            Path to calibration_observations.csv (for station markers)

        Returns
        -------
        Path
            Path to the generated HTML file
        """
        from datetime import datetime

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # ── JSON helper (handles numpy / nan) ─────────────────────────
        def _js(obj):
            """Serialise *obj* to a JSON string safe for embedding."""
            import math
            def _conv(o):
                if isinstance(o, (np.integer,)):
                    return int(o)
                if isinstance(o, (np.floating,)):
                    v = float(o)
                    return None if math.isnan(v) or math.isinf(v) else v
                if isinstance(o, np.ndarray):
                    return [_conv(x) for x in o.tolist()]
                if isinstance(o, float):
                    return None if math.isnan(o) or math.isinf(o) else o
                return o
            if isinstance(obj, (list, tuple)):
                obj = [_conv(x) for x in obj]
            elif isinstance(obj, dict):
                obj = {k: _conv(v) for k, v in obj.items()}
            else:
                obj = _conv(obj)
            return json.dumps(obj)

        # ── Gather data ───────────────────────────────────────────────
        config_summary = config_summary or {}
        now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        try:
            summary = self.get_summary()
        except ValueError:
            summary = None

        obj_stats = self.get_objective_statistics()
        param_stats = self.get_parameter_statistics()

        # ── Per-species metrics ───────────────────────────────────────
        species_metrics = []
        if matched_data is not None and not matched_data.empty:
            from ..objective_functions import ObjectiveFunction
            for sp in sorted(matched_data['species'].unique()):
                sp_data = matched_data[matched_data['species'] == sp]
                obs = sp_data['observed'].values
                sim = sp_data['simulated'].values
                valid = ~(np.isnan(obs) | np.isnan(sim))
                obs_v, sim_v = obs[valid], sim[valid]
                if len(obs_v) >= 2:
                    species_metrics.append({
                        "species": sp, "n": int(len(obs_v)),
                        "KGE": float(ObjectiveFunction.kge(obs_v, sim_v)),
                        "NSE": float(ObjectiveFunction.nse(obs_v, sim_v)),
                        "RMSE": float(ObjectiveFunction.rmse(obs_v, sim_v)),
                        "PBIAS": float(ObjectiveFunction.pbias(obs_v, sim_v)),
                        "R2": float(np.corrcoef(obs_v, sim_v)[0, 1] ** 2),
                        "obs_mean": float(np.mean(obs_v)),
                        "sim_mean": float(np.mean(sim_v)),
                    })

        # ── Detect which sections will appear ─────────────────────────
        has_history = bool(self.history)
        has_correlations = has_history and len(self.parameter_definitions) >= 2
        has_sensitivity = bool(self.sensitivity_results)
        has_matched = matched_data is not None and not matched_data.empty

        # ── Load basin spatial data for map (if available) ────────────
        map_data = {}  # will hold GeoJSON strings + basin stats
        if shapefile_dir and HAS_FIONA:
            shapefile_dir_p = Path(shapefile_dir)
            net_features = []  # keep for station matching
            try:
                # Find GeoPackage files
                hru_files = list(shapefile_dir_p.glob("*_basinHru.gpkg"))
                net_files = list(shapefile_dir_p.glob("*_riverNetwork.gpkg"))

                if hru_files:
                    with fiona.open(str(hru_files[0])) as src:
                        features = list(src)
                        schema = src.schema
                    # Build GeoJSON FeatureCollection
                    hru_fc = {"type": "FeatureCollection",
                              "features": features}
                    map_data["hru_geojson"] = json.dumps(hru_fc)
                    map_data["n_hrus"] = len(features)

                    # Compute area stats & bounds from features
                    props = schema["properties"]
                    has_unitarea = "unitarea" in props
                    has_area = "area" in props
                    total_area = 0.0
                    min_x, min_y = 180.0, 90.0
                    max_x, max_y = -180.0, -90.0
                    for f in features:
                        p = f["properties"]
                        if has_unitarea and p.get("unitarea"):
                            total_area += float(p["unitarea"])
                        elif has_area and p.get("area"):
                            total_area += float(p["area"]) / 1e6
                        # Walk geometry coords for bounds
                        geom = f["geometry"]
                        coords = geom.get("coordinates", [])

                        def _walk(c):
                            """Recursively walk nested coord lists."""
                            if (isinstance(c, (list, tuple))
                                    and len(c) >= 2
                                    and isinstance(c[0], (int, float))):
                                return [c]
                            pts = []
                            for item in c:
                                pts.extend(_walk(item))
                            return pts

                        for pt in _walk(coords):
                            x, y = float(pt[0]), float(pt[1])
                            min_x, min_y = min(min_x, x), min(min_y, y)
                            max_x, max_y = max(max_x, x), max(max_y, y)

                    if total_area > 0:
                        map_data["total_area_km2"] = total_area
                    map_data["center_lat"] = (min_y + max_y) / 2
                    map_data["center_lon"] = (min_x + max_x) / 2
                    map_data["bounds"] = [[min_y, min_x],
                                          [max_y, max_x]]
                    logger.info(f"Loaded {len(features)} HRUs from "
                                f"{hru_files[0].name}")

                if net_files:
                    with fiona.open(str(net_files[0])) as src:
                        net_features = list(src)
                        net_schema = src.schema
                    net_fc = {"type": "FeatureCollection",
                              "features": net_features}
                    map_data["network_geojson"] = json.dumps(net_fc)
                    map_data["n_reaches"] = len(net_features)

                    props = net_schema["properties"]
                    total_len = 0.0
                    max_ord = 0
                    for f in net_features:
                        p = f["properties"]
                        if "Length" in props and p.get("Length"):
                            total_len += float(p["Length"])
                        elif "length" in props and p.get("length"):
                            total_len += float(p["length"])
                        if "order" in props and p.get("order"):
                            max_ord = max(max_ord, int(p["order"]))
                    if total_len > 0:
                        map_data["total_length_km"] = total_len / 1000
                    if max_ord > 0:
                        map_data["max_order"] = max_ord
                    logger.info(f"Loaded {len(net_features)} reaches "
                                f"from {net_files[0].name}")

                # Load observation stations (if available)
                if observation_data_path and net_features:
                    obs_path = Path(observation_data_path)
                    if obs_path.exists():
                        obs_df = pd.read_csv(obs_path)
                        if "reach_id" in obs_df.columns:
                            station_ids = obs_df["reach_id"].unique()
                            stations = []
                            for sid in station_ids:
                                # Match reach by LINKNO
                                for nf in net_features:
                                    if nf["properties"].get("LINKNO") == sid:
                                        # Centroid of first coord pair
                                        coords = nf["geometry"]["coordinates"]
                                        if isinstance(coords[0], (list, tuple)):
                                            mid = len(coords) // 2
                                            pt = coords[mid]
                                        else:
                                            pt = coords
                                        sp_list = obs_df[
                                            obs_df["reach_id"] == sid
                                        ]["species"].unique().tolist()
                                        n_obs = len(obs_df[
                                            obs_df["reach_id"] == sid])
                                        stations.append({
                                            "lat": float(pt[1]),
                                            "lon": float(pt[0]),
                                            "reach_id": int(sid),
                                            "species": sp_list,
                                            "n_obs": n_obs,
                                        })
                                        break
                            if stations:
                                map_data["stations"] = stations
                                logger.info(
                                    f"Loaded {len(stations)} observation "
                                    f"stations")

            except Exception as e:
                logger.warning(f"Could not load spatial data: {e}")
                import traceback
                logger.debug(traceback.format_exc())

        has_map = bool(map_data.get("hru_geojson") or
                       map_data.get("network_geojson"))

        # ── Build TOC items list (reused by sidebar + JS) ─────────────
        toc = [("summary", "Summary")]
        if has_map:
            toc.append(("basin-map", "Basin Map"))
        toc += [("config", "Configuration"),
                ("parameters", "Parameters")]
        if has_history:
            toc += [("convergence", "Convergence"),
                    ("param-evolution", "Param Evolution")]
        if has_correlations:
            toc.append(("correlations", "Correlations"))
        if has_sensitivity:
            toc.append(("sensitivity", "Sensitivity"))
        if species_metrics:
            toc.append(("performance", "Performance"))
        if has_matched:
            toc += [("timeseries", "Time Series"),
                    ("scatter", "Obs vs Sim"),
                    ("residuals", "Residuals")]
        if obj_stats:
            toc.append(("obj-stats", "Obj. Stats"))
        if not param_stats.empty:
            toc.append(("param-stats", "Param Stats"))
        if has_history:
            toc.append(("eval-history", "History"))

        # ── Header metadata ───────────────────────────────────────────
        basin_id = config_summary.get("basin_id", "")
        variant = config_summary.get("variant", "")
        algo = config_summary.get("algorithm", "DDS")
        title_extra = ""
        if basin_id:
            title_extra += f" &mdash; {basin_id}"
        if variant:
            title_extra += f" (Variant {variant})"

        n_evals = summary.n_evaluations if summary else len(self.history)
        best_obj = summary.best_objective if summary else obj_stats.get("best", float('nan'))
        conv_eval = summary.convergence_eval if summary else 0
        runtime_h = summary.runtime_hours if summary else 0

        # Plotly shared config
        PCFG = ('{responsive:true,displayModeBar:true,'
                'modeBarButtonsToRemove:["lasso2d","select2d"],'
                'toImageButtonOptions:{format:"svg",filename:"openwq_plot"}}')

        # ==============================================================
        # Build HTML
        # ==============================================================
        H = []

        # ── <head> ────────────────────────────────────────────────────
        H.append("""<!DOCTYPE html>
<html lang="en" data-theme="light">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>OpenWQ Calibration Report</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&family=JetBrains+Mono:wght@400;500&display=swap" rel="stylesheet">
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js" charset="utf-8"></script>
<link rel="stylesheet" href="https://unpkg.com/leaflet@1.9.4/dist/leaflet.css" crossorigin="">
<script src="https://unpkg.com/leaflet@1.9.4/dist/leaflet.js" crossorigin=""></script>
<style>
/* ─── Theme tokens ──────────────────────────────────── */
:root,[data-theme="light"]{
  --bg:#f0f2f5;--bg2:#fff;--card:#fff;--card-hover:#fff;
  --border:#e0e3e8;--border2:#d0d4da;
  --text:#1a1a2e;--text2:#4a4a68;--text3:#8a8aa0;
  --accent:#6366f1;--accent2:#8b5cf6;--accent-soft:rgba(99,102,241,.08);
  --green:#10b981;--red:#ef4444;--amber:#f59e0b;
  --plot-bg:#fff;--plot-grid:#f0f0f0;--plot-font:#333;
  --shadow:0 1px 3px rgba(0,0,0,.06),0 1px 2px rgba(0,0,0,.04);
  --shadow-lg:0 10px 25px rgba(0,0,0,.07);
  --glass:rgba(255,255,255,.7);--glass-border:rgba(255,255,255,.3);
  --sidebar-bg:rgba(255,255,255,.85);
  --header-from:#6366f1;--header-to:#8b5cf6;
  --th-bg:#6366f1;
  --bar-track:#e5e7eb;--bar-fill:linear-gradient(90deg,#6366f1,#8b5cf6);
  --kpi-icon:#6366f1;
}
[data-theme="dark"]{
  --bg:#0f1117;--bg2:#1a1b26;--card:#1e1f2e;--card-hover:#252640;
  --border:#2a2b3d;--border2:#3a3b4d;
  --text:#e2e2f0;--text2:#a0a0b8;--text3:#6a6a80;
  --accent:#818cf8;--accent2:#a78bfa;--accent-soft:rgba(129,140,248,.12);
  --green:#34d399;--red:#f87171;--amber:#fbbf24;
  --plot-bg:#1e1f2e;--plot-grid:#2d2e42;--plot-font:#ccc;
  --shadow:0 1px 3px rgba(0,0,0,.3);
  --shadow-lg:0 10px 25px rgba(0,0,0,.4);
  --glass:rgba(30,31,46,.75);--glass-border:rgba(255,255,255,.06);
  --sidebar-bg:rgba(15,17,23,.92);
  --header-from:#4f46e5;--header-to:#7c3aed;
  --th-bg:#4338ca;
  --bar-track:#374151;--bar-fill:linear-gradient(90deg,#818cf8,#a78bfa);
  --kpi-icon:#818cf8;
}

/* ─── Reset & base ──────────────────────────────────── */
*{box-sizing:border-box;margin:0;padding:0}
html{scroll-behavior:smooth}
body{
  font-family:'Inter',system-ui,-apple-system,sans-serif;
  background:var(--bg);color:var(--text);line-height:1.65;
  transition:background .3s,color .3s;
}

/* ─── Layout: sidebar + main ────────────────────────── */
.layout{display:flex;min-height:100vh}
.sidebar{
  position:sticky;top:0;height:100vh;width:220px;flex-shrink:0;
  background:var(--sidebar-bg);backdrop-filter:blur(12px);
  border-right:1px solid var(--border);padding:1.2rem .8rem;
  overflow-y:auto;z-index:100;
  transition:background .3s,border .3s;
}
.sidebar .logo{font-weight:700;font-size:.85rem;color:var(--accent);margin-bottom:1rem;letter-spacing:-.3px}
.sidebar nav a{
  display:block;padding:.35rem .7rem;margin-bottom:2px;
  border-radius:6px;font-size:.8rem;color:var(--text2);
  text-decoration:none;transition:all .15s;white-space:nowrap;
  overflow:hidden;text-overflow:ellipsis;
}
.sidebar nav a:hover{background:var(--accent-soft);color:var(--accent)}
.sidebar nav a.active{background:var(--accent-soft);color:var(--accent);font-weight:600}
.sidebar .theme-toggle{
  margin-top:auto;padding-top:1rem;border-top:1px solid var(--border);
}
.theme-btn{
  display:flex;align-items:center;gap:.5rem;width:100%;
  padding:.45rem .7rem;border:1px solid var(--border);border-radius:8px;
  background:var(--card);color:var(--text2);font-size:.78rem;
  cursor:pointer;transition:all .15s;font-family:inherit;
}
.theme-btn:hover{border-color:var(--accent);color:var(--accent)}
.theme-btn svg{width:16px;height:16px;flex-shrink:0}
.main{flex:1;min-width:0;padding:0}

/* ─── Header ────────────────────────────────────────── */
.header{
  background:linear-gradient(135deg,var(--header-from),var(--header-to));
  color:#fff;padding:2.5rem 2.5rem 2rem;
  position:relative;overflow:hidden;
}
.header::before{
  content:'';position:absolute;inset:0;
  background:radial-gradient(circle at 80% 20%,rgba(255,255,255,.12) 0%,transparent 60%);
}
.header h1{font-size:1.7rem;font-weight:700;letter-spacing:-.4px;position:relative}
.header .subtitle{opacity:.8;font-size:.88rem;margin-top:.3rem;position:relative}
.header .meta{
  display:flex;gap:1.5rem;margin-top:.9rem;font-size:.78rem;
  opacity:.7;flex-wrap:wrap;position:relative;
}
.header .meta span{
  background:rgba(255,255,255,.12);padding:.2rem .6rem;border-radius:20px;
}

/* ─── Container ─────────────────────────────────────── */
.container{max-width:1100px;margin:0 auto;padding:2rem 2rem 1rem}

/* ─── Sections ──────────────────────────────────────── */
.section{margin-bottom:2.5rem;opacity:0;transform:translateY(16px);transition:opacity .5s ease,transform .5s ease}
.section.visible{opacity:1;transform:translateY(0)}
.section h2{
  font-size:1.15rem;font-weight:700;color:var(--text);
  margin-bottom:1rem;display:flex;align-items:center;gap:.5rem;
  letter-spacing:-.3px;
}
.section h2::before{
  content:'';width:4px;height:1.2em;border-radius:2px;
  background:linear-gradient(180deg,var(--accent),var(--accent2));flex-shrink:0;
}

/* ─── Cards ─────────────────────────────────────────── */
.card{
  background:var(--card);border:1px solid var(--border);
  border-radius:12px;padding:1.3rem;margin-bottom:1rem;
  box-shadow:var(--shadow);transition:all .3s;
}
.card:hover{box-shadow:var(--shadow-lg);border-color:var(--border2)}

/* ─── KPI grid ──────────────────────────────────────── */
.kpi-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(170px,1fr));gap:.8rem}
.kpi{
  text-align:center;padding:1.1rem .8rem;
  background:var(--card);border:1px solid var(--border);border-radius:12px;
  box-shadow:var(--shadow);transition:all .25s;position:relative;overflow:hidden;
}
.kpi:hover{transform:translateY(-2px);box-shadow:var(--shadow-lg)}
.kpi::before{
  content:'';position:absolute;top:0;left:0;right:0;height:3px;
  background:var(--bar-fill);border-radius:12px 12px 0 0;
}
.kpi .value{font-size:1.5rem;font-weight:700;color:var(--accent);font-family:'JetBrains Mono',monospace}
.kpi .label{font-size:.7rem;color:var(--text3);text-transform:uppercase;letter-spacing:.8px;margin-top:.3rem}

/* ─── Tables ────────────────────────────────────────── */
.table-wrap{overflow-x:auto;border-radius:10px;border:1px solid var(--border)}
table{width:100%;border-collapse:collapse;font-size:.85rem}
th,td{padding:.55rem .75rem;text-align:left;border-bottom:1px solid var(--border)}
th{
  background:var(--th-bg);color:#fff;font-weight:600;font-size:.72rem;
  text-transform:uppercase;letter-spacing:.5px;position:sticky;top:0;
}
tr{transition:background .15s}
tr:nth-child(even){background:var(--accent-soft)}
tr:hover{background:var(--accent-soft)}
td.num{text-align:right;font-family:'JetBrains Mono',monospace;font-size:.8rem}

/* ─── Plot captions ─────────────────────────────────── */
.plot-caption{font-size:.75rem;color:var(--text3);margin-top:.4rem;text-align:center;font-style:italic}

/* ─── Badges ────────────────────────────────────────── */
.badge{
  display:inline-flex;align-items:center;gap:.25rem;
  padding:.2rem .55rem;border-radius:20px;font-size:.7rem;font-weight:600;
}
.badge-good{background:rgba(16,185,129,.12);color:var(--green)}
.badge-fair{background:rgba(245,158,11,.12);color:var(--amber)}
.badge-poor{background:rgba(239,68,68,.12);color:var(--red)}

/* ─── Config grid ───────────────────────────────────── */
.config-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(260px,1fr));gap:.3rem;font-size:.84rem}
.config-item{display:flex;justify-content:space-between;padding:.35rem 0;border-bottom:1px dotted var(--border)}
.config-item .key{color:var(--text3)}
.config-item .val{font-weight:600;color:var(--text);font-family:'JetBrains Mono',monospace;font-size:.8rem}

/* ─── Range bar ─────────────────────────────────────── */
.range-bar{display:inline-flex;align-items:center;gap:.4rem}
.range-track{width:90px;height:8px;background:var(--bar-track);border-radius:4px;overflow:hidden}
.range-fill{height:100%;border-radius:4px;background:var(--bar-fill);transition:width .6s ease}
.range-pct{font-size:.7rem;color:var(--text3);font-family:'JetBrains Mono',monospace}

/* ─── Collapsible ───────────────────────────────────── */
.collapsible-header{
  display:flex;align-items:center;justify-content:space-between;cursor:pointer;
  padding:.6rem 0;user-select:none;
}
.collapsible-header .chevron{transition:transform .2s;color:var(--text3)}
.collapsible-header.open .chevron{transform:rotate(180deg)}
.collapsible-body{overflow:hidden;transition:max-height .3s ease;max-height:0}
.collapsible-body.open{max-height:none}

/* ─── Basin Map ─────────────────────────────────────── */
.map-container{border-radius:10px;overflow:hidden;border:1px solid var(--border)}
#basinMap{height:520px;width:100%;background:var(--bg2);z-index:1}
.map-info-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(140px,1fr));gap:.6rem;margin-top:.8rem}
.map-info-item{text-align:center;padding:.6rem;background:var(--accent-soft);border-radius:8px}
.map-info-item .val{font-size:1.1rem;font-weight:700;color:var(--accent);font-family:'JetBrains Mono',monospace}
.map-info-item .lbl{font-size:.68rem;color:var(--text3);text-transform:uppercase;letter-spacing:.5px;margin-top:.15rem}
.leaflet-popup-content-wrapper{border-radius:8px !important;font-family:'Inter',sans-serif;font-size:.82rem}
.leaflet-popup-content{margin:10px 14px !important}
.leaflet-popup-content b{color:#6366f1}
[data-theme="dark"] .leaflet-tile{filter:brightness(.8) contrast(1.1) saturate(.85)}
[data-theme="dark"] .leaflet-popup-content-wrapper{background:#1e1f2e !important;color:#e2e2f0 !important}
[data-theme="dark"] .leaflet-popup-tip{border-top-color:#1e1f2e !important}
.legend-box{background:var(--card);border:1px solid var(--border);border-radius:8px;padding:8px 12px;font-size:.75rem;line-height:1.6;color:var(--text)}
.legend-box i{width:14px;height:14px;display:inline-block;margin-right:6px;border-radius:3px;vertical-align:middle}

/* ─── Footer ────────────────────────────────────────── */
.footer{
  text-align:center;padding:2rem;color:var(--text3);font-size:.78rem;
  border-top:1px solid var(--border);
}
.footer a{color:var(--accent);text-decoration:none}
.footer a:hover{text-decoration:underline}

/* ─── Responsive ────────────────────────────────────── */
@media(max-width:900px){
  .sidebar{display:none}
  .container{padding:1.2rem}
  .kpi-grid{grid-template-columns:repeat(2,1fr)}
  .header{padding:1.5rem}
}
</style>
</head>
<body>
""")

        # ── Sidebar ───────────────────────────────────────────────────
        H.append('<div class="layout">')
        H.append('<aside class="sidebar" id="sidebar">')
        H.append('<div class="logo">OpenWQ Calibration</div>')
        H.append('<nav id="sidebarNav">')
        for anchor, label in toc:
            H.append(f'<a href="#{anchor}" data-section="{anchor}">{label}</a>')
        H.append('</nav>')
        H.append("""<div class="theme-toggle">
<button class="theme-btn" id="themeToggle" aria-label="Toggle theme">
<svg id="themeIcon" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
<circle cx="12" cy="12" r="5"/><path d="M12 1v2M12 21v2M4.22 4.22l1.42 1.42M18.36 18.36l1.42 1.42M1 12h2M21 12h2M4.22 19.78l1.42-1.42M18.36 5.64l1.42-1.42"/>
</svg>
<span id="themeLabel">Dark mode</span>
</button>
</div>
</aside>
<div class="main">
""")

        # ── Header ────────────────────────────────────────────────────
        H.append(f"""
<div class="header">
  <h1>OpenWQ Calibration Report{title_extra}</h1>
  <div class="subtitle">Interactive diagnostics &mdash; zoom, pan, hover on any chart</div>
  <div class="meta">
    <span>{now}</span><span>{algo}</span>
    <span>{temporal_resolution}</span><span>{n_evals} evals</span>
  </div>
</div>
<div class="container">
""")

        # ── KPI Summary ───────────────────────────────────────────────
        H.append('<div class="section" id="summary"><h2>Summary</h2><div class="kpi-grid">')
        kpis = [
            (str(n_evals), "Evaluations"),
            (f"{best_obj:.4g}", "Best Objective"),
            (str(conv_eval), "Best at Eval #"),
            (str(len(self.parameter_definitions)), "Parameters"),
        ]
        if runtime_h > 0:
            kpis.append((f"{runtime_h:.2f}h", "Runtime"))
        for v, l in kpis:
            H.append(f'<div class="kpi"><div class="value">{v}</div><div class="label">{l}</div></div>')
        H.append('</div></div>')

        # ── Basin Map (Leaflet.js) ────────────────────────────────────
        if has_map:
            H.append('<div class="section" id="basin-map"><h2>Basin Map</h2>')
            H.append('<div class="card"><div class="map-container">')
            H.append('<div id="basinMap"></div></div>')

            # Basin info grid below the map
            info_items = []
            if "n_hrus" in map_data:
                info_items.append(
                    (str(map_data["n_hrus"]), "HRUs / GRUs"))
            if "n_reaches" in map_data:
                info_items.append(
                    (str(map_data["n_reaches"]), "River Reaches"))
            if "total_area_km2" in map_data:
                info_items.append(
                    (f"{map_data['total_area_km2']:.1f}",
                     "Area (km&sup2;)"))
            if "total_length_km" in map_data:
                info_items.append(
                    (f"{map_data['total_length_km']:.1f}",
                     "Network (km)"))
            if "max_order" in map_data:
                info_items.append(
                    (str(map_data["max_order"]),
                     "Max Strahler Order"))
            if map_data.get("stations"):
                info_items.append(
                    (str(len(map_data["stations"])),
                     "Obs. Stations"))
            if info_items:
                H.append('<div class="map-info-grid">')
                for val, lbl in info_items:
                    H.append(
                        f'<div class="map-info-item">'
                        f'<div class="val">{val}</div>'
                        f'<div class="lbl">{lbl}</div></div>')
                H.append('</div>')

            # Leaflet initialization script
            center_lat = map_data.get("center_lat", 45.0)
            center_lon = map_data.get("center_lon", -66.0)
            bounds_js = json.dumps(map_data.get("bounds",
                                                [[44, -68], [46, -65]]))
            hru_geojson = map_data.get("hru_geojson", "null")
            net_geojson = map_data.get("network_geojson", "null")
            stations_js = json.dumps(map_data.get("stations", []))

            # HRU color palette based on area
            H.append(f"""<script>
(function(){{
  var map = L.map('basinMap', {{zoomControl:true}});

  /* ── Tile layers ────────────────────────── */
  var osmLight = L.tileLayer(
    'https://{{s}}.basemaps.cartocdn.com/light_all/{{z}}/{{x}}/{{y}}@2x.png',
    {{attribution:'&copy; <a href="https://carto.com/">CARTO</a> &copy; <a href="https://osm.org/copyright">OSM</a>',maxZoom:18}});
  var osmTopo = L.tileLayer(
    'https://{{s}}.tile.opentopomap.org/{{z}}/{{x}}/{{y}}.png',
    {{attribution:'&copy; OpenTopoMap',maxZoom:17}});
  var satellite = L.tileLayer(
    'https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{{z}}/{{y}}/{{x}}',
    {{attribution:'&copy; Esri',maxZoom:18}});
  osmLight.addTo(map);

  /* ── Fit bounds ─────────────────────────── */
  var bnds = {bounds_js};
  map.fitBounds(bnds, {{padding:[30,30]}});

  /* ── Stream order → color ───────────────── */
  var ordColors = {{1:'#93c5fd',2:'#60a5fa',3:'#3b82f6',4:'#2563eb',5:'#1d4ed8',6:'#1e40af',7:'#1e3a8a'}};
  var ordWidths = {{1:1.5,2:2,3:2.8,4:3.5,5:4.5,6:5.5,7:7}};

  /* ── HRU layer ──────────────────────────── */
  var hruData = {hru_geojson};
  var hruLayer = null;
  if (hruData) {{
    /* color HRUs by area quintiles (indigo palette) */
    var areas = hruData.features.map(function(f){{
      return f.properties.unitarea || f.properties.area || 0;
    }}).filter(function(a){{return a>0}});
    areas.sort(function(a,b){{return a-b}});
    var q = [0,.2,.4,.6,.8,1].map(function(p){{
      return areas[Math.min(Math.floor(p*(areas.length-1)),areas.length-1)] || 0;
    }});
    var hruPal = ['#e0e7ff','#c7d2fe','#a5b4fc','#818cf8','#6366f1'];
    function hruColor(a) {{
      for (var i=0;i<5;i++) if (a<=q[i+1]) return hruPal[i];
      return hruPal[4];
    }}
    hruLayer = L.geoJSON(hruData, {{
      style: function(f) {{
        var a = f.properties.unitarea || f.properties.area || 0;
        return {{
          fillColor: hruColor(a), fillOpacity: 0.35,
          color: '#6366f1', weight: 1.2, opacity: 0.6
        }};
      }},
      onEachFeature: function(f, layer) {{
        var p = f.properties;
        var pop = '<b>HRU / GRU</b><br>';
        if (p.GRU_ID) pop += 'ID: <b>' + p.GRU_ID + '</b><br>';
        if (p.unitarea) pop += 'Area: ' + p.unitarea.toFixed(2) + ' km&sup2;<br>';
        else if (p.area) pop += 'Area: ' + (p.area/1e6).toFixed(2) + ' km&sup2;<br>';
        if (p.avg_subbas) pop += 'Avg subbasin: ' + p.avg_subbas.toFixed(2) + '<br>';
        layer.bindPopup(pop);
        layer.on('mouseover', function(){{ this.setStyle({{fillOpacity:0.6,weight:2.5}}); }});
        layer.on('mouseout', function(){{ hruLayer.resetStyle(this); }});
      }}
    }}).addTo(map);
  }}

  /* ── River network layer ────────────────── */
  var netData = {net_geojson};
  var netLayer = null;
  if (netData) {{
    netLayer = L.geoJSON(netData, {{
      style: function(f) {{
        var ord = f.properties.order || 1;
        return {{
          color: ordColors[Math.min(ord,7)] || '#3b82f6',
          weight: ordWidths[Math.min(ord,7)] || 2,
          opacity: 0.85, lineCap: 'round', lineJoin: 'round'
        }};
      }},
      onEachFeature: function(f, layer) {{
        var p = f.properties;
        var pop = '<b>River Reach</b><br>';
        if (p.LINKNO) pop += 'ID: <b>' + p.LINKNO + '</b><br>';
        if (p.order) pop += 'Strahler order: ' + p.order + '<br>';
        if (p.Length) pop += 'Length: ' + (p.Length/1000).toFixed(2) + ' km<br>';
        if (p.Slope != null) pop += 'Slope: ' + p.Slope.toFixed(6) + '<br>';
        if (p.uparea) pop += 'Upstream area: ' + p.uparea.toFixed(1) + ' km&sup2;<br>';
        if (p.sinuosity) pop += 'Sinuosity: ' + p.sinuosity.toFixed(2) + '<br>';
        if (p.DSLINKNO && p.DSLINKNO > 0)
          pop += 'Downstream: ' + p.DSLINKNO + '<br>';
        layer.bindPopup(pop);
      }}
    }}).addTo(map);
  }}

  /* ── Observation station markers ────────── */
  var stations = {stations_js};
  var stationLayer = L.layerGroup();
  stations.forEach(function(s) {{
    var icon = L.divIcon({{
      className:'',
      html:'<div style="width:14px;height:14px;background:#ef4444;border:2.5px solid #fff;border-radius:50%;box-shadow:0 0 6px rgba(0,0,0,.35)"></div>',
      iconSize:[14,14], iconAnchor:[7,7]
    }});
    var pop = '<b>Observation Station</b><br>'
      + 'Reach ID: <b>' + s.reach_id + '</b><br>'
      + 'Species: ' + s.species.join(', ') + '<br>'
      + 'Observations: ' + s.n_obs;
    L.marker([s.lat,s.lon],{{icon:icon}}).bindPopup(pop).addTo(stationLayer);
  }});
  stationLayer.addTo(map);

  /* ── Layer control ──────────────────────── */
  var baseMaps = {{"Light":osmLight,"Topo":osmTopo,"Satellite":satellite}};
  var overlays = {{}};
  if (hruLayer) overlays["HRU Boundaries"] = hruLayer;
  if (netLayer) overlays["River Network"] = netLayer;
  if (stations.length) overlays["Obs. Stations"] = stationLayer;
  L.control.layers(baseMaps, overlays, {{collapsed:false,position:'topright'}}).addTo(map);

  /* ── Legend ─────────────────────────────── */
  var legend = L.control({{position:'bottomright'}});
  legend.onAdd = function() {{
    var d = L.DomUtil.create('div','legend-box');
    d.innerHTML =
      '<div style="font-weight:600;margin-bottom:4px">Legend</div>'
      + '<i style="background:#a5b4fc"></i> HRU (area shading)<br>'
      + '<i style="background:#93c5fd;width:20px;height:3px;border-radius:2px"></i> Order 1<br>'
      + '<i style="background:#3b82f6;width:20px;height:4px;border-radius:2px"></i> Order 3<br>'
      + '<i style="background:#1e40af;width:20px;height:5px;border-radius:2px"></i> Order 5+<br>'
      + '<i style="background:#ef4444;border-radius:50%"></i> Obs. Station';
    return d;
  }};
  legend.addTo(map);

  /* ── Scale bar ──────────────────────────── */
  L.control.scale({{imperial:false}}).addTo(map);

  /* ── Re-apply tile filter on theme change ─ */
  window._owqMapInvalidate = function() {{
    setTimeout(function(){{ map.invalidateSize(); }}, 100);
  }};
}})();
</script>""")
            H.append('</div></div>')

        # ── Configuration ─────────────────────────────────────────────
        H.append('<div class="section" id="config"><h2>Configuration</h2><div class="card"><div class="config-grid">')
        for key, val in config_summary.items():
            if val:
                dv = ", ".join(f"{k}: {v}" for k, v in val.items()) if isinstance(val, dict) else str(val)
                H.append(f'<div class="config-item"><span class="key">{key}</span><span class="val">{dv}</span></div>')
        H.append('</div></div></div>')

        # ── Best Parameters ───────────────────────────────────────────
        H.append('<div class="section" id="parameters"><h2>Best Parameters</h2><div class="card"><div class="table-wrap"><table>')
        H.append('<tr><th>Parameter</th><th>Best Value</th><th>Initial</th>'
                 '<th>Lower</th><th>Upper</th><th>Transform</th>'
                 '<th>Type</th><th>Position</th></tr>')
        if summary and summary.best_parameters:
            for pdef in self.parameter_definitions:
                name = pdef["name"]
                val = summary.best_parameters.get(name, float('nan'))
                init_v = pdef.get("initial", float('nan'))
                lb, ub = pdef.get("bounds", (0, 1))
                transform = pdef.get("transform", "linear")
                file_type = pdef.get("file_type", "&mdash;")
                if ub > lb:
                    pct = max(0.0, min(100.0, (val - lb) / (ub - lb) * 100))
                    bar = (f'<div class="range-bar"><div class="range-track">'
                           f'<div class="range-fill" style="width:{pct:.0f}%"></div></div>'
                           f'<span class="range-pct">{pct:.0f}%</span></div>')
                else:
                    bar = "&mdash;"
                H.append(f'<tr><td><strong>{name}</strong></td><td class="num">{val:.6g}</td>'
                         f'<td class="num">{init_v:.4g}</td>'
                         f'<td class="num">{lb:.4g}</td><td class="num">{ub:.4g}</td>'
                         f'<td>{transform}</td><td>{file_type}</td><td>{bar}</td></tr>')
        H.append('</table></div></div></div>')

        # ==============================================================
        # INTERACTIVE PLOTLY CHARTS  (theme-aware)
        # ==============================================================
        # Inject the Plotly helper functions BEFORE any chart script runs
        H.append("""<script>
/* ── Plotly helpers (must load before chart scripts) ── */
window._owqPlotIds = [];
function _owqLayout(extra) {
  var t = document.documentElement.getAttribute('data-theme') === 'dark';
  var gc = t ? '#2d2e42' : '#f0f0f0';
  var base = {
    margin:{l:60,r:25,t:25,b:45},
    plot_bgcolor: t ? '#1e1f2e' : '#fff',
    paper_bgcolor: t ? '#1e1f2e' : '#fff',
    font:{color: t ? '#ccc' : '#333', family:'Inter,sans-serif', size:12},
    colorway:['#6366f1','#ef4444','#10b981','#f59e0b','#8b5cf6','#ec4899','#06b6d4']
  };
  var merged = Object.assign({}, base, extra || {});
  for (var k in merged) {
    if (k.startsWith('xaxis') || k.startsWith('yaxis')) {
      merged[k] = merged[k] || {};
      if (!merged[k].gridcolor) merged[k].gridcolor = gc;
      if (!merged[k].zerolinecolor) merged[k].zerolinecolor = gc;
      if (typeof merged[k].title === 'string')
        merged[k].title = {text: merged[k].title, font:{size:11}};
    }
  }
  return merged;
}
function _owqRegister(id) { window._owqPlotIds.push(id); }
</script>""")

        # ── 1. Convergence ────────────────────────────────────────────
        if has_history:
            all_objs = [e.get("objective", None) for e in self.history]
            evals_arr, best_arr = self.get_convergence_data()

            H.append('<div class="section" id="convergence"><h2>Convergence</h2><div class="card">')
            H.append('<div id="plotConvergence" style="width:100%;height:440px"></div>')
            H.append('<p class="plot-caption">Best objective vs evaluation &mdash; hover for details, drag to zoom</p>')
            H.append(f"""<script>
(function(){{
  var allObj={_js(all_objs)};
  var bestX={_js(evals_arr.tolist())};
  var bestY={_js(best_arr.tolist())};
  var evalNums=allObj.map(function(_,i){{return i}});
  var traces=[
    {{x:evalNums,y:allObj,mode:'markers',name:'All evaluations',
      marker:{{color:'rgba(150,150,170,0.5)',size:7}},
      hovertemplate:'Eval %{{x}}<br>Obj: %{{y:.4g}}<extra></extra>'}},
    {{x:bestX,y:bestY,mode:'lines+markers',name:'Best so far',
      line:{{color:'var(--accent)',width:2.5,shape:'hv'}},marker:{{size:5}},
      hovertemplate:'Eval %{{x}}<br>Best: %{{y:.4g}}<extra></extra>'}}
  ];
  if(bestX.length>0)
    traces.push({{x:[bestX[bestX.length-1]],y:[bestY[bestY.length-1]],
      mode:'markers',name:'Final best',
      marker:{{color:'#ef4444',size:14,symbol:'star'}},
      hovertemplate:'Best: %{{y:.4g}}<extra></extra>'}});
  Plotly.newPlot('plotConvergence',traces,_owqLayout({{
    xaxis:{{title:'Evaluation Number'}},
    yaxis:{{title:'Objective Function Value'}},
    legend:{{x:.98,xanchor:'right',y:.98}},hovermode:'closest'
  }}),{PCFG});
  _owqRegister('plotConvergence');
}})();
</script>""")
            H.append('</div></div>')

        # ── 2. Parameter Evolution ────────────────────────────────────
        if has_history:
            pnames = list(self.history[0].get("parameters", {}).keys()) if self.history else []
            n_p = len(pnames)
            if n_p > 0:
                n_cols_p = min(3, n_p)
                n_rows_p = int(np.ceil(n_p / n_cols_p))
                plot_h = 300 * n_rows_p
                param_data_js = {}
                for pn in pnames:
                    ev, vl = self.get_parameter_evolution(pn)
                    param_data_js[pn] = {"x": ev.tolist(), "y": vl.tolist()}
                bounds_map = {}
                for pdef in self.parameter_definitions:
                    bounds_map[pdef["name"]] = list(pdef.get("bounds", [None, None]))

                H.append('<div class="section" id="param-evolution"><h2>Parameter Evolution</h2><div class="card">')
                H.append(f'<div id="plotParamEvol" style="width:100%;height:{plot_h}px"></div>')
                H.append('<p class="plot-caption">Sampled values across evaluations &mdash; red dashed = bounds</p>')
                H.append(f"""<script>
(function(){{
  var pd={_js(param_data_js)},b={_js(bounds_map)},nm={_js(pnames)};
  var nC={n_cols_p},nR={n_rows_p},traces=[];
  nm.forEach(function(pn,i){{
    var xr=i===0?'x':'x'+(i+1),yr=i===0?'y':'y'+(i+1),d=pd[pn];
    traces.push({{x:d.x,y:d.y,mode:'markers',xaxis:xr,yaxis:yr,
      marker:{{size:6,opacity:.7}},showlegend:false,
      hovertemplate:'Eval %{{x}}<br>'+pn+': %{{y:.4g}}<extra></extra>'}});
    var bb=b[pn]||[null,null];
    [0,1].forEach(function(k){{
      if(bb[k]!==null) traces.push({{x:[d.x[0]||0,d.x[d.x.length-1]||0],y:[bb[k],bb[k]],
        mode:'lines',line:{{color:'#ef4444',width:1.5,dash:'dash'}},
        xaxis:xr,yaxis:yr,showlegend:false,hoverinfo:'skip'}});
    }});
  }});
  var lo=_owqLayout({{grid:{{rows:nR,columns:nC,pattern:'independent',ygap:.2,xgap:.08}},showlegend:false}});
  for(var i=0;i<nm.length;i++){{
    var xs=i===0?'xaxis':'xaxis'+(i+1),ys=i===0?'yaxis':'yaxis'+(i+1);
    lo[xs]=lo[xs]||{{}};lo[xs].title='Evaluation';
    lo[ys]=lo[ys]||{{}};lo[ys].title=nm[i];
  }}
  Plotly.newPlot('plotParamEvol',traces,lo,{PCFG});_owqRegister('plotParamEvol');
}})();
</script>""")
                H.append('</div></div>')

        # ── 3. Parameter Correlations ─────────────────────────────────
        if has_correlations:
            apd = {}
            for entry in self.history:
                for name, value in entry.get("parameters", {}).items():
                    apd.setdefault(name, []).append(value)
            df_corr = pd.DataFrame(apd)
            corr = df_corr.corr()
            cn = list(corr.columns)
            cv = corr.values.tolist()
            ct = [[f"{cv[i][j]:.2f}" for j in range(len(cn))] for i in range(len(cn))]

            H.append('<div class="section" id="correlations"><h2>Parameter Correlations</h2><div class="card">')
            H.append('<div id="plotCorr" style="width:100%;height:480px"></div>')
            H.append('<p class="plot-caption">Pearson correlation between sampled parameter values</p>')
            H.append(f"""<script>
(function(){{
  var nm={_js(cn)},z={_js(cv)},tx={_js(ct)};
  Plotly.newPlot('plotCorr',[{{z:z,x:nm,y:nm,type:'heatmap',
    colorscale:'RdBu',reversescale:true,zmin:-1,zmax:1,
    text:tx,texttemplate:'%{{text}}',textfont:{{size:11}},
    hovertemplate:'%{{x}} vs %{{y}}<br>r = %{{z:.3f}}<extra></extra>',
    colorbar:{{title:'r',titleside:'right'}}}}],
    _owqLayout({{xaxis:{{tickangle:-45}},yaxis:{{autorange:'reversed'}}}}),{PCFG});
  _owqRegister('plotCorr');
}})();
</script>""")
            H.append('</div></div>')

        # ── 4. Sensitivity Analysis ───────────────────────────────────
        if has_sensitivity:
            sens = self.sensitivity_results
            sm_type = sens.get("method", "morris")
            sn = sens.get("parameter_names", [])
            H.append('<div class="section" id="sensitivity"><h2>Sensitivity Analysis</h2><div class="card">')
            H.append('<div id="plotSens" style="width:100%;height:440px"></div>')
            if sm_type == "morris" and sn:
                mu_s = sens.get("mu_star", [])
                sig = sens.get("sigma", [])
                si = list(np.argsort(mu_s)[::-1])
                H.append(f"""<script>
(function(){{
  var nm={_js([sn[i] for i in si])},mu={_js([mu_s[i] for i in si])},sg={_js([sig[i] for i in si])};
  Plotly.newPlot('plotSens',[
    {{y:nm,x:mu,type:'bar',orientation:'h',name:'\\u03bc*',marker:{{color:'var(--accent)'}},
      hovertemplate:'%{{y}}<br>\\u03bc* = %{{x:.4g}}<extra></extra>'}},
    {{y:nm,x:sg,type:'bar',orientation:'h',name:'\\u03c3',marker:{{color:'var(--amber)'}},
      hovertemplate:'%{{y}}<br>\\u03c3 = %{{x:.4g}}<extra></extra>'}}
  ],_owqLayout({{barmode:'group',yaxis:{{autorange:'reversed'}}}}),{PCFG});
  _owqRegister('plotSens');
}})();
</script>""")
            elif sm_type == "sobol" and sn:
                S1 = sens.get("S1", [])
                ST = sens.get("ST", [])
                si = list(np.argsort(ST)[::-1])
                H.append(f"""<script>
(function(){{
  var nm={_js([sn[i] for i in si])},s1={_js([S1[i] for i in si])},st={_js([ST[i] for i in si])};
  Plotly.newPlot('plotSens',[
    {{x:nm,y:s1,type:'bar',name:'S1',marker:{{color:'var(--accent)'}}}},
    {{x:nm,y:st,type:'bar',name:'ST',marker:{{color:'#ef4444'}}}}
  ],_owqLayout({{barmode:'group',xaxis:{{tickangle:-45}}}}),{PCFG});
  _owqRegister('plotSens');
}})();
</script>""")
            H.append('<p class="plot-caption">Parameter sensitivity ranking</p>')
            H.append('</div></div>')

        # ── Species Performance Metrics ───────────────────────────────
        if species_metrics:
            H.append('<div class="section" id="performance"><h2>Model Performance by Species</h2><div class="card"><div class="table-wrap"><table>')
            H.append('<tr><th>Species</th><th>N</th><th>KGE</th><th>NSE</th><th>RMSE</th>'
                     '<th>PBIAS%</th><th>R&sup2;</th><th>Obs&nbsp;Mean</th><th>Sim&nbsp;Mean</th><th>Rating</th></tr>')
            for sm in species_metrics:
                k = sm["KGE"]
                b = ('<span class="badge badge-good">Good</span>' if k > .75
                     else '<span class="badge badge-fair">Fair</span>' if k > .5
                     else '<span class="badge badge-poor">Poor</span>')
                H.append(f'<tr><td><strong>{sm["species"]}</strong></td>'
                         f'<td class="num">{sm["n"]}</td><td class="num">{sm["KGE"]:.3f}</td>'
                         f'<td class="num">{sm["NSE"]:.3f}</td><td class="num">{sm["RMSE"]:.3f}</td>'
                         f'<td class="num">{sm["PBIAS"]:.1f}</td><td class="num">{sm["R2"]:.3f}</td>'
                         f'<td class="num">{sm["obs_mean"]:.3g}</td><td class="num">{sm["sim_mean"]:.3g}</td>'
                         f'<td>{b}</td></tr>')
            H.append('</table></div></div></div>')

        # ── 5. Time-series comparison ─────────────────────────────────
        if has_matched:
            data = matched_data.copy()
            data['datetime'] = pd.to_datetime(data['datetime'])
            unique_species = sorted(data['species'].unique())
            n_sp = len(unique_species)
            ts_data = {}
            for sp in unique_species:
                sp_df = data[data['species'] == sp].sort_values('datetime')
                ts_data[sp] = {
                    "dt": sp_df['datetime'].dt.strftime('%Y-%m-%d %H:%M').tolist(),
                    "obs": [None if np.isnan(v) else float(v) for v in sp_df['observed'].values],
                    "sim": [None if np.isnan(v) else float(v) for v in sp_df['simulated'].values],
                }
            H.append(f'<div class="section" id="timeseries"><h2>Time Series Comparison</h2><div class="card">')
            H.append(f'<div id="plotTS" style="width:100%;height:{max(400,340*n_sp)}px"></div>')
            H.append('<p class="plot-caption">Observed vs simulated for best parameter set</p>')
            H.append(f"""<script>
(function(){{
  var td={_js(ts_data)},sp={_js(unique_species)},nS=sp.length,tr=[];
  sp.forEach(function(s,i){{
    var d=td[s],xr=i===0?'x':'x'+(i+1),yr=i===0?'y':'y'+(i+1);
    tr.push({{x:d.dt,y:d.obs,mode:'markers',name:'Obs ('+s+')',
      xaxis:xr,yaxis:yr,marker:{{size:7,symbol:'circle'}},
      hovertemplate:'%{{x}}<br>Obs: %{{y:.4g}}<extra>'+s+'</extra>',legendgroup:s}});
    tr.push({{x:d.dt,y:d.sim,mode:'markers+lines',name:'Sim ('+s+')',
      xaxis:xr,yaxis:yr,marker:{{size:5,symbol:'x'}},
      line:{{width:1,dash:'dot'}},
      hovertemplate:'%{{x}}<br>Sim: %{{y:.4g}}<extra>'+s+'</extra>',legendgroup:s}});
  }});
  var lo=_owqLayout({{grid:{{rows:nS,columns:1,pattern:'independent',ygap:.14}},hovermode:'closest'}});
  for(var i=0;i<nS;i++){{
    var xs=i===0?'xaxis':'xaxis'+(i+1),ys=i===0?'yaxis':'yaxis'+(i+1);
    lo[xs]={{title:'Date',type:'date'}};lo[ys]={{title:sp[i]+' (mg/L)'}};
  }}
  Plotly.newPlot('plotTS',tr,lo,{PCFG});_owqRegister('plotTS');
}})();
</script>""")
            H.append('</div></div>')

        # ── 6. Scatter obs vs sim ─────────────────────────────────────
        if has_matched:
            n_cs = min(3, n_sp)
            n_rs = int(np.ceil(n_sp / n_cs))
            sc_data = {}
            for sp in unique_species:
                sp_df = data[data['species'] == sp]
                o, s = sp_df['observed'].values, sp_df['simulated'].values
                v = ~(np.isnan(o) | np.isnan(s))
                sc_data[sp] = {"obs": [float(x) for x in o[v]], "sim": [float(x) for x in s[v]]}
            H.append('<div class="section" id="scatter"><h2>Observed vs Simulated</h2><div class="card">')
            H.append(f'<div id="plotSc" style="width:100%;height:{max(400,370*n_rs)}px"></div>')
            H.append('<p class="plot-caption">Dashed = 1:1 line, solid red = linear fit</p>')
            H.append(f"""<script>
(function(){{
  var sd={_js(sc_data)},sp={_js(unique_species)},nC={n_cs},nR={n_rs},tr=[];
  sp.forEach(function(s,i){{
    var d=sd[s],xr=i===0?'x':'x'+(i+1),yr=i===0?'y':'y'+(i+1);
    tr.push({{x:d.obs,y:d.sim,mode:'markers',xaxis:xr,yaxis:yr,
      marker:{{size:7,opacity:.7}},showlegend:false,
      hovertemplate:'Obs: %{{x:.4g}}<br>Sim: %{{y:.4g}}<extra>'+s+'</extra>'}});
    if(d.obs.length>0){{
      var all=d.obs.concat(d.sim),mn=Math.min.apply(null,all),mx=Math.max.apply(null,all);
      tr.push({{x:[mn,mx],y:[mn,mx],mode:'lines',line:{{color:'var(--text3)',width:1.5,dash:'dash'}},
        xaxis:xr,yaxis:yr,showlegend:false,hoverinfo:'skip'}});
      if(d.obs.length>=3){{
        var n=d.obs.length,sx=0,sy=0,sxy=0,sxx=0;
        for(var j=0;j<n;j++){{sx+=d.obs[j];sy+=d.sim[j];sxy+=d.obs[j]*d.sim[j];sxx+=d.obs[j]*d.obs[j]}}
        var sl=(n*sxy-sx*sy)/(n*sxx-sx*sx),ic=(sy-sl*sx)/n;
        tr.push({{x:[mn,mx],y:[sl*mn+ic,sl*mx+ic],mode:'lines',line:{{color:'#ef4444',width:1.5}},
          xaxis:xr,yaxis:yr,showlegend:false,
          hovertemplate:'y='+sl.toFixed(3)+'x+'+ic.toFixed(3)+'<extra></extra>'}});
      }}
    }}
  }});
  var lo=_owqLayout({{grid:{{rows:nR,columns:nC,pattern:'independent',ygap:.15,xgap:.08}}}});
  for(var i=0;i<sp.length;i++){{
    var xs=i===0?'xaxis':'xaxis'+(i+1),ys=i===0?'yaxis':'yaxis'+(i+1);
    lo[xs]={{title:'Observed'}};lo[ys]={{title:'Simulated',scaleanchor:xs,scaleratio:1}};
  }}
  Plotly.newPlot('plotSc',tr,lo,{PCFG});_owqRegister('plotSc');
}})();
</script>""")
            H.append('</div></div>')

        # ── 7. Residuals ──────────────────────────────────────────────
        if has_matched:
            rd = {}
            for sp in unique_species:
                sp_df = data[data['species'] == sp].dropna(subset=['observed', 'simulated']).copy()
                sp_df['residual'] = sp_df['simulated'] - sp_df['observed']
                rd[sp] = {
                    "dt": sp_df['datetime'].dt.strftime('%Y-%m-%d %H:%M').tolist(),
                    "r": [float(v) for v in sp_df['residual'].values],
                }
            H.append('<div class="section" id="residuals"><h2>Residual Analysis</h2><div class="card">')
            H.append(f'<div id="plotRes" style="width:100%;height:{max(400,340*n_sp)}px"></div>')
            H.append('<p class="plot-caption">Residuals (sim &minus; obs) over time and distribution</p>')
            H.append(f"""<script>
(function(){{
  var rd={_js(rd)},sp={_js(unique_species)},nS=sp.length,tr=[];
  sp.forEach(function(s,i){{
    var d=rd[s],xr1='x'+(i*2+1===1?'':(i*2+1)),yr1='y'+(i*2+1===1?'':(i*2+1));
    var xr2='x'+(i*2+2),yr2='y'+(i*2+2);
    tr.push({{x:d.dt,y:d.r,mode:'markers',xaxis:xr1,yaxis:yr1,
      marker:{{size:6,opacity:.7}},showlegend:false,
      hovertemplate:'%{{x}}<br>Res: %{{y:.4g}}<extra>'+s+'</extra>'}});
    if(d.dt.length>0) tr.push({{x:[d.dt[0],d.dt[d.dt.length-1]],y:[0,0],mode:'lines',
      line:{{color:'#ef4444',width:1.5,dash:'dash'}},xaxis:xr1,yaxis:yr1,showlegend:false,hoverinfo:'skip'}});
    tr.push({{x:d.r,type:'histogram',xaxis:xr2,yaxis:yr2,marker:{{opacity:.7}},showlegend:false,
      hovertemplate:'Bin: %{{x:.3g}}<br>Count: %{{y}}<extra></extra>'}});
  }});
  var lo=_owqLayout({{grid:{{rows:nS,columns:2,pattern:'independent',ygap:.15,xgap:.08}}}});
  for(var i=0;i<nS;i++){{
    var x1='xaxis'+(i*2+1===1?'':(i*2+1)),y1='yaxis'+(i*2+1===1?'':(i*2+1));
    var x2='xaxis'+(i*2+2),y2='yaxis'+(i*2+2);
    lo[x1]={{title:'Date',type:'date'}};lo[y1]={{title:sp[i]+' Residual'}};
    lo[x2]={{title:'Residual'}};lo[y2]={{title:'Count'}};
  }}
  Plotly.newPlot('plotRes',tr,lo,{PCFG});_owqRegister('plotRes');
}})();
</script>""")
            H.append('</div></div>')

        # ── Objective Statistics ───────────────────────────────────────
        if obj_stats:
            sl = {"best": "Best", "worst": "Worst", "mean": "Mean",
                  "std": "Std Dev", "median": "Median",
                  "q25": "25th Pctl", "q75": "75th Pctl"}
            H.append('<div class="section" id="obj-stats"><h2>Objective Function Statistics</h2><div class="card"><div class="table-wrap"><table>')
            H.append('<tr><th>Statistic</th><th>Value</th></tr>')
            for k, v in obj_stats.items():
                H.append(f'<tr><td>{sl.get(k, k)}</td><td class="num">{v:.6g}</td></tr>')
            H.append('</table></div></div></div>')

        # ── Parameter Statistics ───────────────────────────────────────
        if not param_stats.empty:
            H.append('<div class="section" id="param-stats"><h2>Parameter Statistics</h2><div class="card"><div class="table-wrap"><table>')
            H.append('<tr><th>Parameter</th><th>Min</th><th>Max</th><th>Mean</th><th>Std</th><th>Median</th></tr>')
            for _, r in param_stats.iterrows():
                H.append(f'<tr><td><strong>{r["parameter"]}</strong></td>'
                         f'<td class="num">{r["min"]:.6g}</td><td class="num">{r["max"]:.6g}</td>'
                         f'<td class="num">{r["mean"]:.6g}</td><td class="num">{r["std"]:.6g}</td>'
                         f'<td class="num">{r["median"]:.6g}</td></tr>')
            H.append('</table></div></div></div>')

        # ── Evaluation History (collapsible) ──────────────────────────
        if has_history:
            pnh = list(self.history[0].get("parameters", {}).keys())
            bov = min((e.get("objective", float('inf')) for e in self.history), default=float('inf'))

            H.append('<div class="section" id="eval-history"><h2>Evaluation History</h2><div class="card">')
            H.append(f'<div class="collapsible-header" id="histToggle" onclick="toggleHist()">'
                     f'<span>{len(self.history)} evaluations</span>'
                     f'<svg class="chevron" width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M6 9l6 6 6-6"/></svg></div>')
            H.append('<div class="collapsible-body" id="histBody"><div class="table-wrap"><table>')
            hdr = '<tr><th>#</th><th>Objective</th>'
            for pn in pnh:
                hdr += f'<th>{pn}</th>'
            H.append(hdr + '</tr>')
            dh = self.history[:200]
            for entry in dh:
                eid = entry.get("eval_id", 0)
                obj = entry.get("objective", float('nan'))
                ib = (abs(obj - bov) < 1e-12) if not (isinstance(obj, float) and np.isnan(obj)) else False
                st = ' style="background:rgba(16,185,129,.1);font-weight:600"' if ib else ''
                row = f'<tr{st}><td class="num">{eid}</td><td class="num">{obj:.6g}</td>'
                params = entry.get("parameters", {})
                for pn in pnh:
                    v = params.get(pn, float('nan'))
                    row += f'<td class="num">{v:.6g}</td>'
                H.append(row + '</tr>')
            if len(self.history) > 200:
                H.append(f'<tr><td colspan="{2+len(pnh)}" style="text-align:center;color:var(--text3)">'
                         f'&hellip; {len(self.history)-200} more (see calibration_history.json)</td></tr>')
            H.append('</table></div></div></div></div>')

        # ── Footer ────────────────────────────────────────────────────
        H.append(f"""
</div><!-- /container -->
</div><!-- /main -->
</div><!-- /layout -->

<div class="footer">
  Generated by <strong>OpenWQ Calibration Framework</strong> &mdash; {now}<br>
  <a href="https://github.com/DiogoCostaPT/openwq">github.com/DiogoCostaPT/openwq</a>
</div>
""")

        # ==============================================================
        # JavaScript: theme engine, scroll-spy, animations, Plotly sync
        # ==============================================================
        H.append("""
<script>
/* ── Theme relayout for Plotly charts ────────────────── */
function _owqRelayoutAll() {
  var t = document.documentElement.getAttribute('data-theme') === 'dark';
  var upd = {
    'plot_bgcolor': t ? '#1e1f2e' : '#fff',
    'paper_bgcolor': t ? '#1e1f2e' : '#fff',
    'font.color': t ? '#ccc' : '#333'
  };
  window._owqPlotIds.forEach(function(id) {
    var el = document.getElementById(id);
    if (el && el.data) {
      Plotly.relayout(id, upd);
      // Update grid colors on all axes
      var lo = el.layout || {};
      for (var k in lo) {
        if (k.startsWith('xaxis') || k.startsWith('yaxis')) {
          var gc = t ? '#2d2e42' : '#f0f0f0';
          var u = {}; u[k+'.gridcolor'] = gc; u[k+'.zerolinecolor'] = gc;
          Plotly.relayout(id, u);
        }
      }
    }
  });
}

/* ── Dark / Light toggle ────────────────────────────── */
(function() {
  var stored = localStorage.getItem('owq-theme');
  if (stored) document.documentElement.setAttribute('data-theme', stored);
  var btn = document.getElementById('themeToggle');
  var icon = document.getElementById('themeIcon');
  var label = document.getElementById('themeLabel');
  function update() {
    var dark = document.documentElement.getAttribute('data-theme') === 'dark';
    label.textContent = dark ? 'Light mode' : 'Dark mode';
    icon.innerHTML = dark
      ? '<circle cx="12" cy="12" r="5"/><path d="M12 1v2M12 21v2M4.22 4.22l1.42 1.42M18.36 18.36l1.42 1.42M1 12h2M21 12h2M4.22 19.78l1.42-1.42M18.36 5.64l1.42-1.42"/>'
      : '<path d="M21 12.79A9 9 0 1 1 11.21 3 7 7 0 0 0 21 12.79z"/>';
  }
  update();
  btn.addEventListener('click', function() {
    var cur = document.documentElement.getAttribute('data-theme');
    var next = cur === 'dark' ? 'light' : 'dark';
    document.documentElement.setAttribute('data-theme', next);
    localStorage.setItem('owq-theme', next);
    update();
    setTimeout(_owqRelayoutAll, 50);
    if (window._owqMapInvalidate) window._owqMapInvalidate();
  });
})();

/* ── Scroll-spy (sidebar active highlighting) ─────── */
(function() {
  var sections = document.querySelectorAll('.section[id]');
  var links = document.querySelectorAll('#sidebarNav a');
  if (!sections.length || !links.length) return;
  var observer = new IntersectionObserver(function(entries) {
    entries.forEach(function(entry) {
      if (entry.isIntersecting) {
        links.forEach(function(a) { a.classList.remove('active'); });
        var link = document.querySelector('#sidebarNav a[data-section="'+entry.target.id+'"]');
        if (link) link.classList.add('active');
      }
    });
  }, {rootMargin:'-20% 0px -70% 0px'});
  sections.forEach(function(s) { observer.observe(s); });
})();

/* ── Fade-in on scroll ────────────────────────────── */
(function() {
  var secs = document.querySelectorAll('.section');
  var obs = new IntersectionObserver(function(entries) {
    entries.forEach(function(e) {
      if (e.isIntersecting) { e.target.classList.add('visible'); obs.unobserve(e.target); }
    });
  }, {threshold:0.08});
  secs.forEach(function(s) { obs.observe(s); });
})();

/* ── Collapsible history ─────────────────────────── */
function toggleHist() {
  var h = document.getElementById('histToggle');
  var b = document.getElementById('histBody');
  h.classList.toggle('open');
  b.classList.toggle('open');
}
</script>
</body>
</html>""")

        # ── Write HTML ────────────────────────────────────────────────
        report_path = output_dir / "calibration_report.html"
        with open(report_path, 'w', encoding='utf-8') as f:
            f.write("\n".join(H))

        logger.info(f"HTML calibration report saved to {report_path}")
        return report_path


# ══════════════════════════════════════════════════════════════════════
# Per-Basin Report  — aggregates results across multiple variants
# ══════════════════════════════════════════════════════════════════════

def generate_basin_report(
    basin_id: str,
    variant_results: Dict[str, Dict],
    output_dir: str,
    shapefile_dir: str = None,
    observation_data_path: str = None,
) -> Path:
    """
    Generate a per-basin HTML report that compares calibration results
    across all variants (A/B/C/D…).

    Parameters
    ----------
    basin_id : str
        Basin identifier (e.g. "CAN_01AM001_meso")
    variant_results : Dict[str, Dict]
        Keys are variant labels ("A", "B", …). Each value is a dict with:
          - results_dir : str — path to that variant's results/ dir
          - calibration_parameters : List[Dict] — parameter definitions
          - matched_data_csv : str (optional) — path to matched data CSV
          - config : Dict (optional) — calibration config summary
    output_dir : str
        Where to write the basin report HTML
    shapefile_dir : str, optional
        Path to basin shapefile directory
    observation_data_path : str, optional
        Path to observation CSV

    Returns
    -------
    Path   to the generated HTML file
    """
    from datetime import datetime
    import math

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    def _js(obj):
        """Serialise *obj* to a JSON string safe for JS embedding."""
        def _conv(o):
            if isinstance(o, (np.integer,)):
                return int(o)
            if isinstance(o, (np.floating,)):
                v = float(o)
                return None if math.isnan(v) or math.isinf(v) else v
            if isinstance(o, np.ndarray):
                return [_conv(x) for x in o.tolist()]
            if isinstance(o, float):
                return None if math.isnan(o) or math.isinf(o) else o
            return o
        if isinstance(obj, (list, tuple)):
            obj = [_conv(x) for x in obj]
        elif isinstance(obj, dict):
            obj = {k: _conv(v) for k, v in obj.items()}
        else:
            obj = _conv(obj)
        return json.dumps(obj)

    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # ── Load per-variant results ─────────────────────────────────────
    variants_data = {}
    for vlabel, vinfo in sorted(variant_results.items()):
        rdir = Path(vinfo["results_dir"])
        analyzer = ResultsAnalyzer(
            results_dir=rdir,
            parameter_definitions=vinfo.get("calibration_parameters", []))
        try:
            summary = analyzer.get_summary()
        except ValueError:
            summary = None

        vdata = {
            "analyzer": analyzer,
            "summary": summary,
            "label": vlabel,
            "config": vinfo.get("config", {}),
        }
        # Load matched data if available
        mdata_path = vinfo.get("matched_data_csv")
        if mdata_path and Path(mdata_path).exists():
            vdata["matched_data"] = pd.read_csv(mdata_path)
        else:
            # Try loading from results dir
            candidate = rdir / "matched_data.csv"
            if candidate.exists():
                vdata["matched_data"] = pd.read_csv(candidate)
            else:
                vdata["matched_data"] = None

        variants_data[vlabel] = vdata

    if not variants_data:
        logger.warning("No variant results found — skipping basin report")
        return None

    # ── Load map data using fiona (reuse same logic) ─────────────────
    map_data = {}
    if shapefile_dir and HAS_FIONA:
        shapefile_dir_p = Path(shapefile_dir)
        net_features = []
        try:
            hru_files = list(shapefile_dir_p.glob("*_basinHru.gpkg"))
            net_files = list(shapefile_dir_p.glob("*_riverNetwork.gpkg"))
            if hru_files:
                with fiona.open(str(hru_files[0])) as src:
                    features = list(src)
                hru_fc = {"type": "FeatureCollection", "features": features}
                map_data["hru_geojson"] = json.dumps(hru_fc)
                map_data["n_hrus"] = len(features)
                total_area = 0.0
                min_x, min_y, max_x, max_y = 180, 90, -180, -90
                for f in features:
                    p = f["properties"]
                    if p.get("unitarea"):
                        total_area += float(p["unitarea"])

                    def _walk(c):
                        if (isinstance(c, (list, tuple)) and len(c) >= 2
                                and isinstance(c[0], (int, float))):
                            return [c]
                        pts = []
                        for item in c:
                            pts.extend(_walk(item))
                        return pts

                    for pt in _walk(f["geometry"]["coordinates"]):
                        x, y = float(pt[0]), float(pt[1])
                        min_x, min_y = min(min_x, x), min(min_y, y)
                        max_x, max_y = max(max_x, x), max(max_y, y)
                if total_area > 0:
                    map_data["total_area_km2"] = total_area
                map_data["center_lat"] = (min_y + max_y) / 2
                map_data["center_lon"] = (min_x + max_x) / 2
                map_data["bounds"] = [[min_y, min_x], [max_y, max_x]]

            if net_files:
                with fiona.open(str(net_files[0])) as src:
                    net_features = list(src)
                net_fc = {"type": "FeatureCollection",
                          "features": net_features}
                map_data["network_geojson"] = json.dumps(net_fc)
                map_data["n_reaches"] = len(net_features)
                total_len, max_ord = 0.0, 0
                for f in net_features:
                    p = f["properties"]
                    if p.get("Length"):
                        total_len += float(p["Length"])
                    if p.get("order"):
                        max_ord = max(max_ord, int(p["order"]))
                if total_len > 0:
                    map_data["total_length_km"] = total_len / 1000
                if max_ord > 0:
                    map_data["max_order"] = max_ord

            if observation_data_path and net_features:
                obs_path = Path(observation_data_path)
                if obs_path.exists():
                    obs_df = pd.read_csv(obs_path)
                    if "reach_id" in obs_df.columns:
                        stations = []
                        for sid in obs_df["reach_id"].unique():
                            for nf in net_features:
                                if nf["properties"].get("LINKNO") == sid:
                                    coords = nf["geometry"]["coordinates"]
                                    mid = len(coords) // 2
                                    pt = coords[mid] if isinstance(
                                        coords[0], (list, tuple)) else coords
                                    stations.append({
                                        "lat": float(pt[1]),
                                        "lon": float(pt[0]),
                                        "reach_id": int(sid),
                                        "species": obs_df[
                                            obs_df["reach_id"] == sid
                                        ]["species"].unique().tolist(),
                                        "n_obs": len(obs_df[
                                            obs_df["reach_id"] == sid]),
                                    })
                                    break
                        if stations:
                            map_data["stations"] = stations
        except Exception as e:
            logger.warning(f"Basin report map data error: {e}")

    has_map = bool(map_data.get("hru_geojson") or
                   map_data.get("network_geojson"))

    # ── Build variant comparison table ───────────────────────────────
    comparison_rows = []
    for vlabel, vd in sorted(variants_data.items()):
        s = vd["summary"]
        row = {
            "variant": vlabel,
            "best_obj": s.best_objective if s else float('nan'),
            "n_evals": s.n_evaluations if s else 0,
            "n_params": len(vd["analyzer"].parameter_definitions),
            "conv_eval": s.convergence_eval if s else 0,
        }
        # Per-species KGE from matched data
        md = vd.get("matched_data")
        if md is not None and not md.empty:
            from ..objective_functions import ObjectiveFunction
            for sp in sorted(md["species"].unique()):
                sp_d = md[md["species"] == sp]
                obs = sp_d["observed"].values
                sim = sp_d["simulated"].values
                valid = ~(np.isnan(obs) | np.isnan(sim))
                if valid.sum() >= 2:
                    row[f"KGE_{sp}"] = float(
                        ObjectiveFunction.kge(obs[valid], sim[valid]))
        comparison_rows.append(row)

    # ── Build HTML ───────────────────────────────────────────────────
    n_variants = len(variants_data)
    H = []

    # <head>
    H.append(f"""<!DOCTYPE html>
<html lang="en" data-theme="light">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>OpenWQ Basin Report — {basin_id}</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&family=JetBrains+Mono:wght@400;500&display=swap" rel="stylesheet">
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js" charset="utf-8"></script>
<link rel="stylesheet" href="https://unpkg.com/leaflet@1.9.4/dist/leaflet.css" crossorigin="">
<script src="https://unpkg.com/leaflet@1.9.4/dist/leaflet.js" crossorigin=""></script>
<style>
:root,[data-theme="light"]{{
  --bg:#f0f2f5;--bg2:#fff;--card:#fff;--card-hover:#fff;
  --border:#e0e3e8;--border2:#d0d4da;
  --text:#1a1a2e;--text2:#4a4a68;--text3:#8a8aa0;
  --accent:#6366f1;--accent2:#8b5cf6;--accent-soft:rgba(99,102,241,.08);
  --green:#10b981;--red:#ef4444;--amber:#f59e0b;
  --shadow:0 1px 3px rgba(0,0,0,.06),0 1px 2px rgba(0,0,0,.04);
  --shadow-lg:0 10px 25px rgba(0,0,0,.07);
  --sidebar-bg:rgba(255,255,255,.85);
  --header-from:#6366f1;--header-to:#8b5cf6;
  --th-bg:#6366f1;--bar-track:#e5e7eb;
  --bar-fill:linear-gradient(90deg,#6366f1,#8b5cf6);
}}
[data-theme="dark"]{{
  --bg:#0f1117;--bg2:#1a1b26;--card:#1e1f2e;--card-hover:#252640;
  --border:#2a2b3d;--border2:#3a3b4d;
  --text:#e2e2f0;--text2:#a0a0b8;--text3:#6a6a80;
  --accent:#818cf8;--accent2:#a78bfa;--accent-soft:rgba(129,140,248,.12);
  --green:#34d399;--red:#f87171;--amber:#fbbf24;
  --shadow:0 1px 3px rgba(0,0,0,.3);
  --shadow-lg:0 10px 25px rgba(0,0,0,.4);
  --sidebar-bg:rgba(15,17,23,.92);
  --header-from:#4f46e5;--header-to:#7c3aed;
  --th-bg:#4338ca;--bar-track:#374151;
  --bar-fill:linear-gradient(90deg,#818cf8,#a78bfa);
}}
*{{box-sizing:border-box;margin:0;padding:0}}
html{{scroll-behavior:smooth}}
body{{font-family:'Inter',system-ui,sans-serif;background:var(--bg);color:var(--text);line-height:1.65;transition:background .3s,color .3s}}
.layout{{display:flex;min-height:100vh}}
.sidebar{{position:sticky;top:0;height:100vh;width:220px;flex-shrink:0;background:var(--sidebar-bg);backdrop-filter:blur(12px);border-right:1px solid var(--border);padding:1.2rem .8rem;overflow-y:auto;z-index:100}}
.sidebar .logo{{font-weight:700;font-size:.85rem;color:var(--accent);margin-bottom:1rem}}
.sidebar nav a{{display:block;padding:.35rem .7rem;margin-bottom:2px;border-radius:6px;font-size:.8rem;color:var(--text2);text-decoration:none;transition:all .15s}}
.sidebar nav a:hover,.sidebar nav a.active{{background:var(--accent-soft);color:var(--accent)}}
.sidebar nav a.active{{font-weight:600}}
.theme-toggle{{margin-top:auto;padding-top:1rem;border-top:1px solid var(--border)}}
.theme-btn{{display:flex;align-items:center;gap:.5rem;width:100%;padding:.45rem .7rem;border:1px solid var(--border);border-radius:8px;background:var(--card);color:var(--text2);font-size:.78rem;cursor:pointer;transition:all .15s;font-family:inherit}}
.theme-btn:hover{{border-color:var(--accent);color:var(--accent)}}
.theme-btn svg{{width:16px;height:16px;flex-shrink:0}}
.main{{flex:1;min-width:0}}
.header{{background:linear-gradient(135deg,var(--header-from),var(--header-to));color:#fff;padding:2.5rem 2.5rem 2rem;position:relative;overflow:hidden}}
.header::before{{content:'';position:absolute;inset:0;background:radial-gradient(circle at 80% 20%,rgba(255,255,255,.12) 0%,transparent 60%)}}
.header h1{{font-size:1.7rem;font-weight:700;letter-spacing:-.4px;position:relative}}
.header .subtitle{{opacity:.8;font-size:.88rem;margin-top:.3rem;position:relative}}
.header .meta{{display:flex;gap:1.5rem;margin-top:.9rem;font-size:.78rem;opacity:.7;flex-wrap:wrap;position:relative}}
.header .meta span{{background:rgba(255,255,255,.12);padding:.2rem .6rem;border-radius:20px}}
.container{{max-width:1100px;margin:0 auto;padding:2rem 2rem 1rem}}
.section{{margin-bottom:2.5rem;opacity:0;transform:translateY(16px);transition:opacity .5s ease,transform .5s ease}}
.section.visible{{opacity:1;transform:translateY(0)}}
.section h2{{font-size:1.15rem;font-weight:700;margin-bottom:1rem;display:flex;align-items:center;gap:.5rem;letter-spacing:-.3px}}
.section h2::before{{content:'';width:4px;height:1.2em;border-radius:2px;background:linear-gradient(180deg,var(--accent),var(--accent2));flex-shrink:0}}
.card{{background:var(--card);border:1px solid var(--border);border-radius:12px;padding:1.3rem;margin-bottom:1rem;box-shadow:var(--shadow);transition:all .3s}}
.card:hover{{box-shadow:var(--shadow-lg);border-color:var(--border2)}}
.kpi-grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(140px,1fr));gap:.8rem}}
.kpi{{text-align:center;padding:1.1rem .8rem;background:var(--card);border:1px solid var(--border);border-radius:12px;box-shadow:var(--shadow);transition:all .25s;position:relative;overflow:hidden}}
.kpi:hover{{transform:translateY(-2px);box-shadow:var(--shadow-lg)}}
.kpi::before{{content:'';position:absolute;top:0;left:0;right:0;height:3px;background:var(--bar-fill);border-radius:12px 12px 0 0}}
.kpi .value{{font-size:1.5rem;font-weight:700;color:var(--accent);font-family:'JetBrains Mono',monospace}}
.kpi .label{{font-size:.7rem;color:var(--text3);text-transform:uppercase;letter-spacing:.8px;margin-top:.3rem}}
.table-wrap{{overflow-x:auto;border-radius:10px;border:1px solid var(--border)}}
table{{width:100%;border-collapse:collapse;font-size:.85rem}}
th,td{{padding:.55rem .75rem;text-align:left;border-bottom:1px solid var(--border)}}
th{{background:var(--th-bg);color:#fff;font-weight:600;font-size:.72rem;text-transform:uppercase;letter-spacing:.5px}}
tr:nth-child(even){{background:var(--accent-soft)}}
td.num{{text-align:right;font-family:'JetBrains Mono',monospace;font-size:.8rem}}
.badge{{display:inline-flex;align-items:center;gap:.25rem;padding:.2rem .55rem;border-radius:20px;font-size:.7rem;font-weight:600}}
.badge-good{{background:rgba(16,185,129,.12);color:var(--green)}}
.badge-fair{{background:rgba(245,158,11,.12);color:var(--amber)}}
.badge-poor{{background:rgba(239,68,68,.12);color:var(--red)}}
.map-container{{border-radius:10px;overflow:hidden;border:1px solid var(--border)}}
#basinMap{{height:520px;width:100%;background:var(--bg2);z-index:1}}
.map-info-grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(140px,1fr));gap:.6rem;margin-top:.8rem}}
.map-info-item{{text-align:center;padding:.6rem;background:var(--accent-soft);border-radius:8px}}
.map-info-item .val{{font-size:1.1rem;font-weight:700;color:var(--accent);font-family:'JetBrains Mono',monospace}}
.map-info-item .lbl{{font-size:.68rem;color:var(--text3);text-transform:uppercase;letter-spacing:.5px;margin-top:.15rem}}
.leaflet-popup-content-wrapper{{border-radius:8px !important;font-family:'Inter',sans-serif;font-size:.82rem}}
.leaflet-popup-content{{margin:10px 14px !important}}
.leaflet-popup-content b{{color:#6366f1}}
[data-theme="dark"] .leaflet-tile{{filter:brightness(.8) contrast(1.1) saturate(.85)}}
[data-theme="dark"] .leaflet-popup-content-wrapper{{background:#1e1f2e !important;color:#e2e2f0 !important}}
[data-theme="dark"] .leaflet-popup-tip{{border-top-color:#1e1f2e !important}}
.legend-box{{background:var(--card);border:1px solid var(--border);border-radius:8px;padding:8px 12px;font-size:.75rem;line-height:1.6;color:var(--text)}}
.legend-box i{{width:14px;height:14px;display:inline-block;margin-right:6px;border-radius:3px;vertical-align:middle}}
.variant-card{{background:var(--card);border:1px solid var(--border);border-radius:12px;padding:1.2rem;margin-bottom:.8rem;box-shadow:var(--shadow);transition:all .3s}}
.variant-card:hover{{box-shadow:var(--shadow-lg)}}
.variant-card h3{{font-size:1rem;font-weight:600;color:var(--accent);margin-bottom:.6rem}}
.variant-card .param-grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(200px,1fr));gap:.3rem;font-size:.82rem}}
.variant-card .param-item{{display:flex;justify-content:space-between;padding:.25rem 0;border-bottom:1px dotted var(--border)}}
.range-bar{{display:inline-flex;align-items:center;gap:.4rem}}
.range-track{{width:90px;height:8px;background:var(--bar-track);border-radius:4px;overflow:hidden}}
.range-fill{{height:100%;border-radius:4px;background:var(--bar-fill);transition:width .6s ease}}
.range-pct{{font-size:.7rem;color:var(--text3);font-family:'JetBrains Mono',monospace}}
.footer{{text-align:center;padding:2rem;color:var(--text3);font-size:.78rem;border-top:1px solid var(--border)}}
.footer a{{color:var(--accent);text-decoration:none}}
@media(max-width:900px){{.sidebar{{display:none}}.container{{padding:1.2rem}}.kpi-grid{{grid-template-columns:repeat(2,1fr)}}}}
</style>
</head>
<body>""")

    # ── TOC ──────────────────────────────────────────────────────────
    toc = [("overview", "Overview")]
    if has_map:
        toc.append(("basin-map", "Basin Map"))
    toc.append(("comparison", "Variant Comparison"))
    for vlabel in sorted(variants_data):
        toc.append((f"variant-{vlabel}", f"Variant {vlabel}"))

    H.append('<div class="layout"><aside class="sidebar" id="sidebar">')
    H.append('<div class="logo">OpenWQ Basin Report</div>')
    H.append('<nav id="sidebarNav">')
    for anchor, label in toc:
        H.append(f'<a href="#{anchor}" data-section="{anchor}">{label}</a>')
    H.append('</nav>')
    H.append("""<div class="theme-toggle">
<button class="theme-btn" id="themeToggle" aria-label="Toggle theme">
<svg id="themeIcon" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
<circle cx="12" cy="12" r="5"/><path d="M12 1v2M12 21v2M4.22 4.22l1.42 1.42M18.36 18.36l1.42 1.42M1 12h2M21 12h2M4.22 19.78l1.42-1.42M18.36 5.64l1.42-1.42"/>
</svg><span id="themeLabel">Dark mode</span></button></div></aside><div class="main">""")

    # ── Header ───────────────────────────────────────────────────────
    H.append(f"""
<div class="header">
  <h1>Basin Report &mdash; {basin_id}</h1>
  <div class="subtitle">{n_variants} variant{"s" if n_variants != 1 else ""} compared</div>
  <div class="meta"><span>{now}</span><span>{n_variants} variants</span></div>
</div>
<div class="container">""")

    # ── Overview KPIs ────────────────────────────────────────────────
    H.append('<div class="section" id="overview"><h2>Overview</h2><div class="kpi-grid">')
    kpis = [(str(n_variants), "Variants"),
            (basin_id, "Basin ID")]
    if map_data.get("n_hrus"):
        kpis.append((str(map_data["n_hrus"]), "HRUs"))
    if map_data.get("n_reaches"):
        kpis.append((str(map_data["n_reaches"]), "Reaches"))
    if map_data.get("total_area_km2"):
        kpis.append((f"{map_data['total_area_km2']:.0f}", "Area (km²)"))
    best_overall = min(
        (r["best_obj"] for r in comparison_rows
         if not math.isnan(r["best_obj"])),
        default=float('nan'))
    if not math.isnan(best_overall):
        kpis.append((f"{best_overall:.4g}", "Best Obj (all)"))
    for v, l in kpis:
        H.append(f'<div class="kpi"><div class="value">{v}</div>'
                 f'<div class="label">{l}</div></div>')
    H.append('</div></div>')

    # ── Basin Map (same Leaflet code) ────────────────────────────────
    if has_map:
        bounds_js = json.dumps(map_data.get("bounds", [[44, -68], [46, -65]]))
        hru_geojson = map_data.get("hru_geojson", "null")
        net_geojson = map_data.get("network_geojson", "null")
        stations_js = json.dumps(map_data.get("stations", []))

        H.append('<div class="section" id="basin-map"><h2>Basin Map</h2><div class="card">')
        H.append('<div class="map-container"><div id="basinMap"></div></div>')
        info_items = []
        if map_data.get("n_hrus"):
            info_items.append((str(map_data["n_hrus"]), "HRUs / GRUs"))
        if map_data.get("n_reaches"):
            info_items.append((str(map_data["n_reaches"]), "River Reaches"))
        if map_data.get("total_area_km2"):
            info_items.append(
                (f"{map_data['total_area_km2']:.1f}", "Area (km²)"))
        if map_data.get("total_length_km"):
            info_items.append(
                (f"{map_data['total_length_km']:.1f}", "Network (km)"))
        if map_data.get("max_order"):
            info_items.append((str(map_data["max_order"]), "Max Order"))
        if map_data.get("stations"):
            info_items.append(
                (str(len(map_data["stations"])), "Obs. Stations"))
        if info_items:
            H.append('<div class="map-info-grid">')
            for val, lbl in info_items:
                H.append(f'<div class="map-info-item"><div class="val">{val}'
                         f'</div><div class="lbl">{lbl}</div></div>')
            H.append('</div>')
        H.append(f"""<script>
(function(){{
  var map=L.map('basinMap',{{zoomControl:true}});
  var osmLight=L.tileLayer('https://{{s}}.basemaps.cartocdn.com/light_all/{{z}}/{{x}}/{{y}}@2x.png',
    {{attribution:'&copy; CARTO &copy; OSM',maxZoom:18}});
  var osmTopo=L.tileLayer('https://{{s}}.tile.opentopomap.org/{{z}}/{{x}}/{{y}}.png',
    {{attribution:'&copy; OpenTopoMap',maxZoom:17}});
  var satellite=L.tileLayer('https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{{z}}/{{y}}/{{x}}',
    {{attribution:'&copy; Esri',maxZoom:18}});
  osmLight.addTo(map);
  map.fitBounds({bounds_js},{{padding:[30,30]}});
  var ordColors={{1:'#93c5fd',2:'#60a5fa',3:'#3b82f6',4:'#2563eb',5:'#1d4ed8',6:'#1e40af',7:'#1e3a8a'}};
  var ordWidths={{1:1.5,2:2,3:2.8,4:3.5,5:4.5,6:5.5,7:7}};
  var hruData={hru_geojson};
  var hruLayer=null;
  if(hruData){{
    var areas=hruData.features.map(function(f){{return f.properties.unitarea||f.properties.area||0}}).filter(function(a){{return a>0}});
    areas.sort(function(a,b){{return a-b}});
    var q=[0,.2,.4,.6,.8,1].map(function(p){{return areas[Math.min(Math.floor(p*(areas.length-1)),areas.length-1)]||0}});
    var hruPal=['#e0e7ff','#c7d2fe','#a5b4fc','#818cf8','#6366f1'];
    function hruColor(a){{for(var i=0;i<5;i++)if(a<=q[i+1])return hruPal[i];return hruPal[4]}}
    hruLayer=L.geoJSON(hruData,{{
      style:function(f){{var a=f.properties.unitarea||f.properties.area||0;return{{fillColor:hruColor(a),fillOpacity:.35,color:'#6366f1',weight:1.2,opacity:.6}}}},
      onEachFeature:function(f,layer){{var p=f.properties;var pop='<b>HRU</b><br>';if(p.GRU_ID)pop+='ID: <b>'+p.GRU_ID+'</b><br>';if(p.unitarea)pop+='Area: '+p.unitarea.toFixed(2)+' km²<br>';layer.bindPopup(pop);
        layer.on('mouseover',function(){{this.setStyle({{fillOpacity:.6,weight:2.5}})}});layer.on('mouseout',function(){{hruLayer.resetStyle(this)}})}}
    }}).addTo(map);
  }}
  var netData={net_geojson};var netLayer=null;
  if(netData){{
    netLayer=L.geoJSON(netData,{{
      style:function(f){{var o=f.properties.order||1;return{{color:ordColors[Math.min(o,7)]||'#3b82f6',weight:ordWidths[Math.min(o,7)]||2,opacity:.85,lineCap:'round',lineJoin:'round'}}}},
      onEachFeature:function(f,layer){{var p=f.properties;var pop='<b>Reach</b><br>';if(p.LINKNO)pop+='ID: <b>'+p.LINKNO+'</b><br>';if(p.order)pop+='Order: '+p.order+'<br>';if(p.Length)pop+='Length: '+(p.Length/1000).toFixed(2)+' km<br>';if(p.uparea)pop+='Upstream: '+p.uparea.toFixed(1)+' km²<br>';layer.bindPopup(pop)}}
    }}).addTo(map);
  }}
  var stations={stations_js};var stationLayer=L.layerGroup();
  stations.forEach(function(s){{
    var icon=L.divIcon({{className:'',html:'<div style="width:14px;height:14px;background:#ef4444;border:2.5px solid #fff;border-radius:50%;box-shadow:0 0 6px rgba(0,0,0,.35)"></div>',iconSize:[14,14],iconAnchor:[7,7]}});
    L.marker([s.lat,s.lon],{{icon:icon}}).bindPopup('<b>Obs Station</b><br>Reach: '+s.reach_id+'<br>Species: '+s.species.join(', ')+'<br>N obs: '+s.n_obs).addTo(stationLayer);
  }});stationLayer.addTo(map);
  var baseMaps={{"Light":osmLight,"Topo":osmTopo,"Satellite":satellite}};
  var overlays={{}};if(hruLayer)overlays["HRUs"]=hruLayer;if(netLayer)overlays["Rivers"]=netLayer;if(stations.length)overlays["Stations"]=stationLayer;
  L.control.layers(baseMaps,overlays,{{collapsed:false,position:'topright'}}).addTo(map);
  var legend=L.control({{position:'bottomright'}});
  legend.onAdd=function(){{var d=L.DomUtil.create('div','legend-box');d.innerHTML='<div style="font-weight:600;margin-bottom:4px">Legend</div><i style="background:#a5b4fc"></i> HRU<br><i style="background:#3b82f6;width:20px;height:3px;border-radius:2px"></i> River<br><i style="background:#ef4444;border-radius:50%"></i> Station';return d}};
  legend.addTo(map);L.control.scale({{imperial:false}}).addTo(map);
}})();
</script>""")
        H.append('</div></div>')

    # ── Variant Comparison Table ─────────────────────────────────────
    H.append('<div class="section" id="comparison"><h2>Variant Comparison</h2>')
    H.append('<div class="card"><div class="table-wrap"><table>')

    # Detect species columns
    species_cols = sorted(set(
        k for r in comparison_rows for k in r if k.startswith("KGE_")))

    H.append('<tr><th>Variant</th><th>Best Objective</th><th>Evaluations</th>'
             '<th>Parameters</th><th>Best at Eval</th>')
    for sc in species_cols:
        H.append(f'<th>{sc}</th>')
    H.append('</tr>')

    for row in comparison_rows:
        is_best = (abs(row["best_obj"] - best_overall) < 1e-10
                   if not math.isnan(best_overall) else False)
        style = ' style="background:rgba(16,185,129,.1);font-weight:600"' if is_best else ''
        H.append(f'<tr{style}><td><strong>Variant {row["variant"]}</strong></td>')
        H.append(f'<td class="num">{row["best_obj"]:.4g}</td>')
        H.append(f'<td class="num">{row["n_evals"]}</td>')
        H.append(f'<td class="num">{row["n_params"]}</td>')
        H.append(f'<td class="num">{row["conv_eval"]}</td>')
        for sc in species_cols:
            v = row.get(sc, float('nan'))
            if math.isnan(v):
                H.append('<td>&mdash;</td>')
            else:
                cls = ("badge-good" if v > 0.5 else
                       "badge-fair" if v > 0 else "badge-poor")
                H.append(f'<td><span class="badge {cls}">'
                         f'{v:.3f}</span></td>')
        H.append('</tr>')
    H.append('</table></div></div></div>')

    # ── Per-Variant Detail Cards ─────────────────────────────────────
    for vlabel, vd in sorted(variants_data.items()):
        s = vd["summary"]
        pdefs = vd["analyzer"].parameter_definitions
        H.append(f'<div class="section" id="variant-{vlabel}">')
        H.append(f'<h2>Variant {vlabel}</h2>')
        H.append('<div class="variant-card">')

        if s:
            H.append(f'<h3>Best Objective: {s.best_objective:.6g} '
                     f'(eval #{s.convergence_eval} of '
                     f'{s.n_evaluations})</h3>')

            # Full parameter table with bounds + range bar
            H.append('<div class="table-wrap" style="margin-top:.8rem">'
                     '<table><tr><th>Parameter</th><th>Best Value</th>'
                     '<th>Initial</th><th>Lower</th><th>Upper</th>'
                     '<th>Transform</th><th>Type</th>'
                     '<th>Position in Range</th></tr>')
            for pdef in pdefs:
                pname = pdef["name"]
                pval = s.best_parameters.get(pname, float('nan'))
                init_v = pdef.get("initial", float('nan'))
                lb, ub = pdef.get("bounds", (0, 1))
                transform = pdef.get("transform", "linear")
                file_type = pdef.get("file_type", "&mdash;")
                if ub > lb and not math.isnan(pval):
                    pct = max(0.0, min(100.0,
                              (pval - lb) / (ub - lb) * 100))
                    bar = (f'<div class="range-bar">'
                           f'<div class="range-track">'
                           f'<div class="range-fill" '
                           f'style="width:{pct:.0f}%"></div></div>'
                           f'<span class="range-pct">{pct:.0f}%'
                           f'</span></div>')
                else:
                    bar = "&mdash;"
                init_str = (f'{init_v:.4g}'
                            if not math.isnan(init_v) else "&mdash;")
                H.append(
                    f'<tr><td><strong>{pname}</strong></td>'
                    f'<td class="num">{pval:.6g}</td>'
                    f'<td class="num">{init_str}</td>'
                    f'<td class="num">{lb:.4g}</td>'
                    f'<td class="num">{ub:.4g}</td>'
                    f'<td>{transform}</td>'
                    f'<td>{file_type}</td><td>{bar}</td></tr>')
            H.append('</table></div>')
        elif pdefs:
            # No calibration results yet — just show parameter definitions
            H.append('<p style="color:var(--text3);margin-bottom:.6rem">'
                     'No calibration results yet</p>')
            H.append('<div class="table-wrap"><table>'
                     '<tr><th>Parameter</th><th>Initial</th>'
                     '<th>Lower</th><th>Upper</th>'
                     '<th>Transform</th><th>Type</th></tr>')
            for pdef in pdefs:
                pname = pdef["name"]
                init_v = pdef.get("initial", "—")
                lb, ub = pdef.get("bounds", (0, 1))
                transform = pdef.get("transform", "linear")
                ftype = pdef.get("file_type", "—")
                H.append(
                    f'<tr><td><strong>{pname}</strong></td>'
                    f'<td class="num">{init_v}</td>'
                    f'<td class="num">{lb:.4g}</td>'
                    f'<td class="num">{ub:.4g}</td>'
                    f'<td>{transform}</td><td>{ftype}</td></tr>')
            H.append('</table></div>')
        else:
            H.append('<p style="color:var(--text3)">No results or '
                     'parameter definitions available</p>')

        # Link to detailed variant report
        rdir = Path(variant_results[vlabel]["results_dir"])
        variant_report = rdir / "calibration_report.html"
        if variant_report.exists():
            rel = variant_report.resolve()
            H.append(f'<p style="margin-top:.8rem;font-size:.82rem">'
                     f'<a href="file://{rel}" '
                     f'style="color:var(--accent)">Open detailed variant '
                     f'report &rarr;</a></p>')

        H.append('</div></div>')

    # ── Footer & JS ──────────────────────────────────────────────────
    H.append(f'</div></div></div>')
    H.append(f'<div class="footer">Generated by <strong>OpenWQ Calibration '
             f'Framework</strong> &mdash; {now}<br>'
             f'<a href="https://github.com/DiogoCostaPT/openwq">'
             f'github.com/DiogoCostaPT/openwq</a></div>')
    H.append("""<script>
(function(){
  var stored=localStorage.getItem('owq-theme');
  if(stored)document.documentElement.setAttribute('data-theme',stored);
  var btn=document.getElementById('themeToggle'),icon=document.getElementById('themeIcon'),label=document.getElementById('themeLabel');
  function update(){var dark=document.documentElement.getAttribute('data-theme')==='dark';label.textContent=dark?'Light mode':'Dark mode';
    icon.innerHTML=dark?'<circle cx="12" cy="12" r="5"/><path d="M12 1v2M12 21v2M4.22 4.22l1.42 1.42M18.36 18.36l1.42 1.42M1 12h2M21 12h2M4.22 19.78l1.42-1.42M18.36 5.64l1.42-1.42"/>':'<path d="M21 12.79A9 9 0 1 1 11.21 3 7 7 0 0 0 21 12.79z"/>';}
  update();btn.addEventListener('click',function(){var cur=document.documentElement.getAttribute('data-theme');var next=cur==='dark'?'light':'dark';document.documentElement.setAttribute('data-theme',next);localStorage.setItem('owq-theme',next);update();});
})();
(function(){var secs=document.querySelectorAll('.section[id]'),links=document.querySelectorAll('#sidebarNav a');
  if(!secs.length||!links.length)return;
  var obs=new IntersectionObserver(function(entries){entries.forEach(function(e){if(e.isIntersecting){links.forEach(function(a){a.classList.remove('active');});var lk=document.querySelector('#sidebarNav a[data-section="'+e.target.id+'"]');if(lk)lk.classList.add('active');}});},{rootMargin:'-20% 0px -70% 0px'});
  secs.forEach(function(s){obs.observe(s);});
})();
(function(){var secs=document.querySelectorAll('.section');
  var obs=new IntersectionObserver(function(entries){entries.forEach(function(e){if(e.isIntersecting){e.target.classList.add('visible');obs.unobserve(e.target);}});},{threshold:.08});
  secs.forEach(function(s){obs.observe(s);});
})();
</script>
</body></html>""")

    report_path = out / f"basin_report_{basin_id}.html"
    with open(report_path, 'w', encoding='utf-8') as f:
        f.write("\n".join(H))

    logger.info(f"Basin report saved to {report_path}")
    return report_path
