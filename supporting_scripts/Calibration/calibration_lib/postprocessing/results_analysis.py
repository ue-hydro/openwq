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
