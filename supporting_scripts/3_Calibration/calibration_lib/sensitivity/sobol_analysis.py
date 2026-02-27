# Copyright 2026, Diogo Costa, diogo.costa@uevora.pt
# This file is part of OpenWQ model.

"""
Sobol Sensitivity Analysis
==========================

Variance-based global sensitivity analysis using Sobol indices.

Reference: Sobol, I.M. (2001). "Global sensitivity indices for nonlinear
mathematical models and their Monte Carlo estimates." Mathematics and
Computers in Simulation, 55(1-3), 271-280.

Sobol indices decompose output variance into contributions from individual
parameters (first-order indices S_i) and interactions (total-order indices S_T).
"""

import numpy as np
from typing import List, Tuple, Optional, Dict
from dataclasses import dataclass
import logging
import json
from pathlib import Path

logger = logging.getLogger(__name__)

# Try to import SALib
try:
    from SALib.sample import saltelli
    from SALib.analyze import sobol as sobol_analyze
    SALIB_AVAILABLE = True
except ImportError:
    SALIB_AVAILABLE = False
    logger.warning("SALib not available. Install with: pip install SALib")


@dataclass
class SobolResult:
    """Container for Sobol sensitivity analysis results."""
    S1: np.ndarray           # First-order indices
    S1_conf: np.ndarray      # S1 confidence intervals
    ST: np.ndarray           # Total-order indices
    ST_conf: np.ndarray      # ST confidence intervals
    param_names: List[str]
    influential_params: List[str]  # Parameters above threshold
    rankings: List[Tuple[str, float, float]]  # (name, S1, ST) sorted by ST
    n_evaluations: int

    def to_dict(self) -> Dict:
        """Convert results to dictionary for JSON serialization."""
        return {
            "S1": {name: {"value": float(s1), "conf": float(c)}
                   for name, s1, c in zip(self.param_names, self.S1, self.S1_conf)},
            "ST": {name: {"value": float(st), "conf": float(c)}
                   for name, st, c in zip(self.param_names, self.ST, self.ST_conf)},
            "influential_params": self.influential_params,
            "rankings": [(name, float(s1), float(st)) for name, s1, st in self.rankings],
            "n_evaluations": self.n_evaluations
        }

    def save(self, filepath: Path) -> None:
        """Save results to JSON file."""
        with open(filepath, 'w') as f:
            json.dump(self.to_dict(), f, indent=2)

    def print_summary(self) -> None:
        """Print a formatted summary of results."""
        print("\n" + "=" * 70)
        print("SOBOL SENSITIVITY ANALYSIS RESULTS")
        print("=" * 70)
        print(f"\nTotal evaluations: {self.n_evaluations}")
        print(f"\nParameter Rankings (by Total-order index ST):")
        print("-" * 70)
        print(f"{'Rank':<5} {'Parameter':<30} {'S1 (first)':<15} {'ST (total)':<15}")
        print("-" * 70)
        for i, (name, s1, st) in enumerate(self.rankings, 1):
            influential = "***" if name in self.influential_params else "   "
            print(f"{i:<5} {influential} {name:<27} {s1:>12.4f} {st:>15.4f}")
        print("-" * 70)
        print(f"\nInfluential parameters (marked ***): {len(self.influential_params)}")
        print(f"Non-influential parameters: {len(self.param_names) - len(self.influential_params)}")
        print("\nInterpretation:")
        print("  S1 = First-order (main effect, no interactions)")
        print("  ST = Total-order (main effect + all interactions)")
        print("  Large ST-S1 gap indicates significant interaction effects")
        print("=" * 70 + "\n")


class SobolAnalysis:
    """
    Sobol variance-based sensitivity analysis.

    Computes first-order (S1) and total-order (ST) sensitivity indices.
    More computationally expensive than Morris but provides quantitative
    importance measures.

    Requires N*(2D+2) model evaluations where N is sample size and D is
    number of parameters.
    """

    def __init__(self,
                 param_names: List[str],
                 bounds: List[Tuple[float, float]],
                 num_samples: int = 1024,
                 seed: Optional[int] = None):
        """
        Initialize Sobol analysis.

        Parameters
        ----------
        param_names : List[str]
            Names of parameters
        bounds : List[Tuple[float, float]]
            (min, max) bounds for each parameter
        num_samples : int
            Base sample size (N). Total evaluations will be N*(2D+2).
            Recommended: 512-2048 for robust estimates.
        seed : int
            Random seed for reproducibility
        """
        if not SALIB_AVAILABLE:
            raise ImportError("SALib is required for Sobol analysis. "
                            "Install with: pip install SALib")

        self.param_names = param_names
        self.bounds = bounds
        self.num_samples = num_samples
        self.seed = seed
        self.n_params = len(param_names)

        # SALib problem definition
        self.problem = {
            'num_vars': self.n_params,
            'names': param_names,
            'bounds': bounds
        }

        # Calculate expected number of evaluations
        # Saltelli sampling: N*(2D+2) for computing second-order indices
        self.expected_evals = num_samples * (2 * self.n_params + 2)

        logger.info(f"Sobol analysis: {self.n_params} params, "
                   f"N={num_samples}, "
                   f"~{self.expected_evals} evaluations")

    def generate_samples(self) -> np.ndarray:
        """
        Generate Saltelli sampling design for Sobol analysis.

        Returns
        -------
        np.ndarray
            Sample matrix of shape (n_samples, n_params)
        """
        samples = saltelli.sample(
            self.problem,
            N=self.num_samples,
            calc_second_order=False,  # Reduces evaluations
            seed=self.seed
        )
        logger.info(f"Generated {len(samples)} Sobol samples")
        return samples

    def analyze(self,
                samples: np.ndarray,
                outputs: np.ndarray,
                threshold: float = 0.05) -> SobolResult:
        """
        Analyze Sobol results.

        Parameters
        ----------
        samples : np.ndarray
            Parameter samples from generate_samples()
        outputs : np.ndarray
            Objective function values for each sample
        threshold : float
            Threshold for influential parameters.
            Parameters with ST > threshold are considered influential.

        Returns
        -------
        SobolResult
            Sensitivity analysis results
        """
        si = sobol_analyze.analyze(
            self.problem,
            outputs,
            calc_second_order=False,
            num_resamples=100,
            seed=self.seed
        )

        S1 = np.array(si['S1'])
        S1_conf = np.array(si['S1_conf'])
        ST = np.array(si['ST'])
        ST_conf = np.array(si['ST_conf'])

        # Identify influential parameters (ST > threshold)
        influential = [
            name for name, st in zip(self.param_names, ST)
            if st > threshold
        ]

        # Create rankings sorted by total-order index
        rankings = sorted(
            [(name, s1, st) for name, s1, st in zip(self.param_names, S1, ST)],
            key=lambda x: x[2],
            reverse=True
        )

        return SobolResult(
            S1=S1,
            S1_conf=S1_conf,
            ST=ST,
            ST_conf=ST_conf,
            param_names=self.param_names,
            influential_params=influential,
            rankings=rankings,
            n_evaluations=len(outputs)
        )
