# Copyright 2026, Diogo Costa, diogo.costa@uevora.pt
# This file is part of OpenWQ model.

"""
Morris Screening Method
=======================

Elementary effects-based sensitivity analysis for parameter screening.

Reference: Morris, M.D. (1991). "Factorial sampling plans for preliminary
computational experiments." Technometrics, 33(2), 161-174.

Morris screening is a computationally efficient method to identify
non-influential parameters that can be fixed at default values,
reducing the dimensionality of the calibration problem.
"""

import numpy as np
from typing import List, Tuple, Optional, Dict, Any
from dataclasses import dataclass
import logging
import json
from pathlib import Path

logger = logging.getLogger(__name__)

# Try to import SALib
try:
    from SALib.sample import morris as morris_sample
    from SALib.analyze import morris as morris_analyze
    SALIB_AVAILABLE = True
except ImportError:
    SALIB_AVAILABLE = False
    logger.warning("SALib not available. Install with: pip install SALib")


@dataclass
class MorrisResult:
    """Container for Morris screening results."""
    mu: np.ndarray              # Mean of elementary effects
    mu_star: np.ndarray         # Mean of absolute elementary effects
    sigma: np.ndarray           # Standard deviation of elementary effects
    param_names: List[str]
    influential_params: List[str]  # Parameters above threshold
    rankings: List[Tuple[str, float]]  # Parameters sorted by mu_star
    n_evaluations: int

    def to_dict(self) -> Dict:
        """Convert results to dictionary for JSON serialization."""
        return {
            "mu": {name: float(val) for name, val in zip(self.param_names, self.mu)},
            "mu_star": {name: float(val) for name, val in zip(self.param_names, self.mu_star)},
            "sigma": {name: float(val) for name, val in zip(self.param_names, self.sigma)},
            "influential_params": self.influential_params,
            "rankings": [(name, float(val)) for name, val in self.rankings],
            "n_evaluations": self.n_evaluations
        }

    def save(self, filepath: Path) -> None:
        """Save results to JSON file."""
        with open(filepath, 'w') as f:
            json.dump(self.to_dict(), f, indent=2)

    def print_summary(self) -> None:
        """Print a formatted summary of results."""
        print("\n" + "=" * 60)
        print("MORRIS SCREENING RESULTS")
        print("=" * 60)
        print(f"\nTotal evaluations: {self.n_evaluations}")
        print(f"\nParameter Rankings (by mu*):")
        print("-" * 40)
        for i, (name, mu_star) in enumerate(self.rankings, 1):
            influential = "***" if name in self.influential_params else "   "
            print(f"{i:3d}. {influential} {name:30s} mu*={mu_star:.4f}")
        print("-" * 40)
        print(f"Influential parameters (marked ***): {len(self.influential_params)}")
        print(f"Non-influential parameters: {len(self.param_names) - len(self.influential_params)}")
        print("=" * 60 + "\n")


class MorrisScreening:
    """
    Morris method for screening influential parameters.

    The Morris method evaluates elementary effects by varying one parameter
    at a time while holding others constant. Parameters with high mean
    absolute elementary effect (mu*) are considered influential.

    Parameters with high sigma/mu* ratio may indicate interaction effects.
    """

    def __init__(self,
                 param_names: List[str],
                 bounds: List[Tuple[float, float]],
                 num_trajectories: int = 10,
                 num_levels: int = 4,
                 seed: Optional[int] = None):
        """
        Initialize Morris screening.

        Parameters
        ----------
        param_names : List[str]
            Names of parameters
        bounds : List[Tuple[float, float]]
            (min, max) bounds for each parameter
        num_trajectories : int
            Number of Morris trajectories (r). More trajectories give more
            stable estimates but require more evaluations.
            Recommended: 10-20 for screening, 50+ for quantitative analysis.
        num_levels : int
            Number of levels in the parameter grid (p). Usually 4 or 6.
        seed : int
            Random seed for reproducibility
        """
        self.param_names = param_names
        self.bounds = bounds
        self.num_trajectories = num_trajectories
        self.num_levels = num_levels
        self.seed = seed
        self.n_params = len(param_names)

        # SALib problem definition
        self.problem = {
            'num_vars': self.n_params,
            'names': param_names,
            'bounds': bounds
        }

        # Calculate expected number of evaluations
        # Morris: (p+1) * r where p is num_levels and r is num_trajectories
        # Actually SALib uses (n_params + 1) * num_trajectories
        self.expected_evals = (self.n_params + 1) * num_trajectories

        logger.info(f"Morris screening: {self.n_params} params, "
                   f"{num_trajectories} trajectories, "
                   f"~{self.expected_evals} evaluations")

    def generate_samples(self) -> np.ndarray:
        """
        Generate Morris sampling design.

        Returns
        -------
        np.ndarray
            Sample matrix of shape (n_samples, n_params)
        """
        if not SALIB_AVAILABLE:
            return self._generate_samples_builtin()

        samples = morris_sample.sample(
            self.problem,
            N=self.num_trajectories,
            num_levels=self.num_levels,
            seed=self.seed
        )
        logger.info(f"Generated {len(samples)} Morris samples")
        return samples

    def _generate_samples_builtin(self) -> np.ndarray:
        """
        Built-in Morris sampling without SALib.

        Creates trajectories by starting from random point and
        moving one parameter at a time by delta = p/(2*(p-1)).
        """
        rng = np.random.default_rng(self.seed)
        p = self.num_levels
        delta = p / (2 * (p - 1))

        samples = []

        for _ in range(self.num_trajectories):
            # Random starting point on grid
            base = rng.integers(0, p - 1, size=self.n_params) / (p - 1)

            # Random order of parameter perturbations
            order = rng.permutation(self.n_params)

            trajectory = [base.copy()]
            current = base.copy()

            for j in order:
                # Decide direction
                if current[j] + delta <= 1.0:
                    current[j] += delta
                else:
                    current[j] -= delta
                trajectory.append(current.copy())

            samples.extend(trajectory)

        samples = np.array(samples)

        # Scale to actual bounds
        for j in range(self.n_params):
            lb, ub = self.bounds[j]
            samples[:, j] = lb + samples[:, j] * (ub - lb)

        return samples

    def analyze(self,
                samples: np.ndarray,
                outputs: np.ndarray,
                threshold: float = 0.1) -> MorrisResult:
        """
        Analyze Morris results.

        Parameters
        ----------
        samples : np.ndarray
            Parameter samples from generate_samples()
        outputs : np.ndarray
            Objective function values for each sample
        threshold : float
            Relative threshold for influential parameters.
            Parameters with mu* > threshold * max(mu*) are considered influential.

        Returns
        -------
        MorrisResult
            Sensitivity analysis results
        """
        if SALIB_AVAILABLE:
            return self._analyze_salib(samples, outputs, threshold)
        else:
            return self._analyze_builtin(samples, outputs, threshold)

    def _analyze_salib(self,
                       samples: np.ndarray,
                       outputs: np.ndarray,
                       threshold: float) -> MorrisResult:
        """Analyze using SALib."""
        si = morris_analyze.analyze(
            self.problem,
            samples,
            outputs,
            num_resamples=100,
            seed=self.seed
        )

        mu = np.array(si['mu'])
        mu_star = np.array(si['mu_star'])
        sigma = np.array(si['sigma'])

        # Identify influential parameters
        max_mu_star = np.max(mu_star)
        influential = [
            name for name, ms in zip(self.param_names, mu_star)
            if ms > threshold * max_mu_star
        ]

        # Create rankings
        rankings = sorted(
            zip(self.param_names, mu_star),
            key=lambda x: x[1],
            reverse=True
        )

        return MorrisResult(
            mu=mu,
            mu_star=mu_star,
            sigma=sigma,
            param_names=self.param_names,
            influential_params=influential,
            rankings=rankings,
            n_evaluations=len(outputs)
        )

    def _analyze_builtin(self,
                         samples: np.ndarray,
                         outputs: np.ndarray,
                         threshold: float) -> MorrisResult:
        """
        Built-in Morris analysis without SALib.

        Computes elementary effects from trajectories.
        """
        n_samples = len(outputs)
        trajectory_size = self.n_params + 1
        n_trajectories = n_samples // trajectory_size

        elementary_effects = {j: [] for j in range(self.n_params)}

        for t in range(n_trajectories):
            start_idx = t * trajectory_size
            trajectory_samples = samples[start_idx:start_idx + trajectory_size]
            trajectory_outputs = outputs[start_idx:start_idx + trajectory_size]

            for i in range(self.n_params):
                # Find which parameter changed between steps
                for step in range(trajectory_size - 1):
                    diff = trajectory_samples[step + 1] - trajectory_samples[step]
                    changed_param = np.argmax(np.abs(diff))

                    if changed_param == i and np.abs(diff[i]) > 1e-10:
                        # Compute elementary effect
                        delta_y = trajectory_outputs[step + 1] - trajectory_outputs[step]
                        delta_x = diff[i]
                        ee = delta_y / delta_x * (self.bounds[i][1] - self.bounds[i][0])
                        elementary_effects[i].append(ee)

        # Compute statistics
        mu = np.zeros(self.n_params)
        mu_star = np.zeros(self.n_params)
        sigma = np.zeros(self.n_params)

        for j in range(self.n_params):
            if elementary_effects[j]:
                ee_array = np.array(elementary_effects[j])
                mu[j] = np.mean(ee_array)
                mu_star[j] = np.mean(np.abs(ee_array))
                sigma[j] = np.std(ee_array)

        # Identify influential parameters
        max_mu_star = np.max(mu_star) if np.max(mu_star) > 0 else 1.0
        influential = [
            name for name, ms in zip(self.param_names, mu_star)
            if ms > threshold * max_mu_star
        ]

        # Create rankings
        rankings = sorted(
            zip(self.param_names, mu_star),
            key=lambda x: x[1],
            reverse=True
        )

        return MorrisResult(
            mu=mu,
            mu_star=mu_star,
            sigma=sigma,
            param_names=self.param_names,
            influential_params=influential,
            rankings=rankings,
            n_evaluations=len(outputs)
        )
