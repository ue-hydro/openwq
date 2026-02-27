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
Dynamically Dimensioned Search (DDS) Algorithm
==============================================

Efficient global optimization for expensive model calibration.

Reference: Tolson, B.A. and Shoemaker, C.A. (2007). "Dynamically dimensioned
search algorithm for computationally efficient watershed model calibration."
Water Resources Research, 43(1).

Key features:
- Designed for computationally expensive models (100-1000 evaluations)
- No derivative information required
- Automatically reduces search dimensionality as optimization progresses
- Simple, robust, and efficient for high-dimensional problems
"""

import numpy as np
from typing import Callable, List, Tuple, Optional, Dict, Any
import logging
from dataclasses import dataclass, field
import json
from pathlib import Path

logger = logging.getLogger(__name__)


@dataclass
class DDSResult:
    """Container for DDS optimization results."""
    best_params: np.ndarray
    best_objective: float
    param_names: List[str]
    history: List[Tuple[int, float, np.ndarray]]
    n_evaluations: int
    converged: bool = False
    convergence_reason: str = ""

    def to_dict(self) -> Dict:
        """Convert results to dictionary for JSON serialization."""
        return {
            "best_params": {
                name: float(val)
                for name, val in zip(self.param_names, self.best_params)
            },
            "best_objective": float(self.best_objective),
            "n_evaluations": self.n_evaluations,
            "converged": self.converged,
            "convergence_reason": self.convergence_reason,
            "history": [
                {"eval": int(h[0]), "objective": float(h[1]),
                 "params": [float(p) for p in h[2]]}
                for h in self.history
            ]
        }

    def save(self, filepath: Path) -> None:
        """Save results to JSON file."""
        with open(filepath, 'w') as f:
            json.dump(self.to_dict(), f, indent=2)


class DDS:
    """
    Dynamically Dimensioned Search optimizer.

    DDS is a simple, efficient, and robust algorithm for calibrating
    computationally expensive simulation models. It adaptively narrows
    the search space as the optimization progresses.

    The algorithm probabilistically perturbs fewer dimensions as iterations
    increase, focusing the search on promising regions while maintaining
    global exploration capability.
    """

    def __init__(self,
                 n_params: int,
                 bounds: List[Tuple[float, float]],
                 param_names: List[str],
                 max_evals: int,
                 r: float = 0.2,
                 seed: Optional[int] = None):
        """
        Initialize DDS optimizer.

        Parameters
        ----------
        n_params : int
            Number of parameters to optimize
        bounds : List[Tuple[float, float]]
            (min, max) bounds for each parameter
        param_names : List[str]
            Names of parameters (for reporting)
        max_evals : int
            Maximum number of objective function evaluations
        r : float
            Perturbation parameter (typically 0.2)
            Controls the standard deviation of perturbations:
            sigma_j = r * (upper_j - lower_j)
        seed : int
            Random seed for reproducibility
        """
        self.n_params = n_params
        self.bounds = np.array(bounds)
        self.param_names = param_names
        self.max_evals = max_evals
        self.r = r
        self.rng = np.random.default_rng(seed)

        # Validate inputs
        if len(bounds) != n_params:
            raise ValueError(f"bounds length ({len(bounds)}) != n_params ({n_params})")
        if len(param_names) != n_params:
            raise ValueError(f"param_names length ({len(param_names)}) != n_params ({n_params})")

        # Compute perturbation standard deviations
        self.sigma = self.r * (self.bounds[:, 1] - self.bounds[:, 0])

        logger.info(f"DDS initialized: {n_params} params, {max_evals} max evals, r={r}")

    def optimize(self,
                 objective_func: Callable[[np.ndarray], float],
                 initial_point: np.ndarray = None,
                 callback: Callable[[int, float, np.ndarray], None] = None,
                 early_stop_threshold: float = None,
                 early_stop_patience: int = 50) -> DDSResult:
        """
        Run DDS optimization.

        Parameters
        ----------
        objective_func : Callable
            Function that takes parameter array and returns scalar objective.
            The objective should be minimized (lower is better).
        initial_point : np.ndarray
            Starting point. If None, uses random point within bounds.
        callback : Callable
            Called after each evaluation: callback(eval_num, obj_val, params)
        early_stop_threshold : float
            Stop if objective falls below this value
        early_stop_patience : int
            Stop if no improvement for this many evaluations

        Returns
        -------
        DDSResult
            Optimization results including best parameters and history
        """
        history = []
        no_improvement_count = 0

        # Initialize with provided point or random
        if initial_point is not None:
            x_best = np.clip(initial_point, self.bounds[:, 0], self.bounds[:, 1])
        else:
            x_best = self._random_point()

        # Evaluate initial point
        f_best = objective_func(x_best)
        history.append((0, f_best, x_best.copy()))

        logger.info(f"DDS Eval 0: Initial objective = {f_best:.6f}")

        if callback:
            callback(0, f_best, x_best)

        # Check early stopping on initial point
        if early_stop_threshold is not None and f_best <= early_stop_threshold:
            return DDSResult(
                best_params=x_best,
                best_objective=f_best,
                param_names=self.param_names,
                history=history,
                n_evaluations=1,
                converged=True,
                convergence_reason=f"Objective below threshold ({early_stop_threshold})"
            )

        # Main DDS loop
        for i in range(1, self.max_evals):
            # Calculate probability of perturbing each dimension
            # P(i) = 1 - ln(i) / ln(max_evals)
            # This decreases from 1 to 0 as i goes from 1 to max_evals
            prob = 1.0 - np.log(i) / np.log(self.max_evals)

            # Generate new candidate by perturbing current best
            x_new = self._perturb(x_best, prob)

            # Evaluate candidate
            f_new = objective_func(x_new)

            # Update if better (greedy acceptance)
            if f_new < f_best:
                x_best = x_new
                f_best = f_new
                no_improvement_count = 0
                logger.info(f"DDS Eval {i}: New best = {f_best:.6f}")
            else:
                no_improvement_count += 1

            history.append((i, f_best, x_best.copy()))

            if callback:
                callback(i, f_best, x_best)

            # Check early stopping conditions
            if early_stop_threshold is not None and f_best <= early_stop_threshold:
                return DDSResult(
                    best_params=x_best,
                    best_objective=f_best,
                    param_names=self.param_names,
                    history=history,
                    n_evaluations=i + 1,
                    converged=True,
                    convergence_reason=f"Objective below threshold ({early_stop_threshold})"
                )

            if no_improvement_count >= early_stop_patience:
                return DDSResult(
                    best_params=x_best,
                    best_objective=f_best,
                    param_names=self.param_names,
                    history=history,
                    n_evaluations=i + 1,
                    converged=True,
                    convergence_reason=f"No improvement for {early_stop_patience} evaluations"
                )

        return DDSResult(
            best_params=x_best,
            best_objective=f_best,
            param_names=self.param_names,
            history=history,
            n_evaluations=self.max_evals,
            converged=False,
            convergence_reason="Maximum evaluations reached"
        )

    def _random_point(self) -> np.ndarray:
        """Generate random point within bounds."""
        return self.bounds[:, 0] + self.rng.random(self.n_params) * (
            self.bounds[:, 1] - self.bounds[:, 0])

    def _perturb(self, x: np.ndarray, prob: float) -> np.ndarray:
        """
        Perturb solution with DDS neighborhood function.

        At least one dimension is always perturbed (guaranteed exploration).
        Each dimension has probability 'prob' of being perturbed.
        """
        x_new = x.copy()

        # Decide which dimensions to perturb
        perturb_mask = self.rng.random(self.n_params) < prob

        # Ensure at least one dimension is perturbed
        if not perturb_mask.any():
            perturb_mask[self.rng.integers(self.n_params)] = True

        # Perturb selected dimensions
        for j in np.where(perturb_mask)[0]:
            x_new[j] = self._perturb_dim(x[j], j)

        return x_new

    def _perturb_dim(self, x_j: float, j: int) -> float:
        """
        Perturb a single dimension using reflection at bounds.

        Uses normal distribution with mean x_j and std sigma_j.
        Reflects at bounds to ensure feasibility.
        """
        lb, ub = self.bounds[j]

        # Normal perturbation
        x_new = x_j + self.sigma[j] * self.rng.standard_normal()

        # Reflect at bounds (keeps solution within feasible region)
        if x_new < lb:
            x_new = lb + (lb - x_new)
            # If still outside after reflection, snap to bound
            if x_new > ub:
                x_new = lb
        elif x_new > ub:
            x_new = ub - (x_new - ub)
            # If still outside after reflection, snap to bound
            if x_new < lb:
                x_new = ub

        return x_new

    def get_state(self) -> Dict[str, Any]:
        """Get optimizer state for checkpointing."""
        return {
            "rng_state": self.rng.bit_generator.state,
            "n_params": self.n_params,
            "bounds": self.bounds.tolist(),
            "param_names": self.param_names,
            "max_evals": self.max_evals,
            "r": self.r
        }

    def set_state(self, state: Dict[str, Any]) -> None:
        """Restore optimizer state from checkpoint."""
        self.rng.bit_generator.state = state["rng_state"]


class RandomSearch:
    """
    Simple random search baseline for comparison.
    """

    def __init__(self,
                 n_params: int,
                 bounds: List[Tuple[float, float]],
                 param_names: List[str],
                 max_evals: int,
                 seed: Optional[int] = None):
        """Initialize random search."""
        self.n_params = n_params
        self.bounds = np.array(bounds)
        self.param_names = param_names
        self.max_evals = max_evals
        self.rng = np.random.default_rng(seed)

    def optimize(self,
                 objective_func: Callable[[np.ndarray], float],
                 callback: Callable[[int, float, np.ndarray], None] = None) -> DDSResult:
        """Run random search optimization."""
        history = []

        # Random initial point
        x_best = self._random_point()
        f_best = objective_func(x_best)
        history.append((0, f_best, x_best.copy()))

        if callback:
            callback(0, f_best, x_best)

        for i in range(1, self.max_evals):
            x_new = self._random_point()
            f_new = objective_func(x_new)

            if f_new < f_best:
                x_best = x_new
                f_best = f_new
                logger.info(f"Random Eval {i}: New best = {f_best:.6f}")

            history.append((i, f_best, x_best.copy()))

            if callback:
                callback(i, f_best, x_best)

        return DDSResult(
            best_params=x_best,
            best_objective=f_best,
            param_names=self.param_names,
            history=history,
            n_evaluations=self.max_evals,
            converged=False,
            convergence_reason="Maximum evaluations reached"
        )

    def _random_point(self) -> np.ndarray:
        """Generate random point within bounds."""
        return self.bounds[:, 0] + self.rng.random(self.n_params) * (
            self.bounds[:, 1] - self.bounds[:, 0])
