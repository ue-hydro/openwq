# Copyright 2026, Diogo Costa, diogo.costa@uevora.pt
# This file is part of OpenWQ model.

"""
Base Optimizer Abstract Class
=============================

Defines the interface for all optimization algorithms.
"""

from abc import ABC, abstractmethod
from typing import Callable, List, Tuple, Optional
import numpy as np
from dataclasses import dataclass


@dataclass
class OptimizerResult:
    """Base result class for optimizers."""
    best_params: np.ndarray
    best_objective: float
    param_names: List[str]
    history: List[Tuple[int, float, np.ndarray]]
    n_evaluations: int


class BaseOptimizer(ABC):
    """Abstract base class for optimization algorithms."""

    @abstractmethod
    def optimize(self,
                 objective_func: Callable[[np.ndarray], float],
                 initial_point: Optional[np.ndarray] = None,
                 callback: Optional[Callable[[int, float, np.ndarray], None]] = None) -> OptimizerResult:
        """
        Run optimization.

        Parameters
        ----------
        objective_func : Callable
            Function to minimize
        initial_point : np.ndarray, optional
            Starting point
        callback : Callable, optional
            Called after each evaluation

        Returns
        -------
        OptimizerResult
            Optimization results
        """
        pass
