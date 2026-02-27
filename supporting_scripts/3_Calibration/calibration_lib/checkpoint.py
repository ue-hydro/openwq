# Copyright 2026, Diogo Costa, diogo.costa@uevora.pt
# This file is part of OpenWQ model.

"""
Checkpoint Manager
==================

Handles saving and restoring calibration state for restart capability.
"""

import json
import pickle
from pathlib import Path
from typing import Dict, Any, Optional, List
from datetime import datetime
import logging
import numpy as np

logger = logging.getLogger(__name__)


class CheckpointManager:
    """
    Manages calibration state persistence for restart capability.
    """

    def __init__(self, checkpoint_dir: Path):
        """
        Initialize checkpoint manager.

        Parameters
        ----------
        checkpoint_dir : Path
            Directory to store checkpoint files
        """
        self.checkpoint_dir = Path(checkpoint_dir)
        self.checkpoint_dir.mkdir(parents=True, exist_ok=True)
        self.checkpoint_file = self.checkpoint_dir / "calibration_state.json"
        self.history_file = self.checkpoint_dir / "calibration_history.pkl"

    def save_state(self,
                   current_eval: int,
                   max_evals: int,
                   best_objective: float,
                   best_params: np.ndarray,
                   param_names: List[str],
                   history: List,
                   algorithm: str,
                   rng_state: Optional[Dict] = None,
                   additional_data: Optional[Dict] = None) -> None:
        """
        Save current calibration state to checkpoint.

        Parameters
        ----------
        current_eval : int
            Current evaluation number
        max_evals : int
            Maximum evaluations
        best_objective : float
            Best objective value found so far
        best_params : np.ndarray
            Best parameter values
        param_names : List[str]
            Parameter names
        history : List
            Evaluation history
        algorithm : str
            Algorithm name
        rng_state : Dict, optional
            Random number generator state
        additional_data : Dict, optional
            Any additional data to save
        """
        state = {
            "timestamp": datetime.now().isoformat(),
            "current_eval": current_eval,
            "max_evals": max_evals,
            "best_objective": float(best_objective),
            "best_params": {
                name: float(val)
                for name, val in zip(param_names, best_params)
            },
            "algorithm": algorithm,
            "n_params": len(param_names),
            "param_names": param_names,
        }

        if rng_state is not None:
            # Can't directly serialize numpy random state to JSON
            # Store as pickle separately
            rng_file = self.checkpoint_dir / "rng_state.pkl"
            with open(rng_file, 'wb') as f:
                pickle.dump(rng_state, f)
            state["has_rng_state"] = True

        if additional_data:
            state["additional_data"] = additional_data

        # Save state
        with open(self.checkpoint_file, 'w') as f:
            json.dump(state, f, indent=2)

        # Save history (may be large, use pickle)
        with open(self.history_file, 'wb') as f:
            pickle.dump(history, f)

        logger.debug(f"Checkpoint saved: eval {current_eval}/{max_evals}, "
                    f"best obj = {best_objective:.6f}")

    def load_state(self) -> Optional[Dict[str, Any]]:
        """
        Load calibration state from checkpoint.

        Returns
        -------
        Dict or None
            Loaded state, or None if no checkpoint exists
        """
        if not self.checkpoint_file.exists():
            logger.info("No checkpoint found")
            return None

        try:
            with open(self.checkpoint_file, 'r') as f:
                state = json.load(f)

            # Load history
            if self.history_file.exists():
                with open(self.history_file, 'rb') as f:
                    state['history'] = pickle.load(f)
            else:
                state['history'] = []

            # Load RNG state if exists
            rng_file = self.checkpoint_dir / "rng_state.pkl"
            if state.get("has_rng_state") and rng_file.exists():
                with open(rng_file, 'rb') as f:
                    state['rng_state'] = pickle.load(f)

            # Convert best_params back to array
            state['best_params_array'] = np.array([
                state['best_params'][name]
                for name in state['param_names']
            ])

            logger.info(f"Checkpoint loaded: eval {state['current_eval']}/{state['max_evals']}")
            return state

        except Exception as e:
            logger.error(f"Error loading checkpoint: {e}")
            return None

    def checkpoint_exists(self) -> bool:
        """Check if a checkpoint exists."""
        return self.checkpoint_file.exists()

    def clear_checkpoint(self) -> None:
        """Remove checkpoint files."""
        for f in [self.checkpoint_file, self.history_file,
                  self.checkpoint_dir / "rng_state.pkl"]:
            if f.exists():
                f.unlink()
        logger.info("Checkpoint cleared")
