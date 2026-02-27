# Copyright 2026, Diogo Costa, diogo.costa@uevora.pt
# This file is part of OpenWQ model.

"""
HPC Support
===========

SLURM and PBS job management for HPC calibration.
"""

from .slurm_manager import SLURMManager, SLURMConfig, SLURMJob, JobState
from .job_coordinator import (
    JobCoordinator,
    AdaptiveJobCoordinator,
    CoordinatorMode,
    CoordinatorState,
    EvaluationTask
)
