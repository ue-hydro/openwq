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
Job Coordinator Module
======================

Coordinates asynchronous job execution for HPC calibration.
Supports both batch job arrays (for sensitivity analysis) and
adaptive worker pools (for DDS optimization).
"""

import json
import time
import threading
import queue
from pathlib import Path
from typing import Dict, List, Callable, Optional, Tuple, Any
from dataclasses import dataclass, field
from enum import Enum
import logging

from .slurm_manager import SLURMManager, SLURMConfig, JobState

logger = logging.getLogger(__name__)


class CoordinatorMode(Enum):
    """Job coordination modes."""
    BATCH_ARRAY = "batch_array"      # Submit all at once (sensitivity analysis)
    WORKER_POOL = "worker_pool"      # Adaptive submission (DDS optimization)
    SEQUENTIAL = "sequential"         # One at a time (debugging)


@dataclass
class EvaluationTask:
    """Represents a single evaluation task."""
    eval_id: int
    parameters: Dict[str, float]
    work_dir: Optional[Path] = None
    job_id: Optional[str] = None
    objective: Optional[float] = None
    status: str = "pending"  # pending, running, completed, failed
    start_time: Optional[float] = None
    end_time: Optional[float] = None


@dataclass
class CoordinatorState:
    """State of the job coordinator for checkpointing."""
    mode: CoordinatorMode
    total_evaluations: int
    completed_evaluations: int
    best_objective: float
    best_parameters: Dict[str, float]
    eval_history: List[Dict]
    pending_tasks: List[int]
    running_tasks: List[int]


class JobCoordinator:
    """
    Coordinates job execution for HPC calibration.

    This class manages the lifecycle of evaluation tasks, from parameter
    generation to result collection. It supports different execution modes
    optimized for different calibration phases.
    """

    def __init__(self,
                 slurm_config: SLURMConfig,
                 work_dir: Path,
                 mode: CoordinatorMode = CoordinatorMode.BATCH_ARRAY,
                 n_workers: int = 10):
        """
        Initialize job coordinator.

        Parameters
        ----------
        slurm_config : SLURMConfig
            SLURM configuration
        work_dir : Path
            Working directory for calibration
        mode : CoordinatorMode
            Execution mode
        n_workers : int
            Number of concurrent workers (for WORKER_POOL mode)
        """
        self.slurm_config = slurm_config
        self.work_dir = Path(work_dir)
        self.mode = mode
        self.n_workers = n_workers

        self.slurm_manager = SLURMManager(slurm_config, work_dir)

        # Task tracking
        self.tasks: Dict[int, EvaluationTask] = {}
        self.task_queue: queue.Queue = queue.Queue()
        self.result_queue: queue.Queue = queue.Queue()

        # State
        self.is_running = False
        self.total_submitted = 0
        self.total_completed = 0

        # Callbacks
        self.on_task_complete: Optional[Callable[[EvaluationTask], None]] = None

        # Create evaluations directory
        self.evals_dir = self.work_dir / "evaluations"
        self.evals_dir.mkdir(parents=True, exist_ok=True)

    def set_completion_callback(self, callback: Callable[[EvaluationTask], None]):
        """Set callback function called when a task completes."""
        self.on_task_complete = callback

    def prepare_evaluation_directory(self,
                                     eval_id: int,
                                     base_config_dir: Path) -> Path:
        """
        Prepare working directory for an evaluation.

        Parameters
        ----------
        eval_id : int
            Evaluation ID
        base_config_dir : Path
            Base configuration directory to copy

        Returns
        -------
        Path
            Path to evaluation working directory
        """
        eval_dir = self.evals_dir / f"eval_{eval_id:04d}"
        eval_dir.mkdir(parents=True, exist_ok=True)

        # Copy base configuration (handled by parameter_handler)
        # Just return the directory path here

        return eval_dir

    def submit_batch(self,
                     tasks: List[EvaluationTask],
                     setup_func: Callable[[int, Dict[str, float], Path], Path]
                     ) -> List[str]:
        """
        Submit a batch of evaluation tasks.

        Parameters
        ----------
        tasks : List[EvaluationTask]
            Tasks to submit
        setup_func : Callable
            Function to set up each evaluation directory
            Signature: (eval_id, parameters, base_dir) -> eval_dir

        Returns
        -------
        List[str]
            Job IDs of submitted tasks
        """
        if self.mode == CoordinatorMode.BATCH_ARRAY:
            return self._submit_as_array(tasks, setup_func)
        elif self.mode == CoordinatorMode.WORKER_POOL:
            return self._submit_to_workers(tasks, setup_func)
        else:
            return self._submit_sequential(tasks, setup_func)

    def _submit_as_array(self,
                         tasks: List[EvaluationTask],
                         setup_func: Callable
                         ) -> List[str]:
        """Submit tasks as a SLURM job array."""
        eval_ids = []
        eval_dirs = []
        params_list = []

        for task in tasks:
            # Set up directory
            eval_dir = setup_func(task.eval_id, task.parameters, self.evals_dir)
            task.work_dir = eval_dir
            task.status = "pending"
            task.start_time = time.time()

            self.tasks[task.eval_id] = task

            eval_ids.append(task.eval_id)
            eval_dirs.append(eval_dir)
            params_list.append(task.parameters)

        # Submit array job
        job_id = self.slurm_manager.submit_array_job(
            eval_ids, eval_dirs, params_list
        )

        if job_id:
            # Update tasks with job IDs
            for i, task in enumerate(tasks):
                task.job_id = f"{job_id}_{i}"
                task.status = "running"

            self.total_submitted += len(tasks)
            return [f"{job_id}_{i}" for i in range(len(tasks))]

        return []

    def _submit_to_workers(self,
                           tasks: List[EvaluationTask],
                           setup_func: Callable
                           ) -> List[str]:
        """Submit tasks to worker pool (individual jobs)."""
        job_ids = []

        for task in tasks:
            eval_dir = setup_func(task.eval_id, task.parameters, self.evals_dir)
            task.work_dir = eval_dir
            task.status = "pending"
            task.start_time = time.time()

            self.tasks[task.eval_id] = task

            job_id = self.slurm_manager.submit_single_job(
                task.eval_id, eval_dir, task.parameters
            )

            if job_id:
                task.job_id = job_id
                task.status = "running"
                job_ids.append(job_id)
                self.total_submitted += 1

        return job_ids

    def _submit_sequential(self,
                           tasks: List[EvaluationTask],
                           setup_func: Callable
                           ) -> List[str]:
        """Submit tasks sequentially (for debugging)."""
        job_ids = []

        for task in tasks:
            eval_dir = setup_func(task.eval_id, task.parameters, self.evals_dir)
            task.work_dir = eval_dir
            task.status = "pending"
            task.start_time = time.time()

            self.tasks[task.eval_id] = task

            job_id = self.slurm_manager.submit_single_job(
                task.eval_id, eval_dir, task.parameters
            )

            if job_id:
                task.job_id = job_id
                task.status = "running"
                job_ids.append(job_id)
                self.total_submitted += 1

                # Wait for completion before submitting next
                self.slurm_manager.wait_for_jobs([job_id])

        return job_ids

    def wait_for_completion(self,
                            timeout: Optional[float] = None,
                            poll_interval: int = 30) -> Dict[int, float]:
        """
        Wait for all submitted tasks to complete.

        Parameters
        ----------
        timeout : float, optional
            Maximum wait time in seconds
        poll_interval : int
            Seconds between status checks

        Returns
        -------
        Dict[int, float]
            Mapping of eval_id to objective value
        """
        start_time = time.time()
        results = {}

        while True:
            # Check timeout
            if timeout and (time.time() - start_time) > timeout:
                logger.warning("Timeout waiting for tasks to complete")
                break

            # Update job states
            self.slurm_manager.update_all_job_states()

            # Check for completed tasks
            all_done = True
            for eval_id, task in self.tasks.items():
                if task.status == "running":
                    # Check if job completed
                    if task.job_id in self.slurm_manager.jobs:
                        job = self.slurm_manager.jobs[task.job_id]
                        if job.state == JobState.COMPLETED:
                            task.status = "completed"
                            task.end_time = time.time()
                            self.total_completed += 1

                            # Call completion callback
                            if self.on_task_complete:
                                self.on_task_complete(task)

                        elif job.state in [JobState.FAILED, JobState.CANCELLED,
                                           JobState.TIMEOUT]:
                            task.status = "failed"
                            task.end_time = time.time()
                            logger.warning(f"Task {eval_id} failed: {job.state}")

                    # Also check for completion marker file
                    if task.work_dir and (task.work_dir / "COMPLETED").exists():
                        if task.status == "running":
                            task.status = "completed"
                            task.end_time = time.time()
                            self.total_completed += 1

                            if self.on_task_complete:
                                self.on_task_complete(task)

                if task.status == "running":
                    all_done = False
                elif task.status == "completed":
                    if task.objective is not None:
                        results[eval_id] = task.objective

            if all_done:
                break

            logger.info(f"Progress: {self.total_completed}/{self.total_submitted} tasks complete")
            time.sleep(poll_interval)

        return results

    def collect_results(self,
                        objective_func: Callable[[Path], float]
                        ) -> Dict[int, float]:
        """
        Collect objective function values from completed tasks.

        Parameters
        ----------
        objective_func : Callable
            Function to compute objective from output directory
            Signature: (output_dir) -> objective_value

        Returns
        -------
        Dict[int, float]
            Mapping of eval_id to objective value
        """
        results = {}

        for eval_id, task in self.tasks.items():
            if task.status == "completed" and task.work_dir:
                output_dir = task.work_dir / "openwq_out"
                if output_dir.exists():
                    try:
                        objective = objective_func(output_dir)
                        task.objective = objective
                        results[eval_id] = objective
                        logger.debug(f"Eval {eval_id}: objective = {objective:.6f}")
                    except Exception as e:
                        logger.warning(f"Failed to compute objective for eval {eval_id}: {e}")
                        task.status = "failed"

        return results

    def get_best_result(self) -> Optional[Tuple[int, Dict[str, float], float]]:
        """
        Get the best result so far.

        Returns
        -------
        Tuple[int, Dict, float] or None
            (eval_id, parameters, objective) or None if no results
        """
        best_eval_id = None
        best_params = None
        best_objective = float('inf')

        for eval_id, task in self.tasks.items():
            if task.objective is not None and task.objective < best_objective:
                best_objective = task.objective
                best_eval_id = eval_id
                best_params = task.parameters

        if best_eval_id is not None:
            return (best_eval_id, best_params, best_objective)
        return None

    def get_state(self) -> CoordinatorState:
        """Get current coordinator state for checkpointing."""
        eval_history = []
        for eval_id, task in self.tasks.items():
            eval_history.append({
                "eval_id": eval_id,
                "parameters": task.parameters,
                "objective": task.objective,
                "status": task.status,
                "start_time": task.start_time,
                "end_time": task.end_time,
            })

        best = self.get_best_result()
        best_obj = best[2] if best else float('inf')
        best_params = best[1] if best else {}

        pending = [eid for eid, t in self.tasks.items() if t.status == "pending"]
        running = [eid for eid, t in self.tasks.items() if t.status == "running"]

        return CoordinatorState(
            mode=self.mode,
            total_evaluations=self.total_submitted,
            completed_evaluations=self.total_completed,
            best_objective=best_obj,
            best_parameters=best_params,
            eval_history=eval_history,
            pending_tasks=pending,
            running_tasks=running,
        )

    def restore_state(self, state: CoordinatorState):
        """Restore coordinator state from checkpoint."""
        self.mode = state.mode
        self.total_submitted = state.total_evaluations
        self.total_completed = state.completed_evaluations

        for eval_data in state.eval_history:
            task = EvaluationTask(
                eval_id=eval_data["eval_id"],
                parameters=eval_data["parameters"],
                objective=eval_data.get("objective"),
                status=eval_data["status"],
                start_time=eval_data.get("start_time"),
                end_time=eval_data.get("end_time"),
            )
            # Reconstruct work_dir
            task.work_dir = self.evals_dir / f"eval_{task.eval_id:04d}"
            self.tasks[task.eval_id] = task

    def cancel_all(self):
        """Cancel all pending and running tasks."""
        cancelled = self.slurm_manager.cancel_all_jobs()
        logger.info(f"Cancelled {cancelled} jobs")

        for task in self.tasks.values():
            if task.status in ["pending", "running"]:
                task.status = "cancelled"
                task.end_time = time.time()

    def get_statistics(self) -> Dict[str, Any]:
        """Get coordinator statistics."""
        stats = {
            "mode": self.mode.value,
            "total_submitted": self.total_submitted,
            "total_completed": self.total_completed,
            "pending": sum(1 for t in self.tasks.values() if t.status == "pending"),
            "running": sum(1 for t in self.tasks.values() if t.status == "running"),
            "completed": sum(1 for t in self.tasks.values() if t.status == "completed"),
            "failed": sum(1 for t in self.tasks.values() if t.status == "failed"),
        }

        # Add runtime statistics
        runtimes = []
        for task in self.tasks.values():
            if task.start_time and task.end_time:
                runtimes.append(task.end_time - task.start_time)

        if runtimes:
            stats["avg_runtime_seconds"] = sum(runtimes) / len(runtimes)
            stats["total_runtime_seconds"] = sum(runtimes)
            stats["min_runtime_seconds"] = min(runtimes)
            stats["max_runtime_seconds"] = max(runtimes)

        # Add objective statistics
        objectives = [t.objective for t in self.tasks.values() if t.objective is not None]
        if objectives:
            stats["best_objective"] = min(objectives)
            stats["worst_objective"] = max(objectives)
            stats["mean_objective"] = sum(objectives) / len(objectives)

        return stats

    def save_results(self, output_file: Path):
        """
        Save all results to a JSON file.

        Parameters
        ----------
        output_file : Path
            Output file path
        """
        results = {
            "statistics": self.get_statistics(),
            "evaluations": []
        }

        for eval_id, task in sorted(self.tasks.items()):
            results["evaluations"].append({
                "eval_id": eval_id,
                "parameters": task.parameters,
                "objective": task.objective,
                "status": task.status,
                "runtime": (task.end_time - task.start_time)
                          if task.start_time and task.end_time else None,
            })

        with open(output_file, 'w') as f:
            json.dump(results, f, indent=2)

        logger.info(f"Results saved to {output_file}")


class AdaptiveJobCoordinator(JobCoordinator):
    """
    Adaptive job coordinator for DDS optimization.

    This coordinator submits jobs adaptively, using results from
    completed evaluations to guide parameter selection for new ones.
    """

    def __init__(self,
                 slurm_config: SLURMConfig,
                 work_dir: Path,
                 n_workers: int = 10,
                 min_pending: int = 5):
        """
        Initialize adaptive coordinator.

        Parameters
        ----------
        slurm_config : SLURMConfig
            SLURM configuration
        work_dir : Path
            Working directory
        n_workers : int
            Number of concurrent workers
        min_pending : int
            Minimum number of pending jobs to maintain
        """
        super().__init__(
            slurm_config, work_dir,
            mode=CoordinatorMode.WORKER_POOL,
            n_workers=n_workers
        )
        self.min_pending = min_pending

    def run_adaptive_optimization(self,
                                  generate_params: Callable[[int], Dict[str, float]],
                                  setup_func: Callable[[int, Dict, Path], Path],
                                  objective_func: Callable[[Path], float],
                                  max_evaluations: int,
                                  checkpoint_func: Optional[Callable] = None
                                  ) -> Tuple[Dict[str, float], float]:
        """
        Run adaptive optimization loop.

        Parameters
        ----------
        generate_params : Callable
            Function to generate parameters for next evaluation
            Signature: (eval_id) -> parameters_dict
        setup_func : Callable
            Function to set up evaluation directory
        objective_func : Callable
            Function to compute objective from output
        max_evaluations : int
            Maximum number of evaluations
        checkpoint_func : Callable, optional
            Function to save checkpoint after each evaluation

        Returns
        -------
        Tuple[Dict, float]
            Best parameters and objective value
        """
        eval_id = len(self.tasks)
        pending_count = 0
        running_count = 0

        while self.total_completed < max_evaluations:
            # Count pending and running
            pending_count = sum(1 for t in self.tasks.values() if t.status == "pending")
            running_count = sum(1 for t in self.tasks.values() if t.status == "running")

            # Submit new tasks if needed
            while (pending_count + running_count < self.n_workers and
                   eval_id < max_evaluations):
                params = generate_params(eval_id)
                task = EvaluationTask(eval_id=eval_id, parameters=params)

                eval_dir = setup_func(eval_id, params, self.evals_dir)
                task.work_dir = eval_dir
                task.start_time = time.time()

                job_id = self.slurm_manager.submit_single_job(
                    eval_id, eval_dir, params
                )

                if job_id:
                    task.job_id = job_id
                    task.status = "running"
                    self.tasks[eval_id] = task
                    self.total_submitted += 1
                    eval_id += 1
                    running_count += 1

            # Update job states
            self.slurm_manager.update_all_job_states()

            # Process completed jobs
            for eid, task in self.tasks.items():
                if task.status == "running":
                    # Check SLURM status
                    if task.job_id in self.slurm_manager.jobs:
                        job = self.slurm_manager.jobs[task.job_id]
                        if job.state == JobState.COMPLETED:
                            task.status = "completed"
                            task.end_time = time.time()

                            # Compute objective
                            if task.work_dir:
                                output_dir = task.work_dir / "openwq_out"
                                if output_dir.exists():
                                    try:
                                        task.objective = objective_func(output_dir)
                                        self.total_completed += 1

                                        logger.info(
                                            f"Eval {eid}: objective = {task.objective:.6f}"
                                        )

                                        if checkpoint_func:
                                            checkpoint_func()

                                    except Exception as e:
                                        logger.warning(f"Eval {eid} objective failed: {e}")
                                        task.status = "failed"

                        elif job.state in [JobState.FAILED, JobState.TIMEOUT]:
                            task.status = "failed"
                            task.end_time = time.time()
                            logger.warning(f"Eval {eid} failed: {job.state}")

            # Brief sleep
            time.sleep(5)

        # Return best result
        best = self.get_best_result()
        if best:
            return best[1], best[2]
        return {}, float('inf')
