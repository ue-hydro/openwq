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
SLURM Manager Module
====================

Manages SLURM job submission, monitoring, and result collection for HPC calibration.
"""

import subprocess
import time
import re
import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union
from dataclasses import dataclass, field
from enum import Enum
import logging

logger = logging.getLogger(__name__)


class JobState(Enum):
    """SLURM job states."""
    PENDING = "PENDING"
    RUNNING = "RUNNING"
    COMPLETED = "COMPLETED"
    FAILED = "FAILED"
    CANCELLED = "CANCELLED"
    TIMEOUT = "TIMEOUT"
    NODE_FAIL = "NODE_FAIL"
    UNKNOWN = "UNKNOWN"


@dataclass
class SLURMJob:
    """Represents a SLURM job."""
    job_id: str
    eval_id: int
    state: JobState = JobState.PENDING
    exit_code: Optional[int] = None
    start_time: Optional[float] = None
    end_time: Optional[float] = None
    work_dir: Optional[Path] = None

    @property
    def runtime(self) -> Optional[float]:
        """Get job runtime in seconds."""
        if self.start_time and self.end_time:
            return self.end_time - self.start_time
        return None


@dataclass
class SLURMConfig:
    """SLURM configuration settings."""
    partition: str = "standard"
    walltime: str = "04:00:00"
    nodes: int = 1
    tasks_per_node: int = 1
    memory: str = "8G"
    account: Optional[str] = None
    qos: Optional[str] = None
    max_concurrent_jobs: int = 50
    poll_interval: int = 30  # seconds

    # Apptainer settings
    sif_path: str = ""
    bind_path: str = ""
    executable_name: str = "mizuroute_lakes_cslm_openwq_fast"
    executable_args: str = "-g 1 1"
    file_manager_path: str = ""


class SLURMManager:
    """
    Manages SLURM job submission and monitoring.

    Supports both job arrays (for sensitivity analysis) and
    individual job submission (for DDS optimization).
    """

    def __init__(self,
                 config: SLURMConfig,
                 work_dir: Path,
                 template_dir: Optional[Path] = None):
        """
        Initialize SLURM manager.

        Parameters
        ----------
        config : SLURMConfig
            SLURM configuration settings
        work_dir : Path
            Working directory for calibration
        template_dir : Path, optional
            Directory containing job templates
        """
        self.config = config
        self.work_dir = Path(work_dir)
        self.template_dir = template_dir or (
            Path(__file__).parent.parent.parent / "hpc_templates"
        )

        self.jobs: Dict[str, SLURMJob] = {}
        self.job_history: List[SLURMJob] = []

        # Create job scripts directory
        self.scripts_dir = self.work_dir / "slurm_scripts"
        self.scripts_dir.mkdir(parents=True, exist_ok=True)

        # Create logs directory
        self.logs_dir = self.work_dir / "slurm_logs"
        self.logs_dir.mkdir(parents=True, exist_ok=True)

    def submit_single_job(self,
                          eval_id: int,
                          eval_dir: Path,
                          parameters: Dict[str, float]) -> Optional[str]:
        """
        Submit a single evaluation job.

        Parameters
        ----------
        eval_id : int
            Evaluation ID
        eval_dir : Path
            Evaluation working directory
        parameters : Dict[str, float]
            Parameter values for this evaluation

        Returns
        -------
        str or None
            Job ID if successful, None otherwise
        """
        # Write parameters file
        params_file = eval_dir / "parameters.json"
        with open(params_file, 'w') as f:
            json.dump(parameters, f, indent=2)

        # Generate job script
        script_content = self._generate_single_job_script(eval_id, eval_dir)
        script_path = self.scripts_dir / f"eval_{eval_id:04d}.sh"

        with open(script_path, 'w') as f:
            f.write(script_content)

        # Submit job
        try:
            result = subprocess.run(
                ["sbatch", str(script_path)],
                capture_output=True,
                text=True,
                check=True
            )

            # Parse job ID from output: "Submitted batch job 12345"
            match = re.search(r'Submitted batch job (\d+)', result.stdout)
            if match:
                job_id = match.group(1)

                job = SLURMJob(
                    job_id=job_id,
                    eval_id=eval_id,
                    state=JobState.PENDING,
                    work_dir=eval_dir
                )
                self.jobs[job_id] = job

                logger.info(f"Submitted job {job_id} for evaluation {eval_id}")
                return job_id

        except subprocess.CalledProcessError as e:
            logger.error(f"Failed to submit job for eval {eval_id}: {e.stderr}")
        except FileNotFoundError:
            logger.error("sbatch command not found - is SLURM available?")

        return None

    def submit_array_job(self,
                         eval_ids: List[int],
                         eval_dirs: List[Path],
                         parameters_list: List[Dict[str, float]]) -> Optional[str]:
        """
        Submit a SLURM array job for multiple evaluations.

        Parameters
        ----------
        eval_ids : List[int]
            List of evaluation IDs
        eval_dirs : List[Path]
            List of evaluation directories
        parameters_list : List[Dict]
            List of parameter dictionaries

        Returns
        -------
        str or None
            Array job ID if successful
        """
        n_evals = len(eval_ids)
        if n_evals == 0:
            return None

        # Write parameters files
        for eval_id, eval_dir, params in zip(eval_ids, eval_dirs, parameters_list):
            params_file = eval_dir / "parameters.json"
            with open(params_file, 'w') as f:
                json.dump({
                    "eval_id": eval_id,
                    "parameters": params
                }, f, indent=2)

        # Write array index mapping
        index_map = {i: {"eval_id": eval_ids[i], "dir": str(eval_dirs[i])}
                     for i in range(n_evals)}
        index_file = self.work_dir / "array_index_map.json"
        with open(index_file, 'w') as f:
            json.dump(index_map, f, indent=2)

        # Generate array job script
        max_concurrent = min(self.config.max_concurrent_jobs, n_evals)
        script_content = self._generate_array_job_script(
            n_evals, max_concurrent, index_file
        )
        script_path = self.scripts_dir / "array_job.sh"

        with open(script_path, 'w') as f:
            f.write(script_content)

        # Submit job
        try:
            result = subprocess.run(
                ["sbatch", str(script_path)],
                capture_output=True,
                text=True,
                check=True
            )

            match = re.search(r'Submitted batch job (\d+)', result.stdout)
            if match:
                job_id = match.group(1)

                # Create job objects for each array task
                for i, eval_id in enumerate(eval_ids):
                    task_job_id = f"{job_id}_{i}"
                    job = SLURMJob(
                        job_id=task_job_id,
                        eval_id=eval_id,
                        state=JobState.PENDING,
                        work_dir=eval_dirs[i]
                    )
                    self.jobs[task_job_id] = job

                logger.info(f"Submitted array job {job_id} with {n_evals} tasks")
                return job_id

        except subprocess.CalledProcessError as e:
            logger.error(f"Failed to submit array job: {e.stderr}")
        except FileNotFoundError:
            logger.error("sbatch command not found - is SLURM available?")

        return None

    def check_job_status(self, job_id: str) -> JobState:
        """
        Check the status of a SLURM job.

        Parameters
        ----------
        job_id : str
            Job ID to check

        Returns
        -------
        JobState
            Current job state
        """
        try:
            result = subprocess.run(
                ["sacct", "-j", job_id, "--format=State,ExitCode",
                 "--noheader", "--parsable2"],
                capture_output=True,
                text=True,
                check=True
            )

            lines = result.stdout.strip().split('\n')
            if lines and lines[0]:
                parts = lines[0].split('|')
                state_str = parts[0].strip()

                state_map = {
                    "PENDING": JobState.PENDING,
                    "RUNNING": JobState.RUNNING,
                    "COMPLETED": JobState.COMPLETED,
                    "FAILED": JobState.FAILED,
                    "CANCELLED": JobState.CANCELLED,
                    "TIMEOUT": JobState.TIMEOUT,
                    "NODE_FAIL": JobState.NODE_FAIL,
                }

                return state_map.get(state_str, JobState.UNKNOWN)

        except subprocess.CalledProcessError:
            pass
        except FileNotFoundError:
            logger.warning("sacct command not found")

        return JobState.UNKNOWN

    def update_all_job_states(self) -> Dict[str, JobState]:
        """
        Update states for all tracked jobs.

        Returns
        -------
        Dict[str, JobState]
            Mapping of job IDs to their current states
        """
        states = {}

        for job_id, job in self.jobs.items():
            if job.state in [JobState.PENDING, JobState.RUNNING]:
                new_state = self.check_job_status(job_id.split('_')[0])

                if new_state != job.state:
                    job.state = new_state

                    if new_state == JobState.RUNNING and job.start_time is None:
                        job.start_time = time.time()
                    elif new_state in [JobState.COMPLETED, JobState.FAILED,
                                       JobState.CANCELLED, JobState.TIMEOUT]:
                        job.end_time = time.time()

                        # Also check for completion marker file
                        if job.work_dir:
                            completed_file = job.work_dir / "COMPLETED"
                            if completed_file.exists():
                                job.state = JobState.COMPLETED

            states[job_id] = job.state

        return states

    def wait_for_jobs(self,
                      job_ids: List[str],
                      timeout: Optional[float] = None) -> Dict[str, JobState]:
        """
        Wait for specified jobs to complete.

        Parameters
        ----------
        job_ids : List[str]
            Job IDs to wait for
        timeout : float, optional
            Maximum wait time in seconds

        Returns
        -------
        Dict[str, JobState]
            Final states of all jobs
        """
        start_time = time.time()
        pending_jobs = set(job_ids)

        while pending_jobs:
            if timeout and (time.time() - start_time) > timeout:
                logger.warning("Timeout waiting for jobs")
                break

            self.update_all_job_states()

            completed = set()
            for job_id in pending_jobs:
                if job_id in self.jobs:
                    if self.jobs[job_id].state not in [JobState.PENDING, JobState.RUNNING]:
                        completed.add(job_id)

            pending_jobs -= completed

            if pending_jobs:
                logger.debug(f"Waiting for {len(pending_jobs)} jobs...")
                time.sleep(self.config.poll_interval)

        return {jid: self.jobs[jid].state for jid in job_ids if jid in self.jobs}

    def get_completed_evaluations(self) -> List[Tuple[int, Path]]:
        """
        Get list of completed evaluation IDs and directories.

        Returns
        -------
        List[Tuple[int, Path]]
            List of (eval_id, work_dir) tuples
        """
        completed = []

        for job in self.jobs.values():
            if job.state == JobState.COMPLETED and job.work_dir:
                # Verify completion marker exists
                if (job.work_dir / "COMPLETED").exists():
                    completed.append((job.eval_id, job.work_dir))

        return completed

    def cancel_job(self, job_id: str) -> bool:
        """
        Cancel a SLURM job.

        Parameters
        ----------
        job_id : str
            Job ID to cancel

        Returns
        -------
        bool
            True if cancellation was successful
        """
        try:
            subprocess.run(
                ["scancel", job_id.split('_')[0]],  # Remove array index
                check=True
            )

            if job_id in self.jobs:
                self.jobs[job_id].state = JobState.CANCELLED
                self.jobs[job_id].end_time = time.time()

            logger.info(f"Cancelled job {job_id}")
            return True

        except subprocess.CalledProcessError as e:
            logger.error(f"Failed to cancel job {job_id}: {e}")
        except FileNotFoundError:
            logger.error("scancel command not found")

        return False

    def cancel_all_jobs(self) -> int:
        """
        Cancel all tracked jobs.

        Returns
        -------
        int
            Number of jobs cancelled
        """
        cancelled = 0

        for job_id, job in self.jobs.items():
            if job.state in [JobState.PENDING, JobState.RUNNING]:
                if self.cancel_job(job_id):
                    cancelled += 1

        return cancelled

    def _generate_single_job_script(self, eval_id: int, eval_dir: Path) -> str:
        """Generate SLURM script for a single evaluation."""
        script = f"""#!/bin/bash
#SBATCH --job-name=openwq_eval_{eval_id:04d}
#SBATCH --partition={self.config.partition}
#SBATCH --time={self.config.walltime}
#SBATCH --nodes={self.config.nodes}
#SBATCH --ntasks-per-node={self.config.tasks_per_node}
#SBATCH --mem={self.config.memory}
#SBATCH --output={self.logs_dir}/eval_{eval_id:04d}_%j.out
#SBATCH --error={self.logs_dir}/eval_{eval_id:04d}_%j.err
"""

        if self.config.account:
            script += f"#SBATCH --account={self.config.account}\n"
        if self.config.qos:
            script += f"#SBATCH --qos={self.config.qos}\n"

        script += f"""
# Set working directory
WORK_DIR="{eval_dir}"
cd $WORK_DIR

echo "Starting evaluation {eval_id} at $(date)"
echo "Working directory: $WORK_DIR"

# Run model with Apptainer
apptainer exec --bind {self.config.bind_path} \\
    --env master_json=$WORK_DIR/openWQ_master.json \\
    {self.config.sif_path} \\
    ./{self.config.executable_name} {self.config.executable_args} \\
    -m {self.config.file_manager_path}

EXIT_CODE=$?

echo "Model completed with exit code: $EXIT_CODE"

# Create completion marker
if [ $EXIT_CODE -eq 0 ]; then
    touch $WORK_DIR/COMPLETED
    echo "SUCCESS" > $WORK_DIR/STATUS
else
    echo "FAILED" > $WORK_DIR/STATUS
fi

echo "Finished evaluation {eval_id} at $(date)"
exit $EXIT_CODE
"""
        return script

    def _generate_array_job_script(self,
                                   n_tasks: int,
                                   max_concurrent: int,
                                   index_file: Path) -> str:
        """Generate SLURM array job script."""
        script = f"""#!/bin/bash
#SBATCH --job-name=openwq_array
#SBATCH --array=0-{n_tasks - 1}%{max_concurrent}
#SBATCH --partition={self.config.partition}
#SBATCH --time={self.config.walltime}
#SBATCH --nodes={self.config.nodes}
#SBATCH --ntasks-per-node={self.config.tasks_per_node}
#SBATCH --mem={self.config.memory}
#SBATCH --output={self.logs_dir}/array_%A_%a.out
#SBATCH --error={self.logs_dir}/array_%A_%a.err
"""

        if self.config.account:
            script += f"#SBATCH --account={self.config.account}\n"
        if self.config.qos:
            script += f"#SBATCH --qos={self.config.qos}\n"

        script += f"""
# Read task info from index file
INDEX_FILE="{index_file}"
TASK_ID=$SLURM_ARRAY_TASK_ID

# Parse work directory from JSON (requires jq)
if command -v jq &> /dev/null; then
    WORK_DIR=$(jq -r ".[\\"$TASK_ID\\"].dir" $INDEX_FILE)
    EVAL_ID=$(jq -r ".[\\"$TASK_ID\\"].eval_id" $INDEX_FILE)
else
    # Fallback: use python
    WORK_DIR=$(python3 -c "import json; d=json.load(open('$INDEX_FILE')); print(d['$TASK_ID']['dir'])")
    EVAL_ID=$(python3 -c "import json; d=json.load(open('$INDEX_FILE')); print(d['$TASK_ID']['eval_id'])")
fi

cd $WORK_DIR

echo "Starting array task $TASK_ID (eval $EVAL_ID) at $(date)"
echo "Working directory: $WORK_DIR"

# Wait for parameters file (in case coordinator is still writing)
PARAMS_FILE="$WORK_DIR/parameters.json"
while [ ! -f "$PARAMS_FILE" ]; do
    echo "Waiting for parameters file..."
    sleep 5
done

# Run model with Apptainer
apptainer exec --bind {self.config.bind_path} \\
    --env master_json=$WORK_DIR/openWQ_master.json \\
    {self.config.sif_path} \\
    ./{self.config.executable_name} {self.config.executable_args} \\
    -m {self.config.file_manager_path}

EXIT_CODE=$?

echo "Model completed with exit code: $EXIT_CODE"

# Create completion marker
if [ $EXIT_CODE -eq 0 ]; then
    touch $WORK_DIR/COMPLETED
    echo "SUCCESS" > $WORK_DIR/STATUS
else
    echo "FAILED" > $WORK_DIR/STATUS
fi

echo "Finished array task $TASK_ID at $(date)"
exit $EXIT_CODE
"""
        return script

    def get_job_statistics(self) -> Dict:
        """
        Get statistics about submitted jobs.

        Returns
        -------
        Dict
            Job statistics
        """
        stats = {
            "total": len(self.jobs),
            "pending": 0,
            "running": 0,
            "completed": 0,
            "failed": 0,
            "cancelled": 0,
            "other": 0,
            "total_runtime": 0.0,
            "avg_runtime": 0.0,
        }

        runtimes = []

        for job in self.jobs.values():
            if job.state == JobState.PENDING:
                stats["pending"] += 1
            elif job.state == JobState.RUNNING:
                stats["running"] += 1
            elif job.state == JobState.COMPLETED:
                stats["completed"] += 1
                if job.runtime:
                    runtimes.append(job.runtime)
            elif job.state == JobState.FAILED:
                stats["failed"] += 1
            elif job.state == JobState.CANCELLED:
                stats["cancelled"] += 1
            else:
                stats["other"] += 1

        if runtimes:
            stats["total_runtime"] = sum(runtimes)
            stats["avg_runtime"] = sum(runtimes) / len(runtimes)

        return stats
