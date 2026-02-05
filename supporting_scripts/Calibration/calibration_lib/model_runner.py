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
Model Runner Module
===================

Executes OpenWQ model in Docker or Apptainer containers.
"""

import subprocess
import os
import time
from pathlib import Path
from typing import Optional, Dict, Tuple, List
from concurrent.futures import ProcessPoolExecutor, as_completed
import logging
import shutil

logger = logging.getLogger(__name__)


class ModelRunner:
    """
    Handles model execution in containerized environments.
    """

    def __init__(self,
                 runtime: str,
                 docker_container_name: str = None,
                 docker_compose_path: str = None,
                 apptainer_sif_path: str = None,
                 apptainer_bind_path: str = None,
                 executable_name: str = "mizuroute_lakes_cslm_openwq_fast",
                 executable_args: str = "-g 1 1",
                 file_manager_path: str = None,
                 timeout_seconds: int = 7200):
        """
        Initialize model runner.

        Parameters
        ----------
        runtime : str
            "docker" or "apptainer"
        docker_container_name : str
            Name of the Docker container (e.g., "docker_openwq")
        docker_compose_path : str
            Path to docker-compose.yml
        apptainer_sif_path : str
            Path to Apptainer .sif file
        apptainer_bind_path : str
            Bind mount specification "host_path:container_path"
        executable_name : str
            Name of the executable
        executable_args : str
            Additional command-line arguments
        file_manager_path : str
            Path to mizuRoute file manager (inside container)
        timeout_seconds : int
            Maximum execution time
        """
        self.runtime = runtime
        self.docker_container_name = docker_container_name
        self.docker_compose_path = docker_compose_path
        self.apptainer_sif_path = apptainer_sif_path
        self.apptainer_bind_path = apptainer_bind_path
        self.executable_name = executable_name
        self.executable_args = executable_args
        self.file_manager_path = file_manager_path
        self.timeout_seconds = timeout_seconds

        # Parse Docker compose to get volume mapping
        self.docker_host_path = None
        self.docker_container_path = "/code"
        if docker_compose_path:
            self._parse_docker_compose()

    def _parse_docker_compose(self):
        """Parse docker-compose.yml to extract volume mapping."""
        try:
            with open(self.docker_compose_path, 'r') as f:
                content = f.read()
            # Simple parsing for "volumes: - host:container"
            import re
            match = re.search(r'volumes:\s*\n\s*-\s*([^:]+):([^:\s]+)', content)
            if match:
                host_rel = match.group(1).strip()
                container = match.group(2).strip().rstrip(':Z')
                # Convert relative path to absolute
                compose_dir = Path(self.docker_compose_path).parent
                self.docker_host_path = str((compose_dir / host_rel).resolve())
                self.docker_container_path = container
                logger.debug(f"Docker volume: {self.docker_host_path} -> {self.docker_container_path}")
        except Exception as e:
            logger.warning(f"Could not parse docker-compose.yml: {e}")

    def run_single_evaluation(self,
                              eval_dir: Path,
                              master_json_path: str,
                              eval_id: int) -> Tuple[bool, float, str]:
        """
        Run a single model evaluation.

        Parameters
        ----------
        eval_dir : Path
            Working directory for this evaluation
        master_json_path : str
            Path to master JSON file (inside container)
        eval_id : int
            Evaluation identifier

        Returns
        -------
        Tuple[bool, float, str]
            (success, runtime_seconds, error_message)
        """
        start_time = time.time()

        try:
            if self.runtime == "docker":
                success, error = self._run_docker(eval_dir, master_json_path)
            elif self.runtime == "apptainer":
                success, error = self._run_apptainer(eval_dir, master_json_path)
            else:
                return False, 0.0, f"Unknown runtime: {self.runtime}"
        except Exception as e:
            return False, time.time() - start_time, str(e)

        elapsed = time.time() - start_time

        # Save runtime info
        runtime_file = eval_dir / "runtime.txt"
        with open(runtime_file, 'w') as f:
            f.write(f"eval_id: {eval_id}\n")
            f.write(f"runtime_seconds: {elapsed:.2f}\n")
            f.write(f"success: {success}\n")
            if error:
                f.write(f"error: {error}\n")

        return success, elapsed, error

    def _run_docker(self,
                    eval_dir: Path,
                    master_json_path: str) -> Tuple[bool, str]:
        """
        Execute model in existing Docker container.

        Uses docker exec to run the model. Expects container is already running.
        """
        # Convert eval_dir to container path
        eval_dir_abs = str(eval_dir.resolve())
        if self.docker_host_path:
            container_eval_dir = eval_dir_abs.replace(
                self.docker_host_path, self.docker_container_path)
        else:
            # Fallback: assume /code mapping
            container_eval_dir = eval_dir_abs.replace(
                str(Path.home()), "/code")

        container_master_json = f"{container_eval_dir}/openWQ_master.json"

        # Build the command
        # Find executable path - assume it's in the bin directory
        exec_path = f"{self.docker_container_path}/openwq_code/6_mizuroute_cslm_openwq/route/build/openwq/openwq/bin/{self.executable_name}"

        cmd = [
            "docker", "exec",
            "-e", f"master_json={container_master_json}",
            self.docker_container_name,
            "/bin/bash", "-c",
            f"cd {container_eval_dir} && {exec_path} {self.executable_args} -m {self.file_manager_path}"
        ]

        logger.debug(f"Docker command: {' '.join(cmd)}")

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=self.timeout_seconds
            )

            # Save output
            log_file = eval_dir / "model_output.log"
            with open(log_file, 'w') as f:
                f.write("=== STDOUT ===\n")
                f.write(result.stdout)
                f.write("\n=== STDERR ===\n")
                f.write(result.stderr)

            if result.returncode != 0:
                return False, f"Exit code {result.returncode}: {result.stderr[-500:]}"

            # Check if output files were created
            output_dir = eval_dir / "openwq_out" / "HDF5"
            h5_files = list(output_dir.glob("*.h5"))
            if not h5_files:
                return False, "No HDF5 output files generated"

            return True, ""

        except subprocess.TimeoutExpired:
            return False, f"Timeout after {self.timeout_seconds} seconds"
        except Exception as e:
            return False, str(e)

    def _run_apptainer(self,
                       eval_dir: Path,
                       master_json_path: str) -> Tuple[bool, str]:
        """
        Execute model with Apptainer.

        Creates a new container instance per evaluation to enable parallelism.
        """
        # Parse bind path
        if ":" in self.apptainer_bind_path:
            host_path, container_path = self.apptainer_bind_path.split(":", 1)
        else:
            host_path = self.apptainer_bind_path
            container_path = "/code"

        # Convert eval_dir to container path
        eval_dir_abs = str(eval_dir.resolve())
        container_eval_dir = eval_dir_abs.replace(host_path, container_path)
        container_master_json = f"{container_eval_dir}/openWQ_master.json"

        # Build executable path
        exec_dir = f"{container_path}/route/build/openwq/openwq/bin"

        cmd = [
            "apptainer", "exec",
            "--bind", self.apptainer_bind_path,
            "--pwd", exec_dir,
            "--env", f"master_json={container_master_json}",
            self.apptainer_sif_path,
            f"./{self.executable_name}",
            *self.executable_args.split(),
            "-m", self.file_manager_path
        ]

        logger.debug(f"Apptainer command: {' '.join(cmd)}")

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=self.timeout_seconds
            )

            # Save output
            log_file = eval_dir / "model_output.log"
            with open(log_file, 'w') as f:
                f.write("=== STDOUT ===\n")
                f.write(result.stdout)
                f.write("\n=== STDERR ===\n")
                f.write(result.stderr)

            if result.returncode != 0:
                return False, f"Exit code {result.returncode}: {result.stderr[-500:]}"

            # Check if output files were created
            output_dir = eval_dir / "openwq_out" / "HDF5"
            h5_files = list(output_dir.glob("*.h5"))
            if not h5_files:
                return False, "No HDF5 output files generated"

            return True, ""

        except subprocess.TimeoutExpired:
            return False, f"Timeout after {self.timeout_seconds} seconds"
        except Exception as e:
            return False, str(e)

    def run_parallel_evaluations(self,
                                 eval_configs: List[Dict],
                                 n_parallel: int) -> List[Tuple[int, bool, float, str]]:
        """
        Run multiple evaluations in parallel.

        Only works properly with Apptainer (Docker uses single container).

        Parameters
        ----------
        eval_configs : List[Dict]
            List of {"eval_id": int, "eval_dir": Path, "master_json": str}
        n_parallel : int
            Number of parallel workers

        Returns
        -------
        List[Tuple[int, bool, float, str]]
            List of (eval_id, success, runtime, error_msg)
        """
        results = []

        if self.runtime == "docker":
            # Sequential execution for Docker
            logger.warning("Docker mode: running evaluations sequentially")
            for config in eval_configs:
                success, runtime, error = self.run_single_evaluation(
                    config['eval_dir'],
                    config['master_json'],
                    config['eval_id']
                )
                results.append((config['eval_id'], success, runtime, error))
        else:
            # Parallel execution for Apptainer
            with ProcessPoolExecutor(max_workers=n_parallel) as executor:
                futures = {}
                for config in eval_configs:
                    future = executor.submit(
                        self.run_single_evaluation,
                        config['eval_dir'],
                        config['master_json'],
                        config['eval_id']
                    )
                    futures[future] = config['eval_id']

                for future in as_completed(futures):
                    eval_id = futures[future]
                    try:
                        success, runtime, error = future.result()
                        results.append((eval_id, success, runtime, error))
                    except Exception as e:
                        results.append((eval_id, False, 0.0, str(e)))

        return results


class HPCJobGenerator:
    """
    Generates HPC batch job scripts for SLURM or PBS.
    """

    def __init__(self, scheduler: str = "slurm"):
        """
        Initialize job generator.

        Parameters
        ----------
        scheduler : str
            "slurm" or "pbs"
        """
        self.scheduler = scheduler

    def generate_array_job_script(self,
                                  script_path: Path,
                                  calibration_work_dir: str,
                                  apptainer_sif_path: str,
                                  apptainer_bind_path: str,
                                  executable_name: str,
                                  executable_args: str,
                                  file_manager_path: str,
                                  n_evaluations: int,
                                  partition: str = "standard",
                                  walltime: str = "04:00:00",
                                  nodes: int = 1,
                                  tasks_per_node: int = 1,
                                  memory: str = "8G",
                                  max_concurrent: int = 50) -> str:
        """
        Generate a batch job array script for HPC submission.

        Returns the path to the generated script.
        """
        if self.scheduler == "slurm":
            content = self._generate_slurm_array_script(
                calibration_work_dir=calibration_work_dir,
                apptainer_sif_path=apptainer_sif_path,
                apptainer_bind_path=apptainer_bind_path,
                executable_name=executable_name,
                executable_args=executable_args,
                file_manager_path=file_manager_path,
                n_evaluations=n_evaluations,
                partition=partition,
                walltime=walltime,
                nodes=nodes,
                tasks_per_node=tasks_per_node,
                memory=memory,
                max_concurrent=max_concurrent
            )
        elif self.scheduler == "pbs":
            content = self._generate_pbs_array_script(
                calibration_work_dir=calibration_work_dir,
                apptainer_sif_path=apptainer_sif_path,
                apptainer_bind_path=apptainer_bind_path,
                executable_name=executable_name,
                executable_args=executable_args,
                file_manager_path=file_manager_path,
                n_evaluations=n_evaluations,
                partition=partition,
                walltime=walltime,
                nodes=nodes,
                tasks_per_node=tasks_per_node,
                memory=memory,
                max_concurrent=max_concurrent
            )
        else:
            raise ValueError(f"Unknown scheduler: {self.scheduler}")

        with open(script_path, 'w') as f:
            f.write(content)

        # Make executable
        os.chmod(script_path, 0o755)

        return str(script_path)

    def _generate_slurm_array_script(self, **kwargs) -> str:
        """Generate SLURM job array script content."""
        n_evals = kwargs['n_evaluations']
        max_concurrent = kwargs['max_concurrent']

        return f'''#!/bin/bash
#SBATCH --job-name=openwq_calib
#SBATCH --array=0-{n_evals - 1}%{max_concurrent}
#SBATCH --partition={kwargs['partition']}
#SBATCH --time={kwargs['walltime']}
#SBATCH --nodes={kwargs['nodes']}
#SBATCH --ntasks-per-node={kwargs['tasks_per_node']}
#SBATCH --mem={kwargs['memory']}
#SBATCH --output={kwargs['calibration_work_dir']}/logs/eval_%a.out
#SBATCH --error={kwargs['calibration_work_dir']}/logs/eval_%a.err

# OpenWQ Calibration - SLURM Array Job
# Generated by OpenWQ Calibration Framework

# Configuration
CALIB_ROOT="{kwargs['calibration_work_dir']}"
SIF_PATH="{kwargs['apptainer_sif_path']}"
BIND_PATH="{kwargs['apptainer_bind_path']}"
EXECUTABLE="{kwargs['executable_name']}"
EXEC_ARGS="{kwargs['executable_args']}"
FILE_MANAGER="{kwargs['file_manager_path']}"

# Evaluation ID from array index
EVAL_ID=$SLURM_ARRAY_TASK_ID
EVAL_DIR=$CALIB_ROOT/evaluations/eval_$(printf "%04d" $EVAL_ID)

# Parse bind path
HOST_PATH=$(echo $BIND_PATH | cut -d':' -f1)
CONTAINER_PATH=$(echo $BIND_PATH | cut -d':' -f2)

# Convert paths
CONTAINER_EVAL_DIR="${{EVAL_DIR/$HOST_PATH/$CONTAINER_PATH}}"
MASTER_JSON="$CONTAINER_EVAL_DIR/openWQ_master.json"
EXEC_DIR="${{CONTAINER_PATH}}/route/build/openwq/openwq/bin"

# Wait for parameters file (created by coordinator)
PARAMS_FILE="$EVAL_DIR/parameters.json"
echo "Waiting for parameters file: $PARAMS_FILE"
while [ ! -f "$PARAMS_FILE" ]; do
    sleep 5
done

echo "Starting evaluation $EVAL_ID at $(date)"
echo "Working directory: $EVAL_DIR"

# Run the model
apptainer exec \\
    --bind $BIND_PATH \\
    --pwd $EXEC_DIR \\
    --env master_json=$MASTER_JSON \\
    $SIF_PATH \\
    ./$EXECUTABLE $EXEC_ARGS -m $FILE_MANAGER

EXIT_CODE=$?

# Record completion
echo "exit_code: $EXIT_CODE" >> $EVAL_DIR/runtime.txt
echo "end_time: $(date)" >> $EVAL_DIR/runtime.txt

# Signal completion
if [ $EXIT_CODE -eq 0 ]; then
    touch $EVAL_DIR/COMPLETED
else
    touch $EVAL_DIR/FAILED
fi

echo "Evaluation $EVAL_ID completed with exit code $EXIT_CODE at $(date)"
'''

    def _generate_pbs_array_script(self, **kwargs) -> str:
        """Generate PBS job array script content."""
        n_evals = kwargs['n_evaluations']
        max_concurrent = kwargs['max_concurrent']

        return f'''#!/bin/bash
#PBS -N openwq_calib
#PBS -J 0-{n_evals - 1}%{max_concurrent}
#PBS -q {kwargs['partition']}
#PBS -l walltime={kwargs['walltime']}
#PBS -l select={kwargs['nodes']}:ncpus={kwargs['tasks_per_node']}:mem={kwargs['memory']}
#PBS -o {kwargs['calibration_work_dir']}/logs/
#PBS -e {kwargs['calibration_work_dir']}/logs/

# OpenWQ Calibration - PBS Array Job
# Generated by OpenWQ Calibration Framework

# Configuration
CALIB_ROOT="{kwargs['calibration_work_dir']}"
SIF_PATH="{kwargs['apptainer_sif_path']}"
BIND_PATH="{kwargs['apptainer_bind_path']}"
EXECUTABLE="{kwargs['executable_name']}"
EXEC_ARGS="{kwargs['executable_args']}"
FILE_MANAGER="{kwargs['file_manager_path']}"

# Evaluation ID from array index
EVAL_ID=$PBS_ARRAY_INDEX
EVAL_DIR=$CALIB_ROOT/evaluations/eval_$(printf "%04d" $EVAL_ID)

# Parse bind path
HOST_PATH=$(echo $BIND_PATH | cut -d':' -f1)
CONTAINER_PATH=$(echo $BIND_PATH | cut -d':' -f2)

# Convert paths
CONTAINER_EVAL_DIR="${{EVAL_DIR/$HOST_PATH/$CONTAINER_PATH}}"
MASTER_JSON="$CONTAINER_EVAL_DIR/openWQ_master.json"
EXEC_DIR="${{CONTAINER_PATH}}/route/build/openwq/openwq/bin"

# Wait for parameters file
PARAMS_FILE="$EVAL_DIR/parameters.json"
echo "Waiting for parameters file: $PARAMS_FILE"
while [ ! -f "$PARAMS_FILE" ]; do
    sleep 5
done

echo "Starting evaluation $EVAL_ID at $(date)"

# Run the model
apptainer exec \\
    --bind $BIND_PATH \\
    --pwd $EXEC_DIR \\
    --env master_json=$MASTER_JSON \\
    $SIF_PATH \\
    ./$EXECUTABLE $EXEC_ARGS -m $FILE_MANAGER

EXIT_CODE=$?

# Signal completion
if [ $EXIT_CODE -eq 0 ]; then
    touch $EVAL_DIR/COMPLETED
else
    touch $EVAL_DIR/FAILED
fi

echo "Evaluation $EVAL_ID completed with exit code $EXIT_CODE"
'''
