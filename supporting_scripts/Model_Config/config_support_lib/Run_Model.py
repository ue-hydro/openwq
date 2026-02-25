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
Run_Model.py — Execute OpenWQ-coupled model inside Docker or Apptainer containers.

Handles volume mapping, command construction, execution, and output logging.
"""

import os
import sys
import subprocess
import time


def run_model_in_container(
        dir2save_input_files,
        container_runtime,
        executable_path,
        file_manager_path,
        mpi_np,
        docker_container_name="docker_openwq",
        apptainer_sif_path=None,
        apptainer_bind_path=None,
        timeout_seconds=7200,
):
    """Run the model inside a Docker or Apptainer container.

    Parameters:
        dir2save_input_files: Host directory where openWQ_master.json lives (CWD for model)
        container_runtime: "docker" or "apptainer"
        executable_path: Absolute path to executable INSIDE the container
        file_manager_path: Absolute path to mizuRoute control file INSIDE the container
        mpi_np: Number of MPI processes (minimum 2 for mizuRoute)
        docker_container_name: Docker container name (only for docker runtime)
        apptainer_sif_path: Path to .sif image (only for apptainer runtime)
        apptainer_bind_path: Bind mount "host:container" (only for apptainer runtime)
        timeout_seconds: Max execution time in seconds (default 2 hours)

    Returns:
        float: Elapsed time in seconds

    Raises:
        SystemExit on failure
    """
    import Gen_Input_Driver as gJSON_lib

    print("\n" + "=" * 60)
    print("RUNNING MODEL")
    print("=" * 60)

    if container_runtime == "docker":
        cmd, container_work_dir = _build_docker_cmd(
            dir2save_input_files, docker_container_name,
            executable_path, file_manager_path, mpi_np, gJSON_lib)

        print(f"  Runtime:    Docker ({docker_container_name})")
        print(f"  Work dir:   {container_work_dir}")

    elif container_runtime == "apptainer":
        cmd, container_work_dir = _build_apptainer_cmd(
            dir2save_input_files, apptainer_sif_path, apptainer_bind_path,
            executable_path, file_manager_path, mpi_np)

        print(f"  Runtime:    Apptainer ({apptainer_sif_path})")
        print(f"  Bind:       {apptainer_bind_path}")
        print(f"  Work dir:   {container_work_dir}")

    else:
        print(f"ERROR: Unknown container_runtime '{container_runtime}'. "
              f"Use 'docker' or 'apptainer'.")
        sys.exit(1)

    print(f"  Executable: {executable_path}")
    print(f"  Control:    {file_manager_path}")
    print(f"  MPI procs:  {mpi_np}")
    print("-" * 60)

    elapsed = _execute(cmd, dir2save_input_files, timeout_seconds)
    return elapsed


def _build_docker_cmd(dir2save_input_files, docker_container_name,
                      executable_path, file_manager_path, mpi_np, gJSON_lib):
    """Build the docker exec command."""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    docker_compose = os.path.normpath(
        os.path.join(script_dir, '..', '..', '..', 'containers', 'docker-compose.yml'))
    host_root, container_root = gJSON_lib._parse_docker_volume_mount(docker_compose)

    if not host_root or not container_root:
        print("ERROR: Could not determine Docker volume mapping.")
        print("       Check docker-compose.yml path and volume mounts.")
        sys.exit(1)

    abs_save_dir = os.path.abspath(dir2save_input_files)
    if not abs_save_dir.endswith('/'):
        abs_save_dir += '/'
    container_work_dir = gJSON_lib._correct_path_for_docker(
        abs_save_dir, host_root, container_root).rstrip('/')

    shell_cmd = (
        f"cd {container_work_dir} && "
        f"mpirun --allow-run-as-root -np {mpi_np} "
        f"{executable_path} {file_manager_path}"
    )

    cmd = [
        "docker", "exec",
        docker_container_name,
        "/bin/bash", "-c",
        shell_cmd
    ]
    return cmd, container_work_dir


def _build_apptainer_cmd(dir2save_input_files, apptainer_sif_path,
                         apptainer_bind_path, executable_path,
                         file_manager_path, mpi_np):
    """Build the apptainer exec command."""
    if apptainer_bind_path and ":" in apptainer_bind_path:
        host_bind, container_bind = apptainer_bind_path.split(":", 1)
    else:
        host_bind = apptainer_bind_path or ""
        container_bind = "/code"

    abs_save_dir = os.path.abspath(dir2save_input_files)
    container_work_dir = abs_save_dir.replace(host_bind, container_bind)

    shell_cmd = (
        f"cd {container_work_dir} && "
        f"mpirun --allow-run-as-root -np {mpi_np} "
        f"{executable_path} {file_manager_path}"
    )

    cmd = [
        "apptainer", "exec",
        "--bind", apptainer_bind_path,
        apptainer_sif_path,
        "/bin/bash", "-c",
        shell_cmd
    ]
    return cmd, container_work_dir


def _execute(cmd, dir2save_input_files, timeout_seconds):
    """Execute the command and handle output/errors."""
    start_time = time.time()
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout_seconds
        )
        elapsed = time.time() - start_time

        # Save output log
        log_path = os.path.join(dir2save_input_files, "model_output.log")
        with open(log_path, 'w') as f:
            f.write("=== STDOUT ===\n")
            f.write(result.stdout)
            f.write("\n=== STDERR ===\n")
            f.write(result.stderr)

        if result.returncode == 0:
            print(f"\nMODEL RUN SUCCESSFUL ({elapsed:.1f}s)")
            print(f"  Log: {log_path}")

            h5_dir = os.path.join(dir2save_input_files, "openwq_out", "HDF5")
            if os.path.isdir(h5_dir):
                h5_files = [f for f in os.listdir(h5_dir) if f.endswith('.h5')]
                print(f"  HDF5 outputs: {len(h5_files)} files in {h5_dir}")
            else:
                print(f"  WARNING: No HDF5 output directory found at {h5_dir}")
        else:
            print(f"\nMODEL RUN FAILED (exit code {result.returncode}, {elapsed:.1f}s)")
            print(f"  Log: {log_path}")
            stderr_lines = result.stderr.strip().split('\n')
            if stderr_lines:
                print("  Last error lines:")
                for line in stderr_lines[-20:]:
                    print(f"    {line}")
            sys.exit(1)

    except subprocess.TimeoutExpired:
        elapsed = time.time() - start_time
        print(f"\nMODEL RUN TIMED OUT ({timeout_seconds}s)")
        sys.exit(1)
    except FileNotFoundError as e:
        elapsed = time.time() - start_time
        print(f"\nERROR: Command not found: {e}")
        print(f"  Is the container runtime installed and accessible?")
        sys.exit(1)

    return elapsed
