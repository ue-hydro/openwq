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
import threading
from pathlib import Path
from typing import Optional, Dict, Tuple, List
from concurrent.futures import ThreadPoolExecutor, as_completed
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
                 executable_name: str = "mizuroute_lakes_openwq_Release",
                 executable_args: str = "",
                 file_manager_path: str = None,
                 executable_full_path: str = None,
                 command_template: str = None,
                 hostmodel: str = "",
                 calibration_work_dir: str = None,
                 calibration_period: Optional[Tuple[str, str]] = None,
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
        executable_full_path : str, optional
            Full container path to executable. If provided, overrides
            default path construction from executable_name.
        command_template : str, optional
            Custom command template using placeholders:
            {eval_dir}, {exec_path}, {master_json}, {file_manager}, {args}.
            Default: "cd {eval_dir} && mpirun --allow-run-as-root -np 2 -x master_json {exec_path} {file_manager}"
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
        self.executable_full_path = executable_full_path
        self.command_template = command_template
        # Host model governs how the control file is passed to the executable:
        #   SUMMA      → `-m <fileManager>`   (master-file flag)
        #   mizuRoute  → `<control_file>`     (positional argument)
        self.hostmodel = (hostmodel or "").lower()
        # Optional (start, end) calibration window.  When set, the per-eval
        # control file is rewritten so each evaluation simulates ONLY this
        # window instead of the model's full period — keeps runtime + memory
        # low and avoids OOM kills (openWQ holds its output in RAM).
        self.calibration_period = calibration_period
        self.timeout_seconds = timeout_seconds

        # Reference size (bytes) of a COMPLETED evaluation's openWQ output,
        # used to estimate live % progress.  openWQ writes the SAME set of
        # HDF5 outputs (one record per output timestep) regardless of host
        # model, and every eval simulates the same period → the same final
        # size.  So this is host-agnostic.  We seed it up front (from a
        # persisted value or a prior completed eval on disk) so % is shown
        # from the very first chunk, and refresh it the first time an eval
        # completes in this run.
        self._work_dir = calibration_work_dir
        self._ref_file = (Path(calibration_work_dir) / ".openwq_output_reference.txt"
                          if calibration_work_dir else None)
        self._ref_output_bytes = None
        self._init_reference()

        # Parse Docker compose to get volume mapping
        self.docker_host_path = None
        self.docker_container_path = "/code"
        if docker_compose_path:
            self._parse_docker_compose()

        # ── Automatic memory guard ──────────────────────────────────────
        # openWQ's in-RAM footprint is far larger than its on-disk output
        # (~10-20×) and balloons early, so running too many evals at once can
        # exhaust the container and trigger OS OOM-kills.  We detect the
        # container memory ceiling and a per-eval RAM estimate (persisted
        # across runs, refined by live sampling) and cap concurrency so
        # n_parallel × per-eval-RAM stays under the limit — auto-reducing
        # n_parallel with a clear log line instead of dying mid-run.
        self._mem_safety = 0.85   # fraction of the limit we allow to be used
        self._mem_limit_gb = self._detect_mem_limit_gb()
        self._mem_file = (Path(calibration_work_dir) / ".openwq_mem_per_eval_gb"
                          if calibration_work_dir else None)
        self._mem_per_eval_gb = self._load_mem_estimate()
        if self._mem_limit_gb:
            logger.debug(f"Memory guard: container limit ~"
                         f"{self._mem_limit_gb:.1f} GB, per-eval estimate "
                         f"{self._mem_per_eval_gb}.")

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
        # Host→container path mapper, matching the docker-compose volume
        # mount (e.g. /Users/me/Documents -> /code).  EVERY host path handed
        # to the container (eval dir, executable, control file) must be
        # rewritten, otherwise the container can't find them.
        if self.docker_host_path:
            _hroot, _croot = self.docker_host_path, self.docker_container_path
        else:
            # Fallback: assume the user's home dir maps to /code.
            _hroot, _croot = str(Path.home()), "/code"

        def _to_container(p):
            return str(p).replace(_hroot, _croot) if p else p

        eval_dir_abs = str(eval_dir.resolve())
        container_eval_dir = _to_container(eval_dir_abs)
        container_master_json = f"{container_eval_dir}/openWQ_master.json"

        # Executable path inside the container.
        if self.executable_full_path:
            exec_path = _to_container(self.executable_full_path)
        else:
            # Default: assume it's in the standard bin directory
            exec_path = f"{self.docker_container_path}/openwq_code/6_mizuroute_cslm_openwq/route/build/openwq/openwq/bin/{self.executable_name}"

        # Control file (host-model file manager) inside the container.
        container_file_manager = _to_container(self.file_manager_path)

        # Host-model-specific invocation:
        #   SUMMA      → single process, control file behind the `-m` flag.
        #                SUMMA is NOT MPI-domain-decomposed, so launching >1
        #                rank runs duplicate simulations that collide on the
        #                shared output (rank 0 finishes, rank 1 hits STOP 1).
        #                `-s <suffix>` gives each evaluation a UNIQUE SUMMA
        #                output filename in the shared outputPath, so parallel
        #                evals don't fight over the same NetCDF/HDF5 file
        #                (which surfaces as "Permission denied" on create).
        #   mizuRoute  → MPI-parallel across reaches; control file positional.
        if self.hostmodel == "summa":
            n_ranks = 1
            # Per-eval SUMMA output dir (fresh + isolated) avoids stale-file
            # 'Permission denied' on re-create and parallel collisions.
            eval_fm = self._summa_eval_filemanager(eval_dir, container_eval_dir)
            fm_arg = eval_fm or container_file_manager
            model_arg = f"-m {fm_arg}"
        else:
            n_ranks = 2
            # Per-eval control file with the calibration window (if set).
            eval_ctrl = self._mizuroute_eval_control(eval_dir, container_eval_dir)
            model_arg = f"{eval_ctrl or container_file_manager}"

        # Build the shell command using template or default
        if self.command_template:
            shell_cmd = self.command_template.format(
                eval_dir=container_eval_dir,
                exec_path=exec_path,
                master_json=container_master_json,
                file_manager=container_file_manager,
                args=self.executable_args or ""
            )
        else:
            shell_cmd = (
                f"cd {container_eval_dir} && mpirun --allow-run-as-root "
                f"-np {n_ranks} -x master_json {exec_path} {model_arg}"
            )

        cmd = [
            "docker", "exec",
            "-e", f"master_json={container_master_json}",
            self.docker_container_name,
            "/bin/bash", "-c",
            shell_cmd
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
                oom = self._oom_message(result.returncode, result.stderr)
                if oom:
                    return False, oom
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

        # Convert host paths to their in-container equivalents (bind mount).
        def _to_container(p):
            return str(p).replace(host_path, container_path) if p else p

        eval_dir_abs = str(eval_dir.resolve())
        container_eval_dir = _to_container(eval_dir_abs)
        container_master_json = f"{container_eval_dir}/openWQ_master.json"
        container_file_manager = _to_container(self.file_manager_path)

        # Build executable path
        exec_dir = f"{container_path}/route/build/openwq/openwq/bin"

        # SUMMA → `-m <fileManager>` with a per-eval fileManager that writes
        # into this eval's own summa_out (isolated, fresh each run).
        # mizuRoute → control file positionally.
        if self.hostmodel == "summa":
            eval_fm = self._summa_eval_filemanager(eval_dir, container_eval_dir)
            model_args = ["-m", eval_fm or container_file_manager]
        else:
            eval_ctrl = self._mizuroute_eval_control(eval_dir, container_eval_dir)
            model_args = [eval_ctrl or container_file_manager]

        cmd = [
            "apptainer", "exec",
            "--bind", self.apptainer_bind_path,
            "--pwd", exec_dir,
            "--env", f"master_json={container_master_json}",
            self.apptainer_sif_path,
            f"./{self.executable_name}",
            *self.executable_args.split(),
            *model_args
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
                oom = self._oom_message(result.returncode, result.stderr)
                if oom:
                    return False, oom
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

    def _init_reference(self):
        """Seed the reference output size so % progress is available from the
        first chunk.  Tries (1) a persisted value, then (2) the largest
        completed openWQ output among existing eval directories on disk
        (left by a prior run / sensitivity stage)."""
        # 1. persisted reference
        if self._ref_file and self._ref_file.is_file():
            try:
                val = int(self._ref_file.read_text().strip())
                if val > 0:
                    self._ref_output_bytes = val
                    return
            except (ValueError, OSError):
                pass
        # 2. fallback: scan existing eval dirs for the biggest output set
        if self._work_dir:
            evals_dir = Path(self._work_dir) / "evaluations"
            best = 0
            try:
                for d in evals_dir.glob("eval_*"):
                    best = max(best, self._openwq_out_bytes(d))
            except OSError:
                pass
            if best > 0:
                self._ref_output_bytes = best

    def _save_reference(self, nbytes: int):
        """Persist the reference output size for future runs."""
        if self._ref_file and nbytes > 0:
            try:
                self._ref_file.write_text(str(int(nbytes)))
            except OSError:
                pass

    def _summa_eval_filemanager(self, eval_dir, container_eval_dir) -> str:
        """Create a per-evaluation SUMMA fileManager so each eval writes its
        SUMMA output into its OWN ``<eval_dir>/summa_out/`` (fresh every run),
        instead of a shared directory.  This avoids 'Permission denied' on
        re-creating leftover output files and any cross-eval collision when
        running in parallel.  Returns the CONTAINER path of the per-eval
        fileManager (falls back to the baseline path on any problem).

        Only ``outputPath`` is rewritten; ``settingsPath`` / ``forcingPath``
        stay pointed at the shared (read-only) inputs.  When a calibration
        window is set, ``simStartTime`` / ``simEndTime`` are also rewritten so
        the eval simulates only that (short) window.
        """
        import re
        try:
            host_out = Path(eval_dir) / "summa_out"
            host_out.mkdir(parents=True, exist_ok=True)
            text = Path(self.file_manager_path).read_text()
            container_out = f"{container_eval_dir.rstrip('/')}/summa_out/"
            text = re.sub(
                r"(outputPath\s+')[^']*(')",
                lambda m: m.group(1) + container_out + m.group(2), text)
            # Restrict the simulation to the calibration window (if set).
            if self.calibration_period:
                cstart, cend = self.calibration_period
                if cstart:
                    text = re.sub(
                        r"(simStartTime\s+')[^']*(')",
                        lambda m: m.group(1) + str(cstart) + m.group(2), text)
                if cend:
                    text = re.sub(
                        r"(simEndTime\s+')[^']*(')",
                        lambda m: m.group(1) + str(cend) + m.group(2), text)
            host_fm = Path(eval_dir) / "fileManager_eval.txt"
            host_fm.write_text(text)
            # Container path of the per-eval fileManager = the eval dir's
            # container path + the file name.
            return f"{container_eval_dir.rstrip('/')}/fileManager_eval.txt"
        except OSError as e:
            logger.warning(
                f"Could not write per-eval SUMMA fileManager ({e}); "
                f"using the shared one.")
            return None

    def _mizuroute_eval_control(self, eval_dir, container_eval_dir) -> Optional[str]:
        """Create a per-eval mizuRoute control file with ``<sim_start>`` /
        ``<sim_end>`` rewritten to the calibration window.  Returns its
        CONTAINER path, or ``None`` when no window is set / on any problem
        (caller then falls back to the baseline control file).

        Only the two date tags change; every ``<*_dir>`` / file path in the
        control file is absolute, so a copy in the eval dir resolves the same
        inputs.
        """
        if not self.calibration_period:
            return None
        import re
        try:
            text = Path(self.file_manager_path).read_text()
            cstart, cend = self.calibration_period

            def _mz(s):
                # mizuRoute expects 'YYYY-MM-DD HH:MM:SS'.
                s = str(s).strip()
                if re.match(r"^\d{4}-\d{2}-\d{2}$", s):
                    return s + " 00:00:00"
                if re.match(r"^\d{4}-\d{2}-\d{2} \d{2}:\d{2}$", s):
                    return s + ":00"
                return s

            if cstart:
                text = re.sub(r"(<sim_start>\s+)[^!\n]*",
                              lambda m: m.group(1) + _mz(cstart) + "    ", text)
            if cend:
                text = re.sub(r"(<sim_end>\s+)[^!\n]*",
                              lambda m: m.group(1) + _mz(cend) + "    ", text)
            host_ctrl = Path(eval_dir) / "control_eval.txt"
            host_ctrl.write_text(text)
            return f"{container_eval_dir.rstrip('/')}/control_eval.txt"
        except OSError as e:
            logger.warning(
                f"Could not write per-eval mizuRoute control file ({e}); "
                f"using the shared one.")
            return None

    @staticmethod
    def _oom_message(returncode, stderr) -> Optional[str]:
        """If a failure looks like an out-of-memory kill, return an
        actionable message; otherwise ``None``.  Exit 137 = 128 + signal 9
        (SIGKILL), which the Linux OOM-killer uses when the container hits its
        memory ceiling."""
        s = stderr or ""
        low = s.lower()
        if (returncode == 137 or "signal 9 (killed)" in low
                or "out of memory" in low or "oom-kill" in low
                or "oomkill" in low):
            return ("Killed by the OS (out of memory: exit 137 / signal 9). "
                    "openWQ holds its output in RAM, so a long simulation "
                    "period or too many parallel runs can exhaust container "
                    "memory. Fixes: shorten the calibration period in the "
                    "setup report, lower n_parallel, or give Docker / the "
                    "container more RAM.")
        return None

    @staticmethod
    def _openwq_out_bytes(eval_dir) -> int:
        """Total bytes of the openWQ HDF5 output written so far for an eval.

        Host-agnostic progress signal: openWQ appends one record per output
        timestep to each ``openwq_out/**/*.h5`` file, so the directory grows
        ~linearly with simulated time (same for SUMMA and mizuRoute).  Uses
        only ``stat`` (never opens the HDF5), so it can't collide with the
        model still writing them.
        """
        out = Path(eval_dir) / "openwq_out"
        total = 0
        try:
            for p in out.rglob("*.h5"):
                try:
                    total += p.stat().st_size
                except OSError:
                    pass
        except OSError:
            pass
        return total

    def _progress_str(self, eval_dir) -> str:
        """Human-readable per-eval progress: 'NN.N MB (PP%)' (or just MB when
        no reference size is known yet, i.e. before the first eval finishes)."""
        b = self._openwq_out_bytes(eval_dir)
        mb = b / 1_000_000.0
        if self._ref_output_bytes and self._ref_output_bytes > 0:
            pct = int(round(100.0 * b / self._ref_output_bytes))
            pct = max(0, min(99, pct))  # cap < 100 until it actually finishes
            return f"{mb:.1f} MB ({pct}%)"
        return f"{mb:.1f} MB"

    # ── Memory-guard helpers ────────────────────────────────────────────
    def _detect_mem_limit_gb(self) -> Optional[float]:
        """Best-effort container memory ceiling in GB (Docker only).

        Uses the container's explicit ``HostConfig.Memory`` if set, else the
        Docker VM total (``docker info``).  Returns ``None`` for Apptainer/HPC
        (memory is managed by the scheduler there) or if it can't be read.
        """
        if self.runtime != "docker" or not self.docker_container_name:
            return None
        try:
            out = subprocess.run(
                ["docker", "inspect", self.docker_container_name,
                 "--format", "{{.HostConfig.Memory}}"],
                capture_output=True, text=True, timeout=15)
            val = int((out.stdout or "0").strip() or 0)
            if val > 0:
                return val / (1024.0 ** 3)
        except Exception:
            pass
        try:
            out = subprocess.run(
                ["docker", "info", "--format", "{{.MemTotal}}"],
                capture_output=True, text=True, timeout=15)
            val = int((out.stdout or "0").strip() or 0)
            if val > 0:
                return val / (1024.0 ** 3)
        except Exception:
            pass
        return None

    @staticmethod
    def _parse_mem_to_gb(s) -> Optional[float]:
        """Parse a docker-style size ('12.54GiB', '880MiB') into GB."""
        import re
        m = re.match(r"\s*([0-9.]+)\s*([KMGT]?i?B)", str(s), re.I)
        if not m:
            return None
        v = float(m.group(1))
        factor = {"b": 1, "kib": 1024, "kb": 1e3, "mib": 1024 ** 2, "mb": 1e6,
                  "gib": 1024 ** 3, "gb": 1e9, "tib": 1024 ** 4,
                  "tb": 1e12}.get(m.group(2).lower(), 1)
        return v * factor / (1024.0 ** 3)

    def _sample_container_mem_gb(self) -> Optional[float]:
        """Current memory used by the container (GB) via ``docker stats``."""
        if self.runtime != "docker" or not self.docker_container_name:
            return None
        try:
            out = subprocess.run(
                ["docker", "stats", "--no-stream", "--format",
                 "{{.MemUsage}}", self.docker_container_name],
                capture_output=True, text=True, timeout=15)
            used = (out.stdout or "").split("/")[0].strip()  # "12.54GiB / 16GiB"
            return self._parse_mem_to_gb(used)
        except Exception:
            return None

    def _load_mem_estimate(self) -> Optional[float]:
        """Load the persisted per-eval RAM estimate (GB), if any."""
        if not getattr(self, "_mem_file", None):
            f = (Path(self._work_dir) / ".openwq_mem_per_eval_gb"
                 if self._work_dir else None)
        else:
            f = self._mem_file
        try:
            if f and f.is_file():
                return float(f.read_text().strip())
        except Exception:
            pass
        return None

    def _save_mem_estimate(self, gb: float):
        """Persist the per-eval RAM estimate so later runs cap from the start."""
        if not self._mem_file:
            return
        try:
            self._mem_file.write_text(f"{gb:.3f}")
        except Exception:
            pass

    def _effective_parallelism(self, n_requested: int, n_available: int) -> int:
        """Cap *n_requested* so ``n × per-eval-RAM`` fits under the container
        memory limit.  Returns it unchanged when the limit or the per-eval
        estimate is unknown (then the probe step measures it)."""
        n = min(max(1, int(n_requested)), max(1, int(n_available)))
        L, E = self._mem_limit_gb, self._mem_per_eval_gb
        if not L or not E or E <= 0:
            return n
        safe = max(1, int((L * self._mem_safety) // E))
        return min(n, safe)

    def run_parallel_evaluations(self,
                                 eval_configs: List[Dict],
                                 n_parallel: int) -> List[Tuple[int, bool, float, str]]:
        """
        Run multiple (independent) evaluations concurrently.

        Works for **both** runtimes:

        * **Docker** — each evaluation is a separate ``docker exec`` into the
          (already-running) shared container; because every eval has its own
          ``eval_dir`` (distinct ``openwq_in``/``openwq_out``) and the master
          JSON is passed per-exec, concurrent runs don't collide.
        * **Apptainer** — each evaluation is its own container instance.

        Uses a thread pool (the work is subprocess-bound — ``docker exec`` /
        ``apptainer exec`` — so the GIL is released while the model runs).
        With ``n_parallel <= 1`` (or a single config) it runs sequentially,
        with no pool overhead.

        Parameters
        ----------
        eval_configs : List[Dict]
            List of {"eval_id": int, "eval_dir": Path, "master_json": str}
        n_parallel : int
            Maximum number of concurrent model runs.

        Returns
        -------
        List[Tuple[int, bool, float, str]]
            List of (eval_id, success, runtime, error_msg)
        """
        results = []
        total = len(eval_configs)

        # ── Memory guard: pick a safe concurrency ───────────────────────
        n_req = max(1, int(n_parallel or 1))
        if (n_req > 1 and self.runtime == "docker" and self._mem_limit_gb
                and not self._mem_per_eval_gb):
            # No per-eval RAM estimate yet → run this chunk one-at-a-time to
            # measure peak memory safely before ever parallelising.
            n_workers = 1
            logger.info(
                f"    Memory guard: no per-eval RAM estimate yet — running "
                f"1-at-a-time to measure peak before parallelising "
                f"(container ~{self._mem_limit_gb:.1f} GB).")
        else:
            n_workers = self._effective_parallelism(n_req, total)
            if n_workers < n_req:
                logger.info(
                    f"    Memory guard: reducing parallel runs {n_req} → "
                    f"{n_workers} (~{self._mem_per_eval_gb:.1f} GB/eval, "
                    f"container ~{self._mem_limit_gb:.1f} GB, "
                    f"{int(self._mem_safety * 100)}% budget).")

        # ── Live progress feedback ──────────────────────────────────────
        # Per-eval completion lines (printed as each model finishes) plus a
        # periodic heartbeat so the user can see the run is alive even while
        # a long model evaluation is still in flight.
        _lock = threading.Lock()
        _done_ids = set()
        _active = set()              # eval_ids actually executing right now
        _mem_peak = {"per_eval": 0.0}  # GB, refined by live sampling
        _t0 = time.time()
        _stop_hb = threading.Event()

        def _heartbeat():
            # Fires only if a chunk takes a while; short chunks never trigger
            # it (no log spam), long ones get a reassuring per-model update of
            # how much openWQ output each running simulation has written.
            while not _stop_hb.wait(15):
                with _lock:
                    done = set(_done_ids)
                    active_n = len(_active)
                running = [c for c in eval_configs if c['eval_id'] not in done]
                if not running:
                    continue
                elapsed = time.time() - _t0
                logger.info(
                    f"    progress | {len(done)}/{total} done | "
                    f"{len(running)} running | elapsed {elapsed:.0f}s")
                # One model per row, for readability.
                for c in running:
                    logger.info(
                        f"        eval {c['eval_id']:<6} {self._progress_str(c['eval_dir'])}")
                # Sample container memory to refine the per-eval RAM estimate
                # (used to auto-cap parallelism on later chunks and runs).
                if self._mem_limit_gb and active_n > 0:
                    used = self._sample_container_mem_gb()
                    if used:
                        cand = used / active_n
                        if cand > _mem_peak["per_eval"]:
                            _mem_peak["per_eval"] = cand

        def _run_one(config):
            with _lock:
                _active.add(config['eval_id'])
            try:
                success, runtime, error = self.run_single_evaluation(
                    config['eval_dir'], config['master_json'], config['eval_id'])
            finally:
                with _lock:
                    _active.discard(config['eval_id'])
            out_bytes = self._openwq_out_bytes(config['eval_dir'])
            with _lock:
                _done_ids.add(config['eval_id'])
                n_done = len(_done_ids)
                # A successful eval gives the true full-run output size; keep
                # the largest seen as the reference and persist it for % est.
                if success and out_bytes > (self._ref_output_bytes or 0):
                    self._ref_output_bytes = out_bytes
                    self._save_reference(out_bytes)
            status = "done " if success else "FAILED"
            logger.info(
                f"    [{n_done}/{total}] eval {config['eval_id']} {status} "
                f"({runtime:.1f}s, {out_bytes / 1_000_000.0:.1f} MB output)")
            return (config['eval_id'], success, runtime, error)

        hb_thread = threading.Thread(target=_heartbeat, daemon=True)
        hb_thread.start()
        try:
            if n_workers <= 1 or total <= 1:
                # Sequential — no pool overhead (still logs each completion).
                for config in eval_configs:
                    results.append(_run_one(config))
            else:
                logger.info(
                    f"Running {total} evaluations with up to "
                    f"{n_workers} in parallel ({self.runtime})...")
                with ThreadPoolExecutor(max_workers=n_workers) as executor:
                    futures = {executor.submit(_run_one, c): c['eval_id']
                               for c in eval_configs}
                    for future in as_completed(futures):
                        eval_id = futures[future]
                        try:
                            results.append(future.result())
                        except Exception as e:
                            results.append((eval_id, False, 0.0, str(e)))
        finally:
            _stop_hb.set()
            # Persist / refine the per-eval RAM estimate so later chunks (and
            # future `_run.py` invocations) auto-cap parallelism from the start.
            newE = _mem_peak["per_eval"]
            if newE > 0 and (not self._mem_per_eval_gb
                             or newE > self._mem_per_eval_gb):
                self._mem_per_eval_gb = newE
                self._save_mem_estimate(newE)
                logger.info(
                    f"    Memory guard: per-eval RAM estimate now "
                    f"~{newE:.1f} GB (container ~{self._mem_limit_gb:.1f} GB).")

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
