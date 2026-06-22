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
Calibration Driver
==================

Main entry point for running OpenWQ calibration.
Orchestrates parameter handling, model execution, and optimization.
"""

import os
import sys
import shutil
import json
import logging
from pathlib import Path
from typing import Dict, List, Any, Optional
from datetime import datetime
import numpy as np

# Setup logging before imports
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Import calibration modules (relative imports since we're inside calibration_lib)
from .parameter_handler import ParameterHandler
from .model_runner import ModelRunner, HPCJobGenerator
from .objective_functions import ObjectiveFunction
from .checkpoint import CheckpointManager
from .optimization.dds import (DDS, DDSResult, RandomSearch,
                               ParallelDDS, ParallelRandom)
from .sensitivity.morris_screening import MorrisScreening
from .sensitivity.sobol_analysis import SobolAnalysis
from .postprocessing.results_analysis import ResultsAnalyzer


# Filename of the machine-readable snapshot the running calibration persists
# (alongside the HTML report) so a *separate* `--report` process can regenerate
# the report from the latest on-disk state while the run is still going.
_LIVE_SNAPSHOT_NAME = "live_report_snapshot.json"


# Human-readable names for the calibration scheme, shown in the per-eval
# banner so it's obvious which algorithm produced each simulation.
_ALGO_DISPLAY = {
    "DDS": "DDS sequential",
    "DDS_PARALLEL": "DDS parallel chains",
    "PARALLEL_DDS": "DDS parallel chains",
    "RANDOM": "RANDOM parallel",
}


def _print_eval_banner(eval_id, max_evals=None, best_obj=None, algorithm=None):
    """Print a prominent, colour-highlighted banner for each calibration
    simulation so its number stands out in the terminal amid the model's own
    (mizuRoute/openWQ) output. ANSI colour is used only when stdout is a TTY,
    so log files / pipes stay clean."""
    try:
        use_color = bool(sys.stdout.isatty())
    except Exception:
        use_color = False
    BOLD = "\033[1;96m" if use_color else ""   # bold bright cyan
    DIM  = "\033[0;90m" if use_color else ""   # grey
    RST  = "\033[0m"    if use_color else ""
    total = f" / {max_evals}" if max_evals else ""
    scheme = ""
    if algorithm:
        scheme = f" [{_ALGO_DISPLAY.get(algorithm, algorithm)}]"
    best = ""
    if best_obj is not None and best_obj < 1e9:
        best = f"{DIM}      best so far: {best_obj:.4g}{RST}"
    bar = "═" * 66
    print(f"\n{BOLD}{bar}{RST}", flush=True)
    print(f"{BOLD}  ▶  SIMULATION {eval_id}{total}{scheme}{RST}{best}", flush=True)
    print(f"{BOLD}{bar}{RST}", flush=True)


def _json_default(o):
    """json.dump default: make numpy scalars / arrays serialisable."""
    try:
        import numpy as _np
        if isinstance(o, _np.integer):
            return int(o)
        if isinstance(o, _np.floating):
            return float(o)
        if isinstance(o, _np.ndarray):
            return o.tolist()
    except Exception:
        pass
    return str(o)


def _write_live_snapshot(results_dir, calibration_results,
                         calibration_settings, in_progress):
    """Persist a JSON snapshot of the current results so a separate
    ``--report`` invocation can rebuild the HTML report from disk mid-run.

    The HTML report itself is throttled; this snapshot is the authoritative,
    machine-readable state (calibration history + settings + in_progress flag).
    Best-effort: never raise into the calibration loop.
    """
    try:
        results_dir = Path(results_dir)
        results_dir.mkdir(parents=True, exist_ok=True)
        snap = {
            "in_progress": bool(in_progress),
            "timestamp": datetime.now().isoformat(),
            "calibration_results": calibration_results,
            "calibration_settings": calibration_settings,
        }
        with open(results_dir / _LIVE_SNAPSHOT_NAME, "w") as f:
            json.dump(snap, f, indent=2, default=_json_default)
    except Exception as _e:
        logger.debug(f"Live snapshot skipped: {_e}")


def _count_sensitivity_evals(work_dir) -> int:
    """Count completed sensitivity model runs by their eval folders.

    Sensitivity samples use eval ids 10000+i (see _run_sensitivity_analysis),
    so their working dirs are evaluations/eval_1XXXX (id >= 10000).  Lets a
    separate ``--report`` process show Morris/Sobol progress while the
    screening is still running (sensitivity_results.json only appears at the
    very end).
    """
    try:
        evdir = Path(work_dir) / "evaluations"
        if not evdir.is_dir():
            return 0
        n = 0
        for p in evdir.glob("eval_*"):
            try:
                if int(p.name.split("_")[-1]) >= 10000:
                    n += 1
            except (ValueError, IndexError):
                continue
        return n
    except Exception:
        return 0


def _validate_observation_csv(path: str, obs_source: str) -> str:
    """Make sure the observation path is a calibration-ready CSV.

    For ``observation_data_source = "grqa"`` the model config gives a path
    to the raw GRQA database *directory* (per-species ``*_GRQA.csv`` files).
    That cannot be fed straight to the objective function, which expects a
    single prepared CSV (``datetime, reach_id, species, value, units,
    source, [is_primary]``).  Fail here with actionable guidance instead of
    letting ``pd.read_csv`` raise a cryptic ``IsADirectoryError`` deep in
    the objective function.
    """
    if obs_source == "skip" or not path:
        return path
    if os.path.isdir(path):
        raise ValueError(
            "Observation path is a directory, not a prepared calibration "
            f"CSV:\n    {path}\n\n"
            "The GRQA database has not been clipped/prepared for this basin "
            "yet.\n"
            "Run your MODEL-CONFIG script once first — it clips GRQA to the "
            "basin and writes:\n"
            "    <dir2save_input_files>/openwq_in/grqa_clipped_data/"
            "grqa_clipped_observations.csv\n"
            "The calibration then reshapes that automatically (no "
            "re-download).\n\n"
            "Alternatively, set in your model config:\n"
            "      observation_data_source = \"user_csv\"\n"
            "      user_observation_csv    = \"<your_observations.csv>\"\n"
            "  with columns: datetime, reach_id, species, value, units, "
            "source, [is_primary]."
        )
    return path


def _maybe_clean_work_dir(work_dir: Path, resume: bool = False,
                          clean=None) -> None:
    """Offer to clear stale evaluation folders / old results before a run.

    A previous run (especially one with MORE evaluations) leaves ``eval_*``
    folders behind that don't belong to the new calibration.  This detects
    them and, on an interactive terminal, asks whether to delete them for a
    clean run.  Never cleans when resuming.

    ``clean`` can force the decision (True/False) for non-interactive use;
    when None it prompts (interactive) or keeps + warns (non-interactive).
    """
    if resume:
        return
    evals = work_dir / "evaluations"
    stale = sorted(evals.glob("eval_*")) if evals.is_dir() else []
    if not stale:
        return

    n = len(stale)
    decision = clean
    if decision is None:
        if sys.stdin.isatty():
            print(f"\n  The calibration work dir already contains {n} "
                  f"evaluation folder(s) from a previous run:")
            print(f"      {evals}")
            print("  These may be stale (e.g. left over from a run with more "
                  "evaluations) and don't belong to this calibration.")
            try:
                resp = input("  Delete them (and old results) for a clean "
                             "run? [y/N]: ").strip().lower()
            except EOFError:
                resp = ""
            decision = resp in ("y", "yes")
        else:
            decision = False
            logger.warning(
                f"{n} stale evaluation folder(s) found in {evals} "
                f"(non-interactive run — left in place; pass "
                f"clean_work_dir=True to auto-clean).")

    if decision:
        # Before deleting, preserve a % progress reference: the largest
        # openWQ output among the eval dirs about to be removed.  Otherwise a
        # clean run loses the only reference (no eval has necessarily
        # completed yet to persist one) and the live % can't be shown.
        try:
            best = max((ModelRunner._openwq_out_bytes(d) for d in stale),
                       default=0)
            ref_file = work_dir / ".openwq_output_reference.txt"
            prev = 0
            if ref_file.is_file():
                try:
                    prev = int(ref_file.read_text().strip())
                except (ValueError, OSError):
                    prev = 0
            if best > prev:
                ref_file.write_text(str(int(best)))
        except (OSError, ValueError):
            pass

        shutil.rmtree(evals, ignore_errors=True)
        evals.mkdir(parents=True, exist_ok=True)
        # Clear old results (plots / history) too so runs don't mix.
        results = work_dir / "results"
        if results.is_dir():
            shutil.rmtree(results, ignore_errors=True)
        results.mkdir(parents=True, exist_ok=True)
        logger.info(
            f"Cleaned {n} stale evaluation folder(s) and old results for a "
            f"fresh run.")
    else:
        logger.info("Keeping existing folders (stale evaluations from a "
                    "previous run may remain).")


def run_calibration(
        # Paths
        calibration_work_dir: str,
        base_model_config_dir: str = None,
        observation_data_path: str = None,
        test_case_dir: str = None,

        # Integrated config (preferred — replaces paths + container args)
        model_config: Optional[Dict[str, Any]] = None,

        # Container runtime
        container_runtime: str = None,
        docker_container_name: str = None,
        docker_compose_path: str = None,
        apptainer_sif_path: str = None,
        apptainer_bind_path: str = None,
        executable_name: str = "mizuroute_lakes_openwq_Release",
        executable_args: str = "",
        file_manager_path: str = None,
        executable_full_path: str = None,
        command_template: str = None,

        # Calibration parameters
        calibration_parameters: List[Dict] = None,

        # Calibration settings
        algorithm: str = "DDS",
        max_evaluations: int = 500,
        n_parallel: int = 1,
        objective_function: str = "KGE",
        objective_weights: Dict[str, float] = None,
        calibration_targets: Dict = None,
        random_seed: int = None,

        # Sensitivity analysis
        run_sensitivity_first: bool = False,
        sensitivity_method: str = "morris",
        sensitivity_morris_trajectories: int = 10,
        sensitivity_morris_levels: int = 4,
        sensitivity_sobol_samples: int = 1024,
        sensitivity_threshold: float = 0.1,

        # HPC settings
        hpc_enabled: bool = False,
        hpc_scheduler: str = "slurm",
        hpc_partition: str = "standard",
        hpc_walltime: str = "24:00:00",
        hpc_nodes: int = 1,
        hpc_tasks_per_node: int = 4,
        hpc_memory: str = "8G",
        hpc_max_concurrent_jobs: int = 50,

        # Control
        resume: bool = False,
        **kwargs
) -> Dict:
    """
    Run complete calibration workflow.

    Returns
    -------
    Dict
        Calibration results including best parameters and diagnostics
    """

    start_time = datetime.now()

    # Setup
    work_dir = Path(calibration_work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    (work_dir / "results").mkdir(exist_ok=True)
    (work_dir / "logs").mkdir(exist_ok=True)

    # Setup file logging
    log_file = work_dir / "calibration.log"
    file_handler = logging.FileHandler(log_file)
    file_handler.setFormatter(logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    ))
    logging.getLogger().addHandler(file_handler)

    # Offer to clear stale evaluation folders / old results from a previous
    # run (so a smaller new run doesn't inherit leftover eval_* folders).
    _maybe_clean_work_dir(work_dir, resume=resume,
                          clean=kwargs.get("clean_work_dir"))

    # ── Resolve settings from model_config (integrated mode) ──
    if model_config is not None:
        from . import config_integration
        _container = config_integration.get_container_config(model_config)
        _obs = config_integration.get_observation_config(model_config)

        # Calibration override for openWQ debug mode.  The model config may
        # set run_mode_debug=True (handy when configuring/testing), but it
        # makes openWQ write 5 extra diagnostic HDF5 files per
        # species/compartment (chemistry/transport/ewf/ic/ss) — a big
        # speed / memory / output drag across many evaluations.  Default OFF;
        # the user can re-enable it via the setup-report checkbox.  Mutating
        # the shared model_config here propagates to every per-eval config
        # (generate_config_for_eval deep-copies this dict).
        _rmd = kwargs.get("run_mode_debug")
        if _rmd is not None:
            model_config["run_mode_debug"] = bool(_rmd)
            logger.info(
                f"openWQ debug mode (run_mode_debug) overridden to "
                f"{bool(_rmd)} for calibration evaluations")

        # Fill in any blanks from model_config
        container_runtime = container_runtime or "docker"
        executable_full_path = executable_full_path or _container.get("executable_path", "")
        file_manager_path = file_manager_path or _container.get("file_manager_path", "")

        # Observation data path — may be supplied directly or derived from
        # the model config the user already set up.  Honours all sources:
        #   grqa (folder or "auto" download) and user_csv.  The paths come
        #   straight from the model config; nothing is re-downloaded.
        if not observation_data_path:
            obs_source = _obs.get("source", "skip")
            if obs_source in ("grqa", "user_csv"):
                _log_adapter = lambda *a, **k: logger.info(
                    " ".join(str(x) for x in a)) if a else None
                observation_data_path = config_integration.prepare_calibration_observations_csv(
                    model_config, str(work_dir), log=_log_adapter)
                if not observation_data_path:
                    # Fall back to the raw config path so the guard below
                    # raises a clear, actionable error.
                    observation_data_path = (
                        _obs.get("grqa_local_data_path", "")
                        if obs_source == "grqa"
                        else _obs.get("user_observation_csv", ""))

        # Make sure we end up with a real CSV (not the raw GRQA directory).
        observation_data_path = _validate_observation_csv(
            observation_data_path, _obs.get("source", "skip"))

        # Preflight: for mizuRoute + Copernicus SS, warn early if the basin
        # shapefile's mapping_key doesn't intersect the reach segIds (the
        # lumped-vs-delineated trap → all-zero SS loads across every eval).
        try:
            _ss_warn = config_integration.validate_ss_reach_mapping(
                model_config, log=logger.debug)
            if _ss_warn:
                for _ln in _ss_warn.splitlines():
                    logger.warning(_ln)
        except Exception:
            pass  # never let a preflight check block the run

    logger.info("=" * 60)
    logger.info("OPENWQ CALIBRATION FRAMEWORK")
    logger.info("=" * 60)
    logger.info(f"Start time: {start_time}")
    logger.info(f"Working directory: {work_dir}")
    logger.info(f"Parameters to calibrate: {len(calibration_parameters)}")
    logger.info(f"Algorithm: {algorithm}")
    logger.info(f"Max evaluations: {max_evaluations}")
    logger.info(f"Container runtime: {container_runtime}")
    if model_config is not None:
        logger.info("Mode: config-integrated (Gen_Input_Driver per eval)")
    else:
        logger.info("Mode: legacy (copy from test_case_dir)")

    # Validate inputs
    if not calibration_parameters:
        raise ValueError("calibration_parameters cannot be empty")
    if not calibration_targets or "species" not in calibration_targets:
        raise ValueError("calibration_targets must include 'species'")

    # Initialize components
    logger.info("Initializing components...")

    param_handler = ParameterHandler(
        calibration_work_dir=calibration_work_dir,
        model_config=model_config,
        base_model_config_dir=base_model_config_dir,
        test_case_dir=test_case_dir,
        running_on_docker=(container_runtime == "docker"),
        calibration_period=kwargs.get("calibration_period"),
    )

    model_runner = ModelRunner(
        runtime=container_runtime,
        docker_container_name=docker_container_name,
        docker_compose_path=docker_compose_path,
        apptainer_sif_path=apptainer_sif_path,
        apptainer_bind_path=apptainer_bind_path,
        executable_name=executable_name,
        executable_args=executable_args,
        file_manager_path=file_manager_path,
        executable_full_path=executable_full_path,
        command_template=command_template,
        hostmodel=(model_config.get("hostmodel", "") if model_config else ""),
        calibration_work_dir=calibration_work_dir,
        # Per-eval calibration window (start, end) | None.  When set, each
        # eval's control file is rewritten to simulate only this window —
        # keeps runtime + memory low and avoids OOM kills.
        calibration_period=kwargs.get("calibration_period"),
        # Total planned evaluations — drives the cumulative-elapsed + ETA
        # progress lines printed after each eval's spinner finishes.
        total_evaluations=max_evaluations,
    )

    # H5 reader path for objective function
    if base_model_config_dir:
        h5_reader_path = str(Path(base_model_config_dir).parent / "2_Read_Outputs" / "hdf5_support_lib")
    else:
        # Derive from this file's location
        _this_dir = Path(__file__).resolve().parent
        h5_reader_path = str(_this_dir.parent.parent / "2_Read_Outputs" / "hdf5_support_lib")

    # Get temporal resolution settings from kwargs
    temporal_resolution = kwargs.get("temporal_resolution", "native")
    aggregation_method = kwargs.get("aggregation_method", "mean")

    logger.info(f"Temporal resolution: {temporal_resolution}")
    logger.info(f"Aggregation method: {aggregation_method}")

    # Resolve host-aware spatial mapping (mirrors Gen_Report convention).
    # Re-import defensively in case the obs-path branch above didn't fire
    # (it's gated on `model_config is not None`).
    from . import config_integration as _cint
    _spatial = _cint.get_spatial_mapping(model_config) \
        if model_config is not None else {
            "hostmodel": "mizuroute", "h5_mapping_key": "reachID"}

    obj_func = ObjectiveFunction(
        observation_path=observation_data_path,
        target_species=calibration_targets["species"],
        target_reaches=calibration_targets.get("reach_ids", "all"),
        compartments=calibration_targets.get("compartments", ["RIVER_NETWORK_REACHES"]),
        metric=objective_function,
        weights=objective_weights,
        h5_reader_path=h5_reader_path,
        temporal_resolution=temporal_resolution,
        aggregation_method=aggregation_method,
        hostmodel=_spatial.get("hostmodel", "mizuroute"),
        h5_mapping_key=_spatial.get("h5_mapping_key"),
        use_primary_only=kwargs.get("use_primary_only", True),
        zone_select=kwargs.get("zone_select"),
    )

    checkpoint_mgr = CheckpointManager(work_dir / "checkpoints")

    # Extract parameter info
    param_names = [p["name"] for p in calibration_parameters]
    bounds = [p["bounds"] for p in calibration_parameters]
    transforms = [p.get("transform", "linear") for p in calibration_parameters]
    initial_values = np.array([p["initial"] for p in calibration_parameters])

    # Transform bounds and initial values to optimization space
    bounds_opt = []
    initial_opt = []
    for i, (b, t) in enumerate(zip(bounds, transforms)):
        if t == "log":
            bounds_opt.append((np.log10(b[0]), np.log10(b[1])))
            initial_opt.append(np.log10(initial_values[i]))
        else:
            bounds_opt.append(b)
            initial_opt.append(initial_values[i])
    initial_opt = np.array(initial_opt)

    logger.info(f"Parameters: {param_names}")
    logger.info(f"Transforms: {transforms}")

    # =========================================================================
    # Sensitivity Analysis (if requested)
    # =========================================================================

    # Workflow mode (from the setup-report slider):
    #   "sensitivity"  → identify influential parameters only (no calibration)
    #   "both"         → sensitivity first, then calibration
    #   "calibration"  → calibration only (no sensitivity)
    # Back-compat: when calibration_mode isn't supplied, derive it from the
    # legacy run_sensitivity_first flag.
    calibration_mode = kwargs.get("calibration_mode")
    if calibration_mode not in ("sensitivity", "both", "calibration"):
        calibration_mode = "both" if run_sensitivity_first else "calibration"
    _do_sensitivity = calibration_mode in ("sensitivity", "both")
    _do_calibration = calibration_mode in ("both", "calibration")
    logger.info(f"Workflow mode: {calibration_mode} "
                f"(sensitivity={_do_sensitivity}, calibration={_do_calibration})")

    # ── Pre-flight ──────────────────────────────────────────────────────
    # Verify the container image + model executable exist BEFORE launching
    # hundreds of evaluations.  Without this, a misconfigured executable path
    # (e.g. a stale .sbatch that never exported OWQ_EXEC_PATH, so the runner
    # falls back to the remapped model-config path) makes EVERY eval die in
    # <1 s with no output — a whole run wasted on cryptic failures.  Fail fast
    # and clearly instead.
    _pf_ok, _pf_msg = model_runner.preflight()
    if not _pf_ok:
        logger.error("=" * 60)
        logger.error("PRE-FLIGHT CHECK FAILED - aborting before any model run.")
        for _ln in _pf_msg.splitlines():
            logger.error("  " + _ln)
        logger.error("=" * 60)
        raise RuntimeError("Pre-flight check failed - " + _pf_msg.splitlines()[0])

    # ── Pre-flight: observation availability ────────────────────────────
    # Even with a valid executable + obs file, a calibration is futile if no
    # observation overlaps BOTH the target reaches/HRUs AND the simulated time
    # window.  Without this guard every evaluation silently returns the failure
    # penalty (1e10) — the user only discovers the temporal mismatch after a
    # whole run (hours).  Resolve the window the runs will actually simulate
    # (explicit calibration window if set, else the model's declared
    # <sim_start>/<sim_end>) and check ONCE, up front.
    _cp = kwargs.get("calibration_period")
    if _cp and len(_cp) == 2 and _cp[0] and _cp[1]:
        _win_start, _win_end = _cp[0], _cp[1]
    else:
        _win_start = _win_end = None
        if model_config is not None:
            try:
                _mp = _cint.get_model_sim_period(model_config, log=logger.debug)
                if _mp:
                    _win_start, _win_end = _mp["sim_start"], _mp["sim_end"]
            except Exception:
                pass  # fall through to a window-less (spatial-only) check
    try:
        _obs_chk = obj_func.check_observation_availability(_win_start, _win_end)
    except Exception as _e:  # never let the guard itself crash the run
        logger.warning(f"Observation-availability check skipped: {_e}")
        _obs_chk = {"ok": True}
    if not _obs_chk.get("ok", True):
        _msg = _obs_chk.get("message", "")
        logger.error("=" * 60)
        logger.error("PRE-FLIGHT CHECK FAILED - no usable observations for "
                     "the performance metric.")
        for _ln in _msg.splitlines():
            logger.error("  " + _ln)
        logger.error("=" * 60)
        # Print to the terminal too, so it's seen even when console logging is
        # quiet (the file handler above always captures it as well).
        _banner = "=" * 72
        print("\n" + _banner, file=sys.stderr)
        print("OPENWQ CALIBRATION ABORTED - no usable observations", file=sys.stderr)
        print(_banner, file=sys.stderr)
        print(_msg, file=sys.stderr)
        print(_banner + "\n", file=sys.stderr)
        raise RuntimeError(
            "No usable observations - " + _obs_chk.get("headline", "see log"))

    # ── Pre-flight: Copernicus land-cover cache covers the simulation years ──
    # The dynamic source/sink loads are built per simulation year from
    # lulc_areas_all.csv.  A stale single-year cache silently maps EVERY sim
    # year onto that one year (the "all → 1993" proxy trap, with land-cover
    # change over time ignored).  Verify coverage and, interactively, offer to
    # load the missing years before launching the run.
    if model_config is not None:
        try:
            _mp_full = _cint.get_model_sim_period(model_config, log=logger.debug)
            _lulc = _cint.get_lulc_cache_status(
                model_config,
                _mp_full.get("sim_start") if _mp_full else None,
                _mp_full.get("sim_end") if _mp_full else None,
                log=logger.debug)
        except Exception as _e:
            logger.warning(f"LULC cache check skipped: {_e}")
            _lulc = None
        if _lulc and not _lulc.get("ok", True):
            _cy = _lulc["cache_years"]
            _av = _lulc["available_years"]
            _bn = "=" * 72
            _expl = (
                f"{_bn}\n"
                "COPERNICUS LAND-COVER CACHE IS INCOMPLETE\n"
                f"{_bn}\n"
                f"  Cache file              : {_lulc['cache_path']}\n"
                f"  Years currently in cache: "
                f"{_cy if _cy else 'none'}\n"
                f"  Years available (source): "
                f"{(str(_av[0]) + '-' + str(_av[-1])) if _av else 'none'}\n"
                f"  Missing years           : {_lulc['missing_years']}\n\n"
                "Because the cache is missing these years, the dynamic source/sink\n"
                "loads map EVERY simulation year onto the nearest cached year "
                f"(e.g. all\n  years → "
                f"{_cy[0] if _cy else 'the single cached year'}), so land-cover "
                "change over time is ignored.\n" + _bn
            )
            print("\n" + _expl, file=sys.stderr)
            for _ln in _expl.splitlines():
                logger.warning(_ln)
            _do = False
            if getattr(sys.stdin, "isatty", lambda: False)():
                try:
                    _ans = input("Load the missing year(s) now and update "
                                 "lulc_areas_all.csv? [y/N] ").strip().lower()
                    _do = _ans in ("y", "yes")
                except EOFError:
                    _do = False
            else:
                print("(Non-interactive session — not prompting. Re-run "
                      "interactively to load the missing years, or proceed "
                      "with the proxy mapping.)", file=sys.stderr)
            if _do:
                _ok = False
                try:
                    _ok = _cint.regenerate_lulc_cache(
                        model_config,
                        _mp_full.get("sim_start") if _mp_full else None,
                        _mp_full.get("sim_end") if _mp_full else None,
                        log=lambda *a, **k: print(*a, file=sys.stderr))
                except Exception as _e:
                    logger.error(f"LULC regeneration failed: {_e}")
                if _ok:
                    # Drop stale per-eval copies so they re-seed from the fixed
                    # baseline cache on the next evaluation.
                    import glob as _glob
                    for _f in _glob.glob(os.path.join(
                            str(work_dir), "evaluations", "eval_*", "openwq_in",
                            "ss_copernicus_files", "lulc_areas_all.csv")):
                        try:
                            os.remove(_f)
                        except OSError:
                            pass
                    print("✓ lulc_areas_all.csv updated with all available "
                          "years. Continuing with the run.\n", file=sys.stderr)
                else:
                    print("⚠ Could not update the cache — continuing with "
                          "the current (proxy) land-cover data.\n", file=sys.stderr)
            else:
                print("Continuing — each simulation year will use the nearest "
                      "available land-cover year.\n", file=sys.stderr)

    # ── Interim (partial) results report ────────────────────────────────
    # Lets the user open the <stem>_results_report.html WHILE the run is still
    # going — it renders whatever is available so far (influential parameters
    # and/or partial calibration history) with an "in progress" banner.  The
    # generated run-script writes the final (complete) report at the end.
    report_stem = kwargs.get("report_stem")
    _live_history = []          # appended after each optimization evaluation
    _last_report_t = [0.0]

    def _partial_report(in_progress=True, throttle=0.0):
        """Write the user-facing results report from results available so far."""
        if not report_stem or model_config is None:
            return
        now = datetime.now().timestamp()
        if throttle and (now - _last_report_t[0]) < throttle:
            return
        _last_report_t[0] = now
        try:
            from . import Gen_Calibration_Results_Report as _GRR
            _hist = list(_live_history)
            _objs = [h["objective"] for h in _hist] or [float("nan")]
            _best = min(_objs)
            _bestp = next((h["parameters"] for h in _hist
                           if h["objective"] == _best), {})
            _cr = {
                "calibration_mode": calibration_mode,
                "calibration_ran": bool(_hist),
                "best_params": _bestp,
                "best_objective": _best,
                "n_evaluations": len(_hist),
                "history": _hist,
                "converged": False,
                "convergence_reason": "in progress",
            }
            _sr = None
            _sp = work_dir / "results" / "sensitivity_results.json"
            if _sp.exists():
                with open(_sp) as _f:
                    _sr = json.load(_f)
            _cs = {
                "algorithm": algorithm, "max_evaluations": max_evaluations,
                "objective_function": objective_function,
                "temporal_resolution": kwargs.get("temporal_resolution", "native"),
                "aggregation_method": kwargs.get("aggregation_method", "mean"),
                "calibration_targets": calibration_targets or {},
                "random_seed": random_seed,
                "calibration_mode": calibration_mode,
                "use_primary_only": kwargs.get("use_primary_only", True),
                "zone_select": kwargs.get("zone_select"),
            }
            try:
                _md = obj_func.get_matched_data()
            except Exception:
                _md = None
            # Persist matched data so a separate `--report` process can render
            # the performance/obs-map sections from disk mid-run.
            _results_dir = work_dir / "results"
            try:
                if _md is not None and not _md.empty:
                    _results_dir.mkdir(parents=True, exist_ok=True)
                    _md.to_csv(_results_dir / "matched_data.csv", index=False)
            except Exception:
                pass
            # Persist the FULL simulated series (all objective species at the
            # target reaches, incl. species with no obs) so the report can plot
            # a simulated-only time series for them.
            try:
                _sd = obj_func.get_simulated_data()
                if _sd is not None and not _sd.empty:
                    _results_dir.mkdir(parents=True, exist_ok=True)
                    _sd.to_csv(_results_dir / "simulated_data.csv", index=False)
            except Exception:
                pass
            # Persist the machine-readable snapshot (authoritative state for
            # on-demand regeneration via the run-script's --report flag).
            _write_live_snapshot(_results_dir, _cr, _cs, in_progress)
            _GRR.generate_results_report(
                output_dir=str(work_dir), model_config=model_config,
                calibration_parameters=calibration_parameters,
                calibration_settings=_cs, calibration_results=_cr,
                sensitivity_results=_sr, matched_data=_md,
                report_stem=report_stem, in_progress=in_progress)
        except Exception as _e:
            logger.debug(f"Interim report skipped: {_e}")

    sensitive_params = None
    if _do_sensitivity:
        logger.info("=" * 60)
        logger.info("RUNNING SENSITIVITY ANALYSIS")
        logger.info("=" * 60)

        sa_result = _run_sensitivity_analysis(
            method=sensitivity_method,
            param_handler=param_handler,
            model_runner=model_runner,
            obj_func=obj_func,
            calibration_parameters=calibration_parameters,
            bounds_opt=bounds_opt,
            transforms=transforms,
            morris_trajectories=sensitivity_morris_trajectories,
            morris_levels=sensitivity_morris_levels,
            sobol_samples=sensitivity_sobol_samples,
            threshold=sensitivity_threshold,
            work_dir=work_dir,
            random_seed=random_seed,
            n_parallel=n_parallel
        )

        sensitive_params = sa_result.influential_params
        sa_result.print_summary()
        sa_result.save(work_dir / "results" / "sensitivity_results.json")

        logger.info(f"Influential parameters: {sensitive_params}")
        logger.info(f"Non-influential parameters will be fixed at initial values")

        # TODO: Filter calibration_parameters to only include sensitive ones

        # Interim report so influential parameters are viewable immediately
        # (and while calibration runs, in 'both' mode).
        _partial_report(in_progress=_do_calibration)

    # =========================================================================
    # Sensitivity-only workflow → skip calibration entirely and return.
    # =========================================================================
    if not _do_calibration:
        logger.info("=" * 60)
        logger.info("WORKFLOW: influential-parameters only — calibration skipped")
        logger.info("=" * 60)
        results_dir = work_dir / "results"
        results_dir.mkdir(parents=True, exist_ok=True)
        sens_only = {
            "calibration_mode": calibration_mode,
            "calibration_ran": False,
            "best_params": {},
            "best_objective": float("nan"),
            "n_evaluations": 0,
            "converged": False,
            "convergence_reason": "Calibration not run (sensitivity-only workflow)",
            "history": [],
            "sensitive_params": sensitive_params,
            "results_dir": str(results_dir),
        }
        sa_path = results_dir / "sensitivity_results.json"
        if sa_path.exists():
            with open(sa_path) as f:
                sens_only["sensitivity_results"] = json.load(f)
        return sens_only

    # =========================================================================
    # Check for resume
    # =========================================================================

    start_eval = 0
    history = []

    if resume and checkpoint_mgr.checkpoint_exists():
        state = checkpoint_mgr.load_state()
        if state:
            start_eval = state['current_eval'] + 1
            history = state.get('history', [])
            initial_opt = state['best_params_array']
            logger.info(f"Resuming from evaluation {start_eval}")

    # =========================================================================
    # Create objective wrapper
    # =========================================================================

    eval_counter = [start_eval]

    def _append_all_simulated(eval_id):
        """Append this eval's simulated series (target reach(es), objective
        species) to results/all_simulated.csv so the report can draw a light-
        gray 'all simulations' overlay behind the best-fit time series.
        Best-effort — never breaks an evaluation."""
        try:
            _sim = obj_func.get_simulated_data()
            if _sim is None or _sim.empty:
                return
            _sim = _sim.copy()
            _sim["eval_id"] = int(eval_id)
            _asd = work_dir / "results" / "all_simulated.csv"
            _asd.parent.mkdir(parents=True, exist_ok=True)
            _sim.to_csv(_asd, mode="a", index=False, header=not _asd.exists())
        except Exception:
            pass

    def objective_wrapper(params_opt: np.ndarray) -> float:
        """Wrapper that handles parameter transformation and model execution."""
        eval_id = eval_counter[0]
        eval_counter[0] += 1

        # Prominent, highlighted banner so this simulation's number is easy to
        # spot in the terminal amid the model's own output.
        _best_so_far = min((h["objective"] for h in _live_history), default=None)
        _print_eval_banner(eval_id, max_evaluations, _best_so_far, algorithm)
        logger.info(f"--- Evaluation {eval_id} ---")

        # Transform from optimization space to real space
        params_real = np.array([
            ParameterHandler.transform_to_real(p, t)
            for p, t in zip(params_opt, transforms)
        ])

        logger.debug(f"Parameters (real): {dict(zip(param_names, params_real))}")

        # Setup working directory.  Pass parameters+values so generation-time
        # params (e.g. dynamic SS climate-response params) are baked into the
        # config before the SS JSON is generated.
        eval_dir = param_handler.setup_working_directory(
            eval_id, calibration_parameters, params_real)

        # Apply parameters (post-hoc edits to the generated config files)
        param_handler.apply_parameters(eval_dir, calibration_parameters, params_real)

        # Save parameters for reference
        params_file = eval_dir / "parameters.json"
        with open(params_file, 'w') as f:
            json.dump({name: float(val) for name, val in zip(param_names, params_real)}, f, indent=2)

        # Run model
        master_json = str(eval_dir / "openWQ_master.json")
        success, runtime, error = model_runner.run_single_evaluation(
            eval_dir, master_json, eval_id
        )

        if not success:
            logger.warning(f"Evaluation {eval_id} failed: {error}")
            obj_val = 1e10  # Penalty for failed run
        else:
            # Compute objective
            output_dir = eval_dir / "openwq_out"
            obj_val = obj_func.compute(output_dir)

        # Save objective
        obj_file = eval_dir / "objective.txt"
        with open(obj_file, 'w') as f:
            f.write(f"objective: {obj_val}\n")
            f.write(f"success: {success}\n")

        logger.info(f"Evaluation {eval_id}: obj = {obj_val:.6f}, runtime = {runtime:.1f}s")

        # Accumulate live history and refresh the partial report (throttled),
        # so the user can open the results report mid-run and see progress.
        _live_history.append({
            "eval_id": int(eval_id),
            "objective": float(obj_val),
            "parameters": {n: float(v) for n, v in zip(param_names, params_real)},
        })
        if success:
            _append_all_simulated(eval_id)
        _partial_report(in_progress=True, throttle=20.0)

        return obj_val

    def batch_objective(points_opt):
        """Evaluate a LIST of optimization-space parameter vectors, running up
        to ``n_parallel`` model evaluations concurrently.  Same per-eval work
        as ``objective_wrapper`` (setup → run → score → live history) but
        batched + parallel; used by the parallel algorithms (DDS_PARALLEL and
        RANDOM).  Returns objective values in the SAME order as the inputs."""
        n = len(points_opt)
        objs = [1e10] * n
        n_par = max(1, int(n_parallel or 1))
        # For DDS parallel chains the optimizer hands us exactly one candidate
        # per chain per call, in chain order → the position IS the chain id.
        # Tag it so mid-run partial reports can also draw one path per chain.
        _is_chains = algorithm in ("DDS_PARALLEL", "PARALLEL_DDS")

        for chunk_start in range(0, n, n_par):
            chunk = list(range(chunk_start, min(chunk_start + n_par, n)))

            # ── Setup each eval in the chunk (sequential, fast) ──
            eval_configs = []
            real_by_idx = {}
            for idx in chunk:
                params_opt = np.asarray(points_opt[idx])
                eval_id = eval_counter[0]
                eval_counter[0] += 1
                params_real = np.array([
                    ParameterHandler.transform_to_real(p, t)
                    for p, t in zip(params_opt, transforms)
                ])
                real_by_idx[idx] = (eval_id, params_real)
                eval_dir = param_handler.setup_working_directory(
                    eval_id, calibration_parameters, params_real)
                param_handler.apply_parameters(
                    eval_dir, calibration_parameters, params_real)
                with open(eval_dir / "parameters.json", 'w') as f:
                    json.dump({nm: float(v) for nm, v
                               in zip(param_names, params_real)}, f, indent=2)
                eval_configs.append({
                    "eval_id": eval_id,
                    "eval_dir": eval_dir,
                    "master_json": str(eval_dir / "openWQ_master.json"),
                    "_idx": idx,
                })

            # Banner spanning the whole concurrent chunk so the simulation
            # numbers stand out even when several run at once.
            first_id = eval_configs[0]["eval_id"]
            last_id = eval_configs[-1]["eval_id"]
            _best_so_far = min((h["objective"] for h in _live_history),
                               default=None)
            if first_id == last_id:
                _print_eval_banner(first_id, max_evaluations, _best_so_far,
                                   algorithm)
            else:
                _print_eval_banner(f"{first_id}-{last_id}  ({len(eval_configs)} "
                                   f"in parallel)", max_evaluations,
                                   _best_so_far, algorithm)
            logger.info(f"--- Evaluations {first_id}-{last_id} "
                        f"(n_parallel={n_par}) ---")

            # ── Run the chunk (parallel when n_par > 1) ──
            run_results = model_runner.run_parallel_evaluations(
                eval_configs, n_par)
            by_id = {eid: (ok, rt, err)
                     for (eid, ok, rt, err) in run_results}

            # ── Objective for each (sequential; reads HDF5 output) ──
            for cfg in eval_configs:
                idx = cfg["_idx"]
                eval_id, params_real = real_by_idx[idx]
                ok, rt, err = by_id.get(cfg["eval_id"], (False, 0.0, "no result"))
                if ok:
                    obj_val = obj_func.compute(cfg["eval_dir"] / "openwq_out")
                    _append_all_simulated(eval_id)
                else:
                    logger.warning(f"Evaluation {eval_id} failed: {err}")
                    obj_val = 1e10
                objs[idx] = obj_val

                with open(cfg["eval_dir"] / "objective.txt", 'w') as f:
                    f.write(f"objective: {obj_val}\n")
                    f.write(f"success: {ok}\n")
                logger.info(f"Evaluation {eval_id}: obj = {obj_val:.6f}, "
                            f"runtime = {rt:.1f}s")

                _live_history.append({
                    "eval_id": int(eval_id),
                    "objective": float(obj_val),
                    "parameters": {nm: float(v) for nm, v
                                   in zip(param_names, params_real)},
                    "chain": (int(idx) if _is_chains else None),
                })
            _partial_report(in_progress=True, throttle=20.0)

        return objs

    # =========================================================================
    # Run Optimization
    # =========================================================================

    logger.info("=" * 60)
    logger.info(f"RUNNING {algorithm} OPTIMIZATION")
    logger.info("=" * 60)

    # Checkpoint callback
    def checkpoint_callback(eval_num: int, obj_val: float, params: np.ndarray):
        # Transform params back to real space for storage
        params_real = np.array([
            ParameterHandler.transform_to_real(p, t)
            for p, t in zip(params, transforms)
        ])
        checkpoint_mgr.save_state(
            current_eval=eval_num,
            max_evals=max_evaluations,
            best_objective=obj_val,
            best_params=params_real,
            param_names=param_names,
            history=history,
            algorithm=algorithm
        )

    if algorithm == "DDS":
        optimizer = DDS(
            n_params=len(calibration_parameters),
            bounds=bounds_opt,
            param_names=param_names,
            max_evals=max_evaluations - start_eval,
            seed=random_seed
        )
        result = optimizer.optimize(
            objective_wrapper,
            initial_point=initial_opt,
            callback=checkpoint_callback
        )
    elif algorithm == "RANDOM":
        # RANDOM (parallel): independent samples evaluated up to n_parallel at
        # a time via the batch objective.  With n_parallel == 1 this is the
        # plain sequential random baseline.
        optimizer = ParallelRandom(
            n_params=len(calibration_parameters),
            bounds=bounds_opt,
            param_names=param_names,
            max_evals=max_evaluations - start_eval,
            batch=n_parallel,
            seed=random_seed
        )
        result = optimizer.optimize(
            batch_objective,
            callback=checkpoint_callback
        )
    elif algorithm in ("DDS_PARALLEL", "PARALLEL_DDS"):
        # Parallel DDS chains: n_parallel independent DDS chains advanced
        # together, one candidate per chain per round evaluated concurrently —
        # keeps DDS's sample efficiency while using n_parallel model runs.
        optimizer = ParallelDDS(
            n_params=len(calibration_parameters),
            bounds=bounds_opt,
            param_names=param_names,
            max_evals=max_evaluations - start_eval,
            n_chains=n_parallel,
            seed=random_seed
        )
        result = optimizer.optimize(
            batch_objective,
            initial_point=initial_opt,
            callback=checkpoint_callback
        )
    else:
        raise ValueError(f"Unknown algorithm: {algorithm}")

    # =========================================================================
    # Post-processing
    # =========================================================================

    logger.info("=" * 60)
    logger.info("POST-PROCESSING RESULTS")
    logger.info("=" * 60)

    # Transform best parameters back to real space
    best_params_real = {
        name: ParameterHandler.transform_to_real(val, transforms[i])
        for i, (name, val) in enumerate(zip(param_names, result.best_params))
    }

    # Save results
    results_dir = work_dir / "results"
    result.save(results_dir / "optimization_results.json")

    # Save best parameters
    with open(results_dir / "best_parameters.json", 'w') as f:
        json.dump(best_params_real, f, indent=2)

    # Save calibration history in the format expected by ResultsAnalyzer.
    # Entries are (eval_num, objective, params_array) for DDS/RANDOM, or
    # (eval_num, objective, params_array, chain_id) for ParallelDDS — the
    # trailing chain id (None otherwise) lets the report draw one path/chain.
    calibration_history = []
    for entry in result.history:
        eval_num, objective, params_array = entry[0], entry[1], entry[2]
        chain = entry[3] if len(entry) > 3 else None
        hist_entry = {
            "eval_id": int(eval_num),
            "objective": float(objective),
            "parameters": {
                name: ParameterHandler.transform_to_real(val, transforms[i])
                for i, (name, val) in enumerate(
                    zip(param_names, params_array))
            },
            "timestamp": 0,
            "chain": (int(chain) if chain is not None else None),
        }
        calibration_history.append(hist_entry)

    with open(results_dir / "calibration_history.json", 'w') as f:
        json.dump(calibration_history, f, indent=2)

    # Save parameter definitions (for basin report aggregation)
    _param_defs_serializable = []
    for pdef in calibration_parameters:
        d = dict(pdef)
        # Convert tuples to lists for JSON
        if "bounds" in d and isinstance(d["bounds"], tuple):
            d["bounds"] = list(d["bounds"])
        # Convert Path-like objects to strings
        if "path" in d and not isinstance(d["path"], (str, list, dict)):
            d["path"] = str(d["path"])
        _param_defs_serializable.append(d)
    with open(results_dir / "parameter_definitions.json", 'w') as f:
        json.dump(_param_defs_serializable, f, indent=2)

    # Resolve shapefile directory from base_model_config_dir
    # Convention: ../shapefiles/ relative to variant dir
    _shapefile_dir = kwargs.get("shapefile_dir", "")
    if not _shapefile_dir and base_model_config_dir:
        _base_cfg = Path(base_model_config_dir)
        if _base_cfg.exists():
            _candidate = _base_cfg.parent / "shapefiles"
            if _candidate.is_dir():
                _shapefile_dir = str(_candidate)

    # Generate plots and HTML report
    try:
        analyzer = ResultsAnalyzer(
            results_dir,
            parameter_definitions=calibration_parameters
        )
        analyzer.plot_convergence(output_path=results_dir / "convergence.png")
        analyzer.plot_parameter_evolution(
            param_names=param_names,
            output_path=results_dir / "parameter_evolution.png"
        )
        analyzer.plot_parameter_correlations(
            output_path=results_dir / "parameter_correlations.png"
        )

        # Generate performance plots if matched data is available
        matched_data = obj_func.get_matched_data()
        if matched_data is not None and not matched_data.empty:
            logger.info("Generating performance plots...")
            analyzer.generate_performance_plots(
                matched_data=matched_data,
                output_dir=results_dir / "performance",
                temporal_resolution=temporal_resolution
            )
            # Save matched data for basin report aggregation
            matched_data.to_csv(
                results_dir / "matched_data.csv", index=False)

        # Save the full simulated series (all objective species at the target
        # reaches, including those with no observations) so the report can plot
        # a simulated-only time series for objective species that lack obs.
        try:
            _sim_all = obj_func.get_simulated_data()
            if _sim_all is not None and not _sim_all.empty:
                _sim_all.to_csv(results_dir / "simulated_data.csv", index=False)
        except Exception as _e:
            logger.debug(f"Could not persist simulated_data.csv: {_e}")

        # Auto-validation: re-run the best parameters over the validation
        # period (held-out data) and write the combined calibration+validation
        # dataset for the report's "Best Sim (Calibration & Validation)"
        # section.  Idempotent: skipped if no validation_period is given OR the
        # validation data already exists — so a --resume run produces it only if
        # a prior (interrupted) run never got to it, and never redundantly.
        _vp = kwargs.get("validation_period")
        if _vp:
            _vfile = results_dir / "validation_matched_data.csv"
            if _vfile.exists() and _vfile.stat().st_size > 0:
                logger.info("Validation data already present (%s) — skipping "
                            "the validation re-run.", _vfile.name)
            else:
                try:
                    logger.info("Re-running best parameters over the validation "
                                "period %s for the calibration+validation "
                                "report (runs the model once; may take a "
                                "while)...", tuple(_vp))
                    _run_validation_and_combine(
                        work_dir=work_dir, results_dir=results_dir,
                        validation_period=_vp,
                        best_params_real=best_params_real,
                        param_names=param_names,
                        calibration_parameters=calibration_parameters,
                        param_handler=param_handler, model_runner=model_runner,
                        obj_func=obj_func)
                except Exception as _e:
                    logger.warning(f"Validation re-run skipped: {_e}")

        # Generate comprehensive HTML report
        logger.info("Generating HTML calibration report...")

        html_report_path = analyzer.generate_html_report(
            output_dir=results_dir,
            matched_data=matched_data if matched_data is not None else None,
            config_summary={
                "basin_id": kwargs.get("basin_id", ""),
                "variant": kwargs.get("variant", ""),
                "algorithm": algorithm,
                "max_evaluations": max_evaluations,
                "objective_function": objective_function,
                "temporal_resolution": temporal_resolution,
                "aggregation_method": kwargs.get("aggregation_method", "mean"),
                "n_parameters": len(calibration_parameters),
                "calibration_targets": kwargs.get("calibration_targets", {}),
                "random_seed": random_seed,
                "container_runtime": kwargs.get("container_runtime", ""),
            },
            temporal_resolution=temporal_resolution,
            shapefile_dir=_shapefile_dir,
            observation_data_path=observation_data_path,
        )
        logger.info(f"HTML report saved to: {html_report_path}")
    except Exception as e:
        logger.warning(f"Could not generate plots/report: {e}")
        import traceback
        logger.debug(traceback.format_exc())

    # =========================================================================
    # Auto-generate per-basin report (aggregates all completed variants)
    # =========================================================================
    _basin_id = kwargs.get("basin_id", "")
    if _basin_id:
        try:
            from .postprocessing.results_analysis import generate_basin_report

            logger.info("Generating per-basin report...")

            # Discover sibling variant workspaces
            # Convention: workspace_{basin_id}_{variant} in the same parent
            ws_parent = work_dir.parent
            prefix = f"workspace_{_basin_id}_"
            variant_results = {}

            for ws_dir in sorted(ws_parent.iterdir()):
                if not ws_dir.is_dir():
                    continue
                if not ws_dir.name.startswith(prefix):
                    continue

                # Extract variant label from directory name
                v_label = ws_dir.name[len(prefix):]
                v_results_dir = ws_dir / "results"
                v_history = v_results_dir / "calibration_history.json"

                if not v_history.exists():
                    continue  # not yet completed

                # Load parameter definitions if saved
                v_pdefs_path = v_results_dir / "parameter_definitions.json"
                v_pdefs = []
                if v_pdefs_path.exists():
                    with open(v_pdefs_path) as f:
                        v_pdefs = json.load(f)
                    # Convert bounds lists back to tuples
                    for pd in v_pdefs:
                        if "bounds" in pd and isinstance(pd["bounds"], list):
                            pd["bounds"] = tuple(pd["bounds"])

                variant_results[v_label] = {
                    "results_dir": str(v_results_dir),
                    "calibration_parameters": v_pdefs,
                }

            if variant_results:
                # Resolve paths for map data
                _obs_path = observation_data_path or ""
                basin_report_dir = ws_parent / "basin_reports"

                basin_report_path = generate_basin_report(
                    basin_id=_basin_id,
                    variant_results=variant_results,
                    output_dir=str(basin_report_dir),
                    shapefile_dir=_shapefile_dir,
                    observation_data_path=_obs_path,
                )
                if basin_report_path:
                    logger.info(f"Basin report saved to: "
                                f"{basin_report_path}")
            else:
                logger.info("No completed variant workspaces found "
                            "for basin report")
        except Exception as e:
            logger.warning(f"Could not generate basin report: {e}")
            import traceback
            logger.debug(traceback.format_exc())

    # Compute runtime
    end_time = datetime.now()
    runtime_hours = (end_time - start_time).total_seconds() / 3600.0

    # Find convergence eval (evaluation that found the best)
    convergence_eval = 0
    if calibration_history:
        best_obj_val = result.best_objective
        for entry in calibration_history:
            if abs(entry["objective"] - best_obj_val) < 1e-12:
                convergence_eval = entry["eval_id"]
                break

    # Summary
    logger.info("=" * 60)
    logger.info("CALIBRATION COMPLETE")
    logger.info("=" * 60)
    logger.info(f"Best objective: {result.best_objective:.6f}")
    logger.info(f"Total evaluations: {result.n_evaluations + start_eval}")
    logger.info(f"Runtime: {runtime_hours:.2f} hours")
    logger.info(f"Best parameters:")
    for name, val in best_params_real.items():
        logger.info(f"  {name}: {val}")
    logger.info(f"Results saved to: {results_dir}")

    # Build the comprehensive results dict
    # (compatible with Gen_Calibration_Results_Report.generate_results_report)
    calibration_results = {
        "calibration_mode": calibration_mode,
        "calibration_ran": True,
        "best_params": best_params_real,
        "best_objective": result.best_objective,
        "n_evaluations": result.n_evaluations + start_eval,
        "converged": result.converged,
        "convergence_reason": result.convergence_reason,
        "convergence_eval": convergence_eval,
        "runtime_hours": runtime_hours,
        "history": calibration_history,
        "sensitive_params": sensitive_params,
        "results_dir": str(results_dir),
    }

    # Load sensitivity results if available
    sa_results_path = results_dir / "sensitivity_results.json"
    if sa_results_path.exists():
        with open(sa_results_path) as f:
            calibration_results["sensitivity_results"] = json.load(f)

    # Attach matched data / performance metrics for report generation
    try:
        matched_data = obj_func.get_matched_data()
        if matched_data is not None and not matched_data.empty:
            calibration_results["matched_data"] = matched_data
            # Compute per-species metrics
            from .objective_functions import compute_all_metrics
            perf_metrics = compute_all_metrics(
                matched_data, calibration_targets.get("species", [])
            )
            calibration_results["performance_metrics"] = perf_metrics
    except Exception as e:
        logger.debug(f"Could not attach performance data: {e}")

    # Final snapshot marking the run complete (in_progress=False).  A later
    # `--report` invocation reads this so it knows the run finished and drops
    # the "in progress" banner.  We strip the in-memory DataFrame first (the
    # matched data is already persisted to results/matched_data.csv).
    _snap_cr = {k: v for k, v in calibration_results.items()
                if k not in ("matched_data", "performance_metrics")}
    _snap_cs = {
        "algorithm": algorithm, "max_evaluations": max_evaluations,
        "objective_function": objective_function,
        "temporal_resolution": kwargs.get("temporal_resolution", "native"),
        "aggregation_method": kwargs.get("aggregation_method", "mean"),
        "calibration_targets": calibration_targets or {},
        "random_seed": random_seed,
        "calibration_mode": calibration_mode,
        "use_primary_only": kwargs.get("use_primary_only", True),
        "zone_select": kwargs.get("zone_select"),
    }
    _write_live_snapshot(results_dir, _snap_cr, _snap_cs, in_progress=False)

    return calibration_results


# Failure penalty an evaluation gets when it cannot be scored (must match
# objective_functions._FAIL value of 1e10).
_REPORT_PENALTY = 1e10


def _reconstruct_history_from_evals(work_dir):
    """Rebuild the calibration history from the per-evaluation folders.

    Each completed calibration evaluation writes ``evaluations/eval_<id>/
    objective.txt`` (``objective: <loss>`` + ``success: <bool>``) and
    ``parameters.json`` (flat ``{name: value}``) AS IT FINISHES.  These are the
    ground truth and survive an interrupted run, an empty/zero-length history
    pickle, or a stale ``results/*.json`` that a re-run or an HPC fetch left
    behind.  Returns a list sorted by eval id:
    ``[{eval_id, objective, parameters, success}, ...]`` (empty if none).
    """
    import re as _re
    out = []
    evdir = Path(work_dir) / "evaluations"
    if not evdir.is_dir():
        return out
    for d in sorted(evdir.glob("eval_*")):
        m = _re.search(r"eval_(\d+)", d.name)
        if not m:
            continue
        objp = d / "objective.txt"
        if not objp.exists():       # sensitivity-screening evals have none
            continue
        try:
            txt = objp.read_text()
        except Exception:
            continue
        om = _re.search(r"objective:\s*([-+0-9.eE]+)", txt)
        if not om:
            continue
        try:
            obj = float(om.group(1))
        except ValueError:
            continue
        params = {}
        pp = d / "parameters.json"
        if pp.exists():
            try:
                with open(pp) as f:
                    pj = json.load(f)
                if isinstance(pj, dict):
                    params = pj.get("parameters", pj)
            except Exception:
                params = {}
        out.append({"eval_id": int(m.group(1)), "objective": obj,
                    "parameters": params,
                    "success": ("success: True" in txt)})
    out.sort(key=lambda e: e["eval_id"])
    return out


def _load_checkpoint_best(work_dir):
    """Return ``(best_objective, best_params, current_eval)`` from the optimizer
    checkpoint (``checkpoints/calibration_state.json``), which is updated every
    evaluation and therefore holds the authoritative best even when the run was
    interrupted before the final ``results/`` files were written.  Returns
    ``(None, None, None)`` if unavailable."""
    p = Path(work_dir) / "checkpoints" / "calibration_state.json"
    if not p.exists():
        return None, None, None
    try:
        with open(p) as f:
            st = json.load(f)
    except Exception:
        return None, None, None
    bo = st.get("best_objective")
    bo = bo if isinstance(bo, (int, float)) and bo == bo else None
    ce = st.get("current_eval")
    ce = ce if isinstance(ce, int) else None
    return bo, (st.get("best_params") or None), ce


def _n_valid_evals(history):
    """Count history entries with a finite, non-penalty objective."""
    n = 0
    for e in history or []:
        o = e.get("objective") if isinstance(e, dict) else None
        if isinstance(o, (int, float)) and o == o and o < _REPORT_PENALTY:
            n += 1
    return n


def _best_eval_dir(work_dir):
    """Path to the GLOBAL-best evaluation folder — the lowest successful,
    finite objective with HDF5 output — scanned across every
    ``evaluations/eval_*`` folder (they accumulate across the original run and
    every ``--resume``).  Returns the folder path (str) or ``None``."""
    import glob as _glob
    import re as _re
    evdir = Path(work_dir) / "evaluations"
    if not evdir.is_dir():
        return None
    best = None
    for d in _glob.glob(str(evdir / "eval_*")):
        f = os.path.join(d, "objective.txt")
        if not os.path.isfile(f):
            continue
        try:
            t = open(f).read()
        except Exception:
            continue
        if "success: True" not in t:
            continue
        m = _re.search(r"objective:\s*([-0-9eE.+]+)", t)
        if not m:
            continue
        try:
            o = float(m.group(1))
        except ValueError:
            continue
        if o >= _REPORT_PENALTY:
            continue
        if not _glob.glob(os.path.join(d, "openwq_out", "HDF5", "*.h5")):
            continue
        if best is None or o < best[1]:
            best = (d, o)
    return best[0] if best else None


def _best_params_from_eval_dir(eval_dir, param_names):
    """Read the best eval's ``parameters.json`` (real-space values keyed by
    name) and return them as an array in ``param_names`` order, or ``None``."""
    import json as _json
    if not eval_dir:
        return None
    pj = os.path.join(eval_dir, "parameters.json")
    if not os.path.isfile(pj):
        return None
    try:
        d = _json.load(open(pj))
        return np.array([float(d[n]) for n in param_names])
    except Exception:
        return None


def _run_validation_and_combine(*, work_dir, results_dir, validation_period,
                                best_params_real, param_names,
                                calibration_parameters, param_handler,
                                model_runner, obj_func):
    """Re-run the best parameters over the VALIDATION period, then write
    ``results/validation_matched_data.csv`` = calibration-best + validation
    obs/sim pairs, period-tagged, for the report's combined cal/val section.

    Heavy (one model run over the validation window) — called once when a
    calibration completes.  Best-effort: any failure just leaves the section's
    placeholder in place.  The validation eval folder gets a distinct id and no
    ``objective.txt``, so it never affects the best-eval scan."""
    import pandas as pd
    _VAL_ID = 999999

    def _mk_of():
        # Fresh ObjectiveFunction cloned from obj_func's config.  Using a fresh
        # instance per period (one compute each) keeps get_matched_data() bound
        # to THAT period's pairs (no cross-period best-tracking).  The h5 reader
        # path is already on sys.path from obj_func's init.
        return ObjectiveFunction(
            observation_path=obj_func.observation_path,
            target_species=obj_func.target_species,
            target_reaches=obj_func.target_reaches,
            compartments=obj_func.compartments,
            metric=obj_func.metric, weights=obj_func.weights,
            temporal_resolution=obj_func.temporal_resolution,
            aggregation_method=obj_func.aggregation_method,
            hostmodel=obj_func.hostmodel,
            h5_mapping_key=obj_func.h5_mapping_key,
            use_primary_only=getattr(obj_func, "use_primary_only", True),
            zone_select=getattr(obj_func, "zone_select", None))

    # Calibration half + best parameters from the GLOBAL-best eval folder
    # (ground truth across all runs/resumes) — independent of obj_func's
    # in-memory state, so this is correct even on a resume that adds no evals.
    best_dir = _best_eval_dir(work_dir)
    cal_md, bp_arr = None, None
    if best_dir:
        try:
            _cof = _mk_of()
            _cof.compute(Path(best_dir) / "openwq_out")
            cal_md = _cof.get_matched_data()
        except Exception as _e:
            logger.debug(f"Validation: cal half from {best_dir} failed: {_e}")
        bp_arr = _best_params_from_eval_dir(best_dir, param_names)
    if cal_md is None or cal_md.empty:
        try:
            cal_md = obj_func.get_matched_data()
        except Exception:
            cal_md = None
    if bp_arr is None:
        try:
            bp_arr = np.array([best_params_real[n] for n in param_names])
        except Exception:
            logger.warning("Validation: could not resolve best parameters; "
                           "skipping validation re-run.")
            return

    # Run the best parameters over the VALIDATION window.
    eval_dir = param_handler.setup_working_directory(
        _VAL_ID, calibration_parameters, bp_arr)
    param_handler.apply_parameters(
        eval_dir, calibration_parameters, bp_arr)
    master_json = str(Path(eval_dir) / "openWQ_master.json")
    # Temporarily point the runner's window at the validation period.
    _saved = model_runner.calibration_period
    model_runner.calibration_period = tuple(validation_period)
    try:
        ok, _rt, err = model_runner.run_single_evaluation(
            eval_dir, master_json, _VAL_ID)
    finally:
        model_runner.calibration_period = _saved
    if not ok:
        logger.warning(f"Validation re-run did not complete: {err}")
        return
    val_md = None
    try:
        _vof = _mk_of()
        _vof.compute(Path(eval_dir) / "openwq_out")
        val_md = _vof.get_matched_data()
    except Exception as _e:
        logger.warning(f"Validation: matching the validation run failed: {_e}")
    parts = []
    if cal_md is not None and not cal_md.empty:
        cal_md = cal_md.copy()
        cal_md["period"] = "calibration"
        parts.append(cal_md)
    if val_md is not None and not val_md.empty:
        val_md = val_md.copy()
        val_md["period"] = "validation"
        parts.append(val_md)
    if not parts:
        logger.warning("Validation: no matched obs/sim pairs in either period.")
        return
    combined = pd.concat(parts, ignore_index=True)
    results_dir = Path(results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)
    combined.to_csv(results_dir / "validation_matched_data.csv", index=False)
    logger.info("Wrote validation_matched_data.csv "
                "(best fit over calibration + validation).")


def _rebuild_best_fit_from_evals(work_dir, calibration_settings, model_config):
    """Find the GLOBAL-best evaluation by scanning every ``evaluations/eval_*``
    folder (the ground truth — they accumulate across the original run AND
    every ``--resume``, since eval numbering continues), then recompute that
    eval's matched + simulated series via the ObjectiveFunction.

    This makes the report's "best fit" independent of the (churning, possibly
    fragmented or stale) ``results/*.csv`` caches: no matter how many times a
    calibration is interrupted and resumed, the report always reflects the true
    best simulation for ANY basin.  Returns ``(matched_df, sim_df, best_label)``
    or ``(None, None, None)`` when nothing usable is found.
    """
    work_dir = Path(work_dir)
    # Lowest successful, finite objective with HDF5 output across ALL eval
    # folders (handles original run + every --resume).
    best_dir = _best_eval_dir(work_dir)
    if not best_dir:
        return None, None, None

    cs = calibration_settings or {}
    tg = cs.get("calibration_targets") or {}
    species = tg.get("species") or []
    if not species:
        return None, None, None
    # Observation file lives in the calibration work dir.
    obs_path = work_dir / "calibration_observations.csv"
    if not obs_path.is_file():
        return None, None, None
    # h5 reader path (parent of the hdf5_support_lib package).
    _this = Path(__file__).resolve().parent
    h5_reader_path = str(_this.parent.parent / "2_Read_Outputs" / "hdf5_support_lib")
    hostmodel, mapkey = "mizuroute", None
    try:
        from . import config_integration as _cint
        _sm = (_cint.get_spatial_mapping(model_config) if model_config else {})
        hostmodel = _sm.get("hostmodel", "mizuroute")
        mapkey = _sm.get("h5_mapping_key")
    except Exception:
        pass
    try:
        _of = ObjectiveFunction(
            observation_path=str(obs_path),
            target_species=species,
            target_reaches=tg.get("reach_ids", "all"),
            compartments=tg.get("compartments") or ["RIVER_NETWORK_REACHES"],
            metric=cs.get("objective_function", "KGE"),
            weights=cs.get("objective_weights"),
            h5_reader_path=h5_reader_path,
            temporal_resolution=cs.get("temporal_resolution", "native"),
            aggregation_method=cs.get("aggregation_method", "mean"),
            hostmodel=hostmodel, h5_mapping_key=mapkey,
            use_primary_only=cs.get("use_primary_only", True),
            zone_select=cs.get("zone_select"))
        _of.compute(Path(best_dir) / "openwq_out")
        return (_of.get_matched_data(), _of.get_simulated_data(),
                os.path.basename(best_dir))
    except Exception as _e:
        logger.debug(f"best-fit rebuild from eval folders failed: {_e}")
        return None, None, None


def regenerate_results_report(*, calibration_work_dir, model_config,
                              calibration_parameters, calibration_settings,
                              report_stem):
    """Rebuild the calibration results report from whatever state is on disk.

    Safe to call WHILE a calibration / sensitivity run is still going (e.g.
    from a second terminal via the run-script's ``--report`` flag): it reads
    the persisted live snapshot + result files and regenerates the HTML,
    keeping an "in progress" banner until the run has finished.  Returns the
    report path (or None).
    """
    from . import Gen_Calibration_Results_Report as _GRR
    work_dir = Path(calibration_work_dir)
    results_dir = work_dir / "results"
    cs = dict(calibration_settings or {})
    mode = cs.get("calibration_mode") or "both"

    cr = {
        "calibration_mode": mode,
        "calibration_ran": False,
        "best_params": {},
        "best_objective": float("nan"),
        "n_evaluations": 0,
        "history": [],
        "converged": False,
        "convergence_reason": "in progress",
    }
    in_progress = True

    # 1) Live snapshot — authoritative interim/final state written by the run.
    snap_path = results_dir / _LIVE_SNAPSHOT_NAME
    if snap_path.exists():
        try:
            with open(snap_path) as f:
                snap = json.load(f)
            if isinstance(snap.get("calibration_results"), dict):
                cr.update(snap["calibration_results"])
            if isinstance(snap.get("calibration_settings"), dict):
                for k, v in snap["calibration_settings"].items():
                    cs.setdefault(k, v)   # caller's settings win
            in_progress = bool(snap.get("in_progress", True))
        except Exception as _e:
            logger.debug(f"Snapshot read failed: {_e}")

    # 2) Final result files override once the run has completed.
    hist_path = results_dir / "calibration_history.json"
    if hist_path.exists():
        try:
            with open(hist_path) as f:
                hist = json.load(f)
            if hist:
                cr["history"] = hist
                cr["n_evaluations"] = len(hist)
                cr["calibration_ran"] = True
        except Exception:
            pass
    bp_path = results_dir / "best_parameters.json"
    if bp_path.exists():
        try:
            with open(bp_path) as f:
                cr["best_params"] = json.load(f)
        except Exception:
            pass
    opt_path = results_dir / "optimization_results.json"
    if opt_path.exists():
        try:
            with open(opt_path) as f:
                _opt = json.load(f)
            cr["best_objective"] = _opt.get("best_objective", cr["best_objective"])
            cr["converged"] = _opt.get("converged", cr["converged"])
            cr["convergence_reason"] = _opt.get("convergence_reason",
                                                cr["convergence_reason"])
        except Exception:
            pass

    # ── Robustness: trust the per-eval folders + checkpoint over stale aggregates ──
    # results/*.json can be stale, empty, or left over from a DIFFERENT run (an
    # interrupted run never wrote its finals; the live history pickle can be
    # empty; a re-run or an HPC fetch can overwrite results/ with another run's
    # files).  The per-evaluation ``objective.txt`` folders and the optimizer
    # checkpoint are written continuously and are the ground truth, so reconcile
    # against them: if they reveal MORE successful evaluations — or a BETTER
    # best — than the loaded aggregates, use them.  This keeps the report honest
    # for ANY user without depending on a pristine results/ directory.  In a
    # cleanly-finished run all sources agree, so this is a no-op.
    recon = _reconstruct_history_from_evals(work_dir)
    if _n_valid_evals(recon) > _n_valid_evals(cr.get("history")):
        logger.info(
            "Report: results/ history had %d valid evals; reconstructed %d from "
            "evaluations/ folders - using the reconstruction.",
            _n_valid_evals(cr.get("history")), _n_valid_evals(recon))
        cr["history"] = recon
        cr["n_evaluations"] = max(cr.get("n_evaluations", 0), len(recon))
        cr["calibration_ran"] = True

    ck_best, ck_params, ck_eval = _load_checkpoint_best(work_dir)

    # Best objective = lowest (best) finite, non-penalty value across every
    # source: the loaded aggregate, the checkpoint, and the reconstructed evals.
    _cands = []   # list of (objective, params)
    _bo = cr.get("best_objective")
    if isinstance(_bo, (int, float)) and _bo == _bo:
        _cands.append((_bo, cr.get("best_params") or {}))
    if ck_best is not None:
        _cands.append((ck_best, ck_params or {}))
    _rb = min((e for e in cr.get("history", [])
               if isinstance(e.get("objective"), (int, float))),
              key=lambda e: e["objective"], default=None)
    if _rb is not None:
        _cands.append((_rb["objective"], _rb.get("parameters") or {}))
    _finite = [(o, p) for (o, p) in _cands if o < _REPORT_PENALTY]
    if _finite:
        _best_o, _best_p = min(_finite, key=lambda x: x[0])
        cr["best_objective"] = _best_o
        if _best_p:
            cr["best_params"] = _best_p
    if ck_eval is not None:
        cr["n_evaluations"] = max(cr.get("n_evaluations", 0), ck_eval + 1)

    # ── Sensitivity state ──
    sr = None
    sp = results_dir / "sensitivity_results.json"
    if sp.exists():
        try:
            with open(sp) as f:
                sr = json.load(f)
        except Exception:
            sr = None

    # Sensitivity expected but not finished → surface live screening progress
    # (count of completed model runs) so the banner isn't empty during Morris.
    if sr is None and mode in ("sensitivity", "both"):
        _done = _count_sensitivity_evals(work_dir)
        if _done:
            _traj = cs.get("sensitivity_morris_trajectories")
            try:
                _total = (int(_traj) * (len(calibration_parameters) + 1)
                          if _traj else None)
            except (TypeError, ValueError):
                _total = None
            cr["sensitivity_progress"] = (f"{_done}/{_total}" if _total
                                          else str(_done))

    # ── Matched data (performance / observation map) ──
    md = None
    md_path = results_dir / "matched_data.csv"
    if md_path.exists():
        try:
            import pandas as _pd
            md = _pd.read_csv(md_path)
        except Exception:
            md = None

    # Bulletproof best fit: rebuild matched + simulated from the GLOBAL-best
    # eval folder (ground truth across all runs/resumes), overriding the
    # possibly-stale/last-eval CSV.  This guarantees the report's best fit is
    # the true best for any basin, however many times it was interrupted.
    try:
        _bmd, _bsd, _beid = _rebuild_best_fit_from_evals(
            work_dir, cs, model_config)
        if _bmd is not None and not _bmd.empty:
            md = _bmd
            try:
                md.to_csv(md_path, index=False)            # refresh cache
                if _bsd is not None and not _bsd.empty:
                    _bsd.to_csv(results_dir / "simulated_data.csv", index=False)
            except Exception:
                pass
            logger.info("Report best fit taken from %s (global best across all "
                        "eval folders).", _beid)
    except Exception as _e:
        logger.debug(f"best-fit rebuild skipped: {_e}")

    # ── Performance metrics (same best-effort path as the driver) ──
    perf = None
    if md is not None and not md.empty:
        try:
            from .objective_functions import compute_all_metrics
            perf = compute_all_metrics(
                md, (cs.get("calibration_targets") or {}).get("species", []))
        except Exception:
            perf = None

    # Nothing computed yet → definitely still in progress.
    if not cr.get("history") and sr is None:
        in_progress = True

    return _GRR.generate_results_report(
        output_dir=str(work_dir), model_config=model_config or {},
        calibration_parameters=calibration_parameters,
        calibration_settings=cs, calibration_results=cr,
        sensitivity_results=sr, performance_metrics=perf,
        matched_data=md, report_stem=report_stem,
        in_progress=in_progress)


def _run_sensitivity_analysis(
        method: str,
        param_handler: ParameterHandler,
        model_runner: ModelRunner,
        obj_func: ObjectiveFunction,
        calibration_parameters: List[Dict],
        bounds_opt: List,
        transforms: List[str],
        morris_trajectories: int,
        morris_levels: int,
        sobol_samples: int,
        threshold: float,
        work_dir: Path,
        random_seed: int,
        n_parallel: int = 1
):
    """Run sensitivity analysis.

    Sensitivity samples are independent, so their (slow) model runs are
    executed up to ``n_parallel`` at a time.  Per-sample setup and the
    objective computation stay on the main thread; only the model execution
    is parallelised.
    """

    param_names = [p["name"] for p in calibration_parameters]

    if method == "morris":
        sa = MorrisScreening(
            param_names=param_names,
            bounds=bounds_opt,
            num_trajectories=morris_trajectories,
            num_levels=morris_levels,
            seed=random_seed
        )
    elif method == "sobol":
        sa = SobolAnalysis(
            param_names=param_names,
            bounds=bounds_opt,
            num_samples=sobol_samples,
            seed=random_seed
        )
    else:
        raise ValueError(f"Unknown sensitivity method: {method}")

    # Generate samples
    samples = sa.generate_samples()
    logger.info(f"Generated {len(samples)} samples for {method} analysis")

    # Evaluate samples.  Model runs are independent → run up to n_parallel
    # at a time (works for Docker and Apptainer).  Setup + objective compute
    # stay on the main thread; only the model execution is parallelised.
    n_samples = len(samples)
    n_par = max(1, int(n_parallel or 1))
    outputs = [1e10] * n_samples

    for chunk_start in range(0, n_samples, n_par):
        chunk = list(range(chunk_start, min(chunk_start + n_par, n_samples)))

        # ── Setup each eval in the chunk (sequential, fast) ──
        eval_configs = []
        for i in chunk:
            params_real = np.array([
                ParameterHandler.transform_to_real(p, t)
                for p, t in zip(samples[i], transforms)
            ])
            # Setup (bake generation-time params before generation) + patch
            eval_dir = param_handler.setup_working_directory(
                10000 + i, calibration_parameters, params_real)
            param_handler.apply_parameters(
                eval_dir, calibration_parameters, params_real)

            params_file = eval_dir / "parameters.json"
            with open(params_file, 'w') as f:
                json.dump({name: float(val)
                           for name, val in zip(param_names, params_real)},
                          f, indent=2)

            eval_configs.append({
                "eval_id": 10000 + i,
                "eval_dir": eval_dir,
                "master_json": str(eval_dir / "openWQ_master.json"),
                "_idx": i,
            })

        logger.info(
            f"Sensitivity samples {chunk[0] + 1}-{chunk[-1] + 1}/{n_samples} "
            f"(running {len(eval_configs)} with n_parallel={n_par})")

        # ── Run the chunk (parallel when n_par > 1) ──
        run_results = model_runner.run_parallel_evaluations(eval_configs, n_par)
        success_by_id = {eid: ok for (eid, ok, _rt, _err) in run_results}

        # ── Objective for each (sequential; reads HDF5 output) ──
        for cfg in eval_configs:
            if success_by_id.get(cfg["eval_id"], False):
                outputs[cfg["_idx"]] = obj_func.compute(
                    cfg["eval_dir"] / "openwq_out")
            else:
                outputs[cfg["_idx"]] = 1e10

    outputs = np.array(outputs)

    # Analyze
    result = sa.analyze(samples, outputs, threshold=threshold)
    return result


def run_sensitivity_analysis(**kwargs) -> Dict:
    """
    Run only sensitivity analysis (without optimization).
    """
    # Extract relevant kwargs
    work_dir = Path(kwargs['calibration_work_dir'])
    work_dir.mkdir(parents=True, exist_ok=True)
    (work_dir / "results").mkdir(exist_ok=True)

    _model_config = kwargs.get('model_config')
    _container_runtime = kwargs.get('container_runtime', 'docker')

    # Calibration override for openWQ debug mode (see run_calibration) — keeps
    # evaluations from writing the 5 extra diagnostic HDF5 files per
    # species/compartment.  Default OFF; user-overridable via the report.
    _rmd = kwargs.get("run_mode_debug")
    if _rmd is not None and _model_config is not None:
        _model_config["run_mode_debug"] = bool(_rmd)

    # Initialize components (same as run_calibration)
    param_handler = ParameterHandler(
        calibration_work_dir=kwargs['calibration_work_dir'],
        model_config=_model_config,
        base_model_config_dir=kwargs.get('base_model_config_dir'),
        test_case_dir=kwargs.get('test_case_dir'),
        running_on_docker=(_container_runtime == "docker"),
        calibration_period=kwargs.get("calibration_period"),
    )

    model_runner = ModelRunner(
        runtime=_container_runtime,
        docker_container_name=kwargs.get('docker_container_name'),
        docker_compose_path=kwargs.get('docker_compose_path'),
        apptainer_sif_path=kwargs.get('apptainer_sif_path'),
        apptainer_bind_path=kwargs.get('apptainer_bind_path'),
        executable_name=kwargs.get('executable_name', "mizuroute_lakes_openwq_Release"),
        executable_args=kwargs.get('executable_args', ""),
        file_manager_path=kwargs.get('file_manager_path'),
        executable_full_path=kwargs.get('executable_full_path'),
        command_template=kwargs.get('command_template'),
        hostmodel=(_model_config.get("hostmodel", "") if _model_config else ""),
        calibration_work_dir=kwargs.get('calibration_work_dir'),
        calibration_period=kwargs.get("calibration_period"),
    )

    if kwargs.get('base_model_config_dir'):
        h5_reader_path = str(Path(kwargs['base_model_config_dir']).parent / "2_Read_Outputs" / "hdf5_support_lib")
    else:
        _this_dir = Path(__file__).resolve().parent
        h5_reader_path = str(_this_dir.parent.parent / "2_Read_Outputs" / "hdf5_support_lib")

    # Get temporal resolution settings
    temporal_resolution = kwargs.get("temporal_resolution", "native")
    aggregation_method = kwargs.get("aggregation_method", "mean")

    # Resolve observation data path
    _obs_path = kwargs.get('observation_data_path', '')
    if not _obs_path and _model_config:
        from . import config_integration
        _obs = config_integration.get_observation_config(_model_config)
        obs_source = _obs.get("source", "skip")
        if obs_source in ("grqa", "user_csv"):
            # Reuse the observations the model config points to (GRQA folder/
            # "auto" download, or user_csv), reshaped to the objective CSV
            # format (no re-download).
            _log_adapter = lambda *a, **k: logger.info(
                " ".join(str(x) for x in a)) if a else None
            _obs_path = config_integration.prepare_calibration_observations_csv(
                _model_config, str(work_dir), log=_log_adapter)
            if not _obs_path:
                _obs_path = (
                    _obs.get("grqa_local_data_path", "")
                    if obs_source == "grqa"
                    else _obs.get("user_observation_csv", ""))
        _obs_path = _validate_observation_csv(_obs_path, obs_source)

    # Re-import defensively so this branch works even when the obs-path
    # resolution above (which is the only other importer of
    # config_integration in this function) didn't fire.
    from . import config_integration as _cint
    _spatial = _cint.get_spatial_mapping(_model_config) \
        if _model_config is not None else {
            "hostmodel": "mizuroute", "h5_mapping_key": "reachID"}
    obj_func = ObjectiveFunction(
        observation_path=_obs_path,
        target_species=kwargs['calibration_targets']["species"],
        target_reaches=kwargs['calibration_targets'].get("reach_ids", "all"),
        compartments=kwargs['calibration_targets'].get("compartments", ["RIVER_NETWORK_REACHES"]),
        metric=kwargs.get('objective_function', "KGE"),
        weights=kwargs.get('objective_weights'),
        h5_reader_path=h5_reader_path,
        temporal_resolution=temporal_resolution,
        aggregation_method=aggregation_method,
        hostmodel=_spatial.get("hostmodel", "mizuroute"),
        h5_mapping_key=_spatial.get("h5_mapping_key"),
        use_primary_only=kwargs.get("use_primary_only", True),
    )

    # Parameter info
    calibration_parameters = kwargs['calibration_parameters']
    param_names = [p["name"] for p in calibration_parameters]
    bounds = [p["bounds"] for p in calibration_parameters]
    transforms = [p.get("transform", "linear") for p in calibration_parameters]

    bounds_opt = []
    for b, t in zip(bounds, transforms):
        if t == "log":
            bounds_opt.append((np.log10(b[0]), np.log10(b[1])))
        else:
            bounds_opt.append(b)

    result = _run_sensitivity_analysis(
        method=kwargs.get('sensitivity_method', 'morris'),
        param_handler=param_handler,
        model_runner=model_runner,
        obj_func=obj_func,
        calibration_parameters=calibration_parameters,
        bounds_opt=bounds_opt,
        transforms=transforms,
        morris_trajectories=kwargs.get('sensitivity_morris_trajectories', 10),
        morris_levels=kwargs.get('sensitivity_morris_levels', 4),
        sobol_samples=kwargs.get('sensitivity_sobol_samples', 1024),
        threshold=kwargs.get('sensitivity_threshold', 0.1),
        work_dir=work_dir,
        random_seed=kwargs.get('random_seed'),
        n_parallel=kwargs.get('n_parallel', 1)
    )

    result.print_summary()
    result.save(work_dir / "results" / "sensitivity_results.json")

    return {
        "influential_params": result.influential_params,
        "rankings": result.rankings,
        "n_evaluations": result.n_evaluations
    }
