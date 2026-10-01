# Copyright 2026, Diogo Costa, diogo.costa@uevora.pt
# This file is part of OpenWQ model.
#
# This program, openWQ, is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) later version.
# ...GPLv3 (see COPYING).

"""
Chained-model calibration orchestration.
========================================

Runs a SEQUENTIAL chain of host_model-openWQ models per calibration evaluation
(e.g. SUMMA-openWQ -> mizuRoute-openWQ). The chain is an ordered list; index 0
is the most UPSTREAM model, the last entry is the DOWNSTREAM model that is
VALIDATED against observations (its output drives the objective).

Design — zero changes to the DDS / objective / checkpoint machinery:
    `ChainParameterHandler` and `ChainModelRunner` implement the SAME interface
    the single-model `ParameterHandler` / `ModelRunner` expose to the driver
    (`setup_working_directory`, `apply_parameters`, `run_single_evaluation`,
    `preflight`, `calibration_period`). Internally they orchestrate N models.

Per-eval directory layout (so the objective reads `eval_dir/openwq_out` as
today — the LAST model runs directly in `eval_dir`):

    evaluations/eval_0000/
        _chain/m0/           <- upstream model 0 (e.g. SUMMA): config + output
            openwq_out/HDF5/  ...the flux h5 the next model ingests
        _chain/m1/           <- intermediate model 1 (if any)
        openWQ_master.json   <- LAST model (validated) runs here
        openwq_out/HDF5/     <- scored against observations

The downstream EWF is rewired automatically to each eval's upstream folder
(the template only fixes the run order + declares ewf_method="from_openwq_hdf5").
Parameters are namespaced by chain position (each param dict carries
`model_index`); apply routes each to its model's run directory.
"""

import copy
import shutil
import time
import logging
from pathlib import Path

from .parameter_handler import ParameterHandler
from .model_runner import ModelRunner
from . import config_integration

logger = logging.getLogger(__name__)


class ChainParameterHandler:
    """Per-eval config generation + parameter application for a model chain.

    Mirrors ``ParameterHandler`` (the driver calls ``setup_working_directory``
    then ``apply_parameters``), but generates one config per chain model into
    its own run directory and rewires the downstream EWF.
    """

    def __init__(self, chain_configs, calibration_work_dir,
                 calibration_period=None, running_on_docker=True,
                 ml_closures=None, ml_runtime=None):
        self.chain_configs = list(chain_configs)
        self.n = len(self.chain_configs)
        self.calibration_work_dir = Path(calibration_work_dir)
        self.calibration_period = calibration_period

        # Hybrid physics-ML: route each closure / runtime-NN entry to its model
        # by the ``model_index`` the report baked into it (default 0).
        def _for_model(d, i):
            return {k: v for k, v in (d or {}).items()
                    if int(v.get("model_index", 0)) == i}

        # One ParameterHandler per model — reused only for its (eval_dir-based)
        # apply_parameters(), which edits the already-generated config JSONs.
        self._handlers = [
            ParameterHandler(
                calibration_work_dir=str(calibration_work_dir),
                model_config=cfg,
                running_on_docker=running_on_docker,
                calibration_period=calibration_period,
                ml_closures=_for_model(ml_closures, i),
                ml_runtime=_for_model(ml_runtime, i),
            )
            for i, cfg in enumerate(self.chain_configs)
        ]

    # ---- directory / path helpers -------------------------------------------
    def _rundir(self, eval_dir, i):
        """Run directory for model i. The LAST model runs in eval_dir itself
        (so its output lands at eval_dir/openwq_out where the objective reads);
        upstream models run in eval_dir/_chain/m{i}."""
        eval_dir = Path(eval_dir)
        return eval_dir if i == self.n - 1 else eval_dir / "_chain" / f"m{i}"

    def _ewf_src_rel(self, i):
        """Relative EWF source folder (upstream model i-1's HDF5 output) as
        seen from model i's run CWD."""
        if i == self.n - 1:
            # last model runs in eval_dir; upstream lives under _chain/
            return f"_chain/m{i - 1}/openwq_out/HDF5"
        # intermediate model runs in _chain/m{i}; sibling is _chain/m{i-1}
        return f"../m{i - 1}/openwq_out/HDF5"

    def _split(self, calibration_parameters, params_real, i):
        """Return (params, values) belonging to model i (by ``model_index``)."""
        # NB: params_real is a numpy array — use explicit None checks, never
        # `x or []` (ambiguous truth value on arrays).
        _cp = calibration_parameters if calibration_parameters is not None else []
        _pr = params_real if params_real is not None else []
        pi, vi = [], []
        for p, v in zip(_cp, _pr):
            if int(p.get("model_index", 0)) == i:
                pi.append(p)
                vi.append(v)
        return pi, vi

    # ---- interface used by the driver ---------------------------------------
    def setup_working_directory(self, eval_id, calibration_parameters=None,
                                params_real=None):
        eval_dir = self.calibration_work_dir / "evaluations" / f"eval_{eval_id:04d}"
        if eval_dir.exists():
            shutil.rmtree(eval_dir)
        eval_dir.mkdir(parents=True)
        (eval_dir / "openwq_out" / "HDF5").mkdir(parents=True, exist_ok=True)

        for i, cfg in enumerate(self.chain_configs):
            rundir = self._rundir(eval_dir, i)
            rundir.mkdir(parents=True, exist_ok=True)
            (rundir / "openwq_out" / "HDF5").mkdir(parents=True, exist_ok=True)

            eval_cfg = copy.deepcopy(cfg)
            # Rewire a downstream from_openwq_hdf5 EWF to THIS eval's upstream
            # model output (the template's own ewf_h5_source_folder pointed at
            # the standalone domain layout).
            if i > 0 and str(eval_cfg.get("ewf_method", "")).lower() == "from_openwq_hdf5":
                eval_cfg["ewf_h5_source_folder"] = self._ewf_src_rel(i)

            config_integration.generate_config_for_eval(
                eval_cfg, str(rundir), suppress_report=True,
                calibration_period=self.calibration_period,
            )
        return eval_dir

    def apply_parameters(self, eval_dir, calibration_parameters, params_real):
        for i in range(self.n):
            pi, vi = self._split(calibration_parameters, params_real, i)
            if not pi:
                continue
            self._handlers[i].apply_parameters(self._rundir(eval_dir, i), pi, vi)


class ChainModelRunner:
    """Runs each model in the chain, in order, into its own run directory.

    Mirrors ``ModelRunner.run_single_evaluation`` — the objective then reads
    the LAST model's output at ``eval_dir/openwq_out``.
    """

    def __init__(self, runners, n_models):
        self._runners = list(runners)
        self.n = n_models
        self._calibration_period = getattr(self._runners[-1], "calibration_period", None)
        # attributes the driver reads/sets on a ModelRunner
        self._calib_total = 0
        self.total_evaluations = getattr(self._runners[-1], "total_evaluations", 0)

    # calibration_period is set by the driver (e.g. the validation re-run);
    # propagate it to every model so the whole chain uses the same window.
    @property
    def calibration_period(self):
        return self._calibration_period

    @calibration_period.setter
    def calibration_period(self, value):
        self._calibration_period = value
        for r in self._runners:
            try:
                r.calibration_period = value
            except Exception:
                pass

    def _rundir(self, eval_dir, i):
        eval_dir = Path(eval_dir)
        return eval_dir if i == self.n - 1 else eval_dir / "_chain" / f"m{i}"

    def preflight(self):
        for idx, r in enumerate(self._runners):
            ok, msg = r.preflight()
            if not ok:
                return False, f"chain model {idx}: {msg}"
        return True, "chain preflight ok"

    def run_single_evaluation(self, eval_dir, master_json_path, eval_id):
        t0 = time.time()
        _hms = [(getattr(r, "hostmodel", "") or ("model%d" % k))
                for k, r in enumerate(self._runners)]
        logger.info("  ⛓ chain: %s  (m%d output is scored vs observations)",
                    " → ".join("m%d:%s-openWQ" % (k, h)
                                    for k, h in enumerate(_hms)), self.n - 1)
        for i in range(self.n):
            rundir = self._rundir(eval_dir, i)   # Path (ModelRunner calls .resolve())
            _role = ("\U0001f3af target — scored vs observations"
                     if i == self.n - 1
                     else "upstream — output feeds the next model via EWF")
            logger.info("  ⛓ chain step %d/%d: m%d = %s-openWQ  (%s)  [run dir: %s]",
                        i + 1, self.n, i, _hms[i], _role, rundir.name)
            mj = str(rundir / "openWQ_master.json")
            ok, _rt, err = self._runners[i].run_single_evaluation(
                rundir, mj, eval_id)
            if not ok:
                return False, time.time() - t0, \
                    f"chain model {i} (m{i}:{_hms[i]}, {rundir.name}) failed: {err}"
        return True, time.time() - t0, ""

    def run_parallel_evaluations(self, eval_configs, n_parallel):
        """Chain-aware parallel evaluation — mirrors
        ``ModelRunner.run_parallel_evaluations`` so the parallel algorithms
        (DDS_PARALLEL, RANDOM) and the sensitivity stage work for a chain too,
        not just plain sequential DDS.

        Each eval's chain (m0 → … → m{n-1}) runs SEQUENTIALLY inside itself —
        the EWF coupling needs the upstream model's HDF5 output before the
        downstream model runs — but INDEPENDENT evaluations run concurrently:
        every eval has its own ``eval_dir`` (and its own ``_chain/m{i}``), so
        they never collide.  Returns ``[(eval_id, ok, runtime, error), …]`` in
        the SAME order as ``eval_configs`` (what the driver's batch objective
        expects).
        """
        import concurrent.futures as _cf

        def _one(cfg):
            eid = int(cfg["eval_id"])
            ok, rt, err = self.run_single_evaluation(
                Path(cfg["eval_dir"]), cfg.get("master_json", ""), eid)
            return (eid, ok, rt, err)

        n_par = max(1, int(n_parallel or 1))
        if n_par <= 1 or len(eval_configs) <= 1:
            return [_one(c) for c in eval_configs]

        results = [None] * len(eval_configs)
        with _cf.ThreadPoolExecutor(
                max_workers=min(n_par, len(eval_configs))) as ex:
            _fut_to_i = {ex.submit(_one, c): i
                         for i, c in enumerate(eval_configs)}
            for _fut in _cf.as_completed(_fut_to_i):
                results[_fut_to_i[_fut]] = _fut.result()
        return results


def build_chain(chain_configs, calibration_work_dir, container_runtime,
                calibration_period, max_evaluations,
                docker_container_name=None, docker_compose_path=None,
                apptainer_sif_path=None, apptainer_bind_path=None,
                command_template=None, executable_args="",
                ml_closures=None, ml_runtime=None):
    """Build a (ChainParameterHandler, ChainModelRunner) pair from an ordered
    list of loaded model configs. Each model's executable / control file /
    hostmodel are taken from its OWN container config, so a SUMMA model and a
    mizuRoute model in the same chain each run with the right binary."""
    param_handler = ChainParameterHandler(
        chain_configs, calibration_work_dir,
        calibration_period=calibration_period,
        running_on_docker=(container_runtime == "docker"),
        ml_closures=ml_closures,
        ml_runtime=ml_runtime,
    )

    runners = []
    for cfg in chain_configs:
        c = config_integration.get_container_config(cfg)
        runners.append(ModelRunner(
            runtime=container_runtime,
            docker_container_name=docker_container_name,
            docker_compose_path=docker_compose_path or c.get("docker_compose_path"),
            apptainer_sif_path=apptainer_sif_path,
            apptainer_bind_path=apptainer_bind_path,
            executable_full_path=c.get("executable_path", ""),
            executable_args=executable_args,
            file_manager_path=c.get("file_manager_path", ""),
            command_template=command_template,
            hostmodel=cfg.get("hostmodel", ""),
            calibration_work_dir=str(calibration_work_dir),
            calibration_period=calibration_period,
            total_evaluations=max_evaluations,
        ))
    model_runner = ChainModelRunner(runners, len(chain_configs))
    return param_handler, model_runner
