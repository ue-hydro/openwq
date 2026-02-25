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
from .optimization.dds import DDS, DDSResult, RandomSearch
from .sensitivity.morris_screening import MorrisScreening
from .sensitivity.sobol_analysis import SobolAnalysis
from .postprocessing.results_analysis import ResultsAnalyzer


def run_calibration(
        # Paths
        base_model_config_dir: str,
        calibration_work_dir: str,
        observation_data_path: str,
        test_case_dir: str,

        # Container runtime
        container_runtime: str,
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

    logger.info("=" * 60)
    logger.info("OPENWQ CALIBRATION FRAMEWORK")
    logger.info("=" * 60)
    logger.info(f"Start time: {datetime.now()}")
    logger.info(f"Working directory: {work_dir}")
    logger.info(f"Parameters to calibrate: {len(calibration_parameters)}")
    logger.info(f"Algorithm: {algorithm}")
    logger.info(f"Max evaluations: {max_evaluations}")
    logger.info(f"Container runtime: {container_runtime}")

    # Validate inputs
    if not calibration_parameters:
        raise ValueError("calibration_parameters cannot be empty")
    if not calibration_targets or "species" not in calibration_targets:
        raise ValueError("calibration_targets must include 'species'")

    # Initialize components
    logger.info("Initializing components...")

    param_handler = ParameterHandler(
        base_model_config_dir=base_model_config_dir,
        test_case_dir=test_case_dir,
        calibration_work_dir=calibration_work_dir,
        running_on_docker=(container_runtime == "docker")
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
        command_template=command_template
    )

    # H5 reader path for objective function
    h5_reader_path = str(Path(base_model_config_dir).parent / "Read_Outputs" / "hdf5_support_lib")

    # Get temporal resolution settings from kwargs
    temporal_resolution = kwargs.get("temporal_resolution", "native")
    aggregation_method = kwargs.get("aggregation_method", "mean")

    logger.info(f"Temporal resolution: {temporal_resolution}")
    logger.info(f"Aggregation method: {aggregation_method}")

    obj_func = ObjectiveFunction(
        observation_path=observation_data_path,
        target_species=calibration_targets["species"],
        target_reaches=calibration_targets.get("reach_ids", "all"),
        compartments=calibration_targets.get("compartments", ["RIVER_NETWORK_REACHES"]),
        metric=objective_function,
        weights=objective_weights,
        h5_reader_path=h5_reader_path,
        temporal_resolution=temporal_resolution,
        aggregation_method=aggregation_method
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

    sensitive_params = None
    if run_sensitivity_first:
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
            random_seed=random_seed
        )

        sensitive_params = sa_result.influential_params
        sa_result.print_summary()
        sa_result.save(work_dir / "results" / "sensitivity_results.json")

        logger.info(f"Influential parameters: {sensitive_params}")
        logger.info(f"Non-influential parameters will be fixed at initial values")

        # TODO: Filter calibration_parameters to only include sensitive ones

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

    def objective_wrapper(params_opt: np.ndarray) -> float:
        """Wrapper that handles parameter transformation and model execution."""
        eval_id = eval_counter[0]
        eval_counter[0] += 1

        logger.info(f"--- Evaluation {eval_id} ---")

        # Transform from optimization space to real space
        params_real = np.array([
            ParameterHandler.transform_to_real(p, t)
            for p, t in zip(params_opt, transforms)
        ])

        logger.debug(f"Parameters (real): {dict(zip(param_names, params_real))}")

        # Setup working directory
        eval_dir = param_handler.setup_working_directory(eval_id)

        # Apply parameters
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
        return obj_val

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
        optimizer = RandomSearch(
            n_params=len(calibration_parameters),
            bounds=bounds_opt,
            param_names=param_names,
            max_evals=max_evaluations - start_eval,
            seed=random_seed
        )
        result = optimizer.optimize(
            objective_wrapper,
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

    # Save calibration history in the format expected by ResultsAnalyzer
    # result.history entries are tuples: (eval_num, objective, params_array)
    calibration_history = []
    for entry in result.history:
        eval_num, objective, params_array = entry  # Unpack tuple
        hist_entry = {
            "eval_id": int(eval_num),
            "objective": float(objective),
            "parameters": {
                name: ParameterHandler.transform_to_real(val, transforms[i])
                for i, (name, val) in enumerate(
                    zip(param_names, params_array))
            },
            "timestamp": 0,
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
    _base_cfg = Path(base_model_config_dir)
    _shapefile_dir = kwargs.get("shapefile_dir", "")
    if not _shapefile_dir and _base_cfg.exists():
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

    # Summary
    logger.info("=" * 60)
    logger.info("CALIBRATION COMPLETE")
    logger.info("=" * 60)
    logger.info(f"Best objective: {result.best_objective:.6f}")
    logger.info(f"Total evaluations: {result.n_evaluations + start_eval}")
    logger.info(f"Best parameters:")
    for name, val in best_params_real.items():
        logger.info(f"  {name}: {val}")
    logger.info(f"Results saved to: {results_dir}")

    return {
        "best_params": best_params_real,
        "best_objective": result.best_objective,
        "n_evaluations": result.n_evaluations + start_eval,
        "sensitive_params": sensitive_params,
        "results_dir": str(results_dir)
    }


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
        random_seed: int
):
    """Run sensitivity analysis."""

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

    # Evaluate samples
    outputs = []
    for i, sample in enumerate(samples):
        logger.info(f"Sensitivity sample {i+1}/{len(samples)}")

        # Transform to real space
        params_real = np.array([
            ParameterHandler.transform_to_real(p, t)
            for p, t in zip(sample, transforms)
        ])

        # Setup and run
        eval_dir = param_handler.setup_working_directory(10000 + i)
        param_handler.apply_parameters(eval_dir, calibration_parameters, params_real)

        # Save parameters
        params_file = eval_dir / "parameters.json"
        with open(params_file, 'w') as f:
            json.dump({name: float(val) for name, val in zip(param_names, params_real)}, f, indent=2)

        master_json = str(eval_dir / "openWQ_master.json")
        success, runtime, error = model_runner.run_single_evaluation(eval_dir, master_json, 10000 + i)

        if success:
            output_dir = eval_dir / "openwq_out"
            obj_val = obj_func.compute(output_dir)
        else:
            obj_val = 1e10

        outputs.append(obj_val)

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

    # Initialize components (same as run_calibration)
    param_handler = ParameterHandler(
        base_model_config_dir=kwargs['base_model_config_dir'],
        test_case_dir=kwargs['test_case_dir'],
        calibration_work_dir=kwargs['calibration_work_dir'],
        running_on_docker=(kwargs['container_runtime'] == "docker")
    )

    model_runner = ModelRunner(
        runtime=kwargs['container_runtime'],
        docker_container_name=kwargs.get('docker_container_name'),
        docker_compose_path=kwargs.get('docker_compose_path'),
        apptainer_sif_path=kwargs.get('apptainer_sif_path'),
        apptainer_bind_path=kwargs.get('apptainer_bind_path'),
        executable_name=kwargs.get('executable_name', "mizuroute_lakes_openwq_Release"),
        executable_args=kwargs.get('executable_args', ""),
        file_manager_path=kwargs.get('file_manager_path'),
        executable_full_path=kwargs.get('executable_full_path'),
        command_template=kwargs.get('command_template')
    )

    h5_reader_path = str(Path(kwargs['base_model_config_dir']).parent / "Read_Outputs" / "hdf5_support_lib")

    # Get temporal resolution settings
    temporal_resolution = kwargs.get("temporal_resolution", "native")
    aggregation_method = kwargs.get("aggregation_method", "mean")

    obj_func = ObjectiveFunction(
        observation_path=kwargs['observation_data_path'],
        target_species=kwargs['calibration_targets']["species"],
        target_reaches=kwargs['calibration_targets'].get("reach_ids", "all"),
        compartments=kwargs['calibration_targets'].get("compartments", ["RIVER_NETWORK_REACHES"]),
        metric=kwargs.get('objective_function', "KGE"),
        weights=kwargs.get('objective_weights'),
        h5_reader_path=h5_reader_path,
        temporal_resolution=temporal_resolution,
        aggregation_method=aggregation_method
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
        random_seed=kwargs.get('random_seed')
    )

    result.print_summary()
    result.save(work_dir / "results" / "sensitivity_results.json")

    return {
        "influential_params": result.influential_params,
        "rankings": result.rankings,
        "n_evaluations": result.n_evaluations
    }
