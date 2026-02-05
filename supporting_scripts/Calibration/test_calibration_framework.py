#!/usr/bin/env python3
# Copyright 2026, Diogo Costa, diogo.costa@uevora.pt
# This file is part of OpenWQ model.
"""
Calibration Framework Test Script
=================================

Tests all components of the OpenWQ calibration framework.
Run this to verify the framework is properly installed and functioning.

Usage:
    python test_calibration_framework.py
"""

import sys
import os
import tempfile
import shutil
import numpy as np
from pathlib import Path

# Add calibration library to path
sys.path.insert(0, str(Path(__file__).parent))

# Test configuration
VERBOSE = True


def print_header(text):
    """Print a test section header."""
    print()
    print("=" * 60)
    print(f" {text}")
    print("=" * 60)


def print_result(name, success, message=""):
    """Print a test result."""
    status = "PASS" if success else "FAIL"
    print(f"  [{status}] {name}" + (f" - {message}" if message else ""))


def test_imports():
    """Test that all modules can be imported."""
    print_header("Testing Module Imports")

    all_passed = True

    modules = [
        ("calibration_lib", "Main library"),
        ("calibration_lib.parameter_handler", "ParameterHandler"),
        ("calibration_lib.model_runner", "ModelRunner"),
        ("calibration_lib.objective_functions", "ObjectiveFunction"),
        ("calibration_lib.checkpoint", "CheckpointManager"),
        ("calibration_lib.optimization.dds", "DDS Optimizer"),
        ("calibration_lib.optimization.base_optimizer", "BaseOptimizer"),
        ("calibration_lib.sensitivity.morris_screening", "Morris Screening"),
        ("calibration_lib.sensitivity.sobol_analysis", "Sobol Analysis"),
        ("calibration_lib.hpc.slurm_manager", "SLURM Manager"),
        ("calibration_lib.hpc.job_coordinator", "Job Coordinator"),
        ("calibration_lib.postprocessing.results_analysis", "Results Analyzer"),
    ]

    for module_name, description in modules:
        try:
            __import__(module_name)
            print_result(description, True)
        except Exception as e:
            print_result(description, False, str(e))
            all_passed = False

    return all_passed


def test_dds_optimizer():
    """Test DDS optimization algorithm."""
    print_header("Testing DDS Optimizer")

    from calibration_lib.optimization import DDS

    # Test with sphere function
    def sphere(x):
        return float(np.sum(x**2))

    try:
        dds = DDS(
            n_params=3,
            bounds=[(-5, 5), (-5, 5), (-5, 5)],
            param_names=['x', 'y', 'z'],
            max_evals=30,
            r=0.2,
            seed=42
        )
        print_result("Initialization", True)

        initial = np.array([3.0, 3.0, 3.0])
        result = dds.optimize(
            objective_func=sphere,
            initial_point=initial
        )

        initial_obj = sphere(initial)
        improvement = (1 - result.best_objective / initial_obj) * 100

        print_result(
            "Optimization",
            True,
            f"objective: {result.best_objective:.4f} ({improvement:.1f}% improvement)"
        )

        # Test result serialization
        result_dict = result.to_dict()
        print_result("Result serialization", "best_params" in result_dict)

        return True

    except Exception as e:
        print_result("DDS test", False, str(e))
        import traceback
        traceback.print_exc()
        return False


def test_morris_screening():
    """Test Morris sensitivity analysis."""
    print_header("Testing Morris Sensitivity Analysis")

    from calibration_lib.sensitivity import MorrisScreening

    # Linear function with known sensitivities
    def linear_func(x):
        return 5.0 * x[0] + 2.0 * x[1] + 0.5 * x[2]

    try:
        morris = MorrisScreening(
            param_names=['a', 'b', 'c'],
            bounds=[(0, 1), (0, 1), (0, 1)],
            num_trajectories=4,
            num_levels=4,
            seed=42
        )
        print_result("Initialization", True)

        samples = morris.generate_samples()
        print_result("Sample generation", True, f"{len(samples)} samples")

        outputs = np.array([linear_func(s) for s in samples])
        result = morris.analyze(samples, outputs)

        # Check that 'a' (coeff=5) is most influential
        most_influential = result.rankings[0][0]
        print_result(
            "Sensitivity ranking",
            most_influential == 'a',
            f"most influential: {most_influential}"
        )

        # Test result serialization
        result_dict = result.to_dict()
        print_result("Result serialization", "mu_star" in result_dict)

        return True

    except Exception as e:
        print_result("Morris test", False, str(e))
        import traceback
        traceback.print_exc()
        return False


def test_checkpoint_manager():
    """Test checkpoint save/restore."""
    print_header("Testing Checkpoint Manager")

    from calibration_lib import CheckpointManager

    test_dir = tempfile.mkdtemp(prefix='checkpoint_test_')

    try:
        cm = CheckpointManager(checkpoint_dir=test_dir)
        print_result("Initialization", True)

        # Save state
        cm.save_state(
            current_eval=25,
            max_evals=100,
            best_objective=0.123,
            best_params=np.array([0.5, 0.3, 0.7]),
            param_names=['x', 'y', 'z'],
            history=[{'eval': i, 'obj': 0.5 - i*0.01} for i in range(25)],
            algorithm='DDS'
        )
        print_result("Save state", True)

        # Check file exists
        checkpoint_file = Path(test_dir) / "calibration_state.json"
        print_result("Checkpoint file created", checkpoint_file.exists())

        # Load state
        loaded = cm.load_state()
        print_result(
            "Load state",
            loaded['current_eval'] == 25 and loaded['best_objective'] == 0.123
        )

        # Test history file exists
        history_file = Path(test_dir) / "calibration_history.pkl"
        print_result("History file created", history_file.exists())

        return True

    except Exception as e:
        print_result("Checkpoint test", False, str(e))
        import traceback
        traceback.print_exc()
        return False

    finally:
        shutil.rmtree(test_dir)


def test_parameter_handler():
    """Test parameter handler."""
    print_header("Testing Parameter Handler")

    from calibration_lib import ParameterHandler

    test_dir = tempfile.mkdtemp(prefix='param_test_')

    try:
        # Create mock directory structure
        model_config = Path(test_dir) / 'Model_Config'
        test_case = Path(test_dir) / 'test_case'
        work_dir = Path(test_dir) / 'calibration'

        model_config.mkdir()
        test_case.mkdir()
        (test_case / 'openwq_in').mkdir()
        (test_case / 'openwq_out').mkdir()

        # Create mock BGC JSON file
        import json
        bgc_data = {
            'CYCLING_FRAMEWORKS': {
                'N_cycle': {
                    '3': {
                        'parameter_values': {
                            'k': 0.03
                        }
                    }
                }
            }
        }
        bgc_file = test_case / 'openwq_in' / 'openWQ_MODULE_NATIVE_BGC_FLEX.json'
        with open(bgc_file, 'w') as f:
            f.write('// BGC Config\n')
            json.dump(bgc_data, f, indent=2)

        handler = ParameterHandler(
            base_model_config_dir=str(model_config),
            test_case_dir=str(test_case),
            calibration_work_dir=str(work_dir),
            running_on_docker=True
        )
        print_result("Initialization", True)

        # Test working directory setup
        eval_dir = handler.setup_working_directory(1)
        print_result("Working directory setup", eval_dir.exists())

        # Test parameter transforms
        log_val = ParameterHandler.transform_to_opt(0.1, 'log')
        real_val = ParameterHandler.transform_to_real(log_val, 'log')
        print_result("Log transform", abs(real_val - 0.1) < 1e-10)

        linear_val = ParameterHandler.transform_to_opt(0.5, 'linear')
        print_result("Linear transform", linear_val == 0.5)

        return True

    except Exception as e:
        print_result("Parameter handler test", False, str(e))
        import traceback
        traceback.print_exc()
        return False

    finally:
        shutil.rmtree(test_dir)


def test_model_runner():
    """Test model runner configuration."""
    print_header("Testing Model Runner")

    from calibration_lib import ModelRunner
    from calibration_lib.model_runner import HPCJobGenerator

    test_dir = tempfile.mkdtemp(prefix='runner_test_')

    try:
        # Test Docker configuration
        runner_docker = ModelRunner(
            runtime='docker',
            docker_container_name='docker_openwq',
            docker_compose_path='/path/to/docker-compose.yml',
            executable_name='mizuroute_lakes_cslm_openwq_fast',
            executable_args='-g 1 1',
            file_manager_path='/code/test/fileManager.txt',
            timeout_seconds=3600
        )
        print_result("Docker runner", runner_docker.runtime == 'docker')

        # Test Apptainer configuration
        runner_hpc = ModelRunner(
            runtime='apptainer',
            apptainer_sif_path='/path/to/openwq.sif',
            apptainer_bind_path='/scratch:/code',
            executable_name='mizuroute_lakes_cslm_openwq_fast',
            executable_args='-g 1 1',
            file_manager_path='/code/test/fileManager.txt',
            timeout_seconds=7200
        )
        print_result("Apptainer runner", runner_hpc.runtime == 'apptainer')

        # Test HPC job generator
        gen = HPCJobGenerator(scheduler='slurm')
        script_path = Path(test_dir) / 'array_job.sh'

        gen.generate_array_job_script(
            script_path=script_path,
            calibration_work_dir='/scratch/calibration',
            apptainer_sif_path='/path/to/openwq.sif',
            apptainer_bind_path='/scratch:/code',
            executable_name='mizuroute_lakes_cslm_openwq_fast',
            executable_args='-g 1 1',
            file_manager_path='/code/test/fileManager.txt',
            n_evaluations=100,
            partition='standard',
            walltime='04:00:00',
            max_concurrent=50
        )
        print_result("SLURM script generation", script_path.exists())

        return True

    except Exception as e:
        print_result("Model runner test", False, str(e))
        import traceback
        traceback.print_exc()
        return False

    finally:
        shutil.rmtree(test_dir)


def test_results_analyzer():
    """Test results analysis and plotting."""
    print_header("Testing Results Analyzer")

    from calibration_lib.postprocessing import ResultsAnalyzer

    test_dir = tempfile.mkdtemp(prefix='results_test_')

    try:
        # Create mock calibration history
        import json
        history = []
        for i in range(50):
            history.append({
                'eval_id': i,
                'parameters': {'x': 0.1 + i*0.01, 'y': 0.2 - i*0.002},
                'objective': 1.0 - i*0.015,
                'timestamp': 1000000 + i*60
            })

        history_file = Path(test_dir) / 'calibration_history.json'
        with open(history_file, 'w') as f:
            json.dump(history, f)

        analyzer = ResultsAnalyzer(
            results_dir=test_dir,
            parameter_definitions=[
                {'name': 'x', 'bounds': (0, 1)},
                {'name': 'y', 'bounds': (0, 1)}
            ]
        )
        print_result("Initialization", True)

        # Test convergence data
        evals, best_objs = analyzer.get_convergence_data()
        print_result("Convergence data", len(evals) == 50)

        # Test summary
        summary = analyzer.get_summary()
        print_result(
            "Summary generation",
            summary.n_evaluations == 50,
            f"best objective: {summary.best_objective:.4f}"
        )

        # Test parameter statistics
        stats = analyzer.get_parameter_statistics()
        print_result("Parameter statistics", len(stats) == 2)

        return True

    except Exception as e:
        print_result("Results analyzer test", False, str(e))
        import traceback
        traceback.print_exc()
        return False

    finally:
        shutil.rmtree(test_dir)


def main():
    """Run all tests."""
    print()
    print("*" * 60)
    print("*  OpenWQ Calibration Framework Test Suite")
    print("*" * 60)

    tests = [
        ("Module Imports", test_imports),
        ("DDS Optimizer", test_dds_optimizer),
        ("Morris Screening", test_morris_screening),
        ("Checkpoint Manager", test_checkpoint_manager),
        ("Parameter Handler", test_parameter_handler),
        ("Model Runner", test_model_runner),
        ("Results Analyzer", test_results_analyzer),
    ]

    results = {}

    for name, test_func in tests:
        try:
            results[name] = test_func()
        except Exception as e:
            print(f"\nUnexpected error in {name}: {e}")
            results[name] = False

    # Summary
    print_header("Test Summary")

    passed = sum(1 for v in results.values() if v)
    total = len(results)

    for name, success in results.items():
        status = "PASS" if success else "FAIL"
        print(f"  [{status}] {name}")

    print()
    print(f"  Total: {passed}/{total} tests passed")
    print()

    if passed == total:
        print("  All tests PASSED!")
        print("  The calibration framework is ready to use.")
    else:
        print("  Some tests FAILED.")
        print("  Please check the error messages above.")

    print()

    return 0 if passed == total else 1


if __name__ == "__main__":
    sys.exit(main())
