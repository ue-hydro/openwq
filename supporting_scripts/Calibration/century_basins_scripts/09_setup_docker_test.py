#!/usr/bin/env python3
# Copyright 2026, Diogo Costa, diogo.costa@uevora.pt
# This file is part of OpenWQ model.
#
# This program, openWQ, is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) later version.

"""
09_setup_docker_test.py
========================

Build a Docker test workspace for 3 selected basins x 4 variants.
Creates symlinks to static data and rewrites config files with Docker paths.

The existing config files have absolute HPC paths (/projects/F20240...)
that don't work inside the Docker container. This script:
  1. Creates a docker_test/ workspace with symlinks to static data
  2. Rewrites openWQ_master.json with RELATIVE paths (works because
     param_handler copies both master.json and openwq_in/ to eval_dir)
  3. Rewrites fileManager.txt + mizuroute.control with Docker ABSOLUTE paths
  4. Generates Docker-ready calibration configs

Docker volume mount: /Users/diogocosta/Documents -> /code

USAGE:
  python 09_setup_docker_test.py \\
      [--prepared-dir /path/to/century_basins_prepared] \\
      [--docker-test-dir /path/to/docker_test] \\
      [--max-evaluations 5]
"""

import sys
import os
import json
import re
import shutil
import argparse
from pathlib import Path
from typing import List, Dict

SCRIPT_DIR = Path(__file__).resolve().parent

# ─────────────────────────────────────────────────────────────────────
# Docker path constants
# ─────────────────────────────────────────────────────────────────────

DOCKER_HOST_BASE = "/Users/diogocosta/Documents"
DOCKER_CONTAINER_BASE = "/code"

# Test basin selection (headwater, meso, macro)
TEST_BASINS = [
    "CAN_01DG003_headwater",
    "CAN_01AM001_meso",
    "CAN_02BA003_macro",
]

VARIANTS = {
    'A': 'variant_A_nitrogen_full',
    'B': 'variant_B_swat',
    'C': 'variant_C_thermodynamic',
    'D': 'variant_D_phreeqc',
}

# Executable info
EXECUTABLE_DIR = "/code/openwq_code/mizuRoute-OpenWQ/route/build/openwq/openwq/bin"
EXECUTABLE_NAME = "mizuroute_lakes_openwq_Release"

# Calibration lib on host
CALIBRATION_LIB_HOST = str(SCRIPT_DIR.parent)


# ─────────────────────────────────────────────────────────────────────
# Helper: convert host path to Docker container path
# ─────────────────────────────────────────────────────────────────────

def host_to_docker(host_path: str) -> str:
    """Convert a host absolute path to Docker container path."""
    return host_path.replace(DOCKER_HOST_BASE, DOCKER_CONTAINER_BASE)


def make_symlink(src: str, dst: str):
    """Create a RELATIVE symlink so it resolves in both host and Docker container.

    Absolute symlinks break inside Docker because the host path prefix
    (/Users/diogocosta/Documents) doesn't exist at the same location
    inside the container (/code).  Relative symlinks are path-prefix
    agnostic, so they work in both environments.
    """
    dst_path = Path(dst)
    if dst_path.is_symlink() or dst_path.exists():
        if dst_path.is_dir() and not dst_path.is_symlink():
            shutil.rmtree(dst)
        else:
            dst_path.unlink()
    # Compute relative path from the symlink's parent directory to the target
    rel_src = os.path.relpath(src, os.path.dirname(dst))
    os.symlink(rel_src, dst)


# ─────────────────────────────────────────────────────────────────────
# Rewrite openWQ_master.json with relative paths
# ─────────────────────────────────────────────────────────────────────

def rewrite_master_json_relative(src_path: str, dst_path: str):
    """
    Read openWQ_master.json and rewrite all paths to be relative.
    This works because param_handler copies openWQ_master.json + openwq_in/
    to the eval_dir, so relative paths like ./openwq_in/... resolve correctly.
    """
    # Read with comment stripping (JSON doesn't officially support comments)
    with open(src_path, 'r') as f:
        content = f.read()
    # Strip // comments
    lines = content.split('\n')
    clean_lines = [l for l in lines if not l.strip().startswith('//')]
    data = json.loads('\n'.join(clean_lines))

    # Rewrite OPENWQ_INPUT paths
    if "OPENWQ_INPUT" in data:
        inp = data["OPENWQ_INPUT"]
        if "CONFIG_FILEPATH" in inp:
            inp["CONFIG_FILEPATH"] = "./openwq_in/openWQ_config.json"
        if "SINK_SOURCE" in inp:
            for key, ss in inp["SINK_SOURCE"].items():
                if "FILEPATH" in ss:
                    # Keep just the filename
                    basename = os.path.basename(ss["FILEPATH"])
                    ss["FILEPATH"] = f"./openwq_in/{basename}"

    # Rewrite MODULE paths
    if "MODULES" in data:
        for mod_name, mod_cfg in data["MODULES"].items():
            if "MODULE_CONFIG_FILEPATH" in mod_cfg:
                basename = os.path.basename(mod_cfg["MODULE_CONFIG_FILEPATH"])
                mod_cfg["MODULE_CONFIG_FILEPATH"] = f"./openwq_in/{basename}"

    # Rewrite OUTPUT path
    if "OPENWQ_OUTPUT" in data:
        data["OPENWQ_OUTPUT"]["RESULTS_FOLDERPATH"] = "./openwq_out/"

    with open(dst_path, 'w') as f:
        json.dump(data, f, indent=4)

    return data


# ─────────────────────────────────────────────────────────────────────
# Rewrite text config files with Docker absolute paths
# ─────────────────────────────────────────────────────────────────────

def rewrite_text_config(src_path: str, dst_path: str, replacements: Dict[str, str]):
    """
    Read a text config file and apply string replacements.
    Used for fileManager.txt and mizuroute.control.
    Also strips unsupported tags like <calendar>.
    """
    with open(src_path, 'r') as f:
        content = f.read()

    for old, new in replacements.items():
        content = content.replace(old, new)

    # Remove unsupported tags (e.g., <calendar> not recognized by mizuroute_lakes_openwq)
    content = re.sub(r'<calendar>.*\n', '', content)

    with open(dst_path, 'w') as f:
        f.write(content)


def build_path_replacements(
    basin_id: str,
    variant_dir_name: str,
    docker_basin_dir: str,
    hpc_base: str = "/projects/F202407877CPCAA1/diogo/openwq/century_calibration",
) -> Dict[str, str]:
    """
    Build a dict of old->new path replacements for text config files.
    Replaces both HPC and local host paths with Docker container paths
    pointing into the docker_test workspace (not century_basins_prepared).
    """
    hpc_basin = f"{hpc_base}/basins/basin_{basin_id}/{variant_dir_name}"
    docker_variant = f"{docker_basin_dir}/{variant_dir_name}"

    # Local host paths in century_basins_prepared → Docker docker_test paths
    local_prepared = f"/Users/diogocosta/Documents/openwq_code/century_basins_prepared/basin_{basin_id}/{variant_dir_name}"

    return {
        # HPC variant path → Docker docker_test variant path
        hpc_basin: docker_variant,
        # Local host variant path → Docker docker_test variant path
        local_prepared: docker_variant,
        # HPC forcing path → Docker docker_test forcing path
        f"{hpc_base}/basins/basin_{basin_id}/forcing": f"{docker_basin_dir}/forcing",
        # Local host forcing path → Docker docker_test forcing path
        f"/Users/diogocosta/Documents/openwq_code/century_basins_prepared/basin_{basin_id}/forcing": f"{docker_basin_dir}/forcing",
    }


# ─────────────────────────────────────────────────────────────────────
# Generate Docker calibration config
# ─────────────────────────────────────────────────────────────────────

def generate_docker_calibration_config(
    basin_id: str,
    variant_key: str,
    variant_dir_name: str,
    docker_test_dir: str,
    prepared_dir: str,
    max_evaluations: int = 5,
) -> str:
    """Generate a Docker-ready calibration config by reading the HPC one
    and patching the runtime/path settings."""

    # Read the original calibration config to extract parameter definitions
    orig_config_path = os.path.join(
        prepared_dir, f"basin_{basin_id}", "calibration",
        f"calibration_config_{variant_key}.py"
    )

    if not os.path.exists(orig_config_path):
        print(f"  WARNING: No calibration config found: {orig_config_path}")
        return ""

    with open(orig_config_path, 'r') as f:
        orig_content = f.read()

    # Extract parameter list (between "calibration_parameters = [" and the closing "]")
    params_match = re.search(
        r'(calibration_parameters\s*=\s*\[.*?\])',
        orig_content, re.DOTALL
    )
    params_block = params_match.group(1) if params_match else "calibration_parameters = []"

    # Extract objective weights
    weights_match = re.search(r'objective_weights\s*=\s*(\{[^}]+\})', orig_content)
    weights = weights_match.group(1) if weights_match else '{"NO3-N": 1.0, "NH4-N": 0.5}'

    # Extract targets species
    targets_match = re.search(r'"species":\s*(\[[^\]]+\])', orig_content)
    targets = targets_match.group(1) if targets_match else '["NO3-N", "NH4-N"]'

    # Docker paths
    docker_basin_dir = host_to_docker(
        os.path.join(docker_test_dir, "basins", f"basin_{basin_id}")
    )
    docker_variant_dir = f"{docker_basin_dir}/{variant_dir_name}"
    docker_fm_path = f"{docker_variant_dir}/mizuroute_in/fileManager.txt"

    # Host paths (calibration driver runs on host Python)
    host_variant_dir = os.path.join(
        docker_test_dir, "basins", f"basin_{basin_id}", variant_dir_name
    )
    host_workspace = os.path.join(
        docker_test_dir, "calibration",
        f"workspace_{basin_id}_{variant_key}"
    )
    host_obs_path = os.path.join(
        docker_test_dir, "basins", f"basin_{basin_id}",
        "observations", "calibration_observations.csv"
    )

    docker_compose_path = str(
        SCRIPT_DIR.parents[2] / "containers" / "docker-compose.yml"
    )

    config_content = f'''#!/usr/bin/env python3
# Auto-generated Docker test calibration config
# Basin: {basin_id} - Variant {variant_key} ({variant_dir_name})
# Generated by 09_setup_docker_test.py

import sys
import os
import argparse

# Add calibration library to path
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CALIBRATION_LIB = "{CALIBRATION_LIB_HOST}"
if not os.path.isdir(CALIBRATION_LIB):
    CALIBRATION_LIB = os.environ.get("OPENWQ_CALIBRATION_LIB", CALIBRATION_LIB)
sys.path.insert(0, CALIBRATION_LIB)

from calibration_lib.calibration_driver import run_calibration, run_sensitivity_analysis
from calibration_lib.parameter_defaults import validate_parameter_bounds

# =============================================================================
# PATH CONFIGURATION (Docker local test)
# =============================================================================

BASIN_ID = "{basin_id}"
VARIANT = "{variant_key}"
VARIANT_LABEL = "{variant_dir_name}"

# Host paths (Python runs on host)
base_model_config_dir = "{host_variant_dir}"
calibration_work_dir = "{host_workspace}"
observation_data_path = "{host_obs_path}"
test_case_dir = "{host_variant_dir}"

# =============================================================================
# CONTAINER RUNTIME (Docker for local testing)
# =============================================================================

container_runtime = "docker"
docker_container_name = "docker_openwq"
docker_compose_path = "{docker_compose_path}"
executable_name = "{EXECUTABLE_NAME}"
executable_args = "-g 1 1"
file_manager_path = "{docker_fm_path}"

# Apptainer settings (unused in Docker mode)
apptainer_sif_path = ""
apptainer_bind_path = ""

# =============================================================================
# CALIBRATION PARAMETERS
# =============================================================================

{params_block}

# =============================================================================
# CALIBRATION SETTINGS (quick test: {max_evaluations} evaluations)
# =============================================================================

algorithm = "DDS"
max_evaluations = {max_evaluations}
n_parallel = 1
objective_function = "KGE"
temporal_resolution = "monthly"
aggregation_method = "mean"

objective_weights = {weights}

calibration_targets = {{
    "species": {targets},
    "reach_ids": "all",
    "compartments": ["RIVER_NETWORK_REACHES"],
}}

random_seed = 42

# =============================================================================
# SENSITIVITY ANALYSIS (skip for quick Docker test)
# =============================================================================

run_sensitivity_first = False
sensitivity_method = "morris"
sensitivity_morris_trajectories = 5
sensitivity_morris_levels = 4
sensitivity_sobol_samples = 256
sensitivity_threshold = 0.1

# =============================================================================
# HPC SETTINGS (disabled for Docker)
# =============================================================================

hpc_enabled = False
hpc_scheduler = "slurm"
hpc_partition = "standard"
hpc_walltime = "01:00:00"
hpc_nodes = 1
hpc_tasks_per_node = 4
hpc_memory = "8G"
hpc_max_concurrent_jobs = 1

# =============================================================================
# OBSERVATION DATA
# =============================================================================

observation_data_source = "csv"

# =============================================================================
# RUN CALIBRATION
# =============================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=f"Docker Test: {{BASIN_ID}} - Variant {{VARIANT}}"
    )
    parser.add_argument("--resume", action="store_true", help="Resume from checkpoint")
    parser.add_argument("--sensitivity-only", action="store_true", help="SA only")
    parser.add_argument("--dry-run", action="store_true", help="Validate config")
    args = parser.parse_args()

    config = {{
        "base_model_config_dir": base_model_config_dir,
        "calibration_work_dir": calibration_work_dir,
        "observation_data_path": observation_data_path,
        "test_case_dir": test_case_dir,
        "container_runtime": container_runtime,
        "docker_container_name": docker_container_name,
        "docker_compose_path": docker_compose_path,
        "apptainer_sif_path": apptainer_sif_path,
        "apptainer_bind_path": apptainer_bind_path,
        "executable_name": executable_name,
        "executable_args": executable_args,
        "file_manager_path": file_manager_path,
        "calibration_parameters": calibration_parameters,
        "algorithm": algorithm,
        "max_evaluations": max_evaluations,
        "n_parallel": n_parallel,
        "objective_function": objective_function,
        "objective_weights": objective_weights,
        "calibration_targets": calibration_targets,
        "random_seed": random_seed,
        "temporal_resolution": temporal_resolution,
        "aggregation_method": aggregation_method,
        "run_sensitivity_first": run_sensitivity_first,
        "sensitivity_method": sensitivity_method,
        "sensitivity_morris_trajectories": sensitivity_morris_trajectories,
        "sensitivity_morris_levels": sensitivity_morris_levels,
        "sensitivity_sobol_samples": sensitivity_sobol_samples,
        "sensitivity_threshold": sensitivity_threshold,
        "hpc_enabled": hpc_enabled,
        "hpc_scheduler": hpc_scheduler,
        "hpc_partition": hpc_partition,
        "hpc_walltime": hpc_walltime,
        "hpc_nodes": hpc_nodes,
        "hpc_tasks_per_node": hpc_tasks_per_node,
        "hpc_memory": hpc_memory,
        "hpc_max_concurrent_jobs": hpc_max_concurrent_jobs,
        "observation_data_source": observation_data_source,
    }}

    if args.dry_run:
        print("=" * 60)
        print(f"DRY RUN: {{BASIN_ID}} - Variant {{VARIANT}} ({{VARIANT_LABEL}})")
        print("=" * 60)
        print(f"Parameters: {{len(calibration_parameters)}}")
        for p in calibration_parameters:
            print(f"  - {{p['name']}}: {{p['bounds']}} ({{p['transform']}})")
        print(f"Container: Docker ({{docker_container_name}})")
        print(f"File manager: {{file_manager_path}}")
        print(f"Max evaluations: {{max_evaluations}}")
        sys.exit(0)

    if args.sensitivity_only:
        run_sensitivity_analysis(config)
    else:
        run_calibration(config, resume=args.resume)
'''

    # Write config
    output_path = os.path.join(
        docker_test_dir, "calibration",
        f"calibration_config_docker_{basin_id}_{variant_key}.py"
    )
    with open(output_path, 'w') as f:
        f.write(config_content)

    return output_path


# ─────────────────────────────────────────────────────────────────────
# Main: build Docker test workspace
# ─────────────────────────────────────────────────────────────────────

def setup_docker_test(
    prepared_dir: str,
    docker_test_dir: str,
    max_evaluations: int = 5,
):
    """Build the complete Docker test workspace."""

    print("=" * 70)
    print("DOCKER TEST WORKSPACE SETUP")
    print("=" * 70)
    print(f"Source:        {prepared_dir}")
    print(f"Docker test:   {docker_test_dir}")
    print(f"Test basins:   {', '.join(TEST_BASINS)}")
    print(f"Variants:      {', '.join(VARIANTS.keys())}")
    print(f"Max evals:     {max_evaluations}")
    print()

    os.makedirs(docker_test_dir, exist_ok=True)
    os.makedirs(os.path.join(docker_test_dir, "calibration"), exist_ok=True)

    total_configs = 0
    total_symlinks = 0
    total_rewrites = 0

    for basin_id in TEST_BASINS:
        src_basin_dir = os.path.join(prepared_dir, f"basin_{basin_id}")
        dst_basin_dir = os.path.join(docker_test_dir, "basins", f"basin_{basin_id}")

        if not os.path.isdir(src_basin_dir):
            print(f"  WARNING: Basin not found: {src_basin_dir}")
            continue

        print(f"\n[{basin_id}]")
        os.makedirs(dst_basin_dir, exist_ok=True)

        # Docker container path for this basin
        docker_basin_dir = host_to_docker(dst_basin_dir)

        # ── Symlink shared directories ──
        for shared_dir in ["forcing", "observations", "shapefiles"]:
            src = os.path.join(src_basin_dir, shared_dir)
            dst = os.path.join(dst_basin_dir, shared_dir)
            if os.path.exists(src):
                make_symlink(src, dst)
                total_symlinks += 1

        # ── Process each variant ──
        for v_key, v_dir_name in VARIANTS.items():
            src_variant = os.path.join(src_basin_dir, v_dir_name)
            dst_variant = os.path.join(dst_basin_dir, v_dir_name)
            docker_variant = f"{docker_basin_dir}/{v_dir_name}"

            if not os.path.isdir(src_variant):
                print(f"  [{v_key}] SKIP: {v_dir_name} not found")
                continue

            os.makedirs(dst_variant, exist_ok=True)

            # ── Symlink openwq_in/ (param_handler will copy to eval_dir) ──
            src_openwq_in = os.path.join(src_variant, "openwq_in")
            dst_openwq_in = os.path.join(dst_variant, "openwq_in")
            if os.path.exists(src_openwq_in):
                make_symlink(src_openwq_in, dst_openwq_in)
                total_symlinks += 1

            # ── Rewrite openWQ_master.json with relative paths ──
            src_master = os.path.join(src_variant, "openWQ_master.json")
            dst_master = os.path.join(dst_variant, "openWQ_master.json")
            if os.path.exists(src_master):
                rewrite_master_json_relative(src_master, dst_master)
                total_rewrites += 1

            # ── Copy and rewrite mizuroute_in/ ──
            src_miz_in = os.path.join(src_variant, "mizuroute_in")
            dst_miz_in = os.path.join(dst_variant, "mizuroute_in")
            if os.path.isdir(src_miz_in):
                # Create structure
                os.makedirs(dst_miz_in, exist_ok=True)
                os.makedirs(os.path.join(dst_miz_in, "ancillary_data"), exist_ok=True)
                os.makedirs(os.path.join(dst_miz_in, "settings"), exist_ok=True)
                os.makedirs(os.path.join(dst_miz_in, "input"), exist_ok=True)

                # Build path replacements
                replacements = build_path_replacements(
                    basin_id, v_dir_name, docker_basin_dir
                )

                # Rewrite fileManager.txt
                src_fm = os.path.join(src_miz_in, "fileManager.txt")
                dst_fm = os.path.join(dst_miz_in, "fileManager.txt")
                if os.path.exists(src_fm):
                    rewrite_text_config(src_fm, dst_fm, replacements)
                    total_rewrites += 1

                # Rewrite mizuroute.control
                src_ctrl = os.path.join(src_miz_in, "settings", "mizuroute.control")
                dst_ctrl = os.path.join(dst_miz_in, "settings", "mizuroute.control")
                if os.path.exists(src_ctrl):
                    rewrite_text_config(src_ctrl, dst_ctrl, replacements)
                    total_rewrites += 1

                # Symlink ancillary data files
                src_anc = os.path.join(src_miz_in, "ancillary_data")
                if os.path.isdir(src_anc):
                    for f in os.listdir(src_anc):
                        sf = os.path.join(src_anc, f)
                        df = os.path.join(dst_miz_in, "ancillary_data", f)
                        if os.path.isfile(sf):
                            make_symlink(sf, df)
                            total_symlinks += 1

                # Symlink input files
                src_inp = os.path.join(src_miz_in, "input")
                if os.path.isdir(src_inp):
                    for f in os.listdir(src_inp):
                        sf = os.path.join(src_inp, f)
                        df = os.path.join(dst_miz_in, "input", f)
                        if os.path.isfile(sf):
                            make_symlink(sf, df)
                            total_symlinks += 1

            # ── Copy and rewrite summa_in/ ──
            src_summa_in = os.path.join(src_variant, "summa_in")
            dst_summa_in = os.path.join(dst_variant, "summa_in")
            if os.path.isdir(src_summa_in):
                os.makedirs(dst_summa_in, exist_ok=True)

                replacements = build_path_replacements(
                    basin_id, v_dir_name, docker_basin_dir
                )

                # Rewrite SUMMA fileManager.txt
                src_sfm = os.path.join(src_summa_in, "fileManager.txt")
                dst_sfm = os.path.join(dst_summa_in, "fileManager.txt")
                if os.path.exists(src_sfm):
                    rewrite_text_config(src_sfm, dst_sfm, replacements)
                    total_rewrites += 1

                # Symlink all other SUMMA files (.nc, .txt, .TBL, etc.)
                for f in os.listdir(src_summa_in):
                    if f == "fileManager.txt":
                        continue  # Already handled
                    sf = os.path.join(src_summa_in, f)
                    df = os.path.join(dst_summa_in, f)
                    if os.path.isfile(sf):
                        make_symlink(sf, df)
                        total_symlinks += 1

            # ── Create simulation output directories ──
            sim_base = os.path.join(dst_variant, "simulations", "run_1")
            os.makedirs(os.path.join(sim_base, "SUMMA"), exist_ok=True)
            os.makedirs(os.path.join(sim_base, "mizuRoute"), exist_ok=True)

            # ── Create openwq_out/ directory ──
            os.makedirs(os.path.join(dst_variant, "openwq_out"), exist_ok=True)

            print(f"  [{v_key}] {v_dir_name} OK")

        # ── Generate Docker calibration configs ──
        for v_key, v_dir_name in VARIANTS.items():
            src_variant = os.path.join(src_basin_dir, v_dir_name)
            if not os.path.isdir(src_variant):
                continue
            cfg_path = generate_docker_calibration_config(
                basin_id=basin_id,
                variant_key=v_key,
                variant_dir_name=v_dir_name,
                docker_test_dir=docker_test_dir,
                prepared_dir=prepared_dir,
                max_evaluations=max_evaluations,
            )
            if cfg_path:
                total_configs += 1

    # ── Summary ──
    print("\n" + "=" * 70)
    print("DOCKER TEST WORKSPACE SUMMARY")
    print("=" * 70)
    print(f"Symlinks created:      {total_symlinks}")
    print(f"Config files rewritten: {total_rewrites}")
    print(f"Calibration configs:   {total_configs}")
    print(f"Workspace:             {docker_test_dir}")
    print()
    print("Next steps:")
    print(f"  1. cd {SCRIPT_DIR.parents[2]}/containers && docker compose up -d")
    print(f"  2. docker exec docker_openwq ls {EXECUTABLE_DIR}/")
    print(f"  3. Run a forward test:")
    print(f"     docker exec -e master_json={host_to_docker(docker_test_dir)}/basins/basin_{TEST_BASINS[0]}/variant_A_nitrogen_full/openWQ_master.json \\")
    print(f"       docker_openwq /bin/bash -c \\")
    print(f"       \"cd {host_to_docker(docker_test_dir)}/basins/basin_{TEST_BASINS[0]}/variant_A_nitrogen_full && \\")
    print(f"       {EXECUTABLE_DIR}/{EXECUTABLE_NAME} -g 1 1 \\")
    print(f"       -m {host_to_docker(docker_test_dir)}/basins/basin_{TEST_BASINS[0]}/variant_A_nitrogen_full/mizuroute_in/fileManager.txt\"")


def main():
    parser = argparse.ArgumentParser(
        description="Set up Docker test workspace for century basin calibration"
    )
    parser.add_argument(
        '--prepared-dir',
        default='/Users/diogocosta/Documents/openwq_code/century_basins_prepared',
        help='Source prepared basins directory'
    )
    parser.add_argument(
        '--docker-test-dir',
        default='/Users/diogocosta/Documents/openwq_code/docker_test',
        help='Output Docker test workspace directory'
    )
    parser.add_argument(
        '--max-evaluations', type=int, default=5,
        help='Max calibration evaluations for quick test (default: 5)'
    )

    args = parser.parse_args()

    setup_docker_test(
        prepared_dir=args.prepared_dir,
        docker_test_dir=args.docker_test_dir,
        max_evaluations=args.max_evaluations,
    )


if __name__ == "__main__":
    main()
