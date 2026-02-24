#!/bin/bash
# =============================================================================
# run_pilot.sh — Pilot calibration: 5 basins x 3 variants = 15 calibrations
# =============================================================================
#
# Run Morris sensitivity analysis + DDS calibration for 5 diverse basins,
# each with all 3 N-cycle variants, to validate setup and estimate runtimes.
#
# USAGE:
#   sbatch run_pilot.sh
#
# PREREQUISITES:
#   1. Upload prepared basins to HPC:
#      rsync -avz century_basins_prepared/ hpc:/scratch/user/century_calibration/basins/
#   2. Upload Apptainer image:
#      rsync -avz openwq.sif hpc:/scratch/user/century_calibration/
#   3. Upload scripts:
#      rsync -avz century_basins_scripts/ hpc:/scratch/user/century_calibration/scripts/
#
# =============================================================================

#SBATCH --job-name=openwq_pilot
#SBATCH --array=0-14%15
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --partition=standard
#SBATCH --output=slurm_logs/pilot_%A_%a.out
#SBATCH --error=slurm_logs/pilot_%A_%a.err

# =============================================================================
# CONFIGURATION — Edit these paths for your HPC setup
# =============================================================================

# Base directory on HPC
HPC_BASE="/scratch/${USER}/century_calibration"

# Apptainer image
SIF_PATH="${HPC_BASE}/openwq.sif"

# Python environment (must have pandas, numpy, geopandas)
PYTHON_ENV="${HPC_BASE}/env/bin/activate"

# Prepared basins directory
BASINS_DIR="${HPC_BASE}/basins"

# Scripts directory
SCRIPTS_DIR="${HPC_BASE}/scripts"

# =============================================================================
# PILOT BASIN SELECTION — 5 diverse basins
# Choose: 1 headwater, 2 meso, 2 macro; mix of CAN and USA
# =============================================================================

# Edit these basin IDs based on your feasibility_report.csv (basins with GRQA data)
PILOT_BASINS=(
    "CAN_01DG003_headwater"
    "CAN_02GC010_meso"
    "USA_01013500_meso"
    "CAN_05BB001_macro"
    "USA_02011400_macro"
)

# Variants
VARIANTS=("A" "B" "C")

# =============================================================================
# COMPUTE BASIN AND VARIANT FROM ARRAY INDEX
# =============================================================================

N_BASINS=${#PILOT_BASINS[@]}
N_VARIANTS=${#VARIANTS[@]}

BASIN_IDX=$((SLURM_ARRAY_TASK_ID / N_VARIANTS))
VARIANT_IDX=$((SLURM_ARRAY_TASK_ID % N_VARIANTS))

BASIN_ID=${PILOT_BASINS[$BASIN_IDX]}
VARIANT=${VARIANTS[$VARIANT_IDX]}

echo "========================================"
echo "PILOT CALIBRATION"
echo "========================================"
echo "Array task:  ${SLURM_ARRAY_TASK_ID}"
echo "Basin:       ${BASIN_ID}"
echo "Variant:     ${VARIANT}"
echo "Start time:  $(date)"
echo "Node:        $(hostname)"
echo "========================================"

# Create log directory
mkdir -p slurm_logs

# Activate Python environment
if [ -f "${PYTHON_ENV}" ]; then
    source "${PYTHON_ENV}"
fi

# =============================================================================
# RUN CALIBRATION
# =============================================================================

BASIN_DIR="${BASINS_DIR}/basin_${BASIN_ID}"
CALIB_CONFIG="${BASIN_DIR}/calibration/calibration_config_${VARIANT}.py"

if [ ! -f "${CALIB_CONFIG}" ]; then
    echo "ERROR: Calibration config not found: ${CALIB_CONFIG}"
    exit 1
fi

# Run with sensitivity analysis first (Morris), then DDS
echo ""
echo "Running: python ${CALIB_CONFIG}"
python "${CALIB_CONFIG}"
EXIT_CODE=$?

echo ""
echo "========================================"
echo "COMPLETED"
echo "========================================"
echo "Exit code:   ${EXIT_CODE}"
echo "End time:    $(date)"
echo "========================================"

# Create completion marker
if [ ${EXIT_CODE} -eq 0 ]; then
    touch "${BASIN_DIR}/calibration_workspace_${VARIANT}/COMPLETED"
else
    echo "${EXIT_CODE}" > "${BASIN_DIR}/calibration_workspace_${VARIANT}/FAILED"
fi

exit ${EXIT_CODE}
