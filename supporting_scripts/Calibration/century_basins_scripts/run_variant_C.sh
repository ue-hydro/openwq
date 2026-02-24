#!/bin/bash
# =============================================================================
# run_variant_C.sh — Production calibration: Variant C (Thermodynamic N Cycle)
# =============================================================================
#
# 111 basins x 450 DDS evaluations
# 16 parameters (9 BGC + 7 SS)
# Template: Combined nitrification + denitrification_thermodynamic + anammox
# Uses F_T (thermodynamic potential factor) from Jin & Bethke (2003, 2007)
#
# USAGE:
#   sbatch run_variant_C.sh
#   sbatch --array=0-50 run_variant_C.sh   # Only first 51 basins
#
# =============================================================================

#SBATCH --job-name=openwq_C
#SBATCH --array=0-110%50
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=60:00:00
#SBATCH --partition=standard
#SBATCH --output=slurm_logs/variant_C_%A_%a.out
#SBATCH --error=slurm_logs/variant_C_%A_%a.err

# =============================================================================
# CONFIGURATION
# =============================================================================

HPC_BASE="/scratch/${USER}/century_calibration"
SIF_PATH="${HPC_BASE}/openwq.sif"
PYTHON_ENV="${HPC_BASE}/env/bin/activate"
BASINS_DIR="${HPC_BASE}/basins"
VARIANT="C"

# =============================================================================
# DETERMINE BASIN FROM ARRAY INDEX
# =============================================================================

mapfile -t BASIN_LIST < <(ls -d "${BASINS_DIR}"/basin_* 2>/dev/null | sort | while read d; do
    basename "$d" | sed 's/^basin_//'
done)

if [ ${#BASIN_LIST[@]} -eq 0 ]; then
    echo "ERROR: No basin directories found in ${BASINS_DIR}"
    exit 1
fi

if [ ${SLURM_ARRAY_TASK_ID} -ge ${#BASIN_LIST[@]} ]; then
    echo "Array task ${SLURM_ARRAY_TASK_ID} exceeds number of basins (${#BASIN_LIST[@]})"
    exit 0
fi

BASIN_ID=${BASIN_LIST[$SLURM_ARRAY_TASK_ID]}

echo "========================================"
echo "VARIANT C: Thermodynamic N Cycle"
echo "========================================"
echo "Array task:  ${SLURM_ARRAY_TASK_ID}"
echo "Basin:       ${BASIN_ID}"
echo "Max evals:   450"
echo "Parameters:  16"
echo "Start time:  $(date)"
echo "Node:        $(hostname)"
echo "Memory:      ${SLURM_MEM_PER_NODE}"
echo "CPUs:        ${SLURM_CPUS_PER_TASK}"
echo "========================================"

mkdir -p slurm_logs

if [ -f "${PYTHON_ENV}" ]; then
    source "${PYTHON_ENV}"
fi

# =============================================================================
# CHECK COMPLETION / RESUME
# =============================================================================

BASIN_DIR="${BASINS_DIR}/basin_${BASIN_ID}"
WORKSPACE="${BASIN_DIR}/calibration_workspace_${VARIANT}"

if [ -f "${WORKSPACE}/COMPLETED" ]; then
    echo "Basin ${BASIN_ID} Variant ${VARIANT} already completed. Skipping."
    exit 0
fi

OBS_FILE="${BASIN_DIR}/observations/calibration_observations.csv"
if [ ! -f "${OBS_FILE}" ]; then
    echo "WARNING: No observations file for ${BASIN_ID}. Skipping."
    mkdir -p "${WORKSPACE}"
    echo "NO_OBSERVATIONS" > "${WORKSPACE}/SKIPPED"
    exit 0
fi

# =============================================================================
# RUN CALIBRATION
# =============================================================================

CALIB_CONFIG="${BASIN_DIR}/calibration/calibration_config_${VARIANT}.py"

if [ ! -f "${CALIB_CONFIG}" ]; then
    echo "ERROR: Config not found: ${CALIB_CONFIG}"
    exit 1
fi

RESUME_FLAG=""
if [ -f "${WORKSPACE}/calibration_state.json" ]; then
    echo "Found checkpoint, resuming..."
    RESUME_FLAG="--resume"
fi

echo ""
echo "Running: python ${CALIB_CONFIG} ${RESUME_FLAG}"
python "${CALIB_CONFIG}" ${RESUME_FLAG}
EXIT_CODE=$?

echo ""
echo "========================================"
echo "COMPLETED"
echo "========================================"
echo "Exit code:   ${EXIT_CODE}"
echo "End time:    $(date)"
echo "========================================"

if [ ${EXIT_CODE} -eq 0 ]; then
    touch "${WORKSPACE}/COMPLETED"
else
    echo "${EXIT_CODE}" > "${WORKSPACE}/FAILED"
fi

exit ${EXIT_CODE}
