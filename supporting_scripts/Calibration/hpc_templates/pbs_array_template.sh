#!/bin/bash
# =============================================================================
# OpenWQ Calibration - PBS Array Job Template
# =============================================================================
# Use this template for PBS/Torque job schedulers.
# Similar to SLURM array template but with PBS directives.
#
# Variables to set:
#   $CALIB_ROOT     - Root calibration workspace directory
#   $INDEX_FILE     - JSON file mapping array index to eval directory
#   $SIF_PATH       - Path to Apptainer/Singularity image
#   $BIND_PATH      - Bind mount path (host:container)
#   $EXECUTABLE     - Model executable name
#   $EXEC_ARGS      - Executable arguments
#   $FILE_MANAGER   - Path to fileManager.txt
# =============================================================================

#PBS -N openwq_calib
#PBS -J 0-__N_TASKS__
#PBS -q __QUEUE__
#PBS -l walltime=__WALLTIME__
#PBS -l nodes=1:ppn=__PPN__
#PBS -l mem=__MEMORY__
#PBS -o __LOG_DIR__/pbs_^array_index^.out
#PBS -e __LOG_DIR__/pbs_^array_index^.err
# Uncomment if needed:
# #PBS -A __ACCOUNT__

# =============================================================================
# Environment Setup
# =============================================================================
echo "=========================================="
echo "OpenWQ Calibration PBS Array Job"
echo "=========================================="
echo "Job ID: $PBS_JOBID"
echo "Array Index: $PBS_ARRAY_INDEX"
echo "Node: $(hostname)"
echo "Start Time: $(date)"
echo ""

# Load modules (adjust for your system)
# module load singularity
# module load python/3.9

# Set paths
CALIB_ROOT="__CALIB_ROOT__"
INDEX_FILE="__INDEX_FILE__"
SIF_PATH="__SIF_PATH__"
BIND_PATH="__BIND_PATH__"
EXECUTABLE="__EXECUTABLE__"
EXEC_ARGS="__EXEC_ARGS__"
FILE_MANAGER="__FILE_MANAGER__"

# =============================================================================
# Parse Task Information
# =============================================================================
TASK_ID=$PBS_ARRAY_INDEX

# Get work directory from index file
if command -v jq &> /dev/null; then
    WORK_DIR=$(jq -r ".\"$TASK_ID\".dir" "$INDEX_FILE")
    EVAL_ID=$(jq -r ".\"$TASK_ID\".eval_id" "$INDEX_FILE")
else
    WORK_DIR=$(python3 -c "import json; d=json.load(open('$INDEX_FILE')); print(d['$TASK_ID']['dir'])")
    EVAL_ID=$(python3 -c "import json; d=json.load(open('$INDEX_FILE')); print(d['$TASK_ID']['eval_id'])")
fi

echo "Evaluation ID: $EVAL_ID"
echo "Work Directory: $WORK_DIR"

# Validate work directory
if [ ! -d "$WORK_DIR" ]; then
    echo "ERROR: Work directory does not exist: $WORK_DIR"
    exit 1
fi

cd "$WORK_DIR"

# =============================================================================
# Wait for Parameters
# =============================================================================
PARAMS_FILE="$WORK_DIR/parameters.json"
MAX_WAIT=300
WAITED=0

while [ ! -f "$PARAMS_FILE" ]; do
    echo "Waiting for parameters file..."
    sleep 10
    WAITED=$((WAITED + 10))
    if [ $WAITED -ge $MAX_WAIT ]; then
        echo "ERROR: Timeout waiting for parameters file"
        exit 1
    fi
done

echo "Parameters file found"

# =============================================================================
# Run Model
# =============================================================================
echo "=========================================="
echo "Starting Model Execution"
echo "=========================================="

export master_json="$WORK_DIR/openWQ_master.json"

START_TIME=$(date +%s)

# Use singularity if apptainer not available
if command -v apptainer &> /dev/null; then
    apptainer exec \
        --bind "$BIND_PATH" \
        --env master_json="$master_json" \
        "$SIF_PATH" \
        ./"$EXECUTABLE" $EXEC_ARGS -m "$FILE_MANAGER"
else
    singularity exec \
        --bind "$BIND_PATH" \
        --env master_json="$master_json" \
        "$SIF_PATH" \
        ./"$EXECUTABLE" $EXEC_ARGS -m "$FILE_MANAGER"
fi

EXIT_CODE=$?

END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))

echo ""
echo "Model completed with exit code: $EXIT_CODE"
echo "Runtime: $RUNTIME seconds"

# =============================================================================
# Post-processing
# =============================================================================
if [ $EXIT_CODE -eq 0 ]; then
    echo "SUCCESS" > "$WORK_DIR/STATUS"
    touch "$WORK_DIR/COMPLETED"
    echo "{\"eval_id\": $EVAL_ID, \"exit_code\": 0, \"runtime_seconds\": $RUNTIME}" > "$WORK_DIR/completion.json"
    echo "Evaluation $EVAL_ID completed successfully"
else
    echo "FAILED" > "$WORK_DIR/STATUS"
    echo "{\"eval_id\": $EVAL_ID, \"exit_code\": $EXIT_CODE, \"runtime_seconds\": $RUNTIME}" > "$WORK_DIR/completion.json"
    echo "ERROR: Evaluation $EVAL_ID failed"
fi

echo ""
echo "=========================================="
echo "End Time: $(date)"
echo "=========================================="

exit $EXIT_CODE
