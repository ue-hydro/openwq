#!/bin/bash
# =============================================================================
# OpenWQ Calibration - SLURM Array Job Template
# =============================================================================
# Use this template for sensitivity analysis (Morris, Sobol) where all
# evaluations can run independently in parallel.
#
# Variables to set (passed as environment variables or replaced by coordinator):
#   $CALIB_ROOT     - Root calibration workspace directory
#   $INDEX_FILE     - JSON file mapping array index to eval directory
#   $SIF_PATH       - Path to Apptainer/Singularity image
#   $BIND_PATH      - Bind mount path (host:container)
#   $EXECUTABLE     - Model executable name
#   $EXEC_ARGS      - Executable arguments
#   $FILE_MANAGER   - Path to fileManager.txt
# =============================================================================

#SBATCH --job-name=openwq_calib_array
#SBATCH --array=0-__N_TASKS__%__MAX_CONCURRENT__
#SBATCH --partition=__PARTITION__
#SBATCH --time=__WALLTIME__
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=__TASKS_PER_NODE__
#SBATCH --mem=__MEMORY__
#SBATCH --output=__LOG_DIR__/array_%A_%a.out
#SBATCH --error=__LOG_DIR__/array_%A_%a.err
# Uncomment if needed:
# #SBATCH --account=__ACCOUNT__
# #SBATCH --qos=__QOS__

# =============================================================================
# Environment Setup
# =============================================================================
echo "=========================================="
echo "OpenWQ Calibration Array Job"
echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Node: $SLURMD_NODENAME"
echo "Start Time: $(date)"
echo ""

# Load required modules (adjust for your HPC system)
# module load singularity  # or apptainer
# module load python/3.9

# Set paths (these will be replaced by the coordinator)
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
TASK_ID=$SLURM_ARRAY_TASK_ID

# Get work directory from index file
# Requires jq or python to parse JSON
if command -v jq &> /dev/null; then
    WORK_DIR=$(jq -r ".\"$TASK_ID\".dir" "$INDEX_FILE")
    EVAL_ID=$(jq -r ".\"$TASK_ID\".eval_id" "$INDEX_FILE")
else
    # Fallback to python
    WORK_DIR=$(python3 -c "import json; d=json.load(open('$INDEX_FILE')); print(d['$TASK_ID']['dir'])")
    EVAL_ID=$(python3 -c "import json; d=json.load(open('$INDEX_FILE')); print(d['$TASK_ID']['eval_id'])")
fi

echo "Evaluation ID: $EVAL_ID"
echo "Work Directory: $WORK_DIR"

# Validate work directory exists
if [ ! -d "$WORK_DIR" ]; then
    echo "ERROR: Work directory does not exist: $WORK_DIR"
    exit 1
fi

cd "$WORK_DIR"

# =============================================================================
# Wait for Parameters (in case coordinator is still writing)
# =============================================================================
PARAMS_FILE="$WORK_DIR/parameters.json"
MAX_WAIT=300  # 5 minutes
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
cat "$PARAMS_FILE"
echo ""

# =============================================================================
# Run Model
# =============================================================================
echo "=========================================="
echo "Starting Model Execution"
echo "=========================================="

# Set OpenWQ environment variables
export master_json="$WORK_DIR/openWQ_master.json"

# Run with Apptainer
START_TIME=$(date +%s)

apptainer exec \
    --bind "$BIND_PATH" \
    --env master_json="$master_json" \
    "$SIF_PATH" \
    ./"$EXECUTABLE" $EXEC_ARGS -m "$FILE_MANAGER"

EXIT_CODE=$?

END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))

echo ""
echo "Model completed with exit code: $EXIT_CODE"
echo "Runtime: $RUNTIME seconds"

# =============================================================================
# Post-processing
# =============================================================================
echo "=========================================="
echo "Post-processing"
echo "=========================================="

# Create status file
if [ $EXIT_CODE -eq 0 ]; then
    echo "SUCCESS" > "$WORK_DIR/STATUS"
    touch "$WORK_DIR/COMPLETED"

    # Log completion time
    echo "{\"eval_id\": $EVAL_ID, \"exit_code\": 0, \"runtime_seconds\": $RUNTIME}" > "$WORK_DIR/completion.json"

    echo "Evaluation $EVAL_ID completed successfully"
else
    echo "FAILED" > "$WORK_DIR/STATUS"
    echo "{\"eval_id\": $EVAL_ID, \"exit_code\": $EXIT_CODE, \"runtime_seconds\": $RUNTIME}" > "$WORK_DIR/completion.json"

    echo "ERROR: Evaluation $EVAL_ID failed with exit code $EXIT_CODE"
fi

# List output files
echo ""
echo "Output files:"
ls -la "$WORK_DIR/openwq_out/" 2>/dev/null || echo "No output directory found"

echo ""
echo "=========================================="
echo "End Time: $(date)"
echo "=========================================="

exit $EXIT_CODE
