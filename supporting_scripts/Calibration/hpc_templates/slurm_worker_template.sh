#!/bin/bash
# =============================================================================
# OpenWQ Calibration - SLURM Worker Template
# =============================================================================
# Use this template for DDS optimization where a worker pool processes
# tasks adaptively based on optimization algorithm decisions.
#
# Workers poll a task queue directory for new work, execute model runs,
# and report results back to the coordinator.
#
# Variables to set:
#   $CALIB_ROOT     - Root calibration workspace directory
#   $WORKER_ID      - Unique worker identifier
#   $TASK_QUEUE     - Directory where coordinator posts new tasks
#   $SIF_PATH       - Path to Apptainer/Singularity image
#   $BIND_PATH      - Bind mount path (host:container)
#   $EXECUTABLE     - Model executable name
#   $EXEC_ARGS      - Executable arguments
#   $FILE_MANAGER   - Path to fileManager.txt
# =============================================================================

#SBATCH --job-name=openwq_worker
#SBATCH --partition=__PARTITION__
#SBATCH --time=__WALLTIME__
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=__TASKS_PER_NODE__
#SBATCH --mem=__MEMORY__
#SBATCH --output=__LOG_DIR__/worker_%j.out
#SBATCH --error=__LOG_DIR__/worker_%j.err
# Uncomment if needed:
# #SBATCH --account=__ACCOUNT__
# #SBATCH --qos=__QOS__

# =============================================================================
# Environment Setup
# =============================================================================
echo "=========================================="
echo "OpenWQ Calibration Worker"
echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURMD_NODENAME"
echo "Start Time: $(date)"
echo ""

# Load required modules
# module load singularity
# module load python/3.9

# Set paths
CALIB_ROOT="__CALIB_ROOT__"
TASK_QUEUE="$CALIB_ROOT/task_queue"
RESULTS_DIR="$CALIB_ROOT/results_queue"
SIF_PATH="__SIF_PATH__"
BIND_PATH="__BIND_PATH__"
EXECUTABLE="__EXECUTABLE__"
EXEC_ARGS="__EXEC_ARGS__"
FILE_MANAGER="__FILE_MANAGER__"

# Worker identification
WORKER_ID="worker_${SLURM_JOB_ID}_$(hostname)"

# Polling settings
POLL_INTERVAL=10        # seconds between task checks
MAX_IDLE_TIME=1800      # 30 minutes idle before shutdown
HEARTBEAT_INTERVAL=60   # seconds between heartbeats

# Create worker directories
mkdir -p "$TASK_QUEUE/claimed"
mkdir -p "$RESULTS_DIR"

# Register worker
WORKER_FILE="$CALIB_ROOT/workers/$WORKER_ID"
mkdir -p "$(dirname $WORKER_FILE)"
echo "{\"worker_id\": \"$WORKER_ID\", \"job_id\": \"$SLURM_JOB_ID\", \"node\": \"$SLURMD_NODENAME\", \"start_time\": $(date +%s)}" > "$WORKER_FILE"

echo "Worker ID: $WORKER_ID"
echo "Task Queue: $TASK_QUEUE"
echo ""

# =============================================================================
# Worker Functions
# =============================================================================

send_heartbeat() {
    echo "{\"worker_id\": \"$WORKER_ID\", \"timestamp\": $(date +%s), \"status\": \"$1\"}" > "$WORKER_FILE"
}

claim_task() {
    # Try to atomically claim a task
    local task_file="$1"
    local claimed_file="$TASK_QUEUE/claimed/$(basename $task_file)"

    # Use mv for atomic operation
    if mv "$task_file" "$claimed_file" 2>/dev/null; then
        echo "$claimed_file"
        return 0
    fi
    return 1
}

run_evaluation() {
    local task_file="$1"

    # Parse task
    if command -v jq &> /dev/null; then
        EVAL_ID=$(jq -r '.eval_id' "$task_file")
        WORK_DIR=$(jq -r '.work_dir' "$task_file")
    else
        EVAL_ID=$(python3 -c "import json; print(json.load(open('$task_file'))['eval_id'])")
        WORK_DIR=$(python3 -c "import json; print(json.load(open('$task_file'))['work_dir'])")
    fi

    echo "Running evaluation $EVAL_ID in $WORK_DIR"

    cd "$WORK_DIR"

    # Set environment
    export master_json="$WORK_DIR/openWQ_master.json"

    # Run model
    START_TIME=$(date +%s)

    apptainer exec \
        --bind "$BIND_PATH" \
        --env master_json="$master_json" \
        "$SIF_PATH" \
        ./"$EXECUTABLE" $EXEC_ARGS -m "$FILE_MANAGER"

    EXIT_CODE=$?

    END_TIME=$(date +%s)
    RUNTIME=$((END_TIME - START_TIME))

    # Report result
    if [ $EXIT_CODE -eq 0 ]; then
        echo "SUCCESS" > "$WORK_DIR/STATUS"
        touch "$WORK_DIR/COMPLETED"
        STATUS="completed"
    else
        echo "FAILED" > "$WORK_DIR/STATUS"
        STATUS="failed"
    fi

    # Write result file
    RESULT_FILE="$RESULTS_DIR/result_${EVAL_ID}.json"
    echo "{\"eval_id\": $EVAL_ID, \"worker_id\": \"$WORKER_ID\", \"status\": \"$STATUS\", \"exit_code\": $EXIT_CODE, \"runtime_seconds\": $RUNTIME}" > "$RESULT_FILE"

    # Clean up claimed task
    rm -f "$task_file"

    return $EXIT_CODE
}

# =============================================================================
# Main Worker Loop
# =============================================================================
IDLE_TIME=0
LAST_HEARTBEAT=$(date +%s)
TASKS_COMPLETED=0

echo "Starting worker loop..."
echo ""

while true; do
    # Check for shutdown signal
    if [ -f "$CALIB_ROOT/SHUTDOWN" ]; then
        echo "Shutdown signal received"
        break
    fi

    # Check for calibration complete signal
    if [ -f "$CALIB_ROOT/CALIBRATION_COMPLETE" ]; then
        echo "Calibration complete"
        break
    fi

    # Send heartbeat if needed
    NOW=$(date +%s)
    if [ $((NOW - LAST_HEARTBEAT)) -ge $HEARTBEAT_INTERVAL ]; then
        send_heartbeat "idle"
        LAST_HEARTBEAT=$NOW
    fi

    # Look for new tasks
    TASK_FILES=$(ls "$TASK_QUEUE"/task_*.json 2>/dev/null | head -1)

    if [ -n "$TASK_FILES" ]; then
        # Try to claim a task
        CLAIMED=$(claim_task "$TASK_FILES")

        if [ -n "$CLAIMED" ]; then
            IDLE_TIME=0
            send_heartbeat "running"

            echo "----------------------------------------"
            echo "Claimed task: $CLAIMED"

            run_evaluation "$CLAIMED"
            TASK_EXIT=$?

            TASKS_COMPLETED=$((TASKS_COMPLETED + 1))

            if [ $TASK_EXIT -eq 0 ]; then
                echo "Task completed successfully"
            else
                echo "Task failed with exit code $TASK_EXIT"
            fi

            echo "Total tasks completed: $TASKS_COMPLETED"
            echo ""
        fi
    else
        # No tasks available
        sleep $POLL_INTERVAL
        IDLE_TIME=$((IDLE_TIME + POLL_INTERVAL))

        if [ $IDLE_TIME -ge $MAX_IDLE_TIME ]; then
            echo "Max idle time reached, shutting down"
            break
        fi
    fi
done

# =============================================================================
# Cleanup
# =============================================================================
echo ""
echo "=========================================="
echo "Worker Shutdown"
echo "=========================================="
echo "Tasks Completed: $TASKS_COMPLETED"
echo "End Time: $(date)"

# Unregister worker
send_heartbeat "stopped"

echo "=========================================="
exit 0
