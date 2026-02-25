#!/bin/bash

#  Author: Anna Su
#  Date: 2026-02-20

#SBATCH --partition=pi_gerstein_gpu
#SBATCH --job-name=boltz_guided_gen
#SBATCH -c 8                  # Default CPUs; overridden by submit_all_jobs.sh --cpus-per-task
#SBATCH --mem=16G             # Default memory; overridden by submit_all_jobs.sh --mem
#SBATCH --gres=gpu:1
#SBATCH --output=/dev/null    # Default for standalone use; overridden by submit_all_jobs.sh --output
#SBATCH --error=/dev/null     # Default for standalone use; overridden by submit_all_jobs.sh --error
#SBATCH --mail-user=4752279178@vtext.com
#SBATCH --mail-type=ALL

# Note: --output and --error are set by submit_all_jobs.sh to point to the run log directory.
#       The #SBATCH defaults above (/dev/null) prevent stray slurm-*.out files when using standalone sbatch.

# ====================
# Parse Command Line Arguments
# ====================
# These are passed via sbatch --export or command line
# Usage: sbatch --export=ALL,ESM_API_TOKEN='your_token',SMILES='...',PDB_ID='...',CHAIN_ID='...' guided_generation.sh
# Or parameters can be passed as arguments: sbatch guided_generation.sh "SMILES" "PDB_ID" "CHAIN_ID" [MASKING_PCT] [NUM_DECODING_STEPS] [NUM_SAMPLES] [ESM_MODEL] [CACHE_DIR] [SCORE_TYPE] [USE_MSA_SERVER]

# Generate timestamp for this run
# RUN_DIR may already be set by submit_all_jobs.sh; if so, reuse it
# to keep .out/.err/.txt in the same directory.
if [ -z "$RUN_DIR" ]; then
    RUN_TIMESTAMP=$(date +%Y%m%d_%H%M%S)
else
    # Extract timestamp from existing RUN_DIR for other uses
    RUN_TIMESTAMP=$(echo "$RUN_DIR" | grep -oP '\d{8}_\d{6}$')
fi

# Check for command line arguments first, then environment variables
if [ $# -ge 3 ]; then
    SMILES="$1"
    PDB_ID="$2"
    CHAIN_ID="$3"
    MASKING_PCT="${4:-0.4}"
    NUM_DECODING_STEPS="${5:-32}"
    NUM_SAMPLES="${6:-10}"
    ESM_MODEL="${7:-esm3-medium-2024-08}"
    CACHE_DIR="${8:-}"
    SCORE_TYPE="${9:-affinity}"
    USE_MSA_SERVER="${10:-false}"
    TIME_LIMIT="${11:-15:00:00}"
elif [ -n "$SMILES" ] && [ -n "$PDB_ID" ] && [ -n "$CHAIN_ID" ]; then
    # Use environment variables (set via --export)
    MASKING_PCT="${MASKING_PCT:-0.4}"
    NUM_DECODING_STEPS="${NUM_DECODING_STEPS:-32}"
    NUM_SAMPLES="${NUM_SAMPLES:-10}"
    ESM_MODEL="${ESM_MODEL:-esm3-medium-2024-08}"
    CACHE_DIR="${CACHE_DIR:-}"
    SCORE_TYPE="${SCORE_TYPE:-affinity}"
    USE_MSA_SERVER="${USE_MSA_SERVER:-false}"
    TIME_LIMIT="${TIME_LIMIT:-15:00:00}"
else
    echo "ERROR: Missing required parameters!"
    echo "Usage (environment variables):"
    echo "  sbatch --export=ALL,ESM_API_TOKEN='your_token',SMILES='...',PDB_ID='...',CHAIN_ID='...' guided_generation.sh"
    echo ""
    echo "Usage (command line arguments):"
    echo "  sbatch guided_generation.sh SMILES PDB_ID CHAIN_ID [MASKING_PCT] [NUM_DECODING_STEPS] [NUM_SAMPLES] [ESM_MODEL] [CACHE_DIR] [SCORE_TYPE] [USE_MSA_SERVER] [TIME_LIMIT]"
    echo ""
    echo "Required: SMILES, PDB_ID, CHAIN_ID, ESM_API_TOKEN (via --export)"
    echo "Optional: MASKING_PCT (0.4), NUM_DECODING_STEPS (32), NUM_SAMPLES (10), ESM_MODEL (esm3-medium-2024-08), CACHE_DIR (auto), SCORE_TYPE (affinity), USE_MSA_SERVER (false), TIME_LIMIT (15:00:00)"
    exit 1
fi

# Now that PDB_ID is known, finalize RUN_DIR if not already set
if [ -z "$RUN_DIR" ]; then
    RUN_DIR="run_${PDB_ID}_${RUN_TIMESTAMP}"
fi

echo "Job Time Limit: $TIME_LIMIT"

# ====================
# Environment Setup
# ====================
# Strip whitespace from ESM_API_TOKEN to prevent "Illegal header value" errors
if [ -n "$ESM_API_TOKEN" ]; then
    ESM_API_TOKEN=$(echo "$ESM_API_TOKEN" | tr -d '[:space:]')
fi

# Validate ESM_API_TOKEN is set (required for Forge API)
if [ "$ESM_MODEL" != "local" ] && [ -z "$ESM_API_TOKEN" ]; then
    echo "ERROR: ESM_API_TOKEN is required for model: $ESM_MODEL"
    echo "Please set it via: sbatch --export=ALL,ESM_API_TOKEN='your_token',..."
    exit 1
fi

module load miniconda

# Source conda for non-interactive shell before activating
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate boltz_esm

# ====================
# Directory Setup
# ====================
# Use SLURM_SUBMIT_DIR (where sbatch was called from)
SCRIPT_DIR="${SLURM_SUBMIT_DIR}"
PROJECT_DIR="$(dirname "$(dirname "$SCRIPT_DIR")")"

# Create run-specific directories
RUN_LOG_DIR="${PROJECT_DIR}/logs/boltz/${RUN_DIR}"
mkdir -p "$RUN_LOG_DIR"

# Set output and error log paths in the run directory
OUTPUT_LOG="${RUN_LOG_DIR}/boltz_${PDB_ID}_${CHAIN_ID}_${SLURM_JOB_ID}.out"
ERROR_LOG="${RUN_LOG_DIR}/boltz_${PDB_ID}_${CHAIN_ID}_${SLURM_JOB_ID}.err"

# Redirect stdout and stderr to run-specific logs
# Use >> (append) so that if SLURM's --output already points to the same file,
# early output is preserved rather than truncated.
exec 1>>"$OUTPUT_LOG"
exec 2>>"$ERROR_LOG"

# ====================
# GPU Logging Setup
# ====================
GPU_LOG="${RUN_LOG_DIR}/gpu_usage_${SLURM_JOB_ID}.csv"

# ====================
# CPU Logging Setup
# ====================
CPU_LOG="${RUN_LOG_DIR}/cpu_usage_${SLURM_JOB_ID}.csv"

# Print job info
echo "=============================================="
echo "SLURM Job ID: $SLURM_JOB_ID"
echo "Node: $SLURMD_NODENAME"
echo "Start Time: $(date)"
echo "Run Directory: $RUN_DIR"
echo "Working Directory: $SCRIPT_DIR"
echo "Output Log: $OUTPUT_LOG"
echo "Error Log: $ERROR_LOG"
echo "GPU Log: $GPU_LOG"
echo "CPU Log: $CPU_LOG"
echo "=============================================="
echo "Job Parameters:"
echo "  PDB ID: $PDB_ID"
echo "  Chain ID: $CHAIN_ID"
echo "  SMILES: $SMILES"
echo "  Masking %: $MASKING_PCT"
echo "  Decoding Steps: $NUM_DECODING_STEPS"
echo "  Samples per Step: $NUM_SAMPLES"
echo "  ESM Model: $ESM_MODEL"
echo "  Cache Dir: ${CACHE_DIR:-auto}"
echo "  Score Type: $SCORE_TYPE"
echo "  Use MSA Server: $USE_MSA_SERVER"
echo "=============================================="

# Print GPU info
echo "GPU Information:"
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader
echo "=============================================="

# Start GPU monitoring in background
# Logs: timestamp, GPU utilization %, memory used MB, memory total MB, temperature
echo "timestamp,gpu_util_percent,memory_used_mb,memory_total_mb,temperature_c,power_draw_w" > "$GPU_LOG"

nvidia-smi --query-gpu=timestamp,utilization.gpu,memory.used,memory.total,temperature.gpu,power.draw \
    --format=csv,noheader,nounits -l 30 >> "$GPU_LOG" &
GPU_MONITOR_PID=$!

echo "Started GPU monitoring (PID: $GPU_MONITOR_PID)"

# Start CPU monitoring in background
# Logs: timestamp, CPU utilization %, memory used MB, memory total MB, load average (1min)
echo "timestamp,cpu_util_percent,memory_used_mb,memory_total_mb,load_avg_1min" > "$CPU_LOG"

# Monitor CPU every 30 seconds using top and /proc/meminfo
(
    while true; do
        TIMESTAMP=$(date +"%Y-%m-%d %H:%M:%S")
        # Get CPU usage (100 - idle%)
        CPU_IDLE=$(top -bn2 -d 0.5 | grep "Cpu(s)" | tail -1 | awk '{print $8}' | cut -d'%' -f1)
        CPU_UTIL=$(awk "BEGIN {printf \"%.1f\", 100 - $CPU_IDLE}")
        # Get memory info in MB
        MEM_INFO=$(free -m | grep Mem:)
        MEM_TOTAL=$(echo $MEM_INFO | awk '{print $2}')
        MEM_USED=$(echo $MEM_INFO | awk '{print $3}')
        # Get load average (1 minute)
        LOAD_AVG=$(uptime | awk -F'load average:' '{print $2}' | awk -F',' '{print $1}' | xargs)
        
        echo "$TIMESTAMP,$CPU_UTIL,$MEM_USED,$MEM_TOTAL,$LOAD_AVG" >> "$CPU_LOG"
        sleep 30
    done
) &
CPU_MONITOR_PID=$!

echo "Started CPU monitoring (PID: $CPU_MONITOR_PID)"

# ====================
# Main Execution
# ====================
echo ""
echo "Starting Guided Generation..."
echo ""

# Build python command with required arguments
# Pass --run_dir so that main.py places .txt log in the same directory as .out/.err
PYTHON_CMD="python ${SCRIPT_DIR}/main.py \\
    --smiles \"$SMILES\" \\
    --wildtype \"$PDB_ID\" \"$CHAIN_ID\" \\
    --masking_percentage \"$MASKING_PCT\" \\
    --num_decoding_steps \"$NUM_DECODING_STEPS\" \\
    --num_samples_per_step \"$NUM_SAMPLES\" \\
    --model \"$ESM_MODEL\" \\
    --score_type \"$SCORE_TYPE\" \\
    --run_dir \"$RUN_DIR\""

# Add forge token if using API model
if [ "$ESM_MODEL" != "local" ]; then
    PYTHON_CMD="$PYTHON_CMD \\
    --forge_token \"$ESM_API_TOKEN\""
fi

if [ -n "$CACHE_DIR" ]; then
    PYTHON_CMD="$PYTHON_CMD \\
    --cache_dir \"$CACHE_DIR\""
fi

if [ "$USE_MSA_SERVER" = "true" ]; then
    PYTHON_CMD="$PYTHON_CMD \\
    --use_msa_server"
fi

echo "Running command:"
echo "$PYTHON_CMD"
echo ""

# Execute the command
eval "$PYTHON_CMD"

PYTHON_EXIT_CODE=$?

# ====================
# Cleanup & Summary
# ====================
# Stop monitoring processes
kill $GPU_MONITOR_PID 2>/dev/null
wait $GPU_MONITOR_PID 2>/dev/null

kill $CPU_MONITOR_PID 2>/dev/null
wait $CPU_MONITOR_PID 2>/dev/null

echo ""
echo "=============================================="
echo "Job Complete"
echo "End Time: $(date)"
echo "Python Exit Code: $PYTHON_EXIT_CODE"
echo "=============================================="

# ====================
# Generate Resource Usage Plots
# ====================
echo ""
echo "Generating resource usage plots..."

RESOURCE_PLOTS_DIR="${PROJECT_DIR}/results/boltz/${RUN_DIR}"
mkdir -p "$RESOURCE_PLOTS_DIR"

# Run plotting script
python "${SCRIPT_DIR}/plot_resource_usage.py" "$GPU_LOG" "$CPU_LOG" "$RESOURCE_PLOTS_DIR"
PLOT_EXIT_CODE=$?

if [ $PLOT_EXIT_CODE -eq 0 ]; then
    echo "✓ Resource plots saved to: $RESOURCE_PLOTS_DIR"
else
    echo "✗ Failed to generate resource plots (exit code: $PLOT_EXIT_CODE)"
fi

# ====================
# Resource Usage Summary
# ===================="

# Generate GPU usage summary
if [ -f "$GPU_LOG" ]; then
    echo ""
    echo "GPU Usage Summary:"
    echo "------------------"
    # Skip header, calculate stats
    tail -n +2 "$GPU_LOG" | awk -F',' '
    BEGIN { 
        count=0; sum_util=0; sum_mem=0; max_util=0; max_mem=0 
    }
    {
        count++
        sum_util += $2
        sum_mem += $3
        if ($2 > max_util) max_util = $2
        if ($3 > max_mem) max_mem = $3
    }
    END {
        if (count > 0) {
            printf "  Samples: %d\n", count
            printf "  Avg GPU Utilization: %.1f%%\n", sum_util/count
            printf "  Max GPU Utilization: %.1f%%\n", max_util
            printf "  Avg Memory Used: %.0f MB\n", sum_mem/count
            printf "  Max Memory Used: %.0f MB\n", max_mem
        }
    }'
    echo "  Full log: $GPU_LOG"
fi

# Generate CPU usage summary
if [ -f "$CPU_LOG" ]; then
    echo ""
    echo "CPU Usage Summary:"
    echo "------------------"
    # Skip header, calculate stats
    tail -n +2 "$CPU_LOG" | awk -F',' '
    BEGIN { 
        count=0; sum_cpu=0; sum_mem=0; max_cpu=0; max_mem=0; sum_load=0; max_load=0
    }
    {
        count++
        sum_cpu += $2
        sum_mem += $3
        sum_load += $5
        if ($2 > max_cpu) max_cpu = $2
        if ($3 > max_mem) max_mem = $3
        if ($5 > max_load) max_load = $5
    }
    END {
        if (count > 0) {
            printf "  Samples: %d\n", count
            printf "  Avg CPU Utilization: %.1f%%\n", sum_cpu/count
            printf "  Max CPU Utilization: %.1f%%\n", max_cpu
            printf "  Avg Memory Used: %.0f MB\n", sum_mem/count
            printf "  Max Memory Used: %.0f MB\n", max_mem
            printf "  Avg Load (1min): %.2f\n", sum_load/count
            printf "  Max Load (1min): %.2f\n", max_load
        }
    }'
    echo "  Full log: $CPU_LOG"
fi

exit $PYTHON_EXIT_CODE
