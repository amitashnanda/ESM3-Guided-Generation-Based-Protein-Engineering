#!/bin/bash
#  Author: Anna Su
#  Date: 2026-02-20
#  Description: This script reads protein-ligand pairs from a config file and submits
#               a separate SLURM job for each pair.
#  Usage: ./submit_all_jobs.sh [config_file]
#  Default config file: protein_ligand_pairs.conf

set -e

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$(dirname "$SCRIPT_DIR")")"
CONFIG_FILE="${1:-${SCRIPT_DIR}/protein_ligand_pairs.conf}"
JOB_SCRIPT="${SCRIPT_DIR}/guided_generation.sh"
BOLTZ_LOG_DIR="${PROJECT_DIR}/logs/boltz"

# Check if ESM_API_TOKEN is set
if [ -z "$ESM_API_TOKEN" ]; then
    echo "WARNING: ESM_API_TOKEN is not set. Jobs using Forge API models will fail."
    echo "Set it via: export ESM_API_TOKEN='your_token' before running this script."
    echo ""
    read -p "Continue anyway? (y/n) " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        exit 1
    fi
fi

# Check if config file exists
if [ ! -f "$CONFIG_FILE" ]; then
    echo "ERROR: Config file not found: $CONFIG_FILE"
    echo ""
    echo "Please create a config file with the following format:"
    echo "# Comments start with #"
    echo "# Format: PDB_ID,CHAIN_ID,SMILES,MASKING_PCT,NUM_DECODING_STEPS,NUM_SAMPLES,ESM_MODEL,CACHE_DIR,SCORE_TYPE,USE_MSA_SERVER,TIME_LIMIT,CPUS,MEM"
    echo "1RNT,A,NC1=Nc2n(cnc2C(=O)N1)[C@@H]3O[C@H](CO)[C@@H](O)[C@H]3O[P](O)(O)=O,0.4,32,10,esm3-medium-2024-08,./boltz_cache/1rnt,affinity,true,15:00:00,8,16G"
    echo "1DB1,A,CCO,0.3,64,20,local,,affinity,false,24:00:00,12,32G"
    exit 1
fi

# Check if job script exists
if [ ! -f "$JOB_SCRIPT" ]; then
    echo "ERROR: Job script not found: $JOB_SCRIPT"
    exit 1
fi

echo "=============================================="
echo "Submitting Multiple Guided Generation Jobs"
echo "=============================================="
echo "Config File: $CONFIG_FILE"
echo "Job Script: $JOB_SCRIPT"
echo ""

# Counter for jobs
JOB_COUNT=0
SUBMITTED_JOBS=()

# Read config file and submit jobs
while IFS= read -r line || [ -n "$line" ]; do
    # Skip empty lines and comments
    [[ -z "$line" || "$line" =~ ^[[:space:]]*# ]] && continue
    
    # Parse the line (CSV format)
    IFS=',' read -r PDB_ID CHAIN_ID SMILES MASKING_PCT NUM_DECODING_STEPS NUM_SAMPLES ESM_MODEL CACHE_DIR SCORE_TYPE USE_MSA_SERVER TIME_LIMIT CPUS MEM <<< "$line"
    
    # Trim whitespace and set defaults
    PDB_ID=$(echo "$PDB_ID" | xargs)
    CHAIN_ID=$(echo "$CHAIN_ID" | xargs)
    SMILES=$(echo "$SMILES" | xargs)
    MASKING_PCT=$(echo "${MASKING_PCT:-0.4}" | xargs)
    NUM_DECODING_STEPS=$(echo "${NUM_DECODING_STEPS:-32}" | xargs)
    NUM_SAMPLES=$(echo "${NUM_SAMPLES:-10}" | xargs)
    ESM_MODEL=$(echo "${ESM_MODEL:-esm3-medium-2024-08}" | xargs)
    CACHE_DIR=$(echo "$CACHE_DIR" | xargs)
    SCORE_TYPE=$(echo "${SCORE_TYPE:-affinity}" | xargs)
    USE_MSA_SERVER=$(echo "${USE_MSA_SERVER:-false}" | xargs)
    TIME_LIMIT=$(echo "${TIME_LIMIT:-15:00:00}" | xargs)
    CPUS=$(echo "${CPUS:-8}" | xargs)
    MEM=$(echo "${MEM:-16G}" | xargs)
    
    # Validate required fields
    if [ -z "$PDB_ID" ] || [ -z "$CHAIN_ID" ] || [ -z "$SMILES" ]; then
        echo "WARNING: Skipping invalid line: $line"
        continue
    fi
    
    # Generate per-job run directory: run_<PDBID>_<timestamp>
    JOB_TIMESTAMP=$(date +%Y%m%d_%H%M%S)
    JOB_RUN_DIR="run_${PDB_ID}_${JOB_TIMESTAMP}"
    JOB_LOG_DIR="${BOLTZ_LOG_DIR}/${JOB_RUN_DIR}"
    mkdir -p "$JOB_LOG_DIR"
    
    echo "Submitting job for ${PDB_ID}_${CHAIN_ID} (model: ${ESM_MODEL}, time: ${TIME_LIMIT}, cpus: ${CPUS}, mem: ${MEM})..."
    echo "  Log dir: ${JOB_RUN_DIR}"
    
    # Submit job with exported variables (including ESM_API_TOKEN from environment)
    # --output/--error: SLURM writes to the per-job run directory (no stray slurm-*.out)
    # --cpus-per-task and --mem override the #SBATCH defaults in guided_generation.sh
    # RUN_DIR is exported so guided_generation.sh reuses the same directory
    JOB_ID=$(sbatch \
        --job-name="boltz_${PDB_ID}_${CHAIN_ID}" \
        --time="$TIME_LIMIT" \
        --cpus-per-task="$CPUS" \
        --mem="$MEM" \
        --output="${JOB_LOG_DIR}/boltz_${PDB_ID}_${CHAIN_ID}_%j.out" \
        --error="${JOB_LOG_DIR}/boltz_${PDB_ID}_${CHAIN_ID}_%j.err" \
        --export=ALL,ESM_API_TOKEN="$ESM_API_TOKEN",SMILES="$SMILES",PDB_ID="$PDB_ID",CHAIN_ID="$CHAIN_ID",MASKING_PCT="$MASKING_PCT",NUM_DECODING_STEPS="$NUM_DECODING_STEPS",NUM_SAMPLES="$NUM_SAMPLES",ESM_MODEL="$ESM_MODEL",CACHE_DIR="$CACHE_DIR",SCORE_TYPE="$SCORE_TYPE",USE_MSA_SERVER="$USE_MSA_SERVER",TIME_LIMIT="$TIME_LIMIT",RUN_DIR="$JOB_RUN_DIR" \
        "$JOB_SCRIPT" | awk '{print $4}')
    
    # Brief delay to ensure unique timestamps between jobs
    sleep 1
    
    echo "  -> Submitted Job ID: $JOB_ID"
    SUBMITTED_JOBS+=("$JOB_ID:${PDB_ID}_${CHAIN_ID}")
    JOB_COUNT=$((JOB_COUNT + 1))
    
done < "$CONFIG_FILE"

echo ""
echo "=============================================="
echo "Submission Summary"
echo "=============================================="
echo "Total jobs submitted: $JOB_COUNT"
echo ""
echo "Submitted Jobs:"
for job_info in "${SUBMITTED_JOBS[@]}"; do
    IFS=':' read -r job_id job_name <<< "$job_info"
    echo "  Job $job_id: $job_name"
done
echo ""
echo "Monitor jobs with: squeue -u $USER"
echo "Cancel all jobs with: scancel ${SUBMITTED_JOBS[*]%%:*}"
