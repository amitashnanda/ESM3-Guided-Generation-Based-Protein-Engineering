#!/bin/bash
#SBATCH -N 1                      
#SBATCH -C gpu                    
#SBATCH -q regular                
#SBATCH -J esm3-foldx-parallel    
#SBATCH -t 06:00:00              
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --account=""
#SBATCH --mail-type=end,fail
#SBATCH --mail-user=""

#SBATCH --array=0-1 

#SBATCH -o foldx_run_%j.log

# PDB_FILES=("1rnt.pdb" "1hii.pdb" "1fbm.pdb" "1db1.pdb" "6q30.pdb")
# CHAIN_IDS=("A" "A" "A" "A" "A")

# PDB_FILES=("1hii.pdb" "1fbm.pdb" "1db1.pdb")
PDB_FILES=("1db1.pdb" "6q30.pdb")
CHAIN_IDS=("A" "A")

# MASK_RATES=(0.3 0.3 0.4 0.20 0.20)
MASK_RATES=(0.2 0.2)

export SLURM_CPU_BIND=cores

CURRENT_PDB=${PDB_FILES[$SLURM_ARRAY_TASK_ID]}
CURRENT_CHAIN=${CHAIN_IDS[$SLURM_ARRAY_TASK_ID]}
CURRENT_MASK=${MASK_RATES[$SLURM_ARRAY_TASK_ID]}

CURRENT_STEPS=${DECODING_STEPS[$SLURM_ARRAY_TASK_ID]}
CURRENT_SAMPLES=${SAMPLES_PER_STEP[$SLURM_ARRAY_TASK_ID]}
CURRENT_WORKERS=${NUM_WORKERS[$SLURM_ARRAY_TASK_ID]}

echo "========================================================"
echo "STARTING JOB: ${SLURM_JOB_ID}, ARRAY_TASK: ${SLURM_ARRAY_TASK_ID}"
echo "Protein: ${CURRENT_PDB}, Chain: ${CURRENT_CHAIN}, Mask: ${CURRENT_MASK}"
echo "Steps: ${CURRENT_STEPS}, Samples: ${CURRENT_SAMPLES}, Workers: ${CURRENT_WORKERS}"
echo "========================================================"

module load conda

conda activate /global/common/software/nintern/ananda/proteinenv

cd /pscratch/sd/a/ananda/protein_ai/src/


srun python -m esm_foldx_guidedgeneration.main --pdb_filename "$CURRENT_PDB" --chain_id "$CURRENT_CHAIN" --masking_percentage "$CURRENT_MASK" --num_decoding_steps "$CURRENT_STEPS" --num_samples_per_step "$CURRENT_SAMPLES" --num_workers "$CURRENT_WORKERS"

echo "========================================================"
echo "FINISHED JOB: ${SLURM_JOB_ID}, ARRAY_TASK: ${SLURM_ARRAY_TASK_ID}"
echo "========================================================"


