#!/bin/bash
#SBATCH --job-name=bootstrap_pathways
#SBATCH --killable
#SBATCH --requeue
#SBATCH --time=48:00:00
#SBATCH --mem=8G
#SBATCH --ntasks=1
###########SBATCH --array=5000-5999
#SBATCH --output=slurm_out/%x_%A_%a.log
###########SBATCH --gres=gg:g4:1

source /cs/labs/dina/ophirmil12/miniforge3/etc/profile.d/conda.sh

# Activate your specific Conda environment
conda activate project_env

umask 003

# Calculate the actual ID: e.g., 5000 + 0, 5000 + 1...
REAL_ID=$((OFFSET + SLURM_ARRAY_TASK_ID))

echo "Slurm Index: $SLURM_ARRAY_TASK_ID"
echo "Calculated Real ID: $REAL_ID"

python -u p7_bootstrap_pathways.py "$REAL_ID"

echo "Job completed."
