#!/bin/bash
#SBATCH --job-name=merge_pathways
#SBATCH --killable
#SBATCH --requeue
#SBATCH --time=48:00:00
#SBATCH --mem=16G
#SBATCH --ntasks=1
#SBATCH --array=0-0
#SBATCH --output=slurm_out/%x_%A_%a.log
###########SBATCH --gres=gg:g4:1

source /cs/labs/dina/ophirmil12/miniforge3/etc/profile.d/conda.sh

# Activate your specific Conda environment
conda activate project_env

umask 003

echo "Script Starting..."

python -u p2_merge_studies_by_cancer.py

echo "Job completed."
