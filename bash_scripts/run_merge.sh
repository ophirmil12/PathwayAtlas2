#!/bin/bash
#SBATCH --job-name=cyto
#SBATCH --killable
#SBATCH --requeue
#SBATCH --time=12:00:00
#SBATCH --array=29-29
#SBATCH --mem=16G
#SBATCH --ntasks=1
#SBATCH --output=slurm_out/%x_%A_%a.log
###########SBATCH --gres=gg:g4:1

source /cs/labs/dina/ophirmil12/miniforge3/etc/profile.d/conda.sh

# Activate your specific Conda environment
conda activate project_env
umask 003

echo "Script Starting..."
python -u p16_cytoscape_gene_pathogenicity_visualization.py $SLURM_ARRAY_TASK_ID
# python -u p16_get_interface_in_chimera.py
# python -u scatter_s2_vs_delta_means.py
# python -u p16_gene_level_cox.py
# python -u p13_cox_km_analysis.py $SLURM_ARRAY_TASK_ID
# python -u p13_survival_permutation_testing.py $SLURM_ARRAY_TASK_ID
# python -u p13_create_patient_scores_csvs.py $SLURM_ARRAY_TASK_ID

echo "Job completed."
