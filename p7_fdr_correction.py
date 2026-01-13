# Calculate the q_values of cancer-pathway distances after bootstrapping using bh FDR correction
# The q values are calculated per cancer

import os
from os.path import join as pjoin
import pandas as pd
import sys
import glob
from definitions import *
from scipy.stats import false_discovery_control

def perform_fdr_correction(cancer_results_file: str):
    if not os.path.exists(cancer_results_file):
        print(f"    ERROR: {cancer_results_file} not found.")

    cancer_results_df = pd.read_csv(cancer_results_file)
    if not cancer_results_df:
        print(f"    ERROR: could not read {cancer_results_file}.")

    p_values_dw = cancer_results_df['p_value_dw']

    if not p_values_dw:
        print(f"    ERROR: p values not found, please run bootstrap first.")

    try:
        q_values_dw = false_discovery_control(p_values_dw, method='bh')
    except Exception as e:
        print(f"    ERROR during FDR correction: {e}")
        return

    cancer_results_df['q_value_dw'] = q_values_dw
    cancer_results_df.to_csv(cancer_results_file, index=False)
    print(f"    FDR correction completed and saved to {cancer_results_file}.")


if __name__ == '__main__':
    # Get the cancer results file based on SLURM_ARRAY_TASK_ID
    args = sys.argv[1:]
    if len(args) < 1:
        print("Usage: python -u p7_fdr_correction.py '$SLURM_ARRAY_TASK_ID'")
        sys.exit(1)

    cancer_results_files = glob.glob(pjoin(RESULTS_DISTANCES_P, "*.csv"))  # should be 37 files
    if not cancer_results_files:
        print(f"No cancer distance results CSV files found in the specified directory: {RESULTS_DISTANCES_P}")
        sys.exit(1)

    index = int(args[0])
    if index < 0 or index >= len(cancer_results_files):
        print(f"Index {index} is out of range. There are only {len(cancer_results_files)} mutation study files.")
        sys.exit(1)

    cancer_results_file = sorted(cancer_results_files)[index]

    print(f"----- Performing FDR correction for file {cancer_results_file} -----")
    perform_fdr_correction(cancer_results_file)
