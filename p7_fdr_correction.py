# Calculate the q_values of cancer-pathway distances after bootstrapping using bh FDR correction
# The q values are calculated per cancer

import os
from os.path import join as pjoin
import pandas as pd
import sys
import glob
from definitions import *
from scipy.stats import false_discovery_control
import numpy as np
import pickle 

def perform_fdr_correction(cancer_results_df: pd.DataFrame, output_path: str, recalc: bool = False):
    
    if 'p_value' not in cancer_results_df.columns or cancer_results_df['p_value'].isnull().all():
        print(f"    ERROR: p values not found, please run bootstrap first.")

    if (recalc == False) and 'q_value' in cancer_results_df.columns and not cancer_results_df['q_value'].isna().all():
        print(f"    Q values already exist in {cancer_results_file}, skipping FDR correction.")
        print(f"    Percent significant (q < 0.05): {(np.sum(cancer_results_df['q_value'] < 0.05) / (len(cancer_results_df) )) * 100:.2f}%")
        return
    
    # check if there are any p values that are null
    if cancer_results_df['p_value'].isnull().any():
        print(f"    WARNING: Found {cancer_results_df['p_value'].isnull().sum()} null p values in {cancer_results_file}, consider rerunning bootstrap.")
    else:
        cancer_results_df = cancer_results_df.dropna(subset=['p_value'])

        p_values = cancer_results_df['p_value'].tolist()

        try:
            q_values = false_discovery_control(p_values, method='bh')
        except Exception as e:
            print(f"    ERROR during FDR correction: {e}")
            return

        cancer_results_df.loc[:, 'q_value'] = q_values
        cancer_results_df.to_csv(output_path, index=False)
        print(f"    FDR correction completed and saved to {output_path}.")
        print(f"    Percent significant (q < 0.05): {(np.sum(cancer_results_df['q_value'] < 0.05) / (len(cancer_results_df) )) * 100:.2f}%")


def perform_fdr_correction_filtered_pathways(cancer_results_file: str, output_path: str):
    if not os.path.exists(cancer_results_file):
        print(f"    ERROR: {cancer_results_file} not found.")

    cancer_results_df = pd.read_csv(cancer_results_file)

    with open(FILTERED_KEGG_PATHWAY_METADATA_FILE, 'rb') as f:
        pathway_metadata = pickle.load(f)

    pathways_to_keep = pathway_metadata.keys()

    filtered_cancer_results_df = cancer_results_df[
        cancer_results_df["pathway"].isin(pathways_to_keep)
    ]

    print(f"Number of pathways before filtering: {len(cancer_results_df)}, after: {len(filtered_cancer_results_df)}")
    

    perform_fdr_correction(filtered_cancer_results_df, output_path, recalc=True)


def find_differences_in_significance(distances_before_filtering_df, distances_after_filtering_df):
    ALPHA = 0.05

    before = distances_before_filtering_df[['pathway', 'q_value']].copy()
    after = distances_after_filtering_df[['pathway', 'q_value']].copy()

    before_sig = set(before.loc[before['q_value'] < ALPHA, 'pathway'])
    after_sig = set(after.loc[after['q_value'] < ALPHA, 'pathway'])
    after_all = set(after['pathway'])

    # Was significant before, no longer significant (but still present)
    lost_significance = before_sig & after_all - after_sig
    # Was significant before, removed entirely after filtering
    removed = before_sig - after_all
    # Was not significant before, now significant
    gained_significance = after_sig - before_sig

    print(f"\n===== Significance Differences (alpha={ALPHA}) =====")
    print(f"  Significant before: {len(before_sig)}, after: {len(after_sig)}")

    if removed:
        print(f"\n  Significant pathways removed by filtering ({len(removed)}):")
        for p in sorted(removed):
            pathway_name = distances_before_filtering_df.loc[
                distances_before_filtering_df['pathway'] == p,
                'pathway_name'
            ]
            q = before.loc[before['pathway'] == p, 'q_value'].values[0]
            print(f"    - {p}: {pathway_name} (q={q:.4f})")
    else:
        print("\n  No significant pathways were removed by filtering.")

    if lost_significance:
        print(f"\n  Pathways that lost significance after filtering ({len(lost_significance)}):")
        for p in sorted(lost_significance):
            pathway_name = distances_before_filtering_df.loc[
                distances_before_filtering_df['pathway'] == p,
                'pathway_name'
            ]
            q_before = before.loc[before['pathway'] == p, 'q_value'].values[0]
            q_after = after.loc[after['pathway'] == p, 'q_value'].values[0]
            print(f"    - {p} : {pathway_name}  (q before={q_before:.4f} -> q after={q_after:.4f})")
    else:
        print("\n  No pathways lost significance after filtering.")

    if gained_significance:
        print(f"\n  Pathways that gained significance after filtering ({len(gained_significance)}):")
        for p in sorted(gained_significance):
            pathway_name = distances_before_filtering_df.loc[
                distances_before_filtering_df['pathway'] == p,
                'pathway_name'
            ]
            q_after = after.loc[after['pathway'] == p, 'q_value'].values[0]
            q_before_val = before.loc[before['pathway'] == p, 'q_value'].values
            q_before_str = f"{q_before_val[0]:.4f}" if len(q_before_val) > 0 else "not present"
            print(f"    - {p}: {pathway_name} (q before={q_before_str} -> q after={q_after:.4f})")
    else:
        print("\n  No pathways gained significance after filtering.")

    print("=" * 50)



if __name__ == '__main__':
    # Get the cancer results file based on SLURM_ARRAY_TASK_ID
    args = sys.argv[1:]
    if len(args) < 1:
        print("Usage: python -u p7_fdr_correction.py '$SLURM_ARRAY_TASK_ID'")
        sys.exit(1)

    cancer_results_files = glob.glob(pjoin(RESULTS_DISTANCES_P, "*.csv"))  # should be 45 files
    if not cancer_results_files:
        print(f"No cancer distance results CSV files found in the specified directory: {RESULTS_DISTANCES_P}")
        sys.exit(1)

    index = int(args[0])
    if index < 0 or index >= len(cancer_results_files):
        print(f"Index {index} is out of range. There are only {len(cancer_results_files)} cancer results files.")
        sys.exit(1)

    cancer_results_file = sorted(cancer_results_files)[index]

    print(f"----- Performing FDR correction for file {cancer_results_file} -----")

    if not os.path.exists(cancer_results_file):
        print(f"    ERROR: {cancer_results_file} not found.")

    cancer_results_df = pd.read_csv(cancer_results_file)
    perform_fdr_correction(cancer_results_df, cancer_results_file)

    cancer_type = os.path.basename(cancer_results_file).replace('.csv', '')  # get cancer type from file name
    output_path = pjoin(FILTERED_RESULTS_DISTANCES_P, f"{cancer_type}.csv")
    perform_fdr_correction_filtered_pathways(cancer_results_file, output_path)

    filtered_cancer_results_df = pd.read_csv(output_path)

    find_differences_in_significance(cancer_results_df, filtered_cancer_results_df)
