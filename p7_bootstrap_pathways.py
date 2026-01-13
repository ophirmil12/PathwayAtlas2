# Bootstrap each pathway for each cancer type

import sys
import os
import glob
from ensurepip import bootstrap

from definitions import *
import pandas as pd
from kegg_api import KeggNetwork
import numpy as np

def calculate_p_value(observed_distance, bootstrap_distances):
    """
    Calculates the p-value based on the observed distance and the distribution
    of distances from permutations.

    Args:
        observed_distance (float): The distance calculated from the actual data.
        bootstrap_distances (list): List of distances from permuted datasets.

    Returns:
        float: The p-value indicating the significance of the observed distance.
    """
    if not bootstrap_distances:
        return 1.0

    # Convert to numpy array for efficiency
    dist_array = np.array(bootstrap_distances)

    # P-value = (Number of random dists >= observed dist) / Total permutations
    # We add 1 to numerator and denominator for pseudo-count smoothing (prevents p=0)
    n_extreme = np.sum(np.abs(dist_array) >= np.abs(observed_distance))
    p_value = (n_extreme + 1) / (len(dist_array) + 1)

    return p_value if p_value else 1.0

def bootstrap_dw(pathway_scores_df: pd.DataFrame, num_samples_dict: dict) -> list:
    pass

def bootstrap_pathway_for_cancer(pathway_file: str, cancer_distances_file: str) -> None:
    pathway_name = os.path.splitext(os.path.basename(pathway_file))[0]
    cancer_name = os.path.splitext(os.path.basename(cancer_distances_file))[0]

    print(f"----- Bootstrapping pathway {pathway_name} for cancer type {cancer_name}... -----")

    pathway_scores_df = pd.read_csv(pathway_file)
    cancer_distances_df = pd.read_csv(cancer_distances_file)

    kegg_net = KeggNetwork(pathway_name, 'pathway' if pathway_name.startswith('hsa') else 'module')

    num_samples_dict = dict(kegg_net.get_genes_num_cancer_samples(cancer_distances_file))

    dw_distances = bootstrap_dw(pathway_scores_df, num_samples_dict)

    observed_dw = cancer_distances_df.loc[cancer_distances_df['pathway'] == pathway_name, 'dw_distance'].iloc[0]

    p_value_dw = calculate_p_value(observed_dw, dw_distances)

    cancer_distances_df = pd.read_csv(cancer_distances_file)  # Reload
    cancer_distances_df.loc[cancer_distances_df['pathway'] == pathway_name, 'p_value_dw'] = p_value_dw


if __name__ == '__main__':
    args = sys.argv[1:]
    if len(args) < 1:
        print("Usage: python -u p7_bootstrap_pathways.py '$SLURM_ARRAY_TASK_ID'")
        sys.exit(1)

    index = int(args[0])

    all_pathway_files = glob.glob(os.path.join(KEGG_PATHWAY_SCORES_P, f"*.csv"))
    all_cancer_distances_files = glob.glob(os.path.join(RESULTS_DISTANCES_P, f"*.csv"))

    if not all_pathway_files or not all_cancer_distances_files:
        print("No pathway  scores or cancer distances files found.")
        sys.exit(1)

    num_cancer_types = len(all_cancer_distances_files)
    pathway_index = index // num_cancer_types
    cancer_index = index % num_cancer_types

    if pathway_index >= len(all_pathway_files) or cancer_index >= len(all_cancer_distances_files):
        print(f"Index {index} is out of range. There are only {len(all_pathway_files)} pathway files and {len(all_cancer_distances_files)} cancer types.")
        sys.exit(1)

    pathway_file = sorted(all_pathway_files)[pathway_index]
    cancer_distances_file = sorted(all_cancer_distances_files)[cancer_index]