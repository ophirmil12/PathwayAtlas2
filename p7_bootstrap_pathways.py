# Bootstrap each pathway for each cancer type

import sys
import os
import glob
from collections import defaultdict
from definitions import *
import pandas as pd
from kegg_api import KeggNetwork
import numpy as np
from typing import Dict

def calculate_p_value(observed_distance: float, bootstrap_distances: np.array) -> float:
    """
    Calculates the p-value based on the observed distance and the distribution
    of distances from permutations.

    Args:
        observed_distance (float): The distance calculated from the actual data.
        bootstrap_distances (list): List of distances from permuted datasets.

    Returns:
        float: The p-value indicating the significance of the observed distance.
    """
    if not bootstrap_distances or len(bootstrap_distances) == 0:
        print("    WARNING: No bootstrap distances provided for p-value calculation, returning p-value 1.0.")
        return 1.0

    # P-value = (Number of random dists >= observed dist) / Total permutations
    # We add 1 to numerator and denominator for pseudo-count smoothing (prevents p=0)
    n_extreme = np.sum(np.abs(bootstrap_distances) >= np.abs(observed_distance))
    p_value = (n_extreme + 1) / (len(bootstrap_distances) + 1)

    return p_value if p_value else 1.0


def get_pathway_scores_background(pathway_df: pd.DataFrame) -> dict:
    """
    Collects mutation scores from the pathway's genes, categorized by
    mutation type (e.g., 'A>C') and then by score type.
    Returns:
        dict: A nested dictionary: {mut_type: {score_type: [scores]}}.
        Example: {"A>C": {"esm_log_probs": [0.1, 0.2], ...}, ...}
    """
    if pathway_df.empty:
        return {}

    background_scores = defaultdict(lambda: defaultdict(list))

    # Check for essential columns in the pre-aggregated CSV
    required_cols = {"Ref", "Alt"}
    if not required_cols.issubset(pathway_df.columns):
        print("Warning: 'Ref' or 'Alt' columns missing in pathway scores file")
        return {}

    score_types = ['esm_log_probs', 'clinvar_reg_dis_ordered_prob', 'clinvar_reg_global_prob']

    # Efficiently iterate over the single, preloaded DataFrame
    for _, row in pathway_df.iterrows():
        mut_type = f"{str(row['Ref']).upper()}>{str(row['Alt']).upper()}"

        for score_type in score_types:
            # Check if score exists and is not NaN
            if score_type in row and pd.notna(row[score_type]):
                background_scores[mut_type][score_type].append(float(row[score_type]))

    # Convert default dicts to regular dicts for a clean return value
    return {mut: dict(scores) for mut, scores in background_scores.items()}


def sample_pathway_by_proteins(bg_scores_df: pd.DataFrame, num_samples_dict: Dict[str, int]) -> pd.DataFrame:
    """
    Samples the background scores DataFrame based on the number of samples per protein in the cancer df.
    Returns:
        pd.DataFrame: A DataFrame containing the sampled rows.
    """
    sampled_rows = []
    groups = dict(tuple(bg_scores_df.groupby('Protein')))

    for protein, num_samples in num_samples_dict.items():
        if protein not in groups:
            continue

        sampled = groups[protein].sample(
            n=num_samples, replace=True, random_state=None
        )
        sampled_rows.append(sampled)

    if sampled_rows:
        return pd.concat(sampled_rows, ignore_index=True)
    else:
        return pd.DataFrame()


def bootstrap_dw(pathway_scores_df: pd.DataFrame, num_samples_dict: dict, n_iters = 2000) -> np.array:
    """
    Bootstraps the Wasserstein distances for a given pathway based on the number of samples per protein.
    Args:
        pathway_scores_df (pd.DataFrame): DataFrame containing the pathway scores.
        num_samples_dict (dict): Dictionary with the number of samples per protein in the cancer results.
        n_iters (int): Number of bootstrap iterations. Default is 2000.
    Returns:
        np.array: Array of bootstrapped Wasserstein distances.
    """
    distances = np.array

    # Prepare reference dictionary of background scores
    ref_scores_dict = get_pathway_scores_background(pathway_scores_df)

    # Create weighted histogram for background scores (PSSM)
    ref_hist, bin_edges = create_joint_distribution(
        ref_scores_dict, MICHAL_HN1_PSSM, use_pssm=True
    )

    # Check for empty distributions
    if np.sum(ref_hist) == 0:
        print(f"    ERROR: background scores histogram sum = 0. Cannot continue bootstrap.")
        return distances

    for b in range (n_iters):
        # Sample scores for each mutation type based on the number of samples
        sampled_df = sample_pathway_by_proteins(pathway_scores_df, num_samples_dict)
        if sampled_df.empty:
            print(f"    ERROR: No proteins were sampled at time {b + 1}/{n_iters}. Skipping iteration.")
            continue

        # Prepare sampled dictionary of scores
        sampled_scores_dict = get_pathway_scores_background(sampled_df)

        # Create weighted histogram for sampled scores (PSSM)
        sampled_hist, _ = create_joint_distribution(
            sampled_scores_dict, MICHAL_HN1_PSSM, use_pssm=True
        )

        # Check for empty distributions
        if np.sum(sampled_hist) == 0:
            print(f"    WARNING: samples histogram sum = 0 at time {b + 1}/{n_iters}. Skipping iteration.")
            continue

        # Calculate Wasserstein distance between sampled and reference histograms
        _, dw_shift = DistributionDistances.directional_wasserstein_from_hist(
            ref_hist, sampled_hist, bin_edges
        )
        distances = np.append(distances, dw_shift)

        if b % 100 == 0 and b > 0:
            print(f"    Iteration {b + 1}/{n_iters}...")

    print(f"    Completed {n_iters} bootstrap iterations.")

    return distances


def bootstrap_pathway_for_cancer(pathway_scores_file: str, cancer_distances_file: str) -> None:
    """
    Bootstraps the given pathway for the specified cancer type by the directional wasserstein distance and updates
    the cancer distances file with the calculated p-value.
    Args:
        pathway_scores_file (str): Path to the CSV file containing pathway background scores.
        cancer_distances_file (str): Path to the CSV file containing cancer distances (dw distance).
    """
    pathway_name = os.path.splitext(os.path.basename(pathway_scores_file))[0]
    cancer_name = os.path.splitext(os.path.basename(cancer_distances_file))[0]

    print(f"----- Bootstrapping pathway {pathway_name} for cancer type {cancer_name}... -----")

    pathway_scores_df = pd.read_csv(pathway_scores_file)
    cancer_distances_df = pd.read_csv(cancer_distances_file)

    kegg_net = KeggNetwork(pathway_name, 'pathway' if pathway_name.startswith('hsa') else 'module')

    num_samples_dict = dict(kegg_net.get_genes_num_cancer_samples(cancer_distances_file))

    dw_distances = bootstrap_dw(pathway_scores_df, num_samples_dict)

    observed_dw = cancer_distances_df.loc[cancer_distances_df['pathway'] == pathway_name, 'dw_distance'].iloc[0]

    p_value_dw = calculate_p_value(observed_dw, dw_distances)

    cancer_distances_df = pd.read_csv(cancer_distances_file)  # Reload
    cancer_distances_df.loc[cancer_distances_df['pathway'] == pathway_name, 'p_value_dw'] = p_value_dw


if __name__ == '__main__':
    # Get the specific pathway and cancer type based on SLURM_ARRAY_TASK_ID
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

    # Calculate pathway and cancer indices
    num_cancer_types = len(all_cancer_distances_files)
    pathway_index = index // num_cancer_types
    cancer_index = index % num_cancer_types

    # Check index bounds
    if pathway_index >= len(all_pathway_files) or cancer_index >= len(all_cancer_distances_files):
        print(f"Index {index} is out of range. There are only {len(all_pathway_files)} pathway files and {len(all_cancer_distances_files)} cancer types.")
        sys.exit(1)

    # Get the specific files
    pathway_scores_file = sorted(all_pathway_files)[pathway_index]
    cancer_distances_file = sorted(all_cancer_distances_files)[cancer_index]

    print(f"----- Bootstrapping the dw distance in {os.path.basename(cancer_distances_file)} by pathway {os.path.basename(pathway_scores_file)} -----")
    bootstrap_pathway_for_cancer(pathway_scores_file, cancer_distances_file)
