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


def bootstrap_dw(pathway_scores_df: pd.DataFrame, ref_hist, num_samples_dict: dict, n_iters = 2000) -> np.array:
    """
    Bootstraps the Wasserstein distances for a given pathway based on the number of samples per protein.
    Args:
        pathway_scores_df (pd.DataFrame): DataFrame containing the pathway scores.
        ref_hist (np.array): Reference histogram of background scores.
        num_samples_dict (dict): Dictionary with the number of samples per protein in the cancer results.
        n_iters (int): Number of bootstrap iterations. Default is 2000.
    Returns:
        np.array: Array of bootstrapped Wasserstein distances.
    """
    distances = np.array

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


def bootstrap_pathway_for_cancer(pathway_scores_file: str, cancer_scoress_file: str) -> None:
    """
    Creates histogram of the background scores for the given pathway,
    calculates the distance between the background and the cancer scores,
    and bootstraps the distance by sampling from the background scores.
    Creates new file in directory RESULTS_DISTANCES_P with the p-values and the distance added.
    Args:
        pathway_scores_file (str): Path to the CSV file containing pathway background scores.
        cancer_scoress_file (str): Path to the CSV file containing cancer scores.
    """
    pathway_name = os.path.splitext(os.path.basename(pathway_scores_file))[0]
    cancer_name = os.path.splitext(os.path.basename(cancer_scoress_file))[0]

    print(f"----- Bootstrapping pathway {pathway_name} for cancer type {cancer_name}... -----")

    pathway_scores_df = pd.read_csv(pathway_scores_file)
    cancer_scores_df = pd.read_csv(cancer_scoress_file)

    # Create weighted histogram for background and cancer scores
    bg_hist, bin_edges = get_bg_histogram_after_pssm(pathway_name)
    cancer_hist, _ =
    # Calculate the distance between cancer scores and background scores
    observed_distance =

    # Bootstrap the distances
    kegg_net = KeggNetwork(pathway_name, 'pathway' if pathway_name.startswith('hsa') else 'module')
    num_samples_dict = dict(kegg_net.get_genes_num_cancer_samples(cancer_scoress_file))
    boot_distances = bootstrap_dw(pathway_scores_df, num_samples_dict)

    p_value = calculate_p_value(observed_distance, boot_distances)

    cancer_distances_df = pd.read_csv(cancer_scoress_file)  # Reload
    cancer_distances_df.loc[cancer_distances_df['pathway'] == pathway_name, 'p_value'] = p_value


if __name__ == '__main__':
    # Get the specific pathway and cancer type based on SLURM_ARRAY_TASK_ID
    args = sys.argv[1:]
    if len(args) < 1:
        print("Usage: python -u p7_bootstrap_pathways.py '$SLURM_ARRAY_TASK_ID'")
        sys.exit(1)

    index = int(args[0])

    all_pathway_files = glob.glob(os.path.join(KEGG_PATHWAY_SCORES_P, f"*.csv"))
    all_cancer_files = glob.glob(os.path.join(CBIO_CANCER_MUTATIONS_P, f"*.csv"))

    if not all_pathway_files or not all_cancer_files:
        print("No pathway  scores or cancer distances files found.")
        sys.exit(1)

    # Calculate pathway and cancer indices
    num_cancer_types = len(all_cancer_files)
    pathway_index = index // num_cancer_types
    cancer_index = index % num_cancer_types

    # Check index bounds
    if pathway_index >= len(all_pathway_files) or cancer_index >= len(all_cancer_files):
        print(f"Index {index} is out of range. There are only {len(all_pathway_files)} pathway files and {len(all_cancer_files)} cancer types.")
        sys.exit(1)

    # Get the specific files
    pathway_scores_file = sorted(all_pathway_files)[pathway_index]
    cancer_scoress_file = sorted(all_cancer_files)[cancer_index]

    print(f"----- Bootstrapping distance in {os.path.basename(cancer_scoress_file)} by pathway {os.path.basename(pathway_scores_file)} -----")
    bootstrap_pathway_for_cancer(pathway_scores_file, cancer_scoress_file)
