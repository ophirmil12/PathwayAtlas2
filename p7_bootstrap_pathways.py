# Bootstrap each pathway for each cancer type

import sys
import glob
import pandas as pd
from kegg_api import KeggNetwork
from typing import Dict, List
from distance_utils import *

def calculate_p_value(observed_distance: float, bootstrap_distances: List) -> float:
    """
    Calculates the p-value based on the observed distance and the distribution
    of distances from permutations.

    Args:
        observed_distance (float): The distance calculated from the actual data.
        bootstrap_distances (list): List of distances from permuted datasets.

    Returns:
        float: The p-value indicating the significance of the observed distance.
    """

    distances_array = np.array(bootstrap_distances)

    # P-value = (Number of random dists >= observed dist) / Total permutations
    # We add 1 to numerator and denominator for pseudo-count smoothing (prevents p=0)
    n_extreme = np.sum(distances_array >= observed_distance)
    p_value = (n_extreme + 1) / (len(bootstrap_distances) + 1)

    if p_value == 1.0:
        print("    WARNING: P-value calculated as 1.0, indicating no extreme distances found.")

    return p_value


def get_sampled_hist(bg_scores_df: pd.DataFrame, num_samples_dict: Dict[str, int]) -> np.ndarray:
    """
    Samples from the background histogram correctly using PSSM weights.
    """
    sampled_protein_dfs = []

    # Create the group dictionary once
    # Ensure we use 'KeggId' as established
    groups = dict(tuple(bg_scores_df.groupby('KeggId')))

    for kegg_gene, num_samples in num_samples_dict.items():
        if num_samples == 0:
            continue

        if kegg_gene not in groups:
            continue

        # 1. Get the histogram (probabilities) and edges
        protein_hist, bin_edges = get_bg_histogram_after_pssm(groups[kegg_gene])

        # 2. Calculate the values to sample FROM (The X-axis / Bin Centers)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

        # If sum is 0 (empty histogram), skip
        total_prob = protein_hist.sum()
        if total_prob == 0:
            continue
        p_norm = protein_hist / total_prob

        # 3. Sample from the CENTERS, using the HISTOGRAM as the PROBABILITY
        sampled_scores = np.random.choice(
            bin_centers,
            size=num_samples,
            replace=True,
            p=p_norm
        )

        # 4. Create a DataFrame for the sampled scores
        # We perform this step to match your existing pipeline flow
        sampled_protein_df = pd.DataFrame({
            'pathogenic_prob': sampled_scores,
            'KeggId': [kegg_gene] * num_samples
        })
        sampled_protein_dfs.append(sampled_protein_df)

    if sampled_protein_dfs:
        # 5. Concatenate all sampled protein DataFrames
        sampled_pathway_df = pd.concat(sampled_protein_dfs, ignore_index=True)
        # Return the histogram of this new sampled population
        return create_histogram(sampled_pathway_df)
    else:
        # print("    ERROR: No proteins were sampled.")
        return np.array([])


def bootstrap(pathway_scores_df: pd.DataFrame, ref_hist, bin_edges, num_samples_dict: dict, n_iters = 1000) -> List:
    """
    Bootstraps the Wasserstein distances for a given pathway based on the number of samples per protein.
    Args:
        pathway_scores_df (pd.DataFrame): DataFrame containing the pathway scores.
        ref_hist (np.array): Reference histogram of background scores.
        bin_edges (np.array): Edges of the histogram bins.
        num_samples_dict (dict): Dictionary with the number of samples per protein in the cancer results.
        n_iters (int): Number of bootstrap iterations. Default is 2000.
    Returns:
        np.array: Array of bootstrapped Wasserstein distances.
    """
    distances = []

    # Check for empty distributions
    if np.sum(ref_hist) == 0:
        print(f"    ERROR: background scores histogram sum = 0. Cannot continue bootstrap.")
        return distances

    print("    Starting bootstrap iterations...")

    for b in range (n_iters):
        # Sample scores for each mutation type based on the number of samples
        sampled_hist = get_sampled_hist(pathway_scores_df, num_samples_dict)

        # Check for empty distributions
        if np.sum(sampled_hist) == 0:
            print(f"    WARNING: samples histogram sum = 0 at time {b}/{n_iters}. Skipping iteration.")
            continue

        # Calculate Wasserstein distance between sampled and reference histograms
        w_distance = get_wasserstein_distance(ref_hist, sampled_hist, bin_edges)
        distances.append(w_distance)

        if b % 50 == 0 and b > 0:
            print(f"    Iteration {b}/{n_iters}...")

    print(f"    Completed {n_iters} bootstrap iterations.")

    return distances


def bootstrap_pathway_for_cancer(pathway_scores_file: str, cancer_scores_file: str) -> None:
    """
    Creates histogram of the background scores for the given pathway,
    calculates the distance between the background and the cancer scores,
    and bootstraps the distance by sampling from the background scores.
    Creates new file in directory RESULTS_DISTANCES_P with the p-values and the distance added.
    Args:
        pathway_scores_file (str): Path to the CSV file containing pathway background scores.
        cancer_scores_file (str): Path to the CSV file containing cancer scores.
    """
    pathway_name = os.path.splitext(os.path.basename(pathway_scores_file))[0]

    # Create weighted histogram for background and cancer scores
    pathway_scores_df = read_scores_file(pathway_name)
    bg_hist, bin_edges = get_bg_histogram_after_pssm(pathway_scores_df)
    cancer_hist = get_cancer_histogram(pathway_name, cancer_scores_file)

    if np.sum(cancer_hist) == 0 or np.sum(bg_hist) == 0:
        print(f"    WARNING: One of the histograms is empty for pathway {pathway_name}. Skipping bootstrap.")
        return

    # Calculate the distance between cancer scores and background scores
    w_distance = get_wasserstein_distance(bg_hist, cancer_hist, bin_edges)
    delta_means = get_delta_means(bg_hist, cancer_hist, bin_edges)

    # Bootstrap the distances
    kegg_net = KeggNetwork(pathway_name, 'pathway' if pathway_name.startswith('hsa') else 'module')
    num_samples_dict = dict(kegg_net.get_genes_num_cancer_samples(cancer_scores_file))
    if sum(num_samples_dict.values()) == 0:
        print(f"    WARNING: No samples found for pathway {pathway_name} in cancer file {cancer_scores_file}. Skipping bootstrap.")
        return

    boot_distances = bootstrap(pathway_scores_df, bg_hist, bin_edges, num_samples_dict)

    if len(boot_distances) == 0:
        print("    WARNING: No bootstrap distances provided for p-value calculation.")
        return

    p_value = calculate_p_value(w_distance, boot_distances)
    print(f"    Pathway: {pathway_name}, P-value: {p_value}")

    cancer_distances_file = pjoin(RESULTS_DISTANCES_P, os.path.basename(cancer_scores_file))

    print(f"    Updating distances file: {cancer_distances_file}...")
    cancer_distances_df = pd.read_csv(cancer_distances_file)
    cancer_distances_df.loc[cancer_distances_df['pathway'] == pathway_name, 'wasserstein_distance'] = w_distance
    cancer_distances_df.loc[cancer_distances_df['pathway'] == pathway_name, 'delta_means'] = delta_means
    cancer_distances_df.loc[cancer_distances_df['pathway'] == pathway_name, 'p_value'] = p_value
    cancer_distances_df.to_csv(cancer_distances_file, index=False)


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
    cancer_scores_file = sorted(all_cancer_files)[cancer_index]

    cancer_distances_file = pjoin(RESULTS_DISTANCES_P, os.path.basename(cancer_scores_file))

    if not os.path.exists(cancer_distances_file):
        print(f"Creating new distances file for {os.path.basename(cancer_scores_file)}...")
        cancer_distances_df = pd.DataFrame({'pathway': [os.path.splitext(os.path.basename(f))[0] for f in sorted(all_pathway_files)],
            'wasserstein_distance': [np.nan] * len(all_pathway_files),
            'delta_means': [np.nan] * len(all_pathway_files),
            'p_value': [np.nan] * len(all_pathway_files)
        })
        cancer_distances_df.to_csv(cancer_distances_file, index=False)

    else:
        print(f"Bootstrapping already completed for {os.path.basename(cancer_scores_file)} and pathway {os.path.basename(pathway_scores_file)}")
        sys.exit(0)

    print(f"----- Bootstrapping distance in {os.path.basename(cancer_scores_file)} by pathway {os.path.basename(pathway_scores_file)} -----")
    bootstrap_pathway_for_cancer(pathway_scores_file, cancer_scores_file)
