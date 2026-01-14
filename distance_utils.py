# This file contains utils for calculating distances in the end of the PathwayAtlas analysis
# 1. Wasserstein between 2 histograms
# 2. Delta of means between 2 histograms
# 3. PSSM-normalized histogram of background pathway


from definitions import *

import pandas as pd
import numpy as np
import os
import pickle


def get_wasserstein_distance(hist_bg: np.ndarray, hist_cancer: np.ndarray,
                             bin_edges: np.ndarray) -> float:
    """
    Calculates the 1-Wasserstein distance (Earth Mover's Distance)
    between two histograms. Always returns a value >= 0.
    """
    bin_widths = np.diff(bin_edges)

    # Since hist sums to 1, the CDF is just the cumsum of the heights
    # We normalize here just in case the input wasn't perfectly 1.0
    cdf_bg = np.cumsum(hist_bg) / np.sum(hist_bg)
    cdf_cancer = np.cumsum(hist_cancer) / np.sum(hist_cancer)

    # W1 = Integral of |F1 - F2| dx
    return np.sum(np.abs(cdf_bg - cdf_cancer) * bin_widths)


def get_delta_means(hist_bg: np.ndarray, hist_cancer: np.ndarray,
                    bin_edges: np.ndarray) -> float:
    """
    Calculates the difference between the expected values (means)
    of the cancer and background distributions.

    Returns:
        float: E[X_cancer] - E[X_bg]
               (+) Cancer is more pathogenic
               (-) Cancer is more benign
    """
    bin_widths = np.diff(bin_edges)
    bin_centers = bin_edges[:-1] + bin_widths / 2

    # Since heights sum to 1, we don't multiply by bin_widths here
    # We just ensure normalization for safety
    mean_bg = np.sum(bin_centers * hist_bg) / np.sum(hist_bg)
    mean_cancer = np.sum(bin_centers * hist_cancer) / np.sum(hist_cancer)

    return mean_cancer - mean_bg


def get_bg_histogram_after_pssm(pathway_id: str, pssm_matrix=MICHAL_HN1_PSSM, bins=NUMBER_OF_BINS):
    """
    Takes a pathway, finds the scores file (in KEGG_PATHWAY_SCORES_P), takes the "Ref", "Alt" and "pathogenic_prob",
    and turn each of the 12 mutation type to its own histogram, use the PSSM to correct the weight of the mutation,
    and sum them all together (and normalize the histogram heights to sum up to 1).
    :param bins:Number of bins between 0 and 1 TODO check how many to use
    :param pathway_id: KEGG pathway ID, e.g. hsa_00130 or M00873
    :param pssm_matrix: a pssm matrix in dictionary format (see definitions.py)
    :return: np.ndarray
    """
    #  read the pathway file
    #  get the columns
    #  split to 12 groups by the mutation type
    #  for each group, make into histogram
    #    and use pssm to change the mass
    #  sum together
    #  normalize and return
    p_id = pathway_id.replace(":", "_")
    file_path = os.path.join(KEGG_PATHWAY_SCORES_P, f"{p_id}.csv")

    if not os.path.exists(file_path):
        print(f"[Warning] Background file not found: {file_path}")
        return np.zeros(bins)

    # Load only necessary columns to save memory
    df = pd.read_csv(file_path, usecols=["Ref", "Alt", "pathogenic_prob"])

    # Pre-process mutation types
    # Ensure uppercase and format as 'A>C'
    df['mut_type'] = df['Ref'].str.upper() + '>' + df['Alt'].str.upper()

    # Define common bin edges for all histograms
    bin_edges = np.linspace(0, 1, bins + 1)
    final_histogram = np.zeros(bins)

    # Group by mutation type and aggregate
    # This is more efficient than filtering the dataframe 12 times
    grouped = df.groupby('mut_type')['pathogenic_prob']

    total_pssm_weight_found = 0.0

    for mut_type, pssm_weight in pssm_matrix.items():
        if mut_type in grouped.groups:
            # Get all scores for this specific mutation type (e.g., all 'C>T' scores)
            scores = grouped.get_group(mut_type).dropna().values

            if len(scores) == 0:
                continue

            # We take raw counts.
            # This implicitly weights by len(scores), which is the number of
            # available sites (opportunity) for this mutation type in the pathway.
            counts, _ = np.histogram(scores, bins=bin_edges)

            # Apply PSSM weight to the raw counts
            final_histogram += (counts.astype(float) * pssm_weight)

    # Final Normalization
    # If the pathway was missing some mutation types present in the PSSM,
    # we re-normalize so the total probability mass equals 1.
    if final_histogram.sum() > 0:
        final_histogram = final_histogram / final_histogram.sum()

    return final_histogram


def get_cancer_histogram(pathway_id: str, cancer_file: str, bins=NUMBER_OF_BINS):
    """
    Creates a histogram of pathogenic_prob for mutations observed in a specific
    cancer cohort that fall within the genes of a specific pathway.

    :param pathway_id: KEGG pathway ID (e.g., 'hsa00010')
    :param cancer_file: Name of the cancer CSV file (e.g., 'brca_tcga.csv')
    :param bins: Number of bins between 0 and 1 TODO check how many bins to use
    :return: np.ndarray (normalized histogram)
    """
    #  get pathway's genes kegg ids
    #  load columns KeggId and pathogenic_prob from cancer file
    #  take only rows where KeggId is in the pathway's genes
    #  make into histogram
    # 1. Get the list of genes belonging to this pathway
    if not os.path.exists(KEGG_PATHWAY_METADATA_P):
        print(f"[Error] Metadata file not found: {KEGG_PATHWAY_METADATA_P}")
        return np.zeros(bins)

    with open(KEGG_PATHWAY_METADATA_P, "rb") as f:
        pathway_metadata = pickle.load(f)

    if pathway_id not in pathway_metadata:
        print(f"[Warning] Pathway {pathway_id} not found in metadata.")
        return np.zeros(bins)

    # We use a set for O(1) lookup speed
    pathway_genes = set(pathway_metadata[pathway_id]['genes_ids'])

    # 2. Load the cancer mutation file
    cancer_path = os.path.join(CBIO_CANCER_MUTATIONS, cancer_file)
    if not os.path.exists(cancer_path):
        print(f"[Warning] Cancer file not found: {cancer_path}")
        return np.zeros(bins)

    # Load only necessary columns
    # Note: KeggId in cancer files might contain multiple IDs separated by commas
    df = pd.read_csv(cancer_path, usecols=["KeggId", "pathogenic_prob"])

    # 3. Filter rows: only keep mutations in genes belonging to this pathway
    def is_gene_in_pathway(kegg_id_cell):
        if pd.isna(kegg_id_cell):
            return False
        id = str(kegg_id_cell)
        return id.strip() in pathway_genes

    mask = df['KeggId'].apply(is_gene_in_pathway)
    pathway_mutations = df[mask].copy()

    # 4. Create the histogram
    # Extract scores and remove NaNs
    scores = pathway_mutations['pathogenic_prob'].dropna().values

    if len(scores) == 0:
        return np.zeros(bins)

    bin_edges = np.linspace(0, 1, bins + 1)
    counts, _ = np.histogram(scores, bins=bin_edges)

    # 5. Normalize: The sum of the histogram heights should be 1
    total_observed = counts.sum()
    if total_observed > 0:
        return counts.astype(float) / total_observed

    return np.zeros(bins)