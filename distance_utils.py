# This file contains utils for calculating distances in the end of the PathwayAtlas analysis
# 1. Wasserstein between 2 histograms
# 2. Delta of means between 2 histograms
# 3. PSSM-normalized histogram of background pathway


from definitions import *

import pandas as pd
import numpy as np
import os



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

    file_path = os.path.join(KEGG_PATHWAY_SCORES_P, f"{pathway_id}.csv")

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