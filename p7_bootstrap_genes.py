import sys
import os
import glob
import pandas as pd
import numpy as np
from multiprocessing import Pool, cpu_count
from functools import partial
from tqdm import tqdm

from definitions import *
from distance_utils import *
from p7_bootstrap_pathways import bootstrap, calculate_p_value


def process_single_gene(kegg_id, gene_cancer_scores, n_iters, bins):
    """
    Worker function.
    kegg_id: e.g., 'hsa:287'
    gene_cancer_scores: List of observed pathogenic_probs (floats)
    """
    safe_id = kegg_id.replace(":", "_")
    bg_file = pjoin(KEGG_GENE_SCORES_P, f"{safe_id}.csv")

    if not os.path.exists(bg_file):
        return None

    try:
        # Load Background
        bg_df = pd.read_csv(bg_file, usecols=["KeggId", "Ref", "Alt", "pathogenic_prob"])

        # --- FIX: Strict Type Casting ---
        # Convert to numeric, turn errors into NaN, then drop NaNs
        bg_df['pathogenic_prob'] = pd.to_numeric(bg_df['pathogenic_prob'], errors='coerce')
        bg_df = bg_df.dropna(subset=['pathogenic_prob'])

        # 1. Generate Histograms
        bg_hist, bin_edges = get_bg_histogram_after_pssm(bg_df, bins=bins)

        # Observed data is already sanitized in the main loop, but we wrap it in a DF
        obs_df = pd.DataFrame({'pathogenic_prob': gene_cancer_scores})
        obs_hist = create_histogram(obs_df, bins=bins)

        if np.sum(obs_hist) == 0 or np.sum(bg_hist) == 0:
            return None

        # 2. Calculate Distances
        w_dist = get_wasserstein_distance(bg_hist, obs_hist, bin_edges)
        d_means = get_delta_means(bg_hist, obs_hist, bin_edges)

        # 3. Bootstrap
        num_obs = len(gene_cancer_scores)
        num_samples_dict = {kegg_id: num_obs}

        boot_distances = bootstrap(
            bg_df,
            bg_hist,
            bin_edges,
            num_samples_dict,
            n_iters=n_iters,
            bins=bins
        )

        if not boot_distances:
            return None

        p_val = calculate_p_value(w_dist, boot_distances)

        return {
            'gene_id': kegg_id,
            'wasserstein_distance': w_dist,
            'delta_means': d_means,
            'p_value': p_val,
            'num_mutations': num_obs
        }
    except Exception as e:
        sys.stderr.write(f"Error in {kegg_id}: {str(e)}\n")
        return None


def run_gene_bootstrapping(cancer_type):
    n_iters = 1000
    bins = 100

    output_dir = RESULTS_GENE_LEVEL_P
    os.makedirs(output_dir, exist_ok=True)
    output_file = pjoin(output_dir, f"{cancer_type}_gene_distances.csv")

    search_path = pjoin(CBIO_CANCER_MUTATIONS_P, f"{cancer_type}*.csv")
    cancer_files = glob.glob(search_path)
    if not cancer_files:
        print(f"CRITICAL ERROR: No cancer files found for {cancer_type}")
        return

    cancer_path = cancer_files[0]
    print(f"--- Processing {cancer_type} ---")

    # --- FIX: Sanitize Cancer Data ---
    df = pd.read_csv(cancer_path, usecols=["KeggId", "pathogenic_prob"])
    # Convert to numeric, force 'str' probabilities to NaN, then drop
    df['pathogenic_prob'] = pd.to_numeric(df['pathogenic_prob'], errors='coerce')
    df = df.dropna(subset=['KeggId', 'pathogenic_prob'])

    gene_groups = df.groupby('KeggId')['pathogenic_prob'].apply(list).to_dict()
    unique_genes = list(gene_groups.keys())

    # Pre-flight Check
    valid_tasks = []
    for gid in unique_genes:
        safe_id = gid.replace(":", "_")
        if os.path.exists(pjoin(KEGG_GENE_SCORES_P, f"{safe_id}.csv")):
            valid_tasks.append((gid, gene_groups[gid], n_iters, bins))

    print(f"Ready to process: {len(valid_tasks)} genes.")

    if not valid_tasks:
        return

    # Multiprocessing
    num_processes = max(1, cpu_count() - 1)
    with Pool(processes=num_processes) as pool:
        results = []
        for result in tqdm(pool.starmap(process_single_gene, valid_tasks), total=len(valid_tasks)):
            if result is not None:
                results.append(result)

    if results:
        pd.DataFrame(results).to_csv(output_file, index=False)
        print(f"SUCCESS: Created {output_file}")


if __name__ == "__main__":
    run_gene_bootstrapping(sys.argv[1])