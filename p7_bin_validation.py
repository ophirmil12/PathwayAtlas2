import os
import glob
import pandas as pd
import numpy as np
from tqdm import tqdm
from multiprocessing import Pool, cpu_count
from functools import partial

from definitions import *
from distance_utils import *
from kegg_api import KeggNetwork
from p7_bootstrap_pathways import bootstrap, calculate_p_value


def process_single_pathway(pathway_file, cancer_scores_file, cancer_filename, num_bins, n_bootstrap_iters):
    """
    Worker function to process a single pathway.
    Returns a dictionary of results or None if skipped.
    """
    try:
        pathway_name = os.path.splitext(os.path.basename(pathway_file))[0]

        # Load Background Data
        pathway_scores_df = read_scores_file(pathway_name, bins=num_bins)
        if pathway_scores_df is None or (isinstance(pathway_scores_df, np.ndarray) and not pathway_scores_df.any()):
            return None

        # Calculate Histograms
        bg_hist, bin_edges = get_bg_histogram_after_pssm(pathway_scores_df, bins=num_bins)
        cancer_hist = get_cancer_histogram(pathway_name, cancer_filename, bins=num_bins)

        if np.sum(cancer_hist) == 0 or np.sum(bg_hist) == 0:
            return None

        # Setup Bootstrapping Metadata
        kegg_net = KeggNetwork(pathway_name, 'pathway' if pathway_name.startswith('hsa') else 'module')
        num_samples_dict = dict(kegg_net.get_genes_num_cancer_samples(cancer_scores_file))

        if sum(num_samples_dict.values()) == 0:
            return None

        # Calculate Observed Metrics
        w_dist = get_wasserstein_distance(bg_hist, cancer_hist, bin_edges)
        d_means = get_delta_means(bg_hist, cancer_hist, bin_edges)

        # Run Bootstrap
        # Note: bootstrap() contains print statements; in multiprocessing,
        #  these might overlap. You might want to suppress them in bootstrap()
        #  for a cleaner console.
        boot_distances = bootstrap(
            pathway_scores_df,
            bg_hist,
            bin_edges,
            num_samples_dict,
            n_iters=n_bootstrap_iters,
            bins=num_bins
        )

        if len(boot_distances) > 0:
            p_val = calculate_p_value(w_dist, boot_distances)
            return {
                'pathway': pathway_name,
                'bins': num_bins,
                'wasserstein_distance': w_dist,
                'delta_means': d_means,
                'p_value': p_val,
                'total_mutations': sum(num_samples_dict.values()),
                'num_bootstrap_iters': len(boot_distances)
            }
    except Exception as e:
        # Catching exceptions to prevent one pathway from crashing the whole pool
        print(f"Error processing pathway {pathway_file}: {e}")
        return None

    return None


def run_bin_validation():
    # 1. Configuration
    target_cancer = "pan_cancer"
    bin_sizes = [10, 25, 50, 75, 100]
    n_bootstrap_iters = 1000

    # Use most available CPUs, leaving 1 or 2 free for system stability
    num_processes = max(1, cpu_count() - 2)
    print(f"Number of processes: {num_processes} CPU count")

    # Create output directory
    output_dir = pjoin(RESULTS_P, "p7_bin_validation_analysis")
    os.makedirs(output_dir, exist_ok=True)

    # Locate the pan_cancer mutation file
    cancer_files = glob.glob(pjoin(CBIO_CANCER_MUTATIONS_P, f"{target_cancer}*.csv"))
    if not cancer_files:
        print(f"Error: Could not find mutation file for {target_cancer} in {CBIO_CANCER_MUTATIONS_P}")
        return
    cancer_scores_file = cancer_files[0]
    cancer_filename = os.path.basename(cancer_scores_file)

    # Get list of all pathway files
    all_pathway_files = sorted(glob.glob(pjoin(KEGG_PATHWAY_SCORES_P, "*.csv")))

    print(f"Starting bin validation for {target_cancer} using {num_processes} cores...")
    print(f"Pathways found: {len(all_pathway_files)}")

    for num_bins in bin_sizes:
        print(f"\n------ Processing Bin Size: {num_bins} ------")

        # Prepare the function with fixed arguments (partial)
        worker_func = partial(
            process_single_pathway,
            cancer_scores_file=cancer_scores_file,
            cancer_filename=cancer_filename,
            num_bins=num_bins,
            n_bootstrap_iters=n_bootstrap_iters
        )

        results = []

        # Initialize Pool
        with Pool(processes=num_processes) as pool:
            # imap_unordered is often faster and better for memory than map
            # We wrap it in tqdm to visualize progress
            for result in tqdm(pool.imap_unordered(worker_func, all_pathway_files),
                               total=len(all_pathway_files),
                               desc=f"Bins {num_bins}"):
                if result is not None:
                    results.append(result)

        # Save results for this specific bin size
        if results:
            df_results = pd.DataFrame(results)
            output_path = pjoin(output_dir, f"{target_cancer}_bins_{num_bins}.csv")
            df_results.to_csv(output_path, index=False)
            print(f"Saved {len(results)} results for bin {num_bins} to {output_path}")


if __name__ == "__main__":
    run_bin_validation()