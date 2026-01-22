# To provide visual proof of statistical findings, we will create overlay plots where the PSSM-weighted background
# (the expectation) is compared directly against the observed cancer distribution (the reality).

from plot_boot import *
boot_plot_folder()


import os
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pickle
import textwrap
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm

# Set non-interactive backend for server-side plotting (Crucial for multiprocessing)
import matplotlib

matplotlib.use('Agg')

from definitions import (
    RESULTS_DISTANCES_P,
    KEGG_PATHWAY_METADATA_FILE,
    PLOTS_P,
    COLOR_MAP,
    NUMBER_OF_BINS,
    MICHAL_HN1_PSSM, KEGG_PATHWAY_SCORES_P,
)
from distance_utils import get_bg_histogram_after_pssm, get_cancer_histogram


def plot_single_pathway(args):
    (pw_id, pw_name, q_val, w1_dist, d_mean, n_mut, cancer_name, cancer_plot_dir) = args
    bin_edges = np.linspace(0, 1, NUMBER_OF_BINS + 1)

    try:
        # 1. Generate Histograms
        p_id = pw_id.replace(":", "_")
        file_path = os.path.join(KEGG_PATHWAY_SCORES_P, f"{p_id}.csv")

        # FIX: Proper early return if background file is missing
        if not os.path.exists(file_path):
            return False

        df = pd.read_csv(file_path, usecols=["Ref", "Alt", "pathogenic_prob"])

        # FIX: Unpack the tuple returned by distance_utils
        hist_bg, _ = get_bg_histogram_after_pssm(df, pssm_matrix=MICHAL_HN1_PSSM)

        hist_cancer = get_cancer_histogram(pw_id, f"{cancer_name}.csv")

        if np.sum(hist_cancer) == 0:
            return False

        # 2. Plotting
        plt.figure(figsize=(10, 6.5))

        # Expected BG (PDF Steps)
        plt.stairs(hist_bg, bin_edges,
                   color=COLOR_MAP["dark-blue"], fill=True, alpha=0.25,
                   label='Expected (PSSM-Weighted BG)')
        plt.stairs(hist_bg, bin_edges,
                   color=COLOR_MAP["dark-blue"], linewidth=1, alpha=0.6)

        # Observed Cancer (PDF Steps)
        # Use .get() or index safety for COLOR_MAP
        path_color = COLOR_MAP.get("pathogenic")
        benign_color = COLOR_MAP.get("benign")
        color_choice = path_color if d_mean > 0 else benign_color

        plt.stairs(hist_cancer, bin_edges,
                   color=color_choice, fill=True, alpha=0.4,
                   label=f'Observed ({cancer_name.upper()})')
        plt.stairs(hist_cancer, bin_edges,
                   color=color_choice, linewidth=2.5)

        # 3. Stats Box & Title (Same as your logic)
        sig_status = "SIGNIFICANT" if q_val < 0.05 else "Not Significant"
        stats_text = (f"Q-value: {q_val:.2e}\nStatus:  {sig_status}\n"
                      f"W1 Dist: {w1_dist:.4f}\nΔ Mean:  {d_mean:+.4f}\n"
                      f"Mutation #: {n_mut}")

        plt.annotate(stats_text, xy=(0.02, 0.96), xycoords='axes fraction',
                     va='top', ha='left', fontsize=10, family='monospace',
                     bbox=dict(boxstyle="round,pad=0.5", fc="white", ec="black", alpha=0.8))

        wrapped_name = "\n".join(textwrap.wrap(pw_name.split(' - Homo sapiens')[0], width=50))
        plt.title(f"{wrapped_name}\n({pw_id})", fontsize=15, fontweight='bold', pad=15)

        plt.xlabel("Pathogenicity Probability", fontsize=12)
        plt.ylabel("Probability Mass", fontsize=12)
        plt.xlim(0, 1)

        # FIX: Now that hist_bg is an array, .max() works!
        y_max = max(hist_bg.max(), hist_cancer.max()) * 1.3
        plt.ylim(0, y_max)
        plt.legend(loc='upper right', frameon=True)
        plt.grid(axis='y', linestyle=':', alpha=0.4)

        # 4. Save
        clean_filename = pw_id.replace(":", "_")
        save_path = os.path.join(cancer_plot_dir, f"{clean_filename}_overlay.png")
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close()
        return True
    except Exception as e:
        print(f"Error plotting {pw_id}: {e}")
        return False


def plot_all_overlays_parallel():
    # 1. Load Metadata
    with open(KEGG_PATHWAY_METADATA_FILE, 'rb') as f:
        pathway_metadata = pickle.load(f)

    # 2. Identify result files
    result_files = glob.glob(os.path.join(RESULTS_DISTANCES_P, "*.csv"))

    # Define Column Names
    ID_COL, Q_COL, W1_COL, DM_COL, N_MUT = 'pathway', 'q_value', 'wasserstein_distance', 'delta_means', 'num_mutations'

    for res_path in result_files:
        cancer_name = os.path.basename(res_path).replace(".csv", "")
        cancer_plot_dir = os.path.join(PLOTS_P, "p7_histogram_overlays", cancer_name)
        os.makedirs(cancer_plot_dir, exist_ok=True)

        results_df = pd.read_csv(res_path).dropna(subset=[Q_COL])

        # Filter for interesting ones to save massive amount of time/disk
        interesting_df = results_df[results_df[Q_COL] <= 0.2]

        print(f"\nCohort: {cancer_name.upper()} | Plotting {len(interesting_df)} pathways...")

        # Prepare tasks for workers
        tasks = []
        for _, row in interesting_df.iterrows():
            pw_id = str(row[ID_COL])
            pw_name = pathway_metadata.get(pw_id, {}).get('name', pw_id)
            tasks.append((
                pw_id, pw_name, row[Q_COL], row[W1_COL], row[DM_COL], row[N_MUT],
                cancer_name, cancer_plot_dir
            ))

        # 3. Multiprocessing Execution
        # We use a context manager to ensure proper cleanup of workers
        # max_workers=None defaults to CPU count
        with ProcessPoolExecutor() as executor:
            list(tqdm(executor.map(plot_single_pathway, tasks),
                      total=len(tasks),
                      desc=f"CPUs working on {cancer_name}"))


if __name__ == "__main__":
    plot_all_overlays_parallel()