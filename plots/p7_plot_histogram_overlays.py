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
    MICHAL_HN1_PSSM,
)
from distance_utils import get_bg_histogram_after_pssm, get_cancer_histogram


def plot_single_pathway(args):
    """
    Worker function to plot a single pathway.
    'args' is a tuple to support ProcessPoolExecutor.submit.
    """
    (pw_id, pw_name, q_val, w1_dist, d_mean, cancer_name, cancer_plot_dir) = args

    # Define bins and centers locally in worker
    bin_edges = np.linspace(0, 1, NUMBER_OF_BINS + 1)
    bin_centers = bin_edges[:-1] + np.diff(bin_edges) / 2

    try:
        # 1. Generate Histograms
        hist_bg = get_bg_histogram_after_pssm(pw_id, pssm_matrix=MICHAL_HN1_PSSM)

        hist_cancer = get_cancer_histogram(pw_id, f"{cancer_name}.csv")

        if np.sum(hist_cancer) == 0:
            return False

        # 2. Plotting
        plt.figure(figsize=(10, 6.5))

        # --- BACKGROUND (Expected) ---
        # fill=True creates the area-under-the-curve look
        plt.stairs(hist_bg, bin_edges,
                   color=COLOR_MAP["dark blue"], fill=True, alpha=0.25,
                   label='Expected (PSSM-Weighted BG)')

        # Add a thin outline to make the steps clear
        plt.stairs(hist_bg, bin_edges,
                   color=COLOR_MAP["dark blue"], linewidth=1, alpha=0.6)

        # --- OBSERVED CANCER ---
        color_choice = COLOR_MAP["pathogenic"] if d_mean > 0 else COLOR_MAP["benign"]

        # Plot with fill=True for the "Binned PDF" look
        plt.stairs(hist_cancer, bin_edges,
                   color=color_choice, fill=True, alpha=0.4,
                   label=f'Observed ({cancer_name.upper()})')

        # Add a thicker outline to emphasize the cancer distribution steps
        plt.stairs(hist_cancer, bin_edges,
                   color=color_choice, linewidth=2.5)

        # 3. Formatting & Stats Box
        sig_status = "SIGNIFICANT" if q_val < 0.05 else "Not Significant"
        stats_text = (
            f"Q-value: {q_val:.2e}\n"
            f"Status:  {sig_status}\n"
            f"W1 Dist: {w1_dist:.4f}\n"
            f"Δ Mean:  {d_mean:+.4f}"
        )

        plt.annotate(stats_text, xy=(0.02, 0.96), xycoords='axes fraction',
                     va='top', ha='left', fontsize=10, family='monospace',
                     bbox=dict(boxstyle="round,pad=0.5", fc="white", ec="black", alpha=0.8))

        # Text wrapping for title
        wrapped_name = "\n".join(textwrap.wrap(pw_name.split(' - Homo sapiens')[0], width=50))
        plt.title(f"{wrapped_name}\n({pw_id})", fontsize=15, fontweight='bold', pad=15)

        plt.xlabel("Pathogenicity Probability", fontsize=12)
        plt.ylabel("Probability Mass", fontsize=12)
        plt.xlim(0, 1)

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
    ID_COL, Q_COL, W1_COL, DM_COL = 'pathway', 'q_value', 'wasserstein_distance', 'delta_means'

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
                pw_id, pw_name, row[Q_COL], row[W1_COL], row[DM_COL],
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