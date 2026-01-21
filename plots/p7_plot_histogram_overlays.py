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
from tqdm import tqdm

from definitions import (
    RESULTS_DISTANCES_P,
    KEGG_PATHWAY_METADATA_P,
    PLOTS_P,
    MY_PALETTE,
    NUMBER_OF_BINS,
    MICHAL_HN1_PSSM,
    set_paper_palette
)
from distance_utils import get_bg_histogram_after_pssm, get_cancer_histogram


def plot_all_overlays():
    """
    Generates overlay plots for all pathways in all cancer types.
    Uses metadata to provide descriptive titles.
    """
    set_paper_palette()

    # 1. Load Pathway Metadata for Name Lookups
    if not os.path.exists(KEGG_PATHWAY_METADATA_P):
        print(f"Error: Metadata not found at {KEGG_PATHWAY_METADATA_P}")
        return

    with open(KEGG_PATHWAY_METADATA_P, 'rb') as f:
        pathway_metadata = pickle.load(f)

    # 2. Identify result files
    result_files = glob.glob(os.path.join(RESULTS_DISTANCES_P, "*.csv"))
    if not result_files:
        print(f"No results found in {RESULTS_DISTANCES_P}")
        return

    # Define bins
    bin_edges = np.linspace(0, 1, NUMBER_OF_BINS + 1)
    bin_centers = bin_edges[:-1] + np.diff(bin_edges) / 2

    for res_path in result_files:
        cancer_name = os.path.basename(res_path).replace(".csv", "")
        # The cancer mutation filename should match the result filename
        cancer_mut_file = f"{cancer_name}.csv"

        cancer_plot_dir = os.path.join(PLOTS_P, "p7_histogram_overlays", cancer_name)
        os.makedirs(cancer_plot_dir, exist_ok=True)

        results_df = pd.read_csv(res_path)

        ID_COL = 'pathway'
        Q_COL = 'q_value'
        W1_COL = 'wasserstein_distance'
        DM_COL = 'delta_means'

        print(f"\nProcessing Cohort: {cancer_name.upper()}")

        for _, row in tqdm(results_df.iterrows(), total=len(results_df), desc=f"Plotting {cancer_name}"):
            pw_id = str(row[ID_COL])
            q_val = row[Q_COL]
            w1_dist = row[W1_COL]
            d_mean = row[DM_COL]

            # 3. Lookup Real Name from Metadata
            # Fallback to ID if name is missing
            raw_name = pathway_metadata.get(pw_id, {}).get('name', pw_id)
            # Remove the ' - Homo sapiens (human)' suffix if present to keep titles clean
            clean_name = raw_name.split(' - Homo sapiens')[0]
            # Wrap long names so they don't run off the plot
            wrapped_name = "\n".join(textwrap.wrap(clean_name, width=50))

            # 4. Generate Histograms
            hist_bg = get_bg_histogram_after_pssm(pw_id, pssm_matrix=MICHAL_HN1_PSSM)
            hist_cancer = get_cancer_histogram(pw_id, cancer_mut_file)

            if np.sum(hist_cancer) == 0:
                continue

                # 5. Plotting
            plt.figure(figsize=(10, 6.5))

            # Background (PSSM-Expected)
            plt.fill_between(bin_centers, hist_bg, color=MY_PALETTE[3], alpha=0.3, label='Expected (PSSM-Weighted BG)')
            plt.plot(bin_centers, hist_bg, color=MY_PALETTE[3], linewidth=1.5, alpha=0.7)

            # Observed Cancer
            # Color logic: Red if pathogenic shift (>0), Green if benign (<0)
            color_choice = MY_PALETTE[0] if d_mean > 0 else MY_PALETTE[1]
            plt.fill_between(bin_centers, hist_cancer, color=color_choice, alpha=0.4,
                             label=f'Observed ({cancer_name.upper()})')
            plt.plot(bin_centers, hist_cancer, color=color_choice, linewidth=2.5)

            # 6. Formatting & Stats Box
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

            plt.title(f"{wrapped_name}\n({pw_id})", fontsize=15, fontweight='bold', pad=15)
            plt.xlabel("Pathogenicity Probability", fontsize=12)
            plt.ylabel("Probability Mass", fontsize=12)
            plt.xlim(0, 1)

            # Set Y-max dynamically with some headroom
            y_max = max(hist_bg.max(), hist_cancer.max()) * 1.3
            plt.ylim(0, y_max)

            plt.legend(loc='upper right', frameon=True)
            plt.grid(axis='y', linestyle=':', alpha=0.4)

            # 7. Save
            clean_filename = pw_id.replace(":", "_")
            save_path = os.path.join(cancer_plot_dir, f"{clean_filename}_overlay.png")
            plt.savefig(save_path, dpi=200, bbox_inches='tight')
            plt.close()


if __name__ == "__main__":
    plot_all_overlays()