# Volcano plot of the results

from plot_boot import *
boot_plot_folder()

import os
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

from definitions import RESULTS_DISTANCES_P, PLOTS_P, MY_PALETTE, set_paper_palette


def plot_volcano_matched_pairs():
    """
    Generates a pair of matched plots for each cancer result file.
    Plot 1: Delta Mean vs -log10(Q)
    Plot 2: Wasserstein Distance (W1) vs -log10(Q)
    Both plots are colored by the sign of Delta Mean.
    """
    set_paper_palette()

    # 1. Setup paths
    result_files = glob.glob(os.path.join(RESULTS_DISTANCES_P, "*.csv"))
    output_dir = os.path.join(PLOTS_P, "p7_volcano_plots_matched")
    os.makedirs(output_dir, exist_ok=True)

    if not result_files:
        print(f"No results found in {RESULTS_DISTANCES_P}")
        return

    # Column Mapping
    X1_COL = 'delta_means'  # Directional shift
    X2_COL = 'wasserstein_distance'  # Absolute magnitude (Standard Wasserstein)
    Q_COL = 'q_value'
    LABEL_COL = 'pathway_name'

    NUM_LABELS_FOR_TOP_K_HITS = 10
    Q_THRESHOLD = 0.05
    NEG_LOG_Q_THRESHOLD = -np.log10(Q_THRESHOLD)

    for file_path in tqdm(result_files, desc="Generating Matched Volcano Plots"):
        cancer_name = os.path.basename(file_path).replace(".csv", "").upper()
        df = pd.read_csv(file_path)
        
        # filter out rows with missing q_values
        df = df.dropna(subset=[Q_COL])

        # Check if necessary columns exist
        required = [X1_COL, X2_COL, Q_COL, LABEL_COL]
        if any(col not in df.columns for col in required):
            print(f"Skipping {cancer_name}: Missing required columns in {required}")
            continue

        # 2. Prepare Data
        # Calculate -log10(Q) with safety epsilon
        df['neg_log_q'] = -np.log10(df[Q_COL].replace(0, 1e-300))

        # Define color categories
        # Non-significant: grey
        # Significant + Positive Delta: Red (#CB7673)
        # Significant + Negative Delta: Green (#447D68)

        is_sig = df[Q_COL] < Q_THRESHOLD
        is_pos = df[X1_COL] > 0

        # 3. Create the matched subplots
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 8), sharey=True)

        # Helper plotting function to ensure consistency
        def apply_scatter(ax, x_col):
            # Non-significant
            mask_ns = ~is_sig
            ax.scatter(df.loc[mask_ns, x_col], df.loc[mask_ns, 'neg_log_q'],
                       alpha=0.2, color='grey', label='Not Significant')

            # Significant Pathogenic (Positive shift)
            mask_patho = is_sig & is_pos
            ax.scatter(df.loc[mask_patho, x_col], df.loc[mask_patho, 'neg_log_q'],
                       alpha=0.8, color=MY_PALETTE[0], label='Significant Pathogenic ($\Delta\mu > 0$)')

            # Significant Benign (Negative shift)
            mask_benign = is_sig & (~is_pos)
            ax.scatter(df.loc[mask_benign, x_col], df.loc[mask_benign, 'neg_log_q'],
                       alpha=0.8, color=MY_PALETTE[1], label='Significant Benign ($\Delta\mu < 0$)')

            # Labels for top NUM_LABELS_FOR_TOP_K_HITS hits
            # top_hits = df.nsmallest(NUM_LABELS_FOR_TOP_K_HITS, Q_COL)
            # for _, row in top_hits.iterrows():
            #     if row[Q_COL] < Q_THRESHOLD:
            #         ax.text(row[x_col], row['neg_log_q'] + 0.1, row[LABEL_COL],
            #                 fontsize=8, fontweight='bold', ha='center', va='bottom', alpha=0.7)

            # Threshold lines
            ax.axhline(y=NEG_LOG_Q_THRESHOLD, color='black', linestyle='--', alpha=0.4)
            ax.grid(True, which='both', linestyle=':', alpha=0.3)

        # Plot A: Delta Mean
        apply_scatter(ax1, X1_COL)
        ax1.axvline(x=0, color='black', linestyle='-', alpha=0.2)
        ax1.set_title(f"Directional Shift ($\Delta\mu$)", fontsize=14)
        ax1.set_xlabel("Delta Mean (Cancer $\mu$ - Background $\mu$)", fontsize=12)
        ax1.set_ylabel("$-log_{10}(Q\text{-value})$", fontsize=12)

        # Plot B: Wasserstein Magnitude
        apply_scatter(ax2, X2_COL)
        ax2.set_title(f"Absolute Distance ($W_1$)", fontsize=14)
        ax2.set_xlabel("Wasserstein Distance (Total Distribution Change)", fontsize=12)
        ax2.legend(loc='upper right', fontsize=9)

        plt.suptitle(f"Matched Pathway Pathogenicity Analysis: {cancer_name} cancer", fontsize=18, y=1.02)
        plt.tight_layout()

        # Save
        save_path = os.path.join(output_dir, f"matched_volcano_{cancer_name.lower()}.png")
        plt.savefig(save_path, dpi=600, bbox_inches='tight')
        plt.close()

    print(f"\nAll matched volcano plots saved to: {output_dir}")


if __name__ == "__main__":
    plot_volcano_matched_pairs()