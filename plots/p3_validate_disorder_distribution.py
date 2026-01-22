# Histogram of Disorder_score for all human genes.
# should see a bimodal distribution (one peak for ordered, one for disordered).
# Verify that our DISORDERED_THRESHOLD (0.7) sits in a logical "valley."

from plot_boot import *
boot_plot_folder()

import os
import glob
import random
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm

from definitions import (
    KEGG_GENE_SCORES_P,
    PLOTS_P,
    DISORDERED_THRESHOLD,
    COLOR_MAP
)


def plot_disorder_landscape():
    """
    Plots the global distribution of disorder scores across sampled human genes.
    Validates the bimodal nature of protein structure and the 0.7 threshold.
    """
    # 1. Sample Gene Files
    # 2000 genes provide millions of residues, more than enough for a smooth distribution
    all_files = glob.glob(os.path.join(KEGG_GENE_SCORES_P, "*.csv"))
    if not all_files:
        print(f"Error: No files found in {KEGG_GENE_SCORES_P}")
        return

    sample_files = random.sample(all_files, min(2000, len(all_files)))

    disorder_scores = []

    print(f"Reading disorder scores from {len(sample_files)} genes...")
    for f in tqdm(sample_files):
        try:
            # We only need AA_index (to deduplicate) and the score
            df = pd.read_csv(f, usecols=['AA_index', 'Disorder_score'])

            if 'Disorder_score' not in df.columns:
                continue

            # IMPORTANT: Deduplicate by AA_index.
            # One residue has one disorder score, but multiple potential SNVs.
            # We want the distribution of residues.
            gene_residue_scores = df.drop_duplicates('AA_index')['Disorder_score'].dropna()
            disorder_scores.extend(gene_residue_scores.tolist())
        except Exception:
            continue

    if not disorder_scores:
        print("No valid disorder scores found. Ensure p3_disorder_v3_... has run.")
        return

    # 2. Plotting
    plt.figure(figsize=(10, 6))

    # Use COLOR_MAP['pink'] (Pink-ish #F3B8BA)
    main_color = COLOR_MAP['pink']

    sns.histplot(disorder_scores, bins=50, kde=True, color=main_color, edgecolor='white', alpha=0.7)

    # 3. Add Threshold Annotation
    plt.axvline(x=DISORDERED_THRESHOLD, color=COLOR_MAP['red'], linestyle='--', linewidth=2.5,
                label=f'Disorder Threshold ({DISORDERED_THRESHOLD})')

    # Shade the areas
    plt.axvspan(0, DISORDERED_THRESHOLD, alpha=0.1, color='blue', label='Ordered Region')
    plt.axvspan(DISORDERED_THRESHOLD, 1.0, alpha=0.1, color='orange', label='Disordered Region')

    # 4. Formatting
    plt.title("Global Proteome Disorder Distribution", fontsize=15, fontweight='bold')
    plt.xlabel("Metapredict V3 Disorder Score", fontsize=12)
    plt.ylabel("Count of Residues", fontsize=12)
    plt.xlim(0, 1)
    plt.grid(axis='y', linestyle=':', alpha=0.5)
    plt.legend()

    # 5. Save result
    os.makedirs(PLOTS_P, exist_ok=True)
    save_path = os.path.join(PLOTS_P, "p3_validate_disorder_distribution.png")
    plt.savefig(save_path, dpi=600, bbox_inches='tight')

    print(f"\nDisorder distribution plot saved to: {save_path}")

    # Print some stats for verification
    scores_array = np.array(disorder_scores)
    perc_disordered = (scores_array >= DISORDERED_THRESHOLD).sum() / len(scores_array) * 100
    print(f"Total residues sampled: {len(scores_array):,}")
    print(f"Percentage of residues considered disordered: {perc_disordered:.2f}%")

    plt.show()


if __name__ == "__main__":
    plot_disorder_landscape()