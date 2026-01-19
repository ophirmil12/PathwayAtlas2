# ESM Score Separation (ClinVar Ground Truth)
#     Plot: Density plot of esm_log_probs grouped by ClinVar labels (Benign vs. Pathogenic).
#     Insight: If the two curves overlap completely, ESM isn't helping. If they are well-separated, your foundation is strong.

from plot_boot import *
boot_plot_folder()

import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import roc_auc_score

from definitions import (
    CLINVAR_DATA_TABLE_P,
    PLOTS_P,
    MY_PALETTE,
    set_paper_palette
)


def validate_esm_on_clinvar():
    """
    Plots the distribution of ESM scores for Benign vs Pathogenic ClinVar mutations.
    Separates by Ordered vs Disordered regions.
    """
    set_paper_palette()

    # 1. Load ClinVar Data
    if not os.path.exists(CLINVAR_DATA_TABLE_P):
        print(f"Error: ClinVar data not found at {CLINVAR_DATA_TABLE_P}")
        return

    df = pd.read_csv(CLINVAR_DATA_TABLE_P)

    # Column names based on p5_clinvar_reggresors.py
    esm_col = 'wt_not_nadav_marginals_base_wt_score'
    label_col = 'binary_label'  # 0=Benign, 1=Pathogenic
    disorder_col = 'is_disordered'

    # Clean data
    df = df.dropna(subset=[esm_col, label_col])

    # 2. Setup Plot
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True)
    titles = {0: "Ordered Regions (is_disordered=0)", 1: "Disordered Regions (is_disordered=1)"}

    # Define colors for Benign and Pathogenic
    # Benign: Green (#447D68), Pathogenic: Red (#CB7673)
    color_map = {0: MY_PALETTE[1], 1: MY_PALETTE[0]}
    label_map = {0: "Benign", 1: "Pathogenic"}

    for flag in [0, 1]:
        ax = axes[flag]
        subset = df[df[disorder_col] == flag]

        if subset.empty:
            continue

        # Calculate AUC for this subset to show in title
        auc = roc_auc_score(subset[label_col], -subset[esm_col])  # Negate because lower score = more pathogenic

        # Plot Distributions
        for label_val in [0, 1]:
            data = subset[subset[label_col] == label_val][esm_col]
            sns.kdeplot(
                data,
                ax=ax,
                fill=True,
                color=color_map[label_val],
                label=f"{label_map[label_val]} (N={len(data)})",
                linewidth=2,
                alpha=0.4
            )

        ax.set_title(f"{titles[flag]}\nSeparation AUC: {auc:.3f}", fontsize=13, fontweight='bold')
        ax.set_xlabel("ESM Log-Prob (Lower = More Surprising/Pathogenic)", fontsize=11)
        ax.set_ylabel("Density", fontsize=11)
        ax.legend()
        ax.grid(axis='y', linestyle='--', alpha=0.3)

    plt.suptitle("ESM Score Separation: Benign vs Pathogenic (ClinVar)", fontsize=16, y=1.02)
    plt.tight_layout()

    # 3. Save result
    os.makedirs(PLOTS_P, exist_ok=True)
    save_path = os.path.join(PLOTS_P, "p5_clinvar_validate_esm_usage.png")
    plt.savefig(save_path, dpi=600, bbox_inches='tight')

    print(f"\nValidation plot saved to: {save_path}")
    print(f"Total mutations analyzed: {len(df)}")
    plt.show()


if __name__ == "__main__":
    validate_esm_on_clinvar()


