from pathlib import Path
import sys
import os
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

PROJECT_ROOT = Path(__file__).parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from definitions import RESULTS_DISTANCES_P, FIGURES_P, COLOR_MAP, set_paper_palette

set_paper_palette()

def plot_pathogenic_enrichment():
    print("Collecting files...")
    all_files = glob.glob(os.path.join(RESULTS_DISTANCES_P, "*.csv"))
    df_list = []

    for f in all_files:
        df = pd.read_csv(f)
        if df.empty:
            continue
        df_list.append(df)

    if not df_list:
        print("Error: No valid data found.")
        return

    # Combine all distance dataframes
    summary_df = pd.concat(df_list, ignore_index=True)

    # 1. Prepare the Data
    q_threshold = 0.05

    # Handle q-values of 0 for log transformation
    min_q = summary_df[summary_df['q_value'] > 0]['q_value'].min()
    summary_df['q_value_safe'] = summary_df['q_value'].replace(0, min_q / 10)
    summary_df['neg_log_q'] = -np.log10(summary_df['q_value_safe'])

    # Categorize pathways for coloring
    def categorize(row):
        if row['q_value'] < q_threshold:
            return 'Pathogenic (Sig)' if row['delta_means'] > 0 else 'Benign (Sig)'
        else:
            return 'Not Significant'

    summary_df['Significance'] = summary_df.apply(categorize, axis=1)

    # Calculate the exact counts for the plot annotations
    pathogenic_sig_count = len(summary_df[summary_df['Significance'] == 'Pathogenic (Sig)'])
    benign_sig_count = len(summary_df[summary_df['Significance'] == 'Benign (Sig)'])

    # Define strict colors mapping based on your definitions
    palette = {
        'Pathogenic (Sig)': COLOR_MAP['pathogenic'],
        'Benign (Sig)': COLOR_MAP['benign'],
        'Not Significant': COLOR_MAP['grey']
    }

    # 2. Setup the Figure (1x2 Subplots)
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle("Pathways Tends Towards Pathogenicity", fontsize=16, y=1.05)

    # --- PANEL A: Density / Histogram ---
    sns.histplot(
        data=summary_df,
        x='delta_means',
        hue='Significance',
        palette=palette,
        multiple='stack',
        bins=50,
        edgecolor='white',
        ax=axes[0]
    )
    axes[0].axvline(x=0, color='black', linestyle='--', linewidth=1.5)
    axes[0].set_title("Distribution of Delta Means Scores", fontsize=14)
    axes[0].set_xlabel("Delta Means", fontsize=12)
    axes[0].set_ylabel("Number of Pathways", fontsize=12)
    axes[0].grid(True, linestyle=':', alpha=0.5)

    # --- PANEL B: Volcano Plot ---
    sns.scatterplot(
        data=summary_df,
        x='delta_means',
        y='neg_log_q',
        hue='Significance',
        palette=palette,
        alpha=0.7,
        edgecolor=None,
        ax=axes[1]
    )

    # Threshold lines for the volcano plot
    axes[1].axvline(x=0, color='black', linestyle='--', linewidth=1.5)
    axes[1].axhline(y=-np.log10(q_threshold), color=COLOR_MAP['grey'], linestyle='--', label='q = 0.05')

    axes[1].set_title("Volcano Plot of Pathway Enrichment", fontsize=14)
    axes[1].set_xlabel("Delta Means", fontsize=12)
    axes[1].set_ylabel("-log10(q-value)", fontsize=12)
    axes[1].grid(True, linestyle=':', alpha=0.5)

    # ---> NEW: Add Text Annotations to the Volcano Plot <---

    # Top-Left: Benign Count
    axes[1].text(0.05, 0.95, f"n = {benign_sig_count}",
                 transform=axes[1].transAxes, fontsize=12, fontweight='bold',
                 color=COLOR_MAP['benign'], ha='left', va='top',
                 bbox=dict(facecolor='white', alpha=0.8, edgecolor=COLOR_MAP['grey'], boxstyle='round,pad=0.4'))

    # Top-Right: Pathogenic Count
    axes[1].text(0.95, 0.95, f"n = {pathogenic_sig_count}",
                 transform=axes[1].transAxes, fontsize=12, fontweight='bold',
                 color=COLOR_MAP['pathogenic'], ha='right', va='top',
                 bbox=dict(facecolor='white', alpha=0.8, edgecolor=COLOR_MAP['grey'], boxstyle='round,pad=0.4'))

    plt.tight_layout()

    # Save the figure
    save_dir = os.path.join(FIGURES_P, "main_results")
    os.makedirs(save_dir, exist_ok=True)
    save_path = os.path.join(save_dir, "pathogenic_enrichment_distribution.png")

    plt.savefig(save_path, dpi=600, bbox_inches='tight')
    print(f"Plot saved to: {save_path}")

    plt.show()


if __name__ == "__main__":
    plot_pathogenic_enrichment()