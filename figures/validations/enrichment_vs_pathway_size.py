from pathlib import Path
import sys
import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

PROJECT_ROOT = Path(__file__).parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from definitions import RESULTS_DISTANCES_P, COLOR_MAP, set_paper_palette, FIGURES_P

set_paper_palette()


def plot_coverage_independence():
    # 1. Load all result files and combine them
    all_files = glob.glob(os.path.join(RESULTS_DISTANCES_P, "*.csv"))
    df_list = []

    for f in all_files:
        df = pd.read_csv(f)
        if df.empty:
            continue
        df_list.append(df)

    if not df_list:
        print("Error: No valid data found in RESULTS_DISTANCES_P.")
        return

    # Combine all cancer/distance dataframes into one master dataframe
    summary_df = pd.concat(df_list, ignore_index=True)

    # 2. Calculate Coverage
    # Filter out any pathways with 0 genes to avoid division by zero
    summary_df = summary_df[(summary_df['num_genes'] > 0) & (summary_df['num_genes'] < 1500)].copy()

    # 3. Plotting
    plt.figure(figsize=(10, 6))
    sns.regplot(
        data=summary_df,
        x='num_genes',
        y='q_value',
        fit_reg=True,
        scatter_kws={'alpha': 0.5, 'color': COLOR_MAP['dark-blue']},
    )

    # cut the y axis at 0..1
    plt.ylim(-0.05, 1.05)

    plt.title("Validation: Enrichment vs. Pathway's Number of Genes", fontsize=14)
    plt.xlabel("Number of Genes", fontsize=12)
    plt.ylabel("Enrichment (Q-Value)", fontsize=12)
    plt.legend()
    plt.grid(True, linestyle=':', alpha=0.6)

    plt.tight_layout()

    # Ensure directory exists before saving
    save_dir = os.path.join(FIGURES_P, "validations")
    os.makedirs(save_dir, exist_ok=True)

    plt.savefig(os.path.join(save_dir, "enrichment_vs_pathway_size.png"), dpi=600)


if __name__ == "__main__":
    plot_coverage_independence()