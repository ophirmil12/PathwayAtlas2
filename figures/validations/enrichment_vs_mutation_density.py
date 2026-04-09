from pathlib import Path
import sys
import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

PROJECT_ROOT = Path(__file__).parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from definitions import RESULTS_DISTANCES_P, FIGURES_P, COLOR_MAP, set_paper_palette

set_paper_palette()


def plot_mutation_density():
    print("collecting files...")

    # 1. Load all result files from distances (these contain q_value AND num_mutations)
    all_files = glob.glob(os.path.join(RESULTS_DISTANCES_P, "*.csv"))
    df_list = []

    for f in all_files:
        df = pd.read_csv(f)
        # Skip empty files to prevent errors
        if df.empty:
            continue
        df_list.append(df)

    if not df_list:
        print("Error: No valid data found in RESULTS_DISTANCES_P.")
        return

    # Combine all distance dataframes into one master dataframe
    summary_df = pd.concat(df_list, ignore_index=True)

    # 2. Filter and Calculate
    # Prevent DivisionByZero errors
    summary_df = summary_df[summary_df['num_covered_genes'] > 0].copy()

    # Calculate Mutation Density: Observed Cancer Mutations / Covered Genes
    summary_df['density'] = summary_df['num_mutations'] / summary_df['num_covered_genes']

    print("plotting...")
    plt.figure(figsize=(10, 6))

    sns.regplot(
        data=summary_df,
        x='density',
        y='q_value',
        fit_reg=True,
        scatter_kws={'alpha': 0.5, 'color': COLOR_MAP['dark-blue']},
        line_kws={'color': COLOR_MAP['red']}
    )

    # cut the y axis at 0..1
    plt.ylim(-0.05, 1.05)

    plt.title("Validation: Enrichment vs. Mutation Density", fontsize=14)
    plt.xlabel("Mutation Density (Mutations per Covered Gene)", fontsize=12)
    plt.ylabel("Enrichment (Q-Value)", fontsize=12)
    plt.grid(True, linestyle=':', alpha=0.6)

    plt.tight_layout()

    # Safely handle directory creation and save
    save_dir = os.path.join(FIGURES_P, "validations")
    os.makedirs(save_dir, exist_ok=True)

    save_path = os.path.join(save_dir, "enrichment_vs_mutation_density.png")
    plt.savefig(save_path, dpi=600)
    print(f"Plot saved to: {save_path}")


if __name__ == "__main__":
    plot_mutation_density()