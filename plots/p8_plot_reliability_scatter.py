# This script uses the coverage statistics calculated in p8 to check for systematic bias.
# It plots the "Coverage Ratio" (well-covered genes / total genes) against the Pathogenicity Shift (delta means)

from plot_boot import *
boot_plot_folder()

import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from definitions import RESULTS_DISTANCES_P, PLOTS_P, MY_PALETTE, set_paper_palette


def plot_reliability_analysis():
    set_paper_palette()

    # 1. Aggregate all results with coverage data
    all_dfs = []
    result_files = glob.glob(os.path.join(RESULTS_DISTANCES_P, "*.csv"))

    for f in result_files:
        df = pd.read_csv(f)
        if 'num_covered_genes' in df.columns:
            all_dfs.append(df)

    if not all_dfs:
        print("Error: No coverage data found. Run p8_calc_coverage.py first.")
        return

    master_df = pd.concat(all_dfs, ignore_index=True)

    # 2. Calculate Coverage Ratio
    # % of genes in the pathway that meet the 10-mut/1%-density threshold
    master_df['coverage_ratio'] = master_df['num_covered_genes'] / master_df['num_genes']
    master_df['is_significant'] = (master_df['q_value'] < 0.05).astype(str)

    # 3. Plotting
    plt.figure(figsize=(10, 7))

    # Plot non-significant points in light grey
    sns.scatterplot(data=master_df[master_df['is_significant'] == 'False'],
                    x='delta_mean', y='coverage_ratio',
                    alpha=0.15, color='grey', label='Non-Significant')

    # Plot significant hits in color
    sns.scatterplot(data=master_df[master_df['is_significant'] == 'True'],
                    x='delta_mean', y='coverage_ratio',
                    alpha=0.6, color=MY_PALETTE[0], label='Significant Hit (Q < 0.05)')

    # Add reference lines
    plt.axvline(0, color='black', linestyle='--', alpha=0.3)
    plt.axhline(0.5, color='black', linestyle=':', alpha=0.2)  # 50% coverage line

    # 4. Formatting
    plt.title("Reliability Analysis: Data Coverage vs. Pathogenicity Shift", fontsize=15)
    plt.xlabel("Delta Mean (Effect Size)", fontsize=12)
    plt.ylabel("Coverage Ratio (Covered Genes / Total Genes)", fontsize=12)

    # Annotation for Bias Check
    plt.text(-0.15, 0.05, "Low Data / High Noise Area", fontsize=9, color='red', alpha=0.6)
    plt.text(0.1, 0.9, "High Reliability Area", fontsize=9, color='green', alpha=0.6)

    plt.grid(True, linestyle=':', alpha=0.4)
    plt.legend(frameon=True)

    save_path = os.path.join(PLOTS_P, "p8_reliability_scatter.png")
    plt.savefig(save_path, dpi=600, bbox_inches='tight')
    print(f"Reliability plot saved to: {save_path}")
    plt.show()


if __name__ == "__main__":
    plot_reliability_analysis()


