# This script uses the coverage statistics calculated in p8 to check for systematic bias.
# It plots the "Coverage Ratio" (well-covered genes / total genes) against the Pathogenicity Shift (delta means)

from plot_boot import *
boot_plot_folder()

import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from definitions import RESULTS_DISTANCES_P, PLOTS_P, COLOR_MAP


def plot_reliability_analysis():
    all_dfs = []
    result_files = glob.glob(os.path.join(RESULTS_DISTANCES_P, "*.csv"))

    for f in result_files:
        df = pd.read_csv(f)
        df = df.dropna(subset=['q_value'])
        if 'num_covered_genes' in df.columns:
            all_dfs.append(df)

    if not all_dfs:
        print("Error: No coverage data found. Run p8_calc_coverage.py first.")
        return None

    master_df = pd.concat(all_dfs, ignore_index=True)

    # 2. Calculate Coverage Ratio
    master_df['coverage_ratio'] = master_df['num_covered_genes'] / master_df['num_genes']

    master_df['is_significant'] = master_df['q_value'] < 0.05

    # 3. Plotting Original Scatter
    plt.figure(figsize=(10, 7))

    # Masking using Boolean logic
    sns.scatterplot(data=master_df[~master_df['is_significant']],
                    x='delta_means', y='coverage_ratio',
                    alpha=0.17, color='grey', label='Non-Significant')

    sns.scatterplot(data=master_df[master_df['is_significant']],
                    x='delta_means', y='coverage_ratio',
                    alpha=0.6, color=COLOR_MAP['pathogenic'], label='Significant Hit (Q < 0.05)')

    plt.axvline(0, color='black', linestyle='--', alpha=0.35)
    plt.axhline(0.5, color='black', linestyle=':', alpha=0.15)
    plt.title("Reliability Analysis: Data Coverage vs. Pathogenicity Shift", fontsize=15)
    plt.xlabel("Delta Mean (Effect Size)", fontsize=12)
    plt.ylabel("Coverage Ratio (Covered Genes / Total Genes)", fontsize=12)
    plt.grid(True, linestyle=':', alpha=0.4)
    plt.legend(frameon=True)

    save_path = os.path.join(PLOTS_P, "p8_reliability_scatter.png")
    plt.savefig(save_path, dpi=600, bbox_inches='tight')
    print(f"Reliability plot saved to: {save_path}")

    return master_df


def plot_reliability_density_comparison(df):
    """
    Plots two 2D Density (KDE) maps side-by-side to compare the distribution
    of significant vs non-significant pathways.
    """
    if df is None:
        return

    fig, axes = plt.subplots(1, 2, figsize=(16, 7), sharex=True, sharey=True)

    titles = ["Non-Significant Pathways", "Significant Pathways (Q < 0.05)"]
    masks = [~df['is_significant'], df['is_significant']]
    colors = [COLOR_MAP['dark blue'], COLOR_MAP['pathogenic']]

    for i in range(2):
        ax = axes[i]
        subset = df[masks[i]]

        if subset.empty or len(subset) < 2:
            ax.text(0.5, 0.5, "Insufficient Data for Density", ha='center')
            continue

        # --- KDEPLOT ---
        sns.kdeplot(
            data=subset, x='delta_means', y='coverage_ratio',
            fill=True,
            thresh=0.05,
            levels=10,
            cut=0,
            cmap=sns.light_palette(colors[i], as_cmap=True),
            ax=ax
        )

        ax.scatter(subset['delta_means'], subset['coverage_ratio'],
                   s=1.2, alpha=0.15, color=colors[i])

        ax.set_title(titles[i], fontsize=14, fontweight='bold')
        ax.set_xlabel("Delta Mean", fontsize=12)

        # --- EXPLICIT AXIS LIMITS ---
        ax.set_ylim(0, 1.05)  # Coverage ratio is between 0 and 1

        ax.axvline(0, color='black', linestyle='--', alpha=0.3)
        ax.grid(True, linestyle=':', alpha=0.3)

    axes[0].set_ylabel("Coverage Ratio", fontsize=12)
    plt.suptitle("Density Distribution of Reliability Metrics", fontsize=18, y=1.02)
    plt.tight_layout()

    save_path_density = os.path.join(PLOTS_P, "p8_reliability_density_comparison.png")
    plt.savefig(save_path_density, dpi=600, bbox_inches='tight')
    print(f"Density comparison plot saved to: {save_path_density}")


if __name__ == "__main__":
    os.makedirs(PLOTS_P, exist_ok=True)
    master_df = plot_reliability_analysis()
    plot_reliability_density_comparison(master_df)