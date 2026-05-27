from plot_boot import *
boot_plot_folder()

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import os
import re
from scipy.stats import pearsonr
import matplotlib.gridspec as gridspec
from definitions import *


def plot_bin_validation(results_dir):
    # 1. Load and merge data
    all_files = glob.glob(os.path.join(results_dir, "pan_cancer_bins_*.csv"))
    if not all_files:
        print(f"No files found in {results_dir}")
        return

    dfs = []
    for f in all_files:
        df = pd.read_csv(f)
        if 'bins' not in df.columns:
            match = re.search(r'bins_(\d+)', f)
            if match:
                df['bins'] = int(match.group(1))
        df['bins'] = df['bins'].astype(int)
        dfs.append(df)

    full_df = pd.concat(dfs, ignore_index=True)
    dist_pivot = full_df.pivot(index='pathway', columns='bins', values='wasserstein_distance')
    q_pivot = full_df.pivot(index='pathway', columns='bins', values='q_value')

    dist_pivot = dist_pivot.reindex(sorted(dist_pivot.columns), axis=1)
    q_pivot = q_pivot.reindex(sorted(q_pivot.columns), axis=1)

    sns.set_style("ticks")

    # --- FIGURE 1: Wasserstein Distance Stability ---
    if 100 in dist_pivot.columns:
        fig1, ax_dist = plt.subplots(figsize=(10, 8))
        top_pathways = dist_pivot.sort_values(by=100, ascending=False).head(12).index
        colors = sns.color_palette("husl", len(top_pathways))

        for i, pathway in enumerate(top_pathways):
            ax_dist.plot(dist_pivot.columns, dist_pivot.loc[pathway],
                         marker='o', linewidth=2, alpha=0.8, color=colors[i], label=pathway)

        ax_dist.set_title("Wasserstein Distance Stability", fontsize=18, fontweight='bold', pad=20)
        ax_dist.set_xlabel("Number of Bins", fontsize=14)
        ax_dist.set_ylabel("Distance ($W_1$)", fontsize=14)
        ax_dist.set_xticks(dist_pivot.columns)
        ax_dist.legend(title="Top Pathways", bbox_to_anchor=(0.5, -0.15), loc='upper center', ncol=2, fontsize=9)

        sns.despine()
        plt.savefig(pjoin(PLOTS_P, "p7_bin_validation_wasserstein.png"), dpi=300, bbox_inches='tight')
        plt.close(fig1)

    # --- FIGURE 2: Q-value Comparisons (Grid) ---
    if 100 in q_pivot.columns:
        fig2, axes = plt.subplots(2, 2, figsize=(12, 10), gridspec_kw={'hspace': 0.3, 'wspace': 0.3})
        comparison_bins = [10, 25, 50, 75]
        target_q = -np.log10(q_pivot[100] + 1e-4)

        for i, b_size in enumerate(comparison_bins):
            ax_sub = axes[i // 2, i % 2]
            if b_size in q_pivot.columns:
                current_q = -np.log10(q_pivot[b_size] + 1e-4)
                mask = q_pivot[100].notna() & q_pivot[b_size].notna()
                r, _ = pearsonr(current_q[mask], target_q[mask])

                sns.regplot(x=target_q, y=current_q, ax=ax_sub,
                            scatter_kws={'alpha': 0.3, 's': 20, 'color': COLOR_MAP['dark-blue']},
                            line_kws={'color': COLOR_MAP['red'], 'label': f'R={r:.3f}'})

                ax_sub.set_title(f"{b_size} vs 100 Bins", fontsize=12, fontweight='bold')
                ax_sub.set_xlabel("$-log_{10}$ Q (100 Bins)", fontsize=10)
                ax_sub.set_ylabel("$-log_{10}$ Q", fontsize=10)
                ax_sub.legend(loc='lower right', fontsize=10, frameon=False)
                ax_sub.set_xlim(-0.1, 4.2)
                ax_sub.set_ylim(-0.1, 4.2)

        fig2.suptitle("Q-Value Stability Comparison", fontsize=18, fontweight='bold', y=0.95)
        sns.despine()
        plt.savefig(pjoin(PLOTS_P, "p7_bin_validation_qvalue.png"), dpi=300, bbox_inches='tight')
        plt.close(fig2)

    # --- FIGURE 3: Pathway Discovery Count ---
    fig3, ax_count = plt.subplots(figsize=(10, 8))
    sig_05 = (q_pivot < 0.05).sum()
    sig_01 = (q_pivot < 0.01).sum()
    x_indices = np.arange(len(sig_05))
    width = 0.35

    ax_count.bar(x_indices - width / 2, sig_05, width, label='q < 0.05', color=COLOR_MAP['light-blue'],
                 edgecolor='black')
    ax_count.bar(x_indices + width / 2, sig_01, width, label='q < 0.01', color=COLOR_MAP['significant'],
                 edgecolor='black')

    ax_count.set_title("Discovery Sensitivity", fontsize=18, fontweight='bold', pad=20)
    ax_count.set_xlabel("Number of Bins", fontsize=14)
    ax_count.set_ylabel("Significant Pathways Count", fontsize=14)
    ax_count.set_xticks(x_indices)
    ax_count.set_xticklabels(sig_05.index)
    ax_count.legend(loc='upper left', frameon=False)

    for i, v in enumerate(sig_05):
        ax_count.text(i - width / 2, v + 1, str(v), ha='center', va='bottom', fontsize=11, fontweight='bold')
    for i, v in enumerate(sig_01):
        ax_count.text(i + width / 2, v + 1, str(v), ha='center', va='bottom', fontsize=11)

    sns.despine()
    plt.savefig(pjoin(PLOTS_P, "p7_bin_validation_discovery.png"), dpi=300, bbox_inches='tight')
    plt.close(fig3)

    print(f"Saved 3 individual validation plots to {PLOTS_P}")


if __name__ == "__main__":
    plot_bin_validation(pjoin(RESULTS_P, "p7_bin_validation_analysis"))