# Plots the distributions of q-values and delta means (effect sizes)
# across pathways to visualize statistical significance and result trends.

from plot_boot import *
boot_plot_folder()

import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from definitions import RESULTS_DISTANCES_P, PLOTS_P, COLOR_MAP


def plot_q_val_and_delta_distributions():
    # 1. Load all CSV files from the results directory
    all_files = glob.glob(os.path.join(RESULTS_DISTANCES_P, "*.csv"))

    if not all_files:
        print(f"No CSV files found in {RESULTS_DISTANCES_P}")
        return

    df_list = []
    for f in all_files:
        try:
            temp_df = pd.read_csv(f)
            df_list.append(temp_df)
        except Exception as e:
            print(f"Error reading {f}: {e}")

    df = pd.concat(df_list, ignore_index=True)

    # Clean data: drop NaNs in critical columns if any
    df = df.dropna(subset=['q_value', 'delta_means'])

    # 2. Identify significant pathways
    significant_df = df[df['q_value'] < 0.05]

    # 3. Create Plotting Figure
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    plt.subplots_adjust(wspace=0.3)

    # --- Plot A: Q-Value Distribution ---
    sns.histplot(df['q_value'], bins=50, ax=axes[0], color=COLOR_MAP['dark-blue'], kde=True)
    axes[0].axvline(0.05, color=COLOR_MAP['pathogenic'], linestyle='--', label='Threshold (0.05)')
    axes[0].set_title('Distribution of Q-Values')
    axes[0].set_xlabel('Q-Value')
    axes[0].set_ylabel('Frequency')
    axes[0].legend()

    # --- Plot B: All Delta Means Distribution ---
    # Split data into positive and negative for distinct coloring
    pos_delta = df[df['delta_means'] > 0]['delta_means']
    neg_delta = df[df['delta_means'] <= 0]['delta_means']

    sns.histplot(neg_delta, bins=25, ax=axes[1], color=COLOR_MAP['benign'], alpha=0.8)
    sns.histplot(pos_delta, bins=25, ax=axes[1], color=COLOR_MAP['pathogenic'], alpha=0.8)

    axes[1].axvline(0, color='grey', linewidth=1, linestyle='--')
    axes[1].set_title('Distribution of All Delta Means Values')
    axes[1].set_xlabel('Delta Means')
    axes[1].set_ylabel('Frequency')

    # --- Plot C: Significant vs All Delta Means ---
    # Overlaying the significant ones over the background of all
    sns.histplot(df['delta_means'], bins=50, ax=axes[2], color=COLOR_MAP['grey'],
                 label='All Pathsways', alpha=0.3)
    sns.histplot(significant_df['delta_means'], bins=50, ax=axes[2], color=COLOR_MAP['significant'],
                 label='Significant (q < 0.05)', alpha=0.8)

    axes[2].set_title('Delta Means: Significant Rows vs All')
    axes[2].set_xlabel('Delta Means')
    axes[2].set_ylabel('Frequency')
    axes[2].legend()

    # 4. Save the plot
    output_path = os.path.join(PLOTS_P, 'p9_q_val_delta_means_dist.png')
    os.makedirs(PLOTS_P, exist_ok=True)
    plt.savefig(output_path, dpi=600, bbox_inches='tight')
    print(f"Plot saved to: {output_path}")

    # Optional: Print summary statistics
    print(f"Total entries: {len(df)}")
    print(f"Significant entries (q < 0.05): {len(significant_df)} ({len(significant_df) / len(df) * 100:.2f}%)")


if __name__ == "__main__":
    plot_q_val_and_delta_distributions()