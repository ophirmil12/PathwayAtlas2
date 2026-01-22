# This script aggregates all cancer results,
# identifies the top 50 pathways based on how frequently they are significant across different cancers,
# and creates a diverging heatmap.

from plot_boot import *
boot_plot_folder()

import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from definitions import RESULTS_DISTANCES_P, PLOTS_P, set_paper_palette


def plot_pancancer_heatmap():
    # Identify columns
    PW_COL = 'pathway'
    X_COL = 'delta_means'
    Q_COL = 'q_value'

    # 1. Aggregate all results
    all_dfs = []
    result_files = glob.glob(os.path.join(RESULTS_DISTANCES_P, "*.csv"))

    for f in result_files:
        cancer_name = os.path.basename(f).replace(".csv", "").upper()
        df = pd.read_csv(f)

        # Filter and clean
        df = df.dropna(subset=[Q_COL])
        df['cancer'] = cancer_name
        all_dfs.append(df[[PW_COL, X_COL, Q_COL, 'cancer']])

    master_df = pd.concat(all_dfs, ignore_index=True)

    # Pre-emptive cleanup of exact duplicates
    master_df = master_df.drop_duplicates(subset=[PW_COL, 'cancer'])

    # 2. Select top 50 pathways
    TOP_K = 50
    sig_counts = master_df[master_df[Q_COL] < 0.05][PW_COL].value_counts()
    top_pw_names = sig_counts.head(TOP_K).index.tolist()

    if len(top_pw_names) < TOP_K:
        remaining = TOP_K - len(top_pw_names)
        other_pws = master_df[~master_df[PW_COL].isin(top_pw_names)].groupby(PW_COL)[X_COL].apply(
            lambda x: x.abs().mean()).sort_values(ascending=False)
        top_pw_names.extend(other_pws.head(remaining).index.tolist())

    # 3. Pivot data for Heatmap (Using pivot_table for safety)
    plot_df = master_df[master_df[PW_COL].isin(top_pw_names)]

    heatmap_val = plot_df.pivot_table(index=PW_COL, columns='cancer', values=X_COL, aggfunc='mean').fillna(0)
    heatmap_sig = plot_df.pivot_table(index=PW_COL, columns='cancer', values=Q_COL, aggfunc='mean').fillna(1.0)

    # 4. Plotting
    plt.figure(figsize=(22, 16))  # Slightly larger for readability

    ax = sns.heatmap(heatmap_val, cmap="RdBu_r", center=0,
                     cbar_kws={'label': 'Delta Mean ($\Delta\mu$)'},
                     linewidths=.5, linecolor='gray')

    plt.title(f"Pan-Cancer Pathway Pathogenicity: Top {TOP_K} Frequent Hits", fontsize=20, pad=20)
    plt.xlabel("Cancer Type", fontsize=14)
    plt.ylabel("Pathway Name", fontsize=14)

    os.makedirs(PLOTS_P, exist_ok=True)
    save_path = os.path.join(PLOTS_P, "p7_pancancer_heatmap.png")
    plt.savefig(save_path, dpi=600, bbox_inches='tight')
    plt.close()
    print(f"Heatmap successfully saved to {save_path}")


if __name__ == "__main__":
    plot_pancancer_heatmap()

