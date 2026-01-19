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
    set_paper_palette()

    # 1. Aggregate all results
    all_dfs = []
    result_files = glob.glob(os.path.join(RESULTS_DISTANCES_P, "*.csv"))

    for f in result_files:
        cancer_name = os.path.basename(f).replace(".csv", "").upper()
        df = pd.read_csv(f)
        df['cancer'] = cancer_name
        all_dfs.append(df)

    master_df = pd.concat(all_dfs, ignore_index=True)

    # Required columns # TODO check column names
    X_COL = 'delta_mean'
    Q_COL = 'q_value_dw'
    PW_COL = 'pathway_name'

    # 2. Select top 50 pathways
    TOP_K = 50
    # We select pathways that appear as significant (Q < 0.05) in the most cancers
    sig_counts = master_df[master_df[Q_COL] < 0.05][PW_COL].value_counts()
    top_pw_names = sig_counts.head(TOP_K).index.tolist()

    # If not enough significant pathways, fill with those having largest absolute delta_mean
    if len(top_pw_names) < TOP_K:
        remaining = TOP_K - len(top_pw_names)
        other_pws = master_df[~master_df[PW_COL].isin(top_pw_names)].groupby(PW_COL)[X_COL].apply(
            lambda x: x.abs().mean()).sort_values(ascending=False)
        top_pw_names.extend(other_pws.head(remaining).index.tolist())

    # 3. Pivot data for Heatmap
    plot_df = master_df[master_df[PW_COL].isin(top_pw_names)]
    heatmap_val = plot_df.pivot(index=PW_COL, columns='cancer', values=X_COL).fillna(0)
    heatmap_sig = plot_df.pivot(index=PW_COL, columns='cancer', values=Q_COL).fillna(1.0)

    # 4. Plotting
    plt.figure(figsize=(20, 15))

    # Diverging colormap: Blue (Benign) -> White -> Red (Pathogenic)
    # center=0 ensures 0 delta-mean is white
    ax = sns.heatmap(heatmap_val, cmap="RdBu_r", center=0,
                     cbar_kws={'label': 'Delta Mean ($\Delta\mu$)'},
                     linewidths=.5, linecolor='gray')

    # 5. Add stars for significance (Q < 0.05)
    for i in range(len(heatmap_val.index)):
        for j in range(len(heatmap_val.columns)):
            q_val = heatmap_sig.iloc[i, j]
            if q_val < 0.05:
                # Add a small white dot/star in the center of significant cells
                ax.text(j + 0.5, i + 0.5, '●',
                        ha='center', va='center', color='white', fontsize=8)

    plt.title(f"Pan-Cancer Pathway Pathogenicity: Top {TOP_K} Frequent Hits", fontsize=20, pad=20)
    plt.xlabel("Cancer Type", fontsize=14)
    plt.ylabel("Pathway Name", fontsize=14)

    os.makedirs(PLOTS_P, exist_ok=True)
    plt.savefig(os.path.join(PLOTS_P, "p7_pancancer_heatmap.png"), dpi=600, bbox_inches='tight')
    print(f"Heatmap saved to {PLOTS_P}")
    plt.show()


if __name__ == "__main__":
    plot_pancancer_heatmap()

