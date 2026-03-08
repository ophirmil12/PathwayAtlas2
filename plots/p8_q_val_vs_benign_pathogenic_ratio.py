from plot_boot import *
boot_plot_folder()      # coloring scheme using cycler
import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from definitions import *

# Output folder
output_folder = os.path.join(PLOTS_P, "p8_q_val_vs_benign_pathogenic_ratio_plots")
os.makedirs(output_folder, exist_ok=True)

# Q-value thresholds: logarithmic scale from ~0.001 to 0.1, 30 points
q_thresholds = np.logspace(np.log10(0.001), np.log10(1), 30)

# Load all cancer CSVs
csv_files = glob.glob(os.path.join(RESULTS_DISTANCES_P, "*.csv"))

for csv_path in csv_files:
    cancer_type = os.path.splitext(os.path.basename(csv_path))[0]

    df = pd.read_csv(csv_path)

    ratios = []
    valid_thresholds = []

    for q in q_thresholds:
        filtered = df[df['q_value'] <= q]

        n_benign = (filtered['delta_means'] <= 0).sum()
        n_pathogenic = (filtered['delta_means'] > 0).sum()

        if n_pathogenic == 0:
            continue  # Skip / NaN

        ratio = n_benign / n_pathogenic
        ratios.append(ratio)
        valid_thresholds.append(q)

    if not valid_thresholds:
        print(f"Skipping {cancer_type}: no valid thresholds (always 0 pathogenic pathways).")
        continue

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(valid_thresholds, ratios, marker='o', linewidth=2, markersize=6)
    ax.set_xscale('log')
    ax.set_xlabel("Q-value threshold (log scale)", fontsize=13)
    ax.set_ylabel("Benign / Pathogenic ratio", fontsize=13)
    ax.set_title(f"{cancer_type} - Benign/Pathogenic ratio vs Q-value threshold", fontsize=14)
    ax.axvline(0.05, color='grey', linestyle='--', linewidth=1, label='Q=0.05')
    ax.legend(fontsize=11)
    sns.despine(ax=ax)
    plt.tight_layout()

    save_path = os.path.join(output_folder, f"{cancer_type}_q_val_vs_bp_ratio.png")
    plt.savefig(save_path, dpi=300)
    plt.close()
    print(f"Saved: {save_path}")

print("Done.")