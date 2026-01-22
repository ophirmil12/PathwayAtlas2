from plot_boot import *
boot_plot_folder()

import os
import glob
import random
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from tqdm import tqdm
from statsmodels.stats.multicomp import pairwise_tukeyhsd

from definitions import KEGG_GENE_SCORES_P, COLOR_MAP


def validate_codon_positions():
    # 1. Get a sample of files
    all_files = glob.glob(os.path.join(KEGG_GENE_SCORES_P, "*.csv"))
    # Sampling 2000 files provides high statistical power
    sample_files = random.sample(all_files, min(2000, len(all_files)))

    collector = []
    print(f"Reading {len(sample_files)} files...")

    for f in tqdm(sample_files, desc="Processing CSVs"):
        df = pd.read_csv(f)
        if 'esm_log_probs' in df.columns and 'NT_index' in df.columns:
            collector.append(df[['NT_index', 'esm_log_probs']].dropna())

    if not collector:
        print("No data found. Ensure the SNV pipeline has run.")
        return

    full_df = pd.concat(collector)

    # 2. Calculate position in codon (1, 2, or 3)
    full_df['Pos_in_Codon'] = (full_df['NT_index'] % 3) + 1

    # 3. Statistical Testing
    print("\n--- Statistical Analysis ---")

    # Split data into groups for testing
    group1 = full_df[full_df['Pos_in_Codon'] == 1]['esm_log_probs']
    group2 = full_df[full_df['Pos_in_Codon'] == 2]['esm_log_probs']
    group3 = full_df[full_df['Pos_in_Codon'] == 3]['esm_log_probs']

    # One-way ANOVA
    f_stat, p_anova = stats.f_oneway(group1, group2, group3)
    print(f"One-way ANOVA: F={f_stat:.2f}, p={p_anova:.2e}")

    # Tukey HSD Post-hoc (includes correction for multiple comparisons)
    tukey = pairwise_tukeyhsd(endog=full_df['esm_log_probs'],
                              groups=full_df['Pos_in_Codon'],
                              alpha=0.05)
    print("\nTukey HSD Pairwise Comparisons (Corrected P-values):")
    print(tukey)

    # 4. Group means and errors for plotting
    stats_summary = full_df.groupby('Pos_in_Codon')['esm_log_probs'].agg(['mean', 'sem'])
    print("\nSummary Stats:")
    print(stats_summary)

    # 5. Plotting
    plt.figure(figsize=(10, 6))

    # Use the first 3 colors from custom palette
    colors = COLOR_MAP.values().tolist()[:3]

    # Plot bars with error bars representing the Standard Error of the Mean (SEM)
    bars = plt.bar(stats_summary.index, stats_summary['mean'],
                   yerr=stats_summary['sem'],
                   color=colors, capsize=10, edgecolor='black', alpha=0.7)

    plt.title("ESM Pathogenicity Score by Codon Position\n(Biological Sanity Check)", fontsize=14)
    plt.xlabel("Nucleotide Position in Codon", fontsize=12)
    plt.ylabel("Mean ESM Log-Prob (Higher = More Benign)", fontsize=12)
    plt.xticks([1, 2, 3])
    plt.grid(axis='y', linestyle='--', alpha=0.5)

    # Add p-value annotation to the plot
    plt.annotate(f"ANOVA p: {p_anova:.2e}\nPost-hoc: all pairs p < 0.001",
                 xy=(0.05, 0.05), xycoords='axes fraction',
                 bbox=dict(boxstyle="round", fc="white", alpha=0.8))

    plt.tight_layout()
    plt.savefig("p4_plot_mean_3_codon_position_stats.png", dpi=300)
    print("\nPlot saved as p4_plot_mean_3_codon_position_stats.png")
    plt.show()


if __name__ == "__main__":
    validate_codon_positions()
