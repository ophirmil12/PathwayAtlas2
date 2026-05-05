from plot_boot import *
boot_plot_folder()

import glob
import pickle
import pandas as pd
import matplotlib.pyplot as plt
import os
from definitions import (
    KEGG_PATHWAY_METADATA_FILE,
    RESULTS_DISTANCES_P,
    PLOTS_P,
    COLOR_MAP,
    set_paper_palette,
)

OUTPUT_DIR = os.path.join(PLOTS_P, "p6_pathway_size_distributions")
os.makedirs(OUTPUT_DIR, exist_ok=True)


def get_significant_pathway_ids(distances_p: str, q_threshold: float = 0.05) -> set:
    """Collect pathway IDs that are significant (q < threshold) in any cancer CSV."""
    sig_ids = set()
    for fpath in glob.glob(os.path.join(distances_p, '*.csv')):
        try:
            df = pd.read_csv(fpath)
            sig = df[df['q_value'] < q_threshold]['pathway']
            sig_ids.update(sig.tolist())
        except Exception as e:
            print(f'[WARNING] Could not read {fpath}: {e}')
    return sig_ids


def plot_num_genes_hist(num_genes, title, out_path, color, xlim=None):
    fig, ax = plt.subplots(figsize=(7, 4))
    data = [n for n in num_genes if xlim is None or n <= xlim]
    ax.hist(data, bins=40, color=color, edgecolor='white', linewidth=0.5)
    ax.set_title(title, fontsize=13, fontweight='bold')
    ax.set_xlabel('Number of Genes', fontsize=11)
    ax.set_ylabel('Number of Pathways', fontsize=11)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    if xlim is not None:
        ax.set_xlim(0, xlim)
    median_val = sorted(data)[len(data) // 2]
    ax.axvline(median_val, color=COLOR_MAP['dark-red'], linestyle='--',
               linewidth=1.4, label=f'Median: {median_val}')
    ax.legend(fontsize=10)
    plt.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f'Saved: {out_path}')


def main():
    set_paper_palette()

    with open(KEGG_PATHWAY_METADATA_FILE, 'rb') as f:
        metadata = pickle.load(f)

    # --- All hsa pathways ---
    num_genes_all = [len(v['genes_ids']) for k, v in metadata.items() if k.startswith('hsa')]

    plot_num_genes_hist(
        num_genes_all,
        title='Distribution of Number of Genes per Pathway',
        out_path=os.path.join(OUTPUT_DIR, 'p6_distribution_num_genes.png'),
        color=COLOR_MAP['dark-blue'],
    )
    plot_num_genes_hist(
        num_genes_all,
        title='Distribution of Number of Genes per Pathway (0–250)',
        out_path=os.path.join(OUTPUT_DIR, 'p6_distribution_num_genes_0_250.png'),
        color=COLOR_MAP['dark-blue'],
        xlim=250,
    )

    # --- Significant pathways only (q < 0.05 in any cancer) ---
    sig_ids = get_significant_pathway_ids(RESULTS_DISTANCES_P)
    print(f'Significant pathways (any cancer, q<0.05): {len(sig_ids)}')

    num_genes_sig = [
        len(v['genes_ids']) for k, v in metadata.items()
        if k.startswith('hsa') and k in sig_ids
    ]

    plot_num_genes_hist(
        num_genes_sig,
        title='Distribution of Number of Genes per Pathway\n(Significant, q<0.05)',
        out_path=os.path.join(OUTPUT_DIR, 'p6_distribution_num_genes_significant.png'),
        color=COLOR_MAP['significant'],
    )
    plot_num_genes_hist(
        num_genes_sig,
        title='Distribution of Number of Genes per Pathway (0–250)\n(Significant, q<0.05)',
        out_path=os.path.join(OUTPUT_DIR, 'p6_distribution_num_genes_significant_0_250.png'),
        color=COLOR_MAP['significant'],
        xlim=250,
    )

    # print the top 20 pathways with the most genes
    top_20 = sorted(metadata.items(), key=lambda x: len(x[1]['genes_ids']), reverse=True)[:20]
    print("Top 20 Pathways with the Most Genes:")
    for pathway_id, pathway_info in top_20:
        print(f"{pathway_id}: {len(pathway_info['genes_ids'])} genes; Name: {pathway_info['name']}")


if __name__ == '__main__':
    main()