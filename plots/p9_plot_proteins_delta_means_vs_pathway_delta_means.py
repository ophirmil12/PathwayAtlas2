from plot_boot import *
boot_plot_folder()

import os
import pickle
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from definitions import *
from distance_utils import *
from kegg_api import KeggNetwork

# Configuration
cancer_type = 'pan_cancer'


def get_gene_color(q_val):
    """Returns color based on significance thresholds."""
    # Force float conversion just in case pandas read it as a string
    try:
        q_val = float(q_val)
    except (ValueError, TypeError):
        print("???")
        return "#a8a8a8" # Grey fallback if data is weird

    if pd.isna(q_val):
        return "#a8a8a8"
    elif q_val > 0.05:
        return COLOR_MAP['non-significant']  # Orange
    else:
        return COLOR_MAP['significant']  # Purple


def plot_pathway_gene_distributions():
    # 1. Setup Directories
    output_dir = pjoin(PLOTS_P, f'p9_pathway_protein_delta_means_{cancer_type}')
    os.makedirs(output_dir, exist_ok=True)

    # 2. Load Results
    pathway_results_file = pjoin(RESULTS_DISTANCES_P, f'{cancer_type}.csv')
    # Unfiltered gene results (all genes, including low-coverage ones)
    gene_results_unfiltered_file = pjoin(f"{RESULTS_GENE_LEVEL_P}_unfiltered", f'{cancer_type}_gene_distances.csv')
    # Filtered gene results (only genes that passed coverage threshold)
    gene_results_filtered_file = pjoin(RESULTS_GENE_LEVEL_P, f'{cancer_type}_gene_distances.csv')

    if not os.path.exists(pathway_results_file) or not os.path.exists(gene_results_unfiltered_file):
        print(f"Error: Results missing for {cancer_type}. Ensure p7 scripts ran successfully.")
        return

    pathway_df = pd.read_csv(pathway_results_file)
    gene_df_unfiltered = pd.read_csv(gene_results_unfiltered_file).set_index('gene_id')

    # Load filtered dataframe to keep both IDs and q-values
    filtered_df = pd.read_csv(gene_results_filtered_file)
    filtered_df.columns = filtered_df.columns.str.strip()
    if 'gene_id' in filtered_df.columns:
        filtered_df = filtered_df.set_index('gene_id')

    # GUARANTEE q_values are floats, not strings or objects
    if 'q_value' in filtered_df.columns:
        filtered_df['q_value'] = pd.to_numeric(filtered_df['q_value'], errors='coerce')
    else:
        print(f"CRITICAL ERROR: 'q_value' column missing! Columns found: {filtered_df.columns.tolist()}")

    filtered_gene_ids = set(filtered_df.index)

    print(f"Loaded {len(pathway_df)} pathways.")
    print(f"Loaded {len(gene_df_unfiltered)} genes (unfiltered), "
          f"{len(filtered_gene_ids)} passed min-mutations filter.")

    # Load Metadata for titles
    pathway_metadata = {}
    if os.path.exists(KEGG_PATHWAY_METADATA_FILE):
        with open(KEGG_PATHWAY_METADATA_FILE, 'rb') as f:
            pathway_metadata = pickle.load(f)

    set_paper_palette()

    # 3. Iterate through pathways
    for _, row in pathway_df[::-1].iterrows():
        p_id = row['pathway']
        p_delta = row['delta_means']
        p_q = row['q_value']

        # Get genes belonging to this pathway
        if p_id not in pathway_metadata:
            print(f"No metadata for {p_id}")
            continue
        meta = pathway_metadata[p_id]
        pathway_genes = meta.get('genes_ids', [])  # in 'hsa:XXXX' format

        # Split genes into: sufficient data vs insufficient data
        all_pathway_genes = gene_df_unfiltered[gene_df_unfiltered.index.isin(pathway_genes)].copy()

        if all_pathway_genes.empty:
            print(f"current_pathway_genes is empty for {p_id}")
            continue

        sufficient_genes = all_pathway_genes[all_pathway_genes.index.isin(filtered_gene_ids)].copy()
        insufficient_genes = all_pathway_genes[~all_pathway_genes.index.isin(filtered_gene_ids)].copy()

        # 4. Prepare Plotting Data for sufficient-data genes
        sufficient_genes['dot_size'] = np.log2(2 + sufficient_genes['num_mutations']) * 40

        # Pull the q-values from the FILTERED dataframe, matching by index (gene_id)
        sufficient_genes['q_value'] = sufficient_genes.index.map(filtered_df['q_value'])

        if not sufficient_genes.empty:
            n_sig = (sufficient_genes['q_value'] < 0.05).sum()
            print(f"\nPathway {p_id}: {len(sufficient_genes)} genes plotted. {n_sig} have q < 0.05")

            if n_sig == 0:
                print("   -> No significant genes in THIS pathway. Here are the 3 lowest q-values it found:")
                # Print the lowest 3 to see what the script is actually looking at
                lowest_genes = sufficient_genes['q_value'].nsmallest(3)
                print(lowest_genes)
            else:
                print("   -> Significant genes FOUND! They should be purple/blue.")

        sufficient_genes['dot_color'] = sufficient_genes['q_value'].apply(get_gene_color)

        # Insufficient-data genes: fixed small size (num_mutations may be very low/zero)
        insufficient_genes['dot_size'] = np.log2(2 + insufficient_genes.get(
            'num_mutations', pd.Series(0, index=insufficient_genes.index))) * 40

        # 5. Plotting
        plt.figure(figsize=(9, 7))

        # Violin is drawn only over the sufficient-data genes (the "real" distribution)
        if not sufficient_genes.empty:
            sns.violinplot(
                y=sufficient_genes['delta_means'],
                color=COLOR_MAP['light-green'],
                inner=None,
                alpha=0.2
            )

        # --- Sufficient-data genes (filled circles) ---
        if not sufficient_genes.empty:
            jitter = np.random.uniform(-0.15, 0.15, size=len(sufficient_genes))
            plt.scatter(
                x=jitter,
                y=sufficient_genes['delta_means'],
                s=sufficient_genes['dot_size'],
                c=sufficient_genes['dot_color'],
                alpha=0.8,
                edgecolors='white',
                linewidths=0.5,
                zorder=5,
            )

        # --- Insufficient-data genes (hollow circles, faded) ---
        if not insufficient_genes.empty:
            jitter_insuf = np.random.uniform(-0.15, 0.15, size=len(insufficient_genes))
            plt.scatter(
                x=jitter_insuf,
                y=insufficient_genes['delta_means'],
                s=insufficient_genes['dot_size'],
                facecolors='none',           # hollow
                edgecolors=COLOR_MAP['grey'],
                linewidths=1.0,
                alpha=0.6,
                zorder=4,                    # behind sufficient-data dots
                marker='o',
            )

        # --- Pathway star ---
        plt.scatter(
            [0], [p_delta],
            color=COLOR_MAP['red'],
            s=400,
            marker='*',
            edgecolors='black',
            zorder=10,
            label=f'Pathway $\Delta\mu$ (q={p_q:.2e})' if not pd.isna(p_q) else 'Pathway $\Delta\mu$ (q=N/A)'
        )

        # Formatting
        display_name = p_id
        if p_id in pathway_metadata:
            meta = pathway_metadata[p_id]
            display_name = meta.get('name', p_id) if isinstance(meta, dict) else getattr(meta, 'name', p_id)

        plt.axhline(0, color='black', linestyle='--', alpha=0.3)
        plt.title(
            f"{display_name}\nGene-level Significance Distribution ({cancer_type})"
            f"\n({len(sufficient_genes)} genes, {len(insufficient_genes)} insufficient data)",
            fontsize=13
        )
        plt.ylabel("Delta Mean ($\Delta\mu$)", fontsize=12)
        plt.xticks([])

        # Custom Legend
        from matplotlib.lines import Line2D
        legend_elements = [
            Line2D([0], [0], marker='o', color='w', label='q > 0.05',
                   markerfacecolor=COLOR_MAP['non-significant'], markersize=10),
            Line2D([0], [0], marker='o', color='w', label='q < 0.05',
                   markerfacecolor=COLOR_MAP['light-blue'], markersize=10),
            Line2D([0], [0], marker='o', color='w', label='q < 0.01',
                   markerfacecolor=COLOR_MAP['significant'], markersize=10),
            Line2D([0], [0], marker='o', color='w', label='Insufficient data',
                   markerfacecolor='none', markeredgecolor=COLOR_MAP['grey'],
                   markeredgewidth=1.0, markersize=8),
            Line2D([0], [0], marker='*', color='w', label='Pathway Mean',
                   markerfacecolor=COLOR_MAP['red'], markersize=15, markeredgecolor='black'),
        ]
        plt.legend(handles=legend_elements, loc='best', title="Significance")

        plt.tight_layout()

        safe_name = p_id.replace(":", "_")
        plt.savefig(pjoin(output_dir, f"{safe_name}.png"), dpi=300)
        plt.close()


if __name__ == '__main__':
    plot_pathway_gene_distributions()
    print(f"Finished generating plots for {cancer_type}.")