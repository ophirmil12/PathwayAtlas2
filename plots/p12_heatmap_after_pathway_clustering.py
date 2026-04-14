# After clustering the pathways anf merging, plotting the cancer similarity matrix
# and the cancers*pathways merged matrix

from plot_boot import *
boot_plot_folder()  # coloring scheme using cycler
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from scipy.cluster.hierarchy import linkage, leaves_list, dendrogram
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

# Import your project definitions
from definitions import RESULTS_P, PLOTS_P

# --- Setup ---
boot_plot_folder()
output_folder = PLOTS_P
os.makedirs(output_folder, exist_ok=True)


def hierarchical_order(df, axis=0, method='ward', metric='euclidean'):
    """Return reordered index/columns using Ward hierarchical clustering."""
    data = df.fillna(0).values if axis == 0 else df.fillna(0).T.values
    if data.shape[0] < 2:
        return df.index if axis == 0 else df.columns
    Z = linkage(data, method=method, metric=metric)
    order = leaves_list(Z)
    return (df.index if axis == 0 else df.columns)[order]


def plot_clusters(db_type="pathway"):
    print(f"\nProcessing Plots for {db_type.upper()} Clusters...")

    # 1. Load Data
    file_path = os.path.join(RESULTS_P, f"p12_{db_type}_cancer_cluster_enrichment_matrix.csv")
    if not os.path.exists(file_path):
        print(f"File not found: {file_path}")
        return

    df = pd.read_csv(file_path)

    # Create descriptive label for the clusters
    df['label'] = "C" + df['cluster'].astype(str) + ": " + df['short_name'] + " (" + df['count'].astype(str) + ")"

    metadata_cols = ['cluster', 'short_name', 'count', 'contained', 'label']
    cancer_cols = [c for c in df.columns if c not in metadata_cols]

    # Orientation: Cancers as Rows, Clusters as Columns
    # We transpose the input so cancers are index, clusters are columns
    delta_filt = df.set_index('label')[cancer_cols].T

    # --- Plot A: Clustered Heatmap (Normal Orientation) ---
    print(f"Plotting A: {db_type} Cluster Heatmap...")
    row_order = hierarchical_order(delta_filt, axis=0)  # Cancers
    col_order = hierarchical_order(delta_filt, axis=1)  # Clusters

    delta_ordered = delta_filt.loc[row_order, col_order]
    delta_display = delta_ordered.fillna(0)

    # Use 98th percentile for robust color scaling
    vmax = np.nanpercentile(np.abs(delta_filt.values), 98)
    vmin = -vmax

    # Dynamic scaling for "Normal" orientation:
    # Height scales with number of Cancers, Width scales with number of Clusters
    fig_width = max(12, delta_display.shape[1] * 0.18)
    fig_height = max(8, delta_display.shape[0] * 0.35)

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    # Plotting the "Normal way" (No .T on delta_display)
    im = ax.imshow(delta_display.values, aspect='auto', cmap='RdBu_r',
                   vmin=vmin, vmax=vmax, interpolation='nearest')

    # Formatting
    ax.set_yticks(range(delta_display.shape[0]))
    ax.set_yticklabels(delta_display.index, fontsize=9)

    ax.set_xticks(range(delta_display.shape[1]))
    ax.set_xticklabels(delta_display.columns, rotation=90, fontsize=6)

    plt.colorbar(im, ax=ax, fraction=0.02, pad=0.01, label='Median Delta Means')
    ax.set_title(
        f'Functional Cluster Enrichment: {db_type.upper()}\n({delta_display.shape[1]} Clusters, {delta_display.shape[0]} Cancers)',
        fontsize=16)

    plt.tight_layout()
    plt.savefig(os.path.join(output_folder, f'p12_{db_type}_cluster_heatmap.png'), dpi=600, bbox_inches='tight')
    plt.close()

    # --- Plot B: Cancer-Cancer Similarity with Left Dendrogram ---
    print(f"Plotting B: {db_type} Cancer Similarity + Dendrogram...")
    delta_filled = delta_filt.fillna(0)

    # 1. Calculate Similarity
    denom = np.linalg.norm(delta_filled.values, axis=1, keepdims=True) + 1e-12
    normed = delta_filled.values / denom
    sim_mat = pd.DataFrame(normed @ normed.T, index=delta_filled.index, columns=delta_filled.index)

    # 2. Calculate Linkage
    Z = linkage(sim_mat, method='ward', metric='euclidean')

    # 3. Get the order from clustering
    order = leaves_list(Z)
    sim_ordered = sim_mat.iloc[order, order]

    # 4. Create Figure with two columns: [Dendrogram | Heatmap]
    # width_ratios=[1, 10] means the dendrogram takes 10% of the width
    fig, (ax_dendro, ax_matrix) = plt.subplots(1, 2, figsize=(12, 10),
                                               gridspec_kw={'width_ratios': [0.15, 1], 'wspace': 0.2})

    # A. Plot the Dendrogram on the left axis
    with plt.rc_context({'lines.linewidth': 1.5}):
        dendro = dendrogram(Z, orientation='left', ax=ax_dendro,
                            color_threshold=0, above_threshold_color='black',
                            no_labels=True)  # Hide the numeric IDs

    ax_dendro.set_axis_off()  # Remove the box and background from the dendrogram

    # B. Plot the Matrix on the right axis
    im = ax_matrix.imshow(sim_ordered.values, aspect='auto', cmap='RdBu_r', vmin=-1, vmax=1)

    # C. Formatting labels (Attached to the matrix)
    ax_matrix.set_yticks(range(len(sim_ordered)))
    ax_matrix.set_yticklabels(sim_ordered.index, fontsize=9)

    ax_matrix.set_xticks(range(len(sim_ordered)))
    ax_matrix.set_xticklabels(sim_ordered.columns, rotation=90, fontsize=9)

    # D. Colorbar and Title
    plt.colorbar(im, ax=ax_matrix, fraction=0.046, pad=0.04, label='Cosine Similarity')
    ax_matrix.set_title(f'Cancer Similarity: {db_type.upper()} Clusters', fontsize=14, pad=20)

    plt.savefig(os.path.join(output_folder, f'p12_{db_type}_cancer_similarity.png'), dpi=600, bbox_inches='tight')
    plt.close()

    # --- Plot D: PCA ---
    print(f"Plotting D: {db_type} PCA...")
    X = delta_filt.fillna(0)
    X_scaled = StandardScaler().fit_transform(X.values)

    pca = PCA(n_components=min(10, X_scaled.shape[0] - 1))
    coords = pca.fit_transform(X_scaled)
    var_exp = pca.explained_variance_ratio_ * 100

    plt.figure(figsize=(12, 9))
    # Using the project cycler via color palette
    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = prop_cycle.by_key()['color']

    for i, cancer in enumerate(X.index):
        plt.scatter(coords[i, 0], coords[i, 1], s=120, edgecolors='black',
                    alpha=0.8, color=colors[i % len(colors)], label=None)
        plt.annotate(cancer, (coords[i, 0], coords[i, 1]), fontsize=9,
                     xytext=(0, 6), textcoords='offset points', ha='center', weight='bold')

    plt.axhline(0, color='grey', lw=1, ls='--')
    plt.axvline(0, color='grey', lw=1, ls='--')
    plt.xlabel(f'PC1 ({var_exp[0]:.1f}%)', fontsize=12)
    plt.ylabel(f'PC2 ({var_exp[1]:.1f}%)', fontsize=12)
    plt.title(f'PCA: Cancers in {db_type.upper()} Functional Space', fontsize=14)
    sns.despine()
    plt.savefig(os.path.join(output_folder, f'p12_{db_type}_pca.png'), dpi=600, bbox_inches='tight')
    plt.close()


if __name__ == "__main__":
    for db in ["pathway", "module"]:
        plot_clusters(db)
    print(f"\nDone. All p12 cluster plots saved to: {output_folder}")