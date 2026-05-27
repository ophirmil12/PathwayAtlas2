"""
p12_pathway_clustering.py
Cluster KEGG pathways by gene-set similarity using Jaccard distance.

Method:
  1. Load pathway -> gene sets from KEGG_PATHWAY_METADATA_FILE.
     Each entry has a 'genes_ids' list — that is the gene set.
  2. Build a binary (one-hot) matrix: rows = pathways, cols = all unique genes.
  3. Compute pairwise Jaccard similarity between every pair of pathways.
  4. Cluster with average-linkage hierarchical clustering, cutting at
     Jaccard DISTANCE threshold = 1 - PATHWAY_JACCARD_SIMILARITY_THRESHOLD.
     (i.e. two pathways merge if Jaccard similarity >= 0.80)
  5. Save results to KEGG_PATHWAY_CLUSTERING_P.
"""

import os
import pickle
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, fcluster
from definitions import (
    KEGG_PATHWAY_METADATA_FILE,
    KEGG_PATHWAY_CLUSTERING_P,
    PATHWAY_JACCARD_SIMILARITY_THRESHOLD,
)


# ── 1. Load pathway gene sets from metadata ───────────────────────────────────

def load_pathway_gene_sets() -> dict[str, set]:
    """
    Return {pathway_id: set_of_gene_ids} from KEGG_PATHWAY_METADATA_FILE.
    Each metadata entry is expected to have a 'genes_ids' key.
    Pathways with no genes are skipped.
    """
    with open(KEGG_PATHWAY_METADATA_FILE, 'rb') as f:
        metadata: dict = pickle.load(f)

    pathway_gene_sets = {}
    for pathway_id, entry in metadata.items():
        genes = entry.get('genes_ids', [])
        if genes:
            pathway_gene_sets[pathway_id] = set(genes)

    return pathway_gene_sets


# ── 2. Build one-hot matrix ───────────────────────────────────────────────────

def build_onehot_matrix(pathway_gene_sets: dict[str, set]) -> tuple[pd.DataFrame, list]:
    """
    Returns:
        onehot_df   - DataFrame of shape (n_pathways, n_genes), dtype uint8
        pathway_ids - ordered list matching row indices
    """
    pathway_ids = sorted(pathway_gene_sets.keys())
    all_genes   = sorted({g for genes in pathway_gene_sets.values() for g in genes})

    # remove all modules (M...)
    pathway_ids = [pid for pid in pathway_ids if not pid.startswith('M')]

    gene_index = {g: i for i, g in enumerate(all_genes)}
    matrix = np.zeros((len(pathway_ids), len(all_genes)), dtype=np.uint8)

    for row_idx, pid in enumerate(pathway_ids):
        for gene in pathway_gene_sets[pid]:
            matrix[row_idx, gene_index[gene]] = 1

    onehot_df = pd.DataFrame(matrix, index=pathway_ids, columns=all_genes)
    return onehot_df, pathway_ids


# ── 3. Jaccard distance + hierarchical clustering ────────────────────────────

def cluster_pathways(
    onehot_df: pd.DataFrame,
    similarity_threshold: float,
) -> tuple[pd.DataFrame, np.ndarray]:
    """
    Compute pairwise Jaccard distances, run average-linkage clustering,
    cut at distance = 1 - similarity_threshold.

    Returns (result_df, condensed_dist_array).
    """
    distance_threshold = 1.0 - similarity_threshold
    pathway_ids = list(onehot_df.index)

    condensed_dist = pdist(onehot_df.values, metric='jaccard')
    Z = linkage(condensed_dist, method='average')
    cluster_labels = fcluster(Z, t=distance_threshold, criterion='distance')

    result_df = pd.DataFrame({
        'pathway_id': pathway_ids,
        'cluster_id': cluster_labels,
    }).sort_values('cluster_id').reset_index(drop=True)

    return result_df, condensed_dist


# ── 4. Annotate with pathway names ────────────────────────────────────────────

def annotate_with_names(result_df: pd.DataFrame) -> pd.DataFrame:
    """Adds a 'pathway_name' column from the metadata file."""
    with open(KEGG_PATHWAY_METADATA_FILE, 'rb') as f:
        metadata: dict = pickle.load(f)

    def get_name(pid):
        entry = metadata.get(pid, {})
        return entry.get('name', pid) if isinstance(entry, dict) else str(entry)

    result_df = result_df.copy()
    result_df.insert(1, 'pathway_name', result_df['pathway_id'].map(get_name))
    return result_df


# ── 5. Save results ───────────────────────────────────────────────────────────

def save_results(
    result_df: pd.DataFrame,
    onehot_df: pd.DataFrame,
    condensed_dist: np.ndarray,
    output_dir: str,
):
    os.makedirs(output_dir, exist_ok=True)
    pathway_ids = list(onehot_df.index)

    # Cluster assignments
    clusters_path = os.path.join(output_dir, 'pathway_clusters.csv')
    result_df.to_csv(clusters_path, index=False)
    print(f"[OK] Cluster assignments       -> {clusters_path}")

    # One-hot matrix
    onehot_path = os.path.join(output_dir, 'pathway_onehot_matrix.csv')
    onehot_df.to_csv(onehot_path)
    print(f"[OK] One-hot matrix            -> {onehot_path}")

    # Full Jaccard similarity matrix
    sim_matrix = 1.0 - squareform(condensed_dist)
    sim_df = pd.DataFrame(sim_matrix, index=pathway_ids, columns=pathway_ids)
    sim_path = os.path.join(output_dir, 'pathway_jaccard_similarity.csv')
    sim_df.to_csv(sim_path)
    print(f"[OK] Jaccard similarity matrix -> {sim_path}")

    # Cluster summary: size + member list per cluster
    summary = (
        result_df.groupby('cluster_id')['pathway_id']
        .agg(cluster_size='count', members=lambda x: '|'.join(sorted(x)))
        .reset_index()
        .sort_values('cluster_size', ascending=False)
    )
    summary_path = os.path.join(output_dir, 'cluster_summary.csv')
    summary.to_csv(summary_path, index=False)
    print(f"[OK] Cluster summary           -> {summary_path}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    print("=== Pathway Jaccard Clustering ===")
    print(f"Similarity threshold : {PATHWAY_JACCARD_SIMILARITY_THRESHOLD}")

    print("\n[1/4] Loading pathway gene sets from metadata ...")
    pathway_gene_sets = load_pathway_gene_sets()
    print(f"      Loaded {len(pathway_gene_sets)} pathways.")

    print("[2/4] Building one-hot matrix ...")
    onehot_df, pathway_ids = build_onehot_matrix(pathway_gene_sets)
    print(f"      Matrix shape: {onehot_df.shape}  (pathways x genes)")

    print("[3/4] Computing Jaccard distances and clustering ...")
    result_df, condensed_dist = cluster_pathways(onehot_df, PATHWAY_JACCARD_SIMILARITY_THRESHOLD)
    result_df = annotate_with_names(result_df)
    n_clusters = result_df['cluster_id'].nunique()
    print(f"      Found {n_clusters} clusters from {len(pathway_ids)} pathways.")

    print("[4/4] Saving results ...")
    save_results(result_df, onehot_df, condensed_dist, KEGG_PATHWAY_CLUSTERING_P)

    print("\n=== Done ===")


if __name__ == '__main__':
    main()