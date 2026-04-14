import os
import glob
import re
import pandas as pd
import numpy as np
from collections import Counter
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, fcluster
import matplotlib.pyplot as plt
import seaborn as sns

# Import the class from your provided file
from kegg_api import KeggApi
from definitions import RESULTS_DISTANCES_P, RESULTS_P, KEGG_PATHWAY_METADATA_FILE

# Configuration
CUTOFF = 0.85  # Distance threshold for clustering (1 - Jaccard)


def clean_biological_name(name):
    """Strips biological noise words (mimics gseapy label cleaning)."""
    name = name.split(' - Homo sapiens')[0]
    # Regex to remove redundant terms
    noise = [
        r'\bpathway\b', r'\bsignaling\b', r'\bsystem\b', r'\bmetabolism\b',
        r'\bhuman\b', r'\bprotein\b', r'\bmediated\b', r'\bregulation\b',
        r'\borganismal\b', r'\bcellular\b'
    ]
    short_name = name.lower()
    for pattern in noise:
        short_name = re.sub(pattern, '', short_name)
    # Remove non-alphanumeric, collapse spaces, and Title Case
    short_name = re.sub(r'[^a-zA-Z0-9\s]', ' ', short_name)
    return " ".join(short_name.split()[:4]).title()


def get_cluster_short_name(cluster_ids, ids_list, sim_matrix, names_map):
    """Finds the 'Medoid' (most central pathway) of the cluster to use as the name."""
    if len(cluster_ids) == 1:
        return clean_biological_name(names_map[cluster_ids[0]])

    # Get indices of cluster members in the global similarity matrix
    indices = [ids_list.index(cid) for cid in cluster_ids]
    sub_sim = sim_matrix[np.ix_(indices, indices)]

    # Medoid = member with highest average similarity to all other members
    medoid_idx_in_sub = np.argmax(sub_sim.mean(axis=1))
    medoid_id = cluster_ids[medoid_idx_in_sub]

    return clean_biological_name(names_map[medoid_id])

def generate_smart_name(names_list):
    """
    Summarizes a cluster by identifying the most frequent meaningful keywords.
    Mimics gseapy's approach to redundant term naming.
    """
    stop_words = {
        'and', 'the', 'of', 'in', 'to', 'with', 'by', 'pathway', 'signaling',
        'human', 'metabolism', 'hsa', 'organism', 'pathways', 'system'
    }
    words = []
    for name in names_list:
        # Remove punctuation and split
        tokens = re.findall(r'\b\w{3,}\b', name.lower())
        words.extend([w for w in tokens if w not in stop_words])

    if not words:
        return ""

    # Pick the top 3 most common keywords
    most_common = [word for word, count in Counter(words).most_common(3)]
    return " ".join(most_common).title()


def run_clustering(db_type):
    """
    db_type: "pathway" or "module"
    """
    print(f"\n{'=' * 40}\nSTARTING PIPELINE FOR: {db_type.upper()}\n{'=' * 40}")
    api = KeggApi()

    # 1. Load significant pathways/modules from distances CSVs
    files = glob.glob(os.path.join(RESULTS_DISTANCES_P, "*.csv"))
    if not files:
        print(f"No files found in {RESULTS_DISTANCES_P}")
        return

    all_dfs = []
    for f in files:
        cancer_name = os.path.basename(f).replace(".csv", "")
        df = pd.read_csv(f)
        # Filter for significance
        df = df[df['q_value'] < 0.05][['pathway', 'delta_means']]
        df['cancer'] = cancer_name
        all_dfs.append(df)

    if not all_dfs:
        print("No significant pathways found across files.")
        return

    df_main = pd.concat(all_dfs)

    # Pivot to matrix (Pathways x Cancers)
    # Using median (as in R script) to aggregate duplicate pathways across cancers
    # TODO: it might be worth changing median() to max(): in here just once, regenerating the heatmaps,
    #  and seeing if those "pathogenic clusters" visually pop out much clearer than they do with the median
    matrix_df = df_main.pivot_table(index='pathway', columns='cancer',
                                    values='delta_means', aggfunc='median').fillna(0)

    pathway_ids_in_data = matrix_df.index.tolist()

    # 2. Load Metadata from Pickle (Replacing API)
    print(f"Loading {db_type} metadata from local pickle...")
    full_metadata = pd.read_pickle(KEGG_PATHWAY_METADATA_FILE)

    gene_sets = {}
    names_map = {}

    for kid, info in full_metadata.items():
        if kid in pathway_ids_in_data:
            is_module = bool(re.match(r'^M\d+', kid))
            if (db_type == "module" and is_module) or (db_type == "pathway" and not is_module):
                gene_sets[kid] = set(info['genes_ids'])
                names_map[kid] = info['name']

    if not gene_sets:
        print(f"No matching {db_type}s found in metadata for the significant IDs.")
        return

    # 3. Calculate Jaccard Similarity Matrix
    ids = list(gene_sets.keys())
    n = len(ids)
    sim_matrix = np.zeros((n, n))

    print(f"Calculating Jaccard similarity for {n} entries...")
    for i in range(n):
        sim_matrix[i, i] = 1.0
        for j in range(i + 1, n):
            set_i = gene_sets[ids[i]]
            set_j = gene_sets[ids[j]]
            intersection = len(set_i.intersection(set_j))
            union = len(set_i.union(set_j))
            sim = intersection / union if union > 0 else 0
            sim_matrix[i, j] = sim_matrix[j, i] = sim

    # 4. Hierarchical Clustering
    # Distance = 1 - Similarity
    dist_matrix = 1 - sim_matrix
    # Complete linkage is robust for finding compact clusters
    Z = linkage(squareform(dist_matrix), method='complete')

    # fcluster uses the CUTOFF to define flat clusters
    clusters = fcluster(Z, t=CUTOFF, criterion='distance')

    # 5. Create results DataFrames
    df_clusters = pd.DataFrame({'id': ids, 'name': [names_map[i] for i in ids], 'cluster': clusters})

    annotations = []
    for cl_id, group in df_clusters.groupby('cluster'):
        cl_ids = group['id'].tolist()
        # Use Medoid Logic here
        short_name = get_cluster_short_name(cl_ids, ids, sim_matrix, names_map)

        annotations.append({
            'cluster': cl_id,
            'short_name': short_name,
            'count': len(group),
            'contained': " | ".join(group['name'].tolist())
        })

    df_annotations = pd.DataFrame(annotations).sort_values('count', ascending=False)
    # Merge short_name back into the main summary file for use in heatmaps
    df_clusters = df_clusters.merge(df_annotations[['cluster', 'short_name']], on='cluster')

    # 6. Save Results (Maintaining R script filenames)
    df_clusters.to_csv(os.path.join(RESULTS_P, f"p12_{db_type}_clusters_summary.csv"), index=False)
    df_annotations.to_csv(os.path.join(RESULTS_P, f"p12_{db_type}_cluster_annotations.csv"), index=False)

    # 7. Heatmap Visualization
    plt.figure(figsize=(16, 12))
    # Sort matrix by cluster ID for visual grouping
    order = df_clusters.sort_values('cluster').index
    sorted_sim = sim_matrix[np.ix_(order, order)]
    sorted_labels = df_clusters.iloc[order]['name'].values

    sns.heatmap(sorted_sim,
                xticklabels=False,
                yticklabels=sorted_labels,
                cmap="YlGnBu",
                cbar_kws={'label': 'Jaccard Similarity'})

    plt.title(f"KEGG Functional Clusters ({db_type.upper()}) - Similarity Heatmap")
    plt.tight_layout()

    # Maintaining R script filename
    plt.savefig(os.path.join(RESULTS_P, f"p12_simplified_kegg_{db_type}_heatmap.pdf"))
    plt.close()

    print(f"Finished {db_type} pipeline. Results saved to {RESULTS_P}")


if __name__ == "__main__":
    # Run for both pathways and modules
    for target in ["pathway", "module"]:
        run_clustering(target)