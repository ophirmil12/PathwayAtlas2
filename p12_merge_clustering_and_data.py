# For each cluster created by the R script, for each cancer type, take all significant pathways in it,
# and take the median of the delta-means, and create the enrichment matrix [cancers*clusters]

import pandas as pd
import os
from definitions import RESULTS_DISTANCES_P, RESULTS_P

def merge_clusters_to_cancers(db_type="pathway"):
    # --- 1. Load the dynamic R-generated mappings ---
    summary_path = os.path.join(RESULTS_P, f"p12_{db_type}_clusters_summary.csv")
    annotations_path = os.path.join(RESULTS_P, f"p12_{db_type}_cluster_annotations.csv")

    if not os.path.exists(summary_path):
        print(f"Skipping {db_type.upper()}: {summary_path} not found.")
        return

    df_summary = pd.read_csv(summary_path)
    df_annotations = pd.read_csv(annotations_path)

    # --- 2. Load the 38 cancer CSVs ---
    print(f"Loading raw cancer distance matrices for {db_type.upper()}...")
    all_data = []
    
    for file in os.listdir(RESULTS_DISTANCES_P):
        if file.endswith(".csv"):
            cancer_type = file.replace(".csv", "")
            df_cancer = pd.read_csv(os.path.join(RESULTS_DISTANCES_P, file))
            
            df_sig = df_cancer[df_cancer['q_value'] < 0.05].copy()
            df_sig['cancer'] = cancer_type
            
            # The column in your CSVs is always named 'pathway', even if it holds module IDs
            all_data.append(df_sig[['pathway', 'delta_means', 'cancer']])

    df_cancers = pd.concat(all_data)

    # --- 3. Map to Clusters ---
    # The new R script outputs 'id' instead of 'pathway_id'
    df_merged = df_cancers.merge(df_summary, left_on='pathway', right_on='id', how='inner')

    # --- 4. Aggregate & Pivot ---
    # Merge summary (which now has short_name) with cancer data
    df_merged = df_cancers.merge(df_summary, left_on='pathway', right_on='id', how='inner')
    df_grouped = df_merged.groupby(['cancer', 'cluster', 'short_name'])['delta_means'].median().reset_index()

    df_pivot = df_grouped.pivot(index=['cluster', 'short_name'], columns='cancer', values='delta_means').reset_index()

    # --- 5. Attach Contained List (from annotations) ---
    df_final = df_pivot.merge(df_annotations[['cluster', 'count', 'contained']], on='cluster', how='left')

    # Reorder columns: Cluster Info first, then Cancers
    cancer_cols = [c for c in df_final.columns if c not in ['cluster', 'short_name', 'count', 'contained']]
    df_final = df_final[['cluster', 'short_name', 'count', 'contained'] + cancer_cols]

    # --- 6. Export ---
    output_path = os.path.join(RESULTS_P, f"p12_{db_type}_cancer_cluster_enrichment_matrix.csv")
    df_final.to_csv(output_path, index=False)
    print(f"Saved {db_type.upper()} matrix to: {output_path}\n")

if __name__ == "__main__":
    for db in ["pathway", "module"]:
        merge_clusters_to_cancers(db)