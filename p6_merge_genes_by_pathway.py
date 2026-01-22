# Merge gene snvs tables by pathway


import os
import pickle
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
from definitions import KEGG_PATHWAY_METADATA_FILE, KEGG_PATHWAY_SCORES_P, KEGG_GENE_SCORES_P


def merge_single_pathway(pathway_id, gene_ids):
    """
    Worker function to merge all genes for a single pathway.
    """
    output_path = os.path.join(KEGG_PATHWAY_SCORES_P, f"{pathway_id}.csv")

    # Optional: Skip if already exists to allow resuming
    # if os.path.exists(output_path):
    #    return f"Skipped {pathway_id}"

    df_list = []

    for gene_id in gene_ids:
        g_id = gene_id.replace(":", "_")
        gene_snv_file = os.path.join(KEGG_GENE_SCORES_P, f"{g_id}.csv")

        if os.path.exists(gene_snv_file):
            try:
                gene_df = pd.read_csv(gene_snv_file)
                if not gene_df.empty:
                    df_list.append(gene_df)
            except Exception as e:
                print(f"Error reading {gene_snv_file}: {e}")

    if not df_list:
        return f"No data for {pathway_id}"

    # Concatenate all genes at once
    pathway_df = pd.concat(df_list, ignore_index=True)

    # Write to disk exactly once
    pathway_df.to_csv(output_path, index=False)
    return f"Merged {pathway_id}"


def merge_genes_to_pathway():
    # 1. Load Metadata
    if not os.path.exists(KEGG_PATHWAY_METADATA_FILE):
        print(f"ERROR: Metadata file not found at {KEGG_PATHWAY_METADATA_FILE}")
        return

    with open(KEGG_PATHWAY_METADATA_FILE, 'rb') as f:
        pathway_metadata = pickle.load(f)

    os.makedirs(KEGG_PATHWAY_SCORES_P, exist_ok=True)

    pathways_to_process = list(pathway_metadata.items())
    print(f"----- Merging SNVs for {len(pathways_to_process)} pathways... -----")

    # 2. Run in Parallel
    # Using ProcessPoolExecutor to bypass Python's GIL for CSV parsing
    with ProcessPoolExecutor(max_workers=os.cpu_count()) as executor:
        futures = {
            executor.submit(merge_single_pathway, pid, meta['genes_ids']): pid
            for pid, meta in pathways_to_process
        }

        for future in tqdm(as_completed(futures), total=len(futures), desc="Merging Pathways"):
            res = future.result()
            # print(res) # Optional: uncomment for verbose logging


if __name__ == '__main__':
    merge_genes_to_pathway()