# Merge gene snvs tables by pathway

import os
from definitions import *
import pickle
import pandas as pd

def merge_genes_to_pathway():
    if os.path.exists(KEGG_PATHWAY_METADATA_P):
        with open(KEGG_PATHWAY_METADATA_P, 'rb') as f:
            pathway_metadata = pickle.load(f)
    else:
        print(f"    ERROR: KEGG pathway metadata file not found at {KEGG_PATHWAY_METADATA_P}")
        return

    print("----- Merging gene SNVs to pathway SNVs... -----")

    for pathway_id, metadata in pathway_metadata.items():

        print(f"    Merging gene SNVs for pathway: {pathway_id}")
        merged_df = pd.DataFrame()
        output_path = os.path.join(KEGG_PATHWAY_SCORES_P, f"{pathway_id}.csv")

        for gene_id in metadata['genes_ids']:
            gene_snv_file = os.path.join(KEGG_GENE_SCORES_P, f"{gene_id}.csv")

            if os.path.exists(gene_snv_file):
                gene_df = pd.read_csv(gene_snv_file)
                merged_df = pd.concat([merged_df, gene_df], ignore_index=True)
                merged_df.to_csv(output_path, index=False)
            else:
                print(f"    SNV file for gene {gene_id} not found, skipping.")

        print(f"    Saved merged gene SNVs for pathway {pathway_id} to {output_path}")


if __name__ == '__main__':
    # Merge gene SNVs to pathway SNVs
    merge_genes_to_pathway()