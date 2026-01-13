# Go through all genes, for each create V3 version, predict the disorder scores

import os
import pandas as pd
from tqdm import tqdm
import metapredict as meta
import numpy as np


from definitions import V3_version_letters, KEGG_GENE_SCORES_P, DISORDERED_THRESHOLD
from kegg_api import KeggApi, KeggGene


# 1. Init KeggAPI
api = KeggApi()

# 2. Get all human genes
print("(p3) Fetching master gene list from KEGG...")
all_genes = list(api.get_all_genes().keys())


sequences_to_predict = {}
for gene_id in tqdm(all_genes, desc="Loading and cleaning gene sequences"):
    try:
        gene = KeggGene(gene_id, kegg_api=api)
        aa_seq = gene.aa_seq

        if aa_seq and gene.coding_type == "CDS":
            # Clean sequence based on metapredict V3 rules
            processed_seq = V3_version_letters(aa_seq)

            if processed_seq and len(processed_seq) == len(aa_seq):
                sequences_to_predict[gene_id] = processed_seq
            else:
                print(f"   Gene {gene_id} error - V3 sequence length mismatch.")
    except Exception as e:
        print(f"[Error] Loading gene {gene_id}: {e}.")


# Batch Predict Disorder
# Metapredict is highly optimized for dictionary input
print(f"\nPredicting disorder for {len(sequences_to_predict)} sequences...")
disorder_predictions = meta.predict_disorder(sequences_to_predict)

# Merge scores into SNV files
for gene_id, result_list in tqdm(disorder_predictions.items(),
                                 desc="Merging disorder scores into SNV files"):

    # Format path (matching the hsa_XXXX.csv naming convention)
    file_name = gene_id.replace(":", "_") + ".csv"
    snv_file_path = os.path.join(KEGG_GENE_SCORES_P, file_name)

    if not os.path.exists(snv_file_path):
        # This is expected if the gene was not processed in the SNV generation step
        continue

    try:
        # Load the existing SNV table
        snv_df = pd.read_csv(snv_file_path)

        if snv_df.empty:
            continue

        # Metapredict results: index 0 is sequence, index 1 is the np.array of scores
        scores = result_list[1]

        # AA_index in the CSV corresponds to the index in the 'scores' array.
        # We use .map() or .apply() to assign the score based on the AA_index column.
        # We ensure AA_index is treated as integer for indexing.
        snv_df['Disorder_score'] = snv_df['AA_index'].map(
            lambda idx: scores[int(idx)] if int(idx) < len(scores) else np.nan)

        # Assign binary disordered flag (1 if >= threshold, else 0)
        # We handle NaNs by keeping them as NaN or 0 (here we use 0 via comparison)
        snv_df['Is_disordered'] = (snv_df['Disorder_score'] >= DISORDERED_THRESHOLD).astype(int)

        # Save the updated CSV (overwriting the old one)
        snv_df.to_csv(snv_file_path, index=False)

    except Exception as e:
        print(f"[Error] Merging scores for {gene_id}: {e}")

print("\nDisorder score integration complete.")



