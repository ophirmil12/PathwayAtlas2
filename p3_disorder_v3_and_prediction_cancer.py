# Go through all genes, for each create V3 version, predict the disorder scores (for the cancer data)

import os
import re
import sys
import pandas as pd
import numpy as np
import metapredict as meta
from tqdm import tqdm

from definitions import V3_version_letters, CBIO_MUTATION_STUDIES_WITH_SEQUENCES_P, DISORDERED_THRESHOLD


def extract_aa_index(variant_str):
    """
    Parses variant string (e.g., 'L145Q') to return 0-based index.
    Returns None if parsing fails.
    """
    try:
        match = re.search(r'\d+', str(variant_str))
        if match:
            return int(match.group()) - 1
    except:
        pass
    return None


def process_cancer_file(file_idx):
    # 1. Setup Folder
    cancer_files = sorted([f for f in os.listdir(CBIO_MUTATION_STUDIES_WITH_SEQUENCES_P) if f.endswith('.csv')])     # TODO check that is the right folder

    if file_idx >= len(cancer_files):
        print(f"Index {file_idx} out of range.")
        return

    target_file = cancer_files[file_idx]
    file_path = os.path.join(CBIO_MUTATION_STUDIES_WITH_SEQUENCES_P, target_file)
    print(f"Processing file {file_idx}: {target_file}")

    # 2. Load Data
    df = pd.read_csv(file_path)
    if df.empty or 'Sequence' not in df.columns:
        print(f"Skipping {target_file}: Empty or missing 'Sequence' column.")
        return

    # 3. Collect Unique Sequences directly from the DataFrame
    # Create a mapping of KeggId -> Sequence for unique genes in this file
    unique_genes = df[['KeggId', 'Sequence']].dropna().drop_duplicates('KeggId')

    sequences_to_predict = {}
    for _, row in tqdm(unique_genes.iterrows(), total=len(unique_genes), desc="Preparing sequences"):
        kegg_id = row['KeggId']
        raw_seq = row['Sequence']

        # Clean the sequence using V3 rules
        processed_seq = V3_version_letters(raw_seq)

        # Ensure we don't have empty strings and lengths match
        if processed_seq and len(processed_seq) == len(raw_seq):
            sequences_to_predict[kegg_id] = processed_seq

    # 4. Batch Predict Disorder
    # Metapredict returns {kegg_id: [sequence, scores_array]}
    print(f"Predicting disorder for {len(sequences_to_predict)} unique proteins...")
    disorder_results = meta.predict_disorder(sequences_to_predict)

    # 5. Map scores back to the DataFrame rows
    disorder_scores = []
    is_disordered_flags = []

    for _, row in tqdm(df.iterrows(), total=len(df), desc="Mapping scores to mutations"):
        kegg_id = row['KeggId']
        variant = row['Variant']

        score = np.nan
        is_disordered = 0

        if kegg_id in disorder_results:
            aa_idx = extract_aa_index(variant)
            scores_array = disorder_results[kegg_id][1]

            if aa_idx is not None and 0 <= aa_idx < len(scores_array):
                score = scores_array[aa_idx]
                is_disordered = 1 if score >= DISORDERED_THRESHOLD else 0
            # Optional: handle cases where Variant index > Sequence length

        disorder_scores.append(score)
        is_disordered_flags.append(is_disordered)

    # 6. Update DataFrame and Save
    df['Disorder_score'] = disorder_scores
    df['Is_disordered'] = is_disordered_flags

    df.to_csv(file_path, index=False)
    print(f"Successfully updated {target_file}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python p3_disorder_v3_and_prediction_cancer.py <file_index>")
    else:
        idx = int(sys.argv[1])
        process_cancer_file(idx)