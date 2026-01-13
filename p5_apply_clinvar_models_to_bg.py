# Use the trained ClinVar regressors on the SNVs CSVs data


import os
import pickle
import pandas as pd
import numpy as np
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed

from definitions import *
from kegg_api import KeggApi
from p5_clinvar_reggresors import get_regression_over_clinvar


def apply_model_to_gene(kegg_id, models):
    """
    Worker function: Loads a gene's SNV CSV and applies the
    corresponding logistic regression models.
    """
    file_name = kegg_id.replace(":", "_") + ".csv"
    csv_path = os.path.join(KEGG_GENE_SCORES_P, file_name)

    if not os.path.exists(csv_path):
        return None

    try:
        df = pd.read_csv(csv_path)

        # We need both esm_log_probs and Is_disordered to proceed
        if 'esm_log_probs' not in df.columns or 'Is_disordered' not in df.columns:
            return f"Missing columns in {kegg_id}"

        # Initialize result column
        df['pathogenic_prob'] = np.nan

        # Process Ordered and Disordered parts separately
        for flag in [0, 1]:
            if flag not in models:
                continue

            # Select subset
            mask = (df['Is_disordered'] == flag) & (df['esm_log_probs'].notna())
            subset = df[mask]

            if subset.empty:
                continue

            # Reshape features for sklearn [n_samples, 1]
            X = subset[['esm_log_probs']].values

            # predict_proba returns [prob_0, prob_1]. We want prob_1 (Pathogenic)
            probs = models[flag].predict_proba(X)[:, 1]

            # Map results back to the main dataframe
            df.loc[mask, 'pathogenic_prob'] = probs

        # Save updated CSV
        df.to_csv(csv_path, index=False)
        return kegg_id

    except Exception as e:
        return f"Error processing {kegg_id}: {str(e)}"


def run_model_application():
    # 1. Load pre-trained models
    # This calls the function from your attached p5_clinvar_reggresors.py
    models = get_regression_over_clinvar()

    if 0 not in models or 1 not in models:
        print("Error: Models for both ordered (0) and disordered (1) regions are required.")
        return

    # 2. Get gene list
    api = KeggApi()
    print("Fetching gene list...")
    gene_ids = list(api.get_all_genes().keys())

    print(f"Applying ClinVar regressors to {len(gene_ids)} gene files...")

    # 3. Multithreaded processing
    with ThreadPoolExecutor(max_workers=8) as executor:
        futures = {executor.submit(apply_model_to_gene, gid, models): gid for gid in gene_ids}

        success_count = 0
        for future in tqdm(as_completed(futures), total=len(gene_ids), desc="Predicting Pathogenicity"):
            res = future.result()
            if res:
                if res.startswith("Error") or res.startswith("Missing"):
                    # Optional: print(res)
                    pass
                else:
                    success_count += 1

    print(f"\nFinished! Applied models to {success_count} gene files.")
    print(f"Results saved in {KEGG_GENE_SCORES_P}")


if __name__ == "__main__":
    run_model_application()