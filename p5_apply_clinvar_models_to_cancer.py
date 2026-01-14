# Use the trained ClinVar regressors on the cancer CSVs data

import os
import sys
import pandas as pd
import numpy as np
from definitions import *
from p5_clinvar_reggresors import get_regression_over_clinvar



def process_cancer_file_pathogenicity(file_index):
    """
    Applies ClinVar regression models to a single cancer mutation file.
    Uses esm_log_probs and Is_disordered to compute pathogenic_prob.
    """
    # 1. Load pre-trained models
    # This automatically loads from CLINVAR_MODELS_P or trains if missing
    models = get_regression_over_clinvar()

    if 0 not in models or 1 not in models:
        print("Error: Models for both ordered (0) and disordered (1) regions are required.")
        sys.exit(1)

    # 2. Identify target file
    cancer_files = sorted([f for f in os.listdir(CBIO_CANCER_MUTATIONS_UNMERGED) if f.endswith('.csv')])

    if file_index < 0 or file_index >= len(cancer_files):
        print(f"Error: File index {file_index} out of range.")
        sys.exit(1)

    filename = cancer_files[file_index]
    csv_path = os.path.join(CBIO_CANCER_MUTATIONS_UNMERGED, filename)
    print(f"[{file_index}] Applying models to: {filename}")

    # 3. Load Data
    try:
        df = pd.read_csv(csv_path)
    except Exception as e:
        print(f"Error reading CSV: {e}")
        sys.exit(1)

    if df.empty:
        print("File is empty, skipping.")
        return

    # 4. Check for required columns (generated in p3 and p4 steps)
    required = ['esm_log_probs', 'Is_disordered']
    for col in required:
        if col not in df.columns:
            print(f"Error: Missing column '{col}' in {filename}. Ensure p3 and p4 scripts have run.")
            sys.exit(1)

    # 5. Application Logic
    # Initialize/Reset the pathogenic_prob column
    df['pathogenic_prob'] = np.nan

    # Process Ordered (flag=0) and Disordered (flag=1) mutations separately
    for flag in [0, 1]:
        # Create a mask for rows belonging to this disorder category with valid ESM scores
        mask = (df['Is_disordered'] == flag) & (df['esm_log_probs'].notna())

        if not mask.any():
            continue

        # Extract features (sklearn expects 2D array [n_samples, 1])
        X = df.loc[mask, ['esm_log_probs']].values

        # Get pathogenicity probability (class 1)
        # predict_proba returns [prob_benign, prob_pathogenic]
        probs = models[flag].predict_proba(X)[:, 1]

        # Map probabilities back to the main DataFrame
        df.loc[mask, 'pathogenic_prob'] = probs

    # 6. Save results
    try:
        df.to_csv(csv_path, index=False)
        print(f"Successfully updated {filename} with pathogenic_prob.")
    except Exception as e:
        print(f"Error saving updated CSV: {e}")
        sys.exit(1)


if __name__ == "__main__":
    # Expecting Slurm Array Task ID as the first argument
    if len(sys.argv) < 2:
        print("Usage: python p5_apply_clinvar_models_to_cancer.py <file_index>")
        sys.exit(1)

    try:
        target_idx = int(sys.argv[1])
    except ValueError:
        print("Error: file_index must be an integer.")
        sys.exit(1)

    process_cancer_file_pathogenicity(target_idx)