# Training basic regression models over the CLinVar data (CLINVAR_DATA_TABLE_P)
# 2 models, one for ordered regions, one for disordered regions


import os
import pickle
import pandas as pd
from sklearn.linear_model import LogisticRegression


from definitions import *



def get_regression_over_clinvar() -> dict[int, LogisticRegression]:
    """
    Train logistic regression classifiers on ClinVar mutations:
    - one for ordered regions (is_disordered=0)
    - one for disordered regions (is_disordered=1)
    models_path (str, optional): Full path to save/load the models pickle file.
    """
    # --- Check if models are already trained and saved ---
    if os.path.exists(CLINVAR_MODELS_P):
        print(f"Loading pre-trained ClinVar regression models from: {CLINVAR_MODELS_P}")
        with open(CLINVAR_MODELS_P, 'rb') as f:
            models = pickle.load(f)
        return models

    print("Training simple regressors over clinvar...")

    df = pd.read_csv(CLINVAR_DATA_TABLE_P)
    df['is_disordered'] = df['is_disordered']
    df = df.dropna(subset=['binary_label'])

    models = {}
    esm_log_probs = 'wt_not_nadav_marginals_base_wt_score'

    # ---- Train per-flag models ----
    for disorder_flag in [0, 1]:
        subset = df[df['is_disordered'] == disorder_flag].copy()
        subset = subset.dropna(subset=[esm_log_probs])

        if subset.empty:
            continue

        # remove extreme outliers
        lower, upper = subset[esm_log_probs].quantile([0.001, 0.999])
        subset = subset[(subset[esm_log_probs] >= lower) &
                        (subset[esm_log_probs] <= upper)]

        X = subset[[esm_log_probs]].values
        y = subset['binary_label'].astype(int).values

        model = LogisticRegression(solver="lbfgs", max_iter=1000)
        model.fit(X, y)
        models[disorder_flag] = model

    # --- Save the newly trained models to the specified path ---
    print(f"Saving trained models to: {CLINVAR_MODELS_P}")
    with open(CLINVAR_MODELS_P, 'wb') as f:
        pickle.dump(models, f)

    return models

