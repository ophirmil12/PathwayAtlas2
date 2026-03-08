import glob
import os
from definitions import *
import pandas as pd
from os.path import join as pjoin
import numpy as np

if __name__ == '__main__':
    print("----- Creating results CSV files for all cancer types and pathways -----")

    all_pathway_files = glob.glob(os.path.join(KEGG_PATHWAY_SCORES_P, f"*.csv"))
    all_pathway_names = [os.path.splitext(os.path.basename(f))[0] for f in sorted(all_pathway_files)]
    all_cancer_files = glob.glob(os.path.join(CBIO_CANCER_MUTATIONS_P, f"*.csv"))

    for cancer_scores_file in sorted(all_cancer_files):

        cancer_distances_path = pjoin(RESULTS_DISTANCES_P, os.path.basename(cancer_scores_file))
        if not os.path.exists(cancer_distances_path):
            print(f"Creating new distances file for {os.path.basename(cancer_scores_file)}...")
            cancer_distances_df = pd.DataFrame({'pathway': all_pathway_names,
                'wasserstein_distance': [np.nan] * len(all_pathway_names),
                'delta_means': [np.nan] * len(all_pathway_names),
                'p_value': [np.nan] * len(all_pathway_names)
            })
            cancer_distances_df.to_csv(cancer_distances_path, index=False)
        else:
            print(f"Distances file already exists for {os.path.basename(cancer_scores_file)}, skipping creation...")