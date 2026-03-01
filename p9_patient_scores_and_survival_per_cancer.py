# this creates a csv per cancer of patients and their aggregated scores per pathway
# survival rates of each patient are added to conduct Kaplan-Meier survival analysis on the patients with high vs low scores for each pathway

import os
import numpy as np
import pandas as pd
from cbio_api import CbioApi
import pickle
import sys
from definitions import KEGG_PATHWAY_METADATA_FILE, CANCER_PATIENT_SURVIVAL_P, CBIO_CANCER_MUTATIONS_P
from os.path import join as pjoin
import glob


# open the metadata pickle dictionary to get the mapping from pathway KEGG ID to the pathway's metadata
with open(KEGG_PATHWAY_METADATA_FILE, 'rb') as f:
    pathway_id_to_metadata = pickle.load(f)

def create_patient_scores_per_pathway_df(cancer_file: str):
    cancer_df = pd.read_csv(cancer_file)
    df_rows = []  # list of dictionaries with patient_id and aggregated pathway scores

    # iterate over unique patients
    for patient_id in cancer_df['PatientId'].unique():
        patient_df = cancer_df[cancer_df['PatientId'] == patient_id]

        # start row with patient_id
        patient_row = {"PatientId": patient_id}

        for pathway in pathway_id_to_metadata.keys():
            pathway_genes = pathway_id_to_metadata[pathway]['genes_ids']

            # check if any of the patient's mutated genes are in the pathway
            if not any(gene in pathway_genes for gene in patient_df['KeggId']):
                patient_mean_score = np.nan
            else:
                patient_mean_score = (
                    patient_df[patient_df['KeggId'].isin(pathway_genes)]
                    ['pathogenic_prob']
                    .mean()
                )

            # add value under column named as pathway
            patient_row[pathway] = patient_mean_score

        df_rows.append(patient_row)

    cancer_patients_df = pd.DataFrame(df_rows)
    # delete columns of pathways with all NaN values
    cancer_patients_df.dropna(axis=1, how='all', inplace=True)
    return cancer_patients_df


def add_patient_survival_data(cancer_patients_df: pd.DataFrame, cancer_type: str):
    # get the survival data for the patients in the cancer type
    cbio = CbioApi()
    survival_data_df = cbio.get_cancer_survival_data(cancer_type)
    cancer_patients_survival_df = pd.merge(cancer_patients_df, survival_data_df, on='PatientId', how='left')
    cancer_patients_survival_df.to_csv(pjoin(CANCER_PATIENT_SURVIVAL_P, f"{cancer_type}.csv"), index=False)


if __name__ == '__main__':
    # Get the cancer mutations file based on SLURM_ARRAY_TASK_ID
    args = sys.argv[1:]
    if len(args) < 1:
        print("Usage: python -u p9_patient_scores_and_survival_per_cancer.py '$SLURM_ARRAY_TASK_ID'")
        sys.exit(1)

    cancer_results_files = glob.glob(pjoin(CBIO_CANCER_MUTATIONS_P, "*.csv"))  # should be 38 files
    if not cancer_results_files:
        print(f"No cancer mutations CSV files found in the specified directory: {CBIO_CANCER_MUTATIONS_P}")
        sys.exit(1)

    index = int(args[0])
    if index < 0 or index >= len(cancer_results_files):
        print(f"Index {index} is out of range. There are only {len(cancer_results_files)} cancer mutation files.")
        sys.exit(1)

    cancer_mutations_file = sorted(cancer_results_files)[index]

    print(f"----- Summarizing patient scores per pathway and adding survival data for cancer {cancer_mutations_file} -----")
    cancer_patients_df = create_patient_scores_per_pathway_df(cancer_mutations_file)
    cancer_type = os.path.basename(cancer_mutations_file).split('.')[0]  # get cancer type from file name
    add_patient_survival_data(cancer_patients_df, cancer_type)
