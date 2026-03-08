# this creates a csv per cancer of patients and their aggregated scores per pathway
# survival rates of each patient are added to conduct Kaplan-Meier survival analysis on the patients with high vs low scores for each pathway

import os
import numpy as np
import pandas as pd
import pickle
import sys
from definitions import KEGG_PATHWAY_METADATA_FILE, CANCER_PATIENT_SURVIVAL_P, CBIO_CANCER_MUTATIONS_P, CBIO_PATIENT_CLINICAL_STUDIES_P
from os.path import join as pjoin
import glob


# open the metadata pickle dictionary to get the mapping from pathway KEGG ID to the pathway's metadata
with open(KEGG_PATHWAY_METADATA_FILE, 'rb') as f:
    pathway_id_to_metadata = pickle.load(f)

def create_patient_scores_per_pathway_df(cancer_file: str):
    cancer_df = pd.read_csv(
        cancer_file,
        dtype={'StudyId': str, 'PatientId': str},  # Force these to strings
        low_memory=False
    )
    n = (
        cancer_df.groupby('PatientId')['StudyId']
        .nunique()
        .gt(1)
        .sum()
    )

    print(n)
    print(cancer_df.groupby('PatientId')['StudyId'].nunique().value_counts())

    df_rows = []  # list of dictionaries with patient_id and aggregated pathway scores

    # iterate over unique patient, study tuples
    for idx, (patient_id, study_id) in enumerate(
            cancer_df[['PatientId', 'StudyId']]
                    .drop_duplicates()
                    .itertuples(index=False)
    ):
        patient_df = cancer_df[
            (cancer_df['PatientId'] == patient_id) &
            (cancer_df['StudyId'] == study_id)
            ]
        # start row with patient_id and study_id
        patient_row = {"PatientId": patient_id, "StudyId": study_id}

        for pathway in pathway_id_to_metadata.keys():
            # use modules only
            if "hsa" in pathway:
                # TODO: why is that? Don't we want hsa pathways?
                continue

            pathway_genes = pathway_id_to_metadata[pathway]['genes_ids']

            # check if any of the patient's mutated genes are in the pathway
            if not any(gene in pathway_genes for gene in patient_df['KeggId']):
                patient_score = 0
            else:
                pathway_mutations = (
                    patient_df[patient_df['KeggId'].isin(pathway_genes)]
                    ['esm_log_probs']
                )

                # patient_score = pathway_mutations.mean()
                # patient_score = (pathway_mutations > 0.5).sum()
                patient_score = int((pathway_mutations < -15).any())

            # add value under column named as pathway
            patient_row[pathway] = patient_score

        df_rows.append(patient_row)

        if idx % 100 == 0 and idx > 0:
            print(f"    Processed {idx}/{len(cancer_df[['PatientId', 'StudyId']].drop_duplicates())} patients...")

    return pd.DataFrame(df_rows)

def add_patient_survival_data(cancer_patients_df: pd.DataFrame, cancer_type: str):
    study_dfs = []
    for study_id in cancer_patients_df['StudyId'].unique():

        # subset patients of this study only
        study_patients_df = cancer_patients_df[
            cancer_patients_df['StudyId'] == study_id
        ]

        # get survival data for this study
        patient_data_path = pjoin(CBIO_PATIENT_CLINICAL_STUDIES_P, f"{study_id}.csv")
        if not os.path.exists(patient_data_path):
            print(f"Warning: Patient clinical data file not found for study {study_id} at path {patient_data_path}. Skipping survival data merge for this study.")
            continue
        clinical_data_df = pd.read_csv(patient_data_path)
        survival_data_df = add_metastasis_status(clinical_data_df)

        # merge only this study’s patients
        cancer_patients_survival_df = pd.merge(
            study_patients_df,
            survival_data_df,
            on=["PatientId", "StudyId"],
            how="left"
        )

        # append to list
        study_dfs.append(cancer_patients_survival_df)

    if not study_dfs:
        print(f"Error: No survival data was found for any studies in {cancer_type}.")
        return
    merged_df = pd.concat(study_dfs, ignore_index=True)
    # count how many patients have nan metastatic status
    print(f"   {merged_df['Metastatic'].isna().sum()}/{len(merged_df)} patients have no documents metastatic status.")
    # save to csv
    print(f"    Saving merged patient scores and survival data for cancer {cancer_type} to CSV.")
    merged_df.to_csv(pjoin(CANCER_PATIENT_SURVIVAL_P, f"{cancer_type}.csv"), index=False)


def add_metastasis_status(clinical_data_df: pd.DataFrame) -> pd.DataFrame:
    """
    Returns a DataFrame with PatientId and a binary Metastatic column.
    Priority: PATH_M_STAGE > AJCC_PATHOLOGIC_TUMOR_STAGE > not available
    """
    if 'PATH_M_STAGE' in clinical_data_df.columns:
        # M1 = metastatic, M0 = not, MX = unknown
        clinical_data_df['Metastatic'] = clinical_data_df['PATH_M_STAGE'].map(
            lambda x: 1 if str(x).startswith('M1') else (0 if str(x).startswith('M0') else np.nan)
        )

    elif 'AJCC_PATHOLOGIC_TUMOR_STAGE' in clinical_data_df.columns:
        # Stage IV = metastatic
        clinical_data_df['Metastatic'] = clinical_data_df['AJCC_PATHOLOGIC_TUMOR_STAGE'].map(
            lambda x: 1 if 'IV' in str(x).upper() else (
                0 if any(s in str(x).upper() for s in ['I', 'II', 'III']) else np.nan)
        )

    else:
        clinical_data_df['Metastatic'] = np.nan

    # only keep relevant columns if they exist, and drop duplicates
    survival_cols = ['PatientId', 'StudyId', 'Metastatic', 'OS_MONTHS', 'OS_STATUS', 'DFS_MONTHS', 'DFS_STATUS']
    existing_cols = [col for col in survival_cols if col in clinical_data_df.columns]
    return clinical_data_df[existing_cols].drop_duplicates()



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
