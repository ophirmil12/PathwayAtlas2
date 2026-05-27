# this creates a csv per cancer of patients and their aggregated scores per pathway
# survival rates of each patient are added to conduct Kaplan-Meier survival analysis on the patients with high vs low scores for each pathway

import os
import numpy as np
import pandas as pd
import pickle
import sys
from definitions import *
from os.path import join as pjoin
import glob


# open the metadata pickle dictionary to get the mapping from pathway KEGG ID to the pathway's metadata
with open(FILTERED_KEGG_PATHWAY_METADATA_FILE, 'rb') as f:
    pathway_id_to_metadata = pickle.load(f)

def create_patient_scores_per_pathway_df(cancer_file: str):
    cancer_df = pd.read_csv(
        cancer_file,
        dtype={'StudyId': str, 'PatientId': str},
        low_memory=False
    )

    # Drop duplicate mutations for the same patient across overlapping study versions.
    # DUPLICATE_EXCLUSION_COLUMNS = [Chr, Ref, Alt, Protein, Variant, PatientId]
    cancer_df = cancer_df.drop_duplicates(subset=DUPLICATE_EXCLUSION_COLUMNS)

    unique_patients = cancer_df["PatientId"].unique()
    print(f"For cancer {os.path.basename(cancer_file)}: "
          f"{len(unique_patients)} patients, {len(cancer_df)} mutations (after dedup)")

    df_rows = []
    pathways_list = list(pathway_id_to_metadata.keys())

    for idx, patient_id in enumerate(unique_patients):
        patient_df = cancer_df[cancer_df['PatientId'] == patient_id]

        for pathway in pathways_list:
            pathway_genes = pathway_id_to_metadata[pathway]['genes_ids']

            if not any(gene in pathway_genes for gene in patient_df['KeggId']):
                min_esm_log_probs_score = np.nan
                path_count_score = np.nan
            else:
                patient_pathway_df = patient_df[patient_df['KeggId'].isin(pathway_genes)]
                min_esm_log_probs_score = patient_pathway_df['esm_log_probs'].min()
                path_count_score = patient_pathway_df['pathogenic_prob'].gt(0.5).fillna(False).sum()

            pathway_description = pathway_id_to_metadata.get(
                str(pathway), {}
            ).get('name', 'Unknown').split(' - Homo')[0]

            df_rows.append({
                "PatientId": patient_id,
                "pathway_id": pathway,
                "pathway_name": pathway_description,
                "min_esm_log_probs": min_esm_log_probs_score,
                "path_count": path_count_score,
            })

        if idx % 100 == 0 and idx > 0:
            print(f"    Processed {idx}/{len(unique_patients)} patients...")

    return pd.DataFrame(df_rows)

def add_patient_survival_data(cancer_patients_df: pd.DataFrame, cancer_type: str,
                              study_ids: list, output_path: str):
    # Load clinical data from every study for this cancer type and concatenate.
    # The same patient often appears in multiple overlapping study versions; we
    # keep the row with the most complete survival data for each patient.
    all_clinical = []
    for study_id in study_ids:
        patient_data_path = pjoin(CBIO_PATIENT_CLINICAL_STUDIES_P, f"{study_id}.csv")
        if not os.path.exists(patient_data_path):
            print(f"Warning: clinical data not found for study {study_id}, skipping.")
            continue
        clinical_data_df = pd.read_csv(patient_data_path, low_memory=False)
        all_clinical.append(add_metastasis_status(clinical_data_df))

    if not all_clinical:
        print(f"Error: No clinical data found for any study in {cancer_type}.")
        return

    combined = pd.concat(all_clinical, ignore_index=True)

    # Prefer rows with more non-null survival fields, then deduplicate per patient.
    key_cols = ['Metastatic', 'OS_MONTHS', 'OS_STATUS']
    combined['_completeness'] = combined[[c for c in key_cols if c in combined.columns]].notna().sum(axis=1)
    combined = (combined
                .sort_values('_completeness', ascending=False)
                .drop_duplicates(subset=['PatientId'])
                .drop(columns=['_completeness']))

    merged_df = pd.merge(cancer_patients_df, combined, on='PatientId', how='left')
    merged_df.dropna(subset=['Metastatic', 'OS_MONTHS', 'OS_STATUS'], inplace=True)

    print(f"   {merged_df['Metastatic'].isna().sum()}/{len(merged_df)} patients have no documented metastatic status.")
    print(f"   Saving merged patient scores and survival data for cancer {cancer_type} to CSV.")
    merged_df.to_csv(output_path, index=False)


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
    survival_cols = ['PatientId', 'Metastatic', 'OS_MONTHS', 'OS_STATUS']
    existing_cols = [col for col in survival_cols if col in clinical_data_df.columns]
    return clinical_data_df[existing_cols].drop_duplicates()



if __name__ == '__main__':
    # Get the cancer mutations file based on SLURM_ARRAY_TASK_ID
    args = sys.argv[1:]
    if len(args) < 1:
        print("Usage: python -u p9_patient_scores_and_survival_per_cancer.py '$SLURM_ARRAY_TASK_ID'")
        sys.exit(1)

    cancer_results_files = glob.glob(pjoin(CBIO_CANCER_MUTATIONS_P, "*.csv"))
    if not cancer_results_files:
        print(f"No cancer mutations CSV files found in the specified directory: {CBIO_CANCER_MUTATIONS_P}")
        sys.exit(1)

    index = int(args[0])
    if index < 0 or index >= len(cancer_results_files):
        print(f"Index {index} is out of range. There are only {len(cancer_results_files)} cancer mutation files.")
        sys.exit(1)

    cancer_mutations_file = sorted(cancer_results_files)[index]

    cancer_type = os.path.basename(cancer_mutations_file).split('.')[0]  # get cancer type from file name

    patient_scores_csv_outpath = pjoin(CANCER_PATIENT_SURVIVAL_P, f"{cancer_type}_dedup.csv")

    print(f"----- Summarizing patient scores per pathway and adding survival data for cancer {cancer_mutations_file} -----")

    cancer_patients_df = create_patient_scores_per_pathway_df(cancer_mutations_file)

    if cancer_type == 'pan_cancer':
        study_ids = [s for studies in CANCER_STUDIES_DICT.values() for s in studies]
    else:
        study_ids = CANCER_STUDIES_DICT.get(cancer_type, [])
    if not study_ids:
        print(f"Warning: no study list found in CANCER_STUDIES_DICT for '{cancer_type}'. "
              f"Clinical data merge may be incomplete.")
    add_patient_survival_data(cancer_patients_df, cancer_type, study_ids, patient_scores_csv_outpath)


