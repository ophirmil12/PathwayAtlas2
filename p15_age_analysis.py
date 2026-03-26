# This file performs the analysis of cancer signatures in younger population (<50)

import os

from definitions import *
import pandas as pd
import numpy as np
import glob
from os.path import join as pjoin
import sys

def create_age_group_csv(cancer_scores_file: str) -> (pd.DataFrame, pd.DataFrame):
    # Read the cancer scores file
    df = pd.read_csv(cancer_scores_file)
    df["age"] = np.nan  # Initialize age column with NaN values
    
    studies = df['StudyId'].unique()
    for study_id in studies:
        # Get clinical data for the patients in the cancer scores file
        clinical_data_df = pd.read_csv(pjoin(CBIO_PATIENT_CLINICAL_STUDIES_P, f"{study_id}.csv"))
        if "AGE" not in clinical_data_df.columns:
            print(f"AGE column not found in clinical data for study {study_id}. Skipping this study.")
            continue

        for _, row in clinical_data_df.iterrows():
            age = row['AGE']
            patient_id = row['PatientId']
            if pd.isna(age):
                continue
            age = int(age)
            # Update the age column in the cancer mutations dataframe for the corresponding patient
            df.loc[(df['PatientId'] == patient_id) & (df['StudyId'] == study_id), 'age'] = age

    # Group patients into age groups
    df_under_50 = df[df['age'] < 50]
    df_50_and_over = df[df['age'] >= 50]

    # Save the age group dataframes to new CSV files
    under_50_file = pjoin(AGE_ANALYSIS_P, f"{os.path.basename(cancer_scores_file).replace('.csv', '_under_50.csv')}")
    over_50_file = pjoin(AGE_ANALYSIS_P, f"{os.path.basename(cancer_scores_file).replace('.csv', '_50_and_over.csv')}")
    df_under_50.to_csv(under_50_file, index=False)
    df_50_and_over.to_csv(over_50_file, index=False)
    print(f"Created age group CSV files: {under_50_file}, {over_50_file}")

    return df_under_50, df_50_and_over
            
                

if __name__ == '__main__':
    # Get the cancer mutations file based on SLURM_ARRAY_TASK_ID
    args = sys.argv[1:]
    if len(args) < 1:
        print("Usage: python -u p7_bootstrap_pathways.py '$SLURM_ARRAY_TASK_ID'")
        sys.exit(1)
    
    index = int(args[0])
    
    all_pathway_files = glob.glob(os.path.join(KEGG_PATHWAY_SCORES_P, f"*.csv"))
    cancer_types_for_analysis = ["colorectal", "brca", "uterine", "liver", "lung"]
    
    # Calculate pathway and cancer indices
    num_cancer_types = len(cancer_types_for_analysis)
    pathway_index = index // num_cancer_types
    cancer_index = index % num_cancer_types
    
    # Check index bounds
    if pathway_index >= len(all_pathway_files) or cancer_index >= len(cancer_types_for_analysis):
        print(f"Index {index} is out of range. There are only {len(all_pathway_files)} pathway files and {len(cancer_types_for_analysis)} cancer types.")
        sys.exit(1)
    
    # Get the specific files
    pathway_scores_file = sorted(all_pathway_files)[pathway_index]
    cancer_scores_file = sorted(glob.glob(pjoin(CBIO_CANCER_MUTATIONS_P, f"{cancer_types_for_analysis[index]}.csv")))[0]
    pathway_name = os.path.splitext(os.path.basename(pathway_scores_file))[0]

    under_50_file = pjoin(AGE_ANALYSIS_P, f"{os.path.basename(cancer_scores_file).replace('.csv', '_under_50.csv')}")
    over_50_file = pjoin(AGE_ANALYSIS_P, f"{os.path.basename(cancer_scores_file).replace('.csv', '_50_and_over.csv')}")

    if not os.path.exists(under_50_file) or not os.path.exists(over_50_file):
        print(f"Age group CSV files not found for {os.path.basename(cancer_scores_file)}. Creating age group CSV files...")
        df_under_50, df_50_and_over = create_age_group_csv(cancer_scores_file)
    else:
        print(f"Age group CSV files already exist for {os.path.basename(cancer_scores_file)}. Skipping age group CSV file creation.")
        df_under_50 = pd.read_csv(under_50_file)
        df_50_and_over = pd.read_csv(over_50_file)
