# This file performs the analysis of cancer signatures in younger population (<50)

import os

from definitions import *
import pandas as pd
import numpy as np
import glob
from os.path import join as pjoin
import sys
from p7_bootstrap_pathways import *

def create_age_group_csv(cancer_scores_file: str, outpath_under_50: str, outpath_50_and_over: str) -> (pd.DataFrame, pd.DataFrame):
    # Read the cancer scores file
    df = pd.read_csv(cancer_scores_file)
    df["age"] = np.nan  # Initialize age column with NaN values
    #df["early_stage"] = np.nan
    studies = df['StudyId'].unique()

    for study_id in studies:
        # Get clinical data for the patients in the cancer scores file
        clinical_data_df = pd.read_csv(pjoin(CBIO_PATIENT_CLINICAL_STUDIES_P, f"{study_id}.csv"))
        # Filter by stage so only early stage patients are analyzed
        # early_stage_df = get_early_stage_study(clinical_data_df)
        # if early_stage_df.empty:
        #     print(f"    Stage column not found in clinical data for study {study_id}.")

        if "AGE" not in clinical_data_df.columns:
            print(f"    Age column not found in clinical data for study {study_id}.")
            print(clinical_data_df.columns)
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
    df_under_50.to_csv(outpath_under_50, index=False)
    df_50_and_over.to_csv(outpath_50_and_over, index=False)
    print(f"Created age group CSV files: Under 50 with {len(df_under_50)} rows,\n Over 50 with {len(df_50_and_over)} rows.")

    return df_under_50, df_50_and_over
    

def get_early_stage_study(clinical_data_df: pd.DataFrame):
    early_stage_regex = r'\bI{1,2}(?![I])\b'
    for col in STAGE_COLUMNS:
        if col in clinical_data_df.columns:
            return clinical_data_df[clinical_data_df[col].str.contains(early_stage_regex, case=False, na=False, regex=True)]
    print(list(clinical_data_df.columns))
    return pd.DataFrame()



if __name__ == '__main__':
 
    # Get the cancer mutations file based on SLURM_ARRAY_TASK_ID
    args = sys.argv[1:]
    if len(args) < 1:
        print("Usage: python -u p7_bootstrap_pathways.py '$SLURM_ARRAY_TASK_ID'")
        sys.exit(1)
    
    index = int(args[0])
    
    all_pathway_files = glob.glob(os.path.join(KEGG_PATHWAY_SCORES_P, f"*.csv"))
    all_pathway_files = [file for file in all_pathway_files if os.path.splitext(os.path.basename(file))[0] not in EXCLUDED_HSA]
    all_pathway_names = [os.path.splitext(os.path.basename(f))[0] for f in sorted(all_pathway_files)]
    cancer_types_for_analysis = ["colorectal", "brca", "liver", "pan_cancer", "gynecological"]
    
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
    cancer_scores_file = sorted(glob.glob(pjoin(CBIO_CANCER_MUTATIONS_P, f"{cancer_types_for_analysis[cancer_index]}.csv")))[0]
    pathway_name = os.path.splitext(os.path.basename(pathway_scores_file))[0]
    cancer_name = os.path.splitext(os.path.basename(cancer_scores_file))[0]

    under_50_file = pjoin(AGE_ANALYSIS_P, f"{cancer_name}_under_50.csv")
    over_50_file = pjoin(AGE_ANALYSIS_P, f"{cancer_name}_50_and_over.csv")

    if not os.path.exists(under_50_file) or not os.path.exists(over_50_file):
        print(f"Age group CSV files not found for {cancer_name}. Creating age group CSV files...")
        df_under_50, df_50_and_over = create_age_group_csv(cancer_scores_file, under_50_file, over_50_file)
    else:
        print(f"Age group CSV files already exist for {cancer_name}. Skipping age group CSV file creation.")

    # Perform the bootstrap pipeline (distance calculation + p_value calculation) for both age groups
    # For under 50
    print(f"Performing bootstrap analysis for patients under 50 for {cancer_name}...")
    under_50_distances_outpath = pjoin(AGE_ANALYSIS_DISTANCES_P, f"{cancer_name}_under_50.csv")
    if not os.path.exists(under_50_distances_outpath):
        under_50_distances_df = pd.DataFrame({'pathway': all_pathway_names,
                'wasserstein_distance': [np.nan] * len(all_pathway_names),
                'delta_means': [np.nan] * len(all_pathway_names),
                'p_value': [np.nan] * len(all_pathway_names)
            })
        under_50_distances_df.to_csv(under_50_distances_outpath, index=False)
        
    bootstrap_pathway_for_cancer(pathway_scores_file, under_50_file, under_50_distances_outpath)
        

    # For 50 and over
    print(f"Performing bootstrap analysis for patients 50 and over for {cancer_name}...")
    over_50_distances_outpath = pjoin(AGE_ANALYSIS_DISTANCES_P, f"{cancer_name}_over_50.csv")
    if not os.path.exists(over_50_distances_outpath):
        over_50_distances_df = pd.DataFrame({'pathway': all_pathway_names,
                'wasserstein_distance': [np.nan] * len(all_pathway_names),
                'delta_means': [np.nan] * len(all_pathway_names),
                'p_value': [np.nan] * len(all_pathway_names)
            })
        over_50_distances_df.to_csv(over_50_distances_outpath, index=False)

    bootstrap_pathway_for_cancer(pathway_scores_file, over_50_file, over_50_distances_outpath)
