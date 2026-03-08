# Merges all TCGA studies by cancer type and saves them as CSV files.

from definitions import *
import pandas as pd
import glob
import os
from os.path import join as pjoin
from cbio_api import CbioApi

def merge_studies_by_cancer():
    """
    Merge all TCGA studies by cancer type and save them as CSV files in path defined by CBIO_CANCER_MUTATIONS.
    """
    mutation_studies_files = glob.glob(pjoin(CBIO_MUTATION_STUDIES_WITH_SEQUENCES_P, "*.csv"))
    if not mutation_studies_files:
        print(f"    No cancer mutation studies CSV files found in the specified directory: {CBIO_MUTATION_STUDIES_WITH_SEQUENCES_P}")
        return

    studies_dfs_dict = {os.path.splitext(os.path.basename(file))[0] : pd.read_csv(file, low_memory=False) for file in mutation_studies_files}

    pan_cancer_dfs_dict = {}

    for cancer_shortname in CANCER_FULLNAME.keys():
        # find all tcga studies of this cancer type
        output_path = pjoin(CBIO_CANCER_MUTATIONS_P, f"{cancer_shortname}.csv")

        if os.path.exists(output_path):
            print(f"    Merged mutations for {cancer_shortname} already exists at {output_path}, skipping...")
            pan_cancer_dfs_dict[cancer_shortname] = pd.read_csv(output_path, low_memory=False)
            continue

        if cancer_shortname not in CANCER_STUDIES_DICT.keys():
            print(f"    No tcga studies found for cancer type {cancer_shortname}")
        else:
            print(f"    Merging studies for cancer type: {cancer_shortname}...")
            merged_studies_df = pd.concat([studies_dfs_dict[study_id] for study_id in CANCER_STUDIES_DICT[cancer_shortname] if study_id in studies_dfs_dict], ignore_index=True)
            print(f"    Merged studies for {cancer_shortname} has {len(merged_studies_df)} rows before dropping duplicates...")
            merged_studies_df.drop_duplicates(keep='first', inplace=True, ignore_index=True,
                                              subset=DUPLICATE_EXCLUSION_COLUMNS)

            # Save merged DataFrame to CSV
            pan_cancer_dfs_dict[cancer_shortname] = merged_studies_df
            merged_studies_df.to_csv(output_path, index=False)
            print(f"    Saved merged mutations for {cancer_shortname} with {len(merged_studies_df)} rows to {output_path}")


    # Additionally, create a pan-cancer merged file
    print("    Merging all studies into a pan-cancer mutations file...")
    pan_cancer_output_path = pjoin(CBIO_CANCER_MUTATIONS_P, "pan_cancer.csv")
    if not os.path.exists(pan_cancer_output_path):
        pan_cancer_df = pd.concat(pan_cancer_dfs_dict.values(), ignore_index=True)
        print(f"    Pan-cancer merged mutations has {len(pan_cancer_df)} rows before dropping duplicates...")
        pan_cancer_df.drop_duplicates(keep='first', inplace=True, ignore_index=True,
                                  subset=DUPLICATE_EXCLUSION_COLUMNS)
        pan_cancer_df.to_csv(pan_cancer_output_path, index=False)
        print(f"    Saved pan-cancer merged mutations with {len(pan_cancer_df)} rows to {pan_cancer_output_path}")

    # Merge cancer subtypes into their parent cancer type
    for parent_cancer, subtypes in CANCER_SUBTYPES.items():
        parent_output_path = pjoin(CBIO_CANCER_MUTATIONS_P, f"{parent_cancer}.csv")
        if os.path.exists(parent_output_path):
            print(f"    Merged mutations for {parent_cancer} already exists at {parent_output_path}, skipping...")
            continue

        print(f"    Merging subtypes for parent cancer type: {parent_cancer}...")
        parent_df = pd.concat([pan_cancer_dfs_dict[subtype] for subtype in subtypes if subtype in pan_cancer_dfs_dict], ignore_index=True)
        print(f"    Merged subtypes for {parent_cancer} has {len(parent_df)} rows before dropping duplicates...")
        parent_df.drop_duplicates(keep='first', inplace=True, ignore_index=True,
                                    subset=DUPLICATE_EXCLUSION_COLUMNS)

        # Save merged DataFrame to CSV
        parent_df.to_csv(parent_output_path, index=False)
        print(f"    Saved merged mutations for parent cancer type {parent_cancer} with {len(parent_df)} rows to {parent_output_path}")


if __name__ == '__main__':
    print(f"----- Merging processed mutation studies by cancer type -----")
    # Merge studies by cancer type
    merge_studies_by_cancer()