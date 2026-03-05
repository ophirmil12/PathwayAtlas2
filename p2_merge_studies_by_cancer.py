# Merges all TCGA studies by cancer type and saves them as CSV files.

from definitions import *
import pandas as pd
import glob
import os
from os.path import join as pjoin
from cbio_api import CbioApi

def merge_studies_by_cancer():
    """Merge all TCGA studies by cancer type and save them as CSV files in path defined by CBIO_CANCER_MUTATIONS.
        @param cbio: CbioApi object
        @param studies_dfs_dict: dictionary of study id to DataFrame of that study.
    """
    mutation_studies_files = glob.glob(pjoin(CBIO_MUTATION_STUDIES_WITH_SEQUENCES_P, "*.csv"))
    if not mutation_studies_files:
        print(f"    No cancer mutation studies CSV files found in the specified directory: {CBIO_MUTATION_STUDIES_WITH_SEQUENCES_P}")
        return

    studies_dfs_dict = {os.path.splitext(os.path.basename(file))[0] : pd.read_csv(file, low_memory=False) for file in mutation_studies_files}

    pan_cancer_dfs_list = []

    for cancer_shortname in CANCER_FULLNAME.keys():
        # find all tcga studies of this cancer type
        output_path = pjoin(CBIO_CANCER_MUTATIONS_P + "_no_dups", f"{cancer_shortname}.csv")

        if os.path.exists(output_path):
            print(f"    Merged mutations for {cancer_shortname} already exists at {output_path}, skipping...")
            pan_cancer_dfs_list.append(pd.read_csv(output_path, low_memory=False))
            continue

        if cancer_shortname not in CANCER_STUDIES.keys():
            print(f"    No tcga studies found for cancer type {cancer_shortname}")
        else:
            print(f"    Merging studies for cancer type: {cancer_shortname}...")

            merged_df = pd.DataFrame()

            for study_id in CANCER_STUDIES[cancer_shortname]:
                if study_id in studies_dfs_dict:
                    study_df = studies_dfs_dict[study_id]
                    print(f"    Merging study {study_id} with {len(study_df)} rows...")
                    merged_df = pd.concat([merged_df, study_df], ignore_index=True)
                    merged_df.drop_duplicates(keep='first', inplace=True, ignore_index=True,
                                              subset=DUPLICATE_EXCLUSION_COLUMNS)
                else:
                    print(f"    Study file {study_id} not found in the loaded dataframes.")

            # Save merged DataFrame to CSV
            pan_cancer_dfs_list.append(merged_df)
            merged_df.to_csv(output_path, index=False)
            print(f"    Saved merged mutations for {cancer_shortname} to {output_path}")

    # Additionally, create a pan-cancer merged file
    print("    Merging all studies into a pan-cancer mutations file...")
    pan_cancer_output_path = pjoin(CBIO_CANCER_MUTATIONS_P + "_no_dups", "pan_cancer.csv")
    if not os.path.exists(pan_cancer_output_path):
        pan_cancer_df = pd.concat(pan_cancer_dfs_list, ignore_index=True)
        pan_cancer_df.to_csv(pan_cancer_output_path, index=False)
        print(f"    Saved pan-cancer merged mutations to {pan_cancer_output_path}")


if __name__ == '__main__':
    print(f"----- Merging processed mutation studies by cancer type -----")
    # Merge studies by cancer type
    merge_studies_by_cancer()
