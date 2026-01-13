# Merges all TCGA studies by cancer type and saves them as CSV files.

from definitions import *
import pandas as pd
import glob
import os
from os.path import join as pjoin
from cbio_api import CbioApi

def merge_studies_by_cancer(cbio: CbioApi):
    """Merge all TCGA studies by cancer type and save them as CSV files in path defined by CBIO_CANCER_MUTATIONS.
        @param cbio: CbioApi object
        @param studies_dfs_dict: dictionary of study id to DataFrame of that study.
    """
    mutation_studies_files = glob.glob(pjoin(CBIO_MUTATION_STUDIES_WITH_SEQUENCES_P, "*.csv"))
    if not mutation_studies_files:
        print(f"    No cancer mutation studies CSV files found in the specified directory: {CBIO_MUTATION_STUDIES_WITH_SEQUENCES_P}")
        return

    studies_dfs_dict = {os.path.splitext(os.path.basename(file))[0] : pd.read_csv(file) for file in mutation_studies_files}

    cancer_name_dict = cbio.cancer_types_dict()

    for cancer_shortname in cancer_name_dict.values():
        # find all tcga studies of this cancer type
        cancer_shortname = cancer_shortname.lower()
        output_path = pjoin(CBIO_CANCER_MUTATIONS_P, f"{cancer_shortname}.csv")

        if os.path.exists(output_path):
            print(f"    Merged mutations for {cancer_shortname} already exists at {output_path}, skipping...")
            continue

         # Get all studies for this cancer type
        study_ids, study_names = cbio.all_studies_by_keyword(cancer_shortname)
        tcga_study_ids = [id for id in study_ids if id in studies_dfs_dict.keys()]

        if not tcga_study_ids:
            print(f"    No tcga studies found for cancer type {cancer_shortname}")
        else:
            print(f"    Merging studies for cancer type: {cancer_shortname}...")

            merged_df = pd.DataFrame()

            for study_id in tcga_study_ids:
                if study_id in studies_dfs_dict:
                    study_df = studies_dfs_dict[study_id]
                    merged_df = pd.concat([merged_df, study_df], ignore_index=True)
                    merged_df.drop_duplicates(keep='first', inplace=True, ignore_index=True,
                                              subset=DUPLICATE_EXCLUSION_COLUMNS)
                else:
                    print(f"    Study file {study_id} not found in the loaded dataframes.")

            # Save merged DataFrame to CSV
            merged_df.to_csv(output_path, index=False)
            print(f"    Saved merged mutations for {cancer_shortname} to {output_path}")


if __name__ == '__main__':
    print(f"----- Merging processed mutation studies by cancer type -----")

    # Initialize cBioPortal API
    cbio = CbioApi()

    # Merge studies by cancer type
    merge_studies_by_cancer(cbio)
