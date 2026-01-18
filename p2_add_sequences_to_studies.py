# Add sequences to the downloaded tcga mutation studies and merge them by cancer type

import numpy as np
import pandas as pd
from definitions import *
from cbio_api import CbioApi
import glob
import os
from os.path import join as pjoin
from kegg_api import KeggApi, KeggGene
import re
from typing import Dict
import sys

# global dictionary of hugo id to a list of kegg hsa ids
hugo_to_hsa = {}

def is_valid_sequence(sequence: str, variant: str) -> bool:
    """Check if the sequence is valid for the given variant.
        @param sequence: gene sequence
        @param variant: mutation variant information (format XnY)
        @return: True if valid, False otherwise
    """
    if pd.isnull(variant) or not variant:
        return False

    match = re.match(VARIATION_REGEX, variant)
    if not match:
        print(f"    ERROR: Invalid mutation format: {variant}")
        return False

    orig_aa, loc, mut_aa = match.groups()
    idx = int(loc) - 1
    if idx >= len(sequence):
        return False

    return sequence[idx] == orig_aa

def gene_to_sequence(gene: str, variant: str) -> (str, str):
    """Get the sequence for a given gene using the cBioPortal API.
        @param gene: gene symbol
        @param variant: mutation variant information (format XnY)
        @return: hsa id of gene, gene sequence as a string
    """
    if not gene or not variant:
        return np.nan, np.nan

    hsa_ids = []

    if gene in hugo_to_hsa.keys():
        hsa_ids = hugo_to_hsa[gene]
    else:
        try:
            hsa_ids = KeggApi.hugo_to_kegg_hsa(gene)
        except Exception as e:
            print(f"    ERROR: Could not retrieve hsa id for {gene}: {e}")
            return np.nan, np.nan
        hugo_to_hsa[gene] = hsa_ids

    for hsa_id in hsa_ids:

        try:
            kegg_gene = KeggGene(hsa_id)
            sequence = kegg_gene.aa_seq
        except Exception as e:
            print(f"    ERROR: Could not retrieve sequence for {hsa_id}: {e}")
            continue

        if sequence and is_valid_sequence(sequence, variant):
            return hsa_id, sequence

    print(f"    No valid sequence found for gene: {gene} with variant: {variant}")
    return np.nan, np.nan

def add_sequences_to_study(file: str) -> None:
    """Add sequences to all genes in the tcga mutation studies in the given filename.
        @param file: path to the mutation study CSV file
    """
    output_path = pjoin(CBIO_MUTATION_STUDIES_WITH_SEQUENCES_P, os.path.basename(file))
    if os.path.exists(output_path):
        print(f"----- Sequences already added to mutations in file: {file}. Skipping... -----")
        return

    df = pd.read_csv(file)
    print(f"----- Adding sequences to mutations in file: {file}... -----")

    df_with_sequences = df.copy()

    for idx, row in df.iterrows():
        ref_name = row['Protein']
        ref_mut = row['Variant']

        hsa_id, sequence = gene_to_sequence(ref_name, ref_mut)
        df_with_sequences.at[idx, 'KeggId'] = hsa_id
        df_with_sequences.at[idx, 'Sequence'] = sequence

        if idx % 100 == 0:
            print(f"    Processed {idx} / {len(df)} rows.")

    df_with_sequences.to_csv(output_path, index=False)


if __name__ == '__main__':
    args = sys.argv[1:]
    if len(args) < 1:
        print("Usage: python -u p2_add_sequences_to_studies.py '$SLURM_ARRAY_TASK_ID'")
        sys.exit(1)

    mutation_study_files = glob.glob(pjoin(CBIO_MUTATION_STUDIES, "*.csv"))
    if not mutation_study_files:
        print(f"No cancer mutation studies CSV files found in the specified directory: {CBIO_MUTATION_STUDIES}")
        sys.exit(1)

    index = int(args[0])
    if index < 0 or index >= len(mutation_study_files):
        print(f"Index {index} is out of range. There are only {len(mutation_study_files)} mutation study files.")
        sys.exit(1)

    mutation_study_file = sorted(mutation_study_files)[index]

    # Initialize cBioPortal API
    cbio = CbioApi()

    # Add sequences to mutation studies
    add_sequences_to_study(mutation_study_file)
