import os
import glob
import pandas as pd
from scipy.stats import false_discovery_control
from definitions import *


def filter_results(min_mutations: int):
    """
    Saves a copy of the original results after bootstrapping to a backup folder,
    and filters the working files by num_mutations >= min_mutations.
    """
    gene_results_dir = RESULTS_GENE_LEVEL_P
    backup_dir = f"{gene_results_dir}_unfiltered"

    # Ensure backup directory exists
    os.makedirs(backup_dir, exist_ok=True)

    cancer_gene_files = sorted(glob.glob(pjoin(gene_results_dir, "*_gene_distances.csv")))

    if not cancer_gene_files:
        print(f"No result files found in {gene_results_dir}")
        return

    for file_path in cancer_gene_files:
        base_name = os.path.basename(file_path)

        # 1. Read the original data
        df = pd.read_csv(file_path)

        # 2. Save the backup (unfiltered) version
        backup_path = pjoin(backup_dir, base_name)
        df.to_csv(backup_path, index=False)

        # 3. Filter the dataframe
        if "num_mutations" in df.columns:
            initial_count = len(df)
            df = df[df["num_mutations"] >= min_mutations]
            final_count = len(df)

            # 4. Overwrite original file with filtered data
            df.to_csv(file_path, index=False)

            print(f"Processed {base_name}:")
            print(f"  - Kept {final_count}/{initial_count} genes (min_mutations >= {min_mutations})")
            print(f"  - Backup saved to {backup_path}")
        else:
            print(f"  - Skip filtering for {base_name}: 'num_mutations' column not found.")


def perform_gene_fdr_correction(cancer_gene_file: str):
    """
    Reads the gene-level results for a specific cancer,
    calculates q-values from p-values, and updates the CSV.
    """
    try:
        df = pd.read_csv(cancer_gene_file)

        if 'p_value' not in df.columns or df['p_value'].isnull().all():
            print(f"Skipping {os.path.basename(cancer_gene_file)}: No p-values found.")
            return

        # Filter out NaNs for the calculation to prevent errors
        mask = df['p_value'].notna()
        p_values = df.loc[mask, 'p_value'].values

        if len(p_values) == 0:
            return

        # Calculate q-values using Benjamini-Hochberg (BH)
        q_values = false_discovery_control(p_values, method='bh')

        # Map them back to the original dataframe
        df.loc[mask, 'q_value'] = q_values

        # Save updated file
        df.to_csv(cancer_gene_file, index=False)

        sig_count = (q_values < 0.05).sum()
        print(f"Processed {os.path.basename(cancer_gene_file)}: {sig_count} significant genes (q < 0.05).")

    except Exception as e:
        print(f"Error processing {cancer_gene_file}: {e}")


def run_fdr_on_all_results():
    # 1. Define the directory where gene-level results are stored
    gene_results_dir = RESULTS_GENE_LEVEL_P

    # 2. Find all relevant CSV files
    cancer_gene_files = sorted(glob.glob(pjoin(gene_results_dir, "*_gene_distances.csv")))

    if not cancer_gene_files:
        print(f"No result files found in {gene_results_dir}")
        return

    print(f"Starting FDR correction for {len(cancer_gene_files)} files...")

    # 3. Loop through files one by one
    for file_path in cancer_gene_files:
        perform_gene_fdr_correction(file_path)

    print("\nAll FDR corrections completed.")


if __name__ == '__main__':
    filter_results(min_mutations=10)
    run_fdr_on_all_results()