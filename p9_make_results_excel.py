# Prepare an Excel file with a summary of the results

#  We want the file to be as followed:
#  For each sheet, named {cancer type}:
#       Have a table with columns:
#           Pathway ID
#           Pathway name (description from metadata)
#           Wasserstein distance
#           Delta of expectations
#           Q-value
#           Coverage
#           More statistics...
# Rows are sorted by Q-value

import os
import glob
import pickle
import pandas as pd
import numpy as np
from definitions import (
    RESULTS_DISTANCES_P,
    KEGG_PATHWAY_METADATA_P,
    RESULTS_P
)


def create_summary_excel():
    """
    Aggregates per-cancer result CSVs into a single multi-sheet Excel file.
    Renames columns to remove internal shorthand and sorts by significance.
    """

    # 1. Load Pathway Metadata for descriptions
    with open(KEGG_PATHWAY_METADATA_P, 'rb') as f:
        pathway_metadata = pickle.load(f)

    # 2. Identify Files
    result_files = glob.glob(os.path.join(RESULTS_DISTANCES_P, "*.csv"))
    if not result_files:
        print(f"No result files found in {RESULTS_DISTANCES_P}.")
        return

    output_path = os.path.join(RESULTS_P, "PathwayAtlas2_Final_Results.xlsx")

    # 3. Setup Writer
    writer = pd.ExcelWriter(output_path, engine='xlsxwriter')

    # Column mapping to clean up the headers for the final report
    # This maps internal CSV column names to "Safe" Excel headers
    RENAME_DICT = {
        'pathway': 'Pathway ID',
        'pathway_name': 'Pathway Description',
        'wasserstein_distance': 'Wasserstein Distance (W1)',
        'delta_means': 'Pathogenic Shift (Delta)',
        'q_value': 'Significance (Q-value)',
        'p_value': 'Raw P-value',
        'num_mutations': 'Observed Mutations',
        'total_aa_length': 'Pathway AA Length',
        'num_genes': 'Gene Count',
        'num_covered_genes': 'Reliably Covered Genes',
        'coverage_ratio': 'Coverage Percentage'
    }

    print(f"Generating summary for {len(result_files)} cohorts...")

    for file_path in sorted(result_files):
        # Cancer type becomes the tab name
        tab_name = os.path.basename(file_path).replace(".csv", "").upper()[:31]

        try:
            df = pd.read_csv(file_path)
            
            # filter out rows with missing q_values
            df = df.dropna(subset=['q_value'])
            
            if df.empty:
                continue

            # 4. Cleanup & Merge Metadata
            # Find the ID column regardless of exact name
            id_col = 'pathway'  # Based on your code

            df['pathway_name'] = df[id_col].apply(
                lambda x: pathway_metadata.get(str(x), {}).get('name', 'Unknown').split(' - Homo')[0]
            )

            # =Handle division by zero and convert to 0-100%
            num_cov = pd.to_numeric(df['num_covered_genes'], errors='coerce').fillna(0)
            num_tot = pd.to_numeric(df['num_genes'], errors='coerce').replace(0, np.nan)

            # Calculate ratio and round to 2 decimal places
            df['coverage_ratio'] = (num_cov / num_tot).fillna(0) * 100

            # 5. Sorting
            df = df.sort_values(by=['q_value', 'delta_means'], ascending=[True, False])

            # 6. Rename and Reorder
            df = df.rename(columns=RENAME_DICT)

            # Select and order primary columns
            primary_cols = [
                'Pathway ID', 'Pathway Description',
                'Significance (Q-value)', 'Pathogenic Shift (Delta)',
                'Wasserstein Distance (W1)', 'Coverage Percentage'
            ]

            # Keep any other columns (coverage stats, etc.) to the right
            existing_primary = [c for c in primary_cols if c in df.columns]
            other_cols = [c for c in df.columns if c not in existing_primary]
            df = df[existing_primary + other_cols]

            # 7. Export to Excel
            df.to_excel(writer, sheet_name=tab_name, index=False)

            # 8. Visual Formatting
            worksheet = writer.sheets[tab_name]
            worksheet.freeze_panes(1, 0)  # Freeze header
            worksheet.set_column('A:A', 12)  # ID width
            worksheet.set_column('B:B', 50)  # Description width
            worksheet.set_column('C:E', 18)  # Primary metrics width
            worksheet.set_column('F:Z', 14)  # Stats width

        except Exception as e:
            print(f"Skipping {tab_name} due to error: {e}")

    writer.close()
    print(f"\nReport generated: {output_path}")


if __name__ == "__main__":
    create_summary_excel()