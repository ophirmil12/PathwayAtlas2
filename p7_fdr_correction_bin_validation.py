# For each bin-size results file produced by p7_bootstrap_bin_validation.py:
#   1. Copy coverage columns (num_mutations, total_aa_length, num_genes, num_covered_genes,
#      pathway_name) from the canonical pan_cancer distances file in RESULTS_DISTANCES_P
#   2. Save an unfiltered backup copy of the enriched bin results
#   3. Filter out pathways whose gene-coverage percentage is below PATHWAY_COVERAGE_THRESHOLD
#   4. Apply BH FDR correction to the surviving p-values and write q_value back to the file

import os
import glob
import pandas as pd
from scipy.stats import false_discovery_control
from definitions import *

# ------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------
BIN_VALIDATION_DIR = pjoin(RESULTS_P, "p7_bin_validation_analysis")
BACKUP_DIR = pjoin(BIN_VALIDATION_DIR, "backup_with_coverage")

# The canonical pan_cancer distances file already has coverage columns
# written by p6A_calc_coverage_and_filter.py.
COVERAGE_COLUMNS = [
    'num_mutations',
    'total_aa_length',
    'num_genes',
    'num_covered_genes',
    'pathway_name',
]


def load_coverage_source(cancer_name: str = "pan_cancer") -> pd.DataFrame | None:
    """
    Returns the filtered pan_cancer distances DataFrame that contains the
    per-pathway coverage statistics added by p6A.
    The 'pathway' column is used as the join key (matches the bin-results files).
    """
    pattern = pjoin(RESULTS_DISTANCES_P, f"{cancer_name}*.csv")
    matches = glob.glob(pattern)
    if not matches:
        print(f"ERROR: No distances file found for '{cancer_name}' in {RESULTS_DISTANCES_P}")
        return None

    src_path = sorted(matches)[0]
    df = pd.read_csv(src_path)

    missing = [c for c in COVERAGE_COLUMNS if c not in df.columns]
    if missing:
        print(f"ERROR: Coverage source file is missing columns: {missing}")
        print(f"       Please run p6A_calc_coverage_and_filter.py first.")
        return None

    # Normalise the join-key column name to 'pathway'
    if 'pathway' not in df.columns:
        candidate = df.columns[0]
        print(f"WARNING: 'pathway' column not found in source; using '{candidate}' as join key.")
        df = df.rename(columns={candidate: 'pathway'})

    print(f"Loaded coverage source: {src_path}  ({len(df)} pathways after p6A filtering)")
    return df[['pathway'] + COVERAGE_COLUMNS]


def add_coverage(bin_df: pd.DataFrame, coverage_df: pd.DataFrame) -> pd.DataFrame:
    """
    Left-joins coverage columns onto the bin results.
    Pathways absent from the coverage source get NaN (they will be filtered out).
    """
    # Drop stale coverage columns if they exist
    bin_df = bin_df.drop(columns=COVERAGE_COLUMNS, errors='ignore')
    merged = bin_df.merge(coverage_df, on='pathway', how='left')
    return merged


def filter_by_coverage(df: pd.DataFrame) -> pd.DataFrame:
    """
    Keeps only pathways where (num_covered_genes / num_genes) * 100
    >= PATHWAY_COVERAGE_THRESHOLD (40 %).
    Rows with missing coverage info are dropped.
    """
    df = df.copy()
    df = df.dropna(subset=['num_genes', 'num_covered_genes'])

    coverage_pct = (
        df['num_covered_genes']
        / df['num_genes'].replace(0, pd.NA)
    ) * 100
    coverage_pct = coverage_pct.fillna(0)

    mask = coverage_pct >= PATHWAY_COVERAGE_THRESHOLD
    return df[mask].copy()


def apply_fdr(df: pd.DataFrame, bin_label: str) -> pd.DataFrame:
    """
    Applies BH FDR correction to the p_value column.
    Adds a q_value column and returns the updated DataFrame.
    """
    if 'p_value' not in df.columns or df['p_value'].isnull().all():
        print(f"  [bins={bin_label}] ERROR: p_value column missing or all-NaN; skipping FDR.")
        return df

    df = df.dropna(subset=['p_value']).copy()

    try:
        q_values = false_discovery_control(df['p_value'].tolist(), method='bh')
    except Exception as e:
        print(f"  [bins={bin_label}] ERROR during FDR correction: {e}")
        return df

    df['q_value'] = q_values
    return df


def run_fdr_bin_validation():
    os.makedirs(BACKUP_DIR, exist_ok=True)

    # 1. Load coverage reference once
    coverage_df = load_coverage_source(cancer_name="pan_cancer")
    if coverage_df is None:
        return

    # 2. Find all bin-results files  (e.g. pan_cancer_bins_10.csv … pan_cancer_bins_100.csv)
    bin_files = sorted(glob.glob(pjoin(BIN_VALIDATION_DIR, "pan_cancer_bins_*.csv")))
    if not bin_files:
        print(f"No bin-validation result files found in {BIN_VALIDATION_DIR}")
        return

    print(f"\nFound {len(bin_files)} bin-result file(s) to process.\n")

    for bin_path in bin_files:
        filename = os.path.basename(bin_path)
        bin_label = filename.replace("pan_cancer_bins_", "").replace(".csv", "")
        print(f"--- Processing: {filename}  (bins={bin_label}) ---")

        bin_df = pd.read_csv(bin_path)
        print(f"  Loaded {len(bin_df)} pathway results.")

        # Step 1: Attach coverage columns
        enriched_df = add_coverage(bin_df, coverage_df)
        n_with_coverage = enriched_df[COVERAGE_COLUMNS[0]].notna().sum()
        print(f"  Coverage matched for {n_with_coverage}/{len(enriched_df)} pathways.")

        # Step 2: Save unfiltered backup
        backup_path = pjoin(BACKUP_DIR, filename)
        enriched_df.to_csv(backup_path, index=False)
        print(f"  Backup (unfiltered + coverage) saved to: {backup_path}")

        # Step 3: Filter by 40 % coverage threshold
        filtered_df = filter_by_coverage(enriched_df)
        dropped = len(enriched_df) - len(filtered_df)
        print(f"  Coverage filter: kept {len(filtered_df)} pathways, dropped {dropped}.")

        # Step 4: BH FDR correction
        result_df = apply_fdr(filtered_df, bin_label)
        significant = (result_df.get('q_value', pd.Series(dtype=float)) < 0.05).sum()
        print(f"  FDR done. Significant pathways (q<0.05): {significant}/{len(result_df)}")

        # Overwrite the original file with the filtered + FDR-corrected data
        result_df.to_csv(bin_path, index=False)
        print(f"  Saved final results to: {bin_path}\n")

    print("All bin-validation files processed.")


if __name__ == "__main__":
    run_fdr_bin_validation()