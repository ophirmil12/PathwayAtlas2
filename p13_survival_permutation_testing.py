"""
p13_survival_permutation_testing.py
====================================
Permutation test for pathway survival signal.

For each significant pathway S from the Cox analysis:
  1. Compute the real Cox z-score (high-burden vs wildtype) on S's actual genes.
  2. Repeat n_permutations times: sample |S| random genes from the gene pool,
     compute the same Cox z-score on that random gene set.
  3. Empirical p-value = fraction of permutations where |z_null| >= |z_real|.

Gene pool (a): union of all genes in pathways significant in the enrichment
               analysis (FILTERED_RESULTS_DISTANCES_P, q_value < 0.05).

Run per cancer via SLURM array:
    python -u p13_survival_permutation_testing.py $SLURM_ARRAY_TASK_ID [n_permutations]
Default n_permutations = 50.
"""

import glob
import os
import pickle
import random
import sys
import warnings
from collections import defaultdict
from os.path import join as pjoin
from statsmodels.stats.multitest import fdrcorrection

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.exceptions import ConvergenceWarning

from definitions import *
from p13_cox_km_analysis import assign_burden_groups, WILDTYPE, LOW_BURDEN, HIGH_BURDEN

SCORING_TYPE = "min_esm_log_probs"
Q_THRESHOLD = 0.05


# ─────────────────────────────────────────────────────────────────────────────
# Gene pool and target pathway helpers
# ─────────────────────────────────────────────────────────────────────────────

def get_genes_in_significant_pathways() -> list:
    """Return all genes found in pathways that are enrichment-significant (q < 0.05)."""
    genes = set()
    with open(KEGG_PATHWAY_METADATA_FILE, 'rb') as f:
        pathway_metadata = pickle.load(f)

    for file in glob.glob(pjoin(FILTERED_RESULTS_DISTANCES_P, "*.csv")):
        df = pd.read_csv(file)
        if 'q_value' not in df.columns:
            continue
        for pathway in df[df['q_value'] < Q_THRESHOLD]['pathway']:
            pathway_genes = pathway_metadata.get(str(pathway), {}).get('genes_ids', [])
            genes.update(pathway_genes)

    print(f"Gene pool: {len(genes)} genes from enrichment-significant pathways")
    return list(genes)


def get_significant_pathways_in_cox(cancer_type: str) -> pd.DataFrame:
    """Return rows for pathways Cox-significant (q_high_vs_wildtype < 0.05) for one cancer."""
    cox_file = pjoin(KAPLAN_MEIER_P, cancer_type, f"{SCORING_TYPE}.csv")
    if not os.path.exists(cox_file):
        print(f"No Cox results found for {cancer_type} at {cox_file}")
        return pd.DataFrame()

    df = pd.read_csv(cox_file)
    sig = df[df['q_high_vs_wildtype'] < Q_THRESHOLD].copy()
    print(f"Found {len(sig)} Cox-significant pathways for {cancer_type}")
    return sig


# ─────────────────────────────────────────────────────────────────────────────
# Efficient patient-score computation
# ─────────────────────────────────────────────────────────────────────────────

def build_patient_gene_scores(mutations_df: pd.DataFrame) -> dict:
    """
    Pre-compute the minimum ESM log prob per (patient, gene) pair.
    Returns: {(PatientId, StudyId): {gene: min_esm_log_probs}}
    """
    patient_gene_lists = defaultdict(lambda: defaultdict(list))
    for _, row in mutations_df.iterrows():
        key = (row['PatientId'], row['StudyId'])
        patient_gene_lists[key][row['KeggId']].append(row['esm_log_probs'])

    return {
        key: {gene: min(scores) for gene, scores in gene_dict.items()}
        for key, gene_dict in patient_gene_lists.items()
    }


def get_scores_for_gene_set(patient_gene_min: dict, gene_set: set) -> pd.DataFrame:
    """
    For each patient, compute min_esm_log_probs across mutations in gene_set.
    NaN for patients with no mutations in gene_set (wildtype).
    """
    records = []
    for (patient_id, study_id), gene_scores in patient_gene_min.items():
        relevant = [score for gene, score in gene_scores.items() if gene in gene_set]
        records.append({
            'PatientId': patient_id,
            'StudyId': study_id,
            SCORING_TYPE: min(relevant) if relevant else np.nan,
        })
    return pd.DataFrame(records)


# ─────────────────────────────────────────────────────────────────────────────
# Cox fitting (mirrors p13_cox_km_analysis._fit_cox)
# ─────────────────────────────────────────────────────────────────────────────

def _fit_cox_z_score(analysis_df: pd.DataFrame) -> float:
    """
    Fit a 3-group Cox model and return the z-score for high_burden vs wildtype.
    Returns NaN if quality filters fail or the model does not converge.
    """
    if len(analysis_df) < 60:
        return np.nan

    group_counts = analysis_df.groupby('burden_group').size()
    group_events = analysis_df.groupby('burden_group')['OS_STATUS'].sum()

    present_groups = analysis_df['burden_group'].unique()
    if (group_counts[present_groups] < 15).any() or (group_events[present_groups] < 10).any():
        return np.nan

    if analysis_df['burden_group'].std() < 0.1:
        return np.nan

    df = analysis_df.copy()
    df['is_low_burden']  = (df['burden_group'] == LOW_BURDEN).astype(int)
    df['is_high_burden'] = (df['burden_group'] == HIGH_BURDEN).astype(int)

    try:
        cph = CoxPHFitter(penalizer=0.1)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", ConvergenceWarning)
            warnings.simplefilter("ignore", RuntimeWarning)
            cph.fit(
                df[['OS_MONTHS', 'OS_STATUS', 'is_low_burden', 'is_high_burden']],
                duration_col='OS_MONTHS',
                event_col='OS_STATUS',
            )
        if cph.summary['coef'].isna().any() or cph.summary['p'].isna().any():
            return np.nan
        return float(cph.summary.loc['is_high_burden', 'z'])
    except Exception:
        return np.nan


def build_analysis_df(all_patients: pd.DataFrame, patient_scores: pd.DataFrame) -> pd.DataFrame:
    """
    Merge survival metadata with patient pathway scores, assign burden groups,
    restrict to non-metastatic patients with complete OS data.
    """
    full_df = all_patients.merge(patient_scores, on=['PatientId', 'StudyId'], how='left')
    # Negate ESM: higher value → more pathogenic (same as run_survival_analysis)
    full_df[SCORING_TYPE] = -full_df[SCORING_TYPE]
    full_df = assign_burden_groups(full_df, SCORING_TYPE)
    return (
        full_df[full_df['Metastatic'] == 0][['OS_MONTHS', 'OS_STATUS', 'burden_group']]
        .dropna()
    )


# ─────────────────────────────────────────────────────────────────────────────
# Permutation test
# ─────────────────────────────────────────────────────────────────────────────

def run_permutation_test(
    cancer_type: str,
    significant_pathways_df: pd.DataFrame,
    gene_pool: list,
    patient_gene_min: dict,
    survival_df: pd.DataFrame,
    n_permutations: int = 100,
) -> pd.DataFrame:
    """
    For each significant pathway, run a permutation test against random gene sets.
    Empirical p-value (two-sided): fraction of |z_null| >= |z_real|.
    """
    with open(KEGG_PATHWAY_METADATA_FILE, 'rb') as f:
        pathway_metadata = pickle.load(f)

    all_patients = (
        survival_df[['PatientId', 'StudyId', 'Metastatic', 'OS_MONTHS', 'OS_STATUS']]
        .drop_duplicates(subset=['PatientId', 'StudyId'])
    )

    results = []

    for _, pathway_row in significant_pathways_df.iterrows():
        pathway_id   = str(pathway_row['pathway'])
        pathway_name = pathway_row.get('pathway_name', pathway_id)
        pathway_genes = set(pathway_metadata.get(pathway_id, {}).get('genes_ids', []))
        n_genes = len(pathway_genes)

        if n_genes == 0:
            print(f"  Skipping {pathway_id}: no genes in metadata")
            continue
        if n_genes > len(gene_pool):
            print(f"  Skipping {pathway_id}: pathway size {n_genes} > gene pool {len(gene_pool)}")
            continue

        # ── Real pathway z-score ──────────────────────────────────────────────
        real_scores   = get_scores_for_gene_set(patient_gene_min, pathway_genes)
        real_analysis = build_analysis_df(all_patients, real_scores)
        z_real        = _fit_cox_z_score(real_analysis)

        if np.isnan(z_real):
            print(f"  Skipping {pathway_id} ({pathway_name}): real Cox model failed filters")
            continue

        print(f"  {pathway_name} (n={n_genes} genes): z_real={z_real:.3f}")

        # ── Null distribution ─────────────────────────────────────────────────
        z_nulls = []
        for _ in range(n_permutations):
            random_genes  = set(random.sample(gene_pool, n_genes))
            perm_scores   = get_scores_for_gene_set(patient_gene_min, random_genes)
            perm_analysis = build_analysis_df(all_patients, perm_scores)
            z_perm        = _fit_cox_z_score(perm_analysis)
            if not np.isnan(z_perm):
                z_nulls.append(z_perm)

        n_valid = len(z_nulls)
        if n_valid == 0:
            print(f"    All permutations failed Cox — skipping")
            continue

        empirical_p = sum(abs(z) >= abs(z_real) for z in z_nulls) / n_valid
        print(f"    empirical_p={empirical_p:.3f}  ({n_valid}/{n_permutations} valid permutations)")

        results.append({
            'cancer_type':        cancer_type,
            'pathway':            pathway_id,
            'pathway_name':       pathway_name,
            'pathway_size':       n_genes,
            'z_real':             z_real,
            'z_null_mean':        np.mean(z_nulls),
            'z_null_std':         np.std(z_nulls),
            'empirical_p':        empirical_p,
            'n_valid_permutations': n_valid,
            'q_high_vs_wildtype': pathway_row.get('q_high_vs_wildtype', np.nan),
        })

    return pd.DataFrame(results)


# ─────────────────────────────────────────────────────────────────────────────
# Entry point
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    args = sys.argv[1:]
    if len(args) < 1:
        print("Usage: python -u p13_survival_permutation_testing.py $SLURM_ARRAY_TASK_ID [n_permutations]")
        sys.exit(1)

    cancer_patients_files = sorted(glob.glob(pjoin(CANCER_PATIENT_SURVIVAL_P, "*.csv")))
    if not cancer_patients_files:
        print(f"No patient survival files found in {CANCER_PATIENT_SURVIVAL_P}")
        sys.exit(1)

    index = int(args[0])
    if index < 0 or index >= len(cancer_patients_files):
        print(f"Index {index} out of range ({len(cancer_patients_files)} files).")
        sys.exit(1)

    n_permutations = int(args[1]) if len(args) > 1 else 100

    cancer_patients_file = cancer_patients_files[index]
    cancer_type = os.path.basename(cancer_patients_file).replace('.csv', '')

    print(f"----- Permutation testing: {cancer_type} ({n_permutations} permutations) -----")

    # ── Gene pool (a) ─────────────────────────────────────────────────────────
    gene_pool = get_genes_in_significant_pathways()
    if not gene_pool:
        print("Gene pool is empty — no enrichment-significant pathways found.")
        sys.exit(1)

    # ── Cox-significant pathways for this cancer ──────────────────────────────
    significant_pathways_df = get_significant_pathways_in_cox(cancer_type)
    if significant_pathways_df.empty:
        print(f"No Cox-significant pathways for {cancer_type}.")
        sys.exit(0)

    # ── Raw mutation data for this cancer ─────────────────────────────────────
    cancer_mutations_files = sorted(glob.glob(pjoin(CBIO_CANCER_MUTATIONS_P, "*.csv")))
    cancer_mutations_file = next(
        (f for f in cancer_mutations_files
         if os.path.splitext(os.path.basename(f))[0] == cancer_type),
        None,
    )
    if cancer_mutations_file is None:
        print(f"No mutations file found for {cancer_type} in {CBIO_CANCER_MUTATIONS_P}")
        sys.exit(1)

    print(f"Loading mutations from {cancer_mutations_file}...")
    mutations_df = pd.read_csv(
        cancer_mutations_file,
        dtype={'StudyId': str, 'PatientId': str},
        low_memory=False,
        usecols=['PatientId', 'StudyId', 'KeggId', 'esm_log_probs'],
    )
    mutations_df = mutations_df.dropna(subset=['KeggId', 'esm_log_probs'])

    print("Building patient-gene score index...")
    patient_gene_min = build_patient_gene_scores(mutations_df)

    # ── Survival metadata ─────────────────────────────────────────────────────
    survival_df = pd.read_csv(cancer_patients_file)

    # ── Run permutation tests ─────────────────────────────────────────────────
    results_df = run_permutation_test(
        cancer_type=cancer_type,
        significant_pathways_df=significant_pathways_df,
        gene_pool=gene_pool,
        patient_gene_min=patient_gene_min,
        survival_df=survival_df,
        n_permutations=n_permutations,
    )

    if results_df.empty:
        print(f"No permutation results for {cancer_type}.")
        sys.exit(0)

    # fdr correction for empirical p-values
    results_df['corrected_empirical_p'] = fdrcorrection(results_df['empirical_p'])[1]

    out_dir = pjoin(KAPLAN_MEIER_P, cancer_type)
    os.makedirs(out_dir, exist_ok=True)
    out_path = pjoin(out_dir, f"permutation_results_{n_permutations}.csv")
    results_df.to_csv(out_path, index=False)
    print(f"Saved permutation results -> {out_path}")

    # for each cancer type, print all pathways with significant empirical p-values after correction
    sig_empirical = results_df[results_df['corrected_empirical_p'] < 0.05]
    if not sig_empirical.empty:
        print(f"\nSignificant pathways for {cancer_type} after FDR correction (q < 0.05):")
        for _, row in sig_empirical.iterrows():
            print(f"  {row['pathway_name']} (empirical p={row['empirical_p']:.3e}, corrected q={row['corrected_empirical_p']:.3e})")
    else:
        print(f"\nNo pathways with significant empirical p-values after FDR correction for {cancer_type}.")
