"""
p16_gene_level_cox.py – Per-gene multivariate Cox regression for pathway disruption
====================================================================================

Design rationale
----------------
For a given pathway, every sufficiently-mutated gene is encoded as one
continuous covariate in a *single* multivariate Cox model:

  Not mutated  → score = 0  (wildtype reference)
  Mutated      → score = –min_esm_log_probs across the patient's mutations
                 in that gene  (negated so higher = more pathogenic, because
                 raw ESM log-probs are negative and more negative = worse)

All genes enter the model simultaneously.  Each gene's hazard ratio is
therefore adjusted for co-occurring mutations in every other pathway gene.
This is the correct treatment of the "co-occurring mutations" problem: the
coefficient for gene A is estimated while conditioning on gene B, C, … scores,
exactly as in any other multivariate regression.

Validity as an explanatory test
---------------------------------
Yes – this is a valid exploratory approach for pinpointing which genes within a
pathway drive the survival signal, because:
  1. Multivariate conditioning adjusts for co-mutations (see above).
  2. Using the continuous ESM score (not just binary mutated/not) gives more
     power and a biologically meaningful gradient (mild → severe mutation).
  3. FDR correction across genes prevents inflation from the many-gene search.

Main caveats:
  - Power per gene is low; pan-cancer pooling (default) helps.
  - The model does not adjust for clinical confounders (stage, treatment).
  - L2 penalization (penalizer=0.1) shrinks estimates toward 0; treat small
    effects with caution.  Higher penalizer = more shrinkage = more stability.

Usage
-----
  python p16_gene_level_cox.py [--cancer CANCER] [--pathways P1 P2 ...]

  --cancer        single cancer type (e.g. "luad"); default = pan-cancer
  --pathways      space-separated KEGG pathway IDs
                  (default: hsa04080 hsa04722 hsa04216 hsa04137)
  --min-mutated   minimum patients mutated per gene to include it (default: 10)
  --penalizer     L2 penalizer for CoxPH (default: 0.1)

Output
------
  Per-pathway CSV  → <KAPLAN_MEIER_P>/gene_level_cox/<label>/<pathway_id>_gene_cox.csv
  Combined CSV     → <KAPLAN_MEIER_P>/gene_level_cox/<label>/all_pathways_gene_cox.csv

  Columns: pathway_id, pathway_name, kegg_id, coef, hr, hr_lower95, hr_upper95,
           p, q, n_mutated, n_patients, n_events
"""

import argparse
import glob
import os
import pickle
import sys
import warnings
from os.path import join as pjoin

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.exceptions import ConvergenceWarning
from statsmodels.stats.multitest import fdrcorrection

from definitions import (
    CANCER_PATIENT_SURVIVAL_P,
    CBIO_CANCER_MUTATIONS_P,
    KEGG_PATHWAY_METADATA_FILE,
    KAPLAN_MEIER_P,
    CRITICAL_EDGES_P,
)


# ─────────────────────────────────────────────────────────────────────────────
# Data loading
# ─────────────────────────────────────────────────────────────────────────────

def load_clinical(cancer_types: list) -> pd.DataFrame:
    """One OS row per (PatientId, StudyId), non-metastatic filter NOT applied yet."""
    frames = []
    required = {"PatientId", "StudyId", "Metastatic", "OS_MONTHS", "OS_STATUS"}
    for ct in cancer_types:
        path = pjoin(CANCER_PATIENT_SURVIVAL_P, f"{ct}.csv")
        if not os.path.exists(path):
            print(f"  [warn] no clinical file for {ct}, skipping")
            continue
        header = pd.read_csv(path, nrows=0).columns.tolist()
        if not required.issubset(header):
            print(f"  [warn] {ct}.csv missing OS columns, skipping")
            continue
        df = pd.read_csv(path, usecols=list(required))
        frames.append(df.drop_duplicates(subset=["PatientId", "StudyId"]))
    if not frames:
        raise FileNotFoundError("No clinical files found.")
    combined = pd.concat(frames, ignore_index=True)
    # Patients may appear in multiple cancer-type files if studies overlap;
    # keep the first occurrence.
    return combined.drop_duplicates(subset=["PatientId", "StudyId"])


def load_mutations(cancer_types: list, gene_ids: set) -> pd.DataFrame:
    """All mutation rows for genes in gene_ids across the given cancer types.

    Includes the Protein column so callers can build a gene-name → KeggId map.
    """
    frames = []
    for ct in cancer_types:
        path = pjoin(CBIO_CANCER_MUTATIONS_P, f"{ct}.csv")
        if not os.path.exists(path):
            continue
        df = pd.read_csv(path, usecols=[
            "PatientId", "StudyId", "KeggId", "Protein", "esm_log_probs"
        ])
        frames.append(df[df["KeggId"].isin(gene_ids)])
    if not frames:
        return pd.DataFrame(
            columns=["PatientId", "StudyId", "KeggId", "Protein", "esm_log_probs"]
        )
    return pd.concat(frames, ignore_index=True)


def load_bottleneck_kegg_ids(pathway_id: str, protein_to_kegg: dict) -> set:
    """
    Read the saved bottleneck CSV for pathway_id and return a set of KeggIds.

    The bottleneck CSVs store gene names (Protein column).  protein_to_kegg
    is a dict built from the mutation data: {protein_name: kegg_id}.
    Genes not found in protein_to_kegg are warned and skipped.
    """
    bottleneck_path = pjoin(CRITICAL_EDGES_P, f"{pathway_id}_bottlenecks.csv")
    if not os.path.exists(bottleneck_path):
        print(f"  [warn] no bottleneck file at {bottleneck_path}; using all pathway genes")
        return set()

    bn_df = pd.read_csv(bottleneck_path)
    gene_names = bn_df["gene"].tolist()
    kegg_ids = set()
    missing = []
    for name in gene_names:
        kid = protein_to_kegg.get(name)
        if kid:
            kegg_ids.add(kid)
        else:
            missing.append(name)
    if missing:
        print(f"  [warn] {len(missing)} bottleneck genes not found in mutation data "
              f"(no observed mutations?): {missing[:10]}")
    print(f"  Bottleneck genes: {len(gene_names)} named → {len(kegg_ids)} with KeggId")
    return kegg_ids


# ─────────────────────────────────────────────────────────────────────────────
# Matrix construction
# ─────────────────────────────────────────────────────────────────────────────

def build_gene_matrix(
    clinical_df: pd.DataFrame,
    mutations_df: pd.DataFrame,
    min_mutated: int,
) -> tuple:
    """
    Returns (matrix_df, gene_cols).

    matrix_df has OS columns + one column per gene:
      0           – wildtype (not mutated)
      positive    – –min_esm_log_probs  (negated; higher = more pathogenic)

    Genes with < min_mutated patients are dropped.
    """
    # Per (patient, gene): most pathogenic mutation = minimum esm_log_probs
    per_patient_gene = (
        mutations_df
        .groupby(["PatientId", "StudyId", "KeggId"])["esm_log_probs"]
        .min()
        .reset_index()
        .rename(columns={"esm_log_probs": "min_esm"})
    )

    # Wide: rows = patients, columns = gene KeggIds
    wide = per_patient_gene.pivot_table(
        index=["PatientId", "StudyId"],
        columns="KeggId",
        values="min_esm",
        aggfunc="min",
    )
    wide.columns.name = None
    wide = wide.reset_index()

    merged = clinical_df.merge(wide, on=["PatientId", "StudyId"], how="left")
    gene_cols = [c for c in merged.columns if c.startswith("hsa:")]

    for col in gene_cols:
        merged[col] = -merged[col]   # negate: higher = more pathogenic
        merged[col] = merged[col].fillna(0.0)  # 0 = not mutated (wildtype)

    # Filter by mutation frequency
    n_mutated = (merged[gene_cols] > 0).sum()
    keep = n_mutated[n_mutated >= min_mutated].index.tolist()
    dropped = len(gene_cols) - len(keep)
    if dropped:
        print(f"    Dropped {dropped} genes with < {min_mutated} mutated patients; "
              f"{len(keep)} remain.")

    keep_cols = ["PatientId", "StudyId", "Metastatic", "OS_MONTHS", "OS_STATUS"] + keep
    return merged[keep_cols], keep


# ─────────────────────────────────────────────────────────────────────────────
# Cox regression
# ─────────────────────────────────────────────────────────────────────────────

def run_cox(
    matrix_df: pd.DataFrame,
    gene_cols: list,
    pathway_id: str,
    penalizer: float = 0.1,
) -> pd.DataFrame:
    """
    Fit one multivariate Cox model with all gene_cols as covariates.
    Non-metastatic patients only.  Returns tidy per-gene DataFrame or None.
    """
    analysis = (
        matrix_df[matrix_df["Metastatic"] == 0]
        [["OS_MONTHS", "OS_STATUS"] + gene_cols]
        .dropna(subset=["OS_MONTHS", "OS_STATUS"])
    )

    n_patients = len(analysis)
    n_events = int(analysis["OS_STATUS"].sum())
    print(f"    {n_patients} non-metastatic patients, {n_events} events")

    if n_events < 20:
        print(f"    Skipping: too few events ({n_events} < 20)")
        return None

    # Drop genes that end up with < 3 mutated patients after non-metastatic filter
    active = [c for c in gene_cols if (analysis[c] > 0).sum() >= 3]
    inactive = len(gene_cols) - len(active)
    if inactive:
        print(f"    Dropped {inactive} more genes with < 3 mutated non-metastatic patients")
    if not active:
        print("    No genes remain after post-filter, skipping.")
        return None

    cph = CoxPHFitter(penalizer=penalizer)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", ConvergenceWarning)
            warnings.simplefilter("ignore", RuntimeWarning)
            cph.fit(
                analysis[["OS_MONTHS", "OS_STATUS"] + active],
                duration_col="OS_MONTHS",
                event_col="OS_STATUS",
            )
    except Exception as e:
        print(f"    Cox fitting failed: {e}")
        return None

    if cph.summary["coef"].isna().any() or np.isnan(cph.log_likelihood_):
        print("    Skipping: Cox did not converge (NaN coefficients).")
        return None

    summary = cph.summary[
        ["coef", "exp(coef)", "exp(coef) lower 95%", "exp(coef) upper 95%", "p"]
    ].copy()
    summary.index.name = "kegg_id"
    summary = summary.reset_index()
    summary.columns = ["kegg_id", "coef", "hr", "hr_lower95", "hr_upper95", "p"]

    # Attach per-gene mutation counts
    n_mutated_series = (analysis[active] > 0).sum()
    n_mutated_df = n_mutated_series.reset_index()
    n_mutated_df.columns = ["kegg_id", "n_mutated"]
    summary = summary.merge(n_mutated_df, on="kegg_id", how="left")

    summary["pathway_id"] = pathway_id
    summary["n_patients"] = n_patients
    summary["n_events"] = n_events
    return summary


def attach_gene_names(result_df: pd.DataFrame, mutations_df: pd.DataFrame) -> pd.DataFrame:
    """Add a human-readable gene_name column from the Protein field in mutations."""
    name_map = (
        mutations_df[["KeggId", "Protein"]]
        .drop_duplicates()
        .rename(columns={"KeggId": "kegg_id", "Protein": "gene_name"})
    )
    return result_df.merge(name_map, on="kegg_id", how="left")


# ─────────────────────────────────────────────────────────────────────────────
# Entry point
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Per-gene multivariate Cox regression for KEGG pathway disruption"
    )
    parser.add_argument("--cancer", default=None,
                        help="Single cancer type (e.g. 'luad'); omit for pan-cancer")
    parser.add_argument("--pathways", nargs="+",
                        default=["hsa04080", "hsa04722", "hsa04216", "hsa04137"])
    parser.add_argument("--min-mutated", type=int, default=10,
                        help="Min patients mutated in a gene to include it (default: 10)")
    parser.add_argument("--penalizer", type=float, default=0.1,
                        help="L2 penalizer for CoxPH (default: 0.1)")
    args = parser.parse_args()

    all_cancer_csvs = sorted(glob.glob(pjoin(CANCER_PATIENT_SURVIVAL_P, "*.csv")))
    all_cancers = [os.path.basename(p).replace(".csv", "") for p in all_cancer_csvs]

    if args.cancer:
        cancer_types = [args.cancer]
        label = args.cancer
    else:
        cancer_types = [c for c in all_cancers if c != "pan_cancer"]
        label = "pan_cancer"

    print(f"Cancer scope : {label}  ({len(cancer_types)} types)")
    print(f"Pathways     : {args.pathways}")
    print(f"min-mutated  : {args.min_mutated}")
    print(f"penalizer    : {args.penalizer}")

    with open(KEGG_PATHWAY_METADATA_FILE, "rb") as f:
        pathway_meta = pickle.load(f)

    print("\nLoading clinical OS data…")
    clinical_df = load_clinical(cancer_types)
    print(f"  {len(clinical_df)} unique patients")

    out_dir = pjoin(KAPLAN_MEIER_P, "gene_level_cox", label)
    os.makedirs(out_dir, exist_ok=True)

    all_results = []

    for pathway_id in args.pathways:
        if pathway_id not in pathway_meta:
            print(f"\n[warn] {pathway_id} not in KEGG metadata, skipping.")
            continue

        pathway_name = pathway_meta[pathway_id].get("name", "?").split(" - Homo")[0]
        all_gene_ids = set(pathway_meta[pathway_id].get("genes_ids", []))
        print(f"\n── {pathway_id}: {pathway_name} ({len(all_gene_ids)} pathway genes) ──")

        # Load mutations for all pathway genes first to build the name→KeggId map
        print("  Loading mutations for all pathway genes…")
        mutations_df = load_mutations(cancer_types, all_gene_ids)
        if mutations_df.empty:
            print("  No mutations found, skipping.")
            continue

        # Build protein-name → KeggId mapping from observed mutations
        protein_to_kegg = (
            mutations_df[["Protein", "KeggId"]]
            .drop_duplicates()
            .set_index("Protein")["KeggId"]
            .to_dict()
        )

        # Restrict to bottleneck genes
        bottleneck_ids = load_bottleneck_kegg_ids(pathway_id, protein_to_kegg)
        if bottleneck_ids:
            mutations_df = mutations_df[mutations_df["KeggId"].isin(bottleneck_ids)]
        else:
            print("  Falling back to all pathway genes.")

        print(f"  {len(mutations_df)} mutation rows, "
              f"{mutations_df['PatientId'].nunique()} mutated patients "
              f"in {mutations_df['KeggId'].nunique()} bottleneck genes")

        print(f"  Building gene matrix (min_mutated={args.min_mutated})…")
        matrix_df, gene_cols = build_gene_matrix(
            clinical_df, mutations_df, args.min_mutated
        )
        if not gene_cols:
            print("  No genes passed the mutation frequency filter, skipping.")
            continue

        print(f"  Fitting multivariate Cox ({len(gene_cols)} covariates)…")
        result = run_cox(matrix_df, gene_cols, pathway_id, args.penalizer)
        if result is None:
            continue

        reject, q_values = fdrcorrection(result["p"])
        result["q"] = q_values
        result["pathway_name"] = pathway_name
        result = attach_gene_names(result, mutations_df)
        result = result.sort_values("p").reset_index(drop=True)

        out_path = pjoin(out_dir, f"{pathway_id}_gene_cox.csv")
        result.to_csv(out_path, index=False)
        print(f"  Saved {len(result)} gene results → {out_path}")

        sig = result[result["q"] < 0.05]
        print(f"  Significant genes (q < 0.05): {len(sig)}")
        if not sig.empty:
            cols = ["gene_name", "kegg_id", "hr", "hr_lower95", "hr_upper95", "p", "q", "n_mutated"]
            print(sig[cols].to_string(index=False))

        all_results.append(result)

    if not all_results:
        print("\nNo results produced.")
        return

    combined = pd.concat(all_results, ignore_index=True)
    combined_path = pjoin(out_dir, "all_pathways_gene_cox.csv")
    combined.to_csv(combined_path, index=False)
    print(f"\nAll results → {combined_path}")

    print("\n=== Top genes by p-value (q < 0.20) ===")
    top = combined[combined["q"] < 0.20].sort_values("p")
    if top.empty:
        print("  None at q<0.20.  Top 20 by raw p:")
        top = combined.sort_values("p").head(20)
    print(top[["pathway_id", "gene_name", "hr", "p", "q", "n_mutated"]].to_string(index=False))


if __name__ == "__main__":
    main()
