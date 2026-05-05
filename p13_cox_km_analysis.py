"""
p13_cox_km_analysis.py  –  Three-group pathway burden survival analysis
=======================================================================

Design rationale
----------------
Patients are stratified into THREE groups per pathway, not two:

  0 – Wildtype   : no mutation in the pathway at all
  1 – Low burden : mutated, but score above the median of mutated patients
                   (i.e. less pathogenic by ESM / lower path_count)
  2 – High burden: mutated, score at or below the median of mutated patients
                   (i.e. more pathogenic)

Why three groups?
  • Patients with zero pathway mutations are biologically distinct from
    patients with a neutral missense variant (raw ESM score ≈ 0).  Collapsing
    them into a single "low" bin mixes wildtype with benign-mutant and dilutes
    any true signal.
  • Including wildtype patients as the reference group massively increases the
    effective sample size (all patients, not just the mutated minority), which
    is the main reason per-pathway analyses were failing the minimum-events
    filter.
  • This design is standard in TCGA mutation-survival studies that distinguish
    mutant vs. wildtype before further stratifying by severity.

Cox model setup
  • The three-group ordinal variable (0/1/2) is passed as a single continuous
    covariate.  This imposes the constraint that each step up carries the same
    log-HR increment, which is appropriate here because the ordering is
    biologically meaningful (wildtype < low burden < high burden).
  • A single coefficient is extracted, so the model has 1 degree of freedom and
    maximum power given the sample sizes typical in per-cancer TCGA studies.
  • Penalizer=0.1 is kept for numerical stability.

Minimum events filter
  • Changed from (counts < 5).all()  →  (counts < 10).any()
    The old logic skipped only if BOTH groups were event-sparse, meaning one
    arm could have just 1 event and the model would still run, producing
    unstable estimates.  The new threshold (10 per group) follows the standard
    "10 events per variable" rule of thumb for Cox models.
"""

from lifelines import KaplanMeierFitter
from lifelines import CoxPHFitter
from lifelines.statistics import multivariate_logrank_test
from statsmodels.stats.multitest import fdrcorrection
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os
from os.path import join as pjoin
import pandas as pd
from definitions import *
import sys
import glob
import pickle
import numpy as np
import warnings
from lifelines.exceptions import ConvergenceWarning
import scipy.stats

# ─────────────────────────────────────────────────────────────────────────────
# Group label constants  (used in both Cox and KM functions)
# ─────────────────────────────────────────────────────────────────────────────
WILDTYPE    = 0   # no mutation in this pathway
LOW_BURDEN  = 1   # mutated, score above median among mutated (less pathogenic)
HIGH_BURDEN = 2   # mutated, score at/below median among mutated (more pathogenic)

GROUP_COLORS = {
    WILDTYPE:    COLOR_MAP["benign"],   # blue
    LOW_BURDEN:  COLOR_MAP["non-significant"],   # orange
    HIGH_BURDEN: COLOR_MAP["pathogenic"],   # red
}
GROUP_LABELS = {
    WILDTYPE:    'Wildtype (no mutation)',
    LOW_BURDEN:  'Low burden',
    HIGH_BURDEN: 'High burden',
}


# ─────────────────────────────────────────────────────────────────────────────
# Helper: assign three-group labels for one pathway
# ─────────────────────────────────────────────────────────────────────────────
def assign_burden_groups(df: pd.DataFrame, scoring_type: str) -> pd.DataFrame:
    """
    Adds a 'burden_group' column to df.

    Rows where scoring_type is NaN are patients with no pathway mutation
    → wildtype (group 0).

    Among patients WITH a mutation (non-NaN score), the median of their scores
    is used to split them into low (group 1) and high (group 2) burden.

    Parameters
    ----------
    df           : DataFrame that must contain the scoring_type column.
    scoring_type : 'min_esm_log_probs' or 'path_count'.  For ESM scores,
                   caller should have already negated the column so that
                   higher values = more pathogenic.

    Returns
    -------
    df with a new integer column 'burden_group' (0 / 1 / 2).
    """
    df = df.copy()
    df['burden_group'] = WILDTYPE  # default: wildtype

    mutated_mask = df[scoring_type].notna()
    mutated_scores = df.loc[mutated_mask, scoring_type]

    if mutated_scores.empty:
        return df  # no mutated patients in this pathway for this subset

    median_score = mutated_scores.median()

    # High burden: score > median  (more pathogenic because we negated ESM)
    df.loc[mutated_mask & (df[scoring_type] >  median_score), 'burden_group'] = HIGH_BURDEN
    # Low burden:  score <= median
    df.loc[mutated_mask & (df[scoring_type] <= median_score), 'burden_group'] = LOW_BURDEN

    return df


# ─────────────────────────────────────────────────────────────────────────────
# Main analysis loop
# ─────────────────────────────────────────────────────────────────────────────
def run_survival_analysis(
    cancer_patients_file: str,
    output_path: str,
    min_patients: int = 60,
    min_events_per_group: int = 10,
    scoring_type: str = "min_esm_log_probs",
) -> pd.DataFrame:
    """
    Run three-group Cox survival analysis for every pathway in a cancer type.

    For each pathway:
      1. All patients (including those with no pathway mutation) are assigned
         to one of three burden groups.
      2. Non-metastatic patients are selected.
      3. The combined group size must reach min_patients.
      4. Every group must have at least min_events_per_group observed deaths.
      5. A Cox model is fit with burden_group as an ordinal covariate (0/1/2).

    Returns a DataFrame of per-pathway Cox results, FDR-corrected.
    """
    df = pd.read_csv(cancer_patients_file)
    cancer_type = os.path.basename(cancer_patients_file).replace('.csv', '')

    # Negate ESM scores so that higher value = more pathogenic
    # (raw ESM: more negative = worse, so we flip the sign)
    if scoring_type == "min_esm_log_probs":
        df[scoring_type] = -df[scoring_type]

    results = []

    all_pathways = df['pathway_id'].unique()
    print(f"  Running Cox analysis for {len(all_pathways)} pathways in {cancer_type}...")

    for pathway in all_pathways:

        # ── Isolate this pathway's rows (patients WITH a mutation here) ──────
        # Patients absent from this pathway's rows are wildtype for it.
        pathway_mutated_df = df[df['pathway_id'] == pathway].copy()

        # ── Build a full-cohort frame for this pathway ───────────────────────
        # Start from all unique patients in the file (any row gives OS data).
        all_patients = (
            df[['PatientId', 'StudyId', 'Metastatic', 'OS_MONTHS', 'OS_STATUS']]
            .drop_duplicates(subset=['PatientId', 'StudyId'])
        )

        # Merge pathway score onto all patients (NaN for wildtype)
        pathway_scores = pathway_mutated_df[
            ['PatientId', 'StudyId', scoring_type]
        ].drop_duplicates(subset=['PatientId', 'StudyId'])

        full_df = all_patients.merge(
            pathway_scores, on=['PatientId', 'StudyId'], how='left'
        )

        # ── Assign three burden groups ────────────────────────────────────────
        full_df = assign_burden_groups(full_df, scoring_type)

        # ── Restrict to non-metastatic patients ──────────────────────────────
        analysis_df = full_df[full_df['Metastatic'] == 0].copy()
        analysis_df = analysis_df[['OS_MONTHS', 'OS_STATUS', 'burden_group']].dropna()

        # ── Minimum total patients filter ────────────────────────────────────
        if len(analysis_df) < min_patients:
            continue

        # ── Minimum events per group filter (10 events-per-variable rule) ────
        events_per_group = analysis_df.groupby('burden_group')['OS_STATUS'].sum()
        if (events_per_group < min_events_per_group).any():
            # It's acceptable if a group simply doesn't exist (e.g. all wildtype),
            # but any group that does exist must meet the threshold.
            present_groups = events_per_group[events_per_group.index.isin(analysis_df['burden_group'].unique())]
            if (present_groups < min_events_per_group).any():
                continue

        # ── Per-group size check (before Cox) ────────────────────────────────────
        group_counts = analysis_df.groupby('burden_group').size()
        group_events = analysis_df.groupby('burden_group')['OS_STATUS'].sum()

        # Every group that is present must have enough patients AND enough events
        if (group_counts < 15).any() or (group_events < 10).any():
            continue

        # Variance guard: if burden_group is nearly constant, Cox will fail
        group_std = analysis_df['burden_group'].std()
        if group_std < 0.1:
            continue

        # ── Fit Cox model ─────────────────────────────────────────────────────
        result = _fit_cox(analysis_df, pathway)
        if result is not None:
            results.append(result)

    if not results:
        print(f"  No pathways passed filters for {cancer_type}.")
        return pd.DataFrame()

    results_df = pd.concat(results, ignore_index=True)

    # ── FDR correction across all pathways ───────────────────────────────────
    results_df['q_high_vs_wildtype'] = fdrcorrection(results_df['p_high_vs_wildtype'])[1]
    results_df.sort_values('q_high_vs_wildtype', inplace=True)

    results_df['q_high_vs_low'] = fdrcorrection(results_df['p_high_vs_low'])[1]
    results_df.sort_values('q_high_vs_low', inplace=True)

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    results_df.to_csv(output_path, index=False)
    print(f"  Saved {len(results_df)} pathway results to {output_path}")

    return results_df


def _fit_cox(analysis_df, pathway):
    try:
        df = analysis_df.copy()

        with open(KEGG_PATHWAY_METADATA_FILE, 'rb') as f:
            pathway_metadata = pickle.load(f)
        pathway_name = (
            pathway_metadata.get(str(pathway), {})
            .get('name', 'Unknown')
            .split(' - Homo')[0]
        )
        pathway_size = len(pathway_metadata['genes_ids'])
        
        # Dummy encode: wildtype (0) is the reference category
        df['is_low_burden']  = (df['burden_group'] == LOW_BURDEN).astype(int)
        df['is_high_burden'] = (df['burden_group'] == HIGH_BURDEN).astype(int)

        cph = CoxPHFitter(penalizer=0.1)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", ConvergenceWarning)
            warnings.simplefilter("ignore", RuntimeWarning)
            cph.fit(
                df[['OS_MONTHS', 'OS_STATUS', 'is_low_burden', 'is_high_burden']],
                duration_col='OS_MONTHS',
                event_col='OS_STATUS'
            )

        # Check if any coefficients are NaN — indicates failed convergence
        if cph.summary['coef'].isna().any():
            print(f"    Skipping pathway {pathway}: NaN coefficients, convergence failed.")
            return None

        if cph.summary['p'].isna().any():
            print(f"    Skipping pathway {pathway}: NaN p-values, convergence failed.")
            return None

        if np.isnan(cph.log_likelihood_):
            print(f"    Skipping pathway {pathway}: Cox did not converge.")
            return None

        # Extract the two HRs independently
        high_vs_wt = cph.summary.loc['is_high_burden']
        low_vs_wt  = cph.summary.loc['is_low_burden']

        # ── Derive high vs low from the single fit ────────────────────────────
        # HR(high vs low) = exp(coef_high - coef_low)
        # p-value uses the variance-covariance matrix to account for
        # the covariance between the two coefficient estimates
        coef_diff      = high_vs_wt['coef'] - low_vs_wt['coef']
        hr_high_vs_low = np.exp(coef_diff)

        cov            = cph.variance_matrix_
        se_diff        = np.sqrt(
            cov.loc['is_high_burden', 'is_high_burden']
            + cov.loc['is_low_burden',  'is_low_burden']
            - 2 * cov.loc['is_high_burden', 'is_low_burden']
        )
        z_high_vs_low  = coef_diff / se_diff
        p_high_vs_low  = 2 * (1 - scipy.stats.norm.cdf(abs(z_high_vs_low)))

        return pd.DataFrame({
            'pathway':             [pathway],
            'pathway_name':        [pathway_name],
            'pathway_size':        [pathway_size],
            # HRs
            'hr_high_vs_wildtype': [high_vs_wt['exp(coef)']],
            'hr_low_vs_wildtype':  [low_vs_wt['exp(coef)']],
            'hr_high_vs_low':      [hr_high_vs_low],
            # Supporting stats
            'p_high_vs_wildtype':  [high_vs_wt['p']],
            'p_low_vs_wildtype':   [low_vs_wt['p']],
            'p_high_vs_low':       [p_high_vs_low],
            'z_high_vs_wildtype':  [high_vs_wt['z']],
            'z_high_vs_low':       [z_high_vs_low],
            # Group sizes
            'n_total':             [len(df)],
            'n_wildtype':          [(df['burden_group'] == WILDTYPE).sum()],
            'n_low':               [(df['burden_group'] == LOW_BURDEN).sum()],
            'n_high':              [(df['burden_group'] == HIGH_BURDEN).sum()],
            'events_wildtype':     [df.loc[df['burden_group'] == WILDTYPE,    'OS_STATUS'].sum()],
            'events_low':          [df.loc[df['burden_group'] == LOW_BURDEN,  'OS_STATUS'].sum()],
            'events_high':         [df.loc[df['burden_group'] == HIGH_BURDEN, 'OS_STATUS'].sum()],
        })
    except Exception as e:
        print(f"    Cox fitting failed for pathway {pathway}: {e}")
        return None

# Find the time when the smallest group drops below 5 at-risk
def get_truncation_time(group_df, min_at_risk=5):
    # Sort by time, compute rolling at-risk count
    sorted_times = group_df['OS_MONTHS'].sort_values()
    if len(sorted_times) <= min_at_risk:
        return sorted_times.iloc[-1]
    return sorted_times.iloc[-min_at_risk]  # time of 5th-from-last patient

# ─────────────────────────────────────────────────────────────────────────────
# Plotting
# ─────────────────────────────────────────────────────────────────────────────
def plot_significant_pathways(
    results_df: pd.DataFrame,
    cancer_type: str,
    scoring_type: str = "min_esm_log_probs",
    q_threshold: float = 0.05,
):
    """Plot KM curves for every pathway that passes the q-value threshold."""
    sig_df = results_df[(results_df['q_high_vs_wildtype'] <= q_threshold) | (results_df['q_high_vs_low'] <= q_threshold)]
    if sig_df.empty:
        print(f"  No pathways passed q ≤ {q_threshold} for {cancer_type}.")
        return

    print(f"  Plotting {len(sig_df)} significant pathways for {cancer_type}...")
    for _, row in sig_df.iterrows():
        plot_km_curve(row, cancer_type, scoring_type)


def plot_km_curve(result_row: pd.Series, cancer_type: str, scoring_type: str):
    """
    Plot a three-group Kaplan-Meier curve for a single pathway.

    The three curves are:
      Wildtype (no pathway mutation)
      Low burden (mutated, less pathogenic half)
      High burden (mutated, more pathogenic half)

    The plot title includes the q-value, the hazard ratio per group step,
    and the per-group patient/event counts as a legend annotation.

    A multivariate log-rank p-value (testing all three curves simultaneously)
    is computed and shown separately from the Cox HR, which comes from the
    ordinal model.
    """
    pathway = result_row['pathway']

    # ── Load survival data ────────────────────────────────────────────────────
    survival_df = pd.read_csv(pjoin(CANCER_PATIENT_SURVIVAL_P, f"{cancer_type}.csv"))

    # Negate ESM if needed (mirror what run_survival_analysis does)
    if scoring_type == "min_esm_log_probs" and scoring_type in survival_df.columns:
        survival_df[scoring_type] = -survival_df[scoring_type]

    # ── Build full-cohort frame with burden groups ────────────────────────────
    all_patients = (
        survival_df[['PatientId', 'StudyId', 'Metastatic', 'OS_MONTHS', 'OS_STATUS']]
        .drop_duplicates(subset=['PatientId', 'StudyId'])
    )
    pathway_scores = (
        survival_df[survival_df['pathway_id'] == pathway]
        [['PatientId', 'StudyId', scoring_type]]
        .drop_duplicates(subset=['PatientId', 'StudyId'])
    )
    full_df = all_patients.merge(pathway_scores, on=['PatientId', 'StudyId'], how='left')
    full_df = assign_burden_groups(full_df, scoring_type)

    # Restrict to non-metastatic with complete OS data
    plot_df = full_df[full_df['Metastatic'] == 0][
        ['OS_MONTHS', 'OS_STATUS', 'burden_group']
    ].dropna()

    # ── Multivariate log-rank test (all three groups simultaneously) ──────────
    mlr = multivariate_logrank_test(
        plot_df['OS_MONTHS'],
        plot_df['burden_group'],
        plot_df['OS_STATUS'],
    )
    logrank_p = mlr.p_value

    # ── Build KM figure ───────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(9, 6))
    kmf = KaplanMeierFitter()
    legend_patches = []

    for group_id in [WILDTYPE, LOW_BURDEN, HIGH_BURDEN]:
        grp = plot_df[plot_df['burden_group'] == group_id]
        if grp.empty:
            continue

        n       = len(grp)
        events  = int(grp['OS_STATUS'].sum())
        color   = GROUP_COLORS[group_id]
        label   = GROUP_LABELS[group_id]

        kmf.fit(
            durations=grp['OS_MONTHS'],
            event_observed=grp['OS_STATUS'],
            label=f'{label} (n={n}, events={events})',
        )
        kmf.plot_survival_function(ax=ax, color=color, ci_show=True)

        legend_patches.append(
            mpatches.Patch(color=color, label=f'{label}\nn={n}, events={events}')
        )

    # ── Annotations ───────────────────────────────────────────────────────────

    pathway_name = result_row.get('pathway_name', str(result_row['pathway']))

    hr_high_vs_wt  = result_row.get('hr_high_vs_wildtype', float('nan'))
    hr_high_vs_low = result_row.get('hr_high_vs_low',      float('nan'))
    q_hvw          = result_row.get('q_high_vs_wildtype',  float('nan'))
    q_hvl          = result_row.get('q_high_vs_low',       float('nan'))

    table_data = [
        ['',              'Hazard Ratio',                   'P-Value'],
        ['High burden vs Wildtype', f'{hr_high_vs_wt:.2f}', f'{q_hvw:.2e}'],
        ['High burden vs Low',      f'{hr_high_vs_low:.2f}', f'{q_hvl:.2e}'],
    ]

    table = ax.table(
        cellText=table_data[1:],
        colLabels=table_data[0],
        cellLoc='center',
        loc='lower left',
        bbox=[0.02, 0.02, 0.50, 0.24],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8)

    # Header
    for col in range(3):
        #table[0, col].set_facecolor('#FFFFFF')
        table[0, col].set_text_props(color='black', fontweight='bold')

    # Row colors matching the KM curve colors
    row_colors = ['#FFFFFF', '#FFFFFF']
    for row, color in enumerate(row_colors, start=1):
        for col in range(2):
            table[row, col].set_facecolor(color)

    title = (
        f"{CANCER_FULLNAME.get(cancer_type, cancer_type)}  -  {pathway_name}\n"
    )
    ax.set_title(title, fontsize=10, pad=12)
    ax.set_xlabel("Time (months)")
    ax.set_ylabel("Survival probability")
    ax.legend(handles=legend_patches, loc='upper right', fontsize=8, framealpha=0.8)
    ax.set_ylim(0, 1.05)
    trunc_times = [
        get_truncation_time(plot_df[plot_df['burden_group'] == g])
            for g in [WILDTYPE, LOW_BURDEN, HIGH_BURDEN]
            if not plot_df[plot_df['burden_group'] == g].empty
    ]
    ax.set_xlim(0, min(trunc_times))

    plt.tight_layout()

    out_dir = pjoin(KAPLAN_MEIER_P, cancer_type)
    os.makedirs(out_dir, exist_ok=True)
    out_path = pjoin(out_dir, f"{pathway}_{scoring_type}_km_curve.png")
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"    Saved KM plot → {out_path}")


# ─────────────────────────────────────────────────────────────────────────────
# Utility: merge per-cancer CSVs into a pan-cancer file
# ─────────────────────────────────────────────────────────────────────────────
def merge_patient_scores_csvs(pan_cancer_patient_csv_outpath: str):
    patient_scores_csvs = sorted(glob.glob(pjoin(CANCER_PATIENT_SURVIVAL_P, "*.csv")))
    all_dfs = [pd.read_csv(f) for f in patient_scores_csvs]
    pan_cancer_df = pd.concat(all_dfs, ignore_index=True)
    pan_cancer_df.drop_duplicates(
        keep='first', inplace=True, ignore_index=True,
        subset=DUPLICATE_EXCLUSION_COLUMNS,
    )
    pan_cancer_df.to_csv(pan_cancer_patient_csv_outpath, index=False)


if __name__ == '__main__':
    args = sys.argv[1:]
    if len(args) < 1:
        print("Usage: python -u p13_cox_km_analysis.py '$SLURM_ARRAY_TASK_ID'")
        sys.exit(1)

    cancer_patients_files = sorted(glob.glob(pjoin(CANCER_PATIENT_SURVIVAL_P, "*.csv")))
    if not cancer_patients_files:
        print(f"No cancer patients CSV files found in {CANCER_PATIENT_SURVIVAL_P}")
        sys.exit(1)

    index = int(args[0])
    if index < 0 or index >= len(cancer_patients_files):
        print(f"Index {index} out of range ({len(cancer_patients_files)} files).")
        sys.exit(1)

    cancer_patients_file = cancer_patients_files[index]
    cancer_type = os.path.basename(cancer_patients_file).replace('.csv', '')

    print(f"----- Three-group Cox/KM analysis: {cancer_type} -----")

    #scoring_types = ["min_esm_log_probs", "path_count"]
    scoring_types = ["min_esm_log_probs"]
    for scoring_type in scoring_types:
        output_path = pjoin(KAPLAN_MEIER_P, f"{cancer_type}/{scoring_type}.csv")

        # if os.path.exists(output_path):
        #     print(f"  Results already exist, loading from {output_path}")
        #     results_df = pd.read_csv(output_path)
        # else:
        #     results_df = run_survival_analysis(cancer_patients_file, output_path)

        results_df = run_survival_analysis(cancer_patients_file, output_path)

        if results_df is not None and not results_df.empty:
            plot_significant_pathways(results_df, cancer_type, scoring_type)
        else:
            print(f"  No significant pathways to plot for {cancer_type} with {scoring_types}.")