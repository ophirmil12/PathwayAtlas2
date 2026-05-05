from lifelines import KaplanMeierFitter
from lifelines import CoxPHFitter
from lifelines.statistics import logrank_test
from statsmodels.stats.multitest import fdrcorrection
from scipy.stats import false_discovery_control
import matplotlib.pyplot as plt
import os
from os.path import join as pjoin
import pandas as pd
from definitions import *
import sys
import glob
import pickle
import numpy as np

def run_survival_analysis(cancer_patients_file: str, output_path: str, min_patients: int = 30):
    """
    For a given cancer type, run Kaplan-Meier survival analysis for each pathway.
    Stratifies patients into high/low burden groups by median score.
    Returns a DataFrame of results sorted by q-value.
    """
    # load the CSV we built in the previous step
    df = pd.read_csv(cancer_patients_file)
    cancer_type = os.path.basename(cancer_patients_file).replace('.csv', '')

    # identify pathway columns (everything that's not patient/survival metadata)
    non_pathway_cols = {'PatientId', 'StudyId', 'Metastatic', 'OS_MONTHS', 'OS_STATUS', 'DFS_MONTHS', 'DFS_STATUS'}  # adjust to your column names
    pathway_cols = [c for c in df.columns if c not in non_pathway_cols]

    results = []

    for pathway in pathway_cols:

        for percentile in range(25, 75 + 1):
            pathway_df = df[['Metastatic', 'OS_MONTHS', 'OS_STATUS', pathway]].dropna()

            # Filter for non-metastatic patients only
            non_metastatic_df = pathway_df[pathway_df['Metastatic'] == 0]

            # skip pathway if too few patients have mutations in it
            if len(non_metastatic_df) < min_patients:
                continue

            # cox test for non-metastatic patients
            non_metastatic_result = calculate_patient_scores(non_metastatic_df, pathway, percentile)
            if non_metastatic_result is not None:
                results.append(non_metastatic_result)

    results_df = pd.concat(results, ignore_index=True) if results else pd.DataFrame()

    if results_df.empty or results_df.empty:
        print(f"No pathways passed the minimum patient threshold for {cancer_type}")
        return pd.DataFrame()

    # multiple testing correction across all pathways and all percentiles
    try:
        results_df['q_value'] = fdrcorrection(results_df['p_value'])[1]
    except Exception as e:
        print(f"    ERROR during FDR correction: {e}")
        return


    # save summary table
    output_dir = os.path.dirname(output_path)
    os.makedirs(output_dir, exist_ok=True)  # Won't error if it already exists
    results_df.to_csv(output_path, index=False)

    return results_df

def calculate_patient_scores(pathway_df: pd.DataFrame, pathway: str, percentile: int, min_patients: int = 10) -> pd.DataFrame:

    # stratify by percentile
    percentile_score = pathway_df[pathway].quantile(percentile / 100)
    high_group = pathway_df[pathway_df[pathway] > percentile_score]
    low_group = pathway_df[pathway_df[pathway] <= percentile_score]

    # log-rank test
    result = logrank_test(
        high_group['OS_MONTHS'], low_group['OS_MONTHS'],
        event_observed_A=high_group['OS_STATUS'],
        event_observed_B=low_group['OS_STATUS']
    )

    return pd.DataFrame({
        'pathway': pathway,
        'n_total': len(pathway_df),
        'n_high': len(high_group),
        'n_low': len(low_group),
        'percentile_score': percentile_score,
        'p_value': result.p_value,
        'test_statistic': result.test_statistic,
    }, index=[0])

def select_best_cutoff(results_df: pd.DataFrame, pathway: str, tol=1e-4):
    pathway_results = results_df[results_df['pathway'] == pathway]
    # best p-value
    best_p = pathway_results['p_value'].min()

    # keep near-best p-values
    candidates = pathway_results[
        pathway_results['p_value'] <= best_p + tol
    ].copy()

    # compute effect size
    candidates['effect_size'] = np.abs(
        np.log(candidates['hazard_ratio'])
    )

    # pick strongest effect
    best_row = candidates.loc[
        candidates['effect_size'].idxmax()
    ]

    print(f"    Selected best cutoff for pathway {pathway} among {len(candidates)} candidates.")

    return best_row

def plot_significant_pathways(results_df: pd.DataFrame, cancer_type: str, 
                      q_threshold: float = 0.05):
    """Plots KM curves for the top_n most significant pathways passing the q_threshold."""
    # for each pathway, get the best cutoff across percentiles
    results_df = results_df[results_df['q_value'] <= q_threshold]
    if results_df.empty:
        print(f"    No pathways passed the q-value threshold of {q_threshold} for {cancer_type}.")
        return

    for pathway in results_df['pathway'].unique():
        best_row = select_best_cutoff(results_df, pathway)
        print(f"    Plotting pathway {pathway} with q-value {best_row['q_value']:.2e}...")
        plot_km_curve(best_row, cancer_type, pathway)



def plot_km_curve(best_row: pd.Series, cancer_type: str, pathway: str):
    """
    Plot KM curves for a single pathway using the selected percentile cutoff.
    """

    # Load pathway metadata
    with open(KEGG_PATHWAY_METADATA_FILE, 'rb') as f:
        pathway_metadata = pickle.load(f)

    pathway_description = pathway_metadata.get(
        str(pathway), {}
    ).get('name', 'Unknown').split(' - Homo')[0]

    # Load survival data
    survival_df = pd.read_csv(
        pjoin(CANCER_PATIENT_SURVIVAL_P, f"{cancer_type}.csv")
    )

    # Filter for non-metastatic patients
    pathway_df = survival_df[
        survival_df['Metastatic'] == 0
    ][['OS_MONTHS', 'OS_STATUS', pathway]].dropna()

    # Get percentile from best_row
    percentile = best_row['percentile']
    percentile_score = pathway_df[pathway].quantile(percentile / 100)
    # Split groups
    high_group = pathway_df[pathway_df[pathway] > percentile_score]
    low_group = pathway_df[pathway_df[pathway] <= percentile_score]

    # Initialize plot
    fig, ax = plt.subplots(figsize=(8, 6))
    kmf = KaplanMeierFitter()

    # Plot high group
    kmf.fit(
        durations=high_group['OS_MONTHS'],
        event_observed=high_group['OS_STATUS'],
        label=f'Low burden (n={len(high_group)})'
    )
    kmf.plot(ax=ax, color=COLOR_MAP['benign'])

    # Plot LOW group
    kmf.fit(
        durations=low_group['OS_MONTHS'],
        event_observed=low_group['OS_STATUS'],
        label=f'High burden (n={len(low_group)})'
    )
    kmf.plot(ax=ax, color=COLOR_MAP['pathogenic'])

    # Extract stats
    p_value = best_row.get('p_value', None)
    q_value = best_row.get('q_value', None)
    hr = best_row.get('hazard_ratio', None)

    # Title
    title = f"{cancer_type} - {pathway}: {pathway_description}\n"

    if q_value is not None:
        title += f"P = {q_value:.2e}"
    if hr is not None:
        title += f", HR = {hr:.2f}"

    ax.set_title(title)

    # Labels
    ax.set_xlabel("Time (months)")
    ax.set_ylabel("Survival probability")

    plt.tight_layout()
    plt.savefig(pjoin(KAPLAN_MEIER_P, f"{cancer_type}/{pathway}_km_curve.png"))
    plt.close()

if __name__ == '__main__':
    # Get the cancer patients file based on SLURM_ARRAY_TASK_ID
    args = sys.argv[1:]
    if len(args) < 1:
        print("Usage: python -u p9_kaplan_meier_analysis.py '$SLURM_ARRAY_TASK_ID'")
        sys.exit(1)

    cancer_patients_files = glob.glob(pjoin(CANCER_PATIENT_SURVIVAL_P, "*.csv"))  # should be 45 files
    if not cancer_patients_files:
        print(f"No cancer patients CSV files found in the specified directory: {CANCER_PATIENT_SURVIVAL_P}")
        sys.exit(1)

    index = int(args[0])
    if index < 0 or index >= len(cancer_patients_files):
        print(f"Index {index} is out of range. There are only {len(cancer_patients_files)} cancer mutation files.")
        sys.exit(1)

    cancer_patients_file = sorted(cancer_patients_files)[index]
    cancer_type = os.path.basename(cancer_patients_file).replace('.csv', '')  # get cancer type from file name

    print(f"----- Running Kaplan Meier analysis for cancer {cancer_patients_file} -----")

    output_path = pjoin(KAPLAN_MEIER_P, f"{cancer_type}/cox_results.csv")
    if os.path.exists(output_path):
        print(f"    Results already exist for {cancer_type}, skipping csv creation.")
        results_df = pd.read_csv(output_path)
    else:
        results_df = run_survival_analysis(cancer_patients_file, output_path)

    if not results_df.empty:
        print(f"    Plotting all significant pathways for {cancer_type}...")
        plot_significant_pathways(results_df, cancer_type)
    else:
        print(f"    No significant pathways to plot for {cancer_type}.")
