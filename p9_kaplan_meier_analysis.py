from lifelines import KaplanMeierFitter
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


def run_survival_analysis(cancer_patients_file: str, percentile: int, output_path: str, min_patients: int = 10) -> pd.DataFrame:
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
        pathway_df = df[['Metastatic', 'OS_MONTHS', 'OS_STATUS', pathway]].dropna()

        # Separate to patient with metastases and without metastases
        metastatic_df = pathway_df[pathway_df['Metastatic'] == 1]
        non_metastatic_df = pathway_df[pathway_df['Metastatic'] == 0]

        # skip pathway if too few patients have mutations in it
        if len(metastatic_df) < min_patients:
            continue

        if len(non_metastatic_df) < min_patients:
            continue

        # log-rank test for metastatic patients
        metastatic_result = calculate_patient_scores(metastatic_df, pathway, percentile)
        if metastatic_result is not None:
            metastatic_result['Metastatic'] = 1
            results.append(metastatic_result)

        # log-rank test for non-metastatic patients
        non_metastatic_result = calculate_patient_scores(non_metastatic_df, pathway, percentile)
        if non_metastatic_result is not None:
            non_metastatic_result['Metastatic'] = 0
            results.append(non_metastatic_result)

    results_df = pd.concat(results, ignore_index=True) if results else pd.DataFrame()

    if results_df.empty:
        print(f"No pathways passed the minimum patient threshold for {cancer_type}")
        return results_df

    # multiple testing correction across all pathways
    try:
        results_df['q_value'] = fdrcorrection(results_df['p_value'])[1]
    except Exception as e:
        print(f"    ERROR during FDR correction: {e}")
        return

    results_df.sort_values('q_value', inplace=True)

    # save summary table
    results_df.to_csv(output_path, index=False)

    return results_df


def calculate_patient_scores(pathway_df: pd.DataFrame, pathway: str, percentile: int, min_patients: int = 10) -> pd.DataFrame:

    # stratify by percentile
    percentile_score = pathway_df[pathway].quantile(percentile / 100)
    high_group = pathway_df[pathway_df[pathway] > percentile_score]
    low_group = pathway_df[pathway_df[pathway] <= percentile_score]

    # Optional but recommended: skip if the mutated group is too small
    if len(high_group) < min_patients or len(low_group) < min_patients:
        return None

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

def plot_km_curve(cancer_type: str, pathway: str, results_df: pd.DataFrame, percentile: int):
    """
    Plot KM curves for a single pathway, separated by metastatic status.
    """
    # Get pathway description for title
    with open(KEGG_PATHWAY_METADATA_FILE, 'rb') as f:
        pathway_metadata = pickle.load(f)
    if "cluster" in pathway:
        pathway_cluster_annotations = pd.read_csv(pjoin(RESULTS_P, 'p12_pathway_cluster_annotations.csv'))
        cluster_number = pathway.split('_')[-1]
        pathway_description = pathway_cluster_annotations[pathway_cluster_annotations['cluster'] == int(cluster_number)]['short_name'].values[0]
    else:
        pathway_description = pathway_metadata.get(str(pathway), {}).get('name', 'Unknown').split(' - Homo')[0]

    df = pd.read_csv(pjoin(CANCER_PATIENT_SURVIVAL_P, f"{cancer_type}.csv"))
    
    # look up q-value for title
    q_val = results_df.loc[results_df['pathway'] == pathway, 'q_value'].values[0]

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))  # side-by-side plots
    kmf = KaplanMeierFitter()
    
    for idx, metastatic_status in enumerate([0, 1]):
        status_label = "Non-Metastatic" if metastatic_status == 0 else "Metastatic"
        
        # Filter by metastatic status
        pathway_df = df[(df['Metastatic'] == metastatic_status)][['OS_MONTHS', 'OS_STATUS', pathway]].dropna()
        
        percentile_score = pathway_df[pathway].quantile(percentile / 100)
        high_group = pathway_df[pathway_df[pathway] > percentile_score]
        low_group = pathway_df[pathway_df[pathway] <= percentile_score]
        
        ax = axes[idx]
        
        kmf.fit(high_group['OS_MONTHS'], event_observed=high_group['OS_STATUS'], 
                label=f'Low burden (n={len(high_group)})')
        kmf.plot_survival_function(ax=ax, ci_show=True, color=COLOR_MAP['benign'])
        
        kmf.fit(low_group['OS_MONTHS'], event_observed=low_group['OS_STATUS'], 
                label=f'High burden (n={len(low_group)})')
        kmf.plot_survival_function(ax=ax, ci_show=True, color=COLOR_MAP['pathogenic'])
        
        ax.set_title(f'{status_label} Patients')
        ax.set_xlabel('Time (months)')
        ax.set_ylabel('Survival probability')
        ax.legend()

    fig.suptitle(f'{cancer_type} — {pathway}: {pathway_description}\nq = {q_val:.3e} (percentile cutoff: {percentile}th)', 
                 fontsize=14)
    plt.tight_layout()
    
    output_path = pjoin(KAPLAN_MEIER_P, f"{cancer_type}/{pathway}_{percentile}_percentile_km.png")
    try:
        plt.savefig(output_path, dpi=150)
        plt.close()
    except Exception as e:
        print(f"    ERROR saving plot for {cancer_type} — {pathway}: {e}")

def plot_top_pathways(results_df: pd.DataFrame, cancer_type: str, percentile: int, top_n: int = 10,
                      q_threshold: float = 0.05):
    """Plots KM curves for the top_n most significant pathways passing the q_threshold."""
    significant = results_df[results_df['q_value'] < q_threshold].head(top_n)

    if significant.empty:
        print(f"No significant pathways found for {cancer_type} at q < {q_threshold}")
        return

    for _, row in significant.iterrows():
        plot_km_curve(cancer_type, row['pathway'], results_df, percentile)

    print(f"Saved {len(significant)} KM plots for {cancer_type}")

if __name__ == '__main__':
    # Get the cancer patients file based on SLURM_ARRAY_TASK_ID
    args = sys.argv[1:]
    if len(args) < 1:
        print("Usage: python -u p9_kaplan_meier_analysis.py '$SLURM_ARRAY_TASK_ID'")
        sys.exit(1)

    cancer_results_files = glob.glob(pjoin(CANCER_PATIENT_SURVIVAL_P, "*.csv"))  # should be 45 files
    if not cancer_results_files:
        print(f"No cancer patients CSV files found in the specified directory: {CANCER_PATIENT_SURVIVAL_P}")
        sys.exit(1)

    index = int(args[0])
    if index < 0 or index >= len(cancer_results_files):
        print(f"Index {index} is out of range. There are only {len(cancer_results_files)} cancer mutation files.")
        sys.exit(1)

    cancer_mutations_file = sorted(cancer_results_files)[index]

    print(f"----- Running Kaplan Meier analysis for cancer {cancer_mutations_file} -----")
    for percentile in range(10):  # you can adjust these percentiles as needed
        print(f"    Analyzing pathways at {percentile}th percentile...")
        output_dir = pjoin(KAPLAN_MEIER_P, f"{cancer_type}")
        os.makedirs(output_dir, exist_ok=True)  # Won't error if it already exists
        output_path = pjoin(output_dir, f"{percentile}_percentile.csv")

        if os.path.exists(output_path):
            print(f"    Results already exist for {cancer_type} at {percentile}th percentile, skipping csv creation.")
            results_df = pd.read_csv(output_path)
        else:
            results_df = run_survival_analysis(cancer_mutations_file, percentile, output_path)

        if not results_df.empty:
            cancer_type = os.path.basename(cancer_mutations_file).split('.')[0]  # get cancer type from file name

            print(f"    Plotting top pathways for {cancer_type}...")
            plot_top_pathways(results_df, cancer_type, percentile)
