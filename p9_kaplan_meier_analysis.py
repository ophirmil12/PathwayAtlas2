from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from statsmodels.stats.multitest import fdrcorrection
import matplotlib.pyplot as plt
import os


def run_survival_analysis(cancer_patient_file: str, min_patients: int = 20, output_dir: str = None):
    """
    For a given cancer type, run Kaplan-Meier survival analysis for each pathway.
    Stratifies patients into high/low burden groups by median score.
    Returns a DataFrame of results sorted by q-value.
    """
    # load the CSV we built in the previous step
    df = pd.read_csv(cancer_patient_file)

    # identify pathway columns (everything that's not patient/survival metadata)
    non_pathway_cols = {'PatientId', 'OS_MONTHS', 'OS_STATUS'}  # adjust to your column names
    pathway_cols = [c for c in df.columns if c not in non_pathway_cols]

    results = []

    for pathway in pathway_cols:
        # only keep patients who have a score for this pathway (i.e. had mutations in it)
        pathway_df = df[['OS_MONTHS', 'OS_STATUS', pathway]].dropna()

        # skip pathway if too few patients have mutations in it
        if len(pathway_df) < min_patients:
            continue

        # stratify by median
        median_score = pathway_df[pathway].median()
        high_group = pathway_df[pathway_df[pathway] >= median_score]
        low_group = pathway_df[pathway_df[pathway] < median_score]

        # log-rank test
        result = logrank_test(
            high_group['OS_MONTHS'], low_group['OS_MONTHS'],
            event_observed_A=high_group['OS_STATUS'],
            event_observed_B=low_group['OS_STATUS']
        )

        results.append({
            'pathway': pathway,
            'pathway_name': pathway_id_to_metadata[pathway]['name'],  # adjust to your metadata structure
            'n_total': len(pathway_df),
            'n_high': len(high_group),
            'n_low': len(low_group),
            'median_score': median_score,
            'p_value': result.p_value,
            'test_statistic': result.test_statistic,
        })

    results_df = pd.DataFrame(results)

    if results_df.empty:
        print(f"No pathways passed the minimum patient threshold for {cancer_type}")
        return results_df

    # multiple testing correction across all pathways
    results_df['q_value'] = fdrcorrection(results_df['p_value'])[1]
    results_df.sort_values('q_value', inplace=True)

    # save summary table
    results_df.to_csv(pjoin(output_dir, f"{cancer_type}_survival_results.csv"), index=False)

    return results_df


def plot_km_curve(cancer_type: str, pathway: str, results_df: pd.DataFrame, output_dir: str):
    """
    Plot and save a Kaplan-Meier curve for a single pathway in a single cancer type.
    Call this for your top significant hits.
    """
    df = pd.read_csv(pjoin(CANCER_PATIENT_SURVIVAL_P, f"{cancer_type}.csv"))
    pathway_df = df[['OS_MONTHS', 'OS_STATUS', pathway]].dropna()

    median_score = pathway_df[pathway].median()
    high_group = pathway_df[pathway_df[pathway] >= median_score]
    low_group = pathway_df[pathway_df[pathway] < median_score]

    # look up q-value for title
    q_val = results_df.loc[results_df['pathway'] == pathway, 'q_value'].values[0]
    pathway_name = pathway_id_to_metadata[pathway]['name']

    fig, ax = plt.subplots(figsize=(8, 5))
    kmf = KaplanMeierFitter()

    kmf.fit(high_group['OS_MONTHS'], event_observed=high_group['OS_STATUS'], label=f'High burden (n={len(high_group)})')
    kmf.plot_survival_function(ax=ax, ci_show=True)

    kmf.fit(low_group['OS_MONTHS'], event_observed=low_group['OS_STATUS'], label=f'Low burden (n={len(low_group)})')
    kmf.plot_survival_function(ax=ax, ci_show=True)

    ax.set_title(f'{cancer_type} — {pathway_name}\nq = {q_val:.3e}')
    ax.set_xlabel('Time (months)')
    ax.set_ylabel('Survival probability')
    ax.legend()

    plt.tight_layout()
    plt.savefig(pjoin(output_dir, f"{cancer_type}_{pathway}_km.png"), dpi=150)
    plt.close()


def plot_top_pathways(cancer_type: str, results_df: pd.DataFrame, output_dir: str, top_n: int = 10,
                      q_threshold: float = 0.05):
    """Plots KM curves for the top_n most significant pathways passing the q_threshold."""
    significant = results_df[results_df['q_value'] < q_threshold].head(top_n)

    if significant.empty:
        print(f"No significant pathways found for {cancer_type} at q < {q_threshold}")
        return

    for _, row in significant.iterrows():
        plot_km_curve(cancer_type, row['pathway'], results_df, output_dir)

    print(f"Saved {len(significant)} KM plots for {cancer_type}")