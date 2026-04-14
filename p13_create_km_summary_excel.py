# This files checks for each novel hit in the novel hits excel whether theres a km plot showing significant survival effects
import argparse
import glob
import os
from os.path import join as pjoin
import pandas as pd
from definitions import *
import numpy as np
from openpyxl import Workbook
from openpyxl.drawing.image import Image
import pickle

def create_km_summary_excel(summary_type: str = "all"):
    """
    Creates an excel file summarizing the Kaplan-Meier analysis results for each cancer type and pathway.
    For each pathway with a KM plot, it extracts the q-values for metastatic and non-metastatic patients, the pathway description, and includes the KM plot in the excel.
    Saves the excel file in the KAPLAN_MEIER_P directory.
    """
    wb = Workbook()
    ws = wb.active

    base_cols = ['cancer_type', 'pathway_id', 'pathway_name', 'delta_means']
    if summary_type == "novel":
        novel_hits_df = pd.read_excel(pjoin(RESULTS_P, "novel_pathway_hits.xlsx"))
        novel_hits_benign_df = pd.read_excel(pjoin(RESULTS_P, "novel_pathway_hits_benign.xlsx"))
        novel_pathways = set(novel_hits_df['pathway_id'])
        novel_benign_pathways = set(novel_hits_benign_df['pathway_id'])
        ws.append(base_cols + ['q_value', 'confidence_score', 'biological_relevance', 'suggested_mechanism', 'KM_plot'])
        km_plot_col = 'I'
    else:
        novel_hits_df = None
        novel_pathways = set()
        ws.append(base_cols + ['KM_plot'])
        km_plot_col = 'E'

    km_cancer_directories = sorted(glob.glob(pjoin(KAPLAN_MEIER_P, '*')))

    for km_cancer_directory in km_cancer_directories:
        cancer_type = km_cancer_directory.split('/')[-1]

        km_plots = sorted(glob.glob(pjoin(km_cancer_directory, '*.png')))

        if not km_plots:
            print(f"    No km plots were generated for {cancer_type}, skipping...")
            continue

        # for each km plot, add to the excel a row with cancer type, pathway id and full name, q_value_met, q_value_non_met, percentile, and km plot
        for km_plot in km_plots:
            pathway_id = os.path.basename(km_plot).split('_')[0]  # assuming pathway_id is the first part of the filename before the first underscore
            pathway_name = get_pathway_description(pathway_id)
            # get delta means from distances result csv of that cancer in the pathway_id row
            distances_csv_path = pjoin(RESULTS_DISTANCES_P, f"{cancer_type}.csv")
            distances_df = pd.read_csv(distances_csv_path)
            delta_means_series = distances_df[distances_df['pathway'] == pathway_id]['delta_means']
            delta_means = delta_means_series.values[0] if not delta_means_series.empty else np.nan

            if summary_type == "novel":
                if pathway_id not in novel_pathways and pathway_id not in novel_benign_pathways:
                    continue
                else:
                    pathway_cancer_row = novel_hits_df[(novel_hits_df['pathway_id'] == pathway_id) & (novel_hits_df['cancer_type'] == cancer_type)]
                    if pathway_cancer_row.empty:
                        pathway_cancer_row = novel_hits_benign_df[(novel_hits_benign_df['pathway_id'] == pathway_id) & (novel_hits_benign_df['cancer_type'] == cancer_type)]
                    if not pathway_cancer_row.empty:
                        q_value = pathway_cancer_row['q_value'].values[0]
                        print(f"    Adding {pathway_id} for {cancer_type} with q_value {q_value} to the summary excel")
                        confidence_score = pathway_cancer_row['confidence_score'].values[0]
                        relevance = pathway_cancer_row['biological_relevance'].values[0]
                        mechanism = pathway_cancer_row['suggested_mechanism'].values[0]
                    else:
                        q_value = np.nan
                        confidence_score = np.nan
                        relevance = np.nan
                        mechanism = np.nan
                    ws.append([cancer_type, pathway_id, pathway_name, delta_means, f"{q_value:.3f}", f"{confidence_score:.3f}", relevance, mechanism])
            else:
                ws.append([cancer_type, pathway_id, pathway_name, delta_means])

            ws.row_dimensions[ws.max_row].height = 225  # Set height for image
            ws.column_dimensions[km_plot_col].width = 53  # Set width for image column
            img = Image(km_plot)
            img.width = 400
            img.height = 300
            ws.add_image(img, f'{km_plot_col}{ws.max_row}')

    # Save excel at KAPLAN_MEIER_P
    wb.save(pjoin(KAPLAN_MEIER_P, f'km_{summary_type}_pathways_summary.xlsx'))


def get_pathway_description(pathway_id: str) -> str:
    with open(KEGG_PATHWAY_METADATA_FILE, 'rb') as f:
        pathway_metadata = pickle.load(f)
    if "cluster" in pathway_id:
        pathway_cluster_annotations = pd.read_csv(pjoin(RESULTS_P, 'p12_pathway_cluster_annotations.csv'))
        cluster_number = pathway_id.split('_')[-1]
        return pathway_cluster_annotations[pathway_cluster_annotations['cluster'] == int(cluster_number)]['short_name'].values[0]
    else:
        return pathway_metadata.get(str(pathway_id), {}).get('name', 'Unknown')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Summarize Kaplan-Meier analysis results in Excel")

    parser.add_argument(
        "--summary_type",
        type=str,
        choices=["all", "novel"],
        help="Summary type: all significant pathways or novel pathways"
    )
    args = parser.parse_args()
    summary_type = args.summary_type
    print(f"----- Summarizing KM results for {summary_type} pathways -----")

    create_km_summary_excel(summary_type)
    