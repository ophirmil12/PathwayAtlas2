# Here we want for each cancer-results-file to write for each of the pathways some
# statistics about the amount of data we used to get the result
# 1. Number of mutations (rows of the pathway in the cancer SNVs files)
# 2. Total length of the genes in the pathway
# 3. Number of genes in the pathway
# 4. "Covered" genes:
#       Number of genes in the pathway that has record (in cancer SNVs files) of 10 or more mutations
#       or at least 1 mutation per 100 \ 200 bp
# 5. pathway name (metadata)

# Read data from RESULTS_DISTANCES_P folder (for cancer name and pathway ID),
# KEGG_PATHWAY_METADATA_P (for pathway's genes and metadata ["name"],
# KeggGenes objects for "is CDS" and aa_seq and more,
# and CBIO_CANCER_MUTATIONS for the mutations record

# Write data back to RESULTS_DISTANCES_P new columns


import os
import pickle
import pandas as pd
from tqdm import tqdm
from definitions import *



def get_gene_info_cache():
    """
    Scans the KEGG_GENES_P directory once to cache AA lengths.
    This avoids re-loading thousands of pickles for every cancer cohort.
    """
    cache = {}
    pickle_files = [f for f in os.listdir(KEGG_GENES_P) if f.endswith('.pickle')]

    for f in tqdm(pickle_files, desc="Caching gene lengths"):
        path = os.path.join(KEGG_GENES_P, f)
        try:
            with open(path, 'rb') as pf:
                gene_obj = pickle.load(pf)
                # Get AA length (0 if not CDS/missing)
                aa_len = len(gene_obj.aa_seq) if (hasattr(gene_obj, 'aa_seq') and gene_obj.aa_seq) else 0
                cache[gene_obj.kegg_id] = aa_len
        except Exception:
            continue
    return cache


def run_coverage_calculation():
    # 1. Load Pathway Metadata (Names and Gene Lists)
    if not os.path.exists(KEGG_PATHWAY_METADATA_P):
        print(f"Error: Metadata not found at {KEGG_PATHWAY_METADATA_P}")
        return

    with open(KEGG_PATHWAY_METADATA_P, 'rb') as f:
        pathway_metadata = pickle.load(f)

    # 2. Cache Gene Lengths
    gene_lengths = get_gene_info_cache()

    # 3. Iterate through result files
    result_files = [f for f in os.listdir(RESULTS_DISTANCES_P) if f.endswith('.csv')]

    if not result_files:
        print(f"No result files found in {RESULTS_DISTANCES_P}")
        return

    for res_file in result_files:
        print(f"            Updating {res_file}...")
        res_path = os.path.join(RESULTS_DISTANCES_P, res_file)

        # Load the results dataframe
        results_df = pd.read_csv(res_path)
        if results_df.empty:
            continue

        # Attempt to find the pathway ID column
        # Results usually have a column named 'pathway_id' or the first column is the ID
        pw_col = 'pathway_id' if 'pathway_id' in results_df.columns else results_df.columns[0]

        # 4. Load corresponding cancer mutations
        cancer_mut_path = os.path.join(CBIO_CANCER_MUTATIONS_P, res_file)
        if not os.path.exists(cancer_mut_path):
            print(f"Warning: Mutation file not found for {res_file}, skipping.")
            continue

        # Optimized Loading: Single KeggId per row
        mut_df = pd.read_csv(cancer_mut_path, usecols=['KeggId'])
        gene_mutation_counts = mut_df['KeggId'].value_counts().to_dict()

        # 5. Calculate Metrics
        coverage_rows = []

        for _, row in tqdm(results_df.iterrows(), total=len(results_df), desc=f"Calculating coverage for {res_file}"):
            pw_id = row[pw_col]

            if pw_id not in pathway_metadata:
                coverage_rows.append([0, 0, 0, 0, "Unknown"])
                continue

            pw_info = pathway_metadata[pw_id]
            pw_genes = pw_info['genes_ids']
            pw_name = pw_info.get('name', '')

            total_mutations = 0
            total_aa_length = 0
            num_genes = len(pw_genes)
            num_covered_genes = 0

            for g_id in pw_genes:
                # Mutation count for this gene in this cancer
                m_count = gene_mutation_counts.get(g_id, 0)
                total_mutations += m_count

                # AA length from cache
                g_len = gene_lengths.get(g_id, 0)
                total_aa_length += g_len

                # TODO: Lotem, read this and check the rule is correct
                #  Also, lets think this rule through or very long/short sequences
                
                # Coverage Rule: 10+ mutations OR density >= 1% (COVERAGE_PERCENTAGE_THRESHOLD)
                covered = False
                if m_count >= ABSOLUTE_COUNT_THRESHOLD:
                    covered = True
                elif g_len > 0:
                    # (mutations / length) * 100
                    if (m_count / g_len) * 100 >= COVERAGE_PERCENTAGE_THRESHOLD:
                        covered = True

                if covered:
                    num_covered_genes += 1

            coverage_rows.append([
                total_mutations,
                total_aa_length,
                num_genes,
                num_covered_genes,
                pw_name
            ])

        # 6. Merge and Update
        cols = ['num_mutations', 'total_aa_length', 'num_genes', 'num_covered_genes', 'pathway_name']
        stats_df = pd.DataFrame(coverage_rows, columns=cols)

        # Remove existing versions of these columns so they are "fresh"
        results_df = results_df.drop(columns=cols, errors='ignore')

        for col in cols:
            results_df[col] = stats_df[col]

        results_df.to_csv(res_path, index=False)
        print(f"            Successfully updated {res_file}")


if __name__ == "__main__":
    run_coverage_calculation()
    print("\n\nFinished.")