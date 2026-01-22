# Calculate the log-probs (WT-marginals) and put into the CSVs


import os
import torch
import esm
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed

from definitions import *
from kegg_api import KeggApi
from p4_esm_emb_and_log_probs import ScoringCalculator


def process_gene_scores(kegg_id, calculator):
    """
    Worker function to update a single gene's SNV CSV with ESM scores.
    """
    try:
        # 1. Construct path to the existing SNV CSV
        # Format: hsa:249 -> hsa_249.csv
        file_name = kegg_id.replace(":", "_") + ".csv"
        csv_path = os.path.join(KEGG_GENE_SCORES_P, file_name)

        # 2. Skip if the SNV file doesn't exist
        # (This happens for non-CDS genes or genes that failed earlier steps)
        if not os.path.exists(csv_path):
            return None

        # 3. Use the ScoringCalculator to compute scores and update the CSV
        # This method:
        #  - Loads logits from ESM_EMBEDDINGS_P/{kegg_id}.pt
        #  - Calculates log P(mut) - log P(wt)
        #  - Updates the CSV at csv_path
        calculator.save_mutation_scores_to_csv(kegg_id, csv_path)

        return kegg_id

    except Exception as e:
        return f"Error with {kegg_id}: {str(e)}"


def run_log_prob_calculation():
    # 1. Initialize API and Environment
    api = KeggApi()
    device = 'cuda' if torch.cuda.is_available() else 'cpu'

    # 2. Load ESM Model and Alphabet (Needed for ScoringCalculator)
    print(f"Loading ESM model ({ESM1B_MODEL}) to initialize ScoringCalculator...")
    model, alphabet = esm.pretrained.load_model_and_alphabet(ESM1B_MODEL)
    calculator = ScoringCalculator(model, alphabet)

    # 3. Get all human genes
    print("(p4B) Fetching master gene list from KEGG...")
    all_genes_dict = api.get_all_genes(species=KEGG_HOMO_SAPIENS)
    gene_ids = list(all_genes_dict.keys())

    print(f"Starting log-prob calculations for {len(gene_ids)} genes...")

    # 4. Use ThreadPoolExecutor
    # Calculation is done on the GPU/CPU using Tensors,
    # but the bottleneck is often reading/writing CSVs.
    # We use a moderate worker count.
    with ThreadPoolExecutor(max_workers=4) as executor:
        futures = {executor.submit(process_gene_scores, gene_id, calculator): gene_id for gene_id in gene_ids}

        success_count = 0
        for future in tqdm(as_completed(futures), total=len(gene_ids), desc="Updating SNV CSVs"):
            res = future.result()
            if res:
                if res.startswith("Error"):
                    print(f"\n[!] {res}")
                else:
                    success_count += 1

    print(f"\nCalculation Complete!")
    print(f"SNV CSVs updated: {success_count} / 20532 CDS")


if __name__ == "__main__":
    # no_grad is vital to keep memory usage low
    with torch.no_grad():
        run_log_prob_calculation()
