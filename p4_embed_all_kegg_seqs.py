# Embedding all Kegg bg sequences, and only saving the raw logits


import os
import torch
import esm  # pip install fair-esm
from tqdm import tqdm

from definitions import *
from kegg_api import KeggApi, KeggGene
from p4_esm_emb_and_log_probs import ScoringCalculator


def precompute_all_logits():
    # 1. Verify Environment for Slurm
    torch_home = os.environ.get('TORCH_HOME', 'Not Set')
    print(f"CUDA Available: {torch.cuda.is_available()}")
    print(f"TORCH_HOME is: {torch_home}")

    # Ensure the cache directory exists (relative to where the slurm script runs)
    if torch_home != 'Not Set':
        os.makedirs(torch_home, exist_ok=True)

    # 2. Initialize API
    # Create the session once to be shared across all KeggGene objects
    api = KeggApi()

    # 3. ESM Model Loading
    print(f"Loading ESM model: {ESM1B_MODEL}...")
    # This will respect the TORCH_HOME environment variable automatically
    model, alphabet = esm.pretrained.load_model_and_alphabet(ESM1B_MODEL)

    # Initialize your wrapper
    esm_tool = ScoringCalculator(model, alphabet)

    # 4. Get all human genes
    print("Fetching master gene list from KEGG...")
    all_genes_dict = api.get_all_genes(species=KEGG_HOMO_SAPIENS)
    gene_ids = list(all_genes_dict.keys())

    # Ensure output directory exists
    os.makedirs(ESM_EMBEDDINGS_P, exist_ok=True)

    cds_count = 0
    skipped_count = 0
    already_exists_count = 0

    # 5. Iterate and Compute
    # with torch.no_grad() is critical to avoid VRAM accumulation
    with torch.no_grad():
        for kegg_id in tqdm(gene_ids, desc="Precomputing ESM Logits"):
            try:
                # Check if file exists first to allow resuming if Slurm job times out
                file_path = os.path.join(ESM_EMBEDDINGS_P, f"{kegg_id}.pt")
                if os.path.exists(file_path):
                    already_exists_count += 1
                    continue

                # Load gene object (using shared api session for speed/safety)
                gene = KeggGene(kegg_id, kegg_api=api)

                # Skip non-coding genes
                if gene.coding_type != "CDS" or not gene.aa_seq:
                    skipped_count += 1
                    continue

                # Compute and Save (handles long proteins internally)
                esm_tool.get_or_compute_logits(kegg_id)
                cds_count += 1

            except Exception as e:
                print(f"\n[Error] Failed processing {kegg_id}: {e}")
                continue

    print("\n" + "=" * 40)
    print("PRECOMPUTATION SUMMARY")
    print(f"Newly Computed:     {cds_count}")
    print(f"Already Existed:    {already_exists_count}")
    print(f"Non-CDS (Skipped):  {skipped_count}")
    print(f"Logits Directory:   {ESM_EMBEDDINGS_P}")
    print("=" * 40)


if __name__ == "__main__":
    precompute_all_logits()