import os
import torch
from tqdm import tqdm

from definitions import KEGG_HOMO_SAPIENS, ESM_EMBEDDINGS_P
from kegg_api import KeggApi, KeggGene
from p4_esm_emb_and_log_probs import ScoringCalculator


def precompute_all_logits():
    # 1. Setup
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"Starting precomputation on {device}")

    # Initialize your ESM model wrapper
    esm_tool = ScoringCalculator(device=device)
    api = KeggApi()

    # 2. Get all genes
    print("Fetching master gene list from KEGG...")
    all_genes_dict = api.get_all_genes(species=KEGG_HOMO_SAPIENS)
    gene_ids = list(all_genes_dict.keys())

    os.makedirs(ESM_EMBEDDINGS_P, exist_ok=True)

    # 3. Iterate and Compute
    # We use a standard loop because GPU inference doesn't benefit from Python multithreading
    # and we want to avoid OOM errors.

    cds_count = 0
    skipped_count = 0

    for kegg_id in tqdm(gene_ids, desc="Precomputing ESM Logits"):
        try:
            # check if exists first to allow resuming an interrupted run
            file_path = os.path.join(ESM_EMBEDDINGS_P, f"{kegg_id}.pt")

            # 1. Load the gene object (cached locally)
            gene = KeggGene(kegg_id, kegg_api=api)

            # 2. Filter: Only CDS genes have AA sequences for ESM
            if gene.coding_type != "CDS" or not gene.aa_seq:
                skipped_count += 1
                continue

            # 3. Check if already computed
            if os.path.exists(file_path):
                # Optional: Add logic here to verify file integrity if needed
                continue

            # 4. Compute and Save
            # This uses your existing get_or_compute_logits logic
            esm_tool.get_or_compute_logits(kegg_id)
            cds_count += 1

        except Exception as e:
            print(f"\n[Error] Failed processing {kegg_id}: {e}")
            continue

    print("\nPrecomputation Complete!")
    print(f"Logits saved to: {ESM_EMBEDDINGS_P}")
    print(f"CDS processed: {cds_count}")
    print(f"Non-CDS skipped: {skipped_count}")


if __name__ == "__main__":
    # Ensure no other heavy GPU processes are running
    with torch.no_grad():
        precompute_all_logits()