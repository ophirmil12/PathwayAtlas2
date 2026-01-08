from definitions import *
from kegg_api import *
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed

# Recommendations:
# KEGG is generally okay with ~10-15 concurrent requests.
# Don't set this too high (e.g., 100) or you might get IP blocked.
MAX_WORKERS = 10


def download_single_gene(kegg_id, api):
    """Worker function for a single thread."""
    try:
        # redownload=True forces KeggGene to call the API and overwrite the local pickle
        gene_obj = KeggGene(kegg_id, redownload=False, kegg_api=api)

        if not gene_obj.aa_seq:
            return f"Warning: {kegg_id} has no AA sequence (Type: {gene_obj.coding_type})"
        return None  # Success
    except Exception as e:
        return f"Error processing {kegg_id}: {str(e)}"


def recreate_all_genes():
    api = KeggApi()

    print("Fetching list of all human genes from KEGG...")
    all_genes_dict = api.get_all_genes(species=KEGG_HOMO_SAPIENS)
    gene_ids = list(all_genes_dict.keys())

    print(f"Found {len(gene_ids)} genes. Starting multithreaded recreation (Workers: {MAX_WORKERS})...")

    # Use ThreadPoolExecutor for I/O operations
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        # Map the IDs to the worker function
        future_to_gene = {executor.submit(download_single_gene, gid, api): gid for gid in gene_ids}

        # Process as they complete for the progress bar
        for future in tqdm(as_completed(future_to_gene), total=len(gene_ids), desc="Downloading Genes"):
            result = future.result()
            if result:  # If it returned a warning/error string
                print(result)


if __name__ == "__main__":
    recreate_all_genes()