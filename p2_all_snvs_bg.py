import os
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed

from definitions import (
    KEGG_HOMO_SAPIENS,
    KEGG_GENE_SCORES_P,
    KEGG_API_RECOMMENDED_WORKERS
)
from kegg_api import KeggApi, KeggGene


def process_gene_by_id(kegg_id, api):
    """
    Worker function:
    1. Loads/Downloads KeggGene
    2. Checks if CDS
    3. Generates SNV table
    """
    try:
        # This will load from pickle if it exists, otherwise it fetches from API
        # We pass the shared 'api' session to avoid creating new ones in every thread
        gene = KeggGene(kegg_id, redownload=False, kegg_api=api)

        # Only process Protein Coding Sequences
        if gene.coding_type != "CDS":
            return None  # Skip non-coding

        # Format filename: 'hsa:249' -> 'hsa_249.csv'
        file_name = kegg_id.replace(":", "_") + ".csv"
        out_path = os.path.join(KEGG_GENE_SCORES_P, file_name)

        # Generate and save the SNVs
        gene.all_snvs(outpath=out_path, index=True)

        return kegg_id

    except Exception as e:
        return f"Error with {kegg_id}: {str(e)}"


def run_snv_generation_pipeline():
    # 1. Initialize the API
    api = KeggApi()

    # 2. Ensure output directory exists
    os.makedirs(KEGG_GENE_SCORES_P, exist_ok=True)

    # 3. Get all human gene IDs from the API
    print(f"Fetching full gene list for {KEGG_HOMO_SAPIENS}...")
    all_genes_dict = api.get_all_genes(species=KEGG_HOMO_SAPIENS)
    gene_ids = list(all_genes_dict.keys())

    print(f"Starting SNV generation for {len(gene_ids)} genes...")

    # 4. Process using ThreadPoolExecutor
    # We use the recommended worker count to stay within KEGG's politeness limits
    with ThreadPoolExecutor(max_workers=KEGG_API_RECOMMENDED_WORKERS) as executor:
        # Submit all tasks
        futures = {executor.submit(process_gene_by_id, gid, api): gid for gid in gene_ids}

        cds_count = 0
        for future in tqdm(as_completed(futures), total=len(gene_ids), desc="Processing Genes"):
            res = future.result()

            if res:
                if res.startswith("Error"):
                    print(f"\n[!] {res}")
                else:
                    cds_count += 1

    print(f"\nPipeline Complete!")
    print(f"Total Genes Checked: {len(gene_ids)}")
    print(f"CDS SNV Tables Generated: {cds_count}")
    print(f"Output Directory: {KEGG_GENE_SCORES_P}")


if __name__ == "__main__":
    run_snv_generation_pipeline()