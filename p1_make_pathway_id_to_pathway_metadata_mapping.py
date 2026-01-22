# Creating a dictionary and saving as KEGG_PATHWAY_METADATA_FILE:
# {"pathway id":
#   {"genes_ids": ["id 1", "id 2"...],
#    "name": "pathway description:,
#    "type": "module"/"pathway"
# ...
# }


import os
import pickle
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm
from kegg_api import KeggApi
from definitions import *


def fetch_network_metadata(kegg_id, name, network_type, api):
    """Worker function to fetch gene list and format metadata for one network."""
    try:
        # get_gene_list handles both pathways (hsa) and modules (M)
        genes = api.get_gene_list(kegg_id)

        return kegg_id, {
            "genes_ids": list(genes),
            "name": name,
            "type": network_type,
            # You can add more fields here in the future
        }
    except Exception as e:
        print(f"Error fetching metadata for {kegg_id}: {e}")
        return None


def create_pathway_metadata_file():
    api = KeggApi()
    metadata_dict = {}

    print("Fetching lists of all pathways and modules...")
    pathways = api.get_all_pathways(species=KEGG_HOMO_SAPIENS)  # {id: desc}
    modules = api.get_all_modules()  # {id: desc}

    # Prepare task list: (id, name, type)
    tasks = []
    for pid, name in pathways.items():
        tasks.append((pid, name, "pathway"))
    for mid, name in modules.items():
        tasks.append((mid, name, "module"))

    print(f"Processing {len(tasks)} networks...")

    # Use ThreadPoolExecutor for faster API fetching
    with ThreadPoolExecutor(max_workers=10) as executor:
        future_to_id = {
            executor.submit(fetch_network_metadata, t[0], t[1], t[2], api): t[0]
            for t in tasks
        }

        for future in tqdm(as_completed(future_to_id), total=len(tasks), desc="Building Metadata"):
            result = future.result()
            if result:
                kegg_id, info = result
                metadata_dict[kegg_id] = info

    # Ensure the directory exists
    os.makedirs(os.path.dirname(KEGG_PATHWAY_METADATA_FILE), exist_ok=True)

    # Save as pickle (consistent with your other data)
    with open(KEGG_PATHWAY_METADATA_FILE, 'wb') as f:
        pickle.dump(metadata_dict, f)

    print(f"Successfully saved metadata for {len(metadata_dict)} entries to {KEGG_PATHWAY_METADATA_FILE}.")


if __name__ == "__main__":
    create_pathway_metadata_file()