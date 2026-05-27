from matplotlib import lines
import requests
import pandas as pd
import numpy as np
from definitions import *
import os
from uniprot_api import UniprotApi
from Bio import PDB
import io
import networkx as nx
import glob

BRIDGES_DIR = os.path.join(BASE_P, 'results', 'p16_bottleneck_analysis')

def pathway_interaction_matrix(pathway_id: str, score_threshold: int = 700):
    """
    Two API calls: KEGG gene list → STRING physical interaction matrix.
    Example: pathway_interaction_matrix("hsa04110")
    """

    # ── Call 1: KEGG link API → gene list ────────────────────────────────────
    # This endpoint directly returns all genes in a pathway, no KGML parsing needed
    kegg_resp = requests.get(
        f"https://rest.kegg.jp/link/hsa/{pathway_id}"
    )
    kegg_resp.raise_for_status()

    # Response is TSV: pathway_id \t hsa:ENTREZ_ID
    entrez_ids = [
        line.split("\t")[1].strip()          # e.g. "hsa:1017"
        for line in kegg_resp.text.strip().splitlines()
        if line
    ]

    # ── Call 2: STRING network API → interactions ─────────────────────────────
    string_resp = requests.post(
        "https://string-db.org/api/tsv/network",
        data={
            "identifiers":    "\r".join(entrez_ids),
            "species":        9606,
            "network_type":   "physical",
            "required_score": score_threshold,
            "caller_identity": "pathway_atlas"
        }
    )
    string_resp.raise_for_status()

    lines = string_resp.text.strip().splitlines()
    if len(lines) <= 1:
        print("No interactions found — try lowering score_threshold to 400.")
        return pd.DataFrame()

    # Parse TSV → keep only gene symbols and score
    header = lines[0].split("\t")
    df = pd.DataFrame(
        [l.split("\t") for l in lines[1:]],
        columns=header
    )[["preferredName_A", "preferredName_B", "score"]]
    df["score"] = df["score"].astype(float)

    # ── Build binary adjacency matrix ─────────────────────────────────────────
    genes = sorted(set(df["preferredName_A"]) | set(df["preferredName_B"]))
    mat = pd.DataFrame(0, index=genes, columns=genes)

    for _, row in df.iterrows():
        mat.loc[row["preferredName_A"], row["preferredName_B"]] = 1
        mat.loc[row["preferredName_B"], row["preferredName_A"]] = 1  # symmetric

    return mat


def find_bottlenecks(mat: pd.DataFrame) -> pd.DataFrame:

    G = nx.from_pandas_adjacency(mat)

    bc  = nx.betweenness_centrality(G, normalized=True)
    ap  = set(nx.articulation_points(G))
    deg = dict(G.degree())

    results = pd.DataFrame({
        "gene":                    list(bc.keys()),
        "betweenness_centrality":  list(bc.values()),
        "degree":                  [deg[g] for g in bc.keys()],
        "is_articulation_point":   [g in ap for g in bc.keys()],
    }).sort_values("betweenness_centrality", ascending=False).reset_index(drop=True)

    return results

def find_bridges(adj: pd.DataFrame, bottlenecks: pd.DataFrame) -> pd.DataFrame:
    """
    Ranks edges by structural criticality:
    1. Bridges first — removing this single edge disconnects the graph
    2. Among non-bridges, rank by geometric mean BC of endpoints
       (the heuristic proxy, only as a fallback for non-bridge edges)
    """
    G   = nx.from_pandas_adjacency(adj)
    bc  = bottlenecks.set_index("gene")["betweenness_centrality"].to_dict()
    deg = bottlenecks.set_index("gene")["degree"].to_dict()
    
    bridge_edges = set(nx.bridges(G))  # set of (u, v) tuples

    pairs = []
    for u, v in G.edges():
        is_bridge = (u, v) in bridge_edges or (v, u) in bridge_edges
        bc_geo    = np.sqrt(bc.get(u, 0) * bc.get(v, 0))

        pairs.append({
            "gene_a":    u,
            "gene_b":    v,
            "bc_a":      bc.get(u, 0),
            "bc_b":      bc.get(v, 0),
            "degree_a":  deg.get(u, 0),
            "degree_b":  deg.get(v, 0),
            "is_bridge": is_bridge,
            "bc_geo":    bc_geo,
        })

    df = pd.DataFrame(pairs)

    # Sort: bridges first, then by BC geometric mean within each group
    df = (df.sort_values(["is_bridge", "bc_geo"], ascending=[False, False])
            .reset_index(drop=True))

    return df



def find_pdb_complex(uniprot_a: str, uniprot_b: str) -> list[str]:
    """
    Search RCSB for structures containing BOTH UniProt accessions
    (i.e. a co-crystal of the two proteins).
    Returns list of PDB IDs, best resolution first.
    """
    query = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_polymer_entity_container_identifiers"
                                     ".reference_sequence_identifiers.database_accession",
                        "operator": "exact_match",
                        "value": uniprot_a
                    }
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_polymer_entity_container_identifiers"
                                     ".reference_sequence_identifiers.database_accession",
                        "operator": "exact_match",
                        "value": uniprot_b
                    }
                }
            ]
        },
        "return_type": "entry",
        "request_options": {
            "sort": [{"sort_by": "rcsb_entry_info.resolution_combined", "direction": "asc"}],
            "results_verbosity": "compact"
        }
    }

    resp = requests.post("https://search.rcsb.org/rcsbsearch/v2/query", json=query)
    resp.raise_for_status()
    if resp.status_code == 204 or not resp.content:
        return []
    results = resp.json().get("result_set", [])
    return results


def assign_structure_strategy(pairs: pd.DataFrame,
                              uniprot_map: dict) -> pd.DataFrame:
    """
    For each pair, checks PDB for an existing co-crystal.
    Assigns 'pdb' or 'af3' strategy accordingly.
    uniprot_map: {gene_symbol: uniprot_accession}
    """
    strategies = []
    for _, row in pairs.iterrows():
        up_a = uniprot_map.get(row["gene_a"])
        up_b = uniprot_map.get(row["gene_b"])

        if not up_a or not up_b:
            strategies.append({"strategy": "no_uniprot", "pdb_id": None})
            continue

        pdb_hits = find_pdb_complex(up_a, up_b)  # from previous step

        if pdb_hits:
            strategies.append({"strategy": "pdb", "pdb_id": pdb_hits[0]})
        else:
            strategies.append({"strategy": "af3", "pdb_id": None})

    return pairs.join(pd.DataFrame(strategies))


# Get UniProt IDs for all genes via UniProt API
def get_uniprot_ids(genes: list[str]) -> dict:
    """
    Looks up UniProt accessions for human gene symbols.
    Uses the UniProt REST API with gene_exact field.
    Falls back to a gene name search if exact match fails.
    """
    uniprot_map = {}
    uniprot_api = UniprotApi()

    for gene in genes:
 
        uniprot_map[gene] = uniprot_api.uid_from_name(gene)['reviewed'][0]  # take first reviewed entry
        print(f"Mapped {gene} → {uniprot_map[gene]}")
    return uniprot_map



# Genes absent from all cancer mutation files, identified by scanning pathway SNV files.
# KEGG IDs verified against KEGG gene pickle ref_names.
MISSING_GENE_KEGG_IDS = {
    'APELA':  'hsa:100506013',
    'BECN2':  'hsa:441925',
    'NFKBIE': 'hsa:4794',
    'RAB7B':  'hsa:338382',
}


def build_sequence_map():
    id_map = {
        'APELA':  'hsa:100506013',
        'BECN2':  'hsa:441925',
        'NFKBIE': 'hsa:4794',
        'RAB7B':  'hsa:338382',
    }
    for f in glob.glob(os.path.join(CBIO_CANCER_MUTATIONS_P, '*.csv')):
        df = pd.read_csv(f, usecols=['Protein', 'KeggId']).dropna(subset=['KeggId'])
        for protein, id in df.drop_duplicates('Protein')[['Protein', 'KeggId']].values:
            if protein not in id_map:
                id_map[protein] = id


    print(f"KeggId map: {len(id_map)} genes")
    return id_map


def add_ids_to_bridges(id_map):
    for f in sorted(glob.glob(os.path.join(BRIDGES_DIR, '*_bridges.csv'))):
        df = pd.read_csv(f)
        df['keggid_a'] = df['gene_a'].map(id_map)
        df['keggid_b'] = df['gene_b'].map(id_map)
        missing = (set(df.loc[df['seq_a'].isna(), 'gene_a']) |
                   set(df.loc[df['seq_b'].isna(), 'gene_b']))
        if missing:
            print(f"WARNING {os.path.basename(f)}: no id for {missing}")
        print("Saving df...")
        df.to_csv(f, index=False)
        print(f"Updated {os.path.basename(f)}  ({len(df)} bridges, "
              f"{df['keggid_a'].notna().sum()} with keggid_a, {df['keggid_b'].notna().sum()} with keggid_b)")
        return df


if __name__ == '__main__':
    
    pathway_ids = ['hsa04080', 'hsa04722', 'hsa04216', 'hsa04137']

    for pathway_id in pathway_ids:
        print(f"\nProcessing pathway {pathway_id}...")

        print("Generating interaction matrix...")
        matrix_output_path = pjoin(STRING_P, f"{pathway_id}.csv")
        if os.path.exists(matrix_output_path):
            print(f"Adjacency matrix for pathway {pathway_id} already saved.")
            mat = pd.read_csv(matrix_output_path, index_col=0)
        else:
            mat = pathway_interaction_matrix(pathway_id)
            if not mat.empty:
                mat.to_csv(matrix_output_path)
        
        bottleneck_output_path = pjoin(BOTTLENECKS_P, f"{pathway_id}_bottlenecks.csv")
        if os.path.exists(bottleneck_output_path):
            print(f"Bottleneck analysis for pathway {pathway_id} already saved.")
            bottlenecks = pd.read_csv(bottleneck_output_path)
        else:
            bottlenecks = find_bottlenecks(mat)
            if not bottlenecks.empty:
                bottlenecks.to_csv(bottleneck_output_path)

        print("Finding bridges and assigning structurue strategy...")

        bridges_output_path = pjoin(BOTTLENECKS_P, f"{pathway_id}_bridges.csv")
        if os.path.exists(bridges_output_path):
            print(f"Bridge analysis for pathway {pathway_id} already saved.")
            bridges = pd.read_csv(bridges_output_path)
        else:
            bridges = find_bridges(mat, bottlenecks)
            uniprot_map = get_uniprot_ids(bottlenecks["gene"].tolist())
            ranked_pairs = assign_structure_strategy(bridges, uniprot_map)

            if not ranked_pairs.empty:
                ranked_pairs.to_csv(bridges_output_path)

        

        

        

        

