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


def pathway_interaction_matrix(pathway_id: str, score_threshold: int = 700):
    kegg_resp = requests.get(f"https://rest.kegg.jp/link/hsa/{pathway_id}")
    kegg_resp.raise_for_status()

    entrez_ids = [
        line.split("\t")[1].strip()
        for line in kegg_resp.text.strip().splitlines()
        if line
    ]

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

    header = lines[0].split("\t")
    df = pd.DataFrame(
        [l.split("\t") for l in lines[1:]],
        columns=header
    )[["preferredName_A", "preferredName_B", "score"]]
    df["score"] = df["score"].astype(float)

    genes = sorted(set(df["preferredName_A"]) | set(df["preferredName_B"]))
    mat = pd.DataFrame(0, index=genes, columns=genes)

    for _, row in df.iterrows():
        mat.loc[row["preferredName_A"], row["preferredName_B"]] = 1
        mat.loc[row["preferredName_B"], row["preferredName_A"]] = 1

    return mat


def find_bottlenecks(mat: pd.DataFrame) -> pd.DataFrame:
    G   = nx.from_pandas_adjacency(mat)
    bc  = nx.betweenness_centrality(G, normalized=True)
    ap  = set(nx.articulation_points(G))
    deg = dict(G.degree())

    return (pd.DataFrame({
        "gene":                  list(bc.keys()),
        "betweenness_centrality": list(bc.values()),
        "degree":                [deg[g] for g in bc.keys()],
        "is_articulation_point": [g in ap for g in bc.keys()],
    })
    .sort_values("betweenness_centrality", ascending=False)
    .reset_index(drop=True))


def find_critical_edges(adj: pd.DataFrame) -> pd.DataFrame:
    """
    Ranks all edges by edge betweenness centrality (EBC) —
    the fraction of all shortest paths between every pair of nodes
    that pass through each edge.

    High EBC = many pairwise communications in the pathway route
    through this specific interaction. Breaking it disrupts the most paths.
    """
    G   = nx.from_pandas_adjacency(adj)
    ebc = nx.edge_betweenness_centrality(G, normalized=True)
    deg = dict(G.degree())

    return (pd.DataFrame([
                {
                    "gene_a":   u,
                    "gene_b":   v,
                    "edge_bc":  score,
                    "degree_a": deg[u],
                    "degree_b": deg[v],
                }
                for (u, v), score in ebc.items()
            ])
            .sort_values("edge_bc", ascending=False)
            .reset_index(drop=True))


def find_pdb_complex(uniprot_a: str, uniprot_b: str) -> list[str]:
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
    return resp.json().get("result_set", [])


def assign_structure_strategy(pairs: pd.DataFrame, uniprot_map: dict) -> pd.DataFrame:
    strategies = []
    for _, row in pairs.iterrows():
        up_a = uniprot_map.get(row["gene_a"])
        up_b = uniprot_map.get(row["gene_b"])

        if not up_a or not up_b:
            strategies.append({"strategy": "no_uniprot", "pdb_id": None,
                                "uniprot_a": up_a, "uniprot_b": up_b})
            continue

        pdb_hits = find_pdb_complex(up_a, up_b)

        if pdb_hits:
            strategies.append({"strategy": "pdb", "pdb_id": pdb_hits[0],
                                "uniprot_a": up_a, "uniprot_b": up_b})
        else:
            strategies.append({"strategy": "af3", "pdb_id": None,
                                "uniprot_a": up_a, "uniprot_b": up_b})

    return pairs.join(pd.DataFrame(strategies))


def get_uniprot_ids(genes: list[str]) -> dict:
    import pickle

    if os.path.exists(UNIPROT_CACHE_P):
        with open(UNIPROT_CACHE_P, 'rb') as f:
            uniprot_map = pickle.load(f)
        print(f"Loaded UniProt cache ({len(uniprot_map)} entries) from {UNIPROT_CACHE_P}")
    else:
        uniprot_map = {}

    missing = [g for g in genes if g not in uniprot_map]
    if missing:
        uniprot_api = UniprotApi()
        for gene in missing:
            uniprot_map[gene] = uniprot_api.uid_from_name(gene)['reviewed'][0]
            print(f"Mapped {gene} → {uniprot_map[gene]}")
        with open(UNIPROT_CACHE_P, 'wb') as f:
            pickle.dump(uniprot_map, f)
        print(f"Saved UniProt cache ({len(uniprot_map)} entries) → {UNIPROT_CACHE_P}")
    else:
        print(f"All {len(genes)} genes already in UniProt cache — skipping API calls")

    return uniprot_map


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


def add_ids_to_critical_edges(id_map):
    for f in sorted(glob.glob(os.path.join(BRIDGES_DIR, '*_critical_edges.csv'))):
        df = pd.read_csv(f)
        df['keggid_a'] = df['gene_a'].map(id_map)
        df['keggid_b'] = df['gene_b'].map(id_map)
        missing = (set(df.loc[df['keggid_a'].isna(), 'gene_a']) |
                   set(df.loc[df['keggid_b'].isna(), 'gene_b']))
        if missing:
            print(f"WARNING {os.path.basename(f)}: no id for {missing}")
        df.to_csv(f, index=False)
        print(f"Updated {os.path.basename(f)}  ({len(df)} critical edges, "
              f"{df['keggid_a'].notna().sum()} with keggid_a, "
              f"{df['keggid_b'].notna().sum()} with keggid_b)")
    return df

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns

def plot_edge_bc_distributions(pathway_ids: list[str], save_path: str = None):
    """
    For each pathway, loads its critical_edges CSV and plots the distribution
    of edge betweenness centrality values, so you can visually identify
    which values are genuinely extreme vs background noise.

    Shows:
      - Per-pathway KDE + rug plot (one panel per pathway)
      - Shared reference lines at 90th and 95th percentile across all pathways
      - Top 5 edges labeled by name on each panel
    """
    # ── Load all data ─────────────────────────────────────────────────────────
    pathway_data = {}
    all_scores   = []

    for pathway_id in pathway_ids:
        path = pjoin(CRITICAL_EDGES_P, f"{pathway_id}_critical_edges.csv")
        if not os.path.exists(path):
            print(f"WARNING: no critical edges file for {pathway_id}, skipping.")
            continue
        df = pd.read_csv(path)
        if "edge_bc" not in df.columns or df.empty:
            continue
        pathway_data[pathway_id] = df
        all_scores.extend(df["edge_bc"].tolist())

    if not pathway_data:
        print("No data found — run the main pipeline first.")
        return

    all_scores   = np.array(all_scores)
    p90_global   = np.percentile(all_scores, 90)
    p95_global   = np.percentile(all_scores, 95)
    p99_global   = np.percentile(all_scores, 99)

    # ── Layout ────────────────────────────────────────────────────────────────
    n            = len(pathway_data)
    ncols        = 2
    nrows        = (n + 1) // ncols + 1   # +1 for the global panel at the top
    fig          = plt.figure(figsize=(14, 4 * nrows))
    gs           = gridspec.GridSpec(nrows, ncols, figure=fig, hspace=0.55, wspace=0.35)

    palette = sns.color_palette("tab10", n)

    # ── Top panel: all pathways overlaid ──────────────────────────────────────
    ax_global = fig.add_subplot(gs[0, :])

    for i, (pathway_id, df) in enumerate(pathway_data.items()):
        scores = df["edge_bc"]
        sns.kdeplot(scores, ax=ax_global, label=pathway_id,
                    color=palette[i], linewidth=1.8, fill=False)

    ax_global.axvline(p90_global, color="orange", linestyle="--", linewidth=1.2,
                      label=f"90th pct (global) = {p90_global:.4f}")
    ax_global.axvline(p95_global, color="tomato", linestyle="--", linewidth=1.2,
                      label=f"95th pct (global) = {p95_global:.4f}")
    ax_global.axvline(p99_global, color="darkred", linestyle="--", linewidth=1.2,
                      label=f"99th pct (global) = {p99_global:.4f}")

    ax_global.set_title("Edge Betweenness Centrality — all pathways overlaid",
                         fontsize=13, fontweight="bold")
    ax_global.set_xlabel("Edge Betweenness Centrality")
    ax_global.set_ylabel("Density")
    ax_global.legend(fontsize=8, ncol=2)
    sns.despine(ax=ax_global)

    # ── Per-pathway panels ────────────────────────────────────────────────────
    for i, (pathway_id, df) in enumerate(pathway_data.items()):
        row  = (i // ncols) + 1   # offset by 1 for global panel
        col  = i  % ncols
        ax   = fig.add_subplot(gs[row, col])
        scores = df["edge_bc"]

        # KDE
        sns.kdeplot(scores, ax=ax, color=palette[i], linewidth=1.8, fill=True, alpha=0.15)

        # Rug plot — each tick is one edge
        sns.rugplot(scores, ax=ax, color=palette[i], alpha=0.4, height=0.06)

        # Global percentile reference lines
        ax.axvline(p90_global, color="orange", linestyle="--", linewidth=1.0, alpha=0.8)
        ax.axvline(p95_global, color="tomato", linestyle="--", linewidth=1.0, alpha=0.8)
        ax.axvline(p99_global, color="darkred", linestyle="--", linewidth=1.0, alpha=0.8)

        # Local 95th percentile
        p95_local = scores.quantile(0.95)
        ax.axvline(p95_local, color=palette[i], linestyle=":", linewidth=1.2,
                   label=f"local 95th = {p95_local:.4f}")

        # Label top 5 edges by name on the x-axis as vertical ticks
        top5 = df.nlargest(5, "edge_bc")
        for _, edge in top5.iterrows():
            label = f"{edge['gene_a']}–{edge['gene_b']}"
            ax.annotate(label,
                        xy=(edge["edge_bc"], 0),
                        xytext=(edge["edge_bc"], ax.get_ylim()[1] * 0.5),
                        fontsize=5.5, rotation=90, ha="center", color=palette[i],
                        arrowprops=dict(arrowstyle="-", color=palette[i], lw=0.5))

        ax.set_title(pathway_id, fontsize=10, fontweight="bold")
        ax.set_xlabel("Edge BC", fontsize=8)
        ax.set_ylabel("Density", fontsize=8)
        ax.legend(fontsize=7)
        sns.despine(ax=ax)

    # ── Summary stats printed below ───────────────────────────────────────────
    print("\nGlobal EBC thresholds across all pathways:")
    print(f"  90th percentile: {p90_global:.4f}")
    print(f"  95th percentile: {p95_global:.4f}")
    print(f"  99th percentile: {p99_global:.4f}")
    print("\nPer-pathway top edge:")
    for pathway_id, df in pathway_data.items():
        top = df.iloc[0]
        print(f"  {pathway_id}: {top['gene_a']}–{top['gene_b']}  EBC={top['edge_bc']:.4f}")

    fig.suptitle("Edge Betweenness Centrality Distributions\n"
                 "Dashed lines: global 90th / 95th / 99th percentiles",
                 fontsize=14, fontweight="bold", y=1.01)

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"\nPlot saved to {save_path}")
    else:
        plt.show()

    return p90_global, p95_global, p99_global

if __name__ == '__main__':

    pathway_ids = ['hsa04080', 'hsa04722', 'hsa04216', 'hsa04137']

    for pathway_id in pathway_ids:
        print(f"\nProcessing pathway {pathway_id}...")

        matrix_output_path = pjoin(STRING_P, f"{pathway_id}.csv")
        if os.path.exists(matrix_output_path):
            print(f"Adjacency matrix for pathway {pathway_id} already saved.")
            mat = pd.read_csv(matrix_output_path, index_col=0)
        else:
            print("Generating interaction matrix...")
            mat = pathway_interaction_matrix(pathway_id)
            if not mat.empty:
                mat.to_csv(matrix_output_path)

        # bottleneck_output_path = pjoin(CRITICAL_EDGES_P, f"{pathway_id}_bottlenecks.csv")
        # if os.path.exists(bottleneck_output_path):
        #     print(f"Bottleneck analysis for pathway {pathway_id} already saved.")
        #     bottlenecks = pd.read_csv(bottleneck_output_path)
        # else:
        #     bottlenecks = find_bottlenecks(mat)
        #     if not bottlenecks.empty:
        #         bottlenecks.to_csv(bottleneck_output_path, index=False)

        critical_edges_output_path = pjoin(CRITICAL_EDGES_P, f"{pathway_id}_critical_edges.csv")
        if os.path.exists(critical_edges_output_path):
            print(f"Critical edge analysis for pathway {pathway_id} already saved.")
            critical_edges = pd.read_csv(critical_edges_output_path)
        else:
            print("Finding critical edges and assigning structure strategy...")
            critical_edges = find_critical_edges(mat)
            genes          = list(set(critical_edges["gene_a"]) | set(critical_edges["gene_b"]))
            uniprot_map    = get_uniprot_ids(genes)
            critical_edges = assign_structure_strategy(critical_edges, uniprot_map)
            if not critical_edges.empty:
                critical_edges.to_csv(critical_edges_output_path, index=False)

        # Plot after all pathways are processed
    p90, p95, p99 = plot_edge_bc_distributions(
        pathway_ids,
        save_path=pjoin(CRITICAL_EDGES_P, "edge_bc_distributions.png")
    )

    # Use the global 95th percentile as your threshold for marking interactions
    print(f"\nRecommended threshold for 'extreme' edges: {p95:.4f}")