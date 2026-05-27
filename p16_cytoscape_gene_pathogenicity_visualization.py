"""p16_cytoscape_gene_pathogenicity_visualization.py
===================================================
Produces publication-quality KEGG pathway network PNGs.

For each pathway two PNGs are written:
  <id>_network.png              – annotated version with title & full legend
  <id>_network_presentation.png – clean version: no title, bigger fonts,
                                   minimal legend (for slides)

Node encoding
-------------
  Color  – per-gene delta_means (green=benign → cream → mauve=pathogenic);
            only FDR-significant genes are colored; others are grey
  Size   – betweenness centrality computed from the KGML graph
  Stars  – druggability tier appended to gene label
             ★★★  Tier 1  clinically validated
             ★★   Tier 2  strong inferred
             ★    Tier 3A/3B  speculative

Output: plots/p16_pathway_networks/<pathway_id>_network[_presentation].png
"""

import os
import pickle
import time
from os.path import join as pjoin
from xml.etree import ElementTree as ET
import sys
import glob
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import requests

from definitions import KEGG_PATHWAY_METADATA_FILE, PLOTS_P, COLOR_MAP, DATA_P

# ─────────────────────────────────────────────────────────────────────────────
# Configuration
# ─────────────────────────────────────────────────────────────────────────────

BASE_P = "/cs/labs/dina/ophirmil12/PathwayAtlas2"
GENE_LEVEL_P = pjoin(BASE_P, "results", "gene_level_results")
DRUGGABILITY_FILE = pjoin(DATA_P, "papers", "NIHMS80906-supplement-Table_S1.xlsx")

# TARGET_PATHWAYS = ["hsa04722", "hsa04216", "hsa04137", "hsa04080", "hsa00470", "hsa01522", "hsa04950"]
TARGET_PATHWAYS = ["hsa04722"]

Q_THRESHOLD = 0.05

DPI = 600

# Custom diverging colormap: benign (green) → warm cream → pathogenic (mauve)
PATHWAY_CMAP = mcolors.LinearSegmentedColormap.from_list(
    "pathway_atlas",
    [COLOR_MAP["benign"], "#F0EAD6", COLOR_MAP["pathogenic"]],
    N=256,
)

MIN_NODE_SIZE = 1000
MAX_NODE_SIZE = 10000
NONSIG_NODE_SIZE = 80
NONSIG_COLOR = "#636363"

TIER_STARS = {
    "Tier 1":  "★★★",
    "Tier 2":  "★★",
    "Tier 3A": "★",
    "Tier 3B": "★",
}


# ─────────────────────────────────────────────────────────────────────────────
# KGML fetch and parse
# ─────────────────────────────────────────────────────────────────────────────

def fetch_kgml(pathway_id: str) -> ET.Element:
    url = f"https://rest.kegg.jp/get/{pathway_id}/kgml"
    r = requests.get(url, timeout=30)
    r.raise_for_status()
    return ET.fromstring(r.text)


def parse_kgml(root: ET.Element):
    """
    Returns:
      entries  – dict[entry_id -> {kegg_ids, label, x, y, type}]
      relations – list[(entry_id_a, entry_id_b, relation_type)]
    """
    entries = {}
    for entry in root.findall("entry"):
        eid = entry.get("id")
        etype = entry.get("type")
        name_raw = entry.get("name", "")
        kegg_ids = name_raw.split()

        gfx = entry.find("graphics")
        if gfx is None:
            continue
        raw_label = gfx.get("name", "")
        label = raw_label.split(",")[0].split("...")[0].strip()
        try:
            x = float(gfx.get("x", 0))
            y = float(gfx.get("y", 0))
        except (TypeError, ValueError):
            x, y = 0.0, 0.0

        entries[eid] = {
            "kegg_ids": kegg_ids,
            "label": label,
            "x": x,
            "y": -y,
            "type": etype,
        }

    relations = []
    for rel in root.findall("relation"):
        relations.append((rel.get("entry1"), rel.get("entry2"), rel.get("type", "")))

    return entries, relations


def build_graph(entries: dict, relations: list) -> nx.DiGraph:
    G = nx.DiGraph()
    for eid, info in entries.items():
        if info["type"] == "gene":
            G.add_node(eid, **info)
    gene_ids = set(G.nodes())
    for e1, e2, rtype in relations:
        if e1 in gene_ids and e2 in gene_ids:
            G.add_edge(e1, e2, rtype=rtype)
    return G


# ─────────────────────────────────────────────────────────────────────────────
# Data loading
# ─────────────────────────────────────────────────────────────────────────────

def load_gene_results(path: str ) -> pd.DataFrame:
    
    df = pd.read_csv(path)
    df["gene_id"] = df["gene_id"].str.strip()
    return df.set_index("gene_id")


def load_druggability() -> dict:
    """Returns {hgnc_name -> star_string} for all druggable genes."""
    df = pd.read_excel(DRUGGABILITY_FILE, sheet_name="Data",
                       usecols=["hgnc_names", "druggability_tier"])
    df = df.dropna(subset=["hgnc_names", "druggability_tier"])
    result = {}
    for _, row in df.iterrows():
        name  = str(row["hgnc_names"]).strip()
        tier  = str(row["druggability_tier"]).strip()
        stars = TIER_STARS.get(tier, "★")
        # keep the best tier if a gene appears more than once
        if name not in result or len(stars) > len(result[name]):
            result[name] = stars
    return result


# ─────────────────────────────────────────────────────────────────────────────
# Betweenness centrality
# ─────────────────────────────────────────────────────────────────────────────

def compute_betweenness(G: nx.DiGraph) -> dict:
    if len(G) == 0:
        return {}
    ug = G.to_undirected()
    try:
        return nx.betweenness_centrality(ug, normalized=True)
    except Exception:
        return {n: 0.0 for n in G.nodes()}


def scale_sizes(bc_dict: dict, node_list: list) -> np.ndarray:
    vals = np.array([bc_dict.get(n, 0.0) for n in node_list])
    vmax = vals.max()
    if vmax < 1e-9:
        return np.full(len(node_list), (MIN_NODE_SIZE + MAX_NODE_SIZE) / 2)
    norm = vals / vmax
    return MIN_NODE_SIZE + norm * (MAX_NODE_SIZE - MIN_NODE_SIZE)


# ─────────────────────────────────────────────────────────────────────────────
# Plotting
# ─────────────────────────────────────────────────────────────────────────────

def plot_pathway(
    pathway_id: str,
    pathway_name: str,
    cancer_type: str,
    G: nx.DiGraph,
    bc_dict: dict,
    gene_results: pd.DataFrame,
    druggability: dict,
    out_dir: str,
    q_threshold: float = Q_THRESHOLD,
    dpi: int = DPI,
    simple: bool = False,
):
    n_nodes = len(G)
    if n_nodes == 0:
        print(f"  [skip] empty graph for {pathway_id}")
        return

    # ── assign delta_means to each node ──────────────────────────────────────
    node_delta = {}
    for node in G.nodes():
        kegg_ids = G.nodes[node]["kegg_ids"]
        best = None
        for kid in kegg_ids:
            if kid in gene_results.index:
                row = gene_results.loc[kid]
                q  = row["q_value"]     if np.ndim(row["q_value"])    == 0 else float(row["q_value"].iloc[0])
                dm = row["delta_means"] if np.ndim(row["delta_means"]) == 0 else float(row["delta_means"].iloc[0])
                if best is None or q < best[0]:
                    best = (q, dm)
        if best is not None and best[0] < q_threshold:
            node_delta[node] = best[1]

    sig_nodes    = [n for n in G.nodes() if n in node_delta]
    nonsig_nodes = [n for n in G.nodes() if n not in node_delta]

    # ── color map ─────────────────────────────────────────────────────────────
    cmap = PATHWAY_CMAP
    delta_vals = np.array(list(node_delta.values())) if node_delta else np.array([0.0])
    vabs = max(abs(delta_vals).max(), 1e-6)
    norm = mcolors.Normalize(vmin=-vabs, vmax=vabs)
    sig_colors = [cmap(norm(node_delta[n])) for n in sig_nodes]

    # ── node sizes ────────────────────────────────────────────────────────────
    all_nodes = list(G.nodes())
    all_sizes_arr = scale_sizes(bc_dict, all_nodes)
    size_map = dict(zip(all_nodes, all_sizes_arr))

    sig_sizes    = np.array([size_map[n] for n in sig_nodes])    if sig_nodes    else np.array([])
    nonsig_sizes = np.array([size_map[n] for n in nonsig_nodes]) if nonsig_nodes else np.array([])

    # ── labels: sig nodes + top-BC nodes, with druggability stars ─────────────
    bc_vals   = np.array([bc_dict.get(n, 0.0) for n in G.nodes()])
    bc_cutoff = float(np.percentile(bc_vals, 80)) if len(bc_vals) else 0.0
    label_set = set(sig_nodes) | {n for n in G.nodes() if bc_dict.get(n, 0.0) >= bc_cutoff}

    def make_label(node):
        base  = G.nodes[node]["label"]
        stars = druggability.get(base, "")
        # return f"{base}\n{stars}" if stars else base
        return base

    labels = {n: make_label(n) for n in label_set if G.nodes[n]["label"]}

    # ── positions: normalise KGML pixel coords to [0,1]×[0,1] ───────────────
    # KGML uses absolute pixel coords (~800-unit range); without normalisation
    # nodes occupy only a tiny fraction of the axes and appear very small.
    # We scale x and y independently so the full canvas is always utilised
    # while preserving the relative layout.
    raw_xs = np.array([G.nodes[n]["x"] for n in G.nodes()])
    raw_ys = np.array([G.nodes[n]["y"] for n in G.nodes()])
    x_span = max(raw_xs.max() - raw_xs.min(), 1.0)
    y_span = max(raw_ys.max() - raw_ys.min(), 1.0)
    pos = {
        n: (
            (G.nodes[n]["x"] - raw_xs.min()) / x_span,
            (G.nodes[n]["y"] - raw_ys.min()) / y_span,
        )
        for n in G.nodes()
    }

    # ── figure: fixed large size so nodes fill the space ─────────────────────
    aspect = x_span / y_span
    fig_h = 18
    fig_w = max(20, min(32, fig_h * aspect)) + 3   # +3 for colorbar margin
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.set_facecolor("#F5F5F5")
    fig.patch.set_facecolor("white")

    # draw edges — bottom layer (zorder set on returned collection)
    edge_coll = nx.draw_networkx_edges(
        G, pos, ax=ax,
        edge_color="#636363", alpha=0.6,
        arrows=False, arrowsize=18, width=2,
        connectionstyle="arc3,rad=0.05",
    )
    if edge_coll is not None:
        edge_coll.set_zorder(1)

    # non-significant nodes — above edges, faded so colored nodes dominate
    if nonsig_nodes:
        nc = nx.draw_networkx_nodes(
            G, pos, nodelist=nonsig_nodes, ax=ax,
            node_color=NONSIG_COLOR,
            node_size=nonsig_sizes,
            alpha=0.3,
        )
        nc.set_zorder(2)

    # significant nodes — topmost layer, full opacity, white border
    if sig_nodes:
        nc = nx.draw_networkx_nodes(
            G, pos, nodelist=sig_nodes, ax=ax,
            node_color=sig_colors,
            node_size=sig_sizes,
            alpha=1.0,
            edgecolors="white",
            linewidths=1.5,
        )
        nc.set_zorder(3)

    # labels
    if labels:
        fontsize = max(16.0, min(14.0, 160.0 / max(len(labels), 1)))
        if simple:
            fontsize = max(16.0, min(16.0, 200.0 / max(len(labels), 1)))
        nx.draw_networkx_labels(
            G, pos, labels=labels, ax=ax,
            font_size=fontsize, font_color="#111111",
            font_weight="bold",
        )

    # ── colorbar ──────────────────────────────────────────────────────────────
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, fraction=0.025, pad=0.02, aspect=25)
    cbar.set_label("𝔼(Cancer)-𝔼(Background)", fontsize=25 if simple else 16)
    cbar.ax.tick_params(labelsize=14 if simple else 12)

    # ── legend ────────────────────────────────────────────────────────────────
    leg_fs       = 25 if simple else 18
    leg_title_fs = 25 if simple else 18

    legend_handles = []

    if not simple:
        max_bc = max(bc_dict.values()) if bc_dict else 1.0
        for frac, lbl in [(0.25, "25%"), (0.5, "50%"), (1.0, "100%")]:
            s = MIN_NODE_SIZE + frac * (MAX_NODE_SIZE - MIN_NODE_SIZE)
            r = np.sqrt(s / np.pi) * 0.012


    # druggability star key (text-only entries)
    legend_handles += [
        mpatches.Patch(color="none", label="★★★ Clinically validated druggability"),
        mpatches.Patch(color="none", label="★★ Strong inferred druggability"),
        mpatches.Patch(color="none", label="★ Speculative druggability"),
    ]

    # ax.legend(
    #     handles=legend_handles,
    #     title="Node size ∝ betweenness centrality",
    #     loc="lower right", fontsize=leg_fs, title_fontsize=leg_title_fs,
    #     framealpha=0.88, ncol=1, borderpad=1.2,
    #     handlelength=2.0, handleheight=1.4,
    # )

    if not simple:
        short_name = pathway_name.split(" - Homo")[0]
        ax.set_title(
            f"{pathway_id} · {short_name}\n"
            f"{len(sig_nodes)} / {n_nodes} genes FDR-significant ({cancer_type})",
            fontsize=16, fontweight="bold", pad=16,
        )

    ax.axis("off")
    plt.tight_layout()

    os.makedirs(out_dir, exist_ok=True)
    suffix   = "_presentation" if simple else ""
    # out_path = pjoin(out_dir, f"{pathway_id}_network{suffix}_no_legend.png")
    out_path = pjoin(out_dir, f"{pathway_id}_network{suffix}_no_stars.png")
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved → {out_path}  [{len(sig_nodes)} sig / {n_nodes} total]")


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    args = sys.argv[1:]
    if len(args) < 1:
        print("Usage: python -u p9_patient_scores_and_survival_per_cancer.py '$SLURM_ARRAY_TASK_ID'")
        sys.exit(1)

    cancer_results_files = glob.glob(pjoin(GENE_LEVEL_P, "*.csv"))
    if not cancer_results_files:
        print(f"No cancer mutations CSV files found in the specified directory: {GENE_LEVEL_P}")
        sys.exit(1)

    index = int(args[0])
    if index < 0 or index >= len(cancer_results_files):
        print(f"Index {index} is out of range. There are only {len(cancer_results_files)} cancer mutation files.")
        sys.exit(1)

    cancer_gene_distances_file = sorted(cancer_results_files)[index]

    # the file name is in the format "{cancer_type}_gene_distances.csv"
    cancer_type = os.path.basename(cancer_gene_distances_file).split('_gene_distances.csv')[0]
    
    out_dir = pjoin(PLOTS_P, f"p16_pathway_networks/{cancer_type}")
    os.makedirs(out_dir, exist_ok=True)

    print(f"Loading {cancer_type} gene results…")
    gene_results = load_gene_results(cancer_gene_distances_file)
    n_sig = int((gene_results["q_value"] < Q_THRESHOLD).sum())
    print(f"  {len(gene_results):,} genes total; {n_sig:,} with q < {Q_THRESHOLD}")

    print("Loading druggability data…")
    druggability = load_druggability()
    print(f"  {len(druggability):,} druggable genes loaded")

    with open(KEGG_PATHWAY_METADATA_FILE, "rb") as f:
        pathway_meta = pickle.load(f)

    for pathway_id in TARGET_PATHWAYS:
        meta = pathway_meta.get(pathway_id, {})
        pathway_name = meta.get("name", pathway_id)
        short = pathway_name.split(" - Homo")[0]
        print(f"\n{'─'*60}")
        print(f"{pathway_id}: {short}")
        print(f"{'─'*60}")

        print("  Fetching KGML from KEGG…")
        try:
            root = fetch_kgml(pathway_id)
        except Exception as exc:
            print(f"  [error] KGML fetch failed: {exc}")
            continue
        time.sleep(0.4)

        entries, relations = parse_kgml(root)
        gene_count = sum(1 for v in entries.values() if v["type"] == "gene")
        print(f"  {gene_count} gene entries, {len(relations)} relations in KGML")

        G = build_graph(entries, relations)
        print(f"  Graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

        print("  Computing betweenness centrality…")
        bc_dict = compute_betweenness(G)
        top_bc = sorted(bc_dict.items(), key=lambda x: -x[1])[:3]
        for nid, bval in top_bc:
            print(f"    {G.nodes[nid]['label']:14s}  BC={bval:.4f}")

        pathway_labels = {G.nodes[n]["label"] for n in G.nodes()}
        n_drug = sum(1 for lbl in pathway_labels if lbl in druggability)
        print(f"  Druggable genes in pathway: {n_drug}")

        shared = dict(
            pathway_id=pathway_id, pathway_name=pathway_name, cancer_type=cancer_type,
            G=G, bc_dict=bc_dict, gene_results=gene_results,
            druggability=druggability, out_dir=out_dir,
            q_threshold=Q_THRESHOLD, dpi=DPI,
        )
        # plot_pathway(**shared, simple=False)
        plot_pathway(**shared, simple=True)

    print(f"\nDone. PNGs written to {out_dir}/")


if __name__ == "__main__":
    main()
