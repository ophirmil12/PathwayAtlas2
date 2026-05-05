from __future__ import annotations
import numpy as np
from definitions import *
from collections import Counter
import pickle
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import FancyBboxPatch
import seaborn as sns
from collections import Counter
from pathlib import Path
import csv
from kegg_api import KeggGene
from kegg_api import KeggApi
import sys
from figures.consensus.mappings import s2_pathway_to_hsa


 
def compute_overlap_score(pathways: dict[str, set[str]]) -> float:
    """
    Compute the mean number of pathways each gene appears in.
    A score of 1.0 means perfectly non-redundant. Higher = more redundant.
 
    This is the primary metric from Stoney et al. (2018).
    """
    gene_counts = Counter()
    for genes in pathways.values():
        for g in genes:
            gene_counts[g] += 1
 
    if len(gene_counts) == 0:
        return 0.0
    return sum(gene_counts.values()) / len(gene_counts)


def hitting_set_cover_with_provenance(
    pathway_metadata_dict: dict[str, dict],
    consensus_pathways: set[str] | None = None,
    gene_coverage: float = 1.0,
) -> tuple[list[str], dict[str, dict]]:
    """
    Hitting set cover that tracks provenance for omitted pathways.
 
    Parameters
    ----------
    pathways : dict
        pathway_id -> set of gene identifiers.
    protected_pathways : set[str], optional
        Pathway IDs that are forced into the cover set before the greedy
        algorithm runs. Their genes are marked as covered immediately.
        Use this for pathways you know a priori are biologically important
        (e.g. core cancer signaling pathways) and must not be dropped
        regardless of their overlap profile.
    gene_coverage : float
        Fraction of genes to cover (GC parameter).
 
    Returns
    -------
    cover_set : list[str]
        Selected pathway IDs (same as original hitting_set_cover).
    omitted_info : dict
        For each omitted pathway, details on what replaced it.
        Keys: pathway IDs that were NOT selected.
        Values: dict with "total_genes", "covered_by", "uncovered_genes".
    """
    if consensus_pathways is None:
        consensus_pathways = set()

    # Validate protected pathways exist in the input
    missing = consensus_pathways - set(pathways.keys())
    if missing:
        print(f"Warning: {len(missing)} protected pathway IDs not found in input "
              f"(already removed by size/category filter?): {missing}")
        consensus_pathways = consensus_pathways - missing
 
    universe = set()
    for metdata in pathway_metadata_dict.values():
        genes = set(metdata["genes_ids"])
        universe |= genes
 
    target = int(np.ceil(len(universe) * gene_coverage))
 
    # Gene frequencies and values
    gene_freq = Counter()
    for metdata in pathway_metadata_dict.values():
        genes = metdata["genes_ids"]
        for g in genes:
            gene_freq[g] += 1
    gene_value = {g: 1.0 / f for g, f in gene_freq.items()}
 
    uncovered = set(universe)
    covered = set()
    cover_set = []
    cover_set_lookup = set()
 
    # Track which cover_set pathway first covered each gene
    gene_covered_by = {}  # gene -> pathway_id that first covered it

    # ---- Seed protected pathways first ----
    for pid in sorted(consensus_pathways):
        cover_set.append(pid)
        cover_set_lookup.add(pid)
        newly_covered = pathways[pid] & uncovered
        for g in newly_covered:
            gene_covered_by[g] = pid
        covered |= newly_covered
        uncovered -= newly_covered
 
    # ---- Greedy algorithm over remaining pathways ----
    while len(covered) < target and uncovered:
        best_pid = None
        best_score = -1.0
 
        for pid, metdata in pathway_metadata_dict.items():
            genes = set(metdata["genes_ids"])
            if pid in cover_set_lookup:
                continue
            size = len(genes)
            if size == 0:
                continue
            uncov_genes = genes & uncovered
            if not uncov_genes:
                continue
            value_sum = sum(gene_value[g] for g in uncov_genes)
            score = value_sum / size
            if score > best_score:
                best_score = score
                best_pid = pid
 
        if best_pid is None:
            break
 
        cover_set.append(best_pid)
        cover_set_lookup.add(best_pid)
 
        # Record which genes this pathway newly covers
        newly_covered = set(pathway_metadata_dict[best_pid]["genes_ids"]) & uncovered
        for g in newly_covered:
            gene_covered_by[g] = best_pid
        covered |= newly_covered
        uncovered -= newly_covered
 
    # Now build provenance for every omitted pathway
    omitted_info = {}
    for pid, metdata in pathway_metadata_dict.items():
        genes = metdata["genes_ids"]
        if pid in cover_set_lookup:
            continue
 
        # For each gene in this omitted pathway, find which cover_set
        # pathway is responsible for it
        covered_by_map = {}  # cover_pid -> set of shared genes
        uncovered_genes = set()
 
        for g in genes:
            if g in gene_covered_by:
                cpid = gene_covered_by[g]
                covered_by_map.setdefault(cpid, set()).add(g)
            else:
                # Gene wasn't covered at all (GC < 1.0)
                uncovered_genes.add(g)
 
        # Build the "covered_by" summary, sorted by contribution
        covered_by = {}
        for cpid, shared in sorted(
            covered_by_map.items(), key=lambda x: -len(x[1])
        ):
            covered_by[cpid] = {
                "name": pathway_metadata_dict.get(str(cpid), {}).get('name', 'Unknown').split(' - Homo')[0],
                "shared_genes": len(shared),
                "percent_of_omitted": round(len(shared) / len(genes) * 100, 1),
                "gene_ids": shared,
            }
 
        omitted_info[pid] = {
            "name": pathway_metadata_dict.get(str(pid), {}).get('name', 'Unknown').split(' - Homo')[0],
            "total_genes": len(genes),
            "genes_accounted_for": len(genes) - len(uncovered_genes),
            "percent_accounted_for": round(
                (len(genes) - len(uncovered_genes)) / max(len(genes), 1) * 100, 1
            ),
            "n_replacing_pathways": len(covered_by),
            "covered_by": covered_by,
        }
 
    return cover_set, omitted_info


 
def save_omitted_pathways_report_csv(
    omitted_info: dict[str, dict],
    save_path: str,
):
    """
    Save the provenance data as a CSV for further analysis.
 
    Columns: omitted_id, omitted_name, omitted_size, replacer_id,
             replacer_name, shared_genes, percent_of_omitted
 
    Each row is one (omitted_pathway, replacer_pathway) pair.
    """ 
    with open(save_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "omitted_id", "omitted_name", "omitted_size",
            "percent_accounted_for", "replacer_id", "replacer_name",
            "shared_genes", "percent_of_omitted", 
        ])
 
        for pid, info in sorted(omitted_info.items(), key=lambda x: -x[1]["total_genes"]):
            for cpid, cinfo in info["covered_by"].items():
                writer.writerow([
                    pid, info["name"], info["total_genes"],
                    info["percent_accounted_for"],
                    cpid, cinfo["name"],
                    cinfo["shared_genes"], cinfo["percent_of_omitted"],
                ])
 
    print(f"Saved provenance CSV to {save_path}")

def save_omitted_genes_report_csv(
    pathways: dict[str, set[str]],
    cover_set: list[str],
    omitted_info: dict[str, dict],
    save_path: str,
):
    """
    Save a CSV of genes that were dropped by the hitting set cover algorithm —
    i.e. genes that do not appear in any selected pathway.

    Columns: gene_id, n_omitted_pathways, omitted_pathway_ids, omitted_pathway_names

    Each row is one dropped gene.
    """
    # Genes present in at least one selected pathway
    covered_genes: set[str] = set()
    for pid in cover_set:
        covered_genes |= pathways[pid]

    # Reverse map: dropped gene -> list of (pid, name) from omitted pathways
    gene_to_omitted: dict[str, list[tuple[str, str]]] = {}
    for pid, info in omitted_info.items():
        for gene_id in pathways[pid]:
            if gene_id not in covered_genes:
                gene_to_omitted.setdefault(gene_id, []).append((pid, info["name"]))

    with open(save_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "gene_id", "gene_name", "n_omitted_pathways",
            "omitted_pathway_ids", "omitted_pathway_names",
        ])
        for gene_id, entries in sorted(gene_to_omitted.items(), key=lambda x: -len(x[1])):
            kegg_gene = KeggGene(gene_id)
            gene_names = kegg_gene.ref_names
            if not gene_names:
                print(f"Warning: no gene names available for {gene_id} gene.")
                gene_name = np.nan
            else:
                gene_name = gene_names[0]
            pids = ";".join(pid for pid, _ in entries)
            names = ";".join(name for _, name in entries)
            writer.writerow([gene_id, gene_name, len(entries), pids, names])

    print(f"Saved dropped-genes CSV to {save_path} ({len(gene_to_omitted)} genes)")


def _jaccard_matrix(pathways: dict[str, set[str]]) -> np.ndarray:
    """Compute pairwise Jaccard matrix as numpy array."""
    pids = list(pathways.keys())
    n = len(pids)
    mat = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            a = pathways[pids[i]]
            b = pathways[pids[j]]
            inter = len(a & b)
            union = len(a | b)
            jac = inter / union if union > 0 else 0.0
            mat[i, j] = jac
            mat[j, i] = jac
    return mat
 
 
def _upper_triangle(mat: np.ndarray) -> np.ndarray:
    """Extract upper triangle values (excluding diagonal)."""
    mask = np.triu(np.ones_like(mat, dtype=bool), k=1)
    return mat[mask]
 

 
# ============================================================================
# Color scheme
# ============================================================================
 
# Custom colormap: light gray -> warm amber -> red -> dark red
JACCARD_CMAP = LinearSegmentedColormap.from_list(
    "jaccard",
    ["#FFFFFF", "#FAC775", "#EF9F27", "#E24B4A", "#791F1F"],
    N=256,
)
 
COLOR_BEFORE = COLOR_MAP['teal']
COLOR_AFTER = COLOR_MAP['dark-blue']
COLOR_MUTED = COLOR_MAP['grey']
COLOR_BG = COLOR_MAP['bg']
 
# ============================================================================
# 2. Jaccard distribution histogram
# ============================================================================
 

def plot_redundancy_summary(
    pathways_before: dict[str, set[str]],
    pathways_after: dict[str, set[str]],
    stats: dict,
    target_gene_coverage: float = 1,
    save_path: str | Path | None = None,
    figsize: tuple = (14, 8),
    dpi: int = 300,
    vmax: float = 0.8,
) -> plt.Figure:
    """
    Combined 4-panel figure:
      - Top left:  Jaccard heatmap BEFORE
      - Top right: Jaccard heatmap AFTER
      - Bottom left: Jaccard distribution histogram
      - Bottom right: Summary statistics table
 
    This is the main figure for your paper/thesis.
 
    Parameters
    ----------
    pathways_before : dict
        Pathways before set cover (after category pre-filtering).
    pathways_after : dict
        Pathways after set cover.
    stats : dict, optional
        Output of reduce_kegg_redundancy(). If None, computed here.
    save_path : str or Path, optional
        If provided, saves the figure. Supports .pdf, .png, .svg.
    figsize : tuple
        Figure size in inches.
    dpi : int
        Resolution for raster output.
    vmax : float
        Colormap ceiling for heatmaps.
 
    Returns
    -------
    matplotlib.figure.Figure
    """
    # ---- Compute matrices ----
    mat_before = _jaccard_matrix(pathways_before)
    mat_after = _jaccard_matrix(pathways_after)
    np.fill_diagonal(mat_before, 0)
    np.fill_diagonal(mat_after, 0)
 
    vals_before = _upper_triangle(mat_before + np.diag(np.ones(len(mat_before))))
    vals_after = _upper_triangle(mat_after + np.diag(np.ones(len(mat_after))))
    # Recompute without diagonal
    vals_before = _upper_triangle(mat_before)
    vals_after = _upper_triangle(mat_after)

    # ---- Create figure ----
    fig = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(
        2, 3, figure=fig,
        width_ratios=[1, 1, 0.05],
        height_ratios=[1, 0.85],
        hspace=0.35, wspace=0.3,
    )
 
    # ---- Panel A: Heatmap before ----
    ax_hm1 = fig.add_subplot(gs[0, 0])
    im1 = ax_hm1.imshow(mat_before, cmap=JACCARD_CMAP, vmin=0, vmax=vmax,
                          aspect="auto", interpolation="nearest")
    ax_hm1.set_title(f"Pairwise Jaccard Similiarity - Before",
                     fontsize=11, fontweight=500, pad=8)
    ax_hm1.set_xlabel("Pathway index", fontsize=9)
    ax_hm1.set_ylabel("Pathway index", fontsize=9)
    ax_hm1.tick_params(labelsize=7)
    ax_hm1.text(-0.12, 1.08, "A", transform=ax_hm1.transAxes,
                fontsize=16, fontweight=700, va="top")
 
    # ---- Panel B: Heatmap after ----
    ax_hm2 = fig.add_subplot(gs[0, 1])
    im2 = ax_hm2.imshow(mat_after, cmap=JACCARD_CMAP, vmin=0, vmax=vmax,
                          aspect="auto", interpolation="nearest")
    ax_hm2.set_title(f"Pairwise Jaccard Similiarity - After",
                     fontsize=11, fontweight=500, pad=8)
    ax_hm2.set_xlabel("Pathway index", fontsize=9)
    ax_hm2.set_ylabel("Pathway index", fontsize=9)
    ax_hm2.tick_params(labelsize=7)
    ax_hm2.text(-0.12, 1.08, "B", transform=ax_hm2.transAxes,
                fontsize=16, fontweight=700, va="top")
 
    # ---- Shared colorbar ----
    cax = fig.add_subplot(gs[0, 2])
    cbar = fig.colorbar(im2, cax=cax)
    cbar.set_label("Jaccard similarity", fontsize=9)
    cbar.ax.tick_params(labelsize=7)
 
    # ---- Panel C: Histogram ----
    ax_hist = fig.add_subplot(gs[1, 0])
    bins = np.linspace(0, 1, 41)
 
    bin_width = bins[1] - bins[0]
    bar_width = bin_width * 0.45
    bin_centers = (bins[:-1] + bins[1:]) / 2

    counts_before, _ = np.histogram(vals_before, bins=bins, density=True)
    counts_after, _ = np.histogram(vals_after, bins=bins, density=True)

    ax_hist.bar(bin_centers - bar_width / 2, counts_before, width=bar_width,
                alpha=0.8, color=COLOR_BEFORE, edgecolor="white", linewidth=0.4,
                label=f"Before ({len(pathways_before)})")
    ax_hist.bar(bin_centers + bar_width / 2, counts_after, width=bar_width,
                alpha=0.8, color=COLOR_AFTER, edgecolor="white", linewidth=0.4,
                label=f"After ({len(pathways_after)})")
 
    ax_hist.set_xlabel("Pairwise Jaccard similarity", fontsize=10)
    ax_hist.set_ylabel("Density", fontsize=10)
    ax_hist.legend(fontsize=9, frameon=False, loc="upper right")
    ax_hist.spines["top"].set_visible(False)
    ax_hist.spines["right"].set_visible(False)
    ax_hist.tick_params(labelsize=8)
    ax_hist.text(-0.12, 1.08, "C", transform=ax_hist.transAxes,
                 fontsize=16, fontweight=700, va="top")
 
    # ---- Panel D: Summary stats table ----
    ax_stats = fig.add_subplot(gs[1, 1])
    ax_stats.axis("off")
    ax_stats.text(-0.12, 1.08, "D", transform=ax_stats.transAxes,
                  fontsize=16, fontweight=700, va="top")
 
    n_high_before = int(np.sum(vals_before > 0.3))
    n_high_after = int(np.sum(vals_after > 0.3))
    max_jac_before = float(np.max(vals_before)) if len(vals_before) > 0 else 0
    max_jac_after = float(np.max(vals_after)) if len(vals_after) > 0 else 0
    mean_jac_before = float(np.mean(vals_before)) if len(vals_before) > 0 else 0
    mean_jac_after = float(np.mean(vals_after)) if len(vals_after) > 0 else 0
 
    table_data = [
        ["Pathways", str(stats["after_prefilter"]), str(stats["after_set_cover"])],
        ["Unique genes", str(stats["genes_before"]), str(stats["genes_after"])],
        ["Gene coverage", "100%", f"{stats['gene_coverage_achieved']}%"],
        ["Overlap score", f"{stats['overlap_before']:.2f}", f"{stats['overlap_after']:.2f}"],
        ["Mean Jaccard", f"{mean_jac_before:.4f}", f"{mean_jac_after:.4f}"],
        ["Max Jaccard", f"{max_jac_before:.3f}", f"{max_jac_after:.3f}"],
        ["Pairs J > 0.3", f"{n_high_before:,}", f"{n_high_after:,}"],
        ["Redundancy reduction", "", f"{stats['redundancy_reduction_pct']}%"],
    ]
 
    table = ax_stats.table(
        cellText=table_data,
        colLabels=["Metric", "Before", "After"],
        loc="center",
        cellLoc="center",
        colColours=[COLOR_BG, COLOR_BG, COLOR_BG],
    )
 
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1.0, 1.6)
 
    # Style the table
    for (row, col), cell in table.get_celld().items():
        cell.set_edgecolor("#D3D1C7")
        cell.set_linewidth(0.5)
        if row == 0:
            cell.set_text_props(fontweight=500)
            cell.set_facecolor(COLOR_BG)
        else:
            cell.set_facecolor("white")
        # Highlight the "After" column values
        if col == 2 and row > 0:
            cell.set_text_props(color=COLOR_AFTER, fontweight=500)
 
    fig.suptitle(
        f"KEGG Pathway Dataset Redundancy Reduction via Hitting Set Cover (gene coverage target = {target_gene_coverage * 100}%)",
        fontsize=13, fontweight=500, y=0.98,
    )
 
    if save_path:
        fig.savefig(save_path, dpi=dpi, bbox_inches="tight")
        print(f"Saved summary figure to {save_path}")

def plot_pathway_sizes(
    pathways_before: dict[str, set[str]],
    pathways_after: dict[str, set[str]],
    save_path: str | Path | None = None,
    figsize: tuple = (7, 4),
    dpi: int = 300,
) -> plt.Figure:
    """
    Boxplot + strip plot of pathway sizes before and after.
    Shows that set cover doesn't just keep the biggest pathways.
    """
    sizes_before = [len(g) for g in pathways_before.values()]
    sizes_after = [len(g) for g in pathways_after.values()]
 
    df = pd.concat([
        pd.DataFrame({"Pathway size (genes)": sizes_before,
                       "Set": f"Before (n={len(sizes_before)})"}),
        pd.DataFrame({"Pathway size (genes)": sizes_after,
                       "Set": f"After (n={len(sizes_after)})"}),
    ])
 
    fig, ax = plt.subplots(figsize=figsize)
 
    sns.boxplot(
        data=df, x="Set", y="Pathway size (genes)", ax=ax,
        palette=[COLOR_BEFORE, COLOR_AFTER],
        width=0.4, linewidth=0.8, fliersize=3,
        boxprops=dict(alpha=0.6),
    )
    sns.stripplot(
        data=df, x="Set", y="Pathway size (genes)", ax=ax,
        palette=[COLOR_BEFORE, COLOR_AFTER],
        size=3, alpha=0.4, jitter=0.15,
    )
 
    ax.set_xlabel("")
    ax.set_ylabel("Pathway size (number of genes)", fontsize=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(labelsize=9)
 
    # Add median annotations
    for i, sizes in enumerate([sizes_before, sizes_after]):
        med = np.median(sizes)
        ax.annotate(
            f"median = {med:.0f}",
            xy=(i, med), xytext=(i + 0.3, med + max(sizes) * 0.05),
            fontsize=8, color=COLOR_MUTED,
            arrowprops=dict(arrowstyle="-", color=COLOR_MUTED, lw=0.5),
        )
 
    fig.tight_layout()
 
    if save_path:
        fig.savefig(save_path, dpi=dpi, bbox_inches="tight")
        print(f"Saved size distribution to {save_path}")


if __name__ == '__main__':

    #target_gene_coverage = 1.0  
    target_gene_coverage = 0.95

    # Load pathway metadata
    with open(KEGG_PATHWAY_METADATA_FILE, 'rb') as f:
        pathway_metadata = pickle.load(f)

    pathways = {pid: set(d["genes_ids"]) for pid, d in pathway_metadata.items()}

    print(f"Filtering by pathway size, number of pathways before: {len(pathway_metadata)}")
    for pid, d in list(pathway_metadata.items()):
        genes = d["genes_ids"]
        pathway_size = len(genes)
        if pathway_size > MAX_PATHWAY_SIZE or pathway_size < MIN_PATHWAY_SIZE:
            pathway_metadata.pop(pid)
    print(f"Finished filtering by size, number of pathways after: {len(pathway_metadata)}")

    # print("Filtering by categories...")
    # # Fetch the KEGG pathway BRITE hierarchy using the Kegg API function 
    # brite = KeggApi.fetch_kegg_pathway_brite(KEGG_BRITE_JSON_FILE)

    # remove the ids from the metadata dict
    for pid in EXCLUDE_KEGG_IDS:
        if pid in pathway_metadata.keys():
            pathway_metadata.pop(pid)

    filtered = {pid: set(d["genes_ids"]) for pid, d in pathway_metadata.items()}
    overlap_before = compute_overlap_score(filtered)

    print(f"Finished filtering by categories, number of pathways after: {len(pathway_metadata)}")

    print("Filtering using hitting set cover algorithm...")

    # run set cover algorithm

    consensus_s2_csv_file = pjoin(FIGURES_P, "consensus/Table S2.csv")
    consensus_s2_df = pd.read_csv(consensus_s2_csv_file)
    consensus_pathways = set()

    for pathway in consensus_s2_df["Pathway"].unique():
        if not pathway or pathway not in s2_pathway_to_hsa.keys():
            continue
        pathway_ids = s2_pathway_to_hsa[pathway]
        if not pathway_ids:
            continue
        consensus_pathways |= set([pid for pid in pathway_ids if pid in filtered.keys()])

    with open(pjoin(KEGG_P, "cancer_hallmark_kegg_pathways.csv")) as f:
        for row in csv.DictReader(f):
            consensus_pathways.add(row["kegg_id"])

    print(f"Defined {len(consensus_pathways)} consensus pathways.")

    cover_set, omitted_info = hitting_set_cover_with_provenance(pathway_metadata, consensus_pathways=consensus_pathways, gene_coverage=target_gene_coverage)
    save_omitted_pathways_report_csv(omitted_info, pjoin(PATHWAY_FILTERING_P, f"omitted_pathways_info_{target_gene_coverage}_coverage.csv"))
    save_omitted_genes_report_csv(filtered, cover_set, omitted_info, pjoin(PATHWAY_FILTERING_P, f"omitted_genes_{target_gene_coverage}_coverage.csv"))

    reduced_metadata = {pid: d for pid ,d in pathway_metadata.items() if pid in cover_set}

    print(f"Finished filtering by hitting set cover algorithm, number of pathways after: {len(reduced_metadata)}")

    # save updated metadata to pickle
    with open(FILTERED_KEGG_PATHWAY_METADATA_FILE, 'wb') as f:
        pickle.dump(reduced_metadata, f)
    
    print(f"Updated pathway metadata dictionary saved to {FILTERED_KEGG_PATHWAY_METADATA_FILE}.")

    print("Plotting results...")

    reduced = {pid: set(d["genes_ids"]) for pid, d in reduced_metadata.items()}
    overlap_after = compute_overlap_score(reduced)

    universe_before = set()
    for g in filtered.values():
        universe_before |= g
    universe_after = set()
    for g in reduced.values():
        universe_after |= g
    actual_coverage = len(universe_after) / len(universe_before) if universe_before else 0

    stats = {
        "original_pathways": len(pathways),
        "after_prefilter": len(filtered),
        "after_set_cover": len(reduced),
        "pathways_removed": len(filtered) - len(reduced),
        "overlap_before": round(overlap_before, 2),
        "overlap_after": round(overlap_after, 2),
        "redundancy_reduction_pct": round(
            (1 - (overlap_after - 1) / max(overlap_before - 1, 1e-9)) * 100, 1
        ),
        "genes_before": len(universe_before),
        "genes_after": len(universe_after),
        "gene_coverage_achieved": round(actual_coverage * 100, 1),
        "gene_coverage_target": target_gene_coverage * 100,
    }

    plot_redundancy_summary(
        filtered, reduced, stats, target_gene_coverage,
        save_path=pjoin(PATHWAY_FILTERING_P, f"redundancy_summary_{target_gene_coverage}_coverage.png"),
    )

    plot_pathway_sizes(
        pathways, reduced,
        save_path=pjoin(PATHWAY_FILTERING_P, f"pathway_sizes_{target_gene_coverage}_coverage.png"),
    )
