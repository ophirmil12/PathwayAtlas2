"""
scatter_s2_vs_delta_means.py
────────────────────────────
Generates the two key validation plots:
  Plot 1: Scatter plot of Consensus Score vs. Mean Δ-means
  Plot 2: Violin plot of Driver-mapped vs. Other pathway clusters (Max Δ-means)

Features:
- Expands "driver-mapped" pathways to include all pathways in their cluster.
- Scatter plot uses the MEAN Δ-means for the relevant pathways.
- Violin plot uses the MAXIMUM Δ-means per cluster.
- Adds jitter to the scatter plot for better point density visualization.
- Exports only 600 DPI PNGs.

Working directory: ./figures/consensus
Run as: python figures/consensus/scatter_s2_vs_delta_means.py
"""

from pathlib import Path
import sys
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats

# ── project root ──────────────────────────────────────────────────────────────
PROJECT_ROOT = Path(__file__).parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from definitions import (
    RESULTS_DISTANCES_P, 
    CANCER_FULLNAME, 
    KEGG_PATHWAY_CLUSTERING_P
)

# ── local (same folder as this script) ───────────────────────────────────────
HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
from mappings import s2_cancer_to_mine, s2_pathway_to_hsa

# ─────────────────────────────────────────────────────────────────────────────
# 0.  Configuration
# ─────────────────────────────────────────────────────────────────────────────
S2_CSV      = HERE / "Table S2.csv"
S1_CSV      = HERE / "Table S1.csv"
DIST_DIR    = Path(RESULTS_DISTANCES_P)
CLUSTER_CSV = Path(KEGG_PATHWAY_CLUSTERING_P) / "pathway_clusters.csv"
OUT_DIR     = HERE / "output_plots"
OUT_DIR.mkdir(parents=True, exist_ok=True)

DELTA_MEANS_COL = "delta_means"
PATHWAY_ID_COL  = "pathway"

# Pathway category → colour (one colour per S2 category)
CATEGORY_COLOURS = {
    "Apoptosis":                       "#E07B6A",
    "Cell cycle":                      "#6A9FCC",
    "Chromatin SWI/SNF complex":       "#9B7BB8",
    "Chromatin histone modifiers":     "#B89B7B",
    "Chromatin other":                 "#D4A5C9",
    "Epigenetics DNA modifiers":       "#7BB89B",
    "Genome integrity":                "#CC9B6A",
    "Histone modification":            "#6AB8CC",
    "Immune signaling":                "#B86A7B",
    "MAPK signaling":                  "#9BB87B",
    "Metabolism":                      "#B8B86A",
    "NFKB signaling":                  "#6A7BB8",
    "NOTCH signaling":                 "#CC6A9B",
    "Other":                           "#AAAAAA",
    "Other signaling":                 "#CCCCCC",
    "PI3K signaling":                  "#6ACC9B",
    "Protein homeostasis/ubiquitination": "#CC7B6A",
    "RNA abundance":                   "#7B9BCC",
    "RTK signaling":                   "#CC9B9B",
    "Splicing":                        "#9BCCB8",
    "TGFB signaling":                  "#B8CC9B",
    "TOR signaling":                   "#CC6A6A",
    "Transcription factor":            "#9B9B9B",
    "Wnt/B-catenin signaling":         "#6AB8B8",
}


# ─────────────────────────────────────────────────────────────────────────────
# 1.  Load Clusters & Expand the Pathway Mappings
# ─────────────────────────────────────────────────────────────────────────────
if not CLUSTER_CSV.exists():
    raise FileNotFoundError(f"Cluster CSV not found: {CLUSTER_CSV}. Please run the clustering script first.")

cluster_df = pd.read_csv(CLUSTER_CSV)
p2c = dict(zip(cluster_df['pathway_id'], cluster_df['cluster_id']))
c2p = cluster_df.groupby('cluster_id')['pathway_id'].apply(list).to_dict()

expanded_s2_pathway_to_hsa = {}
global_expanded_s2_hsa_ids = set()

for s2_cat, hsa_list in s2_pathway_to_hsa.items():
    expanded = set()
    if hsa_list:  # Safely handle None values
        for hsa in hsa_list:
            expanded.add(hsa)
            cid = p2c.get(hsa)
            if cid is not None:
                expanded.update(c2p[cid])
            
    expanded_s2_pathway_to_hsa[s2_cat] = expanded
    global_expanded_s2_hsa_ids.update(expanded)

print(f"Loaded clusters: mapped S2 drivers to {len(global_expanded_s2_hsa_ids)} total KEGG pathways via cluster expansion.")


# ─────────────────────────────────────────────────────────────────────────────
# 2.  Load Table S1 & Table S2
# ─────────────────────────────────────────────────────────────────────────────
def load_s2(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, header=2, skiprows=[0, 1])
    df.columns = ["Gene", "Cancer_type", "Pathway"]
    return df.dropna(how="all").dropna(subset=["Gene", "Cancer_type", "Pathway"]).reset_index(drop=True)

s2 = load_s2(S2_CSV)

s1 = pd.read_csv(S1_CSV)
s1["Consensus Score"] = pd.to_numeric(s1["Consensus Score"], errors="coerce")           # TODO column for consensus score here
s1 = s1.dropna(subset=["Consensus Score"])
s1_lookup = s1.set_index(["Gene", "Cancer"])["Consensus Score"].to_dict()

# ─────────────────────────────────────────────────────────────────────────────
# 3.  Helper: Load cancer distance CSVs
# ─────────────────────────────────────────────────────────────────────────────
_csv_cache: dict[str, pd.DataFrame | None] = {}

def get_csv(cancer_name: str) -> pd.DataFrame | None:
    if cancer_name not in _csv_cache:
        fpath = DIST_DIR / f"{cancer_name}.csv"
        if not fpath.exists():
            _csv_cache[cancer_name] = None
        else:
            df = pd.read_csv(fpath)
            if DELTA_MEANS_COL not in df.columns or PATHWAY_ID_COL not in df.columns:
                print(f"  [WARN] {cancer_name}.csv missing expected columns.")
                _csv_cache[cancer_name] = None
            else:
                _csv_cache[cancer_name] = df
    return _csv_cache[cancer_name]

def pearson_annotation(x, y):
    mask = np.isfinite(x) & np.isfinite(y)
    r, p = stats.pearsonr(x[mask], y[mask])
    p_str = f"p = {p:.2e}" if p >= 1e-4 else f"p < 1e-4"
    return f"r = {r:.3f},  {p_str}  (n = {mask.sum()})"


# ─────────────────────────────────────────────────────────────────────────────
# 4.  Plot 1: Scatter Consensus Score vs. MEAN Δ-means
# ─────────────────────────────────────────────────────────────────────────────
gene_records = []

for _, row in s2.iterrows():
    gene        = row["Gene"]
    s2_cancer   = row["Cancer_type"]
    s2_pathway  = row["Pathway"]

    # Exact (gene, cancer) match only
    score = s1_lookup.get((gene, s2_cancer))
    if score is None: continue

    my_cancer = s2_cancer_to_mine.get(s2_cancer)
    if my_cancer is None: continue

    hsa_ids = expanded_s2_pathway_to_hsa.get(s2_pathway)
    if not hsa_ids: continue

    df_cancer = get_csv(my_cancer)
    if df_cancer is None: continue

    # Extract relevant pathways and take the MEAN delta_means (reverted per user request)
    relevant = df_cancer[df_cancer[PATHWAY_ID_COL].isin(hsa_ids)]
    if relevant.empty: continue

    mean_delta = relevant[DELTA_MEANS_COL].mean()

    gene_records.append({
        "gene":             gene,
        "s2_cancer":        s2_cancer,
        "my_cancer":        my_cancer,
        "s2_pathway":       s2_pathway,
        "consensus_score":  float(score),
        "delta_means":      mean_delta,
        "cancer_label":     CANCER_FULLNAME.get(my_cancer, my_cancer),
    })

gene_data = pd.DataFrame(gene_records)
print(f"Assembled {len(gene_data)} gene-level data points for Scatter Plot.")

# Wider figure to fit the legend on the side
fig1, ax1 = plt.subplots(figsize=(11, 7))
#fig1.patch.set_facecolor("#F7F5F2")
#ax1.set_facecolor("#F7F5F2")

np.random.seed(42)  # For reproducible jitter

for cat, grp in gene_data.groupby("s2_pathway"):
    col = CATEGORY_COLOURS.get(cat, "#888888")
    
    # Add horizontal jitter so discrete consensus scores don't overlap perfectly
    jitter = np.random.uniform(-0.07, 0.07, size=len(grp))
    x_jittered = grp["consensus_score"] + jitter                        # TODO jitter here
    
    ax1.scatter(
        x_jittered, grp["delta_means"],
        color=col, alpha=0.75, s=45, linewidths=0.4,
        edgecolors="white", label=cat, zorder=3,
    )

x_g = gene_data["consensus_score"].values.astype(float)
y_g = gene_data["delta_means"].values.astype(float)
mask_g = np.isfinite(x_g) & np.isfinite(y_g)

if mask_g.sum() > 2 and x_g[mask_g].std() > 0:
    sl_g, ic_g, *_ = stats.linregress(x_g[mask_g], y_g[mask_g])
    x_ln_g = np.linspace(x_g[mask_g].min(), x_g[mask_g].max(), 200)
    ax1.plot(x_ln_g, sl_g * x_ln_g + ic_g,
             color="#333333", linewidth=1.4, linestyle="--", alpha=0.7, zorder=5)

annot_g = pearson_annotation(x_g, y_g)
ax1.text(0.97, 0.03, annot_g, transform=ax1.transAxes,
         fontsize=10, ha="right", va="bottom",
         bbox=dict(boxstyle="round,pad=0.35", fc="white", alpha=0.85, ec="#CCCCCC"))

ax1.set_xlabel("Consensus Score  (Bailey et al. 2018, Table S1)", fontsize=12)
ax1.set_ylabel("Mean Δ-means  (cluster-expanded pathway distance score)", fontsize=12)
ax1.set_title("S1 Consensus Score  vs.  Pathway Δ-means\n"
              "(one point per gene × cancer × pathway category)", fontsize=13, pad=14)
ax1.spines[["top", "right"]].set_visible(False)
ax1.grid(axis="y", alpha=0.25, linestyle="--")

handles_g = [
    mpatches.Patch(facecolor=CATEGORY_COLOURS.get(c, "#888"), label=c, linewidth=0)
    for c in sorted(CATEGORY_COLOURS) if c in gene_data["s2_pathway"].values
]
# Move legend entirely outside the plot to prevent covering data
ax1.legend(handles=handles_g, fontsize=8, ncol=1, 
           bbox_to_anchor=(1.02, 1.0), loc="upper left",
           framealpha=0.85, edgecolor="#CCCCCC",
           title="Pathway category", title_fontsize=9)

plt.tight_layout()
out_path_1 = OUT_DIR / "scatter_consensus_score_vs_delta_means.png"
fig1.savefig(out_path_1, dpi=600, bbox_inches="tight")
print(f"Saved → {out_path_1}")
plt.close(fig1)


# ─────────────────────────────────────────────────────────────────────────────
# 5.  Plot 2: Violin (Driver-Mapped Clusters vs Other Clusters, MAX Δ-means)
# ─────────────────────────────────────────────────────────────────────────────
def violin_with_strip(ax, groups: dict[str, np.ndarray], colours: dict[str, str],
                      title: str, ylabel: str) -> None:
    labels = list(groups.keys())
    arrays = [groups[l] for l in labels]

    parts = ax.violinplot(arrays, positions=range(len(labels)),
                          showmedians=True, showextrema=False, widths=0.6)
    for i, (pc, lbl) in enumerate(zip(parts["bodies"], labels)):
        pc.set_facecolor(colours[lbl])
        pc.set_alpha(0.65)
        pc.set_edgecolor("#444444")
        pc.set_linewidth(1.0)
    parts["cmedians"].set_color("#222222")
    parts["cmedians"].set_linewidth(1.8)
    
    for i, arr in enumerate(arrays):
        median_val = np.median(arr[np.isfinite(arr)])
        ax.text(i, median_val, f"{median_val:.2f}", 
                ha="center", va="bottom", fontsize=8, fontweight="bold",
                color="#222222")

    rng = np.random.default_rng(42)
    for i, (arr, lbl) in enumerate(zip(arrays, labels)):
        jitter = rng.uniform(-0.15, 0.15, size=len(arr))
        ax.scatter(i + jitter, arr, color=colours[lbl],
                   alpha=0.40, s=9, linewidths=0, zorder=0)

    if len(arrays) == 2 and len(arrays[0]) > 0 and len(arrays[1]) > 0:
        u_stat, p_mw = stats.mannwhitneyu(arrays[0], arrays[1], alternative="two-sided")
        p_str = f"p = {p_mw:.2e}" if p_mw >= 1e-4 else "p < 1e-4"
        n0, n1 = len(arrays[0]), len(arrays[1])
        ax.text(0.97, 0.97, f"Mann-Whitney  {p_str}\nClusters (N) = {n0} vs {n1}",
                transform=ax.transAxes, fontsize=9, ha="right", va="top",
                bbox=dict(boxstyle="round,pad=0.4", fc="white", alpha=0.85, ec="#CCCCCC"))

    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, fontsize=10)
    ax.set_ylabel(ylabel, fontsize=11)
    ax.set_title(title, fontsize=12, pad=10)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", alpha=0.3, linestyle="--")
    ax.axhline(0, color="#888888", linewidth=1.0, linestyle=":", alpha=0.8)


all_my_cancers = sorted({v for v in s2_cancer_to_mine.values() if v is not None})
all_rows = []

for cancer in all_my_cancers:
    df_full = get_csv(cancer)
    if df_full is None: continue
    df_full = df_full.copy()
    
    # Map each KEGG ID to its Cluster ID (fallback to KEGG ID if somehow unclustered)
    df_full['cluster_id'] = df_full[PATHWAY_ID_COL].apply(lambda x: p2c.get(x, x))
    
    # A cluster is "driver-mapped" if ANY of its pathways belong to the expanded S2 set
    df_full['is_driver'] = df_full[PATHWAY_ID_COL].isin(global_expanded_s2_hsa_ids)
    
    # Group by Cluster ID and take the MAX delta_means for the entire cluster
    cluster_agg = df_full.groupby('cluster_id').agg(
        max_delta_means=(DELTA_MEANS_COL, 'max'),
        has_driver=('is_driver', 'any')
    ).reset_index()
    
    cluster_agg["my_cancer"] = cancer
    all_rows.append(cluster_agg)

all_cancer_df = pd.concat(all_rows, ignore_index=True) if all_rows else pd.DataFrame()

if not all_cancer_df.empty:
    fig2, ax2 = plt.subplots(figsize=(7, 6))
    #fig2.patch.set_facecolor("#F7F5F2")
    #ax2.set_facecolor("#F7F5F2")

    grp_driver = all_cancer_df.loc[all_cancer_df["has_driver"],  "max_delta_means"].dropna().values
    grp_other  = all_cancer_df.loc[~all_cancer_df["has_driver"], "max_delta_means"].dropna().values

    # High-contrast colors
    VIOLIN_COLORS = {
        "Driver-mapped\npathways": "#D65F5F", # Strong Red
        "Other\npathways":         "#6C8EBF"  # Steel Blue
    }

    violin_with_strip(
        ax2,
        groups  = {"Driver-mapped\npathways": grp_driver, "Other\npathways": grp_other},
        colours = VIOLIN_COLORS,
        title   = "Max Δ-means: Driver-mapped pathway clusters vs. others\n(global flag, pooled across cancers)",
        ylabel  = "Max Δ-means per Cluster",
    )

    plt.tight_layout()
    out_path_2 = OUT_DIR / "violin_driver_vs_other_global.png"
    fig2.savefig(out_path_2, dpi=600, bbox_inches="tight")
    print(f"Saved → {out_path_2}")
    plt.close(fig2)

# Save data backend
gene_data.to_csv(OUT_DIR / "scatter_data_consensus_vs_delta.csv", index=False)
print("\nDone.")