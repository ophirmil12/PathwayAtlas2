"""
Pathway Co-enrichment Matrix Analysis
--------------------------------------
For each pair of pathways, counts how many cancers they are both significantly
enriched in (q_value < 0.05). Generates three matrices:
  1. Pathogenic-Pathogenic  (delta_means > 0 in both)
  2. Benign-Benign          (delta_means < 0 in both)
  3. Pathogenic-Benign      (delta_means > 0 in one, < 0 in the other)

Pairs that share > 50% of proteins (Jaccard-style: |A ∩ B| / min(|A|,|B|) > 0.5)
are masked (shown as NaN / white in the heatmap).
"""

from plot_boot import *
boot_plot_folder()      # coloring scheme using cycler, and folders traveling

import os
import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
from itertools import combinations
from pathlib import Path
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform

from definitions import BASE_P, RESULTS_DISTANCES_P, KEGG_PATHWAY_METADATA_FILE, PLOTS_P, COLOR_MAP

OUTPUT_DIR        = os.path.join(PLOTS_P, "p9_coenrichment")

Q_VALUE_THRESHOLD  = 0.05
OVERLAP_THRESHOLD  = 0.50   # mask pairs sharing > 50% proteins
MIN_CANCERS        = 5       # only show pathways significant in at least this many cancers
MIN_PAIR_CANCERS   = 2       # row/col shown only if ≥1 pair co-enriched in ≥ this many cancers


print(f"____ Q_VALUE_THRESHOLD: {Q_VALUE_THRESHOLD}"
        f"____ OVERLAP_THRESHOLD: {OVERLAP_THRESHOLD}"
        f"____ MIN_CANCERS: {MIN_CANCERS}"
        f"____ MIN_PAIR_CANCERS: {MIN_PAIR_CANCERS}")


os.makedirs(OUTPUT_DIR, exist_ok=True)


def cluster_order(mat):
    """Reorder rows/cols by hierarchical clustering."""
    filled = np.nan_to_num(mat, nan=0.0)
    dist = np.max(filled) - filled
    np.fill_diagonal(dist, 0)
    try:
        condensed = squareform(dist, checks=False)
        Z = linkage(condensed, method="average")
        return leaves_list(Z)
    except Exception:
        return np.arange(mat.shape[0])

# ── 1. LOAD METADATA ─────────────────────────────────────────────────────────
with open(KEGG_PATHWAY_METADATA_FILE, "rb") as f:
    pathway_metadata = pickle.load(f)

def get_gene_set(pathway_id: str) -> set:
    """Return the set of hsa:XXXX gene IDs for a pathway / module."""
    meta = pathway_metadata.get(pathway_id, {})
    return set(meta.get("genes_ids", []))

def get_pathway_name(pathway_id: str) -> str:
    meta = pathway_metadata.get(pathway_id, {})
    name = meta.get("name", pathway_id)
    # Truncate long names for display
    return name[:55] + "…" if len(name) > 55 else name


# ── 2. LOAD ALL CANCER CSVs ──────────────────────────────────────────────────
csv_files = sorted(Path(RESULTS_DISTANCES_P).glob("*.csv"))
if not csv_files:
    raise FileNotFoundError(f"No CSV files found in {RESULTS_DISTANCES_P}")

print(f"Found {len(csv_files)} cancer CSV files.")

# Build a dict: { pathway_id -> { cancer -> (is_sig, delta_means) } }
pathway_cancer_data: dict[str, dict[str, tuple]] = {}

for csv_path in csv_files:
    cancer = csv_path.stem          # filename without extension = cancer name
    df = pd.read_csv(csv_path)

    required = {"pathway", "q_value", "delta_means"}
    if not required.issubset(df.columns):
        print(f"  Skipping {cancer}: missing columns {required - set(df.columns)}")
        continue

    for _, row in df.iterrows():
        pid  = str(row["pathway"])
        # Skip  modules
        if pid.startswith("M"):
            continue

        sig  = row["q_value"] < Q_VALUE_THRESHOLD
        dm   = row["delta_means"]
        if pd.isna(dm):
            continue
        pathway_cancer_data.setdefault(pid, {})[cancer] = (sig, dm)

print(f"Loaded {len(pathway_cancer_data)} unique pathways across all cancers.")


# ── 3. CLASSIFY EACH PATHWAY PER CANCER ──────────────────────────────────────
# Label: 'pathogenic' (sig & dm>0), 'benign' (sig & dm<0), None (not sig)

def classify(sig: bool, dm: float) -> str | None:
    if not sig:
        return None
    return "pathogenic" if dm > 0 else "benign"

# For each pathway, collect (cancer, label) pairs where it IS significant
pathway_sig_labels: dict[str, dict[str, str]] = {}
for pid, cancer_dict in pathway_cancer_data.items():
    labels = {}
    for cancer, (sig, dm) in cancer_dict.items():
        lbl = classify(sig, dm)
        if lbl:
            labels[cancer] = lbl
    if labels:
        pathway_sig_labels[pid] = labels

# Keep only pathways significant in at least MIN_CANCERS
pathway_sig_labels = {
    pid: labels
    for pid, labels in pathway_sig_labels.items()
    if len(labels) >= MIN_CANCERS
}

pathways = sorted(pathway_sig_labels.keys())
print(f"Pathways significant in ≥{MIN_CANCERS} cancers: {len(pathways)}")


# ── 4. COMPUTE PAIRWISE OVERLAP MASK ─────────────────────────────────────────
print("Computing pairwise protein overlaps …")
gene_sets = {pid: get_gene_set(pid) for pid in pathways}

def overlap_fraction(s1: set, s2: set) -> float:
    if not s1 or not s2:
        return 0.0
    return len(s1 & s2) / min(len(s1), len(s2))

n = len(pathways)
pid_to_idx = {pid: i for i, pid in enumerate(pathways)}

overlap_mask = np.zeros((n, n), dtype=bool)   # True = should be masked
for i, j in combinations(range(n), 2):
    ov = overlap_fraction(gene_sets[pathways[i]], gene_sets[pathways[j]])
    if ov > OVERLAP_THRESHOLD:
        overlap_mask[i, j] = True
        overlap_mask[j, i] = True
np.fill_diagonal(overlap_mask, True)   # mask self-pairs too

print(f"  Masked {overlap_mask.sum() // 2} high-overlap pairs "
      f"(>{OVERLAP_THRESHOLD*100:.0f}% shared proteins).")


# ── 5. BUILD CO-ENRICHMENT COUNT MATRICES ────────────────────────────────────
# For each cancer, collect sets of pathogenic / benign pathway indices
all_cancers = set()
for labels in pathway_sig_labels.values():
    all_cancers.update(labels.keys())
all_cancers = sorted(all_cancers)

pathogenic_per_cancer: dict[str, set] = {c: set() for c in all_cancers}
benign_per_cancer:     dict[str, set] = {c: set() for c in all_cancers}

for pid, labels in pathway_sig_labels.items():
    idx = pid_to_idx[pid]
    for cancer, lbl in labels.items():
        if lbl == "pathogenic":
            pathogenic_per_cancer[cancer].add(idx)
        else:
            benign_per_cancer[cancer].add(idx)

pp_matrix = np.zeros((n, n), dtype=int)   # pathogenic–pathogenic
bb_matrix = np.zeros((n, n), dtype=int)   # benign–benign
pb_matrix = np.zeros((n, n), dtype=int)   # pathogenic–benign

for cancer in all_cancers:
    path_set   = pathogenic_per_cancer[cancer]
    benign_set = benign_per_cancer[cancer]

    # pathogenic–pathogenic
    for i, j in combinations(path_set, 2):
        pp_matrix[i, j] += 1
        pp_matrix[j, i] += 1

    # benign–benign
    for i, j in combinations(benign_set, 2):
        bb_matrix[i, j] += 1
        bb_matrix[j, i] += 1

    # pathogenic–benign (asymmetric: row=pathogenic, col=benign)
    for i in path_set:
        for j in benign_set:
            pb_matrix[i, j] += 1

print(f"Co-enrichment matrices built for {len(all_cancers)} cancers.")


# ── 6. APPLY MASK & CLIP ROWS/COLS WITH NO SIGNAL ────────────────────────────
def apply_mask(mat: np.ndarray) -> np.ndarray:
    masked = mat.astype(float)
    masked[overlap_mask] = np.nan
    return masked

pp_masked = apply_mask(pp_matrix)
bb_masked = apply_mask(bb_matrix)
pb_masked = apply_mask(pb_matrix)

def active_indices(mat: np.ndarray, min_val: int = MIN_PAIR_CANCERS) -> np.ndarray:
    """Rows/cols that have at least one non-NaN entry >= min_val."""
    with np.errstate(invalid="ignore"):
        row_max = np.nanmax(mat, axis=1)
    return np.where(row_max >= min_val)[0]

# Union of active indices across all three matrices
active = np.union1d(
    np.union1d(active_indices(pp_masked), active_indices(bb_masked)),
    active_indices(pb_masked)
)
active = np.array(sorted(active))

print(f"Pathways with ≥1 co-enriched partner: {len(active)}")

pp_sub = pp_masked[np.ix_(active, active)]
bb_sub = bb_masked[np.ix_(active, active)]
pb_sub = pb_masked[np.ix_(active, active)]

# Cluster each matrix independently
pp_ord = cluster_order(pp_sub)
bb_ord = cluster_order(bb_sub)
pp_sub = pp_sub[np.ix_(pp_ord, pp_ord)]
bb_sub = bb_sub[np.ix_(bb_ord, bb_ord)]

pp_labels = [get_pathway_name(pathways[active[i]]) for i in pp_ord]
bb_labels = [get_pathway_name(pathways[active[i]]) for i in bb_ord]
labels_sub = [get_pathway_name(pathways[i]) for i in active]  # kept for pb


# ── 7. PLOTTING ───────────────────────────────────────────────────────────────
def make_cmap(base_hex: str) -> mcolors.LinearSegmentedColormap:
    """White → base_hex colormap."""
    return mcolors.LinearSegmentedColormap.from_list(
        "custom", ["#FFFFFF", base_hex]
    )

def plot_matrix(mat: np.ndarray,
                row_labels: list[str],
                col_labels: list[str],
                title: str,
                base_color: str,
                filename: str):
    """Plot a co-enrichment heatmap and save to OUTPUT_DIR."""
    m, n2 = mat.shape
    # Dynamic figure size: cap at 60 inches
    cell = max(0.25, min(0.55, 40 / max(m, n2)))
    fig_w = min(60, n2 * cell + 4)
    fig_h = min(60, m  * cell + 4)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    cmap = make_cmap(base_color)
    cmap.set_bad(color="#EEEEEE")   # masked cells → light grey

    # Per-matrix vmax at 95th percentile so colours don't wash out
    vals = mat[~np.isnan(mat) & (mat > 0)]
    vmax = max(int(np.percentile(vals, 95)), 1) if len(vals) else 1

    norm = mcolors.PowerNorm(gamma=0.5, vmin=0, vmax=vmax)
    im = ax.imshow(mat, cmap=cmap, norm=norm, aspect="auto",
                interpolation="nearest")

    # Always show labels, font size adapts to count
    fs = max(4, min(8, int(200 / max(m, n2))))
    ax.set_yticks(range(m))
    ax.set_yticklabels(row_labels, fontsize=fs)
    ax.set_xticks(range(n2))
    ax.set_xticklabels(col_labels, rotation=90, fontsize=fs)

    cbar = fig.colorbar(im, ax=ax, shrink=0.5, pad=0.02)
    cbar.set_label("# cancers co-enriched", fontsize=10)

    ax.set_title(title, fontsize=13, fontweight="bold", pad=12)
    ax.set_xlabel("Pathway", fontsize=10)
    ax.set_ylabel("Pathway", fontsize=10)

    plt.tight_layout()
    out_path = os.path.join(OUTPUT_DIR, filename)
    fig.savefig(out_path, dpi=180, bbox_inches="tight")             # TODO change dpi to 300 for final version
    plt.close(fig)
    print(f"  Saved: {out_path}")
    return out_path


print("\nGenerating heatmaps …")

plot_matrix(
    pp_sub, pp_labels, pp_labels,
    title=f"Pathogenic–Pathogenic Co-enrichment\n"
          f"(q<{Q_VALUE_THRESHOLD}, Delta-Means>0 in both; "
          f"grey = >{OVERLAP_THRESHOLD*100:.0f}% shared proteins)",
    base_color=COLOR_MAP["pathogenic"],
    filename="coenrichment_pp.png",
)

plot_matrix(
    bb_sub, bb_labels, bb_labels,
    title=f"Benign–Benign Co-enrichment\n"
          f"(q<{Q_VALUE_THRESHOLD}, Delta-Means<0 in both; "
          f"grey = >{OVERLAP_THRESHOLD*100:.0f}% shared proteins)",
    base_color=COLOR_MAP["benign"],
    filename="coenrichment_bb.png",
)

# For pathogenic–benign the axes differ: rows = pathogenic, cols = benign
# Build index subsets
path_active  = np.array(sorted({i for i in active
                                 if any(pathway_sig_labels[pathways[i]].get(c) == "pathogenic"
                                        for c in all_cancers)}))
benign_active = np.array(sorted({i for i in active
                                   if any(pathway_sig_labels[pathways[i]].get(c) == "benign"
                                          for c in all_cancers)}))

if len(path_active) and len(benign_active):
    pb_rect = pb_masked[np.ix_(path_active, benign_active)]
    row_lbl = [get_pathway_name(pathways[i]) for i in path_active]
    col_lbl = [get_pathway_name(pathways[i]) for i in benign_active]
    plot_matrix(
        pb_rect, row_lbl, col_lbl,
        title=f"Pathogenic–Benign Co-enrichment\n"
              f"(rows: Delta-Means>0, cols: Delta-Means<0; q<{Q_VALUE_THRESHOLD}; "
              f"grey = >{OVERLAP_THRESHOLD*100:.0f}% shared proteins)",
        base_color=COLOR_MAP["light-blue"],
        filename="coenrichment_pb.png",
    )
else:
    print("  Skipping pathogenic–benign: insufficient data.")


# ── 8. EXPORT SUMMARY TABLES ──────────────────────────────────────────────────
def matrix_to_df(mat: np.ndarray, row_ids, col_ids) -> pd.DataFrame:
    df = pd.DataFrame(mat,
                      index=[pathways[i] for i in row_ids],
                      columns=[pathways[i] for i in col_ids])
    df.index.name   = "pathway_row"
    df.columns.name = "pathway_col"
    return df

pp_df = matrix_to_df(pp_masked[np.ix_(active, active)], active, active)
bb_df = matrix_to_df(bb_masked[np.ix_(active, active)], active, active)

pp_df.to_csv(os.path.join(OUTPUT_DIR, "coenrichment_pp.csv"))
bb_df.to_csv(os.path.join(OUTPUT_DIR, "coenrichment_bb.csv"))
if len(path_active) and len(benign_active):
    pb_df = matrix_to_df(pb_masked[np.ix_(path_active, benign_active)],
                         path_active, benign_active)
    pb_df.to_csv(os.path.join(OUTPUT_DIR, "coenrichment_pb.csv"))

print("\nDone. All outputs written to:", OUTPUT_DIR)