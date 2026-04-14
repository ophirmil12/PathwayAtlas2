"""
Recreates panels A and B of Figure 2.

Output:
- A1: p53 pathway pathogenic_prob distributions per substitution type (4x4 grid).
- A2: PSSM Matrix simple heatmap (using MICHAL_HN1_PSSM).
- B1: PSSM-weighted background vs. observed cancer distribution overlay.
     Generated for 3 pathways:
       - hsa04115 (p53 signaling)         → pathogenic color
       - hsa04960 (Aldosterone/sodium)    → benign color
       - hsa05231 (Choline metabolism)    → grey color
"""

import sys
import os
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import seaborn as sns
import pickle

# -- Setup: Nested code location --
PROJECT_ROOT = Path(__file__).parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from definitions import (
    MICHAL_HN1_PSSM,
    COLOR_MAP,
    KEGG_PATHWAY_SCORES_P,
    KEGG_PATHWAY_METADATA_FILE,
    CBIO_CANCER_MUTATIONS_P,
)

# -- 0. Configuration --
OUTPUT_DIR  = Path(__file__).parents[0]
P53_PATHWAY = "hsa04115"
CANCER_NAME = "brca"

NUMBER_OF_BINS = 25   # local, independent of definitions.py

BASES = ['A', 'C', 'G', 'T']
SUBSTITUTION_PAIRS = [f"{i}>{j}" for i in BASES for j in BASES if i != j]

# Pathway definitions: (pathway_id, short_label, nature, bg_color, cancer_color, output_suffix)
B1_PATHWAYS = [
    (
        "hsa04115",
        "p53 Signaling",
        "pathogenic",
        COLOR_MAP.get("light-blue",  "#B7CADB"),
        COLOR_MAP.get("pathogenic",  "#CB7673"),
        "p53_signaling",
    ),
    (
        "hsa04960",
        "Aldosterone-Regulated Sodium Reabsorption",
        "benign",
        COLOR_MAP.get("light-blue",  "#B7CADB"),
        COLOR_MAP.get("benign",      "#447D68"),
        "aldosterone_sodium",
    ),
    (
        "hsa05231",
        "Choline Metabolism in Cancer",
        "natural",
        COLOR_MAP.get("light-blue",  "#B7CADB"),
        COLOR_MAP.get("grey",        "#787A91"),
        "choline_metabolism",
    ),
]

# Publication-quality style
STYLE = {
    "font.family": "DejaVu Sans",
    "font.size": 7,
    "axes.titlesize": 7,
    "axes.labelsize": 6.5,
    "axes.linewidth": 0.6,
    "xtick.labelsize": 5.5,
    "ytick.labelsize": 5.5,
    "xtick.major.width": 0.5,
    "ytick.major.width": 0.5,
    "xtick.major.size": 2.5,
    "ytick.major.size": 2.5,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.05,
}
plt.rcParams.update(STYLE)

COL_DIAG  = "#F5F5F5"
COL_SPINE = "#AAAAAA"


# -- Local histogram utilities (independent of distance_utils.py) -------------

def make_bin_edges(n_bins=NUMBER_OF_BINS):
    return np.linspace(0, 1, n_bins + 1)


def get_bg_histogram_after_pssm(scores_df: pd.DataFrame,
                                 pssm_matrix=MICHAL_HN1_PSSM,
                                 bins=NUMBER_OF_BINS):
    """
    PSSM-weighted background histogram. Replicates distance_utils logic,
    but uses the local NUMBER_OF_BINS.
    """
    scores_df = scores_df.copy()
    scores_df['mut_type'] = scores_df['Ref'].str.upper() + '>' + scores_df['Alt'].str.upper()

    bin_edges = make_bin_edges(bins)
    final_histogram = np.zeros(bins)

    grouped = scores_df.groupby('mut_type')['pathogenic_prob']

    for mut_type, pssm_weight in pssm_matrix.items():
        if mut_type in grouped.groups:
            scores = grouped.get_group(mut_type).dropna().values
            if len(scores) == 0:
                continue
            counts, _ = np.histogram(scores, bins=bin_edges)
            final_histogram += counts.astype(float) * pssm_weight

    if final_histogram.sum() > 0:
        final_histogram /= final_histogram.sum()

    return final_histogram, bin_edges


def get_cancer_histogram(pathway_id: str, cancer_file: str, bins=NUMBER_OF_BINS):
    """
    Normalized histogram of observed cancer mutations for a pathway.
    Replicates distance_utils logic, but uses the local NUMBER_OF_BINS.
    """
    bin_edges = make_bin_edges(bins)

    if not os.path.exists(KEGG_PATHWAY_METADATA_FILE):
        return np.zeros(bins), bin_edges

    with open(KEGG_PATHWAY_METADATA_FILE, "rb") as f:
        pathway_metadata = pickle.load(f)

    entry = pathway_metadata.get(pathway_id)
    if not isinstance(entry, dict):
        return np.zeros(bins), bin_edges

    pathway_genes = set(entry.get('genes_ids', []))
    if not pathway_genes:
        return np.zeros(bins), bin_edges

    cancer_path = os.path.join(CBIO_CANCER_MUTATIONS_P, cancer_file)
    if not os.path.exists(cancer_path):
        return np.zeros(bins), bin_edges

    df = pd.read_csv(cancer_path, usecols=["KeggId", "pathogenic_prob"])
    mask = df['KeggId'].apply(lambda x: str(x).strip() in pathway_genes if pd.notna(x) else False)
    scores = df.loc[mask, 'pathogenic_prob'].dropna().values

    if len(scores) == 0:
        return np.zeros(bins), bin_edges

    counts, _ = np.histogram(scores, bins=bin_edges)
    total = counts.sum()
    hist = counts.astype(float) / total if total > 0 else np.zeros(bins)

    return hist, bin_edges


def weighted_mean_std(hist, centers):
    total = hist.sum()
    if total == 0:
        return 0.0, 0.0
    w = hist / total
    mu = np.sum(w * centers)
    sigma = np.sqrt(np.sum(w * (centers - mu) ** 2))
    return mu, sigma


# -- 1. Load p53 pathway CSV --------------------------------------------------
print(f"Loading p53 pathway scores ({P53_PATHWAY})...")

p53_file = os.path.join(KEGG_PATHWAY_SCORES_P, f"{P53_PATHWAY}.csv")
df_p53 = pd.read_csv(p53_file, usecols=["Ref", "Alt", "pathogenic_prob"])
df_p53 = df_p53.dropna(subset=["Ref", "Alt", "pathogenic_prob"])
df_p53["pair"] = df_p53["Ref"].str.strip().str.upper() + ">" + df_p53["Alt"].str.strip().str.upper()

scores_by_pair = {
    pair: df_p53.loc[df_p53["pair"] == pair, "pathogenic_prob"].values
    for pair in SUBSTITUTION_PAIRS
}
print(f"  Loaded {sum(len(v) for v in scores_by_pair.values()):,} mutations.")


# ── Panel A1: 4×4 pathogenic_prob histogram grid ─────────────────────────────
print("Generating A1: pathogenic_prob distributions by substitution type...")

fig_a1, axes = plt.subplots(
    len(BASES), len(BASES),
    figsize=(5.5, 4.8),
    sharex=True, sharey=True,
)
fig_a1.suptitle("Pathogenic Probability Distributions by Substitution Type\n(p53 Signaling Pathway)",
                fontsize=8, fontweight="bold", y=1.01)

for i, ref_base in enumerate(BASES):
    for j, alt_base in enumerate(BASES):
        ax = axes[i, j]
        if ref_base == alt_base:
            ax.set_facecolor(COL_DIAG)
            ax.text(0.5, 0.5, ref_base, ha="center", va="center",
                    transform=ax.transAxes, fontsize=10, color="#999999", fontweight="bold")
            for spine in ax.spines.values():
                spine.set_visible(False)
            ax.set_xticks([])
            ax.set_yticks([])
            continue

        pair = f"{ref_base}>{alt_base}"
        data = scores_by_pair.get(pair, np.array([]))

        if len(data) > 0:
            ax.hist(data, bins=NUMBER_OF_BINS, range=(0, 1),
                    color=COLOR_MAP.get("dark-blue", "#5b7d87"), alpha=0.82,
                    density=True, linewidth=0)
            mu, sigma = np.mean(data), np.std(data)
            ax.text(0.97, 0.95, f"$\\mu$={mu:.2f}\n$\\sigma$={sigma:.2f}",
                    transform=ax.transAxes, ha="right", va="top",
                    fontsize=7, color="#333333", linespacing=1.4)

        ax.set_title(pair, pad=2, fontsize=6.5, color="#333333")
        ax.set_xlim(0, 1)
        for spine in ax.spines.values():
            spine.set_color(COL_SPINE)
            spine.set_linewidth(0.5)
        ax.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
        ax.yaxis.set_major_locator(ticker.MaxNLocator(3, prune="both"))

fig_a1.text(0.5, -0.01, "Pathogenic Prob.", ha="center", va="top", fontsize=12)
fig_a1.text(-0.01, 0.5, "Density", ha="right", va="center", fontsize=12, rotation="vertical")
fig_a1.tight_layout(h_pad=0.5, w_pad=0.4)
fig_a1.savefig(OUTPUT_DIR / "figure_2_a1_p53_llrs.png")
plt.close(fig_a1)
print("  Saved figure_2_a1_p53_llrs.png")


# ── Panel A2: PSSM heatmap ────────────────────────────────────────────────────
print("Generating A2: PSSM matrix heatmap (normalized to %)...")

pssm_raw = np.zeros((4, 4))
for i, ref in enumerate(BASES):
    for j, alt in enumerate(BASES):
        if ref != alt:
            pssm_raw[i, j] = MICHAL_HN1_PSSM.get(f"{ref}>{alt}", 0)

off_diag_mask = ~np.eye(4, dtype=bool)
total_pssm = pssm_raw[off_diag_mask].sum()
pssm_pct = np.where(off_diag_mask, pssm_raw / total_pssm * 100, 0.0)

fig_a2, ax = plt.subplots(figsize=(3.2, 2.8))
fig_a2.suptitle("PSSM Substitution Weights (%)", fontsize=8, fontweight="bold")

annot_array = np.where(
    off_diag_mask,
    np.vectorize(lambda x: f"{x:.1f}%")(pssm_pct),
    "–"
)
mask_diag = np.eye(4, dtype=bool)
sns.heatmap(
    pssm_pct, annot=annot_array, fmt="", cmap="Blues",
    mask=mask_diag, linewidths=0.4, linecolor="#DDDDDD",
    cbar_kws={"shrink": 0.75, "label": "Weight (%)"},
    ax=ax, square=True, xticklabels=BASES, yticklabels=BASES,
    annot_kws={"size": 6}, vmin=0,
)
ax.set_xlabel("Alt Base", labelpad=4)
ax.set_ylabel("Ref Base", labelpad=4)
ax.tick_params(left=False, bottom=False)
for k in range(4):
    ax.add_patch(plt.Rectangle((k, k), 1, 1, fill=True, color="#EEEEEE", lw=0, zorder=3))
    ax.text(k + 0.5, k + 0.5, "–", ha="center", va="center",
            color="#AAAAAA", fontsize=8, zorder=4)
fig_a2.tight_layout()
fig_a2.savefig(OUTPUT_DIR / "figure_2_a2_p53_pssm_heatmap.png")
plt.close(fig_a2)
print("  Saved figure_2_a2_p53_pssm_heatmap.png")


# ── Panel B1: 3 separate figures, one per pathway ────────────────────────────
print("Generating B1 plots...")

# Load each pathway's scores file upfront
def load_pathway_scores(pathway_id):
    fpath = os.path.join(KEGG_PATHWAY_SCORES_P, f"{pathway_id}.csv")
    df = pd.read_csv(fpath, usecols=["Ref", "Alt", "pathogenic_prob"])
    return df.dropna(subset=["Ref", "Alt", "pathogenic_prob"])

for pathway_id, pathway_label, nature, col_bg, col_cancer, suffix in B1_PATHWAYS:
    print(f"  Processing {pathway_id} ({nature})...")

    df_scores = load_pathway_scores(pathway_id)
    hist_bg, bin_edges = get_bg_histogram_after_pssm(df_scores, bins=NUMBER_OF_BINS)
    hist_cancer, _     = get_cancer_histogram(pathway_id, f"{CANCER_NAME}.csv", bins=NUMBER_OF_BINS)

    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    mu_bg,     sigma_bg     = weighted_mean_std(hist_bg,     bin_centers)
    mu_cancer, sigma_cancer = weighted_mean_std(hist_cancer, bin_centers)

    fig, ax = plt.subplots(figsize=(3.8, 2.8))
    fig.suptitle(f"PSSM-Weighted Background vs. Cancer Distribution\n{pathway_label}",
                 fontsize=8, fontweight="bold")

    ax.stairs(hist_bg, bin_edges, color=col_bg, fill=True, alpha=0.5,
              label="Expected (PSSM-Weighted BG)")
    ax.stairs(hist_bg, bin_edges, color=col_bg, linewidth=1.5, alpha=0.8)

    if hist_cancer.sum() > 0:
        ax.stairs(hist_cancer, bin_edges, color=col_cancer, fill=True, alpha=0.3,
                  label=f"Observed ({CANCER_NAME.upper()})")
        ax.stairs(hist_cancer, bin_edges, color=col_cancer, linewidth=1.5)
    else:
        print(f"    Warning: no cancer mutations found for {pathway_id} / {CANCER_NAME}")

    ax.set_xlabel("Pathogenicity Probability", labelpad=3)
    ax.set_ylabel("Probability Mass", labelpad=3)
    ax.set_xlim(0, 1)
    y_max = max(hist_bg.max(), hist_cancer.max() if hist_cancer.sum() > 0 else 0) * 1.3
    ax.set_ylim(0, y_max)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(0.2))
    ax.yaxis.set_major_locator(ticker.MaxNLocator(4, prune="both"))
    ax.legend(fontsize=6, frameon=False, loc="upper right")

    # Dashed vertical mean lines
    ax.axvline(mu_bg,     color=col_bg,     linestyle='--', linewidth=1.2, alpha=0.9)
    ax.axvline(mu_cancer, color=col_cancer, linestyle='--', linewidth=1.2, alpha=0.9)

    # Double-headed arrow + Delta mu label above it
    arrow_y = y_max * 0.80
    delta_mu = mu_cancer - mu_bg

    if nature != "natural":
        ax.annotate(
            "", xy=(mu_cancer, arrow_y), xytext=(mu_bg, arrow_y),
            arrowprops=dict(arrowstyle='<->', color='#444444', lw=1.2),
        )

    ax.text(
        (mu_bg + mu_cancer) / 2, arrow_y * 1.04,
        f"$\\Delta\\mu$ = {delta_mu:+.2f}",
        ha='center', va='bottom', fontsize=9, color='#444444',
    )

    fig.tight_layout()
    out_name = f"figure_2_b1_{suffix}.png"
    fig.savefig(OUTPUT_DIR / out_name)
    plt.close(fig)
    print(f"    Saved {out_name}")

print(f"\nAll outputs saved to: {OUTPUT_DIR}")