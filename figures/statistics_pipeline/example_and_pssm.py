"""
Recreates panels A and B of Figure 2.

Output:
- A1: p53 pathway pathogenic_prob distributions per substitution type (4x4 grid).
- A2: PSSM Matrix simple heatmap (using MICHAL_HN1_PSSM).
- B1: PSSM-weighted background vs. observed cancer distribution overlay.

Data source: KEGG pathway hsa04115 (p53 signaling), real pathogenic_prob scores only.
No synthetic data is generated anywhere in this script.
"""

import sys
import os
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import seaborn as sns

# -- Setup: Nested code location --
PROJECT_ROOT = Path(__file__).parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from definitions import (
    MICHAL_HN1_PSSM,
    COLOR_MAP,
    NUMBER_OF_BINS,
    KEGG_PATHWAY_SCORES_P,
)
from distance_utils import get_bg_histogram_after_pssm, get_cancer_histogram

# -- 0. Configuration --
OUTPUT_DIR  = Path(__file__).parents[0]
P53_PATHWAY = "hsa04115"   # KEGG: TP53 signaling pathway
CANCER_NAME = "brca"       # ← set your cancer cohort name here

BASES = ['A', 'C', 'G', 'T']
SUBSTITUTION_PAIRS = [f"{i}>{j}" for i in BASES for j in BASES if i != j]

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

# Color palette
COL_HIST   = COLOR_MAP.get("dark-blue",  "#2B5C8A")
COL_BG     = COLOR_MAP.get("light-blue", "#6BAED6")
COL_CANCER = COLOR_MAP.get("pathogenic", "#D94F3D")
COL_DIAG   = "#F5F5F5"
COL_SPINE  = "#AAAAAA"


# -- 1. Load p53 pathway CSV (single file, same pattern as p7) ----------------
print(f"Loading p53 pathway scores ({P53_PATHWAY})...")

p53_file = os.path.join(KEGG_PATHWAY_SCORES_P, f"{P53_PATHWAY.replace(':', '_')}.csv")
df_p53 = pd.read_csv(p53_file, usecols=["Ref", "Alt", "pathogenic_prob"])
df_p53 = df_p53.dropna(subset=["Ref", "Alt", "pathogenic_prob"])
df_p53["pair"] = df_p53["Ref"].str.strip().str.upper() + ">" + df_p53["Alt"].str.strip().str.upper()

scores_by_pair = {
    pair: df_p53.loc[df_p53["pair"] == pair, "pathogenic_prob"].values
    for pair in SUBSTITUTION_PAIRS
}

total = sum(len(v) for v in scores_by_pair.values())
print(f"  Loaded {total:,} mutations from {P53_PATHWAY}.")


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
            ax.text(0.5, 0.5, ref_base,
                    ha="center", va="center",
                    transform=ax.transAxes,
                    fontsize=10, color="#999999", fontweight="bold")
            for spine in ax.spines.values():
                spine.set_visible(False)
            ax.set_xticks([])
            ax.set_yticks([])
            continue

        pair = f"{ref_base}>{alt_base}"
        data = scores_by_pair.get(pair, np.array([]))

        if len(data) > 0:
            ax.hist(data, bins=25, range=(0, 1),
                    color=COL_HIST, alpha=0.82,
                    density=True, linewidth=0)

            mu = np.mean(data)
            sigma = np.std(data)
            ax.text(
                0.97, 0.95,
                f"$\\mu$={mu:.2f}\n$\\sigma$={sigma:.2f}",
                transform=ax.transAxes,
                ha="right", va="top",
                fontsize=7, color="#333333",
                linespacing=1.4,
            )

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


# ── Panel A2: PSSM heatmap — normalized to percentages ───────────────────────
print("Generating A2: PSSM matrix heatmap (normalized to %)...")

pssm_raw = np.zeros((4, 4))
for i, ref in enumerate(BASES):
    for j, alt in enumerate(BASES):
        if ref != alt:
            pssm_raw[i, j] = MICHAL_HN1_PSSM.get(f"{ref}>{alt}", 0)

# Normalize: off-diagonal values sum to 100%
off_diag_mask = ~np.eye(4, dtype=bool)
total = pssm_raw[off_diag_mask].sum()
pssm_pct = np.where(off_diag_mask, pssm_raw / total * 100, 0.0)

fig_a2, ax = plt.subplots(figsize=(3.2, 2.8))
fig_a2.suptitle("PSSM Substitution Weights (%)", fontsize=8, fontweight="bold")

annot_array = np.where(
    off_diag_mask,
    np.vectorize(lambda x: f"{x:.1f}%")(pssm_pct),
    "–"
)

mask = np.eye(4, dtype=bool)
sns.heatmap(
    pssm_pct,
    annot=annot_array, fmt="",          # one decimal place, no scientific notation
    cmap="Blues",
    mask=mask,
    linewidths=0.4, linecolor="#DDDDDD",
    cbar_kws={"shrink": 0.75, "label": "Weight (%)"},
    ax=ax, square=True,
    xticklabels=BASES, yticklabels=BASES,
    annot_kws={"size": 6},
    vmin=0,
)
ax.set_xlabel("Alt Base", labelpad=4)
ax.set_ylabel("Ref Base", labelpad=4)
ax.tick_params(left=False, bottom=False)

for k in range(4):
    ax.add_patch(plt.Rectangle(
        (k, k), 1, 1,
        fill=True, color="#EEEEEE", lw=0, zorder=3,
    ))
    ax.text(k + 0.5, k + 0.5, "–",
            ha="center", va="center",
            color="#AAAAAA", fontsize=8, zorder=4)

fig_a2.tight_layout()
fig_a2.savefig(OUTPUT_DIR / "figure_2_a2_p53_pssm_heatmap.png")
plt.close(fig_a2)
print("  Saved figure_2_a2_p53_pssm_heatmap.png")


# ── Panel B1: PSSM-weighted background + cancer overlay ──────────────────────
print("Generating B1: PSSM-weighted background vs. cancer overlay...")

bin_edges = np.linspace(0, 1, NUMBER_OF_BINS + 1)
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

hist_bg, _ = get_bg_histogram_after_pssm(df_p53, pssm_matrix=MICHAL_HN1_PSSM)
hist_cancer = get_cancer_histogram(P53_PATHWAY, f"{CANCER_NAME}.csv")

# Compute μ and σ from bin_centers weighted by histogram mass
def weighted_mean_std(hist, centers):
    w = hist / hist.sum() if hist.sum() > 0 else hist
    mu = np.sum(w * centers)
    sigma = np.sqrt(np.sum(w * (centers - mu) ** 2))
    return mu, sigma

mu_bg,     sigma_bg     = weighted_mean_std(hist_bg,     bin_centers)
mu_cancer, sigma_cancer = weighted_mean_std(hist_cancer, bin_centers)

fig_b1, ax = plt.subplots(figsize=(3.8, 2.8))
fig_b1.suptitle("PSSM-Weighted Background vs. Cancer Distribution",
                fontsize=10, fontweight="bold")

# Background
ax.stairs(hist_bg, bin_edges,
          color=COL_BG, fill=True, alpha=0.5,
          label=(f"Expected (PSSM-Weighted BG)\n"
                 f"$\\mu$={mu_bg:.2f}, $\\sigma$={sigma_bg:.2f}"))
ax.stairs(hist_bg, bin_edges,
          color=COL_BG, linewidth=1.5, alpha=0.8)

# Cancer overlay
if np.sum(hist_cancer) > 0:
    ax.stairs(hist_cancer, bin_edges,
              color=COL_CANCER, fill=True, alpha=0.3,
              label=(f"Observed ({CANCER_NAME.upper()})\n"
                     f"$\\mu$={mu_cancer:.2f}, $\\sigma$={sigma_cancer:.2f}"))
    ax.stairs(hist_cancer, bin_edges,
              color=COL_CANCER, linewidth=1.5)
else:
    print(f"  Warning: no cancer mutations found for {P53_PATHWAY} / {CANCER_NAME}")

ax.set_xlabel("Pathogenicity Probability", labelpad=3)
ax.set_ylabel("Probability Mass",          labelpad=3)
ax.set_xlim(0, 1)

y_max = max(
    hist_bg.max(),
    hist_cancer.max() if np.sum(hist_cancer) > 0 else 0,
) * 1.3
ax.set_ylim(0, y_max)

ax.xaxis.set_major_locator(ticker.MultipleLocator(0.2))
ax.yaxis.set_major_locator(ticker.MaxNLocator(4, prune="both"))
ax.legend(fontsize=6, frameon=False, loc="upper right")

fig_b1.tight_layout()
fig_b1.savefig(OUTPUT_DIR / "figure_2_b1_final_background.png")
plt.close(fig_b1)
print("  Saved figure_2_b1_final_background.png")

print(f"\nAll outputs saved to: {OUTPUT_DIR}")