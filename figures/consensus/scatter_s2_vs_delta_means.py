"""
scatter_pancan_delta_vs_qvalue.py
──────────────────────────────────
Runs for two cancer types: PANCAN and BRCA.

For each:
  Plot 1 – Scatter: delta_means (x) vs −log10(q_value) (y)
  Plot 2 – Violin: S2-mapped vs non-mapped pathways, delta_means distribution
"""

from pathlib import Path
import sys
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
from scipy import stats

# ── paths ─────────────────────────────────────────────────────────────────────
PROJECT_ROOT = Path(__file__).parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from definitions import RESULTS_DISTANCES_P
from mappings import s2_pathway_to_hsa, s2_cancer_to_mine

HERE   = Path(__file__).parent
S2_CSV = HERE / "Table_S2.csv"
S1_CSV = HERE / "Table_S1.csv"

DELTA_COL   = "delta_means"
PATHWAY_COL = "pathway"
QVAL_COL    = "q_value"

# ── cancer-type configs ───────────────────────────────────────────────────────
RUNS = [
    {
        "cancer_label": s2_code,
        "s2_filter":    s2_code,
        "s1_filter":    s2_code,
        "pan_csv":      Path(RESULTS_DISTANCES_P) / f"{file_name}.csv",
        "out_suffix":   s2_code.lower(),
    }
    for s2_code, file_name in s2_cancer_to_mine.items()
    if file_name is not None  # skips PCPG
]

# ── load S2 and S1 once ───────────────────────────────────────────────────────
s2 = pd.read_csv(S2_CSV)
s2.columns = [c.strip() for c in s2.columns]
if "Cancer type" in s2.columns:
    s2 = s2.rename(columns={"Cancer type": "Cancer_type"})

s1 = pd.read_csv(S1_CSV)
s1.columns = [c.strip() for c in s1.columns]
s1["Consensus Score"] = pd.to_numeric(s1["Consensus Score"], errors="coerce")

# ═══════════════════════════════════════════════════════════════════════════════
# MAIN LOOP
# ═══════════════════════════════════════════════════════════════════════════════
for run in RUNS:
    cancer_label = run["cancer_label"]
    print(f"\n{'='*60}")
    print(f"  Running: {cancer_label}")
    print(f"{'='*60}")

    OUT_DIR = HERE / "output_plots" / run["out_suffix"]
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # ── 1. Load cancer-specific CSV ───────────────────────────────────────────
    pan = pd.read_csv(run["pan_csv"])
    pan[DELTA_COL] = pd.to_numeric(pan[DELTA_COL], errors="coerce")
    pan[QVAL_COL]  = pd.to_numeric(pan[QVAL_COL],  errors="coerce")
    pan = pan.dropna(subset=[DELTA_COL, QVAL_COL]).copy()
    pan["neg_log10_q"] = -np.log10(pan[QVAL_COL].clip(lower=1e-300))
    print(f"{run['pan_csv'].name}: {len(pan)} pathways after dropping NaN")

    if len(pan) == 0:
        print(f"Skipping {cancer_label} — no valid rows in {run['pan_csv'].name}")
        continue

    # ── 2. Filter S2 for this cancer type ─────────────────────────────────────
    s2_ct = s2[s2["Cancer_type"].str.upper() == run["s2_filter"].upper()].copy()
    print(f"S2 {cancer_label} rows: {len(s2_ct)}")

    pancan_hsa_ids: set[str] = set()
    cat_to_hsa: dict[str, list[str]] = {}
    for cat in s2_ct["Pathway"].dropna().unique():
        hsa_list = s2_pathway_to_hsa.get(cat)
        if hsa_list:
            cat_to_hsa[cat] = hsa_list
            pancan_hsa_ids.update(hsa_list)
    print(f"KEGG IDs mapped from S2 {cancer_label} categories: {len(pancan_hsa_ids)}")

    # ── 3. Filter S1 for this cancer type ─────────────────────────────────────
    s1_ct = s1[s1["Cancer"].str.upper() == run["s1_filter"].upper()].dropna(subset=["Consensus Score"])
    gene_to_score: dict[str, float] = dict(zip(s1_ct["Gene"], s1_ct["Consensus Score"]))
    print(f"S1 {cancer_label} genes with consensus score: {len(gene_to_score)}")

    # ── 4. Assign consensus score to each S2-mapped pathway ───────────────────
    kegg_to_cats: dict[str, list[str]] = {}
    for cat, hsa_list in cat_to_hsa.items():
        for hid in hsa_list:
            kegg_to_cats.setdefault(hid, []).append(cat)

    cat_to_genes: dict[str, list[str]] = (
        s2_ct.groupby("Pathway")["Gene"].apply(list).to_dict()
    )

    def consensus_score_for_pathway(hid: str) -> float:
        cats   = kegg_to_cats.get(hid, [])
        genes  = [g for cat in cats for g in cat_to_genes.get(cat, [])]
        scores = [gene_to_score[g] for g in genes if g in gene_to_score]
        return float(np.mean(scores)) if scores else np.nan

    pan["is_s2_mapped"]    = pan[PATHWAY_COL].isin(pancan_hsa_ids)
    pan["consensus_score"] = pan[PATHWAY_COL].apply(
        lambda hid: consensus_score_for_pathway(hid) if hid in pancan_hsa_ids else np.nan
    )

    n_mapped   = pan["is_s2_mapped"].sum()
    n_unmapped = (~pan["is_s2_mapped"]).sum()
    print(f"Pathways to plot: {len(pan)}  |  S2-mapped: {n_mapped}  |  Non-mapped: {n_unmapped}")

    # ── 5. Colour setup ───────────────────────────────────────────────────────
    CMAP_NAME  = "YlOrRd"
    score_vals = pan.loc[pan["is_s2_mapped"], "consensus_score"].dropna()
    vmin = score_vals.min() if len(score_vals) else 0
    vmax = score_vals.max() if len(score_vals) else 1
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

    # ── Plot 1: Scatter ───────────────────────────────────────────────────────
    fig1, ax1 = plt.subplots(figsize=(9, 7))

    grey_df = pan[~pan["is_s2_mapped"]]
    ax1.scatter(
        grey_df[DELTA_COL], grey_df["neg_log10_q"],
        color="#BBBBBB", s=50, alpha=0.6,
        linewidths=0.3, edgecolors="none",
        zorder=2,
    )

    col_df = pan[pan["is_s2_mapped"]].copy()
    col_df = col_df.sort_values("consensus_score", ascending=True, na_position="first")
    sc = ax1.scatter(
        col_df[DELTA_COL], col_df["neg_log10_q"],
        c=col_df["consensus_score"].values,
        cmap=CMAP_NAME, norm=norm,
        s=100, alpha=0.90,
        linewidths=0.5, edgecolors="dimgray",
        zorder=3,
    )

    ax1.set_xlabel("Δ$\\mu$  (pathways)", fontsize=12)
    ax1.set_ylabel("−log₁₀(q-value)", fontsize=12)
    ax1.set_title(
        f"{cancer_label} Pathways:  Δ$\\mu$  vs.  Significance\n"
        f"(Bailey et al. 2018 {cancer_label} driver pathways coloured by mean S1 consensus score)",
        fontsize=12, pad=14,
    )
    ax1.spines[["top", "right"]].set_visible(False)
    ax1.grid(alpha=0.18, linestyle="--")

    cbar = fig1.colorbar(sc, ax=ax1, pad=0.02, fraction=0.04)
    cbar.set_label("Mean Consensus Score", fontsize=10)

    ax1.axhline(-np.log10(0.05), color="#888888", linewidth=1.0, linestyle=":", alpha=0.7)
    ax1.text(
        pan[DELTA_COL].max() - 0.02, -np.log10(0.05) +
        0.08, "q = 0.05", color="#888888", fontsize=9, ha="left", va="bottom",
        bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha
        =0.85, ec="#CCCCCC")
    )

    grey_patch = mpatches.Patch(
        facecolor="#BBBBBB",
        label=f"Not in Bailey et al. 2018 {cancer_label}  (n={len(grey_df)})"
    )
    ax1.legend(handles=[grey_patch], fontsize=9, loc="lower right",
               framealpha=0.85, edgecolor="#CCCCCC")

    plt.tight_layout()
    out1 = OUT_DIR / f"scatter_{run['out_suffix']}_delta_vs_qvalue.png"
    fig1.savefig(out1, dpi=600, bbox_inches="tight")
    print(f"Saved scatter → {out1}")
    plt.close(fig1)


    # --- Simplified version for presentations -----
    fig1b, ax1b = plt.subplots(figsize=(12, 7))

    ax1b.scatter(
        grey_df[DELTA_COL], grey_df["neg_log10_q"],
        color="#BBBBBB", s=80, alpha=0.6,
        linewidths=0, edgecolors="none",
        zorder=2,
    )

    sc2 = ax1b.scatter(
        col_df[DELTA_COL], col_df["neg_log10_q"],
        c=col_df["consensus_score"].values,
        cmap=CMAP_NAME, norm=norm,
        s=160, alpha=0.90,
        linewidths=0.5, edgecolors="dimgray",
        zorder=3,
    )

    ax1b.set_xlabel("Δ$\\mu$  (pathways)", fontsize=14)
    ax1b.set_ylabel("−log₁₀(q-value)", fontsize=14)
    ax1b.tick_params(labelsize=12)
    ax1b.spines[["top", "right"]].set_visible(False)
    ax1b.grid(alpha=0.18, linestyle="--")

    cbar2 = fig1b.colorbar(sc2, ax=ax1b, pad=0.02, fraction=0.04)
    cbar2.set_label("Mean Consensus Score", fontsize=11)
    cbar2.ax.tick_params(labelsize=10)

    ax1b.axhline(-np.log10(0.05), color="#555555", linewidth=2.0, linestyle="--", alpha=0.85)
    ax1b.text(
        pan[DELTA_COL].max() - 0.02, -np.log10(0.05) + 0.08,
        "q = 0.05", color="#555555", fontsize=11, ha="left", va="bottom",
        bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.85, ec="#CCCCCC")
    )

    grey_patch2 = mpatches.Patch(
        facecolor="#BBBBBB",
        label=f"Not in Bailey et al. 2018 {cancer_label}  (n={len(grey_df)})"
    )
    ax1b.legend(handles=[grey_patch2], fontsize=10, loc="lower right",
                framealpha=0.85, edgecolor="#CCCCCC")

    plt.tight_layout()
    out1b = OUT_DIR / f"scatter_{run['out_suffix']}_delta_vs_qvalue_presentation.png"
    fig1b.savefig(out1b, dpi=300, bbox_inches="tight")
    print(f"Saved scatter (presentation) → {out1b}")
    plt.close(fig1b)


    # ── Plot 2: Violin ────────────────────────────────────────────────────────
    grp_mapped = pan.loc[ pan["is_s2_mapped"], DELTA_COL].dropna().values
    grp_other  = pan.loc[~pan["is_s2_mapped"], DELTA_COL].dropna().values

    if len(grp_mapped) == 0 or len(grp_other) == 0:
        print(f"Skipping violin for {cancer_label} — no data")
        pan.to_csv(OUT_DIR / f"scatter_{run['out_suffix']}_delta_vs_qvalue_data.csv", index=False)
        continue

    VIOLIN_COLORS = {
        f"Bailey et al. 2018\n{cancer_label} mapped": "#D65F5F",
        "Non-mapped\npathways":                        "#6C8EBF",
    }
    labels = list(VIOLIN_COLORS.keys())
    arrays = [grp_mapped, grp_other]

    fig2, ax2 = plt.subplots(figsize=(6, 6))

    parts = ax2.violinplot(arrays, positions=[0, 1],
                           showmedians=True, showextrema=False, widths=0.55)
    for pc, lbl in zip(parts["bodies"], labels):
        pc.set_facecolor(VIOLIN_COLORS[lbl])
        pc.set_alpha(0.65)
        pc.set_edgecolor("#444444")
        pc.set_linewidth(1.0)
    parts["cmedians"].set_color("#222222")
    parts["cmedians"].set_linewidth(2.0)

    for i, arr in enumerate(arrays):
        med = np.median(arr[np.isfinite(arr)])
        ax2.text(i, med, f"{med:.3f}",
                 ha="center", va="bottom", fontsize=9, fontweight="bold", color="#222222")

    rng = np.random.default_rng(42)
    for i, (arr, lbl) in enumerate(zip(arrays, labels)):
        jitter = rng.uniform(-0.14, 0.14, size=len(arr))
        ax2.scatter(i + jitter, arr,
                    color=VIOLIN_COLORS[lbl], alpha=0.5, s=15,
                    linewidths=0, zorder=0)

    if len(grp_mapped) > 0 and len(grp_other) > 0:
        _, p_mw = stats.mannwhitneyu(grp_mapped, grp_other, alternative="two-sided")
        p_str = f"p = {p_mw:.2e}" if p_mw >= 1e-4 else "p < 1e-4"
        ax2.text(0.97, 0.97,
                 f"Mann-Whitney  {p_str}\nn = {len(grp_mapped)} vs {len(grp_other)}",
                 transform=ax2.transAxes, fontsize=9, ha="right", va="top",
                 bbox=dict(boxstyle="round,pad=0.4", fc="white", alpha=0.85, ec="#CCCCCC"))

    ax2.set_xticks([0, 1])
    ax2.set_xticklabels(labels, fontsize=11)
    ax2.set_ylabel("Δ$\\mu$  (pathway distance score)", fontsize=11)
    ax2.set_title(
        f"Δ$\\mu$  Distribution  –  {cancer_label}\n"
        f"Bailey et al. 2018 {cancer_label}-mapped vs. non-mapped pathways",
        fontsize=11, pad=10,
    )
    ax2.spines[["top", "right"]].set_visible(False)
    ax2.grid(axis="y", alpha=0.25, linestyle="--")
    ax2.axhline(0, color="#888888", linewidth=1.0, linestyle=":", alpha=0.7)

    plt.tight_layout()
    out2 = OUT_DIR / f"violin_{run['out_suffix']}_mapped_vs_other.png"
    fig2.savefig(out2, dpi=600, bbox_inches="tight")
    print(f"Saved violin  → {out2}")
    plt.close(fig2)

    # ── save backing data ─────────────────────────────────────────────────────
    pan.to_csv(OUT_DIR / f"scatter_{run['out_suffix']}_delta_vs_qvalue_data.csv", index=False)

print("\nDone.")