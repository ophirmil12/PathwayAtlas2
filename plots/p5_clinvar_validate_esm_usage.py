# ESM Score Separation (ClinVar Ground Truth)
from sympy.abc import alpha

from plot_boot import *
boot_plot_folder()

import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import roc_auc_score

from definitions import (
    CLINVAR_DATA_TABLE_P,
    PLOTS_P,
    COLOR_MAP,
    set_paper_palette
)


def validate_esm_on_clinvar():
    """
    Plots the distribution of ESM scores for Benign vs Pathogenic ClinVar mutations.
    Separates by Ordered vs Disordered regions, plus a combined panel.
    Also saves a clean single-panel version (all mutations) for presentations.
    """
    set_paper_palette()

    # 1. Load ClinVar Data
    if not os.path.exists(CLINVAR_DATA_TABLE_P):
        print(f"Error: ClinVar data not found at {CLINVAR_DATA_TABLE_P}")
        return

    df = pd.read_csv(CLINVAR_DATA_TABLE_P)

    esm_col      = 'wt_not_nadav_marginals_base_wt_score'
    label_col    = 'binary_label'   # 0=Benign, 1=Pathogenic
    disorder_col = 'is_disordered'

    df = df.dropna(subset=[esm_col, label_col])

    color_map = {0: COLOR_MAP['benign'], 1: COLOR_MAP['pathogenic']}
    label_map = {0: "Benign", 1: "Pathogenic"}

    # ── helper ────────────────────────────────────────────────────────────────
    def _draw_kde_panel(ax, subset, show_n=True):
        """Draw Benign/Pathogenic KDE curves onto *ax* and return the AUC."""
        auc = roc_auc_score(subset[label_col], -subset[esm_col])
        for label_val in [0, 1]:
            data = subset[subset[label_col] == label_val][esm_col]
            label = f"{label_map[label_val]} (N={len(data)})" if show_n else label_map[label_val]
            sns.kdeplot(
                data,
                ax=ax,
                fill=True,
                color=color_map[label_val],
                label=label,
                linewidth=2,
                alpha=0.4,
            )
        ax.set_xlabel("ESM Log-Prob (Lower = More Surprising/Pathogenic)", fontsize=11)
        ax.grid(axis='y', linestyle='--', alpha=0.3)
        ax.legend()
        return auc

    # ── 3-panel figure ────────────────────────────────────────────────────────
    fig, axes = plt.subplots(1, 3, figsize=(20, 6), sharey=True)

    panel_configs = [
        (0, "Ordered Regions (is_disordered=0)",   df[df[disorder_col] == 0]),
        (1, "Disordered Regions (is_disordered=1)", df[df[disorder_col] == 1]),
        (2, "All Mutations (Ordered + Disordered)", df),
    ]

    for ax_idx, title, subset in panel_configs:
        ax = axes[ax_idx]
        if subset.empty:
            ax.set_visible(False)
            continue

        auc = _draw_kde_panel(ax, subset)
        ax.set_title(f"{title}\nSeparation AUC: {auc:.3f}", fontsize=13, fontweight='bold')
        ax.set_ylabel("Density" if ax_idx == 0 else "", fontsize=11)

    plt.suptitle("ESM Score Separation: Benign vs Pathogenic (ClinVar)", fontsize=16, y=1.02)
    plt.tight_layout()

    os.makedirs(PLOTS_P, exist_ok=True)
    save_path = os.path.join(PLOTS_P, "p5_clinvar_validate_esm_usage.png")
    plt.savefig(save_path, dpi=600, bbox_inches='tight')
    print(f"3-panel plot saved to: {save_path}")

    # ── single-panel presentation figure (all mutations, no suptitle) ─────────
    fig2, ax2 = plt.subplots(figsize=(8, 6))

    auc_all = _draw_kde_panel(ax2, df, show_n=False)
    ax2.set_title(None)
    ax2.set_xlabel("ESM LLR", fontsize=22)
    ax2.set_ylabel("Density", fontsize=22)
    ax2.legend(fontsize=22)
    plt.tight_layout()

    save_path_simple = os.path.join(PLOTS_P, "p5_clinvar_validate_esm_usage_simple.png")
    plt.savefig(save_path_simple, dpi=600, bbox_inches='tight')
    print(f"Simple plot saved to:  {save_path_simple}")

    print(f"Total mutations analyzed: {len(df)}")

    plt.close()

    # ── sigmoid presentation figure (one for all ordered/disordered) ──────────────────────────────
    from sklearn.linear_model import LogisticRegression
    import numpy as np

    rng = np.random.default_rng(42)

    X = df[[esm_col]].values
    y = df[label_col].values

    clf = LogisticRegression()
    clf.fit(X, y)

    fig3, ax3 = plt.subplots(figsize=(10, 6))

    x_range = np.linspace(df[esm_col].min(), df[esm_col].max(), 300).reshape(-1, 1)
    y_prob = clf.predict_proba(x_range)[:, 1]

    ax3.plot(x_range, y_prob, color='steelblue', linewidth=2.5, label='Logistic fit', zorder=3, alpha=0.7)

    for label_val in [0, 1]:
        subset = df[df[label_col] == label_val]
        sampled = subset.sample(n=min(200, len(subset)), random_state=42)
        probs = clf.predict_proba(sampled[[esm_col]].values)[:, 1]
        jitter = rng.uniform(-0.05, 0.05, size=len(sampled))
        x_vals = sampled[esm_col].values

        if label_val == 1:  # pathogenic (red)
            for i, (mask, z) in enumerate([(x_vals < -10, 5), (x_vals >= -10, 2)]):
                ax3.scatter(
                    x_vals[mask],
                    probs[mask] + jitter[mask],
                    color=color_map[label_val],
                    alpha=0.8,
                    s=50,
                    zorder=z,
                    label=label_map[label_val] if i == 0 else '_nolegend_',
                )
        else:  # benign (green)
            ax3.scatter(
                x_vals,
                probs + jitter,
                color=color_map[label_val],
                alpha=0.8,
                s=60,
                label=label_map[label_val],
                zorder=4,
            )

    ax3.set_xlabel("ESM LLR", fontsize=18)
    ax3.set_ylabel("Pathogenicity probability", fontsize=18)
    ax3.set_ylim(-0.15, 1.15)
    ax3.set_xlim(-20, 0)
    ax3.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax3.grid(axis='y', linestyle='--', alpha=0.4)
    ax3.legend(fontsize=22)
    plt.tight_layout()

    save_path_sigmoid = os.path.join(PLOTS_P, "p5_clinvar_logistic_sigmoid.png")
    plt.savefig(save_path_sigmoid, dpi=600, bbox_inches='tight')
    print(f"Sigmoid plot saved to: {save_path_sigmoid}")

    plt.close()


    # ── sigmoid presentation figure (one for each, ordered/disordered) ──────────────────────────────
    for disorder_val, disorder_label in [(0, "Ordered"), (1, "Disordered")]:
        subset = df[df[disorder_col] == disorder_val]

        if subset.empty:
            print(f"No data for {disorder_label} regions, skipping sigmoid plot.")
            continue

        X_sub = subset[[esm_col]].values
        y_sub = subset[label_col].values

        clf_sub = LogisticRegression()
        clf_sub.fit(X_sub, y_sub)

        fig_sub, ax_sub = plt.subplots(figsize=(10, 6))

        x_range_sub = np.linspace(subset[esm_col].min(), subset[esm_col].max(), 300).reshape(-1, 1)
        y_prob_sub = clf_sub.predict_proba(x_range_sub)[:, 1]

        ax_sub.plot(x_range_sub, y_prob_sub, color='steelblue', linewidth=2.5, label='Logistic fit', zorder=3, alpha=0.7)

        for label_val in [0, 1]:
            subset_label = subset[subset[label_col] == label_val]
            sampled_label = subset_label.sample(n=min(200, len(subset_label)), random_state=42)
            probs_label = clf_sub.predict_proba(sampled_label[[esm_col]].values)[:, 1]
            jitter_label = rng.uniform(-0.05, 0.05, size=len(sampled_label))
            x_vals_label = sampled_label[esm_col].values

            if label_val == 1:
                for i, (mask, z) in enumerate([(x_vals_label < -10, 5), (x_vals_label >= -10, 2)]):
                    ax_sub.scatter(
                        x_vals_label[mask],
                        probs_label[mask] + jitter_label[mask],
                        color=color_map[label_val],
                        alpha=0.8,
                        s=50,
                        zorder=z,
                        label=label_map[label_val] if i == 0 else '_nolegend_',
                    )
            else:
                ax_sub.scatter(
                    x_vals_label,
                    probs_label + jitter_label,
                    color=color_map[label_val],
                    alpha=0.8,
                    s=60,
                    label=label_map[label_val],
                    zorder=4,
                )

        ax_sub.set_xlabel("ESM LLR", fontsize=18)
        ax_sub.set_ylabel("Pathogenicity probability", fontsize=18)
        ax_sub.set_ylim(-0.15, 1.15)
        ax_sub.set_xlim(-20, 0)
        ax_sub.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
        ax_sub.grid(axis='y', linestyle='--', alpha=0.4)
        ax_sub.legend(fontsize=22)
        plt.tight_layout()
        save_path_sub = os.path.join(PLOTS_P, f"p5_clinvar_logistic_sigmoid_{disorder_label.lower()}.png")
        plt.savefig(save_path_sub, dpi=600, bbox_inches='tight')
        print(f"Sigmoid plot for {disorder_label} regions saved to: {save_path_sub}")
        plt.close()  # ← moved inside the loop

    # ── side-by-side version: ordered | disordered, shared y-axis & legend ───
    fig_sb, axes_sb = plt.subplots(1, 2, figsize=(18, 6), sharey=True)

    for ax_sb, (disorder_val, disorder_label) in zip(axes_sb, [(0, "Ordered"), (1, "Disordered")]):
        subset = df[df[disorder_col] == disorder_val]
        if subset.empty:
            ax_sb.set_visible(False)
            continue

        clf_sb = LogisticRegression()
        clf_sb.fit(subset[[esm_col]].values, subset[label_col].values)

        x_range_sb = np.linspace(subset[esm_col].min(), subset[esm_col].max(), 300).reshape(-1, 1)
        y_prob_sb  = clf_sb.predict_proba(x_range_sb)[:, 1]

        ax_sb.plot(x_range_sb, y_prob_sb, color='steelblue', linewidth=2.5,
                   label='Logistic fit', zorder=3, alpha=0.7)

        for label_val in [0, 1]:
            sub_lbl   = subset[subset[label_col] == label_val]
            sampled   = sub_lbl.sample(n=min(200, len(sub_lbl)), random_state=42)
            probs     = clf_sb.predict_proba(sampled[[esm_col]].values)[:, 1]
            jitter    = rng.uniform(-0.05, 0.05, size=len(sampled))
            x_vals    = sampled[esm_col].values

            if label_val == 1:
                for i, (mask, z) in enumerate([(x_vals < -10, 5), (x_vals >= -10, 2)]):
                    ax_sb.scatter(x_vals[mask], probs[mask] + jitter[mask],
                                  color=color_map[label_val], alpha=0.8, s=50, zorder=z,
                                  label=label_map[label_val] if i == 0 else '_nolegend_')
            else:
                ax_sb.scatter(x_vals, probs + jitter,
                              color=color_map[label_val], alpha=0.8, s=60,
                              label=label_map[label_val], zorder=4)

        ax_sb.set_title(None)
        ax_sb.set_xlabel("")
        ax_sb.set_xlim(-20, 0)
        ax_sb.set_xticks(np.arange(-20, 1, 5))
        ax_sb.set_xticklabels([])
        ax_sb.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
        ax_sb.set_yticklabels([])
        ax_sb.tick_params(left=False, bottom=False)
        ax_sb.grid(axis='y', linestyle='--', alpha=0.4)
        for spine in ax_sb.spines.values():
            spine.set_visible(False)

    axes_sb[0].set_ylabel("Pathogenicity probability", fontsize=22)

    # Single x-axis label centred across both panels
    fig_sb.text(0.5, -0.02, "ESM LLR", ha='center', va='top', fontsize=22)

    # Legend on the upper-right of the right panel only
    handles, labels = axes_sb[1].get_legend_handles_labels()
    axes_sb[1].legend(handles, labels, fontsize=23, loc='upper right', frameon=False)

    plt.tight_layout()
    save_path_sb = os.path.join(PLOTS_P, "p5_clinvar_logistic_sigmoid_side_by_side.png")
    plt.savefig(save_path_sb, dpi=600, bbox_inches='tight')
    plt.close()
    print(f"Side-by-side sigmoid plot saved to: {save_path_sb}")



if __name__ == "__main__":
    validate_esm_on_clinvar()