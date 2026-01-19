# To visualize how a specific cancer signature (the PSSM)
# alters our expectation of pathogenicity for a pathway,
# we will overlay the Uniform Background (every possible mutation is equally likely)
# and the PSSM-Weighted Background (mutations are weighted by their cancer-specific frequency).
#
# This plot essentially answers:
# "Does this specific cancer's mutation profile make this pathway
# look more pathogenic even before we look at actual patient data?"

from plot_boot import *
boot_plot_folder()

import glob
import pandas as pd
import numpy as np
from tqdm import tqdm
from definitions import *
from distance_utils import get_bg_histogram_after_pssm


def get_raw_background_histogram(pathway_id: str, bins=NUMBER_OF_BINS):
    """
    Calculates the distribution of pathogenicity scores assuming
    all mutation types are equally likely (No PSSM).
    """
    p_id = pathway_id.replace(":", "_")
    file_path = os.path.join(KEGG_PATHWAY_SCORES_P, f"{p_id}.csv")

    if not os.path.exists(file_path):
        return None

    try:
        # Load only necessary columns
        df = pd.read_csv(file_path, usecols=["pathogenic_prob"])
        scores = df["pathogenic_prob"].dropna().values

        if len(scores) == 0:
            return np.zeros(bins)

        bin_edges = np.linspace(0, 1, bins + 1)
        counts, _ = np.histogram(scores, bins=bin_edges)

        # Normalize to sum to 1
        return counts.astype(float) / counts.sum()
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return None


def plot_pssm_impact(pathway_id: str, output_dir: str):
    """
    Plots the comparison between Raw and PSSM-weighted background distributions.
    """
    # 1. Get both histograms
    hist_raw = get_raw_background_histogram(pathway_id)
    # Using the PSSM defined in definitions.py
    hist_pssm = get_bg_histogram_after_pssm(pathway_id, pssm_matrix=MICHAL_HN1_PSSM)

    if hist_raw is None or np.sum(hist_pssm) == 0:
        return

    # 2. Calculate Means for annotation
    bin_edges = np.linspace(0, 1, NUMBER_OF_BINS + 1)
    bin_centers = bin_edges[:-1] + np.diff(bin_edges) / 2

    mean_raw = np.sum(bin_centers * hist_raw)
    mean_pssm = np.sum(bin_centers * hist_pssm)
    shift = mean_pssm - mean_raw

    # 3. Plotting
    plt.figure(figsize=(10, 6))

    # Use index 3 (Slate) and index 2 (Orange) from your palette
    plt.plot(bin_centers, hist_raw, label="Uniform Background (Equal Probabilities)",
             color=MY_PALETTE[3], linewidth=2, alpha=0.8)
    plt.fill_between(bin_centers, hist_raw, color=MY_PALETTE[3], alpha=0.1)

    plt.plot(bin_centers, hist_pssm, label="PSSM-Weighted Background (Cancer Signature)",
             color=MY_PALETTE[2], linewidth=2.5)
    plt.fill_between(bin_centers, hist_pssm, color=MY_PALETTE[2], alpha=0.2)

    # Add vertical lines for means
    plt.axvline(mean_raw, color=MY_PALETTE[3], linestyle='--', alpha=0.6)
    plt.axvline(mean_pssm, color=MY_PALETTE[2], linestyle='--', alpha=0.6)

    # 4. Formatting
    plt.title(f"Impact of Cancer Signature on Expected Pathogenicity\nPathway: {pathway_id}", fontsize=14)
    plt.xlabel("Pathogenicity Probability", fontsize=12)
    plt.ylabel("Probability Mass", fontsize=12)

    # Text annotation of the shift
    direction = "More Pathogenic" if shift > 0 else "More Benign"
    plt.annotate(f"Signature Shift: {shift:+.4f}\n({direction})",
                 xy=(0.05, 0.8), xycoords='axes fraction',
                 bbox=dict(boxstyle="round", fc="white", alpha=0.8),
                 fontsize=10, fontweight='bold')

    plt.legend(loc='upper right')
    plt.grid(axis='y', linestyle=':', alpha=0.4)
    plt.xlim(0, 1)

    # Save and Close
    save_path = os.path.join(output_dir, f"{pathway_id.replace(':', '_')}.png")
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()  # CRITICAL: Close the plot to free memory in the loop


if __name__ == "__main__":
    set_paper_palette()

    # Define output directory
    impact_plots_dir = os.path.join(PLOTS_P, "p6_pssm_impact_all")
    os.makedirs(impact_plots_dir, exist_ok=True)

    # Find all pathway files generated in previous steps
    pathway_files = glob.glob(os.path.join(KEGG_PATHWAY_SCORES_P, "*.csv"))

    print(f"Starting PSSM impact analysis for {len(pathway_files)} pathways...")

    for file_path in tqdm(pathway_files, desc="Generating PSSM Impact Plots"):
        # Extract ID (e.g., 'hsa00010' from 'hsa00010.csv')
        pw_id = os.path.basename(file_path).replace(".csv", "")

        try:
            plot_pssm_impact(pw_id, impact_plots_dir)
        except Exception as e:
            print(f"\n[Error] Failed to plot {pw_id}: {e}")
            continue

    print(f"\nDone! All plots saved to: {impact_plots_dir}")