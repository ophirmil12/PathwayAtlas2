import os
import glob
import random
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

# Import paths and palette from your definitions
from definitions import KEGG_GENE_SCORES_P, MY_PALETTE, set_paper_palette

# Hardcoded Evolutionary Substitution Matrix (Dataset S3-revised)
# Format: AA_SUBSTITUTION_MATRIX[RefAA][AltAA] = Probability
AA_SUBSTITUTION_MATRIX = {
    'F': {'F': 0.999437, 'L': 0.000268, 'M': 0.0, 'I': 2.05e-05, 'V': 2.91e-05, 'S': 0.000195, 'P': 0.0, 'T': 0.0, 'A': 0.0, 'Y': 2.05e-05, 'H': 0.0, 'Q': 0.0, 'N': 0.0, 'K': 0.0, 'D': 0.0, 'E': 0.0, 'C': 2.91e-05, 'W': 0.0, 'R': 0.0, 'G': 0.0},
    'L': {'F': 9.53e-05, 'L': 0.999589, 'M': 1.73e-05, 'I': 1.75e-05, 'V': 4.45e-05, 'S': 5.88e-05, 'P': 0.000136, 'T': 0.0, 'A': 0.0, 'Y': 0.0, 'H': 5.86e-06, 'Q': 8.46e-06, 'N': 0.0, 'K': 0.0, 'D': 0.0, 'E': 0.0, 'C': 0.0, 'W': 5.82e-06, 'R': 2.03e-05, 'G': 0.0},
    'M': {'F': 0.0, 'L': 5.34e-05, 'M': 0.999203, 'I': 0.000321, 'V': 0.000176, 'S': 0.0, 'P': 0.0, 'T': 0.000195, 'A': 0.0, 'Y': 0.0, 'H': 0.0, 'Q': 0.0, 'N': 0.0, 'K': 2.05e-05, 'D': 0.0, 'E': 0.0, 'C': 0.0, 'W': 0.0, 'R': 2.91e-05, 'G': 0.0},
    'I': {'F': 1.29e-05, 'L': 4.05e-05, 'M': 8.27e-05, 'I': 0.999442, 'V': 0.000176, 'S': 2.02e-05, 'P': 0.0, 'T': 0.000195, 'A': 0.0, 'Y': 0.0, 'H': 0.0, 'Q': 0.0, 'N': 1.42e-05, 'K': 6.29e-06, 'D': 0.0, 'E': 0.0, 'C': 0.0, 'W': 0.0, 'R': 8.94e-06, 'G': 0.0},
    'V': {'F': 1.23e-05, 'L': 6.55e-05, 'M': 0.000113, 'I': 0.000129, 'V': 0.999433, 'S': 0.0, 'P': 0.0, 'T': 0.0, 'A': 0.000195, 'Y': 0.0, 'H': 0.0, 'Q': 0.0, 'N': 0.0, 'K': 0.0, 'D': 7.82e-06, 'E': 1.27e-05, 'C': 0.0, 'W': 0.0, 'R': 0.0, 'G': 2.91e-05},
    'S': {'F': 9.08e-05, 'L': 7.26e-05, 'M': 0.0, 'I': 1.2e-05, 'V': 0.0, 'S': 0.999413, 'P': 0.000122, 'T': 2.98e-05, 'A': 1.83e-05, 'Y': 1.43e-05, 'H': 0.0, 'Q': 0.0, 'N': 9.05e-05, 'K': 0.0, 'D': 0.0, 'E': 0.0, 'C': 2.48e-05, 'W': 3.34e-06, 'R': 4.15e-05, 'G': 6.56e-05},
    'P': {'F': 0.0, 'L': 0.000259, 'M': 0.0, 'I': 0.0, 'V': 0.0, 'S': 0.000259, 'P': 0.999295, 'T': 4.1e-05, 'A': 5.11e-05, 'Y': 0.0, 'H': 2.14e-05, 'Q': 1.97e-05, 'N': 0.0, 'K': 0.0, 'D': 0.0, 'E': 0.0, 'C': 0.0, 'W': 0.0, 'R': 5.11e-05, 'G': 0.0},
    'T': {'F': 0.0, 'L': 0.0, 'M': 3.57e-05, 'I': 0.000224, 'V': 0.0, 'S': 4.53e-05, 'P': 3.48e-05, 'T': 0.999417, 'A': 0.000176, 'Y': 0.0, 'H': 0.0, 'Q': 0.0, 'N': 2.14e-05, 'K': 1.96e-05, 'D': 0.0, 'E': 0.0, 'C': 0.0, 'W': 0.0, 'R': 2.44e-05, 'G': 0.0},
    'A': {'F': 0.0, 'L': 0.0, 'M': 0.0, 'I': 0.0, 'V': 0.000259, 'S': 3.23e-05, 'P': 4.55e-05, 'T': 0.000243, 'A': 0.999326, 'Y': 0.0, 'H': 0.0, 'Q': 0.0, 'N': 0.0, 'K': 0.0, 'D': 2.33e-05, 'E': 1.77e-05, 'C': 0.0, 'W': 0.0, 'R': 0.0, 'G': 5.11e-05},
    'Y': {'F': 1.87e-05, 'L': 0.0, 'M': 0.0, 'I': 0.0, 'V': 0.0, 'S': 3.48e-05, 'P': 0.0, 'T': 0.0, 'A': 0.0, 'Y': 0.999525, 'H': 0.000195, 'Q': 0.0, 'N': 2.05e-05, 'K': 0.0, 'D': 2.91e-05, 'E': 0.0, 'C': 0.000176, 'W': 0.0, 'R': 0.0, 'G': 0.0},
    'H': {'F': 0.0, 'L': 1.87e-05, 'M': 0.0, 'I': 0.0, 'V': 0.0, 'S': 0.0, 'P': 3.48e-05, 'T': 0.0, 'A': 0.0, 'Y': 0.000259, 'H': 0.999343, 'Q': 7.43e-05, 'N': 4.1e-05, 'K': 0.0, 'D': 5.11e-05, 'E': 0.0, 'C': 0.0, 'W': 0.0, 'R': 0.000176, 'G': 0.0},
    'Q': {'F': 0.0, 'L': 1.87e-05, 'M': 0.0, 'I': 0.0, 'V': 0.0, 'S': 0.0, 'P': 3.48e-05, 'T': 0.0, 'A': 0.0, 'Y': 0.0, 'H': 7.02e-05, 'Q': 0.999607, 'N': 0.0, 'K': 4.1e-05, 'D': 0.0, 'E': 5.11e-05, 'C': 0.0, 'W': 0.0, 'R': 0.000176, 'G': 0.0},
    'N': {'F': 0.0, 'L': 0.0, 'M': 0.0, 'I': 1.87e-05, 'V': 0.0, 'S': 0.000176, 'P': 0.0, 'T': 3.48e-05, 'A': 0.0, 'Y': 1.87e-05, 'H': 3.48e-05, 'Q': 0.0, 'N': 0.999466, 'K': 7.3e-05, 'D': 0.000176, 'E': 0.0, 'C': 0.0, 'W': 0.0, 'R': 0.0, 'G': 0.0},
    'K': {'F': 0.0, 'L': 0.0, 'M': 1.05e-05, 'I': 8.17e-06, 'V': 0.0, 'S': 0.0, 'P': 0.0, 'T': 3.48e-05, 'A': 0.0, 'Y': 0.0, 'H': 0.0, 'Q': 3.48e-05, 'N': 6.72e-05, 'K': 0.999491, 'D': 0.0, 'E': 0.000176, 'C': 0.0, 'W': 0.0, 'R': 0.000176, 'G': 0.0},
    'D': {'F': 0.0, 'L': 0.0, 'M': 0.0, 'I': 0.0, 'V': 1.87e-05, 'S': 0.0, 'P': 0.0, 'T': 0.0, 'A': 3.48e-05, 'Y': 3.23e-05, 'H': 4.55e-05, 'Q': 0.0, 'N': 0.000243, 'K': 0.0, 'D': 0.999374, 'E': 7.33e-05, 'C': 0.0, 'W': 0.0, 'R': 0.0, 'G': 0.000176},
    'E': {'F': 0.0, 'L': 0.0, 'M': 0.0, 'I': 0.0, 'V': 1.87e-05, 'S': 0.0, 'P': 0.0, 'T': 0.0, 'A': 3.48e-05, 'Y': 0.0, 'H': 0.0, 'Q': 4.56e-05, 'N': 0.0, 'K': 0.000243, 'E': 0.999413, 'D': 6.72e-05, 'C': 0.0, 'W': 0.0, 'R': 0.0, 'G': 0.000176},
    'C': {'F': 3.23e-05, 'L': 0.0, 'M': 0.0, 'I': 0.0, 'V': 0.0, 'S': 6.6e-05, 'P': 0.0, 'T': 0.0, 'A': 0.0, 'Y': 0.000243, 'H': 0.0, 'Q': 0.0, 'N': 0.0, 'K': 0.0, 'D': 0.0, 'E': 0.0, 'C': 0.99939, 'W': 4.28e-05, 'R': 0.000195, 'G': 2.91e-05},
    'W': {'F': 0.0, 'L': 3.23e-05, 'M': 0.0, 'I': 0.0, 'V': 0.0, 'S': 4.55e-05, 'P': 0.0, 'T': 0.0, 'A': 0.0, 'Y': 0.0, 'H': 0.0, 'Q': 0.0, 'N': 0.0, 'K': 0.0, 'D': 0.0, 'E': 0.0, 'C': 7.79e-05, 'W': 0.999599, 'R': 0.000215, 'G': 2.91e-05},
    'R': {'F': 0.0, 'L': 1.08e-05, 'M': 1.1e-05, 'I': 1.05e-05, 'V': 0.0, 'S': 5e-05, 'P': 1.52e-05, 'T': 3.03e-05, 'A': 0.0, 'Y': 0.0, 'H': 3.65e-05, 'Q': 4.5e-05, 'N': 0.0, 'K': 0.000162, 'D': 0.0, 'E': 0.0, 'C': 3.89e-05, 'W': 4.11e-05, 'R': 0.999413, 'G': 0.000134},
    'G': {'F': 0.0, 'L': 0.0, 'M': 0.0, 'I': 0.0, 'V': 3.23e-05, 'S': 0.000102, 'P': 0.0, 'T': 0.0, 'A': 4.55e-05, 'Y': 0.0, 'H': 0.0, 'Q': 0.0, 'N': 0.0, 'K': 0.0, 'D': 0.000102, 'E': 0.000141, 'C': 1.36e-05, 'W': 9.59e-06, 'R': 0.000187, 'G': 0.999365}
}

def run_pathogenic_prob_biological_validation():
    """
    Validates pathogenic_prob scores against an evolutionary substitution matrix.
    Checks if higher evolutionary probability correlates with more benign (higher) pathogenic_prob scores.
    """
    set_paper_palette()

    # 1. Collect Data Pairs from a sample of gene files
    all_files = glob.glob(os.path.join(KEGG_GENE_SCORES_P, "*.csv"))
    # Sampling 1000 genes to get a statistically significant population
    sample_files = random.sample(all_files, min(1000, len(all_files)))

    matrix_probs = []
    pathogenic_prob_scores = []

    print(f"Sampling {len(sample_files)} genes for biological validation...")

    for f in sample_files:
        try:
            df = pd.read_csv(f)
            # Ensure columns generated by previous pipeline steps exist
            if 'Variant' not in df.columns or 'pathogenic_prob' not in df.columns:
                continue

            # Drop rows where pathogenic_prob score computation failed
            df = df.dropna(subset=['pathogenic_prob'])

            for _, row in df.iterrows():
                variant = str(row['Variant'])
                # Variant format: [RefAA][Index][AltAA] (e.g., M1L)
                ref_aa = variant[0]
                alt_aa = variant[-1]

                # Dictionary Lookup
                prob = AA_SUBSTITUTION_MATRIX.get(ref_aa, {}).get(alt_aa, 0.0)

                # We analyze non-synonymous mutations (Ref != Alt)
                # and filter out 0 probabilities to allow log-scaling in the plot
                if ref_aa != alt_aa and prob > 0:
                    matrix_probs.append(prob)
                    pathogenic_prob_scores.append(row['pathogenic_prob'])
        except Exception:
            continue

    if not matrix_probs:
        print("Error: No valid mutation data found in the sampled gene files.")
        return

    # 2. Calculate Spearman Rank Correlation
    # Spearman is used because the relationship is expected to be monotonic
    # (more common -> less pathogenic) but not necessarily linear in scale.
    corr, p_val = spearmanr(matrix_probs, pathogenic_prob_scores)

    print("\n" + "="*40)
    print("BIOLOGICAL VALIDATION RESULTS")
    print(f"Spearman Correlation (r): {corr:.4f}")
    print(f"P-value:                  {p_val:.2e}")
    print(f"Mutations Analyzed:       {len(matrix_probs)}")
    print("="*40)

    # 3. Plotting
    plt.figure(figsize=(10, 6))

    # Scatter plot using 'The Green' from MY_PALETTE (#447D68)
    plt.scatter(matrix_probs, pathogenic_prob_scores, alpha=0.1, color=MY_PALETTE[1], s=8)

    # Log scale for X-axis (Mutation probabilities span several orders of magnitude)
    plt.xscale('log')

    plt.title("Biological Validation: pathogenic_prob Scores vs. Evolutionary Substitution Matrix", fontsize=14)
    plt.xlabel("AA Substitution Probability (Evolutionary Matrix)", fontsize=12)
    plt.ylabel("pathogenic_prob (WT-Marginal Score after ClinVar)", fontsize=12)

    # Calculate and plot log-linear trend line
    try:
        log_probs = np.log10(matrix_probs)
        m, b = np.polyfit(log_probs, pathogenic_prob_scores, 1)
        plt.plot(matrix_probs, m * log_probs + b, color=MY_PALETTE[0],
                 label=f'Trend (r={corr:.2f})', linewidth=2)
    except Exception as e:
        print(f"Warning: Could not plot trend line: {e}")

    plt.grid(True, which="both", ls="-", alpha=0.2)
    plt.legend()
    plt.tight_layout()

    plt.savefig("p5_plot_pathogenic_prob_vs_evolutionary_matrix.png", dpi=300)
    print("\nSuccess: Plot saved as p4_plot_pathogenic_prob_vs_evolutionary_matrix.png")
    plt.show()

if __name__ == "__main__":
    run_pathogenic_prob_biological_validation()