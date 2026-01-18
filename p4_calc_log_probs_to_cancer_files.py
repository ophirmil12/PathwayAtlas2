# Calculate the log-probs (WT-marginals) and put into the CSVs of the cancer data


import os
import sys
import torch
import esm

from definitions import *
from p4_esm_emb_and_log_probs import ScoringCalculator


def run_cancer_log_prob_calculation(file_index):
    """
    Applies ESM scoring to a single cancer mutation file based on the provided index.
    Designed for Slurm Array tasks.
    """
    # 1. Setup Environment
    device = 'cuda' if torch.cuda.is_available() else 'cpu'

    # 2. Identify target file
    cancer_files = sorted([f for f in os.listdir(CBIO_MUTATION_STUDIES_WITH_SEQUENCES_P) if f.endswith('.csv')])

    if file_index < 0 or file_index >= len(cancer_files):
        print(f"Error: File index {file_index} is out of range (Total files: {len(cancer_files)}).")
        sys.exit(1)

    filename = cancer_files[file_index]
    csv_path = os.path.join(CBIO_MUTATION_STUDIES_WITH_SEQUENCES_P, filename)

    # 3. Load ESM Model and Alphabet
    print(f"[{file_index}] Processing: {filename}")
    print(f"Loading ESM model ({ESM1B_MODEL})...")
    model, alphabet = esm.pretrained.load_model_and_alphabet(ESM1B_MODEL)
    calculator = ScoringCalculator(model, alphabet)

    # 4. Process the file
    try:
        # This method groups by KeggId, loads .pt files,
        # calculates log P(mut) - log P(wt), and updates the CSV.
        calculator.handle_cancer_csv_optimized(csv_path, recalc_scores=False)       # HERE IT HAPPENS
        print(f"Successfully processed {filename}")
    except Exception as e:
        print(f"Failed to process {filename}: {e}")
        sys.exit(1)


if __name__ == "__main__":
    # Check for required index argument
    if len(sys.argv) < 2:
        print("Usage: python p4_calc_log_probs_to_cancer_files.py <file_index>")
        sys.exit(1)

    try:
        target_idx = int(sys.argv[1])
    except ValueError:
        print("Error: file_index must be an integer.")
        sys.exit(1)

    # no_grad is vital to keep VRAM usage low
    with torch.no_grad():
        run_cancer_log_prob_calculation(target_idx)