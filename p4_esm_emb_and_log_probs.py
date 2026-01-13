# This part only makes the embeddings, save them, and calculate the log-probs, save them to snv/cancer csv


import os
import time
import numpy as np
import pandas as pd
from torch.nn.functional import log_softmax
import torch
import re

from definitions import *
from kegg_api import KeggGene


class ScoringCalculator:
    def __init__(self, model, alphabet, testing=False):
        """
        Initialize the calculator with the model and alphabet.

        Args:
            model - the ESM model
            alphabet - the alphabet, do:
            model, alphabet = pretrained.load_model_and_alphabet(ESM1B_MODEL)
            testing - test mode: without using esm or alphabet, only statics
        """
        if not testing:
            self.device = 'cuda' if torch.cuda.is_available() else 'cpu'
            self.model = model.to(self.device)
            self.alphabet = alphabet
            self.tokenizer = alphabet.get_batch_converter()
            self.model.eval()

            # Map the 20 standard AAs (indices 4-23 in ESM) to indices 0-19
            self.aa_to_idx = {
                ESM_1B_LOGITS_INDEX_TO_AA[i]: i - 4
                for i in range(4, 24)
            }

            print("Model initialized on: ", self.device)
            print("Alphabet indexes from model:       ", [(i, alphabet.get_tok(i)) for i in range(4, 24)])
            print("Alphabet indexes from definitions: ", [(i, ESM_1B_LOGITS_INDEX_TO_AA[i]) for i in range(4, 24)])

    def _get_batch_tokens_from_seq(self, sequence: str, name='WT') -> torch.Tensor:
        """Tokenizes the input sequence using the ESM tokenizer."""
        batch = [(name, sequence)]
        _, _, tokens = self.tokenizer(batch)
        return tokens.to(self.device)

    def _get_batch_tokens_from_kegg_id(self, kegg_id: str, name='WT') -> torch.Tensor:
        """Fetches sequence from KeggGene and tokenizes it."""
        kegg_gene = KeggGene(kegg_id)
        sequence = kegg_gene.aa_seq
        return self._get_batch_tokens_from_seq(sequence, name=kegg_id)

    def _get_logits(self, batch_tokens: torch.Tensor) -> torch.Tensor:
        """Returns raw logits from the model."""
        self.model.eval()
        with torch.no_grad():
            out = self.model(
                batch_tokens,
                repr_layers=[33],
                return_contacts=False,
            )
            logits = out['logits']
        return logits

    def _get_aa_logits(self, batch_tokens: torch.Tensor) -> torch.Tensor:
        """Returns the part of the AA in the logits (removes start/end tokens)."""
        # [batch, seq_len, vocab] -> [seq_len, 20]
        return self._get_logits(batch_tokens)[0, 1:-1, 4:24]

    def handle_long_protein(self, sequence: str) -> torch.Tensor:
        """
        Handles proteins longer than ESM_MAX_LENGTH using a sliding window.
        Each residue is assigned logits from exactly one window:
          - First window: from start to (midpoint + half stride)
          - Middle windows: exactly `stride` residues each
          - Final window (when it reaches the end): take everything until the end
        """
        seq_len = len(sequence)
        window_size = ESM_MAX_LENGTH
        stride = SLIDING_WINDOW_STRIDE
        half_window = window_size // 2
        half_stride = stride // 2

        logits_full = torch.zeros(seq_len, 20, device=self.device)

        pos = 0
        window_idx = 0

        for start in range(0, seq_len, stride):
            end = min(start + window_size, seq_len)
            subseq = sequence[start:end]  # subsequence for this window

            # Run model on this window
            batch_tokens = self._get_batch_tokens_from_seq(subseq, f"WT_window_{window_idx}")
            window_logits = self._get_aa_logits(batch_tokens)  # [window_len, 20]
            window_len = window_logits.shape[0]

            if window_idx == 0:
                # First window: take first half_window + half_stride residues
                take_len = min(half_window + half_stride, window_len)
                logits_full[0:take_len] = window_logits[0:take_len]
                pos = take_len

            elif end == seq_len:
                # Final window: take remaining residues to the end
                take_len = seq_len - pos
                logits_full[pos:] = window_logits[-take_len:]
                pos = seq_len
                break

            else:
                # Middle windows: take central `stride` residues
                center = window_len // 2
                half_stride = stride // 2
                start_idx = max(0, center - half_stride)
                end_idx = min(window_len, center + half_stride)
                take_len = end_idx - start_idx

                logits_full[pos:pos + take_len] = window_logits[start_idx:end_idx]
                pos += take_len

            window_idx += 1

        return logits_full

    def save_logits(self, kegg_id: str, logits: torch.Tensor):
        """Saves the computed logits tensor to disk."""
        file_path = os.path.join(ESM_EMBEDDINGS_P, f"{kegg_id}.pt")
        # Ensure tensor is on CPU before saving to save GPU memory/compatibility
        torch.save(logits.cpu(), file_path)
        # print(f"Saved logits for {kegg_id} to {file_path}")

    def get_or_compute_logits(self, kegg_id: str) -> tuple[torch.Tensor, str]:
        """
        1. Checks if {kegg_id}.pt exists in ESM_EMBEDDINGS_P.
        2. If yes, loads and returns it.
        3. If no, computes it (handling long sequences), saves it, and returns it.

        Returns:
            logits (Tensor): [seq_len, 20]
            sequence (str): The amino acid sequence corresponding to the logits
        """
        file_path = os.path.join(ESM_EMBEDDINGS_P, f"{kegg_id}.pt")
        kegg_gene = KeggGene(kegg_id)
        sequence = kegg_gene.aa_seq

        # 1. Try Loading
        if os.path.exists(file_path):
            try:
                # Load map_location='cpu' to avoid OOM if loading many files, move to device later
                logits = torch.load(file_path, map_location=self.device)

                # Sanity check: length match
                if logits.shape[0] != len(sequence):
                    print(
                        f"[Warning] Loaded logits length ({logits.shape[0]}) != sequence length ({len(sequence)}) for {kegg_id}. Recomputing.")
                    # Fall through to computation logic
                else:
                    return logits, sequence
            except Exception as e:
                print(f"[Error] Failed to load existing logits for {kegg_id}: {e}. Recomputing.")

        # 2. Compute
        # print(f"Computing logits for {kegg_id} (len {len(sequence)})...")
        if len(sequence) > ESM_MAX_LENGTH:
            logits = self.handle_long_protein(sequence)
        else:
            batch_tokens = self._get_batch_tokens_from_seq(sequence, kegg_id)
            logits = self._get_aa_logits(batch_tokens)

        # 3. Save
        self.save_logits(kegg_id, logits)

        return logits, sequence

    def score_all_mutations(self, kegg_id: str) -> torch.Tensor:
        """
        Scores all single amino acid substitutions for a given sequence.

        Returns:
            log_probs: Tensor of shape [seq_len, 20] representing log P(mut) - log P(wt).
        """
        # Get logits (either loaded from file or computed and saved)
        logits, sequence = self.get_or_compute_logits(kegg_id)

        # Ensure logits are on the correct device for calculation
        logits = logits.to(self.device)

        # Map sequence AA characters to ESM indices
        # Use .get(aa, -1) to handle unknown amino acids
        wt_indices_list = [self.aa_to_idx.get(aa, -1) for aa in sequence]
        wt_indices = torch.tensor(wt_indices_list, dtype=torch.long, device=self.device)

        # Filter out invalid AAs (-1)
        valid_mask = wt_indices >= 0
        valid_positions = valid_mask.nonzero(as_tuple=True)[0]

        if len(valid_positions) == 0:
            print(f"INFO: No valid amino acids found in {kegg_id}.")
            return torch.empty(0, 20, device=self.device)

        # Calculate Log Softmax of RAW logits -> gives log P(x)
        log_probs_all = log_softmax(logits[valid_positions], dim=-1)  # [valid_len, 20]

        # Extract log P(wt)
        wt_indices_valid = wt_indices[valid_positions]
        wt_log_probs = log_probs_all.gather(1, wt_indices_valid.unsqueeze(1))  # [valid_len, 1]

        # Score = log P(mut) - log P(wt)
        mutation_scores = log_probs_all - wt_log_probs

        return mutation_scores

    def save_mutation_scores_to_csv(self, kegg_id: str, csv_path: str):
        """
        Computes mutation scores for a single protein and adds them to an existing SNVs CSV file.
        The scores are added to a new column named 'esm_log_probs'.

        Args:
            csv_path: Path to the existing CSV file with variant information.
            :param kegg_id: kegg id of gene
        """
        # Compute scores for the given sequence
        try:
            score_matrix = self.score_all_mutations(kegg_id)
            # Move to CPU for easy indexing with Pandas
            score_matrix = score_matrix.cpu()
        except Exception as e:
            print(f"Failed to score mutations for {kegg_id}: {e}")
            return

        # Load the existing CSV file
        try:
            df = pd.read_csv(csv_path)
        except FileNotFoundError:
            raise f"Error: Input CSV file not found at {csv_path}"
        except Exception as e:
            raise f"Error reading CSV file: {e}"

        # Initialize a new 'esm_log_probs' column with nan
        if 'esm_log_probs' not in df.columns:
            df['esm_log_probs'] = np.nan

        # Iterate through the DataFrame and add scores
        for idx, row in df.iterrows():
            try:
                aa_pos = int(row['AA_index'])  # amino acid index (start in 0)

                variant = str(row['Variant'])  # e.g., M0L, M2R...
                mut_aa = variant[-1]  # Extract mutant amino acid

                if mut_aa == STOP_AA:  # stop codon, so score is stays np.nan
                    continue

                mut_idx = self.aa_to_idx.get(mut_aa)

                if mut_idx is not None and aa_pos < score_matrix.shape[0]:
                    score = score_matrix[aa_pos, mut_idx].item()
                    df.at[idx, 'esm_log_probs'] = score

            except Exception as e:
                print(f"[Warning] Could not assign score for row {idx} in {csv_path}: {e}")

        # Save the updated DataFrame to a new CSV file
        df.to_csv(csv_path, index=False)

    # def handle_cancer_row(self, row: pd.Series) -> pd.Series:
    #     """
    #     Takes a row with KeggId and RefrenceSeq, and produce all esm_log_probs
    #     """
    #     output_columns = ['esm_log_probs']
    #     results = pd.Series([np.nan] * len(output_columns), index=output_columns)
    #
    #     try:
    #         kegg_id = row['KeggId']
    #         mutation = row['Variant']
    #
    #         if pd.isna(kegg_id) or pd.isna(mutation):
    #             print("Missing KeggId or Variant")
    #             return results
    #
    #         # Parse the mutation string (e.g., 'p.V600E')
    #         # This regex captures the wild type AA, position, and mutant AA.
    #         match = re.match(r'(?:p\.)?([A-Z*])(\d+)([A-Z*])', str(mutation))
    #         if not match or not kegg_id or kegg_id == np.nan:
    #             print(f"No match to {mutation} or kegg_id: {kegg_id} problem.")
    #             return results  # Return NaNs if the format is not recognized
    #
    #         wt_aa, pos_str, mut_aa = match.groups()
    #         position = int(pos_str)  # 1-based position
    #
    #         if mut_aa == STOP_AA:
    #             print("The mutated AA is a stop codon, skipping...")
    #             return results
    #
    #         mut_aa_idx = self.aa_to_idx.get(mut_aa)
    #         if mut_aa_idx is None:
    #             return results
    #
    #         # Load and extract the ESM log probability
    #         esm_file = os.path.join(ESM_EMBEDDINGS_P, f"{kegg_id}.pt")
    #
    #         if not os.path.exists(esm_file):
    #             print(f"\nCould not find: {esm_file}")
    #             return results
    #
    #         esm_tensor = torch.load(esm_file, map_location='cpu')  # Expected format: [Length, 20]
    #
    #         idx_0 = position - 1
    #
    #         if idx_0 < esm_tensor.shape[0]:
    #             # We need to perform the log-softmax logic here on the raw logits
    #             logits_at_pos = esm_tensor[idx_0]  # [20]
    #             log_probs = log_softmax(logits_at_pos, dim=0)
    #
    #             wt_idx = self.aa_to_idx.get(wt_aa)
    #             if wt_idx is not None:
    #                 wt_log_prob = log_probs[wt_idx]
    #                 mut_log_prob = log_probs[mut_aa_idx]
    #                 score = mut_log_prob - wt_log_prob
    #                 results['esm_log_probs'] = score.item()
    #             else:
    #                 # Fallback if WT mismatch, just take mut log prob (unnormalized) or NaN
    #                 results['esm_log_probs'] = log_probs[mut_aa_idx].item()
    #
    #     except Exception as e:
    #         # In case of any error during processing, return the series of NaNs
    #         print(f"row exception: {e}")
    #
    #     return results
    #
    # def handle_cancer_csv(self, csv_path: str, recalc_scores=False):
    #     """
    #     EFFICIENT VERSION: Reads a CSV, pre-filters for mutations with available data,
    #     calculates scores, and saves the enriched DataFrame. Reports total time taken.
    #     """
    #     file_basename = os.path.basename(csv_path)
    #     print(f"--- Starting Cancer CSV Processing: {file_basename} ---")
    #     print(f"Loading mutation data from {file_basename}...")
    #     try:
    #         df = pd.read_csv(csv_path)
    #     except FileNotFoundError:
    #         print(f"[Error] Input file not found: {csv_path}")
    #         return
    #
    #     #### This part is for skipping files that are ready
    #     # Check if already processed
    #     if not recalc_scores and 'esm_log_probs' in df.columns:
    #         # Check if column is mostly full or completely empty
    #         if df['esm_log_probs'].notna().sum() > 0:
    #             print(f"Skipping {file_basename}, already scored.")
    #             return
    #
    #     if df.empty:
    #         return df
    #
    #     print(f"Scoring {len(df)} rows...")
    #     start_time = time.time()
    #
    #     scores_df = df.apply(self.handle_cancer_row, axis=1)
    #
    #     # Merge results
    #     df['esm_log_probs'] = scores_df['esm_log_probs']
    #
    #     df.to_csv(csv_path, index=False)
    #     print(f"Done. Time: {time.time() - start_time:.2f}s.\n")
    #     return df

    def handle_cancer_csv_optimized(self, csv_path: str, recalc_scores=False):
        """
        Highly optimized version of cancer CSV processing.
        Uses get_or_compute_logits to ensure data exists, then processes
        all associated mutations in memory.
        """
        file_basename = os.path.basename(csv_path)
        print(f"--- Starting Optimized Processing: {file_basename} ---")

        try:
            df = pd.read_csv(csv_path)
        except Exception as e:
            print(f"[Error] Could not read CSV: {e}")
            return

        # 1. Check if processing is needed
        if not recalc_scores and 'esm_log_probs' in df.columns:
            if df['esm_log_probs'].notna().sum() > 0:
                print(f"Skipping {file_basename}, already scored.")
                return df

        if df.empty:
            return df

        if 'esm_log_probs' not in df.columns:
            df['esm_log_probs'] = np.nan

        # 2. Pre-compile Regex
        mut_regex = re.compile(r'(?:p\.)?([A-Z*])(\d+)([A-Z*])')
        start_time = time.time()

        # 3. Group by KeggId
        grouped = df.groupby('KeggId')
        total_genes = len(grouped)

        print(f"Processing {len(df)} rows across {total_genes} unique genes...")

        for i, (kegg_id, group_indices) in enumerate(grouped.groups.items(), 1):
            if pd.isna(kegg_id) or str(kegg_id).lower() == 'nan':
                continue

            try:
                # REPLACED MANUAL LOADING WITH CLASS METHOD
                # This automatically loads if exists, or computes if missing.
                logits, _ = self.get_or_compute_logits(str(kegg_id))

                # Ensure logits are on the current device for math, then apply log_softmax
                logits = logits.to(self.device)
                log_probs_matrix = log_softmax(logits, dim=-1)

                for idx in group_indices:
                    variant_str = str(df.at[idx, 'Variant'])
                    match = mut_regex.match(variant_str)

                    if not match:
                        continue

                    wt_aa, pos_str, mut_aa = match.groups()
                    if mut_aa == STOP_AA:
                        continue

                    pos_idx = int(pos_str) - 1
                    mut_aa_idx = self.aa_to_idx.get(mut_aa)
                    wt_aa_idx = self.aa_to_idx.get(wt_aa)

                    if pos_idx < 0 or pos_idx >= log_probs_matrix.shape[0] or mut_aa_idx is None:
                        continue

                    lp_mut = log_probs_matrix[pos_idx, mut_aa_idx].item()
                    if wt_aa_idx is not None:
                        lp_wt = log_probs_matrix[pos_idx, wt_aa_idx].item()
                        df.at[idx, 'esm_log_probs'] = lp_mut - lp_wt
                    else:
                        df.at[idx, 'esm_log_probs'] = lp_mut

            except Exception as e:
                print(f"Error processing gene {kegg_id}: {e}")

            if i % 50 == 0:
                print(f"Progress: {i}/{total_genes} genes processed...")

        # 4. Save results
        df.to_csv(csv_path, index=False)
        print(f"Done! Processed {total_genes} genes in {time.time() - start_time:.2f}s.\n")
        return df