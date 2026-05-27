"""
Adds protein sequences to the p16 bridge CSVs, and generates AF3 input JSONs
for all bridge pairs that require AF3 (strategy == 'af3').

Sequence source:
  Primary   — cancer mutation files (Protein -> Sequence column, covers 496/500 genes)
  Fallback  — KEGG gene pickle aa_seq for the 4 genes absent from cancer mutations:
              APELA (hsa:100506013), BECN2 (hsa:441925),
              NFKBIE (hsa:4794), RAB7B (hsa:338382)

Outputs:
  - Updated *_bridges.csv files (seq_a, seq_b columns added)
  - data/af3_output/<GENE_A>_<GENE_B>/step1.json  for every AF3 pair
  - data/af3_output/pair_list.txt                  list of pair directories (for SLURM array)
"""

import os
import glob
import json
import pickle

import pandas as pd

from definitions import CBIO_CANCER_MUTATIONS_P, KEGG_GENES_P, BASE_P

BRIDGES_DIR = os.path.join(BASE_P, 'results', 'p16_bottleneck_analysis')
AF3_OUTPUT_P = os.path.join(BASE_P, 'data', 'af3_output')

# Genes absent from all cancer mutation files, identified by scanning pathway SNV files.
# KEGG IDs verified against KEGG gene pickle ref_names.
MISSING_GENE_KEGG_IDS = {
    'APELA':  'hsa:100506013',
    'BECN2':  'hsa:441925',
    'NFKBIE': 'hsa:4794',
    'RAB7B':  'hsa:338382',
}


def build_sequence_map():
    seq_map = {}
    for f in glob.glob(os.path.join(CBIO_CANCER_MUTATIONS_P, '*.csv')):
        df = pd.read_csv(f, usecols=['Protein', 'Sequence']).dropna(subset=['Sequence'])
        for protein, seq in df.drop_duplicates('Protein')[['Protein', 'Sequence']].values:
            if protein not in seq_map:
                seq_map[protein] = seq

    for gene, kegg_id in MISSING_GENE_KEGG_IDS.items():
        if gene not in seq_map:
            pickle_path = os.path.join(KEGG_GENES_P, f"hsa_{kegg_id.replace('hsa:', '')}.pickle")
            obj = pickle.load(open(pickle_path, 'rb'))
            seq_map[gene] = obj['aa_seq']
            print(f"  Resolved {gene} from KEGG pickle ({kegg_id}): {len(obj['aa_seq'])} aa")

    print(f"Sequence map: {len(seq_map)} genes")
    return seq_map


def add_sequences_to_bridges(seq_map):
    for f in sorted(glob.glob(os.path.join(BRIDGES_DIR, '*_bridges.csv'))):
        df = pd.read_csv(f)
        df['seq_a'] = df['gene_a'].map(seq_map)
        df['seq_b'] = df['gene_b'].map(seq_map)
        missing = (set(df.loc[df['seq_a'].isna(), 'gene_a']) |
                   set(df.loc[df['seq_b'].isna(), 'gene_b']))
        if missing:
            print(f"WARNING {os.path.basename(f)}: no sequence for {missing}")
        df.to_csv(f, index=False)
        print(f"Updated {os.path.basename(f)}  ({len(df)} bridges, "
              f"{df['seq_a'].notna().sum()} with seq_a, {df['seq_b'].notna().sum()} with seq_b)")


def create_af3_inputs(seq_map):
    os.makedirs(AF3_OUTPUT_P, exist_ok=True)
    pair_dirs = []

    for f in sorted(glob.glob(os.path.join(BRIDGES_DIR, '*_bridges.csv'))):
        df = pd.read_csv(f)
        af3_pairs = df[df['strategy'] == 'af3']

        for _, row in af3_pairs.iterrows():
            gene_a, gene_b = row['gene_a'], row['gene_b']
            seq_a = seq_map.get(gene_a)
            seq_b = seq_map.get(gene_b)
            if not seq_a or not seq_b:
                print(f"  Skipping {gene_a}-{gene_b}: missing sequence")
                continue

            pair_name = f"{gene_a}_{gene_b}"
            pair_dir = os.path.join(AF3_OUTPUT_P, pair_name)
            os.makedirs(pair_dir, exist_ok=True)

            step1 = {
                "name": pair_name,
                "modelSeeds": [1],
                "sequences": [
                    {"protein": {"id": "A", "sequence": seq_a}},
                    {"protein": {"id": "B", "sequence": seq_b}},
                ],
                "dialect": "alphafold3",
                "version": 1,
            }
            with open(os.path.join(pair_dir, 'step1.json'), 'w') as jf:
                json.dump(step1, jf, indent=2)
            pair_dirs.append(pair_dir)

    pair_list_path = os.path.join(AF3_OUTPUT_P, 'pair_list.txt')
    with open(pair_list_path, 'w') as pf:
        pf.write('\n'.join(pair_dirs))

    print(f"\nCreated {len(pair_dirs)} AF3 input JSONs")
    print(f"Pair list: {pair_list_path}")
    print(f"Submit with: sbatch --array=0-{len(pair_dirs) - 1} p16_run_msa_for_af3.slurm")


if __name__ == '__main__':
    seq_map = build_sequence_map()
    add_sequences_to_bridges(seq_map)
    create_af3_inputs(seq_map)
