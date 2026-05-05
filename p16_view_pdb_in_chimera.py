
import subprocess
import glob
import io
import os
import numpy as np
import pandas as pd
import requests
from Bio import PDB

from definitions import *

BRIDGES_DIR = os.path.join(BASE_P, 'results', 'p16_bottleneck_analysis')


def get_interface_residues_from_text(pdb_text: str, chain_a: str, chain_b: str,
                                     cutoff: float = 4.5) -> pd.DataFrame:
    """
    Returns a DataFrame of interface residues with one row per residue,
    including chain, residue number, residue name, and which atoms make contact.

    Output columns:
        chain        : chain ID (chain_a or chain_b)
        res_number   : PDB residue sequence number (used by ChimeraX/PyMOL)
        res_name     : three-letter amino acid code
        n_contacts   : number of atomic contacts with the partner chain
        contact_atoms: list of atom names making contact
    """
    parser = PDB.PDBParser(QUIET=True)
    struct = parser.get_structure("x", io.StringIO(pdb_text))
    model  = struct[0]

    def heavy_atoms(chain_id):
        """Returns list of (res_number, res_name, atom) for all heavy atoms."""
        return [
            (res.get_id()[1], res.get_resname(), atom)
            for res in model[chain_id]
            if res.get_id()[0] == " "  # exclude HETATMs (ligands, waters)
            for atom in res.get_atoms()
            if atom.element != "H"
        ]

    atoms_a = heavy_atoms(chain_a)
    atoms_b = heavy_atoms(chain_b)

    # Track contacts per residue
    contacts_a = {}  # {res_number: {"res_name": str, "atoms": set()}}
    contacts_b = {}

    for rid_a, rname_a, at_a in atoms_a:
        for rid_b, rname_b, at_b in atoms_b:
            if (at_a - at_b) <= cutoff:
                if rid_a not in contacts_a:
                    contacts_a[rid_a] = {"res_name": rname_a, "atoms": set()}
                if rid_b not in contacts_b:
                    contacts_b[rid_b] = {"res_name": rname_b, "atoms": set()}

                contacts_a[rid_a]["atoms"].add(at_a.get_name())
                contacts_b[rid_b]["atoms"].add(at_b.get_name())

    def to_df(contacts, chain_id):
        rows = []
        for res_num, info in sorted(contacts.items()):
            rows.append({
                "chain":         chain_id,
                "res_number":    res_num,
                "res_name":      info["res_name"],
                "n_contacts":    len(info["atoms"]),
                "contact_atoms": sorted(info["atoms"])
            })
        return pd.DataFrame(rows)

    return pd.concat([to_df(contacts_a, chain_a),
                      to_df(contacts_b, chain_b)], ignore_index=True)


def get_interface_from_binary(pdb_path: str, chain_a: str, chain_b: str,
                               cutoff: float = 5,
                               af3_json: str = None) -> pd.DataFrame:
    """
    Calls the interface binary and parses its output into a DataFrame.
    """
    cmd = ["/cs/staff/dina/utils/interface", pdb_path, chain_a, chain_b,
           str(cutoff), "-c"]  # -c for residue-residue contacts

    if af3_json:
        cmd += ["-j", af3_json]

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"Error: {result.stderr}")
        return pd.DataFrame()

    # Parse output: format is "res_number_a  chain_a - res_number_b  chain_b"
    # e.g. "199  A - 99  B"
    rows = []
    for line in result.stdout.strip().splitlines():
        parts = line.split()
        # data lines have exactly 5 tokens with '-' in the middle
        if len(parts) != 5 or parts[2] != "-":
            continue
        rows.append({
            "chain_a":      parts[1],
            "res_number_a": int(parts[0]),
            "chain_b":      parts[4],
            "res_number_b": int(parts[3]),
        })

    return pd.DataFrame(rows)


def _score_to_hex(score: float, light: tuple, dark: tuple) -> str:
    """Interpolate between light (score=0) and dark (score=1) RGB tuples."""
    r = int(light[0] + score * (dark[0] - light[0]))
    g = int(light[1] + score * (dark[1] - light[1]))
    b = int(light[2] + score * (dark[2] - light[2]))
    return f"#{r:02x}{g:02x}{b:02x}"

# Purple scale for chain_a mutations: lavender → deep purple
_PURPLE_LIGHT = (220, 180, 255)
_PURPLE_DARK  = (70,   0, 130)

# Orange scale for chain_b mutations: light peach → dark burnt orange
_ORANGE_LIGHT = (255, 220, 160)
_ORANGE_DARK  = (160,  50,   0)


def view_in_chimerax(pdb_id: str, chain_a: str, chain_b: str,
                     interface_a: list, interface_b: list,
                     mut_positions_a: dict = None, mut_positions_b: dict = None,
                     script_path: str = "view_interface.cxc",
                     file_path: str = None):
    """
    Generates a ChimeraX .cxc script.

    Visualization style:
      - Whole protein: molecular surface, 80% transparent (ghostly envelope)
      - Interface residues: ribbon, opaque, light grey (different shade per chain)
      - Mutations: surface recolored by score — purple scale for chain_a,
                   orange scale for chain_b (score 0=light, 1=dark)

    mut_positions_a / mut_positions_b: dict {residue_number: score (0-1)} or None
    If file_path is provided, opens that local file instead of fetching from RCSB.
    """
    resi_a = ",".join(map(str, interface_a))
    resi_b = ",".join(map(str, interface_b))

    open_cmd = f"open {file_path}" if file_path else f"open {pdb_id}"

    lines = [
        open_cmd,
        "hide all",

        # Slightly transparent surface for each chain — distinct light greys
        f"surface #1/{chain_a}",
        f"surface #1/{chain_b}",
        f"color #1/{chain_a} #dcdcdc target s",   # gainsboro (lighter)
        f"color #1/{chain_b} #a0a0a0 target s",   # darker grey
        f"transparency #1 25 surfaces",

        # Interface residues: opaque ribbon, matching grey per chain
        f"show #1/{chain_a}:{resi_a} cartoons",
        f"show #1/{chain_b}:{resi_b} cartoons",
        f"color #1/{chain_a}:{resi_a} #dcdcdc target r",
        f"color #1/{chain_b}:{resi_b} #a0a0a0 target r",
        f"transparency #1/{chain_a}:{resi_a} 0 cartoons",
        f"transparency #1/{chain_b}:{resi_b} 0 cartoons",

        # Focus on interface
        f"view #1/{chain_a}:{resi_a} | #1/{chain_b}:{resi_b}",

        "lighting soft",
        "set bgColor white",
    ]

    # Mutation coloring on the surface (no style change)
    if mut_positions_a:
        for res_num, score in mut_positions_a.items():
            hex_color = _score_to_hex(score, _PURPLE_LIGHT, _PURPLE_DARK)
            lines.append(f"color #1/{chain_a}:{res_num} {hex_color} target s")

    if mut_positions_b:
        for res_num, score in mut_positions_b.items():
            hex_color = _score_to_hex(score, _ORANGE_LIGHT, _ORANGE_DARK)
            lines.append(f"color #1/{chain_b}:{res_num} {hex_color} target s")

    script = "\n".join(lines)

    with open(script_path, "w") as f:
        f.write(script)
    print(f"        Script written to {script_path}")


import json
import glob
import pandas as pd

PAIR_LIST = "/cs/labs/dina/ophirmil12/PathwayAtlas2/data/af3_output/pair_list.txt"

def get_actifptm(pair_dir: str) -> float | None:
    """Reads actifptm from the summary_confidences.json in a pair's step3_out."""
    matches = glob.glob(f"{pair_dir}/step3_out/*/*_summary_confidences.json")
    if not matches:
        return None
    with open(matches[0]) as f:
        return json.load(f).get("actifptm")


def load_all_actifptms(pair_dirs: list[str]) -> pd.DataFrame:

    rows = []
    for pair_dir in pair_dirs:
        pair_name = pair_dir.rstrip("/").split("/")[-1]
        actifptm      = get_actifptm(pair_dir)
        if actifptm is None:
            tier = "missing"
        elif actifptm >= 0.75:
            tier = "high_confidence"
        elif actifptm >= 0.5:
            tier = "low_confidence"
        else:
            tier = "unreliable"
    

        rows.append({
            "pair_name": pair_name,
            "pair_dir":  pair_dir,
            "actifptm":      actifptm,
            "tier": tier,
        })

    return pd.DataFrame(rows)

def residue_score(cancer_mutations_df):
    probs  = cancer_mutations_df["pathogenic_prob"]
    counts = (cancer_mutations_df["count"] if "count" in cancer_mutations_df.columns
              else pd.Series(1, index=cancer_mutations_df.index))

    max_prob        = probs.max()
    weighted_burden = (probs * counts).sum()

    # Normalize burden by total patients to get a frequency
    total_patients  = counts.sum()
    burden_freq     = weighted_burden / total_patients

    # Geometric mean rewards positions that are BOTH
    # highly pathogenic AND recurrently mutated
    combined = np.sqrt(max_prob * burden_freq)

    return pd.Series({
        "max_pathogenic_prob": max_prob,
        "n_unique_mutations":  len(cancer_mutations_df),
        "total_patients":      total_patients,
        "burden_freq":         burden_freq,
        "residue_score":       combined,
    })


def compute_residue_scores(mutations_df: pd.DataFrame) -> pd.DataFrame:
    """
    From a merged cancer mutations DataFrame with columns [Protein, Variant, pathogenic_prob],
    extracts residue positions from Variant (e.g. 'V382M' -> 382), computes per-residue
    scores via residue_score(), and normalizes to [0, 1] globally.

    Returns DataFrame with columns [Protein, protein_position, residue_score, ...].
    """
    df = mutations_df.dropna(subset=["Variant", "pathogenic_prob"]).copy()
    df["protein_position"] = df["Variant"].str.extract(r'(\d+)')[0].astype(int)

    scores = (df
              .groupby(["Protein", "protein_position"], group_keys=False)
              .apply(residue_score)
              .reset_index())

    max_score = scores["residue_score"].max()
    if max_score > 0:
        scores["residue_score"] /= max_score

    return scores


def mutation_dicts_for_genes(residue_scores_df: pd.DataFrame,
                             gene_a: str, gene_b: str) -> tuple[dict, dict]:
    """
    Returns (dict_a, dict_b) mapping residue_number -> normalized_score (0-1)
    for each gene. Returns None for a gene with no scored mutations.
    """
    def gene_dict(gene):
        sub = residue_scores_df[residue_scores_df["Protein"] == gene]
        if sub.empty:
            print(f"    WARNING: Protein {gene} not found in data")
            return None, None
        return dict(zip(sub["protein_position"], sub["residue_score"])) 

    return gene_dict(gene_a), gene_dict(gene_b)

if __name__ == '__main__':

    # Load all cancer mutations once and compute per-residue scores
    print("Loading cancer mutations and computing residue scores...")
    all_mutations = pd.read_csv(pjoin(CBIO_CANCER_MUTATIONS_P, "pan_cancer.csv"), usecols=["Protein", "Variant", "pathogenic_prob"])
    residue_scores_df = compute_residue_scores(all_mutations)
    print(f"  Scored {len(residue_scores_df)} residues across "
          f"{residue_scores_df['Protein'].nunique()} genes")

    for f in sorted(glob.glob(os.path.join(BRIDGES_DIR, '*_bridges.csv'))):
        pathway_name = os.path.basename(f).split('_')[0]

        print(f"Locating known pdb structures from bridges in pathway {pathway_name}...")
        bridges_df = pd.read_csv(f)
        pdb_rows = bridges_df[bridges_df['strategy'] == 'pdb']

        print(f"    Found {len(pdb_rows)} pdb structures...")

        for _, bridge_row in pdb_rows.iterrows():
            pdb_id = bridge_row['pdb_id']
            gene_a = bridge_row['gene_a']
            gene_b = bridge_row['gene_b']
            try:
                pdb_text = requests.get(f"https://files.rcsb.org/download/{pdb_id}.pdb").text
                pdb_path = pjoin(BRIDGES_DIR, f"pdb_files/{pdb_id}.pdb")
                with open(pdb_path, "w") as pdb_file:
                    pdb_file.write(pdb_text)

                interface_path = pjoin(BRIDGES_DIR, f"interface_csvs/pdb_structures/{pathway_name}_{pdb_id}.csv")
                if os.path.exists(interface_path):
                    print(f"    Interface residues csv already exists for {pathway_name}, {pdb_id}")
                    interface = pd.read_csv(interface_path)
                else:
                    print(f"    Creating interface csv for {pdb_id}...")
                    interface = get_interface_from_binary(pdb_path, "A", "B")
                    if not interface.empty:
                        interface.to_csv(interface_path)

                if interface.empty:
                    print(f"    No interface found for {pdb_id}, skipping.")
                    continue

                iface_a = sorted(interface["res_number_a"].unique().tolist())
                iface_b = sorted(interface["res_number_b"].unique().tolist())

                mut_a, mut_b = mutation_dicts_for_genes(residue_scores_df, gene_a, gene_b)

                script_path = pjoin(BRIDGES_DIR, f"chimera_scripts/pdb_structures/{pathway_name}_{pdb_id}.cxc")
                if os.path.exists(script_path):
                    print(f"    Chimera script already exists for {pathway_name}, {pdb_id}")
                else:
                    print(f"    Creating Chimera script for complex {pdb_id} ({gene_a}/{gene_b})...")
                    view_in_chimerax(
                        pdb_id          = pdb_id,
                        chain_a         = "A",
                        chain_b         = "B",
                        interface_a     = iface_a,
                        interface_b     = iface_b,
                        mut_positions_a = mut_a,
                        mut_positions_b = mut_b,
                        script_path     = script_path,
                    )
            except Exception as e:
                print(f"    Caught exception: {e}")

    # AF3 pairs from pair_list.txt
    pair_list_path = pjoin(DATA_P, "af3_output/pair_list.txt")
    print(f"\nProcessing AF3 pairs from {pair_list_path}...")

    with open(pair_list_path) as fh:
        pair_dirs = [line.strip() for line in fh if line.strip()]

    # Load and filter
    df = load_all_actifptms(pair_dirs)

    print(f"Total pairs:           {len(df)}")
    print(f"Missing output:        {df['actifptm'].isna().sum()}")
    print(f"actifptm >= 0.75 (high):   {(df['actifptm'] >= 0.75).sum()}")
    print(f"actifptm 0.5-0.75 (low):   {df['actifptm'].between(0.5, 0.75).sum()}")
    print(f"actifptm < 0.5 (unreliable):{(df['actifptm'] < 0.5).sum()}")

    # Split into tiers
    high_confidence = df[df["actifptm"] >= 0.75].sort_values("actifptm", ascending=False)

    print("\nHigh confidence pairs:")
    print(high_confidence[["pair_name", "actifptm"]].to_string(index=False))

    for row in df.itertuples():
        pair_dir = row.pair_dir
        pair_name = row.pair_name
        tier = row.tier
        gene_a, gene_b = pair_name.split("_")
        cif_path = os.path.join(pair_dir, "step3_out", pair_name.lower(), f"{pair_name.lower()}_model.cif")

        if not os.path.exists(cif_path):
            print(f"    CIF not found, skipping: {cif_path}")
            continue

        print(f"    Processing AF3 pair {pair_name}...")
        try:
            interface_path = pjoin(BRIDGES_DIR, f"interface_csvs/{tier}/af3_{pair_name}.csv")
            if os.path.exists(interface_path):
                print(f"    Interface residues csv already exists for {pair_name}")
                interface = pd.read_csv(interface_path)
            else:
                print(f"    Creating interface csv for {pair_name}...")
                interface = get_interface_from_binary(cif_path, "A", "B")
                if not interface.empty:
                    interface.to_csv(interface_path)

            if interface.empty:
                print(f"    No interface found for {pair_name}, skipping.")
                continue

            iface_a = sorted(interface["res_number_a"].unique().tolist())
            iface_b = sorted(interface["res_number_b"].unique().tolist())

            mut_a, mut_b = mutation_dicts_for_genes(residue_scores_df, gene_a, gene_b)

            script_path = pjoin(BRIDGES_DIR, f"chimera_scripts/{tier}/af3_{pair_name}.cxc")
            if os.path.exists(script_path):
                print(f"    Chimera script already exists for {pair_name}")
            else:
                print(f"    Creating Chimera script for {pair_name}...")
                view_in_chimerax(
                    pdb_id          = pair_name,
                    chain_a         = "A",
                    chain_b         = "B",
                    interface_a     = iface_a,
                    interface_b     = iface_b,
                    mut_positions_a = mut_a,
                    mut_positions_b = mut_b,
                    script_path     = script_path,
                    file_path       = cif_path,
                )
        except Exception as e:
            print(f"    Caught exception: {e}")
