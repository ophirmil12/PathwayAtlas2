
import subprocess
import glob
import io
import os
import json
import numpy as np
import pandas as pd
import requests
from Bio import PDB
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

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

# Ribbon colors: lightest tier
_RIBBON_A = "#A9D0E8"   # light blue
_RIBBON_B = "#F0B8BF"   # light pink

# Interface surface patch colors: mid tier
_BLUE_LIGHT = (169, 208, 232)
_BLUE_DARK  = ( 92, 138, 178)

_PINK_LIGHT = (240, 184, 191)
_PINK_DARK  = (216, 111, 122)

# Mutated interface residues: darkest tier
_BLUE_MUT  = "#1A3A5C"
_PINK_MUT  = "#7B1010"


def view_in_chimerax(pdb_id: str, chain_a: str, chain_b: str,
                     interface_a: list, interface_b: list,
                     mut_positions_a: dict = None, mut_positions_b: dict = None,
                     non_interface_transparency: int = 50,
                     script_path: str = "view_interface.cxc",
                     file_path: str = None):
    """
    Generates a ChimeraX .cxc script.

    Visualization style:
            - Full ribbon for both chains: light blue (A) / light pink (B)
      - Non-interface residues: ribbon only, no extra marking
            - Interface residues: molecular surface patch in a darker blue/pink tone
            - Mutated interface residues: darkest blue/pink tone regardless of score

    mut_positions_a / mut_positions_b: dict {residue_number: score (0-1)} or None
    non_interface_transparency: transparency percent for the cartoons (0 = opaque, 100 = invisible)
    If file_path is provided, opens that local file instead of fetching from RCSB.
    """
    resi_a = ",".join(map(str, interface_a))
    resi_b = ",".join(map(str, interface_b))
    open_cmd = f"open {os.path.basename(file_path)}" if file_path else f"open {pdb_id}"

    _default_a = _score_to_hex(0.35, _BLUE_LIGHT, _BLUE_DARK)
    _default_b = _score_to_hex(0.35, _PINK_LIGHT, _PINK_DARK)

    lines = [
        open_cmd,
        "hide all",

        # Full ribbon for both proteins
        f"show #1/{chain_a} cartoons",
        f"show #1/{chain_b} cartoons",
        f"color #1/{chain_a} {_RIBBON_A} target r",
        f"color #1/{chain_b} {_RIBBON_B} target r",
        f"transparency #1/{chain_a} {non_interface_transparency} target r",
        f"transparency #1/{chain_b} {non_interface_transparency} target r",

        # Interface: molecular surface patch (score=0 default, overridden per residue below)
        f"surface #1/{chain_a}:{resi_a}",
        f"surface #1/{chain_b}:{resi_b}",
        f"color #1/{chain_a}:{resi_a} {_default_a} target s",
        f"color #1/{chain_b}:{resi_b} {_default_b} target s",

        # Focus on interface
        f"view #1/{chain_a}:{resi_a} | #1/{chain_b}:{resi_b}",

        "lighting soft",
        "set bgColor white",
    ]

    # Per-residue score coloring on the surface patch (only for scored interface residues)
    iface_set_a = set(interface_a)
    iface_set_b = set(interface_b)

    if mut_positions_a:
        for res_num, score in mut_positions_a.items():
            if res_num in iface_set_a:
                lines.append(f"color #1/{chain_a}:{res_num} {_BLUE_MUT} target s")

    if mut_positions_b:
        for res_num, score in mut_positions_b.items():
            if res_num in iface_set_b:
                lines.append(f"color #1/{chain_b}:{res_num} {_PINK_MUT} target s")

    # Mutation coloring on ribbons (all mutations, both interface and non-interface)
    if mut_positions_a:
        for res_num, score in mut_positions_a.items():
            lines.append(f"color #1/{chain_a}:{res_num} {_BLUE_MUT} target r")

    if mut_positions_b:
        for res_num, score in mut_positions_b.items():
            lines.append(f"color #1/{chain_b}:{res_num} {_PINK_MUT} target r")

    with open(script_path, "w") as f:
        f.write("\n".join(lines))
    print(f"        Script written to {script_path}")


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

def run_pair_stats(gene_a: str, gene_b: str,
                   interface_df: pd.DataFrame,
                   residue_scores_df: pd.DataFrame) -> dict | None:
    """
    Mann-Whitney U test (one-sided, greater): do interface residues have
    higher max pathogenicity scores than non-interface residues?

    Uses max_pathogenic_prob per residue (already computed in residue_scores_df).
    Combines both proteins in the pair. Returns None if either group has < 2 points.
    """
    iface_set_a = set(interface_df["res_number_a"].unique())
    iface_set_b = set(interface_df["res_number_b"].unique())

    iface_scores, non_iface_scores = [], []

    for gene, iface_set in [(gene_a, iface_set_a), (gene_b, iface_set_b)]:
        sub = residue_scores_df[residue_scores_df["Protein"] == gene]
        for _, r in sub.iterrows():
            bucket = iface_scores if r["protein_position"] in iface_set else non_iface_scores
            bucket.append(r["max_pathogenic_prob"])

    if len(iface_scores) < 2 or len(non_iface_scores) < 2:
        return None

    stat, p = mannwhitneyu(iface_scores, non_iface_scores, alternative='greater')
    return {
        "n_interface":             len(iface_scores),
        "n_non_interface":         len(non_iface_scores),
        "mean_interface_score":    float(np.mean(iface_scores)),
        "mean_non_interface_score":float(np.mean(non_iface_scores)),
        "U_statistic":             float(stat),
        "p_value":                 float(p),
        "q_value":                 None,  # filled after FDR
    }


def _apply_fdr(stats_df: pd.DataFrame, mask: pd.Series) -> pd.DataFrame:
    valid = mask & stats_df["p_value"].notna()
    if valid.sum() > 1:
        _, q_vals, _, _ = multipletests(stats_df.loc[valid, "p_value"], method='fdr_bh')
        stats_df.loc[valid, "q_value"] = q_vals
    return stats_df


if __name__ == '__main__':
    os.makedirs(BRIDGES_DIR, exist_ok=True)
    os.makedirs(pjoin(BRIDGES_DIR, "interface_csvs"), exist_ok=True)
    os.makedirs(pjoin(BRIDGES_DIR, "chimera_scripts"), exist_ok=True)

    # Create output directories
    for sub in ["pdb_files",
                "interface_csvs/pdb_structures", "chimera_scripts/pdb_structures",
                "interface_csvs/high_confidence", "chimera_scripts/high_confidence"]:
        os.makedirs(pjoin(BRIDGES_DIR, sub), exist_ok=True)

    # Load all cancer mutations once and compute per-residue scores
    print("Loading cancer mutations and computing residue scores...")
    all_mutations = pd.read_csv(pjoin(CBIO_CANCER_MUTATIONS_P, "pan_cancer.csv"),
                                usecols=["Protein", "Variant", "pathogenic_prob"])
    residue_scores_df = compute_residue_scores(all_mutations)
    print(f"  Scored {len(residue_scores_df)} residues across "
          f"{residue_scores_df['Protein'].nunique()} genes")

    # Build gene-pair lookup from bridges CSVs (avoids fragile string split on pair_name)
    pair_to_genes = {}
    for f in glob.glob(os.path.join(BRIDGES_DIR, '*_bridges.csv')):
        pathway_name = os.path.basename(f).split('_')[0]
        for _, row in pd.read_csv(f).iterrows():
            pair_to_genes[f"{row['gene_a']}_{row['gene_b']}"] = (row['gene_a'], row['gene_b'], pathway_name)

    stats_rows = []

    # ── PDB pairs ─────────────────────────────────────────────────────────────
    for f in sorted(glob.glob(os.path.join(BRIDGES_DIR, '*_bridges.csv'))):
        pathway_name = os.path.basename(f).split('_')[0]
        bridges_df   = pd.read_csv(f)
        pdb_rows     = bridges_df[bridges_df['strategy'] == 'pdb']
        print(f"\n[PDB] Pathway {pathway_name}: {len(pdb_rows)} pairs")

        for _, bridge_row in pdb_rows.iterrows():
            pdb_id = bridge_row['pdb_id']
            gene_a = bridge_row['gene_a']
            gene_b = bridge_row['gene_b']
            print(f"  {gene_a}–{gene_b} ({pdb_id})")
            try:
                pdb_text = requests.get(f"https://files.rcsb.org/download/{pdb_id}.pdb").text
                pdb_path = pjoin(BRIDGES_DIR, f"pdb_files/{pdb_id}.pdb")
                with open(pdb_path, "w") as pf:
                    pf.write(pdb_text)

                interface = get_interface_from_binary(pdb_path, "A", "B")
                pair_name = f"{gene_a}_{gene_b}"
                if not interface.empty:
                    interface.to_csv(pjoin(BRIDGES_DIR, f"interface_csvs/pdb_structures/pdb_{pair_name}.csv"), index=False)
                else:
                    print(f"    No interface found, skipping.")
                    continue

                iface_a  = sorted(interface["res_number_a"].unique().tolist())
                iface_b  = sorted(interface["res_number_b"].unique().tolist())
                mut_a, mut_b = mutation_dicts_for_genes(residue_scores_df, gene_a, gene_b)

                cxc_dir = pjoin(BRIDGES_DIR, "chimera_scripts/pdb_structures")
                link_name = f"{gene_a}_{gene_b}_{pdb_id}.pdb"
                link = pjoin(cxc_dir, link_name)
                if not os.path.exists(link):
                    os.symlink(os.path.abspath(pdb_path), link)

                view_in_chimerax(
                    pdb_id          = pdb_id,
                    chain_a         = "A", chain_b = "B",
                    interface_a     = iface_a, interface_b = iface_b,
                    mut_positions_a = mut_a,  mut_positions_b = mut_b,
                    script_path     = pjoin(cxc_dir, f"pdb_{pair_name}.cxc"),
                    file_path       = link,
                )

                result = run_pair_stats(gene_a, gene_b, interface, residue_scores_df)
                stats_rows.append({
                    "pair_name": pair_name, "gene_a": gene_a, "gene_b": gene_b,
                    "strategy": "pdb", "pathway": pathway_name, "actifptm": None,
                    **(result or {"n_interface": None, "n_non_interface": None,
                                  "mean_interface_score": None, "mean_non_interface_score": None,
                                  "U_statistic": None, "p_value": None, "q_value": None}),
                })

            except Exception as e:
                print(f"    Exception: {e}")

    # ── AF3 pairs ──────────────────────────────────────────────────────────────
    pair_list_path = pjoin(DATA_P, "af3_output/pair_list.txt")
    print(f"\n[AF3] Loading pairs from {pair_list_path}...")
    with open(pair_list_path) as fh:
        pair_dirs = [line.strip() for line in fh if line.strip()]

    af3_df = load_all_actifptms(pair_dirs)
    print(f"  Total: {len(af3_df)}  |  high: {(af3_df['actifptm'] >= 0.75).sum()}  "
          f"|  low: {af3_df['actifptm'].between(0.5, 0.75).sum()}  "
          f"|  unreliable: {(af3_df['actifptm'] < 0.5).sum()}  "
          f"|  missing: {af3_df['actifptm'].isna().sum()}")

    # Process only high_confidence AF3 pairs (CXC + stats)
    for row in af3_df.itertuples():
        if row.tier != "high_confidence":
            continue

        pair_dir  = row.pair_dir
        pair_name = row.pair_name
        gene_a, gene_b, pathway_name = pair_to_genes.get(pair_name, (None, None, None))
        if gene_a is None:
            print(f"  {pair_name}: not found in bridges CSV, skipping.")
            continue

        cif_path = pjoin(pair_dir, "step3_out", pair_name.lower(), f"{pair_name.lower()}_model.cif")
        if not os.path.exists(cif_path):
            print(f"  {pair_name}: CIF not found, skipping.")
            continue

        print(f"  {pair_name}")
        try:
            interface = get_interface_from_binary(cif_path, "A", "B")
            if interface.empty:
                print(f"    No interface found, skipping.")
                continue
            interface.to_csv(pjoin(BRIDGES_DIR, f"interface_csvs/high_confidence/af3_{pair_name}.csv"), index=False)

            iface_a  = sorted(interface["res_number_a"].unique().tolist())
            iface_b  = sorted(interface["res_number_b"].unique().tolist())
            mut_a, mut_b = mutation_dicts_for_genes(residue_scores_df, gene_a, gene_b)

            cxc_dir = pjoin(BRIDGES_DIR, "chimera_scripts/high_confidence")
            view_in_chimerax(
                pdb_id          = pair_name,
                chain_a         = "A", chain_b = "B",
                interface_a     = iface_a, interface_b = iface_b,
                mut_positions_a = mut_a,  mut_positions_b = mut_b,
                script_path     = pjoin(cxc_dir, f"af3_{pair_name}.cxc"),
                file_path       = cif_path,
            )
            link = pjoin(cxc_dir, os.path.basename(cif_path))
            if not os.path.exists(link):
                os.symlink(os.path.abspath(cif_path), link)

            result = run_pair_stats(gene_a, gene_b, interface, residue_scores_df)
            stats_rows.append({
                "pair_name": pair_name, "gene_a": gene_a, "gene_b": gene_b,
                "strategy": "af3", "pathway": pathway_name, "actifptm": row.actifptm,
                **(result or {"n_interface": None, "n_non_interface": None,
                              "mean_interface_score": None, "mean_non_interface_score": None,
                              "U_statistic": None, "p_value": None, "q_value": None}),
            })

        except Exception as e:
            print(f"    Exception: {e}")

    # ── FDR correction (separately for PDB and AF3) ───────────────────────────
    stats_df = pd.DataFrame(stats_rows)
    for strategy in ["pdb", "af3"]:
        stats_df = _apply_fdr(stats_df, stats_df["strategy"] == strategy)

    out_path = pjoin(BRIDGES_DIR, "p16_interface_pathogenicity_stats.csv")
    stats_df.to_csv(out_path, index=False)
    print(f"\nSaved stats for {len(stats_df)} pairs to {out_path}")
    sig = stats_df[stats_df["q_value"].notna() & (stats_df["q_value"] < 0.05)]
    print(f"Significant pairs (q < 0.05): {len(sig)} / {stats_df['p_value'].notna().sum()}")
