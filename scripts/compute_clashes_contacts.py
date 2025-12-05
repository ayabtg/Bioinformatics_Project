#!/usr/bin/env python3
"""
Compute steric clashes and residue-residue contacts
for cleaned RNA PDBs in data/pdb_clean/{train,test}.

Outputs:
  data/geometry_summary.csv with columns:
    pdb_id, split, n_chains, total_length,
    n_clashes, clashes_per_100_nt, max_clash_overlap,
    n_contacts, contacts_per_residue
"""

from pathlib import Path
import csv
from Bio.PDB import PDBParser, NeighborSearch

PROJECT_ROOT = Path(__file__).resolve().parents[1]

CLEAN_TRAIN = PROJECT_ROOT / "data" / "pdb_clean" / "train"
CLEAN_TEST = PROJECT_ROOT / "data" / "pdb_clean" / "test"
OUT_CSV = PROJECT_ROOT / "data" / "geometry_summary.csv"

RNA_RES = {"A", "C", "G", "U"}

# Very simple VDW radii in Å
VDW = {
    "C": 1.70,
    "N": 1.55,
    "O": 1.52,
    "P": 1.80,
    "S": 1.80,
}


def get_element(atom):
    """Robust element from Biopython Atom."""
    el = atom.element.strip()
    if el:
        return el
    # Fallback: first letter of atom name
    return atom.get_name().strip()[0]


def analyze_structure(structure):
    """Return (n_chains, total_len, clash_stats, contact_stats)."""
    # --- basic counts ---
    n_chains = 0
    total_len = 0

    rna_residues = []
    for model in structure:
        for chain in model:
            chain_res = [
                res for res in chain
                if res.get_resname().strip() in RNA_RES
            ]
            if not chain_res:
                continue
            n_chains += 1
            total_len += len(chain_res)
            rna_residues.extend(chain_res)

    # Build list of heavy atoms
    atoms = []
    for res in rna_residues:
        for atom in res:
            if get_element(atom) == "H":
                continue
            atoms.append(atom)

    if not atoms or total_len == 0:
        return n_chains, total_len, (0, 0.0, 0.0), (0, 0.0)

    ns = NeighborSearch(atoms)

    # --- clashes ---
    clash_pairs = set()
    max_overlap = 0.0

    # use a loose cutoff for candidate pairs
    for i, atom in enumerate(atoms):
        center = atom.coord
        neighbors = ns.search(center, 3.5)  # Å
        for nb in neighbors:
            if nb is atom:
                continue
            # ensure unique ordering
            if id(nb) < id(atom):
                continue

            res_i = atom.get_parent()
            res_j = nb.get_parent()
            # ignore same residue (mostly bonded)
            if res_i is res_j:
                continue

            el1 = get_element(atom)
            el2 = get_element(nb)
            r1 = VDW.get(el1, 1.7)
            r2 = VDW.get(el2, 1.7)

            diff = atom.coord - nb.coord
            dist = (diff * diff).sum() ** 0.5

            overlap = (r1 + r2) - dist
            if overlap > 0.4:  # clash threshold
                clash_pairs.add((res_i, res_j))
                if overlap > max_overlap:
                    max_overlap = overlap

    n_clashes = len(clash_pairs)
    clashes_per_100 = (n_clashes / total_len) * 100.0 if total_len else 0.0

    # --- contacts ---
    # residue-residue contact if any atom pair < 4.5 Å
    residue_pairs = set()
    for atom1, atom2 in ns.search_all(4.5):
        res1 = atom1.get_parent()
        res2 = atom2.get_parent()
        if res1 is res2:
            continue
        # unique unordered pair
        if id(res1) < id(res2):
            pair = (res1, res2)
        else:
            pair = (res2, res1)
        residue_pairs.add(pair)

    n_contacts = len(residue_pairs)
    contacts_per_res = (n_contacts / total_len) if total_len else 0.0

    clash_stats = (n_clashes, clashes_per_100, max_overlap)
    contact_stats = (n_contacts, contacts_per_res)
    return n_chains, total_len, clash_stats, contact_stats


def process_folder(folder, split, writer, parser):
    for pdb_file in sorted(folder.glob("*.pdb")):
        pdb_id = pdb_file.stem
        structure = parser.get_structure(pdb_id, pdb_file)

        n_chains, total_len, clash_stats, contact_stats = analyze_structure(structure)
        n_clashes, clashes_per_100, max_overlap = clash_stats
        n_contacts, contacts_per_res = contact_stats

        writer.writerow({
            "pdb_id": pdb_id,
            "split": split,
            "n_chains": n_chains,
            "total_length": total_len,
            "n_clashes": n_clashes,
            "clashes_per_100_nt": clashes_per_100,
            "max_clash_overlap": max_overlap,
            "n_contacts": n_contacts,
            "contacts_per_residue": contacts_per_res,
        })

        print(f"[{split}] {pdb_id}: len={total_len}, "
              f"clashes={n_clashes}, contacts={n_contacts}")


def main():
    parser = PDBParser(PERMISSIVE=True, QUIET=True)
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)

    with open(OUT_CSV, "w", newline="") as f:
        fieldnames = [
            "pdb_id", "split", "n_chains", "total_length",
            "n_clashes", "clashes_per_100_nt", "max_clash_overlap",
            "n_contacts", "contacts_per_residue",
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        process_folder(CLEAN_TRAIN, "train", writer, parser)
        process_folder(CLEAN_TEST, "test", writer, parser)


if __name__ == "__main__":
    main()





