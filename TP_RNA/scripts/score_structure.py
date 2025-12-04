#!/usr/bin/env python3
import os
import argparse
import numpy as np
from Bio.PDB import PDBParser

# -------------------------
# Base pairs to consider
# -------------------------
BASE_PAIRS = [
    "AA", "AU", "AC", "AG",
    "UU", "UC", "UG",
    "CC", "CG",
    "GG"
]

def get_base_pair(a, b):
    """Return base pair type (unordered) or None."""
    pair = "".join(sorted([a, b]))
    return pair if pair in BASE_PAIRS else None

def extract_c3_atoms(structure, atom_name):
    """
    Return list of tuples (residue index, base name, atom object).
    """
    atoms = []
    for model in structure:
        for chain in model:
            for res in chain:
                if atom_name in res:
                    base = res.resname.strip()[0]  # A, U, C, G
                    idx = res.get_id()[1]
                    atoms.append((idx, base, res[atom_name]))
            break  # use model 0 only
    return atoms

def load_potentials(directory):
    """
    Read potentials from .potential files.
    Files contain *only the potential values*, one per distance bin.
    Returns dict: { "AU": [values], "CG": [values], ... }
    """
    potentials = {}
    for fname in os.listdir(directory):
        if fname.endswith(".potential"):
            pair = fname.replace(".potential", "")
            values = []
            with open(os.path.join(directory, fname)) as f:
                for line in f:
                    line = line.strip()
                    if line:
                        values.append(float(line))
            potentials[pair] = values
    return potentials

def score_structure(pdb_path, potentials_dir, output_csv, cutoff, min_sep, atom):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("rna", pdb_path)

    # Extract atoms
    atoms = extract_c3_atoms(structure, atom)
    atoms.sort(key=lambda x: x[0])  # sort by residue index

    # Load potentials
    potentials = load_potentials(potentials_dir)

    nbins = len(next(iter(potentials.values())))
    bin_width = cutoff / nbins

    # Prepare output
    residues = [idx for idx, _, _ in atoms]
    scores = [0.0] * len(residues)

    # ----------------------
    # Score each residue
    # ----------------------
    for i in range(len(atoms)):
        idx_i, base_i, atom_i = atoms[i]

        total_score = 0.0

        for j in range(len(atoms)):
            if abs(i - j) < min_sep:
                continue

            idx_j, base_j, atom_j = atoms[j]

            pair = get_base_pair(base_i, base_j)
            if pair is None or pair not in potentials:
                continue

            # Compute distance
            distance = atom_i - atom_j
            dist = np.linalg.norm(distance)

            if dist > cutoff:
                continue

            # Find bin index
            bin_index = int(dist / bin_width)
            if bin_index >= nbins:
                bin_index = nbins - 1

            # Retrieve the score
            vals = potentials[pair]
            total_score += vals[bin_index]

        scores[i] = total_score

    # Write output CSV
    os.makedirs(os.path.dirname(output_csv), exist_ok=True)
    with open(output_csv, "w") as f:
        f.write("residue,score\n")
        for idx, score in zip(residues, scores):
            f.write(f"{idx},{score:.4f}\n")

    print(f"Scores saved to: {output_csv}")


# -------------------------------------------------
# CLI parsing (PARAMETRIZATION)
# -------------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Score RNA structure using learned potentials.")
    parser.add_argument("--input-pdb", required=True, help="Input PDB file.")
    parser.add_argument("--potentials-dir", required=True, help="Directory with *.potential files.")
    parser.add_argument("--output-csv", required=True, help="Where the per-position scores will be saved.")

    # Parameterized options (your contribution)
    parser.add_argument("--cutoff", type=float, default=20.0, help="Maximum distance to consider (Å).")
    parser.add_argument("--min-sep", type=int, default=4, help="Minimum |i-j| separation.")
    parser.add_argument("--atom", type=str, default="C3'", help="Atom used for scoring distances.")

    args = parser.parse_args()

    score_structure(
        pdb_path=args.input_pdb,
        potentials_dir=args.potentials_dir,
        output_csv=args.output_csv,
        cutoff=args.cutoff,
        min_sep=args.min_sep,
        atom=args.atom
    )
