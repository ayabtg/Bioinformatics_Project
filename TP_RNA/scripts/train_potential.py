#!/usr/bin/env python3
import os
import argparse
import numpy as np
from Bio.PDB import PDBParser

# Base pairs to consider
BASE_PAIRS = [
    "AA","AU","AC","AG",
    "UU","UC","UG",
    "CC","CG",
    "GG"
]

def get_base_pair(a, b):
    """Return base pair type (unordered)."""
    pair = "".join(sorted([a, b]))
    return pair if pair in BASE_PAIRS else None

def extract_c3_atoms(structure):
    """Return list of tuples (residue index, base, C3' atom)."""
    atoms = []
    for model in structure:
        for chain in model:
            for res in chain:
                if "C3'" in res:
                    base = res.resname.strip()[0]   # A, U, C, G
                    idx  = res.get_id()[1]
                    atoms.append((idx, base, res["C3'"]))
        break  # Always use model 0
    return atoms

def compute_potentials(data_dir, output_dir):
    bins = np.arange(0, 21, 1)  # 0–20 Å, 1 Å bins
    counts = {bp: np.zeros(20) for bp in BASE_PAIRS}
    ref_counts = np.zeros(20)

    parser = PDBParser(QUIET=True)

    for pdb_file in os.listdir(data_dir):
        if not pdb_file.endswith(".pdb"):
            continue

        structure = parser.get_structure("rna", os.path.join(data_dir, pdb_file))
        atoms = extract_c3_atoms(structure)

        for i in range(len(atoms)):
            idx1, base1, atom1 = atoms[i]
            for j in range(i+4, len(atoms)):     # |i-j| ≥ 4
                idx2, base2, atom2 = atoms[j]

                d = atom1 - atom2
                dist = d

                if dist > 20:
                    continue

                bp = get_base_pair(base1, base2)
                bin_id = int(dist)

                ref_counts[bin_id] += 1
                if bp:
                    counts[bp][bin_id] += 1

    # Compute potentials
    for bp in BASE_PAIRS:
        # observed frequencies
        obs = counts[bp] / (np.sum(counts[bp]) + 1e-8)
        ref = ref_counts / (np.sum(ref_counts) + 1e-8)

        # potential formula
        potential = -np.log((obs + 1e-8) / (ref + 1e-8))
        potential = np.clip(potential, -10, 10)

        # save
        out = os.path.join(output_dir, f"{bp}.potential")
        np.savetxt(out, potential, fmt="%.5f")

    print("Potentials computed and saved.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    compute_potentials(args.data_dir, args.output_dir)
