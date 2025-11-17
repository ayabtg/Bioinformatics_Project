#!/usr/bin/env python3
import os
import argparse
import numpy as np
from Bio.PDB import PDBParser

BASE_PAIRS = ["AA","AU","AC","AG","UU","UC","UG","CC","CG","GG"]

def get_base_pair(a, b):
    pair = "".join(sorted([a, b]))
    return pair if pair in BASE_PAIRS else None

def load_potentials(pot_dir):
    potentials = {}
    for bp in BASE_PAIRS:
        path = os.path.join(pot_dir, f"{bp}.potential")
        potentials[bp] = np.loadtxt(path)
    return potentials

def extract_atoms(structure):
    atoms = []
    for model in structure:
        for chain in model:
            for res in chain:
                if "C3'" in res:
                    atoms.append((res.get_id()[1], res.resname[0], res["C3'"]))
        break
    return atoms

def score_structure(pdb_path, potentials):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("rna", pdb_path)

    atoms = extract_atoms(structure)
    scores = {idx: 0.0 for (idx, _, _) in atoms}

    for i in range(len(atoms)):
        idx1, b1, a1 = atoms[i]
        for j in range(i+4, len(atoms)):
            idx2, b2, a2 = atoms[j]
            bp = get_base_pair(b1, b2)
            if not bp:
                continue

            d = a1 - a2
            dist = int(min(max(d, 0), 19))  # bin 0–19

            scores[idx1] += potentials[bp][dist]
            scores[idx2] += potentials[bp][dist]

    # Convert to sorted list
    return sorted(scores.items())

def save_profile(result, out_csv):
    with open(out_csv, "w") as f:
        f.write("position,score\n")
        for pos, sc in result:
            f.write(f"{pos},{sc}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("pdb_file")
    parser.add_argument("output_csv")
    parser.add_argument("--pot-dir", default="output/potentials")
    args = parser.parse_args()

    pot = load_potentials(args.pot_dir)
    result = score_structure(args.pdb_file, pot)
    save_profile(result, args.output_csv)

    print("Profile saved.")
