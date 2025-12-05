#!/usr/bin/env python3
"""
Score RNA structures using trained distance-based potentials (.potential format).
Fully compatible with the updated training and plotting scripts.
"""

import os
import sys
import argparse
import numpy as np
from Bio import PDB


# ----------------------------------------------------------------------
# FIX: allow importing modules from TP_RNA/scripts
# ----------------------------------------------------------------------
PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.append(PROJECT_ROOT)

from scripts.train_potential import BASE_PAIRS, SORTED_PAIR_TO_CANONICAL
from scripts.plot_potential import load_all_potentials   # now reads .potential files


# ----------------------------------------------------------------------
# ARGUMENTS
# ----------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Score RNA structures using trained distance-based potentials."
    )

    parser.add_argument(
        "input",
        nargs="+",
        help="PDB/mmCIF file(s) or a directory containing PDB/mmCIF files."
    )

    parser.add_argument(
        "--atom-type", "-a",
        type=str, default="C3'",
        help="Atom used for distances (default: C3')."
    )

    parser.add_argument(
        "--potentials-dir", "-p",
        type=str, default="output/potentials",
        help="Directory containing <BP>.potential files."
    )

    parser.add_argument(
        "--csv-output", "-o",
        type=str, default=None,
        help="Save scores to CSV."
    )

    return parser.parse_args()


# ----------------------------------------------------------------------
# EXTRACT ATOMS
# ----------------------------------------------------------------------
def extract_atoms_from_chain(chain, atom_type="C3'"):
    atoms = []
    seq_idx = 0

    for residue in chain:
        name = residue.get_resname().upper()
        if not name or name[0] not in {"A", "U", "C", "G"}:
            continue

        # find C3' or equivalent
        atom = None
        for a in residue:
            if a.id.upper() in {atom_type.upper(), atom_type.replace("'", "*")}:
                atom = a
                break
        if atom is None:
            continue

        atoms.append((seq_idx, name[0], np.array(atom.coord, float)))
        seq_idx += 1

    return atoms


# ----------------------------------------------------------------------
# COMPUTE ALL DISTANCES FOR EACH BASE PAIR
# ----------------------------------------------------------------------
def get_distances(structure, atom_type="C3'"):
    distances = {bp: [] for bp in BASE_PAIRS}

    for model in structure:
        for chain in model:
            atoms = extract_atoms_from_chain(chain, atom_type)
            n = len(atoms)

            for i in range(n):
                _, base1, pos1 = atoms[i]

                # minimum sequence separation = 4 (same logic as training)
                for j in range(i + 4, n):
                    _, base2, pos2 = atoms[j]

                    d = float(np.linalg.norm(pos1 - pos2))
                    if d <= 0:
                        continue

                    bp = "".join(sorted([base1, base2]))

                    if bp in distances:
                        distances[bp].append(d)

        break  # single model only

    return distances


# ----------------------------------------------------------------------
# LINEAR INTERPOLATION
# ----------------------------------------------------------------------
def interpolate(scoreprofile, dist):
    bins, scores = scoreprofile

    if dist <= bins[0]:
        return scores[0]
    if dist >= bins[-1]:
        return scores[-1]

    step = bins[1] - bins[0]
    idx = int((dist - bins[0]) / step)

    if idx >= len(scores) - 1:
        return scores[-1]

    d0 = bins[idx]
    d1 = bins[idx + 1]
    t = (dist - d0) / (d1 - d0)

    return scores[idx] * (1 - t) + scores[idx + 1] * t


# ----------------------------------------------------------------------
# SCORE STRUCTURE
# ----------------------------------------------------------------------
def score_structure(distances, potentials):
    total = 0.0

    for bp in BASE_PAIRS:
        if bp not in distances or bp not in potentials:
            continue

        dist_list = distances[bp]
        profile = potentials[bp]  # (bins, scores)

        for d in dist_list:
            total += interpolate(profile, d)

    return total


# ----------------------------------------------------------------------
# MAIN
# ----------------------------------------------------------------------
def main():
    args = parse_args()

    # Load potentials (distance + score)
    potentials = load_all_potentials(args.potentials_dir)

    PDBparser = PDB.PDBParser(QUIET=True)
    CIFparser = PDB.MMCIFParser(QUIET=True)

    results = {}

    # If input is a directory → score all PDB files inside
    if os.path.isdir(args.input[0]):
        files = [
            os.path.join(args.input[0], f)
            for f in os.listdir(args.input[0])
            if f.endswith((".pdb", ".cif"))
        ]
    else:
        files = args.input

    for filepath in files:
        print(f"[INFO] Scoring {filepath}")

        if filepath.endswith(".pdb"):
            parser = PDBparser
        elif filepath.endswith(".cif"):
            parser = CIFparser
        else:
            print(f"[WARNING] Unsupported file format: {filepath}")
            continue

        structure = parser.get_structure("RNA", filepath)

        # Extract raw distances for each base pair
        dists = get_distances(structure, args.atom_type)

        # Score using potentials
        score = score_structure(dists, potentials)
        results[os.path.basename(filepath)] = score

        print(f" → Score = {score}")

    # Save CSV if requested
    if args.csv_output:
        with open(args.csv_output, "w") as f:
            for name, val in results.items():
                f.write(f"{name},{val}\n")
        print(f"[INFO] Scores written to {args.csv_output}")


if __name__ == "__main__":
    main()
