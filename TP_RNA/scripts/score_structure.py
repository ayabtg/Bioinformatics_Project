#!/usr/bin/env python3
import argparse
import os
from Bio import PDB

import sys, os
sys.path.append(os.path.dirname(__file__))


from Train_potential import (
    extract_c3_atoms,
    calculate_ED,
    BASE_PAIRS,
    SORTED_PAIR_TO_CANONICAL
)

from plot_potential import load_all_potentials




def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Score RNA structures using distance-based statistical potentials."
    )

    parser.add_argument(
        "input",
        nargs="+",
        help="PDB/mmCIF files OR a single directory containing structures."
    )

    parser.add_argument(
        "--atom-type", "-a",
        type=str,
        default="C3'",
        help="Atom type for distance calculations (default: C3')."
    )

    parser.add_argument(
        "--csv-output", "-o",
        type=str,
        default=None,
        help="Save results into a CSV file."
    )

    parser.add_argument(
        "--potentials-dir", "-p",
        type=str,
        default="data/potentials",
        help="Directory containing *.txt potential files."
    )

    return parser.parse_args()




def get_raw_distances(structure, atom_type="C3'"):
    distances = {bp: [] for bp in BASE_PAIRS}

    for model in structure:
        for chain in model:
            atoms = extract_c3_atoms(chain, atom_type)
            n = len(atoms)

            for i in range(n):
                idx1, a1, b1, _ = atoms[i]
                for j in range(i + 1, n):
                    idx2, a2, b2, _ = atoms[j]

                    # Skip very close sequence positions (<4)
                    if idx2 - idx1 < 4:
                        continue

                    dist = calculate_ED(a1, a2)

                    pair = ''.join(sorted((b1, b2)))
                    pair = SORTED_PAIR_TO_CANONICAL.get(pair, None)

                    if pair is not None:
                        distances[pair].append(dist)

    return distances




def score_interpolated(profile, dist):
    bins, scores = profile

    if dist <= bins[0]:
        return scores[0]
    if dist >= bins[-1]:
        return scores[-1]

    bw = bins[1] - bins[0]
    idx = int((dist - bins[0]) // bw)

    fraction = (dist - bins[idx]) / bw
    return scores[idx] + fraction * (scores[idx + 1] - scores[idx])





def compute_total_score(distances, potentials):
    total = 0.0

    for bp, dlist in distances.items():
        profile = potentials[bp]

        for d in dlist:
            total += score_interpolated(profile, d)

    return total





def main():
    args = parse_args()

    potentials = load_all_potentials(args.potentials_dir)

    PDBparser = PDB.PDBParser(QUIET=True)
    CIFparser = PDB.MMCIFParser(QUIET=True)

    all_scores = {}

    # Case 1 → user gave a directory
    if os.path.isdir(args.input[0]):
        folder = args.input[0]
        files = [
            os.path.join(folder, f)
            for f in os.listdir(folder)
            if f.lower().endswith(".pdb") or f.lower().endswith(".cif")
        ]
    else:
        files = args.input

    for path in files:
        print(f"[INFO] Processing {path}")

        parser = CIFparser if path.endswith(".cif") else PDBparser
        structure = parser.get_structure("RNA", path)

        distances = get_raw_distances(structure, args.atom_type)
        score = compute_total_score(distances, potentials)

        print(f"  → Score: {score}")
        all_scores[path] = score

    
    if args.csv_output:
        with open(args.csv_output, "w") as f:
            for k, v in all_scores.items():
                f.write(f"{k},{v}\n")
        print(f"[INFO] Saved CSV to: {args.csv_output}")


if __name__ == "__main__":
    main()
