#!/usr/bin/env python3
import os
import argparse
import numpy as np
from Bio.PDB import PDBParser

# Base pairs considered in the potential files
BASE_PAIRS = ["AA","AU","AC","AG","UU","UC","UG","CC","CG","GG"]

def get_base_pair(a, b):
    """Return base-pair type (unordered)."""
    pair = "".join(sorted([a, b]))
    return pair if pair in BASE_PAIRS else None

def extract_atoms(structure, atom_name):
    """
    Extract atoms as (res_index, chain_id, base, atom_obj).
    Includes altLoc filtering (keeps altLoc ' ' and 'A').
    """
    atoms = []
    for model in structure:
        for chain in model:
            for res in chain:

                if atom_name not in res:
                    continue

                atom = res[atom_name]

                # Filter alternative atom positions
                if atom.get_altloc() not in (" ", "A"):
                    continue

                idx = res.get_id()[1]
                base = res.resname.strip()[0]
                chain_id = chain.id

                atoms.append((idx, chain_id, base, atom))
        break  # only first model
    return atoms

def linear_interpolate(value, table):
    """
    Interpolate score from potential table for a distance value.
    """
    bin0 = int(value)
    if bin0 >= len(table) - 1:
        return None

    bin1 = bin0 + 1
    alpha = value - bin0
    return (1 - alpha) * table[bin0] + alpha * table[bin1]

def score_structure(pdb_path, potentials, atom_name, cutoff, min_sep):
    """
    Compute score of an RNA structure using:
      ✔ Euclidean distances
      ✔ Linear interpolation
      ✔ Intrachain-only pairs
      ✔ Sequence separation ≥ min_sep
      ✔ Steric clash detection (< 2 Å)
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("rna", pdb_path)

    atoms = extract_atoms(structure, atom_name)

    # Initialize residue scores
    residue_scores = {idx: 0.0 for idx, _, _, _ in atoms}
    total_score = 0.0

    n = len(atoms)

    for i in range(n):
        idx1, chain1, b1, a1 = atoms[i]

        for j in range(i + min_sep, n):
            idx2, chain2, b2, a2 = atoms[j]

            # Require same chain (intrachain)
            if chain1 != chain2:
                continue

            bp = get_base_pair(b1, b2)
            if bp is None:
                continue

            dist = np.linalg.norm(a1.coord - a2.coord)

            # Steric clash
            if dist < 2.0:
                print(f"[Warning] Steric clash between residues {idx1} and {idx2}: {dist:.2f} Å")
                continue

            if dist > cutoff:
                continue

            score = linear_interpolate(dist, potentials[bp])
            if score is None:
                continue

            # Add score contribution
            residue_scores[idx1] += score
            residue_scores[idx2] += score
            total_score += score

    # Convert residue scores to sorted list for CSV export
    profile = sorted(residue_scores.items())
    return total_score, profile

def save_csv_profile(profile, total_score, out_csv):
    """Save profile + total score to CSV file."""
    with open(out_csv, "w") as f:
        f.write("position,score\n")
        for pos, val in profile:
            f.write(f"{pos},{val}\n")
        f.write(f"\nTOTAL_SCORE,{total_score}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Score an RNA structure using trained statistical potentials.")

    parser.add_argument("--input-pdb", required=True, help="Path to the PDB file.")
    parser.add_argument("--output-csv", required=True, help="Output CSV path.")
    parser.add_argument("--potentials-dir", required=True, help="Directory containing *.potential files.")

    parser.add_argument("--cutoff", type=float, default=20.0, help="Maximum interaction distance.")
    parser.add_argument("--min-sep", type=int, default=4, help="Minimum sequence separation.")
    parser.add_argument("--atom", type=str, default="C3'", help="Atom name to use for distance calculations.")

    args = parser.parse_args()

    # Load potential files
    potentials = {}
    for bp in BASE_PAIRS:
        fpath = os.path.join(args.potentials_dir, f"{bp}.potential")
        potentials[bp] = np.loadtxt(fpath)

    # Score structure
    total, profile = score_structure(
        args.input_pdb,
        potentials,
        args.atom,
        args.cutoff,
        args.min_sep
    )

    save_csv_profile(profile, total, args.output_csv)

    print(f"[OK] Scoring complete. Total score = {total:.3f}")
    print(f"[OK] Profile saved to {args.output_csv}")
