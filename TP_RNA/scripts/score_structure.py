#!/usr/bin/env python3
import os
import argparse
import numpy as np
from Bio.PDB import PDBParser

# Base pairs considered in the potential files
BASE_PAIRS = ["AA","AU","AC","AG","UU","UC","UG","CC","CG","GG"]

def get_base_pair(a, b):
    """Return base pair type (unordered)."""
    pair = "".join(sorted([a, b]))
    return pair if pair in BASE_PAIRS else None

def extract_atoms(structure):
    """
    Extract atoms as (res_index, chain_id, base, C3' atom).
    Includes altLoc filtering: keep only altLoc 'A' or ' '.
    """
    atoms = []
    for model in structure:       # Only first model is used
        for chain in model:
            for res in chain:

                if "C3'" not in res:
                    continue

                atom = res["C3'"]

                # Filter undesired altLoc (alternative coordinates)
                if atom.get_altloc() not in (" ", "A"):
                    continue

                idx = res.get_id()[1]
                base = res.resname.strip()[0]
                chain_id = chain.id

                atoms.append((idx, chain_id, base, atom))
        break
    return atoms

def score_structure(pdb_path, potentials, max_distance=20.0):
    """
    Compute:
      ✔ Euclidean distance
      ✔ Linear interpolation
      ✔ Intrachain-only distances
      ✔ Ignore steric clashes (<2 Å)
      ✔ Total score + residue-level profile
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("rna", pdb_path)

    atoms = extract_atoms(structure)

    # Initialize residue scores
    scores = {idx: 0.0 for idx, _, _, _ in atoms}
    total_score = 0.0

    n = len(atoms)

    for i in range(n):
        idx1, chain1, b1, a1 = atoms[i]

        for j in range(i + 4, n):   # sequence separation >= 4
            idx2, chain2, b2, a2 = atoms[j]

            # Intrachain filter
            if chain1 != chain2:
                continue

            bp = get_base_pair(b1, b2)
            if bp is None:
                continue

            # Euclidean distance
            dist = np.linalg.norm(a1.coord - a2.coord)

            # Steric clash
            if dist < 2.0:
                print(f"Warning: steric clash between residues {idx1} and {idx2}: {dist:.2f} Å")
                continue

            if dist > max_distance:
                continue

            # Linear interpolation
            bin0 = int(dist)
            if bin0 >= len(potentials[bp]) - 1:
                continue  # outside potential table range

            bin1 = bin0 + 1
            alpha = dist - bin0

            interp_score = (1 - alpha) * potentials[bp][bin0] + alpha * potentials[bp][bin1]

            # Add to residue-level scores
            scores[idx1] += interp_score
            scores[idx2] += interp_score

            # Add to global score
            total_score += interp_score

    # Sorted residue profile
    profile = sorted(scores.items())
    return total_score, profile

def save_profile(profile, total_score, out_csv):
    """Save residue profile + total score."""
    with open(out_csv, "w") as f:
        f.write("position,score\n")
        for pos, sc in profile:
            f.write(f"{pos},{sc}\n")
        f.write(f"\nTOTAL_SCORE,{total_score}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Score an RNA structure using trained statistical potentials.")
    parser.add_argument("pdb_file")
    parser.add_argument("output_csv")
    parser.add_argument("--pot-dir", default="output/potentials")
    parser.add_argument("--max-distance", type=float, default=20.0)
    args = parser.parse_args()

    # Load potential files
    potentials = {}
    for bp in BASE_PAIRS:
        path = os.path.join(args.pot_dir, f"{bp}.potential")
        potentials[bp] = np.loadtxt(path)

    # Score the structure
    total_score, profile = score_structure(args.pdb_file, potentials, args.max_distance)

    # Save results
    save_profile(profile, total_score, args.output_csv)

    print(f"Scoring complete. Total score = {total_score:.3f}")
