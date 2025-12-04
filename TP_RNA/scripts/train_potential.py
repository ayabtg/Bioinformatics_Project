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

def extract_atoms(structure, atom_name="C3'"):
    """Return list of tuples (residue index, base, chosen atom)."""
    atoms = []
    for model in structure:
        for chain in model:
            for res in chain:
                if atom_name in res:
                    base = res.resname.strip()[0]   # A, U, C, G
                    idx  = res.get_id()[1]
                    atoms.append((idx, base, res[atom_name]))
        break  # Always use model 0
    return atoms

def compute_potentials(data_dir, output_dir, n_bins=20, cutoff=20.0, min_sep=4, atom_name="C3'"):
    """Compute knowledge-based potentials from a folder of PDB files."""
    # Distance bins: from 0 to cutoff, divided into n_bins intervals
    bin_edges = np.linspace(0.0, cutoff, n_bins + 1)
    counts = {bp: np.zeros(n_bins) for bp in BASE_PAIRS}
    ref_counts = np.zeros(n_bins)

    parser = PDBParser(QUIET=True)

    for pdb_file in os.listdir(data_dir):
        if not pdb_file.endswith(".pdb"):
            continue

        structure = parser.get_structure("rna", os.path.join(data_dir, pdb_file))
        atoms = extract_atoms(structure, atom_name)

        for i in range(len(atoms)):
            idx1, base1, atom1 = atoms[i]
            for j in range(i + min_sep, len(atoms)):     # |i-j| ≥ min_sep
                idx2, base2, atom2 = atoms[j]

                # Bio.PDB overloads '-' to return the distance between 2 atoms (float, in Å)
                dist = atom1 - atom2

                if dist < 0 or dist > cutoff:
                    continue

                # Map distance to bin index
                bin_id = np.searchsorted(bin_edges, dist, side="right") - 1
                if bin_id < 0 or bin_id >= n_bins:
                    continue

                bp = get_base_pair(base1, base2)

                ref_counts[bin_id] += 1
                if bp:
                    counts[bp][bin_id] += 1

    # Compute potentials
    for bp in BASE_PAIRS:
        # observed and reference frequencies
        obs = counts[bp] / (np.sum(counts[bp]) + 1e-8)
        ref = ref_counts / (np.sum(ref_counts) + 1e-8)

        # potential formula (capped between -10 and 10)
        potential = -np.log((obs + 1e-8) / (ref + 1e-8))
        potential = np.clip(potential, -10, 10)

        # save
        out = os.path.join(output_dir, f"{bp}.potential")
        np.savetxt(out, potential, fmt="%.5f")

    print("Potentials computed and saved.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Train RNA knowledge-based potentials from PDB files."
    )

    parser.add_argument(
        "--data-dir",
        required=True,
        help="Folder containing training PDB files.",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Folder where potentials will be saved.",
    )
    parser.add_argument(
        "--bins",
        type=int,
        default=20,
        help="Number of distance bins (default: 20).",
    )
    parser.add_argument(
        "--cutoff",
        type=float,
        default=20.0,
        help="Maximum distance in Å to consider (default: 20).",
    )
    parser.add_argument(
        "--min-sep",
        type=int,
        default=4,
        help="Minimum |i-j| sequence separation (default: 4).",
    )
    parser.add_argument(
        "--atom",
        type=str,
        default="C3'",
        help="Atom name used for distances (default: C3').",
    )

    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    compute_potentials(
        args.data_dir,
        args.output_dir,
        n_bins=args.bins,
        cutoff=args.cutoff,
        min_sep=args.min_sep,
        atom_name=args.atom,
    )
