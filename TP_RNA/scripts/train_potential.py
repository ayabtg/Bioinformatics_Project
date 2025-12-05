#!/usr/bin/env python3
import os
import argparse
import numpy as np
from Bio.PDB import PDBParser

# Base pairs to consider (unordered)
BASE_PAIRS = [
    "AA","AU","AC","AG",
    "UU","UC","UG",
    "CC","CG",
    "GG"
]

def get_base_pair(a, b):
    """Return unordered base pair type."""
    pair = "".join(sorted([a, b]))
    return pair if pair in BASE_PAIRS else None


def extract_atoms(structure, atom_name="C3'"):
    """
    Extract atoms as (res_index, chain_id, base, atom_object).
    Keeps only altLoc 'A' or ' '.
    """
    atoms = []
    for model in structure:
        for chain in model:
            for res in chain:

                if atom_name not in res:
                    continue

                atom = res[atom_name]

                # altLoc filtering
                if atom.get_altloc() not in (" ", "A"):
                    continue

                idx = res.get_id()[1]
                base = res.resname.strip()[0]   # A/U/C/G
                chain_id = chain.id

                atoms.append((idx, chain_id, base, atom))
        break  # model 0 only
    return atoms


def save_potentials_with_bin_centers(potentials, bin_edges, out_dir):
    """
    Save potentials as 2-column files:
    Distance(Å)   Score
    """
    os.makedirs(out_dir, exist_ok=True)

    # bin centers (e.g. 0.5, 1.5, 2.5…)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.0

    for bp, scores in potentials.items():
        out_path = os.path.join(out_dir, f"{bp}.potential")

        with open(out_path, "w") as f:
            f.write("# Distance(Å) Score\n")
            for d, s in zip(bin_centers, scores):
                f.write(f"{d:.3f} {s:.5f}\n")

        print(f"[INFO] Wrote {out_path}")


def compute_potentials(data_dir, output_dir, bin_width=1.0, cutoff=20.0, min_sep=4, atom_name="C3'"):
    """Compute knowledge-based potentials."""

    # Create distance bins 0–cutoff with step=bin_width
    bin_edges = np.arange(0, cutoff + bin_width, bin_width)
    n_bins = len(bin_edges) - 1

    # Initialize counters
    counts = {bp: np.zeros(n_bins) for bp in BASE_PAIRS}
    ref_counts = np.zeros(n_bins)

    parser = PDBParser(QUIET=True)

    for pdb_file in os.listdir(data_dir):
        if not pdb_file.endswith(".pdb"):
            continue

        structure = parser.get_structure("rna", os.path.join(data_dir, pdb_file))
        atoms = extract_atoms(structure, atom_name)

        for i in range(len(atoms)):
            idx1, chain1, base1, atom1 = atoms[i]

            for j in range(i + min_sep, len(atoms)):
                idx2, chain2, base2, atom2 = atoms[j]

                # Intrachain only
                if chain1 != chain2:
                    continue

                # Euclidean distance
                dist = np.linalg.norm(atom1.coord - atom2.coord)

                if dist <= 0 or dist > cutoff:
                    continue

                # assign to bin
                bin_id = np.searchsorted(bin_edges, dist) - 1
                if bin_id < 0 or bin_id >= n_bins:
                    continue

                bp = get_base_pair(base1, base2)

                ref_counts[bin_id] += 1
                if bp:
                    counts[bp][bin_id] += 1

    # Compute potentials
    potentials = {}
    for bp in BASE_PAIRS:
        obs = counts[bp] / (np.sum(counts[bp]) + 1e-8)
        ref = ref_counts / (np.sum(ref_counts) + 1e-8)

        pot = -np.log((obs + 1e-8) / (ref + 1e-8))
        pot = np.clip(pot, -10, 10)

        potentials[bp] = pot

    # Save using bin centers
    save_potentials_with_bin_centers(potentials, bin_edges, output_dir)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Train RNA knowledge-based potentials from PDB structures."
    )

    parser.add_argument("--data-dir", required=True, help="Folder containing training PDB files.")
    parser.add_argument("--output-dir", required=True, help="Folder where potentials will be saved.")
    parser.add_argument("--bin-width", type=float, default=1.0, help="Bin width in Å (default = 1.0).")
    parser.add_argument("--cutoff", type=float, default=20.0, help="Maximum distance in Å (default = 20).")
    parser.add_argument("--min-sep", type=int, default=4, help="Minimum |i-j| separation (default = 4).")
    parser.add_argument("--atom", type=str, default="C3'", help="Atom used for distance (default = C3').")

    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    compute_potentials(
        data_dir=args.data_dir,
        output_dir=args.output_dir,
        bin_width=args.bin_width,
        cutoff=args.cutoff,
        min_sep=args.min_sep,
        atom_name=args.atom,
    )
