#!/usr/bin/env python3
import os
import sys
import argparse
from math import log
from typing import List

import numpy as np
from Bio import PDB



BASE_PAIRS = ['AA', 'AU', 'AC', 'AG', 'UU', 'UC', 'UG', 'CC', 'CG', 'GG']

RESIDUE_NAME_MAP = {
    'A': 'A', 'ADE': 'A', 'DA': 'A',
    'C': 'C', 'CYT': 'C', 'DC': 'C',
    'G': 'G', 'GUA': 'G', 'DG': 'G',
    'U': 'U', 'URA': 'U', 'DU': 'U'
}

SORTED_PAIR_TO_CANONICAL = {
    'AA': 'AA',
    'AC': 'AC',
    'AG': 'AG',
    'AU': 'AU',
    'CC': 'CC',
    'CG': 'CG',
    'CU': 'UC',
    'GG': 'GG',
    'GU': 'UG',
    'UU': 'UU',
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Train RNA statistical potential.")

    parser.add_argument(
        "inputs",
        nargs="+",
        help="PDB/MMCIF files or directories containing PDB/MMCIF files."
    )

    parser.add_argument(
        "--out-dir", "-o",
        default="potentials",
        help="Output directory (default: potentials)."
    )

    parser.add_argument("--round-decimals", "-r", type=int, default=3)
    parser.add_argument("--min-distance", "-mind", type=float, default=2.0)

    parser.add_argument(
        "--atom-type", "-a",
        type=str,
        default="C3'",
        help="Atom type to use (default: C3')."
    )

    parser.add_argument(
        "--density-method", "-dm",
        type=str,
        choices=["histogram", "kde"],
        default="histogram"
    )

    parser.add_argument("--kde-bandwidth", "-kbw", type=float, default=0.5)
    parser.add_argument("--max-distance", "-maxd", type=float, default=20.0)
    parser.add_argument("--min-seq-sep", "-s", type=int, default=4)
    parser.add_argument("--bin-width", "-w", type=float, default=1.0)
    parser.add_argument("--max-score", "-ms", type=float, default=10.0)

    return parser.parse_args()

def collect_pdb_paths(inputs: List[str]) -> List[str]:
    """Collect all PDB / CIF / ENT files from the given paths."""
    pdb_paths = []
    for path in inputs:
        if os.path.isdir(path):
            for name in os.listdir(path):
                full = os.path.join(path, name)
                if (
                    os.path.isfile(full)
                    and (
                        name.lower().endswith(".pdb")
                        or name.lower().endswith(".cif")
                        or name.lower().endswith(".ent")
                    )
                ):
                    pdb_paths.append(full)
        elif os.path.isfile(path):
            pdb_paths.append(path)
        else:
            print(f"[WARN] Invalid path: {path}", file=sys.stderr)
    return sorted(set(pdb_paths))


def calculate_ED(a1, a2) -> float:
    """Euclidean distance between two atoms."""
    return float(np.linalg.norm(a1.coord - a2.coord))


def extract_c3_atoms(chain, atom_type: str = "C3'"):
    """
    Extract list of (seq_index, atom, base, bfactor) for residues
    that contain the requested atom.
    """
    atoms_list = []
    residues = list(chain)

    for i, residue in enumerate(residues):
        name = residue.resname.strip()
        base = RESIDUE_NAME_MAP.get(name, "")

        if not base:
            continue

        try:
            atom = residue[atom_type]
        except KeyError:
            continue

        bfactor = atom.get_bfactor()
        atoms_list.append((i, atom, base, bfactor))

    return atoms_list


def canonical_pair(r1: str, r2: str) -> str:
    """Return canonical base-pair name (unordered, mapped to our set)."""
    pair = ''.join(sorted((r1, r2)))
    return SORTED_PAIR_TO_CANONICAL.get(pair, "")


def update_distance_counts(
    c3_atoms,
    observed_counts,
    reference_counts,
    distance_bins,
    min_distance,
    max_distance,
    min_seq_sep
):
    n = len(c3_atoms)
    n_bins = len(distance_bins) - 1

    for i in range(n):
        seq1, atom1, res1, bf1 = c3_atoms[i]
        for j in range(i + 1, n):
            seq2, atom2, res2, bf2 = c3_atoms[j]

            # Skip pairs too close along the sequence
            if seq2 - seq1 < min_seq_sep:
                continue

            dist = calculate_ED(atom1, atom2)

            if dist < min_distance:
                continue
            if dist >= max_distance:
                continue

            bin_id = int(np.searchsorted(distance_bins, dist) - 1)
            if not (0 <= bin_id < n_bins):
                continue

            pair = canonical_pair(res1, res2)
            if pair in observed_counts:
                observed_counts[pair][bin_id] += 1

            reference_counts[bin_id] += 1

def compute_frequency(counts, pseudo: float = 1e-12):
    """Convert raw counts to frequencies with a small pseudocount."""
    total = counts.sum()
    if total <= 0:
        return None
    freq = counts / total
    return np.clip(freq, pseudo, None)


def compute_single_score(f_obs: float, f_ref: float, max_score: float) -> float:
    """Single-bin potential score using -log(f_obs/f_ref)."""
    if f_obs > 0 and f_ref > 0:
        return min(-log(f_obs / f_ref), max_score)
    return max_score


def compute_scores(
    observed_counts,
    reference_counts,
    n_bins,
    round_decimals,
    distance_bins,
    min_distance,
    max_score
):
    """Compute statistical potential scores for each base pair and distance bin."""
    scores = {}
    centers = (distance_bins[:-1] + distance_bins[1:]) / 2.0

    f_ref = compute_frequency(reference_counts)
    if f_ref is None:
        f_ref = np.ones(n_bins) / n_bins

    for bp, counts in observed_counts.items():
        f_obs = compute_frequency(counts)
        if f_obs is None:
            scores[bp] = [max_score] * n_bins
            continue

        bp_scores = []
        for r in range(n_bins):
            if centers[r] < min_distance:
                bp_scores.append(round(max_score, round_decimals))
            elif reference_counts[r] == 0:
                bp_scores.append(round(max_score, round_decimals))
            else:
                val = compute_single_score(f_obs[r], f_ref[r], max_score)
                bp_scores.append(round(val, round_decimals))

        scores[bp] = bp_scores

    return scores


def save_scores(scores, distance_bins, out_dir: str):
    """
    Save one file per base pair in the format:

    # Distance(Å)  Score
    0.5  10.0
    1.5  9.3
    ...
    """
    os.makedirs(out_dir, exist_ok=True)
    centers = (distance_bins[:-1] + distance_bins[1:]) / 2.0

    for bp, arr in scores.items():
        path = os.path.join(out_dir, f"{bp}.txt")

        with open(path, "w") as f:
            f.write("# Distance(Å)  Score\n")
            for d, s in zip(centers, arr):
                f.write(f"{d:.1f}  {s}\n")

        print(f"[INFO] Wrote {path}")

# Main execution

def main():
    args = parse_args()

    pdb_paths = collect_pdb_paths(args.inputs)
    if not pdb_paths:
        print("[ERROR] No PDB files found.")
        sys.exit(1)

    # Define distance bins
    num_bins = int(args.max_distance / args.bin_width) + 1
    distance_bins = np.linspace(0, args.max_distance, num_bins)
    n_bins = len(distance_bins) - 1

    # Initialize counts
    observed_counts = {bp: np.zeros(n_bins, dtype=int) for bp in BASE_PAIRS}
    reference_counts = np.zeros(n_bins, dtype=int)

    PDBparser = PDB.PDBParser(QUIET=True)
    CIFparser = PDB.MMCIFParser(QUIET=True)

    # Loop over all structures and accumulate distances
    for pdb_path in pdb_paths:
        if pdb_path.lower().endswith(".cif"):
            parser = CIFparser
        else:
            parser = PDBparser

        structure = parser.get_structure("RNA", pdb_path)

        for model in structure:
            for chain in model:
                atoms = extract_c3_atoms(chain, args.atom_type)
                if atoms:
                    update_distance_counts(
                        atoms,
                        observed_counts,
                        reference_counts,
                        distance_bins,
                        args.min_distance,
                        args.max_distance,
                        args.min_seq_sep,
                    )

    # Compute potentials
    scores = compute_scores(
        observed_counts,
        reference_counts,
        n_bins,
        args.round_decimals,
        distance_bins,
        args.min_distance,
        args.max_score,
    )

    save_scores(scores, distance_bins, args.out_dir)


if __name__ == "__main__":

    main()
