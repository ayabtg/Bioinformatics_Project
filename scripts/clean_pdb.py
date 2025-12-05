#!/usr/bin/env python3
"""
Clean and validate RNA PDB dataset.

- Reads PDBs from data/pdb_train and data/pdb_test
- Keeps only RNA residues (A, C, G, U)
- Removes water and non-RNA residues
- Discards altLoc atoms except ' ' or 'A'
- Checks for missing backbone atoms
- Writes cleaned PDBs to data/pdb_clean/train and data/pdb_clean/test
- Writes a CSV summary to data/clean_summary.csv
"""

import os
import csv
from pathlib import Path

from Bio.PDB import PDBParser, PDBIO, Select

# --- configuration ---
PROJECT_ROOT = Path(__file__).resolve().parents[1]
TRAIN_DIR = PROJECT_ROOT / "data" / "pdb_train"
TEST_DIR = PROJECT_ROOT / "data" / "pdb_test"
OUT_TRAIN_DIR = PROJECT_ROOT / "data" / "pdb_clean" / "train"
OUT_TEST_DIR = PROJECT_ROOT / "data" / "pdb_clean" / "test"
SUMMARY_CSV = PROJECT_ROOT / "data" / "clean_summary.csv"

RNA_RES = {"A", "C", "G", "U"}
BACKBONE_ATOMS = {"P", "O5'", "C5'", "C4'", "C3'", "O3'"}


class RNASelect(Select):
    """Select only RNA residues and good altLoc atoms."""
    def accept_residue(self, residue):
        return residue.get_resname().strip() in RNA_RES

    def accept_atom(self, atom):
        # skip hydrogens
        if atom.element == "H":
            return 0
        # keep only main altLoc (blank or 'A')
        if atom.get_altloc() not in (" ", "A"):
            return 0
        return 1


def analyze_structure(structure):
    """
    Return (n_chains, total_len, missing_backbone_atoms)
    considering only RNA residues.
    """
    n_chains = 0
    total_len = 0
    missing_backbone = 0

    for model in structure:
        for chain in model:
            # collect RNA residues
            rna_residues = [
                res for res in chain
                if res.get_resname().strip() in RNA_RES
            ]
            if not rna_residues:
                continue

            n_chains += 1
            total_len += len(rna_residues)

            for res in rna_residues:
                present = {atom.get_name() for atom in res}
                if not BACKBONE_ATOMS.issubset(present):
                    missing_backbone += 1

    return n_chains, total_len, missing_backbone


def process_folder(in_dir, out_dir, split_name, writer, parser):
    out_dir.mkdir(parents=True, exist_ok=True)

    for fname in sorted(os.listdir(in_dir)):
        if not fname.lower().endswith(".pdb"):
            continue

        pdb_id = fname.split(".")[0]
        pdb_path = in_dir / fname

        structure = parser.get_structure(pdb_id, pdb_path)
        n_chains, total_len, missing_bb = analyze_structure(structure)

        # simple validity rules (you can tweak thresholds)
        valid = True
        reason = ""
        if n_chains == 0:
            valid = False
            reason = "no_rna_chain"
        elif total_len < 10:
            valid = False
            reason = "too_short"
        elif missing_bb > 0:
            valid = False
            reason = "missing_backbone_atoms"

        # write cleaned PDB even if invalid (so you can inspect)
        io = PDBIO()
        io.set_structure(structure)
        io.save(str(out_dir / fname), select=RNASelect())

        writer.writerow({
            "pdb_id": pdb_id,
            "split": split_name,
            "n_chains": n_chains,
            "total_length": total_len,
            "n_residues_missing_backbone": missing_bb,
            "valid": int(valid),
            "reason": reason,
        })

        print(f"[{split_name}] {pdb_id}: chains={n_chains}, len={total_len}, "
              f"missing_bb={missing_bb}, valid={valid}")


def main():
    parser = PDBParser(PERMISSIVE=True, QUIET=True)

    SUMMARY_CSV.parent.mkdir(parents=True, exist_ok=True)
    with open(SUMMARY_CSV, "w", newline="") as f:
        fieldnames = [
            "pdb_id", "split", "n_chains", "total_length",
            "n_residues_missing_backbone", "valid", "reason"
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        process_folder(TRAIN_DIR, OUT_TRAIN_DIR, "train", writer, parser)
        process_folder(TEST_DIR, OUT_TEST_DIR, "test", writer, parser)


if __name__ == "__main__":
    main()

