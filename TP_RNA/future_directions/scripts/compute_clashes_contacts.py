#!/usr/bin/env python
import argparse
from pathlib import Path
import math
import csv

def parse_pdb_atoms(pdb_path: Path):
    atoms = []
    residues = set()
    chains = set()

    with pdb_path.open() as f:
        for line in f:
            if not line.startswith(("ATOM", "HETATM")):
                continue

            chain_id = line[21].strip() or "?"
            resseq = line[22:26].strip()
            try:
                resseq_int = int(resseq)
            except ValueError:
                resseq_int = None

            # coordinates
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue

            atom_name = line[12:16]
            # rough check to skip hydrogens
            if atom_name.strip().startswith("H"):
                continue

            chains.add(chain_id)
            if resseq_int is not None:
                residues.add((chain_id, resseq_int))
            atoms.append((x, y, z))

    return atoms, len(residues), len(chains)

def count_contacts_and_clashes(atoms, contact_cut=4.5, clash_cut=2.0):
    contact_cut2 = contact_cut ** 2
    clash_cut2 = clash_cut ** 2
    n = len(atoms)
    contacts = 0
    clashes = 0

    for i in range(n):
        x1, y1, z1 = atoms[i]
        for j in range(i + 1, n):
            x2, y2, z2 = atoms[j]
            dx = x1 - x2
            dy = y1 - y2
            dz = z1 - z2
            d2 = dx*dx + dy*dy + dz*dz

            if d2 < clash_cut2:
                clashes += 1
            elif d2 < contact_cut2:
                contacts += 1

    return contacts, clashes

def process_dir(pdb_dir: Path, split_name: str, rows: list):
    for pdb_path in sorted(pdb_dir.glob("*.pdb")):
        pdb_id = pdb_path.stem
        atoms, length, n_chains = parse_pdb_atoms(pdb_path)
        contacts, clashes = count_contacts_and_clashes(atoms)
        rows.append({
            "pdb_id": pdb_id,
            "split": split_name,
            "length": length,
            "chains": n_chains,
            "contacts": contacts,
            "clashes": clashes,
        })
        print(f"Processed {pdb_id}: len={length}, chains={n_chains}, "
              f"contacts={contacts}, clashes={clashes}")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--train_dir", type=Path, default=Path("pdb_train_clean"))
    ap.add_argument("--test_dir", type=Path, default=Path("pdb_test"))
    ap.add_argument("--out_csv", type=Path, default=Path("results/pdb_stats.csv"))
    args = ap.parse_args()

    args.out_csv.parent.mkdir(exist_ok=True)

    rows = []
    process_dir(args.train_dir, "train", rows)
    process_dir(args.test_dir, "test", rows)

    with args.out_csv.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["pdb_id", "split", "length", "chains", "contacts", "clashes"],
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nSaved dataset stats to {args.out_csv}")

if __name__ == "__main__":
    main()

