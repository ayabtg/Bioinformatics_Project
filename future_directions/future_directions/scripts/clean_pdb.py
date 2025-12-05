#!/usr/bin/env python
import argparse
from pathlib import Path

def clean_pdb_file(in_path: Path, out_path: Path):
    """Simple PDB cleaner:
    - keep only ATOM/HETATM/TER/END lines
    - drop alternate locations not ' ' or 'A'
    """
    with in_path.open() as f_in, out_path.open("w") as f_out:
        for line in f_in:
            if not line.startswith(("ATOM", "HETATM", "TER", "END")):
                continue

            if line.startswith(("ATOM", "HETATM")) and len(line) >= 17:
                altloc = line[16]
                if altloc not in (" ", "A"):
                    continue

            f_out.write(line)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input_dir", type=Path, default=Path("pdb_train"),
        help="Folder with raw PDBs",
    )
    parser.add_argument(
        "--output_dir", type=Path, default=Path("pdb_train_clean"),
        help="Where to save cleaned PDBs",
    )
    args = parser.parse_args()

    args.output_dir.mkdir(exist_ok=True)

    for in_pdb in sorted(args.input_dir.glob("*.pdb")):
        out_pdb = args.output_dir / in_pdb.name
        clean_pdb_file(in_pdb, out_pdb)
        print(f"Cleaned {in_pdb.name} -> {out_pdb}")

if __name__ == "__main__":
    main()

