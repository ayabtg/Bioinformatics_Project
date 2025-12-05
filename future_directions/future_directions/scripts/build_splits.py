#!/usr/bin/env python
import argparse
from pathlib import Path

def write_list(pdb_dir: Path, out_file: Path, split_name: str):
    with out_file.open("w") as f:
        for pdb in sorted(pdb_dir.glob("*.pdb")):
            pdb_id = pdb.stem  # e.g. 1A3M
            f.write(pdb_id + "\n")
    print(f"Wrote {split_name} list to {out_file}")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--train_dir", type=Path, default=Path("pdb_train_clean"))
    ap.add_argument("--test_dir", type=Path, default=Path("pdb_test"))
    ap.add_argument("--out_dir", type=Path, default=Path("results"))
    args = ap.parse_args()

    args.out_dir.mkdir(exist_ok=True)

    write_list(args.train_dir, args.out_dir / "train_ids.txt", "train")
    write_list(args.test_dir, args.out_dir / "test_ids.txt", "test")

if __name__ == "__main__":
    main()

