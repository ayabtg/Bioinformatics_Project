#!/usr/bin/env python3
"""
Create train_ids.txt and test_ids.txt from pdb_train / pdb_test folders.
"""

from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
TRAIN_DIR = PROJECT_ROOT / "data" / "pdb_train"
TEST_DIR = PROJECT_ROOT / "data" / "pdb_test"
OUT_DIR = PROJECT_ROOT / "data"

def get_ids(folder):
    pdb_ids = []
    for f in sorted(folder.iterdir()):
        if f.suffix.lower() == ".pdb":
            pdb_ids.append(f.stem)  # remove .pdb
    return pdb_ids

def write_list(ids, path):
    with open(path, "w") as f:
        for pid in ids:
            f.write(pid + "\n")

def main():
    train_ids = get_ids(TRAIN_DIR)
    test_ids = get_ids(TEST_DIR)

    OUT_DIR.mkdir(exist_ok=True, parents=True)
    write_list(train_ids, OUT_DIR / "train_ids.txt")
    write_list(test_ids, OUT_DIR / "test_ids.txt")

    print("Train IDs:", train_ids)
    print("Test IDs:", test_ids)

if __name__ == "__main__":
    main()

