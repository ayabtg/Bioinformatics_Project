#!/usr/bin/env python3
"""score_structure.py
Template to score an RNA 3D structure using learned potentials.
Reads a PDB file and produces a profile of energy scores.
"""
import sys


def score_structure(pdb_path, out_path):
    # TODO: implement real scoring using potentials
    print(f"Scoring {pdb_path} -> {out_path}")
    with open(out_path, 'w') as f:
        f.write(f"# placeholder profile for {pdb_path}\n")


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print('Usage: score_structure.py <pdb_path> <out_profile>')
    else:
        score_structure(sys.argv[1], sys.argv[2])
