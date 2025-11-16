#!/usr/bin/env python3
"""score_structure.py
Template pour scorer une structure et produire un fichier de profil (placeholder).
"""
import sys


def score_structure(pdb_path, out_path):
    # TODO: implémenter le calcul de score réel
    print(f"Scoring {pdb_path} -> {out_path}")
    with open(out_path, 'w') as f:
        f.write(f"# profil placeholder pour {pdb_path}\n")


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print('Usage: score_structure.py <pdb_path> <out_profile>')
    else:
        score_structure(sys.argv[1], sys.argv[2])
