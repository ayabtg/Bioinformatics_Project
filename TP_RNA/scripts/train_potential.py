#!/usr/bin/env python3
"""train_potential.py
Template pour entraîner un potentiel (placeholder).
"""
import argparse


def main():
    parser = argparse.ArgumentParser(description="Entraîner le potentiel — template")
    parser.add_argument("--data-dir", default="../data/pdb_train", help="Répertoire des PDB d'entraînement")
    parser.add_argument("--output-dir", default="../output/potentials", help="Répertoire de sortie pour les potentiels")
    args = parser.parse_args()

    # TODO: implémenter l'entraînement réel
    print(f"Placeholder: entraîner sur {args.data_dir}, sauver dans {args.output_dir}")


if __name__ == '__main__':
    main()
