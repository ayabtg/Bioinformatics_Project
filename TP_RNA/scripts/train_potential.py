#!/usr/bin/env python3
"""train_potential.py
Template to train a statistical potential from RNA 3D structures.
Reads PDB files, extracts base-pair distances, bins them, and learns energy profiles.
"""
import argparse


def main():
    parser = argparse.ArgumentParser(description="Train the statistical potential — template")
    parser.add_argument("--data-dir", default="../data/pdb_train", help="Directory with training PDB files")
    parser.add_argument("--output-dir", default="../output/potentials", help="Output directory for potential files")
    args = parser.parse_args()

    # TODO: implement real potential training
    print(f"Placeholder: train on {args.data_dir}, save to {args.output_dir}")


if __name__ == '__main__':
    main()
