#!/usr/bin/env python
import argparse
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", type=Path, default=Path("results/pdb_stats.csv"))
    ap.add_argument("--outdir", type=Path, default=Path("figures"))
    args = ap.parse_args()

    args.outdir.mkdir(exist_ok=True)

    df = pd.read_csv(args.csv)

    # lengths
    plt.figure()
    df["length"].hist(bins=20)
    plt.xlabel("Sequence length (residues)")
    plt.ylabel("Count")
    plt.title("PDB length distribution")
    plt.tight_layout()
    plt.savefig(args.outdir / "length_distribution.png")
    plt.close()

    # chains
    plt.figure()
    df["chains"].hist(bins=range(1, df["chains"].max() + 2))
    plt.xlabel("Number of chains")
    plt.ylabel("Count")
    plt.title("Number of chains per PDB")
    plt.tight_layout()
    plt.savefig(args.outdir / "chains_distribution.png")
    plt.close()

    # contacts
    plt.figure()
    df["contacts"].hist(bins=20)
    plt.xlabel("Contacts per structure")
    plt.ylabel("Count")
    plt.title("Contacts distribution")
    plt.tight_layout()
    plt.savefig(args.outdir / "contacts_distribution.png")
    plt.close()

    # clashes
    plt.figure()
    df["clashes"].hist(bins=20)
    plt.xlabel("Clashes per structure")
    plt.ylabel("Count")
    plt.title("Clashes distribution")
    plt.tight_layout()
    plt.savefig(args.outdir / "clashes_distribution.png")
    plt.close()

    print(f"Saved plots in {args.outdir}")

if __name__ == "__main__":
    main()

