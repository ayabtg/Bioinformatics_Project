#!/usr/bin/env python3
"""
Produce basic dataset-level plots:
- chain length distribution
- number of RNA chains per PDB
using data/clean_summary.csv
"""

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SUMMARY_CSV = PROJECT_ROOT / "data" / "clean_summary.csv"
FIG_DIR = PROJECT_ROOT / "figures"

def main():
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(SUMMARY_CSV)

    # Plot 1: total_length distribution (per PDB)
    plt.figure()
    df["total_length"].hist(bins=20)
    plt.xlabel("Total RNA length (residues)")
    plt.ylabel("Number of PDBs")
    plt.title("RNA chain length distribution")
    plt.tight_layout()
    plt.savefig(FIG_DIR / "length_distribution.png")
    plt.close()

    # Plot 2: number of RNA chains per PDB
    plt.figure()
    df["n_chains"].value_counts().sort_index().plot(kind="bar")
    plt.xlabel("Number of RNA chains per PDB")
    plt.ylabel("Count")
    plt.title("Distribution of RNA chains per structure")
    plt.tight_layout()
    plt.savefig(FIG_DIR / "chains_per_pdb.png")
    plt.close()

    print("Saved figures in", FIG_DIR)

if __name__ == "__main__":
    main()

