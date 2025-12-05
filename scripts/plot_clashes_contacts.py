#!/usr/bin/env python3
"""
Plots for steric clashes and contacts using data/geometry_summary.csv:
 - clashes_per_100_nt histogram
 - contacts_per_residue histogram
 - total_length vs n_contacts scatter
"""

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

PROJECT_ROOT = Path(__file__).resolve().parents[1]
GEOM_CSV = PROJECT_ROOT / "data" / "geometry_summary.csv"
FIG_DIR = PROJECT_ROOT / "figures"

def main():
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(GEOM_CSV)

    # Plot 1: clashes per 100 nt
    plt.figure()
    df["clashes_per_100_nt"].hist(bins=10)
    plt.xlabel("Steric clashes per 100 nt")
    plt.ylabel("Number of PDBs")
    plt.title("Distribution of steric clashes")
    plt.tight_layout()
    plt.savefig(FIG_DIR / "clashes_per_100_nt.png")
    plt.close()

    # Plot 2: contacts per residue
    plt.figure()
    df["contacts_per_residue"].hist(bins=10)
    plt.xlabel("Contacts per residue")
    plt.ylabel("Number of PDBs")
    plt.title("Distribution of contacts per residue")
    plt.tight_layout()
    plt.savefig(FIG_DIR / "contacts_per_residue.png")
    plt.close()

    # Plot 3: length vs total contacts
    plt.figure()
    plt.scatter(df["total_length"], df["n_contacts"])
    plt.xlabel("Total length (nt)")
    plt.ylabel("Number of contacts")
    plt.title("Contacts vs RNA length")
    plt.tight_layout()
    plt.savefig(FIG_DIR / "contacts_vs_length.png")
    plt.close()

    print("Saved clash/contact figures in", FIG_DIR)

if __name__ == "__main__":
    main()

