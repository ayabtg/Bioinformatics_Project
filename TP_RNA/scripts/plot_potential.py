#!/usr/bin/env python3
#Plotting script for RNA statistical potential visualization.
import matplotlib
matplotlib.use("Agg")  # MUST be before pyplot import
import os
import sys
import argparse
from typing import Dict, List, Tuple
import numpy as np
import matplotlib.pyplot as plt

#Argument parsing

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Plot RNA distance-based statistical potentials."
    )
    parser.add_argument(
        "--potentials-dir", "-p",
        type=str,
        default="data/potentials",
        help="Directory containing *.potential files. Default: data/potentials"
    )
    parser.add_argument(
        "--output-dir", "-o",
        type=str,
        default="plots",
        help="Directory to save potential plots. Default: plots"
    )
    return parser.parse_args()

# Load potential data
def load_potential(filepath: str) -> Tuple[List[float], List[float]]:
    """Read one .potential file (distance, score)."""
    distances = []
    scores = []
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            try:
                d, s = line.split()[:2]
                distances.append(float(d))
                scores.append(float(s))
            except Exception:
                print(f"[WARN] Could not parse line: {line}", file=sys.stderr)

    return distances, scores


def load_all_potentials(potentials_dir: str) -> Dict[str, Tuple[List[float], List[float]]]:
    """Load all base-pair potential files into memory."""
    potentials = {}

    for filename in sorted(os.listdir(potentials_dir)):
        if not filename.endswith(".txt"):       # <-- FIXED
            continue

        bp = filename.replace(".txt", "")       # <-- FIXED
        path = os.path.join(potentials_dir, filename)

        distances, scores = load_potential(path)
        if distances:
            potentials[bp] = (distances, scores)
            print(f"[INFO] Loaded {len(distances)} distances for {bp}")

    return potentials

# Plot potential curve
def plot_potential_curve(base_pair, distances, scores, output_path):
    distances = np.array(distances)
    scores = np.array(scores)

    plt.figure(figsize=(10, 6))

    # Line plot
    plt.plot(distances, scores, 'o-', linewidth=2, label=f"{base_pair} potential")

    # Favorable region (scores < 0)
    plt.fill_between(
        distances, scores, 0,
        where=(scores < 0),
        color='green', alpha=0.3,
        label='Favorable (<0)'
    )

    # Unfavorable region (scores > 0)
    plt.fill_between(
        distances, scores, 0,
        where=(scores > 0),
        color='red', alpha=0.3,
        label='Unfavorable (>0)'
    )

    plt.axhline(0, color='black', linewidth=1)

    plt.xlabel("Distance (Å)")
    plt.ylabel("Pseudo-energy Score")
    plt.title(f"Statistical Potential: {base_pair}")
    plt.grid(alpha=0.3)
    plt.legend()

    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()

    print(f"[INFO] Saved potential plot → {output_path}")

# Main
def main():
    args = parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    print("\n[INFO] Loading potentials...")
    potentials = load_all_potentials(args.potentials_dir)

    if not potentials:
        print("[ERROR] No .potential files found!", file=sys.stderr)
        sys.exit(1)

    print("\n[INFO] Generating potential plots...")

    for bp, (distances, scores) in potentials.items():
        output_path = os.path.join(args.output_dir, f"{bp}_potential.png")
        plot_potential_curve(bp, distances, scores, output_path)

    print(f"\n[SUCCESS] All potential plots saved in: {args.output_dir}/")


if __name__ == "__main__":
    main()
