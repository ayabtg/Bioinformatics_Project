# README.md – TP_RNA Project

## Overview

This project implements a small RNA 3D structure analysis pipeline based on the TP instructions.

The goal is to:

- Build a dataset of RNA structures from the Protein Data Bank (PDB)
- Learn a statistical potential based on C3′–C3′ distances between nucleotides
- Use this potential to score and evaluate new RNA structures
- Visualise the resulting per-position pseudo-energy profile

This scoring approach is inspired by RNA knowledge-based potentials used in structural biology.
The pseudo-energy approximates the Gibbs free energy of an RNA conformation.

## Project Structure

```
TP_RNA/
│
├── data/
│   ├── pdb_train/          # PDB files for training the potential
│   └── pdb_test/           # PDB files for scoring (RNA-Puzzles, etc.)
│
├── output/
│   ├── potentials/         # 10 base-pair potentials (20 bins each)
│   ├── profiles/           # scoring profiles (CSV)
│   └── plots/              # profile plots (PNG/PDF)
│
├── scripts/
│   ├── train_potential.py  # training script
│   ├── score_structure.py  # scoring script
│   └── plot_profiles.R     # plotting script
│
└── README.md
```

## What the project does

### 1. Dataset creation (PDB files)

You must download several RNA-only structures from the RCSB PDB:

- only RNA, no DNA, no proteins
- preferably X-ray structures
- small hairpins, duplexes, tRNA fragments, etc.

Place them here:

```
data/pdb_train/   → used to train the potential
data/pdb_test/    → used to evaluate structures
```

### 2. Learning the RNA statistical potential (train_potential.py)

This script follows exactly the TP consignes:

**✔ Extract distances**

- Only C3′ atoms
- Only intrachain distances
- Only residue pairs separated by at least 3 positions: (i, i+4), (i, i+5), (i, i+6), …
- Only distances ≤ 20 Å

**✔ Build 10 base-pair distributions**

For each pair type: AA, AU, AC, AG, UU, UC, UG, CC, CG, GG

Compute a distance histogram with 20 bins between 0 and 20 Å.

**✔ Compute observed frequencies**

```
f_obs(i,j,r) = N(i,j,r) / N(i,j)
```

**✔ Compute reference frequency (all bases = "X")**

```
f_ref(r) = N(X,X,r) / N(X,X)
```

**✔ Compute pseudo-energy (TP formula)**

```
ū(i,j,r) = − log( f_obs(i,j,r) / f_ref(r) )
```

**✔ Output**

The script writes 10 files, one for each base-pair type, each containing 20 values (1 per bin), capped at a maximum value of 10.

Saved to: `output/potentials/`

### 3. Scoring new RNA structures (score_structure.py)

This script evaluates a structure from pdb_test/ using the trained potential.

It repeats the same extraction rules:

- C3′ atoms
- intrachain
- residues i and i+4 or more
- distances ≤ 20 Å

**✔ For each distance:**

- determine base-pair type (AA, AU, …)
- find the two closest bins
- use linear interpolation to estimate score
- sum all scores for each position → per-position profile

**✔ Output:**

```
output/profiles/<name>.profile.csv
```

Example:

```
position,score
1,-0.23
2,1.45
3,0.11
...
```

You should also compute: the total pseudo-energy of the structure

### 4. Plotting the profiles (plot_profiles.R)

This script uses R + ggplot2 to produce scoring profile plots.

**Input:**

```
output/profiles/<name>.profile.csv
```

**Output:**

```
output/plots/<name>.png
```

The plot helps identify:

- Low-energy regions → good / native-like
- High-energy regions → unusual / possibly misfolded

## How to run the project

Train the potential:

```bash
 python TP_RNA/scripts/train_potential.py --data-dir data/pdb_train --output-dir output/potentials
```

Score a structure:

```bash
python scripts/score_structure.py data/pdb_test/<name>.pdb output/profiles/<name>.profile.csv --pot-dir output/potentials
```

Plot the profile:

```bash
Rscript scripts/plot_profiles.R output/profiles/<name>.profile.csv output/plots/<name>.png
```

##  Requirements

- Python 3
- R + ggplot2
- Recommended Python packages:
  - numpy
  - pandas
  - biopython (optional, for easier PDB parsing)

## Quick usage

- Put training PDBs in `data/pdb_train/`
- Put test PDBs in `data/pdb_test/`
- Run `train_potential.py` to compute potentials
- Run `score_structure.py` to evaluate a structure
- Run `plot_profiles.R` to visualise the profile
