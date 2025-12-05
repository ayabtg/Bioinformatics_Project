# 📘 User Manual — RNA 3D Scoring Pipeline  
*(Parameterized Version)*

This document explains how to use the parameterized Python scripts to:

1. **Train statistical potentials** from RNA training structures.  
2. **Score a new RNA structure** using the learned potentials.  
3. **Customize all parameters** using command-line flags (CLI parametrization).

---

## 1️⃣ Requirements

Make sure you have:

- Python 3 installed  
- Required libraries (`numpy`, `Bio.PDB`, etc.)  
- RNA structure files organized as follows:

TP_RNA/
├── data/
│ ├── pdb_train/ ← PDB files for training
│ └── pdb_test/ ← PDB files for scoring
├── scripts/
│ ├── train_potential.py
│ └── score_structure.py
└── output/
├── potentials/ ← generated potentials
└── profiles/ ← generated scoring profiles

yaml
Copier le code

---

## 2️⃣ Training Potentials (`train_potential.py`)

This script learns statistical pairwise potentials (AA, AC, AG, … UU) from RNA 3D structures.

### **General Command**

```bash
python3 scripts/train_potential.py \
    --data-dir <TRAINING_DATA_DIRECTORY> \
    --output-dir <OUTPUT_DIRECTORY> \
    [--bins N] \
    [--cutoff DIST] \
    [--min-sep S] \
    [--atom ATOM_NAME]
Available Parameters
Parameter	Description	Default
--data-dir	Directory containing PDB training structures	required
--output-dir	Where potential files will be saved	required
--bins	Number of distance bins	20
--cutoff	Maximum distance (Å)	20
--min-sep	Minimum	i-j
--atom	Atom used for distance computation	"C3'"

Example
bash
Copier le code
python3 scripts/train_potential.py \
    --data-dir data/pdb_train \
    --output-dir output/potentials \
    --bins 20 \
    --cutoff 20 \
    --min-sep 4 \
    --atom "C3'"
This creates files such as:

python-repl
Copier le code
AA.potential
AC.potential
AG.potential
...
UU.potential
3️⃣ Scoring an RNA Structure (score_structure.py)
This script evaluates an RNA structure using the trained potentials and produces a CSV scoring profile.

General Command
bash
Copier le code
python3 scripts/score_structure.py \
    --input-pdb <PDB_FILE_TO_SCORE> \
    --potentials-dir <POTENTIALS_DIRECTORY> \
    --output-csv <OUTPUT_CSV_FILE> \
    [--cutoff DIST] \
    [--min-sep S] \
    [--atom ATOM_NAME]
Parameters
Parameter	Description	Default
--input-pdb	Structure to score	required
--potentials-dir	Directory containing .potential files	required
--output-csv	Output CSV path	required
--cutoff	Maximum atom-atom distance	20
--min-sep	Minimum	i-j
--atom	Atom to analyze	"C3'"

Example
bash
Copier le code
python3 scripts/score_structure.py \
    --input-pdb data/pdb_test/1EHZ.pdb \
    --potentials-dir output/potentials \
    --output-csv output/profiles/1EHZ_profile.csv \
    --cutoff 20 \
    --min-sep 4 \
    --atom "C3'"
Output file:

bash
Copier le code
output/profiles/1EHZ_profile.csv
The CSV contains, for each nucleotide:

residue index

residue type (A, U, G, C)

total interaction score

individual contributions

4️⃣ Notes & Best Practices
✔️ Changing atoms
The scoring behavior changes depending on which atom is chosen ("C3'", "C4'", "P"…).

✔️ No need to modify the code
All behavior is controlled through CLI flags — this is the purpose of parametrization.

✔️ Reproducibility
Each run documents its own parameters, making comparisons straightforward.

5️⃣ Quick Reference (Cheat Sheet)
Train:
bash
Copier le code
python3 scripts/train_potential.py \
    --data-dir data/pdb_train \
    --output-dir output/potentials
Score:
bash
Copier le code
python3 scripts/score_structure.py \
    --input-pdb data/pdb_test/1XYZ.pdb \
    --potentials-dir output/potentials \
    --output-csv output/profiles/1XYZ_profile.csv
