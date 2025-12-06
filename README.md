## RNA 3D Scoring Pipeline – Full User Manual
This project allows you to:
1. Train a distance-based statistical potential from RNA 3D structures.
2. Use the trained potentials to score new RNA structures.
3. Visualize the learned potentials.

This README is written for a user who can follow the steps exactly in order and will be able to run everything successfully.


## Folder structure:

```
Bioinformatics_Project-main/
├── TP_RNA/
│ ├── data/
│ │ ├── pdb_train/ # Training RNA structures (.pdb)
│ │ └── pdb_test/ # Test RNA structures (.pdb)
│ ├── output/
│ │ ├── potentials/ # Trained potentials (.txt)
│ │ └── profiles/ # Output score files (.csv)
│ └── scripts/
│ ├── train_potential.py
│ ├── score_structure.py
│ └── plot_potential.py
└── README.md
```

## Software Requirements:
- You must have Python installed. Check Python Version (python --version)
- You must have Python 3.8 or newer.

## Install libraries in VScode
```
pip install numpy biopython matplotlib
```
Verify installation: 
```
python -c "import numpy, Bio, matplotlib"
```
If no error appears, the environment is correctly set up.

# Step 1 — Training RNA Statistical Potentials

1. Training Data:
   
Training RNA structures (.pdb format) must be placed in:
```
TP_RNA/data/pdb_train/
```
The provided dataset can be used to verify functionality.

3. Training command:
```
python TP_RNA/scripts/train_potential.py \
    TP_RNA/data/pdb_train \
    --out-dir TP_RNA/output/potentials
```
   
4. Expected Output:
   
The following potential files are generated:
```
AA.txt  AC.txt  AG.txt  AU.txt  
CC.txt  CG.txt  GG.txt  
UC.txt  UG.txt  UU.txt 
```
These files represent trained distance-based statistical energy potentials.

# Step 2 — Verify Scoring with Provided Data (Quick Validation)

This step confirms that the existing trained potentials and test RNA structures work correctly.

Command: 
```
python TP_RNA/scripts/score_structure.py \
    TP_RNA/data/pdb_test \
    --potentials-dir TP_RNA/output/potentials \
    --csv-output TP_RNA/output/profiles/scores.csv
```

Expected Output:
A file is created -
```
TP_RNA/output/profiles/scores.csv
```
If this file is created successfully, the scoring system is working correctly.


# Step 3 — Scoring New RNA Structures
1. Input Data:
   
Place any RNA structures (.pdb format) to be placed in:
```
TP_RNA/data/pdb_test/
```
2. Scoring Command:
```
python TP_RNA/scripts/score_structure.py \
    TP_RNA/data/pdb_test \
    --potentials-dir TP_RNA/output/potentials \
    --csv-output TP_RNA/output/profiles/scores.csv
```
3. Output:
```
TP_RNA/output/profiles/scores.csv
```
Each RNA structure receives a global statistical energy score.

# Step 4 — Plotting the Potentials
1. Command:
```
python TP_RNA/scripts/plot_potential.py \
    --potentials-dir TP_RNA/output/potentials \
    --output-dir TP_RNA/output/plots
```
2. Output:
Plots are saved in:
```
TP_RNA/output/plots/
```
Each plot shows:

X-axis: Distance (Å)

Y-axis: Statistical Energy

Expected plot:
<img width="2968" height="1767" alt="image" src="https://github.com/user-attachments/assets/483c0998-e653-4ff4-825d-7e60d96ae1c2" />

## Additional work: PDB dataset preparation and quality analysis

This section documents the additional work done on the **PDB dataset cleaning and analysis** pipeline in the `future_directions/` folder.

### Overview

The pipeline takes raw PDB structures and:

1. Cleans and validates them (`clean_pdb.py`).
2. Builds final train/test ID lists (`build_splits.py`).
3. Computes simple steric contacts and clashes per structure (`compute_clashes_contacts.py`).
4. Produces dataset-level plots for:
   - sequence length (number of residues)
   - number of chains
   - number of contacts
   - number of clashes

The whole process is reproducible with a few command-line calls.

---

### Directory layout

All paths below are relative to the `future_directions/` folder:

```text
future_directions/
├── pdb_train/           # raw training PDB files
├── pdb_test/            # raw test PDB files
├── pdb_train_clean/     # cleaned training PDBs (created by the pipeline)
├── scripts/
│   ├── clean_pdb.py
│   ├── build_splits.py
│   ├── compute_clashes_contacts.py
│   └── plot_basic_stats.py
├── results/             # CSVs and ID lists
└── figures/             # PNG plots
```
Pipeline Overview:
<img width="2968" height="1767" alt="image" src="https://github.com/ayabtg/Bioinformatics_Project/blob/main/future_directions/figures/pipeline.png" />

We implemented an end-to-end pipeline that (i) cleans and validates the raw RNA PDB structures, (ii) computes steric clashes and residue–residue contacts on the cleaned models, (iii) builds explicit train/test lists, and (iv) generates dataset-level statistics and figures (length, chain count, clash density, contact density).

## Requirements

Python 3.8+

Recommended to use a virtual environment.

Install Python dependencies:

```
pip install pandas matplotlib
```
## Methods

## 1. Cleaning & validating PDB files

This step removes non-structural records and some problematic alternate locations.

From inside future_directions/:
```
python scripts/clean_pdb.py

```
### Input

pdb_train/

### Output

pdb_train_clean/ – cleaned PDB files (same filenames as input)

## 2. Building train/test ID lists

We generate simple lists of PDB IDs (without the .pdb extension) for later use.
```
python scripts/build_splits.py
```
### Outputs (in results/):

train_ids.txt – one PDB ID per line for the cleaned training set

test_ids.txt – one PDB ID per line for the test set

## 3. Computing clashes and contacts

This script parses each PDB, estimates:

number of residues (sequence length),

number of chains,

approximate non-bonded contacts,

approximate steric clashes.

```
python scripts/compute_clashes_contacts.py
```
### Output

results/pdb_stats.csv

Each row corresponds to one structure and contains:

<img width="2968" height="1767" alt="image" src="https://github.com/ayabtg/Bioinformatics_Project/blob/main/future_directions/figures/Screenshot%202025-12-05%20at%2023.41.27.png" />

## 4. Generating dataset-level plots

We visualize the distributions of the basic statistics computed above.
```
python scripts/plot_basic_stats.py
```

### Inputs

results/pdb_stats.csv

### Outputs (in figures/):

length_distribution.png

chains_distribution.png

contacts_distribution.png

clashes_distribution.png

<img width="2968" height="1767" alt="image" src="https://github.com/ayabtg/Bioinformatics_Project/blob/main/future_directions/figures/clashes_distribution.png" />
## 5. End-to-end reproduction

To regenerate all outputs from scratch (assuming pdb_train/ and pdb_test/ already exist):
```
cd future_directions

python scripts/clean_pdb.py
python scripts/build_splits.py
python scripts/compute_clashes_contacts.py
python scripts/plot_basic_stats.py

```
After these commands:

cleaned PDBs are in pdb_train_clean/

ID lists are in results/train_ids.txt and results/test_ids.txt

per-structure statistics are in results/pdb_stats.csv

figures are in figures/

## 6. Summary of generated files
| File / folder                       | Description                               |
| ----------------------------------- | ----------------------------------------- |
| `pdb_train_clean/`                  | Cleaned training PDB files                |
| `results/train_ids.txt`             | Final list of training PDB IDs            |
| `results/test_ids.txt`              | Final list of test PDB IDs                |
| `results/pdb_stats.csv`             | Length, chains, contacts, clashes per PDB |
| `figures/length_distribution.png`   | Histogram of sequence lengths             |
| `figures/chains_distribution.png`   | Histogram of number of chains             |
| `figures/contacts_distribution.png` | Histogram of contacts per structure       |
| `figures/clashes_distribution.png`  | Histogram of clashes per structure        |

## 7. Notes and assumptions

Only ATOM, HETATM, TER, and END records are kept during cleaning.

Atoms with alternate location indicators other than " " or "A" are discarded.

Hydrogens are ignored when computing contacts and clashes.

Contact and clash thresholds are currently fixed at 4.5 Å and 2.0 Å, respectively.
