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

# Step 1 — Verify Scoring with Provided Data (Quick Validation)

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

# Step 2 — Training RNA Statistical Potentials

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

# Step 3 — Scoring New RNA Structures
1. Input Data:
   
Place any RNA structures (.pdb format) to be placed in:
```
TP_RNA/data/pdb_set/
```
2. Scoring Command:
```
python TP_RNA/scripts/score_structure.py \
    TP_RNA/data/pdb_set \
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
