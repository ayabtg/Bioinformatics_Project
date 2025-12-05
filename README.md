RNA 3D Scoring Pipeline – Full User Manual
This project allows you to:
1. Train a distance-based statistical potential from RNA 3D structures.
2. Use the trained potentials to score new RNA structures.
3. Visualize the learned potentials.

This README is written for a user who can follow the steps exactly in order and will be able to run everything successfully.

Folder structure:
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

Software Requirements:
You must have Python installed. Check Python Version (python --version)
You must have Python 3.8 or newer.
