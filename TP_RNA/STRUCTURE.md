# Project structure — TP_RNA

This file describes the directory layout for the RNA practical (TP_RNA) and gives usage examples.

Directory layout

Bioinformatics_Project/
│
├── README.md                           # global documentation for the entire repository
├── .gitignore                          # Git ignore rules
│
└── TP_RNA/                             # RNA practical project (main module)
    │
    ├── data/                           # input RNA 3D structures
    │   ├── pdb_train/                  # PDB files used to train the potential
    │   └── pdb_test/                   # PDB files used for testing/scoring
    │
    ├── output/                         # automatically generated outputs
    │   ├── profiles/                   # score profiles (CSV)
    │   ├── plots/                      # profile plots (PNG/PDF)
    │   └── potentials/                 # learned potentials (AA, AU, AC, …)
    │
    ├── scripts/                        # core code of the project
    │   ├── train_potential.py          # compute statistical potentials from pdb_train/
    │   ├── score_structure.py          # score a structure using learned potentials
    │   └── plot_profiles.R             # plot score profiles
    │
    ├── STRUCTURE.md                    # detailed explanation of folder layout
    └── .gitkeep                        # keeps empty folders tracked by Git



Conventions and expected formats

- `data/pdb_train/`: put PDB files used for training here. Filenames should keep the PDB identifier (e.g. `1ABC.pdb`).
- `data/pdb_test/`: PDB files used for testing or benchmarking.
- `output/potentials/`: learned statistical potentials per base-pair type. Each file contains 20 energy values (one per distance bin, 1–20 Å).
  - Recommended filenames: `AA.potential`, `AU.potential`, `AC.potential`, `AG.potential`, `UU.potential`, `UC.potential`, `UG.potential`, `CC.potential`, `CG.potential`, `GG.potential`.
  - Format: plain text or CSV with one energy value per line (or per distance bin).
- `output/profiles/`: each scored PDB should produce a CSV profile. Recommended filename: `<pdb_id>.profile.csv`.
  - Minimal CSV format expected by `plot_profiles.R`: header `position,score`.
  - Example:

    position,score
    1,0.12
    2,0.45
    3,-0.03

- `output/plots/`: images created from profiles (PNG, PDF). Recommended filename: `<pdb_id>.profile.png`.

Scripts

- `scripts/train_potential.py` (template):
  - CLI options:
    - `--data-dir` : directory with training PDBs (default `../data/pdb_train`).
    - `--output-dir`: directory to save potential files (default `../output/potentials`).
  - Workflow: reads PDB files, extracts C3'–C3' distances for intrachain pairs (|i–j| ≥ 4, distance ≤ 20 Å), bins distances (1 Å), counts base-pair frequencies, and writes 10 potential files (one per base-pair type).
  - Example (PowerShell):

    python .\TP_RNA\scripts\train_potential.py --data-dir .\TP_RNA\data\pdb_train --output-dir .\TP_RNA\output\potentials

- `scripts/score_structure.py` (template):
  - Workflow: reads a PDB file, loads pre-computed potentials from `output/potentials/`, extracts C3'–C3' distances (same criteria as training), looks up energies in the potentials, and sums them to produce a total score saved as a CSV profile.
  - Usage (PowerShell):

    python .\TP_RNA\scripts\score_structure.py .\TP_RNA\data\pdb_test\example.pdb .\TP_RNA\output\profiles\example.profile.csv

Additional notes:

- `train_potential.py` produces potential files (e.g. `AA`, `AU`, etc.) saved in `output/potentials/`.
- `score_structure.py` produces per-structure profile CSVs with `position,score` columns and writes them to `output/profiles/`.

- `scripts/plot_profiles.R` (R script using `ggplot2`):
  - Usage (PowerShell):

    Rscript .\TP_RNA\scripts\plot_profiles.R .\TP_RNA\output\profiles\example.profile.csv .\TP_RNA\output\plots\example.png

PowerShell quick commands

- List the TP_RNA tree recursively:
```powershell
Get-ChildItem -Path .\TP_RNA -Recurse
```

- List `scripts` contents:
```powershell
Get-ChildItem -Path .\TP_RNA\scripts
```

