# RNA Structural Dataset Processing – TP_RNA Project

This project builds an RNA structural dataset from multiple PDB files, performs
cleaning and validation, computes steric clashes and residue–residue contacts,
generates dataset-level statistics, and produces train/test ID lists along with
figures for analysis.

The project follows the instructions of the TP_RNA assignment and uses Python
scripts (Biopython + pandas + matplotlib) to process and analyse RNA structures.

---

## 📁 Repository Structure

```
TP_RNA/
│
├── data/
│   ├── pdb_train/               # Raw training PDB files
│   ├── pdb_test/                # Raw test PDB files
│   ├── pdb_clean/               # Cleaned PDBs (generated)
│   │   ├── train/
│   │   └── test/
│   ├── clean_summary.csv        # Cleaning statistics
│   ├── geometry_summary.csv     # Clash + contact statistics
│   ├── train_ids.txt            # List of training PDB IDs
│   └── test_ids.txt             # List of test PDB IDs
│
├── scripts/
│   ├── clean_pdb.py             # Cleaning + validation
│   ├── build_splits.py          # Generate ID lists
│   ├── compute_clashes_contacts.py # Compute steric clashes & contacts
│   ├── plot_basic_stats.py      # Length + chain count plots
│   └── plot_clashes_contacts.py # Clash & contact plots
│
├── figures/
│   ├── length_distribution.png
│   ├── chains_per_pdb.png
│   ├── clashes_per_100_nt.png
│   ├── contacts_per_residue.png
│   └── contacts_vs_length.png
│
└── README.md
```

---

## 🚀 How to Run the Pipeline

Run all steps from the **repo root (`TP_RNA/`)**.

### 1️⃣ Clean and validate PDB files

```bash
python scripts/clean_pdb.py
```

Outputs:
- Cleaned PDBs → `data/pdb_clean/`
- Cleaning statistics → `data/clean_summary.csv`

---

### 2️⃣ Build train/test ID lists

```bash
python scripts/build_splits.py
```

Outputs:
- `data/train_ids.txt`
- `data/test_ids.txt`

---

### 3️⃣ Compute steric clashes & residue–residue contacts

```bash
python scripts/compute_clashes_contacts.py
```

Outputs:
- `data/geometry_summary.csv`  
  (clashes, contacts, normalised metrics)

---

### 4️⃣ Generate dataset-level plots

```bash
python scripts/plot_basic_stats.py
python scripts/plot_clashes_contacts.py
```

Outputs saved in:
- `figures/`

---

## 🔬 Methods Summary

### **Cleaning & Validation**
Using Biopython, each PDB is processed to:
- keep only RNA residues (A, C, G, U)
- remove water, hetero-atoms, non-RNA entities
- keep only primary altLoc atoms
- check backbone completeness (P, O5’, C5’, C4’, C3’, O3’)
- count RNA chains and total RNA length

Results recorded in `clean_summary.csv`.

---

### **Steric Clashes**
A steric clash is defined as a pair of heavy atoms with:

\[
(r_1 + r_2) - d > 0.4 \, Å
\]

Computed using NeighborSearch.

### **Contacts**
Two residues are considered in contact if *any* atom pair is within **4.5 Å**.

### **Statistics**
We compute:
- clashes per 100 nucleotides  
- contacts per residue  
- total contacts  
- length vs contact scaling  

---

## 📊 Dataset-Level Results

Figures generated include:

- **RNA Length Distribution**  
- **Number of RNA Chains per PDB**
- **Clashes per 100 Nucleotides**
- **Contacts per Residue**
- **Contacts vs Length** (scatter plot)

These figures allow comparison of RNA packing density and interaction patterns across structures.

---

## 🔁 Pipeline Diagram

(Insert your diagram image here, e.g. `pipeline.png`)

---

## 📝 Interpretation Summary

- RNA lengths range from very small RNAs (18–24 nt) up to large complexes (940 nt).  
- Longer RNAs naturally show higher absolute counts of contacts and steric clashes.  
- When normalised (clashes per 100 nt, contacts per residue), RNAs display similar packing densities.  
- The contact vs length plot confirms that larger RNAs form more 3D interactions.

These patterns are consistent with RNA structural behaviour in the PDB.

---

## ✔️ Final Notes

This pipeline is modular, reproducible, and easy to extend
(e.g., add energy potentials, secondary structure parsing, or RMSD comparisons).

All scripts were successfully tested on the provided dataset.

