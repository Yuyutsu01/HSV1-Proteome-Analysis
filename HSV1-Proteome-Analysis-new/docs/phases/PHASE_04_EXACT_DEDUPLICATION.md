# Phase 04: Exact Deduplication

## 1. Objective
Identify and remove exact sequence duplicates arising from the diploid inverted repeat segments of the HSV-1 genome, producing an authoritative single-copy proteome dataset.

## 2. Inputs
- `data/intermediate/qc_passed_proteins.csv` (77 proteins)

## 3. Processing
- Executed `scripts/04_exact_deduplication.py`.
- Computed SHA-256 sequence hashes for all 77 proteins.
- Grouped identical sequences and retained the canonical primary locus (lower genomic coordinate / standard naming).
- Logged all duplicate mappings into `duplicates.csv`.

## 4. Outputs
- `data/processed/unique_proteins.csv` (**74 unique proteins**)
- `data/processed/unique_proteins.fasta` (74 FASTA records)
- `data/intermediate/duplicates.csv` (3 duplicate records)
- `results/tables/dataset_summary.csv` (Summary counts)

## 5. Validation
- **Exact duplicates identified:** 3
- **Unique single-copy proteins:** 74
- **Deduplication Verification:**

| Duplicate Protein ID | Duplicate Gene | Retained Protein ID | Retained Gene | Reason |
| :--- | :--- | :--- | :--- | :--- |
| `YP_009137133.1` | `RL2` | `YP_009137074.1` | `RL2` | Exact sequence match (ICP0) |
| `YP_009137134.1` | `RL1` | `YP_009137073.1` | `RL1` | Exact sequence match (ICP34.5) |
| `YP_009137149.1` | `RS1` | `YP_009137135.1` | `RS1` | Exact sequence match (ICP4) |

## 6. Scientific Decisions
- Removed **only exact 100% identical sequence duplicates**. Homologous gene families with sequence divergence (e.g., US10/US11, glycoprotein homologs) were strictly preserved as distinct unique proteins.

## 7. Important Limitations
- Reduces dataset to 74 unique polypeptides; diploid gene dosage effects in viral replication are not represented in static sequence feature matrices.

## 8. Files Generated
- [`data/processed/unique_proteins.csv`](file:///c:/Users/shiva/OneDrive/Desktop/HSV1-Proteome-Analysis/HSV1-Proteome-Analysis-new/data/processed/unique_proteins.csv)
- [`data/processed/unique_proteins.fasta`](file:///c:/Users/shiva/OneDrive/Desktop/HSV1-Proteome-Analysis/HSV1-Proteome-Analysis-new/data/processed/unique_proteins.fasta)
- [`data/intermediate/duplicates.csv`](file:///c:/Users/shiva/OneDrive/Desktop/HSV1-Proteome-Analysis/HSV1-Proteome-Analysis-new/data/intermediate/duplicates.csv)
- [`results/tables/dataset_summary.csv`](file:///c:/Users/shiva/OneDrive/Desktop/HSV1-Proteome-Analysis/HSV1-Proteome-Analysis-new/results/tables/dataset_summary.csv)

## 9. Status
**COMPLETE & VALIDATED**
