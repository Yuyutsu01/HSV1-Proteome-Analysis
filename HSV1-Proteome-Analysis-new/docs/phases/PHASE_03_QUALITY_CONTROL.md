# Phase 03: Quality Control

## 1. Objective
Perform rigorous sequence-level quality control on all 77 extracted HSV-1 protein sequences to confirm translational validity before deduplication and representation learning.

## 2. Inputs
- `data/intermediate/raw_proteins.csv`

## 3. Processing
- Executed `scripts/03_quality_control.py`.
- Checked each sequence for:
  1. Standard 20 IUPAC amino acid alphabet validity (no ambiguous `X`, `B`, `Z`, `J`).
  2. Absence of internal stop codons (`*`).
  3. Minimum sequence length constraint ($\ge 30$ aa).
  4. Non-null identifiers and descriptions.

## 4. Outputs
- `data/intermediate/qc_passed_proteins.csv` (77 records)
- `data/intermediate/qc_report.csv` (77 records with pass/fail flags and QC metrics)
- `results/logs/qc_summary.json` (Summary metrics)

## 5. Validation
- **Pass Rate:** **77 / 77 proteins passed QC (100.0%)**.
- **Internal Stop Codons:** 0 detected.
- **Ambiguous Amino Acids:** 0 detected.
- **Minimum Length:** 88 aa (`US12`), safely above the 30 aa threshold.

## 6. Scientific Decisions
- Verified that all RefSeq CDS translations represent mature, unambiguous amino acid sequences suitable for direct biochemical and neural network analysis.

## 7. Important Limitations
- QC validates sequence syntax and completeness; biological post-translational modifications (cleavages, phosphorylation, glycosylation) are not modeled in raw primary sequences.

## 8. Files Generated
- [`data/intermediate/qc_passed_proteins.csv`](file:///c:/Users/shiva/OneDrive/Desktop/HSV1-Proteome-Analysis/HSV1-Proteome-Analysis-new/data/intermediate/qc_passed_proteins.csv)
- [`data/intermediate/qc_report.csv`](file:///c:/Users/shiva/OneDrive/Desktop/HSV1-Proteome-Analysis/HSV1-Proteome-Analysis-new/data/intermediate/qc_report.csv)
- [`results/logs/qc_summary.json`](file:///c:/Users/shiva/OneDrive/Desktop/HSV1-Proteome-Analysis/HSV1-Proteome-Analysis-new/results/logs/qc_summary.json)

## 9. Status
**COMPLETE & VALIDATED**
