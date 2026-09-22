# Phase 07: Physicochemical Feature Extraction

## 1. Objective
Extract exactly 25 documented physicochemical and compositional sequence descriptors for all 74 unique HSV-1 proteins, calculate parametric and non-parametric descriptive statistics, and generate publication-quality distribution figures.

## 2. Inputs
- `data/processed/unique_proteins.csv` (74 unique proteins)
- `data/processed/unique_proteins.fasta`
- `data/annotations/temporal_annotations_final.csv`

## 3. Processing
- Executed `scripts/07_physicochemical_features.py` using Biopython `ProtParam.ProteinAnalysis`.
- Extracted exactly 25 documented features:
  1. `length`: Sequence length (aa)
  2. `molecular_weight`: Molecular mass (Da)
  3. `aromaticity`: Relative frequency of Phe, Trp, Tyr
  4. `instability_index`: Guruprasad empirical instability index
  5. `isoelectric_point`: Theoretical isoelectric point (pI)
  6–25. `fraction_A` through `fraction_Y`: Relative frequency of each of the 20 standard amino acids.
- Computed Pearson and Spearman correlation matrices.
- Executed `tests/test_phase7_features.py` (8 automated test cases).

## 4. Outputs
- `data/processed/physicochemical_features.csv` ($74 \times 25$ feature matrix)
- `data/processed/feature_metadata.csv` (25 feature metadata entries)
- `results/tables/physicochemical_statistics.csv` (Parametric/non-parametric statistics)
- `results/tables/pearson_correlation_matrix.csv` ($25 \times 25$)
- `results/tables/spearman_correlation_matrix.csv` ($25 \times 25$)
- `results/figures/feature_distributions.png`
- `results/figures/class_feature_distributions.png`
- `results/figures/amino_acid_composition_profile.png`
- `results/figures/physicochemical_correlation_matrix.png`
- `results/logs/phase7_validation_report.txt`

## 5. Validation and Exact Descriptive Statistics

| Feature | Mean | Standard Deviation | Min | Median | Max |
| :--- | :---: | :---: | :---: | :---: | :---: |
| **Sequence Length (aa)** | 522.95 | 441.77 | 88.00 | 400.00 | 3,139.00 |
| **Molecular Weight (Da)** | 56,458.95 | 47,210.14 | 9,201.70 | 43,622.75 | 333,438.10 |
| **Isoelectric Point (pI)** | 7.76 | 1.80 | 4.85 | 7.42 | 11.80 |
| **Instability Index** | 51.43 | 11.18 | 31.68 | 51.30 | 95.54 |
| **Aromaticity** | 0.0691 | 0.0218 | 0.0062 | 0.0711 | 0.1351 |

- **Validation Suite:** `tests/test_phase7_features.py` (**8/8 tests passed**).
  - Missing/NaN/Infinite values: **0**
  - Sequence length match with FASTA: **100%**
  - Amino acid composition sum: $\sum \text{fraction}_{AA} = 1.000000 \pm 10^{-6}$ for all 74 proteins.
  - Aromaticity valid bounds: $[0, 1]$

## 6. Scientific Decisions
- **Exactly 25 Documented Features:** The original manuscript text claimed 26 features but explicitly defined only 25. The revised pipeline adheres strictly to the 25 documented features and does not invent or infer an undocumented 26th feature.
- **Terminology:** Proteins with instability index $>40$ are described neutrally as "potentially unstable according to the empirical instability index" rather than claiming proven intrinsic structural disorder.

## 7. Important Limitations
- Features describe primary sequence properties and amino acid composition; 3D structural conformations, quaternary assemblies, and membrane topologies are not directly modeled.

## 8. Files Generated
- [`data/processed/physicochemical_features.csv`](file:///c:/Users/shiva/OneDrive/Desktop/HSV1-Proteome-Analysis/HSV1-Proteome-Analysis-new/data/processed/physicochemical_features.csv)
- [`results/tables/physicochemical_statistics.csv`](file:///c:/Users/shiva/OneDrive/Desktop/HSV1-Proteome-Analysis/HSV1-Proteome-Analysis-new/results/tables/physicochemical_statistics.csv)
- [`tests/test_phase7_features.py`](file:///c:/Users/shiva/OneDrive/Desktop/HSV1-Proteome-Analysis/HSV1-Proteome-Analysis-new/tests/test_phase7_features.py)

## 9. Status
**COMPLETE & VALIDATED**
