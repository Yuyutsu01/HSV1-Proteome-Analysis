# Dataset and Biological Annotations

## HSV-1 Proteome Dataset Construction

### 1. Reference Genome Specification
- **Organism:** Human alphaherpesvirus 1 (HSV-1)
- **Strain:** Strain 17 (17+)
- **NCBI RefSeq Accession:** `NC_001806.2`
- **Download Manifest:** `data/raw/ncbi/download_manifest.json`
- **Raw GenBank Record:** `data/raw/ncbi/NC_001806.2.gbk` (SHA256: `aaee9fc50cfaf8f5bf9f2b8b939f8adba061ca81d11ff982f6dbd9fef038b3a0`)

---

### 2. Extraction and Deduplication Cascade

```
NCBI RefSeq NC_001806.2 (77 raw translated CDS records)
  │
  ├── Phase 02: Extraction ─────────────► 77 raw proteins extracted (data/intermediate/raw_proteins.csv)
  │
  ├── Phase 03: Quality Control ────────► 77/77 passed QC (zero internal stop codons, valid IUPAC alphabet)
  │
  └── Phase 04: Deduplication ──────────► 3 exact duplicate sequences removed
                                          │
                                          └──► 74 Unique Proteins (data/processed/unique_proteins.csv)
```

#### Exact Duplicate Mappings
HSV-1 possesses an inverted repeat genomic architecture ($TR_L-U_L-IR_L-IR_S-U_S-TR_S$). In the NCBI RefSeq annotation `NC_001806.2`, three diploid genes residing in the inverted repeat regions are translated identically twice. These duplicate entries were identified by exact sequence hash matching and deduplicated to preserve single-copy proteomic integrity:

| Duplicate Protein ID | Duplicate Gene | Retained Protein ID | Retained Gene | Product Description | Status |
| :--- | :--- | :--- | :--- | :--- | :--- |
| `YP_009137133.1` | `RL2` | `YP_009137074.1` | `RL2` | E3 ubiquitin ligase ICP0 | Exact duplicate removed |
| `YP_009137134.1` | `RL1` | `YP_009137073.1` | `RL1` | Neurovirulence factor ICP34.5 | Exact duplicate removed |
| `YP_009137149.1` | `RS1` | `YP_009137135.1` | `RS1` | Major transcriptional activator ICP4 | Exact duplicate removed |

*Source record:* [`data/intermediate/duplicates.csv`](file:///c:/Users/shiva/OneDrive/Desktop/HSV1-Proteome-Analysis/HSV1-Proteome-Analysis-new/data/intermediate/duplicates.csv)

---

### 3. Protein Sequence Profile ($N=74$)

| Property | Value | Notes / Outliers |
| :--- | :--- | :--- |
| **Total Unique Proteins** | 74 | Single-copy non-redundant dataset |
| **Minimum Length** | 88 aa | `US12` (ICP47, TAP inhibitor) |
| **Maximum Length** | 3,139 aa | `UL36` (Large tegument protein) |
| **Mean Sequence Length** | 522.95 aa | Standard deviation: 441.77 aa |
| **Median Sequence Length** | 400.00 aa | Interquartile range: 247.25 – 633.50 aa |

---

### 4. Biological Temporal Annotations

Biological temporal-expression labels were compiled from authoritative virological literature and independently audited in Phase 6.

#### Temporal Class Breakdown
- **Immediate-Early ($\alpha$ / IE):** 5 proteins (6.76%)
- **Early ($\beta$ / E):** 15 proteins (20.27%)
- **Late ($\gamma$ / L):** 54 proteins (72.97%)
- **Total:** 74 proteins (100.0%)

#### Canonical Immediate-Early Regulators ($N=5$)
1. `YP_009137074.1` — **RL2** / ICP0 (E3 ubiquitin ligase)
2. `YP_009137135.1` — **RS1** / ICP4 (Major transcriptional activator)
3. `YP_009137136.1` — **US1** / ICP22 (Host RNA polymerase II modulator)
4. `YP_009137129.1` — **UL54** / ICP27 (Multifunctional mRNA export/splicing regulator)
5. `YP_009137148.1` — **US12** / ICP47 (TAP peptide transporter inhibitor)

#### Authoritative Handling of UL44 (glycoprotein C / gC)
- **Protein ID:** `YP_009137119.1`
- **Gene Symbol:** `UL44`
- **Product:** Glycoprotein C (gC)
- **Primary Class:** **Late** (Unambiguous biological consensus)
- **Late Subclass:** **Conflicting** ($\gamma_1$ leaky-late kinetics reported in early studies vs. $\gamma_2$ true-late strict DNA-replication dependence reported under phosphonoacetic acid / cycloheximide release assays).
- **Verification Status:** `conflicting` (subclass level only).
- **Resolution:** Retained in the primary 3-class dataset under `temporal_class = Late`. Subclass remains explicitly flagged as `Conflicting` without forced subclass assignment.

#### Annotation Audit Outcome
- **73 / 74** annotations verified with unanimous literature agreement.
- **1 / 74** annotation verified at primary class level (`Late`) with documented subclass ambiguity (`UL44`).
- **0** uncertain or unsupported primary temporal annotations.
- **Authoritative file:** [`data/annotations/temporal_annotations_final.csv`](file:///c:/Users/shiva/OneDrive/Desktop/HSV1-Proteome-Analysis/HSV1-Proteome-Analysis-new/data/annotations/temporal_annotations_final.csv)
