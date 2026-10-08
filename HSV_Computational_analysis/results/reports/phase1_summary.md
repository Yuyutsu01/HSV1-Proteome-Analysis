# Phase 1 Dataset Construction & Raw Data Ingestion Report

**Ingestion Timestamp**: `2026-10-08T12:52:06.846093+00:00`  
**Raw Source**: `C:\Users\shiva\OneDrive\Desktop\HSV1-Proteome-Analysis\HSV_Computational_analysis\raw`  

## 1. Summary Metrics

| Metric | Count | Description |
| :--- | :--- | :--- |
| **RAW_RECORD_COUNT** | `22,710` | Total raw FASTA records ingested from source files |
| **VALID_SEQUENCE_COUNT** | `22,710` | Raw records passing objective sequence quality criteria |
| **EXCLUDED_RECORD_COUNT** | `0` | Raw records failing sequence quality criteria |
| **UNIQUE_SEQUENCE_COUNT** | `22,689` | Distinct non-redundant high-quality sequences |
| **EXACT_DUPLICATE_RECORDS** | `21` | Redundant sequence instances collapsed into canonical entries |

## 2. Ingested Raw Datasets

| Source File | SHA-256 Checksum | File Size (Bytes) | Record Count |
| :--- | :--- | :--- | :--- |
| `HSV2_non_redundant (2).fasta` | `86c9d627286aa4ca...` | 7,095,643 | 10,540 |
| `non_redundant.fasta` | `15001269a835b412...` | 7,954,295 | 12,170 |

## 3. Quality Control (QC) Breakdown

| QC Status | Record Count | Category Description |
| :--- | :--- | :--- |
| `CANONICAL_20AA` | 15,229 | Evaluated by IUPAC standard rule |
| `AMBIGUOUS_RESIDUE` | 7,481 | Evaluated by IUPAC standard rule |

## 4. Provenance and Immutability Guarantee

- Original files in `./raw/` remained strictly read-only and unmutated.
- Immutability verified via SHA-256 cryptographic fingerprints.
- Exact 1-to-1 provenance mapping preserved in `data/intermediate/raw_derived_sequences.tsv`.
- Every canonical sequence references its exact constituent raw records in `data/processed/unique_high_quality_sequences.tsv`.
