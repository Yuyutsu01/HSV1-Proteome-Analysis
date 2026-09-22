# Phase 01: NCBI Data Acquisition

## 1. Objective
Download and cryptographically verify the authoritative reference genome for Human alphaherpesvirus 1 (strain 17) directly from the NCBI RefSeq database.

## 2. Inputs
- **RefSeq Accession ID:** `NC_001806.2`
- **Database:** NCBI Entrez / Nucleotide
- **Format:** Full GenBank flat file (`.gbk`) with complete feature annotations

## 3. Processing
- Executed `scripts/01_download_ncbi.py` to retrieve the complete GenBank record.
- Calculated SHA-256 hash immediately upon download to guarantee file integrity.
- Generated `download_manifest.json` capturing metadata, accession, timestamp, and SHA-256 hash.

## 4. Outputs
- `data/raw/ncbi/NC_001806.2.gbk` (152,222 bp double-stranded linear DNA)
- `data/raw/ncbi/download_manifest.json`

## 5. Validation
- **SHA-256 Hash:** `aaee9fc50cfaf8f5bf9f2b8b939f8adba061ca81d11ff982f6dbd9fef038b3a0`
- Verified valid GenBank structure containing complete CDS annotations with `/translation` tags.

## 6. Scientific Decisions
- **Reference Strain Selection:** HSV-1 strain 17 (`NC_001806.2`) is the standard laboratory reference genome with comprehensive experimental characterization. Multi-strain analysis was not performed in this phase.

## 7. Important Limitations
- Single reference strain dataset; strain-specific sequence polymorphisms and clinical isolate variations are excluded from this baseline reference.

## 8. Files Generated
- [`data/raw/ncbi/NC_001806.2.gbk`](file:///c:/Users/shiva/OneDrive/Desktop/HSV1-Proteome-Analysis/HSV1-Proteome-Analysis-new/data/raw/ncbi/NC_001806.2.gbk)
- [`data/raw/ncbi/download_manifest.json`](file:///c:/Users/shiva/OneDrive/Desktop/HSV1-Proteome-Analysis/HSV1-Proteome-Analysis-new/data/raw/ncbi/download_manifest.json)

## 9. Status
**COMPLETE & VALIDATED**
