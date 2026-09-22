# Phase 02: Protein Extraction

## 1. Objective
Extract all translated protein-coding sequence (CDS) products from the downloaded NCBI GenBank file `NC_001806.2.gbk` and format them into structured tabular and FASTA datasets.

## 2. Inputs
- `data/raw/ncbi/NC_001806.2.gbk`

## 3. Processing
- Executed `scripts/02_extract_proteins.py` using Biopython's `SeqIO` parser.
- Extracted every `CDS` feature containing a valid `/protein_id`, `/gene`, `/product`, and `/translation`.
- Structured the metadata into a tabular CSV and generated a standard FASTA file.

## 4. Outputs
- `data/intermediate/raw_proteins.csv` (77 records)
- `data/raw/fasta/raw_proteins.fasta` (77 sequence records)

## 5. Validation
- Exactly 77 translated CDS features extracted.
- Zero missing identifiers (`protein_id`, `gene`, `product`, `sequence`).
- Sequence length consistency confirmed between FASTA and CSV.

## 6. Scientific Decisions
- Extracted all translated CDS entries annotated in the authoritative RefSeq record without pre-filtering, capturing both single-copy unique genes and diploid inverted-repeat copies.

## 7. Important Limitations
- Includes duplicate translations originating from inverted repeat genomic regions ($TR_L/IR_L$ and $IR_S/TR_S$).

## 8. Files Generated
- [`data/intermediate/raw_proteins.csv`](file:///c:/Users/shiva/OneDrive/Desktop/HSV1-Proteome-Analysis/HSV1-Proteome-Analysis-new/data/intermediate/raw_proteins.csv)
- [`data/raw/fasta/raw_proteins.fasta`](file:///c:/Users/shiva/OneDrive/Desktop/HSV1-Proteome-Analysis/HSV1-Proteome-Analysis-new/data/raw/fasta/raw_proteins.fasta)

## 9. Status
**COMPLETE & VALIDATED**
