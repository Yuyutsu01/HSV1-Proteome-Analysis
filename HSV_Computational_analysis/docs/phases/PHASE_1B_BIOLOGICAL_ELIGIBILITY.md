# Phase 1B: Biological Eligibility & Source-Type Curation

## 1. Objective
Establish objective biological provenance, taxonomy, completeness, and source-type classifications for all 22,689 unique sequences derived from Phase 1A.

## 2. Evidence Hierarchy
1. Explicit FASTA header brackets (`[...]`) and database identifiers.
2. Verified accession prefixes (`sp|`, `tr|`, `prf|`, `pdb|`, `patent`).
3. UniProtKB mnemonic identifiers.
4. No sequence homology or machine learning inference was permitted during Phase 1B.

## 3. Curated Metrics Overview
- **Raw Ingested Records**: `22,710`
- **Unique Canonical Sequences**: `22,689`

### Taxonomy (Virus Species)
- **HSV-1**: `11,802`
- **HSV-2**: `10,263`
- **OTHER_HERPESVIRUS**: `40`
- **NON_HERPESVIRUS**: `371`
- **UNKNOWN / CONFLICTING**: `213`

### Biological Eligibility Breakdown
- **ELIGIBLE_PRIMARY**: `13,320`
- **ELIGIBLE_SEQUENCE_ANALYSIS_ONLY**: `8,692`
- **EXCLUDE_NON_HSV**: `411`
- **EXCLUDE_SYNTHETIC**: `49`
- **EXCLUDE_TECHNICAL**: `146`
- **EXCLUDE_ENGINEERED**: `50`
- **EXCLUDE_UNRESOLVED_IDENTITY**: `14`
- **REVIEW_REQUIRED**: `7`

## 4. Absolute Constraints & Governance
- **Temporal Annotation**: NOT PERFORMED (Strictly reserved for Phase 2).
- **Machine Learning / ProtBERT**: NONE.
- **Physical Features**: NONE.
- **Data Modification**: Raw `./raw/` files remained strictly read-only and unmutated.
