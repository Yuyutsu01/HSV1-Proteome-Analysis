# Dataset Card: Curated HSV Protein Sequence Corpus

## 1. Dataset Summary
- **Corpus Name**: Curated HSV Protein Sequence Dataset / Canonical HSV Protein Sequence Corpus
- **Total Canonical Sequences**: `22,689`
- **Supervised Temporal Dataset**: `16,657` sequences ([final_temporal_supervised_dataset.csv](file:///data/processed/final_temporal_supervised_dataset.csv))
- **Extended Curated Corpus**: `22,013` sequences ([final_curated_hsv_protein_corpus.csv](file:///data/processed/final_curated_hsv_protein_corpus.csv))
- **Viruses Represented**: Human alphaherpesvirus 1 (HSV-1, $N=11,802$) and Human alphaherpesvirus 2 (HSV-2, $N=10,263$).

## 2. Raw Sources & Deduplication
- **Raw Input Files**:
  - `raw/HSV2_non_redundant (2).fasta` ($10,540$ records)
  - `raw/non_redundant.fasta` ($12,170$ records)
- **Total Raw Records**: `22,710`
- **Exact Duplicate Collapsing**: $21$ cross-file identical sequence instances collapsed into deterministic SHA-256 identifiers (`CANONICAL_<hash>`).
- **Residue Modification**: Zero residues altered ($100.00\%$ exact SHA-256 match).

## 3. Biological Eligibility & Tiering
- **Tier A (`TIER_A_HIGH_CONFIDENCE_GROUND_TRUTH`)**: `10,482` sequences (Full-length reference & isolate proteins with Level 1 experimental temporal evidence).
- **Tier B (`TIER_B_VERIFIED_EXTENDED_BIOLOGICAL_CORPUS`)**: `11,531` sequences (Contains $9,980$ full-length sequences and $1,551$ partial/fragment sequences; $6,175$ labeled + $5,356$ uncharacterized).
- **Tier C (`TIER_C_CONFLICTING_QUARANTINED`)**: `5` sequences (Quarantined cross-file taxonomic conflicts).
- **Tier D (`TIER_D_EXCLUDED`)**: `671` sequences (Non-HSV hosts [411], technical PDB chains [147], recombinant mutants [50], synthetic patents [49], unassigned [14]).

## 4. Temporal Kinetic Classes
- **IMMEDIATE_EARLY ($lpha$)**: `1,552` sequences
- **EARLY ($eta$)**: `4,140` sequences
- **LATE ($\gamma_1 / \gamma_2$)**: `10,965` sequences
- **UNKNOWN**: `5,356` sequences in curated corpus ($6,027$ across all tiers)
- **CONFLICTING**: `5` sequences in manual review quarantine

## 5. Provenance & Authenticity Verification
- **Traceable Accessions**: $22,689 / 22,689$ ($100.00\%$) across NCBI GenBank ($22,292$), PDB ($146$), NCBI RefSeq ($100$), UniProtKB ($98$), Patent ($49$), and PRF ($4$).
- **Sequence Fidelity**: $22,689 / 22,689$ ($100.00\%$) exact SHA-256 match against source records.
- **Target-Species Verification**: $28$ core viral genes validated by independent HSV-2 experimental time-course literature.

## 6. Limitations & Appropriate Use
- Contains diverse laboratory strains and clinical isolate polymorphism variants.
- Supervised temporal classification tasks should utilize `final_temporal_supervised_dataset.csv` ($N=16,657$).
- Self-supervised representation learning and language model embedding tasks may utilize the entire `final_curated_hsv_protein_corpus.csv` ($N=22,013$).
