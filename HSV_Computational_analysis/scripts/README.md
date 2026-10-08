# Scripts Directory: Pipeline Reproducibility Index

This directory contains deterministic, standalone scripts required to reproduce
the entire biological preprocessing, eligibility curation, temporal expression
annotation, and dataset validation pipeline.

| Script | Phase | Purpose | Inputs | Outputs | Status |
| :--- | :--- | :--- | :--- | :--- | :--- |
| `phase1_ingestion.py` | Phase 1A | Raw FASTA ingestion, validation, and SHA-256 deduplication | `raw/*.fasta` | `unique_high_quality_sequences.tsv` | Frozen |
| `phase1_audit.py` | Phase 1A | Computational QC and non-canonical residue composition audit | `unique_high_quality_sequences.tsv` | `results/reports/phase1_*` | Frozen |
| `phase1b_curation.py` | Phase 1B | Biological eligibility, source-type curation, and host filtering | `unique_high_quality_sequences.tsv` | `phase1b_biological_dataset.csv` | Frozen |
| `phase1b_reconciliation.py` | Phase 1B | Biological eligibility reconciliation and cross-tabulation | `phase1b_biological_dataset.csv` | `results/tables/phase1b_*` | Frozen |
| `phase1b1_conflict_resolution.py` | Phase 1B.1 | Manual review and quarantine of cross-file taxonomic conflicts | `phase1b_conflicts.csv` | `phase1b_final_conflict_resolution.csv` | Frozen |
| `phase2_annotation.py` | Phase 2 | Evidence-grounded biological annotation & kinetic curation | `phase1b_biological_dataset.csv` | `phase2_biological_annotation.csv` | Frozen |
| `phase2_1_reconciliation.py` | Phase 2.1 | Decoupling temporal annotation from supervised ground truth | `phase2_biological_annotation.csv` | `phase2_1_reconciled_annotation.csv` | Frozen |
| `phase2_1_1_provenance_audit.py` | Phase 2.1.1 | Temporal evidence provenance and transfer detection audit | `phase2_1_reconciled_annotation.csv` | `phase2_1_1_final_annotation.csv` | Frozen |
| `phase2_1_2_cross_species_audit.py` | Phase 2.1.2 | Cross-species orthology transfer audit | `phase2_1_1_final_annotation.csv` | `phase2_1_2_final_annotation.csv` | Frozen |
| `phase2_1_2A_target_species_audit.py` | Phase 2.1.2A | Evidence-level target-species temporal verification matrix | `phase2_1_2_final_annotation.csv` | `phase2_1_2A_final_annotation.csv` | Frozen |
| `phase2_final_dataset_freeze.py` | Phase 2 Freeze | Final supervised temporal dataset extraction & freeze | `phase2_1_2A_final_annotation.csv` | `final_temporal_supervised_dataset.csv` | Frozen |
| `phase2_5_sequence_authenticity_audit.py` | Phase 2.5 | Sequence authenticity, provenance, and dataset integrity | `raw/*.fasta`, `phase2_1_2A_final_annotation.csv` | `phase2_5_*` | Frozen |
| `phase2_cleanup_and_freeze.py` | Pre-Phase-3 | Repository structure cleanup, inventory, and manifest freeze | All workspace artifacts | `PROJECT_MANIFEST.csv`, `PROJECT_STATUS.md` | Frozen |
