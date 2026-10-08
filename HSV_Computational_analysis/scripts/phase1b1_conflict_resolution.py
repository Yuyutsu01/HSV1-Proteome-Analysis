#!/usr/bin/env python3
"""
===============================================================================
Phase 1B.1: Final Manual Review of Biological Conflicts & Quarantine
===============================================================================
This module audits and resolves the 7 review-required records identified during
Phase 1B reconciliation using authoritative evidence rules.

Allowed Resolution Outcomes:
1. RESOLVED: Authoritative biological concordance between standard database records.
2. UNRESOLVED — QUARANTINED: Cross-file taxonomic/species conflict with no definitive external resolution.
3. EXCLUDED — DOCUMENTED_REASON: Degenerate/technical fragment unsuitable for proteome analysis.
===============================================================================
"""

import os
import sys
from pathlib import Path
from typing import Dict, List, Any
import pandas as pd


def resolve_conflicts(base_dir: Path):
    data_dir = base_dir / "data"
    results_dir = base_dir / "results"
    tables_dir = results_dir / "tables"
    logs_dir = results_dir / "logs"
    
    df_dataset = pd.read_csv(data_dir / "processed" / "phase1b_biological_dataset.csv")
    df_raw = pd.read_csv(data_dir / "processed" / "raw_derived_sequences.tsv", sep="\t")
    
    # 7 Review-required sequences
    review_ids = [
        "CANONICAL_01eb55d7fe80",
        "CANONICAL_08f271887ce9",
        "CANONICAL_40dd30c0e1cc",
        "CANONICAL_5082249f592d",
        "CANONICAL_a9996922e85c",
        "CANONICAL_e9476297d514",
        "CANONICAL_ec8d790b87fd"
    ]
    
    resolution_records = []
    
    for seq_id in review_ids:
        row = df_dataset[df_dataset["sequence_id"] == seq_id].iloc[0]
        raw_recs = df_raw[df_raw["sequence_sha256"].str.startswith(seq_id.replace("CANONICAL_", ""))]
        
        r1 = raw_recs.iloc[0]
        r2 = raw_recs.iloc[1] if len(raw_recs) > 1 else raw_recs.iloc[0]
        
        acc1 = r1["accession"]
        acc2 = r2["accession"]
        src1 = r1["raw_record_id"]
        src2 = r2["raw_record_id"]
        
        if seq_id == "CANONICAL_40dd30c0e1cc":
            # NCBI CAB06702.1 (UL54) and UniProt P28276 (ICP27) are exact synonyms for HSV-2 strain HG52 mRNA export factor
            res_status = "RESOLVED"
            ver_species = "HSV-2"
            ver_gene = "UL54 / ICP27"
            ver_protein = "mRNA export factor ICP27"
            ver_strain = "HG52"
            ev_source = "UniProtKB/Swiss-Prot P28276 and NCBI GenBank CAB06702.1"
            res_reason = "Authoritative biological concordance between UniProt (P28276) and NCBI (CAB06702.1) where UL54 and ICP27 are synonymous gene/protein names for HSV-2 HG52 mRNA export factor"
            
            # Update dataset row
            df_dataset.loc[df_dataset["sequence_id"] == seq_id, "virus_species"] = "HSV-2"
            df_dataset.loc[df_dataset["sequence_id"] == seq_id, "organism"] = "Human alphaherpesvirus 2"
            df_dataset.loc[df_dataset["sequence_id"] == seq_id, "strain"] = "HG52"
            df_dataset.loc[df_dataset["sequence_id"] == seq_id, "gene_symbol"] = "UL54"
            df_dataset.loc[df_dataset["sequence_id"] == seq_id, "protein_name"] = "mRNA export factor ICP27"
            df_dataset.loc[df_dataset["sequence_id"] == seq_id, "taxonomy_status"] = "VERIFIED"
            df_dataset.loc[df_dataset["sequence_id"] == seq_id, "protein_identity_status"] = "VERIFIED"
            df_dataset.loc[df_dataset["sequence_id"] == seq_id, "biological_eligibility"] = "ELIGIBLE_PRIMARY"
            df_dataset.loc[df_dataset["sequence_id"] == seq_id, "dataset_group"] = "GROUP A — VERIFIED HSV NATURAL PROTEINS"
            df_dataset.loc[df_dataset["sequence_id"] == seq_id, "manual_review_required"] = False
            df_dataset.loc[df_dataset["sequence_id"] == seq_id, "review_reason"] = "RESOLVED: Authoritative synonym concordance (UniProt P28276 / NCBI CAB06702.1)"
            
        elif seq_id == "CANONICAL_08f271887ce9":
            # Length = 1 amino acid (M) with conflicting gE/US6 and HSV-1/HSV-2 labels
            res_status = "EXCLUDED — DOCUMENTED_REASON"
            ver_species = "UNKNOWN"
            ver_gene = "UNKNOWN"
            ver_protein = "UNRESOLVED_FRAGMENT"
            ver_strain = "UNKNOWN"
            ev_source = "Sequence property analysis (Length=1 residue)"
            res_reason = "Degenerate single-amino-acid fragment (length=1) with conflicting cross-file taxonomy and gene annotations; excluded from proteome analysis"
            
            # Update dataset row
            df_dataset.loc[df_dataset["sequence_id"] == seq_id, "biological_eligibility"] = "EXCLUDE_TECHNICAL"
            df_dataset.loc[df_dataset["sequence_id"] == seq_id, "eligibility_reason"] = "Degenerate 1-residue fragment with contradictory annotations"
            df_dataset.loc[df_dataset["sequence_id"] == seq_id, "dataset_group"] = "GROUP F — ENGINEERED / SYNTHETIC / TECHNICAL"
            df_dataset.loc[df_dataset["sequence_id"] == seq_id, "manual_review_required"] = False
            df_dataset.loc[df_dataset["sequence_id"] == seq_id, "review_reason"] = "EXCLUDED: 1-aa degenerate fragment"
            
        else:
            # 5 Cross-file taxonomy conflicts without independent authoritative resolution: Quarantined
            res_status = "UNRESOLVED — QUARANTINED"
            ver_species = "UNKNOWN"
            ver_gene = r1.get("gene_symbol", "UNKNOWN") if r1.get("gene_symbol") == r2.get("gene_symbol") else "CONFLICTING"
            ver_protein = r1.get("protein_name", "UNKNOWN")
            ver_strain = "UNKNOWN"
            ev_source = "FASTA headers (Conflicting species: HSV-1 vs HSV-2)"
            res_reason = "AUTHORITATIVE RESOLUTION UNAVAILABLE — Conflicting species headers between source datasets without independent ground truth; quarantined in review queue"
            
            # Maintain explicit quarantine in dataset
            df_dataset.loc[df_dataset["sequence_id"] == seq_id, "biological_eligibility"] = "REVIEW_REQUIRED"
            df_dataset.loc[df_dataset["sequence_id"] == seq_id, "dataset_group"] = "GROUP G — UNRESOLVED / CONFLICTING"
            df_dataset.loc[df_dataset["sequence_id"] == seq_id, "manual_review_required"] = True
            df_dataset.loc[df_dataset["sequence_id"] == seq_id, "review_reason"] = "QUARANTINED: Cross-file taxonomy conflict (HSV-1 vs HSV-2)"
            
        resolution_records.append({
            "sequence_id": seq_id,
            "accession_1": acc1,
            "accession_2": acc2,
            "source_1": src1,
            "source_2": src2,
            "verified_species": ver_species,
            "verified_gene": ver_gene,
            "verified_protein": ver_protein,
            "verified_strain": ver_strain,
            "evidence_source": ev_source,
            "resolution_status": res_status,
            "resolution_reason": res_reason
        })
        
    df_resolutions = pd.DataFrame(resolution_records)
    df_resolutions.to_csv(tables_dir / "phase1b_final_conflict_resolution.csv", index=False)
    
    # Save updated dataset
    df_dataset.to_csv(data_dir / "processed" / "phase1b_biological_dataset.csv", index=False)
    
    # Update exclusion audit table
    df_exclusions = pd.read_csv(tables_dir / "phase1b_exclusion_audit.csv")
    if "CANONICAL_08f271887ce9" not in df_exclusions["sequence_id"].values:
        new_excl = pd.DataFrame([{
            "sequence_id": "CANONICAL_08f271887ce9",
            "original_record_ids": "HSV2_non_redundant (2).fasta#rec_009804;non_redundant.fasta#rec_009807",
            "exclusion_category": "EXCLUDE_TECHNICAL",
            "exclusion_reason": "Degenerate 1-residue fragment with contradictory annotations",
            "evidence_source": "Sequence property analysis (Length=1 residue)"
        }])
        df_exclusions = pd.concat([df_exclusions, new_excl], ignore_index=True)
        df_exclusions.to_csv(tables_dir / "phase1b_exclusion_audit.csv", index=False)
    
    # Recalculate reconciliation totals
    ct_group = pd.crosstab(df_dataset["dataset_group"], df_dataset["biological_eligibility"], margins=True, margins_name="Total")
    counts_elig = df_dataset["biological_eligibility"].value_counts().to_dict()
    counts_species = df_dataset["virus_species"].value_counts().to_dict()
    
    report_path = logs_dir / "phase1b1_final_resolution_report.txt"
    with open(report_path, "w", encoding="utf-8") as f:
        f.write("=" * 80 + "\n")
        f.write("PHASE 1B.1: FINAL CONFLICT RESOLUTION & RECONCILIATION REPORT\n")
        f.write("=" * 80 + "\n\n")
        
        f.write("1. RESOLUTION STATUS FOR THE 7 REVIEW-REQUIRED RECORDS\n")
        f.write("-" * 80 + "\n")
        f.write(df_resolutions[["sequence_id", "verified_species", "verified_gene", "resolution_status", "resolution_reason"]].to_string(index=False) + "\n\n")
        
        f.write("2. RECONCILED BIOLOGICAL ELIGIBILITY TOTALS\n")
        f.write("-" * 80 + "\n")
        for k, v in counts_elig.items():
            f.write(f"  {k}: {v:,}\n")
        f.write(f"  Total: {len(df_dataset):,}\n\n")
        
        f.write("3. RECONCILED DATASET GROUP x ELIGIBILITY CROSS-TABULATION\n")
        f.write("-" * 80 + "\n")
        f.write(ct_group.to_string() + "\n\n")
        
        f.write("=" * 80 + "\n")
        f.write("FINAL STATUS: PHASE 1B.1 COMPLETE — ALL CONFLICTS RESOLVED OR QUARANTINED\n")
        f.write("=" * 80 + "\n")
        
    print(f"Phase 1B.1 completed successfully! Report saved to {report_path}")
    return counts_elig


if __name__ == "__main__":
    base_p = Path(__file__).resolve().parent.parent
    res = resolve_conflicts(base_p)
    print(res)
