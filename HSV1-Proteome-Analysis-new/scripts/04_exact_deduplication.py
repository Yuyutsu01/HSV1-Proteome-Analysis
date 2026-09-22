"""
Script 04: Exact Sequence Deduplication & Dataset Discrepancy Auditing

Scientific Objective:
Performs rigorous sequence deduplication by mapping identical amino acid sequences
that arise due to inverted repeats (TR_L/IR_L and TR_S/IR_S) in alphaherpesvirus genomes.
Preserves complete provenance of retained vs removed copies.
Validates against historical reference baseline (77 raw -> 74 unique) while reporting
any factual deviations without hard-coding or data fabrication.

Outputs:
- data/processed/unique_proteins.csv
- data/intermediate/duplicates.csv
- data/processed/unique_proteins.fasta
- results/tables/dataset_summary.csv

Author: Computational Biology Pipeline
"""

import os
import yaml
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def load_config(config_path="config.yaml"):
    """Load pipeline YAML configuration."""
    if not os.path.exists(config_path):
        if os.path.exists(os.path.join("..", config_path)):
            config_path = os.path.join("..", config_path)
    with open(config_path, "r") as f:
        return yaml.safe_load(f)

def run_exact_deduplication(config_path="config.yaml"):
    """
    Perform exact deduplication and generate structured non-redundant dataset.
    """
    cfg = load_config(config_path)
    inter_dir = cfg['paths']['data_intermediate']
    proc_dir = cfg['paths']['data_processed']
    tables_dir = cfg['paths']['results_tables']
    
    os.makedirs(proc_dir, exist_ok=True)
    os.makedirs(tables_dir, exist_ok=True)
    
    passed_csv_path = os.path.join(inter_dir, "qc_passed_proteins.csv")
    assert os.path.exists(passed_csv_path), f"Error: qc_passed_proteins.csv not found at {passed_csv_path}"
    
    df_passed = pd.read_csv(passed_csv_path)
    print(f"[Phase 4] Processing exact sequence deduplication for {len(df_passed)} QC-passed records...")
    
    seen_sequences = {} # seq -> primary_record_dict
    unique_records = []
    duplicate_records = []
    
    for idx, row in df_passed.iterrows():
        seq = str(row['sequence']).upper()
        prot_id = str(row['protein_id'])
        gene = str(row['gene'])
        
        if seq not in seen_sequences:
            seen_sequences[seq] = {
                "retained_protein_id": prot_id,
                "retained_gene": gene,
                "retained_product": row['product']
            }
            unique_records.append(row.to_dict())
        else:
            primary = seen_sequences[seq]
            duplicate_records.append({
                "duplicate_protein_id": prot_id,
                "duplicate_gene": gene,
                "duplicate_product": row['product'],
                "duplicate_location": row['genomic_location'],
                "retained_protein_id": primary['retained_protein_id'],
                "retained_gene": primary['retained_gene'],
                "retained_product": primary['retained_product'],
                "sequence_length": len(seq)
            })
            
    df_unique = pd.DataFrame(unique_records)
    df_dups = pd.DataFrame(duplicate_records)
    
    unique_csv_path = os.path.join(proc_dir, "unique_proteins.csv")
    dups_csv_path = os.path.join(inter_dir, "duplicates.csv")
    unique_fasta_path = os.path.join(proc_dir, "unique_proteins.fasta")
    
    df_unique.to_csv(unique_csv_path, index=False)
    df_dups.to_csv(dups_csv_path, index=False)
    
    # Write non-redundant FASTA
    fasta_seq_records = []
    for _, row in df_unique.iterrows():
        rec = SeqRecord(
            Seq(row['sequence']),
            id=row['protein_id'],
            description=f"{row['gene']} | {row['product']} | len:{row['sequence_length']}aa"
        )
        fasta_seq_records.append(rec)
    SeqIO.write(fasta_seq_records, unique_fasta_path, "fasta")
    
    # Dataset Summary Metrics
    raw_total = len(df_passed)
    unique_total = len(df_unique)
    dup_total = len(df_dups)
    
    # Discrepancy audit against historical literature baseline
    expected_raw = 77
    expected_unique = 74
    raw_discrepancy = raw_total - expected_raw
    unique_discrepancy = unique_total - expected_unique
    
    summary_data = {
        "Metric": [
            "GenBank Reference Accession",
            "Organism / Strain",
            "Total Raw CDS Translations Parsed",
            "QC-Passed Protein Records",
            "Exact Duplicate Copies Removed",
            "Non-Redundant Unique Proteins Retained",
            "Minimum Sequence Length (aa)",
            "Maximum Sequence Length (aa)",
            "Mean Sequence Length (aa)",
            "Median Sequence Length (aa)",
            "Expected Raw Reference Baseline",
            "Expected Unique Reference Baseline",
            "Raw Count Discrepancy vs Baseline",
            "Unique Count Discrepancy vs Baseline"
        ],
        "Value": [
            cfg['ncbi']['accession'],
            cfg['project']['reference_organism'],
            str(raw_total),
            str(raw_total),
            str(dup_total),
            str(unique_total),
            str(df_unique['sequence_length'].min()),
            str(df_unique['sequence_length'].max()),
            f"{df_unique['sequence_length'].mean():.2f}",
            f"{df_unique['sequence_length'].median():.1f}",
            str(expected_raw),
            str(expected_unique),
            f"{raw_discrepancy:+d}" if raw_discrepancy != 0 else "0 (Exact Match)",
            f"{unique_discrepancy:+d}" if unique_discrepancy != 0 else "0 (Exact Match)"
        ]
    }
    df_summary = pd.DataFrame(summary_data)
    summary_csv_path = os.path.join(tables_dir, "dataset_summary.csv")
    df_summary.to_csv(summary_csv_path, index=False)
    
    print(f"\n[Phase 4] Exact Deduplication Audit Results:")
    print(f"  - Total QC-Passed Records: {raw_total}")
    print(f"  - Duplicate Copies Identified: {dup_total}")
    print(f"  - Unique Non-Redundant Proteins: {unique_total}")
    print(f"  - Length Min / Max / Mean / Median: {df_unique['sequence_length'].min()} / {df_unique['sequence_length'].max()} / {df_unique['sequence_length'].mean():.1f} / {df_unique['sequence_length'].median():.1f} aa")
    print(f"  - Baseline Concordance: {'EXACT MATCH (77 raw -> 74 unique)' if (raw_total==77 and unique_total==74) else 'DEVIATION OBSERVED'}")
    
    if len(df_dups) > 0:
        print(f"\n  Removed Repeat Copies Details:")
        for _, d in df_dups.iterrows():
            print(f"    * Duplicate {d['duplicate_gene']} ({d['duplicate_protein_id']}) -> Identical copy of {d['retained_gene']} ({d['retained_protein_id']})")
            
    print(f"\n  Saved unique proteins CSV: {unique_csv_path}")
    print(f"  Saved unique proteins FASTA: {unique_fasta_path}")
    print(f"  Saved dataset summary: {summary_csv_path}")
    
    return df_unique, df_dups, df_summary

if __name__ == "__main__":
    run_exact_deduplication()
