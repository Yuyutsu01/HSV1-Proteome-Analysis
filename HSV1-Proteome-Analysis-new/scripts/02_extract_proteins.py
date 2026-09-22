"""
Script 02: CDS Feature Extraction & Protein Sequence Parsing

Scientific Objective:
Parse the official RefSeq GenBank record (NC_001806.2) using BioPython SeqIO.
Extract all Coding Sequences (CDS) with annotated translations, recording complete
genomic metadata (locus tag, gene name, product description, coordinate spans, strand).
Outputs structured CSV and FASTA files without manual intervention.

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

def extract_proteins(config_path="config.yaml"):
    """
    Extract all translated CDS records from the downloaded GenBank file.
    """
    cfg = load_config(config_path)
    accession = cfg['ncbi']['accession']
    gbk_path = os.path.join(cfg['paths']['data_raw_ncbi'], f"{accession}.gbk")
    fasta_out_dir = cfg['paths']['data_raw_fasta']
    inter_out_dir = cfg['paths']['data_intermediate']
    
    os.makedirs(fasta_out_dir, exist_ok=True)
    os.makedirs(inter_out_dir, exist_ok=True)
    
    fasta_out_path = os.path.join(fasta_out_dir, "raw_proteins.fasta")
    csv_out_path = os.path.join(inter_out_dir, "raw_proteins.csv")
    
    print(f"[Phase 2] Parsing GenBank record: {gbk_path}")
    assert os.path.exists(gbk_path), f"Error: GenBank file not found at {gbk_path}"
    
    record = SeqIO.read(gbk_path, "genbank")
    genome_accession = record.id
    genome_strain = record.annotations.get("source", "Human herpesvirus 1 strain 17")
    
    extracted_records = []
    fasta_records = []
    
    cds_index = 0
    for feature in record.features:
        if feature.type == "CDS":
            qualifiers = feature.qualifiers
            translation = qualifiers.get("translation", [None])[0]
            if not translation:
                continue
            
            cds_index += 1
            protein_id = qualifiers.get("protein_id", [f"CDS_{cds_index:03d}"])[0]
            gene = qualifiers.get("gene", qualifiers.get("locus_tag", [f"CDS_{cds_index:03d}"]))[0]
            locus_tag = qualifiers.get("locus_tag", [""])[0]
            product = qualifiers.get("product", ["Uncharacterized protein"])[0]
            protein_name = product
            
            location_str = str(feature.location)
            strand = "+" if feature.location.strand == 1 else ("-" if feature.location.strand == -1 else ".")
            start_pos = int(feature.location.start)
            end_pos = int(feature.location.end)
            
            clean_seq = translation.strip().upper()
            seq_len = len(clean_seq)
            
            row_dict = {
                "protein_id": protein_id,
                "gene": gene,
                "locus_tag": locus_tag,
                "protein_name": protein_name,
                "product": product,
                "sequence": clean_seq,
                "sequence_length": seq_len,
                "genome_accession": genome_accession,
                "genome_strain": genome_strain,
                "start_pos": start_pos,
                "end_pos": end_pos,
                "strand": strand,
                "genomic_location": location_str
            }
            extracted_records.append(row_dict)
            
            # Create SeqRecord for FASTA output
            desc = f"{gene} | {product} | loc:{start_pos}..{end_pos}({strand})"
            fasta_rec = SeqRecord(Seq(clean_seq), id=protein_id, description=desc)
            fasta_records.append(fasta_rec)
            
    df_raw = pd.DataFrame(extracted_records)
    df_raw.to_csv(csv_out_path, index=False)
    SeqIO.write(fasta_records, fasta_out_path, "fasta")
    
    print(f"[Phase 2] Successfully extracted {len(df_raw)} translated CDS records.")
    print(f"  - Saved structured CSV: {csv_out_path}")
    print(f"  - Saved FASTA: {fasta_out_path}")
    
    return df_raw

if __name__ == "__main__":
    extract_proteins()
