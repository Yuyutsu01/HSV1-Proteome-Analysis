#!/usr/bin/env python3
"""
===============================================================================
Phase 1B: Biological Eligibility & Source-Type Curation Pipeline
===============================================================================
This module executes rigorous, evidence-grounded biological curation for the
22,689 unique sequences derived from Phase 1A.

Strict Rules Followed:
1. No temporal labels (IE, Early, Late) or machine learning inferences.
2. No biological hallucination: Unresolved fields remain UNKNOWN/UNRESOLVED.
3. Strict taxonomy evidence tracking (Header, UniProt, RefSeq, GenBank).
4. Explicit conflict detection and manual review queuing.
5. Surjective provenance preservation across all intermediate and processed tables.
===============================================================================
"""

import os
import sys
import re
import hashlib
import json
from pathlib import Path
from typing import Dict, List, Tuple, Any, Optional
from collections import Counter
import pandas as pd


# Canonical Amino Acid sets
CANONICAL_20 = set("ACDEFGHIKLMNPQRSTVWY")
AMBIGUOUS_CODES = set("BZXJ")
NON_STANDARD_CODES = set("UO")
STOP_CODES = set("*")


def extract_header_fields(raw_header: str) -> Dict[str, str]:
    """
    Extracts explicit fields from a FASTA header line.
    """
    hdr = raw_header.lstrip(">").strip()
    res = {
        "raw_header": hdr,
        "accession": "UNKNOWN",
        "entry_name": "UNKNOWN",
        "protein_name": "UNKNOWN",
        "gene_symbol": "UNKNOWN",
        "organism": "UNRESOLVED",
        "strain": "UNKNOWN",
        "isolate": "UNKNOWN",
        "source_db": "UNKNOWN",
        "is_patent": False,
        "is_pdb": False,
        "is_uniprot": False,
        "is_prf": False,
        "is_ncbi": False,
    }
    
    if not hdr:
        return res

    # Check for organism in brackets e.g. [Human alphaherpesvirus 1]
    org_match = re.search(r"\[(.*?)\]", hdr)
    if org_match:
        res["organism"] = org_match.group(1).strip()
        desc_text = hdr[:org_match.start()] + hdr[org_match.end():]
    else:
        desc_text = hdr

    # Check for strain / isolate in organism or description
    strain_match = re.search(r"\bstrain\s+([A-Za-z0-9_\-\./]+)", hdr, re.IGNORECASE)
    if strain_match:
        res["strain"] = strain_match.group(1).strip()
        
    isolate_match = re.search(r"\bisolate\s+([A-Za-z0-9_\-\./]+)", hdr, re.IGNORECASE)
    if isolate_match:
        res["isolate"] = isolate_match.group(1).strip()

    # Parse database specific formats
    if hdr.startswith("pdb|") or "|pdb|" in hdr:
        res["is_pdb"] = True
        res["source_db"] = "PDB"
        parts = hdr.split("|")
        if len(parts) >= 3:
            res["accession"] = parts[1].strip()
            res["protein_name"] = parts[2].strip()
    elif hdr.startswith("sp|") or hdr.startswith("tr|"):
        res["is_uniprot"] = True
        res["source_db"] = "UniProtKB"
        parts = hdr.split("|")
        if len(parts) >= 3:
            res["accession"] = parts[1].strip()
            rest = parts[2].strip().split(maxsplit=1)
            res["entry_name"] = rest[0].strip()
            if len(rest) > 1:
                res["protein_name"] = rest[1].strip()
    elif hdr.startswith("prf||"):
        res["is_prf"] = True
        res["source_db"] = "PRF"
        parts = hdr.split("||")
        if len(parts) >= 2:
            tokens = parts[1].strip().split(maxsplit=1)
            res["accession"] = tokens[0].strip()
            if len(tokens) > 1:
                res["protein_name"] = tokens[1].strip()
    elif "patent" in hdr.lower():
        res["is_patent"] = True
        res["source_db"] = "PATENT"
        tokens = hdr.split(maxsplit=1)
        res["accession"] = tokens[0].strip()
        if len(tokens) > 1:
            res["protein_name"] = tokens[1].strip()
    else:
        # Standard NCBI / GenBank header: >ACCESSION Description
        res["is_ncbi"] = True
        res["source_db"] = "NCBI_GenBank"
        tokens = desc_text.strip().split(maxsplit=1)
        res["accession"] = tokens[0].strip()
        if len(tokens) > 1:
            res["protein_name"] = tokens[1].strip()

    # Extract Gene Symbol if present (e.g. UL33, US6, RL2, RS1, ICP4, TK, etc.)
    gene_match = re.search(r"\b(UL\d+[A-Za-z]?|US\d+[A-Za-z]?|RL\d+[A-Za-z]?|RS\d+[A-Za-z]?|ICP\d+|g[B-N])\b", hdr, re.IGNORECASE)
    if gene_match:
        res["gene_symbol"] = gene_match.group(1).upper()
        
    return res


def determine_taxonomy(
    headers_info: List[Dict[str, str]]
) -> Tuple[str, str, str, str, bool]:
    """
    Curates taxonomy using strictly explicit evidence from raw headers.
    Returns: (organism, virus_species, taxonomy_status, taxonomy_evidence, has_conflict)
    """
    organisms = list({h["organism"] for h in headers_info if h["organism"] != "UNRESOLVED"})
    
    # Check for conflicting organisms across multi-header duplicates
    if len(organisms) > 1:
        # Check if they belong to different species
        hsv1_terms = ["human alphaherpesvirus 1", "human herpesvirus 1", "synthetic human alphaherpesvirus 1"]
        hsv2_terms = ["human alphaherpesvirus 2", "human herpesvirus 2"]
        
        has_hsv1 = any(any(t in o.lower() for t in hsv1_terms) for o in organisms)
        has_hsv2 = any(any(t in o.lower() for t in hsv2_terms) for o in organisms)
        
        if has_hsv1 and has_hsv2:
            return (
                " / ".join(organisms),
                "UNKNOWN",
                "CONFLICTING",
                f"Conflicting taxonomy between raw headers: {'; '.join(organisms)}",
                True
            )
            
    if organisms:
        org_str = organisms[0]
        org_lower = org_str.lower()
        
        # Explicit HSV-1
        if "human alphaherpesvirus 1" in org_lower or "human herpesvirus 1" in org_lower:
            return org_str, "HSV-1", "VERIFIED", f"Explicit FASTA header [{org_str}]", False
            
        # Explicit HSV-2
        if "human alphaherpesvirus 2" in org_lower or "human herpesvirus 2" in org_lower:
            return org_str, "HSV-2", "VERIFIED", f"Explicit FASTA header [{org_str}]", False
            
        # Other Herpesviruses
        other_herpes = [
            "betaherpesvirus", "alphaherpesvirus", "gammaherpesvirus",
            "bovine alphaherpesvirus", "gallid alphaherpesvirus", "equid alphaherpesvirus",
            "anatid alphaherpesvirus", "macacine alphaherpesvirus"
        ]
        if any(h in org_lower for h in other_herpes):
            return org_str, "OTHER_HERPESVIRUS", "VERIFIED", f"Explicit FASTA header [{org_str}]", False
            
        # Non-Herpesvirus (Hosts, other virus families, vectors)
        return org_str, "NON_HERPESVIRUS", "VERIFIED", f"Explicit FASTA header [{org_str}]", False
        
    # If organism not bracketed in header, check UniProt mnemonic or PDB
    for h in headers_info:
        if h["is_uniprot"]:
            entry = h["entry_name"]
            if "_HHV1" in entry or "_HHV11" in entry or "_HHV1E" in entry:
                return "Human alphaherpesvirus 1", "HSV-1", "SUPPORTED", f"UniProt entry name mnemonic {entry}", False
            if "_HHV2" in entry or "_HHV2H" in entry:
                return "Human alphaherpesvirus 2", "HSV-2", "SUPPORTED", f"UniProt entry name mnemonic {entry}", False
            if "_HUMAN" in entry:
                return "Homo sapiens", "NON_HERPESVIRUS", "SUPPORTED", f"UniProt host mnemonic {entry}", False
            if "_MOUSE" in entry:
                return "Mus musculus", "NON_HERPESVIRUS", "SUPPORTED", f"UniProt host mnemonic {entry}", False

    return "UNRESOLVED", "UNKNOWN", "UNRESOLVED", "No explicit taxonomy in FASTA header", False


def determine_source_type_and_completeness(
    headers_info: List[Dict[str, str]],
    organism: str,
    virus_species: str
) -> Tuple[str, str, str, str, str]:
    """
    Classifies source_type, completeness_status, natural_or_engineered,
    source_type_evidence, and source_type_status.
    """
    raw_hdrs_lower = [h["raw_header"].lower() for h in headers_info]
    combined_hdrs = " ".join(raw_hdrs_lower)
    
    # 1. Patent sequences
    if any(h["is_patent"] for h in headers_info) or "patent" in combined_hdrs:
        return (
            "SYNTHETIC",
            "UNKNOWN",
            "SYNTHETIC",
            "Header explicitly declares patent sequence",
            "VERIFIED"
        )
        
    # 2. PDB structural derived
    if any(h["is_pdb"] for h in headers_info) or "pdb|" in combined_hdrs or "chain " in combined_hdrs:
        return (
            "PDB_DERIVED",
            "PDB_FRAGMENT",
            "NATURAL_OR_RECOMBINANT_CONSTRUCT",
            "Header explicitly indicates PDB coordinate chain",
            "VERIFIED"
        )
        
    # 3. Engineered / Recombinant vectors / Mutants
    recomb_patterns = r"\b(recombinant|mutant|mutated|engineered|expression vector|cloning vector|chimera)\b"
    if re.search(recomb_patterns, combined_hdrs) or "synthetic human alphaherpesvirus" in combined_hdrs:
        return (
            "ENGINEERED_RECOMBINANT",
            "ENGINEERED_FRAGMENT",
            "ENGINEERED",
            "Header explicitly declares recombinant/mutant/vector construct",
            "VERIFIED"
        )
        
    # 4. Partial / Fragments
    fragment_patterns = r"\b(partial|fragment|truncated|subsequence)\b"
    if re.search(fragment_patterns, combined_hdrs):
        return (
            "NATURAL_VIRAL_FRAGMENT",
            "PARTIAL_FRAGMENT",
            "NATURAL",
            "Header explicitly declares partial/fragment status",
            "VERIFIED"
        )
        
    # 5. Natural Viral / Host Proteins
    if virus_species in ["HSV-1", "HSV-2", "OTHER_HERPESVIRUS", "NON_HERPESVIRUS"] and organism != "UNRESOLVED":
        return (
            "NATURAL_VIRAL_PROTEIN" if virus_species in ["HSV-1", "HSV-2", "OTHER_HERPESVIRUS"] else "OTHER_TECHNICAL",
            "FULL_LENGTH",
            "NATURAL",
            f"Standard natural isolate record from {organism}",
            "SUPPORTED"
        )
        
    return (
        "UNKNOWN",
        "UNKNOWN",
        "UNKNOWN",
        "Insufficient explicit metadata to confirm source type",
        "UNRESOLVED"
    )


def determine_sequence_composition(seq: str) -> Tuple[str, int, int]:
    """
    Analyzes exact amino acid composition.
    Returns: (sequence_composition_status, ambiguous_count, noncanonical_count)
    """
    chars = set(seq.strip().upper())
    ambig_count = sum(seq.count(c) for c in chars.intersection(AMBIGUOUS_CODES))
    non_std_count = sum(seq.count(c) for c in chars.intersection(NON_STANDARD_CODES))
    stop_count = sum(seq.count(c) for c in chars.intersection(STOP_CODES))
    non_canon_total = ambig_count + non_std_count + stop_count
    
    distinct_non_canon = len(chars - CANONICAL_20)
    
    if distinct_non_canon > 1:
        status = "MULTIPLE_NON_CANONICAL_TYPES"
    elif stop_count > 0:
        status = "STOP_SYMBOL_PRESENT"
    elif non_std_count > 0:
        status = "NON_STANDARD_RESIDUES"
    elif ambig_count > 0:
        status = "AMBIGUOUS_RESIDUES"
    else:
        status = "CANONICAL_20AA"
        
    return status, ambig_count, non_canon_total


def determine_biological_eligibility(
    virus_species: str,
    source_type: str,
    completeness: str,
    taxonomy_status: str,
    protein_id_status: str,
    has_conflict: bool,
    composition_status: str
) -> Tuple[str, str, bool, str, str]:
    """
    Evaluates biological eligibility and category group according to Part H and I.
    """
    # 1. Conflict handling -> REVIEW_REQUIRED
    if has_conflict or taxonomy_status == "CONFLICTING" or protein_id_status == "CONFLICTING":
        return (
            "REVIEW_REQUIRED",
            "Conflicting biological metadata across duplicate source records",
            True,
            "Conflicting taxonomy or protein identity between source files",
            "GROUP G — UNRESOLVED / CONFLICTING"
        )
        
    # 2. Exclude Non-HSV
    if virus_species == "NON_HERPESVIRUS":
        return (
            "EXCLUDE_NON_HSV",
            "Non-HSV biological organism (Host/other virus/vector)",
            False,
            "NONE",
            "GROUP E — NON-HERPESVIRUS"
        )
        
    # 3. Other Herpesviruses (Retain in Group D)
    if virus_species == "OTHER_HERPESVIRUS":
        return (
            "EXCLUDE_NON_HSV",
            "Non-target herpesvirus species (Retained in Group D)",
            False,
            "NONE",
            "GROUP D — OTHER HERPESVIRUSES"
        )
        
    # 4. Exclude Engineered / Synthetic / Technical
    if source_type == "SYNTHETIC":
        return (
            "EXCLUDE_SYNTHETIC",
            "Synthetic or patent sequence",
            False,
            "NONE",
            "GROUP F — ENGINEERED / SYNTHETIC / TECHNICAL"
        )
    if source_type == "ENGINEERED_RECOMBINANT":
        return (
            "EXCLUDE_ENGINEERED",
            "Engineered recombinant construct or cloning vector",
            False,
            "NONE",
            "GROUP F — ENGINEERED / SYNTHETIC / TECHNICAL"
        )
    if source_type == "PDB_DERIVED":
        return (
            "EXCLUDE_TECHNICAL",
            "PDB coordinate structure chain / crystallization construct",
            False,
            "NONE",
            "GROUP F — ENGINEERED / SYNTHETIC / TECHNICAL"
        )
        
    # 5. Unresolved Identity
    if taxonomy_status == "UNRESOLVED" or virus_species == "UNKNOWN":
        return (
            "EXCLUDE_UNRESOLVED_IDENTITY",
            "Unresolved biological taxonomy and provenance",
            True,
            "Taxonomy cannot be established from explicit source evidence",
            "GROUP G — UNRESOLVED / CONFLICTING"
        )
        
    # 6. Natural HSV Fragments
    if source_type == "NATURAL_VIRAL_FRAGMENT" or completeness == "PARTIAL_FRAGMENT":
        return (
            "ELIGIBLE_SEQUENCE_ANALYSIS_ONLY",
            "Verified natural HSV fragment; eligible for sequence-level analyses",
            False,
            "NONE",
            "GROUP C — HSV FRAGMENTS"
        )
        
    # 7. Sequences with Ambiguous residues (retained for sequence analysis)
    if composition_status in ["AMBIGUOUS_RESIDUES", "MULTIPLE_NON_CANONICAL_TYPES"]:
        return (
            "ELIGIBLE_SEQUENCE_ANALYSIS_ONLY",
            "Natural HSV protein containing ambiguous IUPAC residues (X/B/Z/J)",
            False,
            "NONE",
            "GROUP B — HSV NATURAL PROTEINS WITH INCOMPLETE METADATA"
        )
        
    # 8. Primary Eligible Natural Full-Length Proteins
    if virus_species in ["HSV-1", "HSV-2"] and source_type == "NATURAL_VIRAL_PROTEIN":
        if protein_id_status == "VERIFIED" and taxonomy_status == "VERIFIED":
            return (
                "ELIGIBLE_PRIMARY",
                "Verified natural full-length HSV protein meeting all biological criteria",
                False,
                "NONE",
                "GROUP A — VERIFIED HSV NATURAL PROTEINS"
            )
        else:
            return (
                "ELIGIBLE_PRIMARY",
                "Supported natural full-length HSV protein",
                False,
                "NONE",
                "GROUP B — HSV NATURAL PROTEINS WITH INCOMPLETE METADATA"
            )
            
    return (
        "REVIEW_REQUIRED",
        "Indeterminate eligibility combination",
        True,
        "Requires manual biological review",
        "GROUP G — UNRESOLVED / CONFLICTING"
    )


def run_phase1b_curation(base_dir: Path) -> Dict[str, Any]:
    """
    Executes complete Phase 1B biological eligibility and source-type curation.
    """
    data_dir = base_dir / "data"
    results_dir = base_dir / "results"
    
    unique_tsv_path = data_dir / "processed" / "unique_high_quality_sequences.tsv"
    raw_tsv_path = data_dir / "processed" / "raw_derived_sequences.tsv"
    
    if not unique_tsv_path.exists() or not raw_tsv_path.exists():
        raise FileNotFoundError("Phase 1A processed tables missing. Run Phase 1A ingestion first.")
        
    df_unique = pd.read_csv(unique_tsv_path, sep="\t")
    df_raw = pd.read_csv(raw_tsv_path, sep="\t")
    
    print(f"Loaded {len(df_unique):,} unique sequences and {len(df_raw):,} raw records.")
    
    # Map raw records by sequence_sha256
    raw_by_hash = df_raw.groupby("sequence_sha256")
    
    curated_records = []
    cross_file_overlaps = []
    conflicts = []
    exclusions = []
    review_queue = []
    
    for idx, u_row in df_unique.iterrows():
        seq_id = u_row["canonical_sequence_id"]
        seq_hash = u_row["sequence_sha256"]
        seq = u_row["sequence"]
        seq_len = u_row["sequence_length"]
        
        # Get all raw records for this sequence
        raw_group = raw_by_hash.get_group(seq_hash)
        raw_headers = raw_group["raw_header"].tolist()
        raw_record_ids = raw_group["raw_record_id"].tolist()
        source_files = sorted(raw_group["source_file"].unique().tolist())
        
        # Parse all headers
        headers_info = [extract_header_fields(h) for h in raw_headers]
        
        # 1. Check for cross-file overlap
        if len(source_files) > 1:
            # Multi-file sequence
            # Check consistency of accession and organism
            accs = list({h["accession"] for h in headers_info if h["accession"] != "UNKNOWN"})
            orgs = list({h["organism"] for h in headers_info if h["organism"] != "UNRESOLVED"})
            is_consistent = (len(orgs) <= 1)
            
            overlap_type = "SAME_ACCESSION_IDENTICAL_SEQ" if len(accs) == 1 else "DIFF_ACCESSION_IDENTICAL_SEQ"
            resolution = "UNRESOLVED_FLAGGED_FOR_REVIEW" if not is_consistent else "RETAINED_CANONICAL_REPRESENTATIVE"
            
            cross_file_overlaps.append({
                "sequence_id": seq_id,
                "source_file_1": source_files[0],
                "source_file_2": source_files[1],
                "overlap_type": overlap_type,
                "metadata_consistency": "CONSISTENT" if is_consistent else "CONFLICTING",
                "resolution": resolution,
                "notes": f"Raw headers: {' | '.join(raw_headers)}"
            })
            
        # 2. Taxonomic Curation
        organism, virus_species, tax_status, tax_evidence, has_tax_conflict = determine_taxonomy(headers_info)
        
        if has_tax_conflict:
            conflicts.append({
                "sequence_id": seq_id,
                "field": "virus_species / organism",
                "value_1": headers_info[0]["organism"],
                "value_2": headers_info[1]["organism"],
                "source_1": raw_record_ids[0],
                "source_2": raw_record_ids[1],
                "resolution_status": "REQUIRES_MANUAL_REVIEW"
            })
            
        # 3. Strain / Isolate Curation
        strains = list({h["strain"] for h in headers_info if h["strain"] != "UNKNOWN"})
        if len(strains) > 1:
            strain = "CONFLICTING"
            conflicts.append({
                "sequence_id": seq_id,
                "field": "strain",
                "value_1": strains[0],
                "value_2": strains[1],
                "source_1": raw_record_ids[0],
                "source_2": raw_record_ids[1],
                "resolution_status": "REQUIRES_MANUAL_REVIEW"
            })
        elif len(strains) == 1:
            strain = strains[0]
        else:
            strain = "UNKNOWN"
            
        isolates = list({h["isolate"] for h in headers_info if h["isolate"] != "UNKNOWN"})
        isolate = isolates[0] if isolates else "UNKNOWN"
        
        # 4. Protein / Gene Identity Curation
        prot_names = list({h["protein_name"] for h in headers_info if h["protein_name"] != "UNKNOWN"})
        gene_symbols = list({h["gene_symbol"] for h in headers_info if h["gene_symbol"] != "UNKNOWN"})
        accessions = list({h["accession"] for h in headers_info if h["accession"] != "UNKNOWN"})
        
        primary_prot = prot_names[0] if prot_names else "UNKNOWN"
        primary_gene = gene_symbols[0] if gene_symbols else "UNKNOWN"
        primary_acc = accessions[0] if accessions else "UNKNOWN"
        
        source_db = headers_info[0]["source_db"]
        source_record_type = "FASTA_HEADER_DERIVED"
        
        # Determine protein identity status (5 tiers)
        prot_lower = primary_prot.lower()
        if len(prot_names) > 1 and not (len(prot_names) == 2 and any(p in prot_names[1] for p in prot_names[0].split())):
            # Check for conflict in protein naming across multi-header duplicates
            if any(p.startswith("UL") or p.startswith("US") or p.startswith("glycoprotein") for p in prot_names):
                distinct_genes = {re.findall(r"\b(UL\d+|US\d+|RL\d+|RS\d+|ICP\d+|g[B-N])\b", p, re.I)[0].upper() 
                                  for p in prot_names if re.search(r"\b(UL\d+|US\d+|RL\d+|RS\d+|ICP\d+|g[B-N])\b", p, re.I)}
                if len(distinct_genes) > 1:
                    prot_id_status = "CONFLICTING"
                    prot_id_evidence = f"Conflicting protein gene identities across source headers: {'; '.join(prot_names)}"
                else:
                    prot_id_status = "VERIFIED"
                    prot_id_evidence = f"Explicit protein annotation: {primary_prot}"
            else:
                prot_id_status = "VERIFIED"
                prot_id_evidence = f"Explicit protein annotation: {primary_prot}"
        elif any(term in prot_lower for term in ["hypothetical", "uncharacterized", "predicted", "putative"]):
            prot_id_status = "PREDICTED_OR_HYPOTHETICAL"
            prot_id_evidence = f"Predicted/hypothetical annotation in source header: {primary_prot}"
        elif primary_prot.startswith("Sequence") and "patent" in prot_lower:
            prot_id_status = "UNRESOLVED"
            prot_id_evidence = "Patent sequence number placeholder"
        elif primary_prot == "UNKNOWN":
            prot_id_status = "UNRESOLVED"
            prot_id_evidence = "No protein descriptor available in source header"
        elif any(term in prot_lower for term in ["viral protein", "membrane protein", "structural protein", "tegument protein", "protein"]):
            prot_id_status = "CURATED_SUPPORTED"
            prot_id_evidence = f"Curated / general descriptor: {primary_prot}"
        else:
            prot_id_status = "VERIFIED"
            prot_id_evidence = f"Explicit characterized protein annotation: {primary_prot}"
            
        # 5. Source-Type & Completeness
        source_type, completeness, nat_or_eng, src_evidence, src_type_status = determine_source_type_and_completeness(
            headers_info, organism, virus_species
        )
        
        # 6. Sequence Composition Status
        comp_status, ambig_cnt, noncanon_cnt = determine_sequence_composition(seq)
        
        # 7. Biological Eligibility & Dataset Grouping
        has_any_conflict = (has_tax_conflict or strain == "CONFLICTING")
        eligibility, elig_reason, review_req, review_reason, group_cat = determine_biological_eligibility(
            virus_species, source_type, completeness, tax_status, prot_id_status, has_any_conflict, comp_status
        )
        
        # Record into manual review queue if required
        if review_req:
            review_queue.append({
                "sequence_id": seq_id,
                "virus_species": virus_species,
                "organism": organism,
                "source_type": source_type,
                "taxonomy_status": tax_status,
                "protein_identity_status": prot_id_status,
                "review_reason": review_reason,
                "source_files": ";".join(source_files),
                "raw_headers": " | ".join(raw_headers)
            })
            
        # Record into exclusion audit if excluded from primary natural proteome
        if eligibility.startswith("EXCLUDE_"):
            exclusions.append({
                "sequence_id": seq_id,
                "original_record_ids": ";".join(raw_record_ids),
                "exclusion_category": eligibility,
                "exclusion_reason": elig_reason,
                "evidence_source": tax_evidence if "tax" in elig_reason.lower() else src_evidence
            })
            
        curated_records.append({
            "sequence_id": seq_id,
            "source_file": ";".join(source_files),
            "original_record_ids": ";".join(raw_record_ids),
            "accession": primary_acc,
            "raw_header": " | ".join(raw_headers),
            "sequence": seq,
            "sequence_length": seq_len,
            "organism": organism,
            "virus_species": virus_species,
            "strain": strain,
            "isolate": isolate,
            "gene_symbol": primary_gene,
            "protein_name": primary_prot,
            "protein_accession": primary_acc,
            "source_database": source_db,
            "source_record_type": source_record_type,
            "source_type": source_type,
            "natural_or_engineered": nat_or_eng,
            "completeness_status": completeness,
            "sequence_composition_status": comp_status,
            "ambiguous_character_count": ambig_cnt,
            "noncanonical_character_count": noncanon_cnt,
            "taxonomy_status": tax_status,
            "taxonomy_evidence": tax_evidence,
            "protein_identity_status": prot_id_status,
            "protein_identity_evidence": prot_id_evidence,
            "source_type_status": src_type_status,
            "source_type_evidence": src_evidence,
            "biological_eligibility": eligibility,
            "eligibility_reason": elig_reason,
            "dataset_group": group_cat,
            "manual_review_required": review_req,
            "review_reason": review_reason,
            "provenance_status": "SURJECTIVELY_MAPPED_TO_RAW_FASTA"
        })
        
    df_curated = pd.DataFrame(curated_records)
    
    # Save Intermediate and Processed tables
    intermediate_dir = data_dir / "intermediate"
    processed_dir = data_dir / "processed"
    tables_dir = results_dir / "tables"
    logs_dir = results_dir / "logs"
    docs_dir = base_dir / "docs" / "phases"
    
    for d in [intermediate_dir, processed_dir, tables_dir, logs_dir, docs_dir]:
        d.mkdir(parents=True, exist_ok=True)
        
    df_curated.to_csv(intermediate_dir / "phase1b_biological_curation.csv", index=False)
    
    # Final processed biological dataset (without temporary internal intermediate columns)
    final_cols = [
        "sequence_id", "source_file", "original_record_ids", "accession",
        "raw_header", "sequence", "sequence_length", "organism",
        "virus_species", "strain", "isolate", "gene_symbol", "protein_name",
        "protein_accession", "source_database", "source_record_type",
        "source_type", "natural_or_engineered", "completeness_status",
        "sequence_composition_status", "ambiguous_character_count",
        "noncanonical_character_count", "taxonomy_status", "taxonomy_evidence",
        "protein_identity_status", "protein_identity_evidence",
        "source_type_status", "source_type_evidence", "biological_eligibility",
        "eligibility_reason", "dataset_group", "manual_review_required",
        "review_reason", "provenance_status"
    ]
    df_final = df_curated[final_cols].copy()
    df_final.to_csv(processed_dir / "phase1b_biological_dataset.csv", index=False)
    
    # Export specific audit and review tables
    df_overlaps = pd.DataFrame(cross_file_overlaps)
    df_overlaps.to_csv(tables_dir / "phase1b_cross_file_overlap.csv", index=False)
    
    df_conflicts = pd.DataFrame(conflicts)
    df_conflicts.to_csv(tables_dir / "phase1b_conflicts.csv", index=False)
    
    df_review = pd.DataFrame(review_queue)
    df_review.to_csv(tables_dir / "phase1b_manual_review_queue.csv", index=False)
    
    df_exclusions = pd.DataFrame(exclusions)
    df_exclusions.to_csv(tables_dir / "phase1b_exclusion_audit.csv", index=False)
    
    # Calculate Summary Statistics
    counts_species = df_final["virus_species"].value_counts().to_dict()
    counts_prot_status = df_final["protein_identity_status"].value_counts().to_dict()
    counts_source_type = df_final["source_type"].value_counts().to_dict()
    counts_completeness = df_final["completeness_status"].value_counts().to_dict()
    counts_eligibility = df_final["biological_eligibility"].value_counts().to_dict()
    counts_group = df_final["dataset_group"].value_counts().to_dict()
    
    summary_metrics = {
        "RAW_RECORDS": len(df_raw),
        "UNIQUE_SEQUENCES": len(df_final),
        "VIRUS_SPECIES": counts_species,
        "PROTEIN_IDENTITY_STATUS": counts_prot_status,
        "SOURCE_TYPE": counts_source_type,
        "COMPLETENESS_STATUS": counts_completeness,
        "BIOLOGICAL_ELIGIBILITY": counts_eligibility,
        "DATASET_GROUPS": counts_group,
        "CROSS_FILE_OVERLAPS": len(df_overlaps),
        "DETECTED_CONFLICTS": len(df_conflicts),
        "MANUAL_REVIEW_COUNT": len(df_review),
        "EXCLUSION_COUNT": len(df_exclusions)
    }
    
    # Write summary text
    with open(logs_dir / "phase1b_summary.txt", "w", encoding="utf-8") as f:
        f.write("=" * 80 + "\n")
        f.write("PHASE 1B: BIOLOGICAL ELIGIBILITY & SOURCE-TYPE CURATION SUMMARY\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"RAW RECORDS: {len(df_raw):,}\n")
        f.write(f"UNIQUE SEQUENCES: {len(df_final):,}\n\n")
        
        f.write("1. TAXONOMY DISTRIBUTION (virus_species):\n")
        for k, v in counts_species.items():
            f.write(f"  {k}: {v:,}\n")
            
        f.write("\n2. PROTEIN IDENTITY STATUS:\n")
        for k, v in counts_prot_status.items():
            f.write(f"  {k}: {v:,}\n")
            
        f.write("\n3. SOURCE TYPE:\n")
        for k, v in counts_source_type.items():
            f.write(f"  {k}: {v:,}\n")
            
        f.write("\n4. COMPLETENESS STATUS:\n")
        for k, v in counts_completeness.items():
            f.write(f"  {k}: {v:,}\n")
            
        f.write("\n5. BIOLOGICAL ELIGIBILITY:\n")
        for k, v in counts_eligibility.items():
            f.write(f"  {k}: {v:,}\n")
            
        f.write("\n6. DATASET GROUPS:\n")
        for k, v in counts_group.items():
            f.write(f"  {k}: {v:,}\n")
            
        f.write(f"\n7. CONFLICTS & AUDIT METRICS:\n")
        f.write(f"  Cross-File Overlaps: {len(df_overlaps)}\n")
        f.write(f"  Direct Metadata Conflicts: {len(df_conflicts)}\n")
        f.write(f"  Manual Review Queue Items: {len(df_review):,}\n")
        f.write(f"  Total Exclusions from Primary Proteome: {len(df_exclusions):,}\n\n")
        
        f.write("=" * 80 + "\n")
        f.write("PHASE 1B COMPLETE — BIOLOGICAL ELIGIBILITY CURATED; TEMPORAL ANNOTATION NOT YET PERFORMED.\n")
        f.write("=" * 80 + "\n")
        
    # Write Phase 1B Documentation Markdown
    doc_content = f"""# Phase 1B: Biological Eligibility & Source-Type Curation

## 1. Objective
Establish objective biological provenance, taxonomy, completeness, and source-type classifications for all 22,689 unique sequences derived from Phase 1A.

## 2. Evidence Hierarchy
1. Explicit FASTA header brackets (`[...]`) and database identifiers.
2. Verified accession prefixes (`sp|`, `tr|`, `prf|`, `pdb|`, `patent`).
3. UniProtKB mnemonic identifiers.
4. No sequence homology or machine learning inference was permitted during Phase 1B.

## 3. Curated Metrics Overview
- **Raw Ingested Records**: `{len(df_raw):,}`
- **Unique Canonical Sequences**: `{len(df_final):,}`

### Taxonomy (Virus Species)
- **HSV-1**: `{counts_species.get('HSV-1', 0):,}`
- **HSV-2**: `{counts_species.get('HSV-2', 0):,}`
- **OTHER_HERPESVIRUS**: `{counts_species.get('OTHER_HERPESVIRUS', 0):,}`
- **NON_HERPESVIRUS**: `{counts_species.get('NON_HERPESVIRUS', 0):,}`
- **UNKNOWN / CONFLICTING**: `{counts_species.get('UNKNOWN', 0):,}`

### Biological Eligibility Breakdown
- **ELIGIBLE_PRIMARY**: `{counts_eligibility.get('ELIGIBLE_PRIMARY', 0):,}`
- **ELIGIBLE_SEQUENCE_ANALYSIS_ONLY**: `{counts_eligibility.get('ELIGIBLE_SEQUENCE_ANALYSIS_ONLY', 0):,}`
- **EXCLUDE_NON_HSV**: `{counts_eligibility.get('EXCLUDE_NON_HSV', 0):,}`
- **EXCLUDE_SYNTHETIC**: `{counts_eligibility.get('EXCLUDE_SYNTHETIC', 0):,}`
- **EXCLUDE_TECHNICAL**: `{counts_eligibility.get('EXCLUDE_TECHNICAL', 0):,}`
- **EXCLUDE_ENGINEERED**: `{counts_eligibility.get('EXCLUDE_ENGINEERED', 0):,}`
- **EXCLUDE_UNRESOLVED_IDENTITY**: `{counts_eligibility.get('EXCLUDE_UNRESOLVED_IDENTITY', 0):,}`
- **REVIEW_REQUIRED**: `{counts_eligibility.get('REVIEW_REQUIRED', 0):,}`

## 4. Absolute Constraints & Governance
- **Temporal Annotation**: NOT PERFORMED (Strictly reserved for Phase 2).
- **Machine Learning / ProtBERT**: NONE.
- **Physical Features**: NONE.
- **Data Modification**: Raw `./raw/` files remained strictly read-only and unmutated.
"""
    with open(docs_dir / "PHASE_1B_BIOLOGICAL_ELIGIBILITY.md", "w", encoding="utf-8") as f:
        f.write(doc_content)
        
    print("Phase 1B Curation completed successfully!")
    return summary_metrics


if __name__ == "__main__":
    base_path = Path(__file__).resolve().parent.parent
    metrics = run_phase1b_curation(base_path)
    print(json.dumps(metrics, indent=2))
