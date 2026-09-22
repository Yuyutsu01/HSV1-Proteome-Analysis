"""
Script 08: ProtBERT Protein-Language-Model Embedding Generation

Scientific Objective & Theoretical Concept:
-------------------------------------------
Protein Language Models (pLMs) such as ProtBERT (Elnaggar et al., 2021) leverage transformer
encoder architectures trained on massive evolutionary sequence databases (UniRef100, 216M proteins).
By modeling context-dependent amino acid probabilities through masked language modeling (MLM),
ProtBERT captures implicit biophysical constraints, secondary structure motifs, and evolutionary
signatures beyond static scalar physicochemical features.

Embedding Architecture & Strategy:
- Base Model: Rostlab/prot_bert (1024-dimensional final hidden state representation).
- Tokenization: Amino acid residues separated by single whitespace; enclosed with [CLS] and [SEP].
- Sequence Chunking: ProtBERT has a maximum sequence context of 512 tokens.
  For long viral proteins (e.g., UL36 giant tegument = 3,139 aa), we apply a deterministic
  sliding-window chunking strategy:
  * Chunk Size: 500 amino acids
  * Overlap: 100 amino acids (Step = 400 amino acids)
  * Effective token count per chunk: 500 + 2 special tokens ([CLS], [SEP]) = 502 tokens (<= 512 max limit).
  * Final chunk of any protein spans [max(0, L - 500), L] to ensure 100% complete residue coverage.
- Residue-Level Pooling: Within each chunk, special tokens ([CLS], [SEP], [PAD]) are masked out.
  The final hidden states of the actual amino acid residues are averaged to produce a chunk embedding.
- Protein-Level Aggregation: Multi-chunk embeddings are aggregated via length-weighted average
  based on the number of amino acids in each chunk:
  e_protein = sum(w_k * e_k) / sum(w_k), where w_k = len(chunk_k).
- Label Independence: No biological or temporal labels are exposed during embedding generation.

Author: Computational Biology Pipeline
"""

import os
import sys
import json
import yaml
import hashlib
import platform
import numpy as np
import pandas as pd
import torch
from transformers import BertTokenizer, BertModel

# Set deterministic execution seeds
RANDOM_SEED = 42
torch.manual_seed(RANDOM_SEED)
np.random.seed(RANDOM_SEED)
if torch.cuda.is_available():
    torch.cuda.manual_seed_all(RANDOM_SEED)

def load_config(config_path="config.yaml"):
    """Load pipeline YAML configuration."""
    if not os.path.exists(config_path):
        if os.path.exists(os.path.join("..", config_path)):
            config_path = os.path.join("..", config_path)
    with open(config_path, "r") as f:
        return yaml.safe_load(f)

def compute_sha256(filepath):
    """Compute SHA-256 hash for provenance."""
    hasher = hashlib.sha256()
    with open(filepath, "rb") as f:
        while chunk := f.read(65536):
            hasher.update(chunk)
    return hasher.hexdigest()

def chunk_sequence(sequence, chunk_size=500, overlap=100):
    """
    Deterministically partitions a protein sequence into overlapping chunks.
    
    Parameters:
        sequence (str): Full amino acid sequence.
        chunk_size (int): Maximum residue length per chunk (default 500).
        overlap (int): Overlap between successive chunks (default 100).
        
    Returns:
        list of tuples: (chunk_str, start_idx, end_idx)
    """
    seq_len = len(sequence)
    if seq_len <= chunk_size:
        return [(sequence, 0, seq_len)]
    
    step = chunk_size - overlap
    chunks = []
    start = 0
    while start < seq_len:
        end = min(start + chunk_size, seq_len)
        chunks.append((sequence[start:end], start, end))
        if end == seq_len:
            break
        start += step
        # If remaining residues would form a tiny tail, ensure the final window captures the terminal segment
        if start + chunk_size > seq_len and start < seq_len:
            final_start = max(0, seq_len - chunk_size)
            if final_start > chunks[-1][1]:
                chunks.append((sequence[final_start:seq_len], final_start, seq_len))
            break
            
    return chunks

def extract_protein_embedding(model, tokenizer, sequence, device, chunk_size=500, overlap=100):
    """
    Computes a 1024-dimensional ProtBERT embedding for a complete protein sequence
    using sliding-window chunking and special-token-masked mean pooling.
    
    Returns:
        tuple: (embedding_1024_np_array, num_chunks, covered_residues, coverage_fraction)
    """
    chunks = chunk_sequence(sequence, chunk_size=chunk_size, overlap=overlap)
    seq_len = len(sequence)
    
    # Track covered residue indices for coverage validation
    covered_set = set()
    chunk_embeddings = []
    chunk_weights = []
    
    for chunk_seq, start_idx, end_idx in chunks:
        # Mark covered indices
        for i in range(start_idx, end_idx):
            covered_set.add(i)
            
        # Format as space-separated amino acid tokens
        spaced_seq = " ".join(list(chunk_seq))
        
        # Tokenize sequence
        inputs = tokenizer(
            spaced_seq,
            return_tensors="pt",
            padding=False,
            truncation=True,
            max_length=512
        )
        
        input_ids = inputs["input_ids"].to(device)
        attention_mask = inputs["attention_mask"].to(device)
        
        with torch.no_grad():
            outputs = model(input_ids=input_ids, attention_mask=attention_mask)
            last_hidden_state = outputs.last_hidden_state  # Shape: (1, seq_len_tokens, 1024)
            
        # Mask out special tokens: [CLS] is token 0, [SEP] is token -1 (or last active token)
        # Residue tokens are at indices 1 to (num_residues)
        num_residues = len(chunk_seq)
        residue_tokens_hidden = last_hidden_state[0, 1:1 + num_residues, :]  # Shape: (num_residues, 1024)
        
        # Mean pooling across amino acid residues within this chunk
        chunk_mean = residue_tokens_hidden.mean(dim=0).cpu().numpy()  # Shape: (1024,)
        chunk_embeddings.append(chunk_mean)
        chunk_weights.append(num_residues)
        
    # Aggregate chunks via length-weighted average
    total_weight = sum(chunk_weights)
    aggregated_embedding = np.zeros(1024, dtype=np.float32)
    for c_emb, c_w in zip(chunk_embeddings, chunk_weights):
        aggregated_embedding += (c_w / total_weight) * c_emb
        
    covered_residues = len(covered_set)
    coverage_fraction = float(covered_residues / seq_len)
    
    return aggregated_embedding, len(chunks), covered_residues, coverage_fraction

def generate_protbert_embeddings(config_path="config.yaml"):
    """
    Executes Phase 8: Generates ProtBERT embeddings for all 74 proteins,
    validates dimensions, coverage, missing values, and checks reproducibility.
    """
    cfg = load_config(config_path)
    proc_dir = cfg['paths']['data_processed']
    tables_dir = cfg['paths']['results_tables']
    logs_dir = cfg['paths']['results_logs']
    
    os.makedirs(proc_dir, exist_ok=True)
    os.makedirs(tables_dir, exist_ok=True)
    os.makedirs(logs_dir, exist_ok=True)
    
    unique_csv = os.path.join(proc_dir, "unique_proteins.csv")
    output_emb_path = os.path.join(proc_dir, "protbert_embeddings.npy")
    output_meta_path = os.path.join(proc_dir, "protbert_metadata.csv")
    output_cfg_path = os.path.join(proc_dir, "protbert_config.json")
    
    assert os.path.exists(unique_csv), f"Missing {unique_csv}"
    df_unique = pd.read_csv(unique_csv)
    
    print("=" * 60)
    print("PHASE 8: PROTBERT EMBEDDING GENERATION (Rostlab/prot_bert)")
    print("=" * 60)
    print(f"Total Unique Proteins to Embed: {len(df_unique)}")
    
    # 1. Device selection
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Execution Device: {device} (PyTorch {torch.__version__})")
    
    # 2. Model & Tokenizer Initialization
    model_name = "Rostlab/prot_bert"
    print(f"Loading Tokenizer & Model: {model_name}...")
    tokenizer = BertTokenizer.from_pretrained(model_name, do_lower_case=False)
    model = BertModel.from_pretrained(model_name)
    model.to(device)
    model.eval()
    
    chunk_size = 500
    overlap = 100
    
    # Save ProtBERT configuration metadata
    protbert_config = {
        "model_name": model_name,
        "model_architecture": "BertModel",
        "hidden_dimension": 1024,
        "max_model_tokens": 512,
        "special_tokens": ["[CLS]", "[SEP]", "[PAD]"],
        "effective_amino_acid_capacity": 510,
        "chunk_size": chunk_size,
        "overlap": overlap,
        "pooling_method": "residue_mean_excluding_special_tokens",
        "aggregation_method": "length_weighted_chunk_average",
        "random_seed": RANDOM_SEED,
        "device": str(device),
        "pytorch_version": torch.__version__,
        "transformers_version": getattr(sys.modules.get('transformers'), '__version__', 'unknown'),
        "python_version": platform.python_version()
    }
    with open(output_cfg_path, "w") as f:
        json.dump(protbert_config, f, indent=2)
    print(f"[Saved] ProtBERT Configuration: {output_cfg_path}")
    
    # 3. Embedding Generation Loop
    embeddings_list = []
    metadata_rows = []
    
    for idx, row in df_unique.iterrows():
        prot_id = str(row['protein_id'])
        gene = str(row['gene'])
        seq = str(row['sequence'])
        expected_len = int(row['sequence_length'])
        
        emb, n_chunks, cov_res, cov_frac = extract_protein_embedding(
            model=model,
            tokenizer=tokenizer,
            sequence=seq,
            device=device,
            chunk_size=chunk_size,
            overlap=overlap
        )
        
        if cov_frac < 0.9999:
            raise ValueError(f"CRITICAL ERROR: Protein {prot_id} ({gene}) has incomplete coverage: {cov_frac:.4f}")
            
        embeddings_list.append(emb)
        metadata_rows.append({
            "embedding_index": idx,
            "protein_id": prot_id,
            "gene": gene,
            "sequence_length": expected_len,
            "number_of_chunks": n_chunks,
            "covered_residues": cov_res,
            "coverage_fraction": round(cov_frac, 6)
        })
        
        if (idx + 1) % 15 == 0 or (idx + 1) == len(df_unique):
            print(f"  Processed {idx + 1:2d} / {len(df_unique):2d} proteins (Latest: {gene:<6} Length: {expected_len:<5} Chunks: {n_chunks})")
            
    embedding_matrix = np.array(embeddings_list, dtype=np.float32)  # Shape: (74, 1024)
    df_meta = pd.DataFrame(metadata_rows)
    
    # Save main outputs
    np.save(output_emb_path, embedding_matrix)
    df_meta.to_csv(output_meta_path, index=False)
    print(f"[Saved] ProtBERT Embedding Matrix (.npy): {output_emb_path} (Shape: {embedding_matrix.shape})")
    print(f"[Saved] ProtBERT Metadata: {output_meta_path}")
    
    # 4. Long Sequence Handling Audit Table
    df_seq_handling = df_meta.sort_values(by="sequence_length", ascending=False).copy()
    seq_handling_csv = os.path.join(tables_dir, "protbert_sequence_handling.csv")
    df_seq_handling.to_csv(seq_handling_csv, index=False)
    print(f"[Saved] Long Sequence Handling Audit Table: {seq_handling_csv}")
    
    # 5. Embedding Matrix Statistics Table
    matrix_mean = float(np.mean(embedding_matrix))
    matrix_std = float(np.std(embedding_matrix))
    matrix_min = float(np.min(embedding_matrix))
    matrix_max = float(np.max(embedding_matrix))
    l2_norms = np.linalg.norm(embedding_matrix, axis=1)
    mean_l2 = float(np.mean(l2_norms))
    std_l2 = float(np.std(l2_norms))
    min_l2 = float(np.min(l2_norms))
    max_l2 = float(np.max(l2_norms))
    
    stats_records = [
        {"metric": "Matrix Row Count (Proteins)", "value": str(embedding_matrix.shape[0])},
        {"metric": "Embedding Dimensionality", "value": str(embedding_matrix.shape[1])},
        {"metric": "Global Mean", "value": f"{matrix_mean:.6f}"},
        {"metric": "Global Standard Deviation", "value": f"{matrix_std:.6f}"},
        {"metric": "Global Minimum Value", "value": f"{matrix_min:.6f}"},
        {"metric": "Global Maximum Value", "value": f"{matrix_max:.6f}"},
        {"metric": "Mean L2 Norm per Protein", "value": f"{mean_l2:.4f}"},
        {"metric": "Std L2 Norm per Protein", "value": f"{std_l2:.4f}"},
        {"metric": "Min L2 Norm", "value": f"{min_l2:.4f}"},
        {"metric": "Max L2 Norm", "value": f"{max_l2:.4f}"},
    ]
    df_stats = pd.DataFrame(stats_records)
    stats_csv = os.path.join(tables_dir, "protbert_embedding_statistics.csv")
    df_stats.to_csv(stats_csv, index=False)
    print(f"[Saved] Embedding Matrix Statistics: {stats_csv}")
    
    # 6. Reproducibility Test: Re-embed sample proteins and evaluate numerical delta
    print("\n--- Running Reproducibility Verification ---")
    repro_check_genes = ['UL36', 'UL30', 'UL19', 'US12']
    max_abs_diff = 0.0
    mean_abs_diff_list = []
    
    for g in repro_check_genes:
        row = df_unique[df_unique['gene'] == g].iloc[0]
        seq = str(row['sequence'])
        orig_idx = df_unique[df_unique['gene'] == g].index[0]
        orig_emb = embedding_matrix[orig_idx]
        
        re_emb, _, _, _ = extract_protein_embedding(
            model=model,
            tokenizer=tokenizer,
            sequence=seq,
            device=device,
            chunk_size=chunk_size,
            overlap=overlap
        )
        diff = np.abs(orig_emb - re_emb)
        max_d = float(np.max(diff))
        mean_d = float(np.mean(diff))
        max_abs_diff = max(max_abs_diff, max_d)
        mean_abs_diff_list.append(mean_d)
        
    avg_mean_diff = float(np.mean(mean_abs_diff_list))
    repro_pass = (max_abs_diff < 1e-5)
    print(f"  Reproducibility Check on {repro_check_genes}: Max Abs Diff = {max_abs_diff:.2e}, Mean Abs Diff = {avg_mean_diff:.2e} -> {'PASS' if repro_pass else 'FAIL'}")
    
    # 7. Automated Integrity Checks
    v1_rows = (embedding_matrix.shape[0] == 74)
    v2_dim = (embedding_matrix.shape[1] == 1024)
    v3_shape = (embedding_matrix.shape == (74, 1024))
    v4_nans = (np.isnan(embedding_matrix).sum() == 0)
    v5_infs = (np.isinf(embedding_matrix).sum() == 0)
    v6_unique = (df_meta['protein_id'].nunique() == 74)
    v7_order = all(df_meta['protein_id'].values == df_unique['protein_id'].values)
    v8_lens = all(df_meta['sequence_length'].values == df_unique['sequence_length'].values)
    v9_cov = all(df_meta['coverage_fraction'].values == 1.0)
    
    # Check that long proteins (>500 aa) have > 1 chunks and short have 1 chunk
    long_mask = df_meta['sequence_length'] > 500
    v10_long_chunks = all(df_meta.loc[long_mask, 'number_of_chunks'] > 1)
    short_mask = df_meta['sequence_length'] <= 500
    v11_short_chunks = all(df_meta.loc[short_mask, 'number_of_chunks'] == 1)
    
    all_checks = [
        v1_rows, v2_dim, v3_shape, v4_nans, v5_infs,
        v6_unique, v7_order, v8_lens, v9_cov,
        v10_long_chunks, v11_short_chunks, repro_pass
    ]
    all_passed = all(all_checks)
    
    # 8. Write Provenance & Validation Logs
    raw_hash = compute_sha256(unique_csv)
    emb_hash = compute_sha256(output_emb_path)
    meta_hash = compute_sha256(output_meta_path)
    
    prov_log_lines = [
        "==================================================",
        "PROTBERT GENERATION PROVENANCE LOG",
        "==================================================",
        f"Execution Timestamp:     {pd.Timestamp.now().isoformat()}",
        f"Python Version:          {platform.python_version()}",
        f"PyTorch Version:         {torch.__version__}",
        f"Transformers Version:    {getattr(sys.modules.get('transformers'), '__version__', 'unknown')}",
        f"Model Name:              {model_name}",
        f"Device:                  {device}",
        f"Random Seed:             {RANDOM_SEED}",
        f"Chunk Size:              {chunk_size} aa",
        f"Overlap:                 {overlap} aa",
        f"Pooling Method:          Residue-Level Mean (Excluding [CLS], [SEP], [PAD])",
        f"Aggregation Method:      Length-Weighted Chunk Mean",
        f"Input unique_proteins:   {raw_hash}",
        f"Output protbert.npy:     {emb_hash}",
        f"Output metadata:         {meta_hash}",
        "=================================================="
    ]
    prov_log_path = os.path.join(logs_dir, "protbert_generation.log")
    with open(prov_log_path, "w") as f:
        f.write("\n".join(prov_log_lines))
    print(f"[Saved] Provenance Log: {prov_log_path}")
    
    # Validation Report
    longest_row = df_meta.loc[df_meta['sequence_length'].idxmax()]
    
    report_lines = [
        "==================================================",
        "PHASE 8 VALIDATION REPORT: PROTBERT EMBEDDINGS",
        "==================================================",
        f"Dataset:                  74 proteins",
        f"Embedding matrix shape:   {embedding_matrix.shape[0]} x {embedding_matrix.shape[1]}",
        f"Model:                    {model_name}",
        f"Pooling strategy:         Residue-Level Mean (Excluding Special Tokens)",
        f"Chunking strategy:        Sliding Window (Size={chunk_size} aa, Overlap={overlap} aa)",
        f"Longest protein:          {longest_row['gene']} ({longest_row['protein_id']}, {longest_row['sequence_length']} aa)",
        f"Maximum number of chunks: {df_meta['number_of_chunks'].max()} chunks ({longest_row['gene']})",
        f"Coverage validation:      {'PASS' if v9_cov else 'FAIL'} (100% residues covered across all 74 proteins)",
        f"NaN / Inf validation:     {'PASS' if (v4_nans and v5_infs) else 'FAIL'} (0 NaN, 0 Inf)",
        f"Reproducibility result:   {'PASS' if repro_pass else 'FAIL'} (Max Abs Diff = {max_abs_diff:.2e})",
        f"Overall Phase 8 status:   {'PASS' if all_passed else 'FAIL'}",
        "==================================================",
        "",
        "Detailed Check List:",
        f"  - Number of proteins = 74:            {'PASS' if v1_rows else 'FAIL'}",
        f"  - Embedding dimension = 1024:         {'PASS' if v2_dim else 'FAIL'}",
        f"  - Matrix shape = 74 x 1024:           {'PASS' if v3_shape else 'FAIL'}",
        f"  - No NaN values:                      {'PASS' if v4_nans else 'FAIL'}",
        f"  - No infinite values:                 {'PASS' if v5_infs else 'FAIL'}",
        f"  - Unique protein IDs = 74:            {'PASS' if v6_unique else 'FAIL'}",
        f"  - Metadata order matches embeddings:  {'PASS' if v7_order else 'FAIL'}",
        f"  - Sequence lengths match input:       {'PASS' if v8_lens else 'FAIL'}",
        f"  - 100% Coverage fraction:             {'PASS' if v9_cov else 'FAIL'}",
        f"  - Long proteins chunked (>500 aa):    {'PASS' if v10_long_chunks else 'FAIL'}",
        f"  - Short proteins single chunk:        {'PASS' if v11_short_chunks else 'FAIL'}",
        f"  - Deterministic reproducibility:      {'PASS' if repro_pass else 'FAIL'}",
        "\n[STOP CONDITION REACHED] Phase 8 completed. Standby for Phase 9."
    ]
    report_text = "\n".join(report_lines)
    report_path = os.path.join(logs_dir, "phase8_validation_report.txt")
    with open(report_path, "w") as f:
        f.write(report_text)
    print(f"[Saved] Validation Report: {report_path}")
    print("\n" + report_text)
    
    if not all_passed:
        print("CRITICAL ERROR: Phase 8 validation failed!")
        sys.exit(1)

if __name__ == "__main__":
    generate_protbert_embeddings()
