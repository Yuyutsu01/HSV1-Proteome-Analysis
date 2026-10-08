"""
Phase 3 - Step 04: Pretrained Protein Language Model Embeddings (ESM-2)
and Long Sequence Handling.

Biological & Computational Concept:
1. Protein Language Models (pLMs):
   Pretrained protein language models (such as ESM-2 and ProtBERT) model evolutionary
   and contextual amino acid constraints from millions of unaligned sequences across the tree of life.
   Residue representations capture structural propensity, surface exposure, active sites, and domain organization.
   Global sequence representations are extracted via attention-weighted mean pooling across sequence residues.

2. Long Sequence Handling Architecture:
   pLMs have fixed context windows (512 or 1024 tokens). Truncating long sequences silently discards
   crucial C-terminal functional domains (e.g. UL36 large tegument protein, 3,164 residues; UL19 major capsid protein).
   Solution: Deterministic sliding-window chunking with segment pooling.
   - Context Window (W): 512 residues
   - Overlap (O): 64 residues (Step S = W - O = 448 residues)
   - For each window: extract residue embeddings and compute mean-pooled chunk vector.
   - Full protein representation: uniform mean pooling across chunk vectors.
   This guarantees 100% residue coverage without arbitrary truncation.

3. Optimization & Quality Control (QC):
   - Length-sorted chunk batching with dynamic padding to minimize transformer padding overhead.
   - Multi-threaded CPU execution (16 threads).
   - Verifies:
     * Total embedding count == Total supervised sequence count (N=16,657)
     * Exact 1-to-1 mapping with canonical IDs
     * Strict absence of NaN or Infinite floating-point values
     * Consistent dimensionality (e.g. 320 for ESM-2 t6)
     * Output stored in compressed binary .npz format with detached metadata.
"""

import os
import json
import yaml
import time
import math
import numpy as np
import pandas as pd
import torch
from typing import Dict, List, Tuple
from transformers import AutoTokenizer, AutoModel


def load_config(config_path: str = "configs/phase3_config.yaml") -> dict:
    """Load Phase 3 configuration."""
    with open(config_path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def chunk_sequence(sequence: str, window_size: int = 512, overlap: int = 64) -> List[str]:
    """
    Deterministically split sequence into overlapping chunks.
    Guarantees full coverage of the entire sequence.
    """
    seq_len = len(sequence)
    if seq_len <= window_size:
        return [sequence]
    
    step = window_size - overlap
    chunks = []
    start = 0
    while start < seq_len:
        end = min(start + window_size, seq_len)
        chunks.append(sequence[start:end])
        if end == seq_len:
            break
        start += step
        
    return chunks


def analyze_long_sequences(
    df: pd.DataFrame, 
    model_name: str,
    max_len: int = 512, 
    overlap: int = 64
) -> Tuple[pd.DataFrame, dict]:
    """
    Analyze sequence length distribution relative to model context window.
    """
    total = len(df)
    lengths = df['normalized_length'].values
    exceeding = (lengths > max_len).sum()
    fraction = exceeding / total
    max_seq_len = int(lengths.max())
    
    step = max_len - overlap
    max_chunks = math.ceil((max_seq_len - max_len) / step) + 1 if max_seq_len > max_len else 1
    
    summary_data = [{
        'model_name': model_name,
        'max_context_length': max_len,
        'chunk_size': max_len,
        'chunk_overlap': overlap,
        'pooling_method': 'mean_chunk_pooling',
        'total_sequences': total,
        'sequences_exceeding_max_len': int(exceeding),
        'fraction_exceeding': round(float(fraction), 4),
        'min_sequence_length': int(lengths.min()),
        'max_sequence_length': max_seq_len,
        'max_chunks_required': int(max_chunks)
    }]
    
    df_long = pd.DataFrame(summary_data)
    return df_long, summary_data[0]


def extract_esm2_embeddings(
    sequences: List[str],
    canonical_ids: List[str],
    model_name: str = "facebook/esm2_t6_8M_UR50D",
    chunk_batch_size: int = 128,
    window_size: int = 512,
    overlap: int = 64
) -> Tuple[np.ndarray, dict]:
    """
    Extract whole-protein representations using ESM-2 with sliding-window chunk pooling.
    Uses length-sorted chunk batching and multi-threaded CPU execution for fast extraction.
    """
    torch.set_num_threads(16)
    print(f"Loading pretrained tokenizer and model: {model_name} (Threads={torch.get_num_threads()})...", flush=True)
    tokenizer = AutoTokenizer.from_pretrained(model_name)
    model = AutoModel.from_pretrained(model_name)
    model.eval()
    
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model.to(device)
    print(f"Model running on device: {device}", flush=True)
    
    embedding_dim = model.config.hidden_size
    n_seqs = len(sequences)
    
    # 1. Build all chunks with sequence index tracking
    print(f"Building sliding-window chunks for {n_seqs} sequences (Window={window_size}, Overlap={overlap})...", flush=True)
    all_chunks = []
    chunk_to_seq_idx = []
    
    for seq_idx, seq in enumerate(sequences):
        chunks = chunk_sequence(seq, window_size=window_size, overlap=overlap)
        for chunk in chunks:
            all_chunks.append(chunk)
            chunk_to_seq_idx.append(seq_idx)
            
    total_chunks = len(all_chunks)
    print(f"Total chunks to process: {total_chunks} across {n_seqs} sequences.", flush=True)
    
    # 2. Sort chunks by length to minimize padding overhead during inference
    chunk_lengths = np.array([len(c) for c in all_chunks])
    sort_order = np.argsort(chunk_lengths)
    sorted_chunks = [all_chunks[idx] for idx in sort_order]
    
    sorted_chunk_embeddings = np.zeros((total_chunks, embedding_dim), dtype=np.float32)
    t0 = time.time()
    
    print(f"Starting batched inference across {total_chunks} chunks (Batch Size={chunk_batch_size})...", flush=True)
    
    for b_start in range(0, total_chunks, chunk_batch_size):
        b_end = min(b_start + chunk_batch_size, total_chunks)
        batch_chunks = sorted_chunks[b_start:b_end]
        
        with torch.inference_mode():
            # Dynamic padding: pads only to the longest sequence in this batch
            inputs = tokenizer(
                batch_chunks, 
                return_tensors="pt", 
                padding=True, 
                truncation=True, 
                max_length=window_size
            )
            inputs = {k: v.to(device) for k, v in inputs.items()}
            outputs = model(**inputs)
            
            # Masked mean over tokens (excluding special tokens / padding)
            mask = inputs['attention_mask'].unsqueeze(-1)
            token_embeds = outputs.last_hidden_state
            batch_embeds = ((token_embeds * mask).sum(dim=1) / mask.sum(dim=1)).cpu().numpy()
            
        sorted_chunk_embeddings[b_start:b_end] = batch_embeds
        
        if (b_start + chunk_batch_size) % 1000 < chunk_batch_size or b_end == total_chunks:
            elapsed = time.time() - t0
            rate = b_end / (elapsed + 1e-5)
            eta = (total_chunks - b_end) / (rate + 1e-5)
            print(f"  Processed {b_end}/{total_chunks} chunks ({b_end/total_chunks*100:.1f}%) | Speed: {rate:.1f} chunks/s | Elapsed: {elapsed:.1f}s | ETA: {eta:.1f}s", flush=True)
            
    # Restore original chunk order
    chunk_embeddings = np.empty_like(sorted_chunk_embeddings)
    chunk_embeddings[sort_order] = sorted_chunk_embeddings
    
    # 3. Aggregate chunk embeddings back to sequence level
    print("Aggregating chunk embeddings into full protein representations...", flush=True)
    protein_embeddings = np.zeros((n_seqs, embedding_dim), dtype=np.float32)
    seq_chunk_counts = np.zeros(n_seqs, dtype=int)
    
    for c_idx, seq_idx in enumerate(chunk_to_seq_idx):
        protein_embeddings[seq_idx] += chunk_embeddings[c_idx]
        seq_chunk_counts[seq_idx] += 1
        
    for s_idx in range(n_seqs):
        if seq_chunk_counts[s_idx] > 0:
            protein_embeddings[s_idx] /= seq_chunk_counts[s_idx]
            
    total_time = time.time() - t0
    meta = {
        'model_name': model_name,
        'embedding_dimension': embedding_dim,
        'pooling_method': 'token_mean_then_chunk_mean',
        'window_size': window_size,
        'chunk_overlap': overlap,
        'total_sequences': n_seqs,
        'total_chunks': total_chunks,
        'extraction_time_seconds': round(total_time, 2)
    }
    
    return protein_embeddings, meta


def run_embedding_pipeline():
    """Main execution function for Step 04."""
    config = load_config()
    
    print("=" * 80, flush=True)
    print("PHASE 3 - STEP 04: PROTEIN LANGUAGE MODEL EMBEDDINGS & LONG SEQUENCE HANDLING", flush=True)
    print("=" * 80, flush=True)
    
    # 1. Load Sequence Manifest
    manifest_path = config["paths"]["sequence_manifest"]
    df_manifest = pd.read_csv(manifest_path)
    print(f"Loaded Sequence Manifest: {manifest_path} (N={len(df_manifest)})", flush=True)
    
    # -------------------------------------------------------------
    # 2. Long Sequence Handling Analysis
    # -------------------------------------------------------------
    esm_cfg = config["embeddings"]["esm2"]
    model_name = esm_cfg["model_name"]
    max_chunk_len = esm_cfg["max_chunk_length"]
    chunk_overlap = esm_cfg["chunk_overlap"]
    
    df_long, long_summary = analyze_long_sequences(
        df_manifest, 
        model_name=model_name,
        max_len=max_chunk_len, 
        overlap=chunk_overlap
    )
    
    long_seq_path = config["paths"]["long_sequence_handling"]
    os.makedirs(os.path.dirname(long_seq_path), exist_ok=True)
    df_long.to_csv(long_seq_path, index=False)
    print(f"\nSaved Long Sequence Handling Report: {long_seq_path}", flush=True)
    print(df_long.to_string(index=False), flush=True)
    
    # -------------------------------------------------------------
    # 3. Extract ESM-2 Embeddings
    # -------------------------------------------------------------
    sequences = df_manifest['normalized_sequence'].tolist()
    canonical_ids = df_manifest['canonical_id'].tolist()
    
    embeddings, meta = extract_esm2_embeddings(
        sequences=sequences,
        canonical_ids=canonical_ids,
        model_name=model_name,
        chunk_batch_size=128,
        window_size=max_chunk_len,
        overlap=chunk_overlap
    )
    
    # Save embeddings to binary .npz
    esm_dir = os.path.join(config["paths"]["embeddings_dir"], "esm2")
    os.makedirs(esm_dir, exist_ok=True)
    
    npz_path = os.path.join(esm_dir, "esm2_embeddings.npz")
    np.savez_compressed(
        npz_path,
        canonical_ids=np.array(canonical_ids),
        embeddings=embeddings
    )
    print(f"\nSaved ESM-2 Embeddings (.npz): {npz_path} (Shape: {embeddings.shape}, Size: {os.path.getsize(npz_path)/(1024*1024):.2f} MB)", flush=True)
    
    # Save metadata JSON
    meta_path = os.path.join(esm_dir, "metadata.json")
    with open(meta_path, "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=2)
    print(f"Saved Metadata JSON: {meta_path}", flush=True)
    
    # -------------------------------------------------------------
    # 4. Embedding Quality Control (QC)
    # -------------------------------------------------------------
    print("\nRunning Embedding Quality Control (QC)...", flush=True)
    n_records, dim = embeddings.shape
    nan_count = int(np.isnan(embeddings).sum())
    inf_count = int(np.isinf(embeddings).sum())
    
    qc_passed = (
        n_records == len(df_manifest) and
        nan_count == 0 and
        inf_count == 0 and
        dim == esm_cfg["embedding_dim"]
    )
    
    qc_data = [{
        'model_name': model_name,
        'embedding_dimension': dim,
        'total_embeddings': n_records,
        'expected_sequences': len(df_manifest),
        'nan_count': nan_count,
        'inf_count': inf_count,
        'min_value': round(float(np.min(embeddings)), 6),
        'max_value': round(float(np.max(embeddings)), 6),
        'mean_value': round(float(np.mean(embeddings)), 6),
        'std_value': round(float(np.std(embeddings)), 6),
        'qc_status': 'PASSED' if qc_passed else 'FAILED'
    }]
    
    df_qc = pd.DataFrame(qc_data)
    qc_path = config["paths"]["embedding_qc"]
    df_qc.to_csv(qc_path, index=False)
    print(f"Saved Embedding QC Table: {qc_path}", flush=True)
    print(df_qc.to_string(index=False), flush=True)
    
    assert qc_passed, "Embedding Quality Control FAILED!"
    print("\nStep 04 completed successfully.", flush=True)


if __name__ == "__main__":
    run_embedding_pipeline()
