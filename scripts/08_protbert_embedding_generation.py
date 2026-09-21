"""
Phase 8: ProtBERT Deep Transformer Embeddings Generation

Model & Pipeline Specifications:
- Pretrained Model: Rostlab/prot_bert (HuggingFace Transformers, PyTorch)
- Architecture: 30-layer Transformer, 1024 hidden dimensions, 16 attention heads (420M parameters)
- Tokenization: Amino acid residues separated by single space (" ".join(list(sequence)))
- Max Sequence Handling: Standard BERT has a max length limit (typically 1024 tokens).
  For sequences exceeding the chunk limit (e.g., UL36 with 3139 residues), sliding window chunking
  (chunk_size=512, overlap=64) is applied and residue representations are averaged over the entire protein length.
- Pooling Strategy: Global Mean-Pooling over per-residue hidden states (excluding [CLS] and [SEP] special tokens).
- Output Representation: Fixed 1024-dimensional embedding vector per HSV-1 protein.

Author: Computational Biology Pipeline
"""

import os
import sys
import torch
import numpy as np
import pandas as pd
from Bio import SeqIO
from transformers import BertTokenizer, BertModel

def generate_protbert_embeddings(fasta_path="processed_data/non_redundant.fasta",
                                 output_npy="embeddings/protbert_embeddings.npy",
                                 output_meta="embeddings/protbert_metadata.csv",
                                 model_name="Rostlab/prot_bert",
                                 device=None):
    """
    Generate contextual 1024-dim ProtBERT embeddings for each HSV-1 protein.
    """
    if device is None:
        device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"Using compute device: {device}")
    
    print(f"Loading pretrained ProtBERT ({model_name})...")
    tokenizer = BertTokenizer.from_pretrained(model_name, do_lower_case=False)
    model = BertModel.from_pretrained(model_name)
    model = model.to(device)
    model.eval()
    
    records = list(SeqIO.parse(fasta_path, "fasta"))
    print(f"Processing {len(records)} protein sequences...")
    
    embeddings_list = []
    metadata_list = []
    
    with torch.no_grad():
        for i, rec in enumerate(records):
            prot_id = rec.id
            seq = str(rec.seq).upper()
            seq_len = len(seq)
            
            # ProtBERT expects space-separated amino acid tokens
            # Replace ambiguous residues or non-standard characters with <unk>
            clean_seq = " ".join(list(seq))
            
            # If sequence fits within 1022 tokens (+ CLS, SEP = 1024)
            if seq_len <= 1020:
                inputs = tokenizer(clean_seq, return_tensors="pt", add_special_tokens=True)
                inputs = {k: v.to(device) for k, v in inputs.items()}
                outputs = model(**inputs)
                # Last hidden state: shape (1, seq_len + 2, 1024)
                # Slice [1:-1] to exclude [CLS] and [SEP]
                token_embeddings = outputs.last_hidden_state[0, 1:-1, :]
                mean_embedding = token_embeddings.mean(dim=0).cpu().numpy()
            else:
                # Handle long sequences (e.g. UL36, UL30) via sliding window chunking
                chunk_size = 500
                step = 400
                chunk_embeds = []
                for start in range(0, seq_len, step):
                    end = min(start + chunk_size, seq_len)
                    sub_seq = seq[start:end]
                    sub_seq_str = " ".join(list(sub_seq))
                    inputs = tokenizer(sub_seq_str, return_tensors="pt", add_special_tokens=True)
                    inputs = {k: v.to(device) for k, v in inputs.items()}
                    outputs = model(**inputs)
                    token_embeddings = outputs.last_hidden_state[0, 1:-1, :]
                    chunk_mean = token_embeddings.mean(dim=0).cpu().numpy()
                    chunk_embeds.append(chunk_mean * len(sub_seq))
                    if end >= seq_len:
                        break
                # Weighted average over chunks
                mean_embedding = np.sum(chunk_embeds, axis=0) / seq_len
                
            embeddings_list.append(mean_embedding)
            metadata_list.append({
                'protein_id': prot_id,
                'gene': rec.description.split('|')[0] if '|' in rec.description else prot_id,
                'length': seq_len,
                'embedding_dim': len(mean_embedding)
            })
            
            if (i + 1) % 10 == 0 or (i + 1) == len(records):
                print(f"Processed [{i+1}/{len(records)}] proteins.")
                
    embeddings_arr = np.array(embeddings_list, dtype=np.float32)
    os.makedirs(os.path.dirname(output_npy), exist_ok=True)
    np.save(output_npy, embeddings_arr)
    
    meta_df = pd.DataFrame(metadata_list)
    meta_df.to_csv(output_meta, index=False)
    
    print(f"\nSuccessfully generated ProtBERT embeddings matrix: {embeddings_arr.shape}")
    print(f"Saved to {output_npy} and metadata to {output_meta}")
    return embeddings_arr, meta_df

if __name__ == "__main__":
    generate_protbert_embeddings()
