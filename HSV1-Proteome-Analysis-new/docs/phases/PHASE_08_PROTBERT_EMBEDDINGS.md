# Phase 08: ProtBERT Embedding Generation

## 1. Objective
Generate protein language model representations ($74 \times 1024$) for all 74 unique HSV-1 proteins using `Rostlab/prot_bert`, implementing a sliding-window chunking strategy with length-weighted pooling for long polypeptides, followed by a methodological chunking audit.

## 2. Inputs
- `data/processed/unique_proteins.csv` (74 unique proteins)
- Pretrained Transformer Model: `Rostlab/prot_bert` (Hugging Face)

## 3. Processing
- Executed `scripts/08_protbert_embeddings.py`.
- **Tokenization:** Space-separated amino acid tokens with special boundary tokens `[CLS]` and `[SEP]`.
- **Long Sequence Handling:**
  - Chunk size: $C = 500$ amino acids.
  - Overlap: $O = 100$ amino acids.
  - Step size: $S = 400$ amino acids.
  - Chunk token length: $500 + 2 = 502$ tokens (including `[CLS]` and `[SEP]`).
  - Pooling: Mean-pooling over residue positions per chunk (excluding special tokens), followed by sequence-length-weighted mean pooling across chunks for multi-chunk proteins.
- Executed `scripts/08_audit_chunking.py` to evaluate architectural context limits and aggregation sensitivity.
- Executed `tests/test_phase8_protbert.py` (6 automated test cases).

## 4. Outputs
- `data/processed/protbert_embeddings.npy` ($74 \times 1024$ NumPy array, float32)
- `data/processed/protbert_metadata.csv` (74 rows, aligned metadata)
- `data/processed/protbert_config.json` (Configuration parameters)
- `results/tables/protbert_sequence_handling.csv` (Chunk counts per protein)
- `results/tables/protbert_chunking_audit.csv`
- `results/tables/protbert_aggregation_sensitivity.csv` (Sensitivity metrics)
- `results/tables/protbert_embedding_statistics.csv`
- `results/logs/phase8_chunking_audit.txt`
- `results/logs/phase8_validation_report.txt`

## 5. Validation and Methodological Audit

### Sequence Chunking Breakdown ($N=74$)
- **Single-chunk proteins ($\le 500$ aa):** 45 proteins (60.81%)
- **Multi-chunk proteins ($> 500$ aa):** 29 proteins (39.19%)
- **Longest protein:** `UL36` (3,139 aa) $\rightarrow$ 8 chunks.
- **Coverage:** **100.0% sequence coverage across all 74 proteins (zero truncation)**.

### Model Context Verification
- `model.config.max_position_embeddings = 40000`
- `tokenizer.model_max_length = 1000000000000000019884624838656` (effectively unbounded).
- *Methodological note:* The chosen 500-aa sequence chunks (502 tokens with special tokens) remain strictly within the effective 512-token context window utilized during ProtBERT pretraining, preventing out-of-distribution positional extrapolation.

### Aggregation Sensitivity Benchmark (Method A vs Method B)
Four representative benchmark proteins spanning sequence length classes were audited comparing length-weighted pooling (Method A) vs uniform chunk averaging (Method B):

| Benchmark Protein | Gene | Sequence Length (aa) | Number of Chunks | Cosine Similarity (A vs B) | L2 Norm Delta |
| :--- | :--- | :---: | :---: | :---: | :---: |
| `YP_009137111.1` | `UL36` | 3,139 | 8 | **0.999928** | 0.009384 |
| `YP_009137094.1` | `UL19` | 1,374 | 4 | **0.999870** | 0.013589 |
| `YP_009137105.1` | `UL30` | 1,235 | 3 | **0.999347** | 0.030588 |
| `YP_009137129.1` | `UL54` | 512 | 2 | **0.999991** | 0.003507 |

*Note: This sensitivity analysis evaluates a four-protein benchmark audit and is not generalized to an equivalence claim across all possible pooling variants.*

- **Validation Suite:** `tests/test_phase8_protbert.py` (**6/6 tests passed**). Zero NaNs/Infs, deterministic float32 reproduction ($\text{max } |\Delta| = 0.0$).

## 6. Scientific Decisions
- **Zero Label Leakage:** ProtBERT embeddings were generated in an entirely unsupervised manner without class labels or temporal supervision.
- **Pooling Logic:** Length-weighted averaging was selected as primary to proportionally weight the final partial chunk.

## 7. Important Limitations
- Protein language models capture statistical sequence semantics and evolutionary co-occurrence; viral-host interface dynamics and physical complex assemblies are not explicitly represented.

## 8. Files Generated
- [`data/processed/protbert_embeddings.npy`](file:///c:/Users/shiva/OneDrive/Desktop/HSV1-Proteome-Analysis/HSV1-Proteome-Analysis-new/data/processed/protbert_embeddings.npy)
- [`results/tables/protbert_chunking_audit.csv`](file:///c:/Users/shiva/OneDrive/Desktop/HSV1-Proteome-Analysis/HSV1-Proteome-Analysis-new/results/tables/protbert_chunking_audit.csv)
- [`tests/test_phase8_protbert.py`](file:///c:/Users/shiva/OneDrive/Desktop/HSV1-Proteome-Analysis/HSV1-Proteome-Analysis-new/tests/test_phase8_protbert.py)

## 9. Status
**COMPLETE & VALIDATED**
