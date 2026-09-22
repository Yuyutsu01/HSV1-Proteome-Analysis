# Data Processing Decision Record

## Data Preprocessing, Alignment, and Feature Engineering Rules

This document details the specific algorithmic rules governing data processing, feature computation, and matrix construction.

---

### 1. Exact Deduplication
- **Method:** SHA-256 cryptographic hashing on translated amino acid sequences.
- **Rules:**
  - Sequences with identical hashes are grouped.
  - The canonical locus (standard naming / lower coordinate) is retained.
  - Exactly 3 duplicates identified and removed: `YP_009137133.1` (RL2), `YP_009137134.1` (RL1), `YP_009137149.1` (RS1).
  - All removals are logged in `data/intermediate/duplicates.csv`.

---

### 2. Protein ID Alignment and Master Index
- **Rule:** Every matrix row must correspond to the exact same protein ID in the same canonical order.
- **Verification:**
  - `representation_master_index.csv` defines the ground-truth ordering ($N=74$).
  - Bijective 1:1 row alignment verified across FASTA, CSVs, NPY matrices, and metadata files.

---

### 3. Amino Acid Composition Representation
- **Rule:** 20 canonical amino acid fractions are computed as:
  $$\text{fraction}_{AA} = \frac{\text{count}(AA)}{\text{sequence length}}$$
- **Integrity:** $\sum_{i=1}^{20} \text{fraction}_i = 1.000000 \pm 10^{-6}$.
- Ambiguous amino acids (`B`, `Z`, `X`, `J`) are strictly disallowed and confirmed absent during QC.

---

### 4. ProtBERT Tokenization, Chunking, and Pooling
- **Model:** `Rostlab/prot_bert`
- **Tokenization:** Space-separated residues (`"M A T E..."`) with `[CLS]` at the start and `[SEP]` at the end.
- **Chunking Parameters:**
  - Window capacity: $C = 500$ amino acids.
  - Overlap: $O = 100$ amino acids.
  - Stride/Step: $S = 400$ amino acids.
  - Token sequence length per chunk: $\le 502$ tokens.
- **Pooling Logic:**
  - Intra-chunk: Mean pooling over hidden states corresponding to amino acid residues (excluding special `[CLS]` and `[SEP]` tokens and attention-padded positions).
  - Inter-chunk: Length-weighted mean pooling across all chunk vectors for multi-chunk proteins:
    $$\mathbf{e}_{\text{protein}} = \sum_{k=1}^K \frac{L_k}{\sum_{j=1}^K L_j} \mathbf{e}_k$$
    where $L_k$ is the number of valid residues in chunk $k$.

---

### 5. Multi-Modal Representation Matrix Construction
- **Matrix Specifications:**
  - $X_{\text{physicochemical}} \in \mathbb{R}^{74 \times 25}$ (float64)
  - $X_{\text{protbert}} \in \mathbb{R}^{74 \times 1024}$ (float32)
  - $X_{\text{combined\_raw}} = [X_{\text{physicochemical}} \,\|\, X_{\text{protbert}}] \in \mathbb{R}^{74 \times 1049}$ (float64)
- **Zero Pre-Scaling Rule:** $X_{\text{combined\_raw}}$ is preserved strictly in unscaled form. No global standardization is performed prior to cross-validation or clustering experiments.
