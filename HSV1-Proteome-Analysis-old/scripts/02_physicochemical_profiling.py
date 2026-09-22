"""
Phase 2: Physicochemical Feature Extraction for HSV-1 Proteome

Biological & Computational Context:
Calculates 26 physicochemical features for each unique protein sequence:
- Sequence Length (L): number of residues
- Molecular Weight (MW): in Daltons / converted to kDa
- Isoelectric Point (pI): theoretical pH at which net charge is zero
- Instability Index (II): Guruprasad index (>40 indicates in vitro instability)
- Aromaticity (Ar): relative frequency of Phe + Trp + Tyr
- 20 Canonical Amino Acid Compositions (f_A, f_C, ..., f_Y): molar percentages / fractions

Dimensionality Note:
The full feature matrix is 74 x 26. In some earlier literature, when 20 amino acid fractions
are used, their sum is 1.0 (a linear dependency / simplex constraint). Dropping 1 amino acid
fraction yields 25 linearly independent variables for unregularized covariance PCA.
We document both 26-feature full representation and standardized PCA analysis.

Author: Computational Biology Pipeline
"""

import os
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

def compute_physicochemical_features(fasta_path="processed_data/non_redundant.fasta",
                                     curated_path="processed_data/curated_annotations.csv",
                                     output_path="processed_data/physicochemical_features.csv"):
    """
    Extract 26 physicochemical properties for each deduplicated protein sequence.
    """
    annotations_df = pd.read_csv(curated_path)
    annot_map = {row['protein_id']: row for _, row in annotations_df.iterrows()}
    
    records = list(SeqIO.parse(fasta_path, "fasta"))
    print(f"Calculating physicochemical properties for {len(records)} proteins...")
    
    amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 
                   'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    
    feature_rows = []
    
    for rec in records:
        prot_id = rec.id
        seq_str = str(rec.seq).upper()
        # Clean sequence of any non-standard residues for ProtParam
        clean_seq = "".join([aa for aa in seq_str if aa in amino_acids])
        
        analyser = ProteinAnalysis(clean_seq)
        
        length = len(seq_str)
        mw_da = analyser.molecular_weight()
        mw_kda = mw_da / 1000.0
        aromaticity = analyser.aromaticity()
        instability_index = analyser.instability_index()
        isoelectric_point = analyser.isoelectric_point()
        
        # Amino acid composition (fraction 0.0 - 1.0, or percent)
        if callable(getattr(analyser, 'amino_acids_percent', None)):
            aa_comp = analyser.amino_acids_percent()
        elif hasattr(analyser, 'amino_acids_percent'):
            aa_comp = analyser.amino_acids_percent
        else:
            aa_comp = analyser.get_amino_acids_percent()
        
        annot = annot_map.get(prot_id, {})
        gene = annot.get('gene', rec.description.split('|')[0] if '|' in rec.description else 'UNKNOWN')
        temporal_class = annot.get('temporal_class', 'Unassigned')
        functional_cat = annot.get('functional_category', 'Uncharacterized')
        product = annot.get('product', rec.description)
        
        row_dict = {
            'protein_id': prot_id,
            'gene': gene,
            'product': product,
            'temporal_class': temporal_class,
            'functional_category': functional_cat,
            'length': length,
            'mw_kda': mw_kda,
            'isoelectric_point': isoelectric_point,
            'instability_index': instability_index,
            'aromaticity': aromaticity
        }
        
        # Add 20 amino acid fractions
        for aa in amino_acids:
            row_dict[f'aa_{aa}'] = aa_comp.get(aa, 0.0)
            
        feature_rows.append(row_dict)
        
    df_features = pd.DataFrame(feature_rows)
    df_features.to_csv(output_path, index=False)
    print(f"Saved {len(df_features)} protein feature profiles to {output_path}")
    print(f"Feature matrix dimensions: {len(df_features)} rows x {len([c for c in df_features.columns if c.startswith('aa_') or c in ['length', 'mw_kda', 'isoelectric_point', 'instability_index', 'aromaticity']])} numerical features.")
    
    return df_features

if __name__ == "__main__":
    compute_physicochemical_features()
