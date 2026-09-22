"""
Phase 1: Dataset Acquisition and Exact Redundancy Removal for HSV-1 (NC_001806.2)

Biological Context:
Herpes Simplex Virus Type 1 (HSV-1) strain 17 contains terminal repeats (TR_L, TR_S)
and internal inverted repeats (IR_L, IR_S). Consequently, certain genes located in repeat
regions (e.g., ICP0/RL2, ICP4/RS1, ICP34.5/RL1) are duplicated in the raw GenBank/RefSeq record.
Exact sequence deduplication ensures that statistical analyses and machine learning models
are not biased by identical duplicate proteins.

Author: Computational Biology Pipeline
"""

import os
import sys
import urllib.request
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def ensure_directories():
    """Create project directories if they do not exist."""
    dirs = ['data', 'processed_data', 'embeddings', 'results', 'figures', 'scripts', 'notebooks']
    for d in dirs:
        os.makedirs(d, exist_ok=True)
    print("Project directories verified.")

def fetch_refseq_genbank(accession="NC_001806.2", output_path="data/NC_001806.2.gbk"):
    """
    Download the complete RefSeq GenBank record for HSV-1 strain 17 from NCBI.
    """
    if os.path.exists(output_path) and os.path.getsize(output_path) > 10000:
        print(f"GenBank record already exists at {output_path}")
        return output_path
    
    print(f"Fetching {accession} from NCBI Entrez...")
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={accession}&rettype=gbwithparts&retmode=text"
    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64)'})
    with urllib.request.urlopen(req) as response, open(output_path, 'wb') as out_file:
        out_file.write(response.read())
    print(f"Downloaded {accession} to {output_path} ({os.path.getsize(output_path)} bytes)")
    return output_path

# Standardized functional & temporal dictionary curated from NCBI, UniProt (UP000009294) and Roizman et al. (Fields Virology)
# Gene -> (Temporal Class, Functional Category, Full Protein Name)
HSV1_GENE_ANNOTATIONS = {
    # Immediate Early (IE / alpha)
    'RL2': ('Immediate-Early', 'Regulatory/Transactivator', 'Infected cell protein 0 (ICP0 / IE110)'),
    'RS1': ('Immediate-Early', 'Regulatory/Transactivator', 'Infected cell protein 4 (ICP4 / IE175)'),
    'US1': ('Immediate-Early', 'Regulatory/Transactivator', 'Infected cell protein 22 (ICP22 / IE68)'),
    'US12': ('Immediate-Early', 'Immune Evasion', 'Infected cell protein 47 (ICP47 / TAP inhibitor)'),
    'UL54': ('Immediate-Early', 'RNA Export/Processing', 'Infected cell protein 27 (ICP27 / IE63)'),
    
    # Early (E / beta) - Replication, Nucleotide metabolism, Repair
    'UL2': ('Early', 'DNA Replication/Repair', 'Uracil-DNA glycosylase (UNG)'),
    'UL5': ('Early', 'DNA Replication/Replication Fork', 'Helicase-primase helicase subunit'),
    'UL8': ('Early', 'DNA Replication/Replication Fork', 'Helicase-primase accessory subunit'),
    'UL9': ('Early', 'DNA Replication/Origin Binding', 'Origin-binding protein (OBP)'),
    'UL12': ('Early', 'DNA Recombination/Processing', 'Deoxyribonuclease (Alkaline nuclease)'),
    'UL23': ('Early', 'Nucleotide Metabolism', 'Thymidine kinase (TK)'),
    'UL29': ('Early', 'DNA Replication/SSB', 'Single-stranded DNA-binding protein (ICP8)'),
    'UL30': ('Early', 'DNA Replication/Polymerase', 'DNA polymerase catalytic subunit (Pol)'),
    'UL39': ('Early', 'Nucleotide Metabolism', 'Ribonucleotide reductase large subunit (ICP6 / RR1)'),
    'UL40': ('Early', 'Nucleotide Metabolism', 'Ribonucleotide reductase small subunit (RR2)'),
    'UL42': ('Early', 'DNA Replication/Processivity', 'DNA polymerase processivity factor (Pol accessory)'),
    'UL50': ('Early', 'Nucleotide Metabolism', 'dUTP pyrophosphatase (dUTPase)'),
    'UL52': ('Early', 'DNA Replication/Replication Fork', 'Helicase-primase primase subunit'),
    'UL53': ('Late', 'Glycoprotein/Entry & Fusion', 'Glycoprotein K (gK)'), # gK is often gamma
    'UL13': ('Early', 'Protein Kinase', 'Serine/threonine-protein kinase UL13'),
    'US3': ('Early', 'Protein Kinase/Anti-apoptotic', 'Serine/threonine-protein kinase US3'),
    'RL1': ('Late', 'Neurovirulence/Host Defense', 'Infected cell protein 34.5 (ICP34.5 / gamma-34.5)'),
    
    # Late (L / gamma) - Capsid, Tegument, Glycoproteins, Packaging, Assembly
    'UL1': ('Late', 'Glycoprotein/Entry', 'Glycoprotein L (gL)'),
    'UL3': ('Late', 'Nuclear Egress', 'Nuclear egress protein UL3'),
    'UL4': ('Late', 'Tegument', 'Tegument protein UL4'),
    'UL6': ('Late', 'Capsid Portal', 'Capsid portal protein UL6'),
    'UL7': ('Late', 'Tegument', 'Tegument protein UL7'),
    'UL10': ('Late', 'Glycoprotein/Entry', 'Glycoprotein M (gM)'),
    'UL11': ('Late', 'Tegument/Secondary Envelopment', 'Myristylated tegument protein UL11'),
    'UL14': ('Late', 'Tegument', 'Tegument protein UL14'),
    'UL15': ('Late', 'DNA Packaging/Terminase', 'DNA packaging terminase subunit 1'),
    'UL16': ('Late', 'Tegument', 'Tegument protein UL16'),
    'UL17': ('Late', 'Capsid Vertex/Packaging', 'Capsid vertex-specific component UL17'),
    'UL18': ('Late', 'Capsid/Triplex', 'Capsid triplex subunit VP23'),
    'UL19': ('Late', 'Capsid/Major', 'Major capsid protein VP5'),
    'UL20': ('Late', 'Envelope/Egress', 'Envelope protein UL20'),
    'UL21': ('Late', 'Tegument', 'Tegument protein UL21'),
    'UL22': ('Late', 'Glycoprotein/Entry', 'Glycoprotein H (gH)'),
    'UL24': ('Late', 'Pathogenicity/Egress', 'Protein UL24'),
    'UL25': ('Late', 'Capsid Vertex/Packaging', 'Capsid vertex-specific component UL25'),
    'UL26': ('Late', 'Capsid Scaffolding/Protease', 'Capsid maturation protease / VP24'),
    'UL26.5': ('Late', 'Capsid Scaffolding', 'Small capsomer-interacting protein / VP22a scaffold'),
    'UL27': ('Late', 'Glycoprotein/Entry & Fusion', 'Glycoprotein B (gB)'),
    'UL28': ('Late', 'DNA Packaging/Terminase', 'DNA packaging terminase subunit 2'),
    'UL31': ('Late', 'Nuclear Egress', 'Nuclear egress lamina protein UL31'),
    'UL32': ('Late', 'DNA Packaging', 'Envelope/packaging protein UL32'),
    'UL33': ('Late', 'DNA Packaging/Terminase', 'DNA packaging terminase subunit 3'),
    'UL34': ('Late', 'Nuclear Egress', 'Inner nuclear membrane protein UL34'),
    'UL35': ('Late', 'Capsid/Small Capsid', 'Capsid protein VP26'),
    'UL36': ('Late', 'Tegument/Large Tegument', 'Large inner tegument protein (VP1/2 / UL36)'),
    'UL37': ('Late', 'Tegument/Inner Tegument', 'Inner tegument protein UL37'),
    'UL38': ('Late', 'Capsid/Triplex', 'Capsid triplex subunit VP19C'),
    'UL41': ('Late', 'Host Shutoff', 'Virion host shutoff protein (vhs / UL41)'),
    'UL43': ('Late', 'Envelope/Membrane', 'Membrane protein UL43'),
    'UL44': ('Late', 'Glycoprotein/Attachment', 'Glycoprotein C (gC)'),
    'UL45': ('Late', 'Envelope', 'Envelope protein UL45'),
    'UL46': ('Late', 'Tegument', 'Tegument protein VP11/12'),
    'UL47': ('Late', 'Tegument', 'Tegument protein VP13/14'),
    'UL48': ('Late', 'Tegument/Transactivator', 'Tegument transactivator protein (VP16 / alpha-TIF)'),
    'UL49': ('Late', 'Tegument/Major', 'Major tegument protein (VP22)'),
    'UL49A': ('Late', 'Glycoprotein/Entry', 'Glycoprotein N (gN)'),
    'UL51': ('Late', 'Tegument/Secondary Envelopment', 'Tegument protein UL51'),
    'UL55': ('Late', 'Nuclear Egress', 'Nuclear egress auxiliary protein UL55'),
    'UL56': ('Late', 'Pathogenicity/Egress', 'Pathogenicity protein UL56'),
    'US2': ('Late', 'Tegument/Immune Modulation', 'Tegument protein US2'),
    'US4': ('Late', 'Glycoprotein/Entry', 'Glycoprotein G (gG)'),
    'US5': ('Late', 'Glycoprotein/Immune Evasion', 'Glycoprotein J (gJ)'),
    'US6': ('Late', 'Glycoprotein/Primary Receptor Binding', 'Glycoprotein D (gD)'),
    'US7': ('Late', 'Glycoprotein/Cell-to-Cell Spread', 'Glycoprotein I (gI)'),
    'US8': ('Late', 'Glycoprotein/Cell-to-Cell Spread', 'Glycoprotein E (gE)'),
    'US8A': ('Late', 'Tegument/Membrane', 'Protein US8A'),
    'US9': ('Late', 'Tegument/Anterograde Transport', 'Envelope protein US9'),
    'US10': ('Late', 'Tegument', 'Tegument protein US10'),
    'US11': ('Late', 'RNA-binding/PKR Antagonist', 'Tegument RNA-binding protein (US11)'),
    'LAT': ('Latency', 'Non-coding/Exon', 'Latency-associated transcript protein ORF'),
}

def parse_and_deduplicate(gbk_path="data/NC_001806.2.gbk"):
    """
    Parse CDS features from GenBank record, extract protein sequence,
    and deduplicate by exact amino acid sequence.
    """
    record = SeqIO.read(gbk_path, "genbank")
    raw_proteins = []
    
    for feature in record.features:
        if feature.type == "CDS":
            qualifiers = feature.qualifiers
            translation = qualifiers.get('translation', [None])[0]
            if not translation:
                continue
            
            protein_id = qualifiers.get('protein_id', ['UNKNOWN'])[0]
            gene = qualifiers.get('gene', [qualifiers.get('locus_tag', ['UNKNOWN'])[0]])[0]
            product = qualifiers.get('product', ['Uncharacterized protein'])[0]
            locus_tag = qualifiers.get('locus_tag', [''])[0]
            location = str(feature.location)
            
            raw_proteins.append({
                'protein_id': protein_id,
                'gene': gene,
                'locus_tag': locus_tag,
                'product': product,
                'sequence': translation.upper(),
                'length': len(translation),
                'location': location
            })
    
    raw_df = pd.DataFrame(raw_proteins)
    print(f"Extracted {len(raw_df)} total raw CDS records from GenBank.")
    
    # Save raw fasta
    raw_records = [
        SeqRecord(Seq(row['sequence']), id=row['protein_id'], description=f"{row['gene']}|{row['product']}")
        for _, row in raw_df.iterrows()
    ]
    SeqIO.write(raw_records, "data/raw_proteins.fasta", "fasta")
    
    # Sequence Deduplication (Exact match)
    seen_sequences = {}
    dedup_proteins = []
    duplicates_removed = []
    
    for _, row in raw_df.iterrows():
        seq = row['sequence']
        if seq not in seen_sequences:
            seen_sequences[seq] = row['gene']
            dedup_proteins.append(row)
        else:
            duplicates_removed.append((row['gene'], row['protein_id'], seen_sequences[seq]))
    
    dedup_df = pd.DataFrame(dedup_proteins)
    print(f"Non-redundant unique protein sequences: {len(dedup_df)}")
    print(f"Duplicates removed ({len(duplicates_removed)}): {duplicates_removed}")
    
    # Map curated annotations
    temporal_classes = []
    functional_cats = []
    detailed_names = []
    
    for _, row in dedup_df.iterrows():
        gene = row['gene']
        # Try exact or clean gene name
        clean_gene = gene.replace('HSV1gp', 'UL').strip()
        if clean_gene in HSV1_GENE_ANNOTATIONS:
            t_class, f_cat, d_name = HSV1_GENE_ANNOTATIONS[clean_gene]
        elif gene in HSV1_GENE_ANNOTATIONS:
            t_class, f_cat, d_name = HSV1_GENE_ANNOTATIONS[gene]
        else:
            # Infer from product name
            prod = row['product'].lower()
            if 'glycoprotein' in prod or 'membrane' in prod or 'envelope' in prod:
                t_class = 'Late'
                f_cat = 'Glycoprotein/Membrane'
            elif 'capsid' in prod or 'tegument' in prod:
                t_class = 'Late'
                f_cat = 'Structural/Capsid-Tegument'
            elif 'polymerase' in prod or 'kinase' in prod or 'helicase' in prod or 'uracil' in prod:
                t_class = 'Early'
                f_cat = 'Replication/Enzyme'
            elif 'infected cell protein 0' in prod or 'infected cell protein 4' in prod or 'icp27' in prod:
                t_class = 'Immediate-Early'
                f_cat = 'Regulatory/Transactivator'
            else:
                t_class = 'Unassigned'
                f_cat = 'Other/Uncharacterized'
            d_name = row['product']
            
        temporal_classes.append(t_class)
        functional_cats.append(f_cat)
        detailed_names.append(d_name)
        
    dedup_df['temporal_class'] = temporal_classes
    dedup_df['functional_category'] = functional_cats
    dedup_df['detailed_name'] = detailed_names
    
    # Check ambiguous residues
    ambiguous_counts = []
    for seq in dedup_df['sequence']:
        ambig = sum(1 for aa in seq if aa not in "ACDEFGHIKLMNPQRSTVWY")
        ambiguous_counts.append(ambig)
    dedup_df['ambiguous_residues'] = ambiguous_counts
    
    # Save non-redundant fasta
    nr_records = [
        SeqRecord(Seq(row['sequence']), id=row['protein_id'], description=f"{row['gene']}|{row['temporal_class']}|{row['functional_category']}")
        for _, row in dedup_df.iterrows()
    ]
    SeqIO.write(nr_records, "processed_data/non_redundant.fasta", "fasta")
    
    # Save dataset statistics
    stats = {
        'Metric': [
            'GenBank Accession',
            'Raw CDS count',
            'Non-redundant unique proteins retained',
            'Duplicate sequences removed',
            'Min sequence length (aa)',
            'Max sequence length (aa)',
            'Mean sequence length (aa)',
            'Median sequence length (aa)',
            'Proteins with ambiguous residues',
            'Proteins with assigned temporal class',
            'Proteins with assigned functional category'
        ],
        'Value': [
            'NC_001806.2 (HSV-1 strain 17)',
            str(len(raw_df)),
            str(len(dedup_df)),
            f"{len(duplicates_removed)} (e.g. repeat copies of RL1/ICP34.5, RL2/ICP0, RS1/ICP4)",
            str(dedup_df['length'].min()),
            str(dedup_df['length'].max()),
            f"{dedup_df['length'].mean():.2f}",
            f"{dedup_df['length'].median():.1f}",
            str(sum(dedup_df['ambiguous_residues'] > 0)),
            f"{sum(dedup_df['temporal_class'] != 'Unassigned')}/{len(dedup_df)}",
            f"{sum(dedup_df['functional_category'] != 'Other/Uncharacterized')}/{len(dedup_df)}"
        ]
    }
    stats_df = pd.DataFrame(stats)
    stats_df.to_csv("processed_data/dataset_statistics.csv", index=False)
    dedup_df.to_csv("processed_data/curated_annotations.csv", index=False)
    
    print("\nDataset Statistics Summary:")
    print(stats_df.to_string(index=False))
    return dedup_df, stats_df

if __name__ == "__main__":
    ensure_directories()
    fetch_refseq_genbank()
    parse_and_deduplicate()
