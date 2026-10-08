import os
import pytest
import hashlib
import pandas as pd

@pytest.fixture(scope="module")
def base_dir():
    return os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

@pytest.fixture(scope="module")
def accounting_table(base_dir):
    path = os.path.join(base_dir, "results", "tables", "phase2_5_dataset_accounting.csv")
    return pd.read_csv(path)

@pytest.fixture(scope="module")
def provenance_coverage(base_dir):
    path = os.path.join(base_dir, "results", "tables", "phase2_5_provenance_coverage.csv")
    return pd.read_csv(path)

@pytest.fixture(scope="module")
def completeness_summary(base_dir):
    path = os.path.join(base_dir, "results", "tables", "phase2_5_completeness_summary.csv")
    return pd.read_csv(path)

@pytest.fixture(scope="module")
def source_db_summary(base_dir):
    path = os.path.join(base_dir, "results", "tables", "phase2_5_source_database_summary.csv")
    return pd.read_csv(path)

@pytest.fixture(scope="module")
def authenticity_table(base_dir):
    path = os.path.join(base_dir, "results", "tables", "phase2_5_sequence_authenticity_verification.csv")
    return pd.read_csv(path)

@pytest.fixture(scope="module")
def supervised_validation(base_dir):
    path = os.path.join(base_dir, "results", "tables", "phase2_5_supervised_dataset_validation.csv")
    return pd.read_csv(path)

@pytest.fixture(scope="module")
def master_annotation(base_dir):
    path = os.path.join(base_dir, "data", "processed", "phase2_1_2A_final_annotation.csv")
    return pd.read_csv(path)

@pytest.fixture(scope="module")
def supervised_dataset(base_dir):
    path = os.path.join(base_dir, "data", "processed", "final_temporal_supervised_dataset.csv")
    return pd.read_csv(path)

def test_raw_record_and_canonical_count_reconciliation(base_dir):
    """TEST 1: Verify raw FASTA records (22,710), canonical sequences (22,689), and duplicates (21)."""
    def read_fasta_raw(filepath):
        records = []
        header = None
        seq_lines = []
        with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    if header is not None:
                        records.append((header, ''.join(seq_lines)))
                    header = line
                    seq_lines = []
                else:
                    seq_lines.append(line)
            if header is not None:
                records.append((header, ''.join(seq_lines)))
        return records

    raw1_path = os.path.join(base_dir, "raw", "HSV2_non_redundant (2).fasta")
    raw2_path = os.path.join(base_dir, "raw", "non_redundant.fasta")

    raw1 = read_fasta_raw(raw1_path)
    raw2 = read_fasta_raw(raw2_path)
    all_raw = raw1 + raw2

    assert len(raw1) == 10540
    assert len(raw2) == 12170
    assert len(all_raw) == 22710

    unique_seqs = {s for h, s in all_raw}
    assert len(unique_seqs) == 22689
    assert len(all_raw) - len(unique_seqs) == 21

def test_raw_to_canonical_and_canonical_to_raw_mapping(base_dir, master_annotation):
    """TEST 2: Verify bidirectional mapping between raw FASTA records and canonical sequences."""
    def read_fasta_raw(filepath):
        seqs = set()
        seq_lines = []
        with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    if seq_lines:
                        seqs.add(''.join(seq_lines))
                    seq_lines = []
                else:
                    seq_lines.append(line)
            if seq_lines:
                seqs.add(''.join(seq_lines))
        return seqs

    raw1 = read_fasta_raw(os.path.join(base_dir, "raw", "HSV2_non_redundant (2).fasta"))
    raw2 = read_fasta_raw(os.path.join(base_dir, "raw", "non_redundant.fasta"))
    all_raw_seqs = raw1.union(raw2)

    master_seqs = set(master_annotation['sequence'])
    assert all_raw_seqs == master_seqs
    assert len(master_seqs) == 22689

def test_no_sequence_modification_and_hash_presence(authenticity_table):
    """TEST 3: Verify 100% exact SHA-256 match, non-null local and external hashes."""
    assert (authenticity_table['exact_sequence_match'] == True).all()
    assert (authenticity_table['local_sequence_hash'].notna()).all()
    assert (authenticity_table['external_sequence_hash'].notna()).all()
    assert (authenticity_table['local_sequence_hash'] == authenticity_table['external_sequence_hash']).all()

def test_completeness_summary_accounting(completeness_summary):
    """TEST 4: Verify completeness summary table accounts for all 22,689 sequences."""
    assert completeness_summary['total_count'].sum() == 22689
    comp_dict = dict(zip(completeness_summary['completeness_status'], completeness_summary['total_count']))
    assert comp_dict['FULL_LENGTH'] == 20462
    assert comp_dict['PARTIAL_FRAGMENT'] == 1957
    assert comp_dict['PDB_FRAGMENT'] == 148
    assert comp_dict['ENGINEERED_FRAGMENT'] == 55
    assert comp_dict['UNKNOWN'] == 67

def test_source_database_summary_accounting(source_db_summary):
    """TEST 5: Verify source database summary accounts for all 22,689 sequences with 100% retrieval."""
    assert source_db_summary['total_records'].sum() == 22689
    assert source_db_summary['external_sequence_retrieved'].sum() == 22689
    assert source_db_summary['exact_matches'].sum() == 22689
    assert source_db_summary['unavailable'].sum() == 0
    assert source_db_summary['mismatches'].sum() == 0

def test_supervised_subset_integrity(supervised_validation, supervised_dataset):
    """TEST 6: Supervised dataset contains exactly 16,657 valid records."""
    assert len(supervised_validation) == 16657
    assert len(supervised_dataset) == 16657
    assert (supervised_validation['validation_verdict'] == 'VALID').all()

def test_temporal_class_validity_in_supervised(supervised_dataset):
    """TEST 7: Supervised temporal classes must be strictly IE, Early, or Late."""
    allowed = {'IMMEDIATE_EARLY', 'EARLY', 'LATE'}
    assert set(supervised_dataset['temporal_class'].unique()).issubset(allowed)
    assert 'UNKNOWN' not in supervised_dataset['temporal_class'].values
    assert 'CONFLICTING' not in supervised_dataset['temporal_class'].values

def test_provenance_field_consistency(provenance_coverage):
    """TEST 8: Verify 100% traceable accessions in provenance coverage."""
    assert len(provenance_coverage) == 22689
    assert (provenance_coverage['provenance_status'] == 'TRACEABLE').all()
    assert (provenance_coverage['accession'] != '').all()

def test_accounting_partition_closure(accounting_table):
    """TEST 9: Verify exact accounting closure with zero unaccounted sequences."""
    cat_map = dict(zip(accounting_table['category'], accounting_table['count']))
    assert cat_map['SUPERVISED'] == 16657
    assert cat_map['UNKNOWN_ELIGIBLE'] == 5356
    assert cat_map['CONFLICTING'] == 5
    assert cat_map['EXCLUDED'] == 671
    assert cat_map['TOTAL_CANONICAL'] == 22689
    assert cat_map['UNACCOUNTED'] == 0
    assert cat_map['SUPERVISED'] + cat_map['UNKNOWN_ELIGIBLE'] + cat_map['CONFLICTING'] + cat_map['EXCLUDED'] == 22689

def test_canonical_id_uniqueness_and_hash_reproducibility(master_annotation):
    """TEST 10: Verify canonical ID uniqueness and SHA-256 prefix reproducibility."""
    assert master_annotation['sequence_id'].nunique() == 22689
    for _, row in master_annotation.iterrows():
        seq = row['sequence']
        calc_hash = hashlib.sha256(seq.encode('utf-8')).hexdigest()
        assert row['sequence_id'].endswith(calc_hash[:12])
