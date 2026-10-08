import os
import pytest
import pandas as pd

@pytest.fixture(scope="module")
def base_dir():
    return os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

@pytest.fixture(scope="module")
def supervised_dataset(base_dir):
    path = os.path.join(base_dir, "data", "processed", "final_temporal_supervised_dataset.csv")
    return pd.read_csv(path)

@pytest.fixture(scope="module")
def curated_corpus(base_dir):
    path = os.path.join(base_dir, "data", "processed", "final_curated_hsv_protein_corpus.csv")
    return pd.read_csv(path)

@pytest.fixture(scope="module")
def status_reconciliation(base_dir):
    path = os.path.join(base_dir, "results", "tables", "final_sequence_status_reconciliation.csv")
    return pd.read_csv(path)

@pytest.fixture(scope="module")
def temporal_summary(base_dir):
    path = os.path.join(base_dir, "results", "tables", "final_temporal_class_summary.csv")
    return pd.read_csv(path)

@pytest.fixture(scope="module")
def species_summary(base_dir):
    path = os.path.join(base_dir, "results", "tables", "final_species_temporal_summary.csv")
    return pd.read_csv(path)

@pytest.fixture(scope="module")
def gene_summary(base_dir):
    path = os.path.join(base_dir, "results", "tables", "final_gene_temporal_summary.csv")
    return pd.read_csv(path)

@pytest.fixture(scope="module")
def master_annotation(base_dir):
    path = os.path.join(base_dir, "data", "processed", "phase2_1_2A_final_annotation.csv")
    return pd.read_csv(path)

def test_canonical_count_and_id_uniqueness(master_annotation, status_reconciliation):
    """TEST 1 & 2: Verify canonical count is exactly 22,689 with unique canonical IDs."""
    assert len(master_annotation) == 22689
    assert master_annotation['sequence_id'].nunique() == 22689
    assert len(status_reconciliation) == 22689
    assert status_reconciliation['canonical_id'].nunique() == 22689

def test_no_new_or_missing_canonical_sequences(master_annotation, status_reconciliation):
    """TEST 3 & 4: Verify exact match of canonical sequence IDs across all reconciliation tables."""
    assert set(master_annotation['sequence_id']) == set(status_reconciliation['canonical_id'])

def test_supervised_dataset_temporal_classes_only(supervised_dataset):
    """TEST 5: Supervised dataset must contain ONLY IMMEDIATE_EARLY, EARLY, LATE."""
    allowed_classes = {'IMMEDIATE_EARLY', 'EARLY', 'LATE'}
    actual_classes = set(supervised_dataset['temporal_class'].unique())
    assert actual_classes.issubset(allowed_classes)
    assert len(supervised_dataset) == 16657

def test_unknown_never_enters_supervised_dataset(supervised_dataset):
    """TEST 6: UNKNOWN temporal class must never enter the supervised dataset."""
    assert 'UNKNOWN' not in supervised_dataset['temporal_class'].values

def test_conflicting_never_enters_supervised_dataset(supervised_dataset):
    """TEST 7: CONFLICTING temporal class must never enter the supervised dataset."""
    assert 'CONFLICTING' not in supervised_dataset['temporal_class'].values

def test_excluded_sequences_never_enter_supervised_dataset(supervised_dataset, status_reconciliation):
    """TEST 8: Excluded sequences must never enter the supervised dataset."""
    excluded_ids = set(status_reconciliation[status_reconciliation['final_status'] == 'EXCLUDED']['canonical_id'])
    assert len(excluded_ids) == 671
    assert not set(supervised_dataset['canonical_id']).intersection(excluded_ids)

def test_temporal_class_normalization(master_annotation, supervised_dataset):
    """TEST 9: Temporal classes must be normalized strictly to standard strings."""
    valid_classes = {'IMMEDIATE_EARLY', 'EARLY', 'LATE', 'UNKNOWN', 'CONFLICTING'}
    assert set(master_annotation['temporal_class'].unique()).issubset(valid_classes)
    assert set(supervised_dataset['temporal_class'].unique()).issubset({'IMMEDIATE_EARLY', 'EARLY', 'LATE'})

def test_species_validity(supervised_dataset):
    """TEST 10: Species in supervised dataset must be valid HSV species."""
    assert set(supervised_dataset['species'].unique()).issubset({'HSV-1', 'HSV-2'})

def test_temporal_evidence_presence_for_all_supervised_records(supervised_dataset):
    """TEST 11: Every supervised sequence must have explicit temporal evidence metadata."""
    assert supervised_dataset['temporal_evidence_source'].notna().all()
    assert supervised_dataset['temporal_evidence_type'].notna().all()
    assert (supervised_dataset['temporal_evidence_source'] != '').all()

def test_sequence_content_unchanged(master_annotation, supervised_dataset):
    """TEST 12: Sequence content and lengths match exactly between annotation and supervised dataset."""
    annot_map = dict(zip(master_annotation['sequence_id'], master_annotation['sequence']))
    for _, row in supervised_dataset.iterrows():
        assert annot_map[row['canonical_id']] == row['sequence']
        assert len(row['sequence']) == row['length']

def test_temporal_summary_closure(temporal_summary, master_annotation):
    """TEST 13: Summary counts match master dataset exactly."""
    total_seqs = temporal_summary['total_sequences'].sum()
    assert total_seqs == 22689
    ie_count = temporal_summary[temporal_summary['temporal_class'] == 'IMMEDIATE_EARLY']['total_sequences'].values[0]
    early_count = temporal_summary[temporal_summary['temporal_class'] == 'EARLY']['total_sequences'].values[0]
    late_count = temporal_summary[temporal_summary['temporal_class'] == 'LATE']['total_sequences'].values[0]
    assert ie_count == 1552
    assert early_count == 4140
    assert late_count == 10965

def test_species_temporal_counts_reconciliation(species_summary, master_annotation):
    """TEST 14: Species temporal cross-tabulation sums to HSV-1 (11,802) and HSV-2 (10,263)."""
    hsv1_row = species_summary[species_summary['virus_species'] == 'HSV-1']
    hsv2_row = species_summary[species_summary['virus_species'] == 'HSV-2']
    assert hsv1_row['TOTAL'].values[0] == 11802
    assert hsv2_row['TOTAL'].values[0] == 10263

def test_gene_temporal_counts_reconciliation(gene_summary, master_annotation):
    """TEST 15: Gene temporal table totals reconcile to 22,689."""
    assert gene_summary['TOTAL'].sum() == 22689

def test_final_status_reconciliation_completeness(status_reconciliation):
    """TEST 16: Final status reconciliation partitions all 22,689 records with zero orphans."""
    status_counts = status_reconciliation['final_status'].value_counts()
    assert status_counts['SUPERVISED_TEMPORAL'] == 16657
    assert status_counts['UNKNOWN_TEMPORAL'] == 5356
    assert status_counts['CONFLICTING_TEMPORAL'] == 5
    assert status_counts['EXCLUDED'] == 671
    assert status_counts.sum() == 22689
