import os
import pytest
import pandas as pd

@pytest.fixture(scope="module")
def base_dir():
    return os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

@pytest.fixture(scope="module")
def final_annotation(base_dir):
    path = os.path.join(base_dir, "data", "processed", "phase2_1_2_final_annotation.csv")
    return pd.read_csv(path)

@pytest.fixture(scope="module")
def ground_truth(base_dir):
    path = os.path.join(base_dir, "data", "processed", "phase2_1_2_ground_truth_dataset.csv")
    return pd.read_csv(path)

@pytest.fixture(scope="module")
def extended_corpus(base_dir):
    path = os.path.join(base_dir, "data", "processed", "phase2_1_2_extended_labeled_corpus.csv")
    return pd.read_csv(path)

@pytest.fixture(scope="module")
def cross_species_audit(base_dir):
    path = os.path.join(base_dir, "results", "tables", "phase2_1_2_cross_species_audit.csv")
    return pd.read_csv(path)

@pytest.fixture(scope="module")
def transfer_matrix(base_dir):
    path = os.path.join(base_dir, "results", "tables", "phase2_1_2_transfer_matrix.csv")
    return pd.read_csv(path)

def test_row_count_preservation(final_annotation):
    """TEST 1: Verify total row count remains exactly 22,689."""
    assert len(final_annotation) == 22689

def test_canonical_id_preservation(final_annotation):
    """TEST 2: Verify canonical IDs are unique and well-formatted."""
    assert final_annotation['sequence_id'].nunique() == 22689
    assert final_annotation['canonical_sequence_hash'].nunique() == 22689

def test_no_temporal_label_alteration(final_annotation):
    """TEST 3: Verify original temporal labels match validated temporal labels."""
    assert (final_annotation['original_temporal_class'] == final_annotation['temporal_class']).all()
    assert (final_annotation['validated_temporal_class'] == final_annotation['temporal_class']).all()

def test_cross_species_transfer_population(final_annotation, cross_species_audit):
    """TEST 4: Verify cross-species transfer count is exactly 7,656."""
    cs_records = final_annotation[final_annotation['cross_species_transfer'] == 'YES']
    assert len(cs_records) == 7656
    assert len(cross_species_audit) == 7656

def test_cross_species_validity_categories(final_annotation):
    """TEST 5: Verify allowed values for cross_species_validity_category."""
    allowed = {
        'VALIDATED_TARGET_SPECIES', 'SUPPORTED_ORTHOLOGY_INFERENCE',
        'INSUFFICIENT_TRANSFER_SUPPORT', 'CONFLICTING_TRANSFER',
        'MANUAL_REVIEW_REQUIRED', 'UNKNOWN'
    }
    actual = set(final_annotation['cross_species_validity_category'].unique())
    assert actual.issubset(allowed)

def test_cross_species_ground_truth_statuses(final_annotation):
    """TEST 6: Verify allowed values for cross_species_ground_truth_status."""
    allowed = {
        'PRIMARY_GROUND_TRUTH', 'SECONDARY_SENSITIVITY_ANALYSIS',
        'EXTENDED_CORPUS_ONLY', 'QUARANTINED', 'EXCLUDED_FROM_TEMPORAL_MODELING'
    }
    actual = set(final_annotation['cross_species_ground_truth_status'].unique())
    assert actual.issubset(allowed)

def test_ground_truth_dataset_integrity(ground_truth, final_annotation):
    """TEST 7: Verify ground truth dataset contains exactly 10,482 Tier A records."""
    assert len(ground_truth) == 10482
    tier_a_in_final = final_annotation[final_annotation['annotation_tier'] == 'TIER_A_HIGH_CONFIDENCE_GROUND_TRUTH']
    assert len(tier_a_in_final) == 10482
    assert set(ground_truth['sequence_id']) == set(tier_a_in_final['sequence_id'])

def test_primary_and_sensitivity_partitions(ground_truth):
    """TEST 8: Verify primary ground truth + secondary sensitivity sums to 10,482."""
    n_prim = (ground_truth['cross_species_ground_truth_status'] == 'PRIMARY_GROUND_TRUTH').sum()
    n_sec = (ground_truth['cross_species_ground_truth_status'] == 'SECONDARY_SENSITIVITY_ANALYSIS').sum()
    assert n_prim + n_sec == 10482
    assert n_prim == 9105
    assert n_sec == 1377

def test_quarantined_exclusion_from_ground_truth(ground_truth, final_annotation):
    """TEST 9: Verify quarantined (Tier C) records cannot enter ground truth."""
    quarantined = final_annotation[final_annotation['cross_species_ground_truth_status'] == 'QUARANTINED']
    assert len(quarantined) == 5
    assert not set(quarantined['sequence_id']).intersection(set(ground_truth['sequence_id']))

def test_extended_corpus_integrity(extended_corpus, final_annotation):
    """TEST 10: Verify extended corpus contains exactly 16,657 temporally labeled sequences."""
    assert len(extended_corpus) == 16657
    labeled_in_final = final_annotation[final_annotation['temporal_class'].isin(['IMMEDIATE_EARLY', 'EARLY', 'LATE'])]
    assert len(labeled_in_final) == 16657
    assert set(extended_corpus['sequence_id']) == set(labeled_in_final['sequence_id'])

def test_transfer_matrix_closure(transfer_matrix):
    """TEST 11: Verify transfer matrix totals close to 7,656."""
    hsv1_to_hsv2 = transfer_matrix[transfer_matrix['source_species'] == 'HSV-1']['TOTAL'].values[0]
    hsv2_to_hsv1 = transfer_matrix[transfer_matrix['source_species'] == 'HSV-2']['TOTAL'].values[0]
    assert hsv1_to_hsv2 == 7656
    assert hsv2_to_hsv1 == 0

def test_dataset_closure_across_all_statuses(final_annotation):
    """TEST 12: Verify total count across all ground-truth statuses equals 22,689."""
    counts = final_annotation['cross_species_ground_truth_status'].value_counts()
    total = counts.sum()
    assert total == 22689
