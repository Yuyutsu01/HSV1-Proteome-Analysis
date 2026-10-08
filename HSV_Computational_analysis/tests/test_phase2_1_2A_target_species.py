import os
import pytest
import pandas as pd

@pytest.fixture(scope="module")
def base_dir():
    return os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

@pytest.fixture(scope="module")
def final_annotation(base_dir):
    path = os.path.join(base_dir, "data", "processed", "phase2_1_2A_final_annotation.csv")
    return pd.read_csv(path)

@pytest.fixture(scope="module")
def primary_ground_truth(base_dir):
    path = os.path.join(base_dir, "data", "processed", "phase2_1_2A_ground_truth_dataset.csv")
    return pd.read_csv(path)

@pytest.fixture(scope="module")
def sensitivity_dataset(base_dir):
    path = os.path.join(base_dir, "data", "processed", "phase2_1_2A_sensitivity_dataset.csv")
    return pd.read_csv(path)

@pytest.fixture(scope="module")
def evidence_matrix(base_dir):
    path = os.path.join(base_dir, "results", "tables", "phase2_1_2A_target_species_evidence.csv")
    return pd.read_csv(path)

@pytest.fixture(scope="module")
def gene_summary(base_dir):
    path = os.path.join(base_dir, "results", "tables", "phase2_1_2A_gene_temporal_evidence_summary.csv")
    return pd.read_csv(path)

def test_row_count_and_sequence_preservation(final_annotation):
    """TEST 1: Verify total row count remains exactly 22,689 with unique hashes."""
    assert len(final_annotation) == 22689
    assert final_annotation['sequence_id'].nunique() == 22689
    assert final_annotation['canonical_sequence_hash'].nunique() == 22689

def test_no_temporal_label_alteration(final_annotation):
    """TEST 2: Verify original temporal labels remain preserved."""
    assert (final_annotation['original_temporal_class'] == final_annotation['temporal_class']).all()

def test_evidence_presence_for_confirmed_records(final_annotation, evidence_matrix):
    """TEST 3: Verify every CONFIRMED_TARGET_SPECIES record has explicit evidence."""
    confirmed = final_annotation[final_annotation['target_species_validation_status'] == 'CONFIRMED_TARGET_SPECIES']
    assert (confirmed['target_species_evidence_source'] != 'NONE').all()
    assert (confirmed['target_species_evidence_identifier'] != 'NONE').all()
    assert len(evidence_matrix) == 28

def test_temporal_class_exact_agreement(final_annotation, evidence_matrix):
    """TEST 4: Verify confirmed records match the exact temporal class in the evidence matrix."""
    ev_dict = {row['gene_symbol']: row['temporal_class'] for _, row in evidence_matrix.iterrows()}
    cs_confirmed = final_annotation[
        (final_annotation['cross_species_transfer'] == 'YES') &
        (final_annotation['target_species_validation_status'] == 'CONFIRMED_TARGET_SPECIES')
    ]
    for _, row in cs_confirmed.iterrows():
        c_gene = row['canonical_gene']
        if c_gene in ev_dict:
            assert row['temporal_class'] == ev_dict[c_gene]

def test_no_orthology_in_primary_ground_truth(primary_ground_truth):
    """TEST 5: Verify no ORTHOLOGY_SUPPORTED_ONLY record enters primary ground truth."""
    assert (primary_ground_truth['cross_species_ground_truth_status'] == 'PRIMARY_GROUND_TRUTH').all()
    assert (primary_ground_truth['target_species_validation_status'] == 'CONFIRMED_TARGET_SPECIES').all()

def test_primary_and_sensitivity_partition_closure(primary_ground_truth, sensitivity_dataset, final_annotation):
    """TEST 6: Verify primary ground truth + sensitivity equals Tier A (10,482)."""
    tier_a = final_annotation[final_annotation['annotation_tier'] == 'TIER_A_HIGH_CONFIDENCE_GROUND_TRUTH']
    assert len(tier_a) == 10482
    assert len(primary_ground_truth) + len(sensitivity_dataset) == 10482
    assert len(primary_ground_truth) == 9097
    assert len(sensitivity_dataset) == 1385

def test_quarantined_and_conflicts_excluded_from_ground_truth(primary_ground_truth, final_annotation):
    """TEST 7: Verify conflicting or quarantined records cannot enter primary ground truth."""
    quarantined = final_annotation[final_annotation['cross_species_ground_truth_status'] == 'QUARANTINED']
    assert len(quarantined) == 5
    assert not set(quarantined['sequence_id']).intersection(set(primary_ground_truth['sequence_id']))

def test_gene_summary_kinetic_agreement(gene_summary):
    """TEST 8: Verify all 74 core genes have temporal kinetic agreement between HSV-1 and HSV-2."""
    assert len(gene_summary) == 74
    assert (gene_summary['kinetic_agreement'] == 'AGREEMENT').all()

def test_dataset_closure_across_all_statuses(final_annotation):
    """TEST 9: Verify total across all validation statuses equals 22,689."""
    counts = final_annotation['target_species_validation_status'].value_counts()
    assert counts.sum() == 22689
