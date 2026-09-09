import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse
from core import (aggregate_symbols, assign_global_cell_ids, batch_entropy,
                  called_count_reconciliation, composition_intervals, depth_match,
                  integer_counts, mad_flags, sample_seed, thin_counts, uncertainty)


@pytest.mark.parametrize('bad', [[[-1, 2]], [[np.nan, 2]], [[.3, 1]], [[np.inf, 0]]])
def test_non_count_input_rejected(bad):
    with pytest.raises(ValueError, match='integer'):
        integer_counts(bad)


def test_thinning_preserves_sparse_integer_contract_and_endpoints():
    x = sparse.csr_matrix([[0, 12, 100], [3, 0, 40]])
    y = thin_counts(x, .5, 42)
    assert sparse.issparse(y) and y.dtype == np.int64
    assert (x - y).min() >= 0
    np.testing.assert_array_equal(y.toarray(), thin_counts(x, .5, 42).toarray())
    np.testing.assert_array_equal(x.toarray(), thin_counts(x, 1, 42).toarray())
    assert thin_counts(x, 0, 42).nnz == 0
    with pytest.raises(ValueError): thin_counts(x, 1.1, 42)


def test_kb_filtered_recount_difference_is_audited_not_rejected():
    filtered = sparse.csr_matrix([[3, 0], [0, 4]])
    unfiltered_called = sparse.csr_matrix([[2, 0], [0, 5]])
    result = called_count_reconciliation(filtered, unfiltered_called)
    assert result['discordant_entries'] == 2
    assert result['discordant_barcodes'] == 2
    assert result['filtered_minus_unfiltered_umi'] == 0
    assert result['absolute_umi_difference'] == 2
    with pytest.raises(ValueError, match='different shapes'):
        called_count_reconciliation(filtered, sparse.csr_matrix([[1, 2]]))


def test_sample_seed_is_deterministic_and_signed_int32_safe():
    seed = sample_seed(42, 'Ctrl_P19')
    assert seed == 1054802820
    assert 0 <= seed <= np.iinfo(np.int32).max
    assert seed == sample_seed(42, 'Ctrl_P19')
    assert seed != sample_seed(42, 'KO_P19')


def test_global_cell_ids_preserve_barcode_and_roundtrip(tmp_path):
    barcodes = pd.Index(['AAAC', 'TTTG'], name='barcode')
    a = ad.AnnData(sparse.csr_matrix([[1], [2]]),
                   obs=pd.DataFrame({'barcode': barcodes}, index=barcodes))
    assign_global_cell_ids(a, 'sample1')
    assert a.obs_names.name == 'cell_id'
    assert list(a.obs_names) == ['sample1:AAAC', 'sample1:TTTG']
    assert list(a.obs['barcode']) == ['AAAC', 'TTTG']
    path = tmp_path / 'roundtrip.h5ad'
    a.write_h5ad(path)
    restored = ad.read_h5ad(path)
    assert restored.obs_names.name == 'cell_id'
    assert list(restored.obs['barcode']) == ['AAAC', 'TTTG']


def test_depth_matching_grouped_reproducible_and_row_order_preserved():
    obs = pd.DataFrame({'sample': ['b', 'a', 'd', 'c'] * 200,
                        'age': ['P19', 'P19', 'P28', 'P28'] * 200})
    x = sparse.csr_matrix(np.array([[2000, 0], [200, 0], [0, 1000], [0, 500]] * 200))
    matched, stats = depth_match(x, obs, groupby=['age'], seed=42)
    np.testing.assert_array_equal(matched.toarray(), depth_match(x, obs, groupby=['age'], seed=42)[0].toarray())
    assert (x - matched).min() >= 0
    assert stats.set_index('sample').loc['b', 'keep_probability'] == .1
    assert stats.set_index('sample').loc['d', 'keep_probability'] == .5
    assert abs(stats.set_index('sample').loc['b', 'median_after'] - 200) < 10
    np.testing.assert_array_equal(matched.toarray()[1::4], x.toarray()[1::4])
    assert not matched.toarray()[0::4, 1].any()
    with pytest.raises(ValueError, match='shallowest'):
        depth_match(x, obs, target=600, groupby=['age'])


def test_duplicate_symbols_sum_counts_without_loss():
    x = sparse.csr_matrix([[2, 3, 4], [0, 1, 7]])
    y, symbols = aggregate_symbols(x, ['A', 'A', 'B'])
    assert list(symbols) == ['A', 'B']
    np.testing.assert_array_equal(y.toarray(), [[5, 4], [1, 7]])
    np.testing.assert_array_equal(x.sum(axis=1), y.sum(axis=1))


def test_entropy_and_intervals_keep_uncertain_probability_mass():
    p = np.array([[1, 0], [.5, .5], [0, 1]])
    u = uncertainty(p)
    np.testing.assert_allclose(u['normalized_entropy'], [0, 1, 0], atol=1e-10)
    obs = pd.DataFrame({'sample': ['a'] * 3})
    result = composition_intervals(p, obs, ['A', 'B'], 1000, 42)
    np.testing.assert_allclose(result['expected_count'], [1.5, 1.5])
    assert result['expected_fraction'].sum() == 1
    assert np.all(result['label_interval_low'] <= result['label_interval_high'])
    pd.testing.assert_frame_equal(result, composition_intervals(p, obs, ['A', 'B'], 1000, 42))
    with pytest.raises(ValueError, match='sum to one'):
        uncertainty(np.array([[.8, .8]]))


def test_mad_ties_and_neighbor_entropy_single_batch():
    assert not mad_flags(np.ones(10), 5)[0].any()
    assert batch_entropy(np.ones((10, 2)), ['a'] * 10, 3).shape == (10,)
    assert np.isnan(batch_entropy(np.ones((10, 2)), ['a'] * 10, 3)).all()
    # With two observations, neighbors exclude self.
    np.testing.assert_allclose(batch_entropy(np.array([[0], [1]]), ['a', 'b'], 1), [0, 0])
