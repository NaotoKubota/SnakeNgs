"""Sparse count contracts, depth matching, and uncertainty summaries."""
import hashlib
import numpy as np
import pandas as pd
from scipy import sparse


def integer_counts(matrix, name='counts'):
    x = sparse.csr_matrix(matrix, copy=True)
    x.sum_duplicates()
    if np.any(~np.isfinite(x.data)) or np.any(x.data < 0) or np.any(np.abs(x.data - np.rint(x.data)) > 1e-6):
        raise ValueError(f'{name}: expected nonnegative integer UMI counts; scaled/log data are invalid')
    x.data = np.rint(x.data).astype(np.int64)
    x.eliminate_zeros()
    return x


def called_count_reconciliation(filtered, unfiltered_called):
    """Summarize expected recount differences after kb barcode filtering.

    kb's bustools filter creates an allowlist and corrects nearby barcode
    records to that allowlist before recounting UMIs.  Consequently, its
    filtered matrix is not necessarily an exact row subset of the unfiltered
    matrix, even though every called barcode must occur in both.
    """
    filtered = integer_counts(filtered, 'filtered called counts')
    unfiltered_called = integer_counts(unfiltered_called, 'unfiltered called counts')
    if filtered.shape != unfiltered_called.shape:
        raise ValueError('Filtered and unfiltered called-count matrices have different shapes')
    delta = filtered - unfiltered_called
    delta.eliminate_zeros()
    filtered_umi = int(filtered.sum())
    unfiltered_umi = int(unfiltered_called.sum())
    absolute = int(np.abs(delta.data).sum())
    return {
        'discordant_entries': int(delta.nnz),
        'discordant_barcodes': int(np.count_nonzero(delta.getnnz(axis=1))),
        'filtered_umi': filtered_umi,
        'unfiltered_called_umi': unfiltered_umi,
        'filtered_minus_unfiltered_umi': filtered_umi - unfiltered_umi,
        'absolute_umi_difference': absolute,
        'relative_absolute_difference': absolute / max(1, filtered_umi),
    }


def sample_seed(seed, name):
    # Annoy (used by Scrublet's approximate-neighbor implementation) accepts
    # only a signed 32-bit seed. Keep one deterministic seed contract that is
    # safe for Annoy, NumPy, Scanpy, and igraph.
    modulus = np.iinfo(np.int32).max + 1
    return (int(seed) + int.from_bytes(hashlib.sha256(str(name).encode()).digest()[:4], 'little')) % modulus


def assign_global_cell_ids(adata, sample, barcode_key='barcode'):
    """Prefix cell IDs while retaining the original library barcode safely."""
    barcodes = pd.Index(adata.obs_names.astype(str))
    if not barcodes.is_unique:
        raise ValueError(f'{sample}: duplicate cell barcodes')
    if barcode_key in adata.obs:
        stored = adata.obs[barcode_key].astype(str).to_numpy()
        if not np.array_equal(stored, barcodes.to_numpy()):
            raise ValueError(f'{sample}: {barcode_key} column differs from input cell IDs')
    else:
        adata.obs[barcode_key] = barcodes.to_numpy()
    adata.obs_names = pd.Index([f'{sample}:{barcode}' for barcode in barcodes], name='cell_id')


def mad_flags(values, n_mads, upper=False):
    v = np.asarray(values, dtype=float)
    median = float(np.median(v))
    mad = float(np.median(np.abs(v - median)))
    # Constant/mostly tied distributions provide no useful MAD cutoff.
    flags = np.zeros(len(v), dtype=bool) if mad == 0 else v < median - n_mads * mad
    if upper and mad > 0:
        flags |= v > median + n_mads * mad
    return flags, {'median': median, 'mad': mad, 'lower': median - n_mads * mad,
                   'upper': median + n_mads * mad if upper else None}


def thin_counts(matrix, probability, seed):
    if not 0 <= probability <= 1:
        raise ValueError('Thinning probability must be in [0, 1]; upsampling is not supported')
    x = integer_counts(matrix)
    rng = np.random.default_rng(seed)
    x.data = rng.binomial(x.data, probability).astype(np.int64)
    x.eliminate_zeros()
    return x


def depth_match(matrix, obs, enabled=True, target=None, groupby=(), seed=0):
    x = integer_counts(matrix)
    totals = np.asarray(x.sum(axis=1)).ravel()
    samples = obs['sample'].astype(str).to_numpy()
    records = []
    probabilities = {}
    groupby = list(groupby)
    groups = obs.groupby(groupby, observed=True, sort=True).indices.values() if groupby else [np.arange(len(obs))]
    for indices in groups:
        medians = {s: float(np.median(totals[indices][samples[indices] == s])) for s in sorted(set(samples[indices]))}
        if min(medians.values()) <= 0:
            raise ValueError('Depth matching encountered a sample with zero median counts')
        goal = min(medians.values()) if target is None else float(target)
        if enabled and goal > min(medians.values()):
            raise ValueError('depth.target_median exceeds the shallowest sample in a matching group')
        for sample, median in medians.items():
            if sample in probabilities:
                raise ValueError('A sample belongs to multiple depth-matching groups')
            probabilities[sample] = goal / median if enabled else 1.0
            records.append({'sample': sample, 'median_before': median,
                            'target_median': goal if enabled else median,
                            'keep_probability': probabilities[sample], 'enabled': enabled})
    # Preserve row order regardless of metadata or sample ordering.
    matched = x.copy()
    for sample, probability in probabilities.items():
        rows = np.flatnonzero(samples == sample)
        positions = np.concatenate([np.arange(x.indptr[i], x.indptr[i + 1]) for i in rows])
        rng = np.random.default_rng(sample_seed(seed, sample))
        matched.data[positions] = rng.binomial(x.data[positions], probability)
    matched.eliminate_zeros()
    after = np.asarray(matched.sum(axis=1)).ravel()
    for record in records:
        values = after[samples == record['sample']]
        record['median_after'] = float(np.median(values))
        record['zero_cells_after'] = int(np.sum(values == 0))
    return matched, pd.DataFrame(records)


def aggregate_symbols(matrix, symbols):
    """Aggregate duplicate symbols before total-normalization; never discard UMIs."""
    symbols = pd.Index(symbols).astype(str)
    codes, names = pd.factorize(symbols, sort=False)
    mapping = sparse.csr_matrix((np.ones(len(codes)), (np.arange(len(codes)), codes)),
                                shape=(len(codes), len(names)))
    return (sparse.csr_matrix(matrix) @ mapping).tocsr(), pd.Index(names)


def validate_probabilities(p, min_classes=2):
    p = np.asarray(p, dtype=float)
    # Float32 class aggregation can exceed [0, 1] by a few ulps. Accept only
    # numerical round-off, clip it, and continue to enforce row normalization.
    tolerance = 1e-6
    if (p.ndim != 2 or p.shape[1] < min_classes or np.any(~np.isfinite(p))
            or np.any(p < -tolerance) or np.any(p > 1 + tolerance)):
        raise ValueError(f'Expected a finite cell-by-class probability matrix with at least {min_classes} classes')
    p = np.clip(p, 0, 1)
    if not np.allclose(p.sum(axis=1), 1, atol=1e-5):
        raise ValueError('Class probabilities must sum to one')
    return p / p.sum(axis=1, keepdims=True)


def uncertainty(p):
    p = validate_probabilities(p)
    entropy = -(p * np.log(np.clip(p, 1e-30, 1))).sum(axis=1)
    ordered = np.sort(p, axis=1)
    return {'max_probability': p.max(axis=1), 'entropy': entropy,
            'normalized_entropy': entropy / np.log(p.shape[1]),
            'margin': ordered[:, -1] - ordered[:, -2]}


def composition_intervals(probabilities, obs, classes, draws, seed):
    """Independent conditional label draws; these intervals exclude donor variation."""
    p = validate_probabilities(probabilities, min_classes=1)
    records = []
    for sample in sorted(obs['sample'].astype(str).unique()):
        rows = np.flatnonzero(obs['sample'].astype(str).to_numpy() == sample)
        q = p[rows]
        rng = np.random.default_rng(sample_seed(seed, sample))
        cumulative = q.cumsum(axis=1)
        cumulative[:, -1] = 1.0
        simulations = np.empty((draws, len(classes)), dtype=float)
        for draw in range(draws):
            labels = (rng.random((len(rows), 1)) > cumulative).sum(axis=1)
            simulations[draw] = np.bincount(labels, minlength=len(classes)) / len(rows)
        lo, hi = np.quantile(simulations, [0.025, 0.975], axis=0)
        hard = np.bincount(q.argmax(axis=1), minlength=len(classes))
        for k, label in enumerate(classes):
            records.append({'sample': sample, 'cell_type': label, 'n_cells': len(rows),
                            'hard_count': int(hard[k]), 'expected_count': float(q[:, k].sum()),
                            'expected_fraction': float(q[:, k].mean()),
                            'label_interval_low': float(lo[k]), 'label_interval_high': float(hi[k]),
                            'interval_scope': 'conditional_label_uncertainty_only'})
    return pd.DataFrame(records)


def map_major_cell_types(labels, mapping, default, unknown='Unknown'):
    """Map subtype labels to major cell types while preserving the unknown label."""
    return np.asarray([unknown if str(label) == unknown else mapping.get(str(label), default)
                       for label in labels], dtype=str)


def aggregate_class_probabilities(probabilities, classes, mapping, default):
    """Collapse subtype probabilities into mutually exclusive major cell types."""
    p = validate_probabilities(probabilities)
    classes = np.asarray(classes, dtype=str)
    if p.shape[1] != len(classes):
        raise ValueError('Probability columns and class names differ in length')
    mapped = np.asarray([mapping.get(label, default) for label in classes], dtype=str)
    major_classes = np.asarray(sorted(set(mapped)), dtype=str)
    major = np.column_stack([p[:, mapped == label].sum(axis=1) for label in major_classes])
    return validate_probabilities(major, min_classes=1), major_classes


def batch_entropy(representation, batches, k):
    from sklearn.neighbors import NearestNeighbors
    labels, names = pd.factorize(np.asarray(batches))
    if len(names) <= 1:
        return np.full(len(labels), np.nan)
    k = min(k, len(labels) - 1)
    # X=None excludes each training point itself, including tied observations.
    idx = NearestNeighbors(n_neighbors=k).fit(representation).kneighbors(return_distance=False)
    p = np.column_stack([(labels[idx] == b).mean(axis=1) for b in range(len(names))])
    return -(p * np.log(np.clip(p, 1e-30, 1))).sum(axis=1) / np.log(len(names))
