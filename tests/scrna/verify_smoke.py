"""Verify real-algorithm smoke outputs, including reloadable trained models."""
import argparse
from pathlib import Path
import json
import re
import anndata as ad
import numpy as np
import pandas as pd


def verify(workdir, reload_models=False):
    root = Path(workdir)
    a = ad.read_h5ad(root / 'analysis/final.h5ad')
    assert a.n_obs > 0 and a.obs_names.is_unique and a.var_names.is_unique
    for layer in ['counts_raw', 'counts', 'counts_depthmatched']:
        if layer in a.layers:
            x = a.layers[layer]
            assert np.all(np.isfinite(x.data)) and np.all(x.data >= 0)
            np.testing.assert_allclose(x.data, np.rint(x.data))
    if 'counts_depthmatched' in a.layers:
        assert (a.layers['counts'] - a.layers['counts_depthmatched']).min() >= 0
    assert set(a.uns['scvi_hvg_names']).issubset(a.var_names)
    assert len(a.uns['scvi_hvg_names']) < a.n_vars  # full genes survive HVG modeling
    for embedding in ['X_scVI', 'X_umap_scVI', 'X_umap_uncorrected']:
        assert np.isfinite(a.obsm[embedding]).all()
    if 'scanvi_probabilities' in a.obsm:
        p = a.obsm['scanvi_probabilities']
        np.testing.assert_allclose(p.sum(axis=1), 1, atol=1e-5)
        assert np.all((a.obs['scanvi_normalized_entropy'] >= 0) & (a.obs['scanvi_normalized_entropy'] <= 1.00001))
        composition = pd.read_csv(root / 'analysis/composition.tsv', sep='\t')
        np.testing.assert_allclose(composition.groupby('sample')['expected_fraction'].sum(), 1, atol=1e-5)
    content = (root / 'report/report.html').read_text()
    assert content.startswith('<!doctype html>')
    assert '<html lang="en">' in content
    assert not re.search(r'[\u3040-\u30ff\u3400-\u9fff]', content)
    assert len(re.findall(r'<img ', content)) >= 9
    assert not re.search(r'<(?:script|img|link)[^>]+(?:src|href)=["\']https?://', content)
    assert 'Synthetic smoke test &lt;not biological data&gt;' in content
    for ext in ['pdf', 'svg', 'png']:
        assert len(list((root / 'report/figures').glob(f'*.{ext}'))) >= 9
    if reload_models:
        import scvi
        base = ad.read_h5ad(root / 'analysis/scvi.h5ad')
        hvg = base[:, list(base.uns['scvi_hvg_names'])].copy()
        model = scvi.model.SCVI.load(str(root / 'models/scvi'), adata=hvg)
        assert model.get_latent_representation().shape[0] == a.n_obs
        if 'scanvi_probabilities' in a.obsm:
            hvg = a[:, list(a.uns['scvi_hvg_names'])].copy()
            model = scvi.model.SCANVI.load(str(root / 'models/scanvi'), adata=hvg)
            np.testing.assert_allclose(model.predict(soft=True).sum(axis=1), 1, atol=1e-5)
    print(f'Smoke verification passed: {a.n_obs} cells, {a.n_vars} genes')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('workdir')
    parser.add_argument('--reload-models', action='store_true')
    args = parser.parse_args()
    verify(args.workdir, args.reload_models)
