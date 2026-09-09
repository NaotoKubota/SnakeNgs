import sys
sys.path.insert(0, snakemake.scriptdir)
import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import celltypist
from core import aggregate_symbols
from runtime import logged, versions, write_json

with logged(snakemake):
    cfg = snakemake.params.cfg
    p = cfg['celltypist']
    a = ad.read_h5ad(snakemake.input.adata)
    x, symbols = aggregate_symbols(a.layers[cfg['depth']['model_layer']], a.var['gene_symbol'])
    annotation = ad.AnnData(x.astype(np.float32), obs=a.obs.copy(), var=pd.DataFrame(index=symbols))
    sc.pp.normalize_total(annotation, target_sum=1e4)
    sc.pp.log1p(annotation)
    model = celltypist.models.Model.load(snakemake.input.model)
    overlap = len(symbols.intersection(model.features))
    if overlap < p['min_model_genes']:
        raise ValueError(f'Only {overlap} CellTypist genes overlap; check species, symbols and model')
    if p['majority_voting']:
        sc.tl.leiden(a, resolution=p['over_clustering_resolution'], key_added='celltypist_overcluster',
                     flavor='igraph', n_iterations=2, random_state=cfg['seed'])
    result = celltypist.annotate(annotation, model=model, majority_voting=p['majority_voting'],
        over_clustering=a.obs['celltypist_overcluster'].astype(str).to_numpy() if p['majority_voting'] else None,
        min_prop=p['min_prop'], mode='best match')
    labels = result.predicted_labels.reindex(a.obs_names)
    probabilities = result.probability_matrix.reindex(a.obs_names)
    a.obs['celltypist_label'] = labels['predicted_labels'].astype(str)
    a.obs['celltypist_majority'] = labels['majority_voting'].astype(str) if p['majority_voting'] else a.obs['celltypist_label']
    values = probabilities.to_numpy()
    if values.shape[1] < 2 or np.any(~np.isfinite(values)):
        raise ValueError('CellTypist must return finite scores for at least two labels')
    sorted_scores = np.sort(values, axis=1)
    a.obs['celltypist_probability'] = sorted_scores[:, -1]
    a.obs['celltypist_margin'] = sorted_scores[:, -1] - sorted_scores[:, -2]
    confident = (a.obs['celltypist_probability'] >= p['seed_probability']) & (a.obs['celltypist_margin'] >= p['seed_margin'])
    if p['majority_voting']:
        confident &= a.obs['celltypist_label'] == a.obs['celltypist_majority']
    unknown = cfg['scanvi']['unlabeled_category']
    if unknown in probabilities.columns:
        raise ValueError('scanvi.unlabeled_category conflicts with a CellTypist model class')
    seeds = a.obs['celltypist_label'].where(confident, unknown)
    counts = seeds[seeds != unknown].value_counts()
    supported = counts.index[counts >= cfg['scanvi']['min_seed_cells']]
    seeds = seeds.where(seeds.isin(supported), unknown)
    a.obs['scanvi_seed'] = pd.Categorical(seeds, categories=sorted(supported.tolist()) + [unknown])
    a.obsm['celltypist_probabilities'] = values.astype(np.float32)
    a.uns['celltypist_classes'] = probabilities.columns.to_numpy(dtype=str)
    a.write_h5ad(snakemake.output.adata, compression='gzip')
    write_json(snakemake.output.summary, {'model_gene_overlap': overlap,
        'seed_counts': seeds.value_counts().to_dict(), 'n_seed_classes': len(supported),
        'score_note': 'CellTypist one-vs-rest sigmoid scores do not sum to one.', 'versions': versions()})
