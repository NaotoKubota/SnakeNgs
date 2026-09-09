import sys
sys.path.insert(0, snakemake.scriptdir)
import anndata as ad
import numpy as np
import scanpy as sc
import scvi
from core import batch_entropy, integer_counts
from model_utils import configure_training, save_history
from runtime import logged, versions, write_json


with logged(snakemake):
    cfg = snakemake.params.cfg
    p = cfg['integration']
    configure_training(cfg, snakemake.threads)
    a = ad.read_h5ad(snakemake.input[0])
    layer = cfg['depth']['model_layer']
    a.layers[layer] = integer_counts(a.layers[layer], 'scVI input')
    expressed = np.asarray((a.layers[layer] > 0).sum(axis=0)).ravel() >= 3
    candidate = a[:, expressed].copy()
    if candidate.n_vars < 3:
        raise ValueError('Fewer than three genes detected in three cells')
    # Seurat v3 HVGs are selected from integer counts per batch.
    if candidate.n_vars <= p['n_hvg']:
        candidate.var['highly_variable'] = True
    else:
        sc.pp.highly_variable_genes(candidate, layer=layer, flavor='seurat_v3',
                                    n_top_genes=p['n_hvg'], batch_key=p['batch_key'])
    a.var['highly_variable'] = a.var_names.isin(candidate.var_names[candidate.var['highly_variable']])
    hvg = a[:, a.var['highly_variable']].copy()
    sc.pp.pca(hvg, n_comps=min(50, hvg.n_obs - 1, hvg.n_vars - 1),
              svd_solver='arpack', random_state=cfg['seed'])
    a.obsm['X_pca_uncorrected'] = hvg.obsm['X_pca'].copy()
    sc.pp.neighbors(a, use_rep='X_pca_uncorrected', n_neighbors=p['n_neighbors'],
                     key_added='uncorrected', random_state=cfg['seed'])
    sc.tl.umap(a, neighbors_key='uncorrected', random_state=cfg['seed'])
    a.obsm['X_umap_uncorrected'] = a.obsm['X_umap'].copy()
    scvi.model.SCVI.setup_anndata(hvg, layer=layer, batch_key=p['batch_key'],
        categorical_covariate_keys=p['categorical_covariate_keys'] or None,
        continuous_covariate_keys=p['continuous_covariate_keys'] or None)
    model = scvi.model.SCVI(hvg, n_latent=p['n_latent'], n_layers=p['n_layers'], gene_likelihood='nb')
    model.train(max_epochs=p['max_epochs'], accelerator=p['accelerator'], devices=p['devices'],
                 batch_size=p['batch_size'], early_stopping=False)
    a.obsm['X_scVI'] = model.get_latent_representation()
    sc.pp.neighbors(a, use_rep='X_scVI', n_neighbors=p['n_neighbors'], random_state=cfg['seed'])
    sc.tl.umap(a, random_state=cfg['seed'])
    a.obsm['X_umap_scVI'] = a.obsm['X_umap'].copy()
    sc.tl.leiden(a, resolution=p['leiden_resolution'], flavor='igraph',
                 n_iterations=2, key_added='leiden_scVI', random_state=cfg['seed'])
    a.obs['batch_entropy_before'] = batch_entropy(a.obsm['X_pca_uncorrected'], a.obs[p['batch_key']], p['n_neighbors'])
    a.obs['batch_entropy_scVI'] = batch_entropy(a.obsm['X_scVI'], a.obs[p['batch_key']], p['n_neighbors'])
    a.uns['scvi_hvg_names'] = hvg.var_names.to_numpy(dtype=str)
    model.save(snakemake.output.model, overwrite=True, save_anndata=False)
    save_history(model, snakemake.output.history)
    a.write_h5ad(snakemake.output.adata, compression='gzip')
    write_json(snakemake.output.summary, {'n_hvg': hvg.n_vars, 'n_cells': a.n_obs,
        'batch_key': p['batch_key'], 'n_batches': a.obs[p['batch_key']].nunique(),
        'model_layer': layer, 'versions': versions()})
