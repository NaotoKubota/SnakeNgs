import sys
sys.path.insert(0, snakemake.scriptdir)
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
from core import (aggregate_class_probabilities, composition_intervals,
                  map_major_cell_types, uncertainty, validate_probabilities)
from model_utils import configure_training, save_history
from runtime import logged, versions, write_json


with logged(snakemake):
    cfg = snakemake.params.cfg
    p = cfg['scanvi']
    annotation = cfg['annotation']
    major_map = annotation['major_celltype_map']
    major_default = annotation['unmapped_major_celltype']
    a = ad.read_h5ad(snakemake.input.adata)
    if p['enabled']:
        configure_training(cfg, snakemake.threads)
        hvg = a[:, list(a.uns['scvi_hvg_names'])].copy()
        seed_classes = sorted(set(hvg.obs['scanvi_seed'].astype(str)) - {p['unlabeled_category']})
        if len(seed_classes) < 2:
            raise ValueError('scANVI needs >=2 sufficiently supported seed classes. Inspect CellTypist model/thresholds; disable scanvi explicitly to retain CellTypist-only results.')
        pretrained = scvi.model.SCVI.load(snakemake.input.model, adata=hvg)
        model = scvi.model.SCANVI.from_scvi_model(pretrained, adata=hvg,
                  unlabeled_category=p['unlabeled_category'], labels_key='scanvi_seed')
        model.train(max_epochs=p['max_epochs'], n_samples_per_label=p['n_samples_per_label'],
                    accelerator=cfg['integration']['accelerator'], devices=cfg['integration']['devices'],
                    batch_size=cfg['integration']['batch_size'], early_stopping=False)
        # Integrate latent-posterior uncertainty with independent stochastic draws.
        total, squares, conditional_entropy = None, None, None
        for _ in range(p['posterior_samples']):
            prediction = model.predict(soft=True, use_posterior_mean=False)
            if not isinstance(prediction, pd.DataFrame) or not prediction.index.equals(a.obs_names):
                raise ValueError('Unexpected scANVI prediction index/type')
            q = validate_probabilities(prediction.to_numpy())
            if total is None:
                classes = prediction.columns.to_numpy(dtype=str)
                total = np.zeros_like(q)
                squares = np.zeros_like(q)
                conditional_entropy = np.zeros(q.shape[0])
            elif not np.array_equal(classes, prediction.columns):
                raise ValueError('scANVI class ordering changed between posterior draws')
            total += q
            squares += q * q
            conditional_entropy += uncertainty(q)['entropy']
        probabilities = total / p['posterior_samples']
        sd = np.sqrt(np.maximum(0, squares / p['posterior_samples'] - probabilities ** 2))
        metrics = uncertainty(probabilities)
        for name, values in metrics.items():
            a.obs['scanvi_' + name] = values
        a.obs['scanvi_latent_information'] = np.maximum(0, metrics['entropy'] - conditional_entropy / p['posterior_samples'])
        a.obsm['scanvi_probabilities'] = probabilities.astype(np.float32)
        a.obsm['scanvi_probability_sd'] = sd.astype(np.float32)
        a.uns['scanvi_classes'] = classes
        a.obs['scanvi_label'] = pd.Categorical(classes[probabilities.argmax(axis=1)])
        confident = (metrics['max_probability'] >= p['min_probability']) & (metrics['normalized_entropy'] <= p['max_normalized_entropy'])
        a.obs['annotation_uncertain'] = ~confident
        a.obs['cell_type'] = a.obs['scanvi_label'].astype(str).where(confident, p['unlabeled_category'])
        a.obs['scanvi_major_label'] = pd.Categorical(map_major_cell_types(
            a.obs['scanvi_label'].astype(str), major_map, major_default, p['unlabeled_category']))
        a.obs['major_cell_type'] = pd.Categorical(map_major_cell_types(
            a.obs['cell_type'].astype(str), major_map, major_default, p['unlabeled_category']))
        a.obsm['X_scANVI'] = model.get_latent_representation()
        # Preserve the scVI graph as well as its UMAP before adding the scANVI graph.
        sc.pp.neighbors(a, use_rep='X_scANVI', n_neighbors=cfg['integration']['n_neighbors'],
                         key_added='scanvi', random_state=cfg['seed'])
        sc.tl.umap(a, neighbors_key='scanvi', random_state=cfg['seed'])
        a.obsm['X_umap_scANVI'] = a.obsm['X_umap'].copy()
        composition = composition_intervals(probabilities, a.obs, classes, p['composition_draws'], cfg['seed'])
        major_probabilities, major_classes = aggregate_class_probabilities(
            probabilities, classes, major_map, major_default)
        composition_major = composition_intervals(
            major_probabilities, a.obs, major_classes, p['composition_draws'], cfg['seed'])
        model.save(snakemake.output.model, overwrite=True, save_anndata=False)
        save_history(model, snakemake.output.history)
        summary = {'enabled': True, 'n_classes': len(classes), 'classes': classes.tolist(),
            'uncertain_cells': int((~confident).sum()), 'posterior_samples': p['posterior_samples'],
            'composition_draws': p['composition_draws'],
            'major_classes': major_classes.tolist(),
            'scope': 'Latent posterior and conditional label uncertainty, fixed model weights. Not calibrated error probabilities or biological replicate confidence intervals.'}
    else:
        a.obs['cell_type'] = a.obs['celltypist_majority'].astype(str)
        a.obs['annotation_uncertain'] = a.obs['scanvi_seed'].astype(str) == p['unlabeled_category']
        a.obs['cell_type'] = a.obs['cell_type'].where(~a.obs['annotation_uncertain'], p['unlabeled_category'])
        a.obs['major_cell_type'] = pd.Categorical(map_major_cell_types(
            a.obs['cell_type'].astype(str), major_map, major_default, p['unlabeled_category']))
        composition = a.obs.groupby(['sample', 'cell_type'], observed=True).size().rename('hard_count').reset_index()
        composition['n_cells'] = composition['sample'].map(a.obs.groupby('sample', observed=True).size())
        composition['expected_count'] = composition['hard_count']
        composition['expected_fraction'] = composition['hard_count'] / composition['n_cells']
        composition['label_interval_low'] = np.nan
        composition['label_interval_high'] = np.nan
        composition['interval_scope'] = 'not_computed_scanvi_disabled'
        composition_major = a.obs.groupby(['sample', 'major_cell_type'], observed=True).size().rename('hard_count').reset_index()
        composition_major = composition_major.rename(columns={'major_cell_type': 'cell_type'})
        composition_major['n_cells'] = composition_major['sample'].map(a.obs.groupby('sample', observed=True).size())
        composition_major['expected_count'] = composition_major['hard_count']
        composition_major['expected_fraction'] = composition_major['hard_count'] / composition_major['n_cells']
        composition_major['label_interval_low'] = np.nan
        composition_major['label_interval_high'] = np.nan
        composition_major['interval_scope'] = 'not_computed_scanvi_disabled'
        Path(snakemake.output.model).mkdir(parents=True, exist_ok=True)
        write_json(Path(snakemake.output.model) / 'disabled.json', {'enabled': False})
        pd.DataFrame(columns=['step', 'metric', 'column', 'value']).to_csv(snakemake.output.history, sep='\t', index=False)
        summary = {'enabled': False, 'uncertain_cells': int(a.obs['annotation_uncertain'].sum())}
    for key in ['donor', 'condition', 'batch']:
        metadata = a.obs[['sample', key]].drop_duplicates().set_index('sample')[key]
        composition[key] = composition['sample'].map(metadata)
        composition_major[key] = composition_major['sample'].map(metadata)
    a.uns['uncertainty_scope'] = 'Conditional on reference labels, observed cells and fitted weights; excludes donor variation, mapping, SoupX and reference-model uncertainty.'
    a.obs.to_csv(snakemake.output.cells, sep='\t', index_label='cell_id')
    composition.to_csv(snakemake.output.composition, sep='\t', index=False)
    composition_major.to_csv(snakemake.output.composition_major, sep='\t', index=False)
    a.write_h5ad(snakemake.output.adata, compression='gzip')
    write_json(snakemake.output.summary, dict(summary, versions=versions()))
