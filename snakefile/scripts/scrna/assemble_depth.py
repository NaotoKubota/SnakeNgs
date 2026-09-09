import sys
sys.path.insert(0, snakemake.scriptdir)
import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.io import mmread
from core import assign_global_cell_ids, depth_match, integer_counts
from runtime import logged, versions, write_json


with logged(snakemake):
    cfg = snakemake.params.cfg
    datasets, removed = [], []
    for i, sample in enumerate(snakemake.params.samples):
        a = ad.read_h5ad(snakemake.input.adata[i])
        genes = pd.read_csv(snakemake.input.genes[i], header=None, dtype=str)[0]
        barcodes = pd.read_csv(snakemake.input.barcodes[i], header=None, dtype=str)[0]
        if not np.array_equal(genes, a.var_names) or not np.array_equal(barcodes, a.obs_names):
            raise ValueError(f'{sample}: SoupX barcode/gene order differs from kb input')
        counts = integer_counts(mmread(snakemake.input.matrices[i]).T, 'SoupX output')
        if counts.shape != a.shape:
            raise ValueError(f'{sample}: SoupX shape differs from kb input')
        a.layers['counts'] = counts
        a.X = counts.copy()
        rho = pd.read_csv(snakemake.input.rhos[i], sep='\t', index_col=0)
        if not rho.index.equals(a.obs_names):
            raise ValueError(f'{sample}: SoupX rho barcodes differ from kb input')
        a.obs['soupx_rho'] = rho['rho'].to_numpy()
        a.obs['counts_corrected'] = np.asarray(counts.sum(axis=1)).ravel()
        # Preserve the raw library barcode in obs['barcode'], while using a
        # distinct index name for globally unique sample-prefixed cell IDs.
        assign_global_cell_ids(a, sample)
        zero = a.obs['counts_corrected'].to_numpy() == 0
        removed.extend({'cell_id': c, 'sample': sample, 'reason': 'zero_after_soupx'} for c in a.obs_names[zero])
        a = a[~zero].copy()
        if a.n_obs < cfg['qc']['min_cells_per_sample']:
            raise ValueError(f'{sample}: too few nonempty cells after SoupX')
        if datasets:
            if set(a.var_names) != set(datasets[0].var_names):
                raise ValueError('Sample gene sets differ; use one reference for all libraries')
            a = a[:, datasets[0].var_names].copy()
            if not np.array_equal(a.var['gene_symbol'], datasets[0].var['gene_symbol']):
                raise ValueError('Samples have conflicting gene-symbol mappings')
        datasets.append(a)
    combined = ad.concat(datasets, join='inner', merge='same')
    if not combined.obs_names.is_unique:
        raise ValueError('Cell IDs are not globally unique')
    for key in cfg['integration']['continuous_covariate_keys']:
        combined.obs[key] = pd.to_numeric(combined.obs[key], errors='raise')
    d = cfg['depth']
    matched, table = depth_match(combined.layers['counts'], combined.obs,
                                 d['enabled'], d['target_median'], d['groupby'], cfg['seed'])
    if d['enabled']:
        combined.layers['counts_depthmatched'] = matched
        combined.obs['counts_depthmatched'] = np.asarray(matched.sum(axis=1)).ravel()
    selected = combined.layers[d['model_layer']]
    zero = np.asarray(selected.sum(axis=1)).ravel() == 0
    removed.extend({'cell_id': c, 'sample': combined.obs.loc[c, 'sample'],
                    'reason': 'zero_after_depth_matching'} for c in combined.obs_names[zero])
    combined = combined[~zero].copy()
    retained = combined.obs.groupby('sample', observed=True).size()
    if any(retained.get(s, 0) < cfg['qc']['min_cells_per_sample'] for s in snakemake.params.samples):
        raise ValueError('Depth matching leaves too few cells in at least one sample')
    # All genes retained for CellTypist. scVI will train only on its separate HVG view.
    combined.X = combined.layers[d['model_layer']].astype(np.float32).copy()
    sc.pp.normalize_total(combined, target_sum=1e4)
    sc.pp.log1p(combined)
    combined.uns['count_contract'] = {
        'counts_raw': 'Unmodified kb UMI before ambient correction, retained cells',
        'counts': 'SoupX corrected integer UMI; raw kb if SoupX disabled',
        'counts_depthmatched': 'Binomial thinning of counts, only present if enabled',
        'X': 'log1p(counts / cell_total * 10000) using model_layer',
        'model_layer': d['model_layer'], 'source': 'kb',
    }
    table.to_csv(snakemake.output.depth, sep='\t', index=False)
    pd.DataFrame(removed, columns=['cell_id', 'sample', 'reason']).to_csv(snakemake.output.removed, sep='\t', index=False)
    combined.write_h5ad(snakemake.output.adata, compression='gzip')
    write_json(snakemake.output.summary, {'n_cells': combined.n_obs, 'n_genes': combined.n_vars,
        'removed_empty_cells': len(removed), 'model_layer': d['model_layer'], 'versions': versions()})
