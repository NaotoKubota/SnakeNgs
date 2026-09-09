import sys
sys.path.insert(0, snakemake.scriptdir)
import json
from pathlib import Path
import re
import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.io import mmwrite
from core import called_count_reconciliation, integer_counts, mad_flags, sample_seed
from runtime import logged, versions, write_json


def read_counts(path, cfg, t2g):
    source = ad.read_h5ad(path)
    layer = cfg['count_layer']
    if layer == 'auto':
        if cfg['workflow'] == 'nac':
            required = ['mature', 'nascent', 'ambiguous']
            if any(k not in source.layers for k in required):
                raise ValueError(f'{path}: nac auto counts require layers {required}; set count_layer explicitly for other layouts')
            matrix = sum((integer_counts(source.layers[k], k) for k in required))
        else:
            matrix = source.X
    else:
        matrix = source.X if layer == 'X' else source.layers[layer]
    x = integer_counts(matrix, f'{path}:{layer}')
    if not source.obs_names.is_unique or not source.var_names.is_unique:
        raise ValueError(f'{path}: duplicate barcodes or gene IDs')
    result = ad.AnnData(x, obs=pd.DataFrame(index=source.obs_names.copy()), var=source.var.copy())
    column = cfg['gene_symbol_column']
    if column in source.var:
        symbols = source.var[column].astype(str)
    elif t2g:
        mapping = pd.read_csv(t2g[0], sep='\t', header=None, dtype=str)
        if mapping.shape[1] < 3:
            raise ValueError('t2g requires transcript, gene ID, gene symbol columns')
        genes = mapping[[1, 2]].drop_duplicates()
        if genes[1].duplicated().any():
            raise ValueError('t2g assigns conflicting symbols to the same gene ID')
        symbols = pd.Series(result.var_names.map(genes.set_index(1)[2]), index=result.var_names).fillna('')
    else:
        # Reused kb --gene-names H5ADs may carry symbols directly as var_names.
        if result.var_names.str.match(r'^ENS[A-Z]*G\d+').any():
            raise ValueError('Gene IDs need kb.gene_symbol_column or reference.t2g for QC/CellTypist')
        symbols = pd.Series(result.var_names, index=result.var_names)
    result.var['gene_symbol'] = [s if s and s != 'nan' else g for s, g in zip(symbols, result.var_names)]
    return result


def run(sm):
    cfg = sm.params.cfg
    sample = sm.wildcards.sample
    seed = sample_seed(cfg['seed'], sample)
    sc.settings.n_jobs = sm.threads
    cells = read_counts(sm.input.filtered, cfg['kb'], sm.input.t2g)
    droplets = read_counts(sm.input.unfiltered, cfg['kb'], sm.input.t2g)
    if set(cells.var_names) != set(droplets.var_names):
        raise ValueError('Filtered and unfiltered kb gene sets differ')
    droplets = droplets[:, cells.var_names].copy()
    if not cells.obs_names.isin(droplets.obs_names).all():
        raise ValueError('Called kb barcodes are absent from unfiltered kb matrix')
    # kb's bustools filter corrects nearby barcode records to the called-cell
    # allowlist and then recounts UMIs. The filtered matrix is therefore the
    # authoritative cell matrix, but need not be an exact subset of the
    # unfiltered droplet matrix. Record that expected difference for auditing.
    reconciliation = called_count_reconciliation(cells.X, droplets[cells.obs_names].X)
    raw_totals = np.asarray(droplets.X.sum(axis=1)).ravel()
    rank = np.argsort(-raw_totals, kind='stable')
    pd.DataFrame({'rank': np.arange(1, len(rank) + 1), 'total_counts': raw_totals[rank],
                  'kb_called': droplets.obs_names[rank].isin(cells.obs_names)}).to_csv(sm.output.rank, sep='\t', index=False)
    q = cfg['qc']
    for name, regex in [('mt', q['mt_regex']), ('hb', q['hb_regex'])]:
        cells.var[name] = cells.var['gene_symbol'].str.contains(regex, regex=True).to_numpy()
    sc.pp.calculate_qc_metrics(cells, qc_vars=['mt', 'hb'], percent_top=None, log1p=False, inplace=True)
    cells.obs['barcode'] = cells.obs_names
    excluded = {'R1', 'R2', 'kb_filtered', 'kb_unfiltered', 'cellranger_metrics'}
    for key, value in sm.params.metadata.items():
        if key not in excluded:
            cells.obs[key] = value
    reasons = pd.DataFrame(index=cells.obs_names)
    reasons['low_genes'] = cells.obs['n_genes_by_counts'] < q['min_genes']
    reasons['low_counts'] = cells.obs['total_counts'] < q['min_counts']
    reasons['high_mt'] = cells.obs['pct_counts_mt'] > q['max_pct_mt']
    reasons['high_hb'] = cells.obs['pct_counts_hb'] > q['max_pct_hb']
    for metric, key in [('n_genes_by_counts', 'max_genes'), ('total_counts', 'max_counts')]:
        if q[key] is not None:
            reasons[key] = cells.obs[metric] > q[key]
    thresholds = {}
    for metric in ['total_counts', 'n_genes_by_counts']:
        reasons['mad_' + metric], thresholds[metric] = mad_flags(np.log1p(cells.obs[metric]), q['mad_n'], q['mad_upper'])
    cells.obs['qc_pass'] = ~reasons.any(axis=1)
    cells.obs['qc_reason'] = reasons.apply(lambda r: ';'.join(r.index[r]), axis=1)
    cells.obs['doublet_score'] = np.nan
    cells.obs['predicted_doublet'] = False
    qc = cells[cells.obs['qc_pass']].copy()
    if qc.n_obs < q['min_cells_per_sample']:
        raise ValueError(f'{sample}: only {qc.n_obs} cells pass QC; review sample-specific thresholds')
    doublet_threshold = None
    if cfg['scrublet']['enabled']:
        import scrublet as scr
        pars = cfg['scrublet']
        model = scr.Scrublet(qc.X, expected_doublet_rate=pars['expected_doublet_rate'],
                            sim_doublet_ratio=pars['sim_doublet_ratio'], random_state=seed)
        scores, predicted = model.scrub_doublets(
            n_prin_comps=min(pars['n_prin_comps'], qc.n_obs - 2, qc.n_vars - 1),
            min_counts=pars['min_counts'], min_cells=pars['min_cells'],
            min_gene_variability_pctl=pars['min_gene_variability_pctl'], verbose=True)
        if pars['threshold'] is not None:
            predicted = model.call_doublets(threshold=pars['threshold'])
        if predicted is None or not hasattr(model, 'threshold_') or not np.all(np.isfinite(scores)):
            raise ValueError('Scrublet did not find a valid threshold; inspect log and set scrublet.threshold')
        doublet_threshold = float(model.threshold_)
        cells.obs.loc[qc.obs_names, 'doublet_score'] = scores
        cells.obs.loc[qc.obs_names, 'predicted_doublet'] = predicted.astype(bool)
    cells.obs['retained'] = cells.obs['qc_pass'] & ~cells.obs['predicted_doublet']
    cells.obs.to_csv(sm.output.cells, sep='\t', index_label='cell_id')
    kept = cells[cells.obs['retained']].copy()
    if kept.n_obs < q['min_cells_per_sample']:
        raise ValueError(f'{sample}: too few singlets after doublet removal ({kept.n_obs})')
    kept.layers['counts_raw'] = kept.X.copy()
    kept.obs['soupx_cluster'] = '0'
    if cfg['soupx']['enabled']:
        provisional = kept.copy()
        sc.pp.normalize_total(provisional, target_sum=1e4)
        sc.pp.log1p(provisional)
        n_comps = min(30, provisional.n_obs - 1, provisional.n_vars - 1)
        sc.pp.pca(provisional, n_comps=n_comps, svd_solver='arpack', random_state=seed)
        sc.pp.neighbors(provisional, n_neighbors=min(15, provisional.n_obs - 1), random_state=seed)
        sc.tl.leiden(provisional, resolution=cfg['soupx']['clusters_resolution'],
                     flavor='igraph', n_iterations=2, random_state=seed)
        kept.obs['soupx_cluster'] = provisional.obs['leiden'].astype(str)
    out = Path(sm.output.soupx)
    out.mkdir(parents=True, exist_ok=True)
    mmwrite(out / 'toc.mtx', kept.X.T)
    pd.Series(kept.var_names).to_csv(out / 'genes.tsv', header=False, index=False)
    pd.Series(kept.obs_names).to_csv(out / 'barcodes.tsv', header=False, index=False)
    kept.obs[['soupx_cluster']].to_csv(out / 'clusters.tsv', sep='\t', index_label='barcode')
    lo, hi = cfg['soupx']['soup_range']
    empty = (raw_totals >= lo) & (raw_totals <= hi) & ~droplets.obs_names.isin(cells.obs_names)
    profile = np.asarray(droplets.X[empty].sum(axis=0)).ravel()
    if cfg['soupx']['enabled'] and (empty.sum() < cfg['soupx']['min_empty_droplets'] or profile.sum() == 0):
        raise ValueError(f'{sample}: insufficient kb empty droplets ({empty.sum()}); review soupx.soup_range')
    pd.DataFrame({'gene': kept.var_names, 'counts': profile,
                  'est': profile / max(1, profile.sum())}).to_csv(out / 'soup_profile.tsv', sep='\t', index=False)
    write_json(out / 'settings.json', {'enabled': cfg['soupx']['enabled'],
        'contamination_fraction': cfg['soupx']['contamination_fraction'],
        'seed': seed, 'n_empty': int(empty.sum()), 'sample': sample})
    kept.write_h5ad(sm.output.adata, compression='gzip')
    write_json(sm.output.summary, {'sample': sample, 'versions': versions(), 'kb_called': cells.n_obs,
        'qc_pass': int(cells.obs['qc_pass'].sum()), 'doublets': int(cells.obs['predicted_doublet'].sum()),
        'retained': kept.n_obs, 'empty_droplets': int(empty.sum()),
        'doublet_threshold': doublet_threshold, 'scrublet_enabled': cfg['scrublet']['enabled'],
        'mad_thresholds_log1p': thresholds, 'qc_reasons': reasons.sum().astype(int).to_dict(),
        'median_umi': float(np.median(kept.obs['total_counts'])),
        **{f'kb_recount_{key}': value for key, value in reconciliation.items()},
        'mt_genes_found': int(cells.var['mt'].sum()), 'hb_genes_found': int(cells.var['hb'].sum())})


with logged(snakemake):
    run(snakemake)
