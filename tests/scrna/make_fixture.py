"""Generate small synthetic kb H5ADs and an optional locally trained toy model.

Never use the toy model to annotate biological data.
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import sparse
import anndata as ad
import yaml

ROOT = Path(__file__).resolve().parents[2]


def make_fixture(out, train_model=False):
    out = Path(out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    cfg = yaml.safe_load((ROOT / 'examples/scrna_seq/config.yaml').read_text())
    cfg['workdir'] = str(out / 'results')
    cfg['sample_table'] = str(out / 'samples.tsv')
    cfg['input_mode'] = 'counts'
    cfg['reference']['build'] = False
    cfg['cellranger']['enabled'] = True
    cfg['kb']['count_layer'] = 'auto'
    cfg['qc'].update(min_genes=10, min_counts=20, max_pct_mt=100, max_pct_hb=100, mad_n=20, min_cells_per_sample=30)
    cfg['scrublet'].update(n_prin_comps=5, threshold=.9, min_counts=1, min_cells=3, min_gene_variability_pctl=0)
    cfg['soupx'].update(min_empty_droplets=10, contamination_fraction=.08)
    cfg['integration'].update(n_hvg=100, n_latent=5, n_layers=1, max_epochs=2, n_neighbors=10, batch_size=32)
    cfg['celltypist'].update(model='synthetic_test_only.pkl', model_path=str(out / 'toy.pkl'),
                            min_model_genes=10, seed_probability=.5, seed_margin=.05, majority_voting=False)
    cfg['scanvi'].update(min_seed_cells=5, max_epochs=2, n_samples_per_label=10,
                         posterior_samples=3, composition_draws=50)
    cfg['resources'].update(analysis_threads=1, analysis_mem_mb=4000)
    cfg['report'].update(title='Synthetic smoke test <not biological data>', dpi=100,
                         markers={'A': ['Gene0', 'Gene1'], 'B': ['Gene20', 'Gene21']})
    rng = np.random.default_rng(42)
    genes = pd.Index([f'Gene{i}' for i in range(600)])
    records, training = [], []
    for sample, depth in [('a', 1.0), ('b', .5)]:
        labels = np.array(['TypeA'] * 60 + ['TypeB'] * 60)
        rates = np.tile(rng.gamma(1.5, 1.0, size=600), (120, 1))
        rates[:60, :20] = 5
        rates[60:, 20:40] = 5
        counts = rng.poisson(rates * depth * rng.lognormal(0, .4, size=rates.shape))
        # kb H5ADs name their observation index "barcode".
        barcodes = pd.Index([f'BC{i:04d}' for i in range(120)], name='barcode')
        cells = ad.AnnData(sparse.csr_matrix(counts), obs=pd.DataFrame({'label': labels}, index=barcodes),
                           var=pd.DataFrame({'gene_name': genes}, index=genes))
        # kb NAC X is mature only even when --sum total writes a separate MTX.
        cells.layers['mature'] = sparse.csr_matrix(counts // 2)
        cells.layers['nascent'] = sparse.csr_matrix(counts // 3)
        cells.layers['ambiguous'] = sparse.csr_matrix(counts - counts // 2 - counts // 3)
        cells.X = cells.layers['mature'].copy()
        cells.write_h5ad(out / f'{sample}.filtered.h5ad')
        empty = rng.poisson(.04, size=(100, 600))
        droplets = ad.AnnData(sparse.csr_matrix(np.vstack([counts, empty])),
            obs=pd.DataFrame(index=list(barcodes) + [f'EMPTY{i}' for i in range(100)]), var=cells.var.copy())
        for name, denominator in [('mature', 2), ('nascent', 3)]:
            droplets.layers[name] = sparse.csr_matrix(np.vstack([counts // denominator, empty // denominator]))
        droplets.layers['ambiguous'] = droplets.X - droplets.layers['mature'] - droplets.layers['nascent']
        # Real kb output can differ at called barcodes because bustools corrects
        # nearby barcode records to the filtered allowlist before recounting.
        row, column = np.argwhere(droplets.layers['mature'][:120].toarray() > 0)[0]
        mature = droplets.layers['mature'].tolil()
        mature[row, column] -= 1
        droplets.layers['mature'] = mature.tocsr()
        droplets.X = droplets.layers['mature'].copy()
        droplets.write_h5ad(out / f'{sample}.unfiltered.h5ad')
        pd.DataFrame([{'Estimated Number of Cells': '120', 'Mean Reads per Cell': '10,000',
                       'Valid Barcodes': '95.0%'}]).to_csv(out / f'{sample}.metrics.csv', index=False)
        records.append({'sample': sample, 'batch': sample, 'donor': sample, 'condition': 'control',
            'kb_filtered': str(out / f'{sample}.filtered.h5ad'),
            'kb_unfiltered': str(out / f'{sample}.unfiltered.h5ad'), 'cellranger_metrics': str(out / f'{sample}.metrics.csv')})
        cells.X = sparse.csr_matrix(counts)
        cells.obs_names = sample + ':' + cells.obs_names
        training.append(cells)
    pd.DataFrame(records).to_csv(out / 'samples.tsv', sep='\t', index=False)
    (out / 'config.yaml').write_text(yaml.safe_dump(cfg, sort_keys=False))
    if train_model:
        import scanpy as sc
        import celltypist
        data = ad.concat(training)
        sc.pp.normalize_total(data, target_sum=1e4)
        sc.pp.log1p(data)
        model = celltypist.train(data, labels='label', check_expression=True, n_jobs=1, max_iter=200)
        model.write(str(out / 'toy.pkl'))
    else:
        (out / 'toy.pkl').touch()
    return out / 'config.yaml'


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('out')
    parser.add_argument('--train-model', action='store_true')
    args = parser.parse_args()
    print(make_fixture(args.out, args.train_model))
