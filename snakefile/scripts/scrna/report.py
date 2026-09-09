import sys
sys.path.insert(0, snakemake.scriptdir)
import base64
import html
import json
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from core import aggregate_symbols
from report_utils import design_notes, embedded_file, table
from runtime import logged


def run(sm):
    cfg = sm.params.cfg
    r = cfg['report']
    a = ad.read_h5ad(sm.input.adata)
    out = Path(sm.output.figures)
    out.mkdir(parents=True, exist_ok=True)
    requested = r['font_family']
    available = {f.name for f in font_manager.fontManager.ttflist}
    actual = requested if requested in available else ('Liberation Sans' if 'Liberation Sans' in available else 'DejaVu Sans')
    plt.rcParams.update({'font.family': actual, 'font.size': r['font_size'],
        'axes.labelsize': r['font_size'], 'axes.titlesize': r['font_size'] + 1,
        'xtick.labelsize': r['font_size'] - 1, 'ytick.labelsize': r['font_size'] - 1,
        'legend.fontsize': r['font_size'] - 1, 'axes.spines.top': False,
        'axes.spines.right': False, 'pdf.fonttype': 42, 'ps.fonttype': 42,
        'svg.fonttype': 'none', 'savefig.dpi': r['dpi'], 'axes.linewidth': 0.6})
    width = r['width_mm'] / 25.4
    sections = []
    figures = []
    palette = r['palette']
    sample_colors = {s: palette[i % len(palette)] for i, s in enumerate(sm.params.samples)}
    celltypes = sorted(a.obs['cell_type'].astype(str).unique())
    cell_colors = {s: r['celltype_colors'].get(s, palette[i % len(palette)]) for i, s in enumerate(celltypes)}
    cell_colors[cfg['scanvi']['unlabeled_category']] = '#999999'
    rng = np.random.default_rng(cfg['seed'])
    chosen = rng.choice(a.n_obs, min(a.n_obs, r['max_plot_cells']), replace=False)

    def add(title, explanation, body):
        anchor = f'section-{len(sections)}'
        sections.append((anchor, title, f'<section id="{anchor}"><h2>{html.escape(title)}</h2><p>{explanation}</p>{body}</section>'))

    def save(fig, name, caption):
        fig.tight_layout()
        for suffix in ['png', 'svg', 'pdf']:
            fig.savefig(out / f'{name}.{suffix}', dpi=r['dpi'], facecolor='white')
        plt.close(fig)
        png = base64.b64encode((out / f'{name}.png').read_bytes()).decode()
        figures.append(name)
        downloads = ' · '.join(embedded_file(out / f'{name}.{ext}', ext.upper(), mime)
            for ext, mime in [('png', 'image/png'), ('svg', 'image/svg+xml'), ('pdf', 'application/pdf')])
        return f'<figure><img alt="{html.escape(caption, quote=True)}" src="data:image/png;base64,{png}"><figcaption>{html.escape(caption)} · {downloads}</figcaption></figure>'

    def umap(ax, embedding, key, title):
        points = a.obsm[embedding][chosen]
        labels = a.obs[key].astype(str).to_numpy()[chosen]
        names = sorted(a.obs[key].astype(str).unique())
        colors = sample_colors if key == 'sample' else cell_colors if key == 'cell_type' else {s: palette[i % len(palette)] for i, s in enumerate(names)}
        for name in names:
            keep = labels == name
            ax.scatter(points[keep, 0], points[keep, 1], s=r['point_size'], color=colors[name], label=name, linewidths=0, rasterized=True)
        ax.set(title=title, xlabel='UMAP 1', ylabel='UMAP 2', xticks=[], yticks=[])
        ax.legend(loc='upper left', bbox_to_anchor=(1.01, 1), frameon=False, markerscale=3)

    qc = pd.DataFrame([json.loads(Path(p).read_text()) for p in sm.input.qc])
    summaries = [json.loads(Path(p).read_text()) for p in sm.input.summaries]
    notes = design_notes(a.obs, cfg['integration']['batch_key'])
    if actual != requested:
        notes.append(f'The requested font, {requested}, was unavailable; {actual} was used instead.')
    if len(celltypes) > len(palette):
        notes.append('The number of cell types exceeds the default palette size. Set celltype_colors or split cell types into separate panels for publication.')
    if (qc['mt_genes_found'] == 0).any():
        notes.append('No mitochondrial genes were identified in some samples. Check the gene symbols and mt_regex setting.')
    note_html = '<ul>' + ''.join(f'<li>{html.escape(n)}</li>' for n in notes) + '</ul>' if notes else ''
    add('Analysis overview', 'This primary analysis uses kb-derived UMI counts. Cell Ranger supplies sequencing QC metrics.',
        f'<div class="stats"><b>{a.n_obs:,}</b> retained cells / nuclei · <b>{a.n_vars:,}</b> genes · <b>{a.obs["sample"].nunique()}</b> libraries · <b>{a.obs["donor"].nunique()}</b> donors</div>'
        + note_html + '<p>kb → sample QC / Scrublet → SoupX → library depth normalization → scVI → CellTypist → scANVI → report</p>')
    meta_columns = [k for k in ['sample', 'donor', 'condition', 'batch', 'age', 'sex'] if k in a.obs]
    add('Samples and experimental design', 'Additional sequencing lanes from the same library are combined. Computational correction cannot disentangle batch effects from biological conditions when they are confounded.', table(a.obs[meta_columns].drop_duplicates()))

    fig, axes = plt.subplots(1, 2, figsize=(width, 3.0))
    positions = np.arange(len(qc))
    for column, label, color in [('retained', 'Retained', '#009E73'), ('doublets', 'Doublet', '#D55E00')]:
        bottom = qc['retained'].to_numpy() if column == 'doublets' else np.zeros(len(qc))
        axes[0].bar(positions, qc[column], bottom=bottom, label=label, color=color)
    axes[0].bar(positions, qc['kb_called'] - qc['qc_pass'], bottom=qc['qc_pass'], label='QC excluded', color='#AAAAAA')
    axes[0].set(xticks=positions, xticklabels=qc['sample'], ylabel='Called barcodes', title='Cell retention before SoupX')
    axes[0].tick_params(axis='x', rotation=90); axes[0].legend(frameon=False)
    for path, sample in zip(sm.input.ranks, sm.params.samples):
        rank = pd.read_csv(path, sep='\t')
        keep = np.unique(np.geomspace(1, len(rank), min(len(rank), 2000)).astype(int) - 1)
        axes[1].loglog(rank['rank'].to_numpy()[keep], np.maximum(1, rank['total_counts'].to_numpy()[keep]), label=sample, color=sample_colors[sample])
    axes[1].set(xlabel='Barcode rank', ylabel='kb UMI', title='Unfiltered droplet rank'); axes[1].legend(frameon=False)
    add('QC and cell retention', 'QC uses the kb filtered matrix because bustools corrects nearby barcode records to the called-cell allowlist before recounting UMIs. Small count differences from the corresponding rows of the unfiltered droplet matrix are expected and recorded below. Scrublet removes predicted doublets after detected-gene, UMI, mitochondrial, hemoglobin, and within-sample MAD filtering. Cells with zero total UMI after SoupX are excluded and recorded separately.',
        save(fig, 'qc_retention', 'QC retention and kb droplet rank') + table(qc[['sample', 'kb_called', 'qc_pass', 'doublets', 'retained', 'doublet_threshold', 'kb_recount_filtered_minus_unfiltered_umi', 'kb_recount_relative_absolute_difference']]))
    qc_cells = pd.concat([pd.read_csv(p, sep='\t') for p in sm.input.qc_cells], ignore_index=True)
    fig, axes = plt.subplots(1, 3, figsize=(width, 2.8))
    for ax, metric, label in zip(axes, ['total_counts', 'n_genes_by_counts', 'pct_counts_mt'], ['UMI', 'Detected genes', 'Mitochondrial UMI (%)']):
        for flag, color, name in [(False, '#AAAAAA', 'Excluded'), (True, '#0072B2', 'Retained')]:
            values = qc_cells.loc[qc_cells['retained'] == flag, metric].to_numpy()
            if len(values):
                ax.hist(np.log10(1 + values) if metric != 'pct_counts_mt' else values, bins=50, alpha=0.6, color=color, label=name)
        ax.set(xlabel=('log10(1 + ' + label + ')') if metric != 'pct_counts_mt' else label, ylabel='Cells')
    axes[0].legend(frameon=False)
    add('QC distributions', 'QC thresholds depend on the tissue, cell or nucleus preparation, and experimental conditions. Compare retained and excluded distributions to assess whether cell types with low RNA content are preferentially lost.', save(fig, 'qc_distributions', 'Pre-filter distributions of called cells'))
    fig, axes = plt.subplots(1, 2, figsize=(width, 2.7))
    for path, sample in zip(sm.input.qc_cells, sm.params.samples):
        scores = pd.read_csv(path, sep='\t')['doublet_score'].dropna()
        if len(scores): axes[0].hist(scores, bins=40, histtype='step', label=sample, color=sample_colors[sample])
    axes[0].set(xlabel='Scrublet score', ylabel='QC-passing cells')
    if axes[0].get_legend_handles_labels()[0]: axes[0].legend(frameon=False)
    soup = pd.DataFrame([json.loads(Path(p).read_text()) for p in sm.input.soup])
    axes[1].bar(soup['sample'], soup['removed_fraction'] * 100, color=[sample_colors[s] for s in soup['sample']])
    axes[1].set(ylabel='UMI removed (%)', title='SoupX ambient correction'); axes[1].tick_params(axis='x', rotation=90)
    add('Doublets and ambient RNA', 'Scrublet operates separately on each library. The SoupX ambient profile is estimated from low-UMI droplets in the unfiltered kb matrix, excluding all called barcodes. Failed automatic contamination estimates are never silently replaced with a fixed fraction. Disabled steps are identified in the accompanying tables.', save(fig, 'doublets_ambient', 'Doublet scores and ambient UMI removal') + table(soup))

    depth = pd.read_csv(sm.input.depth, sep='\t')
    fig, ax = plt.subplots(figsize=(width, 2.6))
    x = np.arange(len(depth))
    ax.bar(x - .18, depth['median_before'], .36, label='Full corrected depth', color='#0072B2')
    ax.bar(x + .18, depth['median_after'], .36, label='Matched depth' if cfg['depth']['enabled'] else 'Matching disabled', color='#E69F00')
    ax.set(xticks=x, xticklabels=depth['sample'], ylabel='Median UMI per cell'); ax.tick_params(axis='x', rotation=90); ax.legend(frameon=False)
    add('Library depth normalization', f'For visualization and CellTypist, each cell is normalized to a total of 10,000 counts and log1p-transformed. Optional depth matching uses binomial thinning based on sample median UMI counts. scVI/scANVI receives integer UMI counts from <code>{html.escape(cfg["depth"]["model_layer"])}</code>. Thinning discards information and cannot recover RNA that was undetected at low sequencing depth.', save(fig, 'library_depth', 'Sample median UMI before and after thinning') + table(depth))

    fig, axes = plt.subplots(2, 2, figsize=(width, 6))
    for row, (embedding, label) in enumerate([('X_umap_uncorrected', 'PCA before correction'), ('X_umap_scVI', 'scVI')]):
        umap(axes[row, 0], embedding, 'sample', label + ': sample')
        umap(axes[row, 1], embedding, 'condition', label + ': condition')
    add('scVI batch correction', 'UMAP is computed separately from the uncorrected PCA and scVI latent representations. Axes and distances are not directly comparable across independent embeddings. Improved batch mixing alone does not establish successful correction; also assess biological conditions and preservation of marker expression.', save(fig, 'integration', 'Uncorrected and scVI embeddings'))
    fig, axes = plt.subplots(1, 2, figsize=(width, 2.6))
    for metric, label in [('batch_entropy_before', 'Before'), ('batch_entropy_scVI', 'scVI')]:
        vals = a.obs[metric].dropna()
        if len(vals): axes[0].hist(vals, bins=30, histtype='step', label=label)
    axes[0].set(xlabel='Normalized neighbor batch entropy', ylabel='Cells')
    if axes[0].get_legend_handles_labels()[0]: axes[0].legend(frameon=False)
    for path, label in zip(sm.input.histories, ['scVI', 'scANVI']):
        history = pd.read_csv(path, sep='\t')
        selected = history[history['metric'] == 'elbo_train']
        if len(selected): axes[1].plot(selected['step'], selected['value'], label=label)
    axes[1].set(xlabel='Epoch', ylabel='Training ELBO loss')
    if axes[1].get_legend_handles_labels()[0]: axes[1].legend(frameon=False)
    add('Integration diagnostics', 'Neighbor batch entropy ranges from 0 for separation to 1 for equal mixing across batches. Its expected value depends on batch proportions and cell-type composition. Training curves diagnose convergence; they do not measure accuracy on independent validation data.', save(fig, 'integration_diagnostics', 'Neighbor entropy and training traces'))

    fig, axes = plt.subplots(1, 2, figsize=(width, 3.3))
    embedding = 'X_umap_scANVI' if cfg['scanvi']['enabled'] else 'X_umap_scVI'
    umap(axes[0], embedding, 'cell_type', 'Final cell type')
    if 'scanvi_normalized_entropy' in a.obs:
        scatter = axes[1].scatter(a.obsm[embedding][chosen, 0], a.obsm[embedding][chosen, 1],
            c=a.obs['scanvi_normalized_entropy'].to_numpy()[chosen], s=r['point_size'],
            cmap='cividis', vmin=0, vmax=1, linewidths=0, rasterized=True)
        fig.colorbar(scatter, ax=axes[1], label='Normalized label entropy')
        axes[1].set(title='scANVI uncertainty', xlabel='UMAP 1', ylabel='UMAP 2', xticks=[], yticks=[])
    else:
        axes[1].text(.5, .5, 'scANVI disabled', ha='center'); axes[1].axis('off')
    add('Cell types and uncertainty', 'When enabled, scANVI is trained using sufficiently supported, high-confidence CellTypist labels as seeds. Repeated latent-variable draws yield class probabilities, standard deviations, and entropy. Cells that fail the confidence thresholds receive the configured unknown label. These probabilities do not necessarily detect reference-model errors or cell types absent from training.', save(fig, 'annotation', 'Final annotation and label uncertainty'))
    if cfg['scanvi']['enabled']:
        fig, axes = plt.subplots(1, 2, figsize=(width, 2.7))
        axes[0].hist(a.obs['scanvi_max_probability'], bins=40, color='#0072B2')
        axes[0].axvline(cfg['scanvi']['min_probability'], color='#D55E00', linestyle='--')
        axes[0].set(xlabel='Maximum class probability', ylabel='Cells')
        cross = pd.crosstab(a.obs['celltypist_label'], a.obs['scanvi_label'], normalize='index')
        axes[1].imshow(cross, aspect='auto', cmap='cividis', vmin=0, vmax=1)
        axes[1].set(xticks=np.arange(len(cross.columns)), xticklabels=cross.columns,
                    yticks=np.arange(len(cross.index)), yticklabels=cross.index,
                    xlabel='scANVI', ylabel='CellTypist')
        axes[1].tick_params(axis='x', rotation=90)
        add('Annotation agreement', 'Agreement is influenced by training on CellTypist-derived seeds and is not an independent accuracy assessment. CellTypist one-vs-rest scores do not sum to one across classes and have a different interpretation from scANVI class probabilities.', save(fig, 'annotation_diagnostics', 'scANVI confidence and seed-model label agreement'))

    requested_markers = list(dict.fromkeys(g for genes in r['markers'].values() for g in genes))
    x, symbols = aggregate_symbols(a.layers[cfg['depth']['model_layer']], a.var['gene_symbol'])
    markers = [g for g in requested_markers if g in symbols]
    if markers:
        totals = np.asarray(x.sum(axis=1)).ravel()
        subset = x[:, symbols.get_indexer(markers)].multiply((1e4 / np.maximum(1, totals))[:, None]).tocsr()
        subset.data = np.log1p(subset.data)
        means = np.vstack([np.asarray(subset[a.obs['cell_type'].astype(str).to_numpy() == group].mean(axis=0)).ravel() for group in celltypes])
        fig, ax = plt.subplots(figsize=(width, max(2.5, len(celltypes) * .18)))
        plot = ax.imshow(means, cmap='cividis', aspect='auto')
        ax.set(xticks=np.arange(len(markers)), xticklabels=markers, yticks=np.arange(len(celltypes)), yticklabels=celltypes)
        ax.tick_params(axis='x', rotation=90); fig.colorbar(plot, ax=ax, label='Mean log1p(CP10k)')
        add('Marker expression', 'Configured markers are displayed by cell type using normalized expression from the selected count layer. These values are not scVI-derived corrected expression estimates.', save(fig, 'markers', 'Marker expression by final cell type'))

    composition = pd.read_csv(sm.input.composition, sep='\t')
    fig, ax = plt.subplots(figsize=(width, 3))
    samples = list(sm.params.samples)
    ctypes = sorted(composition['cell_type'].unique())
    positions = np.arange(len(samples))
    for i, label in enumerate(ctypes):
        sub = composition[composition['cell_type'] == label].set_index('sample').reindex(samples)
        offset = (i - (len(ctypes) - 1) / 2) * .7 / max(1, len(ctypes))
        means = sub['expected_fraction'].to_numpy()
        color = cell_colors.get(label, palette[i % len(palette)])
        ax.scatter(positions + offset, means, color=color, label=label, s=12)
        if sub['label_interval_low'].notna().any():
            # Quantiles need not straddle the expectation for finite discrete draws.
            ax.vlines(positions + offset, sub['label_interval_low'], sub['label_interval_high'], colors=color, linewidth=.8)
    ax.set(xticks=positions, xticklabels=samples, ylabel='Cell fraction', ylim=(0, 1)); ax.tick_params(axis='x', rotation=90)
    ax.legend(loc='upper left', bbox_to_anchor=(1.01, 1), frameon=False)
    add('Propagation of uncertainty to cell composition', 'Expected cell counts and fractions are calculated by summing scANVI probabilities. Repeated probabilistic label assignments provide 95% intervals. All probability mass is included, including cells displayed as unknown. These intervals describe label uncertainty conditional on the observed cells and fitted model; they are not biological confidence intervals and exclude donor variation, sampling variation, and SoupX estimation error. When scANVI is disabled, only hard counts are reported.', save(fig, 'composition', 'Expected composition and conditional label-imputation intervals') + table(composition) + embedded_file(sm.input.composition, 'Download composition TSV', 'text/tab-separated-values'))

    metrics_html = ''
    for sample, path in zip(sm.params.samples, sm.input.cr):
        frame = pd.read_csv(path, dtype=str, keep_default_na=False)
        if len(frame) != 1:
            raise ValueError(f'{path}: expected one Cell Ranger metrics summary row')
        metrics_html += f'<details><summary>{html.escape(sample)}</summary>' + table(frame.T.reset_index().set_axis(['Metric', 'Value'], axis=1)) + '</details>'
    add('Mapping QC', 'Cell Ranger metrics describe the same FASTQ data processed with a separate quantification method. Cell and UMI definitions differ from kb and should not be substituted directly.', metrics_html or '<p>Cell Ranger metrics collection is disabled.</p>')
    for sample, path in zip(sm.params.samples, sm.input.kb):
        payload = json.loads(Path(path).read_text())
        sections[-1] = (sections[-1][0], sections[-1][1], sections[-1][2].replace('</section>', f'<details><summary>kb inspect: {html.escape(sample)}</summary><pre>{html.escape(json.dumps(payload, indent=2))}</pre></details></section>'))
    provenance = json.loads(Path(sm.input.provenance).read_text())
    provenance['stage_summaries'] = summaries
    provenance['mapping_versions'] = {str(p): Path(p).read_text() for p in sm.input.mapping_versions}
    provenance['figure_style'] = {'requested_font': requested, 'actual_font': actual,
        'font_size_pt': r['font_size'], 'width_mm': r['width_mm'], 'dpi': r['dpi']}
    add('Reproducibility and outputs', 'Random seeds, configuration, source hashes, model hashes, and software versions are recorded. Figures use white backgrounds and a default palette designed for color-vision accessibility. Editable SVG, font-embedded PDF, and PNG files are provided. Check readability at the final publication size.',
        '<p>final.h5ad contains counts_raw / counts / counts_depthmatched, latent representations, UMAP coordinates, and class probabilities. cells.tsv.gz contains per-cell QC metrics, labels, and uncertainty measures.</p>'
        + f'<details><summary>Configuration and provenance</summary><pre>{html.escape(json.dumps(provenance, indent=2, ensure_ascii=False))}</pre></details>')
    refs = [('kb-python', 'https://github.com/pachterlab/kb_python'),
        ('SoupX', 'https://github.com/constantAmateur/SoupX'),
        ('Scrublet', 'https://github.com/swolock/scrublet'),
        ('scVI/scANVI', 'https://docs.scvi-tools.org/en/1.3.3/'),
        ('CellTypist', 'https://celltypist.readthedocs.io/en/latest/'),
        ('Cell Press figure guidelines', 'https://www.cell.com/figureguidelines')]
    add('Method references', 'Consult the current tool documentation and journal guidelines when selecting analysis settings and preparing final figures.', '<ul>' + ''.join(f'<li><a href="{url}">{name}</a></li>' for name, url in refs) + '</ul>')
    nav = ' · '.join(f'<a href="#{anchor}">{html.escape(title)}</a>' for anchor, title, _ in sections)
    document = '<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">'
    document += f'<title>{html.escape(r["title"])}</title>'
    document += '''<style>body{font:16px/1.65 Arial,"Noto Sans",sans-serif;color:#20252b;background:#fff;max-width:1180px;margin:2rem auto;padding:0 1rem}h1{font-size:2rem}h2{font-size:1.35rem;border-bottom:2px solid #0072b2;padding-bottom:.3rem}section{margin:3rem 0}nav{padding:1rem;background:#f2f6f8}a{color:#00679f}figure{margin:1.5rem 0}img{max-width:100%;height:auto}figcaption{color:#52606a;font-size:.85rem}table{border-collapse:collapse;font-size:.85rem}td,th{padding:.4rem .65rem;border-bottom:1px solid #dce2e6;text-align:left}th{background:#f2f6f8}.table-wrap{overflow:auto}pre{overflow:auto;background:#f6f8fa;padding:1rem;font-size:.8rem}details{margin:.8rem 0}.stats{font-size:1.15rem;background:#eef6fa;padding:1rem}@media print{nav{display:none}section{break-inside:avoid}}</style>'''
    document += f'<body><h1>{html.escape(r["title"])}</h1><p>Self-contained analysis report · kb counts</p><nav>{nav}</nav>'
    document += ''.join(body for _, _, body in sections) + '</body></html>'
    Path(sm.output.html).write_text(document, encoding='utf-8')


with logged(snakemake):
    run(snakemake)
