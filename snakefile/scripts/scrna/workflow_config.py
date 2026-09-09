"""Validate metadata before constructing the DAG (stdlib + PyYAML only)."""
from pathlib import Path
import copy
import csv
import math
import re


def resolve_path(value, base):
    p = Path(value).expanduser()
    return str((base / p).resolve()) if not p.is_absolute() else str(p.resolve())


def prepare_config(config, defaults, base):
    cfg = copy.deepcopy(defaults)

    def merge(dst, src, prefix=""):
        for key, value in src.items():
            if key not in dst:
                raise ValueError(f"Unknown setting: {prefix}{key}")
            if isinstance(dst[key], dict) and key not in ("markers", "celltype_colors"):
                if not isinstance(value, dict):
                    raise ValueError(f"{prefix}{key} must be a mapping")
                merge(dst[key], value, prefix + key + ".")
            else:
                dst[key] = value
    merge(cfg, config)
    if cfg['input_mode'] not in ('fastq', 'counts'):
        raise ValueError('input_mode must be fastq or counts')
    if cfg['kb']['workflow'] not in ('standard', 'nac'):
        raise ValueError('kb.workflow must be standard or nac')
    if cfg['depth']['model_layer'] not in ('counts', 'counts_depthmatched'):
        raise ValueError('depth.model_layer must be counts or counts_depthmatched')
    if not cfg['depth']['enabled'] and cfg['depth']['model_layer'] != 'counts':
        raise ValueError('counts_depthmatched requires depth.enabled=true')
    for section, key in [('reference', 'build'), ('cellranger', 'enabled'),
                         ('cellranger', 'include_introns'), ('cellranger', 'create_bam'),
                         ('scrublet', 'enabled'), ('soupx', 'enabled'),
                         ('depth', 'enabled'), ('scanvi', 'enabled'),
                         ('celltypist', 'majority_voting'), ('qc', 'mad_upper')]:
        if not isinstance(cfg[section][key], bool):
            raise ValueError(f'{section}.{key} must be a YAML boolean')
    numeric = {
        'qc': ['min_genes', 'min_counts', 'mad_n', 'min_cells_per_sample'],
        'scrublet': ['n_prin_comps', 'sim_doublet_ratio', 'min_counts', 'min_cells'],
        'soupx': ['min_empty_droplets', 'clusters_resolution'],
        'integration': ['n_hvg', 'n_latent', 'n_layers', 'max_epochs', 'batch_size',
                        'devices', 'n_neighbors', 'leiden_resolution'],
        'celltypist': ['min_model_genes', 'over_clustering_resolution'],
        'scanvi': ['min_seed_cells', 'max_epochs', 'n_samples_per_label',
                   'posterior_samples', 'composition_draws'],
        'report': ['font_size', 'dpi', 'width_mm', 'point_size', 'max_plot_cells'],
        'resources': ['mapping_threads', 'mapping_mem_mb', 'analysis_threads', 'analysis_mem_mb', 'model_mem_mb'],
    }
    for section, keys in numeric.items():
        for key in keys:
            value = cfg[section][key]
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
                raise ValueError(f'{section}.{key} must be positive')
    for section, keys in {'scrublet': ['expected_doublet_rate', 'threshold'],
                          'soupx': ['contamination_fraction'],
                          'celltypist': ['min_prop', 'seed_probability', 'seed_margin'],
                          'scanvi': ['min_probability', 'max_normalized_entropy']}.items():
        for key in keys:
            v = cfg[section][key]
            if v is not None and (isinstance(v, bool) or not isinstance(v, (int, float)) or not 0 <= v <= 1):
                raise ValueError(f'{section}.{key} must be in [0, 1]')
    if not 0 < cfg['scrublet']['expected_doublet_rate'] < 1:
        raise ValueError('scrublet.expected_doublet_rate must be in (0, 1)')
    if cfg['scanvi']['posterior_samples'] < 2:
        raise ValueError('scanvi.posterior_samples must be >=2')
    integer_settings = {
        'qc': ['min_genes', 'min_counts', 'min_cells_per_sample'],
        'scrublet': ['n_prin_comps', 'min_counts', 'min_cells'], 'soupx': ['min_empty_droplets'],
        'integration': ['n_hvg', 'n_latent', 'n_layers', 'max_epochs', 'batch_size', 'devices', 'n_neighbors'],
        'scanvi': ['min_seed_cells', 'max_epochs', 'n_samples_per_label', 'posterior_samples', 'composition_draws'],
        'celltypist': ['min_model_genes'], 'report': ['max_plot_cells'],
        'resources': ['mapping_threads', 'mapping_mem_mb', 'analysis_threads', 'analysis_mem_mb', 'model_mem_mb', 'gpu'],
    }
    for section, keys in integer_settings.items():
        for key in keys:
            if type(cfg[section][key]) is not int:
                raise ValueError(f'{section}.{key} must be an integer')
    if type(cfg['seed']) is not int or not 0 <= cfg['seed'] < 2**32:
        raise ValueError('seed must be an integer in [0, 2**32)')
    if not cfg['report']['palette'] or any(not re.fullmatch(r'#[0-9a-fA-F]{6}', c) for c in cfg['report']['palette']):
        raise ValueError('report.palette must contain #RRGGBB colors')
    for key in ['max_pct_mt', 'max_pct_hb']:
        if not isinstance(cfg['qc'][key], (int, float)) or not 0 <= cfg['qc'][key] <= 100:
            raise ValueError(f'qc.{key} must be in [0,100]')
    if not 0 <= cfg['scrublet']['min_gene_variability_pctl'] <= 100:
        raise ValueError('scrublet.min_gene_variability_pctl must be in [0,100]')
    for floor, ceiling in [('min_genes', 'max_genes'), ('min_counts', 'max_counts')]:
        if cfg['qc'][ceiling] is not None and cfg['qc'][ceiling] < cfg['qc'][floor]:
            raise ValueError(f'qc.{ceiling} must be >= qc.{floor}')
    if cfg['depth']['target_median'] is not None and cfg['depth']['target_median'] <= 0:
        raise ValueError('depth.target_median must be positive or null')
    lo, hi = cfg['soupx']['soup_range']
    if not 0 <= lo < hi:
        raise ValueError('soupx.soup_range must be [lower, upper], lower < upper')
    for key in ['mt_regex', 'hb_regex']:
        re.compile(cfg['qc'][key])
    if cfg['integration']['accelerator'] not in ('cpu', 'gpu'):
        raise ValueError('integration.accelerator must be cpu or gpu')
    if cfg['integration']['accelerator'] == 'gpu' and cfg['resources']['gpu'] < 1:
        raise ValueError('GPU training requires resources.gpu >= 1')
    if cfg['integration']['devices'] != 1:
        raise ValueError('This workflow supports one CPU/GPU training device per job')
    for key in ['workdir', 'sample_table']:
        cfg[key] = resolve_path(cfg[key], base)
    for section, keys in {'reference': ['fasta', 'gtf', 'index', 't2g', 'cdna_t2c', 'intron_t2c'],
                          'cellranger': ['transcriptome'], 'kb': ['whitelist'],
                          'celltypist': ['model_path']}.items():
        for key in keys:
            if cfg[section][key]:
                cfg[section][key] = resolve_path(cfg[section][key], base)
    for key, value in cfg['containers'].items():
        if not value or not isinstance(value, str):
            raise ValueError(f'containers.{key} must specify an image')
        if '://' not in value:
            cfg['containers'][key] = resolve_path(value, base)
    if cfg['input_mode'] == 'fastq':
        required = ['fasta', 'gtf'] if cfg['reference']['build'] else ['index', 't2g']
        if not cfg['reference']['build'] and cfg['kb']['workflow'] == 'nac':
            required += ['cdna_t2c', 'intron_t2c']
        for key in required:
            if not cfg['reference'][key]:
                raise ValueError(f'reference.{key} is required')
        if cfg['cellranger']['enabled'] and not cfg['cellranger']['transcriptome']:
            raise ValueError('cellranger.transcriptome is required')
    with open(cfg['sample_table'], newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        required = {'sample', 'batch', 'donor', 'condition'}
        required |= {'R1', 'R2'} if cfg['input_mode'] == 'fastq' else {'kb_filtered', 'kb_unfiltered'}
        if cfg['input_mode'] == 'counts' and cfg['cellranger']['enabled']:
            required.add('cellranger_metrics')
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f'Sample table missing columns: {sorted(missing)}')
        reserved = {'barcode', 'cell_type', 'total_counts', 'n_genes_by_counts', 'pct_counts_mt',
                    'pct_counts_hb', 'qc_pass', 'qc_reason', 'retained', 'predicted_doublet', 'doublet_score'}
        if reserved.intersection(reader.fieldnames):
            raise ValueError(f'Sample metadata uses reserved output columns: {sorted(reserved.intersection(reader.fieldnames))}')
        samples = {}
        used_fastqs = set()
        for row in reader:
            if None in row:
                raise ValueError('Sample table has more fields than headers')
            row = {k: (v or '').strip() for k, v in row.items()}
            sample = row['sample']
            if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]*', sample):
                raise ValueError(f'Invalid sample ID: {sample!r}')
            if sample in samples:
                raise ValueError(f'Duplicate library {sample}; combine top-up lanes in one row')
            if any(not row[k] for k in required):
                raise ValueError(f'{sample}: missing required metadata/path')
            for key in ['R1', 'R2']:
                if cfg['input_mode'] == 'fastq':
                    row[key] = [resolve_path(p.strip(), base) for p in row[key].split(',') if p.strip()]
            if cfg['input_mode'] == 'fastq':
                if not row['R1'] or len(row['R1']) != len(row['R2']) or len(row['R1']) > 999:
                    raise ValueError(f'{sample}: R1/R2 must have 1–999 matching lanes')
                for p in row['R1'] + row['R2']:
                    if p in used_fastqs:
                        raise ValueError(f'FASTQ reused across reads/libraries: {p}')
                    if not p.endswith('.fastq.gz') and not p.endswith('.fq.gz'):
                        raise ValueError(f'Gzipped FASTQ required: {p}')
                    used_fastqs.add(p)
            for key in ['kb_filtered', 'kb_unfiltered', 'cellranger_metrics']:
                if row.get(key):
                    row[key] = resolve_path(row[key], base)
            samples[sample] = row
    if not samples:
        raise ValueError('Sample table is empty')
    metadata_keys = [cfg['integration']['batch_key']] + cfg['depth']['groupby']
    metadata_keys += cfg['integration']['categorical_covariate_keys'] + cfg['integration']['continuous_covariate_keys']
    for key in metadata_keys:
        if any(not row.get(key) for row in samples.values()):
            raise ValueError(f'Missing sample metadata for {key}')
    for key in cfg['integration']['continuous_covariate_keys']:
        if any(not math.isfinite(float(row[key])) for row in samples.values()):
            raise ValueError(f'{key} contains non-finite covariates')
    return cfg, samples
