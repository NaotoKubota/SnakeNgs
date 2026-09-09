"""kb → sample QC/Scrublet → SoupX → depth → scVI → CellTypist → scANVI → HTML.

snakemake -s snakefile/scrna_seq.smk --configfile examples/scrna_seq/config.yaml \
    --cores 16 --use-singularity --rerun-incomplete
"""
import sys
import re
import hashlib
from pathlib import Path
import yaml

SOURCE = Path(workflow.basedir).resolve()
sys.path.insert(0, str(SOURCE / 'scripts' / 'scrna'))
from workflow_config import prepare_config
with open(SOURCE.parent / 'examples' / 'scrna_seq' / 'config.yaml') as handle:
    DEFAULTS = yaml.safe_load(handle)
CFG, SAMPLES = prepare_config(config, DEFAULTS, Path.cwd())
# Imported helpers are not automatically included in Snakemake's script checksum.
CFG['_source_sha256'] = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                        for p in sorted((SOURCE / 'scripts' / 'scrna').glob('*')) if p.is_file()}
workdir: CFG['workdir']
C = CFG['containers']
NAC = CFG['kb']['workflow'] == 'nac'
FASTQ = CFG['input_mode'] == 'fastq'
BUILD = FASTQ and CFG['reference']['build']
REF = {k: f'reference/{name}' if BUILD else CFG['reference'][k]
       for k, name in [('index', 'index.idx'), ('t2g', 't2g.txt'),
                       ('cdna_t2c', 'cdna_t2c.txt'), ('intron_t2c', 'intron_t2c.txt')]}

def kb_file(wc, which):
    return f'kb/{wc.sample}/counts_{which}/adata.h5ad' if FASTQ else SAMPLES[wc.sample]['kb_' + which]

def cr_metrics(sample):
    return f'cellranger/{sample}/outs/metrics_summary.csv' if FASTQ else SAMPLES[sample]['cellranger_metrics']

wildcard_constraints:
    sample='|'.join(re.escape(s) for s in SAMPLES)

rule all:
    input:
        'report/report.html',
        'analysis/final.h5ad',
        'analysis/composition.tsv',
        'analysis/cells.tsv.gz'
    container: C['analysis']

rule provenance:
    input:
        samples=CFG['sample_table'],
        sources=([str(p) for p in sorted((SOURCE / 'scripts' / 'scrna').glob('*')) if p.is_file()]
                 + [str(SOURCE.parent / 'containers' / name) for name in
                    ['scrna-analysis.def', 'scrna-soupx.def', 'scrna-requirements.txt']]
                 + [str(SOURCE.parent / 'examples' / 'scrna_seq' / 'config.yaml')]),
        snakefile=str(SOURCE / 'scrna_seq.smk')
    output: 'provenance/run.json'
    params: cfg=CFG, samples=SAMPLES
    container: C['analysis']
    log: 'logs/provenance.log'
    script: 'scripts/scrna/provenance.py'

if BUILD:
    rule kb_reference:
        input: fasta=CFG['reference']['fasta'], gtf=CFG['reference']['gtf']
        output:
            index=REF['index'], t2g=REF['t2g'],
            captures=[REF['cdna_t2c'], REF['intron_t2c']] if NAC else [],
            versions='provenance/kb_reference.json'
        params:
            workflow=CFG['kb']['workflow'],
            captures=['-f2', 'reference/intron.fa', '-c1', REF['cdna_t2c'], '-c2', REF['intron_t2c']] if NAC else []
        threads: CFG['resources']['mapping_threads']
        resources: mem_mb=CFG['resources']['mapping_mem_mb']
        container: C['kb']
        log: 'logs/kb_reference.log'
        benchmark: 'benchmarks/kb_reference.tsv'
        shell:
            r"""
            kb ref -i {output.index:q} -g {output.t2g:q} -f1 reference/cdna.fa \
                -t {threads} --workflow {params.workflow:q} {params.captures:q} \
                {input.fasta:q} {input.gtf:q} > {log:q} 2>&1
            python -c 'import json, platform, importlib.metadata; print(json.dumps(dict(python=platform.python_version(), kb_python=importlib.metadata.version("kb-python"))))' > {output.versions:q}
            """

if FASTQ:
    rule kb_count:
        input:
            fastqs=lambda wc: [p for pair in zip(SAMPLES[wc.sample]['R1'], SAMPLES[wc.sample]['R2']) for p in pair],
            index=REF['index'], t2g=REF['t2g'],
            captures=[REF['cdna_t2c'], REF['intron_t2c']] if NAC else [],
            whitelist=[CFG['kb']['whitelist']] if CFG['kb']['whitelist'] else []
        output:
            filtered='kb/{sample}/counts_filtered/adata.h5ad',
            unfiltered='kb/{sample}/counts_unfiltered/adata.h5ad',
            inspect='kb/{sample}/inspect.json',
            versions='provenance/kb_{sample}.json'
        params:
            out='kb/{sample}',
            technology=CFG['kb']['technology'],
            workflow=CFG['kb']['workflow'],
            memory=f"{max(1, int(CFG['resources']['mapping_mem_mb'] * 0.8))}M",
            captures=['-c1', REF['cdna_t2c'], '-c2', REF['intron_t2c'], '--sum', 'total'] if NAC else [],
            whitelist=['-w', CFG['kb']['whitelist']] if CFG['kb']['whitelist'] else []
        threads: CFG['resources']['mapping_threads']
        resources: mem_mb=CFG['resources']['mapping_mem_mb']
        container: C['kb']
        log: 'logs/kb_{sample}.log'
        benchmark: 'benchmarks/kb_{sample}.tsv'
        shell:
            r"""
            kb count -i {input.index:q} -g {input.t2g:q} -o {params.out:q} \
                -x {params.technology:q} -t {threads} -m {params.memory:q} \
                --workflow {params.workflow:q} --h5ad --filter bustools --overwrite \
                {params.captures:q} {params.whitelist:q} {input.fastqs:q} > {log:q} 2>&1
            python -c 'import json, platform, importlib.metadata; print(json.dumps(dict(python=platform.python_version(), kb_python=importlib.metadata.version("kb-python"))))' > {output.versions:q}
            """

    if CFG['cellranger']['enabled']:
        rule stage_cellranger_fastqs:
            input:
                r1=lambda wc: SAMPLES[wc.sample]['R1'],
                r2=lambda wc: SAMPLES[wc.sample]['R2']
            output: directory('fastqs/{sample}')
            container: C['analysis']
            log: 'logs/stage_{sample}.log'
            script: 'scripts/scrna/stage_fastqs.py'

        rule cellranger_metrics:
            input:
                fastqs='fastqs/{sample}',
                # Track original FASTQs as well as staged directory (top-up changes).
                reads=lambda wc: SAMPLES[wc.sample]['R1'] + SAMPLES[wc.sample]['R2'],
                transcriptome=CFG['cellranger']['transcriptome']
            output:
                metrics='cellranger/{sample}/outs/metrics_summary.csv',
                web='cellranger/{sample}/outs/web_summary.html',
                versions='provenance/cellranger_{sample}.txt'
            params:
                introns=str(CFG['cellranger']['include_introns']).lower(),
                bam=str(CFG['cellranger']['create_bam']).lower(),
                chemistry=CFG['cellranger']['chemistry'],
                expect=['--expect-cells', str(CFG['cellranger']['expect_cells'])] if CFG['cellranger']['expect_cells'] else []
            threads: CFG['resources']['mapping_threads']
            resources:
                mem_mb=CFG['resources']['mapping_mem_mb'],
                cr_mem_gb=max(1, CFG['resources']['mapping_mem_mb'] // 1024)
            container: C['cellranger']
            log: 'logs/cellranger_{sample}.log'
            benchmark: 'benchmarks/cellranger_{sample}.tsv'
            shell:
                r"""
                mkdir -p cellranger provenance
                cellranger --version > {output.versions:q}
                # Snakemake creates declared output parents before execution.
                # Give Cell Ranger a separate, resumable pipestance directory.
                cellranger count --id {wildcards.sample:q} --output-dir cellranger/{wildcards.sample}/run \
                    --transcriptome {input.transcriptome:q} --fastqs "$PWD/{input.fastqs}" \
                    --sample {wildcards.sample:q} --include-introns {params.introns:q} \
                    --create-bam {params.bam:q} --chemistry {params.chemistry:q} \
                    --localcores {threads} --localmem {resources.cr_mem_gb} \
                    {params.expect:q} > {log:q} 2>&1
                cp cellranger/{wildcards.sample}/run/outs/metrics_summary.csv {output.metrics:q}
                cp cellranger/{wildcards.sample}/run/outs/web_summary.html {output.web:q}
                """

rule sample_qc:
    input:
        filtered=lambda wc: kb_file(wc, 'filtered'),
        unfiltered=lambda wc: kb_file(wc, 'unfiltered'),
        t2g=[REF['t2g']] if FASTQ else ([CFG['reference']['t2g']] if CFG['reference']['t2g'] else [])
    output:
        adata='qc/{sample}/singlets.h5ad',
        cells='qc/{sample}/cells.tsv.gz',
        rank='qc/{sample}/barcode_rank.tsv.gz',
        soupx=directory('qc/{sample}/soupx_input'),
        summary='qc/{sample}/summary.json'
    params: cfg=CFG, metadata=lambda wc: SAMPLES[wc.sample]
    threads: CFG['resources']['analysis_threads']
    resources: mem_mb=CFG['resources']['analysis_mem_mb']
    container: C['analysis']
    log: 'logs/qc_{sample}.log'
    benchmark: 'benchmarks/qc_{sample}.tsv'
    script: 'scripts/scrna/sample_qc.py'

rule soupx:
    input: 'qc/{sample}/soupx_input'
    output:
        matrix='soupx/{sample}/corrected.mtx',
        genes='soupx/{sample}/genes.tsv',
        barcodes='soupx/{sample}/barcodes.tsv',
        rho='soupx/{sample}/rho.tsv',
        summary='soupx/{sample}/summary.json'
    params: cfg=CFG
    threads: CFG['resources']['analysis_threads']
    resources: mem_mb=CFG['resources']['analysis_mem_mb']
    container: C['soupx']
    log: 'logs/soupx_{sample}.log'
    benchmark: 'benchmarks/soupx_{sample}.tsv'
    script: 'scripts/scrna/soupx.R'

rule assemble_depth:
    input:
        adata=expand('qc/{sample}/singlets.h5ad', sample=SAMPLES),
        matrices=expand('soupx/{sample}/corrected.mtx', sample=SAMPLES),
        genes=expand('soupx/{sample}/genes.tsv', sample=SAMPLES),
        barcodes=expand('soupx/{sample}/barcodes.tsv', sample=SAMPLES),
        rhos=expand('soupx/{sample}/rho.tsv', sample=SAMPLES)
    output:
        adata='analysis/prepared.h5ad',
        depth='analysis/depth.tsv',
        removed='analysis/empty_after_correction.tsv',
        summary='analysis/prepared.json'
    params: cfg=CFG, samples=list(SAMPLES)
    threads: CFG['resources']['analysis_threads']
    resources: mem_mb=CFG['resources']['analysis_mem_mb']
    container: C['analysis']
    log: 'logs/assemble_depth.log'
    benchmark: 'benchmarks/assemble_depth.tsv'
    script: 'scripts/scrna/assemble_depth.py'

rule scvi:
    input: 'analysis/prepared.h5ad'
    output:
        adata='analysis/scvi.h5ad',
        model=directory('models/scvi'),
        history='analysis/scvi_history.tsv',
        summary='analysis/scvi.json'
    params: cfg=CFG
    threads: CFG['resources']['analysis_threads']
    resources: mem_mb=CFG['resources']['model_mem_mb'], gpu=CFG['resources']['gpu']
    container: C['analysis']
    log: 'logs/scvi.log'
    benchmark: 'benchmarks/scvi.tsv'
    script: 'scripts/scrna/integrate.py'

rule celltypist_model:
    input: [CFG['celltypist']['model_path']] if CFG['celltypist']['model_path'] else []
    output: model='models/celltypist.pkl', metadata='models/celltypist.json'
    params: cfg=CFG
    container: C['analysis']
    log: 'logs/celltypist_model.log'
    script: 'scripts/scrna/celltypist_model.py'

rule celltypist:
    input: adata='analysis/scvi.h5ad', model='models/celltypist.pkl'
    output: adata='analysis/celltypist.h5ad', summary='analysis/celltypist.json'
    params: cfg=CFG
    threads: CFG['resources']['analysis_threads']
    resources: mem_mb=CFG['resources']['model_mem_mb']
    container: C['analysis']
    log: 'logs/celltypist.log'
    benchmark: 'benchmarks/celltypist.tsv'
    script: 'scripts/scrna/annotate.py'

rule scanvi:
    input: adata='analysis/celltypist.h5ad', model='models/scvi'
    output:
        adata='analysis/final.h5ad',
        model=directory('models/scanvi'),
        composition='analysis/composition.tsv',
        cells='analysis/cells.tsv.gz',
        history='analysis/scanvi_history.tsv',
        summary='analysis/scanvi.json'
    params: cfg=CFG
    threads: CFG['resources']['analysis_threads']
    resources: mem_mb=CFG['resources']['model_mem_mb'], gpu=CFG['resources']['gpu']
    container: C['analysis']
    log: 'logs/scanvi.log'
    benchmark: 'benchmarks/scanvi.tsv'
    script: 'scripts/scrna/scanvi.py'

rule report:
    input:
        adata='analysis/final.h5ad',
        composition='analysis/composition.tsv',
        depth='analysis/depth.tsv',
        qc=expand('qc/{sample}/summary.json', sample=SAMPLES),
        qc_cells=expand('qc/{sample}/cells.tsv.gz', sample=SAMPLES),
        ranks=expand('qc/{sample}/barcode_rank.tsv.gz', sample=SAMPLES),
        soup=expand('soupx/{sample}/summary.json', sample=SAMPLES),
        cr=[cr_metrics(s) for s in SAMPLES] if CFG['cellranger']['enabled'] else [],
        kb=expand('kb/{sample}/inspect.json', sample=SAMPLES) if FASTQ else [],
        summaries=['analysis/prepared.json', 'analysis/scvi.json', 'analysis/celltypist.json', 'analysis/scanvi.json', 'models/celltypist.json'],
        histories=['analysis/scvi_history.tsv', 'analysis/scanvi_history.tsv'],
        provenance='provenance/run.json',
        mapping_versions=(expand('provenance/kb_{sample}.json', sample=SAMPLES)
                          + (['provenance/kb_reference.json'] if BUILD else [])
                          + (expand('provenance/cellranger_{sample}.txt', sample=SAMPLES) if CFG['cellranger']['enabled'] else [])) if FASTQ else []
    output: html='report/report.html', figures=directory('report/figures')
    params: cfg=CFG, samples=list(SAMPLES)
    threads: CFG['resources']['analysis_threads']
    resources: mem_mb=CFG['resources']['analysis_mem_mb']
    container: C['analysis']
    log: 'logs/report.log'
    benchmark: 'benchmarks/report.tsv'
    script: 'scripts/scrna/report.py'
