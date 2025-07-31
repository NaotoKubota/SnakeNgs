
'''
Snakefile for ONT sequencing data analysis for plasmid assembly

Usage:
    snakemake -s ONT_plasmid_assembly.smk --configfile <path to config.yaml> --cores <int> --use-singularity
'''

workdir: config["workdir"]
samples_dict = config["samples"]
samples = list(samples_dict.keys())
reference = config["reference"]

rule all:
    input:
        expand("assembly/medaka/{sample}/consensus.fasta", sample=samples),
        expand("assembly/medaka_denovo/{sample}/consensus.fasta", sample=samples),
        expand("sequali/{sample}.fastq.gz.json", sample=samples),
        expand("sequali/{sample}.fastq.gz.html", sample=samples),
        "multiqc/multiqc_report.html"

rule concat_fastq:
    wildcard_constraints:
        sample = "|".join([re.escape(x) for x in samples])
    input:
        input_dir = lambda wildcards: samples_dict[wildcards.sample]
    output:
        fastq = "fastq/{sample}.fastq.gz"
    benchmark:
        "benchmark/concat_fastq_{sample}.txt"
    shell:
        '''
        cat {input.input_dir}/*fastq.gz > {output}
        '''

rule qc:
    wildcard_constraints:
        sample = "|".join([re.escape(x) for x in samples])
    container:
        "docker://quay.io/biocontainers/sequali:0.12.0--py311haab0aaa_1"
    input:
        fastq = "fastq/{sample}.fastq.gz"
    output:
        json = "sequali/{sample}.fastq.gz.json",
        html = "sequali/{sample}.fastq.gz.html"
    threads:
        2
    benchmark:
        "benchmark/sequali_{sample}.txt"
    log:
        "log/sequali_{sample}.log"
    shell:
        """
        sequali \
        --outdir sequali \
        {input.fastq} \
        >& {log}
        """

rule flye:
    wildcard_constraints:
        sample = "|".join([re.escape(x) for x in samples])
    container:
        "docker://staphb/flye:2.9.6"
    input:
        fastq = "fastq/{sample}.fastq.gz"
    output:
        assembly_fasta = "assembly/flye/{sample}/assembly.fasta"
    params:
        output_dir = "assembly/flye/{sample}"
    threads:
        workflow.cores / len(samples)
    benchmark:
        "benchmark/flye_{sample}.txt"
    log:
        "log/flye_{sample}.log"
    shell:
        '''
        flye \
        --nano-hq \
        {input} \
        --out-dir {params.output_dir} \
        --threads {threads} \
        --meta \
        >& {log}
        '''

rule medaka_consensus:
    wildcard_constraints:
        sample = "|".join([re.escape(x) for x in samples])
    container:
        "docker://ontresearch/medaka:shac4e11bfa4e65668b28739ba32edc3af12baf7574"
    input:
        fastq = "fastq/{sample}.fastq.gz",
        reference = reference
    output:
        consensus_fasta = "assembly/medaka/{sample}/consensus.fasta"
    params:
        output_dir = "assembly/medaka/{sample}"
    threads:
        workflow.cores / len(samples)
    benchmark:
        "benchmark/medaka_consensus_{sample}.txt"
    log:
        "log/medaka_consensus_{sample}.log"
    shell:
        '''
        medaka_consensus \
        -i {input.fastq} \
        -d {input.reference} \
        -o {params.output_dir} \
        -t {threads} \
        -g \
        --bacteria \
        >& {log}
        '''

rule medaka_consensus_denovo:
    wildcard_constraints:
        sample = "|".join([re.escape(x) for x in samples])
    container:
        "docker://ontresearch/medaka:shac4e11bfa4e65668b28739ba32edc3af12baf7574"
    input:
        fastq = "fastq/{sample}.fastq.gz",
        assembly_fasta = "assembly/flye/{sample}/assembly.fasta"
    output:
        consensus_fasta = "assembly/medaka_denovo/{sample}/consensus.fasta"
    params:
        output_dir = "assembly/medaka_denovo/{sample}"
    threads:
        workflow.cores / len(samples)
    benchmark:
        "benchmark/medaka_consensus_denovo_{sample}.txt"
    log:
        "log/medaka_consensus_denovo_{sample}.log"
    shell:
        '''
        medaka_consensus \
        -i {input.fastq} \
        -d {input.assembly_fasta} \
        -o {params.output_dir} \
        -t {threads} \
        -g \
        --bacteria \
        >& {log}
        '''

rule multiqc:
    container:
        "docker://multiqc/multiqc:v1.28"
    input:
        json = expand("sequali/{sample}.fastq.gz.json", sample = samples)
    output:
        "multiqc/multiqc_report.html"
    benchmark:
        "benchmark/multiqc.txt"
    log:
        "log/multiqc.log"
    shell:
        """
        rm -rf multiqc && \
        mkdir -p multiqc/log && \
        cp {input.json} multiqc/log && \
        multiqc -o multiqc/ multiqc/log >& {log} && \
        rm -rf multiqc/log
        """
