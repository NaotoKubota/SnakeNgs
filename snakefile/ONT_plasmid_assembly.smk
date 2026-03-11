
'''
Snakefile for ONT sequencing data analysis for plasmid assembly

Usage:
    snakemake -s ONT_plasmid_assembly.smk --configfile <path to config.yaml> --cores <int> --use-singularity
'''

workdir: config["workdir"]
samples_dict = config["samples"]
samples = list(samples_dict.keys())
reference = config["reference"]

include: "common/containers.smk"

wildcard_constraints:
    sample = "|".join([re.escape(x) for x in samples])

rule all:
    input:
        expand("assembly/medaka/{sample}/consensus.fasta", sample=samples),
        expand("assembly/medaka_denovo/{sample}/consensus.fasta", sample=samples),
        expand("sequali/{sample}.fastq.gz.json", sample=samples),
        expand("sequali/{sample}.fastq.gz.html", sample=samples),
        "multiqc/multiqc_report.html"

rule concat_fastq:
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
    container:
        CONTAINERS["sequali"]
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
    container:
        CONTAINERS["flye"]
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
    container:
        CONTAINERS["medaka"]
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
    container:
        CONTAINERS["medaka"]
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
        CONTAINERS["multiqc"]
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
