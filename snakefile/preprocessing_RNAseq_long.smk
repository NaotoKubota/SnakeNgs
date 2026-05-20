
'''
Snakefile for long-read RNA-seq preprocessing

Usage:
    snakemake -s preprocessing_RNAseq_long.smk --configfile <path to config.yaml> --cores <int> --use-singularity
'''

workdir: config["workdir"]

include: "common/containers.smk"
include: "common/functions.smk"

samples, samples_dict = normalize_samples(config["samples"])

wildcard_constraints:
    sample = "|".join([re.escape(x) for x in samples])

# When samples is given as a dict {sample: [run1, run2, ...]}, fastq files for
# each run are concatenated into fastq_merged/ before downstream rules. When
# samples is a plain list, the merge rule is not emitted and rules read
# fastq/ directly.
def _fastq_input_single(wildcards):
    base = "fastq_merged" if samples_dict else "fastq"
    return f"{base}/{wildcards.sample}.fastq.gz"

rule all:
    input:
        bam = expand("minimap2/{sample}/{sample}_Aligned.out.bam", sample = samples),
        bigwig = expand("bigwig/{sample}.bw", sample = samples),
        isoquant_gene_counts = "isoquant/OUT/OUT.transcript_model_grouped_tpm.tsv",
        multiqc = "multiqc/multiqc_report.html"

if samples_dict:

    rule merge_fastq_single:
        input:
            R1 = lambda wc: [f"fastq/{run}.fastq.gz" for run in samples_dict[wc.sample]]
        output:
            R1 = "fastq_merged/{sample}.fastq.gz"
        benchmark:
            "benchmark/merge_fastq_{sample}.txt"
        log:
            "log/merge_fastq_{sample}.log"
        shell:
            "cat {input.R1} > {output.R1} 2> {log}"

rule gtf2bed:
    container:
        CONTAINERS["minimap2"]
    input:
        gtf = config["gtf"]
    output:
        annobed = "annotation/anno.bed"
    threads:
        1
    benchmark:
        "benchmark/gtf2bed.txt"
    log:
        "log/gtf2bed.log"
    shell:
        """
        paftools.js gff2bed {input.gtf} > {output.annobed} 2> {log}
        """

rule qc:
    container:
        CONTAINERS["sequali"]
    input:
        R1 = _fastq_input_single
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
        {input.R1} \
        >& {log}
        """

rule mapping:
    container:
        CONTAINERS["minimap2"]
    input:
        R1 = _fastq_input_single,
        annobed = "annotation/anno.bed"
    output:
        sam = "minimap2/{sample}/{sample}_Aligned.out.sam"
    threads:
        workflow.cores / 2
    benchmark:
        "benchmark/minimap2_{sample}.txt"
    log:
        "log/minimap2_{sample}.log"
    shell:
        """
        minimap2 \
        -t {threads} \
        -ax {config[preset]} \
        --MD \
        --cs \
        -u f \
        --junc-bed {input.annobed} \
        {config[minimap2_options]} \
        {config[genome_fasta]} \
        {input.R1} > {output.sam} 2> {log}
        """

rule sort:
    container:
        CONTAINERS["samtools"]
    input:
        "minimap2/{sample}/{sample}_Aligned.out.sam"
    output:
        "minimap2/{sample}/{sample}_Aligned.out.bam"
    threads:
        workflow.cores / 4
    benchmark:
        "benchmark/samtools_{sample}.txt"
    log:
        "log/samtools_{sample}.log"
    shell:
        "samtools sort -@ {threads} -O bam -o {output} {input} >& {log} && "
        "samtools index {output} && "
        "rm -rf {input}"

rule bigwig:
    container:
        CONTAINERS["deeptools"]
    input:
        "minimap2/{sample}/{sample}_Aligned.out.bam"
    output:
        "bigwig/{sample}.bw"
    threads:
        workflow.cores / 4
    benchmark:
        "benchmark/bamCoverage_{sample}.txt"
    log:
        "log/bamCoverage_{sample}.log"
    shell:
        "bamCoverage -b {input} -o {output} -p {threads} --binSize 1 >& {log}"

rule isoquant:
    container:
        CONTAINERS["isoquant"]
    input:
        bams = expand("minimap2/{sample}/{sample}_Aligned.out.bam", sample = samples),
        gtf = config["gtf"]
    output:
        isoquant = "isoquant/OUT/OUT.transcript_model_grouped_tpm.tsv"
    params:
        labels = " ".join(samples)
    threads:
        workflow.cores
    benchmark:
        "benchmark/isoquant.txt"
    log:
        "log/isoquant.log"
    shell:
        """
        isoquant.py \
        --reference {config[genome_fasta]} \
        --genedb {input.gtf} \
        --bam {input.bams} \
        --data_type {config[data_type]} \
        --output isoquant \
        --threads {threads} \
        --labels {params.labels} \
        {config[isoquant_options]} \
        >& {log}
        """

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
