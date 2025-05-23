
'''
Snakefile for long-read RNA-seq preprocessing

Usage:
    snakemake -s preprocessing_RNAseq_long.smk --configfile <path to config.yaml> --cores <int> --use-singularity
'''

workdir: config["workdir"]
samples = config["samples"]

rule all:
    input:
        bam = expand("minimap2/{sample}/{sample}_Aligned.out.bam", sample = samples),
        bigwig = expand("bigwig/{sample}.bw", sample = samples),
        isoquant_gene_counts = "isoquant/OUT/OUT.transcript_model_grouped_tpm.tsv",
        multiqc = "multiqc/multiqc_report.html"

rule gtf2bed:
    container:
        "docker://quay.io/biocontainers/minimap2:2.28--h577a1d6_4"
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
    wildcard_constraints:
        sample = "|".join([re.escape(x) for x in samples])
    container:
        "docker://quay.io/biocontainers/sequali:0.12.0--py311haab0aaa_1"
    input:
        R1 = "fastq/{sample}.fastq.gz"
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
    wildcard_constraints:
        sample = "|".join([re.escape(x) for x in samples])
    container:
        "docker://quay.io/biocontainers/minimap2:2.28--h577a1d6_4"
    input:
        R1 = "fastq/{sample}.fastq.gz",
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
    wildcard_constraints:
        sample = "|".join([re.escape(x) for x in samples])
    container:
        "docker://quay.io/biocontainers/samtools:1.18--h50ea8bc_1"
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
    wildcard_constraints:
        sample = "|".join([re.escape(x) for x in samples])
    container:
        "docker://quay.io/biocontainers/deeptools:3.5.4--pyhdfd78af_1"
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
        "docker://quay.io/biocontainers/isoquant:3.6.3--hdfd78af_0"
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
        multiqc --config multiqc/multiqc_config.yaml -o multiqc/ multiqc/log >& {log} && \
        rm -rf multiqc/log
        """
