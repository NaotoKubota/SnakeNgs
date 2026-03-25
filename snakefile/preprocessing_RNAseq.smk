
'''
Snakefile for RNA-seq preprocessing (paired-end and single-end)

Usage:
    snakemake -s preprocessing_RNAseq.smk --configfile <path to config.yaml> --cores <int> --use-singularity

config.yaml must contain:
    layout: "paired" or "single"
'''

## DON'T CHANGE BELOW THIS LINE ##

workdir: config["workdir"]
star_index = config["star_index"]
samples = config["samples"]
gtf = config["gtf"]
layout = config["layout"]

wildcard_constraints:
    sample = "|".join([re.escape(x) for x in samples])

include: "common/containers.smk"
include: "common/functions.smk"
include: "common/sort_bam.smk"
include: "common/bigwig.smk"
include: "common/rnaseq_metrics.smk"

if layout == "paired":

    rule all:
        input:
            multiqc = "multiqc/multiqc_report.html",
            bam = expand("star/{sample}/{sample}_Aligned.out.bam", sample = samples),
            bigwig = expand("bigwig/{sample}.bw", sample = samples),
            RnaSeqMetrics = expand("metrics/{sample}.picard.analysis.CollectRnaSeqMetrics", sample = samples),
            InsertSizeMetrics = expand("metrics/{sample}.picard.analysis.CollectInsertSizeMetrics", sample = samples)

    rule qc:
        container:
            CONTAINERS["fastp"]
        input:
            R1 = "fastq/{sample}_1.fastq.gz",
            R2 = "fastq/{sample}_2.fastq.gz"
        output:
            R1 = "fastp/{sample}_1.fastq.gz",
            R2 = "fastp/{sample}_2.fastq.gz",
            json = "fastp/log/{sample}.json"
        params:
            html = "fastp/log/{sample}.html"
        threads:
            8
        benchmark:
            "benchmark/fastp_{sample}.txt"
        log:
            "log/fastp_{sample}.log"
        shell:
            "fastp -i {input.R1} -I {input.R2} "
            "-o {output.R1} -O {output.R2} -w {threads} "
            "-h {params.html} -j {output.json} >& {log}"

    rule mapping:
        container:
            CONTAINERS["star"]
        input:
            R1 = "fastp/{sample}_1.fastq.gz",
            R2 = "fastp/{sample}_2.fastq.gz",
        output:
            sam = "star/{sample}/{sample}_Aligned.out.sam",
            log = "star/{sample}/{sample}_Log.final.out"
        params:
            outdir = "star/{sample}/{sample}_"
        threads:
            workflow.cores
        benchmark:
            "benchmark/star_{sample}.txt"
        log:
            "log/star_{sample}.log"
        shell:
            "STAR --runThreadN {threads} --genomeDir {star_index} "
            "--outFilterMultimapNmax 1 "
            "--readFilesIn {input.R1} {input.R2} --readFilesCommand zcat --outFileNamePrefix {params.outdir} >& {log}"

    rule CollectInsertSizeMetrics:
        container:
            CONTAINERS["picard"]
        input:
            bam = "star/{sample}/{sample}_Aligned.out.bam"
        output:
            InsertSizeMetrics = "metrics/{sample}.picard.analysis.CollectInsertSizeMetrics",
            InsertSizeMetrics_pdf = "metrics/{sample}.picard.analysis.CollectInsertSizeMetrics.pdf"
        threads:
            1
        benchmark:
            "benchmark/picard_CollectInsertSizeMetrics_{sample}.txt"
        log:
            "log/picard_CollectInsertSizeMetrics_{sample}.log"
        shell:
            "picard CollectInsertSizeMetrics -I {input.bam} -O {output.InsertSizeMetrics} --Histogram_FILE {output.InsertSizeMetrics_pdf} >& {log}"

    rule multiqc:
        container:
            CONTAINERS["multiqc"]
        input:
            json = expand("fastp/log/{sample}.json", sample = samples),
            starlog = expand("star/{sample}/{sample}_Log.final.out", sample = samples),
            RnaSeqMetrics = expand("metrics/{sample}.picard.analysis.CollectRnaSeqMetrics", sample = samples),
            InsertSizeMetrics = expand("metrics/{sample}.picard.analysis.CollectInsertSizeMetrics", sample = samples)
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
            cp {input.json} {input.starlog} {input.RnaSeqMetrics} {input.InsertSizeMetrics} multiqc/log && \
            cat /usr/local/lib/python3.12/site-packages/multiqc/config_defaults.yaml | sed -e '$afastp:\\n  s_name_filenames: true' > multiqc/multiqc_config.yaml && \
            multiqc --config multiqc/multiqc_config.yaml -o multiqc/ multiqc/log >& {log} && \
            rm -rf multiqc/log
            """

else:

    rule all:
        input:
            multiqc = "multiqc/multiqc_report.html",
            bam = expand("star/{sample}/{sample}_Aligned.out.bam", sample = samples),
            bigwig = expand("bigwig/{sample}.bw", sample = samples),
            RnaSeqMetrics = expand("metrics/{sample}.picard.analysis.CollectRnaSeqMetrics", sample = samples)

    rule qc:
        container:
            CONTAINERS["fastp"]
        input:
            R1 = "fastq/{sample}.fastq.gz"
        output:
            R1 = "fastp/{sample}.fastq.gz",
            json = "fastp/log/{sample}.json"
        params:
            html = "fastp/log/{sample}.html"
        threads:
            8
        benchmark:
            "benchmark/fastp_{sample}.txt"
        log:
            "log/fastp_{sample}.log"
        shell:
            "fastp -i {input.R1} "
            "-o {output.R1} -w {threads} "
            "-h {params.html} -j {output.json} >& {log}"

    rule mapping:
        container:
            CONTAINERS["star"]
        input:
            R1 = "fastp/{sample}.fastq.gz"
        output:
            sam = "star/{sample}/{sample}_Aligned.out.sam",
            log = "star/{sample}/{sample}_Log.final.out"
        params:
            outdir = "star/{sample}/{sample}_"
        threads:
            workflow.cores
        benchmark:
            "benchmark/star_{sample}.txt"
        log:
            "log/star_{sample}.log"
        shell:
            "STAR --runThreadN {threads} --genomeDir {star_index} "
            "--outFilterMultimapNmax 1 "
            "--readFilesIn {input.R1} --readFilesCommand zcat --outFileNamePrefix {params.outdir} >& {log}"

    rule multiqc:
        container:
            CONTAINERS["multiqc"]
        input:
            json = expand("fastp/log/{sample}.json", sample = samples),
            starlog = expand("star/{sample}/{sample}_Log.final.out", sample = samples),
            RnaSeqMetrics = expand("metrics/{sample}.picard.analysis.CollectRnaSeqMetrics", sample = samples)
        output:
            "multiqc/multiqc_report.html"
        benchmark:
            "benchmark/multiqc.txt"
        log:
            "log/multiqc.log"
        shell:
            "rm -rf multiqc && "
            "mkdir -p multiqc/log && "
            "cp {input.json} {input.starlog} {input.RnaSeqMetrics} multiqc/log && "
            "multiqc -o multiqc/ multiqc/log >& {log} && "
            "rm -rf multiqc/log"

## DON'T CHANGE ABOVE THIS LINE ##
