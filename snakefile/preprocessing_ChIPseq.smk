
'''
Snakefile for ChIP-seq preprocessing (paired-end and single-end)

Usage:
    snakemake -s preprocessing_ChIPseq.smk --configfile <path to config.yaml> --cores <int> --use-singularity

config.yaml must contain:
    layout: "paired" or "single"
'''

## DON'T CHANGE BELOW THIS LINE ##

workdir: config["workdir"]
bowtie2_index = config["bowtie2_index"]
samples = config["samples"]
layout = config["layout"]

wildcard_constraints:
    sample = "|".join([re.escape(x) for x in samples])

include: "common/containers.smk"
include: "common/sort_bowtie2.smk"
include: "common/markdup.smk"
include: "common/index_bam.smk"
include: "common/plotFingerprint.smk"
include: "common/bigwig_chipseq.smk"

rule all:
    input:
        multiqc = "multiqc_preprocessing/multiqc_report.html",
        bam = expand("bowtie2/{sample}.sort.rmdup.bam", sample = samples),
        bai = expand("bowtie2/{sample}.sort.rmdup.bam.bai", sample = samples),
        bigwig = expand("bigwig/{sample}.bw", sample = samples),
        png = "plotFingerprint/fingerprint.png",
        qc = "plotFingerprint/fingerprint.qc.txt",
        tab = "plotFingerprint/fingerprint.tab"

if layout == "paired":

    rule qc:
        container:
            CONTAINERS["fastp"]
        input:
            R1 = "fastq/{sample}_1.fastq.gz",
            R2 = "fastq/{sample}_2.fastq.gz"
        output:
            R1 = "fastp/{sample}_1.fq.gz",
            R2 = "fastp/{sample}_2.fq.gz",
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

    rule mapping_bowtie2:
        container:
            CONTAINERS["bowtie2"]
        input:
            R1 = "fastp/{sample}_1.fq.gz",
            R2 = "fastp/{sample}_2.fq.gz"
        output:
            sam = temp("bowtie2/{sample}.sam")
        params:
            bowtie2_args = config["bowtie2_args"]
        threads:
            1000
        benchmark:
            "benchmark/bowtie2_{sample}.txt"
        log:
            "log/bowtie2/{sample}.log"
        shell:
            "bowtie2 -p {threads} -x {bowtie2_index} "
            "{params.bowtie2_args} "
            "-1 {input.R1} -2 {input.R2} -S {output.sam} >& {log}"

    rule CollectInsertSizeMetrics:
        container:
            CONTAINERS["picard"]
        input:
            bam = "bowtie2/{sample}.sort.rmdup.bam",
            bai = "bowtie2/{sample}.sort.rmdup.bam.bai"
        output:
            InsertSizeMetrics = "metrics/CollectInsertSizeMetrics/{sample}",
            InsertSizeMetrics_pdf = "metrics/CollectInsertSizeMetrics/{sample}.pdf"
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
            bowtie2log = expand("log/bowtie2/{sample}.log", sample = samples),
            picardlog = expand("log/picard/{sample}.log", sample = samples),
            InsertSizeMetrics = expand("metrics/CollectInsertSizeMetrics/{sample}", sample = samples),
            plotFingerprintlog_qc = "plotFingerprint/fingerprint.qc.txt",
            plotFingerprintlog_tab = "plotFingerprint/fingerprint.tab"
        output:
            "multiqc_preprocessing/multiqc_report.html"
        benchmark:
            "benchmark/multiqc.txt"
        log:
            "log/multiqc.log"
        shell:
            """
            rm -rf multiqc_preprocessing && \
            mkdir -p multiqc_preprocessing/log/picard/CollectInsertSizeMetrics && \
            cp {input.json} {input.bowtie2log} {input.plotFingerprintlog_qc} {input.plotFingerprintlog_tab} multiqc_preprocessing/log && \
            cp {input.picardlog} multiqc_preprocessing/log/picard && \
            cp {input.InsertSizeMetrics} multiqc_preprocessing/log/picard/CollectInsertSizeMetrics && \
            cat /usr/local/lib/python3.12/site-packages/multiqc/config_defaults.yaml | \
            sed -e '$afastp:\\n  s_name_filenames: true' -e '$apicard_config:\\n  s_name_filenames: true' \
            > multiqc_preprocessing/multiqc_config.yaml && \
            multiqc --config multiqc_preprocessing/multiqc_config.yaml -o multiqc_preprocessing/ multiqc_preprocessing/log >& {log} && \
            rm -rf multiqc_preprocessing/log
            """

else:

    rule qc:
        container:
            CONTAINERS["fastp"]
        input:
            R1 = "fastq/{sample}.fastq.gz"
        output:
            R1 = "fastp/{sample}.fq.gz",
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

    rule mapping_bowtie2:
        container:
            CONTAINERS["bowtie2"]
        input:
            R1 = "fastp/{sample}.fq.gz"
        output:
            sam = temp("bowtie2/{sample}.sam")
        params:
            bowtie2_args = config["bowtie2_args"]
        threads:
            1000
        benchmark:
            "benchmark/bowtie2_{sample}.txt"
        log:
            "log/bowtie2/{sample}.log"
        shell:
            "bowtie2 -p {threads} -x {bowtie2_index} "
            "{params.bowtie2_args} "
            "-U {input.R1} -S {output.sam} >& {log}"

    rule multiqc:
        container:
            CONTAINERS["multiqc"]
        input:
            json = expand("fastp/log/{sample}.json", sample = samples),
            bowtie2log = expand("log/bowtie2/{sample}.log", sample = samples),
            picardlog = expand("log/picard/{sample}.log", sample = samples),
            plotFingerprintlog_qc = "plotFingerprint/fingerprint.qc.txt",
            plotFingerprintlog_tab = "plotFingerprint/fingerprint.tab"
        output:
            "multiqc_preprocessing/multiqc_report.html"
        benchmark:
            "benchmark/multiqc.txt"
        log:
            "log/multiqc.log"
        shell:
            "rm -rf multiqc_preprocessing && "
            "mkdir -p multiqc_preprocessing/log && "
            "cp {input.json} {input.bowtie2log} {input.picardlog} {input.plotFingerprintlog_qc} {input.plotFingerprintlog_tab} multiqc_preprocessing/log && "
            "multiqc -o multiqc_preprocessing/ multiqc_preprocessing/log >& {log} && "
            "rm -rf multiqc_preprocessing/log"

## DON'T CHANGE ABOVE THIS LINE ##
