
'''
Snakefile for RNA-seq preprocessing

Usage:
    snakemake -s preprocessing_RNAseq_single.smk --configfile <path to config.yaml> --cores <int> --use-singularity
'''

## DON'T CHANGE BELOW THIS LINE ##

workdir: config["workdir"]
star_index = config["star_index"]
samples = config["samples"]

rule all:
    input:
        multiqc = "multiqc/multiqc_report.html",
        bam = expand("star/{sample}/{sample}_Aligned.out.bam", sample = samples),
        bigwig = expand("bigwig/{sample}.bw", sample = samples),
        RnaSeqMetrics = expand("metrics/{sample}.picard.analysis.CollectRnaSeqMetrics", sample = samples)

rule qc:
    wildcard_constraints:
        sample = "|".join([re.escape(x) for x in samples])
    container:
        "docker://quay.io/biocontainers/fastp:0.23.4--hadf994f_2"
    input:
        R1 = "fastq/{sample}.fastq.gz"
    output:
        R1 = "fastp/{sample}.fastq.gz",
        json = "fastp/log/{sample}_fastp.json"
    params:
        html = "fastp/log/{sample}_fastp.html"
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
    wildcard_constraints:
        sample = "|".join([re.escape(x) for x in samples])
    container:
        "docker://quay.io/biocontainers/star:2.7.11a--h0033a41_0"
    input:
        R1 = "fastp/{sample}.fastq.gz"
    output:
        sam = "star/{sample}/{sample}_Aligned.out.sam",
        log = "star/{sample}/{sample}_Log.final.out"
    params:
        outdir = "star/{sample}/{sample}_"
    threads:
        1000
    benchmark:
        "benchmark/star_{sample}.txt"
    log:
        "log/star_{sample}.log"
    shell:
        "STAR --runThreadN {threads} --genomeDir {star_index} "
        "--outFilterMultimapNmax 1 "
        "--readFilesIn {input.R1} --readFilesCommand zcat --outFileNamePrefix {params.outdir} >& {log}"

rule sort:
    wildcard_constraints:
        sample = "|".join([re.escape(x) for x in samples])
    container:
        "docker://quay.io/biocontainers/samtools:1.18--h50ea8bc_1"
    input:
        "star/{sample}/{sample}_Aligned.out.sam"
    output:
        "star/{sample}/{sample}_Aligned.out.bam"
    threads:
        8
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
        "star/{sample}/{sample}_Aligned.out.bam"
    output:
        "bigwig/{sample}.bw"
    threads:
        8
    benchmark:
        "benchmark/bamCoverage_{sample}.txt"
    log:
        "log/bamCoverage_{sample}.log"
    shell:
        "bamCoverage -b {input} -o {output} -p {threads} --binSize 1 >& {log}"

rule makeRefFlat:
    wildcard_constraints:
        sample = "|".join([re.escape(x) for x in samples])
    container:
        "docker://quay.io/biocontainers/ucsc-gtftogenepred:469--h9b8f530_0"
    output:
        refFlat = "refFlat.txt"
    params:
        gtf = config["gtf"]
    threads:
        1
    benchmark:
        "benchmark/refFlat.txt"
    log:
        "log/refFlat.log"
    shell:
        "gtfToGenePred -genePredExt -geneNameAsName2 {params.gtf} refFlat.tmp >& {log} && "
        "cat refFlat.tmp | awk -F'\t' -v OFS='\t' '{{print $12,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10}}' > {output.refFlat}"

rule CollectRnaSeqMetrics:
    wildcard_constraints:
        sample = "|".join([re.escape(x) for x in samples])
    container:
        "docker://quay.io/biocontainers/picard:3.1.1--hdfd78af_0"
    input:
        bam = "star/{sample}/{sample}_Aligned.out.bam",
        refFlat = "refFlat.txt"
    output:
        RnaSeqMetrics = "metrics/{sample}.picard.analysis.CollectRnaSeqMetrics"
    threads:
        1
    benchmark:
        "benchmark/picard_CollectRnaSeqMetrics_{sample}.txt"
    log:
        "log/picard_CollectRnaSeqMetrics_{sample}.log"
    shell:
        "picard CollectRnaSeqMetrics -I {input.bam} -O {output.RnaSeqMetrics} --REF_FLAT {input.refFlat} --STRAND_SPECIFICITY NONE >& {log}"

rule multiqc:
    container:
        "docker://multiqc/multiqc:latest"
    input:
        json = expand("fastp/log/{sample}_fastp.json", sample = samples),
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
