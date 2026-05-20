'''
Snakefile for CLIP-seq preprocessing (HITS-CLIP, iCLIP-seq, etc.)

Usage:
    snakemake -s preprocessing_CLIPseq.smk --configfile <path to config.yaml> --cores <int> --use-singularity
'''

## DON'T CHANGE BELOW THIS LINE ##

workdir: config["workdir"]
star_index = config["star_index"]
fastp_args = config["fastp_args"]
outFilterMultimapNmax = config["outFilterMultimapNmax"]
normalize_bigwig = config.get("normalize_bigwig", False)

include: "common/containers.smk"
include: "common/functions.smk"

samples, samples_dict = normalize_samples(config["samples"])

wildcard_constraints:
    sample = "|".join([re.escape(x) for x in samples])

# When samples is given as a dict {sample: [run1, run2, ...]}, fastq files for
# each run are concatenated into fastq_merged/ before QC. When samples is a
# plain list, the merge rule is not emitted and rule qc reads fastq/ directly.
def _qc_input_single(wildcards):
    base = "fastq_merged" if samples_dict else "fastq"
    return {"R1": f"{base}/{wildcards.sample}.fastq.gz"}

rule all:
    input:
        multiqc = "multiqc_preprocessing/multiqc_report.html",
        bam = expand("star/{sample}/{sample}_Aligned.rmdup.out.bam", sample = samples),
        bigwig = expand("bigwig/{sample}.bw", sample = samples)

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

rule qc:
    container:
        CONTAINERS["fastp"]
    input:
        unpack(_qc_input_single)
    output:
        R1 = "fastp/{sample}.fastq.gz",
        json = "fastp/log/{sample}.json"
    params:
        html = "fastp/log/{sample}.html",
        fastp_args = fastp_args
    threads:
        8
    benchmark:
        "benchmark/fastp_{sample}.txt"
    log:
        "log/fastp_{sample}.log"
    shell:
        "fastp -i {input.R1} "
        "-o {output.R1} -w {threads} {params.fastp_args} "
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
        outdir = "star/{sample}/{sample}_",
        outFilterMultimapNmax = outFilterMultimapNmax
    threads:
        workflow.cores
    benchmark:
        "benchmark/star_{sample}.txt"
    log:
        "log/star_{sample}.log"
    shell:
        "STAR --runThreadN {threads} --genomeDir {star_index} "
        "--outFilterMultimapNmax {params.outFilterMultimapNmax} "
        "--readFilesIn {input.R1} --readFilesCommand zcat --outFileNamePrefix {params.outdir} >& {log}"

rule sort:
    container:
        CONTAINERS["samtools"]
    input:
        "star/{sample}/{sample}_Aligned.out.sam"
    output:
        "star/{sample}/{sample}_Aligned.out.bam"
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

rule markdup:
    container:
        CONTAINERS["picard"]
    input:
        "star/{sample}/{sample}_Aligned.out.bam"
    output:
        "star/{sample}/{sample}_Aligned.rmdup.out.bam"
    threads:
        8
    benchmark:
        "benchmark/picard_{sample}.txt"
    log:
        "log/picard_{sample}.log"
    shell:
        "picard MarkDuplicates "
        "-I {input} -O {output} "
        "-M {log} "
        "--REMOVE_DUPLICATES true "
        "--TMP_DIR star/tmp "
        ">& {log} && "
        "rm -rf star/*.sort.bam.bai"

rule index:
    container:
        CONTAINERS["samtools"]
    input:
        "star/{sample}/{sample}_Aligned.rmdup.out.bam"
    output:
        "star/{sample}/{sample}_Aligned.rmdup.out.bam.bai"
    shell:
        "samtools index {input}"

rule bigwig:
    container:
        CONTAINERS["deeptools"]
    input:
        bam = "star/{sample}/{sample}_Aligned.rmdup.out.bam",
        bai = "star/{sample}/{sample}_Aligned.rmdup.out.bam.bai"
    output:
        "bigwig/{sample}.bw"
    params:
        normalize = "--normalizeUsing CPM" if normalize_bigwig else ""
    threads:
        8
    benchmark:
        "benchmark/bamCoverage_{sample}.txt"
    log:
        "log/bamCoverage_{sample}.log"
    shell:
        "bamCoverage -b {input.bam} -o {output} -p {threads} --binSize 1 {params.normalize} >& {log}"

rule multiqc:
    container:
        CONTAINERS["multiqc"]
    input:
        json = expand("fastp/log/{sample}.json", sample = samples),
        starlog = expand("star/{sample}/{sample}_Log.final.out", sample = samples),
        picardlog = expand("log/picard_{sample}.log", sample = samples)
    output:
        "multiqc_preprocessing/multiqc_report.html"
    benchmark:
        "benchmark/multiqc.txt"
    log:
        "log/multiqc.log"
    shell:
        "rm -rf multiqc_preprocessing && "
        "mkdir -p multiqc_preprocessing/log && "
        "cp {input.json} {input.starlog} {input.picardlog} multiqc_preprocessing/log && "
        "multiqc -o multiqc_preprocessing/ multiqc_preprocessing/log >& {log} && "
        "rm -rf multiqc_preprocessing/log"

## DON'T CHANGE ABOVE THIS LINE ##
