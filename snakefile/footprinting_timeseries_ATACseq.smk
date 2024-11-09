
'''
Snakefile for ATAC-seq footprinting time series analysis

Usage:
    snakemake -s footprinting_timeseries_ATACseq.smk --configfile <path to config.yaml> --cores <int> --use-singularity
'''

## DON'T CHANGE BELOW THIS LINE ##

def load_experiment(file):
    experiment_dict = {}
    with open(file, "r") as f:
        for line in f:
            sample, bam, peak, group = line.strip().split("\t")
            if sample == "sample":
                continue
            experiment_dict[sample] = {"bam": bam, "peak": peak, "group": group}
    return experiment_dict
experiment_dict = load_experiment(config["experiment_table"])

workdir: config["workdir"]
samples = list(experiment_dict.keys())
bams = [experiment_dict[x]["bam"] for x in experiment_dict]
peaks = [experiment_dict[x]["peak"] for x in experiment_dict]

def deduplicate_preserve_order(lst):
    seen = set()
    return [x for x in lst if not (x in seen or seen.add(x))]
groups = deduplicate_preserve_order([experiment_dict[x]["group"] for x in experiment_dict])

rule all:
    input:
        footprints_bw = expand("TOBIAS/FootprintScores/{merge_group}.sort.rmdup_footprints.bw", merge_group = groups),
        bindetect_results = "TOBIAS/BINDetect/bindetect_results.xlsx",
        merge_overview = "TOBIAS/BINDetect/merge_overview.txt"

rule merge_peaks:
    container:
        "docker://quay.io/biocontainers/bedtools:2.24--1"
    input:
        expand("{peak}", peak = peaks)
    output:
        temp("merged_peaks.bed")
    benchmark:
        "benchmark/merge_peaks.txt"
    log:
        "log/merge_peaks.log"
    shell:
        '''
        cat {input} | sort -k1,1 -k2,2n | bedtools merge -i stdin -d 1000 | grep -v chrM > {output}
        '''

rule merge_bam:
    wildcard_constraints:
        merge_group = "|".join([re.escape(x) for x in groups])
    container:
        "docker://quay.io/biocontainers/samtools:1.18--h50ea8bc_1"
    input:
        bams = lambda wildcards: expand("{bam}", bam = [experiment_dict[x]["bam"] for x in experiment_dict if experiment_dict[x]["group"] == str(wildcards.merge_group)])
    output:
        bam = "{merge_group}.sort.rmdup.bam"
    threads:
        workflow.cores / 2
    benchmark:
        "benchmark/TOBIAS/bam_merge/{merge_group}.txt"
    log:
        "log/TOBIAS/bam_merge/{merge_group}.log"
    shell:
        '''
        samtools merge \
        -@ {threads} \
        -o {output.bam} \
        {input.bams} \
        > {log} 2>&1
        '''

rule ATACorrect:
    wildcard_constraints:
        merge_group = "|".join([re.escape(x) for x in groups])
    container:
        "docker://quay.io/biocontainers/tobias:0.16.0--py38h24c8ff8_0"
    input:
        bam = "{merge_group}.sort.rmdup.bam",
        peaks = "merged_peaks.bed"
    output:
        outdir = directory("TOBIAS/ATACorrect/{merge_group}"),
        corrected_bw = "TOBIAS/ATACorrect/{merge_group}/{merge_group}.sort.rmdup_corrected.bw"
    params:
        genome_fasta = config["genome_fasta"]
    threads:
        workflow.cores / 2
    benchmark:
        "benchmark/TOBIAS/ATACorrect/{merge_group}.txt"
    log:
        "log/TOBIAS/ATACorrect/{merge_group}.log"
    shell:
        '''
        TOBIAS ATACorrect \
        --bam {input.bam} \
        --genome {params.genome_fasta} \
        --peaks {input.peaks} \
        --outdir {output.outdir} \
        --cores {threads} \
        > {log} 2>&1
        '''

rule FootprintScores:
    wildcard_constraints:
        merge_group = "|".join([re.escape(x) for x in groups])
    container:
        "docker://quay.io/biocontainers/tobias:0.16.0--py38h24c8ff8_0"
    input:
        corrected_bw = "TOBIAS/ATACorrect/{merge_group}/{merge_group}.sort.rmdup_corrected.bw",
        peaks = "merged_peaks.bed"
    output:
        footprints_bw = "TOBIAS/FootprintScores/{merge_group}.sort.rmdup_footprints.bw"
    threads:
        workflow.cores / 2
    benchmark:
        "benchmark/TOBIAS/FootprintScores/{merge_group}.txt"
    log:
        "log/TOBIAS/FootprintScores/{merge_group}.log"
    shell:
        '''
        TOBIAS FootprintScores \
        --signal {input.corrected_bw} \
        --regions {input.peaks} \
        --output {output.footprints_bw} \
        --cores {threads} \
        > {log} 2>&1
        '''

rule BINDetect:
    container:
        "docker://quay.io/biocontainers/tobias:0.16.0--py38h24c8ff8_0"
    input:
        footprints_bw = expand("TOBIAS/FootprintScores/{merge_group}.sort.rmdup_footprints.bw", merge_group = groups),
        peaks = "merged_peaks.bed"
    output:
        outdir = directory("TOBIAS/BINDetect"),
        bindetect_results = "TOBIAS/BINDetect/bindetect_results.xlsx"
    params:
        cluster_motifs = config["cluster_motifs"],
        genome_fasta = config["genome_fasta"],
        footprints_bw_all = lambda wildcards: "TOBIAS/FootprintScores/{" + ",".join(groups) + "}.sort.rmdup_footprints.bw"
    threads:
        workflow.cores
    benchmark:
        "benchmark/TOBIAS/BINDetect.txt"
    log:
        "log/TOBIAS/BINDetect.log"
    shell:
        '''
        TOBIAS BINDetect \
        --motifs {params.cluster_motifs} \
        --signals {params.footprints_bw_all} \
        --genome {params.genome_fasta} \
        --peaks {input.peaks} \
        --outdir {output.outdir} \
        --time-series \
        --cores {threads} \
        > {log} 2>&1
        '''

rule merge_overview:
    input:
        bindetect_results = "TOBIAS/BINDetect/bindetect_results.xlsx"
    output:
        merge_overview = "TOBIAS/BINDetect/merge_overview.txt"
    threads:
        1
    benchmark:
        "benchmark/TOBIAS/merge_overview.txt"
    log:
        "log/TOBIAS/merge_overview.log"
    shell:
        '''
        cat TOBIAS/BINDetect/*/*_overview.txt | \
        grep -v "TFBS_chr" | \
        sort -k1,1 -k2,2n | \
        awk '!a[$0]++' > {output} 2> {log}
        '''

## DON'T CHANGE ABOVE THIS LINE ##
