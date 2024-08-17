
'''
Snakefile for differential analysis in ATAC-seq

Usage:
    snakemake -s differential_ATACseq.smk --configfile <path to config.yaml> --cores <int> --use-singularity
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
samples = experiment_dict.keys()
bams = [experiment_dict[x]["bam"] for x in experiment_dict]
peaks = [experiment_dict[x]["peak"] for x in experiment_dict]
groups = ["Ref", "Alt"]

rule all:
    input:
        deseq2_results = "deseq2_results.tsv"

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

rule make_saf:
    input:
        "merged_peaks.bed"
    output:
        "merged_peaks.saf"
    benchmark:
        "benchmark/make_saf.txt"
    log:
        "log/make_saf.log"
    shell:
        '''
        cat {input} | \
        awk -F"\t" -v OFS="\t" '{{print $1":"$2"-"$3,$1,$2,$3,"+"}}' > {output}
        '''

rule read_count:
    container:
        "docker://quay.io/biocontainers/subread:2.0.6--he4a0461_0"
    input:
        bam = expand("{bam}", bam = bams),
        saf = "merged_peaks.saf"
    output:
        temp("all_counts.tsv")
    threads:
        16
    benchmark:
        "benchmark/read_count.txt"
    log:
        "log/all_counts.log"
    shell:
        '''
        featureCounts \
        -a {input.saf} \
        -o {output} \
        -F SAF \
        -T {threads} \
        -O \
        -p \
        {input.bam} \
        > {log} 2>&1
        '''

rule make_count_table:
    input:
        "all_counts.tsv"
    output:
        "all_counts_for_DESeq2.tsv"
    benchmark:
        "benchmark/make_count_table.txt"
    log:
        "log/make_count_table.log"
    shell:
        '''
        cat {input} | \
        cut -f 1,7- | \
        sed -e 1d | \
        sed -e 's@Geneid@PeakID@' \
        > {output}
        '''

rule deseq2:
    container:
        "docker://quay.io/biocontainers/bioconductor-deseq2:1.42.0--r43hf17093f_0"
    input:
        experiment_table = config["experiment_table"],
        counts = "all_counts_for_DESeq2.tsv"
    output:
        "deseq2_results.tsv"
    benchmark:
        "benchmark/deseq2.txt"
    log:
        "log/deseq2.log"
    shell:
        '''
        Rscript /rhome/naotok/ZhengLab_pipeline/utils/deseq2.R \
        {input.experiment_table} {input.counts} {output} \
        > {log} 2>&1
        '''

## DON'T CHANGE ABOVE THIS LINE ##
