'''
Snakefile for RNA-seq BAM to BigWig conversion

Usage:
    snakemake -s bam2bw_RNAseq.smk --configfile <path to config.yaml> --cores <int> --use-singularity

config.yaml must contain:
    workdir: /path/to/output
    experiment_table: /path/to/experiment_table.tsv

experiment_table.tsv must contain:
    sample\tbam
'''

## DON'T CHANGE BELOW THIS LINE ##

import csv
import re


def load_experiment(file):
    experiment_dict = {}
    with open(file, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("experiment_table is empty or missing header")

        required_columns = {"sample", "bam"}
        missing_columns = required_columns - set(reader.fieldnames)
        if missing_columns:
            raise ValueError(
                "experiment_table must contain columns: sample and bam. "
                f"Missing: {', '.join(sorted(missing_columns))}"
            )

        for row in reader:
            sample = row["sample"].strip()
            bam = row["bam"].strip()

            if not sample:
                raise ValueError("experiment_table contains an empty sample value")
            if not bam:
                raise ValueError(f"experiment_table contains an empty bam value for sample: {sample}")

            experiment_dict[sample] = {"bam": bam}

    if not experiment_dict:
        raise ValueError("No samples found in experiment_table")

    return experiment_dict


experiment_dict = load_experiment(config["experiment_table"])

workdir: config["workdir"]
samples = experiment_dict.keys()

wildcard_constraints:
    sample = "|".join([re.escape(x) for x in samples])

include: "common/containers.smk"

rule all:
    input:
        bigwig = expand("bigwig/{sample}.bw", sample=samples)

rule bam2bw:
    container:
        CONTAINERS["deeptools"]
    input:
        bam=lambda wildcards: experiment_dict[wildcards.sample]["bam"]
    output:
        "bigwig/{sample}.bw"
    threads:
        workflow.cores / 4
    benchmark:
        "benchmark/bamCoverage_{sample}.txt"
    log:
        "log/bamCoverage_{sample}.log"
    shell:
        "bamCoverage -b {input.bam} -o {output} -p {threads} --binSize 1 >& {log}"

## DON'T CHANGE ABOVE THIS LINE ##
