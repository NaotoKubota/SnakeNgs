
'''
Snakefile for gene count quantification from single-cell/nucleus RNA-seq data by cellranger

Usage:
    snakemake -s count.smk --configfile <path to config.yaml> --cores <int> --use-singularity
'''

import sys
import pandas as pd
workdir: config["workdir"]

# Read experiment table to get sample names and fastq files
# Make string of list of R1 and R2 files for each sample
# ex. sample1_L1_R1 sample1_L1_R2 sample1_L2_R1 sample1_L2_R2...
def read_experiment_table(file):
    df = pd.read_csv(file, sep="\t")
    samples_values = df["sample"].values
    R1_values = df["R1"].values
    R2_values = df["R2"].values
    samples_dict = {}
    for i in range(len(samples_values)):
        sample = str(samples_values[i])
        samples_dict[sample] = {}
        R1 = R1_values[i].split(",")
        R2 = R2_values[i].split(",")
        samples_dict[sample]["R1"] = R1
        samples_dict[sample]["R2"] = R2
    return samples_dict

def fastq_for_cellranger(samples_dict, workdir):
    fastq_dict = {}
    for sample in samples_dict:
        fastq_dict[sample] = {}
        for read in ["R1", "R2"]:
            fastq_dict[sample][read] = {}
            for i in range(len(samples_dict[sample][read])):
                fastq_dict[sample][read][i+1] = {}
                fastq = samples_dict[sample][read][i]
                fastq_converted_name = f"{workdir}/fastq/{sample}_S1_L00{i+1}_{read}_001.fastq.gz"
                fastq_dict[sample][read][i+1]["fastq"] = fastq
                fastq_dict[sample][read][i+1]["fastq_converted_name"] = fastq_converted_name
    return fastq_dict

samples_dict = read_experiment_table(config["experiment_table"])
samples = [str(x) for x in list(samples_dict.keys())]
fastq_dict = fastq_for_cellranger(samples_dict, config["workdir"])
fastq_converted_name_all = [fastq_dict[sample][read][i+1]["fastq_converted_name"] for sample in samples for read in ["R1", "R2"] for i in range(len(fastq_dict[sample][read]))]

include: "common/containers.smk"

wildcard_constraints:
    sample = "|".join([re.escape(str(x)) for x in samples])

rule all:
    input:
        expand("count/{sample}/outs/web_summary.html", sample = samples),
        "multiqc/multiqc_report.html"

rule fastq_naming_conversion:
    input:
        config["experiment_table"]
    output:
        fastq_converted_name_all
    threads:
        1
    benchmark:
        "benchmark/fastq_naming_conversion.txt"
    log:
        "log/fastq_naming_conversion.log"
    run:
        for sample in fastq_dict:
            for read in ["R1", "R2"]:
                for i in range(len(fastq_dict[sample][read])):
                    fastq = fastq_dict[sample][read][i+1]["fastq"]
                    fastq_converted_name = fastq_dict[sample][read][i+1]["fastq_converted_name"]
                    shell(f"ln -sf {fastq} {fastq_converted_name}")

rule count:
    container:
        CONTAINERS["cellranger"]
    input:
        fastq_converted_name = lambda wildcards: [fastq_dict[wildcards.sample][read][i+1]["fastq_converted_name"] for read in ["R1", "R2"] for i in range(len(fastq_dict[wildcards.sample][read]))],
        transcriptome = config["transcriptome"]
    output:
        "count/{sample}/outs/web_summary.html"
    threads:
        workflow.cores
    resources:
        mem_mb=64
    benchmark:
        "benchmark/count_{sample}.txt"
    log:
        "log/count_{sample}.log"
    shell:
        """
        rm -rf count/{wildcards.sample} && \
        cellranger count \
        --output-dir count/{wildcards.sample} \
        --id {wildcards.sample} \
        --transcriptome {input.transcriptome} \
        --fastqs fastq \
        --sample {wildcards.sample} \
        --include-introns true \
        --create-bam {config[create_bam]} \
        --localcores {threads} \
        --localmem {resources.mem_mb} \
        >& {log}
        """

rule multiqc:
    container:
        CONTAINERS["multiqc"]
    input:
        expand("count/{sample}/outs/web_summary.html", sample=samples)
    output:
        "multiqc/multiqc_report.html"
    threads:
        1
    benchmark:
        "benchmark/multiqc.txt"
    log:
        "log/multiqc.log"
    shell:
        """
        multiqc --force -o multiqc/ count/ >& {log}
        """
