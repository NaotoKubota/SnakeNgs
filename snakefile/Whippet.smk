
'''
Snakefile for Whippet

Usage:
    snakemake -s Whippet.smk --configfile config.yaml --cores <int> --use-singularity
'''

import os
import pandas as pd
workdir: config['workdir']
absolute_path_to_workdir = os.path.abspath(config['workdir'])

'''experiment_table
sample	fastq	group
sample_01	/path/to/sample_01_1.fastq.gz,/path/to/sample_01_2.fastq.gz	Ref
sample_02	/path/to/sample_02_1.fastq.gz,/path/to/sample_02_2.fastq.gz	Ref
sample_03	/path/to/sample_03_1.fastq.gz,/path/to/sample_03_2.fastq.gz	Ref
sample_04	/path/to/sample_04_1.fastq.gz,/path/to/sample_04_2.fastq.gz	Alt
sample_05	/path/to/sample_05_1.fastq.gz,/path/to/sample_05_2.fastq.gz	Alt
sample_06	/path/to/sample_06_1.fastq.gz,/path/to/sample_06_2.fastq.gz	Alt
'''

# Read in the experiment table
def read_experiment_table(experiment_table):
    df = pd.read_csv(experiment_table, sep = "\t")
    samples = df['sample'].tolist()
    samples_Ref = df[df['group'] == 'Ref']['sample'].tolist()
    samples_Alt = df[df['group'] == 'Alt']['sample'].tolist()
    # Make dictionary with sample name as key and fastq file as value
    '''
    {
        sample1: {
            R1: /path/to/sample1_R1.fastq.gz,
            R2: /path/to/sample1_R2.fastq.gz
        },
        sample2: {
            R1: /path/to/sample2_R1.fastq.gz,
            R2: /path/to/sample2_R2.fastq.gz
        }
    }
    '''
    sample_fastq_dict = {sample: {'R1': os.path.join(absolute_path_to_workdir, fastq.split(',')[0]), 'R2': os.path.join(absolute_path_to_workdir, fastq.split(',')[1])} for sample, fastq in zip(df['sample'], df['fastq'])}
    return samples, samples_Ref, samples_Alt, sample_fastq_dict

samples, samples_Ref, samples_Alt, sample_fastq_dict = read_experiment_table(config['experiment_table'])

rule all:
    input:
        psi = expand("quant/all/{sample}.psi.gz", sample = samples),
        diff = "delta/group.diff.gz"

rule quant:
    wildcard_constraints:
        sample = "|".join(samples)
    container:
        "docker://cloxd/whippet:1.6.1"
    input:
        R1 = lambda wildcards: sample_fastq_dict[wildcards.sample]['R1'],
        R2 = lambda wildcards: sample_fastq_dict[wildcards.sample]['R2'],
        whippet_index = config['whippet_index']
    output:
        psi = "quant/all/{sample}.psi.gz"
    threads:
        workflow.cores
    benchmark:
        "benchmark/whippet_quant_{sample}.txt"
    log:
        "log/whippet_quant_{sample}.log"
    shell:
        """
        mkdir -p quant/{wildcards.sample}/ quant/all/ && \
        julia /whippet/bin/whippet-quant.jl \
		{input.R1} \
		{input.R2} \
		-x {input.whippet_index} \
		-o quant/{wildcards.sample}/{wildcards.sample} \
        >& {log} && \
        cp quant/{wildcards.sample}/{wildcards.sample}.psi.gz quant/all/
        """

rule delta:
    container:
        "docker://cloxd/whippet:1.6.1"
    input:
        psi_all = expand("quant/all/{sample}.psi.gz", sample = samples)
    output:
        diff = "delta/group.diff.gz"
    params:
        psi_Ref = ",".join([x + ".psi.gz" for x in samples_Ref]),
        psi_Alt = ",".join([x + ".psi.gz" for x in samples_Alt])
    threads:
        1
    benchmark:
        "benchmark/whippet_delta.txt"
    log:
        "log/whippet_delta.log"
    shell:
        """
        mkdir -p delta/ && \
        cd quant/all && \
        julia /whippet/bin/whippet-delta.jl \
        -a {params.psi_Alt}, \
        -b {params.psi_Ref}, \
        -o ../../delta/group \
        -r 10 \
        >& ../../{log}
        """

