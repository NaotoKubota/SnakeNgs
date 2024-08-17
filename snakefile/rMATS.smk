
'''
Snakefile for rMATS

Usage:
    snakemake -s rMATS.smk --configfile config.yaml --cores <int> --use-singularity --singularity-args "--bind $HOME:$HOME" --rerun-incomplete
'''

import os
import pandas as pd
workdir: config['workdir']

'''experiment_table
sample	bam	group
sample_01	/path/to/sample_01.bam	Ref
sample_02	/path/to/sample_02.bam	Ref
sample_03	/path/to/sample_03.bam	Ref
sample_04	/path/to/sample_04.bam	Alt
sample_05	/path/to/sample_05.bam	Alt
sample_06	/path/to/sample_06.bam	Alt
'''

def get_bam_list(experiment_table):
    df = pd.read_csv(experiment_table, sep = '\t')
    Ref_bam = df[df['group'] == 'Ref']['bam'].tolist()
    Alt_bam = df[df['group'] == 'Alt']['bam'].tolist()
    return Ref_bam, Alt_bam

Ref_bam, Alt_bam = get_bam_list(config['experiment_table'])

rule all:
    input:
        summary = "results/summary.txt"

rule makeRmatsInput:
    output:
        bam_Ref_list = "input/bam_Ref.txt",
        bam_Alt_list = "input/bam_Alt.txt"
    params:
        Ref_bam = ",".join(Ref_bam),
        Alt_bam = ",".join(Alt_bam)
    threads:
        1
    benchmark:
        "benchmark/makeRmatsInput.txt"
    shell:
        """
        echo {params.Ref_bam} > {output.bam_Ref_list}
        echo {params.Alt_bam} > {output.bam_Alt_list}
        """

rule rmats:
    container:
        "docker://xinglab/rmats:v4.3.0"
    input:
        bam_Ref_list = "input/bam_Ref.txt",
        bam_Alt_list = "input/bam_Alt.txt"
    output:
        summary = "results/summary.txt"
    threads:
        workflow.cores
    benchmark:
        "benchmark/rmats.txt"
    log:
        "log/rmats.log"
    shell:
        """
        python /rmats/rmats.py \
        --b1 {input.bam_Alt_list} \
        --b2 {input.bam_Ref_list} \
        --gtf {config[gtf]} \
        -t {config[layout]} \
        --readLength {config[readLength]} \
        --novelSS \
        --cstat {config[cstat]} \
        --anchorLength {config[anchorLength]} \
        --mil {config[mil]} \
        --mel {config[mel]} \
        --nthread {threads} \
        --od results/ \
        --tmp results/tmp \
        >& {log}
        """

