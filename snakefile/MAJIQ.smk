
'''
Snakefile for MAJIQ

Usage:
    snakemake -s MAJIQ.smk --configfile config.yaml --cores <int> --use-singularity
'''

import os
from pathlib import Path
import pandas as pd
import configparser
workdir: config["workdir"]

'''experiment_table.tsv
sample	bam	group
sample_01	/path/to/sample_01/sample_01.bam	Ref
sample_02	/path/to/sample_02/sample_02.bam	Ref
sample_03	/path/to/sample_03/sample_03.bam	Ref
sample_04	/path/to/sample_04/sample_04.bam	Alt
sample_05	/path/to/sample_05/sample_05.bam	Alt
sample_06	/path/to/sample_06/sample_06.bam	Alt
'''

'''setting_file.ini
[info]
bamdirs=/path/to/sample_01,/path/to/sample_02,/path/to/sample_03,/path/to/sample_04,/path/to/sample_05,/path/to/sample_06
genome=hg19
strandness=None
[experiments]
Ref=sample_01,sample_02,sample_03
Alt=sample_04,sample_05,sample_06
'''

# Convert experiment_table.tsv to setting_file.ini
def convert_experiment_table_to_setting_file(experiment_table, genome, strandness):
    '''
    Convert experiment_table.tsv to setting_file.ini
    '''
    df = pd.read_csv(experiment_table, sep = '\t')
    bamdirs = ",".join(list(set([os.path.dirname(x) for x in df['bam'].tolist()])))
    Ref_bam = ",".join([os.path.splitext(os.path.basename(x))[0] for x in df[df['group'] == 'Ref']['bam'].tolist()])
    Alt_bam = ",".join([os.path.splitext(os.path.basename(x))[0] for x in df[df['group'] == 'Alt']['bam'].tolist()])
    with open("setting_file.ini", "w") as f:
        f.write(f"[info]\n")
        f.write(f"bamdirs={bamdirs}\n")
        f.write(f"genome={genome}\n")
        f.write(f"strandness={strandness}\n")
        f.write(f"[experiments]\n")
        f.write(f"Ref={Ref_bam}\n")
        f.write(f"Alt={Alt_bam}\n")

# Get samples per group
def get_samples_per_group(experiment_table):
    '''
    Get samples per group from experiment_table.tsv
    '''
    df = pd.read_csv(experiment_table, sep = '\t')
    samples_Ref = [os.path.splitext(os.path.basename(x))[0] for x in df[df['group'] == 'Ref']['bam'].tolist()]
    samples_Alt = [os.path.splitext(os.path.basename(x))[0] for x in df[df['group'] == 'Alt']['bam'].tolist()]
    return samples_Ref, samples_Alt
samples_Ref, samples_Alt = get_samples_per_group(config["experiment_table"])

rule all:
    input:
        voila_tsv = expand("{module}/Alt-Ref.{module}.tsv", module = ["het", "deltapsi"]),
        summary = expand("voila_modulize/{module}/summary.tsv", module = ["het", "deltapsi"])

rule convert_experiment_table_to_setting_file:
    output:
        setting_file = "setting_file.ini",
    threads:
        1
    benchmark:
        "benchmark/convert_experiment_table_to_setting_file.txt"
    log:
        "log/convert_experiment_table_to_setting_file.log"
    run:
        convert_experiment_table_to_setting_file(config["experiment_table"], config["genome"], config["strandness"])

rule build:
    container:
        config["container"]
    input:
        setting_file = "setting_file.ini",
        gff = config["gff"]
    output:
        splicegraph = "build/splicegraph.sql",
        majiq_ref = expand("build/{sample_ref}.majiq", sample_ref = samples_Ref),
        majiq_alt = expand("build/{sample_alt}.majiq", sample_alt = samples_Alt)
    threads:
        workflow.cores
    benchmark:
        "benchmark/build.txt"
    log:
        "log/build.log"
    shell:
        """
        majiq build \
        -c {input.setting_file} \
        {input.gff} \
        -o build \
        -j {threads} \
        > {log} 2>&1
        """

rule heterogen:
    container:
        config["container"]
    input:
        majiq_ref = expand("build/{sample_ref}.majiq", sample_ref = samples_Ref),
        majiq_alt = expand("build/{sample_alt}.majiq", sample_alt = samples_Alt)
    output:
        voila = "het/Alt-Ref.het.voila"
    params:
        majiq_ref = " ".join(["build/" + sample + ".majiq" for sample in samples_Ref]),
        majiq_alt = " ".join(["build/" + sample + ".majiq" for sample in samples_Alt])
    threads:
        workflow.cores / 2
    benchmark:
        "benchmark/heterogen.txt"
    log:
        "log/heterogen.log"
    shell:
        """
        majiq heterogen \
        -j {threads} \
        -o het \
        -n Alt Ref \
        -grp1 {params.majiq_alt} \
        -grp2 {params.majiq_ref} \
        > {log} 2>&1
        """

rule deltapsi:
    container:
        config["container"]
    input:
        majiq_ref = expand("build/{sample_ref}.majiq", sample_ref = samples_Ref),
        majiq_alt = expand("build/{sample_alt}.majiq", sample_alt = samples_Alt)
    output:
        voila = "deltapsi/Alt-Ref.deltapsi.voila"
    params:
        majiq_ref = " ".join(["build/" + sample + ".majiq" for sample in samples_Ref]),
        majiq_alt = " ".join(["build/" + sample + ".majiq" for sample in samples_Alt])
    threads:
        workflow.cores / 2
    benchmark:
        "benchmark/deltapsi.txt"
    log:
        "log/deltapsi.log"
    shell:
        """
        majiq deltapsi \
        -j {threads} \
        -o deltapsi \
        -n Alt Ref \
        -grp1 {params.majiq_alt} \
        -grp2 {params.majiq_ref} \
        > {log} 2>&1
        """

rule voila_tsv:
    container:
        config["container"]
    wildcard_constraints:
        module = "|".join(["het", "deltapsi"])
    input:
        voila = "{module}/Alt-Ref.{module}.voila",
        splicegraph = "build/splicegraph.sql"
    output:
        tsv = "{module}/Alt-Ref.{module}.tsv"
    threads:
        workflow.cores / 2
    benchmark:
        "benchmark/voila_tsv_{module}.txt"
    log:
        "log/voila_tsv_{module}.log"
    shell:
        """
        voila tsv \
        -j {threads} \
        -f {output.tsv} \
        {input.splicegraph} \
        {input.voila} \
        > {log} 2>&1
        """

rule voila_modulize:
    container:
        config["container"]
    wildcard_constraints:
        module = "|".join(["het", "deltapsi"])
    input:
        splicegraph = "build/splicegraph.sql",
        voila = "{module}/Alt-Ref.{module}.voila"
    output:
        summary = "voila_modulize/{module}/summary.tsv"
    params:
        outdir = "voila_modulize/{module}"
    threads:
        workflow.cores / 2
    benchmark:
        "benchmark/voila_modulize_{module}.txt"
    log:
        "log/voila_modulize_{module}.log"
    shell:
        """
        voila modulize \
        -j {threads} \
        -d {params.outdir} \
        --show-all \
        --overwrite \
        {input.splicegraph} \
        {input.voila} \
        > {log} 2>&1
        """
