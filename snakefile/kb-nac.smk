
'''
Snakefile for gene count quantification from single-nuclei RNA-seq data by kb

Usage:
    snakemake -s kb-nac.smk --configfile <path to config.yaml> --cores <int> --use-singularity
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
		R1 = R1_values[i].split(",")
		R2 = R2_values[i].split(",")
		R1_R2 = ""
		for j in range(len(R1)):
			R1_R2 += R1[j] + " " + R2[j] + " "
		samples_dict[sample] = R1_R2.strip()
	return samples_dict
samples_dict = read_experiment_table(config["experiment_table"])
samples = [str(x) for x in list(samples_dict.keys())]

rule all:
	input:
		adata_unfiliterd_h5ad = expand("kb/{sample}/counts_unfiltered/adata.h5ad", sample = samples),
		adata_filtered_h5ad = expand("kb/{sample}/counts_filtered/adata.h5ad", sample = samples),
		inspect = expand("kb/{sample}/inspect.json", sample = samples)

rule kb_ref:
	container:
		"docker://quay.io/biocontainers/kb-python:0.28.2--pyhdfd78af_2"
	input:
		dna_fasta = config["dna_fasta"],
		gtf = config["gtf"]
	output:
		kb_index = "kb_index/index.idx",
		kb_t2g = "kb_index/t2g.txt",
		kb_cdna_t2c = "kb_index/cdna_t2c.txt",
		kb_intron_t2c = "kb_index/intron_t2c.txt"
	threads:
		workflow.cores
	benchmark:
		"benchmark/kb_ref.txt"
	log:
		"log/kb_ref.log"
	shell:
		"""
		kb ref \
		-i {output.kb_index} \
		-g {output.kb_t2g} \
		-f1 kb_index/cdna.fa \
		-f2 kb_index/intron.fa \
		-c1 {output.kb_cdna_t2c} \
		-c2 {output.kb_intron_t2c} \
		--workflow nac \
		-t {threads} \
		--verbose \
		{input.dna_fasta} {input.gtf} >& {log}
		"""

rule kb_count:
	wildcard_constraints:
		sample = "|".join([re.escape(str(x)) for x in samples])
	container:
		"docker://quay.io/biocontainers/kb-python:0.28.2--pyhdfd78af_2"
	input:
		kb_index = "kb_index/index.idx",
		kb_t2g = "kb_index/t2g.txt",
		kb_cdna_t2c = "kb_index/cdna_t2c.txt",
		kb_intron_t2c = "kb_index/intron_t2c.txt"
	output:
		adata_unfiltered_h5ad = "kb/{sample}/counts_unfiltered/adata.h5ad",
		adata_filtered_h5ad = "kb/{sample}/counts_filtered/adata.h5ad",
		inspect = "kb/{sample}/inspect.json"
	params:
		R1_R2 = lambda wildcards: samples_dict[wildcards.sample],
		kb_output = "kb/{sample}"
	threads:
		workflow.cores / 4
	benchmark:
		"benchmark/kb_count_{sample}.txt"
	log:
		"log/kb_count_{sample}.log"
	shell:
		"""
		kb count \
		-i {input.kb_index} \
		-g {input.kb_t2g} \
		-c1 {input.kb_cdna_t2c} \
		-c2 {input.kb_intron_t2c} \
		-x {config[technology]} \
		-o {params.kb_output} \
		-t {threads} \
		--workflow nac \
		--h5ad \
		--gene-names \
		--sum total \
		--filter bustools \
		--overwrite \
		--verbose \
		{params.R1_R2} >& {log}
		"""
