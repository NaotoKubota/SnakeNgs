
'''
Snakefile for gene count quantification from single-cell/nucleus RNA-seq data by STARsolo

Usage:
    snakemake -s STARsolo.smk --configfile <path to config.yaml> --cores <int> --use-singularity
'''

import sys
import pandas as pd
workdir: config["workdir"]
star_index = config["star_index"]
barcode_whitelist = config["barcode_whitelist"]
soloUMIlen = config["soloUMIlen"]
soloCellFilter = config["soloCellFilter"]
clipAdapterType = config["clipAdapterType"]
outFilterScoreMin = config["outFilterScoreMin"]
soloCBmatchWLtype = config["soloCBmatchWLtype"]
soloUMIfiltering = config["soloUMIfiltering"]
soloUMIdedup = config["soloUMIdedup"]
soloBarcodeReadLength = config["soloBarcodeReadLength"]

# experiment_table
'''
sample	R1	R2
sample1	sample1_L1_R1,sample1_L2_R1	sample1_L1_R2,sample1_L2_R2
sample2	sample2_L1_R1,sample2_L2_R1	sample2_L1_R2,sample2_L2_R2
'''
# Make string of list of R1 and R2 files for each sample
# ex. samples_dict[sample1]["R1"] = sample1_L1_R1,sample1_L2_R1
# ex. samples_dict[sample1]["R2"] = sample1_L1_R2,sample1_L2_R2
def read_experiment_table(file):
	df = pd.read_csv(file, sep="\t")
	samples_values = df["sample"].values
	R1_values = df["R1"].values
	R2_values = df["R2"].values
	samples_dict = {}
	for i in range(len(samples_values)):
		sample = str(samples_values[i])
		R1 = R1_values[i]
		R1_split = R1.split(",")
		R2 = R2_values[i]
		R2_split = R2.split(",")
		samples_dict[sample] = {"R1": R1, "R1_split": R1_split, "R2": R2, "R2_split": R2_split}
	return samples_dict
samples_dict = read_experiment_table(config["experiment_table"])
samples = [str(x) for x in list(samples_dict.keys())]

include: "common/containers.smk"

wildcard_constraints:
	sample = "|".join([re.escape(str(x)) for x in samples])

rule all:
	input:
		bam = expand("star/{sample}/Aligned.sortedByCoord.out.bam", sample = samples),
		bai = expand("star/{sample}/Aligned.sortedByCoord.out.bam.bai", sample = samples)

rule STARsolo:
	container:
		CONTAINERS["star"]
	input:
		R1_split = lambda wildcards: samples_dict[wildcards.sample]["R1_split"],
		R2_split = lambda wildcards: samples_dict[wildcards.sample]["R2_split"]
	output:
		bam = "star/{sample}/Aligned.sortedByCoord.out.bam"
	params:
		R1 = lambda wildcards: samples_dict[wildcards.sample]["R1"],
		R2 = lambda wildcards: samples_dict[wildcards.sample]["R2"]
	threads:
		workflow.cores
	benchmark:
		"benchmark/STARsolo_{sample}.txt"
	log:
		"log/STARsolo_{sample}.log"
	shell:
		"""
		ulimit -n 10000 && \
		mkdir -p star/{wildcards.sample} && \
		cd star/{wildcards.sample} && \
		STAR \
		--runThreadN {threads} \
		--soloType CB_UMI_Simple \
		--genomeDir {star_index} \
		--readFilesIn \
		{params.R2} \
		{params.R1} \
		--readFilesCommand zcat \
		--soloCBwhitelist {barcode_whitelist} \
		--soloUMIlen {soloUMIlen} \
		--soloBarcodeReadLength {soloBarcodeReadLength} \
		--soloCellFilter {soloCellFilter} \
		--soloFeatures Gene GeneFull SJ Velocyto \
		--clipAdapterType {clipAdapterType} \
		--outFilterScoreMin {outFilterScoreMin} \
		--soloCBmatchWLtype {soloCBmatchWLtype} \
		--soloUMIfiltering {soloUMIfiltering} \
		--soloUMIdedup {soloUMIdedup} \
		--outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
		--outSAMtype BAM SortedByCoordinate \
		&> ../../{log}
		"""

rule index:
	container:
		CONTAINERS["samtools"]
	input:
		bam = "star/{sample}/Aligned.sortedByCoord.out.bam"
	output:
		bai = "star/{sample}/Aligned.sortedByCoord.out.bam.bai"
	threads:
		workflow.cores / 4
	benchmark:
		"benchmark/samtools_index_{sample}.txt"
	log:
		"log/samtools_index_{sample}.log"
	shell:
		"""
		samtools index {input.bam} &> {log}
		"""
