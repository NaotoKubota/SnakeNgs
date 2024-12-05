
'''
Snakefile for Small RNA-seq processing

Usage:
    snakemake -s SmallRNAseq.smk --configfile <path to config.yaml> --cores <int> --use-singularity
'''

workdir: config["workdir"]
bowtie2_index = config["bowtie2_index"]
samples = config["samples"]
gtf = config["gtf"]

rule all:
	input:
		bigwig = expand("bigwig/{sample}.bw", sample = samples),
		bigbed = expand("bigbed/{sample}.bb", sample = samples),
		featureCounts_results = expand("quantification/{biotype}.txt", biotype = ["miRNA", "snoRNA", "snRNA", "rRNA"]),
		multiqc = "multiqc/multiqc_report.html"

rule qc:
	wildcard_constraints:
		sample = "|".join([re.escape(x) for x in samples])
	container:
		"docker://quay.io/biocontainers/fastp:0.23.4--hadf994f_2"
	input:
		R1 = "fastq/{sample}.fastq.gz"
	output:
		R1 = "fastp/{sample}.fastq.gz",
		json = "fastp/log/{sample}.json"
	params:
		html = "fastp/log/{sample}.html"
	threads:
		8
	benchmark:
		"benchmark/fastp_{sample}.txt"
	log:
		"log/fastp_{sample}.log"
	shell:
		'''
		fastp -i {input.R1} \
		-o {output.R1} -w {threads} \
		--trim_front1 {config[trim_front1]} --trim_tail1 {config[trim_tail1]} \
		--length_required {config[length_required]} --qualified_quality_phred {config[qualified_quality_phred]} \
		-h {params.html} -j {output.json} >& {log}
		'''

rule mapping:
	wildcard_constraints:
		sample = "|".join([re.escape(x) for x in samples])
	container:
		"docker://quay.io/biocontainers/bowtie2:2.5.2--py39h6fed5c7_0"
	input:
		R1 = "fastp/{sample}.fastq.gz"
	output:
		sam = temp("bowtie2/{sample}.sam")
	threads:
		workflow.cores
	benchmark:
		"benchmark/bowtie2_{sample}.txt"
	log:
		"log/bowtie2/{sample}.log"
	shell:
		'''
		bowtie2 -p {threads} -x {bowtie2_index} \
		--very-sensitive-local \
		-U {input.R1} -S {output.sam} >& {log}
		'''

rule sort:
	wildcard_constraints:
		sample = "|".join([re.escape(x) for x in samples])
	container:
		"docker://quay.io/biocontainers/samtools:1.18--h50ea8bc_1"
	input:
		"bowtie2/{sample}.sam"
	output:
		"bowtie2/{sample}.sort.bam"
	threads:
		8
	benchmark:
		"benchmark/samtools_{sample}.txt"
	log:
		"log/samtools_{sample}.log"
	shell:
		'''
		samtools sort -@ {threads} -O bam -o {output} {input} >& {log}
		'''

rule index:
	wildcard_constraints:
		sample = "|".join([re.escape(x) for x in samples])
	container:
		"docker://quay.io/biocontainers/samtools:1.18--h50ea8bc_1"
	input:
		"bowtie2/{sample}.sort.bam"
	output:
		"bowtie2/{sample}.sort.bam.bai"
	shell:
		'''
		samtools index {input}
		'''

rule bigwig:
	wildcard_constraints:
		sample = "|".join([re.escape(x) for x in samples])
	container:
		"docker://quay.io/biocontainers/deeptools:3.5.4--pyhdfd78af_1"
	input:
		bam = "bowtie2/{sample}.sort.bam",
		bai = "bowtie2/{sample}.sort.bam.bai"
	output:
		"bigwig/{sample}.bw"
	threads:
		8
	benchmark:
		"benchmark/bamCoverage_{sample}.txt"
	log:
		"log/bamCoverage_{sample}.log"
	shell:
		'''
		bamCoverage -b {input.bam} -o {output} -p {threads} --binSize 1 --normalizeUsing CPM >& {log}
		'''

rule bed:
	wildcard_constraints:
		sample = "|".join([re.escape(x) for x in samples])
	container:
		"docker://quay.io/biocontainers/bedtools:2.24--1"
	input:
		bam = "bowtie2/{sample}.sort.bam",
		bai = "bowtie2/{sample}.sort.bam.bai"
	output:
		"bed/{sample}.bed"
	threads:
		1
	benchmark:
		"benchmark/bamToBed_{sample}.txt"
	log:
		"log/bamToBed_{sample}.log"
	shell:
		'''
		bamToBed -i {input.bam} > {output} 2> {log}
		'''

rule bedtobigbed:
	wildcard_constraints:
		sample = "|".join([re.escape(x) for x in samples])
	container:
		"docker://quay.io/biocontainers/ucsc-bedtobigbed:473--he8037a5_0"
	input:
		bed = "bed/{sample}.bed"
	output:
		"bigbed/{sample}.bb"
	params:
		chrom_sizes = config["chrom_sizes"]
	benchmark:
		"benchmark/bedToBigBed_{sample}.txt"
	log:
		"log/bedToBigBed_{sample}.log"
	shell:
		'''
		bedToBigBed {input.bed} {params.chrom_sizes} {output} >& {log}
		'''

rule create_gtf_subset:
	wildcard_constraints:
		biotype = "|".join([re.escape(x) for x in ["miRNA", "snoRNA", "snRNA", "rRNA"]])
	input:
		gtf = gtf
	output:
		"gtf/{biotype}.gtf"
	benchmark:
		"benchmark/create_gtf_subset_{biotype}.txt"
	log:
		"log/create_gtf_subset_{biotype}.log"
	shell:
		'''
		grep -e 'gene_biotype \"{wildcards.biotype}\"' {input.gtf} > {output}
		'''

rule quantification:
	container:
		"docker://quay.io/biocontainers/subread:2.0.6--he4a0461_0"
	wildcard_constraints:
		biotype = "|".join([re.escape(x) for x in ["miRNA", "snoRNA", "snRNA", "rRNA"]])
	input:
		bam = expand("bowtie2/{sample}.sort.bam", sample = samples),
		gtf = "gtf/{biotype}.gtf"
	output:
		results = "quantification/{biotype}.txt",
		summary = "quantification/{biotype}.txt.summary"
	threads:
		workflow.cores / 2
	benchmark:
		"benchmark/featureCounts_{biotype}.txt"
	log:
		"log/featureCounts_{biotype}.log"
	shell:
		'''
		featureCounts \
		-t exon \
		-g gene_id \
		-O \
		-s 1 \
		-T {threads} \
		-a {input.gtf} \
		-o {output.results} \
		{input.bam} \
		>& {log}
		'''

rule multiqc:
	container:
		"docker://multiqc/multiqc:v1.25"
	input:
		json = expand("fastp/log/{sample}.json", sample = samples),
		bowtie2log = expand("log/bowtie2/{sample}.log", sample = samples)
	output:
		"multiqc/multiqc_report.html"
	benchmark:
		"benchmark/multiqc.txt"
	log:
		"log/multiqc.log"
	shell:
		'''
		rm -rf multiqc && \
		mkdir -p multiqc/log && \
		cp {input.json} {input.bowtie2log} multiqc/log && \
		multiqc -o multiqc/ multiqc/log >& {log} && \
		rm -rf multiqc/log
		'''
