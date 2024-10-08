
'''
Snakefile for iCLIP-seq preprocessing

Usage:
    snakemake -s preprocessing_iCLIPseq.smk --configfile <path to config.yaml> --cores <int> --use-singularity
'''

## DON'T CHANGE BELOW THIS LINE ##

workdir: config["workdir"]
star_index = config["star_index"]
samples = config["samples"]

rule all:
    input:
        multiqc = "multiqc_preprocessing/multiqc_report.html",
        bam = expand("star/{sample}/{sample}_Aligned.out.bam", sample = samples),
        bigwig = expand("bigwig/{sample}.bw", sample = samples),

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
		"fastp -i {input.R1} "
		"-o {output.R1} -w {threads} "
		"--trim_front1 {config[trim_front1]} --max_len1 {config[max_len1]} --length_required {config[length_required]} "
		"-h {params.html} -j {output.json} >& {log}"

rule mapping:
	wildcard_constraints:
		sample = "|".join([re.escape(x) for x in samples])
	container:
		"docker://quay.io/biocontainers/star:2.7.11a--h0033a41_0"
	input:
		R1 = "fastp/{sample}.fastq.gz"
	output:
		sam = "star/{sample}/{sample}_Aligned.out.sam",
		log = "star/{sample}/{sample}_Log.final.out"
	params:
		outdir = "star/{sample}/{sample}_"
	threads:
		workflow.cores
	benchmark:
		"benchmark/star_{sample}.txt"
	log:
		"log/star_{sample}.log"
	shell:
		"STAR --runThreadN {threads} --genomeDir {star_index} "
		"--outFilterMultimapNmax 1 "
		"--readFilesIn {input.R1} --readFilesCommand zcat --outFileNamePrefix {params.outdir} >& {log}"

rule sort:
	wildcard_constraints:
		sample = "|".join([re.escape(x) for x in samples])
	container:
		"docker://quay.io/biocontainers/samtools:1.18--h50ea8bc_1"
	input:
		"star/{sample}/{sample}_Aligned.out.sam"
	output:
		"star/{sample}/{sample}_Aligned.out.bam"
	threads:
		workflow.cores / 4
	benchmark:
		"benchmark/samtools_{sample}.txt"
	log:
		"log/samtools_{sample}.log"
	shell:
		"samtools sort -@ {threads} -O bam -o {output} {input} >& {log} && "
		"samtools index {output} && "
		"rm -rf {input}"

rule markdup:
	wildcard_constraints:
		sample = "|".join([re.escape(x) for x in samples])
	container:
		"docker://quay.io/biocontainers/picard:3.1.1--hdfd78af_0"
	input:
		"star/{sample}/{sample}_Aligned.out.bam"
	output:
		"star/{sample}/{sample}.sort.rmdup.bam"
	threads:
		8
	benchmark:
		"benchmark/picard_{sample}.txt"
	log:
		"log/picard_{sample}.log"
	shell:
		"picard MarkDuplicates "
		"-I {input} -O {output} "
		"-M {log} "
		"--REMOVE_DUPLICATES true "
		"--TMP_DIR star/tmp && "
		"rm -rf star/*_Aligned.bam.bai"

rule index:
	wildcard_constraints:
		sample = "|".join([re.escape(x) for x in samples])
	container:
		"docker://quay.io/biocontainers/samtools:1.18--h50ea8bc_1"
	input:
		"star/{sample}/{sample}.sort.rmdup.bam"
	output:
		"star/{sample}/{sample}.sort.rmdup.bam.bai"
	shell:
		"samtools index {input}"

rule bigwig:
	wildcard_constraints:
		sample = "|".join([re.escape(x) for x in samples])
	container:
		"docker://quay.io/biocontainers/deeptools:3.5.4--pyhdfd78af_1"
	input:
		bam = "star/{sample}/{sample}.sort.rmdup.bam",
		bai = "star/{sample}/{sample}.sort.rmdup.bam.bai"
	output:
		"bigwig/{sample}.bw"
	threads:
		8
	benchmark:
		"benchmark/bamCoverage_{sample}.txt"
	log:
		"log/bamCoverage_{sample}.log"
	shell:
		"bamCoverage -b {input.bam} -o {output} -p {threads} --binSize 1 --normalizeUsing CPM >& {log}"

rule multiqc:
	container:
		"docker://multiqc/multiqc:v1.25"
	input:
		json = expand("fastp/log/{sample}.json", sample = samples),
		starlog = expand("star/{sample}/{sample}_Log.final.out", sample = samples),
		picardlog = expand("log/picard_{sample}.log", sample = samples)
	output:
		"multiqc_preprocessing/multiqc_report.html"
	benchmark:
		"benchmark/multiqc.txt"
	log:
		"log/multiqc.log"
	shell:
		"rm -rf multiqc_preprocessing && "
		"mkdir -p multiqc_preprocessing/log && "
		"cp {input.json} {input.starlog} {input.picardlog} multiqc_preprocessing/log && "
		"multiqc -o multiqc_preprocessing/ multiqc_preprocessing/log >& {log} && "
		"rm -rf multiqc_preprocessing/log"
