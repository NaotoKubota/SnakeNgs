
'''
Snakefile for bam2cram conversion

Usage:
    snakemake -s bam2cram.smk --config workdir=/path/to/workdir reference=/path/to/reference.fa bam_list=/path/to/bam_list.txt --cores <int> --use-singularity
'''

def get_samples(bam_list):
    samples = []
    with open(bam_list, "r") as f:
        for line in f:
            samples.append(line.strip().replace(".bam", ""))
    return samples

samples = get_samples(config["bam_list"])

rule all:
    input:
        cram = expand("{sample}.cram", sample = samples),
        cram_index = expand("{sample}.cram.crai", sample = samples)

rule bam2cram:
    wildcard_constraints:
        sample = "|".join([re.escape(x) for x in samples])
    container:
        "docker://quay.io/biocontainers/samtools:1.18--h50ea8bc_1"
    input:
        bam = "{sample}.bam"
    output:
        cram = "{sample}.cram"
    threads:
        workflow.cores / 2
    shell:
        """
        samtools view -@ {threads} -T {config[reference]} -C -o {output.cram} {input.bam} && samtools index {output.cram}
        rm -rf {input.bam}
        """

rule cram_index:
    wildcard_constraints:
        sample = "|".join([re.escape(x) for x in samples])
    container:
        "docker://quay.io/biocontainers/samtools:1.18--h50ea8bc_1"
    input:
        cram = "{sample}.cram"
    output:
        cram_index = "{sample}.cram.crai"
    threads:
        workflow.cores / 4
    shell:
        "samtools index -@ {threads} {input.cram}"
