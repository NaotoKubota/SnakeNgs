
'''
Snakefile for ATAC-seq callpeak

Usage:
    snakemake -s callpeak_ATACseq.smk --configfile <path to config.yaml> --cores <int> --use-singularity
'''

## DON'T CHANGE BELOW THIS LINE ##

workdir: config["workdir"]
samples = config["samples"]

include: "common/containers.smk"

wildcard_constraints:
    sample = "|".join([re.escape(x) for x in samples])

rule all:
    input:
        multiqc = "multiqc_callpeak/multiqc_report.html",
        xls = expand("macs2/{sample}/{sample}_peaks.xls", sample = samples)

rule macs2:
    container:
        CONTAINERS["macs2"]
    input:
        bam = "bowtie2/{sample}.sort.rmdup.bam"
    output:
        outdir = directory("macs2/{sample}"),
        xls = "macs2/{sample}/{sample}_peaks.xls",
        bdg = "macs2/{sample}/{sample}_treat_pileup.bdg"
    params:
        genomesize = config["genomesize"],
        fileformat = config["fileformat"],
        qvalue = config["qvalue"],
        broad = config["broad"],
        broad_cutoff = config["broad_cutoff"]
    benchmark:
        "benchmark/macs2/{sample}.txt"
    log:
        "log/macs2/{sample}.log"
    shell:
        "if [ '{params.broad}' == 'True' ]; then "
        "macs2 callpeak "
        "-t {input.bam} "
        "-g {params.genomesize} -f {params.fileformat} "
        "-q {params.qvalue} "
        "--broad --broad-cutoff {params.broad_cutoff} "
        "-B -n {wildcards.sample} --outdir {output.outdir} >& {log}; "
        "else "
        "macs2 callpeak "
        "-t {input.bam} "
        "-g {params.genomesize} -f {params.fileformat} "
        "-q {params.qvalue} "
        "-B -n {wildcards.sample} --outdir {output.outdir} >& {log}; "
        "fi"

rule multiqc:
    container:
        CONTAINERS["multiqc"]
    input:
        macs2log = expand("macs2/{sample}/{sample}_peaks.xls", sample = samples)
    output:
        "multiqc_callpeak/multiqc_report.html"
    benchmark:
        "benchmark/multiqc_callpeak.txt"
    log:
        "log/multiqc_callpeak.log"
    shell:
        "rm -rf multiqc_callpeak && "
        "mkdir -p multiqc_callpeak/log && "
        "cp {input.macs2log} multiqc_callpeak/log && "
        "multiqc -o multiqc_callpeak/ multiqc_callpeak/log >& {log} && "
        "rm -rf multiqc_callpeak/log"

## DON'T CHANGE ABOVE THIS LINE ##
