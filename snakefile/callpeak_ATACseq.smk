
'''
Snakefile for ATAC-seq callpeak

Usage:
    snakemake -s callpeak_ATACseq.smk --configfile <path to config.yaml> --cores <int> --use-singularity
'''

## DON'T CHANGE BELOW THIS LINE ##

workdir: config["workdir"]
samples = config["samples"]

rule all:
    input:
        multiqc = "multiqc_callpeak/multiqc_report.html",
        xls = expand("macs2/{sample}/{sample}_peaks.xls", sample = samples),
        bigwig = expand("macs2/{sample}/{sample}_treat_pileup.bw", sample = samples)

rule macs2:
    wildcard_constraints:
        sample = "|".join([re.escape(x) for x in samples])
    container:
        "docker://quay.io/biocontainers/macs2:2.2.9.1--py39hf95cd2a_0"
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

rule bedgraphtobigwig:
    wildcard_constraints:
        sample = "|".join(samples)
    container:
        "docker://quay.io/biocontainers/ucsc-bedgraphtobigwig:445--h954228d_0"
    input:
        bedgraph = "macs2/{sample}/{sample}_treat_pileup.bdg"
    output:
        bigwig = "macs2/{sample}/{sample}_treat_pileup.bw"
    params:
        assembly = config["assembly"]
    benchmark:
        "benchmark/bedgraphtobigwig/{sample}.txt"
    log:
        "log/bedgraphtobigwig/{sample}.log"
    shell:
        "wget -O {wildcards.sample}_{params.assembly}.chrom.sizes "
        "https://hgdownload.cse.ucsc.edu/goldenpath/{params.assembly}/bigZips/{params.assembly}.chrom.sizes && "
        "bedGraphToBigWig {input.bedgraph} "
        "{wildcards.sample}_{params.assembly}.chrom.sizes {output.bigwig} >& {log} && "
        "rm -rf {wildcards.sample}_{params.assembly}.chrom.sizes"

rule multiqc:
    container:
        "docker://multiqc/multiqc:v1.25"
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
