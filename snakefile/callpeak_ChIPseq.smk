
'''
Snakefile for ChIP-seq callpeak

Usage:
    snakemake -s callpeak_ChIPseq.smk --configfile <path to config.yaml> --cores <int> --use-singularity
'''

## DON'T CHANGE BELOW THIS LINE ##

def load_experiment(file):
    experiment_dict = {}
    with open(file, "r") as f:
        for line in f:
            sample, target, control = line.strip().split("\t")
            if sample == "sample":
                continue
            experiment_dict[sample] = {"target": target, "control": control}
    return experiment_dict
experiment_dict = load_experiment(config["general"]["experiment_table"])

workdir: config["general"]["workdir"]

rule all:
    input:
        multiqc = "multiqc_callpeak/multiqc_report.html",
        xls = expand("macs2/{sample}/{sample}_peaks.xls", sample = experiment_dict),
        bigwig = expand("macs2/{sample}/{sample}_treat_pileup.bw", sample = experiment_dict)

rule macs2:
    wildcard_constraints:
        sample = "|".join(experiment_dict)
    container:
        "docker://quay.io/biocontainers/macs2:2.2.9.1--py39hf95cd2a_0"
    input:
        target = lambda wildcards: experiment_dict[wildcards.sample]["target"],
        control = lambda wildcards: experiment_dict[wildcards.sample]["control"]
    output:
        outdir = directory("macs2/{sample}"),
        xls = "macs2/{sample}/{sample}_peaks.xls",
        bdg = "macs2/{sample}/{sample}_treat_pileup.bdg"
    params:
        genomesize = config["macs2"]["genomesize"],
        fileformat = config["macs2"]["fileformat"],
        qvalue = config["macs2"]["qvalue"],
        broad = config["macs2"]["broad"],
        broad_cutoff = config["macs2"]["broad_cutoff"]
    benchmark:
        "benchmark/macs2/{sample}.txt"
    log:
        "log/macs2/{sample}.log"
    shell:
        "if [ '{params.broad}' == 'True' ]; then "
        "macs2 callpeak "
        "-t {input.target} -c {input.control} "
        "-g {params.genomesize} -f {params.fileformat} "
        "-q {params.qvalue} "
        "--broad --broad-cutoff {params.broad_cutoff} "
        "-B -n {wildcards.sample} --outdir {output.outdir} >& {log}; "
        "else "
        "macs2 callpeak "
        "-t {input.target} -c {input.control} "
        "-g {params.genomesize} -f {params.fileformat} "
        "-q {params.qvalue} "
        "-B -n {wildcards.sample} --outdir {output.outdir} >& {log}; "
        "fi"

rule bedgraphtobigwig:
    wildcard_constraints:
        sample = "|".join(experiment_dict)
    container:
        "docker://quay.io/biocontainers/ucsc-bedgraphtobigwig:445--h954228d_0"
    input:
        bedgraph = "macs2/{sample}/{sample}_treat_pileup.bdg"
    output:
        bigwig = "macs2/{sample}/{sample}_treat_pileup.bw"
    params:
        assembly = config["bedgraphtobigwig"]["assembly"]
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
        macs2log = expand("macs2/{sample}/{sample}_peaks.xls", sample = experiment_dict)
    output:
        "multiqc_callpeak/multiqc_report.html"
    benchmark:
        "benchmark/multiqc_callpeak.txt"
    log:
        "log/multiqc_callpeak.log"
    shell:
        "rm -rf multiqc && "
        "mkdir -p multiqc/log && "
        "cp {input.macs2log} multiqc/log && "
        "multiqc -o multiqc_callpeak/ multiqc/log >& {log} && "
        "rm -rf multiqc/log"

## DON'T CHANGE ABOVE THIS LINE ##
