
'''
Snakefile for ChIP-seq preprocessing

Usage:
    snakemake -s preprocessing_ChIPseq.smk --configfile <path to config.yaml> --cores <int> --use-singularity
'''

## DON'T CHANGE BELOW THIS LINE ##

workdir: config["workdir"]
bowtie2_index = config["bowtie2_index"]
samples = config["samples"]

rule all:
    input:
        multiqc = "multiqc_preprocessing/multiqc_report.html",
        bam = expand("bowtie2/{sample}.sort.rmdup.bam", sample = samples),
        bai = expand("bowtie2/{sample}.sort.rmdup.bam.bai", sample = samples),
        bigwig = expand("bigwig/{sample}.bw", sample = samples),
        png = "plotFingerprint/fingerprint.png",
        qc = "plotFingerprint/fingerprint.qc.txt",
        tab = "plotFingerprint/fingerprint.tab"

rule qc:
    wildcard_constraints:
        sample = "|".join([re.escape(x) for x in samples])
    container:
        "docker://quay.io/biocontainers/fastp:0.23.4--hadf994f_2"
    input:
        R1 = "fastq/{sample}_1.fastq.gz",
        R2 = "fastq/{sample}_2.fastq.gz"
    output:
        R1 = "fastp/{sample}_1.fq.gz",
        R2 = "fastp/{sample}_2.fq.gz",
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
        "fastp -i {input.R1} -I {input.R2} "
        "-o {output.R1} -O {output.R2} -w {threads} "
        "-h {params.html} -j {output.json} >& {log}"

rule mapping_bowtie2:
    wildcard_constraints:
        sample = "|".join([re.escape(x) for x in samples])
    container:
        "docker://quay.io/biocontainers/bowtie2:2.5.2--py39h6fed5c7_0"
    input:
        R1 = "fastp/{sample}_1.fq.gz",
        R2 = "fastp/{sample}_2.fq.gz"
    output:
        sam = temp("bowtie2/{sample}.sam")
    params:
        bowtie2_args = config["bowtie2_args"]
    threads:
        1000
    benchmark:
        "benchmark/bowtie2_{sample}.txt"
    log:
        "log/bowtie2/{sample}.log"
    shell:
        "bowtie2 -p {threads} -x {bowtie2_index} "
        "{params.bowtie2_args} "
        "-1 {input.R1} -2 {input.R2} -S {output.sam} >& {log}"

rule sort:
    wildcard_constraints:
        sample = "|".join([re.escape(x) for x in samples])
    container:
        "docker://quay.io/biocontainers/samtools:1.18--h50ea8bc_1"
    input:
        "bowtie2/{sample}.sam"
    output:
        temp("bowtie2/{sample}.sort.bam")
    threads:
        8
    benchmark:
        "benchmark/samtools_{sample}.txt"
    log:
        "log/samtools_{sample}.log"
    shell:
        "samtools sort -@ {threads} -O bam -o {output} {input} >& {log} && "
        "samtools index {output}"

rule markdup:
    wildcard_constraints:
        sample = "|".join([re.escape(x) for x in samples])
    container:
        "docker://quay.io/biocontainers/picard:3.1.1--hdfd78af_0"
    input:
        "bowtie2/{sample}.sort.bam"
    output:
        "bowtie2/{sample}.sort.rmdup.bam"
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
        "--TMP_DIR bowtie2/tmp && "
        "rm -rf bowtie2/*.sort.bam.bai"

rule index:
    wildcard_constraints:
        sample = "|".join([re.escape(x) for x in samples])
    container:
        "docker://quay.io/biocontainers/samtools:1.18--h50ea8bc_1"
    input:
        "bowtie2/{sample}.sort.rmdup.bam"
    output:
        "bowtie2/{sample}.sort.rmdup.bam.bai"
    shell:
        "samtools index {input}"

rule plotFingerprint:
    wildcard_constraints:
        sample = "|".join([re.escape(x) for x in samples])
    container:
        "docker://quay.io/biocontainers/deeptools:3.5.4--pyhdfd78af_1"
    input:
        bam = expand("bowtie2/{sample}.sort.rmdup.bam", sample = samples),
        bai = expand("bowtie2/{sample}.sort.rmdup.bam.bai", sample = samples)
    output:
        png = "plotFingerprint/fingerprint.png",
        qc = "plotFingerprint/fingerprint.qc.txt",
        tab = "plotFingerprint/fingerprint.tab"
    params:
        labels = expand("{sample}", sample = samples)
    threads:
        8
    benchmark:
        "benchmark/plotFingerprint.txt"
    log:
        "log/plotFingerprint.log"
    shell:
        "plotFingerprint "
        "-b {input.bam} "
        "--labels {params.labels} "
        "--minMappingQuality 30 --skipZeros "
        "--numberOfSamples 50000 "
        "--binSize 10000 "
        "-T Fingerprints "
        "--plotFile {output.png} "
        "--outQualityMetrics {output.qc} "
        "--outRawCounts {output.tab} >& {log}"

rule bigwig:
    wildcard_constraints:
        sample = "|".join([re.escape(x) for x in samples])
    container:
        "docker://quay.io/biocontainers/deeptools:3.5.4--pyhdfd78af_1"
    input:
        bam = "bowtie2/{sample}.sort.rmdup.bam",
        bai = "bowtie2/{sample}.sort.rmdup.bam.bai"
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
        bowtie2log = expand("log/bowtie2/{sample}.log", sample = samples),
        picardlog = expand("log/picard_{sample}.log", sample = samples),
        plotFingerprintlog_qc = "plotFingerprint/fingerprint.qc.txt",
        plotFingerprintlog_tab = "plotFingerprint/fingerprint.tab"
    output:
        "multiqc_preprocessing/multiqc_report.html"
    benchmark:
        "benchmark/multiqc.txt"
    log:
        "log/multiqc.log"
    shell:
        """
        rm -rf multiqc_preprocessing && \
        mkdir -p multiqc_preprocessing/log && \
        cp {input.json} {input.bowtie2log} {input.picardlog} {input.plotFingerprintlog_qc} {input.plotFingerprintlog_tab} multiqc_preprocessing/log && \
        cat /usr/local/lib/python3.12/site-packages/multiqc/config_defaults.yaml | sed -e '$afastp:\\n  s_name_filenames: true' > multiqc_preprocessing/multiqc_config.yaml && \
        multiqc --config multiqc_preprocessing/multiqc_config.yaml -o multiqc_preprocessing/ multiqc_preprocessing/log >& {log} && \
        rm -rf multiqc_preprocessing/log
        """

## DON'T CHANGE ABOVE THIS LINE ##
