rule plotFingerprint:
    container:
        CONTAINERS["deeptools"]
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
