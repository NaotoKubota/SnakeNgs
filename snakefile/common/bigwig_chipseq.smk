rule bigwig:
    container:
        CONTAINERS["deeptools"]
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
