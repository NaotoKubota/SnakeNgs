rule bigwig:
    container:
        CONTAINERS["deeptools"]
    input:
        "star/{sample}/{sample}_Aligned.out.bam"
    output:
        "bigwig/{sample}.bw"
    threads:
        workflow.cores / 4
    benchmark:
        "benchmark/bamCoverage_{sample}.txt"
    log:
        "log/bamCoverage_{sample}.log"
    shell:
        "bamCoverage -b {input} -o {output} -p {threads} --binSize 1 >& {log}"
