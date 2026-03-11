rule markdup:
    container:
        CONTAINERS["picard"]
    input:
        "bowtie2/{sample}.sort.bam"
    output:
        "bowtie2/{sample}.sort.rmdup.bam"
    threads:
        8
    benchmark:
        "benchmark/picard_{sample}.txt"
    log:
        "log/picard/{sample}.log"
    shell:
        "picard MarkDuplicates "
        "-I {input} -O {output} "
        "-M {log} "
        "--REMOVE_DUPLICATES true "
        "--TMP_DIR bowtie2/tmp "
        ">& {log} && "
        "rm -rf bowtie2/*.sort.bam.bai"
